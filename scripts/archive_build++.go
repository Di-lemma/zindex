/*
# archive_build.py Re-Written in Go
### Rewritten By ck4sp3r 
# Something something Fuck python
*/


package main

import (
	"archive/tar"
	"bytes"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"io/fs"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
	"time"
)

type chunkRecord struct {
	ChunkID           int `json:"chunk_id"`
	CompressedOffset  int `json:"compressed_offset"`
	CompressedSize    int `json:"compressed_size"`
	UncompressedStart int `json:"uncompressed_start"`
	UncompressedEnd   int `json:"uncompressed_end"`
}

type chunkIndex struct {
	Chunks []chunkRecord `json:"chunks"`
}

type config struct {
	InputDir  string
	Output    string
	ChunkSize int
	Level     int
	Debug     bool
}

func parseArgs() (config, error) {
	var cfg config

	flag.IntVar(&cfg.ChunkSize, "chunk-size", 1024*1024, "Uncompressed bytes per zstd frame (default: 1048576)")
	flag.IntVar(&cfg.Level, "level", 9, "zstd compression level (default: 9)")
	flag.BoolVar(&cfg.Debug, "debug", false, "Print debug logs to stdout")
	flag.Parse()

	args := flag.Args()
	if len(args) != 2 {
		return cfg, errors.New("usage: archive_build [--chunk-size N] [--level N] [--debug] INPUT_DIR OUTPUT")
	}

	cfg.InputDir = args[0]
	cfg.Output = args[1]
	return cfg, nil
}

func normalizeOutputPath(output string) string {
	if strings.HasSuffix(output, ".tar.zst") {
		return output
	}
	return output + ".tar.zst"
}

func debugLog(enabled bool, message string) {
	if enabled {
		fmt.Printf("[debug] %s\n", message)
	}
}

func renderProgress(doneChunks, totalChunks, processedBytes, totalBytes int, ttyMode bool) {
	percent := 100
	if totalChunks > 0 {
		percent = (doneChunks * 100) / totalChunks
	}

	line := fmt.Sprintf(
		"[progress] chunks %d/%d (%3d%%) | bytes %d/%d",
		doneChunks,
		totalChunks,
		percent,
		processedBytes,
		totalBytes,
	)

	if ttyMode {
		fmt.Printf("\r%s", line)
		if doneChunks == totalChunks {
			fmt.Print("\n")
		}
		return
	}

	fmt.Println(line)
}

func collectTxtFiles(inputDir string) ([]string, error) {
	var files []string
	err := filepath.WalkDir(inputDir, func(path string, d fs.DirEntry, walkErr error) error {
		if walkErr != nil {
			return walkErr
		}
		if d.IsDir() {
			return nil
		}
		if strings.HasSuffix(d.Name(), ".txt") {
			files = append(files, path)
		}
		return nil
	})
	if err != nil {
		return nil, err
	}
	sort.Strings(files)
	if len(files) == 0 {
		return nil, fmt.Errorf("no .txt files found under %s", inputDir)
	}
	return files, nil
}

func buildTarFile(inputDir string, files []string, tempTarPath string) error {
	tf, err := os.Create(tempTarPath)
	if err != nil {
		return err
	}
	defer tf.Close()

	tw := tar.NewWriter(tf)
	defer tw.Close()

	for _, path := range files {
		info, err := os.Stat(path)
		if err != nil {
			return err
		}

		rel, err := filepath.Rel(inputDir, path)
		if err != nil {
			return err
		}
		arcname := filepath.ToSlash(rel)

		hdr := &tar.Header{
			Name:     arcname,
			Mode:     int64(info.Mode().Perm()),
			Size:     info.Size(),
			ModTime:  time.Unix(0, 0),
			Format:   tar.FormatPAX,
			Typeflag: tar.TypeReg,
			Uid:      0,
			Gid:      0,
			Uname:    "",
			Gname:    "",
		}
		if err := tw.WriteHeader(hdr); err != nil {
			return err
		}

		fh, err := os.Open(path)
		if err != nil {
			return err
		}
		_, copyErr := io.Copy(tw, fh)
		closeErr := fh.Close()
		if copyErr != nil {
			return copyErr
		}
		if closeErr != nil {
			return closeErr
		}
	}

	return nil
}

func compressChunk(chunk []byte, level int, zstdPath string) ([]byte, error) {
	cmd := exec.Command(zstdPath, "--compress", "--stdout", "--quiet", fmt.Sprintf("-%d", level))
	cmd.Stdin = bytes.NewReader(chunk)

	var stdout bytes.Buffer
	var stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr

	if err := cmd.Run(); err != nil {
		msg := strings.TrimSpace(stderr.String())
		if msg == "" {
			msg = err.Error()
		}
		return nil, fmt.Errorf("zstd failed: %s", msg)
	}

	return stdout.Bytes(), nil
}

func validateZindexChunkIndex(chunks []chunkRecord, archiveSize int) error {
	if len(chunks) == 0 {
		return errors.New("internal error: chunk index is empty")
	}

	expectedChunkID := 0
	expectedCompressedOffset := 0
	expectedUncompressedStart := 0

	for _, chunk := range chunks {
		if chunk.ChunkID != expectedChunkID {
			return errors.New("internal error: non-sequential chunk_id values")
		}
		if chunk.CompressedOffset != expectedCompressedOffset {
			return errors.New("internal error: non-contiguous compressed offsets")
		}
		if chunk.UncompressedStart != expectedUncompressedStart {
			return errors.New("internal error: non-contiguous uncompressed ranges")
		}
		if chunk.CompressedSize <= 0 {
			return errors.New("internal error: chunk compressed_size must be positive")
		}
		if chunk.UncompressedEnd <= chunk.UncompressedStart {
			return errors.New("internal error: invalid uncompressed range in chunk index")
		}

		expectedChunkID++
		expectedCompressedOffset += chunk.CompressedSize
		expectedUncompressedStart = chunk.UncompressedEnd
	}

	if expectedCompressedOffset != archiveSize {
		return errors.New("internal error: archive size and index offsets do not match")
	}

	return nil
}

func writeArchiveAndIndex(tempTarPath, outputPath string, chunkSize, level int, debug bool) (string, error) {
	indexPath := outputPath + ".json"
	if err := os.MkdirAll(filepath.Dir(outputPath), 0o755); err != nil {
		return "", err
	}

	tarStat, err := os.Stat(tempTarPath)
	if err != nil {
		return "", err
	}
	if tarStat.Size() <= 0 {
		return "", errors.New("internal error: empty tar payload")
	}

	totalBytes := int(tarStat.Size())
	totalChunks := (totalBytes + chunkSize - 1) / chunkSize
	ttyMode := isTerminal(os.Stdout)
	nonTTYInterval := totalChunks / 20
	if nonTTYInterval < 1 {
		nonTTYInterval = 1
	}

	zstdPath, err := exec.LookPath("zstd")
	if err != nil {
		return "", errors.New("zstd CLI not found in THE FUCKING PATH")
	}

	debugLog(debug, fmt.Sprintf("writing archive to %s", outputPath))
	debugLog(debug, fmt.Sprintf("writing index to %s", indexPath))
	debugLog(debug, fmt.Sprintf("tar size=%d bytes; chunks=%d; level=%d", totalBytes, totalChunks, level))

	in, err := os.Open(tempTarPath)
	if err != nil {
		return "", err
	}
	defer in.Close()

	out, err := os.Create(outputPath)
	if err != nil {
		return "", err
	}
	defer out.Close()

	buf := make([]byte, chunkSize)
	chunks := make([]chunkRecord, 0, totalChunks)
	compressedOffset := 0
	processed := 0

	for chunkNo := 1; chunkNo <= totalChunks; chunkNo++ {
		n, readErr := io.ReadFull(in, buf)
		if readErr != nil && !errors.Is(readErr, io.ErrUnexpectedEOF) {
			return "", readErr
		}
		if n == 0 {
			break
		}

		chunkData := buf[:n]
		compressed, err := compressChunk(chunkData, level, zstdPath)
		if err != nil {
			return "", err
		}
		if _, err := out.Write(compressed); err != nil {
			return "", err
		}

		uncompressedStart := processed
		processed += n

		chunks = append(chunks, chunkRecord{
			ChunkID:           chunkNo - 1,
			CompressedOffset:  compressedOffset,
			CompressedSize:    len(compressed),
			UncompressedStart: uncompressedStart,
			UncompressedEnd:   processed,
		})
		compressedOffset += len(compressed)

		shouldRender := ttyMode || chunkNo == totalChunks || (chunkNo%nonTTYInterval == 0)
		if shouldRender {
			renderProgress(chunkNo, totalChunks, processed, totalBytes, ttyMode)
		}

		if errors.Is(readErr, io.ErrUnexpectedEOF) {
			break
		}
	}

	if err := validateZindexChunkIndex(chunks, compressedOffset); err != nil {
		return "", err
	}
	debugLog(debug, "chunk index validation passed")

	idxPayload := chunkIndex{Chunks: chunks}
	idxBytes, err := json.MarshalIndent(idxPayload, "", "  ")
	if err != nil {
		return "", err
	}
	idxBytes = append(idxBytes, '\n')
	if err := os.WriteFile(indexPath, idxBytes, 0o644); err != nil {
		return "", err
	}

	return indexPath, nil
}

func isTerminal(f *os.File) bool {
	info, err := f.Stat()
	if err != nil {
		return false
	}
	return (info.Mode() & os.ModeCharDevice) != 0
}

/*
# Suddenly 

  # # #
#       #			  #######
#       # # # # # # # # # # # #  ###########
#       #			######------
  # # #				######------
  # # #				########### 
#       # # # # # # # # # # # #   #######
#       # 			 
#       #
  # # #

# Penis
*/

func main() {
	cfg, err := parseArgs()
	if err != nil {
		fmt.Fprintln(os.Stderr, err.Error())
		os.Exit(1)
	}

	if cfg.ChunkSize <= 0 {
		fmt.Fprintln(os.Stderr, "--chunk-size must be positive")
		os.Exit(1)
	}

	inputDir, err := filepath.Abs(cfg.InputDir)
	if err != nil {
		fmt.Fprintln(os.Stderr, err.Error())
		os.Exit(1)
	}
	if st, err := os.Stat(inputDir); err != nil || !st.IsDir() {
		fmt.Fprintf(os.Stderr, "input directory does not exist: %s\n", inputDir)
		os.Exit(1)
	}

	outputAbs, err := filepath.Abs(cfg.Output)
	if err != nil {
		fmt.Fprintln(os.Stderr, err.Error())
		os.Exit(1)
	}
	outputAbs = normalizeOutputPath(outputAbs)

	files, err := collectTxtFiles(inputDir)
	if err != nil {
		fmt.Fprintln(os.Stderr, err.Error())
		os.Exit(1)
	}
	debugLog(cfg.Debug, fmt.Sprintf("found %d .txt files under %s", len(files), inputDir))

	tempFile, err := os.CreateTemp("", "archive_build_*.tar")
	if err != nil {
		fmt.Fprintln(os.Stderr, err.Error())
		os.Exit(1)
	}
	tempTarPath := tempFile.Name()
	tempFile.Close()
	defer os.Remove(tempTarPath)

	if err := buildTarFile(inputDir, files, tempTarPath); err != nil {
		fmt.Fprintln(os.Stderr, err.Error())
		os.Exit(1)
	}

	tarInfo, err := os.Stat(tempTarPath)
	if err != nil {
		fmt.Fprintln(os.Stderr, err.Error())
		os.Exit(1)
	}
	debugLog(cfg.Debug, fmt.Sprintf("built tar payload with %d bytes", tarInfo.Size()))

	indexPath, err := writeArchiveAndIndex(tempTarPath, outputAbs, cfg.ChunkSize, cfg.Level, cfg.Debug)
	if err != nil {
		fmt.Fprintln(os.Stderr, err.Error())
		os.Exit(1)
	}

	fmt.Println(outputAbs)
	fmt.Println(indexPath)
}
