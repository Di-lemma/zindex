/*
    zindex.cpp — tiny searcher over chunked .tar.zst archives in 1337 lines.

    Archives are stored as concatenated independent zstd frames. Each
    archive has a sibling JSON file named <archive>.json with per-frame
    offsets and byte ranges under a top-level "chunks" array.
*/

#define _FILE_OFFSET_BITS 64

#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <fnmatch.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <time.h>
#include <unistd.h>
#include <limits.h>

extern "C"
{
#include <zstd.h>
}

#include <functional>
#include <string>
#include <vector>

#if defined(__aarch64__) || defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#define HAVE_NEON 1
#else
#define HAVE_NEON 0
#endif

#define MAX_THREADS 32
#define MAX_QUEUE_SIZE 8192ULL
#define MAX_LINE_PREVIEW 2048
#define MAX_FILENAME_LEN 512
#define ASCII_SET_SIZE 256
#define RESULT_BATCH_SIZE 64
#define MEMCHR_HYBRID_MAX 32

#define ZSTD_IN_CHUNK (1u << 20) /* 1 MiB compressed feed */
#define ZSTD_OUT_WIN (1u << 20)  /* 1 MiB decompressed window */

typedef struct
{
  const char *pat;
  size_t len;
  int bad[ASCII_SET_SIZE];
} PatCtx;

static PatCtx g_pat;

static void build_badchar_table(const char *pattern, size_t patLen, int bad[ASCII_SET_SIZE])
{
  if (!pattern || patLen == 0)
    return;
  for (int i = 0; i < ASCII_SET_SIZE; i++)
    bad[i] = (int)patLen;
  for (size_t i = 0; i + 1 < patLen; i++)
    bad[(unsigned char)pattern[i]] = (int)(patLen - 1 - i);
}

static inline bool bytes_equal(const char *a, const char *b, size_t len)
{
#if HAVE_NEON
  while (len >= 16)
  {
    uint8x16_t va = vld1q_u8((const uint8_t *)a);
    uint8x16_t vb = vld1q_u8((const uint8_t *)b);
    uint8x16_t cmp = vceqq_u8(va, vb);
    uint64x2_t lanes = vreinterpretq_u64_u8(cmp);
    if ((vgetq_lane_u64(lanes, 0) != UINT64_MAX) || (vgetq_lane_u64(lanes, 1) != UINT64_MAX))
      return false;
    a += 16;
    b += 16;
    len -= 16;
  }
  for (size_t i = 0; i < len; ++i)
    if ((unsigned char)a[i] != (unsigned char)b[i])
      return false;
  return true;
#else
  return memcmp(a, b, len) == 0;
#endif
}

static inline ssize_t BMH_find(const char *hay, size_t hayLen,
                               const char *pat, size_t patLen,
                               const int bad[256])
{
  if (!hay || !pat || patLen == 0 || hayLen < patLen)
    return -1;
  if (patLen <= MEMCHR_HYBRID_MAX)
  {
    const char *base = hay;
    const char *limit = hay + (hayLen - patLen + 1);
    const unsigned char first = (unsigned char)pat[0];
    while (hay < limit)
    {
      const void *found = memchr(hay, first, (size_t)(limit - hay));
      if (!found)
        return -1;
      const char *cand = (const char *)found;
      if (patLen == 1 || bytes_equal(cand + 1, pat + 1, patLen - 1))
        return (ssize_t)(cand - base);
      hay = cand + 1;
    }
    return -1;
  }
  size_t i = 0;
  while (i <= hayLen - patLen)
  {
    bool match = false;
    if (patLen <= 16)
    {
      match = bytes_equal(hay + i, pat, patLen);
    }
    else
    {
      if (bytes_equal(hay + i + patLen - 16, pat + patLen - 16, 16))
        match = bytes_equal(hay + i, pat, patLen);
    }
    if (match)
      return (ssize_t)i;
    unsigned char last = (unsigned char)hay[i + patLen - 1];
    int shift = bad[last];
    i += (shift > 0) ? shift : 1;
  }
  return -1;
}

typedef struct
{
  char filename[MAX_FILENAME_LEN];
  char internalFilename[MAX_FILENAME_LEN];
  int64_t offset;
  char preview[MAX_LINE_PREVIEW];
} SearchResult;

typedef struct PrintNode
{
  SearchResult item;
  struct PrintNode *next;
} PrintNode;

typedef struct
{
  PrintNode *head, *tail;
  size_t size, max_size;
  bool done;
  pthread_mutex_t lock;
  pthread_cond_t not_empty;
  pthread_cond_t not_full;
} PrintQueue;

typedef struct
{
  SearchResult items[RESULT_BATCH_SIZE];
  size_t count;
} ResultBatch;

static void printqueue_init(PrintQueue *q)
{
  memset(q, 0, sizeof(*q));
  q->max_size = MAX_QUEUE_SIZE;
  pthread_mutex_init(&q->lock, NULL);
  pthread_cond_init(&q->not_empty, NULL);
  pthread_cond_init(&q->not_full, NULL);
}
static void printqueue_destroy(PrintQueue *q)
{
  pthread_mutex_lock(&q->lock);
  for (PrintNode *cur = q->head; cur;)
  {
    PrintNode *n = cur->next;
    free(cur);
    cur = n;
  }
  q->head = q->tail = NULL;
  q->size = 0;
  pthread_mutex_unlock(&q->lock);
  pthread_mutex_destroy(&q->lock);
  pthread_cond_destroy(&q->not_empty);
  pthread_cond_destroy(&q->not_full);
}
static bool printqueue_push_batch(PrintQueue *q, const SearchResult *items, size_t count)
{
  if (!q || !items || count == 0)
    return true;

  PrintNode *head = NULL;
  PrintNode *tail = NULL;
  for (size_t i = 0; i < count; ++i)
  {
    PrintNode *node = (PrintNode *)malloc(sizeof(*node));
    if (!node)
    {
      for (PrintNode *cur = head; cur;)
      {
        PrintNode *next = cur->next;
        free(cur);
        cur = next;
      }
      return false;
    }
    memcpy(&node->item, &items[i], sizeof(items[i]));
    node->item.filename[sizeof(node->item.filename) - 1] = '\0';
    node->item.preview[sizeof(node->item.preview) - 1] = '\0';
    node->next = NULL;
    if (tail)
      tail->next = node;
    else
      head = node;
    tail = node;
  }

  pthread_mutex_lock(&q->lock);
  while (q->size + count > q->max_size && !q->done)
    pthread_cond_wait(&q->not_full, &q->lock);
  if (q->done)
  {
    pthread_mutex_unlock(&q->lock);
    for (PrintNode *cur = head; cur;)
    {
      PrintNode *next = cur->next;
      free(cur);
      cur = next;
    }
    return false;
  }
  if (q->tail)
    q->tail->next = head;
  else
    q->head = head;
  q->tail = tail;
  q->size += count;
  pthread_cond_signal(&q->not_empty);
  pthread_mutex_unlock(&q->lock);
  return true;
}
static bool printqueue_pop(PrintQueue *q, SearchResult *out)
{
  pthread_mutex_lock(&q->lock);
  while (q->size == 0 && !q->done)
    pthread_cond_wait(&q->not_empty, &q->lock);
  if (q->size == 0 && q->done)
  {
    pthread_mutex_unlock(&q->lock);
    return false;
  }
  PrintNode *node = q->head;
  q->head = node->next;
  if (!q->head)
    q->tail = NULL;
  q->size--;
  memcpy(out, &node->item, sizeof(*out));
  free(node);
  pthread_cond_signal(&q->not_full);
  pthread_mutex_unlock(&q->lock);
  return true;
}
static void printqueue_mark_done(PrintQueue *q)
{
  pthread_mutex_lock(&q->lock);
  q->done = true;
  pthread_cond_broadcast(&q->not_empty);
  pthread_cond_broadcast(&q->not_full);
  pthread_mutex_unlock(&q->lock);
}

static bool result_batch_flush(PrintQueue *q, ResultBatch *batch)
{
  if (!batch || batch->count == 0)
    return true;
  bool ok = printqueue_push_batch(q, batch->items, batch->count);
  batch->count = 0;
  return ok;
}

static bool result_batch_add(PrintQueue *q, ResultBatch *batch, const SearchResult *res)
{
  if (!batch || !res)
    return false;
  batch->items[batch->count++] = *res;
  if (batch->count < RESULT_BATCH_SIZE)
    return true;
  return result_batch_flush(q, batch);
}

static void *printer_thread(void *arg)
{
  PrintQueue *q = (PrintQueue *)arg;
  char line[MAX_LINE_PREVIEW + MAX_FILENAME_LEN + 128];
  SearchResult sr;
  setvbuf(stdout, NULL, _IOFBF, 1 << 20);
  while (printqueue_pop(q, &sr))
  {
    int n = snprintf(line, sizeof(line), "%s\n", sr.preview);
    if (n > 0)
      (void)fwrite(line, 1, (size_t)n, stdout);
  }
  return NULL;
}

static void make_preview(const char *text, size_t textLen,
                         size_t matchPos, size_t /*patLen*/,
                         const char * /*pattern*/, SearchResult *outRes)
{
  if (!text || !outRes || matchPos >= textLen)
  {
    if (outRes)
      outRes->preview[0] = '\0';
    return;
  }
  size_t start = matchPos, end = matchPos;
  while (start > 0 && text[start - 1] != '\n')
    --start;
  while (end < textLen && text[end] != '\n' && text[end] != '\0')
    ++end;
  size_t len = end > start ? (end - start) : 0;
  if (len >= MAX_LINE_PREVIEW - 1)
    len = MAX_LINE_PREVIEW - 1;
  memcpy(outRes->preview, text + start, len);
  outRes->preview[len] = '\0';
  outRes->internalFilename[0] = '\0';
  const char *exts[] = {".txt", ".csv", ".json", ".md", ".log", NULL};
  for (int i = 0; exts[i]; ++i)
  {
    char *pos = strstr(outRes->preview, exts[i]);
    if (!pos)
      continue;
    char *begin = pos;
    while (begin > outRes->preview)
    {
      char c = *(begin - 1);
      if (c == ' ' || c == '/' || c == '\\' || c == ':')
        break;
      --begin;
    }
    size_t fLen = (size_t)(pos - begin + strlen(exts[i]));
    if (fLen && fLen < MAX_FILENAME_LEN)
    {
      strncpy(outRes->internalFilename, begin, fLen);
      outRes->internalFilename[fLen] = '\0';
    }
    break;
  }
}

typedef struct
{
  int64_t chunk_id;
  int64_t compressed_offset;
  int64_t compressed_size;
  int64_t uncompressed_start;
  int64_t uncompressed_end;
} ChunkInfo;

typedef struct
{
  ChunkInfo *chunks;
  size_t count;
} IndexInfo;

typedef struct
{
  int fd;
  const unsigned char *data;
  size_t size;
  bool mapped;
} ArchiveView;

static IndexInfo *parse_index_json(const char *idx_filename)
{
  struct JsonCursor
  {
    const char *p;
    const char *end;
  };

  auto skip_ws = [](JsonCursor &cur) {
    while (cur.p < cur.end && isspace((unsigned char)*cur.p))
      ++cur.p;
  };

  auto parse_string = [&](JsonCursor &cur, std::string &out) -> bool {
    skip_ws(cur);
    if (cur.p >= cur.end || *cur.p != '"')
      return false;
    ++cur.p;
    out.clear();
    while (cur.p < cur.end)
    {
      char ch = *cur.p++;
      if (ch == '"')
        return true;
      if (ch == '\\')
      {
        if (cur.p >= cur.end)
          return false;
        char esc = *cur.p++;
        switch (esc)
        {
        case '"':
        case '\\':
        case '/':
          out.push_back(esc);
          break;
        case 'b':
          out.push_back('\b');
          break;
        case 'f':
          out.push_back('\f');
          break;
        case 'n':
          out.push_back('\n');
          break;
        case 'r':
          out.push_back('\r');
          break;
        case 't':
          out.push_back('\t');
          break;
        default:
          return false;
        }
        continue;
      }
      out.push_back(ch);
    }
    return false;
  };

  auto expect_char = [&](JsonCursor &cur, char want) -> bool {
    skip_ws(cur);
    if (cur.p >= cur.end || *cur.p != want)
      return false;
    ++cur.p;
    return true;
  };

  auto parse_int64 = [&](JsonCursor &cur, int64_t &out) -> bool {
    skip_ws(cur);
    if (cur.p >= cur.end)
      return false;
    char *next = NULL;
    errno = 0;
    long long value = strtoll(cur.p, &next, 10);
    if (next == cur.p || errno != 0)
      return false;
    cur.p = next;
    out = (int64_t)value;
    return true;
  };

  std::function<bool(JsonCursor &)> skip_value;
  skip_value = [&](JsonCursor &cur) -> bool {
    skip_ws(cur);
    if (cur.p >= cur.end)
      return false;
    if (*cur.p == '"')
    {
      std::string tmp;
      return parse_string(cur, tmp);
    }
    if (*cur.p == '{')
    {
      ++cur.p;
      skip_ws(cur);
      if (cur.p < cur.end && *cur.p == '}')
      {
        ++cur.p;
        return true;
      }
      while (cur.p < cur.end)
      {
        std::string key;
        if (!parse_string(cur, key) || !expect_char(cur, ':') || !skip_value(cur))
          return false;
        skip_ws(cur);
        if (cur.p < cur.end && *cur.p == ',')
        {
          ++cur.p;
          continue;
        }
        return expect_char(cur, '}');
      }
      return false;
    }
    if (*cur.p == '[')
    {
      ++cur.p;
      skip_ws(cur);
      if (cur.p < cur.end && *cur.p == ']')
      {
        ++cur.p;
        return true;
      }
      while (cur.p < cur.end)
      {
        if (!skip_value(cur))
          return false;
        skip_ws(cur);
        if (cur.p < cur.end && *cur.p == ',')
        {
          ++cur.p;
          continue;
        }
        return expect_char(cur, ']');
      }
      return false;
    }
    if (cur.end - cur.p >= 4 && !strncmp(cur.p, "true", 4))
    {
      cur.p += 4;
      return true;
    }
    if (cur.end - cur.p >= 5 && !strncmp(cur.p, "false", 5))
    {
      cur.p += 5;
      return true;
    }
    if (cur.end - cur.p >= 4 && !strncmp(cur.p, "null", 4))
    {
      cur.p += 4;
      return true;
    }
    int64_t ignored = 0;
    return parse_int64(cur, ignored);
  };

  auto parse_chunk = [&](JsonCursor &cur, ChunkInfo &chunk) -> bool {
    chunk = {};
    if (!expect_char(cur, '{'))
      return false;
    skip_ws(cur);
    if (cur.p < cur.end && *cur.p == '}')
    {
      ++cur.p;
      return true;
    }
    while (cur.p < cur.end)
    {
      std::string key;
      if (!parse_string(cur, key) || !expect_char(cur, ':'))
        return false;
      if (key == "chunk_id")
      {
        if (!parse_int64(cur, chunk.chunk_id))
          return false;
      }
      else if (key == "compressed_offset")
      {
        if (!parse_int64(cur, chunk.compressed_offset))
          return false;
      }
      else if (key == "compressed_size")
      {
        if (!parse_int64(cur, chunk.compressed_size))
          return false;
      }
      else if (key == "uncompressed_start")
      {
        if (!parse_int64(cur, chunk.uncompressed_start))
          return false;
      }
      else if (key == "uncompressed_end")
      {
        if (!parse_int64(cur, chunk.uncompressed_end))
          return false;
      }
      else if (!skip_value(cur))
      {
        return false;
      }

      skip_ws(cur);
      if (cur.p < cur.end && *cur.p == ',')
      {
        ++cur.p;
        continue;
      }
      return expect_char(cur, '}');
    }
    return false;
  };

  FILE *fp = fopen(idx_filename, "rb");
  if (!fp)
  {
    fprintf(stderr, "open failed: %s\n", idx_filename);
    return NULL;
  }
  if (fseeko(fp, 0, SEEK_END) != 0)
  {
    fclose(fp);
    return NULL;
  }
  off_t sz = ftello(fp);
  if (sz < 0 || fseeko(fp, 0, SEEK_SET) != 0)
  {
    fclose(fp);
    return NULL;
  }

  std::string json;
  json.resize((size_t)sz);
  if (sz > 0 && fread(&json[0], 1, (size_t)sz, fp) != (size_t)sz)
  {
    fclose(fp);
    fprintf(stderr, "short read: %s\n", idx_filename);
    return NULL;
  }
  fclose(fp);

  JsonCursor cur = {json.data(), json.data() + json.size()};
  if (!expect_char(cur, '{'))
  {
    fprintf(stderr, "Index parse error: invalid root object in %s\n", idx_filename);
    return NULL;
  }

  std::vector<ChunkInfo> chunks;
  while (cur.p < cur.end)
  {
    skip_ws(cur);
    if (cur.p < cur.end && *cur.p == '}')
    {
      ++cur.p;
      break;
    }

    std::string key;
    if (!parse_string(cur, key) || !expect_char(cur, ':'))
    {
      fprintf(stderr, "Index parse error: malformed object in %s\n", idx_filename);
      return NULL;
    }

    if (key == "chunks")
    {
      if (!expect_char(cur, '['))
      {
        fprintf(stderr, "Index parse error: malformed chunks array in %s\n", idx_filename);
        return NULL;
      }
      skip_ws(cur);
      if (cur.p < cur.end && *cur.p == ']')
      {
        ++cur.p;
      }
      else
      {
        while (cur.p < cur.end)
        {
          ChunkInfo chunk = {};
          if (!parse_chunk(cur, chunk))
          {
            fprintf(stderr, "Index parse error: malformed chunk in %s\n", idx_filename);
            return NULL;
          }
          chunks.push_back(chunk);
          skip_ws(cur);
          if (cur.p < cur.end && *cur.p == ',')
          {
            ++cur.p;
            continue;
          }
          if (!expect_char(cur, ']'))
          {
            fprintf(stderr, "Index parse error: malformed chunks array in %s\n", idx_filename);
            return NULL;
          }
          break;
        }
      }
    }
    else if (!skip_value(cur))
    {
      fprintf(stderr, "Index parse error: failed to skip key '%s' in %s\n", key.c_str(), idx_filename);
      return NULL;
    }

    skip_ws(cur);
    if (cur.p < cur.end && *cur.p == ',')
    {
      ++cur.p;
      continue;
    }
    if (cur.p < cur.end && *cur.p == '}')
    {
      ++cur.p;
      break;
    }
  }

  if (chunks.empty())
    return NULL;

  IndexInfo *info = (IndexInfo *)calloc(1, sizeof(*info));
  if (!info)
    return NULL;
  info->count = chunks.size();
  info->chunks = (ChunkInfo *)calloc(info->count, sizeof(ChunkInfo));
  if (!info->chunks)
  {
    free(info);
    return NULL;
  }
  memcpy(info->chunks, chunks.data(), info->count * sizeof(ChunkInfo));
  return info;
}
static void free_index_info(IndexInfo *info)
{
  if (!info)
    return;
  if (info->chunks)
    free(info->chunks);
  free(info);
}

static bool map_archive_file(const char *path, ArchiveView *view)
{
  if (!view)
    return false;
  memset(view, 0, sizeof(*view));
  view->fd = -1;

  int fd = open(path, O_RDONLY);
  if (fd < 0)
  {
    fprintf(stderr, "open failed: %s\n", path);
    return false;
  }

  struct stat st;
  if (fstat(fd, &st) != 0)
  {
    fprintf(stderr, "fstat failed on %s: %s\n", path, strerror(errno));
    close(fd);
    return false;
  }

  if (!S_ISREG(st.st_mode) || st.st_size <= 0)
  {
    view->fd = fd;
    return true;
  }

  if ((uint64_t)st.st_size > (uint64_t)SIZE_MAX)
  {
    view->fd = fd;
    return true;
  }

  void *mapping = mmap(NULL, (size_t)st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  if (mapping == MAP_FAILED)
  {
    view->fd = fd;
    return true;
  }

#if defined(MADV_SEQUENTIAL)
  (void)madvise(mapping, (size_t)st.st_size, MADV_SEQUENTIAL);
#endif
#if defined(MADV_WILLNEED)
  (void)madvise(mapping, (size_t)st.st_size, MADV_WILLNEED);
#endif

  view->fd = fd;
  view->data = (const unsigned char *)mapping;
  view->size = (size_t)st.st_size;
  view->mapped = true;
  return true;
}

static void close_archive_view(ArchiveView *view)
{
  if (!view)
    return;
  if (view->mapped && view->data && view->size > 0)
    munmap((void *)view->data, view->size);
  if (view->fd >= 0)
    close(view->fd);
  memset(view, 0, sizeof(*view));
  view->fd = -1;
}

static void search_in_buffer(PrintQueue *pq,
                             ResultBatch *batch,
                             const char *filename,
                             const char *text,
                             size_t textLen,
                             const PatCtx *ctx,
                             int64_t baseOff)
{
  if (!pq || !text || !ctx || !ctx->pat || ctx->len == 0 || textLen < ctx->len)
    return;
  const int *bad = ctx->bad;
  const size_t plen = ctx->len;
  const char *pat = ctx->pat;
  size_t offset = 0;
  size_t hits = 0;
  const size_t MAX_HITS = 1000;
  while (offset + plen <= textLen && hits < MAX_HITS)
  {
    ssize_t rel = BMH_find(text + offset, textLen - offset, pat, plen, bad);
    if (rel < 0)
      break;
    size_t pos = offset + (size_t)rel;
    SearchResult sr = {};
    strncpy(sr.filename, filename, sizeof(sr.filename) - 1);
    sr.offset = baseOff + (int64_t)pos;
    make_preview(text, textLen, pos, plen, pat, &sr);
    if (sr.preview[0] != '\0')
    {
      if (!result_batch_add(pq, batch, &sr))
      {
        fprintf(stderr, "[warn] print queue saturated; dropping rest of this buffer\n");
        break;
      }
      ++hits;
    }
    offset = pos + 1;
  }
  if (hits == MAX_HITS)
    fprintf(stderr, "[warn] hit cap reached for buffer in %s at offset %lld; truncating additional matches\n",
            filename ? filename : "(unknown)", (long long)baseOff);
}

static void stream_decompress_and_search_mapped(PrintQueue *pq,
                                         ResultBatch *batch,
                                         const unsigned char *compressed,
                                         size_t compressed_size,
                                         const ChunkInfo *c,
                                         const PatCtx *ctx,
                                         const char *archive,
                                         char *buf,
                                         ZSTD_DCtx *dctx)
{
  const size_t OVER = (ctx->len > 1) ? (ctx->len - 1) : 0;

  memset(buf, 0, OVER);
  char *dst = buf + OVER;

  int64_t absOff = c->uncompressed_start;
  bool first = true;

  ZSTD_inBuffer in = {compressed, compressed_size, 0};
  while (in.pos < in.size)
  {
    ZSTD_outBuffer out = {dst, ZSTD_OUT_WIN, 0};
    size_t ret = ZSTD_decompressStream(dctx, &out, &in);
    if (ZSTD_isError(ret))
    {
      fprintf(stderr, "zstd error on %s chunk %lld: %s\n", archive, (long long)c->chunk_id, ZSTD_getErrorName(ret));
      /* Reset session if available, then continue with next input slice */
#if defined(ZSTD_VERSION_MAJOR)
      (void)ZSTD_DCtx_reset(dctx, ZSTD_reset_session_only);
#endif
      break;
    }

    if (out.pos > 0)
    {
      char *search_ptr = first ? dst : buf;
      size_t search_len = first ? out.pos : OVER + out.pos;
      int64_t base = first ? absOff : absOff - (int64_t)OVER;

      search_in_buffer(pq, batch, archive, search_ptr, search_len, ctx, base);

      absOff += (int64_t)out.pos;
      first = false;

      if (OVER > 0)
      {
        if (out.pos >= OVER)
        {
          memcpy(buf, dst + out.pos - OVER, OVER);
        }
        else
        {
          memmove(buf, buf + out.pos, OVER - out.pos);
          memcpy(buf + (OVER - out.pos), dst, out.pos);
        }
      }
    }
    /* ret==0 => end of frame; loop continues if more compressed input remains */
  }
}

static void stream_decompress_and_search_pread(PrintQueue *pq,
                                               ResultBatch *batch,
                                               int fd,
                                               const ChunkInfo *c,
                                               const PatCtx *ctx,
                                               const char *archive,
                                               char *buf,
                                               char *inbuf,
                                               ZSTD_DCtx *dctx)
{
  const size_t OVER = (ctx->len > 1) ? (ctx->len - 1) : 0;

  memset(buf, 0, OVER);
  char *dst = buf + OVER;

  int64_t absOff = c->uncompressed_start;
  off_t file_off = (off_t)c->compressed_offset;
  size_t remaining = (size_t)c->compressed_size;
  bool first = true;

  while (remaining > 0)
  {
    size_t to_read = remaining < ZSTD_IN_CHUNK ? remaining : ZSTD_IN_CHUNK;
    ssize_t r = pread(fd, inbuf, to_read, file_off);
    if (r <= 0)
    {
      fprintf(stderr, "pread failed on %s (chunk %lld): %s\n", archive, (long long)c->chunk_id, strerror(errno));
      break;
    }
    file_off += (off_t)r;
    remaining -= (size_t)r;

    ZSTD_inBuffer in = {inbuf, (size_t)r, 0};

    while (in.pos < in.size)
    {
      ZSTD_outBuffer out = {dst, ZSTD_OUT_WIN, 0};
      size_t ret = ZSTD_decompressStream(dctx, &out, &in);
      if (ZSTD_isError(ret))
      {
        fprintf(stderr, "zstd error on %s chunk %lld: %s\n", archive, (long long)c->chunk_id, ZSTD_getErrorName(ret));
#if defined(ZSTD_VERSION_MAJOR)
        (void)ZSTD_DCtx_reset(dctx, ZSTD_reset_session_only);
#endif
        break;
      }

      if (out.pos > 0)
      {
        char *search_ptr = first ? dst : buf;
        size_t search_len = first ? out.pos : OVER + out.pos;
        int64_t base = first ? absOff : absOff - (int64_t)OVER;

        search_in_buffer(pq, batch, archive, search_ptr, search_len, ctx, base);

        absOff += (int64_t)out.pos;
        first = false;

        if (OVER > 0)
        {
          if (out.pos >= OVER)
          {
            memcpy(buf, dst + out.pos - OVER, OVER);
          }
          else
          {
            memmove(buf, buf + out.pos, OVER - out.pos);
            memcpy(buf + (OVER - out.pos), dst, out.pos);
          }
        }
      }
    }
  }
}

static void search_indexed_file(PrintQueue *pq,
                                ResultBatch *batch,
                                const char *zstFile,
                                const char *idxFile,
                                const PatCtx *ctx,
                                char *buf,
                                char *inbuf,
                                ZSTD_DCtx *dctx)
{
  IndexInfo *info = parse_index_json(idxFile);
  if (!info)
    return;

  ArchiveView archive = {};
  if (!map_archive_file(zstFile, &archive))
  {
    free_index_info(info);
    return;
  }

  for (size_t i = 0; i < info->count; ++i)
  {
    const ChunkInfo *C = &info->chunks[i];
    if (archive.mapped)
    {
      if (C->compressed_offset < 0 || C->compressed_size < 0)
        continue;
      uint64_t start = (uint64_t)C->compressed_offset;
      uint64_t size = (uint64_t)C->compressed_size;
      if (start > archive.size || size > archive.size - start)
      {
        fprintf(stderr, "chunk %lld out of bounds in %s\n", (long long)C->chunk_id, zstFile);
        continue;
      }
      stream_decompress_and_search_mapped(
          pq,
          batch,
          archive.data + start,
          (size_t)size,
          C,
          ctx,
          zstFile,
          buf,
          dctx);
    }
    else
    {
      stream_decompress_and_search_pread(pq, batch, archive.fd, C, ctx, zstFile, buf, inbuf, dctx);
    }
    
    // Reset the decompression context between chunks for proper stream processing
#if defined(ZSTD_VERSION_MAJOR)
    ZSTD_DCtx_reset(dctx, ZSTD_reset_session_only);
#endif
  }

  (void)result_batch_flush(pq, batch);

  close_archive_view(&archive);
  free_index_info(info);
}

typedef struct
{
  const char **zstFiles;
  const char **idxFiles;
  int fileCount;
  int start, end;
  PrintQueue *pqueue;
  pthread_t tid;
  int thread_id;
} WorkerArg;

static void *worker_thread(void *arg)
{
  WorkerArg *W = (WorkerArg *)arg;
  ResultBatch batch = {};
  
  // Allocate reusable buffers once per thread
  const size_t OVER = (g_pat.len > 1) ? (g_pat.len - 1) : 0;
  char *buf = (char *)malloc(ZSTD_OUT_WIN + OVER);
  char *inbuf = (char *)malloc(ZSTD_IN_CHUNK);
  ZSTD_DCtx *dctx = ZSTD_createDCtx();
  
  if (!buf || !inbuf || !dctx)
  {
    fprintf(stderr, "Thread %d: Failed to allocate buffers\n", W->thread_id);
    if (buf) free(buf);
    if (inbuf) free(inbuf);
    if (dctx) ZSTD_freeDCtx(dctx);
    return NULL;
  }
  
  for (int i = W->start; i < W->end && i < W->fileCount; ++i)
  {
    if (!W->zstFiles[i] || !W->idxFiles[i])
      continue;
    search_indexed_file(W->pqueue, &batch, W->zstFiles[i], W->idxFiles[i], &g_pat, buf, inbuf, dctx);
  }
  (void)result_batch_flush(W->pqueue, &batch);
  
  // Free buffers once at thread completion
  ZSTD_freeDCtx(dctx);
  free(inbuf);
  free(buf);
  
  return NULL;
}

// Helper function to build full path
static char* build_path(const char *dir, const char *file)
{
  size_t dir_len = strlen(dir);
  size_t file_len = strlen(file);
  bool need_slash = (dir_len > 0 && dir[dir_len - 1] != '/');
  
  size_t total_len = dir_len + (need_slash ? 1 : 0) + file_len + 1;
  char *path = (char *)malloc(total_len);
  if (!path)
    return NULL;
  
  strcpy(path, dir);
  if (need_slash)
    strcat(path, "/");
  strcat(path, file);
  
  return path;
}

int main(int argc, char **argv)
{
  const char *search_dir = ".";
  const char *pattern = NULL;
  
  // Parse command line arguments
  int opt;
  while ((opt = getopt(argc, argv, "d:h")) != -1)
  {
    switch (opt)
    {
    case 'd':
      search_dir = optarg;
      break;
    case 'h':
      fprintf(stderr, "Usage: %s [-d directory] <pattern>\n", argv[0]);
      fprintf(stderr, "  -d directory  Directory to search (default: current directory)\n");
      fprintf(stderr, "  -h           Show this help message\n");
      return 0;
    default:
      fprintf(stderr, "Usage: %s [-d directory] <pattern>\n", argv[0]);
      return 1;
    }
  }
  
  // Get the pattern from remaining arguments
  if (optind >= argc)
  {
    fprintf(stderr, "Error: No search pattern specified\n");
    fprintf(stderr, "Usage: %s [-d directory] <pattern>\n", argv[0]);
    return 1;
  }
  pattern = argv[optind];
  
  const size_t patLen = strlen(pattern);
  if (!patLen)
  {
    fprintf(stderr, "Empty pattern.\n");
    return 1;
  }

  // Check if directory exists
  struct stat dir_stat;
  if (stat(search_dir, &dir_stat) != 0)
  {
    fprintf(stderr, "Error: Cannot access directory '%s': %s\n", search_dir, strerror(errno));
    return 1;
  }
  if (!S_ISDIR(dir_stat.st_mode))
  {
    fprintf(stderr, "Error: '%s' is not a directory\n", search_dir);
    return 1;
  }

  g_pat.pat = pattern;
  g_pat.len = patLen;
  build_badchar_table(g_pat.pat, g_pat.len, g_pat.bad);

  DIR *d = opendir(search_dir);
  if (!d)
  {
    fprintf(stderr, "opendir %s failed: %s\n", search_dir, strerror(errno));
    return 1;
  }

  size_t cap = 64, count = 0;
  char **zst = (char **)malloc(cap * sizeof(*zst));
  char **idx = (char **)malloc(cap * sizeof(*idx));
  if (!zst || !idx)
  {
    fprintf(stderr, "OOM\n");
    if (zst)
      free(zst);
    if (idx)
      free(idx);
    closedir(d);
    return 1;
  }

  struct dirent *de;
  while ((de = readdir(d)) != NULL)
  {
    if (fnmatch("*.tar.zst", de->d_name, 0) == 0)
    {
      // Build full paths for the .zst and matching .json index files
      char *zstPath = build_path(search_dir, de->d_name);
      if (!zstPath)
      {
        fprintf(stderr, "Failed to build path for %s\n", de->d_name);
        continue;
      }
      
      // Build the index filename: <archive>.json
      char idxName[MAX_FILENAME_LEN];
      int idxLen = snprintf(idxName, sizeof(idxName), "%s.json", de->d_name);
      if (idxLen < 0 || (size_t)idxLen >= sizeof(idxName))
      {
        free(zstPath);
        continue;
      }
      
      char *idxPath = build_path(search_dir, idxName);
      if (!idxPath)
      {
        fprintf(stderr, "Failed to build path for %s\n", idxName);
        free(zstPath);
        continue;
      }
      
      // Check if index file exists
      struct stat st;
      if (stat(idxPath, &st) == 0)
      {
        if (count >= cap)
        {
          cap *= 2;
          char **nz = (char **)realloc(zst, cap * sizeof(*nz));
          char **ni = (char **)realloc(idx, cap * sizeof(*ni));
          if (!nz || !ni)
          {
            fprintf(stderr, "realloc failed\n");
            free(zstPath);
            free(idxPath);
            for (size_t k = 0; k < count; k++)
            {
              free(zst[k]);
              free(idx[k]);
            }
            free(zst);
            free(idx);
            closedir(d);
            return 1;
          }
          zst = nz;
          idx = ni;
        }
        zst[count] = zstPath;
        idx[count] = idxPath;
        ++count;
      }
      else
      {
        free(zstPath);
        free(idxPath);
      }
    }
  }
  closedir(d);

  if (count == 0)
  {
    fprintf(stderr, "No *.tar.zst + *.tar.zst.json pairs found in directory: %s\n", search_dir);
    free(zst);
    free(idx);
    return 0;
  }

  PrintQueue pq;
  printqueue_init(&pq);

  pthread_t printerTid;
  if (pthread_create(&printerTid, NULL, printer_thread, &pq) != 0)
  {
    fprintf(stderr, "printer thread failed\n");
    for (size_t i = 0; i < count; i++)
    {
      free(zst[i]);
      free(idx[i]);
    }
    free(zst);
    free(idx);
    printqueue_destroy(&pq);
    return 1;
  }

  int threads = (int)((count < MAX_THREADS) ? count : MAX_THREADS);
  WorkerArg *args = (WorkerArg *)calloc((size_t)threads, sizeof(*args));
  if (!args)
  {
    fprintf(stderr, "OOM (worker args)\n");
    printqueue_mark_done(&pq);
    pthread_join(printerTid, NULL);
    for (size_t i = 0; i < count; i++)
    {
      free(zst[i]);
      free(idx[i]);
    }
    free(zst);
    free(idx);
    printqueue_destroy(&pq);
    return 1;
  }

  int per = (int)count / threads, rem = (int)count % threads, start = 0;
  for (int t = 0; t < threads; ++t)
  {
    int load = per + (t < rem ? 1 : 0);
    args[t].zstFiles = (const char **)zst;
    args[t].idxFiles = (const char **)idx;
    args[t].fileCount = (int)count;
    args[t].start = start;
    args[t].end = start + load;
    args[t].pqueue = &pq;
    args[t].thread_id = t;
    start += load;
    if (pthread_create(&args[t].tid, NULL, worker_thread, &args[t]) != 0)
    {
      fprintf(stderr, "worker %d spawn failed\n", t);
      args[t].tid = 0;
    }
  }

  for (int t = 0; t < threads; ++t)
    if (args[t].tid)
      pthread_join(args[t].tid, NULL);

  printqueue_mark_done(&pq);
  pthread_join(printerTid, NULL);

  for (size_t i = 0; i < count; i++)
  {
    free(zst[i]);
    free(idx[i]);
  }
  free(zst);
  free(idx);
  free(args);
  printqueue_destroy(&pq);
  return 0;
}
