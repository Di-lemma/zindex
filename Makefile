CXX ?= clang++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS +=
LDFLAGS +=

ZSTD_CFLAGS := $(shell pkg-config --cflags libzstd 2>/dev/null)
ZSTD_LIBS := $(shell pkg-config --libs libzstd 2>/dev/null || printf '%s' '-lzstd')

TARGET := zindex
SRC := zindex.cpp

.PHONY: all clean check

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(ZSTD_CFLAGS) -o $@ $< $(LDFLAGS) $(ZSTD_LIBS) -pthread

clean:
	rm -f $(TARGET)
	rm -rf .check-tmp
