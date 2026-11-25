# SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
# SPDX-License-Identifier: GPL-3.0-or-later
# location_pca: Progressive Covering Algorithm for the Traffic Counting Location Problem
#
# Main targets (run `make help` for the full list):
#   make                 optimized build              -> bin/location_pca
#   make debug           -O0 + sanitizers build       -> bin/location_pca-debug
#   make test            unit tests + planar tests (Gurobi NOT needed)
#   make check           'test' + regression runs on instances (needs Gurobi)
#   make clean           remove build products
#
# Variables you may override on the command line, e.g. `make GUROBI_HOME=/opt/gurobi1100/linux64`:
#   GUROBI_HOME  Gurobi installation directory (default: $GUROBI_HOME, else the newest
#                /opt/gurobi*/linux64 or /Library/gurobi*/macos_universal2)
#   GUROBI_LIB   Gurobi library name without 'lib' prefix (default: detected, e.g. gurobi130)
#   NO_GUROBI=1  build without Gurobi: -a planar works, the MIP-based methods do not
#   CC           C compiler (GCC >= 13 or Clang >= 18)
#   NATIVE=1     tune the release build for this machine (-march=native)
#   SANITIZE=    disable the sanitizers of the debug and test builds

PROGRAM   := location_pca
SRC_DIR   := src
TEST_DIR  := tests
BUILD_DIR := build
BIN_DIR   := bin

.DEFAULT_GOAL := all

# ---------------------------------------------------------------------------
# Compiler and C standard
# ---------------------------------------------------------------------------

# Use the newest standard the compiler accepts; the snippet exercises the C23
# features the code relies on (nullptr, bool keyword, attributes).
C23_PROBE := int main(void) { [[maybe_unused]] bool ok = true; void *p = nullptr; return p != nullptr; }
C_STD := $(shell for std in c23 c2x; do \
             echo '$(C23_PROBE)' | $(CC) -std=$$std -x c -o /dev/null -c - >/dev/null 2>&1 \
             && { echo $$std; break; }; done)
ifeq ($(strip $(C_STD)),)
$(error '$(CC)' does not support C23 (nullptr, bool, [[attributes]]). Use GCC >= 13 or Clang >= 18, e.g. 'make CC=gcc-13')
endif

WARNINGS := -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wstrict-prototypes \
            -Wmissing-prototypes -Wformat=2 -Wvla -Wcast-qual -Wundef
CPPFLAGS += -MMD -MP
CFLAGS   ?=
LDLIBS   += -lm

SANITIZE ?= -fsanitize=address,undefined -fno-omit-frame-pointer

RELEASE_FLAGS := -O3 -DNDEBUG $(if $(NATIVE),-march=native)
DEBUG_FLAGS   := -O0 -g3 $(SANITIZE)
TEST_FLAGS    := -O1 -g $(SANITIZE)

# ---------------------------------------------------------------------------
# Gurobi
# ---------------------------------------------------------------------------

GUROBI_HOME ?= $(shell ls -d /opt/gurobi*/linux64 /Library/gurobi*/macos_universal2 2>/dev/null \
                      | sort -V 2>/dev/null | tail -n 1)
GUROBI_LIB_FILE := $(firstword $(filter-out %_light.so %_light.dylib, \
                     $(wildcard $(GUROBI_HOME)/lib/libgurobi[0-9]*.so \
                                $(GUROBI_HOME)/lib/libgurobi[0-9]*.dylib)))
GUROBI_LIB ?= $(patsubst lib%,%,$(basename $(notdir $(GUROBI_LIB_FILE))))

GUROBI_CPPFLAGS := -isystem $(GUROBI_HOME)/include
ifeq ($(NO_GUROBI),)
GUROBI_LDFLAGS  := -L$(GUROBI_HOME)/lib -Wl,-rpath,$(GUROBI_HOME)/lib
GUROBI_LDLIBS   := -l$(GUROBI_LIB)
endif

# Building without Gurobi: exact_stub.c replaces the MIP interface and the
# planar algorithm keeps working; the MIP-based methods stop at run time.
NO_GUROBI ?= $(if $(and $(GUROBI_HOME),$(GUROBI_LIB)),,1)

# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------

# Solver-independent core, shared by the program and the unit tests.
# The planar algorithm (src/planar/) is always compiled in: it never calls
# the MIP solver and needs no Gurobi at run time.
PLANAR_SOURCES := $(wildcard $(SRC_DIR)/planar/*.c)
CORE_SOURCES := $(filter-out $(SRC_DIR)/main.c $(SRC_DIR)/exact_gurobi.c $(SRC_DIR)/exact_stub.c, \
                  $(wildcard $(SRC_DIR)/*.c)) $(PLANAR_SOURCES)
# Without Gurobi, exact_stub.c replaces exact_gurobi.c and -a planar still works.
MIP_SOURCES := $(if $(NO_GUROBI),$(SRC_DIR)/exact_stub.c,$(SRC_DIR)/exact_gurobi.c)
APP_SOURCES  := $(CORE_SOURCES) $(MIP_SOURCES) $(SRC_DIR)/main.c
TEST_SOURCES := $(CORE_SOURCES) $(TEST_DIR)/unit_tests.c $(TEST_DIR)/exact_bruteforce.c

release_objects = $(patsubst %.c,$(BUILD_DIR)/release/%.o,$(APP_SOURCES))
debug_objects   = $(patsubst %.c,$(BUILD_DIR)/debug/%.o,$(APP_SOURCES))
test_objects    = $(patsubst %.c,$(BUILD_DIR)/test/%.o,$(TEST_SOURCES))

# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------

.PHONY: all debug test check clean format help gurobi-check

all: $(BIN_DIR)/$(PROGRAM)

debug: $(BIN_DIR)/$(PROGRAM)-debug

# The planar algorithm is part of the program; its C tests and the Python
# oracle (weighted graphs vs. exhaustive search) run with 'make test'.
test: $(BUILD_DIR)/test/unit_tests $(BUILD_DIR)/test/planar_embed_test $(BIN_DIR)/$(PROGRAM)
	$(BUILD_DIR)/test/unit_tests
	$(BUILD_DIR)/test/planar_embed_test
	python3 $(TEST_DIR)/planar/test_oracle.py $(BIN_DIR)/$(PROGRAM)
	python3 $(TEST_DIR)/planar/test_oracle.py $(BIN_DIR)/$(PROGRAM) --klein-marx-only

check: test
ifeq ($(NO_GUROBI),)
	@if [ -f $(TEST_DIR)/regression.sh ]; then $(TEST_DIR)/regression.sh $(BIN_DIR)/$(PROGRAM); \
	 else echo "check: $(TEST_DIR)/regression.sh is not present; skipping." >&2; fi
	python3 $(TEST_DIR)/planar/test_location.py
else
	@echo "check: skipped the Gurobi regression runs (NO_GUROBI build)." >&2
endif

gurobi-check:
ifeq ($(NO_GUROBI),)
	@if [ ! -f "$(GUROBI_HOME)/include/gurobi_c.h" ] || [ -z "$(GUROBI_LIB)" ]; then \
	    echo "error: Gurobi not found (GUROBI_HOME='$(GUROBI_HOME)')." >&2; \
	    echo "       Install Gurobi and run e.g. 'make GUROBI_HOME=/opt/gurobi1300/linux64'," >&2; \
	    echo "       or build without it with 'make NO_GUROBI=1' (-a planar only)." >&2; \
	    exit 1; \
	fi
else
	@echo "note: building without Gurobi; -a planar works, exact/pca stop at run time." >&2
endif

$(BIN_DIR)/$(PROGRAM): $(release_objects) | gurobi-check
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(RELEASE_FLAGS) $(LDFLAGS) $(GUROBI_LDFLAGS) $^ -o $@ $(GUROBI_LDLIBS) $(LDLIBS)

$(BIN_DIR)/$(PROGRAM)-debug: $(debug_objects) | gurobi-check
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(DEBUG_FLAGS) $(LDFLAGS) $(GUROBI_LDFLAGS) $^ -o $@ $(GUROBI_LDLIBS) $(LDLIBS)

$(BUILD_DIR)/test/unit_tests: $(test_objects)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(TEST_FLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD_DIR)/test/planar_embed_test: $(BUILD_DIR)/test/$(TEST_DIR)/planar/embed_test.o \
$(BUILD_DIR)/test/$(SRC_DIR)/planar/embed.o $(BUILD_DIR)/test/$(SRC_DIR)/planar/graph.o \
$(BUILD_DIR)/test/$(SRC_DIR)/planar/steiner.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(TEST_FLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

# The planar sources use VLAs bounded by piece sizes on purpose.
$(BUILD_DIR)/release/$(SRC_DIR)/planar/%.o $(BUILD_DIR)/debug/$(SRC_DIR)/planar/%.o \
$(BUILD_DIR)/test/$(SRC_DIR)/planar/%.o: WARNINGS := $(filter-out -Wvla,$(WARNINGS))

# Only exact_gurobi.c needs the Gurobi header.
$(BUILD_DIR)/release/$(SRC_DIR)/exact_gurobi.o \
$(BUILD_DIR)/debug/$(SRC_DIR)/exact_gurobi.o: CPPFLAGS += $(GUROBI_CPPFLAGS)
$(BUILD_DIR)/release/$(SRC_DIR)/exact_gurobi.o \
$(BUILD_DIR)/debug/$(SRC_DIR)/exact_gurobi.o: | gurobi-check

$(BUILD_DIR)/release/%.o: %.c
	@mkdir -p $(@D)
	$(CC) -std=$(C_STD) $(WARNINGS) $(CPPFLAGS) $(CFLAGS) $(RELEASE_FLAGS) -c $< -o $@

$(BUILD_DIR)/debug/%.o: %.c
	@mkdir -p $(@D)
	$(CC) -std=$(C_STD) $(WARNINGS) $(CPPFLAGS) $(CFLAGS) $(DEBUG_FLAGS) -c $< -o $@

$(BUILD_DIR)/test/%.o: %.c
	@mkdir -p $(@D)
	$(CC) -std=$(C_STD) $(WARNINGS) $(CPPFLAGS) -I$(SRC_DIR) $(CFLAGS) $(TEST_FLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

format:
	clang-format -i $(SRC_DIR)/*.[ch] $(TEST_DIR)/*.[ch]

help:
	@echo "Targets:"
	@echo "  all (default)  build $(BIN_DIR)/$(PROGRAM) (optimized; NATIVE=1 adds -march=native)"
	@echo "  debug          build $(BIN_DIR)/$(PROGRAM)-debug with AddressSanitizer/UBSan"
	@echo "  test           build and run the unit tests and the planar-algorithm tests"
	@echo "                 (no Gurobi required)"
	@echo "  check          'test' plus regression runs on instances (requires Gurobi)"
	@echo "  format         apply .clang-format to all sources"
	@echo "  clean          remove $(BUILD_DIR)/ and $(BIN_DIR)/"
	@echo ""
	@echo "Detected configuration:"
	@echo "  CC          = $(CC) (-std=$(C_STD))"
	@echo "  GUROBI_HOME = $(GUROBI_HOME)"
	@echo "  GUROBI_LIB  = $(GUROBI_LIB)"

-include $(release_objects:.o=.d) $(debug_objects:.o=.d) $(test_objects:.o=.d)
