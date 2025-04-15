# Makefile for pcst_fast C++ library and Python bindings

CXX = g++
CXXFLAGS_BASE = -std=c++23 -O3 -Wall -Wextra -pedantic -fPIC
CXXFLAGS_RELEASE = $(CXXFLAGS_BASE) -DNDEBUG
CXXFLAGS_DEBUG = $(CXXFLAGS_BASE) -g
LDFLAGS =
LDLIBS = -pthread

SRCDIR = src
INCLUDEDIR = include
TESTDIR = tests
BINDINGSDIR = bindings/python
EXTERNALDIR = external

OBJDIR = obj
OBJDIR_RELEASE = $(OBJDIR)/release
OBJDIR_DEBUG = $(OBJDIR)/debug
BINDIR = bin
LIBDIR = lib
PYTHON_MODULE_DIR = .

PYTHON = python3
PYTHON_CONFIG = $(PYTHON)-config
RM = rm -f
RMDIR = rm -rf
MKDIR_P = mkdir -p

MIN_GCC_MAJOR := 13
MIN_GCC_MINOR := 0
MIN_CLANG_MAJOR := 15
MIN_CLANG_MINOR := 0

INCLUDES = -I$(INCLUDEDIR) \
           -I$(EXTERNALDIR)/googletest/googletest/include \
           -I$(EXTERNALDIR)/pybind11/include

PYTHON_CFLAGS := $(shell $(PYTHON_CONFIG) --cflags)
PYTHON_LDFLAGS := $(shell $(PYTHON_CONFIG) --ldflags --embed || $(PYTHON_CONFIG) --ldflags)

CORE_SRCS = $(wildcard $(SRCDIR)/*.cc $(SRCDIR)/pruning/*.cc)
TEST_NAMES = logger end_to_end pcst_core_algorithm no_pruner simple_pruner gw_pruner strong_pruner
TEST_BASE_SRCS = $(foreach test,$(TEST_NAMES),$(TESTDIR)/$(test)_test.cc)
BINDING_SRC = $(BINDINGSDIR)/pcst_fast_pybind.cc
GTEST_SRCS = $(EXTERNALDIR)/googletest/googletest/src/gtest-all.cc
GTEST_MAIN_SRCS = $(EXTERNALDIR)/googletest/googletest/src/gtest_main.cc

CORE_OBJS_RELEASE = $(patsubst $(SRCDIR)/%.cc,$(OBJDIR_RELEASE)/src/%.o,$(filter %.cc,$(CORE_SRCS)))
BINDING_OBJS_RELEASE = $(patsubst $(BINDINGSDIR)/%.cc,$(OBJDIR_RELEASE)/bindings/%.o,$(filter %.cc,$(BINDING_SRC)))
CORE_OBJS_DEBUG = $(patsubst $(SRCDIR)/%.cc,$(OBJDIR_DEBUG)/src/%.o,$(filter %.cc,$(CORE_SRCS)))
TEST_BASE_OBJS_DEBUG = $(patsubst $(TESTDIR)/%.cc,$(OBJDIR_DEBUG)/tests/%.o,$(filter %.cc,$(TEST_BASE_SRCS)))
GTEST_OBJS_DEBUG = $(patsubst $(EXTERNALDIR)/googletest/googletest/src/%.cc,$(OBJDIR_DEBUG)/gtest/%.o,$(GTEST_SRCS))
GTEST_MAIN_OBJS_DEBUG = $(patsubst $(EXTERNALDIR)/googletest/googletest/src/%.cc,$(OBJDIR_DEBUG)/gtest/%.o,$(GTEST_MAIN_SRCS))

TEST_EXECS = $(patsubst $(TESTDIR)/%.cc,$(BINDIR)/%,$(filter %.cc,$(TEST_BASE_SRCS)))
PYTHON_MODULE_NAME = pcst_fast
_PYTHON_EXT_CMD = $(PYTHON) -c "import sysconfig; suffix = sysconfig.get_config_var('EXT_SUFFIX'); print(suffix if suffix is not None else sysconfig.get_config_var('SO'))"
_PYTHON_EXT_RAW := $(shell $(_PYTHON_EXT_CMD))
ifeq ($(_PYTHON_EXT_RAW),)
    $(warning Failed to get Python EXT_SUFFIX/SO using command: $(_PYTHON_EXT_CMD))
	PYTHON_MODULE_EXT := .so
else
	PYTHON_MODULE_EXT := $(strip $(_PYTHON_EXT_RAW))
endif
ifeq ($(strip $(PYTHON_MODULE_DIR)),.)
    PYTHON_MODULE = $(PYTHON_MODULE_NAME)$(PYTHON_MODULE_EXT)
else
    PYTHON_MODULE = $(PYTHON_MODULE_DIR)/$(PYTHON_MODULE_NAME)$(PYTHON_MODULE_EXT)
endif

DIRS_TO_CREATE = $(filter-out ., $(OBJDIR_RELEASE)/src $(OBJDIR_RELEASE)/src/pruning $(OBJDIR_RELEASE)/bindings \
                   $(OBJDIR_DEBUG)/src $(OBJDIR_DEBUG)/src/pruning $(OBJDIR_DEBUG)/tests $(OBJDIR_DEBUG)/tests/pruning $(OBJDIR_DEBUG)/gtest \
                   $(BINDIR) $(LIBDIR) $(PYTHON_MODULE_DIR))

.PHONY: all test clean python_binding run_tests build_tests check_compiler astyle content

check_compiler:
	@echo "--- Checking Compiler Version ---"
	@{ \
	CXX_TO_CHECK='$(CXX)'; \
	CXX_BASENAME=$$(basename $$CXX_TO_CHECK); \
	CXX_VERSION_STR=$$($$CXX_TO_CHECK --version 2>/dev/null || echo "unknown"); \
	CXX_MAJOR_VER=0; \
	CXX_MINOR_VER=0; \
	COMPILER_TYPE=unknown; \
	IS_GCC=0; IS_CLANG=0; \
	if echo "$$CXX_BASENAME" | grep -q -E 'g\+\+'; then IS_GCC=1; \
	elif echo "$$CXX_BASENAME" | grep -q -E 'gcc'; then IS_GCC=1; \
	elif echo "$$CXX_VERSION_STR" | grep -i -q 'gcc' && ! echo "$$CXX_VERSION_STR" | grep -i -q 'clang'; then IS_GCC=1; fi; \
	if echo "$$CXX_BASENAME" | grep -q -E 'clang\+\+'; then IS_CLANG=1; \
	elif echo "$$CXX_VERSION_STR" | grep -i -q 'clang'; then IS_CLANG=1; fi; \
	if [ $$IS_GCC -eq 1 ]; then \
		COMPILER_TYPE=g++; MIN_MAJOR_REQ=$(MIN_GCC_MAJOR); MIN_MINOR_REQ=$(MIN_GCC_MINOR); \
		VER_PARSED=$$(echo "$$CXX_VERSION_STR" | grep -o -E '[0-9]+\.[0-9]+\.?[0-9]*' | head -n 1); \
		if [ -n "$$VER_PARSED" ]; then \
			CXX_MAJOR_VER=$$(echo $$VER_PARSED | cut -d. -f1); \
			CXX_MINOR_VER=$$(echo $$VER_PARSED | cut -d. -f2); \
		fi; \
	elif [ $$IS_CLANG -eq 1 ]; then \
		COMPILER_TYPE=clang++; MIN_MAJOR_REQ=$(MIN_CLANG_MAJOR); MIN_MINOR_REQ=$(MIN_CLANG_MINOR); \
		VER_PARSED=$$(echo "$$CXX_VERSION_STR" | grep -o -E 'version [0-9]+\.[0-9]+\.?[0-9]*' | head -n 1 | grep -o -E '[0-9]+\.[0-9]+\.?[0-9]*'); \
		if [ -n "$$VER_PARSED" ]; then \
			CXX_MAJOR_VER=$$(echo $$VER_PARSED | cut -d. -f1); \
			CXX_MINOR_VER=$$(echo $$VER_PARSED | cut -d. -f2); \
		fi; \
	fi; \
	echo "Compiler: $$COMPILER_TYPE ($$CXX_TO_CHECK), Version: $$CXX_MAJOR_VER.$$CXX_MINOR_VER"; \
	echo "Minimum Required: $$MIN_MAJOR_REQ.$$MIN_MINOR_REQ"; \
	if [ "$$COMPILER_TYPE" = "unknown" ]; then \
		echo "Warning: Compiler type not recognized. Skipping version check."; \
	elif ! ( [ $$CXX_MAJOR_VER -gt $$MIN_MAJOR_REQ ] || { [ $$CXX_MAJOR_VER -eq $$MIN_MAJOR_REQ ] && [ $$CXX_MINOR_VER -ge $$MIN_MINOR_REQ ]; } ); then \
		echo "Error: Minimum required $$COMPILER_TYPE version is $$MIN_MAJOR_REQ.$$MIN_MINOR_REQ. Found $$CXX_MAJOR_VER.$$CXX_MINOR_VER." >&2; \
		exit 1; \
	else \
		echo "Version check passed."; \
	fi; \
	}
	@echo "--- End Compiler Version Check ---"


all: check_compiler python_binding

astyle:
	@echo "Formatting C++ files"
	astyle --style=google --mode=c --recursive --suffix=none *.cc,*.h

content:
	@echo "Generate repository content"
	python content.py

python_binding: check_compiler $(PYTHON_MODULE)

$(PYTHON_MODULE): | $(PYTHON_MODULE_DIR) $(DIRS_TO_CREATE)

$(PYTHON_MODULE): $(BINDING_OBJS_RELEASE) $(CORE_OBJS_RELEASE)
	@echo "Linking Python module $@"
	$(CXX) $(CXXFLAGS_RELEASE) $(LDFLAGS) $(BINDING_OBJS_RELEASE) $(CORE_OBJS_RELEASE) $(PYTHON_CFLAGS) $(PYTHON_LDFLAGS) -o "$@" -shared $(LDLIBS)

test: check_compiler run_tests

build_tests: check_compiler $(TEST_EXECS)

$(BINDIR)/%_test: | $(BINDIR) $(DIRS_TO_CREATE)

$(BINDIR)/%_test: $(OBJDIR_DEBUG)/tests/%_test.o $(CORE_OBJS_DEBUG) $(GTEST_OBJS_DEBUG) $(GTEST_MAIN_OBJS_DEBUG)
	@echo "Linking test executable $@"
	$(CXX) $(CXXFLAGS_DEBUG) $(LDFLAGS) $^ -o $@ $(LDLIBS)

run_tests: build_tests
	@echo "Running tests..."
	@for testexec in $(TEST_EXECS); do \
		echo "--- Running $$testexec ---"; \
		./$$testexec || exit 1; \
		echo ""; \
	done
	@echo "All tests passed."


$(DIRS_TO_CREATE):
	@$(MKDIR_P) $@

clean:
	@echo "Cleaning build files..."
	$(RMDIR) $(OBJDIR) $(BINDIR) $(LIBDIR)
	@echo "Attempting to remove Python module pattern: $(PYTHON_MODULE)"
	$(RM) $(PYTHON_MODULE) || echo "Python module '$(PYTHON_MODULE)' not found or already removed."
	@echo "Clean complete."

DEPFILES_RELEASE = $(CORE_OBJS_RELEASE:.o=.d) $(BINDING_OBJS_RELEASE:.o=.d)
DEPFILES_DEBUG = $(CORE_OBJS_DEBUG:.o=.d) $(TEST_BASE_OBJS_DEBUG:.o=.d) $(GTEST_OBJS_DEBUG:.o=.d) $(GTEST_MAIN_OBJS_DEBUG:.o=.d)
-include $(DEPFILES_RELEASE) $(DEPFILES_DEBUG)