# Ludwig Schmidt (ludwigschmidt2@gmail.com) 2017
# Modified based on http://make.paulandlesley.org/autodep.html and user updates.

# Compiler and Flags
CXX = g++
# Consider adding -g for debug builds: CXXFLAGS = -std=c++23 -Wall -Wextra -g -O0 -fPIC
CXXFLAGS = -std=c++23 -Wall -Wextra -O3 -fPIC

# Directories
SRCDIR = src
OBJDIR = obj
GTESTDIR = external/googletest/googletest

# Source Files
CORE_SRCS = $(SRCDIR)/pcst_fast.cc $(SRCDIR)/pruning_strategy.cc
TEST_SRCS = $(SRCDIR)/pcst_fast_test.cc
PYBIND_SRCS = $(SRCDIR)/pcst_fast_pybind.cc
GTEST_SRCS_ALL = $(GTESTDIR)/src/gtest-all.cc
GTEST_SRCS_MAIN = $(GTESTDIR)/src/gtest_main.cc

# Object Files (using pattern substitution)
CORE_OBJS = $(CORE_SRCS:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)
TEST_OBJS = $(TEST_SRCS:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)
PYBIND_OBJS = $(PYBIND_SRCS:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)
GTEST_OBJS_ALL = $(GTEST_SRCS_ALL:$(GTESTDIR)/src/%.cc=$(OBJDIR)/%.o)
GTEST_OBJS_MAIN = $(GTEST_SRCS_MAIN:$(GTESTDIR)/src/%.cc=$(OBJDIR)/%.o)
GTEST_OBJS = $(GTEST_OBJS_ALL) $(GTEST_OBJS_MAIN)

# Executable/Library Names
TEST_TARGET = pcst_fast_test
PYBIND_TARGET = pcst_fast.so

# Python Config (adjust python version if needed)
PYTHON_VERSION = python3.11
PYTHON_CFLAGS = $(shell $(PYTHON_VERSION)-config --cflags)
PYTHON_LDFLAGS = $(shell $(PYTHON_VERSION)-config --ldflags)

# Include Paths
INCLUDES = -I $(SRCDIR) -I $(GTESTDIR)/include -I external/pybind11/include

# Phony targets
.PHONY: all clean test run_test pybind run_py_test

# Default target
all: $(TEST_TARGET) $(PYBIND_TARGET)

# --- Compilation Rules ---

# Rule to compile .cc files from SRCDIR to OBJDIR
$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	@mkdir -p $(OBJDIR) # Ensure OBJDIR exists
	@echo "Compiling $< -> $@"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

# Rule to compile gtest-all.cc
$(OBJDIR)/gtest-all.o: $(GTESTDIR)/src/gtest-all.cc
	@mkdir -p $(OBJDIR)
	@echo "Compiling $< -> $@"
	$(CXX) $(CXXFLAGS) -I $(GTESTDIR)/include -I $(GTESTDIR) -c -o $@ $<

# Rule to compile gtest_main.cc
$(OBJDIR)/gtest_main.o: $(GTESTDIR)/src/gtest_main.cc
	@mkdir -p $(OBJDIR)
	@echo "Compiling $< -> $@"
	$(CXX) $(CXXFLAGS) -I $(GTESTDIR)/include -c -o $@ $<

# --- Linking Rules ---

# Link the test executable
$(TEST_TARGET): $(TEST_OBJS) $(CORE_OBJS) $(GTEST_OBJS)
	@echo "Linking $@..."
	$(CXX) $(CXXFLAGS) -o $@ $^ # $^ includes all prerequisites

# Link the Python binding shared library
# $(PYBIND_OBJS) depends on the generic rule above
$(PYBIND_TARGET): $(PYBIND_OBJS) $(CORE_OBJS)
	@echo "Linking shared library $@..."
	$(CXX) $(CXXFLAGS) -shared $(INCLUDES) $(PYTHON_CFLAGS) $(PYTHON_LDFLAGS) -o $@ $^

# --- Convenience Targets ---

test: $(TEST_TARGET)
run_test: test
	@echo "Running tests..."
	./$(TEST_TARGET)

pybind: $(PYBIND_TARGET)
run_py_test: pybind
	@echo "Running Python tests..."
	python -m pytest src/test_pcst_fast.py

# Target previously named 'run_tests'
run_tests: run_test # Changed to match the pattern

# Target previously named 'run_pcst_fast_test'
run_pcst_fast_test: run_test # Changed to match the pattern

# Target previously named 'run_pcst_fast_py_test'
run_pcst_fast_py_test: run_py_test # Changed to match the pattern

# Target previously named 'pcst_fast_py'
pcst_fast_py: $(PYBIND_TARGET) # Changed to depend on the library target

# Target previously named 'pcst_fast_test'
pcst_fast_test: $(TEST_TARGET) # Changed to depend on the executable target

clean:
	@echo "Cleaning up..."
	rm -rf $(OBJDIR)
	rm -f $(TEST_TARGET)
	rm -f $(PYBIND_TARGET)
	rm -f src/pcst_fast.pyc src/__pycache__/* # Clean Python cache too

# Add header dependencies here if needed, or use auto-dependency generation
# Example (manual):
# $(OBJDIR)/pcst_fast.o: $(SRCDIR)/pcst_fast.h $(SRCDIR)/priority_queue.h $(SRCDIR)/pairing_heap.h $(SRCDIR)/pruning_strategy.h
# $(OBJDIR)/pruning_strategy.o: $(SRCDIR)/pruning_strategy.h $(SRCDIR)/pcst_fast.h
# $(OBJDIR)/pcst_fast_test.o: $(SRCDIR)/pcst_fast.h $(SRCDIR)/test_helpers.h
# $(OBJDIR)/pcst_fast_pybind.o: $(SRCDIR)/pcst_fast.h
