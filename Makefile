CXX = g++
CXXFLAGS = -std=c++23 -Wall -Wextra -O3 -fPIC

SRCDIR = src
OBJDIR = obj
GTESTDIR = external/googletest/googletest

CORE_SRCS = \
	$(SRCDIR)/pcst_fast.cc \
	$(SRCDIR)/logger.cc \
	$(SRCDIR)/ipruner.cc \
	$(SRCDIR)/advanced_pruner_base.cc \
	$(SRCDIR)/no_pruner.cc \
	$(SRCDIR)/simple_pruner.cc \
	$(SRCDIR)/gw_pruner.cc \
	$(SRCDIR)/strong_pruner.cc \
	$(SRCDIR)/connect_final_pruner.cc \
	$(SRCDIR)/pruner_factory.cc

TEST_SRCS = $(SRCDIR)/pcst_fast_test.cc
PYBIND_SRCS = $(SRCDIR)/pcst_fast_pybind.cc
GTEST_SRCS_ALL = $(GTESTDIR)/src/gtest-all.cc
GTEST_SRCS_MAIN = $(GTESTDIR)/src/gtest_main.cc

CORE_OBJS = $(CORE_SRCS:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)
TEST_OBJS = $(TEST_SRCS:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)
PYBIND_OBJS = $(PYBIND_SRCS:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)
GTEST_OBJS_ALL = $(GTEST_SRCS_ALL:$(GTESTDIR)/src/%.cc=$(OBJDIR)/%.o)
GTEST_OBJS_MAIN = $(GTEST_SRCS_MAIN:$(GTESTDIR)/src/%.cc=$(OBJDIR)/%.o)
GTEST_OBJS = $(GTEST_OBJS_ALL) $(GTEST_OBJS_MAIN)

TEST_TARGET = pcst_fast_test
PYBIND_TARGET = pcst_fast.so

PYTHON_VERSION = python3.11
PYTHON_CFLAGS = $(shell $(PYTHON_VERSION)-config --cflags)
PYTHON_LDFLAGS = $(shell $(PYTHON_VERSION)-config --ldflags)

$(PYBIND_OBJS): CXXFLAGS += $(PYTHON_CFLAGS)

INCLUDES = -I $(SRCDIR) -I $(GTESTDIR)/include -I external/pybind11/include

.PHONY: all clean test pybind py_test

all: $(TEST_TARGET) $(PYBIND_TARGET)

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	@mkdir -p $(OBJDIR)
	@echo "Compiling $< -> $@"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

$(OBJDIR)/gtest-all.o: $(GTESTDIR)/src/gtest-all.cc
	@mkdir -p $(OBJDIR)
	@echo "Compiling $< -> $@"
	$(CXX) $(CXXFLAGS) -I $(GTESTDIR)/include -I $(GTESTDIR) -c -o $@ $<

$(OBJDIR)/gtest_main.o: $(GTESTDIR)/src/gtest_main.cc
	@mkdir -p $(OBJDIR)
	@echo "Compiling $< -> $@"
	$(CXX) $(CXXFLAGS) -I $(GTESTDIR)/include -c -o $@ $<

$(TEST_TARGET): $(TEST_OBJS) $(CORE_OBJS) $(GTEST_OBJS)
	@echo "Linking $@..."
	$(CXX) $(CXXFLAGS) -o $@ $^

$(PYBIND_TARGET): $(PYBIND_OBJS) $(CORE_OBJS)
	@echo "Linking shared library $@..."
	$(CXX) $(CXXFLAGS) -shared $(PYTHON_LDFLAGS) -o $@ $^

test: $(TEST_TARGET)

run_test: test
	@echo "Running tests..."
	./$(TEST_TARGET)

pybind: $(PYBIND_TARGET)

py_test: pybind
	@echo "Running Python tests..."
	python -m pytest src/test_pcst_fast.py

clean:
	@echo "Cleaning up..."
	rm -rf $(OBJDIR)
	rm -f $(TEST_TARGET)
	rm -f $(PYBIND_TARGET)
	rm -f src/*.pyc src/__pycache__/*
