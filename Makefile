# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Iinclude -std=c++17 -Wall -g

# Directories
SRC_DIR = src
INCLUDE_DIR = include
BIN_DIR = src/bin
BUILD_DIR = build

# Source files
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)

# Object files
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRC_FILES))

# Executable files
EXECUTABLES = $(BIN_DIR)/test_minimal_graph $(BIN_DIR)/test_parametric_minimal

# Default target (builds everything)
all: $(EXECUTABLES)

# Individual targets
test_minimal_graph: $(BIN_DIR)/test_minimal_graph
test_parametric_minimal: $(BIN_DIR)/test_parametric_minimal

# Rule to build test_minimal_graph
$(BIN_DIR)/test_minimal_graph: $(OBJ_FILES) $(BUILD_DIR)/test_minimal_graph.o | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Rule to build test_parametric_minimal
$(BIN_DIR)/test_parametric_minimal: $(OBJ_FILES) $(BUILD_DIR)/test_parametric_minimal.o | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Rule to compile normal source files into object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to compile test files into object files
$(BUILD_DIR)/test_minimal_graph.o: $(BIN_DIR)/test_minimal_graph.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/test_parametric_minimal.o: $(BIN_DIR)/test_parametric_minimal.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Ensure necessary directories exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Clean up build and bin directories
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)/test_minimal_graph $(BIN_DIR)/test_parametric_minimal

.PHONY: all clean test_minimal_graph test_parametric_minimal


