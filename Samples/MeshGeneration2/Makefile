# program executable, program folder, and library folder
TARGET = run.exe
BUILD = build
LIB := ../../GFD
SOURCE = Main.cpp

# c++ compiler and compiler flags
CXX = g++
CXXFLAGS := -g -Wall -Wextra -pedantic -std=c++11#-O3 #
N = 1

# commands for insert and delete folder
MKDIR_P ?= mkdir
RMDIR_R ?= rmdir /q/s

# source files for the program
SOURCES := $(SOURCE) $(wildcard $(LIB)/Types/*.cpp) $(wildcard $(LIB)/Mesh/*.cpp) $(wildcard $(LIB)/Output/*.cpp) $(wildcard $(LIB)/BlockDec/*.cpp) $(wildcard $(LIB)/Discrete/*.cpp)
OBJECTS := $(addprefix $(BUILD)/,$(patsubst %.cpp,%.o,$(patsubst $(LIB)/%,%,$(SOURCES))))

# general part of makefile
all: $(BUILD) $(BUILD)/$(TARGET)
$(BUILD)/$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@
$(BUILD)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
$(BUILD)/%.o: $(LIB)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
.PHONY: clean
clean:
	$(RMDIR_R) $(BUILD)
run: $(BUILD) $(BUILD)/$(TARGET)
	cd $(BUILD) && $(TARGET) && cd ../

$(BUILD):
	-md build\Types
	-md build\Mesh
	-md build\Output
	-md build\BlockDec
	-md build\Discrete