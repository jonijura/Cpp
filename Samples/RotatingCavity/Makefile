# program executable, program folder, and library folder
TARGET = run
BUILD = build
LIB := ../../GFD

# c++ compiler and compiler flags
CXX = g++
CXXFLAGS := -g -Wall -Wextra -pedantic -std=c++11


# commands for insert and delete folder
MKDIR_P ?= md
RMDIR_R ?= rmdir /q/s

# source files for the program
SOURCES := Main.cpp $(filter-out $(LIB)/Types/MpiEasy.cpp, $(wildcard $(LIB)/Types/*.cpp)) $(wildcard $(LIB)/Mesh/*.cpp) $(wildcard $(LIB)/Output/*.cpp)
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
run: $(BUILD)/$(TARGET)
	cd $(BUILD) && $(TARGET) && cd ../

$(BUILD):
	md $(BUILD)\Types
	md $(BUILD)\Mesh
	md $(BUILD)\Output