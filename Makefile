SHELL=bash
CXX := g++
CXXFLAGS := -O4 -std=c++11 -Wall
LIBS := -Isrc/zstr/src
libo := -lz
SRC_DIR := src
BIN_DIR	:= bin
OBJ_FILES := CountsTable.cpp
MAIN_FILES := mergeTagCounts.cpp countTags.cpp

BINARIES 	= $(addprefix $(BIN_DIR)/, $(MAIN_FILES:.cpp=))
OBJECTS		= $(addprefix $(SRC_DIR)/,$(OBJ_FILES:.cpp=.o))
VERSION := $(shell git describe --tags --always | sed 's/countTags\///')

all: $(addprefix $(BIN_DIR)/, $(MAIN_FILES:.cpp=))

clean:
	rm -f $(OBJECTS)
	rm -f $(BINARIES)
	rm -df $(BIN_DIR)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(BIN_DIR)/%: $(SRC_DIR)/%.cpp $(OBJECTS)
	mkdir -p $(BIN_DIR)
	echo "#define VERSION \"$(VERSION)\"" > src/version.h
	$(CXX) $(CXXFLAGS) $(LIBS) -o $@ $^ $(libo)

test:
	cd test && make all

test_init:
	cd test && make init
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: test

cines:
	rsync -av --delete  ./ cines:~/shared/compil/countTags/

changelog:
	gitchangelog ^countTags/0.4.5 HEAD
