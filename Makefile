SHELL=bash
CXX = g++
CXXFLAGS = -O4 -std=c++11
SRC_DIR = src
BIN_DIR	= bin
OBJ_FILES = CountsTable.cpp
MAIN_FILES = mergeTagCounts.cpp countTags.cpp

BINARIES 	= $(addprefix $(BIN_DIR)/, $(MAIN_FILES:.cpp=))
OBJECTS		= $(addprefix $(SRC_DIR)/,$(OBJ_FILES:.cpp=.o))

all: $(addprefix $(BIN_DIR)/, $(MAIN_FILES:.cpp=))

clean:
	rm -f $(OBJECTS)
	rm -f $(BINARIES)
	rm -df $(BIN_DIR)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(BIN_DIR)/%: $(SRC_DIR)/%.cpp $(OBJECTS)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(LIBS) -o $@ $^

test:
	cd test && make test_ok

test_init:
	cd test && make init
	$(CXX) $(CXXFLAGS) -o $@ $^
