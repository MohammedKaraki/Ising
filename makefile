EXECUTABLE = ising.x
SOURCES = main.cpp
OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))

CXX = clang++
CXX_GENERAL_FLAGS = -g -O2  -Wall -Wpedantic -Wextra

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) -o $(EXECUTABLE) $(OBJECTS)

$(OBJECTS): %.o: %.cpp
	$(CXX) $(CXX_GENERAL_FLAGS)  -c $< -o $@

.dep: $(SOURCES)
	rm -f ./.dep
	$(CXX) $(CXX_GENERAL_FLAGS) -MM $^ > .dep

-include .dep

.PHONY: clean

clean:
	rm -rf *.o
	rm -rf .dep
