# Provided by @adavidzh (GitHub
CXX = g++
CXXFLAGS = -O3 -std=c++11
OMPFLAGS = -fopenmp
WARNFLAGS = -Wall -pedantic -Wextra


all: neep

SRCS=$(wildcard *.C)
OBJS=$(SRCS:.C=.o )

%.o: %.C
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) -c $< -o $@

neep: $(OBJS)
	$(CXX) $(OBJS) -o neep $(OMPFLAGS)

clean:
	rm -rf *.o neep
	