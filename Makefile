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
	rm -rf *.o a.out check neep
check:
	@if $(CXX) $(OMPFLAGS) _.c ; then echo success ;  exit 0 ; else echo "Your compiler does not support OpenMP." ; exit 1 ; fi
	@touch check
	
