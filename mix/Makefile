CXXFLAGS=-O2 -Wall -W -std=c++11

all: mix.o mixfunc.o
	$(CXX) $(CXXFLAGS) -o mix mix.o mixfunc.o -lgsl -lgslcblas

mix.o: mix.h
