CXXFLAGS=-O2 -Wall -W -std=c++11

all: mix.o mixfunc.o mix2.o
	$(CXX) $(CXXFLAGS) -o mix mix.o mixfunc.o -lgsl -lgslcblas
	$(CXX) $(CXXFLAGS) -o mix2 mix2.o -lgsl -lgslcblas

mix.o: mix.h
