CXXFLAGS=-O2 -Wall -W -std=c++11

all: mix.o mixfunc.o mixtime.o
	$(CXX) $(CXXFLAGS) -o mix mix.o mixfunc.o -lgsl -lgslcblas
	$(CXX) $(CXXFLAGS) -o mixtime mixtime.o mixfunc.o -lgsl -lgslcblas

mix.o: mix.h
mixtime.o: mix.h
mixfunc.o: mix.h
