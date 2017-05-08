FC=mpifort #gfortran
CC=mpicc #gcc
CXX=mpicxx

all: test_call_C

test_call_C: testic.o test_call_C.o MPItest.o
	$(CXX) -o test_call_C test_call_C.o testic.o MPItest.o -lgfortran

testic.o: testic.c
	$(CC) -c testic.c

test_call_C.o: test_call_C.F90
	$(FC) -c test_call_C.F90

MPItest.o:
	$(CXX) -c MPItest.cpp

clean:
	rm -f *.o
