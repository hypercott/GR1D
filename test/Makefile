FC=gfortran-6
CC=gcc-6

all: test_call_C

test_call_C: testic.o test_call_C.o
	$(FC) -o test_call_C test_call_C.o testic.o

testic.o: testic.c
	$(CC) -c testic.c

test_call_C.o: test_call_C.F90
	$(FC) -c test_call_C.F90

