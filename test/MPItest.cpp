#include <iostream>
#include <string>
#include "mpi.h"

using namespace std;

extern "C"
void mpiinit_(){
  MPI_Init(NULL,NULL);//&argc, &argv);
}

extern "C"
void mpiwork_(){
  int my_rank, n_procs;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
  std::cout << "I am rank " << my_rank+1 <<"/"<< n_procs << std::endl;
}

