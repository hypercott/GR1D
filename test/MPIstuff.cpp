#include <iostream>
#include <string>
#include "mpi.h"

using namespace std;

int MPIinit(int argc, char **argv){
  MPI_Init(&argc, &argv);
}

int MPIwork(){
  int my_rank, n_procs;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
  std::cout << "I am rank " << my_rank+1 <<"/"<< n_procs << std::endl;
}

int main(int argc, char **argv){
  // call the MPI function
  MPIinit(argc, argv);
  MPIwork();
  MPI_Finalize();
}
