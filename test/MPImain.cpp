#include "mpi.h"

int main(int argc, char **argv){
  // call the MPI function
  MPIinit(argc, argv);
  MPIwork();
  MPI_Finalize();
}
