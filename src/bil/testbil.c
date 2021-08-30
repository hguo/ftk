#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <bil.h>

int main(int argc, char **argv)
{
  int np, rank; 

  int required = MPI_THREAD_MULTIPLE, provided; 
  MPI_Init_thread(&argc, &argv, required, &provided); 
  fprintf(stderr, "%d, %d\n", required, provided); 

  // MPI_Init(&argc, &argv); 
  MPI_Comm_size(MPI_COMM_WORLD, &np); 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
  
  BIL_Init(MPI_COMM_WORLD); 
  if (rank == 0) {
    const int st[] = {0, 0, 0, 0}; 
    const int sz[] = {1, 50, 500, 500}; 
    float *buf = (float*)malloc(sizeof(float)*sz[0]*sz[1]*sz[2]*sz[3]); 
    BIL_Add_block_nc(4, st, sz, "isabel01.nc", "QVAPOR", (void**)&buf); 
  } else {
    const int st[] = {0, 50, 0, 0}; 
    const int sz[] = {1, 30, 500, 500}; 
    float *buf = (float*)malloc(sizeof(float)*sz[0]*sz[1]*sz[2]*sz[3]); 
    BIL_Add_block_nc(4, st, sz, "isabel01.nc", "QVAPOR", (void**)&buf); 
  }
  BIL_Read(); 
  BIL_Finalize(); 

  MPI_Finalize(); 

  return 0; 
}
