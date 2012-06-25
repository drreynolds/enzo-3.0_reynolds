
#ifdef USE_MPI
#include "communicators.h"
#endif
#include "preincludes.h"

void Reduce_Times(double time, double *time_array){
  int nprocs, my_rank; 
#ifdef USE_MPI
  MPI_Comm_size(EnzoTopComm, &nprocs);
  MPI_Comm_rank(EnzoTopComm, &my_rank);

  MPI_Gather(&time, 1, MPI_DOUBLE, 
	     time_array, 1, MPI_DOUBLE, 
	     0, EnzoTopComm);
#else
  time_array[0] = time;
#endif

  return;
}
