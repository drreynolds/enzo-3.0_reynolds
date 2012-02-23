/***********************************************************************
/
/  COMMUNICATION ROUTINE: DISTRIBUTE ACTIVE PARTICLES TO PROCESSORS
/
/  written by: John Wise
/  date:       May, 2009
/  modified1   July, 2009 by John Wise: adapted for stars
/  modified2:  December, 2011 by John Wise: adapted for active particles
/
/  PURPOSE: Takes a list of active particles moves and sends/receives
/           them to all processors
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <map>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "ActiveParticle.h"
#include "SortCompareFunctions.h"

void my_exit(int status);

int CommunicationShareActiveParticles(int *NumberToMove, 
				      ActiveParticleType** &SendList,
				      int &NumberOfReceives, 
				      ActiveParticleType** &SharedList)
{

  int i, type, proc;
  int NumberOfSends, NumberOfNewParticles, NumberOfNewParticlesThisProcessor;
  ActiveParticleType_info *ap_info;

  int TotalNumberToMove = 0;
  for (proc = 0; proc < NumberOfProcessors; proc++)
    TotalNumberToMove += NumberToMove[proc];
  //std::sort(SendList, SendList+TotalNumberToMove, cmp_ap_proc());

  SharedList = NULL;

#ifdef USE_MPI
  MPI_Arg Count;
  MPI_Arg SendCount;
  MPI_Arg RecvCount;
  MPI_Arg stat;
#endif /* USE_MPI */

  if (NumberOfProcessors > 1) {

#ifdef USE_MPI

    /* We will have one collective communication per particle type.
       In the future, there might be a way to consolidate all of the
       buffers into one buffer and communicate it as a whole. */

    for (type = 0; type < EnabledActiveParticlesCount; type++) {

      ap_info = EnabledActiveParticles[type];

      /* Create a MPI packed buffer from the active particles */
      
      Eint32 position = 0, position_on_proc = 0;
      int count, header_size, element_size, size;
      int *mpi_buffer_size, *mpi_recv_buffer_size;
      char *mpi_buffer, *mpi_recv_buffer, *temp_buffer;
      mpi_buffer_size = new int[NumberOfProcessors];
      mpi_recv_buffer_size = new int[NumberOfProcessors];

      // First determine the buffer size, then we can fill it.
      header_size = ap_info->buffer_instance->ReturnHeaderSize();
      element_size = ap_info->buffer_instance->ReturnElementSize();
      size = 0;
      for (i = 0; i < TotalNumberToMove; i++)
	if (SendList[i]->ReturnType() == type) size++;
      mpi_buffer = new char[header_size+size*element_size];

      // Pack the buffer, ordered by destination processor
      position = 0;
      temp_buffer = mpi_buffer;
      for (proc = 0; proc < NumberOfProcessors; proc++) {
	position_on_proc = 0;
	ap_info->allocate_buffer(SendList, size,
				 temp_buffer, mpi_buffer_size[proc], 
				 position_on_proc, proc);
	position += position_on_proc;
	temp_buffer = mpi_buffer + position;
      }

      /* Get counts from each processor to allocate buffers. */

      MPI_Arg *MPI_SendListCount = new MPI_Arg[NumberOfProcessors];
      MPI_Arg *MPI_SendListDisplacements = new MPI_Arg[NumberOfProcessors];

      int *RecvListCount = new int[NumberOfProcessors];
      MPI_Arg *MPI_RecvListCount = new MPI_Arg[NumberOfProcessors];
      MPI_Arg *MPI_RecvListDisplacements = new MPI_Arg[NumberOfProcessors];

      int jjj;

      for (jjj = 0; jjj < NumberOfProcessors; jjj++) {
	RecvListCount[jjj] = 0;
	MPI_RecvListCount[jjj] = 0;
	MPI_RecvListDisplacements[jjj] = 0;
      }

      NumberOfSends = 0;
      for (jjj = 0; jjj < NumberOfProcessors; jjj++) {
	MPI_SendListDisplacements[jjj] = NumberOfSends;
	NumberOfSends += mpi_buffer_size[jjj];
	MPI_SendListCount[jjj] = mpi_buffer_size[jjj];
      }

      SendCount = 1;
      RecvCount = 1;

#ifdef MPI_INSTRUMENTATION
      starttime = MPI_Wtime();
#endif /* MPI_INSTRUMENTATION */

      /*****************************************
         Share the active particle type counts
      ******************************************/
    
      stat = MPI_Alltoall(mpi_buffer_size, SendCount, IntDataType,
			  mpi_recv_buffer_size, RecvCount, IntDataType, MPI_COMM_WORLD);
      if (stat != MPI_SUCCESS) my_exit(EXIT_FAILURE);

      /* Allocate buffers and generate displacement list. */

      NumberOfReceives = 0;  
      for (i = 0; i < NumberOfProcessors; i++) {
	MPI_RecvListDisplacements[i] = NumberOfReceives;
	NumberOfReceives += mpi_recv_buffer_size[i];
	MPI_RecvListCount[i] = mpi_recv_buffer_size[i];
      }

      mpi_recv_buffer = new char[NumberOfReceives];

      /********************************
          Share the active particles
      *********************************/

      stat = MPI_Alltoallv(mpi_buffer, MPI_SendListCount, MPI_SendListDisplacements,
			   MPI_PACKED,
			   mpi_recv_buffer, MPI_RecvListCount, MPI_RecvListDisplacements,
			   MPI_PACKED,
			   MPI_COMM_WORLD);
      if (stat != MPI_SUCCESS) my_exit(EXIT_FAILURE);

      /* Unpack the MPI buffers into an array of active particles */

      // Determine how many particles we have received from the buffer
      // size (NumberOfReceives is in bytes)
      NumberOfNewParticles = (NumberOfReceives - NumberOfProcessors*header_size) / element_size; 
      SharedList = new ActiveParticleType*[NumberOfNewParticles];

      // Now convert the MPI buffer
      count = 0;
      for (proc = 0; proc < NumberOfProcessors; proc++) {
	NumberOfNewParticlesThisProcessor = (MPI_RecvListCount[proc] - header_size) / element_size;
	if (NumberOfNewParticlesThisProcessor > 0)
	  ap_info->unpack_buffer(mpi_recv_buffer + MPI_RecvListDisplacements[proc], 
				 MPI_RecvListCount[proc],
				 NumberOfNewParticlesThisProcessor,
				 SharedList, count);
      }

      NumberOfReceives = NumberOfNewParticles;

#ifdef MPI_INSTRUMENTATION
      endtime = MPI_Wtime();
      timer[9] += endtime-starttime;
      counter[9] ++;
      timer[10] += double(NumberOfReceives);
      GlobalCommunication += endtime-starttime;
      CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
      delete [] MPI_SendListCount;
      delete [] MPI_SendListDisplacements;

      delete [] RecvListCount;
      delete [] MPI_RecvListCount;
      delete [] MPI_RecvListDisplacements;

      delete [] mpi_buffer_size;
      delete [] mpi_recv_buffer_size;
      delete [] mpi_buffer;

    } // ENDFOR types

#endif /* USE_MPI */    

  } // ENDIF multi-processor
  else {
    NumberOfReceives = TotalNumberToMove;
    SharedList = SendList;
  }
  
  // First sort the list by destination grid, so the searching for
  // grids is more efficient.
  //qsort(SharedList, NumberOfReceives, star_data_size, compare_star_grid);
  std::sort(SharedList, SharedList+NumberOfReceives, cmp_ap_grid());

  return SUCCESS;

}
