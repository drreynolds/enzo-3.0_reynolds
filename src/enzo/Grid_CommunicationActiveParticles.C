/***********************************************************************
/
/  GRID CLASS (SEND ACTIVE PARTICLES FROM REAL GRID TO 'FAKE' 
/              (REPLICATED) GRID)
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:  John Wise -- re-purposed for active particles
/  date:       December, 2011
/
/  NOTES:  Adapted from grid::CommunicationSendParticles().
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <map>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "communication.h"
#include "CommunicationUtilities.h"
#include "ActiveParticle.h"

/* Send active particle from this grid to ToGrid on processor
   ToProcessor, using FromNumber particles counting from FromStart.
   Place into ToGrid at particle number ToStart. If ToStart = -1, then
   add to end. */

int grid::CommunicationSendActiveParticle(grid *ToGrid, int ToProcessor)
{

  char *buffer;
  int i, j, type, dim, index, TransferSize;
  int header_size, element_size, buffer_size;
  ActiveParticleType_info *ap_info;

  if (CommunicationShouldExit(ProcessorNumber, ToProcessor))
    return SUCCESS;

  if (NumberOfActiveParticles == 0)
    return SUCCESS;

  for (type = 0; type < EnabledActiveParticlesCount; type++) {

      ap_info = EnabledActiveParticles[type];

      /* Allocate buffer in ToProcessor.  This is automatically done
	 in StarListToBuffer in the local processor. */

      // Determine the buffer size
      header_size = ap_info->buffer_instance->ReturnHeaderSize();
      element_size = ap_info->buffer_instance->ReturnElementSize();
      size = 0;
      for (i = 0; i < NumberOfActiveParticles; i++)
	if (ActiveParticles[i]->ReturnType() == type) size++;
      TransferSize = header_size + size*element_size;

#ifdef USE_MPI
      if (CommunicationDirection == COMMUNICATION_RECEIVE)
	buffer = (char*) CommunicationReceiveBuffer[CommunicationReceiveIndex];
#endif

  /* If this is the from processor, pack fields and delete stars. */

  if (MyProcessorNumber == ProcessorNumber) {
    buffer = new char[TransferSize];
    int position = 0;
    ap_info->allocate_buffer(ActiveParticles, NumberOfActiveParticles,
			     buffer, buffer_size, position);
    for (i = 0; i < NumberOfActiveParticles; i++)
      delete ActiveParticles[i];
    delete[] ActiveParticles;
    ActiveParticles = NULL;
  }
    
  /* Send buffer. */

#ifdef USE_MPI

  /* only send if processor numbers are not identical */

  if (ProcessorNumber != ToProcessor) {

    MPI_Status status;
    MPI_Arg PCount, Count = TransferSize;
    MPI_Arg Source = ProcessorNumber;
    MPI_Arg Dest = ToProcessor;
    MPI_Arg stat;

#ifdef MPI_INSTRUMENTATION
    starttime = MPI_Wtime();
#endif
    if (MyProcessorNumber == ProcessorNumber)
      CommunicationBufferedSend(buffer, Count, MPI_PACKED, 
				Dest, MPI_SENDAP_TAG, MPI_COMM_WORLD, 
				BUFFER_IN_PLACE);

    if (MyProcessorNumber == ToProcessor) {

      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
	MPI_Irecv(buffer, Count, MPI_PACKED, Source,
		  MPI_SENDAP_TAG, MPI_COMM_WORLD,
		  CommunicationReceiveMPI_Request+CommunicationReceiveIndex);

	CommunicationReceiveGridOne[CommunicationReceiveIndex] = this;
	CommunicationReceiveGridTwo[CommunicationReceiveIndex] = ToGrid;
	CommunicationReceiveCallType[CommunicationReceiveIndex] = 20;

	CommunicationReceiveBuffer[CommunicationReceiveIndex] = (float *) buffer;
	CommunicationReceiveDependsOn[CommunicationReceiveIndex] = 
	  CommunicationReceiveCurrentDependsOn;
	CommunicationReceiveIndex++;
      }

      if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE)
	MPI_Recv(buffer, Count, MPI_PACKED, Source,
		 MPI_SENDAP_TAG, MPI_COMM_WORLD, &status);

    } // ENDIF MyProcessorNumber == ToProcessor

#ifdef MPI_INSTRUMENTATION
    /* Zhiling Lan's instrumented part */
    endtime = MPI_Wtime();
    timer[7] += endtime-starttime;
    counter[7] ++;
    timer[8] += double(TransferSize);
    timer[28] += double(TransferSize*TransferSize);
    timer[27] += (endtime-starttime)*(endtime-starttime);
#endif /* MPI_INSTRUMENTATION */
  
  } // end: if (ProcessorNumber != ToProcessor)

#endif /* USE_MPI */

  /* If this is the to processor, unpack fields. */

  if (MyProcessorNumber == ToProcessor &&
      (CommunicationDirection == COMMUNICATION_SEND_RECEIVE ||
       CommunicationDirection == COMMUNICATION_RECEIVE)) {
    ActiveParticles = new ActiveParticleType*[NumberOfActiveParticles];
    ap_info->unpack_buffer(mpi_recv_buffer
    
    RecvStars = StarBufferToList(buffer, TransferSize);
    InsertStarAfter(ToGrid->Stars, RecvStars);
    for (cstar = ToGrid->Stars; cstar; cstar = cstar->NextStar)
      cstar->CurrentGrid = ToGrid;
    delete [] buffer;
			  
  } // end: if (MyProcessorNumber...)

  return SUCCESS;
}

