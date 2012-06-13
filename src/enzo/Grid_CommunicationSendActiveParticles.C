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

#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, 
                              int Target, int Tag, MPI_Comm CommWorld, 
			      int BufferSize);
#endif /* USE_MPI */


/* Send active particle from this grid to ToGrid on processor
   ToProcessor, using FromNumber particles counting from FromStart.
   Place into ToGrid at particle number ToStart. If ToStart = -1, then
   add to end. */

int grid::CommunicationSendActiveParticles(grid *ToGrid, int ToProcessor, bool DeleteParticles)
{

  if (CommunicationShouldExit(ProcessorNumber, ToProcessor))
    return SUCCESS;

  if (NumberOfActiveParticles == 0)
    return SUCCESS;

  char *buffer;
  Eint32 position = 0;
  Eint32 TransferSize;
  int npart, i, j, size, type, dim, index,
    NumberOfNewParticles;
  int header_size, element_size, buffer_size, ap_id;
  ActiveParticleType_info *ap_info;
  ActiveParticleType **NewParticles;

  /* Serial case */

  if (NumberOfProcessors == 1) {
    ToGrid->AddActiveParticles(this->ActiveParticles, 
			       this->NumberOfActiveParticles);
    if (ToGrid != this && DeleteParticles == true)
      this->NumberOfActiveParticles = 0;
    delete[] this->ActiveParticles;
    return SUCCESS;
  } // ENDIF serial case

  for (type = 0; type < EnabledActiveParticlesCount; type++) {

      ap_info = EnabledActiveParticles[type];

      /* Allocate buffer in ToProcessor.  This is automatically done
	 in StarListToBuffer in the local processor. */

      // Determine the buffer size (we know the particle types only on
      // the host processor.  For the destination processor, we need
      // to allocate a buffer, so make it the maximum size it can be
      // (NumberOfActiveParticles).  This can be improved if we can
      // determine the number of active particles with this type on
      // the other processor.

      header_size = ap_info->return_header_size();
      element_size = ap_info->return_element_size();
      if (MyProcessorNumber == ProcessorNumber) {
	size = 0;
	for (i = 0; i < NumberOfActiveParticles; i++)
	  if (ActiveParticles[i]->ReturnType() == type) size++;
      } else {
	size = NumberOfActiveParticles;
      }
      TransferSize = header_size + size*element_size;

#ifdef USE_MPI
      if (CommunicationDirection == COMMUNICATION_RECEIVE)
	buffer = (char*) CommunicationReceiveBuffer[CommunicationReceiveIndex];
      else
#endif
	buffer = new char[TransferSize];

  /* If this is the from processor, pack fields and delete stars. */

  if (MyProcessorNumber == ProcessorNumber) {
    position = 0;
    ap_id = ap_info->GetEnabledParticleID();
    ap_info->allocate_buffer(ActiveParticles, NumberOfActiveParticles,
			     buffer, TransferSize, buffer_size, position, ap_id, -1);
    if (DeleteParticles == true) {
      for (i = 0; i < NumberOfActiveParticles; i++)
	if (ActiveParticles[i]->ReturnType() == type)
	  this->RemoveActiveParticle(ActiveParticles[i]->ReturnID());
    }
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

    // Zero out NumberOfActiveParticles in the first pass, so we can
    // use grid::AddActiveParticles.
    if (type == 0) ToGrid->NumberOfActiveParticles = 0;

    // Extract the number of particles in the buffer
#ifdef USE_MPI
    position = 0;
    MPI_Unpack(buffer, TransferSize, &position, &NumberOfNewParticles,
	       1, IntDataType, MPI_COMM_WORLD);
#endif /* USE_MPI */    
    NewParticles = new ActiveParticleType*[NumberOfNewParticles];
    buffer_size = header_size + NumberOfNewParticles*element_size;
    npart = 0;
    ap_info->unpack_buffer(buffer, buffer_size, NumberOfNewParticles,
			   NewParticles, npart);

    for (i = 0; i < NumberOfNewParticles; i++)
      NewParticles[i]->AssignCurrentGrid(ToGrid);
    ToGrid->AddActiveParticles(NewParticles, NumberOfNewParticles);

    delete[] NewParticles;
    delete[] buffer;

  } // end: if (MyProcessorNumber...)

  } // ENDFOR particle types

  return SUCCESS;
}

