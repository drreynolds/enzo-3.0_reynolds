/***********************************************************************
/
/  ALLOCATE MPI PACKED BUFFER FROM ACTIVE PARTICLE ARRAY
/
/  written by: John Wise
/  date:       December, 2011
/  modified1:
/
************************************************************************/

#include "preincludes.h"

#ifdef USE_MPI
#include "communicators.h"
#endif
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"
#include "ActiveParticle.h"

int ParticleBufferHandler::_AllocateBuffer(char *buffer, Eint32 total_buffer_size,
					   int &buffer_size, Eint32 &position)
{

#ifdef USE_MPI
  int i;
  Eint32 header_size;

  /* Calculate the buffer size */

  // Header: Number of particles
  MPI_Pack_size(1, IntDataType, EnzoTopComm, &header_size);
  buffer_size = header_size + this->NumberOfBuffers * this->ElementSizeInBytes;

  /* Allocate buffer and pack the data */

  //position = 0;
  MPI_Pack(&this->NumberOfBuffers, 1, IntDataType, buffer, total_buffer_size, 
	   &position, EnzoTopComm);

  if (this->NumberOfBuffers > 0) {

    for (i = 0; i < MAX_DIMENSION; i++)
      MPI_Pack(this->pos[i], NumberOfBuffers, MY_MPIFLOAT, buffer, total_buffer_size,
	       &position, EnzoTopComm);
    for (i = 0; i < MAX_DIMENSION; i++)
      MPI_Pack(this->vel[i], NumberOfBuffers, FloatDataType, buffer, total_buffer_size,
	       &position, EnzoTopComm);
    MPI_Pack(this->Mass, NumberOfBuffers, MPI_DOUBLE, buffer, total_buffer_size,
	     &position, EnzoTopComm);
    MPI_Pack(this->BirthTime, NumberOfBuffers, FloatDataType, buffer, total_buffer_size,
	     &position, EnzoTopComm);
    MPI_Pack(this->DynamicalTime, NumberOfBuffers, FloatDataType, buffer, total_buffer_size,
	     &position, EnzoTopComm);
    MPI_Pack(this->Metallicity, NumberOfBuffers, FloatDataType, buffer, total_buffer_size,
	     &position, EnzoTopComm);
    MPI_Pack(this->Identifier, NumberOfBuffers, PINTDataType, buffer, total_buffer_size,
	     &position, EnzoTopComm);
    MPI_Pack(this->level, NumberOfBuffers, IntDataType, buffer, total_buffer_size,
	     &position, EnzoTopComm);
    MPI_Pack(this->GridID, NumberOfBuffers, IntDataType, buffer, total_buffer_size,
	     &position, EnzoTopComm);
    MPI_Pack(this->type, NumberOfBuffers, IntDataType, buffer, total_buffer_size,
	     &position, EnzoTopComm);

  } // ENDIF NumberOfBuffers > 0

#endif /* USE_MPI */
  return SUCCESS;

}
