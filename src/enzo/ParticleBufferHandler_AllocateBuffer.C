/***********************************************************************
/
/  ALLOCATE MPI PACKED BUFFER FROM ACTIVE PARTICLE ARRAY
/
/  written by: John Wise
/  date:       December, 2011
/  modified1:
/
************************************************************************/

#include <map>
#include <iostream>
#include <stdexcept>

#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdlib.h>
#include <stdio.h>
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

int ParticleBufferHandler::_AllocateBuffer(char *buffer, int &buffer_size,
					   int &position)
{

#ifdef USE_MPI
  int i, header_size;

  /* Calculate the buffer size */

  // Header: Number of particles
  MPI_Pack_size(1, IntDataType, MPI_COMM_WORLD, &header_size);
  buffer_size = header_size + this->NumberOfBuffers * this->ElementSizeInBytes;

  /* Allocate buffer and pack the data */

  buffer = new char[buffer_size];
  position = 0;
  MPI_Pack(&this->NumberOfBuffers, 1, IntDataType, buffer, buffer_size, 
	   &position, MPI_COMM_WORLD);
  for (i = 0; i < MAX_DIMENSION; i++)
    MPI_Pack(this->pos[i], NumberOfBuffers, MY_MPIFLOAT, buffer, buffer_size,
	     &position, MPI_COMM_WORLD);
  for (i = 0; i < MAX_DIMENSION; i++)
    MPI_Pack(this->vel[i], NumberOfBuffers, FloatDataType, buffer, buffer_size,
	     &position, MPI_COMM_WORLD);
  MPI_Pack(this->Mass, NumberOfBuffers, MPI_DOUBLE, buffer, buffer_size,
	   &position, MPI_COMM_WORLD);
  MPI_Pack(this->BirthTime, NumberOfBuffers, FloatDataType, buffer, buffer_size,
	   &position, MPI_COMM_WORLD);
  MPI_Pack(this->DynamicalTime, NumberOfBuffers, FloatDataType, buffer, buffer_size,
	   &position, MPI_COMM_WORLD);
  MPI_Pack(this->Metallicity, NumberOfBuffers, FloatDataType, buffer, buffer_size,
	   &position, MPI_COMM_WORLD);
  MPI_Pack(this->Identifier, NumberOfBuffers, PINTDataType, buffer, buffer_size,
	   &position, MPI_COMM_WORLD);
  MPI_Pack(this->level, NumberOfBuffers, IntDataType, buffer, buffer_size,
	   &position, MPI_COMM_WORLD);
  MPI_Pack(this->GridID, NumberOfBuffers, IntDataType, buffer, buffer_size,
	   &position, MPI_COMM_WORLD);
  MPI_Pack(this->type, NumberOfBuffers, IntDataType, buffer, buffer_size,
	   &position, MPI_COMM_WORLD);

#endif /* USE_MPI */
  return SUCCESS;

}
