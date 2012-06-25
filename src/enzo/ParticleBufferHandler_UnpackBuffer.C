/***********************************************************************
/
/  UNPACK A MPI PACKED BUFFER INTO A PARTICLE BUFFER
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
#include "communicators.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
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

int ParticleBufferHandler::_UnpackBuffer(char *buffer, int buffer_size,
					 Eint32 &position)
{

  position = 0;

#ifdef USE_MPI
  int i, nbuffers;

  /* Unpack header (number of buffers) */

  MPI_Unpack(buffer, buffer_size, &position, &nbuffers, 1, IntDataType, 
	     EnzoTopComm);
  // check if we're receiving the right number of buffers
  assert(nbuffers == this->NumberOfBuffers);

  /* Unpack data */

  if (this->NumberOfBuffers > 0) {
    for (i = 0; i < MAX_DIMENSION; i++)
      MPI_Unpack(buffer, buffer_size, &position, this->pos[i], 
		 nbuffers, MY_MPIFLOAT, EnzoTopComm);
    for (i = 0; i < MAX_DIMENSION; i++)
      MPI_Unpack(buffer, buffer_size, &position, this->vel[i], 
		 nbuffers, FloatDataType, EnzoTopComm);
    MPI_Unpack(buffer, buffer_size, &position, this->Mass, 
	       nbuffers, MPI_DOUBLE, EnzoTopComm);
    MPI_Unpack(buffer, buffer_size, &position, this->BirthTime, 
	       nbuffers, FloatDataType, EnzoTopComm);
    MPI_Unpack(buffer, buffer_size, &position, this->DynamicalTime, 
	       nbuffers, FloatDataType, EnzoTopComm);
    MPI_Unpack(buffer, buffer_size, &position, this->Metallicity, 
	       nbuffers, FloatDataType, EnzoTopComm);
    MPI_Unpack(buffer, buffer_size, &position, this->Identifier,
	       nbuffers, PINTDataType, EnzoTopComm);
    MPI_Unpack(buffer, buffer_size, &position, this->level,
	       nbuffers, IntDataType, EnzoTopComm);
    MPI_Unpack(buffer, buffer_size, &position, this->GridID,
	       nbuffers, IntDataType, EnzoTopComm);
    MPI_Unpack(buffer, buffer_size, &position, this->type,
	       nbuffers, IntDataType, EnzoTopComm);
  } // ENDIF NumberOfBuffers > 0

#endif /* USE_MPI */
  return SUCCESS;

}
