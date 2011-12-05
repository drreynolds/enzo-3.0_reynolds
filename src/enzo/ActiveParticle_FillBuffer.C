/***********************************************************************
/
/  ACTIVE PARTICLE ROUTINE (FILL BUFFER FOR COMMUNICATION)
/
/  written by: John Wise
/  date:       December, 2011
/  modified1:  
/  date:       
/
/  PURPOSE:
/
************************************************************************/

#include <string.h>
#include <map>
#include <iostream>
#include <stdexcept>
#include <typeinfo>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "hdf5.h"
#include "h5utilities.h"

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
#include "EventHooks.h"
#include "ActiveParticle.h"

ParticleBufferHandler* ActiveParticleType::FillBuffer
(ParticleBufferHandler* buffer, ActiveParticleType **particles, int NumberOfParticles,
 int start)
{

  int i, dim;
  for (i = start; i < NumberOfParticles; i++) {
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      buffer->pos[dim][i] = particles[i]->pos[dim];
      buffer->vel[dim][i] = particles[i]->vel[dim];
    }
    buffer->Mass[i] = particles[i]->Mass;
    buffer->BirthTime[i] = particles[i]->BirthTime;
    buffer->DynamicalTime[i] = particles[i]->DynamicalTime;
    buffer->Identifier[i] = particles[i]->Identifier;
    buffer->level[i] = particles[i]->level;
    buffer->GridID[i] = particles[i]->GridID;
    buffer->type[i] = particles[i]->type;
  } // ENDFOR particles
  return buffer;

}
