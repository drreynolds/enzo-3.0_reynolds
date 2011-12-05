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

void ParticleBufferHandler::FillBuffer(ActiveParticleType *np)
{
  int dim;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    this->pos[dim] = np->pos[dim];
    this->vel[dim] = np->vel[dim];
  }
  this->Mass = np->Mass;
  this->BirthTime = np->BirthTime;
  this->DynamicalTime = np->DynamicalTime;
  this->Identifier = np->Identifier;
  this->level = np->level;
  this->GridID = np->GridID;
  this->type = np->type;
  return;
}
