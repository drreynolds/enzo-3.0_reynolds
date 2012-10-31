/**********************************************************************
/
/  GRID CLASS (UPDATE THE PARTICLE DATA USING THE DATA STORED IN 
/              THE ACTIVE PARTICLE WITH IDENTIFIER == ID)
/
/  written by: Nathan Goldbaum
/  date:       June, 2012
/
/  PURPOSE:
/
************************************************************************/

#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "ActiveParticle.h"

int grid::UpdateParticleWithActiveParticle(PINT ID)
{
  int i,n,nFound,dim;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  nFound = -1;
  for (i = 0; i < NumberOfActiveParticles; i++) {
    for (n = 0; n < NumberOfParticles; n++) {
      if (this->ActiveParticles[i]->Identifier == ParticleNumber[n]) {
	nFound = n;
	break;
      }
    }
    if (nFound == -1)
      ENZO_FAIL("Cannot find active particle in particle list");
    ParticleMass[nFound] = this->ActiveParticles[i]->Mass;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      ParticlePosition[dim][nFound] = this->ActiveParticles[i]->pos[dim];
      ParticleVelocity[dim][nFound] = this->ActiveParticles[i]->vel[dim];
    }
  }

  return SUCCESS;
}
