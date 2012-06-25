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
  int i,n,dim;

  this->SortParticlesByNumber();

  n = NumberOfParticles - NumberOfActiveParticles;
  for (i = 0; i < NumberOfActiveParticles; i++) {
    if (this->ActiveParticles[i]->Identifier != ParticleNumber[n])
      ENZO_FAIL("Particle identifiers are inconsistent!");
    ParticleMass[n] = this->ActiveParticles[i]->Mass;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      ParticlePosition[dim][n] = this->ActiveParticles[i]->pos[dim];
      ParticleVelocity[dim][n] = this->ActiveParticles[i]->vel[dim];
    }
    n++;
  }

  return SUCCESS;
}
