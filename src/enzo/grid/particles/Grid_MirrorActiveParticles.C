/***********************************************************************
/
/  GRID CLASS (MIRROR ACTIVE AND NORMAL PARTICLE INFO)
/
/  written by: John Wise
/  date:       November, 2011
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/
 
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "ActiveParticle.h"

int grid::MirrorActiveParticles(void)
{

  if ((NumberOfActiveParticles == 0) || (NumberOfParticles == 0))
    return SUCCESS;

  int i, n, dim, lastID;

  LCAPERF_START("grid_MirrorActiveParticles");

  this->SortParticlesByNumber();
  this->SortActiveParticlesByNumber();
  
  // Now sorted, mark duplicates.
  lastID = ParticleNumber[0];
  for (i = 1; i < NumberOfParticles; i++) {
    if (ParticleNumber[i] == lastID) {
      ParticleMass[i] = FLOAT_UNDEFINED;
    }
    lastID = ParticleNumber[i];
  }
  
  // Clear them.
  this->CleanUpMovedParticles();

  // Normal -> active particles (position, velocity only!)  Should be
  // used before ActiveParticleHandler to get the correct position and
  // velocity for feedback

  n = NumberOfParticles - NumberOfActiveParticles;
  for (i = 0; i < NumberOfActiveParticles; i++) {
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      if (this->ActiveParticles[i]->Identifier != ParticleNumber[n]) 
	ENZO_FAIL("Particle identifiers are inconsistent!");
      this->ActiveParticles[i]->pos[dim] = ParticlePosition[dim][n];
      this->ActiveParticles[i]->vel[dim] = ParticleVelocity[dim][n];
    }
    n++;
  }

  LCAPERF_STOP("grid_MirrorActiveParticles");
  return SUCCESS;

}
