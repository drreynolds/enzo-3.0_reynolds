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
 
#include <string.h>
#include <map>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <math.h>
#include <assert.h>
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

  int i, n, dim;

  LCAPERF_START("grid_MirrorActiveParticles");

  this->SortParticlesByNumber();

  // Normal -> active particles (position, velocity only!)  Should be
  // used before ActiveParticleHandler to get the correct position and
  // velocity for feedback

  n = NumberOfParticles - NumberOfActiveParticles;
  for (i = 0; i < NumberOfActiveParticles; i++) {
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      this->ActiveParticles[i]->pos[dim] = ParticlePosition[dim][n];
      this->ActiveParticles[i]->vel[dim] = ParticleVelocity[dim][n];
    }
    n++;
  }

  LCAPERF_STOP("grid_MirrorActiveParticles");
  return SUCCESS;

}
