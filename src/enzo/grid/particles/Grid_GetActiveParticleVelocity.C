/***********************************************************************
/
/  GRID CLASS (RETURN ACTIVE PARTICLE VELOCITIES AS A 2D POINTER ARRAY)
/
/  written by: Nathan Goldbaum
/  date:       November, 2012
/  modified1:  
/
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
#include "TopGridData.h"
#include "fortran.def"
#include "CosmologyParameters.h"

#include "ActiveParticle.h"

void grid::GetActiveParticleVelocity(float *ActiveParticleVelocity[]) 
{
  int i, dim;

  for (i = 0; i < NumberOfActiveParticles; i++) {
    FLOAT* vel = ActiveParticles[i]->ReturnVelocity();
    for (dim = 0; dim < GridRank; dim++)
      ActiveParticleVelocity[dim][i] = vel[dim];
  }

}
