/***********************************************************************
/
/  GRID CLASS (APPEND NEW ACTIVE PARTICLE DATA TO PARTICLE ARRAYS)
/
/  written by: John Wise
/  date:       December, 2011
/  modified1:  
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
#include "TopGridData.h"
#include "fortran.def"
#include "CosmologyParameters.h"

#include "ActiveParticle.h"

int grid::AppendNewActiveParticles(ActiveParticleType **NewParticles,
				   int NumberOfNewParticles)
{

  if (NumberOfNewParticles == 0)
    return SUCCESS;

  int OldNumberOfParticles = NumberOfParticles;
  NumberOfParticles += NumberOfNewParticles;

  /* Create new particle arrays */

  int i, index, dim;
  FLOAT *pos[MAX_DIMENSION];
  float *vel[MAX_DIMENSION];
  float *Mass;
  PINT *Number;
  
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim] = new FLOAT[NumberOfParticles];
    vel[dim] = new float[NumberOfParticles];
  }
  Mass = new float[NumberOfParticles];
  Number = new PINT[NumberOfParticles];

  /* Copy existing particles */

  for (i = 0; i < OldNumberOfParticles; i++) {
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      pos[dim][i] = ParticlePosition[dim][i];
      vel[dim][i] = ParticleVelocity[dim][i];
    }
    Mass[i] = ParticleMass[i];
    Number[i] = ParticleNumber[i];
  }

  /* Copy new active particle data */

  for (i = 0, index = OldNumberOfParticles; i < NumberOfNewParticles;
       i++, index++) {
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      pos[dim][index] = NewParticles[i]->pos[dim];
      vel[dim][index] = NewParticles[i]->vel[dim];
    }
    Mass[index] = NewParticles[i]->Mass;
    Number[index] = NewParticles[i]->Identifier;
  }

  /* Delete old particle arrays and copy new ones */

  this->DeleteParticles();
  this->SetParticlePointers(Mass, Number, pos, vel);

  return SUCCESS;

}
