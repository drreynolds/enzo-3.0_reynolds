/***********************************************************************
/
/  GRID CLASS (APPEND ACTIVE PARTICLE DATA TO PARTICLE ARRAYS)
/
/  written by: John Wise
/  date:       December, 2011
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/
#include <map>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <math.h>
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
#include "fortran.def"
#include "CosmologyParameters.h"

#include "ActiveParticle.h"

int grid::AppendActiveParticles(void)
{

  if (NumberOfActiveParticles == 0)
    return SUCCESS;

  if (ProcessorNumber != MyProcessorNumber) {
    this->NumberOfParticles += NumberOfActiveParticles;
    return SUCCESS;
  }

  int OldNumberOfParticles = NumberOfParticles;
  NumberOfParticles += NumberOfActiveParticles;

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

  /* Copy active particle data */

  for (i = 0, index = OldNumberOfParticles; i < NumberOfActiveParticles;
       i++, index++) {
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      pos[dim][index] = ActiveParticles[i]->pos[dim];
      vel[dim][index] = ActiveParticles[i]->vel[dim];
    }
    Mass[index] = ActiveParticles[i]->Mass;
    Number[index] = ActiveParticles[i]->Identifier;
  }

  /* Delete old particle arrays and copy new ones */

  this->DeleteParticles();
  this->SetParticlePointers(Mass, Number, pos, vel);

  return SUCCESS;

}
