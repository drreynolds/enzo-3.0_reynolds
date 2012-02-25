/***********************************************************************
/
/  GRID CLASS (DETACH ACTIVE PARTICLE DATA FROM PARTICLE ARRAYS)
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

int grid::DetachActiveParticles(void)
{

  if (NumberOfActiveParticles == 0)
    return SUCCESS;

  if (MyProcessorNumber != ProcessorNumber) {
    this->NumberOfParticles -= NumberOfActiveParticles;
    return SUCCESS;
  }

  /* Sort the particles by ID, so the active particles are at the end
     of the arrays. */

  this->SortParticlesByNumber();

  int NewNumberOfParticles = NumberOfParticles - NumberOfActiveParticles;

  /* Create new particle arrays */

  int i, index, dim;
  FLOAT *pos[MAX_DIMENSION];
  float *vel[MAX_DIMENSION];
  float *Mass;
  PINT *Number;

  if (NewNumberOfParticles > 0) {
  
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      pos[dim] = new FLOAT[NewNumberOfParticles];
      vel[dim] = new float[NewNumberOfParticles];
    }
    Mass = new float[NewNumberOfParticles];
    Number = new PINT[NewNumberOfParticles];

    /* Copy normal particles.  All active particles are located at the
       end of the arrays. */
    
    for (i = 0; i < NewNumberOfParticles; i++) {
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
	pos[dim][i] = ParticlePosition[dim][i];
	vel[dim][i] = ParticleVelocity[dim][i];
      }
      Mass[i] = ParticleMass[i];
      Number[i] = ParticleNumber[i];
    }
  } // ENDIF NumberOfParticles > 0
  else {
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      pos[dim] = NULL;
      vel[dim] = NULL;
    }
    Mass = NULL;
    Number = NULL;
  }

  /* Copy active particle data in the normal particle arrays to
     ActiveParticle variable */

  
  for (i = 0, index = NewNumberOfParticles; i < NumberOfActiveParticles; 
       i++, index++) {
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      ActiveParticles[i]->pos[dim] = ParticlePosition[dim][index];
      ActiveParticles[i]->vel[dim] = ParticleVelocity[dim][index];
    }
    ActiveParticles[i]->Mass = ParticleMass[index];
  }
       
  this->DeleteParticles();
  this->SetParticlePointers(Mass, Number, pos, vel);

  NumberOfParticles -= NumberOfActiveParticles;

  return SUCCESS;

}
