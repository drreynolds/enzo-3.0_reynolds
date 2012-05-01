/***********************************************************************
/
/  GRID CLASS (ADD PARTICLES BELONGING TO THIS GRID FROM A POINTER ARRAY)
/
/  written by: Peng Wang
/  date:       January, 2009
/  modified1:  Nathan Goldbaum (porting to Active Particles)
/  date:       March, 2012
/
/  PURPOSE:
/
************************************************************************/

#include <string.h>
#include <map>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <iostream>

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

int grid::AddActiveParticle(ActiveParticleType* ThisParticle)
{

  bool IsHere;
  FLOAT* TPpos;
  int i;

  /* Return if this doesn't involve us */
  if (MyProcessorNumber != ProcessorNumber) return SUCCESS;

  if (ThisParticle->ReturnID() == ID) return SUCCESS;

  IsHere = false;
  TPpos = ThisParticle->ReturnPosition();
  if (TPpos[0] > GridLeftEdge[0] &&
      TPpos[0] < GridRightEdge[0] &&
      TPpos[1] > GridLeftEdge[1] &&
      TPpos[1] < GridRightEdge[1] &&
      TPpos[2] > GridLeftEdge[2] &&
      TPpos[2] < GridRightEdge[2]) {
    IsHere = true;
  }
  
  if (!IsHere) {
    return SUCCESS;
  }

  // If this particle is already on the list, we do nothing
  for (i = 0; i < NumberOfActiveParticles; i++) 
    if (ThisParticle->ReturnID() == ActiveParticles[i]->ReturnID())
      return SUCCESS;
  
  /* Copy the old and new active particles to a new list 
     and get rid of the old list */

  ActiveParticleType **OldActiveParticles = ActiveParticles;
  ActiveParticles = new ActiveParticleType*[NumberOfParticles];

  for (i = 0; i < NumberOfActiveParticles; i++) 
    ActiveParticles[i] = OldActiveParticles[i];

  delete [] OldActiveParticles;

  ThisParticle->SetGridID(ID);
  ThisParticle->AssignCurrentGrid(this);
  ActiveParticles[NumberOfActiveParticles++] = ThisParticle;

  int OldNumberOfParticles = NumberOfParticles;
  NumberOfParticles += 1;

    /* Create new particle arrays */

  int index, dim;
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

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim][OldNumberOfParticles] = ThisParticle->pos[dim];
    vel[dim][OldNumberOfParticles] = ThisParticle->vel[dim];
  }
  Mass[OldNumberOfParticles] = ThisParticle->Mass;
  Number[OldNumberOfParticles] = ThisParticle->Identifier;

    /* Delete old particle arrays and copy new ones */

  this->DeleteParticles();
  this->SetParticlePointers(Mass, Number, pos, vel);

  return SUCCESS;
}
