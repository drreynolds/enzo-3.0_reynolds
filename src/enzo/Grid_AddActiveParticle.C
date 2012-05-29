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
  int i,j;

  /* Return if this doesn't involve us */
  if (MyProcessorNumber != ProcessorNumber) return SUCCESS;

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
  
  /* We should already have checked if the particle is on this grid so this should
     never happen */
  if (!IsHere) {
    return SUCCESS;
  }

  /* Copy the old and new active particles to a new list 
     and get rid of the old list */

  /* If this particle is already on the list, it needs to be moved to
     the end of the list. This needs to happen since the copy of the
     particle in the grid list needs to be updated */

  int iskip = -1;

  for (i = 0; i < NumberOfActiveParticles; i++) 
    if (ThisParticle->ReturnID() == ActiveParticles[i]->ReturnID())
	iskip = i;   

  if (iskip != -1) {
    NumberOfActiveParticles--;
  }

  ActiveParticleType **OldActiveParticles = ActiveParticles;
  ActiveParticles = new ActiveParticleType*[NumberOfActiveParticles++];

  j = 0;
  for (i = 0; i < NumberOfActiveParticles; i++) {
    if (i != iskip)
      j++;
    else
      continue;
    ActiveParticles[j] = OldActiveParticles[i];
  }

  delete [] OldActiveParticles;

  ThisParticle->SetGridID(ID);
  ThisParticle->AssignCurrentGrid(this);
  ActiveParticles[NumberOfActiveParticles-1] = ThisParticle;

  /* Update arrays for the non-active particles*/

  /* If the particle is already on the list then overwrite it. */
  int SavedIndex = -1;
  for (i = 0; i < NumberOfParticles; i++) 
    if (ParticleNumber[i] == ThisParticle->Identifier) {
      SavedIndex = i;
    }

  if (SavedIndex != -1)
    NumberOfParticles--;
  
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

  j = 0;
  for (i = 0; i < NumberOfParticles; i++) {
    if (i != SavedIndex)
      j++;
    else
      continue;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      pos[dim][j] = ParticlePosition[dim][i];
      vel[dim][j] = ParticleVelocity[dim][i];
    }
    Mass[j] = ParticleMass[i];
    Number[j] = ParticleNumber[i];
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

  if (NumberOfActiveParticles != NumberOfParticles)
    ENZO_VFAIL("Number of active particles (%d) != Number of particles (%d)",
	       NumberOfActiveParticles, NumberOfParticles);

  return SUCCESS;
}
