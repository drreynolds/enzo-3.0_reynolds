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
  FLOAT* pos;
  int i;

  /* Return if this doesn't involve us */
  if (MyProcessorNumber != ProcessorNumber) return SUCCESS;

  IsHere = false;
  pos = ThisParticle->ReturnPosition();
  if (pos[0] > GridLeftEdge[0] &&
      pos[0] < GridRightEdge[0] &&
      pos[1] > GridLeftEdge[1] &&
      pos[1] < GridRightEdge[1] &&
      pos[2] > GridLeftEdge[2] &&
      pos[2] < GridRightEdge[2]) {
    IsHere = true;
  }
  
  if (!IsHere) {
    return SUCCESS;
  }

  NumberOfActiveParticles += 1;

  /* Copy the old and new ones to a new list 
     and get rid of the old list */

  ActiveParticleType **OldActiveParticles = ActiveParticles;
  ActiveParticles = new ActiveParticleType*[NumberOfParticles];

  for (i = 0; i < NumberOfActiveParticles - 1; i++) 
    ActiveParticles[i] = OldActiveParticles[i];

  delete [] OldActiveParticles;

  ThisParticle->SetGridID(ID);
  ThisParticle->AssignCurrentGrid(this);
  ActiveParticles[NumberOfActiveParticles-1] = ThisParticle;

  return SUCCESS;
}
