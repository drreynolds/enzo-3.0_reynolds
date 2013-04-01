/**********************************************************************
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

int grid::AddActiveParticle(ActiveParticleType* ThisParticle)
{

  bool IsHere;
  FLOAT* TPpos;
  int i,j;

  /* Return if this doesn't involve us */
  if (MyProcessorNumber != ProcessorNumber) return SUCCESS;

  IsHere = false;
  TPpos = ThisParticle->ReturnPosition();
  if (TPpos[0] >= GridLeftEdge[0] &&
      TPpos[0] < GridRightEdge[0] &&
      TPpos[1] >= GridLeftEdge[1] &&
      TPpos[1] < GridRightEdge[1] &&
      TPpos[2] >= GridLeftEdge[2] &&
      TPpos[2] < GridRightEdge[2]) {
    IsHere = true;
  }
  
  /* We should already have checked if the particle is on this grid so this should
     never happen */
  if (!IsHere) {
    return FAIL;
  }

  /* Copy the old and new active particles to a new list 
     and get rid of the old list */

  /* If this particle is already on the list, it needs to be moved to
     the end of the list. This needs to happen since the copy of the
     particle in the grid list needs to be updated */

  int iskip = -1;
  for (i = 0; i < NumberOfActiveParticles; i++) 
    if (ThisParticle->ReturnID() == ActiveParticles[i]->ReturnID()) {
	iskip = i;   
    }

  if (iskip != -1) {
    NumberOfActiveParticles--;
  }

  ActiveParticleType **OldActiveParticles = ActiveParticles;
  ActiveParticles = new ActiveParticleType*[NumberOfActiveParticles+1]();
  
  j = 0;
  if (NumberOfActiveParticles > 0) {
    for (i = 0; i <= NumberOfActiveParticles; i++) {
      if (i == iskip) {
	delete OldActiveParticles[i];
	OldActiveParticles[i] = NULL;
	continue;
      } else {
	ActiveParticles[j] = OldActiveParticles[i];    
	j++;
      }
    }
  } 
  else if (iskip != -1) {
    delete OldActiveParticles[0];
    OldActiveParticles[0] = NULL;
  }

  ThisParticle->SetGridID(ID);
  ThisParticle->AssignCurrentGrid(this);
  ActiveParticles[NumberOfActiveParticles++] = ThisParticle;

  delete [] OldActiveParticles;
  
  return SUCCESS;
}
