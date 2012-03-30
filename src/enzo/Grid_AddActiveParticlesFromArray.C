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

#include <stdio.h>
#include <math.h>
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
#include "EventHooks.h"
#include "ActiveParticle.h"

int grid::AddActiveParticlesFromArray(ActiveParticleType** ActiveParticleList, int nParticles)
{

  int *NewIndex, Count, i;
  FLOAT pos;

  /* Return if this doesn't involve us */
  if (MyProcessorNumber != ProcessorNumber) return SUCCESS;

  if (nParticles < 1) return SUCCESS;
		      
  *NewIndex = new int[Size];
  Count = 0;
  for (i = 0; i < nParticles; i++) {
    pos = ActiveParticleList[i].ReturnPosition();
    if (pos[0] > GridLeftEdge[0] &&
	pos[0] < GridRightEdge[0] &&
	pos[1] > GridLeftEdge[1] &&
	pos[1] < GridRightEdge[1] &&
	pos[2] > GridLeftEdge[2] &&
	pos[2] < GridRightEdge[2]) {
      NewIndex[Count++] = i;
      AddedNewParticleNumber[i] = 1;
    }
  }

  if (Count == 0) {
    delete [] NewIndex;
    return SUCCESS;
  }

  NumberOfActiveParticles += Count;

  /* Copy the old and new ones to a new list 
     and get rid of the old list */

  ActiveParticleType **OldActiveParticles = ActiveParticles;
  ActiveParticles = new ActiveParticleType*[NumberOfParticles];

  for (i = 0; i < NumberOfActiveParticles - Count; i++) 
    ActiveParticles[i] = OldActiveParticles[i];

  delete [] OldActiveParticles;

  for (i = 0; i < Count; i++)
    {
      ActiveParticleList[NewIndex[i]]->SetGridID(ID);
      ActiveParticleList[NewIndex[i]]->AssignCurrentGrid(this);
      ActiveParticles[NumberOfActiveParticles - Count + i] = ActiveParticleList[NewIndex[i]];
    }

  delete [] NewIndex;

  return SUCCESS;
}
