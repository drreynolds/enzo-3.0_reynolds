/***********************************************************************
/
/  DISABLE THE ASSOCIATED ACTIVE PARTICLE (Delete frrom particle list)
/
/  written by: John Wise
/  date:       December, 2009
/  modified1:  Nathan Goldbaum, December 2011 (porting to active particles)
/
************************************************************************/

#include <map>
#include <iostream>
#include <stdexcept>

#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdlib.h>
#include <stdio.h>
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
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"
#include "ActiveParticle.h"

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);

int ActiveParticleType::DisableParticle(LevelHierarchyEntry *LevelArray[])
{

  int i, nPart, NumberOfGrids, changedGrid = INT_UNDEFINED, found = FALSE;
  HierarchyEntry **Grids;
  
  NumberOfGrids = GenerateGridArray(LevelArray, this->level, &Grids);
  for (i = 0; i < NumberOfGrids; i++) {
    found = Grids[i]->GridData->RemoveActiveParticle(this->ReturnID());
    if (found) {
      changedGrid = i;
      break;
    }
  } // ENDFOR grids

  /* Now clean up deleted particle on the local processor and adjust
     NumberOfParticles on others */

#ifdef USE_MPI
  CommunicationAllReduceValues(&changedGrid, 1, MPI_MAX);
#endif

  if (changedGrid == INT_UNDEFINED) {
    if (debug)
      this->PrintInfo();
    ENZO_VFAIL("DisableParticle: WARNING -- "
	       "particle %"ISYM" not found...\n", this->Identifier)
  }

  Grids[changedGrid]->GridData->NumberOfActiveParticles--;

  delete [] Grids;

  return SUCCESS;

}
