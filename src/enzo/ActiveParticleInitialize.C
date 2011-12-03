/***********************************************************************
/
/  ACTIVE PARTICLE INITIALIZATION
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:  November, 2011 (JHW) -- Converted to active particles.
/
/  PURPOSE: Contains all routines to initialize the star particles.
/
************************************************************************/

#include <map>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#include <stdio.h>
#include "performance.h"
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
#include "ActiveParticle.h"

int FindTotalNumberOfParticles(LevelHierarchyEntry *LevelArray[]);
void RecordTotalActiveParticleCount(HierarchyEntry *Grids[], int NumberOfGrids,
				    int TotalActiveParticleCountPrevious[]);

int ActiveParticleInitialize(HierarchyEntry *Grids[], TopGridData *MetaData,
			     int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			     int ThisLevel, int TotalActiveParticleCountPrevious[])
{

  /* Return if this does not concern us */
  if (EnabledActiveParticlesCount == 0) return SUCCESS;

  LCAPERF_START("ActiveParticleInitialize");

  /* Set MetaData->NumberOfParticles and prepare TotalActiveParticleCountPrevious
     these are to be used in CommunicationUpdateActiveParticleCount 
     in ActiveParticleFinalize */  

  MetaData->NumberOfParticles = FindTotalNumberOfParticles(LevelArray);
  NumberOfOtherParticles = MetaData->NumberOfParticles - NumberOfActiveParticles;
  RecordTotalActiveParticleCount(Grids, NumberOfGrids, 
				 TotalActiveParticleCountPrevious);

  /* TODO: Merging */

  /* Active particle initialization
     1. mirror quantities from active to normal particles
  */

#ifdef UNUSED
  int grid_num;
  for (grid_num = 0; grid_num < NumberOfGrids; grid_num++) {
    Grids[grid_num]->GridData->MirrorActiveParticles(COPY_OUT);
  } // ENDFOR grids
#endif

  LCAPERF_STOP("ActiveParticleInitialize");
  return SUCCESS;

}
