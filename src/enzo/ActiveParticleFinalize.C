/***********************************************************************
/
/  ACTIVE PARTICLE INITIALIZATION
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:  November, 2011 (JHW) -- converting to active particles
/
/  PURPOSE: Contains all routines to finalize the star particles.
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif
#include <map>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "performance.h"
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

/* prototypes */

int CommunicationUpdateStarParticleCount(HierarchyEntry *Grids[],
					 TopGridData *MetaData,
					 int NumberOfGrids,
					 int TotalStarParticleCountPrevious[]);


int ActiveParticleFinalize(HierarchyEntry *Grids[], TopGridData *MetaData,
			   int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			   int level, int TotalStarParticleCountPrevious[])
{
  int i;

  if (!StarParticleCreation && !StarParticleFeedback)
    return SUCCESS;

  LCAPERF_START("ActiveParticleFinalize");

  /* Update the star particle counters. */

  // Needs to be implimented

  //CommunicationUpdateActiveParticleCount(Grids, MetaData, NumberOfGrids,
  //                                       TotalStarParticleCountPrevious);

  /* Update position and velocity of star particles from the actual
     particles */

  int grid_num;
  for (grid_num = 0; grid_num < NumberOfGrids; grid_num++) {
    Grids[grid_num]->GridData->MirrorActiveParticles(COPY_IN);
  } // ENDFOR grids
  
  /* Call finalization routines for each active particle type  */

  int ActiveParticleID;

  for (i = 0 ; i < EnabledActiveParticlesCount; i++) {
    
    ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];
    ActiveParticleID = ActiveParticleTypeToEvaluate->GetEnabledParticleID();
    

    ActiveParticleTypeToEvaluate->after_evolvelevel_function(Grids,MetaData,NumberOfGrids,LevelArray, 
							     level,TotalStarParticleCountPrevious,
							     ActiveParticleID);

  }

  LCAPERF_STOP("ActiveParticleFinalize");
  return SUCCESS;

}
