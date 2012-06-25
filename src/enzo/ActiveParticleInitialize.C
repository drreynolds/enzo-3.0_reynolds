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

#include "preincludes.h"
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

  int i;

  /* Return if this does not concern us */
  if (EnabledActiveParticlesCount == 0) return SUCCESS;

  LCAPERF_START("ActiveParticleInitialize");

  /* Set MetaData->NumberOfParticles and prepare TotalActiveParticleCountPrevious
     these are to be used in CommunicationUpdateActiveParticleCount 
     in ActiveParticleFinalize */  

  MetaData->NumberOfParticles = FindTotalNumberOfParticles(LevelArray);
  NumberOfOtherParticles = MetaData->NumberOfParticles;// - NumberOfActiveParticles;
  RecordTotalActiveParticleCount(Grids, NumberOfGrids, 
				 TotalActiveParticleCountPrevious);

  
  /* Active particle initialization
     1. copy quantities from active to normal particles
  */

  int grid_num;
  for (grid_num = 0; grid_num < NumberOfGrids; grid_num++) {
    Grids[grid_num]->GridData->AppendActiveParticles();
  } // ENDFOR grids

  /* 2. Call initialization routines for each active particle type */

  int ActiveParticleID;

  for (i = 0 ; i < EnabledActiveParticlesCount; i++) {
    
    ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];
    ActiveParticleID = ActiveParticleTypeToEvaluate->GetEnabledParticleID();

    ActiveParticleTypeToEvaluate->before_evolvelevel_function(Grids,MetaData,NumberOfGrids,LevelArray, 
							      ThisLevel,TotalActiveParticleCountPrevious,
							      ActiveParticleID);

  }

  LCAPERF_STOP("ActiveParticleInitialize");
  return SUCCESS;

}
