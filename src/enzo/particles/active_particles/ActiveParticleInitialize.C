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

#define NO_DEBUG

int FindTotalNumberOfParticles(LevelHierarchyEntry *LevelArray[]);
void RecordTotalActiveParticleCount(HierarchyEntry *Grids[], int NumberOfGrids,
				    int TotalActiveParticleCountPrevious[]);

int ActiveParticleInitialize(HierarchyEntry *Grids[], TopGridData *MetaData,
			     int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			     int ThisLevel)
{

  int i;

  /* Return if this does not concern us */
  if (EnabledActiveParticlesCount == 0) return SUCCESS;

  LCAPERF_START("ActiveParticleInitialize");

  int *TotalActiveParticleCount = new int[NumberOfGrids]();

  MetaData->NumberOfParticles = FindTotalNumberOfParticles(LevelArray);
  NumberOfOtherParticles = MetaData->NumberOfParticles;// - NumberOfActiveParticles;

  if (NextActiveParticleID == INT_UNDEFINED)
    NextActiveParticleID = NumberOfOtherParticles + NumberOfActiveParticles;
  
  /* Call initialization routines for each active particle type */

  int ActiveParticleID;

  for (i = 0 ; i < EnabledActiveParticlesCount; i++) {
    
    ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];
    ActiveParticleID = ActiveParticleTypeToEvaluate->GetEnabledParticleID();

    ActiveParticleTypeToEvaluate->BeforeEvolveLevel(Grids,MetaData,NumberOfGrids,LevelArray, 
							      ThisLevel,TotalActiveParticleCount,
							      ActiveParticleID);

  }

#ifdef DEBUG

  int nParticles;
  ActiveParticleType** ParticleList = NULL;

  ParticleList = ActiveParticleFindAll(LevelArray, &nParticles, 0);

  if (nParticles > 0) {
    PINT IDList[nParticles];
    for (i = 0; i < nParticles; i++)
      IDList[i] = ParticleList[i]->ReturnID();
    std::sort(IDList, IDList + sizeof(IDList)/sizeof(IDList[0]));
    for (i = 0; i < nParticles-1; i++)
      if (IDList[i] == IDList[i+1]) {
	ENZO_FAIL("Two active particles have identical IDs"); }
  }

  if (NumberOfProcessors > 1)
    for (i = 0; i < nParticles; i++)
      delete ParticleList[i];

  delete [] ParticleList;

#endif

  delete [] TotalActiveParticleCount;

  LCAPERF_STOP("ActiveParticleInitialize");
  return SUCCESS;

}
