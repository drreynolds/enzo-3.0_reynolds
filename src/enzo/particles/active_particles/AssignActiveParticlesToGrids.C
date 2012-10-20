/***********************************************************************
/
/  AssignActiveParticlesToGrids:
/  Assigning a list of active particles to the appropriate grid.
/
/  written by: Nathan Goldbaum
/  date:       June, 2012
/  modified1:
/
/  PURPOSE: Given a list of active particles, this function will 
/           traverse the hierarchy and assign them to the appropriate 
/           grid.  If a particle with the same ID as the one under 
/           consideration is already assigned to the wrong grid, that 
/           particle is deleted.  If the particle is already assigned 
/           to the correct grid, the mirrored particle is updated.
/
************************************************************************/

#ifdef USE_MPI
#include "communicators.h"
#endif 

#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "units.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CommunicationUtilities.h"
#include "phys_constants.h"

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int CommunicationSyncNumberOfParticles(HierarchyEntry *GridHierarchyPointer[],int NumberOfGrids);

int AssignActiveParticlesToGrids(ActiveParticleType** ParticleList, int nParticles, 
				 LevelHierarchyEntry *LevelArray[]) 
{
  int LevelMax, SavedGrid, NumberOfGrids, i, level, NumberOfLevelGrids, gridnum;
  HierarchyEntry **LevelGrids = NULL;
  FLOAT* pos = NULL;
  float mass;
  
  for (i = 0; i<nParticles; i++) {
    // Find the grid and processor this particle lives on
    LevelMax = SavedGrid = -1;
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
      NumberOfLevelGrids = GenerateGridArray(LevelArray, level, &LevelGrids);     
      for (gridnum = 0; gridnum < NumberOfLevelGrids; gridnum++) 
	if (LevelGrids[gridnum]->GridData->ReturnProcessorNumber() == MyProcessorNumber)
	  if (LevelGrids[gridnum]->GridData->PointInGrid(ParticleList[i]->ReturnPosition()) == true &&
	      LevelGrids[gridnum]->GridData->isLocal() == true) { 
	    SavedGrid = gridnum;
	    LevelMax = level;
	  }
      delete [] LevelGrids;
      LevelGrids = NULL;
    }
    
    /* Assign the merged particles to grids. */
    
    if (NumberOfProcessors == 1) {
      grid* OldGrid = ParticleList[i]->ReturnCurrentGrid();
      int ID = ParticleList[i]->ReturnID();
      int foundP = FALSE; 
      int foundAP = FALSE;
      NumberOfGrids = GenerateGridArray(LevelArray, LevelMax, &LevelGrids); 
      
      // If the particle moved we need to add it to the new grid and remove
      // it from the old grid.
      if (OldGrid != LevelGrids[SavedGrid]->GridData) {
	if (LevelGrids[SavedGrid]->GridData->AddActiveParticle(static_cast<ActiveParticleType*>(ParticleList[i])) == FAIL)
	  ENZO_FAIL("Active particle grid assignment failed!\n");
	if (SavedGrid != -1) {
	  foundAP = OldGrid->RemoveActiveParticle(ID,LevelGrids[SavedGrid]->GridData->ReturnProcessorNumber());
	  foundP = OldGrid->RemoveParticle(ID);
	  if ((foundP != TRUE) || (foundAP != TRUE))
	    return FAIL;
	  OldGrid->CleanUpMovedParticles();
	}
      }
      // If the particle didn't change grids, we still need to mirror the AP data to the particle list.
      else {
	LevelGrids[SavedGrid]->GridData->UpdateParticleWithActiveParticle(ParticleList[i]->ReturnID());
      }
    }
    else {
#ifdef USE_MPI
      /* Find the processor which has the maximum value of
	 LevelMax and assign the accreting particle to the
	 SavedGrid on that processor.  */
      struct { Eint32 value; Eint32 rank; } sendbuf, recvbuf;
      MPI_Comm_rank(EnzoTopComm, &sendbuf.rank); 
      sendbuf.value = LevelMax;
      MPI_Allreduce(&sendbuf, &recvbuf, 1, MPI_2INT, MPI_MAXLOC, EnzoTopComm);
      NumberOfGrids = GenerateGridArray(LevelArray, recvbuf.value, &LevelGrids); 
      if (LevelMax == recvbuf.value) {
	if (LevelGrids[SavedGrid]->GridData->AddActiveParticle(static_cast<ActiveParticleType*>(ParticleList[i])) == FAIL) {
	  ENZO_FAIL("Active particle grid assignment failed"); 
	} 
	// Still need to mirror the AP data to the particle list.
	else {
	  LevelGrids[SavedGrid]->GridData->UpdateParticleWithActiveParticle(ParticleList[i]->ReturnID());
	}
      }
      LevelMax = recvbuf.value;
#endif // endif parallel
    }
    
    /* Sync the updated particle counts accross all proccessors */
    
    CommunicationSyncNumberOfParticles(LevelGrids, NumberOfGrids);
    
    delete [] LevelGrids;
    
  }

  return SUCCESS;
}
