/***********************************************************************
/
/  COMMUNICATION ROUTINE: SYNCHRONIZE GLOBAL ACTIVE PARTICLE LIST ACROSS
/                         PROCESSORS
/
/  written by: John Wise
/  date:       March 2009
/  modified1:  Nathan Goldbaum
/              Porting to active particles
/
/
/  PURPOSE: Generates a list of all active particles in the simulation
/           with id = ActiveParticleIDToFind.
/
*************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <map>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "ActiveParticle.h"
#include "SortCompareFunctions.h"

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);

int ActiveParticleFindAll(LevelHierarchyEntry *LevelArray[], ActiveParticleType* &AllActiveParticles, 
			  int ActiveParticleIDToFind)
{
  int i, level, type, ap_id, GridNum, TotalNumberOfActiveParticles, LocalNumberOfActiveParticles,
    header_size, element_size, buffer_size;
  ActiveParticleType *LocalActiveParticles = NULL, *GridActiveParticles = NULL;
  HierarchyEntry **Grids;
  int NumberOfGrids, *NumberOfActiveParticlesInGrids;
  ActiveParticleType_info *ap_info;

#ifdef USE_MPI
  MPI_Arg Count;
  MPI_Arg SendCount;
  MPI_Arg RecvCount;
  MPI_Arg stat;
#endif

  TotalNumberOfActiveParticles = 0;
  LocalNumberOfActiveParticles = 0;

  for (type = 0; type < EnabledActiveParticlesCount; type++) {
      
    ap_info = EnabledActiveParticles[type];
    ap_id = ap_info->GetEnabledParticleID();
    
    if (ap_id == ActiveParticleIDToFind) {
      
      /* Traverse the hierarchy and generate a buffer of all of the active 
	 particles of this type on this processor */
      
      for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
	NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
	NumberOfActiveParticlesInGrids = new int[NumberOfGrids];
	
	for(GridNum = 0; GridNum < NumberOfGrids; GridNum++) {
	  
	  NumberOfActiveParticlesInGrids[GridNum] = 0;
	  LocalNumberOfActiveParticles += Grids[GridNum]->GridData->
	    ReturnNumberOfActiveParticlesOfThisType(ActiveParticleIDToFind);
	  
	} /* ENDFOR grids */
	
	delete [] Grids;
	delete [] NumberOfActiveParticlesInGrids;
	
      } /* ENDFOR level */
      
      break;

    } /* ENDIF id == id to search for

  } /* ENDFOR active particle type */
  
    /**************************************************/
    /*                                                */
    /* Share active particle counts on all processors */
    /*                                                */
    /**************************************************/


    /**************************************************/
    /*                                                */
    /* Gather the active particles on all processors  */
    /*                                                */
    /**************************************************/

  if (NumberOfProcessors > 1) {
    
#ifdef USE_MPI
    header_size = ap_info->buffer_instance->ReturnHeaderSize();
    element_size = ap_info->buffer_instance->ReturnElementSize();

#endif /* USE_MPI */

  } // ENDIF multi-processor
  else {

  } // ENDIF serial

  return SUCCESS;

}
