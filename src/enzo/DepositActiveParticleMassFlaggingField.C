/***********************************************************************
/
/  COMMUNICATION ROUTINE: SET MASS FLAGGING FIELD FOR MUSTREFINE 
/                         ACTIVE PARTICLES
/
/  written by: Nathan Goldbaum
/  date:       January, 2011
/  modified:   
/
/  PURPOSE:
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
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
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "communication.h"
#include "SortCompareFunctions.h"

#include "ActiveParticle.h"


int DepositActiveParticleMassFlaggingField(LevelHierarchyEntry* LevelArray[],
					   int level)
{

  int i,ActiveParticleID;

  for (i = 0 ; i < EnabledActiveParticlesCount ; i++) {
    ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];
    ActiveParticleID = ActiveParticleTypeToEvaluate->GetEnabledParticleID();
    ActiveParticleTypeToEvaluate->flagging_function(LevelArray,level);
  }
  
  return SUCCESS;
    
}
