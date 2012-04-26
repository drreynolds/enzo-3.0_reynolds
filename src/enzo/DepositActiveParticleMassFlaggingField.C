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

#include "ActiveParticle.h"


int DepositActiveParticleMassFlaggingField(LevelHierarchyEntry* LevelArray[],
					   int level, int TopGridDims[] )
{

  int i,ActiveParticleID;

  for (i = 0 ; i < EnabledActiveParticlesCount ; i++) {
    ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];
    ActiveParticleID = ActiveParticleTypeToEvaluate->GetEnabledParticleID();
    ActiveParticleTypeToEvaluate->flagging_function(LevelArray,level,TopGridDims,ActiveParticleID);
  }
  
  return SUCCESS;
    
}
