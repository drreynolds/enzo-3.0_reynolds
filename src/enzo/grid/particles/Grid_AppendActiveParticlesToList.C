/***********************************************************************
/
/  GRID CLASS (APPEND ACTIVE PARTICLE DATA TO AN ACTIVE PARTICLE ARRAY)
/
/  written by: Nathan Goldbaum
/  date:       March, 2012
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/
#include "preincludes.h"

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
#include "fortran.def"
#include "CosmologyParameters.h"

#include "ActiveParticle.h"

int grid::AppendActiveParticlesToList(ActiveParticleType** APArray, int offset, int search_id) {
  
  // Return if this does not concern us
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  int PNum, count=0;

  for (PNum = 0; PNum < NumberOfActiveParticles; PNum++) 
    if (search_id == ActiveParticles[PNum]->ReturnType()) 
      APArray[offset+count++] = ActiveParticles[PNum];
      
  return SUCCESS;
} 
