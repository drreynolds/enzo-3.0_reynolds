/***********************************************************************
/
/  DETACH MIRRORED ACTIVE PARTICLES FROM GRIDS
/
/  written by: John Wise
/  date:       February, 2011
/  modified1:
/
/  PURPOSE: For the gravity solver and mass refinement fields, we
/           mirrored the active particles' positions, velocities, and 
/           masses to normal particles.  Here we remove them.
/
************************************************************************/

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

int DetachActiveParticles(LevelHierarchyEntry *LevelArray[], int level)
{

  if (EnabledActiveParticlesCount == 0)
    return SUCCESS;

  int i, grid1;
  LevelHierarchyEntry *Temp;
  
  for (i = level; i < MAX_DEPTH_OF_HIERARCHY; i++)
    for (Temp = LevelArray[i]; Temp; Temp = Temp->NextGridThisLevel)
      Temp->GridData->DetachActiveParticles();

  return SUCCESS;

}
