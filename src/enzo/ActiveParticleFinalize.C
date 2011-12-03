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

int CommunicationUpdateActiveParticleCount(HierarchyEntry *Grids[],
					 TopGridData *MetaData,
					 int NumberOfGrids,
					 int TotalActiveParticleCountPrevious[]);


int ActiveParticleFinalize(HierarchyEntry *Grids[], TopGridData *MetaData,
			   int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			   int level, int TotalActiveParticleCountPrevious[])
{

  if (EnabledActiveParticlesCount == 0) return SUCCESS;

  LCAPERF_START("ActiveParticleFinalize");

  /* Update the star particle counters. */

  CommunicationUpdateActiveParticleCount(Grids, MetaData, NumberOfGrids,
					 TotalActiveParticleCountPrevious);

  /* Update position and velocity of star particles from the actual
     particles */

#ifdef UNUSED
  int grid_num;
  for (grid_num = 0; grid_num < NumberOfGrids; grid_num++) {
    Grids[grid_num]->GridData->MirrorActiveParticles(COPY_IN);
  } // ENDFOR grids
#endif

  /* TODO Items:
     1. Add feedback spheres
     2. Accretion from grid to particles (or from pre-determined rates)
     3. Output sink particle information
     4. Change particle statuses for stellar births and deaths
     5. Reset accretion rates
     6. Store MBH mass accretion history
     7. Merge star particles.  It is probably not necessary here and
        should be done in ActiveParticleInitialize.

     JHW Thought: this should be done in a "finalize" routine for each
     active particle type.

   */

  LCAPERF_STOP("ActiveParticleFinalize");
  return SUCCESS;

}
