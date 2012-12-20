/***********************************************************************
/
/  ACTIVE PARTICLE HELPER ROUTINE:
/   Sets the number of active particles for each type in this grid.
/
/  written by: Stephen Skory
/  date:       December 2012
/
*************************************************************************/

#ifdef USE_MPI
#include "communicators.h"
#endif /* USE_MPI */

#include "preincludes.h"
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


int grid::SetActiveParticleTypeCounts(int types[MAX_ACTIVE_PARTICLE_TYPES+2]) {

  // Return if this does not concern us
  if (MyProcessorNumber != ProcessorNumber)
    return 0;
  
  for (int j = 0; j<MAX_ACTIVE_PARTICLE_TYPES; j++) {
    ActiveParticleTypeCount[j] = types[j+2];
  }
  return SUCCESS;
}
