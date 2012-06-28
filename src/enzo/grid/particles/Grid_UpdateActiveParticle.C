/***********************************************************************
/
/  GRID CLASS (Update the active particle on this grid with the field 
/              info from id == ThisParticle->identifier)
/
/  written by: Nathan Goldbaum
/  date:       June 2012
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
#include "ActiveParticle.h"

int grid::UpdateActiveParticle(ActiveParticleType* ThisParticle) {
    // Return if this doesn't concern us
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  int i;
  int iFound = -1, pFound = -1;
  float CellVolume = 1.0;

  for (i = 0; i<NumberOfActiveParticles; i++) {
    if (this->ActiveParticles[i]->ReturnID() == ThisParticle->ReturnID()) {
      iFound = i;
      break;
    }
  }

  if (iFound == -1)
    return SUCCESS;

  this->ActiveParticles[iFound] = ThisParticle;

  ENZO_FAIL("NJG: Need to check this in a debugger\n");

  return SUCCESS;
}
