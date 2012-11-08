/***********************************************************************
/
/  REMOVE AN ACTIVE PARTICLE BY ITS ID
/
/  written by: Nathan Goldbaum
/  date:       December 2011
/  modified1:
/
/  RETURNS: 0 = not found; 1 = removed
/
************************************************************************/

#include "preincludes.h"
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "ActiveParticle.h"

int grid::RemoveActiveParticle(PINT ID, int NewProcessorNumber)
{

  int i,j,found = FALSE;

  if (MyProcessorNumber != ProcessorNumber)
    return found;

  if (NumberOfActiveParticles == 0)
    return found;

  for (i=0; i < NumberOfActiveParticles; i++)
    if (this->ActiveParticles[i]->ReturnID() == ID) {
      found = TRUE;
      break;
    }
  
  if (found == FALSE)
    return found;

  if (NumberOfActiveParticles > 1) {
    ActiveParticleType** temp = new ActiveParticleType*[NumberOfActiveParticles-1]();
    
    for (j=0; j < i; j++)
      temp[j] = this->ActiveParticles[j];
    
    for (j=i+1; j < NumberOfActiveParticles; j++)
      temp[j-1] = this->ActiveParticles[j];
    
    // Only free memory if the particle was communicated to another
    // processor, otherwise we will create a dangling pointer in the
    // reference to this active particle in the new grid
    if (ProcessorNumber != NewProcessorNumber) {
      delete this->ActiveParticles[i];
      this->ActiveParticles[i] = NULL;
    }

    delete [] ActiveParticles;
    
    ActiveParticles = temp;
  }  else { // Removing the only AP on the list
    if (ProcessorNumber != NewProcessorNumber) {
      delete this->ActiveParticles[0];
      this->ActiveParticles[0] = NULL;
    }
    delete [] this->ActiveParticles;
    this->ActiveParticles = NULL;
  }

  this->NumberOfActiveParticles--;

  return found;

}
