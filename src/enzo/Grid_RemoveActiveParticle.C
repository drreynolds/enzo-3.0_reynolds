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

#include <string.h>
#include <map>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <math.h>
#include <assert.h>
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

int grid::RemoveActiveParticle(PINT ID)
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

  ActiveParticleType** temp = new ActiveParticleType*[NumberOfActiveParticles-1];

  for (j=0; j < i; j++)
    temp[j] = ActiveParticles[j];
  
  for (j=i+1; j < NumberOfActiveParticles; j++)
    temp[j-1] = ActiveParticles[j];

  delete [] ActiveParticles;

  ActiveParticles = temp;
  
  return found;

}
