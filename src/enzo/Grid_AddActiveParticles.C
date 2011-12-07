/***********************************************************************
/
/  GRID CLASS (ADD ACTIVE PARTICLE DATA TO GRID)
/
/  written by: John Wise
/  date:       December, 2011
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/
#include <map>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

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

int grid::AddActiveParticles(ActiveParticleType **NewParticles,
			     int NumberOfNewParticles, int start)
{

  if (NumberOfNewParticles == 0)
    return SUCCESS;

  int i, index;
  int OldNumberOfActiveParticles = this->NumberOfActiveParticles;
  ActiveParticleType **OldActiveParticles = this->ActiveParticles;

  this->NumberOfActiveParticles += NumberOfNewParticles;
  this->ActiveParticles = new ActiveParticleType*[this->NumberOfActiveParticles];
  for (i = 0; i < OldNumberOfActiveParticles; i++) {
    this->ActiveParticles[i] = OldActiveParticles[i];
  }
  for (i = start, index = OldNumberOfActiveParticles; i < NumberOfNewParticles; 
       i++, index++) {
    this->ActiveParticles[index] = NewParticles[i];
  }
  
  delete[] OldActiveParticles; 

  return SUCCESS;

}
