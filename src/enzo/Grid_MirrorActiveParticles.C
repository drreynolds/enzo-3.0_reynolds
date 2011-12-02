/***********************************************************************
/
/  GRID CLASS (MIRROR ACTIVE AND NORMAL PARTICLE INFO)
/
/  written by: John Wise
/  date:       November, 2011
/  modified1:  
/
/  PURPOSE:
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

int grid::MirrorActiveParticles(int direction)
{

  int i;

  LCAPERF_START("grid_MirrorActiveParticles");

  // Normal -> active particles
  if (direction == COPY_IN) {

    for (i = 0; i < NumberOfActiveParticles; i++) {
      this->ActiveParticles[i]->UpdatePositionVelocity();
    }

  }

  // Active -> normal particles
  else if (direction == COPY_OUT) {

    for (i = 0; i < NumberOfActiveParticles; i++) {
      this->ActiveParticles[i]->MirrorToParticle();
    }

  } 

  else {
    ENZO_FAIL("Bad direction value");
  }

  LCAPERF_STOP("grid_MirrorActiveParticles");
  return SUCCESS;

}
