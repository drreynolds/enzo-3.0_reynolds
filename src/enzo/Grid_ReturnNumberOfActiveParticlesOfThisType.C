/***********************************************************************
/
/  ACTIVE PARTICLE HELPER ROUTINE:
/   Return Number of Active Particles with ID = ActiveParticleIDToFind
/
/  written by: Nathan Goldbaum
/  date:       March 2012
/
*************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <map>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
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


int grid::ReturnNumberOfActiveParticlesOfThisType(int ActiveParticleIDToFind) {
  
  int NumberOfActiveParticlesOfThisType = 0;
  for (int j = 0; j<NumberOfActiveParticles; j++) {
    if (this->ActiveParticles[j]->GetEnabledParticleID() == ActiveParticleIDToFind) {
      NumberOfActiveParticlesOfThisType++;
    }
  }
  return NumberOfActiveParticlesOfThisType;
}
