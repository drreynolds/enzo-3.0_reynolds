/***********************************************************************
/
/  GRID CLASS (Find the accreting particle with id equal to ID on this 
/              grid's active particle list and increase its mass and 
/              momentum by AccretedMass and AccretedMomentum 
/              respectively)
/
/  written by: Nathan Goldbaum
/  date:       April 2012
/
************************************************************************/

#include <string.h>
#include <map>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <iostream>

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

int grid::AddMassAndMomentumToAccretingParticle(float AccretedMass, float AccretedMomentum[], 
						ActiveParticleType* ThisParticle, LevelHierarchyEntry *LevelArray[]) {

  // Return if this doesn't concern us
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  int i;
  int iFound = -1;
  float CellVolume = 1.0;

  for (i = 0; i<NumberOfActiveParticles; i++) {
    if (this->ActiveParticles[i]->ReturnID() == ThisParticle->ReturnID()) {
      iFound = i;
      break;
    }
  }

  if (iFound == -1)
    return SUCCESS;

  for (i = 0; i < GridRank; i++)
    CellVolume+=CellWidth[i][0];

  float VelocityAdjustment[3] = {
    AccretedMomentum[0]/(ThisParticle->Mass*CellVolume),
    AccretedMomentum[1]/(ThisParticle->Mass*CellVolume),
    AccretedMomentum[2]/(ThisParticle->Mass*CellVolume)};

  // Masses are actually densities
  this->ActiveParticles[iFound]->AddMass(AccretedMass);
  this->ActiveParticles[iFound]->AdjustVelocity(VelocityAdjustment);

  return SUCCESS;
}
