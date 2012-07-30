/***********************************************************************
/
/  GET TOTAL GAS MASS ON THIS GRID
/
/  written by: Nathan Goldbaum
/  date:       June, 2012
/  modified1:
/
/
************************************************************************/
#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "phys_constants.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CosmologyParameters.h"
#include "ActiveParticle.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::SumGasMass(float *mass)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Get indices in BaryonField for density, internal energy, thermal energy, velocity */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  float CellVolume = POW(CellWidth[0][0],3);
  float MassOnGrid = 0.0;
  int i,j,k,index;
  
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	MassOnGrid+=BaryonField[DensNum][index];
      }
    }
  }

  for (i = 0; i < NumberOfActiveParticles; i++)
    MassOnGrid+=this->ActiveParticles[i]->ReturnMass();
  
  *mass += MassOnGrid;//*CellVolume;

  return SUCCESS;
}
