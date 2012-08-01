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

float grid::SumGasMassGZ()
{

  if (MyProcessorNumber != ProcessorNumber)
    return 0;

  /* Get indices in BaryonField for density, internal energy, thermal energy, velocity */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  float CellVolume = POW(CellWidth[0][0],3);
  float MassOnGrid = 0.0;
  int size=1,index=0, i;
  
  for (i = 0; i < 3; i++)
    size *= GridDimension[i];

  for (index = 0; index < size; index++)
    MassOnGrid+=BaryonField[DensNum][index];

  for (i = 0; i < NumberOfActiveParticles; i++)
    MassOnGrid+=this->ActiveParticles[i]->ReturnMass();
  
  return MassOnGrid;//*CellVolume;

}
