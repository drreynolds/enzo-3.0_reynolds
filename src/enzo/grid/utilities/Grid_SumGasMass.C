/***********************************************************************
/
/  GRID CLASS (SUMS THE GAS + AP Mass)
/
/  written by: Nathan Goldbaum
/  date:       2012
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
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
#include "fortran.def"
#include "Grid.h"
 
float grid::SumGasMass()
{
  
  if (ProcessorNumber != MyProcessorNumber)
    return 0;
  
  float mass = 0;
  int size = 1;
  int i, j, k, index;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, MetalNum;
  
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
	 Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    ENZO_FAIL("");
  }

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
     for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
       index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
       for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	 mass += BaryonField[DensNum][index];
       }
     }
  }

  for (int i = 0; i < NumberOfActiveParticles; i++)
    mass += ActiveParticles[i]->ReturnMass();

  return mass;

}
