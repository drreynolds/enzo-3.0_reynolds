/***********************************************************************
/
/  GRID CLASS (Search for zones within the accretion radius of an
/              accreting active particle and return weighted densities)
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

int grid::FindAverageDensityInAccretionZone(ActiveParticleType* ThisParticle, FLOAT AccretionRadius, 
					    float *WeightedSum, float *SumOfWeights, int *NumberOfCells,
					    FLOAT BondiHoyleRadius) {
  /* Return if this doesn't involve us */
  if (MyProcessorNumber != ProcessorNumber) 
    return SUCCESS;

  /* Check whether the cube that circumscribes the accretion zone intersects with this grid */

  if (GridLeftEdge[0] > ParticlePosition[0]+AccretionRadius || GridRightEdge[0] < ParticlePosition[0]-AccretionRadius ||
      GridLeftEdge[1] > ParticlePosition[1]+AccretionRadius || GridRightEdge[1] < ParticlePosition[1]-AccretionRadius ||
      GridLeftEdge[2] > ParticlePosition[2]+AccretionRadius || GridRightEdge[2] < ParticlePosition[2]-AccretionRadius)
    return SUCCESS;

  FLOAT *ParticlePosition, radius2, CellSize, KernelRadius;
  int i, j, k, dim, size=1, index;
  
  /* Get indices in BaryonField for density, internal energy, thermal energy, velocity */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  /* Compute cell width and find botoom left and top right corners of grid */

  CellSize = CellWidth[0][0];

  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  ParticlePosition = ThisParticle->ReturnPosition();
  
  if (BondiHoyleRadius < CellSize/4.0)
    KernelRadius = CellSize/4.0;
  else if (BondiHoyleRadius < AccretionRadius/2.0)
    KernelRadius = BondiHoyleRadius;
  else
    KernelRadius = AccretionRadius/2.0;
  
  for (k = 0; k < GridDimension[0]; k++) {
    for (j = 0; j < GridDimension[1]; j++) {
      index = (k*GridDimension[1] + j)*GridDimension[0];
      for (i = 0; i < GridDimension[2]; index++, i++) {
	radius2 = POW((CellLeftEdge[0][i] + 0.5*CellWidth[0][i]) - ParticlePosition[0],2) +
	  POW((CellLeftEdge[1][j] + 0.5*CellWidth[1][j]) - ParticlePosition[1],2) +
	  POW((CellLeftEdge[2][k] + 0.5*CellWidth[2][k]) - ParticlePosition[2],2);   
	if ((AccretionRadius*AccretionRadius) > radius2) {
	  *WeightedSum += BaryonField[DensNum][index]*exp(-radius2/(KernelRadius*KernelRadius)); 
	  *SumOfWeights += exp(-radius2/(KernelRadius*KernelRadius));
	  (*NumberOfCells)++;
	}
      }
    }
  }

  return SUCCESS;

}
