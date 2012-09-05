/*********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED IN THE REFINEMENT 
/             ZONE OF AN ACTIVE PARTICLE)
/
/  written by: Nathan Goldbaum
/  date:       January, 2012
/  modified1: Stephen Skory, Sept 2012, name change.
/
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
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
#include "Grid.h"

int grid::DepositRefinementZone(int level, FLOAT* ParticlePosition,
    FLOAT RefinementRadius)
{
  /* Return if this grid is not on this processor. */

  int dim, method = 0, ParticleMassMethod, i, j, k, NumberOfFlaggedCells = 0, size=1;
  float MustRefineMass;
  FLOAT CellSize, LeftCorner[MAX_DIMENSION], RightCorner[MAX_DIMENSION];
  
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Check whether accretion zone overlaps with the grid Need to
     include ghost zones because the ParticleMassFlaggingField covers
     the ghost zones as well.*/

  CellSize = CellWidth[0][0];

  for (dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
    LeftCorner[dim] = CellLeftEdge[dim][0];
    RightCorner[dim] = LeftCorner[dim] + CellSize*FLOAT((GridDimension[dim]+1));
  }

  if ((LeftCorner[0] > ParticlePosition[0]+RefinementRadius) || (RightCorner[0] < ParticlePosition[0]-RefinementRadius) ||
      (LeftCorner[1] > ParticlePosition[1]+RefinementRadius) || (RightCorner[1] < ParticlePosition[1]-RefinementRadius) ||
      (LeftCorner[2] > ParticlePosition[2]+RefinementRadius) || (RightCorner[2] < ParticlePosition[2]-RefinementRadius))
    return SUCCESS;

  /* Error checks */

  if (ParticleMassFlaggingField == NULL)
    ENZO_FAIL("Particle Mass Flagging Field is undefined!");

  for (method = 0; method < MAX_FLAGGING_METHODS; method++) {
    if (CellFlaggingMethod[method] == 4)
      ParticleMassMethod = method;
  }

  /* Find mass that will trigger refinement */
  
  MustRefineMass = 1.001*MinimumMassForRefinement[ParticleMassMethod] *
    POW(RefineBy, level * MinimumMassForRefinementLevelExponent[ParticleMassMethod]);

  /* Temporarily set the flagging field, then we will increase the
     particle mass flagging field above the minimimum needed to trigger refinemtn */

  FlaggingField = new int[size];
  for (i = 0; i < size; i++)
    FlaggingField[i] = 0;

  /* Loop over all cells and flag the ones that overlap the accretion zone */

  int index;

  for (i = 0; i < GridDimension[0]; i++)
    for (j = 0; j < GridDimension[1]; j++)
      for (k = 0; k < GridDimension[2]; k++) {
	index = (k*GridDimension[1]+j)*GridDimension[0]+i;
	// Need to do this sort of expensive check since a derefined
	// cell could enclose the accretion zone yet still be centered
	// outside the accretion radius
	if (!((CellLeftEdge[0][i] > ParticlePosition[0]+RefinementRadius) || 
	      (CellLeftEdge[1][j] > ParticlePosition[1]+RefinementRadius) || 
	      (CellLeftEdge[2][k] > ParticlePosition[2]+RefinementRadius) || 
	      (CellLeftEdge[0][i]+CellSize < ParticlePosition[0]-RefinementRadius) ||
	      (CellLeftEdge[1][j]+CellSize < ParticlePosition[1]-RefinementRadius) ||
	      (CellLeftEdge[2][k]+CellSize < ParticlePosition[2]-RefinementRadius)))
	      {
	    FlaggingField[index] = 1;
	    NumberOfFlaggedCells++;
	}
      }

  /* Set ParticleMassFlaggingField appropriately */

  for (i = 0; i < size; i++)
    ParticleMassFlaggingField[i] += (FlaggingField[i] > 0) ? MustRefineMass : 0;

  delete [] FlaggingField;
  FlaggingField = NULL;

  return SUCCESS;
}
