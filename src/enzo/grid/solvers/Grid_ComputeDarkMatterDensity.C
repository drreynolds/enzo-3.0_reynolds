/***********************************************************************
/
/  GRID CLASS (COMPUTE THE DARK MATTER DENSITY)
/
/  written by: Stephen Skory
/  date:       August, 2012
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

int grid::ComputeDarkMatterDensity(float *DarkMatterDensity)
{
    /* Return if this doesn't concern us. */
    
    if (ProcessorNumber != MyProcessorNumber)
      return SUCCESS;

    int i, j, k, dim;
    int ActiveDim[MAX_DIMENSION];

    /* Compute the size of the field. */
    
    int size = 1;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];
    
    /* Get grid particle density field */
    float SaveGravityResolution = GravityResolution;
    GravityResolution = 1;
    this->InitializeGravitatingMassFieldParticles(RefineBy);
    this->ClearGravitatingMassFieldParticles();
    this->DepositParticlePositions(this, Time,
                     GRAVITATING_MASS_FIELD_PARTICLES);
    GravityResolution = SaveGravityResolution;
    
    /* If present, write out the GravitatingMassFieldParticles. */
    
    for (dim = 0; dim < 3; dim++)
      ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] +1;
    
    if (GravitatingMassFieldParticles != NULL) {
      
      /* Set dimensions. */
      
      int StartIndex[] = {0,0,0}, EndIndex[] = {0,0,0};
      for (dim = 0; dim < GridRank; dim++) {
        StartIndex[dim] = nint((GridLeftEdge[dim] -
                   GravitatingMassFieldParticlesLeftEdge[dim])/
                   GravitatingMassFieldParticlesCellSize);
        EndIndex[dim] = nint((GridRightEdge[dim] -
                   GravitatingMassFieldParticlesLeftEdge[dim])/
                   GravitatingMassFieldParticlesCellSize) - 1;
      }
        
      /* Copy active part of field into grid */
      
      for (k = StartIndex[2]; k <= EndIndex[2]; k++)
        for (j = StartIndex[1]; j <= EndIndex[1]; j++)
          for (i = StartIndex[0]; i <= EndIndex[0]; i++)
            DarkMatterDensity[(i-StartIndex[0]) +
             (j-StartIndex[1])*ActiveDim[0] +
             (k-StartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
              GravitatingMassFieldParticles[ i +
              j*GravitatingMassFieldParticlesDimension[0] +
              k*GravitatingMassFieldParticlesDimension[0]*
              GravitatingMassFieldParticlesDimension[1]];
    
      /* Clean up if we modified the resolution. */
      
      if (SelfGravity && GravityResolution != 1)
        this->DeleteGravitatingMassFieldParticles();
    
    } // end of (if GravitatingMassFieldParticles != NULL)
    
  return SUCCESS;
}
