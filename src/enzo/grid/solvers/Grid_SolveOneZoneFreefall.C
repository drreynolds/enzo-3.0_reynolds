/***********************************************************************
/
/  GRID CLASS (SOLVE THE ANALYTICAL SOLUTION FOR FREE-FALL COLLAPSE)
/
/  written by: Britton Smith
/  date:       October, 2010
/  modified1:  
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
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
#include "fortran.def"
#include "CosmologyParameters.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);

int grid::SolveOneZoneFreefall()
{

  /* Return if this doesn't concern us. */
  
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  this->DebugCheck("SolveRadiativeCooling");

  /* Declarations */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
    
  /* Find fields: density, total energy, velocity1-3. */

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
        ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  /* Find Multi-species fields. */

  if (MultiSpecies)
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
            ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }

  /* Compute size of the current grid. */

  int i, f, dim, size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  /* Calculate units. */

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  /* Get gamma field for updating densities and energies. */

  float *gamma_field = new float[size];
  if (this->ComputeGammaField(gamma_field) == FAIL) {
    ENZO_FAIL("Error in grid->ComputeGammaField.\n");
  }

  /* Compute pressure field. */
  float *pressure = new float[size];
  if (this->ComputePressure(Time, pressure) == FAIL) {
    ENZO_FAIL("Error in grid->ComputePressure.\n");
  }

  /* Compute ratio of pressure gradient force to graviational force.
     Equation 9 of Omukai et al (2005). */

  float *force_factor = new float[size];
  if (this->ComputeOneZoneCollapseFactor(force_factor) == FAIL) {
    ENZO_FAIL("Error in ComputeOneZoneCollapseFactor.\n");
  }

  /* Metal cooling codes. */

  int MetalNum = 0, SNColourNum = 0;
  int MetalFieldPresent = FALSE;

  // First see if there's a metal field (so we can conserve species in
  // the solver)
  MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);
  SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields);
  MetalFieldPresent = (MetalNum != -1 || SNColourNum != -1);

  /* Calculate new density and energy. */

  float new_density, density_ratio;
  float FreefallTimeConstant = POW(((32 * GravitationalConstant) / 
				    (3 * pi)), 0.5);

  /* Update density and pressure history. */

  if (CollapseHistory[0] == NULL) {
    CollapseHistory[0] = new float*[2];
    CollapseHistory[0][0] = new float[size]; // density
    CollapseHistory[0][1] = new float[size]; // pressure
  }
  else {
    if (CollapseHistory[1] == NULL) {
      CollapseHistory[1] = new float*[2];
      CollapseHistory[1][0] = new float[size]; // density
      CollapseHistory[1][1] = new float[size]; // pressure
    }
    // move t-1 values into t-2
    for (i = 0; i < size; i++) {
      CollapseHistory[1][0][i] = CollapseHistory[0][0][i];
      CollapseHistory[1][1][i] = CollapseHistory[0][1][i];
    }
  }

  // move current values into t-1
  for (i = 0; i < size; i++) {
    CollapseHistory[0][0][i] = BaryonField[DensNum][i];
    CollapseHistory[0][1][i] = pressure[i];
  }

  /* Update all cells. */

  for (i = 0;i < size;i++) {

    /* Modify the equation for free-fall collapse with a factor 
       taking into account the ratio of pressure gradient force 
       to gravity following Equation 9 from Omukai et al. (2005). */

    new_density = POW((POW(BaryonField[DensNum][i], -0.5) - 
		       (0.5 * FreefallTimeConstant * dtFixed *
			POW((1 - force_factor[i]), 0.5))), -2.);
    density_ratio = new_density / BaryonField[DensNum][i];

    /* Update enegy. */

    BaryonField[TENum][i] += (gamma_field[i] - 1) * BaryonField[TENum][i] * 
      FreefallTimeConstant * POW(BaryonField[DensNum][i], 0.5) * dtFixed;
    if (DualEnergyFormalism) {
      BaryonField[GENum][i] = BaryonField[TENum][i];
    }

    /* Update density. */

    BaryonField[DensNum][i] = new_density;

    if (i == 0) {
      fprintf(stderr, "One-zone collapse: rho[0] = %"ESYM" g/cm^3, f = %"FSYM,
	      (BaryonField[DensNum][i] * DensityUnits), force_factor[i]);
    }

    /* Update species fields. */

    if (MultiSpecies) {
      BaryonField[DeNum][i] *= density_ratio;
      BaryonField[HINum][i] *= density_ratio;
      BaryonField[HIINum][i] *= density_ratio;
      BaryonField[HeINum][i] *= density_ratio;
      BaryonField[HeIINum][i] *= density_ratio;
      BaryonField[HeIIINum][i] *= density_ratio;
      if (MultiSpecies > 1) {
	BaryonField[HMNum][i] *= density_ratio;
	BaryonField[H2INum][i] *= density_ratio;
	BaryonField[H2IINum][i] *= density_ratio;
	if (i == 0) {
	  fprintf(stderr, ", f_H2 = %"ESYM,
		  (BaryonField[H2INum][i] /
		   BaryonField[DensNum][i]));
	}
      }
      if (MultiSpecies > 2) {
	BaryonField[DINum][i] *= density_ratio;
	BaryonField[DIINum][i] *= density_ratio;
	BaryonField[HDINum][i] *= density_ratio;
      }
    }

    if (MetalFieldPresent) {
      BaryonField[MetalNum][i] *= density_ratio;
      if (MultiMetals) {
	BaryonField[MetalNum+1][i] *= density_ratio;
	BaryonField[MetalNum+2][i] *= density_ratio;
      }
    }

    if (i == 0) {
      fprintf(stderr, ".\n");
    }

  }

  delete [] gamma_field;
  delete [] pressure;
  delete [] force_factor;

  return SUCCESS;

}
