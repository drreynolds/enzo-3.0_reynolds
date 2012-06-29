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

  int i, j, k, index, dim, size = 1;
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

  int max_rho_index;
  float max_rho = -1.0;
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) { // nothing
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) { // metallicity
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) { // energy
	index = i + j*GridDimension[0] + k*GridDimension[0]*GridDimension[1];
	if (BaryonField[DensNum][index] > max_rho) {
	  max_rho = BaryonField[DensNum][index];
	  max_rho_index = index;
	}
      }
    }
  }

  /* Update all cells. */

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) { // nothing
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) { // metallicity
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) { // energy

	index = i + j*GridDimension[0] + k*GridDimension[0]*GridDimension[1];

	/* Modify the equation for free-fall collapse with a factor 
	   taking into account the ratio of pressure gradient force 
	   to gravity following Equation 9 from Omukai et al. (2005). */

	new_density = POW((POW(BaryonField[DensNum][index], -0.5) - 
			   (0.5 * FreefallTimeConstant * dtFixed *
			    POW((1 - force_factor[index]), 0.5))), -2.);
	density_ratio = new_density / BaryonField[DensNum][index];

	/* Update enegy. */

	BaryonField[TENum][index] += (gamma_field[index] - 1) * 
	  BaryonField[TENum][index] * FreefallTimeConstant * 
	  POW(BaryonField[DensNum][index], 0.5) * dtFixed;
	if (DualEnergyFormalism) {
	  BaryonField[GENum][index] = BaryonField[TENum][index];
	}

	/* Update density. */

	BaryonField[DensNum][index] = new_density;

	if (index == max_rho_index) {
	  fprintf(stderr, "One-zone collapse: rho[%"ISYM", %"ISYM", %"ISYM"] = %"ESYM" g/cm^3, f = %"FSYM,
		  i, j, k, (BaryonField[DensNum][index] * DensityUnits), 
		  force_factor[index]);
	}

	/* Update species fields. */

	if (MultiSpecies) {
	  BaryonField[DeNum][index] *= density_ratio;
	  BaryonField[HINum][index] *= density_ratio;
	  BaryonField[HIINum][index] *= density_ratio;
	  BaryonField[HeINum][index] *= density_ratio;
	  BaryonField[HeIINum][index] *= density_ratio;
	  BaryonField[HeIIINum][index] *= density_ratio;
	  if (MultiSpecies > 1) {
	    BaryonField[HMNum][index] *= density_ratio;
	    BaryonField[H2INum][index] *= density_ratio;
	    BaryonField[H2IINum][index] *= density_ratio;
	    if (index == max_rho_index) {
	      fprintf(stderr, ", f_H2 = %"ESYM,
		      (BaryonField[H2INum][index] /
		       BaryonField[DensNum][index]));
	    }
	  }
	  if (MultiSpecies > 2) {
	    BaryonField[DINum][index] *= density_ratio;
	    BaryonField[DIINum][index] *= density_ratio;
	    BaryonField[HDINum][index] *= density_ratio;
	  }
	}

	if (MetalFieldPresent) {
	  BaryonField[MetalNum][index] *= density_ratio;
	  if (MultiMetals) {
	    BaryonField[MetalNum+1][index] *= density_ratio;
	    BaryonField[MetalNum+2][index] *= density_ratio;
	  }
	}

	if (index == max_rho_index) {
	  fprintf(stderr, ".\n");
	}

      }
    }
  }

  delete [] gamma_field;
  delete [] pressure;
  delete [] force_factor;

  return SUCCESS;

}
