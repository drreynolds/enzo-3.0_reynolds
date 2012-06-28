/***********************************************************************
/
/  GRID CLASS (COMPUTE THE RATIO OF PRESSURE TO GRAVITY FOR ONE-ZONE
/              COLLAPSE)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE: Computes an approximate ratio of pressure to gravitational 
/           force as modifier to the freefall collapse solution.
/           This follows equations 6-9 of Omukai et al (2005).
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
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
 
int grid::ComputeOneZoneCollapseFactor(float *force_factor)
{
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;  

  /* Compute the size of the fields. */
 
  int i, t, size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  if (!TestProblemData.OneZoneFreefallAdjustCollapse) {
    for (i = 0; i < size; i++) {
      force_factor[i] = 0.0;
    }
    return SUCCESS;
  }
 
  /* Check for density and pressure history. */
  for (t = 0; t < 2; t++) {
    if (CollapseHistory[t] == NULL) {
      for (i = 0; i < size; i++) {
	force_factor[i] = 0.0;
      }
      return SUCCESS;
    }
  }

  /* Find Density, if possible. */
 
  int DensNum;  
  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0)
    ENZO_FAIL("Cannot find density.");

  float gamma_eff;
  for (i = 0; i < size; i++) {

    /* Calculate the effective adiabatic index, dlog(p)/dlog(rho). */
    gamma_eff = log10(CollapseHistory[0][1][i] / CollapseHistory[1][1][i]) /
      log10(CollapseHistory[0][0][i] / CollapseHistory[1][0][i]);

    if (gamma_eff < 0.83) {
      force_factor[i] = 0.0;
    }
    else if (gamma_eff < 1.0) {
      force_factor[i] = 0.6 + 2.5 * (gamma_eff - 1) -
	6.0 * POW((gamma_eff - 1.0), 2.);
    }
    else if (gamma_eff < (4./3.)) {
      force_factor[i] = 1.0 + 0.2 * (gamma_eff - (4./3.)) -
    	2.9 * POW((gamma_eff - (4./3.)), 2.);
    }
    force_factor[i] = max(force_factor[i], 0.0);
    force_factor[i] = min(force_factor[i], 0.95);

  }
 
  return SUCCESS;
}
