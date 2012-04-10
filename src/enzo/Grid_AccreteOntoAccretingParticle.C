/***********************************************************************
/
/  GRID CLASS (Calculate the accretion rate and subtract accreted mass 
/              from the grid.)
/
/  written by: Nathan Goldbaum
/  date:       April 2012
/
/  note:       Equation numbers refer to Krumholz McKee & Klein (2004)
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

float bondi_alpha(float x);

int grid::AccreteOntoAccretingParticle(ActiveParticleType_AccretingParticle* ThisParticle, FLOAT AccretionRadius, 
				       float AverageDensity, float &SubtractedMass) {
  /* Return if this doesn't involve us */
  if (MyProcessorNumber != ProcessorNumber) 
    return SUCCESS;
  
  FLOAT *ParticlePosition, LeftCorner[MAX_DIMENSION], RightCorner[MAX_DIMENSION], CellSize, KernelRadius;
  int i, j, k, dim, index, size=1;
  float lambda_c = 0.25*exp(1.5), CellMass, CellVolume = 1;
  
  
  /* Get indices in BaryonField for density, internal energy, thermal energy, velocity */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }
  
  /* Compute cell width and find bottom left and top right corners of grid */

  CellSize = CellWidth[0][0];

  for (dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
    LeftCorner[dim] = CellLeftEdge[dim][0];
    RightCorner[dim] = LeftCorner[dim] + CellSize*(GridDimension[dim]);
    CellVolume*=CellWidth[dim][0];
  }

  ParticlePosition = ThisParticle->ReturnPosition();
  
  /* Check whether the cube that circumscribes the accretion zone intersects with this grid */

  if (LeftCorner[0] > ParticlePosition[0]+AccretionRadius || RightCorner[0] < ParticlePosition[0]-AccretionRadius ||
      LeftCorner[1] > ParticlePosition[1]+AccretionRadius || RightCorner[1] < ParticlePosition[1]-AccretionRadius ||
      LeftCorner[2] > ParticlePosition[2]+AccretionRadius || RightCorner[2] < ParticlePosition[2]-AccretionRadius)
    return SUCCESS;

  // Eqn 13
  if (ThisParticle->BondiHoyleRadius < CellSize/4.0)
    KernelRadius = CellSize/4.0;
  else if (ThisParticle->BondiHoyleRadius < AccretionRadius/2.0)
    KernelRadius = ThisParticle->BondiHoyleRadius;
  else
    KernelRadius = AccretionRadius/2.0;
  
  // Eqn 12
  RhoInfinity = AverageDensity/bondi_alpha(1.2*CellSize/ThisParticle->BondiHoyleRadius);  

  // Eqn 11
  ThisParticle->AccretionRate = (4*pi*RhoInfinity*POW(ThisParticle->BondiHoyleRadius,2)*
				 sqrt(POW(lambda_c*ThisParticle->cInfinity,2) +
				      POW(ThisParticle->vInfinity,2)));

  for (i = 0; i < GridDimension[0]; i++) {
    for (j = 0; j < GridDimension[1]; j++) {
      index = (i*GridDimension[1] + j)*GridDimension[0];
      for (k = 0; k < GridDimension[2]; index++, k++) {
	radius2 = POW((CellLeftEdge[0][i] + 0.5*CellWidth[0][i]) - ParticlePosition[0],2) +
	  POW((CellLeftEdge[1][j] + 0.5*CellWidth[1][j]) - ParticlePosition[1],2) +
	  POW((CellLeftEdge[2][k] + 0.5*CellWidth[2][k]) - ParticlePosition[2],2);   
	if ((AccretionRadius*AccretionRadius) > radius2) {
	  // Subtract mass from this cell
	  CellMass = BaryonField[DensNum][index]*CellVolume;
	  Weight = exp(-radius2/(KernelRadius*KernelRadius))/SumOfWeights;
	  CellMassAccreted =  this->dtFixed * ThisParticle->AccretionRate * Weight;
	  if (CellMassAccreted > 0.25*CellMass) 
	    CellMassAccreted = 0.25*CellMass;
	  SubtractedMass += CellMassAccreted
	  
	  /* The true accretion rate is somewhat less than this due to
	     angular momentum conservation.  Subdivide the cell into
	     NDIV^2 subcells and estimate the reduction assuming
	     ballistic orbits. See the discussion near Eqn 15. */

	  /* Find the components of the momentum vector transverse and
	     perpendicular to the vector connecting the sink and the
	     cell.  Modify the radial component of the momentum vector
	     so that momentum is conserved but leave the transverse
	     component unchanged */
	    
	  // Some conveneint shorthand
	  rhocell = BaryonField[DensNum][index];
	  mcell = rhocell*CellVolume;
	  

	  /* Don't worry about conserving angular momentum if we're
	     accreting no mass from the cell or if we are accreting
	     all of the mass from it.  Note that paccrete and eaccrete
	     are total energy and momentum, not energy an momentum
	     densities. */

	  if ((

	  // Find the radius vector
	  reff[0] = (CellLeftEdge[0][i] + 0.5*CellWidth[0][i]) - ParticlePosition[0];
	  reff[1] = (CellLeftEdge[1][j] + 0.5*CellWidth[1][j]) - ParticlePosition[1];
	  reff[1] = (CellLeftEdge[2][k] + 0.5*CellWidth[2][k]) - ParticlePosition[2];
	  rsqr = reff[0]*reff[0]+reff[1]*reff[1]+reff[2]*reff[2];

	  // Compute the parallel component of the momentum density
	  rdotp = reff[0]*pcell[0]+reff[1]*pcell[1]+rcell[2]*pcell[2];
	  prad[0] = reff[0]*rdotp/rsqr;
	  prad[1] = reff[1]*rdotp/rsqr;
	  prad[2] = reff[2]*rdotp/rsqr;
	    
	  // Compute the transverse component of the momentum density
	  ptrans[0] = pcell[0] - prad[0];
	  ptrans[1] = pcell[1] - prad[1];
	  ptrans[2] = pcell[2] - prad[2];

	  // Compute the new radial component of the momentum density
	  pradnew[0] = prad[0]*rhocell/(rhocell+CellSubtractedMass/CellVolume);
	  pradnew[1] = prad[1]*rhocell/(rhocell+CellSubtractedMass/CellVolume);
	  pradnew[2] = prad[2]*rhocell/(rhocell+CellSubtractedMass/CellVolume);
	
	  // Set the new transverse momentum
	  ptransnew[0] = ptrans[0];
	  ptransnew[1] = ptrans[1];
	  ptransnew[2] = ptrans[2];

	  // Compute the amount of momentum (not momentum density) accreted
	  paccrete[0] = CellVolume*(pcell[0] - pradnew[0] - ptransnew[0]);
	  paccrete[1] = CellVolume*(pcell[1] - pradnew[1] - ptransnew[1]);
	  paccrete[2] = CellVolume*(pcell[2] - pradnew[2] - ptransnew[2]);
	  
	  BaryonField[

	}
      }
    }
  }

  return SUCCESS;
}

/* Routine to return alpha, defined as rho/rho_inf, for a critical
   Bondi accretion solution.  The argument is x = r / rBondiHoyle. 
   Adapted from Orion, courtesy of Mark Krumholz */

float bondi_alpha(float x) {

#define XMIN 0.01
#define XMAX 2.0
#define NTABLE 51

  float alphatable[NTABLE];
  float lamba_c, xtable, xtablep1, alpha_exp;
  int idx;

  /* This is a precomputed table of alpha values.  These correspond to x values
     that run from 0.01 to 2.0 with uniform logarithmic spacing.  The reason for
     this choice of range is that the asymptotic expressions are accurate to
     better than 2% outside this range */

  alphatable = {820.254, 701.882, 600.752, 514.341, 440.497, 377.381, 323.427,                                                  
		277.295, 237.845, 204.1, 175.23, 150.524, 129.377, 111.27, 95.7613,                                             
		82.4745, 71.0869, 61.3237, 52.9498, 45.7644, 39.5963, 34.2989,                                                  
		29.7471, 25.8338, 22.4676, 19.5705, 17.0755, 14.9254, 13.0714,                                                  
		11.4717, 10.0903, 8.89675, 7.86467, 6.97159, 6.19825, 5.52812,                                                  
		4.94699, 4.44279, 4.00497, 3.6246, 3.29395, 3.00637, 2.75612,                                                   
		2.53827, 2.34854, 2.18322, 2.03912, 1.91344, 1.80378, 1.70804,                                                  
		1.62439};

  // A constant that appears in the following formulae.  This hardcoded value is
  // valid for an isothermal gas.
  lambda_c = 0.25*exp(1.5);

  // deal with the off-the-table cases
  if (x < XMIN) 
    bondi_alpha = lambda_x / sqrt(2.*x*x);
  else if (x >= XMAX)
    bondi_alpha = exp(1./x);
  else {
    // we are on the table
    
    idx = floor((NTABLE-1)*log(x/XMIN)/log(XMAX/XMIN));
    xtable = exp(log(XMIN) + idx*log(XMAX/XMIN)/(NTABLE-1));
    xtablep1 = exp(log(XMIN) + (idx+1)*log(XMAX/XMIN)/(NTABLE-1));
    alpha_exp = log(x/xtable) / log(xtablep1/xtable);

    bondi_alpha = alphatable[idx] * POW(alphatable[idx+1]/alphatable[idx],alpha_exp)
  }

#undef NTABLE
#undef XMIN
#undef XMAX

  return bondi_alpha;
}
