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
#include "phys_constants.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "ActiveParticle.h"
#include "phys_constants.h"

#define NO_DEBUG

float bondi_alpha(float x);

int grid::AccreteOntoAccretingParticle(ActiveParticleType* ThisParticle, FLOAT AccretionRadius, 
				       float AverageDensity, float SumOfWeights, float *AccretedMass, 
				       float AccretedMomentum[], bool *SinkIsOnThisGrid, float vInfinity, 
				       float cInfinity, FLOAT BondiHoyleRadius, float *AccretionRate) {

  /* Return if this doesn't involve us */
  if (MyProcessorNumber != ProcessorNumber) 
    return SUCCESS;
  
  /* Check whether the cube that circumscribes the accretion zone intersects with this grid */

  FLOAT *ParticlePosition = ThisParticle->ReturnPosition();

#ifdef DEBUG
  fprintf(stderr,
	  "Left[0] = %"GSYM", Left[1] = %"GSYM", Left[2] = %"GSYM"\n"
	  "Right[0] = %"GSYM", Right[1] = %"GSYM", Right[2] = %"GSYM"\n",
	  ParticlePosition[0]-AccretionRadius,ParticlePosition[1]-AccretionRadius,ParticlePosition[2]-AccretionRadius,
	  ParticlePosition[0]+AccretionRadius,ParticlePosition[1]+AccretionRadius,ParticlePosition[2]+AccretionRadius);
#endif

  if ((GridLeftEdge[0] > ParticlePosition[0]+AccretionRadius) || (GridRightEdge[0] < ParticlePosition[0]-AccretionRadius) ||
      (GridLeftEdge[1] > ParticlePosition[1]+AccretionRadius) || (GridRightEdge[1] < ParticlePosition[1]-AccretionRadius) ||
      (GridLeftEdge[2] > ParticlePosition[2]+AccretionRadius) || (GridRightEdge[2] < ParticlePosition[2]-AccretionRadius))
    return SUCCESS;

  /* Check whether the sink lives on this grid */
  if ((GridLeftEdge[0] < ParticlePosition[0]) && (GridRightEdge[0] > ParticlePosition[0]) &&
      (GridLeftEdge[1] < ParticlePosition[1]) && (GridRightEdge[1] > ParticlePosition[1]) &&
      (GridLeftEdge[2] < ParticlePosition[2]) && (GridRightEdge[2] > ParticlePosition[2])) {
    *SinkIsOnThisGrid = true;
  }

  FLOAT CellSize, KernelRadius, radius2;

  int i, j, k, dim, index;
  float lambda_c = 0.25*exp(1.5), CellMass, CellVolume = 1., SmallRhoFac = 10., 
    SmallEFac = 10.;

  float RhoInfinity, vsink[3], vgas[3], mcell, etot, eint, Weight, maccreted, 
    rhocell, pcell[3], paccrete[3], eaccrete, mnew, rhonew, reff[3], rsqr, 
    rdotp, prad[3], ptrans[3], pradnew[3], ptransnew[3], eintnew, 
    pnew[3], kenew;

  int isub, jsub, ksub, excluded, NDIV = 8;
  
  float xdist, ydist, zdist, dist, jsp[3], jspsqr, esp, rmin, dxmin, ACCRETERADMIN = 2.0, huge = 1.0e30;

  for (i = 0; i < 3; i++) {
    vsink[i] = 0;
    vgas[i] = 0;
    pcell[i] = 0;
    paccrete[i] = 0;
    reff[i] = 0;
    prad[i] = 0;
    ptrans[i] = 0;
    pradnew[i] = 0;
    ptransnew[i] = 0;
    pnew[i] = 0;
  }

  /* Get indices in BaryonField for density, internal energy, thermal energy, velocity */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }
  
  /* Compute cell width and find bottom left and top right corners of grid */

  CellSize = CellWidth[0][0];

  for (dim = 0; dim < GridRank; dim++) {
    CellVolume*=CellWidth[dim][0];
  }
  
  // Eqn 13
  if (BondiHoyleRadius < CellSize/4.0)
    KernelRadius = CellSize/4.0;
  else if (BondiHoyleRadius < AccretionRadius/2.0)
    KernelRadius = BondiHoyleRadius;
  else
    KernelRadius = AccretionRadius/2.0;
  
  // Eqn 12
  RhoInfinity = AverageDensity/bondi_alpha(1.2*CellSize/BondiHoyleRadius);  

  // Eqn 11
  *AccretionRate = (4*pi*RhoInfinity*POW(BondiHoyleRadius,2)*
		   sqrt(POW(lambda_c*cInfinity,2) + POW(vInfinity,2)));

  vsink[0] = ThisParticle->ReturnVelocity()[0];
  vsink[1] = ThisParticle->ReturnVelocity()[1];
  vsink[2] = ThisParticle->ReturnVelocity()[2];

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	radius2 = 
	  POW((CellLeftEdge[0][i] + 0.5*CellWidth[0][i]) - ParticlePosition[0],2) +
	  POW((CellLeftEdge[1][j] + 0.5*CellWidth[1][j]) - ParticlePosition[1],2) +
	  POW((CellLeftEdge[2][k] + 0.5*CellWidth[2][k]) - ParticlePosition[2],2);   
	if ((AccretionRadius*AccretionRadius) > radius2) {
#ifdef DEBUG
	  fprintf(stderr,
		  "CellLeftEdge[0][i] = %"GSYM", CellRightEdge[0][i] =%"GSYM"\n"
		  "CellLeftEdge[1][j] = %"GSYM", CellRightEdge[1][j] =%"GSYM"\n"
		  "CellLeftEdge[2][k] = %"GSYM", CellRightEdge[2][k] =%"GSYM"\n",
		  CellLeftEdge[0][i],CellLeftEdge[0][i]+CellWidth[0][i],CellLeftEdge[1][j],CellLeftEdge[1][j]+CellWidth[1][j],
		  CellLeftEdge[2][k],CellLeftEdge[2][k]+CellWidth[2][k]);
#endif
	  // useful shorthand
	  vgas[0] = BaryonField[Vel1Num][index];
	  vgas[1] = BaryonField[Vel2Num][index];
	  vgas[2] = BaryonField[Vel3Num][index];
	  rhocell = BaryonField[DensNum][index];
	  mcell = rhocell*CellVolume;
	  // These are momentum densities in the frame of the sink.
	  pcell[0] = rhocell*vgas[0] - vsink[0]*rhocell;
	  pcell[1] = rhocell*vgas[1] - vsink[1]*rhocell;
	  pcell[2] = rhocell*vgas[2] - vsink[2]*rhocell;
	  	  
	  // TE and GE are stored per unit mass
	  if (HydroMethod == 0) { // PPM
	    etot = rhocell*BaryonField[TENum][index];
	    if (DualEnergyFormalism)
	      eint = rhocell*BaryonField[GENum][index];
	    else
	      eint = etot - 0.5*rhocell*(vgas[0]*vgas[0] + vgas[1]*vgas[1] + vgas[2]*vgas[2]);
	  } else {  // Zeus hydro (total energy is really internal energy)
	    eint = rhocell*BaryonField[TENum][index];
	    etot = eint + 0.5*rhocell*(vgas[0]*vgas[0] + vgas[1]*vgas[1] + vgas[2]*vgas[2]);
	  }
	  
	  // Calculate mass we need to subtract from this cell
	  Weight = exp(-radius2/(KernelRadius*KernelRadius))/SumOfWeights;
	  maccreted =  this->dtFixed * (*AccretionRate) * Weight;
	  if (maccreted > 0.1*mcell) 
	    maccreted = 0.1*mcell;
	  
	  /* The true accretion rate is somewhat less than this due to
	     angular momentum conservation.  Subdivide the cell into
	     NDIV^2 subcells and estimate the reduction assuming
	     ballistic orbits. See the discussion near Eqn 15. */
	  
	  excluded = 0;
	  for (ksub = 0; ksub < NDIV-1; ksub++) {
	    zdist = CellLeftEdge[2][k] + CellWidth[2][k]*(float(ksub)+0.5)/NDIV - ParticlePosition[2];
	    for (jsub = 0; jsub < NDIV-1; jsub++) {
	      ydist = CellLeftEdge[1][j] + CellWidth[1][j]*(float(jsub)+0.5)/NDIV - ParticlePosition[1];
	      for (isub = 0; isub < NDIV-1; isub++) {
		xdist = CellLeftEdge[0][i] + CellWidth[0][i]*(float(jsub)+0.5)/NDIV - ParticlePosition[1];

		dist = sqrt(xdist*xdist+ydist*ydist+zdist*zdist);
		if (dist == 0.0)
		  dist = CellWidth[0][0]/huge;

		// Compute specific angular momentum
		jsp[0] = ydist*(vgas[2] - vsink[2]) -
		  zdist*(vgas[1] - vsink[1]);
		jsp[1] = zdist*(vgas[0] - vsink[0]) - 
		  xdist*(vgas[2] - vsink[2]);
		jsp[2] = xdist*(vgas[1] - vsink[1]) -
		  ydist*(vgas[0] - vsink[0]);
		
		jspsqr = jsp[0]*jsp[0]+jsp[1]*jsp[1]+jsp[2]*jsp[2];

		// Compute specific kinetic + gravitational energy
		esp = (POW((vgas[0] - vsink[0]),2) + 
		       POW((vgas[1] - vsink[1]),2) +
		       POW((vgas[2] - vsink[2]),2)) / 2.0 - GravConst * mcell/dist;
		
		// Compute distance of closest approach
		if (esp > 0.0)
		  rmin = huge*CellWidth[0][0];
		else
		  rmin = -GravConst*mcell/(2.0*esp) *
		    (1.0 - sqrt(1.0 + 2.0*jspsqr*esp/POW(GravConst*mcell,2)));
		
		dxmin = rmin / CellWidth[0][0];
		if (dxmin > ACCRETERADMIN)
		  excluded+=1;

	      } // ksub
	    } // jsub
	  } // ksub

	  // Scale down maccrete
	  maccreted = maccreted/POW(NDIV,3) *
	    (POW(NDIV,3)-excluded);
	    
	  /* Don't worry about conserving angular momentum if we're
	     accreting no mass from the cell or if we are accreting
	     all of the mass from it.  Note that paccrete and eaccrete
	     are total energy and momentum, not energy and momentum
	     densities or specific energy and momentum */
	  
	  if ((maccreted == 0) || (mcell-maccreted < 2.0*SmallRhoFac*SmallRho*CellVolume)) {
	    paccrete[0] = pcell[0]*CellVolume*(maccreted/mcell);
	    paccrete[1] = pcell[1]*CellVolume*(maccreted/mcell);
	    paccrete[2] = pcell[2]*CellVolume*(maccreted/mcell);
	    
	    eaccrete = etot*CellVolume*(maccreted/mcell);
	  } 
	  
	  /* Find the components of the momentum vector transverse and
	     perpendicular to the vector connecting the sink and the
	     cell.  Modify the radial component of the momentum vector
	     so that momentum is conserved but leave the transverse
	     component unchanged */
	  	  
	  else {

	    // Keep cell mass well above density floor
	    if ((mcell - maccreted)/CellVolume > SmallRhoFac*SmallRho)
	      mnew = mcell - maccreted;
	    else {
	      mnew = SmallRhoFac*SmallRho*CellVolume;
	      maccreted = mcell - mnew;
	    }

	    rhonew = mnew/CellVolume;

	    // Find the radius vector
	    reff[0] = (CellLeftEdge[0][i] + 0.5*CellWidth[0][i]) - ParticlePosition[0];
	    reff[1] = (CellLeftEdge[1][j] + 0.5*CellWidth[1][j]) - ParticlePosition[1];
	    reff[2] = (CellLeftEdge[2][k] + 0.5*CellWidth[2][k]) - ParticlePosition[2];
	    rsqr = reff[0]*reff[0]+reff[1]*reff[1]+reff[2]*reff[2];
	    
	    // Prevent a floating point error if close to central cell
	    if (rsqr <= POW(1e-7*CellWidth[0][0],2)) {
	      reff[0] = 0.0;
	      reff[1] = 0.0;
	      reff[2] = 0.0;
	      rsqr = 1.0;
	    }
	    
	    // Compute the parallel component of the momentum density
	    rdotp = reff[0]*pcell[0]+reff[1]*pcell[1]+reff[2]*pcell[2];
	    prad[0] = reff[0]*rdotp/rsqr;
	    prad[1] = reff[1]*rdotp/rsqr;
	    prad[2] = reff[2]*rdotp/rsqr;
	    
	    // Compute the transverse component of the momentum density
	    ptrans[0] = pcell[0] - prad[0];
	    ptrans[1] = pcell[1] - prad[1];
	    ptrans[2] = pcell[2] - prad[2];
	    
	    // Compute the new radial component of the momentum density
	    pradnew[0] = prad[0]*mnew/mcell;
	    pradnew[1] = prad[1]*mnew/mcell;
	    pradnew[2] = prad[2]*mnew/mcell;
	    
	    // Set the new transverse momentum
	    ptransnew[0] = ptrans[0];
	    ptransnew[1] = ptrans[1];
	    ptransnew[2] = ptrans[2];
	    
	    // Compute the amount of momentum (not momentum density) accreted
	    paccrete[0] = CellVolume*(pcell[0] - pradnew[0] - ptransnew[0]);
	    paccrete[1] = CellVolume*(pcell[1] - pradnew[1] - ptransnew[1]);
	    paccrete[2] = CellVolume*(pcell[2] - pradnew[2] - ptransnew[2]);
	    
	    // Compute new total internal energy (total, not density, not specific)
	    eintnew = eint * rhonew/rhocell * CellVolume;
	    
	    /* Compute the new momentum densities.  Note that we do
	       not use pcell here because we need to do this
	       calculation in the grid frame, not the sink frame. */

	    pnew[0] = rhocell*vgas[0] - paccrete[0]/CellVolume;
	    pnew[1] = rhocell*vgas[1] - paccrete[1]/CellVolume;
	    pnew[2] = rhocell*vgas[2] - paccrete[2]/CellVolume;
 
	    // Compute new total kinetic energy (not density)
	    kenew = (pnew[0]*pnew[0] + pnew[1]*pnew[1] + pnew[2]*pnew[2]) / 
	      (2.0 * rhonew) * CellVolume;
	  
	    // Compute the amount of energy (not energy density) to be accreted
	    eaccrete = etot*CellVolume - (eintnew + kenew);
	    
	    // This is actually a density since particle masses are stored
	    // in density units.
	    *AccretedMass += maccreted/CellVolume;
	    AccretedMomentum[0] += paccrete[0];
	    AccretedMomentum[1] += paccrete[1];
	    AccretedMomentum[2] += paccrete[2];

	    // Update the state variables
	    BaryonField[DensNum][index] -= maccreted/CellVolume;
	    
	    if (HydroMethod == 0) { // PPM
	      BaryonField[TENum][index] -= eaccrete/mnew;
	      if (DualEnergyFormalism)
		BaryonField[GENum][index] = eintnew/mnew;
	    } else // Zeus
	      BaryonField[TENum][index] = eintnew/mnew;

	    BaryonField[Vel1Num][index] = BaryonField[Vel1Num][index]*mcell/mnew - paccrete[0]/mnew;
	    BaryonField[Vel2Num][index] = BaryonField[Vel2Num][index]*mcell/mnew - paccrete[1]/mnew;
	    BaryonField[Vel3Num][index] = BaryonField[Vel3Num][index]*mcell/mnew - paccrete[2]/mnew;

	    // Check if mass or energy is too small, correct if necessary
	    if (BaryonField[DensNum][index] < SmallRhoFac*SmallRho) {
	      BaryonField[DensNum][index] = SmallRhoFac*SmallRho;
	      BaryonField[Vel1Num][index] = vgas[0];
	      BaryonField[Vel1Num][index] = vgas[1];
	      BaryonField[Vel1Num][index] = vgas[2];
	    }
	    
	    if (HydroMethod == 0) { // PPM
	      if (BaryonField[TENum][index] - 
		  0.5 * (BaryonField[Vel1Num][index]*BaryonField[Vel1Num][index] + 
			 BaryonField[Vel2Num][index]*BaryonField[Vel2Num][index] +
			 BaryonField[Vel3Num][index]*BaryonField[Vel3Num][index]) 
		  < SmallEFac*SmallEint) {
		BaryonField[TENum][index] = SmallEFac*SmallEint + 
		  0.5 * (BaryonField[Vel1Num][index]*BaryonField[Vel1Num][index] +
			 BaryonField[Vel2Num][index]*BaryonField[Vel2Num][index] +
			 BaryonField[Vel3Num][index]*BaryonField[Vel3Num][index]);
		if (DualEnergyFormalism)
		  BaryonField[GENum][index] = SmallEFac*SmallEint;
	      }
	      else  // Zeus
		if (BaryonField[TENum][index] < SmallEFac*SmallEint)
		  BaryonField[TENum][index] = SmallEFac*SmallEint;
	    }
	  }
	}
      }
    }
  }

  return SUCCESS;
}

#undef DEBUG

/* Routine to return alpha, defined as rho/rho_inf, for a critical
   Bondi accretion solution.  The argument is x = r / rBondiHoyle. 
   Adapted from Orion, courtesy of Mark Krumholz */

float bondi_alpha(float x) {

#define XMIN 0.01
#define XMAX 2.0
#define NTABLE 51

  float lambda_c, xtable, xtablep1, alpha_exp;
  int idx;

  /* This is a precomputed table of alpha values.  These correspond to x values
     that run from 0.01 to 2.0 with uniform logarithmic spacing.  The reason for
     this choice of range is that the asymptotic expressions are accurate to
     better than 2% outside this range */

  float alphatable[NTABLE] = {820.254, 701.882, 600.752, 514.341, 440.497, 377.381, 323.427,
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
    return lambda_c / sqrt(2.*x*x);
  else if (x >= XMAX)
    return exp(1./x);
  else {
    // we are on the table
    
    idx = floor((NTABLE-1)*log(x/XMIN)/log(XMAX/XMIN));
    xtable = exp(log(XMIN) + idx*log(XMAX/XMIN)/(NTABLE-1));
    xtablep1 = exp(log(XMIN) + (idx+1)*log(XMAX/XMIN)/(NTABLE-1));
    alpha_exp = log(x/xtable) / log(xtablep1/xtable);

    return alphatable[idx] * POW(alphatable[idx+1]/alphatable[idx],alpha_exp);
  }

#undef NTABLE
#undef XMIN
#undef XMAX

}
