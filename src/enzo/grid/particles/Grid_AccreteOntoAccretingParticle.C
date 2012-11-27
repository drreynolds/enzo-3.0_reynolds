/***********************************************************************
/
/  (Calculate the accretion rate and subtract accreted mass from the 
/   grid.)
/
/  written by: Nathan Goldbaum
/  date:       April 2012
/
/  note:       Equation numbers refer to Krumholz McKee & Klein (2004)
************************************************************************/

#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "units.h"
#include "Fluxes.h"
#include "GridList.h"
#include "phys_constants.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "ActiveParticle.h"
#include "phys_constants.h"


#define NO_DEBUG_AP

float bondi_alpha(float x);

int grid::AccreteOntoAccretingParticle(ActiveParticleType** ThisParticle,FLOAT AccretionRadius,
				       float* AccretionRate){

  /* Return if this doesn't involve us */
  if (MyProcessorNumber != ProcessorNumber) 
    return SUCCESS;

  /* Check whether the cube that circumscribes the accretion zone intersects with this grid */

  FLOAT *ParticlePosition = (*ThisParticle)->ReturnPosition();
  
  if ((GridLeftEdge[0] > ParticlePosition[0]+AccretionRadius) || (GridRightEdge[0] < ParticlePosition[0]-AccretionRadius) ||
      (GridLeftEdge[1] > ParticlePosition[1]+AccretionRadius) || (GridRightEdge[1] < ParticlePosition[1]-AccretionRadius) ||
      (GridLeftEdge[2] > ParticlePosition[2]+AccretionRadius) || (GridRightEdge[2] < ParticlePosition[2]-AccretionRadius))
    return SUCCESS;

  FLOAT CellSize, KernelRadius, radius2;

  int i, j, k, dim, index;
  float lambda_c = 0.25*exp(1.5), CellMass, CellVolume = 1., SmallRhoFac = 10., 
    SmallEFac = 10., AccretedMass = 0, AccretedMomentum[3], 
    RhoInfinity, vsink[3], vgas[3], mcell, etot, eint, ke, Weight, maccreted, 
    rhocell, pcell[3], paccrete[3], etotnew, mnew, rhonew, reff[3], rsqr, 
    rdotp, prad[3], ptrans[3], pradnew[3], ptransnew[3], eintnew, 
    pnew[3], kenew, xdist, ydist, zdist, dist, jsp[3], jspsqr, esp, rmin, 
    dxmin, huge = 1.0e30, WeightedSum = 0, 
    SumOfWeights = 0, AverageDensity = 0;

  int isub, jsub, ksub, excluded, NDIV = 8, NumberOfCells=0;

  for (i = 0; i < 3; i++) {
    AccretedMomentum[i] = 0;
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
    jsp[i] = 0;
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

  /* Get grid info and calculate index of the sink host cell */
  int size = this->GetGridSize(), maxexcluded=0;
  int *nexcluded = new int[size]();
  int cindex = (GridEndIndex[0] - GridStartIndex[0])/2 + GridStartIndex[0];
  int cgindex = GRIDINDEX_NOGHOST(cindex,cindex,cindex);

  /* Find the Bondi-Hoyle radius */
  float vInfinity, cInfinity, CellTemperature;
  float velx = BaryonField[Vel1Num][cgindex];
  float vely = BaryonField[Vel2Num][cgindex];
  float velz = BaryonField[Vel3Num][cgindex];
  FLOAT BondiHoyleRadius;
  float *Temperature = new float[size]();
  float msink = (*ThisParticle)->ReturnMass()*CellVolume;

  this->ComputeTemperatureField(Temperature);

  vsink[0] = (*ThisParticle)->ReturnVelocity()[0];
  vsink[1] = (*ThisParticle)->ReturnVelocity()[1];
  vsink[2] = (*ThisParticle)->ReturnVelocity()[2];

  vInfinity = sqrt(pow(vsink[0] - velx,2) + 
		   pow(vsink[1] - vely,2) + 
		   pow(vsink[2] - velz,2));

  CellTemperature = (JeansRefinementColdTemperature > 0) ? JeansRefinementColdTemperature : Temperature[cgindex];
  cInfinity = sqrt(Gamma*kboltz*CellTemperature/(Mu*mh))/GlobalLengthUnits*GlobalTimeUnits;
  BondiHoyleRadius = GravitationalConstant*msink/
    (pow(vInfinity,2) + pow(cInfinity,2));

  delete [] Temperature;

  // Eqn 13
  if (BondiHoyleRadius < CellSize/4.0)
    KernelRadius = CellSize/4.0;
  else if (BondiHoyleRadius < AccretionRadius/2.0)
    KernelRadius = BondiHoyleRadius;
  else
    KernelRadius = AccretionRadius/2.0;

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
     for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
       index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
       for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	 radius2 = 
	   POW((CellLeftEdge[0][i] + 0.5*CellWidth[0][i]) - ParticlePosition[0],2) +
	   POW((CellLeftEdge[1][j] + 0.5*CellWidth[1][j]) - ParticlePosition[1],2) +
	   POW((CellLeftEdge[2][k] + 0.5*CellWidth[2][k]) - ParticlePosition[2],2);   
	 if ((AccretionRadius*AccretionRadius) > radius2) {
	   WeightedSum += BaryonField[DensNum][index]*exp(-radius2/(KernelRadius*KernelRadius)); 
	   SumOfWeights += exp(-radius2/(KernelRadius*KernelRadius));
	   NumberOfCells++;
	   vgas[0] = BaryonField[Vel1Num][index];
	   vgas[1] = BaryonField[Vel2Num][index];
	   vgas[2] = BaryonField[Vel3Num][index];
	   
	   /* The true accretion rate is somewhat less than this due to
	      angular momentum conservation.  Subdivide the cell into
	      NDIV^2 subcells and estimate the reduction assuming
	      ballistic orbits. See the discussion near Eqn 15. */	  
	   for (ksub = 0; ksub < NDIV-1; ksub++) {
	     zdist = CellLeftEdge[2][k] + CellWidth[2][k]*(float(ksub)+0.5)/NDIV - ParticlePosition[2];
	     for (jsub = 0; jsub < NDIV-1; jsub++) {
	       ydist = CellLeftEdge[1][j] + CellWidth[1][j]*(float(jsub)+0.5)/NDIV - ParticlePosition[1];
	       for (isub = 0; isub < NDIV-1; isub++) {
		 xdist = CellLeftEdge[0][i] + CellWidth[0][i]*(float(jsub)+0.5)/NDIV - ParticlePosition[0];
		 
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
			POW((vgas[2] - vsink[2]),2)) / 2.0 - GravitationalConstant * msink/dist;
		 
		 // Compute distance of closest approach
		 if (esp > 0.0)
		   rmin = huge*CellWidth[0][0];
		 else
		   rmin = -GravitationalConstant*msink/(2.0*esp) *
		     (1.0 - sqrt(1.0 + 2.0*jspsqr*esp/POW(GravitationalConstant*msink,2)));
		 
		 dxmin = rmin / CellWidth[0][0];
		 if (dxmin >= 0.25)
		   nexcluded[index]+=1;
		 
	       } // ksub
	     } // jsub
	   } // ksub
	   
	   if (abs(i-cindex) <= 1 && abs(j-cindex) <= 1 && abs(k-cindex) <= 1)
	     maxexcluded = max(nexcluded[index],maxexcluded);
	   
	 }
       }
     }
  }
  
  // Correct the central cell
  if (nexcluded[cgindex] > 0)
    if (KernelRadius / CellWidth[0][0] >= 0.25)
      nexcluded[cgindex] = maxexcluded;
    else
      nexcluded[cgindex] = 0;
    
  AverageDensity = WeightedSum/SumOfWeights;
  
  // Eqn 12
  RhoInfinity = AverageDensity/bondi_alpha(1.2*CellSize / BondiHoyleRadius);  

  // Eqn 11
  *AccretionRate = (4*pi*RhoInfinity*POW(BondiHoyleRadius,2)*
		    sqrt(POW(lambda_c*cInfinity,2) + POW(vInfinity,2)));

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	radius2 = 
	  POW((CellLeftEdge[0][i] + 0.5*CellWidth[0][i]) - ParticlePosition[0],2) +
	  POW((CellLeftEdge[1][j] + 0.5*CellWidth[1][j]) - ParticlePosition[1],2) +
	  POW((CellLeftEdge[2][k] + 0.5*CellWidth[2][k]) - ParticlePosition[2],2);   
	if ((AccretionRadius*AccretionRadius) > radius2) {
	  // useful shorthand
	  vgas[0] = BaryonField[Vel1Num][index];
	  vgas[1] = BaryonField[Vel2Num][index];
	  vgas[2] = BaryonField[Vel3Num][index];
	  rhocell = BaryonField[DensNum][index];
	  mcell = rhocell*CellVolume;
	  // These are momenta in the frame of the sink.
	  pcell[0] = mcell*vgas[0] - mcell*vsink[0];
	  pcell[1] = mcell*vgas[1] - mcell*vsink[1];
	  pcell[2] = mcell*vgas[2] - mcell*vsink[2];
	  	  
	  // TE and GE are stored per unit mass
	  if (HydroMethod == 0) { // PPM
	    etot = mcell*BaryonField[TENum][index];
	    if (DualEnergyFormalism)
	      eint = mcell*BaryonField[GENum][index];
	    else
	      eint = etot - 0.5*mcell*(vgas[0]*vgas[0] + vgas[1]*vgas[1] + vgas[2]*vgas[2]);
	  } else {  // Zeus hydro (total energy is really internal energy)
	    eint = mcell*BaryonField[TENum][index];
	    etot = eint + 0.5*mcell*(vgas[0]*vgas[0] + vgas[1]*vgas[1] + vgas[2]*vgas[2]);
	  }
	  
	  ke = 0.5*mcell*(vgas[0]*vgas[0] + vgas[1]*vgas[1] + vgas[2]*vgas[2]);

	  // Calculate mass we need to subtract from this cell
	  Weight = exp(-radius2/(KernelRadius*KernelRadius))/SumOfWeights;
	  maccreted =  this->dtFixed * (*AccretionRate) * Weight;
	  if (maccreted > 0.25*mcell) 
	    maccreted = 0.25*mcell;
	  
	  // Scale down maccreted
	  maccreted = maccreted/POW(NDIV,3) *
	    (POW(NDIV,3)-nexcluded[index]);

	  /* Don't worry about conserving angular momentum if we're
	     accreting no mass from the cell or if we are accreting
	     all of the mass from it.  */
	  
	  if ((maccreted == 0) || (mcell-maccreted < 2.0*SmallRhoFac*SmallRho*CellVolume)) {
	    paccrete[0] = pcell[0]*CellVolume*(maccreted/mcell);
	    paccrete[1] = pcell[1]*CellVolume*(maccreted/mcell);
	    paccrete[2] = pcell[2]*CellVolume*(maccreted/mcell);
	    
	    etotnew = etot * (1.0 - maccreted/mcell);
	    eintnew = eint * (1.0 - maccreted/mcell);
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
	    
	    // Compute the parallel component of the momentum 
	    rdotp = reff[0]*pcell[0]+reff[1]*pcell[1]+reff[2]*pcell[2];
	    prad[0] = reff[0]*rdotp/rsqr;
	    prad[1] = reff[1]*rdotp/rsqr;
	    prad[2] = reff[2]*rdotp/rsqr;
	    
	    // Compute the transverse component of the momentum 
	    ptrans[0] = pcell[0] - prad[0];
	    ptrans[1] = pcell[1] - prad[1];
	    ptrans[2] = pcell[2] - prad[2];
	    
	    // Compute the new radial component of the momentum 
	    pradnew[0] = prad[0]*mnew/mcell;
	    pradnew[1] = prad[1]*mnew/mcell;
	    pradnew[2] = prad[2]*mnew/mcell;
	    
	    // Set the new transverse momentum
	    ptransnew[0] = ptrans[0];
	    ptransnew[1] = ptrans[1];
	    ptransnew[2] = ptrans[2];
	    
	    // Compute the amount of momentum accreted
	    paccrete[0] = pcell[0] - pradnew[0] - ptransnew[0];
	    paccrete[1] = pcell[1] - pradnew[1] - ptransnew[1];
	    paccrete[2] = pcell[2] - pradnew[2] - ptransnew[2];
	    
	    // Compute new total internal energy. By construction,
	    // this keeps the specific internal energy constant after
	    // accretion
	    eintnew = eint * (1.0 - maccreted/mcell);
  
	    /* Compute the new momentum of the cell.  Note that we do
	       not use pcell here because we need to do this
	       calculation in the grid frame, not the sink frame. */
	    pnew[0] = mcell*vgas[0] - paccrete[0];
	    pnew[1] = mcell*vgas[1] - paccrete[1];
	    pnew[2] = mcell*vgas[2] - paccrete[2];
 
	    // Compute new total kinetic energy (not density)
	    kenew = (pnew[0]*pnew[0] + pnew[1]*pnew[1] + pnew[2]*pnew[2]) / 
	      (2.0 * mnew);
	  
	    // Compute the new total energy
	    etotnew = eintnew + kenew;
	    	    
	    // This is actually a density since particle masses are stored
	    // in density units.
	    AccretedMass += maccreted/CellVolume;
	    AccretedMomentum[0] += paccrete[0];
	    AccretedMomentum[1] += paccrete[1];
	    AccretedMomentum[2] += paccrete[2];

#ifdef DEBUG_AP
	    if (index == cgindex)
	      printf("Sink Density: %"GOUTSYM", Cell Density: %"GOUTSYM", New Density: %"GOUTSYM"\n",
		     maccreted/CellVolume,
		     BaryonField[DensNum][index],
		     BaryonField[DensNum][index]-maccreted/CellVolume);
#endif 

	    // Update the state variables
	    BaryonField[DensNum][index] -= maccreted/CellVolume;
	    
	    if (HydroMethod == 0) { // PPM
	      BaryonField[TENum][index] = etotnew/mnew;
	      if (DualEnergyFormalism)
		BaryonField[GENum][index] = eintnew/mnew;
	    } else // Zeus
	      BaryonField[TENum][index] = eintnew/mnew;

	    BaryonField[Vel1Num][index] = pnew[0]/mnew;
	    BaryonField[Vel2Num][index] = pnew[1]/mnew;
	    BaryonField[Vel3Num][index] = pnew[2]/mnew;

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

  float OldMass = (*ThisParticle)->Mass*CellVolume;
  float *OldVel = (*ThisParticle)->vel;
  
  float NewVelocity[3] = {
    (OldMass*OldVel[0]+AccretedMomentum[0])/(OldMass+AccretedMass*CellVolume),
    (OldMass*OldVel[1]+AccretedMomentum[1])/(OldMass+AccretedMass*CellVolume),
    (OldMass*OldVel[2]+AccretedMomentum[2])/(OldMass+AccretedMass*CellVolume)};

  (*ThisParticle)->AddMass(AccretedMass);
  (*ThisParticle)->SetVelocity(NewVelocity);

  delete [] nexcluded;

  return SUCCESS;
}

#undef DEBUG_AP
