/***********************************************************************
/
/  Cen & Ostriker star formation
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
#include "TopGridData.h"
#include "EventHooks.h"
#include "ActiveParticle.h"

/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

class ActiveParticleType_CenOstriker;
class CenOstrikerBufferHandler;

class CenOstrikerGrid : private grid {
  friend class ActiveParticleType_CenOstrikerGrid;
};

class ActiveParticleType_CenOstriker : public ActiveParticleType
{
public:
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static int EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormtiondata &data);
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);
  static ParticleBufferHandler *AllocateBuffers(int NumberOfParticles);
};

int ActiveParticleType_CenOstriker::EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  CenOstrikerGrid *thisgrid =
    static_cast<SampleParticleGrid *>(thisgrid_orig);
  
  float BaryonMass,VelocityDivergence,TotalDensity,DynamicalTime,
    IsothermalSoundSpeedSquared,JeansMass,StarFraction, RandomNumber;
  float SoundSpeedConstant = 1.3095d8;
  int i, j, k, dim, index, offset_y, offset_z;
  int NumberOfNewParticles = 0;

  /* Define some pointers for readability */

  float *density = thisGrid->BaryonField[data.DensNum];
  float *velx = thisGrid->BaryonField[data.Vel1Num];
  float *vely = thisGrid->BaryonField[data.Vel2Num];
  float *velz = thisGrid->BaryonField[data.Vel3Num];

  FLOAT dx = data.LengthUnits * thisGrid->CellWidth[0][0]

  // Pre-calculate serialized offsets for the 3D data field.  Used for
  // the divergence.
  offset_y = thisGrid->GridDimension[0];
  offset_z = thisGrid->GridDimension[0] * thisGrid->GridDimension[1];

  for (k = thisGrid->GridStartIndex[2]; k <= thisGrid->GridEndIndex[2]; k++) {
    for (j = thisGrid->GridStartIndex[1]; j <= thisGrid->GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(thisGrid->GridStartIndex[0], j, k);
      for (i = thisGrid->GridStartIndex[0]; i <= thisGrid->GridEndIndex[0]; i++, index++) {
	
	// 0. If no more room for particles, quit.
	if (data.NumberOfNewParticles >=
	    data.MaxNumberOfNewParticles)
          continue;
	
	// 1. Finest level of refinement
	if (thisGrid->BaryonField[thisGrid->NumberOfBaryonFields][index] != 0.0) 
	  continue;
	
	// 2. Density greater than threshold
	if (density[index] < StarMakerOverDensityThreshold)
	  continue;
	
	/* 3. Negative divergence: For ZEUS, the velocities are
	   face-centered, and all of the other routines have
	   cell-centered velocities. */
	
	if (HydroMethod == Zeus_Hydro) {
	  VelocityDivergence = velx[index+1] - velx[index] +
	    vely[index+offset_y] - vely[index] + 
	    velz[index+offset_z] - velz[index];
	} else {
	  VelocityDivergence = velx[index+1] - velx[index-1] + 
	    vely[index+offset_y] - vely[index-offset_y] + 
	    velz[index+offset_z] - velz[index-offset_z];
	}

	if (VelocityDivergence > 0.0) continue;

	// 4. t_cool < t_freefall (skip if T > 11000 K)
	TotalDensity = ( density[index] + data.DarkMatterDensity[index] ) * 
	  data.DensityUnits;
	DynamicalTime = sqrt(3.0 * M_PI / 32.0 / GravConst / TotalDensity) / data.TimeUnits;
	if (DynamicalTime < data.CoolingTime[index] && 
	    data.Temperature[index] > 1.1e4)
	  continue;

	// 5. Cell mass is greater than the Jeans Mass

	if (CenOstrikerJeansMassCriterion) {
	  BaryonMass = density[index] * data.DensityUnits * dx / SolarMass;
	  IsothermalSoundSpeedSquared = SoundSpeedConstant * data.Temperature[index];
	  JeansMass = M_PI / (6.0 * sqrt(density[index] * data.DensityUnits)) *
	    POW(M_PI * IsothermalSoundSpeedSquared / GravConst,1.5) / SolarMass;
	  
	  if (BaryonMass < JeansMass)
	    continue;
	}

	// 6) Check to see if star is above threshold (given in units of M_solar)

	StarFraction = min(StarMakerMassEfficiency*thisGrid->ReturnTimeStep()/DynamicalTime, 0.9);
	DynamicalTime = max(DynamicalTime, StarMakerMinimumDynamicalTime*3.156e7/data.TimeUnits);
	
	// 7) If we allow stochastic star formation, make new particles every time the unfulfilled star formation buffer
	//    exceeds the mininimum particle mass

	if (CenOstrikerStochasticStarformation) {
	  if (StarFraction*BaryonMass < StarMakerMinimumMass) {
	    UnfulfilledStarFormationMass += StarFraction*BaryonMass;
	    if (UnfulfilledStarFormationMass < StarMakerMinimumMass) 
	      continue;
	    StarFraction = min(StarMakerMinimumMass/BaryonMass, 0.5);
	    UnfulfilledStarFormationMass -= StarFraction*BaryonMass;
	  } 
	}

	/*
	 * ====================================================================
	 * PARTICLE CREATION
	 * ====================================================================
	 */

	ActiveParticleType_CenOstriker *np = new ActiveParticleType_CenOstriker();
	data.NewParticles[data.NumberOfNewParticles++] = np;

	np->Mass = StarFraction*density[index];
	np->type = CenOstriker;
	np->BirthTime = thisGrid->ReturnTime();
	np->DynamicalTime = DynamicalTime;


	np->pos[0] = thisGrid->CellLeftEdge[0][i] + 0.5*thisGrid->CellWidth[0][i];
	np->pos[0] = thisGrid->CellLeftEdge[1][j] + 0.5*thisGrid->CellWidth[1][j];
	np->pos[0] = thisGrid->CellLeftEdge[2][k] + 0.5*thisGrid->CellWidth[2][k];
	

	double *tvel = thisGrid->AveragedVelocityAtCell(index,data.DensNum,data.Vel1Num);

	np->vel[0] = tvel[0];
	np->vel[1] = tvel[1];
	np->vel[2] = tvel[2];

	if (HasMetalField)
	  np->Metallicity = data.TotalMetals[index];
	else
	  np->Metallicity = 0.0;

	// Remove mass from grid

	density[index] = (1.0 - starfraction)*density[index]

      }
    }
  }
  return 0.;
}

int ActiveParticleType_CenOstriker::EvaluateFeedback(grid *thisGrid_orig, ActiveParticleFormationData &data)
{
  CenOstrikerGrid *thisgrid =
    static_cast<SampleParticleGrid *>(thisgrid_orig);
  
  /* Define some pointers for readability  */
  float *density = thisGrid->BaryonField[data.DensNum];
  float *velx = thisGrid->BaryonField[data.Vel1Num];
  float *vely = thisGrid->BaryonField[data.Vel2Num];
  float *velz = thisGrid->BaryonField[data.Vel3Num];
  float *totalenergy = thisGrid->BaryonField[data.TENum];
  float *gasenergy = thisGrid->BaryonField[data.GENum];
  float dt = thisGrid->dtFixed;
  float dx = float(thisGrid->CellWidth[0][0]);

  /* Find metallicity fields */

  int SNColourNum, MetalNum;

  if (this->IdentifyColourFields(SNColourNum, MetalNum, MBHColourNum,
				 Galaxy1ColourNum, Galaxy2ColourNum) == FAIL) {
    ENZO_FAIL("Error in grid->IdentifyColourFields.\n");
  }

  float *MetalField;
  float *TotalMetals = NULL;
  int MetallicityFlag;

  MetallicityFlag = (MetalNum != -1 || SNColourNum != -1);

  if (MetalNum != -1 && SNColourNum != -1) {
    TotalMetals = new float[size];
    for (i = 0; i < size; i++)
      TotalMetals[i] = BaryonField[MetalNum][i] + BaryonField[SNColourNum][i];
    MetalField = TotalMetals;
  } // ENDIF both metal types                                                                                                                                                                                                                                                 
  else {
    if (MetalNum != -1)
      MetalField = BaryonField[MetalNum];
    else if (SNColourNum != -1)
      MetalField = BaryonField[SNColourNum];
  } // ENDELSE both metal types

  float xv1,xv2,ParticleBirthTime,ParticleDynamicalTimeAtBirth,
    ParticleMass,ParticleInitialMass,ParticleMetalFraction,StarFormationDensityThisTimestep,
    SupernovaEnergyThisTimestep,DensityToAddToEachCell;

  float StellarMassFormedThistimestepOnThisGrid = 0;

  FLOAT xpos,ypos,zpos;

  FLOAT CurrentTime = thisGrid->Time;
  FLOAT xstart = thisGrid->CellLeftEdge[0];
  FLOAT ystart = thisGrid->CellLeftEdge[1];
  FLOAT zstart = thisGrid->CellLeftEdge[2];

  int npart = thisGrid->NumberOfParticles;
  int GridXSize = thisGrid->GridDimension[0];
  int GridYSize = thisGrid->GridDimension[1];
  int GridZSize = thisGrid->GridDimension[2];
  int NumberOfGhostZones = thisGrid->GridStartIndex[0];
  
  int n,i,j,k,index;

  for (n=0;npart-1;n++) {
    if (thisGrid->ActiveParticles[n]->ReturnType() == CenOstriker)
      continue;
  
    xpos = thisGrid->ActiveParticles[n]->pos[0];
    ypos = thisGrid->ActiveParticles[n]->pos[1];
    zpos = thisGrid->ActiveParticles[n]->pos[2];
  
    xvel = thisGrid->ActiveParticles[n]->vel[0];
    yvel = thisGrid->ActiveParticles[n]->vel[1];
    zvel = thisGrid->ActiveParticles[n]->vel[2];

    ParticleBirthTime = thisGrid->ActiveParticles[n]->BirthTime;
    ParticleDynamicalTimeAtBirth = thisGrid->ActiveParticles[n]->DynamicalTime;
    ParticleMass = thisGrid->ActiveParticles[n]->Mass;
    ParticleMetalFraction = thisGrid->ActiveParticles[n]->Metallicity
    
    // Determine how much of a given star particle would have been turned into stars during this timestep.
    // Then, calculate the mass which should have formed during this timestep dt using the integral form
    // of the Cen & Ostriker formula.
    
    xv1 = (CurrentTime - ParticleBirthTime) / ParticleDynamicalTimeAtBirth;
    if (xv1 > 12.0) continue; // current time - creation time >> dynamical time at formation, so ignore
    
    xv2 = (CurrentTime + dt - ParticleBirthTime) / ParticlDynamicalTimeAtBirth;

    // First, calculate the initial mass of the star particle in question
    
    ParticleInitialMass = ParticleMass / (1.0 - StarMassEjectionFraction(1.0 - (1.0 + xv1)*exp(-xv1)));
    
    // Then, calculate the amount of mass that would have formed in this timestep.

    StarFormationDensityThisTimestep = InitialParticleMass * ((1.0 + xv1)*exp(-xv1) - 
							   (1.0 + xv2)*exp(-xv2));
    
    StarFormationDensityThisTimestep = max(min(StarFormationDensityThisTimestep,ParticleMass),0.0);
      
    // Calculate 3D grid indices

    i = int((xpos - xstart)/dx);
    j = int((ypos - ystart)/dy);
    k = int((zpos - zstart)/dz);

    // Check bounds - if star particle is outside of this grid then give a warning and continue
    
    if (i < 0 || i > GridXSize-1 || j < 0 || j > GridYSize-1 || k < 0 || k > GridZSize-1){
      fprintf(stdout, "Particle out of grid; xind, yind, zind, level = %d, $d, $d, $d\n",i,j,k,);
      continue;
    }
      
    // Calculate serial index

    index = GRIDINDEX_NOGHOST(thisGrid->GridStartIndex[0],j,k);

    // skip if very little mass if formed

    if (StarFormationDensityThisTimestep/density[index] < 1e-10 )
      continue;

    // calculate mass added to each cell

    DensityToAddToEachCell = (StarFormationDensityThisTimestep * StarMassEjectionFraction) / StarFeedbackDistTotalCells;

    // If using distributed feedback, check if particle is too close to the boundary and adjust indices accordingly

    if (StarFeedbackDistRadius > 0)
      {
	i = max(1+NumberOfGhostZones+StarFeedbackDistRadius,
		min(GridXSize - NumberOfGhostZones - StarFeedbackDistRadius,i));
	j = max(1+NumberOfGhostZones+StarFeedbackDistRadius,
		min(GridYSize - NumberOfGhostZones - StarFeedbackDistRadius,j));
	k = max(1+NumberOfGhostZones+StarFeedbackDistRadius,
		min(GridZSize - NumberOfGhostZones - StarFeedbackDistRadius,k));	
      }

    // Subtract ejected mass from particle
    
    ParticleMass -= StarFormationDensityThisTimestep*StarMassEjectionFraction;

    // Save particle mass

    thisGrid->ActiveParticles[n]->Mass = ParticleMass

    // Record amount of star formation in this grid

    StellarMassFormedThisTimestepOnThisGrid += StarFormationDensityThisTimestep*dt*POW(dx,3);

    // Calculate supernova energy for this timestep

    SupernovaEnergyThisTimestep = StarEnergyToThermalFeedback * StarFormationDensityThisTimestep * 
      POW(clight/data.VelocityUnits,2) / StarFeedbackDistTotalCells;

#define NO_SHARE_ENERGY
#ifdef SHARE_ENERGY
    SupernovaEnergyThisTimestep *= StarFormationDensityThisTimestep*StarMassEjectionFraction / 
      (StarFormationDensityThisTimestep*StarMassEjectionFraction + ParticleInitialMAss*exp(-xv2)*(1.0+xv2));
#endif /* SHARE_ENERGY */

    // Add energy to the energy field
    for (kc = k - StarFeedbackDistRadius; kc > k + StarFeedbackDistRadius, kc++){
      stepk = abs(kc - k);
      for (jc = j - StarFeedbackDistRadius; jc > j + StarFeedbackDistRadius, jc++){
	stepj = stepk + abs(jc - j);
	for (ic = i - StarFeedbackDistRadiusl ic > i + StarFeedbackDistRadius, ic++){
	  cellstep = stepj + abs(ic - i);
	  DistIndex = GRIDINDEX_NOGHOST(thisGrid->GridStartIndex[0],jc,kc);
	  if (cellstep < StarFeedbackDistCellStep) {
	    DensityRatio = 1.0/(density[DistIndex] + DensityToAddToEachCell);
	    TotalEnergy[DistIndex] = ((TotalEnergy[DistIndex]*density[DistIndex]) + SupernovaEnergyThisTimestep)*DensityRatio;
	    if (DualEnergyFormalism == 1)
	      GasEnergy[DistIndex] = ((GasEnergy[DistIndex]*density[DistIndex]) + SupernovaEnergyThisTimestep)*DensityRatio;

	    // Metal feedback (note that in this function gas metal is a fraction
	    // (rho_metal/rho_gas) rather than a density.  The conversion has 
	    // been done in the handling routine)

	    // The "Cen Method".  This takes into account gas recycling:

	    if (MetallicityFlag == 1)
	      MetalField[DistIndex] = (MetalField[DistIndex]*density[DistIndex]) +
		(StarFormationDensityThisTimestep / StarFeedbackDistTotalCells) *
		(StarMetalYield * (1.0 - ParticleMetalFraction) +
		 StarMassEjectionFraction * ParticleMetalFraction) * DensityRatio;

	    // Mass and momentum feedback

	    velx[DistIndex] = velx[DistIndex]*density[DistIndex] + DensityToAddToEachCell * xvel;
	    vely[DistIndex] = vely[DistIndex]*density[DistIndex] + DensityToAddToEachCell * yvel;
	    velz[DistIndex] = velz[DistIndex]*density[DistIndex] + DensityToAddToEachCell * zvel;
	    density[DistIndex] += DensityToAddToEachCell;
	    velx[DistIndex] /= density[DistIndex];
	    vely[DistIndex] /= density[DistIndex];
	    velz[DistIndex] /= density[DistIndex];
	      

	  }
	}
      }
    }

    
  } // end loop over particles
  
  return SUCCESS;
}

void ActiveParticleType_CenOstriker::DescribeSupplementalData(ActiveParticleFormationDataFlags &flags)
{
  flags.DarkMatterDensity = true;
  flags.CoolingTime = true;
  flags.Temperature = true;
  flags.UnitConversions = true;
  flags.DataFieldNumbers = true;
  flags.MetalField = true;
}

class CenOstrikerParticleBufferHandler : public ParticleBufferHandler
{
  public:
    CenOstrikerBufferHandler(int NumberOfParticles) { }
};

ParticleBufferHandler *ActiveParticleType_CenOstriker::AllocateBuffers(int NumberOfParticles)
{
    CenOstrikerBufferHandler *handler = new CenOstrikerBufferHandler(NumberOfParticles);
    return handler;
}


namespace {
  ActiveParticleType_info *CenOstrikerInfo = new ActiveParticleType_info
    ("CenOstrikerParticle", (&ActiveParticleType_CenOstriker::EvaluateFormation),
     (&ActiveParticleType_CenOstriker::DescribeSupplementalData),
     (&ActiveParticleType_CenOstriker::AllocateBuffers),
     (&ActiveParticleType_CenOstriker::EvaluateFeedback));
}
