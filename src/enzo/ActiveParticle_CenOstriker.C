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

int ActiveParticleType_SampleParticle::EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data)
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

	BaryonMass = density[index] * data.DensityUnits * dx / SolarMass;
	IsothermalSoundSpeedSquared = SoundSpeedConstant * data.Temperature[index];
	JeansMass = M_PI / (6.0 * sqrt(density[index] * data.DensityUnits)) *
	  POW(M_PI * IsothermalSoundSpeedSquared / GravConst,1.5) / SolarMass;
	
	if (BaryonMass < JeansMass)
	  continue;

	// 6) Check to see if star is above threshold (given in units of M_solar)

	StarFraction = min(StarMakerMassEfficiency*thisGrid->dtFixed/tdyn, 0.9);
	
	//RandomNumber = mt_random();
	//if (RandomNumber >= exp(-1.0 * StarMakerMinimumDynamicalTime / tdyn))
	//  continue;

	/*
	 * ====================================================================
	 * PARTICLE CREATION
	 * ====================================================================
	 */

	ActiveParticleType_CenOstriker *np = new ActiveParticleType_CenOstriker();
	data.NewParticles[data.NumberOfNewParticles++] = np;

	np->Mass = StarFraction*density[index];
	np->type = CenOstriker;
	np->BirthTime = thisGrid->Time;

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
  
  FLOAT xstart = thisGrid->CellLeftEdge[0];
  FLOAT ystart = thisGrid->CellLeftEdge[1];
  FLOAT zstart = thisGrid->

  int npart = thisGrid->NumberOfParticles;
  int n;
  

  for (n=0;npart-1;n++) {

    

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

class SampleParticleBufferHandler : public ParticleBufferHandler
{
  public:
    SampleParticleBufferHandler(int NumberOfParticles) { }
};

ParticleBufferHandler *ActiveParticleType_SampleParticle::AllocateBuffers(int NumberOfParticles)
{
    SampleParticleBufferHandler *handler = new SampleParticleBufferHandler(NumberOfParticles);
    return handler;
}


namespace {
  ActiveParticleType_info *CenOstrikerInfo = new ActiveParticleType_info
    ("CenOstrikerParticle", (&ActiveParticleType_CenOstriker::EvaluateFormation),
     (&ActiveParticleType_CenOstriker::DescribeSupplementalData),
     (&ActiveParticleType_CenOstriker::AllocateBuffers),
     (&ActiveParticleType_CenOstriker::EvaluateFeedback));
}
