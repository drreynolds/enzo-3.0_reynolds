/***********************************************************************
/
/  SUBTRACT ACCRETED MASS FROM THE GRID (AND CHANGES TO MOMENTUM)
/
/  written by: Ji-hoon Kim
/  date:       January, 2010
/  modified1: 
/
/ PURPOSE: To subtract mass from the grid, we must consider multiple grids
/          since sometimes the accretion radius often exceeds the grid
/          boundaries.  
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
#include "LevelHierarchy.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int StarParticleSubtractAccretedMass(TopGridData *MetaData, 
				     LevelHierarchyEntry *LevelArray[], int level, 
				     Star *&AllStars)
{

  const double pc = 3.086e18, Msun = 1.989e33, pMass = 1.673e-24, 
    gravConst = 6.673e-8, yr = 3.1557e7, Myr = 3.1557e13;

  Star *cstar;
  int i, l, dim, temp_int, SkipMassRemoval, SphereContained,
      SphereContainedNextLevel, dummy;
  float influenceRadius, RootCellWidth, SNe_dt, mdot;
  float dtForThisStar;
  double EjectaThermalEnergy, EjectaDensity, EjectaMetalDensity, EjectaVolume;
  FLOAT Time;
  LevelHierarchyEntry *Temp;

  if (AllStars == NULL)
    return SUCCESS;

  /* Get time and SNe timestep */

  Temp = LevelArray[level];
  Time = Temp->GridData->ReturnTime();
  if (LastSupernovaTime < 0)
    SNe_dt = 0.0;
  else
    SNe_dt = Time - LastSupernovaTime;
  LastSupernovaTime = Time;
  RootCellWidth = 1.0 / MetaData->TopGridDims[0];

  /* Set the units. */

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);

  for (cstar = AllStars; cstar; cstar = cstar->NextStar) {

//    if (debug)
//      cstar->PrintInfo();

    if (cstar->ReturnType() != BlackHole && 
	ABS(cstar->ReturnType()) != MBH)
      continue;

    if ((ABS(cstar->ReturnType()) == MBH && MBHAccretion <= 0) ||
	cstar->ReturnLastAccretionRate() < tiny_number)
      continue;

    /* Now let us do the job! */

    switch (cstar->ReturnType()) {
    case BlackHole:

      /* If BlackHole, subtract mass from a single cell (old way) so as 
	 not to bother John or others at the moment! */ 

      if (cstar->SubtractAccretedMassFromCell() == FAIL) {
	fprintf(stderr, "Error in star::SubtractAccretedMass.\n");
	ENZO_FAIL("");
      }

      break;

    case MBH:

      /* If MBH, subtract mass from multiple cells 
	 defined by MBHAccretionRadius */

      dtForThisStar = LevelArray[level]->GridData->ReturnTimeStep();

      /* Compute some parameters, similar to Star_CalculateFeedbackParameters */

      cstar->CalculateSubtractionParameters(LevelArray, influenceRadius, RootCellWidth, 
           EjectaDensity, DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
	   VelocityUnits, dtForThisStar);

//      fprintf(stdout, "SPSAM: EjectaDensity=%g, influenceRadius=%g\n", EjectaDensity, influenceRadius); 

      /* Determine if a sphere with influenceRadius is enclosed within grids on this level 
	 we use FindFeedbackSphere function, but we are not doing any "feedback" here; 
	 we simply subtract the mass */

      if (cstar->FindFeedbackSphere(LevelArray, level, influenceRadius, 
				    EjectaDensity, EjectaThermalEnergy, 
				    SphereContained, dummy, DensityUnits, 
				    LengthUnits, TemperatureUnits, TimeUnits, 
				    VelocityUnits, Time) == FAIL) {
	ENZO_FAIL("Error in star::FindFeedbackSphere");
      }

      /* If something weird happens, don't bother. */

      if ( influenceRadius <= tiny_number || 
	   influenceRadius >= RootCellWidth/2 )
	break;

      /* Determine if a sphere is enclosed within the grids on next level
	 If that is the case, we perform SubtractAccretedMass not here, 
	 but in the EvolveLevel of the next level. */
      
      SphereContainedNextLevel = FALSE;

      if (LevelArray[level+1] != NULL) {
	if (cstar->FindFeedbackSphere(LevelArray, level+1, influenceRadius, 
				      EjectaDensity, EjectaThermalEnergy, 
				      SphereContainedNextLevel, dummy, DensityUnits, 
				      LengthUnits, TemperatureUnits, TimeUnits, 
				      VelocityUnits, Time) == FAIL) {
	  fprintf(stderr, "Error in star::FindFeedbackSphere\n");
	  ENZO_FAIL("");
	}
      }

//    fprintf(stdout, "SPSAM: SkipMassRemoval=%d, SphereContained=%d, SphereContainedNextLevel=%d\n", 
//	    SkipMassRemoval, SphereContained, SphereContainedNextLevel);  

      /* Quit this routine when 
	 (1) sphere is not contained, or 
	 (2) sphere is contained, but the next level can contain the sphere, too. */ 
      if ((SphereContained == FALSE) ||
	  (SphereContained == TRUE && SphereContainedNextLevel == TRUE))
	break;
    
      /* Now set cells within the radius to their values after subtraction. */
      
      int CellsModified = 0;
      
      for (l = level; l < MAX_DEPTH_OF_HIERARCHY; l++)
	for (Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel) 
	  Temp->GridData->SubtractAccretedMassFromSphere
	    (cstar, l, influenceRadius, DensityUnits, LengthUnits, 
	     VelocityUnits, TemperatureUnits, TimeUnits, EjectaDensity, 
	     CellsModified);

//    fprintf(stdout, "StarParticleSubtractAccretedMass[%"ISYM"][%"ISYM"]: "
//	    "Radius = %e pc, changed %"ISYM" cells.\n", 
//	    cstar->ReturnID(), level, influenceRadius*LengthUnits/pc, CellsModified); 

#ifdef UNUSED
      temp_int = CellsModified;
      MPI_Reduce(&temp_int, &CellsModified, 1, MPI_INT, MPI_SUM, ROOT_PROCESSOR,
		 MPI_COMM_WORLD);

      if (debug) {
	if (cstar->ReturnFeedbackFlag() != FORMATION)
	  fprintf(stdout, "StarParticleSubtractAccretedMass[%"ISYM"][%"ISYM"]: "
		  "Radius = %"GSYM" pc\n",
		  cstar->ReturnID(), level, influenceRadius*LengthUnits/pc);
	fprintf(stdout, "StarParticleSubtractAccretedMass[%"ISYM"][%"ISYM"]: "
		"changed %"ISYM" cells.\n", 
		cstar->ReturnID(), level, CellsModified);
      }
#endif      

      break;

    } // ENDSWITCH ReturnType()

  } // ENDFOR stars

  return SUCCESS;

}