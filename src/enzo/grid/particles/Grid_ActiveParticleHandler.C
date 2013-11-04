/***********************************************************************
/
/  GRID CLASS (HANDLE THE CREATION AND FEEDBACK OF ACTIVE PARTICLES)
/
/  written by: Greg Bryan
/  date:       March, 1997
/  modified1:  April, 2009 by JHW to have multiple types of star 
/              particles
/  modified2:  May, 2011 by MJT to be in support of active particles
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "ActiveParticle.h"

/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);
 
int grid::ActiveParticleHandler(HierarchyEntry* SubgridPointer, int level,
                                float dtLevelAbove, int &NumberOfNewParticles)
{

  if (EnabledActiveParticlesCount == 0) return SUCCESS;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  /*fprintf(stderr, "G_APH: Currently have %"ISYM"\n",
          this->NumberOfActiveParticles);*/
 
  /* First, set under_subgrid field */
  HierarchyEntry *Subgrid;
  this->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
  for (Subgrid = SubgridPointer; Subgrid; Subgrid = Subgrid->NextGridThisLevel)
    this->ZeroSolutionUnderSubgrid(Subgrid->GridData, ZERO_UNDER_SUBGRID_FIELD);

  /* initialize */

  LCAPERF_START("grid_ActiveParticleHandler");

  /* First we identify the data dependencies */

  struct ActiveParticleFormationDataFlags flags = flags_default;

  int i;
  for (i = 0; i < EnabledActiveParticlesCount; i++)
  {
    ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];
    ActiveParticleTypeToEvaluate->DescribeSupplementalData(flags);
  }

  struct ActiveParticleFormationData supplemental_data = data_default;
  supplemental_data.level = level;
  supplemental_data.GridID = this->ID;

  ActiveParticleType::ConstructData(this, flags, supplemental_data);

  /******************** FORMATION ********************/

  NumberOfNewParticles = 0;
  /* Now we iterate */
  for (i = 0; i < EnabledActiveParticlesCount; i++)
  {
    ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];
    ActiveParticleTypeToEvaluate->EvaluateFormation(
                                this, supplemental_data);
    NumberOfNewParticles += supplemental_data.NumberOfNewParticles;
    
  }

  /* Now we copy the particles from NewParticles into a statically allocated
   * array */

  if (NumberOfNewParticles > 0) {
    this->AddActiveParticles(supplemental_data.NewParticles,
			     NumberOfNewParticles);
    if (debug2)
      printf("Creating %d new active particles\n", NumberOfNewParticles);
  }

  /******************** FEEDBACK ********************/

  /* Now we iterate */
  if (NumberOfActiveParticles > 0)
    for (i = 0; i < EnabledActiveParticlesCount; i++)
      {
	ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];
	ActiveParticleTypeToEvaluate->EvaluateFeedback(this, supplemental_data);
      }
  
  ActiveParticleType::DestroyData(this, supplemental_data);

  //if (debug) printf("StarParticle: end\n");

  LCAPERF_STOP("grid_ActiveParticleHandler");
  return SUCCESS;
}
