/***********************************************************************
/
/  GRID CLASS (UPDATE PARTICLE VELOCITY FROM ACCELERATIONS)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
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
 
/* function prototypes */
 
#define VELOCITY_METHOD3
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
 
int grid::UpdateParticleVelocity(float TimeStep)
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  if ((NumberOfParticles == 0 && NumberOfActiveParticles == 0) || SelfGravity == FALSE) return SUCCESS;
 
  FLOAT a = 1.0, dadt;
#if defined(VELOCITY_METHOD1) || defined(VELOCITY_METHOD2)
  float VelocityMidStep;
#endif
  int i, dim, dim1;
 
  /* If using comoving coordinates, divide by a(t) first. */
 
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(Time + TimeStep, &a, &dadt)
	== FAIL) {
            ENZO_FAIL("Error in CsomologyComputeExpansionFactors.");
    }
 
  /* Loop over dimensions. */
 
  for (dim = 0; dim < GridRank; dim++) {
 
    /* Error check. */
   
    if ((NumberOfParticles > 0) && (ParticleAcceleration[dim] == NULL)) {
            ENZO_FAIL("No ParticleAccleration present.");
    }
 
    if (NumberOfParticles > 0 && ParticleAcceleration[dim] == NULL) {
            ENZO_FAIL("No ParticleAccleration present.");
    }
 

    /* Update velocities.  */
 
    if (ComovingCoordinates) {
      
      FLOAT coef = 0.5*dadt/a*TimeStep;
      FLOAT coef1 = 1.0 - coef;
      FLOAT coef2 = 1.0 / (1.0 + coef);
      
      /* If using comoving coordinates, subtract the (time-centered)
	 drag-like term and add the acceleration. The acceleration has
	 already been divided by a(t). */
      
      for (i = 0; i < NumberOfParticles; i++) {
	
#ifdef VELOCITY_METHOD1
	
        /* i) partially time-centered. */
	
	VelocityMidStep = ParticleVelocity[dim][i] +
	  ParticleAcceleration[dim][i]*0.5*TimeStep;
	
	ParticleVelocity[dim][i] +=
	  (-VelocityMidStep*dadt/a + ParticleAcceleration[dim][i]) * TimeStep;
	
#endif /* VELOCITY_METHOD1 */
	
#ifdef VELOCITY_METHOD2
	
        /* ii) partially backward. */
	
	VelocityMidStep = ParticleVelocity[dim][i] ;
	
	ParticleVelocity[dim][i] +=
	  (-VelocityMidStep*dadt/a + ParticleAcceleration[dim][i]) * TimeStep;
	
#endif /* VELOCITY_METHOD2 */
	
#ifdef VELOCITY_METHOD3
	
        /* iii) Semi-implicit way */
	
        ParticleVelocity[dim][i] = (coef1*ParticleVelocity[dim][i] +
                                    ParticleAcceleration[dim][i]*TimeStep)*coef2;
	
#endif /* VELOCITY_METHOD3 */
	
      }

      if (NumberOfActiveParticles > 0) {
	  
	float **ActiveParticleVelocity = new float*[GridRank];
	for (dim1 = 0; dim1 < GridRank; dim1++)
	  ActiveParticleVelocity[dim1] = new float[NumberOfActiveParticles];
	this->GetActiveParticleVelocity(ActiveParticleVelocity);
	
	for (i = 0; i < NumberOfActiveParticles; i++) {
	
#ifdef VELOCITY_METHOD1
	
          /* i) partially time-centered. */
	
	  VelocityMidStep = ActiveParticleVelocity[dim][i] +
	    ActiveParticleAcceleration[dim][i]*0.5*TimeStep;
	
	  ActiveParticleVelocity[dim][i] +=
	    (-VelocityMidStep*dadt/a + ActiveParticleAcceleration[dim][i]) * TimeStep;
	
#endif /* VELOCITY_METHOD1 */
	  
#ifdef VELOCITY_METHOD2
	  
	  /* ii) partially backward. */
	  
	  VelocityMidStep = ActiveParticleVelocity[dim][i] ;
	  
	  ActiveParticleVelocity[dim][i] +=
	    (-VelocityMidStep*dadt/a + ActiveParticleAcceleration[dim][i]) * TimeStep;
	  
#endif /* VELOCITY_METHOD2 */
	  
#ifdef VELOCITY_METHOD3
	  
	  /* iii) Semi-implicit way */
	  
	  ActiveParticleVelocity[dim][i] = (coef1*ActiveParticleVelocity[dim][i] +
	                         ActiveParticleAcceleration[dim][i]*TimeStep)*coef2;
	  
 
#endif /* VELOCITY_METHOD3 */
	  
        }
	for (dim1 = 0; dim1 < GridRank; dim1++)
	  delete [] ActiveParticleVelocity[dim1];
	delete [] ActiveParticleVelocity;
      }
    }
    else
 
      /* Otherwise, just add the acceleration. */
 
    for (i = 0; i < NumberOfParticles; i++)
      ParticleVelocity[dim][i] += ParticleAcceleration[dim][i] * TimeStep;
    for (i = 0; i < NumberOfActiveParticles; i++) {
      float **ActiveParticleVelocity = new float*[GridRank];
      for (dim1 = 0; dim1 < GridRank; dim1++)
	ActiveParticleVelocity[dim1] = new float[NumberOfActiveParticles];

      this->GetActiveParticleVelocity(ActiveParticleVelocity);

      ActiveParticleVelocity[dim][i] += ActiveParticleAcceleration[dim][i] * TimeStep;

      for (int dim = 0; dim1 < GridRank; dim1++)
	delete [] ActiveParticleVelocity[dim1];
      delete [] ActiveParticleVelocity;
    }
  }


  if (ProblemType == 29)
    for (i = 0; i < NumberOfParticles; i++)
      printf("id=%"PISYM"  %"PSYM" %"PSYM" %"PSYM"\n", ParticleNumber[i],
	     ParticlePosition[0][i], ParticlePosition[1][i], ParticlePosition[2][i]);

 
  return SUCCESS;
}
