/***********************************************************************
/
/  GRID CLASS (BEFORE REBUILDING, REMOVED UNNEEDED ARRAYS)
/
/  written by: Greg Bryan
/  date:       June, 1995
/  modified1:
/
/  PURPOSE:
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
 
 
void grid::CleanUp()
{
 
  int i;
 
  for (i = 0; i < MAX_DIMENSION; i++) {
    delete [] ParticleAcceleration[i];
    delete [] ActiveParticleAcceleration[i];
//    delete [] AccelerationField[i];
 
    ParticleAcceleration[i]       = NULL;
    ActiveParticleAcceleration[i] = NULL;
//    AccelerationField[i]         = NULL;
  }
  delete [] ParticleAcceleration[MAX_DIMENSION];
  ParticleAcceleration[MAX_DIMENSION] = NULL;

  delete [] ActiveParticleAcceleration[MAX_DIMENSION];
  ActiveParticleAcceleration[MAX_DIMENSION] = NULL;
 
  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
    delete [] OldBaryonField[i];
    OldBaryonField[i] = NULL;
  }
 
  delete [] GravitatingMassField;
  delete [] GravitatingMassFieldParticles;
 
  GravitatingMassField          = NULL;
  GravitatingMassFieldParticles = NULL;

#ifdef SAB
  for (i = 0; i < MAX_DIMENSION; i++)
    if (OldAccelerationField[i] != NULL) {
      delete [] OldAccelerationField[i];
      OldAccelerationField[i] = NULL;
    }
#endif

}
