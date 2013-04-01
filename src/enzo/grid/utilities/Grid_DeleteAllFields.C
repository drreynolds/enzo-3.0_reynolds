/***********************************************************************
/
/  GRID CLASS (REMOVE ALL FIELDS)
/
/  written by: Greg Bryan
/  date:       April, 1996
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
 
void grid::DeleteAllFields()
{
 
  int i;
 
  this->DeleteParticles();
  this->DeleteActiveParticles();
 
  for (i = 0; i < MAX_DIMENSION; i++) {
    delete [] ParticleAcceleration[i];
    delete [] ActiveParticleAcceleration[i];
    delete [] AccelerationField[i];
 
    ParticleAcceleration[i]      = NULL;
    ActiveParticleAcceleration[i]      = NULL;
    AccelerationField[i]         = NULL;
  }
  ParticleAcceleration[MAX_DIMENSION] = NULL;
  delete [] ParticleAcceleration[MAX_DIMENSION];
  ActiveParticleAcceleration[MAX_DIMENSION] = NULL;
  delete [] ActiveParticleAcceleration[MAX_DIMENSION];

  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
    delete [] BaryonField[i];
    delete [] OldBaryonField[i];
    BaryonField[i]    = NULL;
    OldBaryonField[i] = NULL;
  }

#ifdef SAB
  for (i = 0; i < MAX_DIMENSION; i++)
    if (OldAccelerationField[i] != NULL) {
      delete [] OldAccelerationField[i];
      OldAccelerationField[i] = NULL;
    }
#endif
 
  delete [] PotentialField;
  delete [] GravitatingMassField;
  delete [] GravitatingMassFieldParticles;
 
  PotentialField                = NULL;
  GravitatingMassField          = NULL;
  GravitatingMassFieldParticles = NULL;
 
}
