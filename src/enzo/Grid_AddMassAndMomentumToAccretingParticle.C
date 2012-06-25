/***********************************************************************
/
/  GRID CLASS (Find the accreting particle with id equal to ID on this 
/              grid's active particle list and increase its mass and 
/              momentum by AccretedMass and AccretedMomentum 
/              respectively)
/
/  written by: Nathan Goldbaum
/  date:       April 2012
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
#include "Hierarchy.h"
#include "ActiveParticle.h"

#define NO_DEBUG

int grid::AddMassAndMomentumToAccretingParticle(float AccretedMass, float AccretedMomentum[], 
						ActiveParticleType* ThisParticle, LevelHierarchyEntry *LevelArray[]) {

  // Return if this doesn't concern us
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  int i;
  int iFound = -1, pFound = -1;
  float CellVolume = 1.0;

  for (i = 0; i<NumberOfActiveParticles; i++) {
    if (this->ActiveParticles[i]->ReturnID() == ThisParticle->ReturnID()) {
      iFound = i;
      break;
    }
  }

  if (iFound == -1)
    return SUCCESS;

  for (i = 0; i < GridRank; i++)
    CellVolume*=CellWidth[i][0];

  float OldMass = ThisParticle->Mass*CellVolume;
  float *OldVel = ThisParticle->vel;

  float NewVelocity[3] = {
    (OldMass*OldVel[0]+AccretedMomentum[0])/(OldMass+AccretedMass*CellVolume),
    (OldMass*OldVel[1]+AccretedMomentum[1])/(OldMass+AccretedMass*CellVolume),
    (OldMass*OldVel[2]+AccretedMomentum[2])/(OldMass+AccretedMass*CellVolume)};

#ifdef DEBUG
  fprintf(stderr,"AccretedMass = %"GSYM"\n",AccretedMass*CellVolume);
  fprintf(stderr,"AccretedMomentum[0] = %"GSYM"\n",AccretedMomentum[0]);
  fprintf(stderr,"AccretedMomentum[1] = %"GSYM"\n",AccretedMomentum[1]);
  fprintf(stderr,"AccretedMomentum[2] = %"GSYM"\n",AccretedMomentum[2]);
#endif

  // Masses are actually densities
  this->ActiveParticles[iFound]->AddMass(AccretedMass);
  this->ActiveParticles[iFound]->SetVelocity(NewVelocity);

#ifdef DEBUG
  float mnew = this->ActiveParticles[iFound]->ReturnMass();
  float* vnew = this->ActiveParticles[iFound]->ReturnVelocity();
  float pnew[3] = {mnew*CellVolume*vnew[0],
		   mnew*CellVolume*vnew[1],
		   mnew*CellVolume*vnew[2]};

  fprintf(stderr,"dm = %"GSYM"\n",mnew*CellVolume - OldMass);
  fprintf(stderr,"dpx = %"GSYM"\n",pnew[0] - OldMass*OldVel[0]);
  fprintf(stderr,"dpy = %"GSYM"\n",pnew[1] - OldMass*OldVel[1]);
  fprintf(stderr,"dpz = %"GSYM"\n",pnew[2] - OldMass*OldVel[2]);
#endif

  // also need to update the mirrored normal particle

  for (i = 0; i<NumberOfParticles; i++) {
    if (this->ParticleNumber[i] == ThisParticle->ReturnID()) {
      pFound = i;
      break;
    }
  }

  // This only happen if the particles aren't mirrored properly
  if (pFound == -1)
    return FAIL;

  this->ParticleMass[pFound] = this->ActiveParticles[iFound]->Mass;
  for (i = 0; i < MAX_DIMENSION; i++)
    this->ParticleVelocity[i][pFound] = this->ActiveParticles[iFound]->vel[i];

  return SUCCESS;
}
