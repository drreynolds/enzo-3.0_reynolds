/***********************************************************************
/
/  GRID CLASS (ADD ACTIVE PARTICLE DATA TO GRID)
/
/  written by: John Wise
/  date:       December, 2011
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "fortran.def"
#include "CosmologyParameters.h"

#include "ActiveParticle.h"

int grid::AddActiveParticles(ActiveParticleType **NewParticles,
			     int NumberOfNewParticles, int start)
{

  if (NumberOfNewParticles == 0)
    return SUCCESS;

  int i, index;
  int OldNumberOfActiveParticles = this->NumberOfActiveParticles;
  ActiveParticleType **OldActiveParticles = this->ActiveParticles;

  this->NumberOfActiveParticles += NumberOfNewParticles;
  this->ActiveParticles = new ActiveParticleType*[this->NumberOfActiveParticles]();
  for (i = 0; i < OldNumberOfActiveParticles; i++) {
    this->ActiveParticles[i] = OldActiveParticles[i];
  }
  for (i = start, index = OldNumberOfActiveParticles; 
       index < this->NumberOfActiveParticles; i++, index++) {
    this->ActiveParticles[index] = NewParticles[i];
    this->ActiveParticles[index]->SetGridID(this->ID);
    this->ActiveParticles[index]->AssignCurrentGrid(this);
  }

#define NO_DEBUG
#ifdef DEBUG
  int dim, inside;
  FLOAT *pos;
  float TotalMass = 0;
  for (i = 0; i < this->NumberOfActiveParticles; i++) {
    pos = this->ActiveParticles[i]->ReturnPosition();
    TotalMass += this->ActiveParticles[i]->ReturnMass();
    inside = this->PointInGrid(pos);
    if (inside == FALSE) {
      fprintf(stdout,"pos[0]: %"PSYM", pos[1]: %"PSYM", pos[2]: %"PSYM"\n",pos[0],pos[1],pos[2]);
      fprintf(stdout,"mass: %"FSYM"\n");
      fprintf(stdout,"GridLeftEdge[0]: %"PSYM", GridLeftEdge[1]: %"PSYM", GridLeftEdge[2]: %"PSYM"\n",
	      GridLeftEdge[0], GridRightEdge[1], GridRightEdge[2]);
      ENZO_FAIL("ActiveParticle outside!\n");
    }
  }
  fprintf(stdout,"AddActiveParticles: Total Mass added to grid = %"FSYM"\n",TotalMass);
#endif /* DEBUG */  

  delete [] OldActiveParticles;

  return SUCCESS;

}
