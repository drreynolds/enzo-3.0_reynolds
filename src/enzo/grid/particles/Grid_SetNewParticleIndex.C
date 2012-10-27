/***********************************************************************
/
/  GRID CLASS (SEARCH FOR ALL STAR PARTICLES AND RETURN HOW MANY)
/
/  written by: Greg Bryan
/  date:       September, 2000
/  modified1:  JHK & JHW (2009)
/  modified2:  John Wise (December, 2011) -- star->active particles
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
#include "ActiveParticle.h"

#define DEBUG

void grid::SetNewParticleIndex(PINT &next_id)
{
  int n, abstype;
  int ori_count = next_id;
  
  for (n = 0; n < NumberOfActiveParticles; n++)
    if (ActiveParticles[n]->Identifier == INT_UNDEFINED) {
      ActiveParticles[n]->Identifier = next_id++;
#ifdef DEBUG
      if (ActiveParticles[n]->Identifier == 276422)
	printf("Debug trap\n!");
      std::cout << "SNPI[" << MyProcessorNumber << "] " << "GridID: " 
		<< this->ID << " APID: " << ActiveParticles[n]->Identifier
		<< std::endl;
#endif
    }
  // Do the same for mirrored particles.  The normal and active new
  // particles are still in the same order as they were created.
  for (n = NumberOfParticles-NumberOfActiveParticles; 
       n < NumberOfParticles; n++)
    if (ParticleNumber[n] == INT_UNDEFINED)
      ParticleNumber[n] = ori_count++;
  return;
}
