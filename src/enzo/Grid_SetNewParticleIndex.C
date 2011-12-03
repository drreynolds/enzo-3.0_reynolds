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

#include <map>
#include <iostream>
#include <stdexcept>
#include <stdlib.h> 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "ActiveParticle.h"

void grid::SetNewParticleIndex(int &NumberCount1, PINT &NumberCount2)
{
  int n, abstype;
  for (n = 0; n < NumberOfActiveParticles; n++)
    if (ActiveParticles[n]->Identifier == INT_UNDEFINED) {
      ActiveParticles[n]->Identifier = NumberCount1 + NumberCount2++;
      printf("New star particle index = %d (%d %d)\n",
	     ActiveParticles[n]->Identifier, NumberCount1, NumberCount2);
    }
  return;
}
