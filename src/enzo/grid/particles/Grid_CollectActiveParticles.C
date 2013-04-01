/***********************************************************************
/
/  GRID CLASS (COLLECT ACTIVE PARTICLES INTO ONE PROCESSOR)
/
/  written by: John Wise
/  date:       May, 2009
/  modified1:  John Wise -- re-purposing for active particles
/  date:       December, 2011
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
#include "ActiveParticle.h"

int grid::CollectActiveParticles(int GridNum, int* &NumberToMove, 
				 int &StartIndex, int &EndIndex, 
				 ActiveParticleType** &List, int CopyDirection)
{
 
  /* Declarations. */

  int i, j, dim, n1, grid, proc;

  /* ----------------------------------------------------------------- */
  /* Copy star out of grid. */

  if (CopyDirection == COPY_OUT) {

    /* If there are no stars to move, we're done. */

    if (NumberOfActiveParticles == 0)
      return SUCCESS;

    /* If this is the correct processor, no copy-outs required. */

    if (MyProcessorNumber == ProcessorNumber)
      return SUCCESS;

    /* Add to the star count to move */

    // NumberOfActiveParticles is still the number of local stars, not the
    // actual total!
    NumberToMove[ProcessorNumber] += NumberOfActiveParticles;
 
    /* Move and delete stars */

    if (ActiveParticles == NULL)
      ENZO_FAIL("ActiveParticles pointer cannot be NULL here.  "
		"NumberOfActiveParticles and pointer are mismatched.");

    for (i = 0, n1 = StartIndex; i < NumberOfActiveParticles; i++, n1++) {
      List[n1] = ActiveParticles[i];
      List[n1]->SetGridID(GridNum);
      List[n1]->SetDestProcessor(ProcessorNumber);
    } // ENDFOR stars

    StartIndex = n1;
    NumberOfActiveParticles = 0;
    delete [] ActiveParticles;
    ActiveParticles = NULL;

  } // end: if (COPY_OUT)
 
  /* ----------------------------------------------------------------- */
  /* Copy stars back into grid. */
 
  else {

    if (MyProcessorNumber != ProcessorNumber)
      return SUCCESS;

    int NumberOfNewActiveParticles = EndIndex - StartIndex;
    this->AddActiveParticles(List, NumberOfNewActiveParticles, StartIndex);
 
  } // end: if (COPY_IN)
 
  return SUCCESS;
}
