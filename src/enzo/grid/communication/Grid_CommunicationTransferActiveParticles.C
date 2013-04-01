/***********************************************************************
/
/  GRID CLASS (COPY ACTIVE PARTICLES INTO OR OUT OF GRID)
/
/  written by: Greg Bryan
/  date:       January, 1999
/  modified1:  Robert Harkness
/  date:       April, 2006
/  modified2:  May, 2009 by John Wise: optimized version to transfer
/                particles in one sweep with collective calls.
/  modified3:  July, 2009 by John Wise: adapted for stars
/  modified4:  February, 2012 by John Wise: adapted for active particles
/
/  PURPOSE:
/
************************************************************************/
#ifdef USE_MPI
#include <communicators.h>
#endif
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

int search_lower_bound(int *arr, int value, int low, int high, 
		       int total);

int grid::CommunicationTransferActiveParticles(grid* Grids[], int NumberOfGrids, 
	       int ThisGridNum, int TopGridDims[], int *&NumberToMove, 
	       int StartIndex, int EndIndex, ActiveParticleType** &List, 
	       int *Layout, int *GStartIndex[], int *GridMap, int CopyDirection)
{
 
  /* Declarations. */
 
  int i, j, k, dim, grid, proc, grid_num, width, bin, CenterIndex;
  int GridPosition[MAX_DIMENSION];
  FLOAT r[MAX_DIMENSION];
  int *ToGrid, *pbin;

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    GridPosition[dim] = 0;
 
  /* ----------------------------------------------------------------- */
  /* Copy active particles out of grid. */
 
  if (CopyDirection == COPY_OUT) {

    /* If there are no active particles to move, we're done. */
 
    if (NumberOfActiveParticles == 0)
      return SUCCESS;

//    if (MyProcessorNumber != ProcessorNumber)
//      return SUCCESS;
 
    /* Count the number of active particles already moved */

    int PreviousTotalToMove = 0;
    for (i = 0; i < NumberOfProcessors; i++)
      PreviousTotalToMove += NumberToMove[i];

    /* Count active particles to move.  Apply perioidic wrap to them. */
 
    ToGrid = new int[NumberOfActiveParticles];

    float DomainWidth[MAX_DIMENSION], DomainWidthInv[MAX_DIMENSION];
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
      DomainWidthInv[dim] = 1.0/DomainWidth[dim];
    }

    // Periodic boundaries
    for (dim = 0; dim < GridRank; dim++) 
      for (i = 0; i < NumberOfActiveParticles; i++) {
	if (ActiveParticles[i]->pos[dim] > DomainRightEdge[dim])
	  ActiveParticles[i]->pos[dim] -= DomainWidth[dim];
	else if (ActiveParticles[i]->pos[dim] < DomainLeftEdge[dim])
	  ActiveParticles[i]->pos[dim] += DomainWidth[dim];
      }


    for (i = 0; i < NumberOfActiveParticles; i++) {

      for (dim = 0; dim < GridRank; dim++) {

	if (Layout[dim] == 1) {
	  GridPosition[dim] = 0;
	} else {

	  CenterIndex = 
	  (int) (TopGridDims[dim] * 
		 (ActiveParticles[i]->pos[dim] - DomainLeftEdge[dim]) *
		 DomainWidthInv[dim]);

	  GridPosition[dim] = 
	    search_lower_bound(GStartIndex[dim], CenterIndex, 0, Layout[dim],
			       Layout[dim]);
	  GridPosition[dim] = min(GridPosition[dim], Layout[dim]-1);

	} // ENDELSE Layout

      } // ENDFOR dim

      grid_num = GridPosition[0] + 
	Layout[0] * (GridPosition[1] + Layout[1]*GridPosition[2]);
      grid = GridMap[grid_num];
      if (grid != ThisGridNum) {
	proc = Grids[grid]->ReturnProcessorNumber();
	NumberToMove[proc]++;
      }

      ToGrid[i] = grid;

    } // ENDFOR active particles

    /* Allocate space. */

    int TotalToMove = 0;
    for (proc = 0; proc < NumberOfProcessors; proc++)
      TotalToMove += NumberToMove[proc];
    int NumberToMoveThisGrid = TotalToMove - PreviousTotalToMove;
    int NumberLeft = NumberOfActiveParticles - NumberToMoveThisGrid;

    if (NumberToMoveThisGrid > 0) {
 
      // Increase the size of the list to include the active particles
      // from this grid

      ActiveParticleType **NewList = new ActiveParticleType*[TotalToMove]();
      for (i = 0; i < PreviousTotalToMove; i++)
	NewList[i] = List[i];
      delete [] List;
      List = NewList;
 
      /* Move active particles into list */

      int n1 = PreviousTotalToMove;
      int index = 0;
      ActiveParticleType **OldActiveParticles = ActiveParticles;
      ActiveParticles = new ActiveParticleType*[NumberLeft]();

      for (i = 0; i < NumberOfActiveParticles; i++) {
	grid = ToGrid[i];
	if (grid != ThisGridNum) {
	  List[n1] = OldActiveParticles[i];
	  List[n1]->SetDestProcessor(MyProcessorNumber);
	  List[n1]->SetGridID(grid);
	  n1++;
	} else {
	  ActiveParticles[index++] = OldActiveParticles[i];
	}
      } // ENDFOR particles

      this->NumberOfActiveParticles = NumberLeft;
      if (NumberLeft == 0) ActiveParticles = NULL;
      
    } // ENDIF NumberToMoveThisProc > 0

    delete [] ToGrid;
 
  } // end: if (COPY_OUT)
 
  /* ----------------------------------------------------------------- */
  /* Copy active particles back into grid. */
 
  else {

    /* Count up total number. */
 
    int TotalNumberOfActiveParticles;
    int NumberOfNewActiveParticles = EndIndex - StartIndex;

    TotalNumberOfActiveParticles = NumberOfActiveParticles + NumberOfNewActiveParticles;
 
    /* Copy active particles from buffer */

    if (NumberOfNewActiveParticles > 0) {
      for (i = 0; i < NumberOfNewActiveParticles; i++) {
        this->AddActiveParticle(List[i]);
      }
    }


  } // end: if (COPY_IN)
 
  return SUCCESS;
}
