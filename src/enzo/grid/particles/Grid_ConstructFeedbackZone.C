/***********************************************************************
/
/  Construct a fake grid for feedback algorithms based on a 
/  list of active particles
/
/  written by: Nathan Goldbaum
/  date:       June 2012
/
************************************************************************/

#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "phys_constants.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "ActiveParticle.h"
#include "phys_constants.h"

grid** ConstructFeedbackZone(ActiveParticleType** ParticleList, int FeedbackRadius, FLOAT dx,
			     HierarchyEntry** Grids)
{

  for (i = 0; i < nParticles; i++) {
    sinkGrid = ParticleList[i]->ReturnCurrentGrid();
    if (sinkGrid == NULL) {
      return FAIL;
    }
  }

  FLOAT* ParticlePosition = ThisParticle->ReturnPosition();

  // This should only happen if the grid pointer is invalid
  if ((GridLeftEdge[0] > ParticlePosition[0]+FeedbackRadius) || 
      (GridLeftEdge[1] > ParticlePosition[1]+FeedbackRadius) || 
      (GridLeftEdge[2] > ParticlePosition[2]+FeedbackRadius) || 
      (GridRightEdge[0] < ParticlePosition[0]-FeedbackRadius) ||
      (GridRightEdge[1] < ParticlePosition[1]-FeedbackRadius) ||
      (GridRightEdge[2] < ParticlePosition[2]-FeedbackRadius))
    ENZO_FAIL("Particle outside own grid!");

  /* Setup grid properties */

  int FeedbackZoneRank = this->GetGridRank();

  int FeedbackZoneDimension[FeedbackZoneRank], size=1;
  FLOAT LeftCellOffset[FeedbackZoneRank], FeedbackZoneLeftEdge[FeedbackZoneRank], 
    FeedbackZoneRightEdge[FeedbackZoneRank];
  FLOAT CellSize = CellWidth[0][0], ncells[FeedbackZoneRank];

  for (int i = 0; i < FeedbackZoneRank; i++) {
    FeedbackZoneDimension[i] = (2*(FeedbackRadius+DEFAULT_GHOST_ZONES)+1);
    size *= FeedbackZoneDimension[i];

    LeftCellOffset[i]        = modf((ParticlePosition[i]-CellLeftEdge[i][0])/CellSize,&ncells[i]);
        
    FeedbackZoneLeftEdge[i]  = CellLeftEdge[0][0]+CellSize*(ncells[i]-FeedbackRadius);
    FeedbackZoneRightEdge[i] = CellLeftEdge[0][0]+CellSize*(ncells[i]+FeedbackRadius+1);
  }

  grid *FeedbackZone = new grid;
  
  FeedbackZone->InheritProperties(this);
  
  FeedbackZone->PrepareGrid(FeedbackZoneRank, FeedbackZoneDimension, 
			    FeedbackZoneLeftEdge,FeedbackZoneRightEdge,0);
  
  FeedbackZone->SetProcessorNumber(MyProcessorNumber);

  FeedbackZone->SetTimeStep(this->ReturnTimeStep());

  

  if (FeedbackZone->AllocateAndZeroBaryonField() == FAIL)
    ENZO_FAIL("FeedbackZone BaryonField allocation failed");

  // Copy zones from this grid (which must overlap the position of the AP).
  // Note, using ZeroVector here will break if a FeedbackZone overlaps with a
  // domain boundary
  float ZeroVector[] = {0,0,0};
  if (FeedbackZone->CopyZonesFromGrid(this,ZeroVector) == FAIL)
    ENZO_FAIL("FeedbackZone copy failed!");

  // if the grid is filled, return
  

  // Next, recursively iterate over the siblings of that grid, copying
  // zones from overlapping grids until the grid is filled

  delete [] FlaggingField;

  return FeedbackZone;
  
}
