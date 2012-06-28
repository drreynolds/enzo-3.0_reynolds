/***********************************************************************
/
/  GRID CLASS (Construct a fake grid for feedback algorithms)
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

int grid::ConstructFeedbackZone(ActiveParticleType* ThisParticle, int FeedbackRadius, FLOAT dx, grid* FeedbackZone)
{
  FLOAT* ParticlePosition = ThisParticle->ReturnPosition();

  // This should only happen if the grid pointer is invalid
  if ((GridLeftEdge[0] > ParticlePosition[0]+FeedbackRadius) || (GridRightEdge[0] < ParticlePosition[0]-FeedbackRadius) ||
      (GridLeftEdge[1] > ParticlePosition[1]+FeedbackRadius) || (GridRightEdge[1] < ParticlePosition[1]-FeedbackRadius) ||
      (GridLeftEdge[2] > ParticlePosition[2]+FeedbackRadius) || (GridRightEdge[2] < ParticlePosition[2]-FeedbackRadius))
    return FAIL;

  /* Setup grid properties */

  int FeedbackZoneRank = FeedbackZone->GetGridRank();

  // Since the grid creation machinery assume we want ghost zones, we need to
  // trick it into giving us a fake grid without ghost zones.  We therefore ask
  // for a cubic grid of dimension (2*(FeedbackRadius-3)+1)^3
  
  int FeedbackZoneDimension[FeedbackZoneRank], size;
  FLOAT LeftCellOffset[FeedbackZoneRank],RightCellOffset[FeedbackZoneRank],
    FeedbackZoneLeftEdge[FeedbackZoneRank], FeedbackZoneRightEdge[FeedbackZoneRank];

  for (int i = 0; i < FeedbackZoneRank; i++) {
    if (FeedbackRadius > DEFAULT_GHOST_ZONES)
      FeedbackZoneDimension[i] = 2*FeedbackRadius+1;
    else
      FeedbackZoneDimension[i] = 2*DEFAULT_GHOST_ZONES+1;
    size *= FeedbackZoneDimension[i];

    LeftCellOffset[i]        = fmod(ParticlePosition[i],dx);
    RightCellOffset[i]       = dx-LeftCellOffset[i];
    
    FeedbackZoneLeftEdge[i]  = ParticlePosition[i]-FeedbackRadius*dx-LeftCellOffset[i];
    FeedbackZoneRightEdge[i] = ParticlePosition[i]+FeedbackRadius*dx+RightCellOffset[i];
  }

  /* Intialize the fake grid */
  FeedbackZone = new grid;
  
  FeedbackZone->InheritProperties(this);
  
  FeedbackZone->PrepareGrid(FeedbackZoneRank, FeedbackZoneDimension, 
			    FeedbackZoneLeftEdge,FeedbackZoneRightEdge,0);
  
  FeedbackZone->SetProcessorNumber(MyProcessorNumber);

  FeedbackZone->AllocateAndZeroBaryonField();
    
  // Allocate flagging field of the same size as BaryonField. If FlaggingField =
  // 0, the corresponding zone in the BaryonField has not been copied yet.
    
  int* FlaggingField = new int[size];
  for (int i = 0; i<size; i++)
    FlaggingField[i] = 0;

  // Copy zones from this grid (which must overlap the position of the AP).
  // Note, using ZeroVector here will break if a FeedbackZone overlaps with a
  // domain boundary
  float ZeroVector[] = {0,0,0};
  FeedbackZone->CopyZonesFromGrid(this,ZeroVector);

  // if the grid is filled, return
  

  // Next, recursively iterate over the siblings of that grid, copying
  // zones from overlapping grids until the grid is filled

  delete [] FlaggingField;

  return SUCCESS;
  
}
