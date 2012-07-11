/***********************************************************************
/
/  GRID CLASS (Copy data from a 'fake' feedback zone grid back to 
/              the real grids)
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

int grid::DistributeFeedbackZone(ActiveParticleType* ThisParticle, FLOAT FeedbackRadius)
{

  FLOAT* ParticlePosition = ThisParticle->ReturnPosition();

  // This should only happen if enzo is corrupted
  if ((GridLeftEdge[0] > ParticlePosition[0]+FeedbackRadius) || (GridRightEdge[0] < ParticlePosition[0]-FeedbackRadius) ||
      (GridLeftEdge[1] > ParticlePosition[1]+FeedbackRadius) || (GridRightEdge[1] < ParticlePosition[1]-FeedbackRadius) ||
      (GridLeftEdge[2] > ParticlePosition[2]+FeedbackRadius) || (GridRightEdge[2] < ParticlePosition[2]-FeedbackRadius))
    return FAIL;

  grid* APGrid = ThisParticle->ReturnCurrentGrid();

  // Copy zones to the AP's grid.  Note, using ZeroVector here will break if a
  // FeedbackZone overlaps with a domain boundary
  float ZeroVector[] = {0,0,0};
  APGrid->CopyZonesFromGrid(this,ZeroVector);

}
