/***********************************************************************
/
/  GRID CLASS (Calculate the accretion rate and subtract accreted mass 
/              from the grid.)
/
/  written by: Nathan Goldbaum
/  date:       April 2012
/
/  note:       Equation numbers refer to Krumholz McKee & Klein (2004)
************************************************************************/

#include <string.h>
#include <map>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <iostream>

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

  FeedbackZone = new grid;
  
  FeedbackZone->InheritProperties(this);
  
  FLOAT LeftCellOffset[3]  = {fmod(ParticlePosition[0],dx),
			      fmod(ParticlePosition[1],dx),
			      fmod(ParticlePosition[2],dx)};
  
  FLOAT RightCellOffset[3] = {dx-LeftCellOffset[0],     
			      dx-LeftCellOffset[1],    
			      dx-LeftCellOffset[2]};
  
  FLOAT FeedbackZoneLeftEdge[3]  = {ParticlePosition[0]-FeedbackRadius*dx-LeftCellOffset[0], 
				    ParticlePosition[1]-FeedbackRadius*dx-LeftCellOffset[1], 
				    ParticlePosition[2]-FeedbackRadius*dx-LeftCellOffset[2]};
  
  FLOAT FeedbackZoneRightEdge[3] = {ParticlePosition[0]+FeedbackRadius*dx+RightCellOffset[0], 
				    ParticlePosition[1]+FeedbackRadius*dx+RightCellOffset[1], 
				    ParticlePosition[2]+FeedbackRadius*dx+RightCellOffset[2]};
  
  int FeedbackZoneRank = FeedbackZone->GetGridRank();
  
  int FeedbackZoneDimension[3] = {FeedbackRadius*2+1, 
				  FeedbackRadius*2+1, 
				  FeedbackRadius*2+1};
  
  FeedbackZone->PrepareGrid(FeedbackZoneRank, FeedbackZoneDimension, 
			    FeedbackZoneLeftEdge,FeedbackZoneRightEdge,0);
  
  FeedbackZone->SetProcessorNumber(MyProcessorNumber);
  
  // First, copy zones from grid that covers the sink particle

  // if the grid is filled, return

  // Next, recursively iterate over the siblings of that grid, copying zones from overlapping grids until the grid is filled

  return SUCCESS;
}
