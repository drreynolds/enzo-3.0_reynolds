/***********************************************************************
/
/  Construct a fake grid for feedback algorithms based on a 
/  list of active particles
/
/  written by: Nathan Goldbaum
/  date:       June 2012
/
************************************************************************/

#ifdef USE_MPI
#include "communicators.h"
#endif

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
#include "CommunicationUtilities.h"
#include "communication.h"

int CommunicationBufferPurge(void);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);



grid** ConstructFeedbackZone(ActiveParticleType** ParticleList, int nParticles, int FeedbackRadius, 
			     FLOAT dx, HierarchyEntry** Grids, int NumberOfGrids)
{
  int i,j,dim,size;
  int FeedbackZoneRank, FeedbackZoneDimension[FeedbackZoneRank];
  FLOAT FeedbackZoneLeftEdge[FeedbackZoneRank], FeedbackZoneRightEdge[FeedbackZoneRank];
  FLOAT CellSize = CellWidth[0][0], ncells[FeedbackZoneRank];

  grid** sinkGrids = new grid*[nParticles]();

  for (i = 0; i < nParticles; i++) {
    sinkGrid[i] = ParticleList[i]->ReturnCurrentGrid();
    if (sinkGrid == NULL)
      ENZO_FAIL("Particle CurrentGrid is invalid!");

    FLOAT* ParticlePosition = ThisParticle->ReturnPosition();
    
    // This should only happen if the grid pointer is invalid
    if ((GridLeftEdge[0] > ParticlePosition[0]+FeedbackRadius) || 
	(GridLeftEdge[1] > ParticlePosition[1]+FeedbackRadius) || 
	(GridLeftEdge[2] > ParticlePosition[2]+FeedbackRadius) || 
	(GridRightEdge[0] < ParticlePosition[0]-FeedbackRadius) ||
	(GridRightEdge[1] < ParticlePosition[1]-FeedbackRadius) ||
	(GridRightEdge[2] < ParticlePosition[2]-FeedbackRadius) ||)
      ENZO_FAIL("Particle outside own grid!");
  }

  /* Setup grid properties */

  grid** FeedbackZones = new grid*[nParticles];
  
  for (i = 0; i < nParticles; i++) {
    FeedbackZoneRank = this->GetGridRank();
    size = 1;

    for (int dim = 0; dim < FeedbackZoneRank; dim++) {
      FeedbackZoneDimension[dim] = (2*(FeedbackRadius+DEFAULT_GHOST_ZONES)+1);
      size *= FeedbackZoneDimension[dim];
      
      FeedbackZoneLeftEdge[dim]  = sinkGrid+CellSize*(ncells[dim]-FeedbackRadius);
      FeedbackZoneRightEdge[dim] = +CellSize*(ncells[dim]+FeedbackRadius+1);
    }
    
    grid *FeedbackZone = new grid;
    
    FeedbackZone->InheritProperties(this);
    
    FeedbackZone->PrepareGrid(FeedbackZoneRank, FeedbackZoneDimension, 
			      FeedbackZoneLeftEdge,FeedbackZoneRightEdge,0);
    
    FeedbackZone->SetProcessorNumber(MyProcessorNumber);
    
    FeedbackZone->SetTimeStep(this->ReturnTimeStep());
        
    if (FeedbackZone->AllocateAndZeroBaryonField() == FAIL)
      ENZO_FAIL("FeedbackZone BaryonField allocation failed");
  }

  /* Post receives */

  CommunicationReceiveIndex = 0;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  CommunicationDirection = COMMUNICATION_POST_RECEIVE;

  for (i = 0; i < nParticles; i++) {
    for (j = 

  }

  // Copy zones from this grid (which must overlap the position of the AP).
  // Note, using ZeroVector here will break if a FeedbackZone overlaps with a
  // domain boundary
  float ZeroVector[] = {0,0,0};
  if (FeedbackZone->CopyZonesFromGrid(this,ZeroVector) == FAIL)
    ENZO_FAIL("FeedbackZone copy failed!");
  
  // if the grid is filled, return
  
  
    // Next, recursively iterate over the siblings of that grid, copying
    // zones from overlapping grids until the grid is filled
  }

  for (i = 0; i < nParticles; i++) {
    FeedbackZones[i] = FeedbackZone;
  }

  return FeedbackZones;
  
}
