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



grid** ConstructFeedbackZones(ActiveParticleType** ParticleList, int nParticles, int FeedbackRadius, 
			     FLOAT dx, HierarchyEntry** Grids, int NumberOfGrids)
{
  int i,j,dim,size;
  int FeedbackZoneRank, FeedbackZoneDimension[FeedbackZoneRank];
  FLOAT LeftCellOffset[FeedbackZoneRank], *ParticlePosition=NULL,
    FeedbackZoneLeftEdge[FeedbackZoneRank], FeedbackZoneRightEdge[FeedbackZoneRank];
  FLOAT CellSize, GridGZLeftEdge, ncells[FeedbackZoneRank];

  /* Build array of sink grids and check for errors */
  grid** sinkGrids = new grid*[nParticles]();

  for (i = 0; i < nParticles; i++) {
    sinkGrids[i] = ParticleList[i]->ReturnCurrentGrid();
    if (sinkGrids[i] == NULL)
      ENZO_FAIL("Particle CurrentGrid is invalid!\n");

    ParticlePosition = ParticleList[i]->ReturnPosition();
    
    // This should only happen if the grid pointer is invalid
    if ((sinkGrids[i]->GetGridLeftEdge(0) > ParticlePosition[0]+FeedbackRadius) || 
	(sinkGrids[i]->GetGridLeftEdge(1) > ParticlePosition[1]+FeedbackRadius) || 
	(sinkGrids[i]->GetGridLeftEdge(2) > ParticlePosition[2]+FeedbackRadius) || 
	(sinkGrids[i]->GetGridRightEdge(0) < ParticlePosition[0]-FeedbackRadius) ||
	(sinkGrids[i]->GetGridRightEdge(1) < ParticlePosition[1]-FeedbackRadius) ||
	(sinkGrids[i]->GetGridRightEdge(2) < ParticlePosition[2]-FeedbackRadius))
      ENZO_FAIL("Particle outside own grid!\n");
  }

  /* Setup Feedback Zones before copying data */

  grid** FeedbackZones = new grid*[nParticles];
  
  for (i = 0; i < nParticles; i++) {
    FeedbackZoneRank = sinkGrids[i]->GetGridRank();
    size = 1;

    for (int dim = 0; dim < FeedbackZoneRank; dim++) {
      FeedbackZoneDimension[dim] = (2*(FeedbackRadius+DEFAULT_GHOST_ZONES)+1);
      size *= FeedbackZoneDimension[dim];
      CellSize = sinkGrids[i]->GetCellWidth(dim,0);
      GridGZLeftEdge = sinkGrids[i]->GetCellLeftEdge(dim,0);
      
      LeftCellOffset[i] = modf((ParticlePosition[i]-GridGZLeftEdge)/CellSize,&ncells[i]);

      FeedbackZoneLeftEdge[dim]  = GridGZLeftEdge + CellSize*(ncells[dim]-FeedbackRadius);
      FeedbackZoneRightEdge[dim] = GridGZLeftEdge + CellSize*(ncells[dim]+FeedbackRadius+1);
    }
    
    grid *FeedbackZone = new grid;
    
    FeedbackZone->InheritProperties(sinkGrids[i]);
    
    FeedbackZone->PrepareGrid(FeedbackZoneRank, FeedbackZoneDimension, 
			      FeedbackZoneLeftEdge,FeedbackZoneRightEdge,0);
    
    FeedbackZone->SetProcessorNumber(MyProcessorNumber);
    
    FeedbackZone->SetTimeStep(sinkGrids[i]->ReturnTimeStep());
        
    // This will only allocate the BaryonField on the host processor
    if (FeedbackZone->AllocateAndZeroBaryonField() == FAIL)
      ENZO_FAIL("FeedbackZone BaryonField allocation failed\n");
  }

  // Copy zones from this grid (which must overlap the position of the AP).
  // Note, using ZeroVector here will break if a FeedbackZone overlaps with a
  // domain boundary
  float ZeroVector[] = {0,0,0};

  /* Post receives */

  CommunicationReceiveIndex = 0;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  CommunicationDirection = COMMUNICATION_POST_RECEIVE;

  for (i = 0; i < nParticles; i++) {
    for (j = 0; j < NumberOfGrids; j++) {
      if (FeedbackZones[i]->CopyZonesFromGrid(Grids[j]->GridData,ZeroVector) == FAIL)
	ENZO_FAIL("FeedbackZone copy failed!\n");
    }
  }

  /* Send data */

  CommunicationDirection = COMMUNICATION_SEND;

  for (i = 0; i < nParticles; i++) {
    for (j = 0; j < NumberOfGrids; j++) {
      if (FeedbackZones[i]->CopyZonesFromGrid(Grids[j]->GridData,ZeroVector) == FAIL)
	ENZO_FAIL("FeedbackZone copy failed!\n");
    }
  }

  /* Receive data */

  if (CommunicationReceiveHandler() == FAIL)
    ENZO_FAIL("CommunicationReceiveHandler() failed!\n");

#ifdef USE_MPI
  CommunicationBufferPurge();
#endif

  return FeedbackZones;
  
}
