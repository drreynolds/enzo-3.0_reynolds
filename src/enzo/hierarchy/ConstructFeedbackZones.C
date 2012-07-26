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
  int FeedbackZoneRank; 
  FLOAT** ParticlePosition = new FLOAT*[nParticles]();

  /* Build array of AP grids and check for errors */
  grid** APGrids = new grid*[nParticles];

  for (i = 0; i < nParticles; i++) {
    APGrids[i] = ParticleList[i]->ReturnCurrentGrid();
    if (APGrids[i] == NULL)
      ENZO_FAIL("Particle CurrentGrid is invalid!\n");

    ParticlePosition[i] = ParticleList[i]->ReturnPosition();
    
    // This should only happen if the grid pointer is invalid
    if ((APGrids[i]->GetGridLeftEdge(0) > ParticlePosition[i][0]+FeedbackRadius) || 
	(APGrids[i]->GetGridLeftEdge(1) > ParticlePosition[i][1]+FeedbackRadius) || 
	(APGrids[i]->GetGridLeftEdge(2) > ParticlePosition[i][2]+FeedbackRadius) || 
	(APGrids[i]->GetGridRightEdge(0) < ParticlePosition[i][0]-FeedbackRadius) ||
	(APGrids[i]->GetGridRightEdge(1) < ParticlePosition[i][1]-FeedbackRadius) ||
	(APGrids[i]->GetGridRightEdge(2) < ParticlePosition[i][2]-FeedbackRadius))
      ENZO_FAIL("Particle outside own grid!\n");
  }

  /* Setup Feedback Zones before copying data */

  grid** FeedbackZones = new grid*[nParticles];
  
  for (i = 0; i < nParticles; i++) {
    FeedbackZoneRank = APGrids[i]->GetGridRank();

    int FeedbackZoneDimension[FeedbackZoneRank];
    FLOAT LeftCellOffset[FeedbackZoneRank],FeedbackZoneLeftEdge[FeedbackZoneRank], 
      FeedbackZoneRightEdge[FeedbackZoneRank];
    FLOAT CellSize, GridGZLeftEdge, ncells[FeedbackZoneRank];

    size = 1;

    for (int dim = 0; dim < FeedbackZoneRank; dim++) {
      FeedbackZoneDimension[dim] = (2*(FeedbackRadius+DEFAULT_GHOST_ZONES)+1);
      size *= FeedbackZoneDimension[dim];
      CellSize = APGrids[i]->GetCellWidth(dim,0);
      GridGZLeftEdge = APGrids[i]->GetCellLeftEdge(dim,0);
      
      LeftCellOffset[dim] = modf((ParticlePosition[i][dim]-GridGZLeftEdge)/CellSize,&ncells[dim]);

      FeedbackZoneLeftEdge[dim]  = GridGZLeftEdge + CellSize*(ncells[dim]-FeedbackRadius);
      FeedbackZoneRightEdge[dim] = GridGZLeftEdge + CellSize*(ncells[dim]+FeedbackRadius+1);
    }
    
    grid *FeedbackZone = new grid;
    
    FeedbackZone->InheritProperties(APGrids[i]);
    
    FeedbackZone->PrepareGrid(FeedbackZoneRank, FeedbackZoneDimension, 
			      FeedbackZoneLeftEdge,FeedbackZoneRightEdge,0);
    
    FeedbackZone->SetProcessorNumber(APGrids[i]->ReturnProcessorNumber());
    
    FeedbackZone->SetTimeStep(APGrids[i]->ReturnTimeStep());
        
    // This will only allocate the BaryonField on the host processor
    if (FeedbackZone->AllocateAndZeroBaryonField() == FAIL)
      ENZO_FAIL("FeedbackZone BaryonField allocation failed\n");

    FeedbackZones[i] = FeedbackZone;
  }

  delete [] APGrids;
  delete [] ParticlePosition; 

  // Copy zones from this grid (which must overlap the position of the AP).
  // Note, using ZeroVector here will break if a FeedbackZone overlaps with a
  // domain boundary
  float ZeroVector[] = {0,0,0};

  /* Post receives */

  CommunicationReceiveIndex = 0;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  CommunicationDirection = COMMUNICATION_POST_RECEIVE;

  for (i = 0; i < nParticles; i++) 
    for (j = 0; j < NumberOfGrids; j++) 
      if (FeedbackZones[i]->CopyZonesFromGrid(Grids[j]->GridData,ZeroVector) == FAIL)
	ENZO_FAIL("FeedbackZone copy failed!\n");
    
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
