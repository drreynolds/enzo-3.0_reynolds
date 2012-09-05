/***********************************************************************
/
/  GALAXY PARTICLE TYPE
/
/  written by: Stephen Skory
/  date:       August, 2012
/
/  PURPOSE:
/
************************************************************************/

#include "preincludes.h"
#include "hdf5.h"
#include "h5utilities.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "EventHooks.h"
#include "ActiveParticle.h"

class ActiveParticleType_GalaxyParticle : public ActiveParticleType
{
public:
  // Constructors
  ActiveParticleType_GalaxyParticle(void) : ActiveParticleType() {
    Radius = 0;
    initialized = 0;
  };
  ActiveParticleType_GalaxyParticle(ActiveParticleType_GalaxyParticle* part) :
    ActiveParticleType(static_cast<ActiveParticleType*>(part)) {
    Radius = part->Radius;
    initialized = part->initialized;
  };

  // Static members
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);
  static int EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data);
  template <class active_particle_class>
    static int BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				 int ThisLevel, int TotalStarParticleCountPrevious[],
				 int GalaxyParticleID);
  template <class active_particle_class>
    static int AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int TotalStarParticleCountPrevious[],
				int GalaxyParticleID);
  static int SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level, int TopGridDims[], int ActiveParticleID);
  static int InitializeParticleType(void);
  
  static std::vector<ParticleAttributeHandler *> AttributeHandlers;

  // Need this to make active particle ID work correctly.
  int GetEnabledParticleID(int myid = -1) {				
    static int ParticleID = -1;						
    if (myid >= 0) {							
      if (ParticleID != -1) ENZO_FAIL("Setting Particle ID Twice!");	
      ParticleID = myid;						
    }									
    return ParticleID;							
  };
  
  // Galaxy Particle specific stuff.
  float Radius;
  int initialized; // Has mass been subtracted from the grid?
};

template <class active_particle_class>
int ActiveParticleType_GalaxyParticle::BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
							    int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
							    int ThisLevel, int TotalStarParticleCountPrevious[],
							    int GalaxyParticleID)
{

  return SUCCESS;
}

template <class active_particle_class>
int ActiveParticleType_GalaxyParticle::AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
							   int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
							   int ThisLevel, int TotalStarParticleCountPrevious[],
							   int GalaxyParticleID)
{

  return SUCCESS;

}
