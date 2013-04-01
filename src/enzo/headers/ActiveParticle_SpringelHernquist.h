/***********************************************************************
/
/  SPRINGEL & HERNQUIST STAR FORMATION
/
/  written by: Stephen Skory
/  date:       October, 2011
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

#include "phys_constants.h"

class ActiveParticleType_SpringelHernquist : public ActiveParticleType
{
public:
  // Constructors
  ActiveParticleType_SpringelHernquist(void) : ActiveParticleType() {};

  // Static members
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);
  static int EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data);
  template <class active_particle_class>
    static int BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				 int ThisLevel, int TotalStarParticleCountPrevious[],
				 int SpringelHernquistID) { return SUCCESS; };
  template <class active_particle_class>
    static int AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int TotalStarParticleCountPrevious[],
				int SpringelHernquistID) { return SUCCESS; };
  template <class active_particle_class>
    static int DepositMass(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int GalaxyParticleID) { return SUCCESS; };
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
  static float OverDensityThreshold, PhysicalDensityThreshold, 
    MinimumDynamicalTime, MinimumMass;
};
