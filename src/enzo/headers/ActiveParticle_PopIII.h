/***********************************************************************
/
/  A particle that represents a Pop III star.
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
#include "units.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CommunicationUtilities.h"
#include "phys_constants.h"
#include "FofLib.h"

class ActiveParticleType_PopIII;


class ActiveParticleType_PopIII : public ActiveParticleType
{
public:
  // Constructors
  ActiveParticleType_PopIII(void) : ActiveParticleType() {
    Lifetime = 0; 
  };
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &supp_data);
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);
  static int EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data);
  template <class active_particle_class>
    static int BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				 int ThisLevel, bool CallEvolvePhotons, 
				 int TotalStarParticleCountPrevious[],
				 int PopIIIParticleID) { return SUCCESS; };
  template <class active_particle_class>
    static int AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int TotalStarParticleCountPrevious[],
				int PopIIIParticleID) {return SUCCESS; };
  template <class active_particle_class>
    static int DepositMass(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int GalaxyParticleID) {return SUCCESS; };
  static int SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level, 
			      int TopGridDims[], int ActiveParticleID);
  static int InitializeParticleType();
  bool IsARadiationSource(FLOAT Time);
  ENABLED_PARTICLE_ID_ACCESSOR

  // Pop III specific active particle parameters
  static float OverDensityThreshold, MetalCriticalFraction, 
    H2CriticalFraction, StarMass;

  float Lifetime;

  static std::vector<ParticleAttributeHandler *> AttributeHandlers;
};
