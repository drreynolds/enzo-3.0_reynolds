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

class ActiveParticleType_Kravtsov;

class KravtsovBufferHandler : public ParticleBufferHandler
{
  public:
  // No extra fields in CenOstriker.  Same base constructor.
  KravtsovBufferHandler(void) : ParticleBufferHandler() {};

  KravtsovBufferHandler(int NumberOfParticles) :
    ParticleBufferHandler(NumberOfParticles) {};

  KravtsovBufferHandler(ActiveParticleType **np, int NumberOfParticles, int type, int proc) :
    ParticleBufferHandler(np, NumberOfParticles, type, proc) {};

  static void AllocateBuffer(ActiveParticleType **np, int NumberOfParticles, char *buffer,
			     Eint32 total_buffer_size, int &buffer_size, Eint32 &position,
			     int type_num, int proc);
  static void UnpackBuffer(char *mpi_buffer, int mpi_buffer_size, int NumberOfParticles,
			   ActiveParticleType **np, int &npart);
  static int ReturnHeaderSize(void) {return HeaderSizeInBytes; }
  static int ReturnElementSize(void) {return ElementSizeInBytes; }

  // Extra fields would go here.  See
  // ActiveParticle_AccretingParticle.C for a particle type that uses
  // extra fields.
};

class ActiveParticleType_Kravtsov : public ActiveParticleType
{
public:
  // Constructors
  ActiveParticleType_Kravtsov(void) : ActiveParticleType() {};

  ActiveParticleType_Kravtsov(KravtsovBufferHandler *buffer, int index) :
    ActiveParticleType(static_cast<ParticleBufferHandler*>(buffer), index) {};
  
  // Static members
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &supp_data);
  static int WriteToOutput(ActiveParticleType **these_particles, int n, int GridRank, hid_t group_id);
  static int ReadFromOutput(ActiveParticleType **&particles_to_read, int &n, int GridRank, hid_t group_id);
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);
  static int EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static int BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
			       int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			       int ThisLevel, int TotalStarParticleCountPrevious[],
			       int SampleParticleID);
  static int AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
			      int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			      int ThisLevel, int TotalStarParticleCountPrevious[],
			      int SampleParticleID);
  static int SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level, int TopGridDims[], int ActiveParticleID);
  static int InitializeParticleType(void);
  
  // Need this to make active particle ID work correctly.
  int GetEnabledParticleID(int myid = -1) {				
    static int ParticleID = -1;						
    if (myid >= 0) {							
      if (ParticleID != -1) ENZO_FAIL("Setting Particle ID Twice!");	
      ParticleID = myid;						
    }									
    return ParticleID;							
  };

  static float DensityThreshold, StarFormationTimeConstant, MinimumStarMass;

};
