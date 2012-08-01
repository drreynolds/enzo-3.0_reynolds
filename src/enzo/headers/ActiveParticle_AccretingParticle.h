/***********************************************************************
/
/ Accreting Particle
/
************************************************************************/

#ifdef USE_MPI
#include "communicators.h"
#endif 

#include "preincludes.h"
#include "hdf5.h"
#include "h5utilities.h"

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
#include "communication.h"
#include "phys_constants.h"
#include "FofLib.h"

class ActiveParticleType_AccretingParticle;

class AccretingParticleBufferHandler : public ParticleBufferHandler
{
public:
  AccretingParticleBufferHandler(void) : ParticleBufferHandler() {
    this->CalculateAccretingParticleElementSize();
};
  AccretingParticleBufferHandler(int NumberOfParticles) : ParticleBufferHandler(NumberOfParticles) {
    // Any extra fields must be added to the buffer
    if (this->NumberOfBuffers > 0) {
      this->AccretionRate = new float[NumberOfParticles];
      this->CalculateAccretingParticleElementSize();
    }
  };
  AccretingParticleBufferHandler(ActiveParticleType **np, int NumberOfParticles, int type, int proc);
  ~AccretingParticleBufferHandler() {
    if (this->NumberOfBuffers > 0) {
      delete [] AccretionRate;
    }
  };
  static void AllocateBuffer(ActiveParticleType **np, int NumberOfParticles, char *buffer, 
			     Eint32 total_buffer_size, int &buffer_size, Eint32 &position, 
			     int type_num, int proc);
  static void UnpackBuffer(char *mpi_buffer, int mpi_buffer_size, int NumberOfParticles,
			   ActiveParticleType **np, int &npart);
  static int ReturnHeaderSize(void) {return HeaderSizeInBytes; }
  static int ReturnElementSize(void) {return ElementSizeInBytes; }
  void CalculateAccretingParticleElementSize(void) {
    Eint32 mpi_flag = 0;
#ifdef USE_MPI
    MPI_Initialized(&mpi_flag);
#endif
    Eint32 size;
    if (mpi_flag == 1) {
#ifdef USE_MPI
      // float: 1 -- AccretionRate
      MPI_Pack_size(1, FloatDataType, EnzoTopComm, &size);
      this->ElementSizeInBytes += size;
#endif
    } else
      this->ElementSizeInBytes += 1*sizeof(float);
  };
  float* AccretionRate;
};

class ActiveParticleType_AccretingParticle : public ActiveParticleType
{
public:
  // Constructors
  ActiveParticleType_AccretingParticle(void) : ActiveParticleType() {
    AccretionRate = 0;
  };
  ActiveParticleType_AccretingParticle(AccretingParticleBufferHandler *buffer, int index) :
    ActiveParticleType(static_cast<ParticleBufferHandler*>(buffer), index) {
    AccretionRate = buffer->AccretionRate[index];
  };
  ActiveParticleType_AccretingParticle(ActiveParticleType_AccretingParticle* part) :
    ActiveParticleType(static_cast<ActiveParticleType*>(part)) {
    AccretionRate = part->AccretionRate;
  };
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static int WriteToOutput(ActiveParticleType **these_particles, int n, int GridRank, hid_t group_id);
  static int ReadFromOutput(ActiveParticleType **&particles_to_read, int &n, int GridRank, hid_t group_id);
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);
  static int EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static int BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
			       int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			       int ThisLevel, int TotalStarParticleCountPrevious[],
			       int AccretingParticleID);
  static int AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
			      int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			      int ThisLevel, int TotalStarParticleCountPrevious[],
			      int AccretingParticleID);
  static int SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level, int TopGridDims[], int ActiveParticleID);
  static int InitializeParticleType();
  
  int GetEnabledParticleID(int myid = -1) {				
    static int ParticleID = -1;						
    if (myid >= 0) {							
      if (ParticleID != -1) ENZO_FAIL("Setting Particle ID Twice!");	
      ParticleID = myid;						
    }									
    return ParticleID;							
  };

  // sink helper routines

  static ActiveParticleType_AccretingParticle**  MergeAccretingParticles
  (int *nParticles, ActiveParticleType** ParticleList, FLOAT LinkingLength, 
   int *ngroups, LevelHierarchyEntry *LevelArray[], int ThisLevel);

  static int Accrete(int nParticles, ActiveParticleType** ParticleList,
		     int AccretionRadius, FLOAT dx, 
		     LevelHierarchyEntry *LevelArray[], int ThisLevel);

  static int AccreteOntoAccretingParticle(grid* AccretionZone, 
					  ActiveParticleType_AccretingParticle* ThisParticle,
					  FLOAT AccretionRadius);

  static float OverflowFactor;
  static int AccretionRadius;   // in units of CellWidth on the maximum refinement level
  static int LinkingLength;     // Should be equal to AccretionRadius
  float AccretionRate;
  static std::vector<ParticleAttributeHandler> AttributeHandlers;
};

