/***********************************************************************
/
/  AN EXAMPLE ACTIVE PARTICLE TYPE
/
/  written by: Matthew Turk
/  date:       May, 2011
/
/  PURPOSE:
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
#include "phys_constants.h"
#include "FofLib.h"

class ActiveParticleType_PopIII;

class PopIIIParticleBufferHandler : public ParticleBufferHandler
{
public:
  PopIIIParticleBufferHandler(void) : ParticleBufferHandler() { this->CalculatePopIIIParticleElementSize(); };

  PopIIIParticleBufferHandler(int NumberOfParticles) : ParticleBufferHandler(NumberOfParticles) {
    if (this->NumberOfBuffers > 0)
      Lifetime = new float[NumberOfParticles];
  };
  PopIIIParticleBufferHandler(ActiveParticleType **np, int NumberOfParticles, int type, int proc);
  ~PopIIIParticleBufferHandler() {
    if (this->NumberOfBuffers > 0)
      delete [] Lifetime;
  };
  static void AllocateBuffer(ActiveParticleType **np, int NumberOfParticles, char *buffer, 
			     Eint32 total_buffer_size, int &buffer_size, Eint32 &position, 
			     int type_num, int proc);
  static void UnpackBuffer(char *mpi_buffer, int mpi_buffer_size, int NumberOfParticles,
			   ActiveParticleType **np, int &npart);
  static int ReturnHeaderSize(void) {return HeaderSizeInBytes; }
  static int ReturnElementSize(void) {return ElementSizeInBytes; }
  void CalculatePopIIIParticleElementSize(void) {
    Eint32 mpi_flag = 0;
#ifdef USE_MPI
    MPI_Initialized(&mpi_flag);
#endif
    Eint32 size;
    if (mpi_flag == 1) {
#ifdef USE_MPI
      // float: 1 -- Lifetime
      MPI_Pack_size(1, FloatDataType, EnzoTopComm, &size);
      this->ElementSizeInBytes += size;
#endif
    }
    else {
      this->ElementSizeInBytes += 1*sizeof(float);
    }
  };
  float *Lifetime;
};

class ActiveParticleType_PopIII : public ActiveParticleType
{
public:
  // Constructors
  ActiveParticleType_PopIII(void) : ActiveParticleType() {
    Lifetime = 0; 
  };
  ActiveParticleType_PopIII(PopIIIParticleBufferHandler *buffer, int index) :
    ActiveParticleType(static_cast<ParticleBufferHandler*>(buffer), index) {
    Lifetime = buffer->Lifetime[index];
  };
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &supp_data);
  static int WriteToOutput(ActiveParticleType **these_particles, int n, int GridRank, hid_t group_id);
  static int ReadFromOutput(ActiveParticleType **&particles_to_read, int &n, int GridRank, hid_t group_id);
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);
  static int EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data);
  template <class active_particle_class>
    static int BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				 int ThisLevel, int TotalStarParticleCountPrevious[],
				 int PopIIIParticleID) { return SUCCESS; };
  template <class active_particle_class>
    static int AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int TotalStarParticleCountPrevious[],
				int PopIIIParticleID) {return SUCCESS; };
  static int SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level, 
			      int TopGridDims[], int ActiveParticleID);
  static int InitializeParticleType();
  ENABLED_PARTICLE_ID_ACCESSOR

  // Pop III specific active particle parameters
  static float OverDensityThreshold, MetalCriticalFraction, 
    H2CriticalFraction, StarMass;

  float Lifetime;
};
