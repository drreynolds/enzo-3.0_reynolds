/***********************************************************************
/
/  Cen & Ostriker star formation
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

class ActiveParticleType_CenOstriker;
class CenOstrikerBufferHandler : public ParticleBufferHandler
{
public:
  // No extra fields in CenOstriker.  Same base constructor.
  CenOstrikerBufferHandler(void) : ParticleBufferHandler() {};
  CenOstrikerBufferHandler(int NumberOfParticles) : ParticleBufferHandler(NumberOfParticles) {
#ifdef EXAMPLE
    this->field = new float[NumberOfParticles];
#endif
  };
  CenOstrikerBufferHandler(ActiveParticleType **np, int NumberOfParticles, int type, int proc) : 
    ParticleBufferHandler(np, NumberOfParticles, type, proc) {
    // Any extra fields must be added to the buffer and this->ElementSizeInBytes
#ifdef EXAMPLE
    this->field = new float[this->NumberOfBuffers];
    index = 0;
    for (i = 0; i < NumberOfParticles; i++)
      if (np[i]->ReturnType() == type && (np[i]->ReturnDestProcessor() == proc || proc==-1)) {
	this->field[index] = np[i]->field;
	index++;
      }
    this->ElementSizeInBytes += 1*sizeof(float);
#endif /* EXAMPLE */
  };
  ~CenOstrikerBufferHandler() {
#ifdef EXAMPLE
    delete[] this->field;
#endif
  };
  static int ReturnHeaderSize(void) {return CenOstrikerBufferHandler::HeaderSizeInBytes; }
  static int ReturnElementSize(void) {return CenOstrikerBufferHandler::ElementSizeInBytes; }
  static void AllocateBuffer(ActiveParticleType **np, int NumberOfParticles, char *buffer, 
			     Eint32 total_buffer_size, int &buffer_size, Eint32 &position, 
			     int type_num, int proc=-1);
  static void UnpackBuffer(char *mpi_buffer, int mpi_buffer_size, int NumberOfParticles,
			   ActiveParticleType **np, int &npart);
};



class ActiveParticleType_CenOstriker : public ActiveParticleType
{
public:
  // Constructors
  ActiveParticleType_CenOstriker(void) : ActiveParticleType() {};
  ActiveParticleType_CenOstriker(CenOstrikerBufferHandler *buffer, int index) :
    ActiveParticleType(static_cast<ParticleBufferHandler*>(buffer), index) {
    // Add any additional fields here
#ifdef EXAMPLE
    field = buffer->field[index];
#endif
  };

  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static int EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);
  static int WriteToOutput(ActiveParticleType **these_particles, int n, int GridRank, hid_t group_id);
  static int ReadFromOutput(ActiveParticleType **&particles_to_read, int &n, int GridRank, hid_t group_id);
  static void UnpackBuffer(ActiveParticleType *np, ParticleBufferHandler **buffer, int place);
  static ParticleBufferHandler *AllocateBuffers(int NumberOfParticles);
  static int BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
			       int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			       int ThisLevel, int TotalActiveParticleCountPrevious[],
			       int CenOstrikerID);
  static int AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
			      int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			      int ThisLevel, int TotalActiveParticleCountPrevious[],
			      int CenOstrikerID);
  static int SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level, int TopGridDims[], int ActiveParticleID);
  static int InitializeParticleType();
  ENABLED_PARTICLE_ID_ACCESSOR
  
  static float OverdensityThreshold, MassEfficiency, MinimumDynamicalTime, 
    MinimumStarMass, MassEjectionFraction, EnergyToThermalFeedback, MetalYield;

  static int FeedbackDistTotalCells, FeedbackDistRadius, FeedbackDistCellStep;

  static bool JeansMassCriterion, StochasticStarFormation, UnigridVelocities, 
    PhysicalOverdensity, dtDependence;

  friend class ParticleBufferHandler;
};

