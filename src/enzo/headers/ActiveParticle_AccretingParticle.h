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

  template <typename APT>
  static APT** MergeAccretingParticles(int *nParticles, ActiveParticleType** ParticleList, FLOAT LinkingLength, 
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

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);

template <typename APT>
APT** ActiveParticleType_AccretingParticle::MergeAccretingParticles
(int *nParticles, ActiveParticleType** ParticleList, FLOAT LinkingLength, 
 int *ngroups, LevelHierarchyEntry *LevelArray[], int ThisLevel)
{
  int i,j;
  int dim;
  int GroupNumberAssignment[*nParticles];
  FLOAT* tempPos = NULL;
  int *groupsize = NULL;
  int **grouplist = NULL;
  APT **MergedParticles = NULL;

  HierarchyEntry** LevelGrids = NULL;

  int NumberOfGrids = GenerateGridArray(LevelArray, ThisLevel, &LevelGrids);

  /* Construct list of sink particle positions to pass to Foflist */
  FLOAT ParticleCoordinates[3*(*nParticles)];
  
  for (i=0; i<(*nParticles); i++) {
    tempPos = ParticleList[i]->ReturnPosition();
    for (dim=0; dim<3; dim++)
      ParticleCoordinates[3*i+dim] = tempPos[dim];
  }

  /* Find mergeable groups using an FOF search */

  *ngroups = FofList((*nParticles), ParticleCoordinates, LinkingLength, GroupNumberAssignment, &groupsize, &grouplist);
  
  MergedParticles = new APT*[*ngroups]();
  
  /* Merge the mergeable groups */

  for (i=0; i<*ngroups; i++) {
    MergedParticles[i] = static_cast<APT*>(ParticleList[grouplist[i][0]]);
    if (groupsize[i] != 1) {
      for (j=1; j<groupsize[i]; j++) {
	MergedParticles[i]->Merge(static_cast<APT*>(ParticleList[grouplist[i][j]]));
	if (ParticleList[grouplist[i][j]]->DisableParticle(LevelArray,MergedParticles[i]->
							   ReturnCurrentGrid()->ReturnProcessorNumber()) == FAIL)
	  ENZO_FAIL("MergeAccretingParticles: DisableParticle failed!\n");
	if (NumberOfProcessors > 1) {
	  delete ParticleList[grouplist[i][j]];
	  ParticleList[grouplist[i][j]] = NULL;
	}
      }
    }
  }

  delete [] groupsize;
  groupsize = NULL;
  for (i=0; i<*ngroups; i++)
    delete [] grouplist[i];
  delete [] grouplist;
  grouplist = NULL;

  /* Loop over the grids and check if any of the merged particles have
     moved. If so, disable the particle on the current grid and assign
     it to the new grid*/

  int NewGrid = -1;

  for (i = 0; i < *ngroups; i++) {
    if (MergedParticles[i]->ReturnCurrentGrid()->PointInGrid(MergedParticles[i]->ReturnPosition()) == false) {
      // Find the grid to transfer to 
      for (j = 0; j < NumberOfGrids; j++) {
	if (LevelGrids[j]->GridData->PointInGrid(MergedParticles[i]->ReturnPosition())) {
	  NewGrid = j;
	  break;
	}
      }
      if (NewGrid == -1)
	ENZO_FAIL("Cannot assign particle to grid after merging!\n");
      int OldProc = MergedParticles[i]->CurrentGrid->ReturnProcessorNumber();
      APT *temp = new APT(MergedParticles[i]);
      MergedParticles[i]->DisableParticle(LevelArray,LevelGrids[NewGrid]->GridData->ReturnProcessorNumber()); 
      if (LevelGrids[NewGrid]->GridData->AddActiveParticle(static_cast<ActiveParticleType*>(temp)) == FAIL)
      	ENZO_FAIL("Active particle grid assignment failed!\n");
      if (MyProcessorNumber == OldProc) {
	delete MergedParticles[i];
	MergedParticles[i] = new APT(temp);
      }
      else if (MyProcessorNumber != temp->CurrentGrid->ReturnProcessorNumber())
	delete temp;
      MergedParticles[i]->AssignCurrentGrid(LevelGrids[NewGrid]->GridData);
    }
  }

  delete [] LevelGrids;

  *nParticles = *ngroups;

  return MergedParticles;
}


