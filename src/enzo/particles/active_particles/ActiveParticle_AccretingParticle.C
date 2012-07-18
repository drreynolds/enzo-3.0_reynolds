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

#ifdef NEW_CONFIG

#include "ParameterControl/ParameterControl.h"
extern Configuration Param;

/* Set default parameter values. */

const char config_sink_particle_defaults[] =
"### SINK PARTICLE DEFAULTS ###\n"
"\n"
"Physics: {\n"
"    ActiveParticles: {\n"
"        AccretingParticle: {\n"
"            OverflowFactor       = 1.01;\n"
"            LinkingLength        = 4;\n   "
"            AccretionRadius      = 4;\n   "
"        };\n"
"    };\n"
"};\n";

#endif

#define NO_DEBUG

/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

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

class AccretingParticleGrid : private grid {
  friend class ActiveParticleType_AccretingParticle;
};

/* Note that we only refer to AccretingParticleGrid here. 
 * Given a grid object, we static cast to get this:
 *
 *    AccretingParticleGrid *thisgrid =
 *      static_cast<AccretingParticleGrid *>(thisgrid_orig); */

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

AccretingParticleBufferHandler::AccretingParticleBufferHandler
(ActiveParticleType **np, int NumberOfParticles, int type, int proc) : 
  ParticleBufferHandler(np, NumberOfParticles, type, proc) 
  {
    // Any extra fields must be added to the buffer and this->ElementSizeInBytes
    int i, index;
    if (this->NumberOfBuffers > 0) {
      this->AccretionRate = new float[this->NumberOfBuffers];
      for (i = 0, index=0; i < NumberOfParticles; i++)
	if (np[i]->ReturnType() == type && (np[i]->ReturnDestProcessor() == proc || proc==-1)) {
	  ActiveParticleType_AccretingParticle* temp = static_cast<ActiveParticleType_AccretingParticle*>(np[i]);
	  this->AccretionRate[index] = temp->AccretionRate;
	index++;
	}
    }
    this->CalculateAccretingParticleElementSize();
  };


float ActiveParticleType_AccretingParticle::OverflowFactor = FLOAT_UNDEFINED;
int ActiveParticleType_AccretingParticle::AccretionRadius = INT_UNDEFINED;
int ActiveParticleType_AccretingParticle::LinkingLength = INT_UNDEFINED;

int ActiveParticleType_AccretingParticle::InitializeParticleType()
{
#ifdef NEW_CONFIG

  // Update the parameter config to include the local defaults.  Note
  // that this does not overwrite any values previously specified.
  Param.Update(config_sink_particle_defaults);

  // Retrieve parameters from Param structure
  Param.GetScalar(OverflowFactor, "Physics.ActiveParticles.AccretingParticle.OverflowFactor");
  Param.GetScalar(LinkingLength, "Physics.ActiveParticles.AccretingParticle.LinkingLength");
  Param.GetScalar(AccretionRadius, "Physics.ActiveParticles.AccretingParticle.AccretionRadius");

#else

  // Leaving these defaults hardcoded for testing. NJG

  OverflowFactor = 1.01;
  LinkingLength = 4;
  AccretionRadius = 4;

#endif

  AccretingParticleBufferHandler *pbuffer = new AccretingParticleBufferHandler();
  delete pbuffer;

  // Need to turn on particle mass flagging if it isn't already turned on.

  bool TurnOnParticleMassRefinement = true;
  int method;
  for (method = 0; method < MAX_FLAGGING_METHODS; method++)
    if (CellFlaggingMethod[method] == 8 || CellFlaggingMethod[method] == 4) {
      TurnOnParticleMassRefinement = false;
      break;
    }
	
  if (TurnOnParticleMassRefinement) {
    method = 0;
    while(CellFlaggingMethod[method] != INT_UNDEFINED)
      method++;
    CellFlaggingMethod[method] = 4;
  }

  /* Add on the Particle Array Handlers */
  typedef ActiveParticleType_AccretingParticle ap;
  AttributeVector &ah = ap::AttributeHandlers;
  ActiveParticleType::SetupBaseParticleAttributes(ah);

  ah.push_back(Handler<ap, float, &ap::AccretionRate>("AccretionRate"));

  return SUCCESS;
}

int ActiveParticleType_AccretingParticle::EvaluateFormation
(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  // No need to do the rest if we're not on the maximum refinement level.
  if (data.level != MaximumRefinementLevel)
    return SUCCESS;

  AccretingParticleGrid *thisGrid =
    static_cast<AccretingParticleGrid *>(thisgrid_orig);
  

  float *tvel;

  int i,j,k,index,method,MassRefinementMethod;

  float *density = thisGrid->BaryonField[data.DensNum];
  float *velx = thisGrid->BaryonField[data.Vel1Num];
  float *vely = thisGrid->BaryonField[data.Vel2Num];
  float *velz = thisGrid->BaryonField[data.Vel3Num];
  float JeansDensityUnitConversion = (Gamma*pi*kboltz) / (Mu*mh*GravConst);
  float CellTemperature = 0;
  float JeansDensity = 0;
  float MassRefinementDensity = 0;
  float DensityThreshold = huge_number;
  float ExtraDensity = 0;

  int GridDimension[3] = {thisGrid->GridDimension[0],
                          thisGrid->GridDimension[1],
                          thisGrid->GridDimension[2]};

  FLOAT dx = thisGrid->CellWidth[0][0];
  
  bool HasMetalField = (data.MetalNum != -1 || data.ColourNum != -1);
  bool JeansRefinement = false;
  bool MassRefinement = false;

  // determine refinement criteria
  for (method = 0; method < MAX_FLAGGING_METHODS; method++) {
    if (CellFlaggingMethod[method] == 2) {
      MassRefinement = true;
      MassRefinementMethod = method;
    }
    if (CellFlaggingMethod[method] == 6) 
      JeansRefinement = true;
  }

  for (k = thisGrid->GridStartIndex[2]; k <= thisGrid->GridEndIndex[2]; k++) {
    for (j = thisGrid->GridStartIndex[1]; j <= thisGrid->GridEndIndex[1]; j++) {
      for (i = thisGrid->GridStartIndex[0]; i <= thisGrid->GridEndIndex[0]; i++) {
	index = GRIDINDEX_NOGHOST(i, j, k);

	// If no more room for particles, throw an ENZO_FAIL
	if (data.NumberOfNewParticles >= data.MaxNumberOfNewParticles)
	  return FAIL;

	// Does this cell violate the Jeans condition?
	DensityThreshold = huge_number;
	if (JeansRefinement) {
	  CellTemperature = (JeansRefinementColdTemperature > 0) ? JeansRefinementColdTemperature : data.Temperature[index];
	  JeansDensity = JeansDensityUnitConversion * OverflowFactor * CellTemperature / 
	    POW(data.LengthUnits*dx*RefineByJeansLengthSafetyFactor,2) / data.DensityUnits;
	  DensityThreshold = min(DensityThreshold,JeansDensity);
	}
	if (DensityThreshold == huge_number)
	  ENZO_FAIL("Error in Accreting Particles: Must refine by jeans length or overdensity!");
	
	if (density[index] <= DensityThreshold) 
	  continue;

	// Passed creation tests, create sink particle


	ActiveParticleType_AccretingParticle *np = new ActiveParticleType_AccretingParticle();
	data.NewParticles[data.NumberOfNewParticles++] = np;

	ExtraDensity = density[index] - DensityThreshold;
	np->Mass = ExtraDensity;   // Particle 'masses' are actually densities
	np->type = np->GetEnabledParticleID();
	np->BirthTime = thisGrid->ReturnTime();

	np->level = data.level;
	np->GridID = data.GridID;
	np->CurrentGrid = thisGrid;
	
	np->pos[0] = thisGrid->CellLeftEdge[0][i] + 0.5*thisGrid->CellWidth[0][i];
	np->pos[1] = thisGrid->CellLeftEdge[1][j] + 0.5*thisGrid->CellWidth[1][j];
	np->pos[2] = thisGrid->CellLeftEdge[2][k] + 0.5*thisGrid->CellWidth[2][k];
	
	tvel = thisGrid->AveragedVelocityAtCell(index,data.DensNum,data.Vel1Num);
	  
	np->vel[0] = tvel[0];
	np->vel[1] = tvel[1];
	np->vel[2] = tvel[2];
	
	if (HasMetalField)
	  np->Metallicity = data.TotalMetals[index];
	else
	  np->Metallicity = 0.0;

	np->AccretionRate = 0.0;
	
	// Remove mass from grid

	density[index] = DensityThreshold;

	// Clean up
	delete [] tvel;

      } // i
    } // j
  } // k
 
  return SUCCESS;
}  

int ActiveParticleType_AccretingParticle::EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  AccretingParticleGrid *thisGrid =
    static_cast<AccretingParticleGrid *>(thisgrid_orig);
  
  int npart = thisGrid->NumberOfActiveParticles;

  for (int n = 0; n < npart; n++) {
    ActiveParticleType_AccretingParticle *ThisParticle = 
      static_cast<ActiveParticleType_AccretingParticle*>(thisGrid->ActiveParticles[n]);
  
    ThisParticle->level = data.level;
  }

  return SUCCESS;
}

void ActiveParticleType_AccretingParticle::DescribeSupplementalData(ActiveParticleFormationDataFlags &flags)
{
  flags.DarkMatterDensity = true;
  flags.Temperature = true;
  flags.UnitConversions = true;
  flags.DataFieldNumbers = true;
  flags.MetalField = true;
}


int ActiveParticleType_AccretingParticle::WriteToOutput(ActiveParticleType **these_particles, int n, int GridRank, hid_t group_id)
{
  /* Create a new subgroup within the active particle group for active particles of type AccretingParticle */
  hid_t AccretingParticleGroupID = H5Gcreate(group_id,"AccretingParticle",0);

  writeScalarAttribute(AccretingParticleGroupID,HDF5_INT,"Number of Accreting Particles",&n);  

  char *ParticlePositionLabel[] =
     {"position_x", "position_y", "position_z"};
  char *ParticleVelocityLabel[] =
     {"velocity_x", "velocity_y", "velocity_z"};

  /* Create temporary buffers to store particle data */

  FLOAT *Position[MAX_DIMENSION];
  float *Velocity[MAX_DIMENSION]; 
  double *Mass = new double[n];
  float *BirthTime = new float[n];
  float *DynamicalTime = new float[n];
  float *Metallicity = new float[n];
  float *AccretionRate = new float[n];
  PINT *ID = new PINT[n];

  int i,dim;

  for (dim = 0; dim < GridRank; dim++) {
    Position[dim] = new FLOAT[n];
    Velocity[dim] = new float[n];
  }

  hsize_t TempInt;
  TempInt = n;
    
  ActiveParticleType_AccretingParticle *ParticleToWrite;
  for (i=0;i<n;i++) {
    ParticleToWrite = static_cast<ActiveParticleType_AccretingParticle*>(these_particles[i]);
    for (dim = 0; dim < GridRank; dim++) {
      Position[dim][i] = ParticleToWrite->pos[dim];
      Velocity[dim][i] = ParticleToWrite->vel[dim];
    }
    Mass[i] = ParticleToWrite->Mass;
    BirthTime[i] = ParticleToWrite->BirthTime;
    DynamicalTime[i] = ParticleToWrite->DynamicalTime;
    Metallicity[i] = ParticleToWrite->Metallicity;
    ID[i] = ParticleToWrite->Identifier;
    AccretionRate[i] = ParticleToWrite->AccretionRate;
  }

  for (dim = 0; dim < GridRank; dim++) {
    WriteDataset(1,&TempInt,ParticlePositionLabel[dim],
		 AccretingParticleGroupID, HDF5_FILE_PREC, (VOIDP) Position[dim]);
  }
  
  for (dim = 0; dim < GridRank; dim++) {
    WriteDataset(1,&TempInt,ParticleVelocityLabel[dim],
		  AccretingParticleGroupID, HDF5_REAL, (VOIDP) Velocity[dim]);
  }
  
  WriteDataset(1,&TempInt,"mass",AccretingParticleGroupID,HDF5_REAL,(VOIDP) Mass);
  WriteDataset(1,&TempInt,"creation_time",AccretingParticleGroupID,HDF5_REAL,(VOIDP) BirthTime);
  WriteDataset(1,&TempInt,"dynamical_time",AccretingParticleGroupID,HDF5_REAL,(VOIDP) DynamicalTime);
  WriteDataset(1,&TempInt,"metallicity_fraction",AccretingParticleGroupID,HDF5_REAL,(VOIDP) Metallicity);
  WriteDataset(1,&TempInt,"identifier",AccretingParticleGroupID,HDF5_PINT,(VOIDP) ID);
  WriteDataset(1,&TempInt,"accretion_rate",AccretingParticleGroupID,HDF5_REAL,(VOIDP) AccretionRate);

  /* Clean up */

  for (dim = 0; dim < GridRank; dim++) {
    delete[] Position[dim];
    delete[] Velocity[dim];
  }
  delete[] Mass;
  delete[] BirthTime;
  delete[] DynamicalTime;
  delete[] Metallicity;
  delete[] AccretionRate;
  H5Gclose(AccretingParticleGroupID);

  return SUCCESS;
}

int ActiveParticleType_AccretingParticle::ReadFromOutput(ActiveParticleType **&particles_to_read, int &n, int GridRank, hid_t group_id)
{

  int i,dim;
  hsize_t TempInt;

  hid_t AccretingParticleGroupID = H5Gopen(group_id,"AccretingParticle");

  readAttribute(AccretingParticleGroupID,HDF5_INT,"Number of Accreting Particles",&n);

  particles_to_read = new ActiveParticleType*[n];

  char *ParticlePositionLabel[] =
     {"position_x", "position_y", "position_z"};
  char *ParticleVelocityLabel[] =
     {"velocity_x", "velocity_y", "velocity_z"};

  FLOAT *Position[MAX_DIMENSION];
  float *Velocity[MAX_DIMENSION];
  double *Mass = new double[n];
  float *BirthTime = new float[n];
  float *DynamicalTime = new float[n];
  float *Metallicity = new float[n];
  float *AccretionRate = new float[n];
  PINT  *ID = new PINT[n];

  for (dim = 0; dim < GridRank; dim++) {
    Position[dim] = new FLOAT[n];
    Velocity[dim] = new float[n];
  }
  
  TempInt = n;
  
  for (dim = 0; dim < GridRank; dim++) {
    ReadDataset(1,&TempInt,ParticlePositionLabel[dim],
		  AccretingParticleGroupID, HDF5_FILE_PREC, (VOIDP) Position[dim]);
  }

  for (dim = 0; dim < GridRank; dim++) {
    ReadDataset(1,&TempInt,ParticleVelocityLabel[dim],
		  AccretingParticleGroupID, HDF5_REAL, (VOIDP) Velocity[dim]);
  }
  ReadDataset(1,&TempInt,"mass",AccretingParticleGroupID,HDF5_R8,(VOIDP) Mass);
  ReadDataset(1,&TempInt,"creation_time",AccretingParticleGroupID,HDF5_REAL,(VOIDP) BirthTime);
  ReadDataset(1,&TempInt,"dynamical_time",AccretingParticleGroupID,HDF5_REAL,(VOIDP) DynamicalTime);
  ReadDataset(1,&TempInt,"metallicity_fraction",AccretingParticleGroupID,HDF5_REAL,(VOIDP) Metallicity);
  ReadDataset(1,&TempInt,"identifier",AccretingParticleGroupID,HDF5_PINT,(VOIDP) ID);
  ReadDataset(1,&TempInt,"accretion_rate",AccretingParticleGroupID,HDF5_REAL,(VOIDP) AccretionRate);


  for (i = 0; i < n; i++) {
    ActiveParticleType_AccretingParticle *np = new ActiveParticleType_AccretingParticle();
    np->Mass = Mass[i];
    np->type = np->GetEnabledParticleID();
    np->BirthTime = BirthTime[i];
    np->DynamicalTime = DynamicalTime[i];
    np->Metallicity = Metallicity[i];
    np->Identifier = ID[i];
    for (dim = 0; dim < GridRank; dim++){
      np->pos[dim] = Position[dim][i];
      np->vel[dim] = Velocity[dim][i];
    }
    np->AccretionRate = AccretionRate[i];
    particles_to_read[i] = static_cast<ActiveParticleType*>(np);
  }

  delete[] Mass;
  delete[] BirthTime;
  delete[] DynamicalTime;
  delete[] Metallicity;
  delete[] ID;
  delete[] AccretionRate;

  for (dim = 0; dim < GridRank; dim++) {
    delete[] Position[dim];
    delete[] Velocity[dim];
  }
  H5Gclose(AccretingParticleGroupID);
  
  return SUCCESS;
}

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);

int ActiveParticleType_AccretingParticle::BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
							    int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
							    int ThisLevel, int TotalStarParticleCountPrevious[],
							    int AccretingParticleID)
{

  return SUCCESS;
}

ActiveParticleType_AccretingParticle** ActiveParticleType_AccretingParticle::MergeAccretingParticles
(int *nParticles, ActiveParticleType** ParticleList, FLOAT LinkingLength, 
 int *ngroups, LevelHierarchyEntry *LevelArray[], int ThisLevel)
{
  int i,j;
  int dim;
  int GroupNumberAssignment[*nParticles];
  FLOAT* tempPos = NULL;
  int *groupsize = NULL;
  int **grouplist = NULL;
  ActiveParticleType_AccretingParticle **MergedParticles = NULL;

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
  
  MergedParticles = new ActiveParticleType_AccretingParticle*[*ngroups]();
  
  /* Merge the mergeable groups */

  for (i=0; i<*ngroups; i++) {
    MergedParticles[i] = static_cast<ActiveParticleType_AccretingParticle*>(ParticleList[grouplist[i][0]]);
    if (groupsize[i] != 1) {
      for (j=1; j<groupsize[i]; j++) {
	MergedParticles[i]->Merge(static_cast<ActiveParticleType_AccretingParticle*>(ParticleList[grouplist[i][j]]));
	if (ParticleList[grouplist[i][j]]->DisableParticle(LevelArray,MergedParticles[i]->
							   ReturnCurrentGrid()->ReturnProcessorNumber()) == FAIL)
	  ENZO_FAIL("MergeAccretingParticles: DisableParticle failed!\n");
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
      MergedParticles[i]->DisableParticle(LevelArray,LevelGrids[NewGrid]->GridData->ReturnProcessorNumber()); 
      if (LevelGrids[j]->GridData->AddActiveParticle(static_cast<ActiveParticleType*>(MergedParticles[i])) == FAIL)
      	ENZO_FAIL("Active particle grid assignment failed!\n");
      MergedParticles[i]->AssignCurrentGrid(LevelGrids[NewGrid]->GridData);
    }
  }

  delete [] LevelGrids;

  *nParticles = *ngroups;

  return MergedParticles;
}

int AssignActiveParticlesToGrids(ActiveParticleType** ParticleList, int nParticles, 
				 LevelHierarchyEntry *LevelArray[]); 

int ActiveParticleType_AccretingParticle::AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
							   int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
							   int ThisLevel, int TotalStarParticleCountPrevious[],
							   int AccretingParticleID)
{

  /* Accreting particles live on the maximum refinement level.  If we are on a lower level, this does not concern us */

  if (ThisLevel == MaximumRefinementLevel)
    {

      /* Generate a list of all sink particles in the simulation box */
      int i,nParticles,NumberOfMergedParticles;
      ActiveParticleType** ParticleList = NULL;

      ParticleList = ActiveParticleFindAll(LevelArray, &nParticles, AccretingParticleID);

      /* Return if there are no accreting particles */
      
      if (nParticles == 0)
	return SUCCESS;

      /* Calculate CellWidth on maximum refinement level */

      // This assumes a cubic box and may not work for simulations with MinimumMassForRefinementLevelExponent
      FLOAT dx = (DomainRightEdge[0] - DomainLeftEdge[0]) /
	(MetaData->TopGridDims[0]*POW(FLOAT(RefineBy),FLOAT(MaximumRefinementLevel)));

      /* Do Merging   */

      ActiveParticleType_AccretingParticle **MergedParticles = NULL;
      
      /* Generate new merged list of sink particles */
      
      MergedParticles = MergeAccretingParticles(&nParticles, ParticleList, LinkingLength*dx,
						&NumberOfMergedParticles,LevelArray,ThisLevel);

      delete [] ParticleList;
   
      // Do merging twice to catch pathological cases where merging
      // leaves multiple sinks inside the same accretion zone.

      ParticleList = new ActiveParticleType*[NumberOfMergedParticles];

      for (i = 0; i<NumberOfMergedParticles; i++)
	ParticleList[i] = static_cast<ActiveParticleType*>(MergedParticles[i]);

      delete [] MergedParticles;

      MergedParticles = MergeAccretingParticles(&nParticles, ParticleList, LinkingLength*dx,
						&NumberOfMergedParticles,LevelArray,ThisLevel);

      delete [] ParticleList;

      if (debug && MyProcessorNumber == 1)
	printf("Number of particles after merging: %"ISYM"\n",NumberOfMergedParticles);

      /* Assign local particles to grids */
 
      ParticleList = new ActiveParticleType*[NumberOfMergedParticles];

      // need to use a bit of redirection because C++ pointer arrays have
      // trouble with polymorphism
      for (i = 0; i<NumberOfMergedParticles; i++)
	ParticleList[i] = static_cast<ActiveParticleType*>(MergedParticles[i]);

      if (AssignActiveParticlesToGrids(ParticleList,NumberOfMergedParticles, LevelArray) == FAIL)
	return FAIL;

      delete [] ParticleList;
      delete [] MergedParticles;
      
      /* Regenerate the global active particle list */
      
      ParticleList = ActiveParticleFindAll(LevelArray, &nParticles, AccretingParticleID);

      /* Do accretion */
      if (Accrete(nParticles,ParticleList,AccretionRadius,dx,LevelArray,ThisLevel) == FAIL)
	ENZO_FAIL("Accreting Particle accretion failed. \n");
     
      delete [] ParticleList;

    }

  return SUCCESS;
}

grid** ConstructFeedbackZones(ActiveParticleType** ParticleList, int nParticles, int FeedbackRadius, 
			     FLOAT dx, HierarchyEntry** Grids, int NumberOfGrids);

int DistributeFeedbackZones(grid** FeedbackZones, int NumberOfFeedbackZones,
			    HierarchyEntry** Grids, int NumberOfGrids);

int ActiveParticleType_AccretingParticle::Accrete(int nParticles, ActiveParticleType** ParticleList,
						  int AccretionRadius, FLOAT dx, 
						  LevelHierarchyEntry *LevelArray[], int ThisLevel)
{
  
  /* Skip accretion if we're not on the maximum refinement level. 
     This should only ever happen right after creation and then
     only in pathological cases where sink creation is happening at 
     the edges of two regions at the maximum refinement level */

  if (ThisLevel < MaximumRefinementLevel)
    return SUCCESS;

  /* For each particle, loop over all of the grids and do accretion 
     if the grid overlaps with the accretion zone                   */
  
  int i, NumberOfGrids;
  HierarchyEntry **Grids = NULL;
  grid *sinkGrid = NULL;
  
  bool SinkIsOnThisProc, SinkIsOnThisGrid;
  
  float SubtractedMass, SubtractedMomentum[3] = {};
  
  NumberOfGrids = GenerateGridArray(LevelArray, ThisLevel, &Grids);
  
  grid** FeedbackZones = ConstructFeedbackZones(ParticleList, nParticles, AccretionRadius, dx, Grids, NumberOfGrids);

  for (i = 0; i < nParticles; i++) {
    grid* FeedbackZone = FeedbackZones[i];
    if (MyProcessorNumber == FeedbackZone->ReturnProcessorNumber()) {
    
      float AccretionRate = 0;
        
      if (FeedbackZone->AccreteOntoAccretingParticle(&ParticleList[i],AccretionRadius*dx,&AccretionRate) == FAIL)
	return FAIL;
  
      // No need to communicate the accretion rate to the other CPUs since this particle is already local.
      static_cast<ActiveParticleType_AccretingParticle*>(ParticleList[i])->AccretionRate = AccretionRate;
    }
  }
  
  DistributeFeedbackZones(FeedbackZones, nParticles, Grids, NumberOfGrids);

  for (i = 0; i < nParticles; i++) {
    delete FeedbackZones[i];    
  }

  delete [] FeedbackZones;

  if (AssignActiveParticlesToGrids(ParticleList, nParticles, LevelArray) == FAIL)
    return FAIL;

  delete [] Grids;
  return SUCCESS;
}

int ActiveParticleType_AccretingParticle::SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level, 
							   int TopGridDims[], int AccretingParticleID)
{
  /* Generate a list of all sink particles in the simulation box */
  int i, nParticles;
  FLOAT *pos = NULL, dx=0;
  ActiveParticleType **AccretingParticleList = NULL ;
  LevelHierarchyEntry *Temp = NULL;
  
  AccretingParticleList = ActiveParticleFindAll(LevelArray, &nParticles, AccretingParticleID);
  
  /* Calculate CellWidth on maximum refinement level */
  
  // this will fail for noncubic boxes or simulations with MinimimMassForRefinementLevelExponent
  dx = (DomainRightEdge[0] - DomainLeftEdge[0]) /
    (TopGridDims[0]*POW(FLOAT(RefineBy),FLOAT(MaximumRefinementLevel)));
  
  for (i=0 ; i<nParticles; i++){
    pos = AccretingParticleList[i]->ReturnPosition();
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      if (Temp->GridData->DepositAccretionZone(level,pos,AccretionRadius*dx) == FAIL) {
	ENZO_FAIL("Error in grid->DepositAccretionZone.\n")
	  }
  }

  return SUCCESS;
}

void AccretingParticleBufferHandler::AllocateBuffer(ActiveParticleType **np, int NumberOfParticles, 
						    char *buffer, Eint32 total_buffer_size, int &buffer_size,
						    Eint32 &position, int type_num, int proc=-1)
{
  AccretingParticleBufferHandler *pbuffer = new AccretingParticleBufferHandler(np, NumberOfParticles, type_num, proc);
  pbuffer->_AllocateBuffer(buffer, total_buffer_size, buffer_size, position);
  // If any extra fields are added in the future, then they would be
  // transferred to the buffer here.
#ifdef USE_MPI
  if (pbuffer->NumberOfBuffers > 0) {
    MPI_Pack(pbuffer->AccretionRate,pbuffer->NumberOfBuffers, FloatDataType, buffer, total_buffer_size,
	     &position, EnzoTopComm);
  }
#endif /* USE_MPI */
  delete pbuffer;
  return;
}

void AccretingParticleBufferHandler::UnpackBuffer(char *mpi_buffer, int mpi_buffer_size, int NumberOfParticles,
						  ActiveParticleType **np, int &npart)
{
  int i;
  Eint32 position = 0;
  AccretingParticleBufferHandler *pbuffer = new AccretingParticleBufferHandler(NumberOfParticles);
  pbuffer->_UnpackBuffer(mpi_buffer, mpi_buffer_size, position);
  // If any extra fields are added in the future, then they would be
  // transferred to the buffer here.
#ifdef USE_MPI
  if (pbuffer->NumberOfBuffers > 0) {
    MPI_Unpack(mpi_buffer, mpi_buffer_size, &position, pbuffer->AccretionRate,
	       pbuffer->NumberOfBuffers, FloatDataType, EnzoTopComm);
  }
#endif
  /* Convert the particle buffer into active particles */
  
  for (i = 0; i < pbuffer->NumberOfBuffers; i++)
    np[npart++] = new ActiveParticleType_AccretingParticle(pbuffer, i);
  delete pbuffer;
  return;
}

namespace {
  ActiveParticleType_info *AccretingParticleInfo = 
    register_ptype <ActiveParticleType_AccretingParticle, AccretingParticleBufferHandler> 
    ("AccretingParticle");
}
std::vector<ParticleAttributeHandler>
  ActiveParticleType_AccretingParticle::AttributeHandlers;

#undef DEBUG
