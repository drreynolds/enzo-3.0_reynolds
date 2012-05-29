/***********************************************************************
/
/ Accreting Particle
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif 

#include <string.h>
#include <map>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <iostream>
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
#include "ActiveParticle.h"
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
    this->AccretionRate = new float[NumberOfParticles];
    this->cInfinity = new float[NumberOfParticles];
    this->vInfinity = new float[NumberOfParticles];
    this->BondiHoyleRadius = new FLOAT[NumberOfParticles];
    this->CalculateAccretingParticleElementSize();
  };
  AccretingParticleBufferHandler(ActiveParticleType **np, int NumberOfParticles, int type, int proc);
  ~AccretingParticleBufferHandler() {
    delete [] AccretionRate;
    delete [] cInfinity;
    delete [] vInfinity;
    delete [] BondiHoyleRadius;
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
      // float: 3 -- AccretionRate, cInfinity, vInfinity
      // FLOAT: 1 -- BondiHoyleRadius
      MPI_Pack_size(3, FloatDataType, MPI_COMM_WORLD, &size);
      this->ElementSizeInBytes += size;
      MPI_Pack_size(1, MY_MPIFLOAT, MPI_COMM_WORLD, &size);
      this->ElementSizeInBytes += size;
#endif
    } else {
      this->ElementSizeInBytes += 3*sizeof(float) + sizeof(FLOAT);
    }
  };
  float* AccretionRate;
  float* cInfinity;
  float* vInfinity;
  FLOAT* BondiHoyleRadius;
};

class AccretingParticleGrid : private grid {
  friend class ActiveParticleType_AccretingParticle;
};

/* Note that we only refer to AccretingParticleGrid here. 
 * Given a grid object, we static case to get this:
 *
 *    AccretingParticleGrid *thisgrid =
 *      static_cast<AccretingParticleGrid *>(thisgrid_orig); */

class ActiveParticleType_AccretingParticle : public ActiveParticleType
{
public:
  // Constructors
  ActiveParticleType_AccretingParticle(void) : ActiveParticleType() {
    
    AccretionRate = 0;
    vInfinity = 0;
    cInfinity = 0;
    BondiHoyleRadius = 0;
  };
  ActiveParticleType_AccretingParticle(AccretingParticleBufferHandler *buffer, int index) :
    ActiveParticleType(static_cast<ParticleBufferHandler*>(buffer), index) {
    AccretionRate = buffer->AccretionRate[index];
    vInfinity = buffer->vInfinity[index];
    cInfinity = buffer->cInfinity[index];
    BondiHoyleRadius = buffer->BondiHoyleRadius[index];
  };
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static int WriteToOutput(ActiveParticleType **these_particles, int n, int GridRank, hid_t group_id);
  static int ReadFromOutput(ActiveParticleType **&particles_to_read, int &n, int GridRank, hid_t group_id);
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);
  static ParticleBufferHandler *AllocateBuffers(int NumberOfParticles);
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
  int AdjustBondiHoyle(grid* CurrentGrid);

  ENABLED_PARTICLE_ID_ACCESSOR

  // sink helper routines

  static ActiveParticleType_AccretingParticle**  MergeAccretingParticles
  (int *nParticles, ActiveParticleType** ParticleList, FLOAT LinkingLength, 
   int *ngroups, LevelHierarchyEntry *LevelArray[]);

  static int Accrete(int nParticles, ActiveParticleType** ParticleList,
		     FLOAT AccretionRadius, LevelHierarchyEntry *LevelArray[], 
		     int ThisLevel);

  static float OverflowFactor;
  static int AccretionRadius;   // in units of CellWidth on the maximum refinement level
  static int LinkingLength;     // Should be equal to AccretionRadius
  float AccretionRate;
  float vInfinity;
  float cInfinity;
  FLOAT BondiHoyleRadius;
};

AccretingParticleBufferHandler::AccretingParticleBufferHandler
(ActiveParticleType **np, int NumberOfParticles, int type, int proc) : 
  ParticleBufferHandler(np, NumberOfParticles, type, proc) 
  {
    // Any extra fields must be added to the buffer and this->ElementSizeInBytes
    int i, index;
    this->AccretionRate = new float[this->NumberOfBuffers];
    this->cInfinity = new float[this->NumberOfBuffers];
    this->vInfinity = new float[this->NumberOfBuffers];
    this->BondiHoyleRadius = new FLOAT[this->NumberOfBuffers];
    for (i = 0, index=0; i < NumberOfParticles; i++)
      if (np[i]->ReturnType() == type && (np[i]->ReturnDestProcessor() == proc || proc==-1)) {
	ActiveParticleType_AccretingParticle* temp = static_cast<ActiveParticleType_AccretingParticle*>(np[i]);
	this->AccretionRate[index] = temp->AccretionRate;
	this->cInfinity[index] = temp->cInfinity;
	this->vInfinity[index] = temp->vInfinity;
	this->BondiHoyleRadius[index] = temp->BondiHoyleRadius;
	index++;
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

  return SUCCESS;
}

int ActiveParticleType_AccretingParticle::EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  // No need to do the rest if we're not on the maximum refinement level.
  if (data.level != MaximumRefinementLevel)
    return SUCCESS;

  AccretingParticleGrid *thisGrid =
    static_cast<AccretingParticleGrid *>(thisgrid_orig);
  

  int GridDimension[MAX_DIMENSION] = {thisGrid->GridDimension[0],
				      thisGrid->GridDimension[1],
				      thisGrid->GridDimension[2]};

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
      index = GRIDINDEX_NOGHOST(thisGrid->GridStartIndex[0], j, k);
      for (i = thisGrid->GridStartIndex[0]; i <= thisGrid->GridEndIndex[0]; i++, index++) {

	// If no more room for particles, throw an ENZO_FAIL
	if (data.NumberOfNewParticles >= data.MaxNumberOfNewParticles)
	  return FAIL;

	// Does this cell violate the Jeans condition?
	if (JeansRefinement) {
	  CellTemperature = (JeansRefinementColdTemperature > 0) ? JeansRefinementColdTemperature : data.Temperature[index];
	  JeansDensity = JeansDensityUnitConversion * OverflowFactor * CellTemperature / 
	    POW(data.LengthUnits*dx*RefineByJeansLengthSafetyFactor,2) / data.DensityUnits;
	  DensityThreshold = min(DensityThreshold,JeansDensity);
	}
	if (MassRefinement) {
	  MassRefinementDensity = MinimumMassForRefinement[MassRefinementMethod]*
	    pow(RefineBy, data.level*MinimumMassForRefinementLevelExponent[MassRefinementMethod])/POW(dx,3);
	  DensityThreshold = min(DensityThreshold,MassRefinementDensity);
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

	np->vInfinity = 0.0;  // Since np->vel = tvel
	np->cInfinity = sqrt(Gamma*kboltz*CellTemperature/(Mu*mh))/data.LengthUnits*data.TimeUnits;
	// Particle "Mass" is actually a density
	np->BondiHoyleRadius = GravitationalConstant*(np->ReturnMass()*POW(dx,3))/
	  (pow(np->vInfinity,2) + pow(np->cInfinity,2));
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
  float *density = thisGrid->BaryonField[data.DensNum];
  float *velx = thisGrid->BaryonField[data.Vel1Num];
  float *vely = thisGrid->BaryonField[data.Vel2Num];
  float *velz = thisGrid->BaryonField[data.Vel3Num];
  float *totalenergy = thisGrid->BaryonField[data.TENum];
  float *gasenergy = thisGrid->BaryonField[data.GENum];

  int GridDimension[3] = {thisGrid->GridDimension[0],
                          thisGrid->GridDimension[1],
                          thisGrid->GridDimension[2]};

  float *vel,CellTemperature;
  int n,i,j,k,index,dim;

  FLOAT *pos, LeftCorner[MAX_DIMENSION];

  FLOAT dx = thisGrid->CellWidth[0][0];

  for (dim = 0; dim < 3; dim++) 
    LeftCorner[dim] = thisGrid->CellLeftEdge[dim][0];

  /* Loop over all of the AccretingParticles on this grid and set the bondi-hoyle radius */

  for (n = 0; n < npart; n++) {
    ActiveParticleType_AccretingParticle *ThisParticle = 
      static_cast<ActiveParticleType_AccretingParticle*>(thisGrid->ActiveParticles[n]);
    pos = ThisParticle->ReturnPosition();
    vel = ThisParticle->ReturnVelocity(); 
    i = int((pos[0] - LeftCorner[0])/thisGrid->CellWidth[0][0]);
    j = int((pos[1] - LeftCorner[1])/thisGrid->CellWidth[0][0]);
    k = int((pos[2] - LeftCorner[2])/thisGrid->CellWidth[0][0]);
    index = GRIDINDEX_NOGHOST(i,j,k);
    ThisParticle->vInfinity = sqrt(pow((vel[0] - velx[index]),2) +
				   pow((vel[1] - vely[index]),2) +
				   pow((vel[2] - velz[index]),2));
    CellTemperature = (JeansRefinementColdTemperature > 0) ? JeansRefinementColdTemperature : data.Temperature[index];
    ThisParticle->cInfinity = sqrt(Gamma*kboltz*CellTemperature/(Mu*mh))/data.LengthUnits*data.TimeUnits;
    // Convert sound speed to enzo internal units.
    ThisParticle->BondiHoyleRadius = GravitationalConstant*(ThisParticle->ReturnMass()*POW(dx,3))/
      (pow(ThisParticle->vInfinity,2) + pow(ThisParticle->cInfinity,2));
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
  float *vInfinity = new float[n];
  float *cInfinity = new float[n];
  FLOAT *BondiHoyleRadius = new FLOAT[n];
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
    vInfinity[i] = ParticleToWrite->vInfinity;
    cInfinity[i] = ParticleToWrite->cInfinity;
    BondiHoyleRadius[i] = ParticleToWrite->BondiHoyleRadius;
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
  WriteDataset(1,&TempInt,"c_infinity",AccretingParticleGroupID,HDF5_REAL,(VOIDP) cInfinity);
  WriteDataset(1,&TempInt,"v_infinity",AccretingParticleGroupID,HDF5_REAL,(VOIDP) vInfinity);
  WriteDataset(1,&TempInt,"bondi_hoyle_radius",AccretingParticleGroupID,HDF5_REAL,(VOIDP) BondiHoyleRadius);
 

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
  delete[] cInfinity;
  delete[] vInfinity;
  delete[] BondiHoyleRadius;
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
  float *vInfinity = new float[n];
  float *cInfinity = new float[n];
  FLOAT *BondiHoyleRadius = new FLOAT[n];
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
  ReadDataset(1,&TempInt,"c_infinity",AccretingParticleGroupID,HDF5_REAL,(VOIDP) cInfinity);
  ReadDataset(1,&TempInt,"v_infinity",AccretingParticleGroupID,HDF5_REAL,(VOIDP) vInfinity);
  ReadDataset(1,&TempInt,"bondi_hoyle_radius",AccretingParticleGroupID,HDF5_REAL,(VOIDP) BondiHoyleRadius);


  for (i = 0; i < n; i++) {
    ActiveParticleType_AccretingParticle *np = new ActiveParticleType_AccretingParticle();
    np->Mass = Mass[i];
    np->type = AccretingParticle;
    np->BirthTime = BirthTime[i];
    np->DynamicalTime = DynamicalTime[i];
    np->Metallicity = Metallicity[i];
    np->Identifier = ID[i];
    for (dim = 0; dim < GridRank; dim++){
      np->pos[dim] = Position[dim][i];
      np->vel[dim] = Velocity[dim][i];
    }
    np->AccretionRate = AccretionRate[i];
    np->cInfinity = cInfinity[i];
    np->vInfinity = vInfinity[i];
    np->BondiHoyleRadius = BondiHoyleRadius[i];
    particles_to_read[i] = static_cast<ActiveParticleType*>(np);
  }

  delete[] Mass;
  delete[] BirthTime;
  delete[] DynamicalTime;
  delete[] Metallicity;
  delete[] ID;
  delete[] AccretionRate;
  delete[] cInfinity;
  delete[] vInfinity;
  delete[] BondiHoyleRadius;

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
(int *nParticles, ActiveParticleType** ParticleList, FLOAT LinkingLength, int *ngroups, LevelHierarchyEntry *LevelArray[])
{
  int i,j;
  int dim;
  int GroupNumberAssignment[*nParticles];
  FLOAT* tempPos = NULL;
  int *groupsize = NULL;
  int **grouplist = NULL;
  ActiveParticleType_AccretingParticle **MergedParticles = NULL;
  bool debug = false;

  /* Construct list of sink particle positions to pass to Foflist */
  FLOAT ParticleCoordinates[3*(*nParticles)];
  
  for (i=0; i<(*nParticles); i++) {
    tempPos = ParticleList[i]->ReturnPosition();
    for (dim=0; dim<3; dim++)
      ParticleCoordinates[3*i+dim] = tempPos[dim];
    if (MyProcessorNumber == 0 && debug)
      fprintf(stderr,"%"ISYM", %"GSYM", %"GSYM", %"GSYM"\n", 
	      i, tempPos[0], tempPos[1], tempPos[2]);
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
	ParticleList[grouplist[i][j]]->DisableParticle(LevelArray);
      }
    }
  }

  *nParticles = *ngroups;

  return MergedParticles;
}

int CommunicationSyncNumberOfParticles(HierarchyEntry *GridHierarchyPointer[],int NumberOfGrids);

int ActiveParticleType_AccretingParticle::AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
						      int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
						      int ThisLevel, int TotalStarParticleCountPrevious[],
						      int AccretingParticleID)
{

  /* Accreting particles live on the maximum refinement level.  If we are on a lower level, this does not concern us */

  if (ThisLevel == MaximumRefinementLevel)
    {

      /* Generate a list of all sink particles in the simulation box */
      int i,grid,nParticles,NumberOfMergedParticles;
      HierarchyEntry **LevelGrids = NULL;
      ActiveParticleType** ParticleList = NULL;
      bool debug = true;

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
						&NumberOfMergedParticles,LevelArray);
      if (MyProcessorNumber == 0 && debug)
	fprintf(stderr,"Number of Accreting Particles After Merging = %"ISYM"\n",NumberOfMergedParticles);

      delete [] ParticleList;
   
      /* Assign local particles to grids */
 
      int level, LevelMax, SavedGrid, NumberOfGrids;
      FLOAT* pos = NULL;
      float mass;

      for (i = 0; i<NumberOfMergedParticles; i++) {
	LevelMax = SavedGrid = -1;
	NumberOfGrids = 0;
	if (MyProcessorNumber == 0 && debug) {
	  pos = MergedParticles[i]->ReturnPosition();
	  mass = MergedParticles[i]->ReturnMass();
	  fprintf(stderr,"%"ISYM", %"GSYM", %"GSYM", %"GSYM", %"FSYM"\n",
		  i,pos[0],pos[1],pos[2],mass);
	}
	for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
	  NumberOfGrids = GenerateGridArray(LevelArray, level, &LevelGrids);     
	  for (grid = 0; grid < NumberOfGrids; grid++) 
	    if (LevelGrids[grid]->GridData->ReturnProcessorNumber() == MyProcessorNumber)
	      if (LevelGrids[grid]->GridData->PointInGrid(MergedParticles[i]->ReturnPosition()) == true &&
		  LevelGrids[grid]->GridData->isLocal() == true) { 
		SavedGrid = grid;
		LevelMax = level;
	      }
	  delete [] LevelGrids;
	  LevelGrids = NULL;
	}
	
	/* Find the processor which has the maximum value of LevelMax
	   The repeated code in the serial and parallel implimentations 
	   kind of sucks - should this be different? */
	
#ifdef USE_MPI
	struct { Eint32 value; Eint32 rank; } sendbuf, recvbuf;
	MPI_Comm_rank(MPI_COMM_WORLD, &sendbuf.rank); 
	sendbuf.value = LevelMax;
	MPI_Allreduce(&sendbuf, &recvbuf, 1, MPI_2INT, MPI_MAXLOC, MPI_COMM_WORLD);
	NumberOfGrids = GenerateGridArray(LevelArray, recvbuf.value, &LevelGrids); 
	if (LevelMax == recvbuf.value) {
	  MergedParticles[i]->AdjustBondiHoyle(LevelGrids[SavedGrid]->GridData);
	  if (LevelGrids[SavedGrid]->GridData->AddActiveParticle(static_cast<ActiveParticleType*>(MergedParticles[i])) == FAIL)
	    ENZO_FAIL("Active particle grid assignment failed"); 
	}
	LevelMax = recvbuf.value;
#else // endif parallel
	NumberOfGrids = GenerateGridArray(LevelArray, LevelMax, &LevelGrids); 
	MergedParticles[i]->AdjustBondiHoyle(LevelGrids[SavedGrid]->GridData);
	if (LevelGrids[SavedGrid]->GridData->AddActiveParticle(static_cast<ActiveParticleType*>(MergedParticles[i])) == FAIL)
	  ENZO_FAIL("Active particle grid assignment failed");
#endif //endif serial
	
	/* Sync the updated particle counts accross all proccessors */

	CommunicationSyncNumberOfParticles(LevelGrids, NumberOfGrids);

	delete [] LevelGrids;

      }

      delete [] MergedParticles;

      /* Regenerate the global active particle list */
      
      ParticleList = ActiveParticleFindAll(LevelArray, &nParticles, AccretingParticleID);

      /* Do accretion */
      
      if (Accrete(nParticles,ParticleList,AccretionRadius*dx,LevelArray,ThisLevel) == FAIL)
	ENZO_FAIL("Accreting Particle accretion failed. \n");
	  
      delete [] ParticleList;
      
    }

  return SUCCESS;
}

int ActiveParticleType_AccretingParticle::Accrete(int nParticles, ActiveParticleType** ParticleList,
						  FLOAT AccretionRadius, LevelHierarchyEntry *LevelArray[],
						  int ThisLevel)
{
  
  /* Skip accretion if we're not on the maximum refinement level. 
     This should only ever happen right after creation and then
     only in pathological cases where sink creation is happening at 
     the edges of two regions at the maximum refinement level */

  if (ThisLevel < MaximumRefinementLevel)
    return SUCCESS;

  /* For each particle, loop over all of the grids and do accretion 
     if the grid overlaps with the accretion zone                   */
  
  int i, grid, NumberOfGrids;
  HierarchyEntry **Grids = NULL;
  HierarchyEntry *sinkGrid = NULL;
  
  bool SinkIsOnThisProc, SinkIsOnThisGrid;
  
  float WeightedSum, SumOfWeights, GlobalWeightedSum, GlobalSumOfWeights, AverageDensity, SubtractedMass, 
    GlobalSubtractedMass, SubtractedMomentum[3], GlobalSubtractedMomentum[3], vInfinity, cInfinity, BondiHoyleRadius, 
    AccretionRate;
  int NumberOfCells = 0;
  
  for (i = 0; i < 3; i++) {
    SubtractedMomentum[i] = 0;
    GlobalSubtractedMomentum[i] = 0;
  }
  
  NumberOfGrids = GenerateGridArray(LevelArray, ThisLevel, &Grids);
  
  for (i = 0; i < nParticles; i++) {
    WeightedSum = SumOfWeights = GlobalWeightedSum = GlobalSumOfWeights = AverageDensity = SubtractedMass = 
      GlobalSubtractedMass = 0, SinkIsOnThisProc = false;
    ActiveParticleType_AccretingParticle* temp = static_cast<ActiveParticleType_AccretingParticle*>(ParticleList[i]);
    vInfinity = temp->vInfinity;
    cInfinity = temp->cInfinity;
    BondiHoyleRadius = temp->BondiHoyleRadius;
    for (grid = 0; grid < NumberOfGrids; grid++) {
      if (Grids[grid]->GridData->
	  FindAverageDensityInAccretionZone(ParticleList[i],AccretionRadius, &WeightedSum, 
					    &SumOfWeights, &NumberOfCells, BondiHoyleRadius) == FAIL)
	return FAIL;
    }
    /* sum up rhobar on root and broadcast result to all processors */
#ifdef USE_MPI
    MPI_Allreduce(&WeightedSum,  &GlobalWeightedSum,  1, FloatDataType, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&SumOfWeights, &GlobalSumOfWeights, 1, FloatDataType, MPI_SUM, MPI_COMM_WORLD);
#else
    GlobalWeightedSum = WeightedSum;
    GlobalSumOfWeights = SumOfWeights;
#endif
    
    AverageDensity = GlobalWeightedSum / GlobalSumOfWeights;
    
    /* Now perform accretion algorithm by modifying the grids locally */
    for (grid = 0; grid < NumberOfGrids; grid++) {
      SinkIsOnThisGrid = false;
      if (Grids[grid]->GridData->
	  AccreteOntoAccretingParticle(ParticleList[i], AccretionRadius, AverageDensity, GlobalSumOfWeights, &SubtractedMass, 
				       SubtractedMomentum, &SinkIsOnThisGrid, vInfinity, cInfinity, 
				       BondiHoyleRadius, &AccretionRate) == FAIL) {
	return FAIL;
      }
      if (SinkIsOnThisGrid) {
	sinkGrid = Grids[grid];
	SinkIsOnThisProc = true;
	break;
      }
    }

#ifdef USE_MPI
    MPI_Allreduce(&SubtractedMass, &GlobalSubtractedMass, 1, FloatDataType, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&SubtractedMomentum, &GlobalSubtractedMomentum, 3, FloatDataType, MPI_SUM, MPI_COMM_WORLD);
#else
    GlobalSubtractedMass = SubtractedMass;
    GlobalSubtractedMomentum = SubtractedMomentum;
#endif
    temp->AccretionRate = AccretionRate;
    
    /* Transfer the mass and momentum to the particle */
    if (SinkIsOnThisProc)
      if (sinkGrid->GridData->AddMassAndMomentumToAccretingParticle(GlobalSubtractedMass, GlobalSubtractedMomentum, 
								    static_cast<ActiveParticleType*>(temp),
								    LevelArray) == FAIL)
	return FAIL;
    
  }
  CommunicationSyncNumberOfParticles(Grids, NumberOfGrids);
  delete [] Grids;
  return SUCCESS;
}

int ActiveParticleType_AccretingParticle::SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level, 
							   int TopGridDims[], int AccretingParticleID)
{
  /* Generate a list of all sink particles in the simulation box */
  int i,dim,nParticles;
  FLOAT *pos = NULL,dx,AccretionRadius;
  ActiveParticleType **AccretingParticleList = NULL ;
  LevelHierarchyEntry *Temp;
  
  AccretingParticleList = ActiveParticleFindAll(LevelArray, &nParticles, AccretingParticleID);
  
  /* Calculate CellWidth on maximum refinement level */
  
  // this will fail for noncubic boxes or simulations with MinimimRefinementLevelExponent
  dx = (DomainRightEdge[0] - DomainLeftEdge[0]) /
    (TopGridDims[0]*POW(FLOAT(RefineBy),FLOAT(MaximumRefinementLevel)));
  
  for (i=0 ; i++ ; i<nParticles){
    pos = AccretingParticleList[i]->ReturnPosition();
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      if (Temp->GridData->DepositAccretionZone(level,pos,AccretionRadius*dx) == FAIL) {
	ENZO_FAIL("Error in grid->DepositAccretionZone.\n")
	  }
  }

  return SUCCESS;
}

int ActiveParticleType_AccretingParticle::AdjustBondiHoyle(grid* CurrentGrid) {
  float *density = CurrentGrid->AccessDensity();
  float *velx = CurrentGrid->AccessVelocity1();
  float *vely = CurrentGrid->AccessVelocity2();
  float *velz = CurrentGrid->AccessVelocity3();
  float *totalenergy = CurrentGrid->AccessTotalEnergy();
  float *gasenergy = CurrentGrid->AccessGasEnergy();

  float *vel = this->ReturnVelocity();
  FLOAT *pos = this->ReturnPosition();
  float CellTemperature;

  int n,i,j,k,index,dim;

  FLOAT LeftCorner[MAX_DIMENSION];

  FLOAT dx = CurrentGrid->GetCellWidth(0);

  int GridStartIndex[3];
  int GridDimension[3];

  for (dim = 0; dim < 3; dim++) {
    LeftCorner[dim] = CurrentGrid->GetGridLeftEdge(dim);
    GridStartIndex[dim] = CurrentGrid->GetGridStartIndex(dim);
    GridDimension[dim] = CurrentGrid->GetGridDimension(dim);
  }

  float *temperature = new float[CurrentGrid->GetGridSize()];
  CurrentGrid->ComputeTemperatureField(temperature);

  i = int((pos[0] - LeftCorner[0])/dx);
  j = int((pos[1] - LeftCorner[1])/dx);
  k = int((pos[2] - LeftCorner[2])/dx);

  index = GRIDINDEX(i,j,k);
  
  this->vInfinity = sqrt(pow((vel[0] - velx[index]),2) +
			 pow((vel[1] - vely[index]),2) +
			 pow((vel[2] - velz[index]),2));
  
  CellTemperature = (JeansRefinementColdTemperature > 0) ? JeansRefinementColdTemperature : temperature[index];
  this->cInfinity = sqrt(Gamma*kboltz*CellTemperature/(Mu*mh))/GlobalLengthUnits*GlobalTimeUnits;
  
  this->BondiHoyleRadius = GravitationalConstant*(this->ReturnMass()*POW(dx,3))/
    (pow(this->vInfinity,2) + pow(this->cInfinity,2));
  
  delete[] temperature;

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
	     &position, MPI_COMM_WORLD);
    MPI_Pack(pbuffer->cInfinity,pbuffer->NumberOfBuffers, FloatDataType, buffer, total_buffer_size,
	     &position, MPI_COMM_WORLD);
    MPI_Pack(pbuffer->vInfinity,pbuffer->NumberOfBuffers, FloatDataType, buffer, total_buffer_size,
	     &position, MPI_COMM_WORLD);
    MPI_Pack(pbuffer->BondiHoyleRadius,pbuffer->NumberOfBuffers, MY_MPIFLOAT, buffer, total_buffer_size,
	     &position, MPI_COMM_WORLD);
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
  if (pbuffer->NumberOfBuffers > 0) {
    MPI_Unpack(mpi_buffer, mpi_buffer_size, &position, pbuffer->AccretionRate,
	       pbuffer->NumberOfBuffers, FloatDataType, MPI_COMM_WORLD);
    MPI_Unpack(mpi_buffer, mpi_buffer_size, &position, pbuffer->cInfinity,
	       pbuffer->NumberOfBuffers, FloatDataType, MPI_COMM_WORLD);
    MPI_Unpack(mpi_buffer, mpi_buffer_size, &position, pbuffer->vInfinity,
	       pbuffer->NumberOfBuffers, FloatDataType, MPI_COMM_WORLD);
    MPI_Unpack(mpi_buffer, mpi_buffer_size, &position, pbuffer->BondiHoyleRadius,
	       pbuffer->NumberOfBuffers, MY_MPIFLOAT, MPI_COMM_WORLD);
  }
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
