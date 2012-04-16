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
class AccretingParticleBufferHandler;

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
  ActiveParticleType_AccretingParticle(void) : ActiveParticleType() {};
  ActiveParticleType_AccretingParticle(ParticleBufferHandler *buffer, int index) :
    ActiveParticleType(buffer, index) {};
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
  static int SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level, int ActiveParticleID);
  static int InitializeParticleType();

  ENABLED_PARTICLE_ID_ACCESSOR

  // sink helper routines

  static int MergeAccretingParticles(int nParticles, ActiveParticleType** AccretingParticleList,
				     ActiveParticleType_AccretingParticle** MergedParticles,
				     FLOAT LinkingLength, int ngroups, LevelHierarchyEntry *LevelArray[]);  

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
  AccretingParticleGrid *thisGrid =
    static_cast<AccretingParticleGrid *>(thisgrid_orig);
  

  int GridDimension[3] = {thisGrid->GridDimension[0],
			  thisGrid->GridDimension[1],
			  thisGrid->GridDimension[2]};

  int i,j,k,index,method;

  float *density = thisGrid->BaryonField[data.DensNum];
  float *velx = thisGrid->BaryonField[data.Vel1Num];
  float *vely = thisGrid->BaryonField[data.Vel2Num];
  float *velz = thisGrid->BaryonField[data.Vel3Num];
  float JLSquared = (Gamma*pi*kboltz)/
    (data.DensityUnits*Mu*mh*GravConst*POW(data.LengthUnits,2));
  float CellTemperature = 0;
  float JeansDensity = 0;
  float MassRefinementDensity = 0;
  float DensityThreshold = huge_number;
  float ExtraDensity = 0;

  FLOAT dx = data.LengthUnits*thisGrid->CellWidth[0][0];
  
  bool HasMetalField = (data.MetalNum != -1 || data.ColourNum != -1);
  bool JeansRefinement = false;
  bool MassRefinement = false;

  // determine refinement criteria
  for (method = 0; method < MAX_FLAGGING_METHODS; method++) {
    if (CellFlaggingMethod[method] == 2)
      MassRefinement = true;
    if (CellFlaggingMethod[method] == 6) 
      JeansRefinement = true;
  }

  for (k = thisGrid->GridStartIndex[2]; k <= thisGrid->GridEndIndex[2]; k++) {
    for (j = thisGrid->GridStartIndex[1]; j <= thisGrid->GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(thisGrid->GridStartIndex[0], j, k);
      for (i = thisGrid->GridStartIndex[0]; i <= thisGrid->GridEndIndex[0]; i++, index++) {

	// 0. If no more room for particles, quit
	if (data.NumberOfNewParticles >= data.MaxNumberOfNewParticles)
	  continue;

	// 1. Is the cell on the maximum refinement level?
	if (data.level != MaximumRefinementLevel)
	  continue;

	// 2. Does cell violate the Jeans condition?
	if (JeansRefinement) {
	  CellTemperature = (JeansRefinementColdTemperature > 0) ? JeansRefinementColdTemperature : data.Temperature[index];
	  JeansDensity = OverflowFactor*POW(RefineByJeansLengthSafetyFactor,2)*JLSquared*CellTemperature/POW(dx,2);
	  DensityThreshold = min(DensityThreshold,JeansDensity);
	}
	else if (MassRefinement) {
	  MassRefinementDensity = MinimumMassForRefinement[method]*
	    pow(RefineBy, data.level*MinimumMassForRefinementLevelExponent[method])/POW(dx,3);
	  DensityThreshold = min(DensityThreshold,MassRefinementDensity);
	}
	else 
	  ENZO_FAIL("Error in Accreting Particles: Must refine by jeans length or overdensity!");
	
	if (density[index] <= DensityThreshold) 
	  continue;

	// Passed creation tests, create sink particle


	ActiveParticleType_AccretingParticle *np = new ActiveParticleType_AccretingParticle();
	data.NewParticles[data.NumberOfNewParticles++] = np;

	ExtraDensity = density[index] - DensityThreshold;
	np->Mass = ExtraDensity*POW(dx,3);
	np->type = AccretingParticle;
	np->BirthTime = thisGrid->ReturnTime();
	
	np->pos[0] = thisGrid->CellLeftEdge[0][i] + 0.5*thisGrid->CellWidth[0][i];
	np->pos[0] = thisGrid->CellLeftEdge[1][j] + 0.5*thisGrid->CellWidth[1][j];
	np->pos[0] = thisGrid->CellLeftEdge[2][k] + 0.5*thisGrid->CellWidth[2][k];
	
	float *tvel = thisGrid->AveragedVelocityAtCell(index,data.DensNum,data.Vel1Num);
	  
	np->vel[0] = tvel[0];
	np->vel[1] = tvel[1];
	np->vel[2] = tvel[2];
	
	if (HasMetalField)
	  np->Metallicity = data.TotalMetals[index];
	else
	  np->Metallicity = 0.0;

	// Remove mass from grid

	density[index] = DensityThreshold;

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

  for (dim = 0; dim < 3; dim++) 
    LeftCorner[dim] = thisGrid->CellLeftEdge[dim][0];

  /* Loop over all of the AccretingParticles on this grid and set the bondi-hoyle radius */

  ActiveParticleType_AccretingParticle *dummy = new ActiveParticleType_AccretingParticle();
  int type_num = dummy->GetEnabledParticleID();
  delete dummy;
     
  for (n = 0; n < npart; n++) {
    ActiveParticleType_AccretingParticle *ThisParticle = 
      static_cast<ActiveParticleType_AccretingParticle*>(thisGrid->ActiveParticles[n]);
    if (ThisParticle->ReturnType() != type_num) continue;
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
    ThisParticle->cInfinity = sqrt(Gamma*kboltz*CellTemperature/(Mu*mh));
    ThisParticle->BondiHoyleRadius = GravConst*ThisParticle->ReturnMass()/
      (pow(ThisParticle->vInfinity,2) + pow(ThisParticle->cInfinity,2));
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

  writeScalarAttribute(AccretingParticleGroupID,HDF5_INT,"Number of Accreting particles",&n);  


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
		  AccretingParticleGroupID, HDF5_FILE_REAL, (VOIDP) Velocity[dim]);
  }
  
  WriteDataset(1,&TempInt,"mass",AccretingParticleGroupID,HDF5_FILE_REAL,(VOIDP) Mass);
  WriteDataset(1,&TempInt,"creation_time",AccretingParticleGroupID,HDF5_FILE_REAL,(VOIDP) BirthTime);
  WriteDataset(1,&TempInt,"dynamical_time",AccretingParticleGroupID,HDF5_FILE_REAL,(VOIDP) DynamicalTime);
  WriteDataset(1,&TempInt,"metallicity_fraction",AccretingParticleGroupID,HDF5_FILE_REAL,(VOIDP) Metallicity);
  WriteDataset(1,&TempInt,"identifier",AccretingParticleGroupID,HDF5_PINT,(VOIDP) ID);
  WriteDataset(1,&TempInt,"accretion_rate",AccretingParticleGroupID,HDF5_FILE_REAL,(VOIDP) AccretionRate);
  WriteDataset(1,&TempInt,"c_infinity",AccretingParticleGroupID,HDF5_FILE_REAL,(VOIDP) cInfinity);
  WriteDataset(1,&TempInt,"v_infinity",AccretingParticleGroupID,HDF5_FILE_REAL,(VOIDP) vInfinity);
  WriteDataset(1,&TempInt,"bondi_hoyle_radius",AccretingParticleGroupID,HDF5_FILE_REAL,(VOIDP) BondiHoyleRadius);
 

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
  ReadDataset(1,&TempInt,"accretion_rate",AccretingParticleGroupID,HDF5_FILE_REAL,(VOIDP) AccretionRate);
  ReadDataset(1,&TempInt,"c_infinity",AccretingParticleGroupID,HDF5_FILE_REAL,(VOIDP) cInfinity);
  ReadDataset(1,&TempInt,"v_infinity",AccretingParticleGroupID,HDF5_FILE_REAL,(VOIDP) vInfinity);
  ReadDataset(1,&TempInt,"bondi_hoyle_radius",AccretingParticleGroupID,HDF5_FILE_REAL,(VOIDP) BondiHoyleRadius);


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

int ActiveParticleType_AccretingParticle::MergeAccretingParticles
(int nParticles, ActiveParticleType** ParticleList, ActiveParticleType_AccretingParticle** MergedParticles, 
FLOAT LinkingLength, int ngroups, LevelHierarchyEntry *LevelArray[])
{
  int i,j;
  int dim;
  int GroupNumberAssignment[nParticles];
  int *groupsize = NULL;
  int **grouplist = NULL;
  FLOAT *pos;
  
  
  /* Construct list of sink particle positions to pass to Foflist */
  FLOAT ParticleCoordinates[3*nParticles];
  
  for (i=0 ; i++ ; i<nParticles) {
    pos = ParticleList[i]->ReturnPosition();
    for (dim=0; dim++; dim<3) { ParticleCoordinates[i*nParticles+dim] = pos[dim]; }
  }
  
  /* Find mergeable groups using an FOF search */

  ngroups = FofList(nParticles, ParticleCoordinates, LinkingLength, GroupNumberAssignment, &groupsize, &grouplist);
  
  MergedParticles = new ActiveParticleType_AccretingParticle*[ngroups];

  /* Merge the mergeable groups */

  for (i=0 ; i++ ; i<ngroups) {
    MergedParticles[i] = static_cast<ActiveParticleType_AccretingParticle*>(ParticleList[grouplist[i][0]]);
    if (groupsize[i] != 1) {
      for (j=1 ; j++ ; j<groupsize[i]) {
	MergedParticles[i]->Merge(static_cast<ActiveParticleType_AccretingParticle*>(ParticleList[grouplist[i][j]]));
	ParticleList[grouplist[i][j]]->DisableParticle(LevelArray);
      }
    }
  }

  nParticles = ngroups;

  return SUCCESS;
}



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
      HierarchyEntry **Grids;
      ActiveParticleType** ParticleList;

      ActiveParticleFindAll(LevelArray, ParticleList, nParticles, AccretingParticleID);

      /* Calculate CellWidth on maximum refinement level */

      // This may not work for simulations with MinimumMassForRefinementLevelExponent
      FLOAT dx = (DomainRightEdge[0] - DomainLeftEdge[0])/(POW(FLOAT(RefineBy),FLOAT(MaximumRefinementLevel)));

      /* Do Merging   */

      ActiveParticleType_AccretingParticle **MergedParticles = NULL;
      
      /* Generate new merged list of sink particles */
      
      if (MergeAccretingParticles(nParticles,ParticleList, MergedParticles, 
				  LinkingLength*dx,NumberOfMergedParticles,LevelArray) == FAIL) 
	ENZO_FAIL("Accreting Particle merging failed");
   
      /* Assign local particles to grids */

      NumberOfGrids = GenerateGridArray(LevelArray, ThisLevel, &Grids);
      
      for (grid = 0; grid < NumberOfGrids; grid++) 
	if (Grids[grid]->GridData->ReturnProcessorNumber() == MyProcessorNumber)
	  for (i = 0; i<NumberOfMergedParticles; i++)
	    if (Grids[grid]->GridData->AddActiveParticle(static_cast<ActiveParticleType*>(MergedParticles[i])) == FAIL) 
	      ENZO_FAIL("Cannot assign accreting particles to grid");
	  
      delete [] Grids;
      delete [] ParticleList;
      delete [] MergedParticles;

      /* Regenerate the global active particle list */
      
      ActiveParticleFindAll(LevelArray, ParticleList, nParticles, AccretingParticleID);

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
  /* For each particle, loop over all of the grids and do accretion 
     if the grid overlaps with the accretion zone                   */
  
  int i, grid, NumberOfGrids;
  HierarchyEntry **Grids;
  HierarchyEntry *sinkGrid;

  bool SinkIsOnThisGrid;

  float WeightedSum = 0, SumOfWeights = 0, GlobalWeightedSum, GlobalSumOfWeights, AverageDensity, SubtractedMass, 
    GlobalSubtractedMass, SubtractedMomentum[3], GlobalSubtractedMomentum[3], vInfinity, cInfinity, BondiHoyleRadius, 
    AccretionRate;
  int NumberOfCells = 0;

  NumberOfGrids = GenerateGridArray(LevelArray, ThisLevel, &Grids);

  for (i = 0; i<nParticles; i++) {
    ActiveParticleType_AccretingParticle* temp = static_cast<ActiveParticleType_AccretingParticle*>(ParticleList[i]);
    vInfinity = temp->vInfinity;
    cInfinity = temp->cInfinity;
    BondiHoyleRadius = temp->BondiHoyleRadius;
    for (grid = 0; grid < NumberOfGrids; grid++) {
      if (Grids[grid]->GridData->
	  FindAverageDensityInAccretionZone(ParticleList[i],AccretionRadius, WeightedSum, 
					    SumOfWeights, NumberOfCells, BondiHoyleRadius) == FAIL) {
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
      
      /* Now perform accretion algorithm by modifying the grids locally */
      if (Grids[grid]->GridData->
	  AccreteOntoAccretingParticle(ParticleList[i], AccretionRadius, AverageDensity, SumOfWeights, SubtractedMass, 
				       SubtractedMomentum, SinkIsOnThisGrid, vInfinity, cInfinity, 
				       BondiHoyleRadius, AccretionRate) == FAIL) {
	return FAIL;
      }
      if (SinkIsOnThisGrid)
	sinkGrid = Grids[grid];
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
    if (sinkGrid->GridData->AddMassAndMomentumToAccretingParticle(GlobalSubtractedMass, GlobalSubtractedMomentum, 
								  static_cast<ActiveParticleType*>(temp),
								  LevelArray) == FAIL)
      return FAIL;
    
  }
  return SUCCESS;
}



class AccretingParticleBufferHandler : public ParticleBufferHandler
{
public:
  AccretingParticleBufferHandler(void) : ParticleBufferHandler() {};
  AccretingParticleBufferHandler(int NumberOfParticles) : ParticleBufferHandler(NumberOfParticles) {
    // Any extra fields must be added to the buffer
    this->AccretionRate = new float[NumberOfParticles];
    this->cInfinity = new float[NumberOfParticles];
    this->vInfinity = new float[NumberOfParticles];
    this->BondiHoyleRadius = new FLOAT[NumberOfParticles];
  };
  AccretingParticleBufferHandler(ActiveParticleType **np, int NumberOfParticles, int type, int proc) : 
    ParticleBufferHandler(np, NumberOfParticles, type, proc) {
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
    this->ElementSizeInBytes += 3*sizeof(float) + sizeof(FLOAT);
  };
  ~AccretingParticleBufferHandler() {};
  static void AllocateBuffer(ActiveParticleType **np, int NumberOfParticles, char *buffer, 
			     Eint32 total_buffer_size, int &buffer_size, Eint32 &position, 
			     int type_num, int proc);
  static void UnpackBuffer(char *mpi_buffer, int mpi_buffer_size, int NumberOfParticles,
			   ActiveParticleType **np, int &npart);
  float* AccretionRate;
  float* cInfinity;
  float* vInfinity;
  FLOAT* BondiHoyleRadius;
};

int ActiveParticleType_AccretingParticle::SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level, int AccretingParticleID)
{
  /* Generate a list of all sink particles in the simulation box */
  int i,dim,nParticles;
  FLOAT *pos,dx,AccretionRadius;
  ActiveParticleType_AccretingParticle **AccretingParticleList ;
  LevelHierarchyEntry *Temp;
  
  //ActiveParticleFindAll(LevelArray, AccretingParticleList, nParticles, AccretingParticleID);
  
  /* Calculate CellWidth on maximum refinement level */
  
  dx = (DomainRightEdge[0] - DomainLeftEdge[0])/(POW(FLOAT(RefineBy),FLOAT(MaximumRefinementLevel)));
  
  for (i=0 ; i++ ; i<nParticles){
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
    MPI_Pack(pbuffer->AccretionRate,pbuffer->NumberOfBuffers, FloatDataType, buffer, buffer_size,
	     &position, MPI_COMM_WORLD);
    MPI_Pack(pbuffer->cInfinity,pbuffer->NumberOfBuffers, FloatDataType, buffer, buffer_size,
	     &position, MPI_COMM_WORLD);
    MPI_Pack(pbuffer->vInfinity,pbuffer->NumberOfBuffers, FloatDataType, buffer, buffer_size,
	     &position, MPI_COMM_WORLD);
    MPI_Pack(pbuffer->BondiHoyleRadius,pbuffer->NumberOfBuffers, MY_MPIFLOAT, buffer, buffer_size,
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
  return;
}

namespace {
  ActiveParticleType_info *AccretingParticleInfo = 
    register_ptype <ActiveParticleType_AccretingParticle, AccretingParticleBufferHandler> 
    ("AccretingParticle");
}
