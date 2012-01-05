/***********************************************************************
/
/ Accreting Sink Particle
/
************************************************************************/

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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "EventHooks.h"
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
"        SinkParticle: {\n"
"            OverflowFactor       = 1.01;\n"
"            LinkingLength        = 4;\n   "
"            AccretionRadius      = 4;\n   "
"        };\n"
"    };\n"
"};\n";

#endif

/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

class ActiveParticleType_SinkParticle;
class SinkParticleBufferHandler;

class SinkParticleGrid : private grid {
  friend class ActiveParticleType_SinkParticle;
};

/* Note that we only refer to SinkParticleGrid here. 
 * Given a grid object, we static case to get this:
 *
 *    SinkParticleGrid *thisgrid =
 *      static_cast<SinkParticleGrid *>(thisgrid_orig); */

class ActiveParticleType_SinkParticle : public ActiveParticleType
{
public:
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static int WriteToOutput(ActiveParticleType *these_particles, int n, int GridRank, hid_t group_id);
  static int ReadFromOutput(ActiveParticleType **particles_to_read, int *n, int GridRank, hid_t group_id);
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);
  static ParticleBufferHandler *AllocateBuffers(int NumberOfParticles);
  static int EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static int BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
			       int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			       int ThisLevel, int TotalStarParticleCountPrevious[],
			       int SinkParticleID);
  static int AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
			      int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			      int ThisLevel, int TotalStarParticleCountPrevious[],
			      int SinkParticleID);
  static int InitializeParticleType();

  ENABLED_PARTICLE_ID_ACCESSOR

  // sink helper routines

  static int MergeSinks(int nParticles, ActiveParticleType_SinkParticle** SinkParticleList, 
			FLOAT LinkingLength, int ngroups, LevelHierarchyEntry *LevelArray[]);  

  static int SinkAccrete(int nParticles, ActiveParticleType_SinkParticle** SinkParticleList,
			 FLOAT AccretionRadius, LevelHierarchyEntry *LevelArray[]);
  

  static float OverflowFactor;
  static int AccretionRadius;   // in units of CellWidth on the maximum refinement level
  static int LinkingLength;     // Should be equal to AccretionRadius
  float AccretionRate;
  
};

float ActiveParticleType_SinkParticle::OverflowFactor = FLOAT_UNDEFINED;
int ActiveParticleType_SinkParticle::AccretionRadius = INT_UNDEFINED;
int ActiveParticleType_SinkParticle::LinkingLength = INT_UNDEFINED;

int ActiveParticleType_SinkParticle::InitializeParticleType()
{
#ifdef NEW_CONFIG

  // Update the parameter config to include the local defaults.  Note
  // that this does not overwrite any values previously specified.
  Param.Update(config_sink_particle_defaults);

  // Retrieve parameters from Param structure
  Param.GetScalar(OverflowFactor, "Physics.ActiveParticles.SinkParticle.OverflowFactor");
  Param.GetScalar(LinkingLength, "Physics.ActiveParticles.SinkParticle.LinkingLength");
  Param.GetScalar(AccretionRadius, "Physics.ActiveParticles.SinkParticle.AccretionRadius");

#else

  // Leaving these defaults hardcoded for testing. NJG

  OverflowFactor = 1.01;
  LinkingLength = 4;
  AccretionRadius = 4;

#endif

  return SUCCESS;
}

int ActiveParticleType_SinkParticle::EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  SinkParticleGrid *thisGrid =
    static_cast<SinkParticleGrid *>(thisgrid_orig);
  

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
	  JeansDensity = POW(RefineByJeansLengthSafetyFactor,2)*JLSquared*CellTemperature/POW(dx,2);
	  DensityThreshold = min(DensityThreshold,JeansDensity);
	}
	else if (MassRefinement) {
	  MassRefinementDensity = MinimumMassForRefinement[method]*
	    pow(RefineBy, data.level*MinimumMassForRefinementLevelExponent[method])/POW(dx,3);
	  DensityThreshold = min(DensityThreshold,MassRefinementDensity);
	}
	else 
	  ENZO_FAIL("Error in Sink Particles: Must refine by jeans length or overdensity!");
	
	if (density[index] <= DensityThreshold) 
	  continue;

	// Passed creation tests, create sink particle


	ActiveParticleType_SinkParticle *np = new ActiveParticleType_SinkParticle();
	data.NewParticles[data.NumberOfNewParticles++] = np;

	ExtraDensity = density[index] - DensityThreshold;
	np->Mass = ExtraDensity*POW(dx,3);
	np->type = SinkParticle;
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

int ActiveParticleType_SinkParticle::EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  return SUCCESS;
}

void ActiveParticleType_SinkParticle::DescribeSupplementalData(ActiveParticleFormationDataFlags &flags)
{
  flags.DarkMatterDensity = true;
  flags.Temperature = true;
  flags.UnitConversions = true;
  flags.DataFieldNumbers = true;
  flags.MetalField = true;
}


int ActiveParticleType_SinkParticle::WriteToOutput(ActiveParticleType *these_particles, int n, int GridRank, hid_t group_id)
{
  hid_t SinkParticleGroupID = H5Gcreate(group_id,"SinkParticle",0);

  writeScalarAttribute(SinkParticleGroupID,HDF5_INT,"number_of_active_particles_of_this_type",&n);  

  ActiveParticleType_SinkParticle *ParticlesToWrite = static_cast<ActiveParticleType_SinkParticle *>(these_particles);

    char *ParticlePositionLabel[] =
     {"position_x", "position_y", "position_z"};
  char *ParticleVelocityLabel[] =
     {"velocity_x", "velocity_y", "velocity_z"};

  /* Create temporary buffers to store particle data */

  FLOAT Position[GridRank][n];
  float Velocity[GridRank][n]; 
  double Mass[n];
  float BirthTime[n];
  float DynamicalTime[n];
  float Metallicity[n];
  
  FLOAT *pos;
  float *vel;

  int i,dim;

  hsize_t TempInt;
  TempInt = n;
    
  for (i=0;i<n;i++) {
    pos = ParticlesToWrite[i].ReturnPosition();
    vel = ParticlesToWrite[i].ReturnVelocity();
        for (dim = 0; dim < GridRank; dim++) {
      Position[dim][i] = pos[dim];
      Velocity[dim][i] = vel[dim];
    }
    Mass[i] = ParticlesToWrite[i].ReturnMass();
    BirthTime[i] = ParticlesToWrite[i].ReturnBirthTime();
    DynamicalTime[i] = ParticlesToWrite[i].ReturnDynamicalTime();
    Metallicity[i] = ParticlesToWrite[i].ReturnMetallicity();
  }

  for (dim = 0; dim < GridRank; dim++) {
    WriteDataset(1,&TempInt,ParticlePositionLabel[dim],
		 SinkParticleGroupID, HDF5_FILE_PREC, (VOIDP) Position[dim]);
  }
  
  for (dim = 0; dim < GridRank; dim++) {
    WriteDataset(1,&TempInt,ParticleVelocityLabel[dim],
		  SinkParticleGroupID, HDF5_FILE_REAL, (VOIDP) Velocity[dim]);
  }
  
  WriteDataset(1,&TempInt,"mass",SinkParticleGroupID,HDF5_FILE_REAL,(VOIDP) Mass);
  WriteDataset(1,&TempInt,"creation_time",SinkParticleGroupID,HDF5_FILE_REAL,(VOIDP) BirthTime);
  WriteDataset(1,&TempInt,"dynamical_time",SinkParticleGroupID,HDF5_FILE_REAL,(VOIDP) DynamicalTime);
  WriteDataset(1,&TempInt,"metallicity_fraction",SinkParticleGroupID,HDF5_FILE_REAL,(VOIDP) Metallicity);

  return SUCCESS;
}

int ActiveParticleType_SinkParticle::ReadFromOutput(ActiveParticleType **particles_to_read, int *n, int GridRank, hid_t group_id)
{

  hid_t SinkParticleGroupID = H5Gopen(group_id,"SinkParticle");

  readAttribute(SinkParticleGroupID,HDF5_INT,"number_of_active_particles_of_this_type",n);

  char *ParticlePositionLabel[] =
     {"position_x", "position_y", "position_z"};
  char *ParticleVelocityLabel[] =
     {"velocity_x", "velocity_y", "velocity_z"};

  FLOAT Position[GridRank][*n];
  float Velocity[GridRank][*n]; 
  double Mass[*n];
  float BirthTime[*n];
  float DynamicalTime[*n];
  float Metallicity[*n];
  
  int i,dim;

  hsize_t TempInt;
  TempInt = *n;
  
  for (dim = 0; dim < GridRank; dim++) {
    ReadDataset(1,&TempInt,ParticlePositionLabel[dim],
		  SinkParticleGroupID, HDF5_FILE_PREC, (VOIDP) Position[dim]);
  }

  for (dim = 0; dim < GridRank; dim++) {
    ReadDataset(1,&TempInt,ParticleVelocityLabel[dim],
		  SinkParticleGroupID, HDF5_FILE_REAL, (VOIDP) Velocity[dim]);
  }
  ReadDataset(1,&TempInt,"mass",SinkParticleGroupID,HDF5_FILE_REAL,(VOIDP) Mass);
  ReadDataset(1,&TempInt,"creation_time",SinkParticleGroupID,HDF5_FILE_REAL,(VOIDP) BirthTime);
  ReadDataset(1,&TempInt,"dynamical_time",SinkParticleGroupID,HDF5_FILE_REAL,(VOIDP) DynamicalTime);
  ReadDataset(1,&TempInt,"metallicity_fraction",SinkParticleGroupID,HDF5_FILE_REAL,(VOIDP) Metallicity);

  for (i = 0; i < *n; i++) {
    ActiveParticleType_SinkParticle *np = new ActiveParticleType_SinkParticle;
    particles_to_read[i] = np;
    np = new ActiveParticleType_SinkParticle();
    np->Mass = Mass[i];
    np->type = SinkParticle;
    np->BirthTime = DynamicalTime[i];
    np->Metallicity = Metallicity[i];
    for (dim = 0; dim < GridRank; dim++){
      np->pos[dim] = Position[dim][i];
      np->vel[dim] = Velocity[dim][i];
    }
  }
  return SUCCESS;
  return SUCCESS;
}

int ActiveParticleType_SinkParticle::BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
						       int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
						       int ThisLevel, int TotalStarParticleCountPrevious[],
						       int SinkParticleID)
{
  /* Sink particles live on the maximum refinement level.  If we are on a lower level, this does not concern us */
  
  if (ThisLevel == MaximumRefinementLevel)
    {
      /* Generate list of all sink particles in the simulation box */
      
      int nParticles,NumberOfMergedParticles; 
      ActiveParticleType_SinkParticle **SinkParticleList;

      // A function that does the communication work (e.g. StarParticleFindAll) should go here.  JHW is working on this
      
      /* Calculate CellWidth on maximum refinement level */
      
      FLOAT dx = (DomainRightEdge[0] - DomainLeftEdge[0])/(POW(FLOAT(RefineBy),FLOAT(MaximumRefinementLevel)));
      
      /* Generate new merged list of sink particles */
      
      if (MergeSinks(nParticles,SinkParticleList,LinkingLength*dx,NumberOfMergedParticles,LevelArray) == FAIL) {
	ENZO_FAIL("Sink Particle merging failed");
      }
   
      /* Broadcast new global list of active particles */

      // A function that does the communication work should go here. JHW is working on this.
        
      
    } // ENDIF: ThisLevel == MaximumRefinementLevel

  return SUCCESS;
}

int ActiveParticleType_SinkParticle::MergeSinks(int nParticles, ActiveParticleType_SinkParticle** SinkParticleList, FLOAT LinkingLength,int ngroups, LevelHierarchyEntry *LevelArray[])
{
  int i,j;
  int dim;
  int GroupNumberAssignment[nParticles];
  int *groupsize = NULL;
  int **grouplist = NULL;
  FLOAT *pos;
  ActiveParticleType_SinkParticle **NewParticles;
  
  
  /* Construct list of sink particle positions to pass to Foflist */
  FLOAT SinkCoordinates[3*nParticles];
  
  for (i=0 ; i++ ; i<nParticles) {
    pos = SinkParticleList[i]->ReturnPosition();
    for (dim=0; dim++; dim<3) { SinkCoordinates[i*nParticles+dim] = pos[dim]; }
  }
  
  /* Find mergeable groups using an FOF search */

  ngroups = FofList(nParticles, SinkCoordinates, LinkingLength, GroupNumberAssignment, &groupsize, &grouplist);
  
  /* Merge the mergeable groups */

  for (i=0 ; i++ ; i<ngroups) {
    NewParticles[i] = SinkParticleList[grouplist[i][0]];
    if (groupsize[i] != 1) {
      for (j=1 ; j++ ; j<groupsize[i]) {
	NewParticles[i]->Merge(SinkParticleList[grouplist[i][j]]);
	SinkParticleList[grouplist[i][j]]->DisableParticle(LevelArray);
      }
    }
  }

  delete [] SinkParticleList;

  SinkParticleList = NewParticles;

  delete [] NewParticles;

  nParticles = ngroups;

  return SUCCESS;
}



int ActiveParticleType_SinkParticle::AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
						      int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
						      int ThisLevel, int TotalStarParticleCountPrevious[],
						      int SinkParticleID)
{

  /* Sink particles live on the maximum refinement level.  If we are on a lower level, this does not concern us */

  if (ThisLevel == MaximumRefinementLevel)
    {

      /* Generate a list of all sink particles in the simulation box */
      int i,nParticles;
      ActiveParticleType_SinkParticle **SinkParticleList ;

      // A function that does the communication work (e.g. StarParticleFindAll) should go here.  JHW is working on this

      /* Calculate CellWidth on maximum refinement level */

      FLOAT dx = (DomainRightEdge[0] - DomainLeftEdge[0])/(POW(FLOAT(RefineBy),FLOAT(MaximumRefinementLevel)));

      /* Do accretion */
      
      if (SinkAccrete(nParticles,SinkParticleList,AccretionRadius*dx,LevelArray) == FAIL){
	ENZO_FAIL("Sink Particle accretion failed")
	  }

    }

  return SUCCESS;
}

int ActiveParticleType_SinkParticle::SinkAccrete(int nParticles, ActiveParticleType_SinkParticle** SinkParticleList,
						 FLOAT AccretionRadius, LevelHierarchyEntry *LevelArray[])
{


  return SUCCESS;
}

class SinkParticleBufferHandler : public ParticleBufferHandler
{
  public:
    SinkParticleBufferHandler(int NumberOfParticles) { }
};

ParticleBufferHandler *ActiveParticleType_SinkParticle::AllocateBuffers(int NumberOfParticles)
{
    SinkParticleBufferHandler *handler = new SinkParticleBufferHandler(NumberOfParticles);
    return handler;
}


namespace {
  ActiveParticleType_info *SinkParticleInfo = 
    register_ptype <ActiveParticleType_SinkParticle> ("SinkParticle");
}
