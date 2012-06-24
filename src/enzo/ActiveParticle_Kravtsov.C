/***********************************************************************
/
/  KRAVTSOV (2003) STAR FORMATION
/
************************************************************************/

#include <string.h>
#include <map>
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

#ifdef NEW_CONFIG

#include "ParameterControl/ParameterControl.h"
extern Configuration Param;

/* Set default parameter values. */

const char config_kravtsov_particle_defaults[] = 
"### KRAVTSOV STAR PARTICLE DEFAULTS ###\n"
"\n"
"Physics: {\n"
"    ActiveParticles: {\n"
"        Kravtsov: {\n"
"            DensityThreshold = 1e6; # [particles per proper cm^3]\n"
"            StarFormationTimeConstant = 4.0e9; # [years]\n"
"            MinimumStarMass = 1.0e9; # [Msun]\n"
"        };\n"
"    };\n"
"};\n";

#endif


/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

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

class KravtsovGrid : private grid {
    friend class ActiveParticleType_Kravtsov;
};

/* Note that we only refer to SampleParticleGrid here. 
 * Given a grid object, we static case to get this:
 *
 *    SampleParticleGrid *thisgrid =
 *      static_cast<SampleParticleGrid *>(thisgrid_orig); */

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

float ActiveParticleType_Kravtsov::DensityThreshold = FLOAT_UNDEFINED;
float ActiveParticleType_Kravtsov::StarFormationTimeConstant = FLOAT_UNDEFINED;
float ActiveParticleType_Kravtsov::MinimumStarMass = FLOAT_UNDEFINED;

int ActiveParticleType_Kravtsov::InitializeParticleType(void) {

#ifdef NEW_CONFIG

  // Update the parameter config to include the local defaults. Note
  // that this does not overwrite any values previously specified.
  Param.Update(config_kravtsov_particle_defaults);

  // Retrieve parameters from Param structure
  Param.GetScalar(DensityThreshold, "Physics.ActiveParticles.Kravtsov.DensityThreshold");
  Param.GetScalar(StarFormationTimeConstant, "Physics.ActiveParticles.Kravtsov.StarFormationTimeConstant");
  Param.GetScalar(MinimumStarMass, "Physics.ActiveParticles.Kravtsov.MinimumStarMass");

#else

  DensityThreshold = StarMakerOverDensityThreshold;
  StarFormationTimeConstant = StarMakerMinimumDynamicalTime;
  MinimumStarMass = StarMakerMinimumMass;
  
#endif

  return SUCCESS;
}


int ActiveParticleType_Kravtsov::EvaluateFormation
(grid *thisgrid_orig, ActiveParticleFormationData &supp_data)
{
  KravtsovGrid *tg =
    static_cast<KravtsovGrid *>(thisgrid_orig);


  /* Make it pretty */

  float *density = tg->BaryonField[supp_data.DensNum];
//  float *velx = tg->BaryonField[supp_data.Vel1Num];
//  float *vely = tg->BaryonField[supp_data.Vel2Num];
//  float *velz = tg->BaryonField[supp_data.Vel3Num];

  float CellWidthTemp = float(tg->CellWidth[0][0]);

  bool HasMetalField = (supp_data.MetalNum != -1 ||
			supp_data.ColourNum != -1);

  int GridDimension[3] = {tg->GridDimension[0],
                          tg->GridDimension[1],
                          tg->GridDimension[2]};

  float gasfrac, starmass, densthresh, timeconstant;
  int i,j,k,index;
  int NumberOfNewParticles = 0;


  // calculate density threshold.  odthresh is in proper particles per
  // cc and (d1/mproton) gives the mean density of the universe in
  // particles/cm^3 (assuming hydrogen is dominant)

  densthresh = DensityThreshold / (supp_data.DensityUnits / mh);


  // calculate time constant for star formation.  This assumes that
  // the user input is in units of years
  
  timeconstant = StarFormationTimeConstant * 3.156e7 / supp_data.TimeUnits;


  // for each zone, : "star" particle is created if the density
  // exceeds some threshold and this is the highest level of
  // refinement.  That's it.

  for (k = tg->GridStartIndex[2]; k <= tg->GridEndIndex[2]; k++) {
    for (j = tg->GridStartIndex[1]; j <= tg->GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(tg->GridStartIndex[0], j, k);
      for (i = tg->GridStartIndex[0]; i <= tg->GridEndIndex[0]; i++, index++) {
	
	// 0. If no more room for particles, quit.
	if (supp_data.NumberOfNewParticles >=
	    supp_data.MaxNumberOfNewParticles)
          continue;
	
	// 1. Finest level of refinement
	if (tg->BaryonField[tg->NumberOfBaryonFields][index] != 0.0) 
	  continue;
	
	// 2. Density greater than threshold
	if (density[index] < densthresh)
	  continue;
	
	/*
	 * ====================================================================
	 * PARTICLE CREATION
	 * ====================================================================
	 */
	
	ActiveParticleType_Kravtsov *np = new ActiveParticleType_Kravtsov();
	supp_data.NewParticles[supp_data.NumberOfNewParticles++] = np;

	// Make sure that we never give put than 90% of the cell's mass into a star particle
	gasfrac = min( 0.9, tg->dtFixed / timeconstant );
	
	// Calculate star mass in solar masses.  If this is less than
	// the user-defined threshold mass, do NOT make a star in this
	// cell.  This is not exactly in keeping with the spirit of
	// the Kravtsov algorithm, and is somewhat degenerate with the
	// density threshold, but we really don't want millions and
	// millions of star particles.

	starmass = gasfrac * density[index] * supp_data.DensityUnits * pow(supp_data.LengthUnits * CellWidthTemp, 3) /  SolarMass;

	// Do not allow stars with mass less than MinimumStarMass
	if (starmass < MinimumStarMass )
	  continue;

	np->Mass =  starmass / supp_data.MassUnits;

	np->type = Kravtsov;
	np->BirthTime = tg->Time;
	
	np->pos[0] = tg->CellLeftEdge[0][i] + 0.5*tg->CellWidth[0][i];
	np->pos[1] = tg->CellLeftEdge[1][j] + 0.5*tg->CellWidth[1][j];
	np->pos[2] = tg->CellLeftEdge[2][k] + 0.5*tg->CellWidth[2][k];

	/*
	  Star velocities averaged over multiple cells to avoid
	  "runaway star particle" phenomenon imethod = 2 is zeus,
	  otherwise PPM
	*/

	float *tvel = tg->AveragedVelocityAtCell(index, supp_data.DensNum,
						  supp_data.Vel1Num);
	np->vel[0] = tvel[0];
	np->vel[1] = tvel[1];
	np->vel[2] = tvel[2];

	/* Set the metallicity */

	if (HasMetalField)
	  np->Metallicity = supp_data.TotalMetals[index];
	else
	  np->Metallicity = 0.0;

      } // ENDFOR i
    } // ENDFOR j
  } // ENDFOR k
  
  return NumberOfNewParticles;
}

// Feedback
int ActiveParticleType_Kravtsov::EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  return SUCCESS;
}

void ActiveParticleType_Kravtsov::DescribeSupplementalData
(ActiveParticleFormationDataFlags &flags)
{
  flags.UnitConversions = true;
  flags.DataFieldNumbers = true;
  flags.MetalField = true;
}

int ActiveParticleType_Kravtsov::WriteToOutput(ActiveParticleType **these_particles, int n, int GridRank, hid_t group_id)
{
    /* Create a new subgroup within the active particle group for active particles of type Kravtsov */
  hid_t KravtsovGroupID = H5Gcreate(group_id,"Kravtsov",0);

  writeScalarAttribute(KravtsovGroupID,HDF5_INT,"Number of Kravtsov Particles",&n);  

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
  PINT *ID = new PINT[n];

  int i,dim;

  for (dim = 0; dim < GridRank; dim++) {
    Position[dim] = new FLOAT[n];
    Velocity[dim] = new float[n];
  }

  hsize_t TempInt;
  TempInt = n;
    
  ActiveParticleType_Kravtsov *ParticleToWrite;
  for (i=0;i<n;i++) {
    ParticleToWrite = static_cast<ActiveParticleType_Kravtsov*>(these_particles[i]);
    for (dim = 0; dim < GridRank; dim++) {
      Position[dim][i] = ParticleToWrite->pos[dim];
      Velocity[dim][i] = ParticleToWrite->vel[dim];
    }
    Mass[i] = ParticleToWrite->Mass;
    BirthTime[i] = ParticleToWrite->BirthTime;
    DynamicalTime[i] = ParticleToWrite->DynamicalTime;
    Metallicity[i] = ParticleToWrite->Metallicity;
    ID[i] = ParticleToWrite->Identifier;
  }

  for (dim = 0; dim < GridRank; dim++) {
    WriteDataset(1,&TempInt,ParticlePositionLabel[dim],
		 KravtsovGroupID, HDF5_FILE_PREC, (VOIDP) Position[dim]);
  }
  
  for (dim = 0; dim < GridRank; dim++) {
    WriteDataset(1,&TempInt,ParticleVelocityLabel[dim],
		  KravtsovGroupID, HDF5_REAL, (VOIDP) Velocity[dim]);
  }
  
  WriteDataset(1,&TempInt,"mass",KravtsovGroupID,HDF5_REAL,(VOIDP) Mass);
  WriteDataset(1,&TempInt,"creation_time",KravtsovGroupID,HDF5_REAL,(VOIDP) BirthTime);
  WriteDataset(1,&TempInt,"dynamical_time",KravtsovGroupID,HDF5_REAL,(VOIDP) DynamicalTime);
  WriteDataset(1,&TempInt,"metallicity_fraction",KravtsovGroupID,HDF5_REAL,(VOIDP) Metallicity);
  WriteDataset(1,&TempInt,"identifier",KravtsovGroupID,HDF5_PINT,(VOIDP) ID);

  /* Clean up */

  for (dim = 0; dim < GridRank; dim++) {
    delete[] Position[dim];
    delete[] Velocity[dim];
  }
  delete[] Mass;
  delete[] BirthTime;
  delete[] DynamicalTime;
  delete[] Metallicity;
  H5Gclose(KravtsovGroupID);

  return SUCCESS;
}

int ActiveParticleType_Kravtsov::ReadFromOutput(ActiveParticleType **&particles_to_read, int &n, int GridRank, hid_t group_id)
{
  int i,dim;
  hsize_t TempInt;
  
  hid_t KravtsovGroupID = H5Gopen(group_id,"Kravtsov");
  
  readAttribute(KravtsovGroupID,HDF5_INT,"Number of Kravtsov Particles",&n);
  
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
  PINT  *ID = new PINT[n];

  for (dim = 0; dim < GridRank; dim++) {
    Position[dim] = new FLOAT[n];
    Velocity[dim] = new float[n];
  }
  
  TempInt = n;
  
  for (dim = 0; dim < GridRank; dim++) {
    ReadDataset(1,&TempInt,ParticlePositionLabel[dim],
		  KravtsovGroupID, HDF5_FILE_PREC, (VOIDP) Position[dim]);
  }

  for (dim = 0; dim < GridRank; dim++) {
    ReadDataset(1,&TempInt,ParticleVelocityLabel[dim],
		  KravtsovGroupID, HDF5_REAL, (VOIDP) Velocity[dim]);
  }
  ReadDataset(1,&TempInt,"mass",KravtsovGroupID,HDF5_R8,(VOIDP) Mass);
  ReadDataset(1,&TempInt,"creation_time",KravtsovGroupID,HDF5_REAL,(VOIDP) BirthTime);
  ReadDataset(1,&TempInt,"dynamical_time",KravtsovGroupID,HDF5_REAL,(VOIDP) DynamicalTime);
  ReadDataset(1,&TempInt,"metallicity_fraction",KravtsovGroupID,HDF5_REAL,(VOIDP) Metallicity);
  ReadDataset(1,&TempInt,"identifier",KravtsovGroupID,HDF5_PINT,(VOIDP) ID);

  for (i = 0; i < n; i++) {
    ActiveParticleType_Kravtsov *np = new ActiveParticleType_Kravtsov();
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
    particles_to_read[i] = static_cast<ActiveParticleType*>(np);
  }

  delete[] Mass;
  delete[] BirthTime;
  delete[] DynamicalTime;
  delete[] Metallicity;
  delete[] ID;

  for (dim = 0; dim < GridRank; dim++) {
    delete[] Position[dim];
    delete[] Velocity[dim];
  }
  H5Gclose(KravtsovGroupID);


  return SUCCESS;
}

int ActiveParticleType_Kravtsov::BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
							    int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
							    int ThisLevel, int TotalStarParticleCountPrevious[],
							    int KravtsovID)
{
  return SUCCESS;
}

int ActiveParticleType_Kravtsov::AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
							int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
							int ThisLevel, int TotalStarParticleCountPrevious[],
							int KravtsovID)
{
  return SUCCESS;
}

int ActiveParticleType_Kravtsov::SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level,
							int TopGridDims[], int KravtsovID)
{
  return SUCCESS;
}

void KravtsovBufferHandler::AllocateBuffer(ActiveParticleType **np, int NumberOfParticles, char *buffer,
						 Eint32 total_buffer_size, int &buffer_size, Eint32 &position,
						 int type_num, int proc)
{
  KravtsovBufferHandler *pbuffer = new KravtsovBufferHandler(np, NumberOfParticles, type_num, proc);
  pbuffer->_AllocateBuffer(buffer, total_buffer_size, buffer_size, position);
#ifdef USE_MPI
  if (pbuffer->NumberOfBuffers > 0) {
    // If any extra fields are added in the future, then they would be
    // transferred to the buffer here.
  }
#endif
  delete pbuffer;
  return;
}

void KravtsovBufferHandler::UnpackBuffer(char *mpi_buffer, int mpi_buffer_size, int NumberOfParticles,
					       ActiveParticleType **np, int &npart)
{
  int i;
  Eint32 position;
  KravtsovBufferHandler *pbuffer = new KravtsovBufferHandler(NumberOfParticles);
  pbuffer->_UnpackBuffer(mpi_buffer, mpi_buffer_size, position);
#ifdef USE_MPI
  if (pbuffer->NumberOfBuffers > 0) {
    // If any extra fields are added in the future, then they would be
    // transferred to the buffer here.
  }
#endif
   /* Convert the particle buffer into active particles */
  
  for (i = 0; i < pbuffer->NumberOfBuffers; i++)
    np[npart++] = new ActiveParticleType_Kravtsov(pbuffer, i);
  delete pbuffer;
  return;
}

namespace {
  ActiveParticleType_info *KravtsovInfo = 
    register_ptype <ActiveParticleType_Kravtsov, KravtsovBufferHandler> ("Kravtsov");
}
