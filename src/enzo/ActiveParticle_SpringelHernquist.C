/***********************************************************************
/
/  SPRINGEL & HERNQUIST STAR FORMATION
/
/  written by: Stephen Skory
/  date:       October, 2011
/
/  PURPOSE:
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

#include <limits.h>
#include "phys_constants.h"

/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

class ActiveParticleType_SpringelHernquist;
class SpringelHernquistBufferHandler : public ParticleBufferHandler
{
public:
  // No extra fields in SpringelHernquist.  Same base constructor.
  SpringelHernquistBufferHandler(void) : ParticleBufferHandler() {};

  SpringelHernquistBufferHandler(int NumberOfParticles) :
    ParticleBufferHandler(NumberOfParticles) {};
  
  SpringelHernquistBufferHandler(ActiveParticleType **np, int NumberOfParticles, int type, int proc) :
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

class SpringelHernquistGrid : private grid {
  friend class ActiveParticleType_SpringelHernquist;
};

/* Note that we only refer to SpringelHernquistGrid here. 
 * Given a grid object, we static case to get this:
 *
 *    SpringelHernquistGrid *thisgrid =
 *      static_cast<SpringelHernquistGrid *>(thisgrid_orig); */

class ActiveParticleType_SpringelHernquist : public ActiveParticleType
{
public:
  // Constructors
  ActiveParticleType_SpringelHernquist(void) : ActiveParticleType() {};

  ActiveParticleType_SpringelHernquist(SpringelHernquistBufferHandler *buffer, int index) :
    ActiveParticleType(static_cast<ParticleBufferHandler*>(buffer), index) {};

  // Static members
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static int WriteToOutput(ActiveParticleType **these_particles, int n, int GridRank, hid_t group_id);
  static int ReadFromOutput(ActiveParticleType **&particles_to_read, int &n, int GridRank, hid_t group_id);
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);
  static int EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static int BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
			       int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			       int ThisLevel, int TotalStarParticleCountPrevious[],
			       int SpringelHernquistID);
  static int AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
			      int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			      int ThisLevel, int TotalStarParticleCountPrevious[],
			      int SpringelHernquistID);
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
  static float OverDensityThreshold, PhysicalDensityThreshold, 
    MinimumDynamicalTime, MinimumMass;
};

float ActiveParticleType_SpringelHernquist::OverDensityThreshold = FLOAT_UNDEFINED;
float ActiveParticleType_SpringelHernquist::PhysicalDensityThreshold = FLOAT_UNDEFINED;
float ActiveParticleType_SpringelHernquist::MinimumDynamicalTime = FLOAT_UNDEFINED;
float ActiveParticleType_SpringelHernquist::MinimumMass = FLOAT_UNDEFINED;

int ActiveParticleType_SpringelHernquist::InitializeParticleType() {
  // get some parameters from the Param object

#ifdef NEW_CONFIG
  Param.GetScalar(OverDensityThreshold, "Physics.ActiveParticles.SpringelHernquist.OverDensityThreshold");
  Param.GetScalar(PhysicalDensityThreshold, "Physics.ActiveParticles.SpringelHernquist.PhysicalDensityThreshold");
  Param.GetScalar(MinimumDynamicalTime, "Physics.ActiveParticles.SpringelHernquist.MinimumDynamicalTime");
  Param.GetScalar(MinimumMass, "Physics.ActiveParticles.SpringelHernquist.MinimumMass");
#else
  OverDensityThreshold = StarMakerOverDensityThreshold;
  PhysicalDensityThreshold = StarMakerSHDensityThreshold;
  MinimumDynamicalTime = StarMakerMinimumDynamicalTime;
  MinimumMass = StarMakerMinimumMass;
#endif

}

int ActiveParticleType_SpringelHernquist::EvaluateFormation
(grid *thisgrid_orig, ActiveParticleFormationData &supp_data)
{
  SpringelHernquistGrid *tg =
    static_cast<SpringelHernquistGrid *>(thisgrid_orig);

  float bmass, tstar;
  int i, j, k, index, offset_y, offset_z;
  int NumberOfNewParticles = 0;
  float r_float, y, x, pstar, starfraction, usn;
  float msolar=SolarMass, mproton=mh, beta=0.1,
    sqrtepssn=2.0e24, kb=kboltz;

  // Define the energy output from supernovae

  usn = (1. - beta) * sqrtepssn / beta / msolar; // Below Eq (2)
  usn = usn * sqrtepssn; // this trick done to avoid floating-point issues (BWO)

  /* Make it pretty */

  float *density = tg->BaryonField[supp_data.DensNum];
  float *velx = tg->BaryonField[supp_data.Vel1Num];
  float *vely = tg->BaryonField[supp_data.Vel2Num];
  float *velz = tg->BaryonField[supp_data.Vel3Num];

  bool HasMetalField = (supp_data.MetalNum != -1 ||
			supp_data.ColourNum != -1);

  int GridDimension[3] = {tg->GridDimension[0],
                          tg->GridDimension[1],
                          tg->GridDimension[2]};

  // Pre-calculate serialized offsets for the 3D data field.  Used for
  // the divergence.
  offset_y = tg->GridDimension[0];
  offset_z = tg->GridDimension[0] * tg->GridDimension[1];

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
	if (density[index] < OverDensityThreshold)
	  continue;

	// Calculate star formation timescale. Eq 21.
	tstar =  31556926.0 * MinimumDynamicalTime * 
	  pow(density[index] * supp_data.DensityUnits / PhysicalDensityThreshold, -0.5);

	/* note: minus sign is because 'coolrate' is actual backwards: when coolrate is 
	 * negative, gas is cooling (which is the opposite of how it's defined in Springel
	 * & Hernquist 2003, MNRAS, 339, 289, eqtn. 16 */
	y = -1.0 * tstar * supp_data.CoolingRate[index] /
	  (density[index] * supp_data.DensityUnits) /
	  (beta * usn - (1.0 - beta) * 1.5 *
	   kb * supp_data.Temperature[index] / 0.6 / mproton);

	// Calculate the fraction of mass in cold clouds.
	if (y <= 0) {
	  x = 0;
	} else {
	  x = 1. + 1./ 2. / y - sqrt(1. / y + 1. / 4. / y / y); // Eq(18)
	}

	// Calculate total baryon mass in the cell.
	bmass = density[index] * supp_data.MassUnits;

	// Calculate a parameter which is compared to a random number.
	pstar = (bmass / MinimumMass) * 
	  (1. - exp(-(1. - beta) * x * tg->dtFixed * 
		    supp_data.TimeUnits / tstar)); // Eq(39)

	// Make the random number.
	r_float = float(rand())/float(RAND_MAX);

	// Finally, if this random number is smaller than pstar, make a star.
	if (r_float > pstar) continue;
	
	/*
	 * ====================================================================
	 * PARTICLE CREATION
	 * ====================================================================
	 */

	starfraction = min(MinimumMass / bmass, 0.5);

	ActiveParticleType_SpringelHernquist *np = new ActiveParticleType_SpringelHernquist();
	supp_data.NewParticles[supp_data.NumberOfNewParticles++] = np;
	//fprintf(stderr, "G_APH: Creating !\n");

	// Give the star mass and then remove the same from the grid.
	np->Mass = starfraction * density[index];
	density[index] -= starfraction * density[index];

	np->type = PARTICLE_TYPE_STAR;
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

// SH star feedback
int ActiveParticleType_SpringelHernquist::EvaluateFeedback
(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  return SUCCESS;
}

void ActiveParticleType_SpringelHernquist::DescribeSupplementalData
(ActiveParticleFormationDataFlags &flags)
{
  flags.CoolingTime = true;
  flags.Temperature = true;
  flags.UnitConversions = true;
  flags.DataFieldNumbers = true;
  flags.MetalField = true;
}

int ActiveParticleType_SpringelHernquist::WriteToOutput(ActiveParticleType **these_particles, int n, int GridRank, hid_t group_id)
{
  /* Create a new subgroup within the active particle group for active particles of type SpringelHernquist */
  hid_t SpringelHernquistGroupID = H5Gcreate(group_id,"SpringelHernquist",0);

  writeScalarAttribute(SpringelHernquistGroupID,HDF5_INT,"Number of Sample Particles",&n);  

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
    
  ActiveParticleType_SpringelHernquist *ParticleToWrite;
  for (i=0;i<n;i++) {
    ParticleToWrite = static_cast<ActiveParticleType_SpringelHernquist*>(these_particles[i]);
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
		 SpringelHernquistGroupID, HDF5_FILE_PREC, (VOIDP) Position[dim]);
  }
  
  for (dim = 0; dim < GridRank; dim++) {
    WriteDataset(1,&TempInt,ParticleVelocityLabel[dim],
		  SpringelHernquistGroupID, HDF5_REAL, (VOIDP) Velocity[dim]);
  }
  
  WriteDataset(1,&TempInt,"mass",SpringelHernquistGroupID,HDF5_REAL,(VOIDP) Mass);
  WriteDataset(1,&TempInt,"creation_time",SpringelHernquistGroupID,HDF5_REAL,(VOIDP) BirthTime);
  WriteDataset(1,&TempInt,"dynamical_time",SpringelHernquistGroupID,HDF5_REAL,(VOIDP) DynamicalTime);
  WriteDataset(1,&TempInt,"metallicity_fraction",SpringelHernquistGroupID,HDF5_REAL,(VOIDP) Metallicity);
  WriteDataset(1,&TempInt,"identifier",SpringelHernquistGroupID,HDF5_PINT,(VOIDP) ID);

  /* Clean up */

  for (dim = 0; dim < GridRank; dim++) {
    delete[] Position[dim];
    delete[] Velocity[dim];
  }
  delete[] Mass;
  delete[] BirthTime;
  delete[] DynamicalTime;
  delete[] Metallicity;
  H5Gclose(SpringelHernquistGroupID);

  return SUCCESS;
}

int ActiveParticleType_SpringelHernquist::ReadFromOutput(ActiveParticleType **&particles_to_read, int &n, int GridRank, hid_t group_id)
{
  int i,dim;
  hsize_t TempInt;

  hid_t SpringelHernquistGroupID = H5Gopen(group_id,"SpringelHernquist");

  readAttribute(SpringelHernquistGroupID,HDF5_INT,"Number of Sample Particles",&n);

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
		  SpringelHernquistGroupID, HDF5_FILE_PREC, (VOIDP) Position[dim]);
  }

  for (dim = 0; dim < GridRank; dim++) {
    ReadDataset(1,&TempInt,ParticleVelocityLabel[dim],
		  SpringelHernquistGroupID, HDF5_REAL, (VOIDP) Velocity[dim]);
  }
  ReadDataset(1,&TempInt,"mass",SpringelHernquistGroupID,HDF5_R8,(VOIDP) Mass);
  ReadDataset(1,&TempInt,"creation_time",SpringelHernquistGroupID,HDF5_REAL,(VOIDP) BirthTime);
  ReadDataset(1,&TempInt,"dynamical_time",SpringelHernquistGroupID,HDF5_REAL,(VOIDP) DynamicalTime);
  ReadDataset(1,&TempInt,"metallicity_fraction",SpringelHernquistGroupID,HDF5_REAL,(VOIDP) Metallicity);
  ReadDataset(1,&TempInt,"identifier",SpringelHernquistGroupID,HDF5_PINT,(VOIDP) ID);

  for (i = 0; i < n; i++) {
    ActiveParticleType_SpringelHernquist *np = new ActiveParticleType_SpringelHernquist();
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
  H5Gclose(SpringelHernquistGroupID);
  
  return SUCCESS;
}

int ActiveParticleType_SpringelHernquist::BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
							    int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
							    int ThisLevel, int TotalStarParticleCountPrevious[],
							    int SpringelHernquistID)
{
  return SUCCESS;
}

int ActiveParticleType_SpringelHernquist::AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
							int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
							int ThisLevel, int TotalStarParticleCountPrevious[],
							int SpringelHernquistID)
{
  return SUCCESS;
}

int ActiveParticleType_SpringelHernquist::SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level,
							int TopGridDims[], int SpringelHernquistID)
{
  return SUCCESS;
}

void SpringelHernquistBufferHandler::AllocateBuffer(ActiveParticleType **np, int NumberOfParticles, char *buffer,
						 Eint32 total_buffer_size, int &buffer_size, Eint32 &position,
						 int type_num, int proc)
{
  SpringelHernquistBufferHandler *pbuffer = new SpringelHernquistBufferHandler(np, NumberOfParticles, type_num, proc);
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

void SpringelHernquistBufferHandler::UnpackBuffer(char *mpi_buffer, int mpi_buffer_size, int NumberOfParticles,
					       ActiveParticleType **np, int &npart)
{
  int i;
  Eint32 position;
  SpringelHernquistBufferHandler *pbuffer = new SpringelHernquistBufferHandler(NumberOfParticles);
  pbuffer->_UnpackBuffer(mpi_buffer, mpi_buffer_size, position);
#ifdef USE_MPI
  if (pbuffer->NumberOfBuffers > 0) {
    // If any extra fields are added in the future, then they would be
    // transferred to the buffer here.
  }
#endif
   /* Convert the particle buffer into active particles */
  
  for (i = 0; i < pbuffer->NumberOfBuffers; i++)
    np[npart++] = new ActiveParticleType_SpringelHernquist(pbuffer, i);
  delete pbuffer;
  return;
}

namespace {
  ActiveParticleType_info *SpringelHernquistInfo = 
    register_ptype <ActiveParticleType_SpringelHernquist, SpringelHernquistBufferHandler> 
    ("SpringelHernquist");
}
