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

#include "ActiveParticle_SampleParticle.h"

/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

class SampleParticleGrid : private grid {
  friend class ActiveParticleType_SampleParticle;
};

/* Note that we only refer to SampleParticleGrid here. 
 * Given a grid object, we static cast to get this:
 *
 *    SampleParticleGrid *thisgrid =
 *      static_cast<SampleParticleGrid *>(thisgrid_orig); */

int ActiveParticleType_SampleParticle::InitializeParticleType(void)
{
  return SUCCESS;
}

int ActiveParticleType_SampleParticle::EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  SampleParticleGrid *thisgrid =
    static_cast<SampleParticleGrid *>(thisgrid_orig);
  fprintf(stderr, "Checking formation of sample particles.\n");
  return 0;
}

int ActiveParticleType_SampleParticle::EvaluateFeedback
(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  return SUCCESS;
}

void ActiveParticleType_SampleParticle::DescribeSupplementalData(ActiveParticleFormationDataFlags &flags)
{
  flags.DarkMatterDensity = true;
}

int ActiveParticleType_SampleParticle::WriteToOutput(ActiveParticleType **these_particles, int n, int GridRank, hid_t group_id)
{
  /* Create a new subgroup within the active particle group for active particles of type SampleParticle */
  hid_t SampleParticleGroupID = H5Gcreate(group_id,"SampleParticle",0);

  writeScalarAttribute(SampleParticleGroupID,HDF5_INT,"Number of Sample Particles",&n);  

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
    
  ActiveParticleType_SampleParticle *ParticleToWrite;
  for (i=0;i<n;i++) {
    ParticleToWrite = static_cast<ActiveParticleType_SampleParticle*>(these_particles[i]);
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
		 SampleParticleGroupID, HDF5_FILE_PREC, (VOIDP) Position[dim]);
  }
  
  for (dim = 0; dim < GridRank; dim++) {
    WriteDataset(1,&TempInt,ParticleVelocityLabel[dim],
		  SampleParticleGroupID, HDF5_REAL, (VOIDP) Velocity[dim]);
  }
  
  WriteDataset(1,&TempInt,"mass",SampleParticleGroupID,HDF5_REAL,(VOIDP) Mass);
  WriteDataset(1,&TempInt,"creation_time",SampleParticleGroupID,HDF5_REAL,(VOIDP) BirthTime);
  WriteDataset(1,&TempInt,"dynamical_time",SampleParticleGroupID,HDF5_REAL,(VOIDP) DynamicalTime);
  WriteDataset(1,&TempInt,"metallicity_fraction",SampleParticleGroupID,HDF5_REAL,(VOIDP) Metallicity);
  WriteDataset(1,&TempInt,"identifier",SampleParticleGroupID,HDF5_PINT,(VOIDP) ID);

  /* Clean up */

  for (dim = 0; dim < GridRank; dim++) {
    delete[] Position[dim];
    delete[] Velocity[dim];
  }
  delete[] Mass;
  delete[] BirthTime;
  delete[] DynamicalTime;
  delete[] Metallicity;
  H5Gclose(SampleParticleGroupID);

  return SUCCESS;
}

int ActiveParticleType_SampleParticle::ReadFromOutput(ActiveParticleType **&particles_to_read, int &n, int GridRank, hid_t group_id)
{
  int i,dim;
  hsize_t TempInt;

  hid_t SampleParticleGroupID = H5Gopen(group_id,"SampleParticle");

  readAttribute(SampleParticleGroupID,HDF5_INT,"Number of Sample Particles",&n);

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
		  SampleParticleGroupID, HDF5_FILE_PREC, (VOIDP) Position[dim]);
  }

  for (dim = 0; dim < GridRank; dim++) {
    ReadDataset(1,&TempInt,ParticleVelocityLabel[dim],
		  SampleParticleGroupID, HDF5_REAL, (VOIDP) Velocity[dim]);
  }
  ReadDataset(1,&TempInt,"mass",SampleParticleGroupID,HDF5_R8,(VOIDP) Mass);
  ReadDataset(1,&TempInt,"creation_time",SampleParticleGroupID,HDF5_REAL,(VOIDP) BirthTime);
  ReadDataset(1,&TempInt,"dynamical_time",SampleParticleGroupID,HDF5_REAL,(VOIDP) DynamicalTime);
  ReadDataset(1,&TempInt,"metallicity_fraction",SampleParticleGroupID,HDF5_REAL,(VOIDP) Metallicity);
  ReadDataset(1,&TempInt,"identifier",SampleParticleGroupID,HDF5_PINT,(VOIDP) ID);

  for (i = 0; i < n; i++) {
    ActiveParticleType_SampleParticle *np = new ActiveParticleType_SampleParticle();
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
  H5Gclose(SampleParticleGroupID);
  
  return SUCCESS;
}

int ActiveParticleType_SampleParticle::SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level,
							int TopGridDims[], int SampleParticleID)
{
  return SUCCESS;
}

void SampleParticleBufferHandler::AllocateBuffer(ActiveParticleType **np, int NumberOfParticles, char *buffer,
						 Eint32 total_buffer_size, int &buffer_size, Eint32 &position,
						 int type_num, int proc)
{
  SampleParticleBufferHandler *pbuffer = new SampleParticleBufferHandler(np, NumberOfParticles, type_num, proc);
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

void SampleParticleBufferHandler::UnpackBuffer(char *mpi_buffer, int mpi_buffer_size, int NumberOfParticles,
					       ActiveParticleType **np, int &npart)
{
  int i;
  Eint32 position;
  SampleParticleBufferHandler *pbuffer = new SampleParticleBufferHandler(NumberOfParticles);
  pbuffer->_UnpackBuffer(mpi_buffer, mpi_buffer_size, position);
#ifdef USE_MPI
  if (pbuffer->NumberOfBuffers > 0) {
    // If any extra fields are added in the future, then they would be
    // transferred to the buffer here.
  }
#endif
   /* Convert the particle buffer into active particles */
  
  for (i = 0; i < pbuffer->NumberOfBuffers; i++)
    np[npart++] = new ActiveParticleType_SampleParticle(pbuffer, i);
  delete pbuffer;
  return;
}

namespace {
  ActiveParticleType_info *SampleParticleInfo = 
    register_ptype <ActiveParticleType_SampleParticle, SampleParticleBufferHandler> 
    ("SampleParticle");
}
