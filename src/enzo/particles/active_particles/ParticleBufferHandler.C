/***********************************************************************
/
/  PARTICLE BUFFER HANDLER ROUTINES
/
/  written by: Matt Turk, John Wise
/  date:       May, 2011
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/
#ifdef USE_MPI
#include "communicators.h"
#endif

#include "preincludes.h"

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
#include "fortran.def"
#include "CosmologyParameters.h"

#include "ActiveParticle.h"


ParticleBufferHandler::ParticleBufferHandler(void)
{
  this->NumberOfBuffers = 0;
  this->CalculateElementSize();
}

ParticleBufferHandler::ParticleBufferHandler(int NumberOfParticles)
{

  this->NumberOfBuffers = NumberOfParticles;
  this->AllocateMemory();
  this->CalculateElementSize();

}

ParticleBufferHandler::ParticleBufferHandler(ActiveParticleType **np,
					     int NumberOfParticles,
					     int type, int proc)
{

  int i, dim, index;

  this->NumberOfBuffers = 0;
  for (i = 0; i < NumberOfParticles; i++)
    if (np[i]->ReturnType() == type &&
	(np[i]->ReturnDestProcessor() == proc || proc==-1)) 
      this->NumberOfBuffers++;
  
  if (this->NumberOfBuffers > 0) {

    this->AllocateMemory();

    for (i = 0, index = 0; i < NumberOfParticles; i++) {
      if (np[i]->ReturnType() == type && 
	  (np[i]->ReturnDestProcessor() == proc || proc==-1)) {
	for (dim = 0; dim < MAX_DIMENSION; dim++) {
	  this->pos[dim][index] = np[i]->pos[dim];
	  this->vel[dim][index] = np[i]->vel[dim];
	}
	this->Mass[index] = np[i]->Mass;
	this->BirthTime[index] = np[i]->BirthTime;
	this->DynamicalTime[index] = np[i]->DynamicalTime;
	this->Metallicity[index] = np[i]->Metallicity;
	this->Identifier[index] = np[i]->Identifier;
	this->level[index] = np[i]->level;
	this->GridID[index] = np[i]->GridID;
	this->type[index] = np[i]->type;
	this->proc[index] = np[i]->dest_processor;
	index++;
      }
    }

  } // ENDIF NumberOfBuffers > 0

  this->CalculateElementSize();

}

ParticleBufferHandler::~ParticleBufferHandler(void)
{
  if (this->NumberOfBuffers > 0) {
    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
      delete[] pos[dim];
      delete[] vel[dim];
    }
    delete[] Mass;
    delete[] BirthTime;
    delete[] DynamicalTime;
    delete[] Metallicity;
    delete[] Identifier;
    delete[] level;
    delete[] GridID;
    delete[] type;
    delete[] proc;
  }
}

void ParticleBufferHandler::AllocateMemory(void)
{
  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    this->pos[dim] = new FLOAT[NumberOfBuffers];
    this->vel[dim] = new float[NumberOfBuffers];
  }
  this->Mass = new double[NumberOfBuffers];
  this->BirthTime = new float[NumberOfBuffers];
  this->DynamicalTime = new float[NumberOfBuffers];
  this->Metallicity = new float[NumberOfBuffers];
  this->Identifier = new PINT[NumberOfBuffers];
  this->level = new int[NumberOfBuffers];
  this->GridID = new int[NumberOfBuffers];
  this->type = new int[NumberOfBuffers];
  this->proc = new int[NumberOfBuffers];
  return;
}

static int ParticleBufferHandler::ReturnHeadersize(void)
{
  return sizeof(int);
  Eint32 mpi_flag = 0;
  this->ElementSizeInBytes = 0;
  this->HeaderSizeInBytes = 0;

#ifdef USE_MPI
  MPI_Initialized(&mpi_flag);
#endif

  Eint32 size;
  if (mpi_flag == 1) {
#ifdef USE_MPI
    // Header:
    // 1. Number of buffers (int)
    MPI_Pack_size(1, IntDataType, EnzoTopComm, &size);
    this->HeaderSizeInBytes += size;
#endif
  } else {
    this->HeaderSizeInBytes = 1*sizeof(int);
  }
  return;
}
