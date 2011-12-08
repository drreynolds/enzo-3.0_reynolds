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
#include "mpi.h"
#endif

#include <map>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

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

void ParticleBufferHandler::CalculateElementSize(void)
{
  // Define the element size in bytes (excluding the processor number)
  // float: 6 -- vel[3], BirthTime, DynamicalTime, Metallicity
  // FLOAT: 3 -- pos[3]
  // double: 1 -- Mass
  // int: 3 -- level, GridID, type
  // PINT: 1 -- Identifier
#ifdef USE_MPI
  int size;
  this->ElementSizeInBytes = 0;
  MPI_Pack_size(6, FloatDataType, MPI_COMM_WORLD, &size);
  this->ElementSizeInBytes += size;
  MPI_Pack_size(3, MY_MPIFLOAT, MPI_COMM_WORLD, &size);
  this->ElementSizeInBytes += size;
  MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &size);
  this->ElementSizeInBytes += size;
  MPI_Pack_size(3, IntDataType, MPI_COMM_WORLD, &size);
  this->ElementSizeInBytes += size;
  MPI_Pack_size(1, PINTDataType, MPI_COMM_WORLD, &size);
  this->ElementSizeInBytes += size;
  
  // Header:
  // 1. Number of buffers (int)
  this->HeaderSizeInBytes = 0;
  MPI_Pack_size(1, IntDataType, MPI_COMM_WORLD, &size);
  this->HeaderSizeInBytes += size;
#else
  this->ElementSizeInBytes = 6*sizeof(float) + 3*sizeof(FLOAT) + 1*sizeof(double) +
    3*sizeof(int) + 1*sizeof(PINT);
  this->HeaderSizeInBytes = 1*sizeof(int);
#endif
  return;
}
