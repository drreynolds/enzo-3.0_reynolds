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

template <class APClass> int ParticleBufferHandler::ElementSize() {
    static int particle_size = 0;
    if (particle_size > 0) return particle_size;
    AttributeVector &handlers = APClass::ParticleAttributeHandlers;
    for(AttributeVector::iterator it = handlers.begin();
        it != handlers.end(); ++it) {
        particle_size += (*it).element_size;
    }
    return particle_size;
}

template <class APClass> void ParticleBufferHandler::Allocate(int Count, char **buffer) {
        
    /* This routine is called for each particle type. */
    /* So we need to re-calculate the element and header size for each. */

    int particle_size = this->ElementSize<APClass>();
    int header_size = APClass::ReturnHeaderSize();

    buffer = &(new char[particle_size * Count + header_size]);

}

template <class APClass> void ParticleBufferHandler::FillBuffer(
        ActiveParticleType **InList_, int InCount, char *buffer) {

    int i;

    if (buffer == NULL) {
        this->Allocate<APClass>(InCount, &buffer);
    }

    AttributeVector &handlers = APClass::ParticleAttributeHandlers;
    APClass **InList = dynamic_cast<APClass**>(InList_);

    for (i = 0; i < InCount; i++) {
        for(AttributeVector::iterator it = handlers.begin();
            it != handlers.end(); ++it) {
            it->GetAttribute(&buffer, InList[i]);
        }
    }
}

template <class APClass> void ParticleBufferHandler::Unpack(
        char *buffer, int offset,
        ActiveParticleType **OutList_, int OutCount) {

    APClass **OutList = dynamic_cast<APClass**>(OutList_);
    AttributeVector &handlers = APClass::ParticleAttributeHandlers;
    APClass *ap;
    int i;

    for (i = 0; i < OutCount; i++) {
        ap = new APClass();
        OutList[i + offset] = ap;
        for(AttributeVector::iterator it = handlers.begin();
            it != handlers.end(); ++it) {
            it->SetAttribute(&buffer, ap);
        }
    }

}

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

  Eint32 mpi_flag = 0;
  this->ElementSizeInBytes = 0;
  this->HeaderSizeInBytes = 0;

#ifdef USE_MPI
  MPI_Initialized(&mpi_flag);
#endif

  Eint32 size;
  if (mpi_flag == 1) {
#ifdef USE_MPI
    MPI_Pack_size(6, FloatDataType, EnzoTopComm, &size);
    this->ElementSizeInBytes += size;
    MPI_Pack_size(3, MY_MPIFLOAT, EnzoTopComm, &size);
    this->ElementSizeInBytes += size;
    MPI_Pack_size(1, MPI_DOUBLE, EnzoTopComm, &size);
    this->ElementSizeInBytes += size;
    MPI_Pack_size(3, IntDataType, EnzoTopComm, &size);
    this->ElementSizeInBytes += size;
    MPI_Pack_size(1, PINTDataType, EnzoTopComm, &size);
    this->ElementSizeInBytes += size;
  
    // Header:
    // 1. Number of buffers (int)
    MPI_Pack_size(1, IntDataType, EnzoTopComm, &size);
    this->HeaderSizeInBytes += size;
#endif
  } else {
    this->ElementSizeInBytes = 6*sizeof(float) + 3*sizeof(FLOAT) + 1*sizeof(double) +
      3*sizeof(int) + 1*sizeof(PINT);
    this->HeaderSizeInBytes = 1*sizeof(int);
  }
  return;
}
