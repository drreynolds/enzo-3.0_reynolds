/*-*-C++-*-*/
/***********************************************************************
/
/  STAR PARTICLE STRUCTURE
/
/  written by: John Wise
/  date:       September, 2005
/  modified1:  John Wise
/  date:       March, 2009 (converted into a class)
/  modified2:  John Wise, Greg Bryan, Britton Smith, Cameron Hummels,
/              Matt Turk
/  date:       May, 2011 (converting from Star to ActiveParticle)
/
/  PURPOSE:
/
************************************************************************/

#ifndef __PARTICLE_ATTRIBUTE_HANDLER_H
#define __PARTICLE_ATTRIBUTE_HANDLER_H

/* http://www.gamedev.net/topic/474803-c-template-pointer-to-member/ */

class ParticleAttributeHandler
{

  public:

    std::string name;
    MPI_Datatype mpitype;
    Eint32 hdf5type;
    int offset;
};

template <class APClass, typename Type, Type APClass::*var>
class Handler : public ParticleAttributeHandler
{
  public:

    Handler(std::string name, int offset = 0) {
        this->name = name;
        this->offset = offset;

        /* Can't use a switch */
        if (typeid(Type) == typeid(int)) {
            this->mpitype = IntDataType;
        } else if (typeid(Type) == typeid(float)) {
            this->mpitype = FloatDataType;
        } else if (typeid(Type) == typeid(double)) {
            this->mpitype = MPI_DOUBLE;
        } else if (typeid(Type) == typeid(FLOAT)) {
            this->mpitype = FLOATDataType;
        } else {
            ENZO_FAIL("Unrecognized data type");
        }
    }

    void UnpackBuffer(char *mpi_buffer, int mpi_buffer_size,
                      int NumberOParticles, int *position,
                      ParticleBufferHandler *PBHInstance) {

        MPI_Unpack(mpi_buffer, mpi_buffer_size, &position,
                   PBHInstance->*var[this->offset],
                   PBHInstance->NumberOfBuffers,
                   mpitype, MPI_COMM_WORLD);
    }

    void *AllocateArray(int n) {
        return (void *) new Type[n];
    }

    void CopyToArray(int pos, void *output_, APClass instance) {
        Type *output = (Type *) output_;
        output[pos] = *(&(instance->*var) + this->offset);
    }

};

template <class APClass, typename Type, int N, Type (APClass::*var)[N]>
class ArrayHandler : public ParticleAttributeHandler
{
  public:

    ArrayHandler(std::string name, int offset = 0) {
        this->name = name;
        this->offset = offset;

        /* Can't use a switch */
        if (typeid(Type) == typeid(int)) {
            this->mpitype = IntDataType;
        } else if (typeid(Type) == typeid(float)) {
            this->mpitype = FloatDataType;
        } else if (typeid(Type) == typeid(double)) {
            this->mpitype = MPI_DOUBLE;
        } else if (typeid(Type) == typeid(FLOAT)) {
            this->mpitype = FLOATDataType;
        } else {
            ENZO_FAIL("Unrecognized data type");
        }
    }

    void UnpackBuffer(char *mpi_buffer, int mpi_buffer_size,
                      int NumberOParticles, int *position,
                      ParticleBufferHandler *PBHInstance) {

        MPI_Unpack(mpi_buffer, mpi_buffer_size, &position,
                   PBHInstance->*var[this->offset],
                   PBHInstance->NumberOfBuffers,
                   mpitype, MPI_COMM_WORLD);
    }

    void *AllocateArray(int n) {
        return (void *) new Type[n];
    }

    void CopyToArray(int pos, void *output_, APClass instance) {
        Type *output = (Type *) output_;
        output[pos] = *(&(instance->*var) + this->offset);
    }

};

typedef std::vector<ParticleAttributeHandler> AttributeVector ;

#endif
