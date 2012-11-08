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

#include "communicators.h"

class ActiveParticleType;

class ParticleAttributeHandler
{

  public:

    std::string name;
    MPI_Datatype mpitype;
    Eint32 hdf5type;
    int element_size;
    int offset;

    virtual void SetAttribute(char **buffer, ActiveParticleType *pp) = 0;

    virtual int GetAttribute(char **buffer, ActiveParticleType *pp)  = 0;

    virtual void PrintAttribute(ActiveParticleType *pp)  = 0;

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
            this->hdf5type = HDF5_INT;
        } else if (typeid(Type) == typeid(float)) {
            this->mpitype = FloatDataType;
            this->hdf5type = HDF5_REAL;
        } else if (typeid(Type) == typeid(double)) {
            this->mpitype = MPI_DOUBLE;
            this->hdf5type = HDF5_R8;
        } else if (typeid(Type) == typeid(FLOAT)) {
            this->mpitype = FLOATDataType;
            this->hdf5type = HDF5_PREC;
        } else {
            ENZO_FAIL("Unrecognized data type");
        }
        this->element_size = sizeof(Type);
    }

    void SetAttribute(char **buffer, ActiveParticleType *pp_) {
        APClass *pp = static_cast<APClass*>(pp_);
        Type *pb = (Type *)(*buffer);
        pp->*var = *(pb++);
        *buffer = (char *) pb;
    }

    int GetAttribute(char **buffer, ActiveParticleType *pp_) {
        APClass *pp = static_cast<APClass*>(pp_);
        Type *pb = (Type *)(*buffer);
        *(pb++) = pp->*var;
        *buffer = (char *) pb;
        return this->element_size;
    }

    void PrintAttribute(ActiveParticleType *pp_) {
        APClass *pp = static_cast<APClass*>(pp_);
        std::cout << std::setprecision(15) << this->name << ": " << pp->*var;
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
            this->hdf5type = HDF5_INT;
        } else if (typeid(Type) == typeid(float)) {
            this->mpitype = FloatDataType;
            this->hdf5type = HDF5_REAL;
        } else if (typeid(Type) == typeid(double)) {
            this->mpitype = MPI_DOUBLE;
            this->hdf5type = HDF5_R8;
        } else if (typeid(Type) == typeid(FLOAT)) {
            this->mpitype = FLOATDataType;
            this->hdf5type = HDF5_PREC;
        } else {
            ENZO_FAIL("Unrecognized data type");
        }
        this->element_size = sizeof(Type);
    }

    void SetAttribute(char **buffer, ActiveParticleType *pp_) {
        APClass *pp = static_cast<APClass*>(pp_);
        Type *pb = (Type *)(*buffer);
        (pp->*var)[this->offset] = *(pb++);
        *buffer = (char *) pb;
    }

    int GetAttribute(char **buffer, ActiveParticleType *pp_) {
        APClass *pp = static_cast<APClass*>(pp_);
        Type *pb = (Type *)(*buffer);
        *(pb++) = (pp->*var)[this->offset];
        *buffer = (char *) pb;
        return this->element_size;
    }

    void PrintAttribute(ActiveParticleType *pp_) {
        APClass *pp = static_cast<APClass*>(pp_);
        std::cout << this->name << ": " << (pp->*var)[this->offset];
    }

};

typedef std::vector<ParticleAttributeHandler*> AttributeVector ;

#endif
