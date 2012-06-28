/*-*-C++-*-*/
/***********************************************************************
/
/  PARTICLE BUFFER HANDLER
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
#ifndef __PARTICLE_BUFFER_HANDLER_H
#define __PARTICLE_BUFFER_HANDLER_H

class ActiveParticleType;

class ParticleBufferHandler
{
public:
  ParticleBufferHandler(void);
  ParticleBufferHandler(int NumberOfParticles);
  ParticleBufferHandler(ActiveParticleType **np, int NumberOfParticles, int type, int proc);
  ~ParticleBufferHandler();
  /*virtual void WriteBuffers(hid_t group);*/
  int _AllocateBuffer(char *buffer, Eint32 total_buffer_size, int &buffer_size, 
		      Eint32 &position); // helper function for derived classes
  void CalculateElementSize(void);
  void AllocateMemory(void);
  int ReturnNumberOfBuffers(void) { return NumberOfBuffers; };
  int _UnpackBuffer(char *buffer, int buffer_size, Eint32 &position);

  int NumberOfBuffers;

  template <class APClass> int ElementSize();

  template <class APClass> int Allocate(int Count, char **buffer);

  template <class APClass> void Unpack(
          char *buffer, int offset,
          ActiveParticleType **OutList_, int OutCount);


  template <class APClass> void FillBuffer(
          ActiveParticleType **InList_, int InCount, char *buffer);

protected:
  static int HeaderSizeInBytes;
  static int ElementSizeInBytes;
  FLOAT	*pos[MAX_DIMENSION];
  float *vel[MAX_DIMENSION];
  double *Mass;		// Msun
  float *BirthTime;
  float *DynamicalTime;      
  float *Metallicity;
  PINT *Identifier;
  int *level;
  int *GridID;
  int *type;
  int *proc;

  friend class ActiveParticleType;

};

#endif
