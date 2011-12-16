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
#ifndef __ACTIVE_PARTICLE_H
#define __ACTIVE_PARTICLE_H

#include "hdf5.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "StarBuffer.h"

struct ActiveParticleFormationData;
struct ActiveParticleFormationDataFlags;
class ParticleBufferHandler;

class ActiveParticleType
{
public:
  void static ConstructData(grid *_grid,
			    ActiveParticleFormationDataFlags &flags,
			    ActiveParticleFormationData &data);
  void static DestroyData(grid *_grid,
			  ActiveParticleFormationData &data);
  int static WriteDataset(int ndims, hsize_t *dims, char *name, hid_t group,
			  hid_t data_type, void *data);
  int static ReadDataset(int ndims, hsize_t *dims, char *name, hid_t group,
			 hid_t data_type, void *read_to);

  /* Several pure virtual functions */
  
  /* This should return the number of new star particles created, and should
   * create them. */
//  ActiveParticleType(){};
//  ~ActiveParticleType(){};
//  ActiveParticleType(ActiveParticleType*& part){};

  ActiveParticleType(void);
  ActiveParticleType(grid *_grid, ActiveParticleFormationData &data);
  ActiveParticleType(grid *_grid, int _id, int _level);
  ActiveParticleType(ActiveParticleType*& part);
  ActiveParticleType(ParticleBufferHandler *buffer, int index);
  ~ActiveParticleType(void);

  void operator=(ActiveParticleType *a);

  template <class active_particle_class> active_particle_class *copy(void);

  int   ReturnID(void) { return Identifier; };
  double ReturnMass(void) { return Mass; };
  float ReturnBirthTime(void) { return BirthTime; };
  float ReturnDynamicalTime(void) { return DynamicalTime; };
  float ReturnMetallicity(void) { return Metallicity; };
  int   ReturnType(void) { return type; };
  int   ReturnLevel(void) { return level; };
  int   ReturnDestProcessor(void) { return dest_processor; };
  int   ReturnGridID(void) { return GridID; };
  grid *ReturnCurrentGrid(void) { return CurrentGrid; };

  void  ReduceLevel(void) { level--; };
  void  IncreaseLevel(void) { level++; };
  void  SetLevel(int i) { level = i; };
  void  SetGridID(int i) { GridID = i; };
  void  SetDestProcessor(int i) { dest_processor = i; };
  void  AssignCurrentGrid(grid *a) { this->CurrentGrid = a; };
  void  AddMass(double dM) { Mass += dM; };
  void  AdjustMassByFactor(double factor) { Mass *= factor; };

  FLOAT *ReturnPosition(void) { return pos; }
  float *ReturnVelocity(void) { return vel; }
  void   ConvertAllMassesToSolar(void);
  void   ConvertMassToSolar(void);
  void   Merge(ActiveParticleType *a);
  //bool  MergableMBH(ActiveParticleType *a);
  float Separation(ActiveParticleType *a);
  float Separation2(ActiveParticleType *a);
  float RelativeVelocity2(ActiveParticleType *a);
  void  UpdatePositionVelocity(void);
  void  MirrorToParticle(void);
  void  CopyFromParticle(grid *_grid, int _id, int _level);
  int   DisableParticle(LevelHierarchyEntry *LevelArray[]);
  int   SphereContained(LevelHierarchyEntry *LevelArray[], int level, 
			float Radius);
  void  PrintInfo(void);
  //void  ActivateNewStar(FLOAT Time, float Timestep);
  //void  DeleteCopyInGrid(void);
  //int   DeleteCopyInGridGlobal(LevelHierarchyEntry *LevelArray[]);

  /* Virtual and pure virtual functions in this base class */

  virtual bool IsARadiationSource(FLOAT Time) { return FALSE; };
  virtual bool Mergable(ActiveParticleType *a);
  virtual int GetEnabledParticleID(int id = -1) = 0;

#ifdef TRANSFER
  RadiationSourceEntry* RadiationSourceInitialize(void);
#endif

protected:
  grid *CurrentGrid;
  FLOAT	pos[MAX_DIMENSION];
  float vel[MAX_DIMENSION];
  double Mass;		// Msun
  float BirthTime;
  float DynamicalTime;      
  float Metallicity;
  PINT Identifier;
  int level;
  int GridID;
  int dest_processor;  // used for communication
  int type;
  
  bool Active;
  int EnabledParticleID;

private: /* Cannot be accessed by subclasses! */
  
  friend class grid;
  friend class ActiveParticleType_info;
  friend class ParticleBufferHandler;

};

/* Comparer functions for sorting particle buffers with std::sort */

struct cmp_ap_grid {
  bool operator()(ActiveParticleType* const& a, ActiveParticleType* const& b) const {
    if (a->ReturnGridID() < b->ReturnGridID()) return true;
    else return false;
  }
};

struct cmp_ap_proc {
  bool operator()(ActiveParticleType* const& a, ActiveParticleType* const& b) const {
    if (a->ReturnDestProcessor() < b->ReturnDestProcessor()) return true;
    else return false;
  }
};

struct cmp_ap_type {
  bool operator()(ActiveParticleType* const& a, ActiveParticleType* const& b) const {
    if (a->ReturnType() < b->ReturnType()) return true;
    else return false;
  }
};


struct ActiveParticleFormationData {
  int NumberOfNewParticles;
  int MaxNumberOfNewParticles;
  ActiveParticleType **NewParticles;
  /* This is where all the pointers that normally get passed into
   * formation routines gets placed. Things like fractional h2, dark
   * matter density, etc etc. Anything that's derived.  It's okay to
   * add to this.  */
  float *DarkMatterDensity;
  float *H2Fraction;
  float *CoolingTime;
  float *CoolingRate;
  float *Temperature;
  float *TotalMetals;
  float DensityUnits;
  float LengthUnits;
  float TemperatureUnits;
  float TimeUnits;
  float VelocityUnits;
  double MassUnits;
  int DensNum;
  int Vel1Num;
  int Vel2Num;
  int Vel3Num;
  int TENum;
  int GENum;
  int MetalNum;
  int MetalIaNum;
  int ColourNum;
  int level;
  int GridID;
};

const struct ActiveParticleFormationData data_default = {
  0,        // NumberOfNewParticles
  0,        // MaxNumberOfNewParticles
  NULL,     // NewParticles
  NULL,     // DarkMatterDensity
  NULL,     // H2Fraction
  NULL,     // CoolingTime
  NULL,     // CoolingRate
  NULL,     // Temperature
  NULL,     // TotalMetals
  0.0,      // DensityUnits
  0.0,      // LengthUnits
  0.0,      // TemperatureUnits
  0.0,      // TimeUnits
  0.0,      // VelocityUnits
  0.0,      //  MassUnits
  -1,       // DensNum
  -1,       // Vel1Num
  -1,       // Vel2Num
  -1,       // Vel3Num
  -1,       // TENum
  -1,       // GENum
  -1,       // MetalNum
  -1,       // MetalIaNum
  -1,       // ColourNum
  -1,       // level
  -1        // GridID
};


struct ActiveParticleFormationDataFlags {
  /* For every entry in the ActiveParticleFormationData struct, we
   * have a bool here. */
  bool DarkMatterDensity;
  bool H2Fraction;
  bool CoolingTime;
  bool CoolingRate;
  bool Temperature;
  bool UnitConversions;
  bool DataFieldNumbers;
  bool MetalField;
};

const struct ActiveParticleFormationDataFlags flags_default = {
  false,    // DarkMatterDensity
  false,    // H2Fraction
  false,    // CoolingTime
  false,    // CoolingRate
  false,    // Temperature
  false,    // UnitConversions
  false,    // DataFieldNumbers
  false     // MetalField
};


//! maps the name of a plug-in to a pointer of the factory pattern
class ActiveParticleType_info;
typedef std::map<std::string, ActiveParticleType_info *> ActiveParticleMap;

ActiveParticleMap &get_active_particle_types();

void EnableActiveParticleType(char *active_particle_type_name);

class ParticleBufferHandler
{
public:
  ParticleBufferHandler(void);
  ParticleBufferHandler(int NumberOfParticles);
  ParticleBufferHandler(ActiveParticleType **np, int NumberOfParticles, int type, int proc);
  ~ParticleBufferHandler();
  /*virtual void WriteBuffers(hid_t group);*/
  int _AllocateBuffer(char *buffer, int &buffer_size, int &position); // helper function for derived classes
  void CalculateElementSize(void);
  void AllocateMemory(void);
  int ReturnHeaderSize(void) { return HeaderSizeInBytes; };
  int ReturnElementSize(void) { return ElementSizeInBytes; };
  int ReturnNumberOfBuffers(void) { return NumberOfBuffers; };
  int _UnpackBuffer(char *buffer, int buffer_size, int &position);

protected:
  int NumberOfBuffers;
  int HeaderSizeInBytes;
  int ElementSizeInBytes;
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

class ActiveParticleType_info
{
public:
       
  /* We will add more functions to this as necessary */
  ActiveParticleType_info
  (std::string this_name,
   int (*ffunc)(grid *thisgrid_orig, ActiveParticleFormationData &data),
   void (*dfunc)(ActiveParticleFormationDataFlags &flags),
   void (*abfunc)(ActiveParticleType **np, int NumberOfParticles, char *buffer, int &buffer_size,
		  int &position, int proc),
   void (*unfunc)(char *mpi_buffer, int mpi_buffer_size, int NumberOfParticles,
		  ActiveParticleType **np, int &npart),
   int (*ifunc)(),
   int (*feedfunc)(grid *thisgrid_orig, ActiveParticleFormationData &data),
   int (*writefunc)(ActiveParticleType *these_particles, int n, int GridRank, hid_t group_id),
   int (*readfunc)(ActiveParticleType **particles_to_read, int *n, int GridRank, hid_t group_id),
   int (*belfunc)(HierarchyEntry *Grids[], TopGridData *MetaData,
		  int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
		  int ThisLevel, int TotalStarParticleCountPrevious[],
		  int ActiveParticleID),
   int (*aelfunc)(HierarchyEntry *Grids[], TopGridData *MetaData,
		  int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
		  int ThisLevel, int TotalStarParticleCountPrevious[],
		  int ActiveParticleID),
   ActiveParticleType *particle,
   ParticleBufferHandler *buffer
   ){
    this->formation_function = ffunc;
    this->describe_data_flags = dfunc;
    this->allocate_buffer = abfunc;
    this->unpack_buffer = unfunc;
    this->particle_instance = particle;
    this->buffer_instance = buffer;
    this->initialize = ifunc;
    this->feedback_function = feedfunc;
    this->write_function = writefunc;
    this->read_function = readfunc;
    this->before_evolvelevel_function = belfunc;
    this->after_evolvelevel_function = aelfunc;
    get_active_particle_types()[this_name] = this;
  }

  static int count(){return get_active_particle_types().size();}
  int GetEnabledParticleID(){return this->MyEnabledParticleID;}

  int Enable(){
    /* 0-indexed */
    this->MyEnabledParticleID = this->TotalEnabledParticleCount++;
    this->particle_instance->GetEnabledParticleID(this->MyEnabledParticleID);
    return this->MyEnabledParticleID;
  }

  int (*initialize)(void);
  int (*formation_function)(grid *thisgrid_orig, ActiveParticleFormationData &data);
  int (*feedback_function)(grid *thisgrid_orig, ActiveParticleFormationData &data);
  int (*write_function)(ActiveParticleType *these_particles, int n, int GridRank, hid_t group_id);
  int (*read_function)(ActiveParticleType **particles_to_read, int *n, int GridRank, hid_t group_id);
  int (*before_evolvelevel_function)(HierarchyEntry *Grids[], TopGridData *MetaData,
				     int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				     int ThisLevel, int TotalStarParticleCountPrevious[],
				     int ActiveParticleID);
  int (*after_evolvelevel_function)(HierarchyEntry *Grids[], TopGridData *MetaData,
				    int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				    int ThisLevel, int TotalStarParticleCountPrevious[],
				    int ActiveParticleID);
  void (*describe_data_flags)(ActiveParticleFormationDataFlags &flags);
  void (*allocate_buffer)(ActiveParticleType **np, int NumberOfParticles, char *buffer, int &buffer_size,
			  int &position, int proc);
  void (*unpack_buffer)(char *mpi_buffer, int mpi_buffer_size, int NumberOfParticles, 
			ActiveParticleType **np, int &npart);
  ActiveParticleType* particle_instance;
  ParticleBufferHandler* buffer_instance;

private:
  /* This is distinct from the global as a redundant error-checking
     pattern */
  static int TotalEnabledParticleCount;
  int MyEnabledParticleID; /* Defaults to 0 */
  int *EnabledParticleIDPointer;

};

template <class active_particle_class, class particle_buffer_handler>
ActiveParticleType_info *register_ptype(std::string name)
{
  active_particle_class *pp = new active_particle_class();
  particle_buffer_handler *buffer = new particle_buffer_handler();
  ActiveParticleType_info *pinfo = new ActiveParticleType_info
    (name,
     (&active_particle_class::EvaluateFormation),
     (&active_particle_class::DescribeSupplementalData),
     (&particle_buffer_handler::AllocateBuffer),
     (&particle_buffer_handler::UnpackBuffer),
     (&active_particle_class::InitializeParticleType),
     (&active_particle_class::EvaluateFeedback),
     (&active_particle_class::WriteToOutput),
     (&active_particle_class::ReadFromOutput),
     (&active_particle_class::BeforeEvolveLevel),
     (&active_particle_class::AfterEvolveLevel),
     pp,
     buffer);
  return pinfo;
}

#define ENABLED_PARTICLE_ID_ACCESSOR					\
  int GetEnabledParticleID(int myid = -1) {				\
    static int ParticleID = -1;						\
    if (myid >= 0) {							\
      if (ParticleID != -1) ENZO_FAIL("Setting Particle ID Twice!");	\
      ParticleID = myid;						\
    }									\
    return ParticleID;							\
  };


#endif

