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

#include <typeinfo>
#include "TopGridData.h"
#include "ParticleAttributeHandler.h"
#include "h5utilities.h"

struct ActiveParticleFormationData;
struct ActiveParticleFormationDataFlags;

class ActiveParticleType
{
public:
  void static ConstructData(grid *_grid,
			    ActiveParticleFormationDataFlags &flags,
			    ActiveParticleFormationData &data);
  void static DestroyData(grid *_grid,
			  ActiveParticleFormationData &data);
  int static WriteDataset(int ndims, hsize_t *dims, const char *name, hid_t group,
			  hid_t data_type, void *data);
  int static ReadDataset(int ndims, hsize_t *dims, const char *name, hid_t group,
			 hid_t data_type, void *read_to);
  void static SetupBaseParticleAttributes(
    std::vector<ParticleAttributeHandler*> &handlers);

  void OutputPositionInformation(void);

  /* Several pure virtual functions */

  /* This should return the number of new star particles created, and should
   * create them. */

  ActiveParticleType(void);
  ActiveParticleType(grid *_grid, ActiveParticleFormationData &data);
  ActiveParticleType(grid *_grid, int _id, int _level);
  ActiveParticleType(ActiveParticleType* part);
  ~ActiveParticleType(void);

  void operator=(ActiveParticleType *a);

  template <class active_particle_class> active_particle_class *copy(void);

  PINT   ReturnID(void) { return Identifier; };
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
  void  ReduceLevel(int x) { level -= x; };
  void  IncreaseLevel(void) { level++; };
  void  IncreaseLevel(int x) { level += x; };
  void  SetLevel(int i) { level = i; };
  void  SetGridID(int i) { GridID = i; };
  void  SetDestProcessor(int i) { dest_processor = i; };
  void  AssignCurrentGrid(grid *a) { this->CurrentGrid = a; };
  void  AddMass(double dM) { Mass += dM; };
  void  AdjustMassByFactor(double factor) { Mass *= factor; };
  void  AdjustVelocity(float VelocityIncrement[]);
  void  SetVelocity(float NewVelocity[]);
  void  SetPosition(FLOAT NewPosition[]);
  void  SetPositionPeriod(FLOAT period[]);

  FLOAT *ReturnPosition(void) { return pos; };
  float *ReturnVelocity(void) { return vel; };
  float ReturnMomentum(int dim) { return Mass*vel[dim]; };
  void   ConvertAllMassesToSolar(void);
  void   ConvertMassToSolar(void);
  void   Merge(ActiveParticleType *a);
  float Separation(ActiveParticleType *a);
  float Separation2(ActiveParticleType *a);
  float RelativeVelocity2(ActiveParticleType *a);
  void  UpdatePositionVelocity(void);
  void  MirrorToParticle(void);
  void  CopyFromParticle(grid *_grid, int _id, int _level);
  int   DisableParticle(LevelHierarchyEntry *LevelArray[], int NewProcessorNumber);
  int   SphereContained(LevelHierarchyEntry *LevelArray[], int level,
			float Radius);
  void  PrintInfo(void);
  //void  ActivateNewStar(FLOAT Time, float Timestep);
  //void  DeleteCopyInGrid(void);
  //int   DeleteCopyInGridGlobal(LevelHierarchyEntry *LevelArray[]);

  /* Virtual and pure virtual functions in this base class */

  virtual bool IsARadiationSource(FLOAT Time) { return false; };
  virtual bool Mergable(ActiveParticleType *a);
  virtual int GetEnabledParticleID(int id = -1) {
    ENZO_FAIL("Not implemented.");
  };

#ifdef TRANSFER
  RadiationSourceEntry* RadiationSourceInitialize(void);
#endif

protected:
  grid *CurrentGrid;
  FLOAT	pos[MAX_DIMENSION];
  float vel[MAX_DIMENSION];
  double Mass;
  float BirthTime;
  float DynamicalTime;
  float Metallicity;
  PINT Identifier;
  int level;
  int GridID;
  int dest_processor;  // used for communication
  int type;

private: /* Cannot be accessed by subclasses! */

  friend class grid;
  friend class ActiveParticleType_info;

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


struct cmp_ap_number {
  bool operator()(ActiveParticleType* const& a, ActiveParticleType* const& b) const {
    if (a->ReturnID() < b->ReturnID()) return true;
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
  FLOAT CellSize;
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
  bool MetalField;
};

const struct ActiveParticleFormationDataFlags flags_default = {
  false,    // DarkMatterDensity
  false,    // H2Fraction
  false,    // CoolingTime
  false,    // CoolingRate
  false,    // Temperature
  false     // MetalField
};

namespace ActiveParticleHelpers {

  template <class APClass> int CalculateElementSize() {
      static int particle_size = 0;
      if (particle_size > 0) return particle_size;
      AttributeVector &handlers = APClass::AttributeHandlers;
      for(AttributeVector::iterator it = handlers.begin();
          it != handlers.end(); ++it) {
          particle_size += (**it).element_size;
      }
      return particle_size;
  }

  template <class APClass> void Allocate(int Count, char **buffer) {

      /* This routine is called for each particle type. */
      /* So we need to re-calculate the element and header size for each. */

      int particle_size = CalculateElementSize<APClass>();
      int header_size = sizeof(int);

      *buffer = new char[particle_size * Count + header_size];

  }

  template <class APClass> void PrintActiveParticle(APClass *ap,
      const std::string prefix = "") {

      AttributeVector &handlers = APClass::AttributeHandlers;
      for(AttributeVector::iterator it = handlers.begin();
        it != handlers.end(); ++it) {
          (*it)->PrintAttribute(ap);
          std::cout << " ";
      }

      std::cout << std::endl;

  }

  template <class APClass> void WriteParticles(
          ActiveParticleType **InList, int ParticleTypeID,
          int TotalParticles,
          const std::string name, hid_t grid_node)
  {
      int i, size = 0, Count = 0;
      char *buffer, *_buffer;
      AttributeVector &handlers = APClass::AttributeHandlers;
      const char *_name = name.c_str();
      hid_t node = H5Gcreate(grid_node, _name, 0);

      /* Now we count up our particles */
      for (i = 0; i < TotalParticles; i++) {
        if (InList[i]->GetEnabledParticleID() == ParticleTypeID) Count++;
      }

      writeScalarAttribute(node, HDF5_INT, "Count", &Count);

      if(Count == 0) {
        H5Gclose(node);
        return;
      }

      int ndims = 1;
      hsize_t dims[1] = {Count};

      APClass *In;

      for (AttributeVector::iterator it = handlers.begin();
          it != handlers.end(); ++it) {
          size = Count * (*it)->element_size;
          _buffer = buffer = new char[size];
          for (i = 0; i < TotalParticles; i++) {
            if (InList[i]->GetEnabledParticleID() != ParticleTypeID)
              continue;
            (*it)->GetAttribute(&_buffer, InList[i]);
          }
          /* Now write it to disk */
          _name = (*it)->name.c_str();
          APClass::WriteDataset(ndims, dims, _name, node,
                                (*it)->hdf5type, buffer);
          delete [] buffer;
      }

      H5Gclose(node);
  }

  template <class APClass> int ReadParticles(
          ActiveParticleType **OutList, int &offset, const std::string name,
          hid_t grid_node)
  {
      int i, size = 0;
      AttributeVector &handlers = APClass::AttributeHandlers;
      const char *_name = name.c_str();
      hid_t node = H5Gopen(grid_node, _name);
      /* We now have a reference to the top level node. */
      int Count;
      readAttribute(node, HDF5_INT, "Count", &Count, false);
      if(Count == 0) {
        H5Gclose(node);
        return offset;
      }

      char *buffer, *_buffer;
      int ndims = 1;
      hsize_t dims[1] = {Count};
      for (i = 0; i < Count; i++) {
          OutList[i+offset] = new APClass();
      }

      for (AttributeVector::iterator it = handlers.begin();
          it != handlers.end(); ++it) {
          size = Count * (*it)->element_size;
          _buffer = buffer = new char[size];
          _name = (*it)->name.c_str();
          APClass::ReadDataset(ndims, dims, _name, node,
                               (*it)->hdf5type, buffer);
          for (i = 0; i < Count; i++) {
              (*it)->SetAttribute(&_buffer, OutList[i+offset]);
          }
          delete [] buffer;
      }

      H5Gclose(node);
      offset += Count;
      return Count;

  }



  template <class APClass> int FillBuffer(
          ActiveParticleType **InList_, int InCount, char *buffer_) {

      int i;
      int size = 0;

      if (buffer_ == NULL) {
          ENZO_FAIL("Buffer not allocated!");
      }
      /* We increment the pointer as we fill */
      char **buffer = &buffer_;

      AttributeVector &handlers = APClass::AttributeHandlers;

      APClass *In;

      for (i = 0; i < InCount; i++) {
          In = dynamic_cast<APClass*>(InList_[i]);
          /* We'll put debugging output here */
          for(AttributeVector::iterator it = handlers.begin();
              it != handlers.end(); ++it) {
              size += (*it)->GetAttribute(buffer, In);
          }
          /*
          std::cout << "APF[" << MyProcessorNumber << "] " << i << " " << size << " ";
          PrintActiveParticle<APClass>(In);
           */
      }
      return size;
  }

  template <class APClass> void Unpack(
          char *buffer_, int offset,
          ActiveParticleType** OutList_, int OutCount) {

      APClass **OutList = reinterpret_cast<APClass**>(OutList_);
      AttributeVector &handlers = APClass::AttributeHandlers;
      APClass *Out;
      int i;
      char *buffer = buffer_;

      for (i = 0; i < OutCount; i++) {
          Out = new APClass();
          OutList[i + offset] = Out;
          for(AttributeVector::iterator it = handlers.begin();
              it != handlers.end(); ++it) {
              (*it)->SetAttribute(&buffer, Out);
          }
          /*
          std::cout << "APU[" << MyProcessorNumber << "] " << i << " ";
          PrintActiveParticle<APClass>(Out);
          */
      }

  }

}

//! maps the name of a plug-in to a pointer of the factory pattern
class ActiveParticleType_info;
typedef std::map<std::string, ActiveParticleType_info *> ActiveParticleMap;

ActiveParticleMap &get_active_particle_types();

void EnableActiveParticleType(char *active_particle_type_name);

ActiveParticleType** ActiveParticleFindAll(LevelHierarchyEntry *LevelArray[], int *GlobalNumberOfActiveParticles,
					   int ActiveParticleIDToFind);

class ActiveParticleType_info
{
public:

  /* We will add more functions to this as necessary */
  ActiveParticleType_info
  (std::string this_name,
   /* These functions hang off the ActiveParticle subclass */
   int (*evaluate_formation)(grid *thisgrid_orig, ActiveParticleFormationData &data),
   void (*describe_data)(ActiveParticleFormationDataFlags &flags),
   int (*initialize)(),
   int (*feedback)(grid *thisgrid_orig, ActiveParticleFormationData &data),
   int (*before_evolvelevel)(HierarchyEntry *Grids[], TopGridData *MetaData,
		  int NumberOfGrids, LevelHierarchyEntry *LevelArray[],
		  int ThisLevel, int TotalStarParticleCountPrevious[],
		  int ActiveParticleID),
   int (*after_evolvelevel)(HierarchyEntry *Grids[], TopGridData *MetaData,
		  int NumberOfGrids, LevelHierarchyEntry *LevelArray[],
		  int ThisLevel, int TotalStarParticleCountPrevious[],
		  int ActiveParticleID),
   int (*deposit_mass)(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[],
				int ThisLevel, int GalaxyParticleID),
   int (*flagfield)(LevelHierarchyEntry *LevelArray[], int level, int TopGridDims[], int ActiveParticleID),
   void (*allocate_buffer)(int Count, char **buffer),
   int (*fill_buffer)(ActiveParticleType **InList_, int InCount, char *buffer),
   void (*unpack_buffer)(char *buffer, int offset, ActiveParticleType** Outlist,
                       int OutCount),
   int (*element_size)(void),
   void (*write_particles)(ActiveParticleType **particles,
                         int type_id, int total_particles,
                         const std::string name, hid_t node),
   int (*read_particles)(ActiveParticleType **particles, int &offset, const
                         std::string name, hid_t node),
   ActiveParticleType *particle)
   {

    this->InitializeParticleType = initialize;
    this->EvaluateFormation = evaluate_formation;
    this->EvaluateFeedback = feedback;
    this->BeforeEvolveLevel = before_evolvelevel;
    this->AfterEvolveLevel = after_evolvelevel;
    this->DepositMass = deposit_mass;
    this->SetFlaggingField = flagfield;
    this->DescribeSupplementalData = describe_data;
    this->FillBuffer = fill_buffer;
    this->AllocateBuffer = allocate_buffer;
    this->UnpackBuffer = unpack_buffer;
    this->ReturnElementSize = element_size;
    this->WriteParticles = write_particles;
    this->ReadParticles = read_particles;
    this->particle_instance = particle;
    this->particle_name = this_name;
    get_active_particle_types()[this_name] = this;
  }

  static int count(){return get_active_particle_types().size();}
  int GetEnabledParticleID(){return this->MyEnabledParticleID;}
  std::string GetEnabledParticleName(){return this->particle_name;}

  int Enable(){
    /* 0-indexed */
    this->MyEnabledParticleID = this->TotalEnabledParticleCount++;
    this->particle_instance->GetEnabledParticleID(this->MyEnabledParticleID);
    return this->MyEnabledParticleID;
  }

  int (*InitializeParticleType)(void);
  int (*EvaluateFormation)(grid *thisgrid_orig, ActiveParticleFormationData &data);
  int (*EvaluateFeedback)(grid *thisgrid_orig, ActiveParticleFormationData &data);
  int (*BeforeEvolveLevel)(HierarchyEntry *Grids[], TopGridData *MetaData,
				     int NumberOfGrids, LevelHierarchyEntry *LevelArray[],
				     int ThisLevel, int TotalStarParticleCountPrevious[],
				     int ActiveParticleID);
  int (*AfterEvolveLevel)(HierarchyEntry *Grids[], TopGridData *MetaData,
				    int NumberOfGrids, LevelHierarchyEntry *LevelArray[],
				    int ThisLevel, int TotalStarParticleCountPrevious[],
				    int ActiveParticleID);
  int (*DepositMass)(HierarchyEntry *Grids[], TopGridData *MetaData,
                    int NumberOfGrids, LevelHierarchyEntry *LevelArray[],
                    int ThisLevel, int ActiveParticleID);
  int (*SetFlaggingField)(LevelHierarchyEntry *LevelArray[], int level, int TopGridDims[], int ActiveParticleID);
  void (*DescribeSupplementalData)(ActiveParticleFormationDataFlags &flags);
  void (*AllocateBuffer)(int Count, char **buffer);
  void (*UnpackBuffer)(char *buffer, int offset, ActiveParticleType **Outlist,
                       int OutCount);
  int (*FillBuffer)(ActiveParticleType **InList, int InCount, char *buffer);
  int (*ReturnElementSize)(void);
  void (*WriteParticles)(ActiveParticleType **InList,
                       int ParticleTypeID, int TotalParticles,
                       const std::string name, hid_t node);
  int (*ReadParticles)(ActiveParticleType **OutList, int &offset, const
                       std::string name, hid_t node);
  ActiveParticleType* particle_instance;
  std::string particle_name;

  // At the moment the communication buffers don't contain a header.
  int ReturnHeaderSize(void) { return 0; }

private:
  /* This is distinct from the global as a redundant error-checking
     pattern */
  static int TotalEnabledParticleCount;
  int MyEnabledParticleID; /* Defaults to 0 */
  int *EnabledParticleIDPointer;

};

template <class APClass>
ActiveParticleType_info *register_ptype(std::string name)
{
  APClass *pp = new APClass();

  ActiveParticleType_info *pinfo = new ActiveParticleType_info
    (name,
     (&APClass::EvaluateFormation),
     (&APClass::DescribeSupplementalData),
     (&APClass::InitializeParticleType),
     (&APClass::EvaluateFeedback),
     (&APClass::template BeforeEvolveLevel<APClass>),
     (&APClass::template AfterEvolveLevel<APClass>),
     (&APClass::template DepositMass<APClass>),
     (&APClass::SetFlaggingField),
     (&ActiveParticleHelpers::Allocate<APClass>),
     (&ActiveParticleHelpers::FillBuffer<APClass>),
     (&ActiveParticleHelpers::Unpack<APClass>),
     (&ActiveParticleHelpers::CalculateElementSize<APClass>),
     (&ActiveParticleHelpers::WriteParticles<APClass>),
     (&ActiveParticleHelpers::ReadParticles<APClass>),
     pp);
  return pinfo;
}

#define ENABLED_PARTICLE_ID_ACCESSOR                                   \
  int GetEnabledParticleID(int myid = -1) {                            \
    static int ParticleID = -1;                                        \
    if (myid >= 0) {                                                   \
      if (ParticleID != -1) ENZO_FAIL("Setting Particle ID Twice!");   \
      ParticleID = myid;                                               \
    }                                                                  \
    return ParticleID;                                                 \
  };


#endif
