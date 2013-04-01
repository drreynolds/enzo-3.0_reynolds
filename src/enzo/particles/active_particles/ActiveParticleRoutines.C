/***********************************************************************
/
/  ROUTINES FOR THE STAR PARTICLE CLASS
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:  November, 2011 (converting into active particles)
/
/  PURPOSE: Instead of restricting star particles to the typical 
/           particle attributes, this class gives more functionality 
/           to them.
/
************************************************************************/
#include "preincludes.h"
#include "performance.h"
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
#include "LevelHierarchy.h"
#include "ActiveParticle.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

/*******************************

   CONSTRUCTORS AND DESTRUCTOR

 *******************************/

ActiveParticleType::ActiveParticleType(void)
{
  int dim;
  CurrentGrid = NULL;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    pos[dim] = vel[dim] = 0.0;
  Mass = BirthTime = DynamicalTime = 0.0;
  level = GridID = type = 0;


  /* The correct indices are assigned in CommunicationUpdateActiveParticleCount 
     in ActiveParticleFinalize.*/
  Identifier = INT_UNDEFINED;
}

ActiveParticleType::ActiveParticleType(ActiveParticleType* part)
{
  int dim;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim] = part->pos[dim];
    vel[dim] = part->vel[dim];
  }
  Mass = part->Mass;
  BirthTime = part->BirthTime;
  DynamicalTime = part->DynamicalTime;
  Metallicity = part->Metallicity;
  Identifier = part->Identifier;
  level = part->level;
  GridID = part->GridID;
  type = part->type;
  CurrentGrid = part->CurrentGrid;
  dest_processor = part->dest_processor;
}

ActiveParticleType::ActiveParticleType(grid *_grid, ActiveParticleFormationData &data)
{
  int dim;
  CurrentGrid = _grid;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    pos[dim] = vel[dim] = 0.0;
  Mass = BirthTime = DynamicalTime = 0.0;
  type = 0;
  level = data.level;
  GridID = data.GridID;

  /* The correct indices are assigned in CommunicationUpdateActiveParticleCount 
     in ActiveParticleFinalize.*/
  Identifier = INT_UNDEFINED;
  dest_processor = -1;
}


ActiveParticleType::ActiveParticleType(grid *_grid, int _id, int _level)
{

  assert(_id < _grid->NumberOfParticles);

  int dim;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim] = _grid->ParticlePosition[dim][_id];
    vel[dim] = _grid->ParticleVelocity[dim][_id];
  }
  CurrentGrid = _grid;
  level = _level;

  GridID = _grid->ID;
  type = _grid->ParticleType[_id];
  Identifier = _grid->ParticleNumber[_id];
  Mass = (double)(_grid->ParticleMass[_id]);


  // No more attributes.  Everything stored in active particles.
//  BirthTime = _grid->ParticleAttribute[0][_id];
//  DynamicalTime = _grid->ParticleAttribute[1][_id];
//  Metallicity = _grid->ParticleAttribute[2][_id];
  this->ConvertMassToSolar();
}

/* No need to delete the accretion arrays because the pointers are
   stored in the copies located in the grid class. */

ActiveParticleType::~ActiveParticleType(void)
{
}

/***************

    OPERATORS

 ***************/

void ActiveParticleType::operator=(ActiveParticleType *a)
{
  int i, dim;
  CurrentGrid = a->CurrentGrid;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim] = a->pos[dim];
    vel[dim] = a->vel[dim];
  }
  Mass = a->Mass;
  BirthTime = a->BirthTime;
  DynamicalTime = a->DynamicalTime;
  Metallicity = a->Metallicity;
  Identifier = a->Identifier;
  level = a->level;
  GridID = a->GridID;
  type = a->type;
  dest_processor = -1;
  return;
}

/**********************

   CONVENIENT ROUTINES

 **********************/

template<class active_particle_class>
active_particle_class *ActiveParticleType::copy(void)
{
  int i, dim;
  active_particle_class *a = new active_particle_class();
  a->CurrentGrid = CurrentGrid;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    a->pos[dim] = pos[dim];
    a->vel[dim] = vel[dim];
  }
  a->Mass = Mass;
  a->BirthTime = BirthTime;
  a->DynamicalTime = DynamicalTime;
  a->Metallicity = Metallicity;
  a->Identifier = Identifier;
  a->level = level;
  a->GridID = GridID;
  a->type = type;
  a->dest_processor = -1;
  return a;
}

void ActiveParticleType::ConvertAllMassesToSolar(void)
{
  const double Msun = 1.989e33;
  double dx;
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits, MassConversion;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, CurrentGrid->Time);
  dx = LengthUnits * CurrentGrid->CellWidth[0][0];
  MassConversion = (float) (dx*dx*dx * double(DensityUnits) / Msun);
  this->Mass *= MassConversion;
  return;
}

void ActiveParticleType::ConvertMassToSolar(void)
{
  const double Msun = 1.989e33;
  double dx;
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits, MassConversion;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, CurrentGrid->Time);
  dx = LengthUnits * CurrentGrid->CellWidth[0][0];
  MassConversion = (float) (dx*dx*dx * double(DensityUnits) / Msun);
  this->Mass *= MassConversion;
  return;
}

void  ActiveParticleType::AdjustVelocity(float VelocityIncrement[])
{ 
  int i;
  for (i = 0; i<3; i++)
    vel[i] += VelocityIncrement[i];
  return;
}

void ActiveParticleType::SetVelocity(float NewVelocity[])
{
  int i;
  for (i = 0; i<3; i++)
    vel[i] = NewVelocity[i];
  return;
}

void ActiveParticleType::SetPosition(FLOAT NewPosition[])
{
  int i;
  for (i = 0; i<3; i++)
    pos[i] = NewPosition[i];
  return;
}

void ActiveParticleType::SetPositionPeriod(FLOAT period[])
{
  int i;
  for (i = 0; i<3; i++) {
    pos[i] = fmod(pos[i], period[i]);
    if (pos[i] < 0) {
      pos[i] += period[i];
    }
  }
  return;
}

void ActiveParticleType::Merge(ActiveParticleType *a)
{
  int dim;
  double ratio1, ratio2;
  ratio1 = Mass / (Mass + a->Mass);
  ratio2 = 1.0 - ratio1;
  Metallicity = ratio1 * Metallicity + ratio2 * a->Metallicity;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim] = ratio1 * pos[dim] + ratio2 * a->pos[dim];
    vel[dim] = ratio1 * vel[dim] + ratio2 * a->vel[dim];
    //accreted_angmom[dim] = ratio1 * accreted_angmom[dim] + ratio2 * a.accreted_angmom[dim];
  }
  Mass += a->Mass;
  //FinalMass += a.FinalMass;
  //DeltaMass += a.DeltaMass;
  //last_accretion_rate += a.last_accretion_rate;
  //NotEjectedMass += a.NotEjectedMass;
  return;
}

bool ActiveParticleType::Mergable(ActiveParticleType *a)
{
  // Only merge yet-to-be born stars
  return type == a->type && type < 0;
}

#ifdef UNUSED
bool ActiveParticleType::MergableMBH(ActiveParticleType *a)
{
  // Merge MBH particle with another 
  return type == a->type && type == MBH;
}
#endif

float ActiveParticleType::Separation2(ActiveParticleType *a)
{
  int dim;
  float dr, result = 0;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    dr = pos[dim] - a->pos[dim];
    result += dr*dr;
  }
  return result;
}

float ActiveParticleType::Separation(ActiveParticleType *a)  { return sqrt(this->Separation2(a)); }

float ActiveParticleType::RelativeVelocity2(ActiveParticleType *a)
{
  int dim;
  float dv, result = 0;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    dv = vel[dim] - a->vel[dim];
    result += dv*dv;
  }
  return result;
}

void ActiveParticleType::UpdatePositionVelocity(void)
{
  LCAPERF_START("star_UpdatePositionVelocity");
  int i, dim;
  int _id = -1;
  if (CurrentGrid != NULL && type >= 0) { // on local processor and active
    // Search for particle
    for (i = 0; i < CurrentGrid->NumberOfParticles; i++)
      if (Identifier == CurrentGrid->ParticleNumber[i]) {
	_id = i;
	break;
      }
    assert(_id >= 0);
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      pos[dim] = CurrentGrid->ParticlePosition[dim][_id];
      vel[dim] = CurrentGrid->ParticleVelocity[dim][_id];
    }
  }
  LCAPERF_STOP("star_UpdatePositionVelocity");
  return;
}

void ActiveParticleType::CopyFromParticle(grid *_grid, int _id, int _level)
{
  int dim;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim] = _grid->ParticlePosition[dim][_id];
    vel[dim] = _grid->ParticleVelocity[dim][_id];
  }
  CurrentGrid = _grid;
  level = _level;
  GridID = _grid->ID;

  // No more attributes.  Everything stored in active particles.
//  BirthTime = _grid->ParticleAttribute[0][_id];
//  DynamicalTime = _grid->ParticleAttribute[1][_id];
//  Metallicity = _grid->ParticleAttribute[2][_id];

  // below is removed because we want to keep Star->Mass as double 
  // during the run - Ji-hoon Kim, Dec.2009
//  Mass = (double)(_grid->ParticleMass[_id]); 
//  this->ConvertMassToSolar();
  return;
}

void ActiveParticleType::PrintInfo(void)
{
  printf("[P%d] ActiveParticle %"ISYM": pos = %"PSYM" %"PSYM" %"PSYM", vel = %"FSYM" %"FSYM" %"FSYM"\n",
	 MyProcessorNumber, Identifier, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]);
  printf("\t birthtime = %"FSYM", tdyn = %"FSYM"\n", BirthTime, DynamicalTime);
  printf("\t Z = %"GSYM"\n", Metallicity);
  printf("\t mass = %"GSYM", type = %"ISYM", grid %"ISYM", lvl %"ISYM"\n", 
	 Mass, type, GridID, level);
  return;
}

#ifdef TRANSFER
RadiationSourceEntry* ActiveParticleType::RadiationSourceInitialize(void)
{
  RadiationSourceEntry *source = new RadiationSourceEntry;
  source->PreviousSource = GlobalRadiationSources;
  source->NextSource     = GlobalRadiationSources->NextSource;
  source->SuperSource    = NULL;  // Define this later (below)
  source->GridID         = GridID;
  source->GridLevel      = level;
  source->Type           = type;
  source->LifeTime       = DynamicalTime;
  source->CreationTime   = BirthTime;
  source->AddedEmissivity = false;
  source->Position       = new FLOAT[3];
  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    if (pos[dim] < DomainLeftEdge[dim])
      source->Position[dim] = pos[dim] + DomainRightEdge[dim] - DomainLeftEdge[dim];
    else if (pos[dim] >= DomainRightEdge[dim])
      source->Position[dim] = pos[dim] - DomainRightEdge[dim] + DomainLeftEdge[dim];
    else
      source->Position[dim] = pos[dim];
  }
  return source;
}
#endif

