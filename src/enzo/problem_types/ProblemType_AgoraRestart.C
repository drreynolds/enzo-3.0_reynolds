/***********************************************************************
/
/  Agora isolated galaxy restart
/
/  written by: Nathan Goldbaum
/  date:       March, 2013
/
/  PURPOSE:
/  https://sites.google.com/site/projectagoraworkspace/metagroup1/group2
/  http://kicp.uchicago.edu/~agertz/disk/readme.txt
/
************************************************************************/

#ifdef NEW_PROBLEM_TYPES
#include <stdio.h>
#include <iostream>
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
#include "ProblemType.h"
#include "EventHooks.h"
#include "phys_constants.h"

#define VCIRC_TABLE_LENGTH 10000

void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

int nlines(const char* fname) {

  FILE* fptr = fopen(fname, "r");
  int ch, n = 0;

  do
  {
    ch = fgetc(fptr);
    if(ch == '\n')
      n++;
  } while (ch != EOF);

  fclose(fptr);

  return n;
}

class ProblemType_AgoraRestart;

class AgoraRestartGrid : private grid
{
  friend class ProblemType_AgoraRestart;
};

class ProblemType_AgoraRestart : public EnzoProblemType
{
private:
  FLOAT SubgridLeft, SubgridRight;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  FLOAT CenterPosition[MAX_DIMENSION];
  FLOAT ScaleLength;
  FLOAT ScaleHeight;
  FLOAT CutOffRadius;
  FLOAT CutOffZ;
  float DiskMass;
  float GasFraction;
  float DiskTemperature;
  float HaloMass;
  float HaloTemperature;
  FLOAT VCircRadius[VCIRC_TABLE_LENGTH];
  float VCircVelocity[VCIRC_TABLE_LENGTH];
  int RefineAtStart;
  int UseMetallicityField;

public:
  ProblemType_AgoraRestart() : EnzoProblemType()
  {
    std::cout << "Creating problem type Agora Restart" << std::endl;
  }

  ~ProblemType_AgoraRestart() {}

  virtual int InitializeFromRestart(
    HierarchyEntry &TopGrid, TopGridData &MetaData)
  {
    return SUCCESS;
  }

  virtual int InitializeSimulation(
    FILE *fptr, FILE *Outfptr,
    HierarchyEntry &TopGrid, TopGridData &MetaData)
  {
    if(debug)
    {
      printf("Entering AgoraRestartInitialize\n");
      fflush(stdout);
    }

    char *DensName = "Density";
    char *TEName   = "TotalEnergy";
    char *GEName   = "GasEnergy";
    char *Vel1Name = "x-velocity";
    char *Vel2Name = "y-velocity";
    char *Vel3Name = "z-velocity";
    char *MetalName = "Metal_Density";

    /* local declarations */

    char line[MAX_LINE_LENGTH];
    int  i, ret, level;

    /* make sure it is 3D */

    if (MetaData.TopGridRank != 3)
    {
      printf("Cannot do AcoraRestart in %"ISYM" dimension(s)\n",
	     MetaData.TopGridRank);
      ENZO_FAIL("Agora Restart simulations must be 3D!");
    }

    for (i=0; i < MAX_DIMENSION; i++)
    {
      this->CenterPosition[i] = 0.5;
    }

    // These come from Oscar's sample output.  The units are:
    // Velocity: km/s
    // Mass: 10^9 Msun
    // Length: kpc
    // Temperature: K
    this->ScaleLength         = .0343218;
    this->ScaleHeight         = .00343218;
    this->CutOffRadius        = .2;
    this->CutOffZ             = 0.03;
    this->DiskMass            = 42.9661;
    this->GasFraction         = 0.2;
    this->DiskTemperature     = 1e4;
    this->HaloMass            = 0.10000;
    this->HaloTemperature     = this->DiskTemperature;
    this->RefineAtStart       = TRUE;
    this->UseMetallicityField = FALSE;

    /* Set no subgrids by default. */

    SubgridLeft = SubgridRight = 0.0;

    /* read input from file */
    while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
    {
      ret = 0;
      ret += sscanf(line, "AgoraRestartCenterPosition = %"PSYM" %"PSYM" %"PSYM,
		    CenterPosition, CenterPosition+1, CenterPosition+2);
      ret += sscanf(line, "AgoraRestartScaleLength = %"PSYM, &ScaleLength);
      ret += sscanf(line, "AgoraRestartScaleHeight = %"PSYM, &ScaleHeight);
      ret += sscanf(line, "AgoraRestartCutOffRadius = %"PSYM, &CutOffRadius);
      ret += sscanf(line, "AgoraRestartCutOffZ = %"PSYM, &CutOffZ);
      ret += sscanf(line, "AgoraRestartDiskMass = %"FSYM, &DiskMass);
      ret += sscanf(line, "AgoraRestartGasFraction = %"FSYM, &GasFraction);
      ret += sscanf(line, "AgoraRestartDiskTemperature = %"FSYM, &DiskTemperature);
      ret += sscanf(line, "AgoraRestartHaloMass = %"FSYM, &HaloMass);
      ret += sscanf(line, "AgoraRestartHaloTemperature = %"FSYM, &HaloTemperature);
      ret += sscanf(line, "AgoraRestartRefineAtStart = %"ISYM, &RefineAtStart);
      ret += sscanf(line, "AgoraRestartUseMetallicityField = %"ISYM, &UseMetallicityField);

      if (ret == 0 && strstr(line, "=") &&
	  (strstr(line, "AgoraRestart") || strstr(line, "TestProblem")) &&
	  line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr,
		"*** warning: the following parameter line was not interpreted:\n%s\n",
		line);

    } // end input from parameter file

    // Read in circular velocity table

    this->ReadInVcircData();

    /* set up top grid */

    float dummy_density = 1.0;
    float dummy_gas_energy = 1.0; // Only used if DualEnergyFormalism is True
    float dummy_total_energy = 1.0;
    float dummy_velocity[3] = {0.0, 0.0, 0.0};
    float dummy_b_field[3] = {0.0, 0.0, 0.0}; // Only set if HydroMethod = mhd_rk

    if (this->InitializeUniformGrid(
	  TopGrid.GridData, dummy_density, dummy_total_energy,
	  dummy_gas_energy, dummy_velocity, dummy_b_field) == FAIL)
    {
      ENZO_FAIL("Error in InitializeUniformGrid");
    }

    this->InitializeGrid(TopGrid.GridData, TopGrid, MetaData);

    this->InitializeParticles(TopGrid.GridData, TopGrid, MetaData);

    /* Convert minimum initial overdensity for refinement to mass
       (unless MinimumMass itself was actually set). */

    if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
      MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
      for (int dim = 0; dim < MetaData.TopGridRank; dim++)
	MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
	  float(MetaData.TopGridDims[dim]);
    }

    /* If requested, refine the grid to the desired level. */

    if (RefineAtStart)
    {
      /* Declare, initialize, and fill out the first level of the LevelArray. */
      LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
      for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
	LevelArray[level] = NULL;
      AddLevel(LevelArray, &TopGrid, 0);

      /* Add levels to the maximum depth or until no new levels are created,
	 and re-initialize the level after it is created. */
      for (level = 0; level < MaximumRefinementLevel; level++) {
	if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
	  fprintf(stderr, "Error in RebuildHierarchy.\n");
	  return FAIL;
	}
	if (LevelArray[level+1] == NULL)
	  break;
	LevelHierarchyEntry *Temp = LevelArray[level+1];
	while (Temp != NULL) {
	  if (this->InitializeGrid(Temp->GridData, TopGrid, MetaData) == FAIL)
	  {
	    ENZO_FAIL("Error in AgoraRestart->InitializeGrid");
	  }
	  Temp = Temp->NextGridThisLevel;
	} // end: loop over grids on this level
      } // end: loop over levels
    }

    /* set up field names and units */
    int count = 0;
    DataLabel[count++] = DensName;
    DataLabel[count++] = TEName;
    if (DualEnergyFormalism)
      DataLabel[count++] = GEName;
    DataLabel[count++] = Vel1Name;
    if(MetaData.TopGridRank > 1)
      DataLabel[count++] = Vel2Name;
    if(MetaData.TopGridRank > 2)
      DataLabel[count++] = Vel3Name;
    if (UseMetallicityField)
      DataLabel[count++] = MetalName;
    for (i = 0; i < count; i++)
      DataUnits[i] = NULL;

    if (MyProcessorNumber == ROOT_PROCESSOR)
    {
      fprintf(Outfptr, "AgoraRestartCenterPosition          = %"PSYM" %"PSYM" %"PSYM"\n",
	      CenterPosition[0], CenterPosition[1], CenterPosition[2]);
      fprintf(Outfptr, "AgoraRestartScaleLength             = %"PSYM"\n", ScaleLength);
      fprintf(Outfptr, "AgoraRestartScaleHeight             = %"PSYM"\n", ScaleHeight);
      fprintf(Outfptr, "AgoraRestartCutOffRadius            = %"PSYM"\n", CutOffRadius);
      fprintf(Outfptr, "AgoraRestartCutOffZ                 = %"PSYM"\n", CutOffZ);
      fprintf(Outfptr, "AgoraRestartDiskMass                = %"FSYM"\n", DiskMass);
      fprintf(Outfptr, "AgoraRestartGasFraction             = %"FSYM"\n", GasFraction);
      fprintf(Outfptr, "AgoraRestartDiskTemperature         = %"FSYM"\n", DiskTemperature);
      fprintf(Outfptr, "AgoraRestartHaloMass                = %"FSYM"\n", HaloMass);
      fprintf(Outfptr, "AgoraRestartHaloTemperature         = %"FSYM"\n", HaloTemperature);
      fprintf(Outfptr, "AgoraRestartRefineAtStart           = %"ISYM"\n", RefineAtStart);
      fprintf(Outfptr, "AgoraRestartUseMetallicityField     = %"ISYM"\n", UseMetallicityField);
    }

    return SUCCESS;

  } // InitializeSimulation

  int InitializeGrid(grid *thisgrid_orig, HierarchyEntry &TopGrid,
		     TopGridData &MetaData)
  {

    AgoraRestartGrid *thisgrid =
      static_cast<AgoraRestartGrid *>(thisgrid_orig);

    if (thisgrid->ProcessorNumber != MyProcessorNumber)
      return SUCCESS;

    /* Get units */
    float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
      TemperatureUnits=1;
    double MassUnits=1;

    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		 &TimeUnits, &VelocityUnits, &MassUnits, thisgrid->Time) == FAIL) {
      ENZO_FAIL("Error in GetUnits.");
    }

    /* Identify physical quantities */
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, MetalNum;

    if (thisgrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					     Vel3Num, TENum) == FAIL) {
      fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
      ENZO_FAIL("");
    }

    int MetallicityField = FALSE;
    if ((MetalNum = FindField(
	   Metallicity, thisgrid->FieldType, thisgrid->NumberOfBaryonFields)
	  ) != -1)
      MetallicityField = TRUE;
    else
      MetalNum = 0;

    int dim, i, j, k, size, index=0;
    float RhoZero, DiskGasEnergy, DiskDensity, HaloGasEnergy, HaloDensity,
      BoxVolume, vcirc;
    FLOAT x, y, z, radius, xy_radius, cellwidth;

    /* Compute size of this grid */
    size = 1;
    for (dim = 0; dim < thisgrid->GridRank; dim++)
      size *= thisgrid->GridDimension[dim];

    cellwidth = thisgrid->CellWidth[0][0];

    /* Compute the size of the box */
    BoxVolume = 1.;
    for (dim = 0; dim < TopGrid.GridData->GetGridRank(); dim++)
      BoxVolume *= (DomainRightEdge[dim] - DomainLeftEdge[dim]);

    /* Find global physical properties */
    RhoZero = this->DiskMass * this->GasFraction / (4.*pi) /
      (POW((this->ScaleLength),2)*(this->ScaleHeight));

    HaloGasEnergy = this->HaloTemperature / Mu / (Gamma - 1) /
      TemperatureUnits;

    HaloDensity = this->HaloMass / BoxVolume;

    DiskGasEnergy = this->DiskTemperature / Mu / (Gamma - 1) /
      TemperatureUnits;

    /* Loop over the mesh. */

    for (k = 0; k < thisgrid->GridDimension[2]; k++)
    {
      for (j = 0; j < thisgrid->GridDimension[1]; j++)
      {
	for (i = 0; i < thisgrid->GridDimension[0]; i++, index++)
	{
	  /* Compute position */

	  x = (thisgrid->CellLeftEdge[0][i] + 0.5*thisgrid->CellWidth[0][i]) *
	    LengthUnits;
	  y = (thisgrid->CellLeftEdge[1][j] + 0.5*thisgrid->CellWidth[1][j]) *
	    LengthUnits;
	  z = (thisgrid->CellLeftEdge[2][k] + 0.5*thisgrid->CellWidth[2][k]) *
	    LengthUnits;

	  x -= this->CenterPosition[0]*LengthUnits;
	  y -= this->CenterPosition[1]*LengthUnits;
	  z -= this->CenterPosition[2]*LengthUnits;

	  radius = sqrt(POW(x, 2) +
			POW(y, 2) +
			POW(z, 2) );

	  xy_radius = sqrt(POW(x, 2) +
			   POW(y, 2) );

	  /* Find disk density, halo density and internal energy */

	  DiskDensity = gauss_mass(RhoZero, x/LengthUnits, y/LengthUnits,
				   z/LengthUnits, cellwidth) / POW(cellwidth, 3);

	  if ( HaloDensity*HaloTemperature > DiskDensity*DiskTemperature )
	  {
	    thisgrid->BaryonField[DensNum][index] = HaloDensity;
	    thisgrid->BaryonField[GENum][index] = HaloGasEnergy;
	  }
	  else // Ok, we're in the disk
	  {
	    thisgrid->BaryonField[DensNum][index] = DiskDensity;
	    thisgrid->BaryonField[GENum][index] = DiskGasEnergy;

	    vcirc = this->InterpolateVcircTable(xy_radius);

	    thisgrid->BaryonField[Vel1Num][index] =
	      -vcirc*y/xy_radius/VelocityUnits;
	    thisgrid->BaryonField[Vel2Num][index] =
	      vcirc*x/xy_radius/VelocityUnits;
	  }

	} // i
      } // j
    } // k

    return SUCCESS;

  } // InitializeGrid

  void InitializeParticles(grid *thisgrid_orig, HierarchyEntry &TopGrid,
			  TopGridData &MetaData)
  {
    AgoraRestartGrid *thisgrid =
      static_cast<AgoraRestartGrid *>(thisgrid_orig);

    // Determine the number of particles of each type
    int nBulge, nDisk, nHalo, nParticles;
    nBulge = nlines("bulge.dat");
    nDisk = nlines("disk.dat");
    nHalo = nlines("halo.dat");
    nParticles = nBulge + nDisk + nHalo;

    // Initialize particle arrays
    PINT *Number = new PINT[nParticles];
    int *Type = new int[nParticles];
    FLOAT *Position[MAX_DIMENSION];
    float *Velocity[MAX_DIMENSION];
    for (int i = 0; i < thisgrid->GridRank; i++)
    {
      Position[i] = new FLOAT[nParticles];
      Velocity[i] = new float[nParticles];
    }
    float *Mass = new float[nParticles];
    float *Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
    for (int i = 0; i < NumberOfParticleAttributes; i++)
    {
      Attribute[i] = new float[nParticles];
      for (int j = 0; j < nParticles; j++)
	Attribute[i][j] = FLOAT_UNDEFINED;
    }

    FLOAT dx = thisgrid->CellWidth[0][0];

    // Read them in and assign them as we go
    int count = 0;
    this->ReadParticlesFromFile(
      Number, Type, Position, Velocity, Mass,
      "bulge.dat", PARTICLE_TYPE_STAR, count, dx);
    this->ReadParticlesFromFile(
      Number, Type, Position, Velocity, Mass,
      "disk.dat", PARTICLE_TYPE_STAR, count, dx);
    this->ReadParticlesFromFile(
      Number, Type, Position, Velocity, Mass,
      "halo.dat", PARTICLE_TYPE_DARK_MATTER, count, dx);

    thisgrid->SetNumberOfParticles(count);
    thisgrid->SetParticlePointers(Mass, Number, Type, Position,
				  Velocity, Attribute);
    MetaData.NumberOfParticles = count;
  }

  float gauss_mass(
    float RhoZero, FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT cellwidth)
  {
    // Computes the total mass in a given cell by integrating the density
    // profile using 5-point Gaussian quadrature.
    // http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
    FLOAT EvaluationPoints [5] = {-0.90617985,-0.53846931,0.0,0.53846931,0.90617985};
    FLOAT Weights [5] = {0.23692689,0.47862867,0.56888889,0.47862867,0.23692689};
    FLOAT xResult [5];
    FLOAT yResult [5];
    FLOAT r, z;
    float Mass = 0;
    int i,j,k;

    for (i=0;i<5;i++)
    {
      xResult[i] = 0.0;
      for (j=0;j<5;j++)
      {
	yResult[j] = 0.0;
	for (k=0;k<5;k++)
	{
	  r = sqrt((POW(xpos+EvaluationPoints[i]*cellwidth/2.0, 2.0) +
		    POW(ypos+EvaluationPoints[j]*cellwidth/2.0, 2.0) ) );
	  z = fabs(zpos+EvaluationPoints[k]*cellwidth/2.0);
	  yResult[j] +=
	    cellwidth/2.0 * Weights[k] * RhoZero *
	    PEXP(-r/this->ScaleLength) *
	    PEXP(-fabs(z)/this->ScaleHeight);
	}
	xResult[i] += cellwidth/2.0*Weights[j]*yResult[j];
      }
      Mass += cellwidth/2.0*Weights[i]*xResult[i];
    }
    return Mass;
  }

  void ReadInVcircData(void)
  {
    FILE *fptr;
    char line[MAX_LINE_LENGTH];
    int i=0, ret;
    float vcirc;
    FLOAT rad;

    fptr = fopen("vcirc.dat" , "r");

    while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
    {
      ret += sscanf(line, "%"PSYM" %"FSYM, &rad, &vcirc);
      this->VCircRadius[i] = rad*3.08567758e21; // 3.08567758e21 = kpc/cm
      this->VCircVelocity[i] = vcirc*1e5; // 1e5 = (km/s)/(cm/s)
      i += 1;
    }

    fclose(fptr);
  } // ReadInVcircData

  float InterpolateVcircTable(FLOAT radius)
  {
    int i;

    for (i = 0; i < VCIRC_TABLE_LENGTH; i++)
      if (radius < this->VCircRadius[i])
	break;

    if (i == 0)
      return (VCircVelocity[i]) * (radius - VCircRadius[0]) / VCircRadius[0];
    else if (i == VCIRC_TABLE_LENGTH)
      ENZO_FAIL("Fell off the circular velocity interpolation table");

    // we know the radius is between i and i-1
    return VCircVelocity[i-1] +
      (VCircVelocity[i] - VCircVelocity[i-1]) *
      (radius - VCircRadius[i-1])  /
      (VCircRadius[i] - VCircRadius[i-1]);
  }

  int ReadParticlesFromFile(PINT *Number, int *Type, FLOAT *Position[],
			    float *Velocity[], float* Mass, const char* fname,
			    Eint32 particle_type, int &c, FLOAT dx)
  {
    FILE *fptr;
    char line[MAX_LINE_LENGTH];
    int ret;
    FLOAT x, y, z;
    float vx, vy, vz;
    double mass;

    float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
      TemperatureUnits=1;
    double MassUnits=1;

    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		 &TimeUnits, &VelocityUnits, &MassUnits, 0) == FAIL) {
      ENZO_FAIL("Error in GetUnits.");
    }

    const FLOAT kpc_cm = 3.08567758e21;
    const float km_cm = 1e5;
    const float msun_g = 1.9891e33;

    fptr = fopen(fname, "r");

    while(fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
    {
      ret +=
	sscanf(line,
	       "%"PSYM" %"PSYM" %"PSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
	       &x, &y, &z, &vx, &vy, &vz, &mass);

      Position[0][c] = x * kpc_cm / LengthUnits + this->CenterPosition[0];
      Position[1][c] = y * kpc_cm / LengthUnits + this->CenterPosition[1];
      Position[2][c] = z * kpc_cm / LengthUnits + this->CenterPosition[2];

      Velocity[0][c] = vx * km_cm / VelocityUnits;
      Velocity[1][c] = vy * km_cm / VelocityUnits;
      Velocity[2][c] = vz * km_cm / VelocityUnits;

      // Particle masses are actually densities.
      Mass[c] = mass * 1e9 * msun_g / MassUnits / dx / dx / dx;
      Type[c] = particle_type;
      Number[c] = c++;
    }

    fclose(fptr);

    return c;
  } // ReadParticlesFromFile

}; // class declaration

//.. register:
namespace {
    EnzoProblemType_creator_concrete<ProblemType_AgoraRestart>
        agora_restart("AgoraRestart");
}

#endif
