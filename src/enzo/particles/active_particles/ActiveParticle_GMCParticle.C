/***********************************************************************
/
/ GMC Particle
/
************************************************************************/

#include "ActiveParticle_GMCParticle.h"

#ifdef NEW_CONFIG

#include "ParameterControl/ParameterControl.h"
extern Configuration Param;

/* Set default parameter values. */

const char config_gmc_particle_defaults[] =
"### GMC PARTICLE DEFAULTS ###\n"
"\n"
"Physics: {\n"
"    ActiveParticles: {\n"
"        GMCParticle: {\n"
"            OverflowFactor       = 1.01;\n"
"            LinkingLength        = 4;\n   "
"            AccretionRadius      = 4;\n   "
"        };\n"
"    };\n"
"};\n";

#endif

void mt_init(unsigned_int seed);
unsigned_long_int mt_random(void);
int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);

ActiveParticleType_GMCParticle::ActiveParticleType_GMCParticle(void) : ActiveParticleType_AccretingParticle() {
  R = M = sigma = 1.0;
  Mdot = MdotAcc = sigmadot = Rdot = Rddot = Mddot = sigmadotAcc = 0.0;
  tau = 0.0;
  Massoc = Mstar = 0.0;
  Tco = 0.0;
  dtau = 1.0e-4;
  dtauSave = 0.0;
  nHIIreg = 0.0;
  MstarRemain = 0.0;
  sigmadot_noacc = Rddot_noacc = Rdot_noacc = Ecl = Ecl_noacc =
    R_noacc = M_noacc = sigma_noacc = Eacc = sigmaISM = 0.0;
  MdotAcc = 0.0;
  dtauOk = HIIregEsc = dissoc = dtauFloor = 0;
  this->CalculateDerivedParameters();
  // This is done twice on purpose.
  this->UpdateDerivedParameters();
  this->UpdateDerivedParameters();
  Ecl = 
    0.5 * aI * M * Rdot * Rdot +
    2.4 * M * sigma * sigma +
    1.5 * M * Mach0 -
    5.0 * (0.6 * aprime * (1 - etaB*etaB) - chi) * M / (R * R * avir0);
}

int ActiveParticleType_GMCParticle::InitializeParticleType()
{
  if (ActiveParticleType_AccretingParticle::InitializeParticleType() == FAIL)
    ENZO_FAIL("AccretingParticle Initialize failesd");
  
  FILE *fpAccTable = fopen("GMCtable.in","r");
  if (fpAccTable == NULL)
    ENZO_FAIL("Error opening GMCtable.in\n");
  
  int i = 0;
  
  float theta, zeta, f, xi, aprime, chi, gamma;

  accTable.zetaLook = new float[NUMZETA];
  accTable.fLook = new float[NUMZETA];
  accTable.xiLook = new float[NUMZETA];
  accTable.aprimeLook = new float[NUMZETA];
  accTable.chiLook = new float[NUMZETA];
  accTable.gammaLook = new float[NUMZETA];

  while (!feof(fpAccTable)){
    fscanf(fpAccTable,"%lf, %lf, %lf, %lf, %lf, %lf",&zeta,&f,&xi,&aprime,&chi,&gamma);
    accTable.zetaLook[i]=zeta;
    accTable.fLook[i]=f;
    accTable.xiLook[i]=xi;
    accTable.aprimeLook[i]=aprime;
    accTable.chiLook[i]=chi;
    accTable.gammaLook[i]=gamma;
    i++;
  }

  fclose(fpAccTable);
  
  i = 0;

  FILE *fpRadSolTable = fopen("RadSolTable","r");

  if (fpRadSolTable == NULL)
      ENZO_FAIL("Error opening RadSolTable.in\n");

  double tau, xShell, xShellPrime;

  radSolTable.tauLook = new float[NUMTAU];
  radSolTable.xShellLook = new float[NUMTAU];
  radSolTable.xPrimeShellLook = new float[NUMTAU];

  fpRadSolTable = fopen("/Users/goldbaum/Documents/RadSolTable","r");
  while (!feof(fpRadSolTable)){
    fscanf(fpRadSolTable,"%lf, %lf, %lf",&tau,&xShell,&xShellPrime);
    radSolTable.tauLook[i] = tau;
    radSolTable.xShellLook[i] = xShell;
    radSolTable.xPrimeShellLook[i] = xShellPrime;
    i++;
  }

  fclose(fpRadSolTable);

  if (GMCParticleRNGSeed == INT_UNDEFINED) {
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      GMCParticleRNGSeed = time(NULL);
    }
    CommunicationBroadcastValue(&GMCParticleRNGSeed, ROOT_PROCESSOR);
  }

  mt_init(GMCParticleRNGSeed);

  // If this is a restart, reset the RNG to the 
  // appropriate position in the RN stream

  unsigned_long_int trash;
  for (int i = 0; i < 1 + GMCParticleRNGCalls; i++)
    trash = mt_random();

  return SUCCESS;
}

#define CII 9.74e5 /* Sound speed in ionized gas, assuming T = 7000 K and mu = 0.61 */

int ActiveParticleType_GMCParticle :: CalculateDerivedParameters() {
  aI     = (3 - krho) / (5 - krho);
  a      = (15 - 5*krho) / (15 - 6.0*krho);
  aprime = a;
  Mach0  = sigma0 / cs;
  avir0  = (5 * sigma0 * sigma0 * R0) / (GravConst * M0);
  etaG   = 3 * aprime * (1.0 - etaB*etaB) / (aI * avir0);
  etaP   = 4 * pi * R0 * R0 * R0 * Pamb / 
    (aI * M0 * sigma0 * sigma0);
  etaE   = (5.0 - krho) / (4 - krho) * 2 * CII / sigma0;
  etaA   = (5.0 - krho) / (4 - krho) * SQRT(10.0/avir0);
  t0     = R0 / sigma0;
  f      = 1.0;
  etaI   = 10/(aI*avir0);
  xi     = 1.11483349347388;

  return SUCCESS;
}

double logInterpolate(int zetaindex, double interpArray[NUMZETA],
                      double zeta, double zetamin, double zetamax){
  double gridArea = log10(zetamax) - log10(zetamin);
  return( (interpArray[zetaindex]*(log10(zetamax) - log10(zeta))
           + interpArray[zetaindex+1]*(log10(zeta) - log10(zetamin)))/gridArea ) ;
}


int ActiveParticleType_GMCParticle :: UpdateDerivedParameters() {
  double rho   = (M*M0) / (f*4./3.*pi*(R*R*R*R0*R0*R0));
  double tff   = SQRT(3*pi / (32*GravConst*rho));
  double tauff = tff / t0;
  double zeta  = MdotAcc*tauff/M;

  if (zeta >= 1000)
    ENZO_FAIL("GMCParticle: Zeta is greater than 1000!\n");
  if (zeta >= 0.001) {
    int zetaindex = (log10(zeta)+3)*10; // Converting to the log space in the zeta lookup table
    double zetamin = accTable.zetaLook[zetaindex];
    double zetamax = accTable.zetaLook[zetaindex+1];

    aprime = logInterpolate(zetaindex,accTable.aprimeLook, zeta, zetamin, zetamax);
    f = logInterpolate(zetaindex,accTable.fLook,zeta,zetamin,zetamax);
    xi = logInterpolate(zetaindex,accTable.xiLook,zeta,zetamin,zetamax);
    chi = logInterpolate(zetaindex,accTable.chiLook,zeta,zetamin,zetamax);
    gamma = logInterpolate(zetaindex,accTable.gammaLook,zeta,zetamin,zetamax);

    etaI = (10 * gamma) / (f * aI * avir0);
  }
  else {
    aprime = a;
    f = 1.0;
    xi = 1.11484;
    chi = 0.0;
    // The above values are the analytic values for krho = 1 and zeta = 0.
    // If using a different value of krho, they will need to be updated.
  }
  etaG = 3 * aprime * (1.0 - etaB * etaB) / (aI * avir0);
  etaA = (5.0 - krho) / (4.0 - krho) * sqrt(xi * xi * 10.0 / (avir0 * f));
  
  return SUCCESS;
}

int ActiveParticleType_GMCParticle::WriteToOutput(ActiveParticleType **these_particles, int n, int GridRank, hid_t group_id)
{
  /* Create a new subgroup within the active particle group for active particles of type AccretingParticle */
  hid_t AccretingParticleGroupID = H5Gcreate(group_id,"GMCParticle",0);

  writeScalarAttribute(AccretingParticleGroupID,HDF5_INT,"Number of GMC Particles",&n);  

  char *ParticlePositionLabel[] =
     {"position_x", "position_y", "position_z"};
  char *ParticleVelocityLabel[] =
     {"velocity_x", "velocity_y", "velocity_z"};

  /* Create temporary buffers to store particle data */

  FLOAT *Position[MAX_DIMENSION];
  float *Velocity[MAX_DIMENSION]; 
  double *Mass         = new double[n];
  float *BirthTime     = new float[n];
  float *DynamicalTime = new float[n];
  float *Metallicity   = new float[n];
  float *AccretionRate = new float[n];
  float *M0            = new float[n];
  float *R0            = new float[n];
  float *sigma0        = new float[n];
  float *R             = new float[n];
  float *Rdot          = new float[n];
  float *M             = new float[n];
  float *Mstar         = new float[n];
  float *Massoc        = new float[n];
  float *MdotStar      = new float[n];
  float *MdotHII       = new float[n];
  float *MstarRemain   = new float[n];
  float *sigma         = new float[n];
  float *tau           = new float[n];
  float *dtau          = new float[n];
  float *Eacc          = new float[n];

  PINT *ID = new PINT[n];

  int i,dim;

  for (dim = 0; dim < GridRank; dim++) {
    Position[dim] = new FLOAT[n];
    Velocity[dim] = new float[n];
  }

  hsize_t TempInt;
  TempInt = n;
    
  ActiveParticleType_GMCParticle *ParticleToWrite;
  for (i=0;i<n;i++) {
    ParticleToWrite = static_cast<ActiveParticleType_GMCParticle*>(these_particles[i]);
    for (dim = 0; dim < GridRank; dim++) {
      Position[dim][i] = ParticleToWrite->pos[dim];
      Velocity[dim][i] = ParticleToWrite->vel[dim];
    }
    Mass[i]               = ParticleToWrite->Mass;
    BirthTime[i]          = ParticleToWrite->BirthTime;
    DynamicalTime[i]      = ParticleToWrite->DynamicalTime;
    Metallicity[i]        = ParticleToWrite->Metallicity;
    ID[i]                 = ParticleToWrite->Identifier;
    AccretionRate[i]      = ParticleToWrite->AccretionRate;
    M0[i]                 = ParticleToWrite->M0;
    R0[i]                 = ParticleToWrite->R0;
    sigma0[i]             = ParticleToWrite->sigma0;
    R[i]                  = ParticleToWrite->R;
    Rdot[i]               = ParticleToWrite->Rdot;
    M[i]                  = ParticleToWrite->M;
    Mstar[i]              = ParticleToWrite->Mstar;
    Massoc[i]             = ParticleToWrite->Massoc;
    MdotStar[i]           = ParticleToWrite->MdotStar;
    MdotHII[i]            = ParticleToWrite->MdotHII;
    MstarRemain[i]        = ParticleToWrite->MstarRemain;
    sigma[i]              = ParticleToWrite->sigma;
    tau[i]                = ParticleToWrite->tau;
    dtau[i]               = ParticleToWrite->dtau;
    Eacc[i]               = ParticleToWrite->Eacc;
  }

  for (dim = 0; dim < GridRank; dim++) {
    WriteDataset(1,&TempInt,ParticlePositionLabel[dim],
		 AccretingParticleGroupID, HDF5_FILE_PREC, (VOIDP) Position[dim]);
  }
  
  for (dim = 0; dim < GridRank; dim++) {
    WriteDataset(1,&TempInt,ParticleVelocityLabel[dim],
		  AccretingParticleGroupID, HDF5_REAL, (VOIDP) Velocity[dim]);
  }
  
  WriteDataset(1,&TempInt,"mass",AccretingParticleGroupID,HDF5_REAL,(VOIDP) Mass);
  WriteDataset(1,&TempInt,"creation_time",AccretingParticleGroupID,HDF5_REAL,(VOIDP) BirthTime);
  WriteDataset(1,&TempInt,"dynamical_time",AccretingParticleGroupID,HDF5_REAL,(VOIDP) DynamicalTime);
  WriteDataset(1,&TempInt,"metallicity_fraction",AccretingParticleGroupID,HDF5_REAL,(VOIDP) Metallicity);
  WriteDataset(1,&TempInt,"identifier",AccretingParticleGroupID,HDF5_PINT,(VOIDP) ID);
  WriteDataset(1,&TempInt,"accretion_rate",AccretingParticleGroupID,HDF5_REAL,(VOIDP) AccretionRate);
  WriteDataset(1,&TempInt,"gmcevol_M0",AccretingParticleGroupID,HDF5_REAL,(VOIDP) M0);
  WriteDataset(1,&TempInt,"gmcevol_R0",AccretingParticleGroupID,HDF5_REAL,(VOIDP) R0);
  WriteDataset(1,&TempInt,"gmcevol_sigma0",AccretingParticleGroupID,HDF5_REAL,(VOIDP) sigma0);
  WriteDataset(1,&TempInt,"gmcevol_R",AccretingParticleGroupID,HDF5_REAL,(VOIDP) R);
  WriteDataset(1,&TempInt,"gmcevol_Rdot",AccretingParticleGroupID,HDF5_REAL,(VOIDP) Rdot);
  WriteDataset(1,&TempInt,"gmcevol_M",AccretingParticleGroupID,HDF5_REAL,(VOIDP) M);
  WriteDataset(1,&TempInt,"gmcevol_Mstar",AccretingParticleGroupID,HDF5_REAL,(VOIDP) Mstar);  
  WriteDataset(1,&TempInt,"gmcevol_Massoc",AccretingParticleGroupID,HDF5_REAL,(VOIDP) Massoc);  
  WriteDataset(1,&TempInt,"gmcevol_MdotStar",AccretingParticleGroupID,HDF5_REAL,(VOIDP) MdotStar);  
  WriteDataset(1,&TempInt,"gmcevol_MdotHII",AccretingParticleGroupID,HDF5_REAL,(VOIDP) MdotHII);  
  WriteDataset(1,&TempInt,"gmcevol_MstarRemain",AccretingParticleGroupID,HDF5_REAL,(VOIDP) MstarRemain);  
  WriteDataset(1,&TempInt,"gmcevol_sigma",AccretingParticleGroupID,HDF5_REAL,(VOIDP) sigma);  
  WriteDataset(1,&TempInt,"gmcevol_tau",AccretingParticleGroupID,HDF5_REAL,(VOIDP) tau);  
  WriteDataset(1,&TempInt,"gmcevol_dtau",AccretingParticleGroupID,HDF5_REAL,(VOIDP) dtau);  
  WriteDataset(1,&TempInt,"gmcevol_Eacc",AccretingParticleGroupID,HDF5_REAL,(VOIDP) Eacc);  
  
  /* Clean up */

  for (dim = 0; dim < GridRank; dim++) {
    delete[] Position[dim];
    delete[] Velocity[dim];
  }
  delete[] Mass;
  delete[] BirthTime;
  delete[] DynamicalTime;
  delete[] Metallicity;
  delete[] AccretionRate;
  delete[] M0;
  delete[] R0;
  delete[] sigma0;
  delete[] R;
  delete[] Rdot;
  delete[] M;
  delete[] Mstar;
  delete[] Massoc;
  delete[] MdotStar;
  delete[] MdotHII;
  delete[] MstarRemain;
  delete[] sigma;
  delete[] tau;
  delete[] dtau;
  delete[] Eacc;
  delete[] ID;
  H5Gclose(AccretingParticleGroupID);

  return SUCCESS;
}

int ActiveParticleType_GMCParticle::ReadFromOutput(ActiveParticleType **&particles_to_read, int &n, 
							 int GridRank, hid_t group_id)
{
  int i,dim;
  hsize_t TempInt;

  hid_t AccretingParticleGroupID = H5Gopen(group_id,"GMCParticle");

  readAttribute(AccretingParticleGroupID,HDF5_INT,"Number of GMC Particles",&n);

  particles_to_read = new ActiveParticleType*[n];

  char *ParticlePositionLabel[] =
     {"position_x", "position_y", "position_z"};
  char *ParticleVelocityLabel[] =
     {"velocity_x", "velocity_y", "velocity_z"};

  FLOAT *Position[MAX_DIMENSION];
  float *Velocity[MAX_DIMENSION];
  double *Mass         = new double[n];
  float *BirthTime     = new float[n];
  float *DynamicalTime = new float[n];
  float *Metallicity   = new float[n];
  float *AccretionRate = new float[n];
  float *M0            = new float[n];
  float *R0            = new float[n];
  float *sigma0        = new float[n];
  float *R             = new float[n];
  float *Rdot          = new float[n];
  float *M             = new float[n];
  float *Mstar         = new float[n];
  float *Massoc        = new float[n];
  float *MdotStar      = new float[n];
  float *MdotHII       = new float[n];
  float *MstarRemain   = new float[n];
  float *sigma         = new float[n];
  float *tau           = new float[n];
  float *dtau          = new float[n];
  float *Eacc          = new float[n];
  PINT  *ID            = new PINT[n];

  for (dim = 0; dim < GridRank; dim++) {
    Position[dim] = new FLOAT[n];
    Velocity[dim] = new float[n];
  }
  
  TempInt = n;
  
  for (dim = 0; dim < GridRank; dim++) {
    ReadDataset(1,&TempInt,ParticlePositionLabel[dim],
		  AccretingParticleGroupID, HDF5_FILE_PREC, (VOIDP) Position[dim]);
  }

  for (dim = 0; dim < GridRank; dim++) {
    ReadDataset(1,&TempInt,ParticleVelocityLabel[dim],
		  AccretingParticleGroupID, HDF5_REAL, (VOIDP) Velocity[dim]);
  }
  ReadDataset(1,&TempInt,"mass",AccretingParticleGroupID,HDF5_R8,(VOIDP) Mass);
  ReadDataset(1,&TempInt,"creation_time",AccretingParticleGroupID,HDF5_REAL,(VOIDP) BirthTime);
  ReadDataset(1,&TempInt,"dynamical_time",AccretingParticleGroupID,HDF5_REAL,(VOIDP) DynamicalTime);
  ReadDataset(1,&TempInt,"metallicity_fraction",AccretingParticleGroupID,HDF5_REAL,(VOIDP) Metallicity);
  ReadDataset(1,&TempInt,"identifier",AccretingParticleGroupID,HDF5_PINT,(VOIDP) ID);
  ReadDataset(1,&TempInt,"accretion_rate",AccretingParticleGroupID,HDF5_REAL,(VOIDP) AccretionRate);
  ReadDataset(1,&TempInt,"gmcevol_M0",AccretingParticleGroupID,HDF5_REAL,(VOIDP) M0);
  ReadDataset(1,&TempInt,"gmcevol_R0",AccretingParticleGroupID,HDF5_REAL,(VOIDP) R0);
  ReadDataset(1,&TempInt,"gmcevol_sigma0",AccretingParticleGroupID,HDF5_REAL,(VOIDP) sigma0);
  ReadDataset(1,&TempInt,"gmcevol_R",AccretingParticleGroupID,HDF5_REAL,(VOIDP) R);
  ReadDataset(1,&TempInt,"gmcevol_Rdot",AccretingParticleGroupID,HDF5_REAL,(VOIDP) Rdot);
  ReadDataset(1,&TempInt,"gmcevol_M",AccretingParticleGroupID,HDF5_REAL,(VOIDP) M);
  ReadDataset(1,&TempInt,"gmcevol_Mstar",AccretingParticleGroupID,HDF5_REAL,(VOIDP) Mstar);  
  ReadDataset(1,&TempInt,"gmcevol_Massoc",AccretingParticleGroupID,HDF5_REAL,(VOIDP) Massoc);  
  ReadDataset(1,&TempInt,"gmcevol_MdotStar",AccretingParticleGroupID,HDF5_REAL,(VOIDP) MdotStar);  
  ReadDataset(1,&TempInt,"gmcevol_MdotHII",AccretingParticleGroupID,HDF5_REAL,(VOIDP) MdotHII);  
  ReadDataset(1,&TempInt,"gmcevol_MstarRemain",AccretingParticleGroupID,HDF5_REAL,(VOIDP) MstarRemain);  
  ReadDataset(1,&TempInt,"gmcevol_sigma",AccretingParticleGroupID,HDF5_REAL,(VOIDP) sigma);  
  ReadDataset(1,&TempInt,"gmcevol_tau",AccretingParticleGroupID,HDF5_REAL,(VOIDP) tau);  
  ReadDataset(1,&TempInt,"gmcevol_dtau",AccretingParticleGroupID,HDF5_REAL,(VOIDP) dtau);  
  ReadDataset(1,&TempInt,"gmcevol_Eacc",AccretingParticleGroupID,HDF5_REAL,(VOIDP) Eacc);  

  for (i = 0; i < n; i++) {
    ActiveParticleType_GMCParticle *np = new ActiveParticleType_GMCParticle();
    np->Mass              = Mass[i];
    np->type              = np->GetEnabledParticleID();
    np->BirthTime         = BirthTime[i];
    np->DynamicalTime     = DynamicalTime[i];
    np->Metallicity       = Metallicity[i];
    np->Identifier        = ID[i];
    for (dim = 0; dim < GridRank; dim++){
      np->pos[dim] = Position[dim][i];
      np->vel[dim] = Velocity[dim][i];
    }
    np->AccretionRate     = AccretionRate[i];
    np->M0                = M0[i];
    np->R0                = R0[i];
    np->sigma0            = sigma0[i];
    np->R                 = R[i];            
    np->Rdot              = Rdot[i];
    np->M                 = M[i];
    np->Mstar             = Mstar[i];
    np->Massoc            = Massoc[i];
    np->MdotStar          = MdotStar[i];
    np->MdotHII           = MdotHII[i];
    np->MstarRemain       = MstarRemain[i];
    np->sigma             = sigma[i];
    np->tau               = tau[i];
    np->dtau              = dtau[i];
    np->Eacc              = Eacc[i];
    particles_to_read[i] = static_cast<ActiveParticleType*>(np);
  }

  delete[] Mass;
  delete[] BirthTime;
  delete[] DynamicalTime;
  delete[] Metallicity;
  delete[] ID;
  delete[] AccretionRate;
  delete[] M0;
  delete[] R0;
  delete[] sigma0;
  delete[] R;
  delete[] Rdot;
  delete[] M;
  delete[] Mstar;
  delete[] Massoc;
  delete[] MdotStar;
  delete[] MdotHII;
  delete[] MstarRemain;
  delete[] sigma;
  delete[] tau;
  delete[] dtau;
  delete[] Eacc;
  for (dim = 0; dim < GridRank; dim++) {
    delete[] Position[dim];
    delete[] Velocity[dim];
  }
  H5Gclose(AccretingParticleGroupID);
  
  return SUCCESS;
}

namespace {
  ActiveParticleType_info *GMCParticleInfo = 
    register_ptype <ActiveParticleType_GMCParticle, GMCParticleBufferHandler> 
    ("GMCParticle");
}
