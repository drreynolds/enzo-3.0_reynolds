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

ActiveParticleType_GMCParticle::ActiveParticleType_GMCParticle(void) : ActiveParticleType_AccretingParticle() {
  // Shouldn't initialize a new particle with this constructor but we
  // set these attributes to prevent undefined behavior
  R = Rdot = M = Mstar = Massoc = 
    MdotStar = MdotHII = MstarRemain = 
    sigma = tau = dtau = Eacc = ReservoirRatio = -1.0;
  R0 = M0 = sigma0 = -1; 
  nHIIreg = 0;
}

ActiveParticleType_GMCParticle::ActiveParticleType_GMCParticle(ActiveParticleType_AccretingParticle *ap,
							       ActiveParticleFormationData &data) 
  : ActiveParticleType_AccretingParticle(ap) {
  // These are scaled dimensionless variables, see Goldbaum et al. 2011.
  R = M = sigma = 1.0;
  // Set initial radius to one cell size
  FLOAT dx = data.CellSize;
  R0 = dx*data.LengthUnits; 
  M0 = ap->ReturnMass()*data.MassUnits;
  sigma0 = SQRT(2.0*GravConst*M0/(5*R0)); // Assuming alpha_vir,0 = 2.0

  Mstar = Massoc = MdotStar = MdotHII = MstarRemain = 0.0;
  ReservoirRatio = 1.0;
  tau = 0.0;
  dtau = 1.0e-4;
  Eacc = 0.0;
  nHIIreg = 0;
}

ActiveParticleType_GMCParticle::ActiveParticleType_GMCParticle
(ActiveParticleType_GMCParticle* part) : ActiveParticleType_AccretingParticle(part) {
  M0             = part->M0;
  R0             = part->R0;
  sigma0         = part->sigma0;
  
  R              = part->R;
  M              = part->M;
  MdotStar       = part->MdotStar;
  MdotHII        = part->MdotHII;
  MstarRemain    = part->MstarRemain;
  Massoc         = part->Massoc;
  sigma          = part->sigma;
  tau            = part->tau;
  dtau           = part->dtau;
  Mstar          = part->Mstar;
  Eacc           = part->Eacc;  
  ReservoirRatio = part->ReservoirRatio;

  nHIIreg        = part->nHIIreg;
}

int ActiveParticleType_GMCParticle::InitializeParticleType()
{
  if (ActiveParticleType_AccretingParticle::InitializeParticleType() == FAIL)
    ENZO_FAIL("AccretingParticle Initialize failed!");
  
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
    fscanf(fpAccTable,"%"FSYM", %"FSYM", %"FSYM", %"FSYM", %"FSYM", %"FSYM,&zeta,&f,&xi,&aprime,&chi,&gamma);
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

  FILE *fpRadSolTable = fopen("RadSolTable.in","r");

  if (fpRadSolTable == NULL)
      ENZO_FAIL("Error opening RadSolTable.in\n");

  double tau, xShell, xShellPrime;

  radSolTable.tauLook = new float[NUMTAU];
  radSolTable.xShellLook = new float[NUMTAU];
  radSolTable.xPrimeShellLook = new float[NUMTAU];

  while (!feof(fpRadSolTable)){
    fscanf(fpRadSolTable,"%lf, %lf, %lf",&tau,&xShell,&xShellPrime);
    radSolTable.tauLook[i] = tau;
    radSolTable.xShellLook[i] = xShell;
    radSolTable.xPrimeShellLook[i] = xShellPrime;
    i++;
  }

  fclose(fpRadSolTable);

  if (GMCParticleRNGSeed == INT_UNDEFINED) {
    GMCParticleRNGSeed = time(NULL);
  }

  mt_init(GMCParticleRNGSeed);

  // If this is a restart, reset the RNG to the 
  // appropriate position in the stream

  unsigned_long_int trash;
  for (int i = 0; i < 1 + GMCParticleRNGCalls; i++)
    trash = mt_random();

  typedef ActiveParticleType_GMCParticle ap;
  AttributeVector &ah = ap::AttributeHandlers;
  ActiveParticleType::SetupBaseParticleAttributes(ah);
  SetupGMCParticleAttributes(ah);

  return SUCCESS;
}

#define CII 9.74e5 /* Sound speed in ionized gas, assuming T = 7000 K and mu = 0.61 */

int ActiveParticleType_GMCParticle::EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data) 
{
  if (ActiveParticleType_AccretingParticle::EvaluateFormation(thisgrid_orig, data) == FAIL)
    return FAIL;

  if (data.NumberOfNewParticles > 0) {
    
    for (int i = 0; i < data.NumberOfNewParticles; i++) {
      ActiveParticleType_GMCParticle* np = 
	new ActiveParticleType_GMCParticle
	(static_cast<ActiveParticleType_AccretingParticle*>(data.NewParticles[i]),data);
      ActiveParticleType* temp = data.NewParticles[i];
      data.NewParticles[i] = np;
      delete temp;
    }
  }

  return SUCCESS;
}

void ActiveParticleType_GMCParticle::SetupGMCParticleAttributes(
		std::vector<ParticleAttributeHandler*> &handlers)
{

  typedef ActiveParticleType_GMCParticle ap;
  
  handlers.push_back(new Handler<ap, float, &ap::M0>("gmcevol_M0"));
  handlers.push_back(new Handler<ap, float, &ap::R0>("gmcevol_R0"));
  handlers.push_back(new Handler<ap, float, &ap::sigma0>("gmcevol_sigma0"));
  
  handlers.push_back(new Handler<ap, float, &ap::R>("gmcevol_R"));
  handlers.push_back(new Handler<ap, float, &ap::Rdot>("gmcevol_Rdot"));
  handlers.push_back(new Handler<ap, float, &ap::M>("gmcevol_M"));
  handlers.push_back(new Handler<ap, float, &ap::Mstar>("gmcevol_Mstar"));
  handlers.push_back(new Handler<ap, float, &ap::MstarRemain>("gmcevol_MstarRemain"));
  handlers.push_back(new Handler<ap, float, &ap::Massoc>("gmcevol_Massoc"));
  handlers.push_back(new Handler<ap, float, &ap::MdotStar>("gmcevol_MdotStar"));
  handlers.push_back(new Handler<ap, float, &ap::MdotHII>("gmcevol_MdotHII"));
  handlers.push_back(new Handler<ap, float, &ap::MstarRemain>("gmcevol_MstarRemain"));
  handlers.push_back(new Handler<ap, float, &ap::sigma>("gmcevol_sigma"));
  handlers.push_back(new Handler<ap, float, &ap::tau>("gmcevol_tau"));
  handlers.push_back(new Handler<ap, float, &ap::dtau>("gmcevol_dtau"));
  handlers.push_back(new Handler<ap, float, &ap::Eacc>("gmcevol_Eacc"));
  handlers.push_back(new Handler<ap, float, &ap::ReservoirRatio>("gmcevol_f"));

  handlers.push_back(new Handler<ap, int, &ap::nHIIreg>("gmcevol_nhIIreg"));
}

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, 
	     FLOAT Time);

int ActiveParticleType_GMCParticle::AdvanceCloudModel(FLOAT Time)
{
  // Get enzo units
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits, MassUnits;
    
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits,
	   &VelocityUnits, &MassUnits, Time); 


  // Declare derived parameters;
  float Mach0 = sigma0 / cs;
  float avir0 = 5*sigma0*sigma0*R0 / (GravConst*M0);
  float t0    = R0 / sigma0;
  float etaP  = 4 * pi * R0 * R0 * R0 * Pamb / 
    (aI * M0 * sigma0 * sigma0);
  // My use of the global units variables will not work for cosmological simulations
  float MdotAcc     = AccretionRate / (MassUnits*M0) * (TimeUnits*t0);
  float rho   = (M*M0) / (ReservoirRatio*4./3.*pi*(R*R*R*R0*R0*R0));
  float tff   = SQRT(3*pi / (32*GravConst*rho));
  float tauff = tff / t0;
  float zeta  = MdotAcc*tauff/M;

  float aprime, xi, chi, gamma, etaI;

  if (zeta > 1000) 
    ENZO_FAIL("GMCParticle: Zeta is greater than 1000!\n");
  if (zeta >= 0.001) {
    int zetaindex = (log10(zeta)+3)*10; // Converting to the log space in the zeta lookup table
    double zetamin = accTable.zetaLook[zetaindex];
    double zetamax = accTable.zetaLook[zetaindex+1];

    aprime = logInterpolate(zetaindex,accTable.aprimeLook, zeta, zetamin, zetamax);
    ReservoirRatio = logInterpolate(zetaindex,accTable.fLook,zeta,zetamin,zetamax);
    xi = logInterpolate(zetaindex,accTable.xiLook,zeta,zetamin,zetamax);
    chi = logInterpolate(zetaindex,accTable.chiLook,zeta,zetamin,zetamax);
    gamma = logInterpolate(zetaindex,accTable.gammaLook,zeta,zetamin,zetamax);

    etaI = (10 * gamma) / (ReservoirRatio * aI * avir0);
  }
  else {
    aprime = a;
    ReservoirRatio = 1.0;
    xi = 1.11484; // Analytic value in the limit M_res = 0
    chi = 0.0;
    etaI = 10/(aI * avir0);
  }

  float etaG  = 3 * aprime * (1.0 - etaB*etaB) / (aI*avir0);
  float etaA  = (5.0 - krho) / (4.0 - krho) * sqrt(xi * xi * 10.0 / (avir0 * ReservoirRatio));

  // Setup temporary state variables;
  float Mdot, sigmadot, Rddot, Mddot, Rdot, Tco, dtauSave, 
    sigmadot_noacc, Rddot_noacc, Rdot_noacc, Ecl, Ecl_noacc, R_noacc,
    M_noacc, sigma_noacc, sigmaISM;

  Mdot = MdotAcc + MdotHII + MdotStar;
  Rdot = sigmadot = Rddot = Mddot = Tco = dtauSave = sigmadot_noacc = 
    Rddot_noacc = Rdot_noacc = Ecl_noacc = R_noacc = M_noacc =
    sigma_noacc = 0;

  Ecl = 0.5*aI*M*Rdot*Rdot + 2.4*M*sigma*sigma + 1.5*M*Mach0 -
    5.0/avir0*(0.6*aprime*(1 - etaB*etaB) - chi) * M/R/R;

  return SUCCESS;
}

namespace {
  ActiveParticleType_info *GMCParticleInfo = 
    register_ptype <ActiveParticleType_GMCParticle> 
    ("GMCParticle");
}

std::vector<ParticleAttributeHandler*>
    ActiveParticleType_GMCParticle::AttributeHandlers;

