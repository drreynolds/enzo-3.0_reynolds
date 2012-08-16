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
  for (int n = 0; n < MAX_NUMBER_OF_HII_REGIONS; n++)
    HIIregions[n].phase = -1;
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
  tau  = this->BirthTime*R0/sigma0;
  dtau = 1.0e-4;
  Eacc = 0.0;
  nHIIreg = 0;
  //for (int n = 0; n < MAX_NUMBER_OF_HII_REGIONS; n++)
  //  HIIregions[n].phase = -1;
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

  for (int n = 0; n < MAX_NUMBER_OF_HII_REGIONS; n++)
    HIIregions[n] = part->HIIregions[n];
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

  setup_HIIregion_output();

  return SUCCESS;
}

#define CII 9.74e5 /* Sound speed in ionized gas, assuming T = 7000 K and mu = 0.61 */
#define STEPMAX 0.001
#define DTINCRMAX 0.1
#define MINSTEP 1.0e-8
#define MU_CLOUD 2.34e-24

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
  handlers.push_back(new HIIregionHandler<ap>);
}

float ActiveParticleType_GMCParticle::sfr() 
{
  float rho, tff, tcr0, mach, avir, sfrff, sfrtot;

  /* Free-fall time */
  rho = (M*M0) / (ReservoirRatio*4./3.*PI*(R*R*R*R0*R0*R0));
  tff = sqrt(3*PI / (32*GravConst*rho));

  /* Crossing time at start of run -- sets time units */
  tcr0 = R0 / sigma0;

  /* Mach number and virial parameter */
  mach = (sigma*sigma0) / cs;
  avir = 5*(sigma*sigma*sigma0*sigma0) * (R*R0) /
    (GravConst * (M*M0));

  /* Star formation rate per free-fall time */
  /* See Krumholz et al. (2005)             */
  sfrff = 0.022*pow(avir/1.3, -0.68)*pow(mach/25.0, -0.33);

  /* Star formation rate */
  sfrtot = sfrff * M / (tff/tcr0);
  return(sfrtot);
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
  float Mach0   = sigma0 / cs;
  float avir0   = 5*sigma0*sigma0*R0 / (GravConst*M0);
  float t0      = R0 / sigma0;
  float etaP    = 4 * pi * R0 * R0 * R0 * Pamb / 
    (aI * M0 * sigma0 * sigma0);
  float etaE    = (5.0 - krho) / (4.0 - krho) * 2 * CII / sigma0;
  float MdotAcc = AccretionRate / (MassUnits*M0) * (TimeUnits*t0);

  float rho, tauMax, tff, tauff, zeta, aprime, xi, chi, gamma, etaI,
    zetamin, zetamax, f;
  int zetaindex;
  bool Converged = false;
  bool HIIregEsc = false;

  while (!Converged) {
    rho     = (M*M0) / (ReservoirRatio*4./3.*pi*(R*R*R*R0*R0*R0));
    tauMax  = Time * TimeUnits / t0;
    tff     = SQRT(3*pi / (32*GravConst*rho));
    tauff   = tff / t0;
    zeta    = MdotAcc*tauff/M;
    if (zeta > 5) 
      ENZO_FAIL("GMCParticle: Zeta is greater than 5!\n");
    if (zeta >= 0.001) {
      zetaindex = floor((log10(zeta)+3)*10); // Converting to the log space of the zeta lookup table
      zetamin = accTable.zetaLook[zetaindex];
      zetamax = accTable.zetaLook[zetaindex+1];

      aprime = logInterpolate(zetaindex,accTable.aprimeLook, zeta, zetamin, zetamax);
      f = logInterpolate(zetaindex,accTable.fLook,zeta,zetamin,zetamax);
      xi = logInterpolate(zetaindex,accTable.xiLook,zeta,zetamin,zetamax);
      chi = logInterpolate(zetaindex,accTable.chiLook,zeta,zetamin,zetamax);
      gamma = logInterpolate(zetaindex,accTable.gammaLook,zeta,zetamin,zetamax);
      
      etaI = (10 * gamma) / (f * aI * avir0);
      if (abs((f - ReservoirRatio)/ReservoirRatio > 1e-5))
	Converged=true;
      ReservoirRatio = f;
    }
    else {
      aprime = a;
      ReservoirRatio = 1.0;
      xi = 1.11484; // Analytic value in the limit M_res = 0
      chi = 0.0;
      etaI = 10/(aI * avir0);
    }
  }

  float etaG  = 3 * aprime * (1.0 - etaB*etaB) / (aI*avir0);
  float etaA  = (5.0 - krho) / (4.0 - krho) * sqrt(xi * xi * 10.0 / (avir0 * ReservoirRatio));

  // Setup temporary state variables;
  float Mdot, sigmadot, Rddot, Rdot, Tco, dtauSave, 
    sigmadot_noacc, Rddot_noacc, Rdot_noacc, Ecl, Ecl_noacc, R_noacc,
    M_noacc, sigma_noacc, sigmaISM, Lambda;

  int i;
  int dtauFloor = 0;

  Rdot = sigmadot = Rddot = Tco = dtauSave = sigmadot_noacc = 
    Rddot_noacc = Rdot_noacc = Ecl_noacc = R_noacc = M_noacc =
    sigma_noacc = 0;

  Ecl = 0.5*aI*M*Rdot*Rdot + 2.4*M*sigma*sigma + 1.5*M*Mach0 -
    5.0/avir0*(0.6*aprime*(1 - etaB*etaB) - chi) * M/R/R;

  /* Main update loop */
  while (tau < tauMax) {
    /* compute current star formation rate */
    MdotStar = -this->sfr();
    
    /* compute current mass loss rate from sufficiently large HII regions */
    for (i=0, MdotHII=0, Tco = 0; i<nHIIreg; i++) {
      MdotHII -= HIIregions[i].mdot;
      Tco += HIIregions[i].Tco;
    }
    
    Mdot = MdotAcc + MdotHII + MdotStar;

    Lambda = (etaV / phiIn) * M * sigma * sigma * sigma / R;
    
    Rddot = (3.9*sigma*sigma + 3.0/(Mach0*Mach0)) / (aI*R) 
      - etaG * M / (R*R)
      - etaP * R*R/M
      - MdotAcc * Rdot / M
      + etaE * MdotHII / M
      - etaA *MdotAcc / SQRT(M*R);
    
    sigmadot = aI / 4.8 * 
      (
       - Rdot * Rddot / sigma
       - etaG * M * Rdot / (R * R * sigma)
       - etaP * R * R * Rdot / (M * sigma)
       + etaE * MdotHII * Rdot / (M * sigma)
       - etaA * MdotAcc * Rdot / (SQRT(M*R)*sigma)
       - MdotAcc * Rdot * Rdot / (M * sigma)
       - (3 - 1.5*phicorr) * MdotAcc * sigma / (aI * M)
       + (3*phicorr)/(2*aI) * MdotAcc/(M*sigma) * (sigmaISM*sigmaISM)/(sigma0*sigma0)
       + etaI * phicorr * MdotAcc / (R*sigma)
       - 1.0/aI * Lambda / (M * sigma)
       );
    
    /* Make sure subcycle step isn't too big.  If so, scale it back */
    if (tau+dtau > tauMax) dtau = tauMax - tau;
    bool dtauOk = 1;
    while((fabs(Rdot*dtau) > STEPMAX*R) ||
	  (fabs(0.5*Rddot*dtau*dtau) > STEPMAX*R) ||
	  (fabs(Mdot*dtau) > STEPMAX*M) ||
	  (fabs(sigmadot*dtau) > STEPMAX*sigma)) {
      dtau /= 2.0;
      dtauOk = 0;
    }
    
    if(dtau/tau < MINSTEP) 
      if (R < 0.01)
	return -3;
      else
	return -4;
    
    /* Update Cloud quantities */

    Rddot_noacc = (3.9*sigma*sigma + 3.0/(Mach0*Mach0)) / (aI*R)
      - etaG * M/(R*R)
      - etaP * R*R/M
      + etaE * (MdotHII) / M;

    Rdot_noacc = Rdot + Rddot_noacc*dtau;

    sigmadot_noacc = aI / 4.8 *
      (
       - Rdot_noacc * Rddot_noacc / sigma
       - etaG * M * Rdot_noacc / (R * R * sigma)
       - etaP * R * R * Rdot_noacc / (M * sigma)
       + etaE * MdotHII * Rdot_noacc / (M * sigma)
       - 1.0/aI * Lambda / (M * sigma)
       	   );

    R_noacc = R + Rdot_noacc*dtau;
    M_noacc = M + MdotHII*dtau;
    sigma_noacc = sigma + sigmadot_noacc*dtau;
    
    R += Rdot*dtau;
    Rdot += Rddot*dtau;
    M += Mdot*dtau;
    sigma += sigmadot*dtau;
    Mstar += -MdotStar*dtau;
    tau += dtau;
      
    Ecl = 0.5*aI*M*Rdot*Rdot + 2.4*M*sigma*sigma + 1.5*M/(Mach0*Mach0) - 
      5.0/avir0*(0.6*aprime*(1-etaB*etaB) - chi)*M*M/R;

    Ecl_noacc = 0.5*aI*M_noacc*Rdot_noacc*Rdot_noacc 
      + 2.4*M_noacc*sigma_noacc*sigma_noacc + 1.5*M_noacc/(Mach0*Mach0) - 
      5.0/avir0*(0.6*aprime*(1-etaB*etaB) - chi)*M_noacc*M_noacc/R_noacc;

    Eacc += Ecl - Ecl_noacc;

    /* Exit if cloud surface density has dropped to the point where it
       should dissociate */
    double aV;
    
    aV = (M*M0)/(pi*R*R*R0*R0)*aVgcm2;
    if (aV < aVmin)
      return -1;

    /* Update state of existing HII regions */
    this->UpdateHIIregions(&HIIregEsc);
    
    /* Exit if an HII region has encompassed the whole cloud */
    if (HIIregEsc)
      return -2;

    /* Check to see if any new HII regions should appear */
    this->CreateHIIregions();
    
    /* Increase time step if we can */
    if (dtauOk) dtau *= (1.0+DTINCRMAX);
        
  }
  return SUCCESS;
}

void ActiveParticleType_GMCParticle::UpdateHIIregions(bool* HIIregEsc) {
  int n;
  double t, p, pdotgas;

  int ActiveCount = 0;

  /* Get physical time */
  t = tau*R0/sigma0;

  for (n=0; n < MAX_NUMBER_OF_HII_REGIONS; n++) {
    
    if (HIIregions[n].phase == -1)
      continue;

    ActiveCount++;

    /* Update based on phase */
    if (HIIregions[n].phase == 0) {

      /* Phase 0, HII region is still driven */

      /* Compute new values of Tco, r, r', and M', first in physical units */
      float tauRad = (t - HIIregions[n].t0)/HIIregions[n].tch;
      int tauIndex = floor(125*log10(tauRad)+625);
      if (tauIndex < 0) tauIndex = 0;
      float tauMin = radSolTable.tauLook[tauIndex];
      float tauMax = radSolTable.tauLook[tauIndex+1];
      float xShell = logInterpolate(tauIndex,radSolTable.xShellLook, tauRad, tauMin, tauMax);
      float xPrimeShell = logInterpolate(tauIndex,radSolTable.xPrimeShellLook, tauRad, tauMin, tauMax);

      HIIregions[n].r = HIIregions[n].rch * xShell;
      HIIregions[n].rdot = HIIregions[n].rch / HIIregions[n].tch * xPrimeShell;

      p = pi/2 * HIIregions[n].nh22*MU_CLOUD*1e22 * POW(HIIregions[n].r,2) * HIIregions[n].rdot;

      pdotgas = 3.33564095e28 * HIIregions[n].L39 * sqrt(xShell);

      HIIregions[n].mdot = pdotgas / (2*CII);
      HIIregions[n].Tco = p*HIIregions[n].rdot/2.0;

      /* Check if we should convert this HII region to phase 1 
	 and store transition values */

      if (t > HIIregions[n].tms) {
	HIIregions[n].phase = 1;
	HIIregions[n].rms = HIIregions[n].r;
	HIIregions[n].rdotms = HIIregions[n].rdot;
	HIIregions[n].Tcoms = HIIregions[n].Tco;
	HIIregions[n].mdot = 0;
	HIIregions[n].pms = p;
      }

    }
    else if (HIIregions[n].phase == 1) {
      
      /* Phase 1, undriven snowplow expansion */
      
      /* Compute new values of Tco, r, r', and M' in physical units */
      HIIregions[n].r = HIIregions[n].rms * 
	POW(3*HIIregions[n].rdotms/HIIregions[n].rms * (t - HIIregions[n].tms) + 1, 1./3.);
      HIIregions[n].rdot = HIIregions[n].rdotms *
	POW(3*HIIregions[n].rdotms/HIIregions[n].rms * (t - HIIregions[n].tms) + 1, -2./3.);
      HIIregions[n].Tco = HIIregions[n].pms*HIIregions[n].rdot/2.0;
      
    }
    else
      ENZO_FAIL("HII region phase unrecognized");

    /* Convert state data to dimensionless gmcevol units */
    HIIregions[n].Tco /= M0*sigma0*sigma0;
    HIIregions[n].mdot *= R0/(sigma0*M0);
    HIIregions[n].r /= R0;
    HIIregions[n].rdot /= sigma0;

  }

  if (ActiveCount != nHIIreg)
    ENZO_FAIL("HII region accounting is inconsistent!");

  /* Turn off HII regions that we no longer need to track */
  int nremove;
  double p_sh;

  *HIIregEsc = 0;

  

}

void ActiveParticleType_GMCParticle::CreateHIIregions() {


}

namespace {
  ActiveParticleType_info *GMCParticleInfo = 
    register_ptype <ActiveParticleType_GMCParticle> 
    ("GMCParticle");
}

std::vector<ParticleAttributeHandler*>
    ActiveParticleType_GMCParticle::AttributeHandlers;

