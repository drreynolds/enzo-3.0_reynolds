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
  //for (intn = 0; n < MAX_NUMBER_OF_HII_REGIONS; n++)
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
  
  float zeta, f, xi, aprime, chi, gamma;

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
  int nremove = 0;

  *HIIregEsc = 0;

  for (n=0; n<MAX_NUMBER_OF_HII_REGIONS; n++) {

    if (HIIregions[n].phase == -1)
      continue;

    /* First check if the his HII region is expanding slower thant he
       cloud velocity dispersion.  If so, remove it. */
    if (HIIregions[n].rdot <= (sigma*sqrt(HIIregions[n].r/R))) {
      HIIregions[n].phase = 3;
      nremove++;
    }
 
    /* Second check if the radius of this HII region exceeds the cloud radius */
    else if (HIIregions[n].r >= R) {
      
      /* Now compare the HII region expansion velocity to the loud
	 escape velocity and set the HII region disruption flag if it
	 is larger */
      
      if (HIIregions[n].rdot*sigma0 > SQRT(2*GravConst*M*M0 / (ReservoirRatio*R*R0))) {
	*HIIregEsc = true;
	HIIregions[n].phase = 2;
	return;
      }

      /* If we're here, the HII region cannot unbind the cloud.  Merge
	 the HII region either if it is a snowplow (phase 1) or if it
	 is driven but has a radius > 1.741 R */
      if ((HIIregions[n].phase == 1) || (HIIregions[n].r > 1.741*R)) { 
	HIIregions[n].phase = 3;
	nremove++;
      }
    }
  }
  
  /* Deal with phase 3 (to-be-merged) HII regions */
  if (nremove != 0) {
    for (n=0; n<MAX_NUMBER_OF_HII_REGIONS; n++) {

      /* Should we destroy this one? */

      if (HIIregions[n].phase == 3) {
	sigma = SQRT(sigma*sigma + 2./3.*HIIregions[n].Tco / M * 
		     SQRT(HIIregions[n].r / (phiIn*R)));
	HIIregions[n].phase = -1;
      }
    }
  }

}

float ranUniform() {
  unsigned_long_int random_int = mt_random();
  const int max_random = (1<<16);
  return (float) (random_int%max_random) / (float) (max_random);
}

float ranpowerlaw(float min, float max, float index) {
  // see http://mathworld.wolfram.com/RandomNumber.html
  float x = ranUniform();
  return POW((POW(max,index+1) - POW(min,index+1))*x + POW(min,index+1),POW(index+1,-1));
}

#define MMINASSOC 20.0

float Massocgen(float Mcl) {
  
  float Mmax = 0.1*Mcl/SolarMass;
  
  return (ranpowerlaw(MMINASSOC, Mmax, -2.0));
}

/* Lookup for ionizing luminosity s and main sequence lifetime tms
   versus m from fit of Parravano et al. (2003) */
float s49Lookup(float m) {

  if (m<5) {
    return(0.0);
  } else if (m<7) {
    return(2.23e-15*POW(m,11.5));
  } else if (m<12) {
    return(3.69e-13*POW(m,8.87));
  } else if (m<20) {
    return(4.8e-12*POW(m,7.85));
  } else if (m<30) {
    return(3.12e-8*POW(m,4.91));
  } else if (m<40) {
    return(2.8e-5*POW(m,2.91));
  } else if (m<60) {
    return(3.49e-4*POW(m,2.23));
  } else {
    return(2.39e-3*POW(m,1.76));
  }
}

float tmsmyrLookup(float m) {

  if (m<3) {
    return(7.65e3*POW(m,-2.8));
  } else if (m<6) {
    return(4.73e3*POW(m,-2.36));
  } else if (m<9) {
    return(2.76e3*POW(m,-2.06));
  } else if (m<12) {
    return(1.59e3*POW(m, -1.81));
  } else {
    return(7.6e2*POW(m,-1.57)+2.3);
  }
}



/*
  Lookup table for v band luminosity, e.g. Parravano et al. (2003).  
  Found using Lejeune spectral library and Geneva evolutionary tracks,
  see Lejeune & Schaerer 2001
 */
float LvLookup(float m) {
  if (m<5) {
    return(0.0);
  } else if (m<7) {
    return(6.431965*POW(m,2.127491));
  } else if (m<10) {
    return(5.978640*POW(m,2.165051));
  } else if (m<12) {
    return(6.261533*POW(m,2.144973));
  } else if (m<15) {
    return(5.437515*POW(m,2.201753));
  } else if (m<20) {
    return(7.649976*POW(m,2.075689));
  } else if (m<25) {
    return(15.810358*POW(m,1.833356));
  } else if (m<40) {
    return(15.057851*POW(m,1.848506));
  } else if (m<60) {
    return(0.175271*POW(m,3.055736));
  } else if (m<85) {
    return(9321.468750*POW(m,0.398047));
  } else {
    return(375.282410*POW(m,1.121127));
  }
}

/*  
  Lookup table for bolometric luminosity, e.g. Parravano et al. (2003).
 */

float LbolLookup(float m) {
  if (m<5) {
    return(0.0);
  } else if (m<7) {
    return(3.332853*POW(m,3.420702));
  } else if (m<10) {
    return(4.362311*POW(m,3.282376));
  } else if (m<12) {
    return(7.052795*POW(m,3.073731));
  } else if (m<15) {
    return(10.108081*POW(m,2.928888));
  } else if (m<20) {
    return(19.547123*POW(m,2.685352));
  } else if (m<25) {
    return(35.611240*POW(m,2.485123));
  } else if (m<40) {
    return(80.401276*POW(m,2.232125));
  } else if (m<60) {
    return(274.946289*POW(m,1.898816));
  } else if (m<85) {
    return(815.505310*POW(m,1.633271));
  } else {
    return(2082.599121*POW(m,1.422234));
  }
}

/* Routine to pick a massive star from the Kroupa (2001) IMF, which is 
   dN / dlog m \propto 
   m^0.7 (0.01 < m < 0.08)
   m^-0.8 (0.08 < m < 0.5)
   m^-1.7 (0.5 < m < 1)
   m^-1.3 (1 < m < 120)
*/

/* Probability of being in each interval of the Kroupa IMF */
#define WGT1 0.4967
#define WGT2 0.4360
#define WGT3 0.0426
#define WGT4 0.0247

float Mstargen() {

  float x, mmin, mmax, gamma;

  /* First decide which power-law interval to draw from. Note we use
     gamma 1 less than the exponent in Kroupa (2001), because we want
     to draw from dN / dm = (dN / dlog m) / m */
  x = ranUniform();
  if (x <= WGT1) {
    mmin = 0.01;
    mmax = 0.08;
    gamma = -0.3;
  } else if (x <= WGT1+WGT2) {
    mmin = 0.08;
    mmax = 0.5;
    gamma = -1.8;
  } else if (x <= WGT1+WGT2+WGT3) {
    mmin = 0.5;
    mmax = 1.0;
    gamma = -2.7;
  } else {
    mmin = 1.0;
    mmax = 120.0;
    gamma = -2.3;
  }

  /* Now draw from a power-law in m and return */
  return(ranpowerlaw(mmin, mmax, gamma));
}

#define LSUN      3.839e33
#define PC        3.09e18
#define YR        (365.25*24.*3600.)
#define MYR       (1.0e6*YR)
#define MSTARMEAN 0.21

struct sort_pfirst {
  bool operator()(const std::pair<float,float> &left, const std::pair<float,float> &right) {
    return left.first < right.first;
  }
};

struct pair_add {
  float operator()(float lhs, const std::pair<float, float>& x) {
    return lhs + x.first;
  }
};

void ActiveParticleType_GMCParticle::CreateHIIregions() {
  float Massocremain;
  int n;
  float m, Lvtot = 0, Lboltot = 0, s49tot = 0, s49sum = 0, tmscut;

  /* Increase the mass of stars available to go into a cluster, stored
     in solar masses. */
  MstarRemain += (-MdotStar)*M0*dtau/SolarMass;

  if (Massoc == 0.0) 
    Massoc = Massocgen(M*M0);

  while (MstarRemain >= 0.5*Massoc) {

    /* Remove mass in stars from the 'waiting to be born' mass */
    MstarRemain -= Massoc;

    /* Allocate memory to hold s49 and tms lists */
    std::vector<std::pair<float,float> > starData;

    /* Generate stars in the association until we have enough to add
       up to the association mass */
    Massocremain = Massoc;
    while (Massocremain > 0) {
      m = Mstargen();

      /* Subtract mass from remaining association mass */
      if (Massocremain > m) {
	Massocremain -= m;
      } else {
	m = Massocremain;
	Massocremain = 0.0;
      }
	
      /* If star is massive enough to bother, store ionizing
	 luminosity, lifetime, Lbol and Lv */
      if (m > 5) {
	starData.push_back(std::make_pair(s49Lookup(m),tmsmyrLookup(m)*MYR));
	Lvtot += LvLookup(m);
	Lboltot += LbolLookup(m);
      }

    }

    /* Generate next association mass */
    Massoc = Massocgen(M*M0);
    
    /* Sort stars by ionizing luminosity */
    std::sort(starData.begin(), starData.end(), sort_pfirst());
    s49tot = std::accumulate(starData.begin(), starData.end(), 0, pair_add());
    
    int starCtr = starData.size();

    for (s49sum=0; starCtr>=0; starCtr--) {
      s49sum += starData[starCtr].first;
      if (s49sum >= s49tot/2.0) break;
    }

    tmscut = starData[starCtr].second;
    
    starData.clear();
    
    /* Set HII region properties */
    
    bool set = false;
    
    for (int n = 0; n < MAX_NUMBER_OF_HII_REGIONS; n++) {
      if (HIIregions[n].phase == -1) {
	HIIregions[n].r = HIIregions[n].rdot = HIIregions[n].mdot = 
	  HIIregions[n].Tco = 0.0;
	
	HIIregions[n].t0 = tau*sigma0/R0;
	HIIregions[n].rms = tau*sigma0/R0+tmscut;
	HIIregions[n].s49 = s49tot;
	HIIregions[n].Lv = Lvtot;
	HIIregions[n].L39 = Lboltot*LSUN/1e39;
	HIIregions[n].mcl = M*M0;
	HIIregions[n].nh22 = M/(pi*R*R) * M0 / (R0*R0) / (1.0e22*MU_CLOUD);
	HIIregions[n].tch = 6456 * YR * sqrt(HIIregions[n].nh22/1.28) *
	  POW(HIIregions[n].L39,5.0/2.0) * POW(s49tot,-3.0/2.0);
	HIIregions[n].rch = .1*PC/(HIIregions[n].s49) *
	  POW((HIIregions[n].L39),2);
	HIIregions[n].phase = 0;
	HIIregions[n].rms = 0;
	HIIregions[n].rdotms = 0;
	HIIregions[n].Tcoms = 0;
	set = true;
	break;
      }
    }
    if (!set)
      ENZO_FAIL("Exceeded MAX_NUMBER_OF_HII_REGIONS!");
  }
}


namespace {
  ActiveParticleType_info *GMCParticleInfo = 
    register_ptype <ActiveParticleType_GMCParticle> 
    ("GMCParticle");
}

std::vector<ParticleAttributeHandler*>
    ActiveParticleType_GMCParticle::AttributeHandlers;

