/***********************************************************************
/
/ GMC Particle
/
************************************************************************/

#include "ActiveParticle_AccretingParticle.h"
#include "StarParticleData.h"

#define NUMZETA 62
#define NUMTAU  1402
#define P_TO_CGS 1.381e-16  /* This is k_b*K*cm^-3 in dyn cm^-2 */

class GMCParticleBufferHandler : public AccretingParticleBufferHandler {};  

class ActiveParticleType_GMCParticle : public ActiveParticleType_AccretingParticle {
public:
  // constructors
  ActiveParticleType_GMCParticle(void);

  // static member functions
  static int InitializeParticleType();
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static int WriteToOutput(ActiveParticleType **these_particles, int n, int GridRank, hid_t group_id);
  static int ReadFromOutput(ActiveParticleType **&particles_to_read, int &n, int GridRank, hid_t group_id);

  // instance member functions
  int CalculateDerivedParameters();
  int UpdateDerivedParameters();

  // See Goldbaum et al. 2011 Table 1
  static const float krho    = 1.0;
  static const float etaB    = 0.5;
  static const float cs      = 0.19;
  static const float etaV    = 1.2;
  static const float phiIn   = 1.0;
  static const float aVgcm2  = 214.3; // Av corresponding to a gas column density of 1 g cm^-2, assuming solar metallicity
  static const float aVmin   = 1.4;
  static const float HIIEff  = 2.0;
  static const float etaIn   = 3.0;
  static const float phicorr = 0.75; // See Goldbaum et al. 2011 section 3
  static const float Pamb    = 3e4 * P_TO_CGS;

  // Scaling factors to physical units.
  float M0, R0, sigma0;

  // State data (in gmcevol units)
  float R, M, sigma, Rdot, Mdot, MdotAcc, sigmadot, sigmadotAcc, 
    Rddot, Mddot, tau, dtau, Mstar, sigmaISM, Tco, MdotStar, MdotHII, 
    MddotHII, Lambda, Gamma, Massoc, MstarRemain, dtauSave,sigmadot_noacc,
    Rddot_noacc, Rdot_noacc, R_noacc, M_noacc, sigma_noacc, Ecl, 
    Ecl_noacc, Eacc;

  /* Derived parameters */
  float aI, a, aprime, Mach0, avir0, etaG, etaP, etaE, etaA, etaI, t0, f, xi, chi, gamma;

  int nHIIreg, dtauOk, HIIregEsc, dissoc, dtauFloor;
};

class HIIregion {
public:

  /* Constants determined at HII region creation and stored in                                                                      
     dimensional numbers */
  float mcl, nh22, s49, t0, tms, rms, rdotms, Tcoms, pms, Lv, L39, tch, rch;
  float r, rdot, Tco, mdot, mddot; /* Externally used quantities, stored in                                                        
				       dimensionless numbers */
  float EHII; /* The energy added to the cloud by the HII region */
  int phase;
  int breakoutFlag;
};

struct accTableStor {
  float* zetaLook;
  float* fLook;
  float* xiLook;
  float* aprimeLook;
  float* chiLook;
  float* gammaLook;
};

accTableStor accTable;

struct radSolTableStor {
  float* tauLook;
  float* xShellLook;
  float* xPrimeShellLook;
};

radSolTableStor radSolTable;
