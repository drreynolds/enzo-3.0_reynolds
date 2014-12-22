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
  ActiveParticleType_GMCParticle(ActiveParticleType_AccretingParticle* ap,
				 ActiveParticleFormationData &data);
  
  ActiveParticleType_GMCParticle(ActiveParticleType_GMCParticle* part);
  
  // static member functions
  static int InitializeParticleType();
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static int WriteToOutput(ActiveParticleType **these_particles, int n, int GridRank, hid_t group_id);
  static int ReadFromOutput(ActiveParticleType **&particles_to_read, int &n, int GridRank, hid_t group_id);
  template <class active_particle_class>
  static int BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
			       int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			       int ThisLevel, bool CallEvolvePhotons,
			       int TotalStarParticleCountPrevious[],
			       int GMCParticleID);
  template <class active_particle_class>
  static int AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int TotalStarParticleCountPrevious[],
				int GMCParticleID);
  template <class active_particle_class>
    static int DepositMass(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int GalaxyParticleID) {return SUCCESS; };

  // instance member functions
  int CalculateDerivedParameters();
  int UpdateDerivedParameters();

  // See Goldbaum et al. 2011 Table 1
  static const float krho;
  static const float etaB;
  static const float cs;
  static const float etaV;
  static const float phiIn;
  static const float aVgcm2;
  static const float aVmin;
  static const float HIIEff;
  static const float etaIn;
  static const float phicorr;
  static const float Pamb;

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

// Static const member variables must be set outside the class definintion.  
// C++11 has constepr for this.

const float ActiveParticleType_GMCParticle::krho    = 1.0;
const float ActiveParticleType_GMCParticle::etaB    = 0.5;
const float ActiveParticleType_GMCParticle::cs      = 0.19;
const float ActiveParticleType_GMCParticle::etaV    = 1.2;
const float ActiveParticleType_GMCParticle::phiIn   = 1.0;
const float ActiveParticleType_GMCParticle::aVgcm2  = 214.3; 
// Av corresponding to a gas column density of 1 g cm^-2, assuming solar metallicity
const float ActiveParticleType_GMCParticle::aVmin   = 1.4;
const float ActiveParticleType_GMCParticle::HIIEff  = 2.0;
const float ActiveParticleType_GMCParticle::etaIn   = 3.0;
const float ActiveParticleType_GMCParticle::phicorr = 0.75; // See Goldbaum et al. 2011 section 3
const float ActiveParticleType_GMCParticle::Pamb    = 3e4 * P_TO_CGS;

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


template <class active_particle_class>
int ActiveParticleType_GMCParticle::BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
						      int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
						      int ThisLevel, bool CallEvolvePhotons,
						      int TotalStarParticleCountPrevious[],
						      int GMCParticleID)
{
  if (ActiveParticleType_AccretingParticle::BeforeEvolveLevel<ActiveParticleType_GMCParticle>
      (Grids, MetaData, NumberOfGrids, LevelArray, ThisLevel, CallEvolvePhotons,
       TotalStarParticleCountPrevious,GMCParticleID) == FAIL)
    ENZO_FAIL("AccretingParticle BeforeEvolveLevel failed!");

  return SUCCESS;
}

template <class active_particle_class>
int ActiveParticleType_GMCParticle::AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
						     int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
						     int ThisLevel, int TotalStarParticleCountPrevious[],
						     int GMCParticleID)
{
  if (ActiveParticleType_AccretingParticle::AfterEvolveLevel<ActiveParticleType_GMCParticle>
      (Grids, MetaData, NumberOfGrids, LevelArray,ThisLevel, TotalStarParticleCountPrevious,GMCParticleID) == FAIL)
    ENZO_FAIL("AccretingParticle AfterEvolveLevel failed!");

  return SUCCESS;
}
