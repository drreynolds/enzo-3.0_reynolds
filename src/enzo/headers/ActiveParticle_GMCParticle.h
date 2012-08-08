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
  static void SetupGMCParticleAttributes(std::vector<ParticleAttributeHandler*> &handlers);
  template <class active_particle_class>
  static int BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				 int ThisLevel, int TotalStarParticleCountPrevious[],
				 int GMCParticleID);
  template <class active_particle_class>
  static int AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int TotalStarParticleCountPrevious[],
				int GMCParticleID);


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
  float R, Rdot, M, Mdot, MdotAcc, MdotStar, MdotHII, MstarRemain, sigma, 
    tau, dtau, Mstar, Massoc, Eacc;

  int nHIIreg;

  /* Attribute handler instance */
  static std::vector<ParticleAttributeHandler *> AttributeHandlers;
};

// Static const member variables must be set outside the class definintion.  

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
  float r, rdot, Tco, mdot, mddot; 
  
  /* Externally used quantities, stored in                                                        
     dimensionless numbers */
  float EHII; /* The energy added to the cloud by the HII region */
  int phase;
  int breakoutFlag;
};

template <class active_particle_class>
int ActiveParticleType_GMCParticle::BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
						      int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
						      int ThisLevel, int TotalStarParticleCountPrevious[],
						      int GMCParticleID)
{
  if (ActiveParticleType_AccretingParticle::BeforeEvolveLevel<active_particle_class>
      (Grids, MetaData, NumberOfGrids, LevelArray,ThisLevel, TotalStarParticleCountPrevious,GMCParticleID) == FAIL)
    ENZO_FAIL("AccretingParticle BeforeEvolveLevel failed!");

  return SUCCESS;
}

template <class active_particle_class>
int ActiveParticleType_GMCParticle::AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
						     int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
						     int ThisLevel, int TotalStarParticleCountPrevious[],
						     int GMCParticleID)
{
  if (ActiveParticleType_AccretingParticle::AfterEvolveLevel<active_particle_class>
      (Grids, MetaData, NumberOfGrids, LevelArray,ThisLevel, TotalStarParticleCountPrevious,GMCParticleID) == FAIL)
    ENZO_FAIL("AccretingParticle AfterEvolveLevel failed!");

  return SUCCESS;
}
