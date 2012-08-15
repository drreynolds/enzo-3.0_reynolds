/***********************************************************************
/
/ GMC Particle
/
************************************************************************/

#include "ActiveParticle_AccretingParticle.h"
#include "StarParticleData.h"
#include "TopGridData.h"

#define NUMZETA 62
#define NUMTAU  1402
#define P_TO_CGS 1.381e-16  /* This is k_b*K*cm^-3 in dyn cm^-2 */
#define MAX_NUMBER_OF_HII_REGIONS 5000

class HIIregion {
public:

  /* Constants determined at HII region creation and stored in         
     dimensional numbers */
  float mcl, nh22, s49, t0, tms, rms, rdotms, Tcoms, pms, Lv, L39, tch, rch;

  /* HII region state variables, stored in                                                        
     dimensionless numbers */
  float r, rdot, Tco, mdot, mddot; 
  
  float EHII; /* The energy added to the cloud by the HII region */
  int phase;
  int breakoutFlag;
};

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

  // Does the grunt work of advancing gmcevol
  int AdvanceCloudModel(FLOAT Time);
  float sfr(void);

   // See Goldbaum et al. 2011 Table 1
  static const float krho;
  static const float aI;
  static const float a;
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
  float R, Rdot, M, MdotStar, MdotHII, MstarRemain, sigma, 
    tau, dtau, Mstar, Massoc, Eacc, ReservoirRatio;

  int nHIIreg;
  HIIregion HIIregions[MAX_NUMBER_OF_HII_REGIONS];

  /* Attribute handler instance */
  static std::vector<ParticleAttributeHandler *> AttributeHandlers;
};

// Static const member variables must be set outside the class definintion.  

const float ActiveParticleType_GMCParticle::krho    = 1.0;
const float ActiveParticleType_GMCParticle::aI      = (3 - krho) / (5 - krho);
const float ActiveParticleType_GMCParticle::a       = (15 - 5*krho) / (15 - 6*krho);
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

hid_t    HIIregion_tid;  /* HDF5 datatype identifier */

void setup_HIIregion_output() {
  HIIregion_tid = H5Tcreate(H5T_COMPOUND, sizeof(HIIregion));
  H5Tinsert(HIIregion_tid, "HIIregion_mcl", HOFFSET(HIIregion, mcl), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_nh22", HOFFSET(HIIregion, nh22), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_s49", HOFFSET(HIIregion, s49), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_t0", HOFFSET(HIIregion, t0), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_tms", HOFFSET(HIIregion, tms), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_rms", HOFFSET(HIIregion, rms), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_rdotms", HOFFSET(HIIregion, rdotms), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_Tcoms", HOFFSET(HIIregion, Tcoms), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_pms", HOFFSET(HIIregion, pms), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_Lv", HOFFSET(HIIregion, Lv), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_L39", HOFFSET(HIIregion, L39), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_tch", HOFFSET(HIIregion, tch), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_tch", HOFFSET(HIIregion, rch), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_r", HOFFSET(HIIregion, r), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_rdot", HOFFSET(HIIregion, rdot), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_Tco", HOFFSET(HIIregion, Tco), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_mdot", HOFFSET(HIIregion, mdot), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_mddot", HOFFSET(HIIregion, mddot), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_EHII", HOFFSET(HIIregion, EHII), H5T_NATIVE_FLOAT);
  H5Tinsert(HIIregion_tid, "HIIregion_phase", HOFFSET(HIIregion, phase), H5T_NATIVE_INT);
  H5Tinsert(HIIregion_tid, "HIIregion_breakoutFlag", HOFFSET(HIIregion, breakoutFlag), H5T_NATIVE_INT);
}

/* Create the memory datatype */

template <class APClass>
class HIIregionHandler : public ParticleAttributeHandler
{
 public:

  HIIregionHandler(int offset = 0) {
    this->name = "gmcevol_HIIregions";
    this->offset = offset;
    this->hdf5type = HIIregion_tid;
    this->element_size = 18*sizeof(float) + 2*sizeof(int);
  }

  void SetAttribute(char **buffer, ActiveParticleType *pp_) {
    APClass *pp = static_cast<APClass*>(pp_);
    float *pbf = (float *)(*buffer);
    pp->HIIregions[this->offset].mcl    = *(pbf++);
    pp->HIIregions[this->offset].nh22   = *(pbf++);
    pp->HIIregions[this->offset].s49    = *(pbf++);
    pp->HIIregions[this->offset].t0     = *(pbf++);
    pp->HIIregions[this->offset].tms    = *(pbf++);
    pp->HIIregions[this->offset].rdotms = *(pbf++);
    pp->HIIregions[this->offset].Tcoms  = *(pbf++);
    pp->HIIregions[this->offset].pms    = *(pbf++);
    pp->HIIregions[this->offset].Lv     = *(pbf++);
    pp->HIIregions[this->offset].L39    = *(pbf++);
    pp->HIIregions[this->offset].tch    = *(pbf++);
    pp->HIIregions[this->offset].rch    = *(pbf++);
    pp->HIIregions[this->offset].r      = *(pbf++);
    pp->HIIregions[this->offset].rdot   = *(pbf++);
    pp->HIIregions[this->offset].Tco    = *(pbf++);
    pp->HIIregions[this->offset].mdot   = *(pbf++);
    pp->HIIregions[this->offset].mddot  = *(pbf++);
    pp->HIIregions[this->offset].EHII   = *(pbf++);
    int *pbi = (int *)(pbf);
    pp->HIIregions[this->offset].phase  = *(pbi++);
    pp->HIIregions[this->offset].breakoutFlag   = *(pbi++);
    *buffer = (char *) pbi;
  }

  int GetAttribute(char **buffer, ActiveParticleType *pp_) {
    APClass *pp = static_cast<APClass*>(pp_);
    this->element_size = 0;
    float *pbf = (float *)(*buffer);
    *(pbf++) = (pp->HIIregions)[this->offset].mcl;
    *(pbf++) = (pp->HIIregions)[this->offset].nh22;
    *(pbf++) = (pp->HIIregions)[this->offset].s49;
    *(pbf++) = (pp->HIIregions)[this->offset].t0;
    *(pbf++) = (pp->HIIregions)[this->offset].tms;
    *(pbf++) = (pp->HIIregions)[this->offset].rdotms;
    *(pbf++) = (pp->HIIregions)[this->offset].Tcoms;
    *(pbf++) = (pp->HIIregions)[this->offset].pms;
    *(pbf++) = (pp->HIIregions)[this->offset].Lv;
    *(pbf++) = (pp->HIIregions)[this->offset].L39;
    *(pbf++) = (pp->HIIregions)[this->offset].tch;
    *(pbf++) = (pp->HIIregions)[this->offset].rch;
    *(pbf++) = (pp->HIIregions)[this->offset].r;
    *(pbf++) = (pp->HIIregions)[this->offset].rdot;
    *(pbf++) = (pp->HIIregions)[this->offset].Tco;
    *(pbf++) = (pp->HIIregions)[this->offset].mdot;
    *(pbf++) = (pp->HIIregions)[this->offset].mddot;
    *(pbf++) = (pp->HIIregions)[this->offset].EHII;
    int *pbi = (int *)(pbf);
    *(pbi++) = (pp->HIIregions)[this->offset].phase;
    *(pbi++) = (pp->HIIregions)[this->offset].breakoutFlag;
    return this->element_size;
  }

  void PrintAttribute(ActiveParticleType *pp_) {
    APClass *pp = static_cast<APClass*>(pp_);
    std::cout << this->name << std::endl;
    std::cout << (pp->HIIregions)[this->offset].mcl << std::endl;
    std::cout << (pp->HIIregions)[this->offset].nh22 << std::endl;
    std::cout << (pp->HIIregions)[this->offset].s49 << std::endl;
    std::cout << (pp->HIIregions)[this->offset].t0 << std::endl;
    std::cout << (pp->HIIregions)[this->offset].tms << std::endl;
    std::cout << (pp->HIIregions)[this->offset].rdotms << std::endl;
    std::cout << (pp->HIIregions)[this->offset].Tcoms << std::endl;
    std::cout << (pp->HIIregions)[this->offset].pms << std::endl;
    std::cout << (pp->HIIregions)[this->offset].Lv << std::endl;
    std::cout << (pp->HIIregions)[this->offset].L39 << std::endl;
    std::cout << (pp->HIIregions)[this->offset].tch << std::endl;
    std::cout << (pp->HIIregions)[this->offset].rch << std::endl;
    std::cout << (pp->HIIregions)[this->offset].r << std::endl;
    std::cout << (pp->HIIregions)[this->offset].rdot << std::endl;
    std::cout << (pp->HIIregions)[this->offset].Tco << std::endl;
    std::cout << (pp->HIIregions)[this->offset].mdot << std::endl;
    std::cout << (pp->HIIregions)[this->offset].mddot << std::endl;
    std::cout << (pp->HIIregions)[this->offset].EHII << std::endl;
    std::cout << (pp->HIIregions)[this->offset].phase << std::endl;
    std::cout << (pp->HIIregions)[this->offset].breakoutFlag << std::endl;
  }

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

double logInterpolate(int zetaindex, double interpArray[NUMZETA],
                      double zeta, double zetamin, double zetamax){
  double gridArea = log10(zetamax) - log10(zetamin);
  return( (interpArray[zetaindex]*(log10(zetamax) - log10(zeta))
           + interpArray[zetaindex+1]*(log10(zeta) - log10(zetamin)))/gridArea ) ;
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

  ActiveParticleType** ParticleList = NULL;
  int nParticles;

  ParticleList = ActiveParticleFindAll(LevelArray, &nParticles, GMCParticleID);

  /* Return if there are no GMC particles */

  if (nParticles == 0)
    return SUCCESS;

  /* Advance the GMC model for each particle */
  
  active_particle_class *GMCParticle = NULL;

  for (int i = 0; i < nParticles; i++) {
    GMCParticle = static_cast<active_particle_class*>(ParticleList[i]);
    GMCParticle->AdvanceCloudModel(MetaData->Time);
  }      
      
  return SUCCESS;
}

