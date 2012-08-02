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

namespace {
  ActiveParticleType_info *GMCParticleInfo = 
    register_ptype <ActiveParticleType_GMCParticle, GMCParticleBufferHandler> 
    ("GMCParticle");
}
