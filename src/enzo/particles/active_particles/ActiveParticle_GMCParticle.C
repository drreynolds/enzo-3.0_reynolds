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

namespace {
  ActiveParticleType_info *GMCParticleInfo = 
    register_ptype <ActiveParticleType_GMCParticle, GMCParticleBufferHandler> 
    ("GMCParticle");
}
