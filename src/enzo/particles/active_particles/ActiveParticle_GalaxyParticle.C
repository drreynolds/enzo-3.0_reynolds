/***********************************************************************
/
/  GALAXY PARTICLE TYPE
/
/  written by: Stephen Skory
/  date:       August, 2012
/
/  PURPOSE:
/
************************************************************************/

#include "ActiveParticle_GalaxyParticle.h"

/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

class GalaxyParticleGrid : private grid {
  friend class ActiveParticleType_GalaxyParticle;
};

/* Note that we only refer to GalaxyParticleGrid here. 
 * Given a grid object, we static cast to get this:
 *
 *    GalaxyParticleGrid *thisgrid =
 *      static_cast<GalaxyParticleGrid *>(thisgrid_orig); */

int ActiveParticleType_GalaxyParticle::InitializeParticleType(void)
{

  ActiveParticleType::SetupBaseParticleAttributes(
    ActiveParticleType_GalaxyParticle::AttributeHandlers);

  return SUCCESS;
}

int ActiveParticleType_GalaxyParticle::EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  GalaxyParticleGrid *thisgrid =
    static_cast<GalaxyParticleGrid *>(thisgrid_orig);
  fprintf(stderr, "Checking formation of galaxy particles.\n");
  return 0;
}

int ActiveParticleType_GalaxyParticle::EvaluateFeedback
(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  return SUCCESS;
}

void ActiveParticleType_GalaxyParticle::DescribeSupplementalData(ActiveParticleFormationDataFlags &flags)
{
  flags.DarkMatterDensity = true;
}

int ActiveParticleType_GalaxyParticle::SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level,
							int TopGridDims[], int GalaxyParticleID)
{
  return SUCCESS;
}

namespace {
  ActiveParticleType_info *GalaxyParticleInfo = 
    register_ptype <ActiveParticleType_GalaxyParticle> 
    ("GalaxyParticle");
}

std::vector<ParticleAttributeHandler*>
  ActiveParticleType_GalaxyParticle::AttributeHandlers;
