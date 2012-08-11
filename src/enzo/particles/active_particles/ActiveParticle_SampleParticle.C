/***********************************************************************
/
/  AN EXAMPLE ACTIVE PARTICLE TYPE
/
/  written by: Matthew Turk
/  date:       May, 2011
/
/  PURPOSE:
/
************************************************************************/

#include "ActiveParticle_SampleParticle.h"

/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

class SampleParticleGrid : private grid {
  friend class ActiveParticleType_SampleParticle;
};

/* Note that we only refer to SampleParticleGrid here. 
 * Given a grid object, we static cast to get this:
 *
 *    SampleParticleGrid *thisgrid =
 *      static_cast<SampleParticleGrid *>(thisgrid_orig); */

int ActiveParticleType_SampleParticle::InitializeParticleType(void)
{

  ActiveParticleType::SetupBaseParticleAttributes(
    ActiveParticleType_SampleParticle::AttributeHandlers);

  return SUCCESS;
}

int ActiveParticleType_SampleParticle::EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  SampleParticleGrid *thisgrid =
    static_cast<SampleParticleGrid *>(thisgrid_orig);
  fprintf(stderr, "Checking formation of sample particles.\n");
  return 0;
}

int ActiveParticleType_SampleParticle::EvaluateFeedback
(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  return SUCCESS;
}

void ActiveParticleType_SampleParticle::DescribeSupplementalData(ActiveParticleFormationDataFlags &flags)
{
  flags.DarkMatterDensity = true;
}

int ActiveParticleType_SampleParticle::SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level,
							int TopGridDims[], int SampleParticleID)
{
  return SUCCESS;
}

namespace {
  ActiveParticleType_info *SampleParticleInfo = 
    register_ptype <ActiveParticleType_SampleParticle> 
    ("SampleParticle");
}

std::vector<ParticleAttributeHandler*>
  ActiveParticleType_SampleParticle::AttributeHandlers;
