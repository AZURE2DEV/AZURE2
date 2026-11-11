// The set-up of a yield on a layer of finite thickness (TargetEffect::IsLayerYield; docs/upgrade/energy_grid.md
// section 5): the refusals that need the stopping evaluated on the energies of the data, with the codes and messages of
// experiment/'s LayerChecks.h.
#ifndef TARGETLAYER_H
#define TARGETLAYER_H

#include <string>

class Config;
class TargetEffect;
class Unit;
class WorkingCopy;

/// A point whose layer has no positive thickness at its beam energy (the stopping there is not positive and finite):
/// prints E_LAYER_STOPPING and returns -1, the set-up's refusal.
int RefuseLayerAtBeamEnergy(const TargetEffect &effect, const Unit &point, const std::string &dataFile, const Config &configure);

/// After the sub-points are committed (WorkingCopy::SetUpTargetEffects): for every point of a yield on a layer, the
/// stopping is positive and finite at every node of its grid (else E_LAYER_STOPPING), and the beam leaves the layer
/// inside the grid at the largest density the working copy may hold, the layout density (else
/// E_LAYER_EXIT_BELOW_MODEL: the grid stops at TargetEffect::minIntegrationEnergy).  E_out is found as the kernel finds
/// it (LayerExit.h, from the kernel's top: the grid's first node without a spread, the point's energy with one).
/// Prints the first refusal and returns false.
bool CheckLayerYields(WorkingCopy &data, const Config &configure);

#endif
