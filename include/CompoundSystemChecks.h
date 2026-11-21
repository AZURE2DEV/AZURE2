#ifndef COMPOUNDSYSTEM_CHECKS_H
#define COMPOUNDSYSTEM_CHECKS_H

#include <vector>

class PPair;
class CompoundSystem;
class LevelInputs;
class ObservedLevelTable;
class Config;

/*!
 * The checks of a compound nucleus and their messages: what is refused at load, and the warnings about channels and
 * levels the calculation handles but the user should know of.  Each function reads the nucleus through const access
 * and writes to the configuration's output stream only.
 */

/// Refuses (returns -1) a channel with L+S odd on an identical pair with spin, after naming each; 0 otherwise.
int RefuseForbiddenIdenticalChannels(const CompoundSystem &compound, const Config &configure);
/// Warns of a channel with L+S odd on an identical spin-0 pair (it is kept: such a model is refused only with spin).
void WarnForbiddenSpinlessChannels(const CompoundSystem &compound, const Config &configure);
/// Warns of a particle channel with l > Config::maxLOrder / 2, whose Legendre terms above maxLOrder are dropped.
void WarnTruncatedLegendreOrders(const CompoundSystem &compound, const Config &configure);
/// Warns where the separation energies and the masses of a nucleus's particle pairs disagree by more than
/// kMassCheckTolerance; returns the number of disagreeing pairs.
int WarnMassSeparationMismatch(const CompoundSystem &compound, const Config &configure);
constexpr double kMassCheckTolerance = 1.0e-3;  ///< MeV
/// Warns of every R-matrix level whose radiative width exceeds 10 % of its particle width, from the observed widths
/// of `table` (CompoundSystem::ObservedWidthsAt).
void WarnLargeRadiativeWidths(const CompoundSystem &compound, const LevelInputs &inputs, const ObservedLevelTable &table, const Config &configure);

/// A bound state or a resonance of a potential that also acts inside the channel radius (f = S of the waves outside):
/// its channel energy (MeV, c.m.; negative: bound) and width (MeV; 0 for a bound state).
struct PotentialState {
  double Ec = 0.0;
  double width = 0.0;
};
/// The states of the potential of (pair, l, 2j) at `radius` with channel energies in [low, high] (MeV), found as
/// the zeros of f - S on a 50 keV grid (above threshold only those narrower than half their energy); empty for a
/// channel without an interior potential.
std::vector<PotentialState> InteriorPotentialStates(const PPair &pair, int l, int twoJ, double radius, double low, double high);
/// Warns of every level with a width in a channel whose potential also acts inside the channel radius that lies within
/// one width (the larger of the level's observed particle width and the state's) of a state of that potential: the
/// potential already produces it (docs/upgrade/interior_potential.md section 3.4, W_POT_LEVEL_DOUBLE).  Returns the
/// number of such levels.
int WarnLevelDoublesPotential(const CompoundSystem &compound, const LevelInputs &inputs, const ObservedLevelTable &table, const Config &configure);
/// Warns of every external-capture final state that gets no capture from the potential's interior although its pair
/// has a potential acting inside the channel radius (docs/upgrade/interior_potential.md section 6.4): a final state,
/// bound or not, in a channel whose own l has no interior boundary constant (W_POT_INTERIOR_FINAL).  Returns the
/// number of such final states.
int WarnInteriorCaptureFinal(const CompoundSystem &compound, const ObservedLevelTable &table, const Config &configure);

#endif
