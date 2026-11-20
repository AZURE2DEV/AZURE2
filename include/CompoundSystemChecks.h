#ifndef COMPOUNDSYSTEM_CHECKS_H
#define COMPOUNDSYSTEM_CHECKS_H

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
/// The `wave` fields of the pairs' potentials against the channels (potential_background.md section 6): an interior
/// background on a channel whose radius lies beyond the matching radius (it would not act) and a `subtract` naming a
/// box state the channel does not have are refused (-1, after the message); a `wave` whose (l, J) is no channel of
/// its pair is warned about (it does nothing).
int CheckPotentialWaves(const CompoundSystem &compound, const Config &configure);
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

#endif
