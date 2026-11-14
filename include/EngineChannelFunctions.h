#ifndef ENGINECHANNELFUNCTIONS_H
#define ENGINECHANNELFUNCTIONS_H

#include <cstddef>
#include <vector>

#include "SystemOfCompound.h"
#include "azure2/physics/LevelParameters.h"

class CompoundSystem;
class Config;
class LevelInputs;
class LevelState;

/*!
 * The engine's side of physics/ level_parameters (docs/upgrade/level_parameters.md): what the physics functions take
 * as input and physics/ cannot compute itself yet.  CompoundSystem::CalcBoundaryConditions, CompoundSystem::TransformIn and
 * CompoundSystem::TransformOut call level_parameters through these (docs/upgrade/cnuc_split.md); so does the test
 * tests/upgrade/level_parameters_of_compound_test.cpp.  Block and channel indices are the System's (0-based), turned
 * into CompoundSystem's through the SystemOfCompoundMap of BuildSystemOfCompound.
 */

class PPair;

/*!
 * The functions of one particle channel of `pair` with orbital l and 2j (PPair::SpinOrbitTwoJ) at the channel energy
 * Ec (MeV, c.m.): S, dS/dE and P, evaluated on one ChannelFunc in that order (P = ChannelFunc::Penetrability above
 * threshold, mu a / (hbar^2 W^2) from WhitFunc below it, 0 at threshold).  `gsl`: GSL's Coulomb functions.  The one
 * implementation of a particle channel's functions: EngineChannelFunctions and the binding's level doors
 * (api/src/pybind_level_doors.cpp) call it.
 */
azure2::physics::ChannelValues ParticleChannelValues(const PPair &pair, int l, int twoJ, double Ec, bool gsl);

/*!
 * Channel functions of the blocks at any level energy E (MeV, E_x), from the engine's ChannelFunc and WhitFunc, with
 * the channel energy E_c = (E - S) - e (MeV, c.m.; S = GetSepE(), e = GetExE()), the order in which the collision
 * matrix forms it at a data point (Unit.cpp) and CompoundSystem::CalcShiftFunctions at a level, so that a level and a data
 * point at the same E_x see the same E_c.  Every call returns all values of the channel: particle channels S, dS/dE
 * and P, evaluated on one ChannelFunc in that order (P = ChannelFunc::Penetrability above threshold and
 * mu a / (hbar^2 W^2) from WhitFunc below it, 0 at threshold); photon channels (E_gamma/hbar c)^(2L+1) (1 with
 * --use-rmc), or the static-moment factor and E2 scale; beta channels 1 (docs/upgrade/cnuc_split.md section 7).
 * The CompoundSystem must outlive the returned function; it is read, not modified.
 */
azure2::physics::ChannelFunctions EngineChannelFunctions(const CompoundSystem &compound, const SystemOfCompoundMap &map,
                                                           const Config &configure);
/// The same, with the two settings it reads from the Config given directly: GSL Coulomb functions, --use-rmc.
azure2::physics::ChannelFunctions EngineChannelFunctions(const CompoundSystem &compound, const SystemOfCompoundMap &map,
                                                           bool gslCoulomb, bool rmc);

/*!
 * The input values of block b (form Observed): E_lambda, the widths (eV), ANCs (fm^-1/2), amplitudes (gammaIsRWA)
 * and their flags, from `inputs` before the transformation (CompoundSystem::Fill's, or a physical vector's): after
 * CompoundSystem::TransformIn they are the transformed amplitudes.
 */
azure2::physics::ProjectLevels ProjectLevelsOfInputs(const CompoundSystem &filled, const LevelInputs &inputs,
                                                     const SystemOfCompoundMap &map, std::size_t b);

/*!
 * External-capture amplitudes of block b's photon channels (ExternalReducedWidth with initial = true), in the
 * layout of the block's gamma; not present unless external capture is on, nor (as in the engine) for an amplitude
 * |gamma| < 1e-8 when zero widths are ignored (IGNORE_ZERO_WIDTHS).  The value depends on the amplitudes the
 * CompoundSystem's levels hold at the time of the call (the level's own particle amplitudes and the final state's): CompoundSystem's
 * TransformIn calls it with the Brune (Park) amplitudes of the level after its first loop.  The amplitudes are those
 * of `inputs`, the final states' N_f^1/2 those of `state` (ExternalReducedWidthOfInputs).
 */
std::vector<azure2::physics::ExternalWidth> ExternalWidthsOfInputs(const CompoundSystem &compound, const LevelState &state,
                                                                   const LevelInputs &inputs, const SystemOfCompoundMap &map,
                                                                   std::size_t b,
                                                                   const Config &configure);

#endif
