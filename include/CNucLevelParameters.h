#ifndef CNUCLEVELPARAMETERS_H
#define CNUCLEVELPARAMETERS_H

#include <cstddef>
#include <vector>

#include "CNucSystem.h"
#include "azure2/physics/LevelParameters.h"

class CNuc;
class Config;

/*!
 * The engine's side of physics/ level_parameters (docs/upgrade/level_parameters.md): what the physics functions take
 * as input and physics/ cannot compute itself yet.  CNuc::CalcBoundaryConditions, CNuc::TransformIn and
 * CNuc::TransformOut call level_parameters through these (docs/upgrade/cnuc_split.md); so does the test
 * tests/upgrade/level_parameters_from_cnuc_test.cpp.  Block and channel indices are the System's (0-based), turned
 * into CNuc's through the CNucSystemMap of BuildSystemFromCNuc.
 */

/*!
 * Channel functions of the blocks at any level energy E (MeV, E_x), from the engine's ChannelFunc and WhitFunc, with
 * the channel energy E_c = (E - S) - e (MeV, c.m.; S = GetSepE(), e = GetExE()), the order in which the collision
 * matrix forms it at a data point (Unit.cpp) and CNuc::CalcShiftFunctions at a level, so that a level and a data
 * point at the same E_x see the same E_c.  Every call returns all values of the channel: particle channels S, dS/dE
 * and P, evaluated on one ChannelFunc in that order (P = ChannelFunc::Penetrability above threshold and
 * mu a / (hbar^2 W^2) from WhitFunc below it, 0 at threshold); photon channels (E_gamma/hbar c)^(2L+1) (1 with
 * --use-rmc), or the static-moment factor and E2 scale; beta channels 1 (docs/upgrade/cnuc_split.md section 7).
 * The CNuc must outlive the returned function; it is read, not modified.
 */
azure2::physics::ChannelFunctions ChannelFunctionsFromCNuc(CNuc &compound, const CNucSystemMap &map,
                                                           const Config &configure);
/// The same, with the two settings it reads from the Config given directly: GSL Coulomb functions, --use-rmc.
azure2::physics::ChannelFunctions ChannelFunctionsFromCNuc(CNuc &compound, const CNucSystemMap &map, bool gslCoulomb,
                                                           bool rmc);

/*!
 * The <levels> values of block b as the file gives them (form Observed): E_lambda, the widths (eV), ANCs (fm^-1/2),
 * amplitudes (gammaIsRWA) and their flags.  The CNuc must be filled (CNuc::Fill) and not initialized: after
 * CNuc::Initialize its levels hold the transformed amplitudes.
 */
azure2::physics::ProjectLevels ProjectLevelsFromCNuc(CNuc &filled, const CNucSystemMap &map, std::size_t b);

/*!
 * External-capture amplitudes of block b's photon channels (ExternalReducedWidth with initial = true), in the
 * layout of the block's gamma; not present unless external capture is on, nor (as in the engine) for an amplitude
 * |gamma| < 1e-8 when zero widths are ignored (IGNORE_ZERO_WIDTHS).  The value depends on the amplitudes the
 * CNuc's levels hold at the time of the call (the level's own particle amplitudes and the final state's): CNuc's
 * TransformIn calls it with the Brune (Park) amplitudes of the level after its first loop.
 */
std::vector<azure2::physics::ExternalWidth> ExternalWidthsFromCNuc(CNuc &compound, const CNucSystemMap &map,
                                                                   std::size_t b, const Config &configure);

#endif
