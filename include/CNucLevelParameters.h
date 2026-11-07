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
 * as input and physics/ cannot compute itself yet.  Nothing in the engine's calculation calls these; the test
 * tests/upgrade/level_parameters_from_cnuc_test.cpp does.  Block and channel indices are the System's (0-based), turned
 * into CNuc's through the CNucSystemMap of BuildSystemFromCNuc.
 */

/*!
 * Channel functions of the blocks at any level energy E (MeV, E_x), from the engine's ChannelFunc and WhitFunc, as
 * CNuc::TransformIn evaluates them: E_c = E - GetExE() - GetSepE(); particle channels S, dS/dE (ChannelFunc), P above
 * threshold (ChannelFunc::Penetrability) and mu a / (hbar^2 W^2) below it (WhitFunc); photon channels
 * (E_gamma/hbar c)^(2L+1) (1 with --use-rmc), or the static-moment factor and E2 scale; beta channels 1.
 * The CNuc must outlive the returned function; it is read, not modified.
 */
azure2::physics::ChannelFunctions ChannelFunctionsFromCNuc(CNuc &compound, const CNucSystemMap &map,
                                                           const Config &configure);

/*!
 * The <levels> values of block b as the file gives them (form Observed): E_lambda, the widths (eV), ANCs (fm^-1/2),
 * amplitudes (gammaIsRWA) and their flags.  The CNuc must be filled (CNuc::Fill) and not initialized: after
 * CNuc::Initialize its levels hold the transformed amplitudes.
 */
azure2::physics::ProjectLevels ProjectLevelsFromCNuc(CNuc &filled, const CNucSystemMap &map, std::size_t b);

/*!
 * External-capture amplitudes of block b's photon channels (CNuc::CalcExternalWidth with isInitial = true), in the
 * layout of the block's gamma; not present unless external capture is on.  The value depends on the amplitudes the
 * CNuc's levels hold at the time of the call (the level's own particle amplitudes and the final state's): CNuc's
 * TransformIn calls it with the Brune (Park) amplitudes of the level after its first loop.
 */
std::vector<azure2::physics::ExternalWidth> ExternalWidthsFromCNuc(CNuc &compound, const CNucSystemMap &map,
                                                                   std::size_t b, const Config &configure);

#endif
