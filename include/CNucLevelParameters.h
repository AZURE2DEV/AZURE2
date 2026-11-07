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
 * Channel functions of the blocks at any level energy E (MeV, E_x), from the engine's ChannelFunc and WhitFunc, as
 * CNuc::TransformIn evaluates them: E_c = E - GetExE() - GetSepE() (in the given order); particle channels S, dS/dE (ChannelFunc), P above
 * threshold (ChannelFunc::Penetrability) and mu a / (hbar^2 W^2) below it (WhitFunc); photon channels
 * (E_gamma/hbar c)^(2L+1) (1 with --use-rmc), or the static-moment factor and E2 scale; beta channels 1.
 * The CNuc must outlive the returned function; it is read, not modified.
 */
/*!
 * The order in which the engine forms the channel energy E_c from the level energy E (MeV, E_x), the excitation
 * energy e = GetExE() of the pair and its separation energy S.  The three orders agree when e = 0 and may differ by
 * one ulp otherwise; each engine caller keeps its own so that its numbers stay bit-identical
 * (docs/upgrade/cnuc_split.md section 5).  They become one when channel_functions replaces this adapter.
 */
enum class ChannelEnergyOrder {
  ExcitationFirst,  ///< (E - e) - S: CNuc::TransformIn
  SumFirst,         ///< E - (S + e): CNuc::CalcBoundaryConditions
  SeparationFirst   ///< (E - S) - e: CNuc::TransformOut
};

/*!
 * Which particle-channel values the returned function computes (the others are left 0): every one, only S_c (what
 * boundary_values, formal_from_brune and brune_from_formal read), or only P_c and dS_c/dE (what
 * amplitudes_from_observed reads).  The engine asks for what each caller reads, so that it evaluates the Coulomb
 * functions it evaluated before the move, in the same order, and no more.
 */
enum class ChannelQuantities { All, Shift, WidthFactors };

azure2::physics::ChannelFunctions ChannelFunctionsFromCNuc(CNuc &compound, const CNucSystemMap &map,
                                                           const Config &configure,
                                                           ChannelEnergyOrder order = ChannelEnergyOrder::ExcitationFirst,
                                                           ChannelQuantities what = ChannelQuantities::All);

/*!
 * The <levels> values of block b as the file gives them (form Observed): E_lambda, the widths (eV), ANCs (fm^-1/2),
 * amplitudes (gammaIsRWA) and their flags.  The CNuc must be filled (CNuc::Fill) and not initialized: after
 * CNuc::Initialize its levels hold the transformed amplitudes.
 */
azure2::physics::ProjectLevels ProjectLevelsFromCNuc(CNuc &filled, const CNucSystemMap &map, std::size_t b);

/*!
 * External-capture amplitudes of block b's photon channels (CNuc::CalcExternalWidth with isInitial = true), in the
 * layout of the block's gamma; not present unless external capture is on, nor (as in the engine) for an amplitude
 * |gamma| < 1e-8 when zero widths are ignored (IGNORE_ZERO_WIDTHS).  The value depends on the amplitudes the
 * CNuc's levels hold at the time of the call (the level's own particle amplitudes and the final state's): CNuc's
 * TransformIn calls it with the Brune (Park) amplitudes of the level after its first loop.
 */
std::vector<azure2::physics::ExternalWidth> ExternalWidthsFromCNuc(CNuc &compound, const CNucSystemMap &map,
                                                                   std::size_t b, const Config &configure);

#endif
