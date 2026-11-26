#ifndef LEVELPOLES_H
#define LEVELPOLES_H

#include <vector>

/*!
 * The poles of the levels of one J^pi group, as an estimate for the energy mesh.
 *
 * The mesh resolves a level with nodes at its width / points per width, the width being the sum of the level's
 * partial widths (AdaptiveIntegrationGrid::LevelResonances).  That is the width of the resonance only for an
 * isolated level.  Two levels of one J^pi that decay through a common channel and lie closer than their widths
 * interfere: the cross section then has one broad structure and one much narrower than either level's partial
 * widths (two levels of widths G1 << G2 at a distance dE << G2 in one channel: about 4 dE^2 G1 G2 / (G1 + G2)^3).
 * A fold across it was under-resolved (22Ne(a,n), 1- levels at E_x 11.8847 and 11.8913 MeV, 3.1 and 97 keV in the
 * n0 channel: a peak 33 eV wide where the mesh took 3.1 keV; a folded point off by 25 %, 2026-11-25).
 *
 * The estimate: the eigenvalues of E_lambda delta_lambda,mu - (i/2) sum_c a_lambda,c a_mu,c, with a_lambda,c the
 * signed square root of the level's partial width in the open particle channel c (the single-pole form of the
 * level matrix at the levels' own energies, without the energy dependence of the penetrabilities and shifts).
 * Their real parts are the pole energies and -2 times their imaginary parts the pole widths.
 */
namespace azure2 {

struct PoleLevel {
  double energy = 0.0;            ///< level energy (MeV, any common frame)
  double width = 0.0;             ///< the sum of its partial widths in the open particle channels (MeV)
  std::vector<double> amplitude;  ///< per channel of the group: signed sqrt(partial width / MeV), 0 where none
};

struct NarrowPole {
  double energy = 0.0;  ///< pole energy, in the frame of the levels' energies (MeV)
  double width = 0.0;   ///< pole width (MeV)
};

/// (Half, not less: a pole 0.5 to 1 times as wide as its levels is resolved by the levels' own anchors within a
/// factor 2 of the nodes, and anchors for such poles move the folded numbers of projects that have no narrow
/// structure at the 1e-5 level, e.g. tests/upgrade T2 by 4e-5; 22Ne(a,n) is the same with 0.5 and 0.8.)
/// The poles that are narrower than half the width of every level with at least a tenth of the pole's
/// weight, in ascending energy: the structures the levels' own widths do not resolve.  Empty for fewer than two
/// levels, when no two levels with a common channel lie closer than the sum of their widths (the levels are then
/// taken as isolated and the eigenvalues are not computed), and when the eigenvalue routine fails.
std::vector<NarrowPole> NarrowPolesOf(const std::vector<PoleLevel> &levels);

}  // namespace azure2

#endif
