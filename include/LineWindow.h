#ifndef AZURE2_LINEWINDOW_H
#define AZURE2_LINEWINDOW_H

#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Constants.h"

class CompoundSystem;
class Config;
class GenMatrixFunc;
class LevelState;
class Site;

namespace azure2 {

struct Results;

/*!
 * Capture into a BROAD unbound final state: the line-shape integration (the field `line_window` of a capture data
 * set or extrapolation; dev's "linewindow W", 2026-11-15).
 *
 * The final state is the energy-normalized scattering state of its one open particle channel at E_f: interior
 * c(E_f) X_lambda, exterior sqrt(2 mu/(pi hbar^2 k_f)) [cos d F + sin d G], d = d_R - phi from the one-level
 * parameters of the final level (observed g~ and E_R, formal gamma^2 = g~^2/(1 - g~^2 S'(E_R)), B = S(E_R)).  Per
 * E_f the cross section is that of the narrow resonance (ECIntegral::IntegrateUnbound) with
 *   - the photon energy E - E_x(E_f) in the internal capture elements and the external integral,
 *   - the internal elements times A(E_f) = c/sqrt(N_f) = sqrt(2 mu/(pi hbar^2 k_f)) (cos d F_a + sin d G_a)/(g~ s),
 *     s = sqrt(2 mu a)/hbar, A^2 = rho(E_f)/N_f (the line shape, unit area for a narrow state),
 *   - the external element sqrt(2 mu/(pi hbar^2 k_f)) [cos d I_F + sin d I_G] (= g~ s A x the integral with the
 *     standing wave normalized at a; written without that division, whose denominator vanishes off resonance),
 * and sigma = int_{E_R-W}^{E_R+W} dE_f sigma(E_f), clipped to E_f > 1 keV and E_gamma > 1 keV; capture into different
 * E_f adds incoherently.  Quadrature: Gauss-Legendre in theta, E_f = E_R + (Gamma/2) tan theta, fixed per point
 * energy and line centre; the weights carry the current parameters.  The line is centred on the final level's
 * observed energy from the current parameters in the Brune parametrizations (a fitted final-level energy moves it),
 * on the energy of the set-up in the formal one.
 */

/// One node of the final-energy quadrature as the T-matrix builders read it (AMatrixFunc / RMatrixFunc::
/// CalculateTMatrix, through GenMatrixFunc::ActiveLineNode) while LineWindowCrossSection runs.  photonEnergy =
/// E - E_x(node); amp = A(E_f), the final-state amplitude that multiplies the internal capture elements;
/// ec[k-1][m-1] = the external-capture T element of that node (it replaces the stored amplitude's element).
struct LineNode {
  double photonEnergy = 0.0;
  double amp = 0.0;
  std::vector<std::vector<complex>> ec;
};

/// One quadrature node, fixed per point energy and line centre: E_f, its weight (dE_f included),
/// sqrt(2 mu / (pi hbar^2 k_f)), F, G, P, S of the final channel at E_f, and the external-capture amplitude split by
/// exterior function (F_f, G_f) for each (spin group, external pathway) of the exit.
struct LineNodeData {
  double ef = 0.0, weight = 0.0, sqrtDens = 0.0, Fa = 0.0, Ga = 0.0, P = 0.0, S = 0.0;
  std::vector<std::vector<std::pair<complex, complex>>> parts;
};

/// The nodes of one site at one line centre, with the terms of 1/N_f at that centre (LineWindow.cpp).  Immutable
/// once kept.
struct LineNodeSet {
  double centre = 0.0;
  std::vector<LineNodeData> nodes;
  std::vector<double> nfTerms;
};

/// The node sets of a run, keyed by what they depend on besides the line centre (entrance and exit keys, the site's
/// c.m. energy, the window): the newest two centres of each are kept, so a fit that moves the centre does not grow
/// the table.  Owned by the run (RunScratch::lineNodes), shared by its threads behind the lock.
struct LineNodeStore {
  std::mutex mutex;
  std::map<std::vector<double>, std::vector<std::shared_ptr<const LineNodeSet>>> sets;
};

/// Why this site's line window cannot be computed ("" if it can): the exit must be a capture pair whose final level
/// lies above exactly one particle threshold, outside the RMC formalism, for a cross section (not angular
/// coefficients, a phase or a polarization observable).
std::string LineWindowProblem(const Site &site, const CompoundSystem *compound, const Config &configure);
/// The same for a data set's reaction before its points exist (`otherObservable`: not a cross section).
std::string LineWindowProblem(int entranceKey, int exitKey, bool otherObservable, const CompoundSystem *compound,
                              const Config &configure);

/// The capture observable of `site` integrated over the final state's line shape, E_f in E_R +- W, with the level
/// matrix already solved in `func`: the site's record is the weighted sum of the nodes' records.  Leaves the T-matrix
/// of the last node in `func`.
void LineWindowCrossSection(GenMatrixFunc &func, Site *site, Results &results);

/// The fixed-node pass of the adjoint over the line shape: node(weight) for each node of the site's quadrature (the
/// nodes of LineWindowCrossSection at the current parameters), that node active in `func` and its T-matrix built, so
/// that d(record)/d(parameters) at fixed nodes is the weighted sum of the nodes' adjoints.  Leaves the T-matrix of
/// the last node in `func`.
void ForEachLineNode(GenMatrixFunc &func, Site *site, const std::function<void(double)> &node);

/// What the line shape itself takes from the parameters, besides the level matrix: the final level's energy (the
/// line centre, in the Brune parametrizations) and its observed amplitude g~, which place the nodes and set the
/// phase and the amplitude of the final scattering state.  The adjoint takes d(record)/d(centre) and
/// d(record)/d(ln g~) by central differences of the line integration itself (LineWindowCrossSectionAt, steps
/// `centreStep` in MeV and `amplitudeStep`), so that they follow the value with its re-placed nodes, and chains
/// them to the widths with dLogAmplitude[ch - 1] = d ln g~ / d gamma_ch (the channels of the level's J-group).
struct LineShapeChain {
  int group = 0, level = 0;  ///< the final level
  bool energy = false;       ///< the centre follows the level's energy parameter
  std::vector<double> dLogAmplitude;
  double centreStep = 0.0, amplitudeStep = 0.0;
};
LineShapeChain LineShapeChainOf(const Site &site, const CompoundSystem *compound, const LevelState &state,
                                const Config &configure);

/// LineWindowCrossSection with the line centre moved by dCentre (MeV; its own nodes, not kept by the run) or, on the
/// nodes of the evaluation, the observed amplitude times (1 + dLogAmplitude).
void LineWindowCrossSectionAt(GenMatrixFunc &func, Site *site, Results &results, double dCentre, double dLogAmplitude);

}  // namespace azure2

#endif
