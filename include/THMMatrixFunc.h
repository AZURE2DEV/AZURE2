#ifndef THMMATRIXFUNC_H
#define THMMATRIXFUNC_H

#include "AMatrixFunc.h"

/// Modified R-matrix (Trojan Horse Method) cross-section calculator.

/*!
 * THM on the shared mesh (docs/upgrade/energy_grid.md section 8, row 10).  sigma_HOES at a node E_n is one more
 * observable of the reaction (SiteObservable.h, kHoes), and the THM fold one more use of the kernels' table.
 *
 * The kernel.  Today's THM fold is no arithmetic of its own: a THM point is folded by the rows of Kernels()
 * (experiment/src/Kernels.cpp) -- the Gaussian resolution in the c.m. frame of x + A, normalised over its window, or a
 * beam profile -- whose supports reach kThmConvolutionRange sigma for a THM point (SupportPoint::thm, set from the data
 * kind).  Nothing of the fold is in this class or in ThmNumerics (Coulomb phases, the bound-state tail, a cubic
 * Lagrange rule: pieces of the vertex at one energy).  The quasi-free window is not an energy fold: it is the average
 * of sigma_HOES over the accepted spectator directions at one energy (ThmSpectatorWindow, 2017 eq. 34), so its weights
 * w_k(E_n) are inside the site's value at the node.  The OES fold of S(E) (THM.md section 2) has no code today: OES
 * points are compared with the on-shell S factor of an ordinary data set; a fold of S would be the diagonal factor
 * sigma -> S per node before the kernel, as the screening factor is (not added here: no formula changes).  The
 * diagonal factors THM has per node (the weight table W(E), the distortion R(E)) are functions of E_n and the data
 * set or experiment, applied at the node inside the site's value, as before.
 *
 * The level matrix.  The HOES amplitude reads individual elements A_lambda,mu of the level matrix of the level-space
 * route (SetLevelSpaceOnly).  A direct reader at the same node is solved by the route its J-group takes (channel
 * space when it has fewer channels than levels).  The energy pass solves a node once for all its readers; a node whose
 * direct readers take the channel route for some J-group is filled a second time in level space for its THM readers
 * (the numbers of both kinds are then those of their own route, bit for bit); the task, its channel rows and the plan
 * are shared either way.
 */

/*!
 * THMMatrixFunc reuses the shared R-matrix interior of AMatrixFunc
 * (ClearMatrices / FillMatrices / InvertMatrices give the level matrix A,
 * unchanged) and replaces only the channel-surface assembly with the
 * half-off-energy-shell (HOES) cross section of the modified R-matrix
 * formalism used for the Trojan Horse Method:
 *
 *   sigma(E) = sum_J (2J+1) (k_f/mu_f) sum_{c_exit} 2 P_{c_exit}
 *                sum_{s_in} | sum_{la,lap} gamma_{la,c_exit} A_{la,lap}
 *                                sum_{c_in in s_in} gamma_{lap,c_in} M_l |^2
 *
 * coherent over levels (through A) and over entrance partial waves sharing a
 * channel spin (each carrying its own transfer form factor M_l, stored on the
 * point), incoherent over distinct entrance channel spins, exit channels and
 * J^pi groups. The entrance penetrability of conventional R-matrix is replaced
 * by M_l; the exit vertex keeps sqrt(2 P). The absolute scale is arbitrary
 * (folded into the segment normalization), as in mrmpy. See
 * docs/THM_IMPLEMENTATION.md and mrmpy MRMModel.unsmeared.
 */

/*!
 * What the HOES adjoint needs beyond the point and the filled compound,
 * computed once per derivative evaluation (BuildThmAdjointContext).
 */
struct ThmAdjointContext {
  /// dS_c/dE at the level energies (BuildShiftDerivTable); Brune only, else null.
  const vector_matrix_r *shiftDeriv = nullptr;
  /// vertex=constant: [j-1][ch-1] dS_c/dE at the lowest level of the J group.
  std::vector<std::vector<double>> lowestShiftDeriv;
  /// lineshape=on: [j-1][la-1] = {dGamma/dE, dGamma/dgamma_1, ...} of the
  /// level's total width (ThmLevelWidth); empty without a line shape.
  std::vector<std::vector<std::vector<double>>> widthDeriv;
  /// The present values of the THM coherent backgrounds (WorkingCopy::ThmCoherentValues), in parameter order.
  const std::vector<double> *coherentValues = nullptr;
};

/// The context for the compound at the levels of `state` (CalcShiftFunctions already run
/// under Brune); `lineshape` asks for the width derivatives.
ThmAdjointContext BuildThmAdjointContext(CNuc *compound, const LevelState &state,
                                         const std::vector<double> &coherentValues,
                                         const Config &configure, const vector_matrix_r *shiftDeriv, bool lineshape);

/*!
 * sigma_HOES of `point` (a site of kind kHoes, SiteObservable.h) into its record, from the level matrix `levels` holds
 * at the point's energy, filled in level space for every J-group of its nucleus (THMMatrixFunc, or an AMatrixFunc with
 * LevelSpaceThroughout(): the energy pass, CalculationPlan.h).  `coherentValues` are the present values of the THM
 * coherent backgrounds, in parameter order (WorkingCopy::ThmCoherentValues; empty without cbackground=).
 */
void HoesCrossSection(AMatrixFunc &levels, Site *point, const std::vector<double> &coherentValues,
                      azure2::Results &results);

class THMMatrixFunc : public AMatrixFunc {
 public:
  THMMatrixFunc(CNuc *compound, const LevelState &state, const Config &configure);
  /// Computes the HOES cross section for the point and stores it as the fit XS.
  void CalculateTHMCrossSection(Site *point, const std::vector<double> &coherentValues, azure2::Results &results);
  /*!
   * Reverse-mode derivative of the HOES cross section of `point` (one point
   * or folding sub-point, before any folding): accumulates fitBar d sigma/dE
   * and d sigma/dgamma into `accum` and, if given, fitBar d sigma/d(coherent
   * background value k) into (*coherentBar)[k] (ctx.coherentValues
   * order).  Assumes ClearMatrices/FillMatrices/InvertMatrices for the point.
   * The cross section itself, recomputed on the way, goes to *sigmaOut (it is
   * the value CalculateTHMCrossSection gives).  Returns false for a
   * configuration it does not cover (RMC; a line shape without widthDeriv).
   */
  bool THMPointAdjoint(Site *point, double fitBar, GradAccum &accum, const ThmAdjointContext &ctx,
                       std::vector<double> *coherentBar = nullptr, double *sigmaOut = nullptr);
};

#endif
