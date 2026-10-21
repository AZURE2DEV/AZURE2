#ifndef AMATRIXFUNC_H
#define AMATRIXFUNC_H

#include "GenMatrixFunc.h"
#include "LevelMatrixSolver.h"
#include "ChannelSpace.h"
#include <memory>

struct GradAccum;
class JGroup;
class Decay;

/*!
 * Cotangents of the T-matrix elements of one decay (AMatrixFunc::PointTBar):
 * t[k-1][m-1] for the internal pathways T(k,m), ec[k-1][m-1] for the
 * external-capture ones.  Same convention as PointAdjoint: zbar = dL/dRe(z) +
 * i dL/dIm(z).
 */
struct TMatrixCotangent {
  std::vector<vector_c> t;
  std::vector<vector_c> ec;
  /// Electron screening: (screening parameter k, d/dU_e in eV) of the points these cotangents came from
  /// (PointTBar), added to the accumulator by TBarAdjoint.
  std::vector<std::pair<int, double>> screen;
};

/*!
 * One term beta * g_ch^T A g_chp of the level-matrix bilinear forms of a
 * J-group: the internal T-matrix depends on the level matrix only through these
 * (AMatrixFunc::CollectBilinearTerms), so its adjoint is a sum over them
 * (AMatrixFunc::BilinearAdjoint).  beta is the cotangent of the bilinear form.
 */
struct BilinearTerm {
  int ch;
  int chp;
  complex beta;
};

/// A function class to calculate the T-Matrix using the A-Matrix

/*!
 * The AMatrixFunc function class calculates the T-Matrix for a given energy point using the compound
 * nucleus object.  The AMatrixFunc class is a child class of GenMatrixFunc, where the cross section is
 * calculated from the T-Matrix.
 *
 * Each J-group takes one of two routes at each energy:
 *
 *   - channel space (the default, ChannelSpace.h): the level matrix is never
 *     formed; the bilinear forms \f$g_c^T A g_{c'}\f$ come from the c x c matrix
 *     \f$(1 - R\Lambda)^{-1} R\f$ of the formal representation of the level
 *     matrix.  Taken when the J-group has fewer coupled channels than levels;
 *   - level space: the level matrix is filled and LU-factored.  Taken otherwise,
 *     near a formal pole or for an ill-conditioned \f$1 - R\Lambda\f$, when the
 *     formal representation cannot be built, with AZURE_LEVEL_SPACE=1, and for
 *     an object switched to it with SetLevelSpaceOnly().
 *
 * Both routes serve the same interface: GetUBilinear, the channel-capture and
 * general bilinear forms, the A g_c solves of the adjoint and, through the
 * Woodbury identity, individual elements of A (GetAMatrixElement).
 */

class AMatrixFunc : public GenMatrixFunc {
 public:
  AMatrixFunc(CNuc *, const Config &configure);
  /*!
   * Returns a pointer to the compound nucleus object.
   */
  CNuc *compound() const { return compound_; };
  /*!
   * Returns a reference to the Config structure.
   */
  const Config &configure() const { return *configure_; };
  /*!
   * Repoints the object at a compound nucleus and configuration so that one
   * instance can serve many energy points without reallocating its scratch
   * storage.  ClearMatrices() must still be called before each point.
   */
  void Reset(CNuc *compound, const Config &configure) {
    compound_ = compound;
    configure_ = &configure;
  };
  /*!
   * Solves every J-group of this object in level space (fill and LU), as
   * before the channel-space route existed.  For callers that want the level
   * matrix itself or many of its elements.
   */
  void SetLevelSpaceOnly(bool levelOnly) { level_only_ = levelOnly; };
  /// True if J-group jGroupNum (1-based) is solved in channel space at the current energy.
  bool IsChannelRoute(int jGroupNum) const {
    return jGroupNum >= 1 && jGroupNum <= (int)route_.size() && route_[jGroupNum - 1] == kRouteChannel;
  };
  /*!
   * Turns the formal-parameter cotangents a gradient pass accumulated in
   * `accum` (GradAccum::fe/fg, see FormalBilinearAdjoint) into cotangents of the
   * level energies and reduced widths, through the adjoint of the generalized
   * eigen-decomposition that defines the formal parameters, and clears them.
   * The parameters held by the compound must be those of the pass.
   */
  void ResolveFormalCotangents(GradAccum &accum, const vector_matrix_r *shiftDeriv);

  /// Drop the stored level matrices.
  void ClearMatrices();
  /// Build \f$A^{-1}\f$ per J-group from the level parameters at this point.
  void FillMatrices(EPoint *);
  /*!
   * Factors the level matrix of every J-group.  The full inverse is no longer
   * formed here: the T-matrix needs only bilinear forms of it, and Inverse() is
   * materialized lazily for the callers that need individual elements.
   */
  void InvertMatrices();
  /// T-matrix per pathway from the inverted level matrix.
  void CalculateTMatrix(EPoint *);
  /*!
   * Instantiated in the parent class.
   */
  void CalculateCrossSection();

  /// Level-matrix element at (J-group, level, level'), all 1-based.
  complex GetAMatrixElement(int, int, int) const;
  const matrix_c &GetAMatrix(int) const;
  matrix_c *GetJSpecAInvMatrix(int);

  /*!
   * Returns \f$ \sum_{\lambda\mu} \gamma_{\lambda c} \gamma_{\mu c'}
   * A_{\lambda\mu} \f$ for a J-group, the only way the A-matrix enters the
   * internal T-matrix.  Reduced widths below 1e-12 are dropped, reproducing the
   * cutoff the explicit double loop over levels applied.
   */
  complex GetUBilinear(int jGroupNum, int chNum, int chpNum);
  /*!
   * The same bilinear form for the channel-capture pathway, where levels whose
   * reduced width in `maskChannel` is negligible are excluded from both indices.
   * Pass maskChannel = 0 to include every level.  No cutoff is applied to the
   * entrance and exit widths themselves.
   */
  complex GetChannelCaptureBilinear(int jGroupNum, int chNum, int chpNum,
                                    int maskChannel);

  /*!
   * Reverse-mode (adjoint) gradient of one energy point's cross section
   * (angle-integrated or differential).  Assumes FillMatrices/InvertMatrices/
   * CalculateTMatrix have already been run for `point`.  `fitBar` is the
   * cotangent dlnL/dmodel; the resulting dlnL/dE and dlnL/dgamma are
   * accumulated into `accum`.
   *
   * Under the Brune formalism the level-energy gradient acquires explicit terms
   * through S(E_lambda); pass `shiftDeriv` (from BuildShiftDerivTable) to have
   * them included analytically.  If `shiftDeriv` is null the energy gradient
   * omits those terms (only valid for non-Brune models).
   *
   * Handles angle-integrated, differential (incl. UPOS) and phase-shift cross
   * sections.  `xsComponent` selects which cross-section component the cotangent
   * refers to (0 = full, 1 = E1, 2 = E2), used for E1/E2 component segments.
   *
   * Returns false (accumulating nothing) when the point's configuration is
   * outside the supported path (RMC / active-level compaction / angular-dist
   * coefficients), so the caller can fall back to finite differences.  See
   * PLAN.md Phases 2-4, 6.
   *
   * `ayCrossSection` (analyzing-power points only): back-propagate the
   * differential cross section the point carries alongside A_y instead of A_y
   * itself.  A target-integrated A_y weights each sub-point's A_y by that cross
   * section, so its adjoint needs both (GradTargetEffectAdjoint).
   */
  bool PointAdjoint(EPoint *point, double fitBar, GradAccum &accum,
                    const vector_matrix_r *shiftDeriv = nullptr,
                    int xsComponent = 0, bool ayCrossSection = false);

  /*!
   * The first half of PointAdjoint: the cotangents of the T-matrix elements
   * held by this object, for `fitBar` times the derivative of `point`'s
   * observable.  The T-matrix may be that of another point at the same energy
   * (a host of point mapping, EData::MapData): only the angle-dependent parts
   * (Legendre polynomials, Coulomb amplitude, A_y coefficients) are read from
   * `point`.  Returns false, with `bar` unusable, where PointAdjoint would.
   */
  bool PointTBar(EPoint *point, double fitBar, TMatrixCotangent &bar,
                 int xsComponent = 0, bool ayCrossSection = false);
  /*!
   * The second half of PointAdjoint: back-propagates T-matrix cotangents to the
   * level energies and reduced widths, at the energy of `energyPoint` (whose
   * T-matrix this object holds).  Linear in `bar`, so the cotangents of all the
   * points sharing one T-matrix may be summed first.
   */
  void TBarAdjoint(EPoint *energyPoint, const TMatrixCotangent &bar, GradAccum &accum,
                   const vector_matrix_r *shiftDeriv = nullptr);
  /*!
   * The internal-pathway cotangents of `bar` as bilinear-form terms, per
   * J-group (terms[j-1]), merged by channel pair.  False if the decay has
   * external capture, whose channel-capture terms are not of this form
   * (TBarAdjoint then takes the dense path).
   */
  bool CollectBilinearTerms(EPoint *energyPoint, const TMatrixCotangent &bar,
                            std::vector<std::vector<BilinearTerm>> &terms);
  /*!
   * Adjoint of the bilinear terms of J-group j: the explicit reduced-width
   * derivatives and, through d(A) = -A d(M) A, the level-matrix ones, at the
   * energy of `energyPoint`.  The cotangent of A is a sum of rank-one terms, so
   * that of the level matrix is formed from the solves A g_c in O(n^2) per
   * term instead of the two O(n^3) products of the dense path.
   */
  void BilinearAdjoint(EPoint *energyPoint, int jNum, const std::vector<BilinearTerm> &terms,
                       GradAccum &accum, const vector_matrix_r *shiftDeriv = nullptr);
  /// Store an \f$A^{-1}\f$ element.
  void AddAInvMatrixElement(int, int, int, complex);
  /*!
   * The single-term form of BilinearAdjoint for the diagonal form B = g_ch^T A g_ch of one channel:
   * accumulates Re(conj(mBar) dB) into accum through the explicit widths of channel chNum and the
   * level matrix (energies and every width of the group, standard, Brune or Park).  The matrices must
   * have been filled and inverted for `point`.  Used for the scattering length
   * (ScatteringLength::Adjoint).  Returns false, accumulating nothing, for the RMC formalism or a group
   * with levels left out of the R-matrix.
   */
  bool BilinearAdjoint(EPoint *point, int jGroupNum, int chNum, complex mBar, GradAccum &accum,
                       const vector_matrix_r *shiftDeriv = nullptr);

 private:
  /// One memoized channel-capture bilinear form, keyed by its pathway.
  struct ChanCapEntry {
    int jGroupNum;
    int chNum;
    int chpNum;
    int maskChannel;
    complex value;
  };

  const vector_r &GetGammaVector(int jGroupNum, int chNum);
  const vector_c &GetChannelSolve(int jGroupNum, int chNum);
  /// A times the reduced-width vector of a channel WITHOUT the 1e-12 cutoff
  /// (the adjoint differentiates the widths as they are), cached per energy.
  const vector_c &GetRawChannelSolve(int jGroupNum, int chNum);
  /// Beam energy in the compound system, as FillMatrices takes it.
  double InEnergy(EPoint *point);
  /// d(level matrix) part of the adjoint of one J-group, given the level-matrix
  /// cotangent Mbar (n x n): structural energy and reduced-width gradients of
  /// the standard, Brune and Park level matrices.  The channel kernel is read
  /// from `energyPoint` (Lo-matrix plus, for Brune and Park particle channels,
  /// the boundary condition) unless `lambda` gives it per channel.
  void LevelMatrixAdjoint(EPoint *energyPoint, int jNum, const matrix_c &Mbar, double inEnergy,
                          GradAccum &accum, const vector_matrix_r *shiftDeriv,
                          const complex *lambda = nullptr);

  enum Route : char { kRouteNone = 0, kRouteLevel = 1, kRouteChannel = 2 };
  /// Per-energy channel-space state of one J-group.
  struct ChannelState {
    std::vector<complex> dinv;   ///< 1/(e_lambda - E), n
    std::vector<complex> lam;    ///< Lambda per channel, nch
    std::vector<complex> lamAct; ///< Lambda over the active channels, m
    std::vector<complex> R;      ///< R over the active channels, m x m
    std::vector<complex> lu;     ///< LU factors of 1 - R Lambda
    std::vector<int> piv;
    std::vector<complex> F;      ///< (1 - R Lambda)^-1 R, column b at F[b*m]
    std::vector<char> fcol;      ///< columns of F computed
    std::vector<vector_c> usolve;  ///< formal channel solves A_f Gamma_c (per channel)
    std::vector<char> uvalid;
    matrix_c inverse;            ///< A, materialized on request
    bool inverseValid = false;
    vector_c w1, w2, w3;         ///< scratch
    std::vector<double> rreal;   ///< scratch of the real R
    std::vector<double> P;       ///< pencil route: LU factors of H - E T (upper Hessenberg)
    std::vector<char> hswap;     ///< pencil route: row interchanges of that LU
    std::vector<double> X;       ///< pencil route: (H - E T)^-1 Qt G over the active channels, n x m
  };
  /// Level matrix of J-group j at this energy, into a_inv_matrices_ (level route).
  void FillLevelMatrix(int j, EPoint *point, double inEnergy);
  /// Channel-space solve of J-group j at this energy; false sends it to the level route.
  bool SolveChannelSpace(int j, EPoint *point, double inEnergy);
  /// R of the pencil route (no formal representation); false sends the J-group to the level route.
  bool SolvePencil(int j, double inEnergy);
  /// R of the formal route; false sends the J-group to the level route.
  bool SolveFormal(int j, double inEnergy);
  /// Column b (position among the active channels) of (1 - R Lambda)^-1 R.
  const complex *FColumn(int j, int b);
  /// A u in channel space (level-space vectors over the active levels).
  void ApplyA(int j, const complex *u, complex *out);
  /// u^T A v in channel space.
  complex ChannelBilinear(int j, const complex *u, const complex *v);
  /// Formal channel solve A_f Gamma_c (zero for a channel without widths), cached per energy.
  const vector_c &FormalSolve(int j, int ch0);
  /*!
   * Adjoint of the bilinear terms of J-group j in the formal representation:
   * O(n c) per term, no n x n matrix.  The cotangents go to the level
   * parameters directly for the standard level matrix (V = 1), and to the
   * formal ones (GradAccum::fe/fg) otherwise.
   */
  void FormalBilinearAdjoint(int j, const std::vector<BilinearTerm> &terms, GradAccum &accum);
  /// Builds tbar_q_ for the T-matrix held now (see the definition).
  bool TBarProducts(Decay *theDecay, bool particleExit);
  /// Cotangent coefficients of the differential sum, [pathway][order], over
  /// Decay::InterferencePathways and its orders.
  std::vector<complex> tbar_q_;
  std::vector<double> tbar_legendre_;
  const Decay *tbar_decay_ = nullptr;
  bool tbar_particle_ = false;
  /// AZURE_CHANNEL_VERIFY: compare the channel-space solve of J-group j with the level matrix.
  void VerifyChannelSpace(int j, EPoint *point, double inEnergy);

  const Config *configure_;
  CNuc *compound_;
  vector_matrix_c a_inv_matrices_;
  /// Factored level matrix of each J-group.
  std::vector<LevelMatrixSolver> solvers_;
  /// Whether the corresponding solver holds a usable factorization.
  std::vector<char> solver_valid_;
  // private:
  std::vector<std::vector<int>> level_active_index_;  // [jGroup-1][origLevel] -> activeIndex (1..k) or 0 if inactive

  // Per-point caches over active levels, indexed [jGroup-1][channel-1].
  /// Reduced widths of a channel, with the 1e-12 cutoff applied.
  std::vector<std::vector<vector_r>> gamma_vectors_;
  std::vector<std::vector<char>> gamma_valid_;
  /// A times the reduced-width vector of a channel.
  std::vector<std::vector<vector_c>> channel_solves_;
  std::vector<std::vector<char>> channel_solve_valid_;
  std::vector<std::vector<vector_c>> raw_solves_;
  std::vector<std::vector<char>> raw_solve_valid_;
  /// Scratch level-matrix cotangent of BilinearAdjoint.
  matrix_c mbar_scratch_;
  /// Scratch of LevelMatrixAdjoint.
  std::vector<double> lma_e_, lma_g_, lma_s_, lma_inv_, lma_resym_, lma_imsym_, lma_de_, lma_sg_, lma_acc_;
  std::vector<char> lma_brunep_;
  std::vector<complex> lma_l_;
  std::vector<ChanCapEntry> chan_cap_cache_;
  vector_c chan_cap_entrance_;
  vector_c chan_cap_exit_;

  // Pre-allocated buffers to prevent memory fragmentation
  mutable std::vector<std::vector<double>> levelGammas_;
  mutable std::vector<double> levelEnergies_;
  mutable std::vector<std::vector<double>> shiftFunctions_;
  mutable int cached_max_levels_;
  mutable int cached_max_channels_;
  /// Compound nucleus (1-based) of the point being evaluated: groups of the
  /// others are skipped by FillMatrices and InvertMatrices.
  int activeNucleus_ = 0;

  /// Route of each J-group at the current energy.
  std::vector<char> route_;
  /// Formal representation of each J-group (per parameter set).
  std::vector<std::shared_ptr<ChannelSpace::FormalLevels>> formal_;
  /// J-group each entry of formal_ belongs to (the object may be Reset to another compound).
  std::vector<const JGroup *> formal_jgroup_;
  /// Channel-space state of each J-group (per energy).
  std::vector<ChannelState> chan_;
  /// Active levels of each J-group (original numbers).
  std::vector<std::vector<int>> active_levels_;
  bool level_only_ = false;
  /// Every level of every J-group in the R-matrix is itself in it (ClearMatrices).
  bool all_levels_active_ = true;
};

#endif
