#ifndef AMATRIXFUNC_H
#define AMATRIXFUNC_H

#include "GenMatrixFunc.h"
#include "LevelMatrixSolver.h"

struct GradAccum;

///A function class to calculate the T-Matrix using the A-Matrix

/*!
 * The AMatrixFunc function class calculates the T-Matrix for a given energy point using the compound
 * nucleus object.  The AMatrixFunc class is a child class of GenMatrixFunc, where the cross section is
 * calculated from the T-Matrix.
 */

class AMatrixFunc : public GenMatrixFunc {
 public:
  AMatrixFunc(CNuc*, const Config &configure);
  /*!
   * Returns a pointer to the compound nucleus object.
   */
  CNuc *compound() const {return compound_;};
  /*!
   * Returns a reference to the Config structure.
   */
  const Config &configure() const {return *configure_;};
  /*!
   * Repoints the object at a compound nucleus and configuration so that one
   * instance can serve many energy points without reallocating its scratch
   * storage.  ClearMatrices() must still be called before each point.
   */
  void Reset(CNuc *compound, const Config &configure)
    {compound_=compound; configure_=&configure;};

  void ClearMatrices();
  void FillMatrices(EPoint*);
  /*!
   * Factors the level matrix of every J-group.  The full inverse is no longer
   * formed here: the T-matrix needs only bilinear forms of it, and Inverse() is
   * materialized lazily for the callers that need individual elements.
   */
  void InvertMatrices();
  void CalculateTMatrix(EPoint*);
  /*!
   * Instantiated in the parent class.
   */
  void CalculateCrossSection();

  complex GetAMatrixElement(int,int,int) const;
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
   */
  bool PointAdjoint(EPoint* point, double fitBar, GradAccum& accum,
                    const vector_matrix_r* shiftDeriv = nullptr,
                    int xsComponent = 0);
  void AddAInvMatrixElement(int,int,int,complex);
 private:
  ///One memoized channel-capture bilinear form, keyed by its pathway.
  struct ChanCapEntry {
    int jGroupNum;
    int chNum;
    int chpNum;
    int maskChannel;
    complex value;
  };

  const vector_r &GetGammaVector(int jGroupNum, int chNum);
  const vector_c &GetChannelSolve(int jGroupNum, int chNum);

  const Config *configure_;
  CNuc *compound_;
  vector_matrix_c a_inv_matrices_;
  ///Factored level matrix of each J-group.
  std::vector<LevelMatrixSolver> solvers_;
  ///Whether the corresponding solver holds a usable factorization.
  std::vector<char> solver_valid_;
  // private:
  std::vector<std::vector<int>> level_active_index_; // [jGroup-1][origLevel] -> activeIndex (1..k) or 0 if inactive

  // Per-point caches over active levels, indexed [jGroup-1][channel-1].
  ///Reduced widths of a channel, with the 1e-12 cutoff applied.
  std::vector<std::vector<vector_r>> gamma_vectors_;
  std::vector<std::vector<char>> gamma_valid_;
  ///A times the reduced-width vector of a channel.
  std::vector<std::vector<vector_c>> channel_solves_;
  std::vector<std::vector<char>> channel_solve_valid_;
  std::vector<ChanCapEntry> chan_cap_cache_;
  vector_c chan_cap_entrance_;
  vector_c chan_cap_exit_;

  // Pre-allocated buffers to prevent memory fragmentation
  mutable std::vector<std::vector<double>> levelGammas_;
  mutable std::vector<double> levelEnergies_;
  mutable std::vector<std::vector<double>> shiftFunctions_;
  mutable int cached_max_levels_;
  mutable int cached_max_channels_;

};

#endif
