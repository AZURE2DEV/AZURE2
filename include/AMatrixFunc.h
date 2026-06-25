#ifndef AMATRIXFUNC_H
#define AMATRIXFUNC_H

#include "GenMatrixFunc.h"

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
  const Config &configure() const {return configure_;};
  
  void ClearMatrices();
  void FillMatrices(EPoint*);
  void InvertMatrices();
  void CalculateTMatrix(EPoint*);
  /*!
   * Instantiated in the parent class.
   */
  void CalculateCrossSection();

  complex GetAMatrixElement(int,int,int) const;
  matrix_c *GetJSpecAInvMatrix(int);

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
  void AddAMatrix(matrix_c);
  void AddAMatrix(matrix_c&&);
 private:
  const Config &configure_;
  CNuc *compound_;
  vector_matrix_c a_inv_matrices_;
  vector_matrix_c a_matrices_;
  int a_matrices_index_;  // Index for thread-safe AddAMatrix calls
  // private:
  std::vector<std::vector<int>> level_active_index_; // [jGroup-1][origLevel] -> activeIndex (1..k) or 0 if inactive
  
  // Pre-allocated buffers to prevent memory fragmentation
  mutable std::vector<std::vector<double>> levelGammas_;
  mutable std::vector<double> levelEnergies_;
  mutable std::vector<std::vector<double>> shiftFunctions_;
  mutable int cached_max_levels_;
  mutable int cached_max_channels_;

};

#endif
