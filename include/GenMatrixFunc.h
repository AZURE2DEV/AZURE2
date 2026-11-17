#ifndef GENMATRIXFUNC_H
#define GENMATRIXFUNC_H

#include "Constants.h"
#include "Config.h"
#include "Site.h"
#include "LevelState.h"
#include "LineWindow.h"
#include "azure2/physics/Observables.h"

#include <cstring>
#include <string>
#include <vector>

class CompoundSystem;
namespace azure2 {
struct Results;
struct EntrancePathways;
struct ExitRow;
}

/// A generalized function class to calculate cross sections

/*!
 * The GenMatrixFunc function class is the general form of the function used
 * to calculate cross section from R-Matrix parameters.  It is the parent class of
 * AMatrixFunc and RMatrixFunc.
 */


class GenMatrixFunc {
 public:
  GenMatrixFunc() {};
  virtual ~GenMatrixFunc() {};
  /*!
   *This virtual function in implemented in the child class.
   */
  virtual void ClearMatrices() = 0;
  /*!
   *This virtual function in implemented in the child class.
   */
  virtual void FillMatrices(Site *) = 0;
  /*!
   *This virtual function in implemented in the child class.
   */
  virtual void InvertMatrices() = 0;
  /*!
   *This virtual function in implemented in the child class.
   */
  virtual void CalculateTMatrix(Site *) = 0;
  /// The T-matrix of another exit pair at the same energy: drops the elements held now (and what was derived from
  /// them) and builds those of `point`'s exit from the solved level matrix, which does not depend on the exit.  Used
  /// for the T-matrix groups of an energy task after the first (SolveAtOneEnergy).
  void RecalculateTMatrix(Site *point) {
    ResetTMatrix();
    CalculateTMatrix(point);
  }
  /// A node of the line-shape integration of a broad capture final state (LineWindow.h), read by CalculateTMatrix
  /// while LineWindowCrossSection runs; null otherwise.
  const azure2::LineNode *ActiveLineNode() const { return lineNode_; }
  void SetActiveLineNode(const azure2::LineNode *node) { lineNode_ = node; }
  /// The exit channel's sqrt(P) for CalculateTMatrix: the point's own; at a node of the line-shape integration, for a
  /// photon channel, that of the photon at E - E_x(E_f) times the final-state amplitude A(E_f), zero below E_x
  /// (LineWindow.cpp).
  double ExitSqrtPenetrability(Site *point, int jNum, int chNum) const;
  /// The external-capture T element (k, m) for CalculateTMatrix: `element`; at a node the node's own, that of the
  /// scattering state at E_f (LineWindow.cpp).
  complex ExternalElementAt(int k, int m, const complex &element) const;
  /// True when `other` holds the same T-matrix bit for bit, internal and external elements (the shared-solve check
  /// of AZURE_PLAN_CHECK, CalculationPlan.h).
  bool SameTMatrix(const GenMatrixFunc &other) const {
    return SameBits(tmatrix_, other.tmatrix_) && SameBits(ec_tmatrix_, other.ec_tmatrix_);
  }
  /// The observable of the point from the T-matrix the subclass built, times the electron-screening factor where the
  /// point carries one (CalculateBareCrossSection, then Unit::ScreeningFactor).
  void CalculateCrossSection(Site *, azure2::Results &results);
  /// The observable without electron screening.
  void CalculateBareCrossSection(Site *, azure2::Results &results);

  /// Store an internal T-matrix element for a resonant pathway.
  void AddTMatrixElement(int, int, complex, int decayNum = 1);
  /// Store an external T-matrix element for an external-capture pathway.
  void AddECTMatrixElement(int, int, complex);
  /// Internal T-matrix element for a resonant pathway.
  complex GetTMatrixElement(int, int, int decayNum = 1) const;
  /// External T-matrix element for an external-capture pathway.
  complex GetECTMatrixElement(int, int) const;
  /// Racah coefficient \f$R_k\f$ of the angular-distribution algebra.
  double GetRk(double, double, double, double, int);
  /// The secondary gamma ray's R_k of a UPOS point for the orders below `numOrders` (mixing ratio included; 0 for
  /// odd k), for a decaying state of spin j2f.
  std::vector<double> SecondaryRk(const Site *point, double j2f, int numOrders);
  /// The refusal of a point whose exit pair no pathway of its entrance pair reaches (pair numbers, 1-based).
  std::string NoPathwayMessage(int entrancePair, int exitPair) const;
  /// The refusal of a capture-A_y point whose terms were not requested at set-up (CompoundSetUpRequests::AddCaptureAy).
  std::string CaptureAyNotBuiltMessage(const Site *point) const;

  /*!
   *This virtual function in implemented in the child class.
   */
  virtual const CompoundSystem *compound() const = 0;
  /// The levels at the evaluated parameter vector (LevelState.h), held by the evaluation that made this object.
  virtual const LevelState &state() const = 0;
  /*!
   *This virtual function in implemented in the child class.
   */
  virtual const Config &configure() const = 0;

 protected:
  /// Vector of internal T-matrix elements accessable to child class.
  std::vector<matrix_c> tmatrix_;
  /// Vector of external T-matrix elements accessable to child class.
  matrix_c ec_tmatrix_;
  static bool SameBits(const vector_c &a, const vector_c &b) {
    return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(complex)) == 0);
  }
  template <class V>
  static bool SameBits(const std::vector<V> &a, const std::vector<V> &b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); i++)
      if (!SameBits(a[i], b[i])) return false;
    return true;
  }
  /// Forget everything derived from the T-matrix (ClearMatrices).
  void InvalidateTMatrixProducts() { rt_decay_ = nullptr; }
  /// Drop the T-matrix elements and their products, keeping the solved matrices (RecalculateTMatrix).
  virtual void ResetTMatrix() {
    tmatrix_.clear();
    ec_tmatrix_.clear();
    InvalidateTMatrixProducts();
  }

 private:
  /// The interference products Z1Z2 T1 conj(T2) of the differential sum,
  /// summed by Legendre order (rt_wl_[L]), for the T-matrix held now: they do
  /// not depend on the angle, so the points that share one T-matrix
  /// (the units of one T-matrix group of an energy task) compute them once, and each point's sum is one term per
  /// order.  Entrance = exit pair, capture, and particle exits other than
  /// UPOS (whose sum keeps only the k-groups of equal channel spins).
  const azure2::ExitRow *rt_decay_ = nullptr;
  std::vector<complex> rt_wl_;
  /// The Coulomb-nuclear interference per (channel spin, l), BuildCoulombInterference.
  azure2::physics::CoulombInterference rt_it_;
  void BuildCoulombInterference(const azure2::EntrancePathways &P, const azure2::ExitRow &exit, bool elastic);
  /// rt_wl_ (and for elastic scattering rt_it_) for the T-matrix held now, once per T-matrix.
  void BuildOrderProducts(const azure2::EntrancePathways &P, const azure2::ExitRow &exit, bool particleExit, bool elastic);
  /// The record of an analysing power at the T-matrix held now (SiteObservable.h, ay_record): the cross section's
  /// coefficients and A_y's -- the amplitude matrix's per exit orbital momentum (Polarization::AmplitudeMatrix::
  /// AddPathwayCoefficients) for a particle exit, Seyler and Weller's per order for capture.
  Site::AnalyzingPowerCoefficients AnalyzingPowerCoefficientsOf(const Site *point, const azure2::EntrancePathways &P,
                                                                const azure2::ExitRow &theDecay, int aa,
                                                                int exitPairNum);
  /// The observables of elastic scattering of `point`'s pair `aa` at its energy (SiteObservable::kElastic).
  Site::ElasticCoefficients ElasticCoefficientsOf(const Site *point, const azure2::EntrancePathways &P,
                                                  const azure2::ExitRow &exit, int aa);
  /// The Legendre coefficients B_L (b/sr) of dsigma/dOmega_c.m. of `point`'s reaction at its energy (for elastic
  /// scattering the nuclear amplitude's part).
  std::vector<double> LegendreCoefficientsOf(const Site *point, const azure2::EntrancePathways &P,
                                             const azure2::ExitRow &exit, int aa, int exitPairNum);
  /// The coherent sums by (J, l, l') of the spin group at hand (azure2/physics/Observables.h): all pathways, and
  /// those of E1 and of E2 radiation kept apart for the reported components.
  azure2::physics::PartialWaveSums partial_sums_;
  azure2::physics::PartialWaveSums partial_sums_E1_;
  azure2::physics::PartialWaveSums partial_sums_E2_;
  /// ActiveLineNode.
  const azure2::LineNode *lineNode_ = nullptr;
};

#endif
