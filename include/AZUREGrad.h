#ifndef AZUREGRAD_H
#define AZUREGRAD_H

#include "Constants.h"
#include <vector>
#include <map>
#include <tuple>
#include <functional>

class CNuc;
class EData;
class Config;
class EPoint;
class ESegment;
class ParamIndexMap;

/*!
 * \brief Reverse-mode (adjoint) gradient support for the AZURE2 forward model.
 *
 * This module implements the analytic adjoint of the log-likelihood with respect
 * to the RWA fit parameters (level energies E_lambda, reduced-width amplitudes
 * gamma_lambda,c, dataset normalizations n_k, and energy shifts).  See PLAN.md for
 * the full derivation and phase structure.
 *
 * Phase 0 provides the parameter-index map: a faithful mirror of the flat
 * ordering produced by CNuc::FillCompoundFromParams / EData::FillNormsFromParams
 * / EData::FillEnergyShiftsFromParams, so that gradient components computed in
 * "physics" coordinates (jGroup, level, channel) can be scattered into the same
 * flat vector the sampler perturbs.
 */

/// Kinds of fit parameter, in their flat-vector ordering.
enum class ParamKind {
  LevelEnergy,  ///< E_lambda  for a (jGroup, level)
  Gamma,        ///< gamma_lambda,c for a (jGroup, level, channel)
  Norm,         ///< dataset normalization n_k for a segment
  EnergyShift,  ///< per-segment energy shift (constant term)
  EnergyShiftSqrt,  ///< per-segment sqrt(E) energy-shift coefficient
  ThmCoherent   ///< Re or Im of a THM coherent-background amplitude (cbackground=)
};

/// Description of one entry in the full (unfiltered) flat parameter vector.
struct ParamDesc {
  ParamKind kind;
  int jGroup;   ///< 1-based JGroup index   (LevelEnergy/Gamma only, else -1)
  int level;    ///< 1-based level index    (LevelEnergy/Gamma only, else -1)
  int channel;  ///< 1-based channel index  (Gamma; ThmCoherent: its 1-based index in that block; else -1)
  int segment;  ///< 1-based segment index  (Norm/EnergyShift/EnergyShiftSqrt only, else -1)
  bool fixed;   ///< whether this parameter is fixed (excluded from packed vector)
};

/*!
 * \brief Bidirectional map between physics coordinates and the flat parameter
 *        vector used by FillCompoundFromParams / FillNormsFromParams.
 *
 * The flat vector layout is, in order:
 *   - for each JGroup j, for each level la:  E_{j,la}, then gamma_{j,la,ch} for ch=1..NumChannels
 *   - the norm parameters (starting at normOffset)
 *   - the energy-shift parameters (starting at energyShiftOffset)
 *   - the sqrt(E) energy-shift coefficients (starting at energyShiftSqrtOffset)
 *   - the THM coherent-background parameters (cbackground=; none without it)
 *
 * "Full" indices run over every entry above (matching the `all_rwa_` vector in
 * AZUREAPI).  The "packed" index is the position within the non-fixed subset
 * (matching the vector the sampler/optimizer actually passes).
 */
class ParamIndexMap {
 public:
  /// Full index (into all_rwa_) of the energy parameter for (jGroup, level).
  int EIndex(int jGroup, int level) const;
  /// Full index of the gamma parameter for (jGroup, level, channel).
  int GammaIndex(int jGroup, int level, int channel) const;

  /// Full index of the norm parameter for a 1-based segment, or -1 if that
  /// segment does not carry a (varying) norm parameter.
  int NormIndex(int segment) const;
  /// Full index of the energy-shift parameter for a 1-based segment, or -1.
  int EnergyShiftIndex(int segment) const;
  /// Full index of the sqrt(E) energy-shift coefficient for a 1-based segment, or -1.
  int EnergyShiftSqrtIndex(int segment) const;

  /// Number of entries in the full (unfiltered) parameter vector.
  int NumFull() const { return (int)desc_.size(); }
  /// Number of non-fixed (packed) parameters.
  int NumPacked() const { return numPacked_; }

  /// Description of full-index i.
  const ParamDesc &Desc(int fullIndex) const { return desc_[fullIndex]; }

  /// Packed index for a full index, or -1 if that parameter is fixed.
  int FullToPacked(int fullIndex) const { return fullToPacked_[fullIndex]; }
  /// Full index for a packed index.
  int PackedToFull(int packedIndex) const { return packedToFull_[packedIndex]; }

  int NormOffset() const { return normOffset_; }
  int EnergyShiftOffset() const { return energyShiftOffset_; }
  int EnergyShiftSqrtOffset() const { return energyShiftSqrtOffset_; }

  friend ParamIndexMap BuildParamIndexMap(CNuc *, EData *, const std::vector<bool> &);

 private:
  std::vector<ParamDesc> desc_;                          ///< per full index
  std::map<std::pair<int, int>, int> eIndex_;            ///< (j,la) -> full idx
  std::map<std::tuple<int, int, int>, int> gammaIndex_;  ///< (j,la,ch) -> full idx
  std::map<int, int> normIndex_;                         ///< segment -> full idx
  std::map<int, int> shiftIndex_;                        ///< segment -> full idx
  std::map<int, int> sqrtShiftIndex_;                    ///< segment -> full idx
  std::vector<int> fullToPacked_;
  std::vector<int> packedToFull_;
  int numPacked_ = 0;
  int normOffset_ = 0;
  int energyShiftOffset_ = 0;
  int energyShiftSqrtOffset_ = 0;
};

/*!
 * \brief Gradient accumulator in "physics" coordinates.
 *
 * The reverse-mode adjoint naturally produces dlnL/dE and dlnL/dgamma indexed by
 * (jGroup, level [, channel]).  These are accumulated here over all energy
 * points, then scattered once into the flat gradient vector via the index map.
 */
struct GradAccum {
  // e[j-1][la-1]            -> dlnL/dE_{j,la}
  std::vector<std::vector<double>> e;
  // gamma[j-1][la-1][ch-1]  -> dlnL/dgamma_{j,la,ch}
  std::vector<std::vector<std::vector<double>>> gamma;

  /// Allocate (and zero) the accumulators to match the compound's dimensions.
  void Init(CNuc *compound);
  /// Reset all accumulators to zero (without reallocating).
  void Zero();
  /// Add another accumulator (same compound dimensions) into this one.  Used to
  /// reduce per-thread partial gradients after a parallel point loop.
  void Add(const GradAccum &o);
  /// Add an energy contribution (1-based indices).
  void AddE(int jGroup, int level, double v) { e[jGroup - 1][level - 1] += v; }
  /// Add a gamma contribution (1-based indices).
  void AddGamma(int jGroup, int level, int channel, double v) {
    gamma[jGroup - 1][level - 1][channel - 1] += v;
  }
  /// Scatter into the flat (full) gradient vector using the index map.
  void Scatter(const ParamIndexMap &map, vector_r &gradFull) const;
};

/*!
 * \brief Build the parameter-index map by mirroring the fill routines.
 *
 * \param compound  the CNuc whose levels/channels define the R-matrix params.
 * \param data      the EData whose segments define norm/energy-shift params.
 * \param fixed     the per-full-index fixed mask (same as AZUREAPI::fixed_).
 *                  If empty, all parameters are treated as non-fixed.
 */
ParamIndexMap BuildParamIndexMap(CNuc *compound, EData *data,
                                 const std::vector<bool> &fixed);

/*!
 * \brief Precompute the energy derivative of the Brune shift function S(E_lambda)
 *        for every (jGroup, level, channel).
 *
 * Returns a table indexed [j-1][la-1][ch-1] holding dS_ch/dE_level evaluated at
 * the level energy, for particle ('P') channels (zero otherwise).  Mirrors the
 * bound/unbound branch of CNuc::CalcShiftFunctions, using ShftFunc (Whittaker,
 * bound) or CoulFunc::PEShift_dE (unbound).  These depend only on level energies
 * and channel radii, so they are computed once per gradient evaluation and
 * reused across all data points.  See PLAN.md Phase 6.
 */
vector_matrix_r BuildShiftDerivTable(CNuc *compound, const Config &configure);

/*!
 * \brief Gradient of CNuc::ParkNormPenalty with respect to the level energies and
 * reduced widths, added to `accum` in physics coordinates.  Zero unless a level
 * has J < 0; a no-op outside the Park formalism.
 */
void AddParkPenaltyGradient(CNuc *compound, const Config &configure, GradAccum &accum);

/*!
 * \brief Shared analytic gradient engine used by both AZUREAPI (log-likelihood)
 *        and the built-in Minuit2 fitter (chi-squared).
 *
 * Accumulates d(objective)/d(E_lambda, gamma_lambda,c) into the energy/gamma
 * entries of `gradFull` (full RWA parameter layout) via the reverse-mode adjoint
 * of the forward model.  The objective enters *only* through the per-point
 * cotangent `fitBar = d(objective)/d(model)` supplied by `fitBarFn`; returning
 * 0 from `fitBarFn` skips that point.  This is the single place the adjoint is
 * wired up, so future changes to the differentiation propagate to every caller.
 *
 * `compound` / `data` must already be filled from the current parameters.
 * Returns false (touching nothing) if any data point is outside the supported
 * analytic path, so the caller can fall back to finite differences.
 *
 * fitBarFn arguments: (segment, 1-based segment index, 0-based point index,
 * combined model value).
 */
using FitBarFn = std::function<double(ESegment *, int, int, double)>;

bool AccumulateEGammaGradient(CNuc *compound, EData *data, const Config &config,
                              const ParamIndexMap &pmap,
                              const vector_matrix_r *shiftDeriv,
                              const FitBarFn &fitBarFn,
                              GradAccum &accum);

/*!
 * \brief Residual vector and its Jacobian for Gauss-Newton / Levenberg-Marquardt.
 *
 * Builds the standardized residuals r_i = (model_i - data_i*n)/(cmErr_i*n) (so
 * that sum r_i^2 = chi^2) and the Jacobian J_{ij} = d r_i / d theta_j, where the
 * columns theta_j are the non-fixed (packed) parameters.  Each row is obtained
 * by one reverse-mode adjoint of the point with cotangent 1/(cmErr_i*n), so the
 * whole Jacobian costs ~2 forwards regardless of the parameter count.
 *
 * `compound`/`data` must already be filled.  `residuals` gets N_res entries;
 * `jacobian` is row-major N_res x nCols (nCols = number of non-fixed params).
 * Energy shifts are not included as analytic columns here; the API adds them by
 * finite differences on the residual vector this function returns (see
 * AZUREAPI::CalculateResidualJacobianRWA).  Returns false if any point is
 * outside the supported analytic path.
 *
 * THM (HOES) segments are outside it: by default their presence returns false.
 * With skipTHM their rows are still assigned (in segment order, like every other
 * row) but left zero, for a caller that fills them itself -- AZUREAPI does, with
 * the profiled norm (see ComputeTHMRows).
 */
bool ComputeResidualJacobian(CNuc *compound, EData *data, const Config &config,
                             const ParamIndexMap &pmap,
                             const vector_matrix_r *shiftDeriv,
                             vector_r &residuals, vector_r &jacobian, int &nCols,
                             bool skipTHM = false);

/*!
 * \brief Per-point d(model)/d(theta) for analytic covariance bands, keyed by
 *        EPoint*.  Columns are the packed parameters of `pmap`; only the R-matrix
 *        (E, gamma) columns are populated (dT/dn = dT/dshift = 0).  Each row is
 *        one adjoint with cotangent 1.  `compound`/`data` must be best-fit-filled.
 *        Returns false (clearing the output) if any point is outside the analytic
 *        path, so the caller can skip the band.  With skipTHM the rows of THM
 *        (HOES) points are left zero for the caller to fill (the adjoint
 *        differentiates the T-matrix observable, not the HOES one).
 */
bool ComputeModelGradients(CNuc *compound, EData *data, const Config &config,
                           const ParamIndexMap &pmap,
                           const vector_matrix_r *shiftDeriv,
                           std::map<EPoint *, vector_r> &gradByPoint,
                           bool skipTHM = false);

/// THM (HOES) rows of one THM segment: model, residuals and their parameter
/// derivatives with the segment norm profiled (see ComputeTHMRows).
struct THMRows {
  int segment = 0;       ///< 1-based segment index
  int firstRow = 0;      ///< row of its first point in the global residual vector
  vector_r m;            ///< model per point
  vector_r r;            ///< standardized residual per point
  vector_r Jm;           ///< d m / d p, row-major nPoints x nCols (E, gamma columns)
  vector_r J;            ///< d r / d p, row-major nPoints x nCols
  double s = 1.0;        ///< scale on the model, 1/n (n* = S_mm/S_md when profiled)
  bool profiled = false; ///< s follows the parameters (profiled, non-degenerate norm)
  vector_r ds;           ///< d s / d p per packed column (zero unless profiled)
  /// Segments of a THM experiment (<thm> experiment[...], profiled jointly
  /// with a shared norm and background): the band rows d q/d p of the curve
  /// the output shows, q = f(p)/s(p0), nPoints x nCols (ThmProfileDerivative).
  /// Empty otherwise, where the band row is J_m + m (ds/dp)/s.  ds is then zero
  /// and J carries the whole derivative through the profile.
  vector_r G;
};

/*!
 * \brief THM (HOES) rows for every THM segment, at the parameters `full`.
 *
 * THM points are outside the analytic adjoint (it differentiates the T-matrix
 * observable), so their model Jacobian J_m is taken by central differences in
 * the level energies and reduced widths (h = 1e-6 (|x| + 1), the step
 * AZURECalc::Gradient uses for its THM part); norm and energy-shift columns are
 * left zero (a caller differences shifts on the whole residual vector).  The
 * residual of a point of a profiled segment is r_i = (s m_i - d_i)/e_i with
 * s = 1/n* = S_md/S_mm, a function of the parameters; its derivative is taken
 * analytically:
 *   d r_i/dp = (s J_mi + m_i ds/dp)/e_i,
 *   ds/dp    = (sum_k d_k J_mk/e_k^2 - 2 s sum_k m_k J_mk/e_k^2) / S_mm.
 * (A fixed THM norm, or a degenerate profile, has ds/dp = 0.)  The segments
 * of a THM experiment are profiled together, norm and background, and their J
 * (and band rows G) come from ThmProfileDerivative.  `full` is the
 * full (Minuit-ordered) parameter vector of `pmap`.  The compound/data are left
 * filled at `full`, each THM point's fit cross section set to its model and a
 * profiled segment's norm to n*.
 */
std::vector<THMRows> ComputeTHMRows(CNuc *compound, EData *data, const Config &config,
                                    const vector_r &full, const ParamIndexMap &pmap);

#endif
