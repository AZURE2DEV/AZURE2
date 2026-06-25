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
  LevelEnergy,   ///< E_lambda  for a (jGroup, level)
  Gamma,         ///< gamma_lambda,c for a (jGroup, level, channel)
  Norm,          ///< dataset normalization n_k for a segment
  EnergyShift    ///< per-segment energy shift
};

/// Description of one entry in the full (unfiltered) flat parameter vector.
struct ParamDesc {
  ParamKind kind;
  int jGroup;   ///< 1-based JGroup index   (LevelEnergy/Gamma only, else -1)
  int level;    ///< 1-based level index    (LevelEnergy/Gamma only, else -1)
  int channel;  ///< 1-based channel index  (Gamma only, else -1)
  int segment;  ///< 1-based segment index  (Norm/EnergyShift only, else -1)
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

  /// Number of entries in the full (unfiltered) parameter vector.
  int NumFull() const { return (int)desc_.size(); }
  /// Number of non-fixed (packed) parameters.
  int NumPacked() const { return numPacked_; }

  /// Description of full-index i.
  const ParamDesc& Desc(int fullIndex) const { return desc_[fullIndex]; }

  /// Packed index for a full index, or -1 if that parameter is fixed.
  int FullToPacked(int fullIndex) const { return fullToPacked_[fullIndex]; }
  /// Full index for a packed index.
  int PackedToFull(int packedIndex) const { return packedToFull_[packedIndex]; }

  int NormOffset() const { return normOffset_; }
  int EnergyShiftOffset() const { return energyShiftOffset_; }

  friend ParamIndexMap BuildParamIndexMap(CNuc*, EData*, const std::vector<bool>&);

 private:
  std::vector<ParamDesc> desc_;            ///< per full index
  std::map<std::pair<int,int>, int> eIndex_;            ///< (j,la) -> full idx
  std::map<std::tuple<int,int,int>, int> gammaIndex_;   ///< (j,la,ch) -> full idx
  std::map<int, int> normIndex_;                        ///< segment -> full idx
  std::map<int, int> shiftIndex_;                       ///< segment -> full idx
  std::vector<int> fullToPacked_;
  std::vector<int> packedToFull_;
  int numPacked_ = 0;
  int normOffset_ = 0;
  int energyShiftOffset_ = 0;
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
  void Init(CNuc* compound);
  /// Reset all accumulators to zero (without reallocating).
  void Zero();
  /// Add an energy contribution (1-based indices).
  void AddE(int jGroup, int level, double v) { e[jGroup-1][level-1] += v; }
  /// Add a gamma contribution (1-based indices).
  void AddGamma(int jGroup, int level, int channel, double v) {
    gamma[jGroup-1][level-1][channel-1] += v;
  }
  /// Scatter into the flat (full) gradient vector using the index map.
  void Scatter(const ParamIndexMap& map, vector_r& gradFull) const;
};

/*!
 * \brief Build the parameter-index map by mirroring the fill routines.
 *
 * \param compound  the CNuc whose levels/channels define the R-matrix params.
 * \param data      the EData whose segments define norm/energy-shift params.
 * \param fixed     the per-full-index fixed mask (same as AZUREAPI::fixed_).
 *                  If empty, all parameters are treated as non-fixed.
 */
ParamIndexMap BuildParamIndexMap(CNuc* compound, EData* data,
                                 const std::vector<bool>& fixed);

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
vector_matrix_r BuildShiftDerivTable(CNuc* compound, const Config& configure);

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
using FitBarFn = std::function<double(ESegment*, int, int, double)>;

bool AccumulateEGammaGradient(CNuc* compound, EData* data, const Config& config,
                              const ParamIndexMap& pmap,
                              const vector_matrix_r* shiftDeriv,
                              const FitBarFn& fitBarFn,
                              GradAccum& accum);

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
 * Energy shifts (rare) are not included as analytic columns here.  Returns false
 * if any point is outside the supported analytic path.
 */
bool ComputeResidualJacobian(CNuc* compound, EData* data, const Config& config,
                             const ParamIndexMap& pmap,
                             const vector_matrix_r* shiftDeriv,
                             vector_r& residuals, vector_r& jacobian, int& nCols);

#endif
