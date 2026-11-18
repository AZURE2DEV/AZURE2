#ifndef COVARIANCEBAND_H
#define COVARIANCEBAND_H

#include "AZUREGrad.h"
#include <string>
#include <vector>
#include <map>

class CompoundSystem;
class WorkingCopy;
class Unit;
class Config;

/// Fitted-parameter covariance, with each column tagged by its parameter
/// identity (cols[a]) so a covariance saved from a fit can be matched onto a
/// later run whose parameters are ordered or fixed differently.  M is dense
/// (nPacked x nPacked) in ParamIndexMap packed order.
struct BandCovariance {
  std::vector<ParamDesc> cols;  ///< empty when loaded from a plain covariance.dat
  std::vector<std::vector<double>> M;
  /// Parameter names of the columns (param.sav names), written as the "#columns" header of covariance.dat when they
  /// are known (a fit of this build); empty otherwise.
  std::vector<std::string> names;
  /// Factor the covariance was multiplied by (the reduced chi^2 of --scale-covariance), 1 = unscaled; written as a
  /// "#scaled" header line when not 1.
  double scale = 1.0;
  bool empty() const { return M.empty(); }
  int size() const { return (int)M.size(); }
};

/// A parameter of a target, a sample or an apparatus that changes the model of a point: an isotopic abundance, an
/// electron-screening potential, an areal density, a Doppler temperature, a sample thickness, a parameter of an
/// apparatus record (a stopping scale, a kernel's width).  The band propagates it (the owner's answer of 2026-11-18);
/// its per-point sensitivity is a central difference of the model (AddApparatusSensitivities), not the adjoint's.
bool IsApparatusBandParameter(ParamKind kind);
/// A parameter the band propagates: the R-matrix parameters (IsRMatrixParameter), the THM coherent background
/// (cbackground=, which enters the HOES amplitude) and the apparatus parameters (IsApparatusBandParameter).
/// Normalizations and energy shifts are not: they scale or move the data, not the model curve.
bool IsBandParameter(ParamKind kind);

/// Packed-column indices of a pmap that the band propagates (IsBandParameter), in packed order.  The band covariance
/// and covariance.dat span just these columns -- the R-matrix ones alone for a project without cbackground= and
/// without a free target or apparatus parameter, as before 2026-11-18.
std::vector<int> BandPackedColumns(const ParamIndexMap &pmap);

/*!
 * The sensitivities of the model of every unit of `data` to the free apparatus band parameters of `pmap`
 * (IsApparatusBandParameter), by central differences at the full parameter vector `full` (step 1e-6 (|x| + 1), as the
 * gradient's differences of these parameters, AZURECalcGradient.cpp), written into the columns (packed indices) of
 * `fullGrad`'s rows (keyed by the units of `data`; a missing row is made, NumPacked() wide).  The model is evaluated
 * on copies of `state` and `data`, which are left as they were.  Nothing is done without such a parameter.
 */
void AddApparatusSensitivities(const CompoundSystem *compound, const LevelState &state, const LevelInputs &inputs,
                               WorkingCopy *data, const Config &config, const ParamIndexMap &pmap, const vector_r &full,
                               std::map<Unit *, vector_r> &fullGrad);

/// Densify the R-matrix sub-block of Minuit's packed lower-triangular covariance
/// (Data()) into a BandCovariance, using pmap for the column identities.  Norm
/// and energy-shift columns are dropped.  Empty on size mismatch.
BandCovariance BuildBandCovarianceFromMinuit(const std::vector<double> &covData,
                                             const ParamIndexMap &pmap);

/// Write/read a BandCovariance to/from a plain text file: the N x N matrix, one row
/// per line, after optional '#' header lines -- "#columns <name> ..." (the column
/// parameters' names, since 2026-11-05) and "#scaled <factor>" -- which numpy.loadtxt
/// and an older build skip; the reader keeps the names (cov.names).  False on I/O
/// error.  A loaded covariance has no column identities (cov.cols is empty), so its
/// columns are taken to be in the current packed parameter order.
bool SaveBandCovariance(const std::string &path, const BandCovariance &cov);
/// Read a covariance matrix from file; false if it cannot be read.
bool LoadBandCovariance(const std::string &path, BandCovariance &cov);

/// Project a saved covariance onto pmap's R-matrix columns (BandPackedColumns
/// order).  If the saved covariance carries column identities (same run) they are
/// matched by identity, leaving unmatched entries zero.  If it was loaded from
/// covariance.dat (no identities) its columns are assumed already in R-matrix
/// order and are copied through when the dimension matches the R-matrix count.
std::vector<std::vector<double>> RemapCovarianceToParamMap(const BandCovariance &saved,
                                                           const ParamIndexMap &pmap);

/// Inputs for the cross-section band: covariance M and per-point sensitivities
/// grad[point], both in R-matrix column order (BandPackedColumns; norms and
/// energy shifts excluded).  dXS(g) = sqrt(g^T M g) (SAMMY Eq. IV E4.2); sum
/// gradient rows first to combine correlated points.
struct BandData {
  const CompoundSystem *compound;
  std::vector<std::vector<double>> M;
  std::map<Unit *, std::vector<double>> grad;
  /// One-sigma cross-section uncertainty from a gradient vector, \f$\sqrt{g^T M g}\f$.
  double dXS(const std::vector<double> &g) const;
};

/// Build a BandData from savedCov for the current (best-fit-filled) level state and data: the caller's `state` (the
/// run's, which the final evaluation filled), whose shift functions it computes, and its `inputs`.  `full`: the
/// full parameter vector of the best fit (the apparatus columns are differenced there); null takes the vector the
/// project's parameters hold.  A covariance loaded from covariance.dat whose size is not the band's but whose
/// "#columns" header names every band column of this run (an extrapolation after a fit with a free target
/// parameter: the target's columns are not in this run) is taken on those columns, by name.
/// False if the analytic sensitivities are unavailable for this model.  THM
/// (HOES) points get g = J_m + m (ds/dp)/s from ComputeTHMRows (central
/// differences of the HOES model plus the profiled scale s = 1/n*): the band of
/// the model as it lies against the data scaled by n*.
/// `fixedMask` (Minuit order) is the fixed set of the fit -- including an
/// external param.par's "fixed" flags; null takes the .azr's alone.
bool BuildBandData(const CompoundSystem *compound, LevelState &state, const LevelInputs &inputs, WorkingCopy *data, azure2::Results &results, const Config &config,
                   const BandCovariance &savedCov, BandData &out,
                   const std::vector<bool> *fixedMask = nullptr, const vector_r *full = nullptr);

#endif
