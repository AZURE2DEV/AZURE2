#ifndef AZURE2_INFERENCE_ABUNDANCEPRIORS_H
#define AZURE2_INFERENCE_ABUNDANCEPRIORS_H

#include <vector>

#include "azure2/inference/GaussianPrior.h"

namespace azure2 {
namespace inference {

/*!
 * The abundance term of the objective (T7; docs/upgrade/objective_arrays.md section 2) as arithmetic on arrays.
 *
 * A project's <abundances> block gives, per sample, one isotopic (or stoichiometric) abundance w per compound
 * nucleus; all are dimensionless.  The abundances are numbered by the ordinal k = (s - 1) N + i (s the sample, i the
 * nucleus, both 1-based, N the number of nuclei), and the arrays below are indexed by k - 1.  The term is, sample by
 * sample,
 *   sum_i ((w_k - w0_k)/sigma_k)^2   over the nuclei whose prior is active, then
 *   ((sum_i w_k) - c_s)^2 / s_s^2    when s_s > 0 (the sum-to-one term; c_s its centre, 1 unless written),
 * added into one sum from 0 in that order.  sum_i w_k is formed from 0 in nucleus order.
 *
 * Which prior is active (sigma > 0 and not replaced by the prior covariance) and which width a sample's sum has (its
 * own "sum" line, else the config's abundanceSumSigma) is the engine's gather (AbundanceBlockOf,
 * include/DiagonalPriors.h): this header holds no rule of eligibility.  The per-nucleus priors are also entries of
 * the list of diagonal priors (DiagonalPriors.h), from the same gather.
 */
struct AbundanceBlock {
  int numSamples = 0;
  int numNuclei = 0;
  std::vector<double> w;       ///< the abundances, per ordinal (k - 1)
  std::vector<double> centre;  ///< the prior's centre w0, per ordinal
  std::vector<double> sigma;   ///< the prior's width, per ordinal; read only where active
  std::vector<char> active;    ///< the nucleus's prior enters the term, per ordinal
  std::vector<double> sumCentre;  ///< c_s, per sample (s - 1)
  std::vector<double> sumSigma;   ///< s_s, per sample; no sum term unless > 0
};

/// sum_i w_k of sample s (1-based), from 0 in nucleus order.
inline double AbundanceSum(const AbundanceBlock &b, int s) {
  double sum = 0.0;
  for (int i = 0; i < b.numNuclei; i++) sum += b.w[(s - 1) * b.numNuclei + i];
  return sum;
}

/// The term from 0: with `withPriors` the per-nucleus priors and the sum terms, sample by sample (the objective's
/// value); without, the sum terms only (the sampler's part: the per-nucleus priors are its log-priors).
inline double AbundanceChi2(const AbundanceBlock &b, bool withPriors) {
  double chi2 = 0.0;
  for (int s = 1; s <= b.numSamples; s++) {
    for (int i = 0; withPriors && i < b.numNuclei; i++) {
      const int k = (s - 1) * b.numNuclei + i;
      if (b.active[k]) chi2 += GaussianPriorChi2(b.w[k], b.centre[k], b.sigma[k]);
    }
    const double sig = b.sumSigma[s - 1];
    if (sig > 0.0) chi2 += GaussianPriorChi2(AbundanceSum(b, s), b.sumCentre[s - 1], sig);
  }
  return chi2;
}

/// d chi2/d w_k of the whole term, per ordinal into dw (resized): the sum term's derivative 2 ((sum w) - c)/s^2
/// first, then the nucleus's prior's 2 (w - w0)/sigma^2 where active.
inline void AbundanceGradient(const AbundanceBlock &b, std::vector<double> &dw) {
  dw.assign(b.w.size(), 0.0);
  for (int k = 0; k < (int)b.w.size(); k++) {
    const int s = k / b.numNuclei + 1;
    const double sig = b.sumSigma[s - 1];
    double g = (sig > 0.0) ? GaussianPriorGradient(AbundanceSum(b, s), b.sumCentre[s - 1], sig) : 0.0;
    if (b.active[k]) g += GaussianPriorGradient(b.w[k], b.centre[k], b.sigma[k]);
    dw[k] = g;
  }
}

/// One least-squares row of a sum term: r = ((sum w) - c)/s, and dr/dw_k = 1/s on every abundance of the sample.
struct AbundanceSumRow {
  int sample = 0;           ///< 1-based
  double residual = 0.0;    ///< ((sum w) - c)/s
  double derivative = 0.0;  ///< 1/s
};

/// The rows of the sum terms, one per sample with s > 0, in sample order (the per-nucleus priors' rows are the
/// diagonal priors').  The squares of the residuals add up to AbundanceChi2(b, false).
inline std::vector<AbundanceSumRow> AbundanceSumRows(const AbundanceBlock &b) {
  std::vector<AbundanceSumRow> rows;
  for (int s = 1; s <= b.numSamples; s++) {
    const double sig = b.sumSigma[s - 1];
    if (!(sig > 0.0)) continue;
    rows.push_back({s, (AbundanceSum(b, s) - b.sumCentre[s - 1]) / sig, 1.0 / sig});
  }
  return rows;
}

}  // namespace inference
}  // namespace azure2

#endif  // AZURE2_INFERENCE_ABUNDANCEPRIORS_H
