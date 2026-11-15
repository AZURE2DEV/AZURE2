// The Levenberg-Marquardt fitter on the residual rows and their analytic Jacobian, and what it shares with the GSL
// trust-region fitter (LeastSquares.h): the free parameters with their priors, the covariance (J^T J + priors)^-1, the
// row driver with its MIGRAD fallback (docs/upgrade/fit.md).
#include "LeastSquares.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "AZURECalc.h"
#include "CompoundSystem.h"
#include "AZUREGrad.h"
#include "AZUREParams.h"
#include "CovarianceBand.h"
#include "Fitters.h"
#include "Objective.h"
#include "ParameterKind.h"
#include "ParameterLabel.h"
#include "WorkingCopy.h"

int PrepareFreeParams(const AZURECalc &self, const vector_r &full, const ROOT::Minuit2::MnUserParameters &mp,
                      vector_r &r, vector_r &J, std::vector<int> &p2f, vector_r &x, std::vector<double> &lo,
                      std::vector<double> &hi, std::vector<double> &pen_nom, std::vector<double> &pen_inv2) {
  // First Jacobian fixes the column <-> parameter mapping.
  if (!self.ResidualJacobian(full, r, J, p2f)) return -1;  // unsupported -> caller falls back
  const int nFree = (int)p2f.size();

  // Free-parameter values, limits, and quadratic-penalty (prior) terms.
  x.assign(nFree, 0.0);
  lo.assign(nFree, -std::numeric_limits<double>::infinity());
  hi.assign(nFree, std::numeric_limits<double>::infinity());
  pen_nom.assign(nFree, 0.0);
  pen_inv2.assign(nFree, 0.0);
  // The diagonal priors by full parameter index (p2f holds full indices of the same map; Objective.h).
  const ParamIndexMap pmap = BuildParamIndexMap(self.compound(), self.data(), std::vector<bool>());
  std::map<int, std::pair<double, double>> prior;  // full index -> (nominal, sigma)
  ForEachDiagonalPrior(self.data(), pmap,
                       [&](int idx, double nominal, double sig) { prior[idx] = {nominal, sig}; });
  for (int a = 0; a < nFree; a++) {
    int f = p2f[a];
    x[a] = full[f];
    const auto &par = mp.Parameter(f);
    if (par.HasLowerLimit()) lo[a] = par.LowerLimit();
    if (par.HasUpperLimit()) hi[a] = par.UpperLimit();
    // Gaussian penalties contribute (1/sigma^2) to the diagonal of J^T J and (x-nominal)/sigma^2 to the gradient.
    auto it = prior.find(f);
    if (it != prior.end()) {
      pen_nom[a] = it->second.first;
      pen_inv2[a] = 1.0 / (it->second.second * it->second.second);
    } else if (f >= pmap.NumFull() || IsLevelSlot(pmap.Desc(f).kind) || pmap.Desc(f).kind == ParamKind::ThmCoherent) {
      // Nuisance priors (levels and THM coherent background only; none without a limits manager).  Looked up by the full-vector index f, not by
      // a: a counts the parameters free in this Jacobian, the limits manager's non-fixed index counts those free
      // when it was applied plus every fixed segment parameter, and the two need not agree.
      double nom, sig;
      if (self.NuisancePriorOf(f, nom, sig)) {
        pen_nom[a] = nom;
        pen_inv2[a] = 1.0 / (sig * sig);
      }
    }
  }
  return nFree;
}

void FinalizeLeastSquaresCovariance(const AZURECalc &self, const vector_r &full, const std::vector<int> &p2f,
                                    const std::vector<double> &pen_inv2, ROOT::Minuit2::MnUserParameters &mp,
                                    LeastSquaresOutcome &out, std::ostream &log) {
  const int nFree = (int)p2f.size();
  vector_r r, J;
  std::vector<int> p2fLocal;
  if (!self.ResidualJacobian(full, r, J, p2fLocal)) return;

  std::vector<double> H(nFree * nFree, 0.0);
  const int nRes = (int)r.size();
  for (int i = 0; i < nRes; i++) {
    const double *Ji = &J[(size_t)i * nFree];
    for (int a = 0; a < nFree; a++) {
      double *Ha = &H[(size_t)a * nFree];
      for (int b = a; b < nFree; b++) Ha[b] += Ji[a] * Ji[b];
    }
  }
  // Pure J^T J diagonal = each parameter's data sensitivity (column-norm^2 in
  // J), captured before priors are folded into H.  Near zero => unconstrained.
  std::vector<double> jtjDiag(nFree);
  for (int a = 0; a < nFree; a++) {
    for (int b = 0; b < a; b++) H[(size_t)a * nFree + b] = H[(size_t)b * nFree + a];
    jtjDiag[a] = H[(size_t)a * nFree + a];
    if (pen_inv2[a] != 0.0) H[(size_t)a * nFree + a] += pen_inv2[a];
  }

  // The Minuit name alone ("width_1_2") does not say which level or channel is
  // meant, so pair it with the level's J^pi and energy and the channel's pair,
  // L and S.
  auto describe = [&](int a) {
    const int minuitIndex = p2f[a];
    std::string text = mp.GetName(minuitIndex);
    const std::string label = AZURELabel::Parameter(self.compound(), self.inputs(), self.data(), minuitIndex);
    if (!label.empty()) text += "  =  " + label;
    return text;
  };

  // Name the free parameters the data barely constrain (near-zero J column).
  auto reportInsensitive = [&]() {
    double mx = 0.0;
    for (int a = 0; a < nFree; a++) mx = std::max(mx, jtjDiag[a]);
    int nInsensitive = 0;
    for (int a = 0; a < nFree; a++) {
      if (mx <= 0.0 || jtjDiag[a] <= 1.e-10 * mx) {
        if (!nInsensitive)
          log << "  Data-insensitive parameter(s) -- fix, constrain, or add "
                                   "a prior to one of these:"
                                << std::endl;
        log << "    " << describe(a) << std::endl;
        nInsensitive++;
      }
    }
    if (!nInsensitive)
      log << "No single parameter is insensitive; the flat direction is a "
                               "degeneracy between two or more correlated parameters."
                            << std::endl;
  };

  // Covariance = (J^T J + priors)^{-1}.  If it is singular (unconstrained or
  // degenerate parameter), retry with a small ridge so the flat direction gets
  // a large finite variance instead of losing the errors and band; warn.
  double maxDiag = 0.0;
  for (int a = 0; a < nFree; a++) maxDiag = std::max(maxDiag, H[(size_t)a * nFree + a]);
  const double ridges[4] = {0.0, maxDiag * 1.e-9, maxDiag * 1.e-6, maxDiag * 1.e-3};
  gsl_matrix *A = gsl_matrix_alloc(nFree, nFree);
  gsl_set_error_handler_off();
  int used = -1;
  for (int t = 0; t < 4; t++) {
    for (int a = 0; a < nFree; a++)
      for (int b = 0; b < nFree; b++)
        gsl_matrix_set(A, a, b, H[(size_t)a * nFree + b] + (a == b ? ridges[t] : 0.0));
    if (gsl_linalg_cholesky_decomp1(A) == 0 && gsl_linalg_cholesky_invert(A) == 0) {
      used = t;
      break;
    }
  }
  if (used < 0) {
    log << "Warning: the fit curvature matrix (J^T J) is singular; parameter errors "
                             "and the uncertainty band are unavailable."
                          << std::endl;
    reportInsensitive();
  } else {
    if (used > 0) {
      log << "Warning: the fit curvature matrix was near-singular; regularized it "
                               "with a small ridge to obtain the covariance. Uncertainties along the "
                               "poorly-constrained direction(s) are large and approximate."
                            << std::endl;
      reportInsensitive();
      // Strongly correlated pairs (|rho| near 1) are the degenerate combinations.
      struct CorrPair {
        double c;
        int a, b;
      };
      std::vector<CorrPair> pairs;
      for (int a = 0; a < nFree; a++) {
        double va = gsl_matrix_get(A, a, a);
        if (va <= 0.0) continue;
        for (int b = a + 1; b < nFree; b++) {
          double vb = gsl_matrix_get(A, b, b);
          if (vb <= 0.0) continue;
          double c = gsl_matrix_get(A, a, b) / std::sqrt(va * vb);
          if (std::fabs(c) > 0.95) pairs.push_back({c, a, b});
        }
      }
      std::sort(pairs.begin(), pairs.end(),
                [](const CorrPair &x, const CorrPair &y) { return std::fabs(x.c) > std::fabs(y.c); });
      for (size_t k = 0; k < pairs.size() && k < 5; k++)
        log << "  Correlated (rho=" << pairs[k].c << "):" << std::endl
                              << "    " << describe(pairs[k].a) << std::endl
                              << "    " << describe(pairs[k].b) << std::endl;
    }
    for (int a = 0; a < nFree; a++) {
      double v = gsl_matrix_get(A, a, a);
      if (v > 0.0) mp.SetError(p2f[a], std::sqrt(v));
    }
    // Name free parameters whose error dwarfs their value (they dominate the band).
    {
      int n = 0;
      for (int a = 0; a < nFree; a++) {
        double v = gsl_matrix_get(A, a, a);
        if (v <= 0.0) continue;
        double val = mp.Value(p2f[a]);
        if (std::sqrt(v) > 10.0 * std::max(std::fabs(val), 1.e-30)) {
          if (!n)
            log << "Note: weakly-determined parameter(s) (error > 1000% of value) "
                                     "dominate the uncertainty band. Constrain or fix these to "
                                     "tighten it:"
                                  << std::endl;
          log << "    " << describe(a) << std::endl;
          n++;
        }
      }
    }
    // The covariance of the free parameters, in the order of p2f.
    out.freeIndices = p2f;
    out.covariance.assign((size_t)nFree * nFree, 0.0);
    for (int a = 0; a < nFree; a++)
      for (int b = 0; b < nFree; b++) out.covariance[(size_t)a * nFree + b] = gsl_matrix_get(A, a, b);
    // Export the full covariance for cross-section bands.  Columns are the
    // packed free parameters, in the same order as p2f, so we tag them with
    // their parameter identities via a matching index map.
    {
      BandCovariance *bandCovOut = &out.band;
      ParamIndexMap pmap = BuildParamIndexMap(self.compound(), self.data(), self.FixedMask());
      if (pmap.NumPacked() == nFree) {
        // Keep only the band columns (R-matrix, THM coherent background): the
        // band is insensitive to norms and energy shifts, so they are dropped
        // from the saved covariance.
        const std::vector<int> rc = BandPackedColumns(pmap);
        const int m = (int)rc.size();
        bandCovOut->cols.resize(m);
        for (int a = 0; a < m; a++) bandCovOut->cols[a] = pmap.Desc(pmap.PackedToFull(rc[a]));
        bandCovOut->M.assign(m, std::vector<double>(m, 0.0));
        for (int a = 0; a < m; a++)
          for (int b = 0; b < m; b++) bandCovOut->M[a][b] = gsl_matrix_get(A, rc[a], rc[b]);
      }
    }
  }
  gsl_matrix_free(A);
}

LeastSquaresOutcome RunLevenbergMarquardt(AZURECalc &self, AZUREParams &params, int maxIter, const char *fitter,
                                          const FitControl &control) {
  std::ostream &log = *control.log;
  LeastSquaresOutcome out;
  ROOT::Minuit2::MnUserParameters &mp = params.GetMinuitParams();
  const int nFull = mp.Params().size();
  // Free exactly what MIGRAD would: the fixed flags of these parameters
  // (param.par "fixed", <parameterSettings>), not a fresh mask from the .azr.
  {
    std::vector<bool> fixed(nFull);
    for (int i = 0; i < nFull; i++) fixed[i] = mp.Parameter(i).IsFixed();
    self.SetFixedMask(fixed);
  }
  vector_r full(nFull);
  for (int i = 0; i < nFull; i++) full[i] = mp.Value(i);

  // First Jacobian fixes the column <-> parameter mapping; also collects the
  // free-parameter start values, projection limits, and Gaussian priors.
  vector_r r, J, x;
  std::vector<int> p2f;
  std::vector<double> lo, hi, pen_nom, pen_inv2;
  const int nFree = PrepareFreeParams(self, full, mp, r, J, p2f, x, lo, hi, pen_nom, pen_inv2);
  out.calls = 1;
  if (nFree < 0) return out;  // analytic Jacobian unsupported -> the caller falls back
  if (nFree == 0) {
    out.objective = self.Chi2Value(full);
    out.calls++;
    return out;
  }

  double cost = self.Chi2Value(full);
  out.calls++;
  const auto stopped = [&]() { return control.stop && control.stop->load(std::memory_order_relaxed); };
  out.status = "iteration_limit";
  double lambda = 1.0e-3;

  std::vector<double> H(nFree * nFree), g(nFree), Hl(nFree * nFree), dx(nFree);
  std::vector<double> gScale(nFree, 1.0);  // see the convergence test
  for (int a = 0; a < nFree; a++)
    if (params.Kind(p2f[a]) == ParamKind::Density) gScale[a] = std::max(std::fabs(x[a]), 1.0);

  for (int iter = 0; iter < maxIter; iter++) {
    if (stopped()) {
      out.status = "stopped";
      break;
    }
    if (iter > 0) out.calls++;
    if (iter > 0 && !self.ResidualJacobian(full, r, J, p2f)) {
      out.status = "failed_jacobian";
      break;
    }
    // This path evaluates through the side-effect-free Chi2Value, so it calls the checkpoint (the program's
    // intermediate files) itself.
    if (iter > 0 && iter % AZURECalc::kIterationOutputInterval == 0 && control.checkpoint) control.checkpoint(full);
    const int nRes = (int)r.size();

    // Gauss-Newton normal equations: H = J^T J (+ penalty diag), g = J^T r (+ penalty grad).
    std::fill(H.begin(), H.end(), 0.0);
    std::fill(g.begin(), g.end(), 0.0);
    for (int i = 0; i < nRes; i++) {
      const double *Ji = &J[(size_t)i * nFree];
      double ri = r[i];
      for (int a = 0; a < nFree; a++) {
        g[a] += Ji[a] * ri;
        double Jia = Ji[a];
        double *Ha = &H[(size_t)a * nFree];
        for (int b = a; b < nFree; b++) Ha[b] += Jia * Ji[b];
      }
    }
    for (int a = 0; a < nFree; a++) {
      for (int b = 0; b < a; b++) H[(size_t)a * nFree + b] = H[(size_t)b * nFree + a];  // symmetrize
      if (pen_inv2[a] != 0.0) {
        H[(size_t)a * nFree + a] += pen_inv2[a];
        g[a] += (x[a] - pen_nom[a]) * pen_inv2[a];
      }
    }

    // Convergence on the (projected) gradient.  A fitted target density (atoms/cm^2, ~1e17) is judged per
    // relative change, |g| N: per atom/cm^2 its gradient is ~1e-16 and the absolute test stopped LM before
    // its first step.  Every other parameter keeps the absolute test.
    double gmax = 0.0;
    for (int a = 0; a < nFree; a++) gmax = std::max(gmax, std::fabs(g[a]) * gScale[a]);
    if (gmax < 1.0e-8) {
      out.status = "converged";
      break;
    }

    // Inner loop: damp until a step decreases the cost.
    bool accepted = false;
    for (int tries = 0; tries < 40 && !accepted && !stopped(); tries++) {
      // (H + lambda*diag(H)) dx = -g
      Hl = H;
      for (int a = 0; a < nFree; a++) Hl[(size_t)a * nFree + a] += lambda * H[(size_t)a * nFree + a];

      gsl_matrix *A = gsl_matrix_alloc(nFree, nFree);
      for (int a = 0; a < nFree; a++)
        for (int b = 0; b < nFree; b++) gsl_matrix_set(A, a, b, Hl[(size_t)a * nFree + b]);
      gsl_vector *bvec = gsl_vector_alloc(nFree);
      for (int a = 0; a < nFree; a++) gsl_vector_set(bvec, a, -g[a]);
      gsl_vector *sol = gsl_vector_alloc(nFree);

      gsl_set_error_handler_off();
      int status = gsl_linalg_cholesky_decomp1(A);
      if (status == 0) status = gsl_linalg_cholesky_solve(A, bvec, sol);

      bool solved = (status == 0);
      if (solved)
        for (int a = 0; a < nFree; a++) dx[a] = gsl_vector_get(sol, a);
      gsl_matrix_free(A);
      gsl_vector_free(bvec);
      gsl_vector_free(sol);

      if (!solved) {
        lambda *= 4.0;
        if (lambda > 1.0e14) {
          accepted = false;
          break;
        }
        continue;
      }

      vector_r xnew(nFree), fullNew = full;
      for (int a = 0; a < nFree; a++) {
        xnew[a] = std::min(std::max(x[a] + dx[a], lo[a]), hi[a]);  // clamp to limits
        fullNew[p2f[a]] = xnew[a];
      }
      double costNew = self.Chi2Value(fullNew);
      out.calls++;

      if (costNew < cost && std::isfinite(costNew)) {
        double rel = (cost - costNew) / std::max(std::fabs(cost), 1.0);
        full = fullNew;
        x = xnew;
        cost = costNew;
        lambda = std::max(lambda * 0.3, 1.0e-12);
        accepted = true;
        if (control.progress) {
          FitProgress step;
          step.kind = FitProgress::Kind::kIteration;
          step.fitter = fitter;
          step.label = "LM iteration";
          step.iteration = iter + 1;
          step.calls = out.calls;
          step.objective = cost;
          step.best = cost;
          step.measure = lambda;  // the damping after the step
          step.x = &full;
          control.progress(step);
        }
        if (rel < 1.0e-9) {
          out.status = "converged";
          iter = maxIter;
        }
      } else {
        lambda *= 4.0;
        if (lambda > 1.0e14) break;
      }
    }
    if (!accepted) {
      if (out.status != "stopped") out.status = stopped() ? "stopped" : "no_decrease";
      break;
    }
  }
  log << std::endl;

  // Write best-fit values back, and parameter errors from (J^T J + penalties)^{-1}.
  for (int i = 0; i < nFull; i++) mp.SetValue(i, full[i]);
  if (out.status != "stopped") FinalizeLeastSquaresCovariance(self, full, p2f, pen_inv2, mp, out, log);
  out.objective = cost;
  return out;
}


FitResult FitByLeastSquares(AZURECalc &objective, AZUREParams &params, const FitControl &control, const char *fitter,
                            LeastSquaresSolver solver) {
  FitResult result;
  result.fitter = fitter;
  const LeastSquaresOutcome ls = solver(objective, params, 200, fitter, control);
  if (ls.objective < 0.0) {
    // The analytic Jacobian is not available for this model -> MIGRAD (numerical; no HESSE, no MINOS).
    const std::string notice = std::string(fitter) +
                               ": the analytic Jacobian is unsupported for this model; falling back to MIGRAD "
                               "(numerical gradient; its covariance, no HESSE, no MINOS).";
    *control.log << "Analytic Jacobian unsupported for this model; "
                    "falling back to MIGRAD."
                 << std::endl;
    result = MigradPass(objective, params.GetMinuitParams(), false, control);
    AdoptMigradEndPoint(params, result);
    result.fitter = fitter;
    result.fellBackToMigrad = true;
    result.notices.insert(result.notices.begin(), notice);
    return result;
  }
  result.objective = ls.objective;
  result.status = ls.status;
  // "no_decrease": no damped step lowers the objective any more, which is how LM ends at a minimum within rounding.
  result.converged = ls.status == "converged" || ls.status == "no_decrease";
  result.calls = ls.calls;
  result.freeIndices = ls.freeIndices;
  result.covariance = ls.covariance;
  result.leastSquaresBand = ls.band;
  result.x = params.GetMinuitParams().Params();
  return result;
}

FitResult FitLevenbergMarquardt(AZURECalc &objective, AZUREParams &params, const FitControl &control) {
  *control.log << "Using Levenberg-Marquardt minimizer with analytic Jacobian." << std::endl;
  return FitByLeastSquares(objective, params, control, "lm", &RunLevenbergMarquardt);
}
