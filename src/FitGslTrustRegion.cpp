// The GSL trust-region fitter (gsl-lm): geodesic-accelerated Levenberg-Marquardt on the residual rows, with the
// analytic Jacobian (LeastSquares.h; docs/upgrade/fit.md).
#include "LeastSquares.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <string>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_vector.h>

#include "AZURECalc.h"
#include "WorkingCopy.h"
#include "CNuc.h"
#include "AZUREParams.h"
#include "Fitters.h"

namespace {

// State shared with the GSL least-squares callbacks; owned by RunGSLNonlinear.
struct GSLNlsData {
  const AZURECalc *self;
  const FitControl *control;
  const char *fitter;
  long calls = 0;  // residual passes
  int nFree;
  int nData;                    // number of data residual rows
  const std::vector<int> *p2f;  // packed free param -> full-vector index
  const std::vector<double> *lo;
  const std::vector<double> *hi;
  const std::vector<double> *pen_nom;
  const std::vector<double> *pen_inv2;
  const std::vector<int> *penRows;  // free-param indices that carry a prior
  vector_r baseFull;                // full vector holding the fixed params
};

// Project x into [lo,hi] and scatter it into the full parameter vector.
void gslExpandFull(const GSLNlsData *d, const gsl_vector *x, vector_r &full) {
  full = d->baseFull;
  for (int a = 0; a < d->nFree; a++) {
    double v = gsl_vector_get(x, a);
    v = std::min(std::max(v, (*d->lo)[a]), (*d->hi)[a]);
    full[(*d->p2f)[a]] = v;
  }
}

// Residual vector: data rows plus Gaussian-prior rows sqrt(1/sigma^2)*(x-nom).
int gslResidualF(const gsl_vector *x, void *params, gsl_vector *f) {
  GSLNlsData *d = static_cast<GSLNlsData *>(params);
  // A stopped fit: the step fails, and the driver returns with the last accepted point.
  if (d->control->stop && d->control->stop->load(std::memory_order_relaxed)) return GSL_EBADFUNC;
  d->calls++;
  vector_r full;
  gslExpandFull(d, x, full);
  vector_r res;
  if (!d->self->ResidualsOnly(full, res) || (int)res.size() != d->nData) return GSL_EBADFUNC;
  for (int i = 0; i < d->nData; i++) gsl_vector_set(f, i, res[i]);
  for (size_t k = 0; k < d->penRows->size(); k++) {
    int a = (*d->penRows)[k];
    double s = std::sqrt((*d->pen_inv2)[a]);
    gsl_vector_set(f, d->nData + (int)k, s * (full[(*d->p2f)[a]] - (*d->pen_nom)[a]));
  }
  return GSL_SUCCESS;
}

// Residual Jacobian: analytic data-row block plus the constant prior rows.
int gslResidualDf(const gsl_vector *x, void *params, gsl_matrix *Jm) {
  GSLNlsData *d = static_cast<GSLNlsData *>(params);
  if (d->control->stop && d->control->stop->load(std::memory_order_relaxed)) return GSL_EBADFUNC;
  d->calls++;
  vector_r full;
  gslExpandFull(d, x, full);
  vector_r r, Jflat;
  std::vector<int> p2f;
  if (!d->self->ResidualJacobian(full, r, Jflat, p2f)) return GSL_EBADFUNC;
  if ((int)p2f.size() != d->nFree || (int)r.size() != d->nData) return GSL_EBADFUNC;
  const int nFree = d->nFree;
  for (int i = 0; i < d->nData; i++)
    for (int a = 0; a < nFree; a++)
      gsl_matrix_set(Jm, i, a, Jflat[(size_t)i * nFree + a]);
  for (size_t k = 0; k < d->penRows->size(); k++) {
    int row = d->nData + (int)k;
    for (int a = 0; a < nFree; a++) gsl_matrix_set(Jm, row, a, 0.0);
    int a = (*d->penRows)[k];
    gsl_matrix_set(Jm, row, a, std::sqrt((*d->pen_inv2)[a]));
  }
  return GSL_SUCCESS;
}

// Per-iteration progress line (chi^2 = ||residual||^2, priors included).
void gslProgress(const size_t iter, void *params,
                 const gsl_multifit_nlinear_workspace *w) {
  const GSLNlsData *d = static_cast<const GSLNlsData *>(params);
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  double chi2 = 0.0;
  gsl_blas_ddot(f, f, &chi2);
  vector_r full = d->baseFull;
  const gsl_vector *x = gsl_multifit_nlinear_position(w);
  for (int a = 0; a < d->nFree; a++) full[(*d->p2f)[a]] = gsl_vector_get(x, a);
  if (d->control->progress) {
    FitProgress step;
    step.kind = FitProgress::Kind::kIteration;
    step.fitter = d->fitter;
    step.label = "GSL-LM iteration";
    step.iteration = (long)iter;
    step.calls = d->calls;
    step.objective = chi2;
    step.best = chi2;
    step.x = &full;
    d->control->progress(step);
  }

  // Intermediate output, so a long GSL fit also leaves usable files behind
  // before it converges rather than only at the end.
  if (iter > 0 && iter % AZURECalc::kIterationOutputInterval == 0 && d->control->checkpoint) d->control->checkpoint(full);
}

}  // namespace

LeastSquaresOutcome RunGSLNonlinear(const AZURECalc &self, AZUREParams &params, int maxIter, const char *fitter,
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

  vector_r r0, J0, x;
  std::vector<int> p2f;
  std::vector<double> lo, hi, pen_nom, pen_inv2;
  const int nFree = PrepareFreeParams(self, full, mp, r0, J0, p2f, x, lo, hi, pen_nom, pen_inv2);
  if (nFree < 0) return out;  // analytic Jacobian unsupported -> the caller falls back
  if (nFree == 0) {
    out.objective = self.Chi2Value(full);
    return out;
  }
  const int nData = (int)r0.size();

  // Gaussian priors on norms / shifts / nuisance params become extra residual rows.
  std::vector<int> penRows;
  for (int a = 0; a < nFree; a++)
    if (pen_inv2[a] > 0.0) penRows.push_back(a);
  const int nRes = nData + (int)penRows.size();

  GSLNlsData d;
  d.self = &self;
  d.control = &control;
  d.fitter = fitter;
  d.calls = 1;  // PrepareFreeParams's pass
  d.nFree = nFree;
  d.nData = nData;
  d.p2f = &p2f;
  d.lo = &lo;
  d.hi = &hi;
  d.pen_nom = &pen_nom;
  d.pen_inv2 = &pen_inv2;
  d.penRows = &penRows;
  d.baseFull = full;

  gsl_multifit_nlinear_fdf fdf;
  fdf.f = gslResidualF;
  fdf.df = gslResidualDf;
  fdf.fvv = nullptr;  // second directional derivative via finite differences of f
  fdf.n = nRes;
  fdf.p = nFree;
  fdf.params = &d;

  // Geodesic-accelerated LM in a trust region; the default solver factorizes J
  // directly (QR) rather than forming J^T J.
  gsl_multifit_nlinear_parameters gp = gsl_multifit_nlinear_default_parameters();
  gp.trs = gsl_multifit_nlinear_trs_lmaccel;

  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *w = gsl_multifit_nlinear_alloc(T, &gp, nRes, nFree);
  if (!w) return out;

  gsl_vector *x0 = gsl_vector_alloc(nFree);
  for (int a = 0; a < nFree; a++) gsl_vector_set(x0, a, x[a]);

  gsl_set_error_handler_off();
  gsl_multifit_nlinear_init(x0, &fdf, w);

  const double xtol = 1.0e-9, gtol = 1.0e-9, ftol = 1.0e-9;
  int info = 0;
  int status = gsl_multifit_nlinear_driver(maxIter, xtol, gtol, ftol,
                                           gslProgress, &d, &info, w);
  log << std::endl;

  std::string reason = "reached the maximum number of iterations";
  out.status = "iteration_limit";
  if (status == GSL_SUCCESS) {
    reason = (info == 1) ? "converged (small step)" : "converged (small gradient)";
    out.status = "converged";
  } else if (control.stop && control.stop->load(std::memory_order_relaxed)) {
    reason = "stopped: the caller's stop request";
    out.status = "stopped";
  } else if (status != GSL_EMAXITER) {
    reason = std::string("stopped: ") + gsl_strerror(status);
    out.status = std::string("failed_") + gsl_strerror(status);
  }
  log << "GSL trust-region solver " << reason << " after "
                        << gsl_multifit_nlinear_niter(w) << " iterations." << std::endl;

  // Read back the (projected) best-fit parameters.
  gsl_vector *xf = gsl_multifit_nlinear_position(w);
  for (int a = 0; a < nFree; a++) {
    double v = gsl_vector_get(xf, a);
    v = std::min(std::max(v, lo[a]), hi[a]);
    full[p2f[a]] = v;
  }
  for (int i = 0; i < nFull; i++) mp.SetValue(i, full[i]);

  gsl_vector_free(x0);
  gsl_multifit_nlinear_free(w);

  double cost = self.Chi2Value(full);
  out.calls = d.calls + 1;

  // Errors and band covariance from (J^T J + priors)^{-1}, as the LM path does.
  if (out.status != "stopped") FinalizeLeastSquaresCovariance(self, full, p2f, pen_inv2, mp, out, log);
  out.objective = cost;
  return out;
}

FitResult FitGslTrustRegion(AZURECalc &objective, AZUREParams &params, const FitControl &control) {
  // GSL trust-region least-squares (geodesic-accelerated); same analytic
  // Jacobian as the LM path, but factorizes J directly (QR).
  *control.log << "Using GSL trust-region least-squares minimizer with analytic Jacobian." << std::endl;
  return FitByLeastSquares(objective, params, control, "gsl-lm", &RunGSLNonlinear);
}
