// inference_prior_covariance: the multivariate Gaussian prior (azure2/inference/PriorCovariance.h) on covariances
// small enough to invert by hand.  chi2 = (p - mu)^T C^-1 (p - mu) against the closed form of a 2 x 2 matrix, the
// gradient against central differences, the whitened rows and their Jacobian L^-1, a fixed parameter marginalised
// (the remaining block's Gaussian), the logarithm of the normalisation, and each refusal of the grammar by its
// message.  An inverse taken of the full matrix before the fixed row is dropped (conditioning instead of
// marginalising) fails the marginal case.
#include "azure2/inference/PriorCovariance.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include "check.h"

using azure2::inference::PriorCovariance;

namespace {

bool Parse(const std::string &text, PriorCovariance::File &file, std::string &err) {
  std::istringstream in(text);
  return PriorCovariance::Parse(in, "prior.cov", file, err);
}

bool Refused(const std::string &text, const std::string &expected) {
  PriorCovariance::File file;
  std::string err;
  if (Parse(text, file, err)) return false;
  if (err.find(expected) == std::string::npos) std::printf("    message: %s\n", err.c_str());
  return err.find(expected) != std::string::npos;
}

bool Near(double a, double b, double rel) { return std::fabs(a - b) <= rel * std::max(1.0, std::fabs(b)); }

}  // namespace

int main() {
  const double pi = 3.14159265358979323846;

  // ---- 2 x 2: C = [[4, 1.2], [1.2, 1]] (sigma 2 and 1, correlation 0.6), det = 2.56,
  //      C^-1 = [[1, -1.2], [-1.2, 4]] / 2.56.
  PriorCovariance::File file;
  std::string err;
  const std::string text2 = "#parametrization 1.000000e+00 0.000000e+00\n"
                            "# a comment\n"
                            "energy_1 width_1_1   # names\n"
                            "2.37 0.5\n"
                            "4.0 1.2\n"
                            "1.2 1.0\n";
  check("the 2 x 2 file parses", Parse(text2, file, err) && file.path == "prior.cov" && file.basis == 1 &&
                                     file.names.size() == 2 && file.names[1] == "width_1_1" && file.mean[0] == 2.37 &&
                                     file.cov[1] == 1.2 && file.cov[2] == 1.2);
  // The run's vector: [energy_1, width_1_1, segment_1_norm]; nothing fixed.
  const std::vector<std::string> names = {"energy_1", "width_1_1", "segment_1_norm"};
  PriorCovariance prior;
  std::vector<std::string> notes;
  check("set up on the run's vector", prior.Setup(file, names, {false, false, false}, 1, err, notes) && prior.Size() == 2 &&
                                          prior.FullIndex()[0] == 0 && prior.FullIndex()[1] == 1 && notes.empty() &&
                                          prior.Path() == "prior.cov" && prior.Marginalized().empty());
  const std::vector<double> p = {2.37 + 1.0, 0.5 - 0.5, 7.0};  // d = (1, -0.5); the norm is not in the prior
  // d^T C^-1 d = (1*1 + 2*(-1.2)*1*(-0.5) + 4*0.25) / 2.56 = 3.2 / 2.56 = 1.25
  check("chi2 = 1.25 by hand (1e-12)", Near(prior.Chi2(p), 1.25, 1e-12));
  check("zero at the mean", prior.Chi2({2.37, 0.5, -3.0}) == 0.0);
  check("C^-1 by hand", Near(prior.Cinv(0, 0), 1.0 / 2.56, 1e-12) && Near(prior.Cinv(0, 1), -1.2 / 2.56, 1e-12) &&
                            Near(prior.Cinv(1, 0), -1.2 / 2.56, 1e-12) && Near(prior.Cinv(1, 1), 4.0 / 2.56, 1e-12));
  // L = [[2, 0], [0.6, 0.8]]: L^-1 = [[0.5, 0], [-0.375, 1.25]]; r = L^-1 d = (0.5, -1.0); r.r = 1.25.
  std::vector<double> r;
  prior.Residuals(p, r);
  check("whitened residuals r = L^-1 (p - mu) by hand", r.size() == 2 && Near(r[0], 0.5, 1e-12) && Near(r[1], -1.0, 1e-12));
  check("L^-1 by hand", Near(prior.Linv(0, 0), 0.5, 1e-12) && prior.Linv(0, 1) == 0.0 && Near(prior.Linv(1, 0), -0.375, 1e-12) &&
                            Near(prior.Linv(1, 1), 1.25, 1e-12));
  std::vector<double> grad(3, 0.0);
  prior.AddGradient(p, grad);
  bool gradient = grad[2] == 0.0;
  for (int i = 0; i < 2; i++) {
    std::vector<double> up = p, down = p;
    const double h = 1e-5;
    up[i] += h;
    down[i] -= h;
    gradient = gradient && Near(grad[i], (prior.Chi2(up) - prior.Chi2(down)) / (2.0 * h), 1e-9);
  }
  // 2 C^-1 d = 2 (1 + 0.6, -1.2 - 2) / 2.56 = (1.25, -2.5)
  check("gradient = central differences (1e-9), nothing on the norm", gradient);
  check("gradient by hand: (1.25, -2.5)", Near(grad[0], 1.25, 1e-12) && Near(grad[1], -2.5, 1e-12));
  std::vector<double> more(3, 10.0);
  prior.AddGradient(p, more);
  check("AddGradient adds", Near(more[0], 11.25, 1e-12) && more[2] == 10.0);
  std::vector<double> g2(3, 0.0);
  prior.AddCurvature(g2);
  check("curvature 2 (C^-1)_ii", Near(g2[0], 2.0 / 2.56, 1e-12) && Near(g2[1], 8.0 / 2.56, 1e-12) && g2[2] == 0.0);
  // The rows: residuals appended, the Jacobian over packed columns (here full index f -> column 2 - f, the norm none).
  std::vector<double> res = {42.0}, jac;
  prior.AppendRows(p, [](int f) { return f < 2 ? 1 - f : -1; }, 2, res, &jac);
  check("AppendRows: residuals after what was there", res.size() == 3 && res[0] == 42.0 && Near(res[1], 0.5, 1e-12) && Near(res[2], -1.0, 1e-12));
  check("AppendRows: Jacobian rows = L^-1 in the caller's columns",
        jac.size() == 4 && jac[0] == 0.0 && Near(jac[1], 0.5, 1e-12) && Near(jac[2], 1.25, 1e-12) && Near(jac[3], -0.375, 1e-12));
  check("ln of the normalisation: -0.5 ln det(2 pi C)", Near(prior.LogNormalization(), -0.5 * (2.0 * std::log(2.0 * pi) + std::log(2.56)), 1e-12));
  // A vector too short for an index: that parameter sits at its mean.
  check("a vector without the second parameter: d = (1, 0)", Near(prior.Chi2({3.37}), 1.0 / 2.56, 1e-12));

  // ---- A fixed parameter is marginalised: the prior on the other is N(mu_2, C_22), not the conditional.
  PriorCovariance marginal;
  notes.clear();
  check("energy_1 fixed: set up with a note", marginal.Setup(file, names, {true, false, false}, 1, err, notes) && marginal.Size() == 1 &&
                                                  marginal.FullIndex()[0] == 1 && marginal.Marginalized().size() == 1 &&
                                                  marginal.Marginalized()[0] == "energy_1" && notes.size() == 1 &&
                                                  notes[0].find("marginalized") != std::string::npos);
  // d_2 = -0.5, C_22 = 1: chi2 = 0.25 (the conditional would give 4 * 0.25 / 2.56 = 0.390625).
  check("the marginal: chi2 = d^2 / C_22 = 0.25", Near(marginal.Chi2(p), 0.25, 1e-12));
  check("the marginal's normalisation", Near(marginal.LogNormalization(), -0.5 * std::log(2.0 * pi * 1.0), 1e-12));

  // ---- 3 x 3, diagonal plus one correlation, against the sum of an independent Gaussian and the 2 x 2 block.
  const std::string text3 = "a b c\n1 2 3\n0.25 0 0\n0 4.0 1.2\n0 1.2 1.0\n";
  PriorCovariance three;
  notes.clear();
  check("the 3 x 3 file", Parse(text3, file, err) && file.basis == -1 && three.Setup(file, {"c", "b", "a"}, {}, 2, err, notes) &&
                              three.Size() == 3 && three.FullIndex()[0] == 2 && three.FullIndex()[2] == 0);
  // The run's order is (c, b, a): a = 1 + 1, b = 2 + 1, c = 3 - 0.5 -> (1/0.5)^2 + 1.25 = 5.25.
  check("chi2 = 4 + 1.25 (1e-12)", Near(three.Chi2({2.5, 3.0, 2.0}), 5.25, 1e-12));
  check("normalisation of the 3 x 3", Near(three.LogNormalization(), -0.5 * (3.0 * std::log(2.0 * pi) + std::log(0.25 * 2.56)), 1e-12));

  // ---- Refusals of the set-up.
  PriorCovariance refused;
  Parse(text2, file, err);
  check("another basis than the run's is refused", !refused.Setup(file, names, {}, 2, err, notes) &&
                                                       err.find("parametrization 1 (Brune)") != std::string::npos &&
                                                       err.find("2 (Park)") != std::string::npos);
  check("a name that is no parameter of the run", !refused.Setup(file, {"energy_1", "segment_1_norm"}, {}, 1, err, notes) &&
                                                      err.find(file.names[1]) != std::string::npos &&
                                                      err.find("not a parameter of this model") != std::string::npos);

  // ---- Refusals of the grammar, each by its message.
  check("too few lines", Refused("a b\n", "expected a line of parameter names"));
  check("a name twice", Refused("a a\n1 2\n1 0\n0 1\n", "listed twice"));
  check("names and means differ in number", Refused("a b\n1\n1 0\n0 1\n", "2 names but 1 means"));
  check("a mean that is not a number", Refused("a b\n1 x\n1 0\n0 1\n", "'x' is not a number"));
  check("nan is not a number here", Refused("a b\n1 nan\n1 0\n0 1\n", "'nan' is not a number"));
  check("too few covariance rows", Refused("a b\n1 2\n1 0\n", "2 parameters but 1 covariance rows"));
  check("a row of the wrong length", Refused("a b\n1 2\n1 0 0\n0 1\n", "covariance row 1 has 3 entries, expected 2"));
  check("an entry that is not a number", Refused("a b\n1 2\n1 0\n0 y\n", "'y' is not a number"));
  check("not symmetric", Refused("a b\n1 2\n1 0.5\n0.4 1\n", "not symmetric"));
  check("not positive definite", Refused("a b\n1 2\n1 2\n2 1\n", "not positive definite"));
  check("a bad basis tag", Refused("#parametrization 3\na b\n1 2\n1 0\n0 1\n", "must be 0, 1 or 2"));
  check("symmetric within 1e-10 is accepted and made exact",
        Parse("a b\n1 2\n1 0.50000000001\n0.5 1\n", file, err) && file.cov[1] == file.cov[2]);
  return result();
}
