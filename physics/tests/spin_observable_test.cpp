// The operator token of a spin observable, its matrices and the trace formula (SpinObservable.h), checked against
// things that do not use the code under test:
//
//  1. the grammar: tokens that parse and tokens that do not, each with the reason; the spin rules (Check);
//  2. the matrices: spin 1/2 vector operators are the Pauli matrices; for spin 1 the Cartesian tensors are traceless
//     in sum (xx + yy + zz = 0) and the spherical tensors are the Madison convention's closed forms in S_x, S_y, S_z
//     (tau_10 = sqrt(3/2) S_z, tau_11 = -(sqrt 3 / 2)(S_x + i S_y), tau_20 = (3 S_z^2 - 2)/sqrt 2,
//     tau_21 = -(sqrt 3 / 2)((S_x + i S_y) S_z + S_z (S_x + i S_y)), tau_22 = (sqrt 3 / 2)(S_x + i S_y)^2), written
//     here element by element; Tr(tau_kq tau_kq^dagger) = 2j + 1 for j = 1/2, 1, 3/2;
//  3. the value for a spin-1/2 particle on a spin-0 target, F = g + h sigma_y: A_y = P_y = 2 Re(g h*)/(|g|^2 + |h|^2),
//     the transfer b:y>e1:y = 1 and b:x>e1:x = b:z>e1:z = (|g|^2 - |h|^2)/(|g|^2 + |h|^2);
//  4. the gradient dO/dF* against central differences of the value in the real and imaginary part of every element
//     of a fixed non-symmetric F, for d + p -> p + d-like spins (1, 1/2 -> 1/2, 1) and eleven tokens (1e-7).
// A wrong normalization, a transposed operator, a wrong sign of a spherical tensor or a missing adjoint term fails at
// least one of these by more than 1e-3.
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

#include "azure2/physics/SpinObservable.h"

using namespace azure2::physics;
typedef std::complex<double> cx;

namespace {

int failures = 0;

void Check(const std::string &name, bool ok, const std::string &detail = "") {
  std::printf("  %s  %s%s%s\n", ok ? "ok  " : "FAIL", name.c_str(), ok || detail.empty() ? "" : "  -- ", ok ? "" : detail.c_str());
  if (!ok) failures++;
}

double MaxDiff(const std::vector<cx> &a, const std::vector<cx> &b) {
  if (a.size() != b.size()) return 1.e300;
  double d = 0.0;
  for (std::size_t i = 0; i < a.size(); i++) d = std::max(d, std::abs(a[i] - b[i]));
  return d;
}

std::vector<cx> Matrix(const std::string &op, double j) {
  SpinObservable o;
  std::string error;
  if (!o.Parse("b:" + op, error)) return {};
  return SpinObservable::TermMatrix(o.terms[0], j);
}

std::vector<cx> Mul(const std::vector<cx> &a, const std::vector<cx> &b, int n) {
  std::vector<cx> c(n * n, cx(0.0, 0.0));
  for (int i = 0; i < n; i++)
    for (int k = 0; k < n; k++)
      for (int j = 0; j < n; j++) c[i * n + j] += a[i * n + k] * b[k * n + j];
  return c;
}

double Value(const std::string &token, double j1, double j2, double j3, double j4, const std::vector<cx> &F,
             std::vector<cx> *gradient = nullptr) {
  SpinObservable o;
  std::string error;
  if (!o.Parse(token, error) || !o.Check(j1, j2, j3, j4, error)) return std::nan("");
  return SpinObservableValue(o, j1, j2, j3, j4, F, gradient);
}

}  // namespace

int main() {
  std::printf("1. the grammar\n");
  for (const char *t : {"b:y", "t:y", "b:y*t:y", "b:yy", "b:xz", "b:iT11", "b:T20", "e1:y", "e2:T22", "b:y>e1:y",
                        "b:x*t:z>e1:x*e2:y"}) {
    SpinObservable o;
    std::string error;
    Check(std::string("parses: ") + t, o.Parse(t, error), error);
  }
  const char *bad[][2] = {{"", "empty side"},
                          {"y", "has no ':'"},
                          {"p:y", "unknown particle"},
                          {"b:w", "unknown operator"},
                          {"b:y*b:x", "appears twice"},
                          {"b:y*e1:y", "entrance and exit particles on one side"},
                          {"e1:y>b:y", "on the entrance side"},
                          {"b:y>t:y", "on the exit side"},
                          {"b:y>e1:y>e2:y", "more than one '>'"},
                          {"b:T00", "needs 1 <= k"},
                          {"b:T12", "needs 1 <= k"},
                          {"b:iT11*t:iT11", "at most one iT_kq"}};
  for (const auto &b : bad) {
    SpinObservable o;
    std::string error;
    const bool parsed = o.Parse(b[0], error);
    Check(std::string("refused: '") + b[0] + "' (" + b[1] + ")", !parsed && error.find(b[1]) != std::string::npos, error);
  }
  {
    SpinObservable o;
    std::string error;
    o.Parse("t:y", error);
    Check("spin 0 has no observable", !o.Check(0.5, 0.0, 0.5, 0.0, error) && error.find("spin 0") != std::string::npos, error);
    o.Parse("b:yy", error);
    Check("a Cartesian tensor needs spin 1", !o.Check(0.5, 0.0, 0.5, 0.0, error) && error.find("spin 1 only") != std::string::npos, error);
    Check("a Cartesian tensor on spin 1", o.Check(1.0, 0.0, 1.0, 0.0, error), error);
    o.Parse("b:T20", error);
    Check("rank 2 exceeds 2j for spin 1/2", !o.Check(0.5, 0.0, 0.5, 0.0, error) && error.find("exceeds 2j") != std::string::npos, error);
  }

  std::printf("\n2. the matrices\n");
  const cx I(0.0, 1.0);
  Check("spin 1/2: x, y, z are the Pauli matrices",
        MaxDiff(Matrix("x", 0.5), {0.0, 1.0, 1.0, 0.0}) < 1e-15 && MaxDiff(Matrix("y", 0.5), {0.0, -I, I, 0.0}) < 1e-15 &&
            MaxDiff(Matrix("z", 0.5), {1.0, 0.0, 0.0, -1.0}) < 1e-15);
  {
    const double r = std::sqrt(0.5);
    // S_x, S_y, S_z of spin 1 in the basis m = 1, 0, -1, written out
    const std::vector<cx> sx = {0, r, 0, r, 0, r, 0, r, 0}, sy = {0, -I * r, 0, I * r, 0, -I * r, 0, I * r, 0},
                          sz = {1, 0, 0, 0, 0, 0, 0, 0, -1};
    Check("spin 1: x, y, z are S_x, S_y, S_z", MaxDiff(Matrix("x", 1.0), sx) < 1e-15 && MaxDiff(Matrix("y", 1.0), sy) < 1e-15 &&
                                                   MaxDiff(Matrix("z", 1.0), sz) < 1e-15);
    std::vector<cx> sum = Matrix("xx", 1.0);
    const std::vector<cx> yy = Matrix("yy", 1.0), zz = Matrix("zz", 1.0);
    for (std::size_t i = 0; i < sum.size(); i++) sum[i] += yy[i] + zz[i];
    Check("spin 1: xx + yy + zz = 0", MaxDiff(sum, std::vector<cx>(9, cx(0.0, 0.0))) < 1e-14);
    std::vector<cx> sp(9);  // S_x + i S_y
    for (int i = 0; i < 9; i++) sp[i] = sx[i] + I * sy[i];
    const double h = std::sqrt(3.0) / 2.0;
    std::vector<cx> t10(9), t11(9), t20 = Mul(sz, sz, 3), t21(9), t22 = Mul(sp, sp, 3);
    const std::vector<cx> a = Mul(sp, sz, 3), b = Mul(sz, sp, 3);
    for (int i = 0; i < 9; i++) {
      t10[i] = std::sqrt(1.5) * sz[i];
      t11[i] = -h * sp[i];
      t20[i] = (3.0 * t20[i] - ((i % 4 == 0) ? 2.0 : 0.0)) / std::sqrt(2.0);
      t21[i] = -h * (a[i] + b[i]);
      t22[i] = h * t22[i];
    }
    Check("spin 1: tau_10 = sqrt(3/2) S_z", MaxDiff(Matrix("T10", 1.0), t10) < 1e-14);
    Check("spin 1: tau_11 = -(sqrt 3 / 2)(S_x + i S_y)", MaxDiff(Matrix("T11", 1.0), t11) < 1e-14);
    Check("spin 1: tau_20 = (3 S_z^2 - 2)/sqrt 2", MaxDiff(Matrix("T20", 1.0), t20) < 1e-14);
    Check("spin 1: tau_21 = -(sqrt 3 / 2)((S_x + i S_y) S_z + S_z (S_x + i S_y))", MaxDiff(Matrix("T21", 1.0), t21) < 1e-14);
    Check("spin 1: tau_22 = (sqrt 3 / 2)(S_x + i S_y)^2", MaxDiff(Matrix("T22", 1.0), t22) < 1e-14);
    std::vector<cx> zzFromT20 = Matrix("T20", 1.0);
    for (cx &x : zzFromT20) x *= std::sqrt(2.0);
    Check("spin 1: zz = sqrt 2 tau_20", MaxDiff(zz, zzFromT20) < 1e-14);
  }
  for (double j : {0.5, 1.0, 1.5})
    for (int k = 1; k <= (int)std::lround(2.0 * j); k++)
      for (int q = 0; q <= k; q++) {
        const std::string op = "T" + std::to_string(k) + std::to_string(q);
        const std::vector<cx> m = Matrix(op, j);
        double trace = 0.0;
        for (const cx &x : m) trace += std::norm(x);
        Check("Tr(tau tau^dagger) = 2j + 1: j = " + std::to_string(j).substr(0, 3) + ", " + op, std::fabs(trace - (2.0 * j + 1.0)) < 1e-13,
              std::to_string(trace));
      }

  std::printf("\n3. spin 1/2 on spin 0: F = g + h sigma_y\n");
  {
    const cx g(0.7, -0.4), h(-0.3, 0.9);
    // sigma_y = [[0, -i], [i, 0]]; F row-major (exit index, entrance index)
    const std::vector<cx> F = {g, -I * h, I * h, g};
    const double D = std::norm(g) + std::norm(h), ay = 2.0 * std::real(g * std::conj(h)) / D, k = (std::norm(g) - std::norm(h)) / D;
    Check("b:y = 2 Re(g h*)/(|g|^2 + |h|^2)", std::fabs(Value("b:y", 0.5, 0, 0.5, 0, F) - ay) < 1e-14);
    Check("e1:y = the same (P_y = A_y)", std::fabs(Value("e1:y", 0.5, 0, 0.5, 0, F) - ay) < 1e-14);
    Check("b:x = b:z = 0", std::fabs(Value("b:x", 0.5, 0, 0.5, 0, F)) < 1e-14 && std::fabs(Value("b:z", 0.5, 0, 0.5, 0, F)) < 1e-14);
    Check("b:y>e1:y = 1", std::fabs(Value("b:y>e1:y", 0.5, 0, 0.5, 0, F) - 1.0) < 1e-14);
    Check("b:x>e1:x = b:z>e1:z = (|g|^2 - |h|^2)/(|g|^2 + |h|^2)",
          std::fabs(Value("b:x>e1:x", 0.5, 0, 0.5, 0, F) - k) < 1e-14 && std::fabs(Value("b:z>e1:z", 0.5, 0, 0.5, 0, F) - k) < 1e-14);
    Check("F = 0 gives 0", Value("b:y", 0.5, 0, 0.5, 0, std::vector<cx>(4, cx(0.0, 0.0))) == 0.0);
  }

  std::printf("\n4. the gradient against central differences (spins 1, 1/2 -> 1/2, 1)\n");
  {
    const int din = 6, dout = 6;
    std::vector<cx> F(din * dout);
    for (int i = 0; i < din * dout; i++) F[i] = cx(std::sin(1.0 + 0.37 * i), std::cos(0.5 + 0.91 * i) * (1.0 + 0.01 * i));
    for (const char *t : {"b:y", "t:y", "b:y*t:y", "b:yy", "b:xz", "b:iT11", "b:T20", "b:T22", "e2:iT21", "b:y>e1:y",
                          "b:xz*t:x>e1:z*e2:T21"}) {
      std::vector<cx> grad;
      const double value = Value(t, 1.0, 0.5, 0.5, 1.0, F, &grad);
      double worst = 0.0;
      const double step = 1e-6;
      for (std::size_t x = 0; x < F.size(); x++) {
        std::vector<cx> up = F, down = F;
        up[x] += step, down[x] -= step;
        const double dRe = (Value(t, 1.0, 0.5, 0.5, 1.0, up) - Value(t, 1.0, 0.5, 0.5, 1.0, down)) / (2.0 * step);
        up = F, down = F;
        up[x] += cx(0.0, step), down[x] -= cx(0.0, step);
        const double dIm = (Value(t, 1.0, 0.5, 0.5, 1.0, up) - Value(t, 1.0, 0.5, 0.5, 1.0, down)) / (2.0 * step);
        // dO/dF* = (dO/dRe + i dO/dIm)/2
        worst = std::max(worst, std::abs(grad[x] - 0.5 * cx(dRe, dIm)));
      }
      Check(std::string(t) + ": value finite, gradient to 1e-7", std::isfinite(value) && worst < 1e-7, "worst " + std::to_string(worst));
    }
  }

  std::printf("\n%s\n", failures ? "FAILED" : "ALL PASSED");
  return failures ? 1 : 0;
}
