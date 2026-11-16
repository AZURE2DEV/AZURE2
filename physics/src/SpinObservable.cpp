// The operator token of a spin observable, its matrices and the trace formula (SpinObservable.h): angular-momentum
// algebra only.
#include "azure2/physics/SpinObservable.h"

#include <cctype>
#include <cmath>
#include <sstream>

#include "azure2/physics/AngularCoupling.h"

namespace azure2 {
namespace physics {

namespace {

typedef std::complex<double> complex;

int Dim(double j) { return (int)std::lround(2.0 * j + 1.0); }

/// Spin matrices S_x, S_y, S_z for spin j, basis m = j, ..., -j.
std::vector<complex> SpinMatrix(int axis, double j) {
  const int n = Dim(j);
  std::vector<complex> s(n * n, complex(0.0, 0.0));
  for (int r = 0; r < n; r++) {
    const double mr = j - r;
    for (int c = 0; c < n; c++) {
      const double mc = j - c;
      if (axis == 2) {
        if (r == c) s[r * n + c] = mr;
        continue;
      }
      // <m+1|S_+|m> = sqrt(j(j+1) - m(m+1)),  <m-1|S_-|m> = sqrt(j(j+1) - m(m-1))
      double plus = 0.0, minus = 0.0;
      if (std::fabs(mr - (mc + 1.0)) < 1.e-9) plus = std::sqrt(j * (j + 1.0) - mc * (mc + 1.0));
      if (std::fabs(mr - (mc - 1.0)) < 1.e-9) minus = std::sqrt(j * (j + 1.0) - mc * (mc - 1.0));
      if (axis == 0)
        s[r * n + c] = complex(0.5 * (plus + minus), 0.0);  // (S+ + S-)/2
      else
        s[r * n + c] = complex(0.0, -0.5 * (plus - minus));  // (S+ - S-)/(2i)
    }
  }
  return s;
}

std::vector<complex> MatMul(const std::vector<complex> &a, const std::vector<complex> &b, int n) {
  std::vector<complex> c(n * n, complex(0.0, 0.0));
  for (int i = 0; i < n; i++)
    for (int k = 0; k < n; k++) {
      const complex aik = a[i * n + k];
      if (aik == complex(0.0, 0.0)) continue;
      for (int j = 0; j < n; j++) c[i * n + j] += aik * b[k * n + j];
    }
  return c;
}

int Axis(char ch) { return ch == 'x' ? 0 : ch == 'y' ? 1 : ch == 'z' ? 2 : -1; }

/// B F A for F (dout x din), A (din x din), B (dout x dout), all row-major.
std::vector<complex> Sandwich(const std::vector<complex> &B, const std::vector<complex> &F,
                              const std::vector<complex> &A, int dout, int din) {
  std::vector<complex> FA(dout * din, complex(0.0, 0.0)), BFA(dout * din, complex(0.0, 0.0));
  for (int o = 0; o < dout; o++)
    for (int i = 0; i < din; i++) {
      const complex f = F[o * din + i];
      if (f == complex(0.0, 0.0)) continue;
      for (int k = 0; k < din; k++) FA[o * din + k] += f * A[i * din + k];
    }
  for (int o = 0; o < dout; o++)
    for (int p = 0; p < dout; p++) {
      const complex b = B[o * dout + p];
      if (b == complex(0.0, 0.0)) continue;
      for (int k = 0; k < din; k++) BFA[o * din + k] += b * FA[p * din + k];
    }
  return BFA;
}

std::vector<complex> Adjoint(const std::vector<complex> &m, int n) {
  std::vector<complex> a(n * n);
  for (int r = 0; r < n; r++)
    for (int c = 0; c < n; c++) a[r * n + c] = std::conj(m[c * n + r]);
  return a;
}

}  // namespace

bool SpinObservable::Parse(const std::string &tok, std::string &error) {
  token = tok;
  terms.clear();
  std::string left = tok, right;
  const std::size_t gt = tok.find('>');
  if (gt != std::string::npos) {
    left = tok.substr(0, gt);
    right = tok.substr(gt + 1);
    if (right.find('>') != std::string::npos) {
      error = "more than one '>'";
      return false;
    }
  }
  auto parseSide = [&](const std::string &side, int sideNum) -> bool {
    if (side.empty()) {
      error = "empty side";
      return false;
    }
    std::stringstream ss(side);
    std::string term;
    while (std::getline(ss, term, '*')) {
      const std::size_t colon = term.find(':');
      if (colon == std::string::npos) {
        error = "term \"" + term + "\" has no ':'";
        return false;
      }
      const std::string who = term.substr(0, colon), op = term.substr(colon + 1);
      SpinTerm t;
      if (who == "b")
        t.who = 1;
      else if (who == "t")
        t.who = 2;
      else if (who == "e1")
        t.who = 3;
      else if (who == "e2")
        t.who = 4;
      else {
        error = "unknown particle \"" + who + "\" (b, t, e1, e2)";
        return false;
      }
      // Without '>' a side may be all entrance or all exit; with '>' the left is entrance, the right exit.
      if (sideNum == 1 && gt != std::string::npos && t.who > 2) {
        error = "exit particle \"" + who + "\" on the entrance side of '>'";
        return false;
      }
      if (sideNum == 2 && t.who < 3) {
        error = "entrance particle \"" + who + "\" on the exit side of '>'";
        return false;
      }
      if (op.size() == 1 && Axis(op[0]) >= 0) {
        t.kind = 1;
        t.a = Axis(op[0]);
      } else if (op.size() == 2 && Axis(op[0]) >= 0 && Axis(op[1]) >= 0) {
        t.kind = 2;
        t.a = Axis(op[0]);
        t.b = Axis(op[1]);
      } else if ((op.size() == 3 && op[0] == 'T') || (op.size() == 4 && op[0] == 'i' && op[1] == 'T')) {
        const std::size_t p = (op[0] == 'i') ? 2 : 1;
        if (!std::isdigit((unsigned char)op[p]) || !std::isdigit((unsigned char)op[p + 1])) {
          error = "bad spherical operator \"" + op + "\"";
          return false;
        }
        t.kind = 3;
        t.k = op[p] - '0';
        t.q = op[p + 1] - '0';
        t.timesI = (op[0] == 'i');
        if (t.k < 1 || t.q > t.k) {
          error = "spherical operator \"" + op + "\" needs 1 <= k and 0 <= q <= k";
          return false;
        }
      } else {
        error = "unknown operator \"" + op + "\"";
        return false;
      }
      for (const SpinTerm &u : terms)
        if (u.who == t.who) {
          error = "particle \"" + who + "\" appears twice";
          return false;
        }
      terms.push_back(t);
    }
    return true;
  };
  if (!parseSide(left, 1)) return false;
  if (gt != std::string::npos && !parseSide(right, 2)) return false;
  bool anyIn = false, anyOut = false;
  for (const SpinTerm &t : terms) (t.who < 3 ? anyIn : anyOut) = true;
  if (gt == std::string::npos && anyIn && anyOut) {
    error = "entrance and exit particles on one side; write <entrance>'>'<exit>";
    return false;
  }
  int nI = 0;
  for (const SpinTerm &t : terms) nI += t.timesI ? 1 : 0;
  if (nI > 1) {
    error = "at most one iT_kq term";
    return false;
  }
  return true;
}

bool SpinObservable::Check(double j1, double j2, double j3, double j4, std::string &error) const {
  const char *const names[] = {"", "beam", "target", "exit particle 1", "exit particle 2"};
  const double spins[] = {0.0, j1, j2, j3, j4};
  for (const SpinTerm &t : terms) {
    const double j = spins[t.who];
    std::ostringstream js;
    js << j;
    if (j < 0.25) {
      error = std::string("the ") + names[t.who] + " has spin 0: no polarization observable";
      return false;
    }
    if (t.kind == 2 && std::fabs(j - 1.0) > 1.e-6) {
      error = std::string("Cartesian tensor operators are defined for spin 1 only (the ") + names[t.who] + " has spin " +
              js.str() + "); use T_kq";
      return false;
    }
    if (t.kind == 3 && t.k > (int)std::lround(2.0 * j)) {
      error = std::string("rank ") + std::to_string(t.k) + " exceeds 2j for the " + names[t.who] + " (spin " + js.str() +
              ")";
      return false;
    }
  }
  return true;
}

SpinObservable::complex SpinObservable::Prefactor() const {
  for (const SpinTerm &t : terms)
    if (t.timesI) return complex(0.0, 1.0);
  return complex(1.0, 0.0);
}

std::vector<SpinObservable::complex> SpinObservable::TermMatrix(const SpinTerm &t, double j) {
  const int n = Dim(j);
  std::vector<complex> m(n * n, complex(0.0, 0.0));
  if (t.kind == 1) {
    m = SpinMatrix(t.a, j);
    for (complex &x : m) x /= j;
  } else if (t.kind == 2) {
    const std::vector<complex> si = SpinMatrix(t.a, j), sj = SpinMatrix(t.b, j);
    const std::vector<complex> ij = MatMul(si, sj, n), ji = MatMul(sj, si, n);
    for (int r = 0; r < n * n; r++) m[r] = 1.5 * (ij[r] + ji[r]);
    if (t.a == t.b)
      for (int r = 0; r < n; r++) m[r * n + r] -= 2.0;
  } else if (t.kind == 3) {
    // <j m'| tau_kq |j m> = sqrt(2k+1) <j m k q | j m'>, so that Tr(tau tau^dagger) = 2j+1.
    const int twoJ = (int)std::lround(2.0 * j);
    for (int r = 0; r < n; r++)
      for (int c = 0; c < n; c++) {
        const int twoMp = twoJ - 2 * r, twoM = twoJ - 2 * c;
        if (twoMp - twoM != 2 * t.q) continue;
        m[r * n + c] = std::sqrt(2.0 * t.k + 1.0) * clebsch_gordan(twoJ, twoM, 2 * t.k, 2 * t.q, twoJ, twoMp);
      }
  }
  return m;
}

std::vector<SpinObservable::complex> SpinObservable::PairOperator(int first, int second, double jFirst,
                                                                    double jSecond) const {
  const int da = Dim(jFirst), db = Dim(jSecond);
  std::vector<complex> ma(da * da, complex(0.0, 0.0)), mb(db * db, complex(0.0, 0.0));
  for (int r = 0; r < da; r++) ma[r * da + r] = 1.0;
  for (int r = 0; r < db; r++) mb[r * db + r] = 1.0;
  for (const SpinTerm &t : terms) {
    if (t.who == first) ma = TermMatrix(t, jFirst);
    if (t.who == second) mb = TermMatrix(t, jSecond);
  }
  const int n = da * db;
  std::vector<complex> m(n * n, complex(0.0, 0.0));
  for (int a1 = 0; a1 < da; a1++)
    for (int b1 = 0; b1 < db; b1++)
      for (int a2 = 0; a2 < da; a2++)
        for (int b2 = 0; b2 < db; b2++) m[(a1 * db + b1) * n + (a2 * db + b2)] = ma[a1 * da + a2] * mb[b1 * db + b2];
  return m;
}

double SpinObservableValue(const SpinObservable &observable, double j1, double j2, double j3, double j4,
                           const std::vector<complex> &F, std::vector<complex> *gradient) {
  const int din = Dim(j1) * Dim(j2), dout = Dim(j3) * Dim(j4);
  if (gradient) gradient->assign(F.size(), complex(0.0, 0.0));
  if ((int)F.size() != din * dout) return 0.0;
  const std::vector<complex> A = observable.PairOperator(1, 2, j1, j2), B = observable.PairOperator(3, 4, j3, j4);
  // N = Tr(F A F^dagger B) = sum_{o,i} conj(F_oi) (B F A)_oi,  D = sum |F_oi|^2.
  const std::vector<complex> BFA = Sandwich(B, F, A, dout, din);
  complex N(0.0, 0.0);
  double D = 0.0;
  for (std::size_t x = 0; x < F.size(); x++) {
    N += std::conj(F[x]) * BFA[x];
    D += std::norm(F[x]);
  }
  if (D <= 0.0) return 0.0;
  const complex c = observable.Prefactor();
  const double O = std::real(c * N) / D;
  if (gradient) {
    // dN/dF* = B F A,  d(N*)/dF* = B^dagger F A^dagger,  dD/dF* = F.
    const std::vector<complex> BhFAh = Sandwich(Adjoint(B, dout), F, Adjoint(A, din), dout, din);
    for (std::size_t x = 0; x < F.size(); x++)
      (*gradient)[x] = (0.5 * (c * BFA[x] + std::conj(c) * BhFAh[x]) - O * F[x]) / D;
  }
  return O;
}

}  // namespace physics
}  // namespace azure2
