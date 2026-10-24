#include "SpinObservable.h"

#include <cctype>
#include <cmath>
#include <sstream>

#include "AngCoeff.h"

namespace {

int Dim(double j) { return (int)std::lround(2.0 * j + 1.0); }

//! Spin matrices S_x, S_y, S_z for spin j, basis m = j, ..., -j.
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
      if (axis == 0) s[r * n + c] = complex(0.5 * (plus + minus), 0.0);   // (S+ + S-)/2
      else s[r * n + c] = complex(0.0, -0.5 * (plus - minus));            // (S+ - S-)/(2i)
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
      if (who == "b") t.who = 1;
      else if (who == "t") t.who = 2;
      else if (who == "e1") t.who = 3;
      else if (who == "e2") t.who = 4;
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

double SpinObservable::SpinOf(const SpinTerm &t, double j1, double j2, double j3, double j4) {
  return t.who == 1 ? j1 : t.who == 2 ? j2 : t.who == 3 ? j3 : j4;
}

bool SpinObservable::Check(double j1, double j2, double j3, double j4, bool, std::string &error) const {
  static const char *names[] = {"", "beam", "target", "exit particle 1", "exit particle 2"};
  for (const SpinTerm &t : terms) {
    const double j = SpinOf(t, j1, j2, j3, j4);
    std::ostringstream js;
    js << j;
    if (j < 0.25) {
      error = std::string(names[t.who]) + " has spin 0: no polarization observable";
      return false;
    }
    if (t.kind == 2 && std::fabs(j - 1.0) > 1.e-6) {
      error = std::string("Cartesian tensor operators are defined for spin 1 only (") + names[t.who] +
              " has spin " + js.str() + "); use T_kq";
      return false;
    }
    if (t.kind == 3 && t.k > (int)std::lround(2.0 * j)) {
      error = std::string("rank ") + std::to_string(t.k) + " exceeds 2j for the " + names[t.who] +
              " (spin " + js.str() + ")";
      return false;
    }
  }
  return true;
}

complex SpinObservable::Prefactor() const {
  for (const SpinTerm &t : terms)
    if (t.timesI) return complex(0.0, 1.0);
  return complex(1.0, 0.0);
}

std::vector<complex> SpinObservable::TermMatrix(const SpinTerm &t, double j) {
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
    for (int r = 0; r < n; r++)
      for (int c = 0; c < n; c++) {
        const double mp = j - r, mm = j - c;
        if (std::fabs(mp - mm - t.q) > 1.e-9) continue;
        m[r * n + c] = std::sqrt(2.0 * t.k + 1.0) * AngCoeff::ClebGord(j, (double)t.k, j, mm, (double)t.q, mp);
      }
  }
  return m;
}
