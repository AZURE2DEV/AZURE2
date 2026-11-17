// ChannelSolve.h against the level route of CollisionMatrix.h, and its LU routines against residuals:
//  - a standard block of six levels and two channels (the formal levels are the levels: V = 1): every column of
//    F = G^T A G, a product A u and a bilinear form u^T A v agree with the factorized level matrix;
//  - at a level energy and 1e-9 MeV from it the level is deflated (kept in the bordered system) and F still agrees;
//  - the route is refused (Route::kPole) for a narrow level exactly at the energy (the bordered solve would lose it);
//  - lu_factor / lu_solve, the real pair and the Hessenberg pair solve A x = b to round-off; a singular matrix
//    returns the pivot ratio 0.
//   physics_channel_solve_test
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

#include "azure2/physics/ChannelSolve.h"
#include "azure2/physics/CollisionMatrix.h"

using namespace azure2::physics;
namespace cs = azure2::physics::channel_space;
typedef std::complex<double> complex;

namespace {
int failures = 0;
void check(const std::string &what, bool ok, double got = 0.0, double want = 0.0) {
  std::printf("  %s  %s", ok ? "ok  " : "FAIL", what.c_str());
  if (!ok) std::printf("  (got %.17g, want %.17g)", got, want);
  std::printf("\n");
  if (!ok) failures++;
}

const int n = 6, nc = 2;
const double kE[n] = {0.8, 1.9, 2.4, 3.7, 5.1, 8.0};
const double kG[n][nc] = {{0.21, 0.02}, {-0.35, 0.11}, {0.08, -0.4}, {0.5, 0.3}, {-0.12, 0.25}, {0.9, -0.7}};

// The formal block of a standard level matrix: the levels themselves.
cs::FormalBlock standard_block(double g00 = kG[0][0]) {
  cs::FormalBlock f;
  f.mode = 0;
  f.n = n;
  f.nch = nc;
  f.usable = true;
  f.identity = true;
  f.e.resize(n);
  f.elo.assign(n, 0.0);
  f.gamma.resize((std::size_t)n * nc);
  f.gcut = f.gamma;
  for (int l = 0; l < n; l++) {
    f.e[l] = kE[l];
    f.emax = std::max(f.emax, std::fabs(kE[l]));
    for (int c = 0; c < nc; c++) f.gamma[(std::size_t)l * nc + c] = (l == 0 && c == 0) ? g00 : kG[l][c];
  }
  f.gcut = f.gamma;
  f.act = {0, 1};
  f.pos = {0, 1};
  f.lamShift.assign(nc, 0.0);
  return f;
}

// max over (a, b) of |F_ab(channel route) - g_a^T A g_b(level route)| / |g_a^T A g_b|, and whether the route was taken.
double compare(double E, const complex *lo, bool &taken, int &deflated) {
  const cs::FormalBlock f = standard_block();
  cs::ChannelSolve st;
  taken = cs::solve_channel_space(f, lo, E, st) == cs::Route::kTaken;
  deflated = (int)st.defl.size();
  if (!taken) return 0.0;
  std::vector<double> shift((std::size_t)n * nc, 0.0), bnd(nc, 0.0);
  const char rad[nc] = {'P', 'P'};
  LevelMatrixTables t;
  t.n = n;
  t.nc = nc;
  t.E = kE;
  t.gamma = &kG[0][0];
  t.shift = shift.data();
  t.lo = lo;
  t.boundary = bnd.data();
  t.rad_type = rad;
  LevelMatrixSolver::matrix_c M;
  fill_level_matrix(t, E, M);
  LevelMatrixSolver solver;
  solver.Decompose(M);
  double worst = 0.0;
  for (int b = 0; b < nc; b++) {
    std::vector<complex> x(n);
    for (int l = 0; l < n; l++) x[l] = kG[l][b];
    solver.Solve(x.data());
    const complex *col = cs::f_column(f, st, b);
    for (int a = 0; a < nc; a++) {
      double g[n];
      for (int l = 0; l < n; l++) g[l] = kG[l][a];
      const complex want = width_bilinear(g, x.data(), n);
      worst = std::max(worst, std::abs(col[a] - want) / std::abs(want));
    }
  }
  // A u and u^T A v for vectors that are not columns of G.
  std::vector<complex> u(n), v(n), au(n), want(n);
  for (int l = 0; l < n; l++) {
    u[l] = complex(0.3 * l - 1.0, 0.1 * l);
    v[l] = complex(1.0 / (l + 1.0), -0.2);
  }
  want = u;
  solver.Solve(want.data());
  cs::apply_a(f, st, u.data(), au.data());
  complex bil(0.0, 0.0), bilWant(0.0, 0.0);
  for (int l = 0; l < n; l++) {
    worst = std::max(worst, std::abs(au[l] - want[l]) / std::abs(want[l]));
    bilWant += v[l] * want[l];
  }
  bil = cs::channel_bilinear(f, st, v.data(), u.data());
  worst = std::max(worst, std::abs(bil - bilWant) / std::abs(bilWant));
  return worst;
}
}  // namespace

int main() {
  const complex lo[nc] = {complex(-0.15, 0.6), complex(0.05, 0.02)};
  bool taken = false;
  int deflated = 0;
  double d = compare(3.0, lo, taken, deflated);
  check("six levels, two channels, between the levels: F, A u and u^T A v as on the level route",
        taken && d < 1e-13, d);
  check("  no level deflated", deflated == 0, deflated, 0);
  d = compare(kE[3], lo, taken, deflated);
  check("at a level energy: the level is deflated, same numbers", taken && deflated == 1 && d < 1e-12, d);
  d = compare(kE[3] + 1.0e-9, lo, taken, deflated);
  check("1e-9 MeV from a level: deflated, same numbers", taken && deflated == 1 && d < 1e-9, d);

  // A deflated level is solved with the channels, which loses a level narrower than the round-off of that solve:
  // gamma^2 |Lambda| = 6e-7 MeV exactly at the energy is refused, and the caller takes the level matrix.
  {
    cs::FormalBlock f = standard_block(1.0e-3);
    f.gamma[1] = 0.0;
    cs::ChannelSolve st;
    check("a narrow level exactly at the energy: Route::kPole",
          cs::solve_channel_space(f, lo, kE[0], st) == cs::Route::kPole);
    check("  and taken 1e-3 MeV away", cs::solve_channel_space(f, lo, kE[0] + 1.0e-3, st) == cs::Route::kTaken);
  }

  // The LU routines.
  {
    const int m = 4;
    std::vector<complex> a(m * m), lu;
    std::vector<double> ar(m * m), lur;
    for (int r = 0; r < m; r++)
      for (int c = 0; c < m; c++) {
        a[r * m + c] = complex(std::sin(1.0 + 3.0 * r + c), std::cos(2.0 * r - c)) + (r == c ? complex(1.5, 0.5) : 0.0);
        ar[r * m + c] = std::sin(1.0 + 3.0 * r + c) + (r == c ? 0.5 : 0.0);
      }
    std::vector<complex> x(m), b(m);
    for (int r = 0; r < m; r++) b[r] = x[r] = complex(r + 1.0, -0.5 * r);
    lu = a;
    std::vector<int> piv(m);
    const double ratio = cs::lu_factor(lu.data(), piv.data(), m);
    cs::lu_solve(lu.data(), piv.data(), m, x.data());
    double res = 0.0;
    for (int r = 0; r < m; r++) {
      complex s(0.0, 0.0);
      for (int c = 0; c < m; c++) s += a[r * m + c] * x[c];
      res = std::max(res, std::abs(s - b[r]));
    }
    check("lu_factor / lu_solve: A x = b", ratio > 0.0 && ratio <= 1.0 && res < 1e-14, res);
    std::vector<double> xr(m * 2), br(m * 2);
    for (int r = 0; r < m; r++) {
      br[r * 2] = xr[r * 2] = r - 1.5;
      br[r * 2 + 1] = xr[r * 2 + 1] = 1.0 / (r + 1.0);
    }
    lur = ar;
    const double rr = cs::lu_factor_real(lur.data(), piv.data(), m);
    cs::lu_solve_real(lur.data(), piv.data(), m, xr.data(), 2);
    res = 0.0;
    for (int r = 0; r < m; r++)
      for (int k = 0; k < 2; k++) {
        double s = 0.0;
        for (int c = 0; c < m; c++) s += ar[r * m + c] * xr[c * 2 + k];
        res = std::max(res, std::fabs(s - br[r * 2 + k]));
      }
    check("lu_factor_real / lu_solve_real: two right-hand sides", rr > 0.0 && res < 1e-14, res);
    // Upper Hessenberg: zero below the first subdiagonal.
    std::vector<double> h = ar;
    for (int r = 0; r < m; r++)
      for (int c = 0; c + 1 < r; c++) h[r * m + c] = 0.0;
    std::vector<double> hl = h;
    std::vector<char> sw(m);
    xr = br;
    const double hr = cs::hessenberg_lu(hl.data(), sw.data(), m);
    cs::hessenberg_solve(hl.data(), sw.data(), m, xr.data(), 2);
    res = 0.0;
    for (int r = 0; r < m; r++)
      for (int k = 0; k < 2; k++) {
        double s = 0.0;
        for (int c = 0; c < m; c++) s += h[r * m + c] * xr[c * 2 + k];
        res = std::max(res, std::fabs(s - br[r * 2 + k]));
      }
    check("hessenberg_lu / hessenberg_solve: two right-hand sides", hr > 0.0 && res < 1e-14, res);
    std::vector<complex> sing(4, complex(1.0, 0.0));
    check("lu_factor of a singular matrix: pivot ratio 0", cs::lu_factor(sing.data(), piv.data(), 2) == 0.0);
    check("lu_factor of an empty matrix: 1", cs::lu_factor(nullptr, nullptr, 0) == 1.0);
  }

  std::printf("%s: %d failure(s)\n", failures ? "FAIL" : "ok", failures);
  return failures ? 1 : 0;
}
