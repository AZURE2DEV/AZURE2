// ChannelSolveAdjoint.h against central differences of forward functions, and its solves against the level route.
// No project file, nothing of the engine.
//  - the standard block (V = 1; six levels, three channels): formal_channel_solve + kernel_column +
//    formal_bilinear_adjoint against differences of solve_channel_space + f_column in every formal energy and width,
//    for L = Re(conj(beta) F_ab);
//  - a transformed block (K v = e N v from generalized_symmetric_eigen, Gamma = V^T G): the channel route's F and
//    level_solve_from_formal against the factorized level matrix M = K - E N - G Lambda G^T; and the whole reverse
//    chain formal_bilinear_adjoint -> resolve_formal_cotangents (Kbar, Nbar) and add_level_widths_from_formal (Gbar)
//    against differences, in every element of K, N and G, of the same L computed from that factorized M (which does
//    not pass through the eigenproblem);
//  - a pencil block (K tridiagonal, N diagonal, so that Qt = Z = 1): pencil_channel_solve against the factorized M.
// Differences: Richardson's (4 D(h/2) - D(h)) / 3 of the central difference D, h = 1e-4; the agreement reached is
// printed with each check.
//   physics_channel_solve_adjoint_test
#include <cmath>
#include <complex>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "azure2/physics/ChannelSolveAdjoint.h"
#include "azure2/physics/SymmetricEigen.h"

using namespace azure2::physics;
namespace cs = azure2::physics::channel_space;
typedef std::complex<double> complex;
typedef LevelMatrixSolver::matrix_c matrix_c;

namespace {
int failures = 0;
void check(const std::string &what, double worst, double tol) {
  const bool ok = worst <= tol;
  std::printf("  %s  %s  (largest relative difference %.2e; limit %.0e)\n", ok ? "ok  " : "FAIL", what.c_str(), worst, tol);
  if (!ok) failures++;
}

const int n = 6, nc = 3;
const double kEin = 3.1;
const complex kLo[nc] = {complex(-0.3, 0.8), complex(0.25, 0.35), complex(-0.7, 1.4)};
const double kE[n] = {0.8, 1.9, 2.4, 3.7, 5.1, 8.0};
const double kG[n][nc] = {{0.21, 0.02, 0.4},  {-0.35, 0.11, -0.2}, {0.08, -0.4, 0.3},
                          {0.5, 0.3, -0.15},  {-0.12, 0.25, 0.33}, {0.9, -0.7, 0.1}};
typedef std::vector<double> Params;

double difference(const std::function<double(const Params &)> &f, const Params &p, int i, double h) {
  auto D = [&](double step) {
    Params a = p, b = p;
    a[i] += step;
    b[i] -= step;
    return (f(a) - f(b)) / (2.0 * step);
  };
  return (4.0 * D(0.5 * h) - D(h)) / 3.0;
}

// --- the standard block: p = (e_k, Gamma_kc) -------------------------------------------------------------------
cs::FormalBlock standard_block(const Params &p) {
  cs::FormalBlock f;
  f.mode = 0;
  f.n = n;
  f.nch = nc;
  f.usable = true;
  f.identity = true;
  f.e.resize(n);
  f.elo.assign(n, 0.0);
  for (int l = 0; l < n; l++) {
    f.e[l] = p[l];
    f.emax = std::max(f.emax, std::fabs(p[l]));
  }
  f.gamma.assign(p.begin() + n, p.end());
  f.gcut = f.gamma;
  f.act = {0, 1, 2};
  f.pos = {0, 1, 2};
  f.lamShift.assign(nc, 0.0);
  return f;
}

// --- the transformed block: q = (K upper triangle, N upper triangle, G) -----------------------------------------
const int nt = n * (n + 1) / 2;
Params transformed_start() {
  Params q;
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++) q.push_back(i == j ? kE[i] : 0.07 * std::sin(1.3 * i + 2.1 * j));
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++) q.push_back(i == j ? 1.0 + 0.05 * i : 0.04 * std::cos(0.7 * i + 1.9 * j));
  for (int l = 0; l < n; l++)
    for (int c = 0; c < nc; c++) q.push_back(kG[l][c]);
  return q;
}
void unpack(const Params &q, std::vector<double> &K, std::vector<double> &N, std::vector<double> &G) {
  K.assign(n * n, 0.0);
  N.assign(n * n, 0.0);
  int t = 0;
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++) K[i * n + j] = K[j * n + i] = q[t++];
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++) N[i * n + j] = N[j * n + i] = q[t++];
  G.assign(q.begin() + 2 * nt, q.end());
}
// M = K - E N - G Lambda G^T, factorized.
void level_route(const std::vector<double> &K, const std::vector<double> &N, const std::vector<double> &G,
                 LevelMatrixSolver &solver) {
  matrix_c M(n, std::vector<complex>(n));
  for (int r = 0; r < n; r++)
    for (int c = 0; c < n; c++) {
      complex s = K[r * n + c] - kEin * N[r * n + c];
      for (int d = 0; d < nc; d++) s -= G[r * nc + d] * G[c * nc + d] * kLo[d];
      M[r][c] = s;
    }
  solver.Decompose(M);
}
complex level_bilinear(const std::vector<double> &K, const std::vector<double> &N, const std::vector<double> &G, int a,
                       int b, std::vector<complex> *ub = nullptr) {
  LevelMatrixSolver solver;
  level_route(K, N, G, solver);
  std::vector<double> ga(n);
  std::vector<complex> x(n);
  for (int l = 0; l < n; l++) {
    ga[l] = G[l * nc + a];
    x[l] = G[l * nc + b];
  }
  solver.Solve(x.data());
  if (ub) *ub = x;
  return width_bilinear(ga.data(), x.data(), n);
}
double relative(const std::vector<complex> &a, const std::vector<complex> &b) {
  double scale = 0.0, diff = 0.0;
  for (std::size_t i = 0; i < a.size(); i++) {
    scale = std::max(scale, std::abs(b[i]));
    diff = std::max(diff, std::abs(a[i] - b[i]));
  }
  return diff / scale;
}
}  // namespace

int main() {
  const double h = 1.0e-4;
  const complex beta(0.7, -1.3);
  const int cha = 0, chb = 2;
  const BilinearTerm term{cha + 1, chb + 1, beta};

  std::printf("the standard block: formal_bilinear_adjoint against differences of solve_channel_space + f_column\n");
  {
    Params p(kE, kE + n);
    for (int l = 0; l < n; l++)
      for (int c = 0; c < nc; c++) p.push_back(kG[l][c]);
    bool taken = true;
    auto L = [&](const Params &q) {
      const cs::FormalBlock f = standard_block(q);
      cs::ChannelSolve st;
      if (cs::solve_channel_space(f, kLo, kEin, st) != cs::Route::kTaken) taken = false;
      return std::real(std::conj(beta) * cs::f_column(f, st, chb)[cha]);
    };
    const cs::FormalBlock f = standard_block(p);
    cs::ChannelSolve st;
    cs::solve_channel_space(f, kLo, kEin, st);
    std::vector<complex> ua, ub, ya(nc), yb(nc);
    cs::formal_channel_solve(f, st, cha, ua);
    cs::formal_channel_solve(f, st, chb, ub);
    cs::kernel_column(f, st, cha, ya);
    cs::kernel_column(f, st, chb, yb);
    LevelCotangents bar(n, nc);
    cs::formal_bilinear_adjoint(term, ua.data(), ub.data(), ya.data(), yb.data(), n, nc, bar);
    double scale = 0.0, diff = 0.0;
    for (int i = 0; i < (int)p.size(); i++) {
      const double fd = difference(L, p, i, h);
      const double an = i < n ? bar.e[i] : bar.gamma[i - n];
      scale = std::max(scale, std::fabs(fd));
      diff = std::max(diff, std::fabs(an - fd));
    }
    check("every formal energy and width", diff / scale, 1.0e-9);
    if (!taken) check("the channel route was taken at every difference point", 1.0, 0.0);
  }

  std::printf("a transformed block (K v = e N v, Gamma = V^T G)\n");
  {
    const Params q0 = transformed_start();
    std::vector<double> K, N, G;
    unpack(q0, K, N, G);
    const EigenResult eig = generalized_symmetric_eigen(K, N, n);
    cs::FormalBlock f;
    f.mode = 1;
    f.n = n;
    f.nch = nc;
    f.usable = true;
    f.identity = false;
    f.adjointOK = true;
    f.e.resize(n);
    f.elo.assign(n, 0.0);
    f.V = eig.vectors;
    f.gcut = G;
    f.gamma.assign(n * nc, 0.0);
    for (int k = 0; k < n; k++) {
      f.e[k] = eig.values[k];
      f.emax = std::max(f.emax, std::fabs(eig.values[k]));
      for (int c = 0; c < nc; c++)
        for (int r = 0; r < n; r++) f.gamma[k * nc + c] += f.V[r * n + k] * G[r * nc + c];
    }
    f.act = {0, 1, 2};
    f.pos = {0, 1, 2};
    f.lamShift.assign(nc, 0.0);
    cs::ChannelSolve st;
    const bool taken = eig.ok && cs::solve_channel_space(f, kLo, kEin, st) == cs::Route::kTaken;
    check("the eigenproblem is solved and the channel route taken", taken ? 0.0 : 1.0, 0.0);

    // forward: F and A g_c against the level route
    std::vector<complex> Fch, Flev, ub_level;
    for (int a = 0; a < nc; a++)
      for (int b = 0; b < nc; b++) {
        Fch.push_back(cs::f_column(f, st, b)[a]);
        Flev.push_back(level_bilinear(K, N, G, a, b));
      }
    check("F = G^T A G: channel route against the factorized level matrix", relative(Fch, Flev), 1.0e-12);
    std::vector<complex> ua, ub, sa(n), sb(n), ya(nc), yb(nc);
    cs::formal_channel_solve(f, st, cha, ua);
    cs::formal_channel_solve(f, st, chb, ub);
    cs::level_solve_from_formal(f, ub, sb);
    level_bilinear(K, N, G, cha, chb, &ub_level);
    check("A g_c = V u_c (level_solve_from_formal) against the factorized level matrix", relative(sb, ub_level), 1.0e-12);

    // reverse: formal cotangents, resolved to K, N and G
    cs::kernel_column(f, st, cha, ya);
    cs::kernel_column(f, st, chb, yb);
    LevelCotangents formal(n, nc);
    cs::formal_bilinear_adjoint(term, ua.data(), ub.data(), ya.data(), yb.data(), n, nc, formal);
    std::vector<std::vector<double>> gbar(n, std::vector<double>(nc));
    for (int k = 0; k < n; k++)
      for (int c = 0; c < nc; c++) gbar[k][c] = formal.gamma[k * nc + c];
    cs::FormalResolve res;
    cs::resolve_formal_cotangents(f, formal.e, gbar, res);
    LevelCotangents level(n, nc);
    cs::add_level_widths_from_formal(f, gbar, level);
    auto L = [&](const Params &q) {
      std::vector<double> k2, n2, g2;
      unpack(q, k2, n2, g2);
      return std::real(std::conj(beta) * level_bilinear(k2, n2, g2, cha, chb));
    };
    double scale = 0.0, dK = 0.0, dN = 0.0, dG = 0.0;
    int t = 0;
    for (int part = 0; part < 2; part++)
      for (int i = 0; i < n; i++)
        for (int j = i; j < n; j++, t++) {
          const double fd = difference(L, q0, t, h);
          // Kbar = X1 + X2, Nbar = -X2; an off-diagonal element moves K_ij and K_ji together.
          auto bar = [&](int r, int c) {
            return part == 0 ? std::real(res.X1[r][c]) + std::real(res.X2[r][c]) : -std::real(res.X2[r][c]);
          };
          const double an = i == j ? bar(i, i) : bar(i, j) + bar(j, i);
          scale = std::max(scale, std::fabs(fd));
          (part == 0 ? dK : dN) = std::max(part == 0 ? dK : dN, std::fabs(an - fd));
        }
    for (int i = 0; i < n * nc; i++, t++) {
      const double fd = difference(L, q0, t, h);
      scale = std::max(scale, std::fabs(fd));
      dG = std::max(dG, std::fabs(level.gamma[i] - fd));
    }
    check("Kbar (resolve_formal_cotangents): every element of K", dK / scale, 1.0e-9);
    check("Nbar (resolve_formal_cotangents): every element of N", dN / scale, 1.0e-9);
    check("Gbar = V Gammabar (add_level_widths_from_formal): every width", dG / scale, 1.0e-9);
  }

  std::printf("a pencil block (K tridiagonal, N diagonal: Qt = Z = 1)\n");
  {
    std::vector<double> K(n * n, 0.0), N(n * n, 0.0), G(n * nc);
    for (int i = 0; i < n; i++) {
      K[i * n + i] = kE[i];
      N[i * n + i] = 1.0 + 0.05 * i;
      if (i + 1 < n) K[i * n + i + 1] = K[(i + 1) * n + i] = 0.06 * (i + 1);
      for (int c = 0; c < nc; c++) G[i * nc + c] = kG[i][c];
    }
    cs::FormalBlock f;
    f.mode = 1;
    f.n = n;
    f.nch = nc;
    f.usable = true;
    f.identity = false;
    f.pencil = true;
    f.H.assign(K.begin(), K.end());
    f.T.assign(N.begin(), N.end());
    f.Qt.assign(n * n, 0.0);
    for (int i = 0; i < n; i++) f.Qt[i * n + i] = 1.0;
    f.Z = f.Qt;
    f.gq = f.gz = f.gamma = f.gcut = G;
    f.act = {0, 1, 2};
    f.pos = {0, 1, 2};
    f.lamShift.assign(nc, 0.0);
    cs::ChannelSolve st;
    const bool taken = cs::solve_channel_space(f, kLo, kEin, st) == cs::Route::kTaken;
    check("the channel route is taken", taken ? 0.0 : 1.0, 0.0);
    double worst = 0.0;
    for (int b = 0; b < nc; b++) {
      std::vector<complex> solve(n), ref;
      cs::pencil_channel_solve(f, st, b, solve);
      level_bilinear(K, N, G, 0, b, &ref);
      worst = std::max(worst, relative(solve, ref));
    }
    check("A g_c (pencil_channel_solve) against the factorized level matrix, every channel", worst, 1.0e-12);
  }

  std::printf("%s\n", failures ? "FAILED" : "all checks passed");
  return failures ? 1 : 0;
}
