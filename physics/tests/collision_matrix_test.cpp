// CollisionMatrix.h against closed forms written out here; no project file, no engine class:
//  - one level, one channel: A = 1 / (E_l - E - gamma^2 (S - B + i P)), U = Omega^2 (1 + 2 i P gamma^2 A), |U| = 1,
//    T = exp(2 i omega) - U;
//  - three levels, three channels: U symmetric and unitary, A the inverse of the filled matrix;
//  - the Brune form: the boundary condition cancels, the diagonal holds S_c(E) - S_c(E_l), the off-diagonal the chord
//    of S_c between the two levels; the Park form: the overlap J on the diagonal and the slope of the chord off it;
//  - R-matrix capture: i gamma^2 on a photon channel's diagonal;
//  - a level without widths at E = E_l: unit diagonal; a width below the cutoff couples nothing; an element that
//    is not finite is left zero.
//   physics_collision_matrix_test
#include <cmath>
#include <complex>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "azure2/physics/CollisionMatrix.h"

using namespace azure2::physics;
typedef std::complex<double> complex;

namespace {
int failures = 0;
void check(const std::string &what, bool ok, double got = 0.0, double want = 0.0) {
  std::printf("  %s  %s", ok ? "ok  " : "FAIL", what.c_str());
  if (!ok) std::printf("  (got %.17g, want %.17g)", got, want);
  std::printf("\n");
  if (!ok) failures++;
}
bool close(complex a, complex b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

// A block's arrays, owned; tables() lends them.
struct Block {
  int n = 0, nc = 0;
  std::vector<double> E, gamma, shift, d1, d2, park, boundary;
  std::vector<complex> lo;
  std::vector<char> rad;
  bool brune = false, isPark = false, rmc = false;
  Block(int levels, int channels) : n(levels), nc(channels) {
    E.assign(n, 0.0);
    gamma.assign((std::size_t)n * nc, 0.0);
    shift.assign((std::size_t)n * nc, 0.0);
    d1.assign((std::size_t)n * nc, 0.0);
    d2.assign((std::size_t)n * nc, 0.0);
    park.assign(n, 1.0);
    boundary.assign(nc, 0.0);
    lo.assign(nc, complex(0.0, 0.0));
    rad.assign(nc, 'P');
  }
  double &g(int l, int c) { return gamma[(std::size_t)l * nc + c]; }
  double &s(int l, int c) { return shift[(std::size_t)l * nc + c]; }
  LevelMatrixTables tables() const {
    LevelMatrixTables t;
    t.n = n;
    t.nc = nc;
    t.E = E.data();
    t.gamma = gamma.data();
    t.shift = shift.data();
    t.shift_d1 = d1.data();
    t.shift_d2 = d2.data();
    t.park_norm = park.data();
    t.lo = lo.data();
    t.boundary = boundary.data();
    t.rad_type = rad.data();
    t.brune = brune;
    t.park = isPark;
    t.rmc = rmc;
    return t;
  }
};
}  // namespace

int main() {
  const complex i(0.0, 1.0);

  // One level, one channel.
  {
    Block b(1, 1);
    const double El = 2.0, g = 0.3, S = 0.25, B = 0.35, P = 0.4, E = 1.7, omega = 0.6, phi = 0.9;
    b.E[0] = El;
    b.g(0, 0) = g;
    b.lo[0] = complex(S - B, P);
    const complex cphase = std::exp(i * omega), hphase = std::exp(-i * phi);
    const double sqrtP = std::sqrt(P);
    ChannelTablesAtE ch;
    ch.coulomb_phase = &cphase;
    ch.hard_sphere_phase = &hphase;
    ch.sqrt_penetrability = &sqrtP;
    const CollisionMatrix cm = collision_matrix(b.tables(), E, ch);
    const complex A = 1.0 / (El - E - g * g * complex(S - B, P));
    const complex U = std::exp(2.0 * i * (omega - phi)) * (1.0 + 2.0 * i * P * g * g * A);
    check("one level, one channel: A = 1 / (E_l - E - gamma^2 (S - B + i P))", close(cm.A[0], A, 1e-15),
          std::abs(cm.A[0]), std::abs(A));
    check("  U = exp(2 i (omega - phi)) (1 + 2 i P gamma^2 A)", close(cm.u(0, 0), U, 4e-16), std::abs(cm.u(0, 0) - U));
    check("  |U| = 1", std::fabs(std::abs(cm.u(0, 0)) - 1.0) < 4e-16, std::abs(cm.u(0, 0)), 1.0);
    check("  T = exp(2 i omega) - U", close(cm.t(0, 0), std::exp(2.0 * i * omega) - U, 1e-15));
    // The phase shift of the level: U = exp(2 i (omega - phi + delta)),
    // tan(delta) = P gamma^2 / (E_l - E - gamma^2 (S - B)).
    const double delta = std::atan2(P * g * g, El - E - g * g * (S - B));
    check("  U = exp(2 i (omega - phi + delta)), tan(delta) = Gamma / 2 / (E_l + Delta - E)",
          close(cm.u(0, 0), std::exp(2.0 * i * (omega - phi + delta)), 1e-15));
  }

  // Three levels, three channels: symmetry, unitarity, A M = 1.
  {
    Block b(3, 3);
    const double Es[3] = {1.2, 2.9, 4.4};
    const double gs[3][3] = {{0.31, -0.12, 0.05}, {0.22, 0.41, -0.3}, {-0.18, 0.09, 0.27}};
    const double S[3] = {-0.3, 0.12, -0.05}, P[3] = {0.7, 0.2, 1.3}, om[3] = {0.4, 0.0, 1.1}, ph[3] = {0.8, 0.3, 0.1};
    std::vector<complex> cp(3), hp(3);
    std::vector<double> sp(3);
    for (int l = 0; l < 3; l++) {
      b.E[l] = Es[l];
      for (int c = 0; c < 3; c++) b.g(l, c) = gs[l][c];
    }
    for (int c = 0; c < 3; c++) {
      b.lo[c] = complex(S[c], P[c]);
      cp[c] = std::exp(i * om[c]);
      hp[c] = std::exp(-i * ph[c]);
      sp[c] = std::sqrt(P[c]);
    }
    ChannelTablesAtE ch;
    ch.coulomb_phase = cp.data();
    ch.hard_sphere_phase = hp.data();
    ch.sqrt_penetrability = sp.data();
    const double E = 2.5;
    const CollisionMatrix cm = collision_matrix(b.tables(), E, ch);
    double asym = 0.0, unit = 0.0, inv = 0.0, tdef = 0.0;
    LevelMatrixSolver::matrix_c M;
    fill_level_matrix(b.tables(), E, M);
    for (int a = 0; a < 3; a++)
      for (int c = 0; c < 3; c++) {
        asym = std::max(asym, std::abs(cm.u(a, c) - cm.u(c, a)));
        complex s(0.0, 0.0), am(0.0, 0.0);
        for (int k = 0; k < 3; k++) {
          s += std::conj(cm.u(k, a)) * cm.u(k, c);
          am += cm.A[(std::size_t)a * 3 + k] * M[k][c];
        }
        unit = std::max(unit, std::abs(s - (a == c ? 1.0 : 0.0)));
        inv = std::max(inv, std::abs(am - (a == c ? 1.0 : 0.0)));
        tdef = std::max(tdef, std::abs(cm.t(a, c) - ((a == c ? cp[a] * cp[a] : complex(0.0, 0.0)) - cm.u(a, c))));
      }
    check("three levels, three channels: U symmetric", asym < 1e-15, asym);
    check("  U unitary (real widths, open channels)", unit < 2e-15, unit);
    check("  A times the filled matrix = 1", inv < 2e-15, inv);
    check("  T = exp(2 i omega_c) delta_cc' - U", tdef < 1e-15, tdef);
    check("  the standard matrix is (E_l - E) delta - sum_c gamma gamma (S - B + i P)",
          close(M[0][1],
                -(gs[0][0] * gs[1][0] * b.lo[0] + gs[0][1] * gs[1][1] * b.lo[1] + gs[0][2] * gs[1][2] * b.lo[2]),
                1e-15) &&
              close(M[2][2], Es[2] - E - (gs[2][0] * gs[2][0] * b.lo[0] + gs[2][1] * gs[2][1] * b.lo[1] +
                                           gs[2][2] * gs[2][2] * b.lo[2]),
                    1e-15));
  }

  // Brune and Park, two levels, one particle channel.
  {
    const double Ea = 1.0, Eb = 3.0, ga = 0.4, gb = -0.25, Sa = -0.2, Sb = 0.1, S = -0.07, P = 0.3, E = 1.9;
    for (double B : {0.0, -1.5}) {
      Block b(2, 1);
      b.brune = true;
      b.E[0] = Ea;
      b.E[1] = Eb;
      b.g(0, 0) = ga;
      b.g(1, 0) = gb;
      b.s(0, 0) = Sa;
      b.s(1, 0) = Sb;
      b.boundary[0] = B;
      b.lo[0] = complex(S - B, P);
      LevelMatrixSolver::matrix_c M;
      fill_level_matrix(b.tables(), E, M);
      const double chord = (Sa * (E - Eb) - Sb * (E - Ea)) / (Ea - Eb);
      const std::string tag = B == 0.0 ? "Brune, B = 0" : "Brune, B = -1.5";
      check(tag + ": diagonal E_l - E - gamma^2 (S(E) - S(E_l) + i P)",
            close(M[0][0], Ea - E - ga * ga * complex(S - Sa, P), 4e-16) &&
                close(M[1][1], Eb - E - gb * gb * complex(S - Sb, P), 4e-16));
      check(tag + ": off-diagonal -gamma_a gamma_b (S(E) - chord(E) + i P)",
            close(M[0][1], -ga * gb * complex(S - chord, P), 1e-13) && close(M[1][0], M[0][1], 1e-13),
            std::abs(M[0][1] + ga * gb * complex(S - chord, P)));
    }
    Block b(2, 1);
    b.isPark = true;
    b.E[0] = Ea;
    b.E[1] = Eb;
    b.g(0, 0) = ga;
    b.g(1, 0) = gb;
    b.s(0, 0) = Sa;
    b.s(1, 0) = Sb;
    b.park[0] = 0.9;
    b.park[1] = 1.2;
    b.boundary[0] = 0.4;
    b.lo[0] = complex(S - 0.4, P);
    LevelMatrixSolver::matrix_c M;
    fill_level_matrix(b.tables(), E, M);
    const double Jab = -ga * gb * (Sa - Sb) / (Ea - Eb);
    check("Park: diagonal (E_l - E) J_ll - gamma^2 (S(E) - S(E_l) + i P)",
          close(M[0][0], (Ea - E) * 0.9 - ga * ga * complex(S - Sa, P), 4e-16) &&
              close(M[1][1], (Eb - E) * 1.2 - gb * gb * complex(S - Sb, P), 4e-16));
    check("Park: off-diagonal (E_a - E) J_ab - gamma_a gamma_b (S(E) - S(E_a) + i P), J_ab = -gamma gamma slope",
          close(M[0][1], (Ea - E) * Jab - ga * gb * complex(S - Sa, P), 1e-15) &&
              close(M[1][0], (Eb - E) * Jab - ga * gb * complex(S - Sb, P), 1e-15));
  }

  // R-matrix capture, the decoupled level, the cutoff, a non-finite element.
  {
    Block b(2, 2);
    b.rmc = true;
    b.rad[1] = 'E';
    b.E[0] = 1.0;
    b.E[1] = 2.0;
    b.g(0, 0) = 0.5;
    b.g(0, 1) = 0.01;
    b.g(1, 0) = 0.5e-12;  // below the cutoff
    b.lo[0] = complex(-0.1, 0.2);
    b.lo[1] = complex(0.0, 1.0e-4);
    LevelMatrixSolver::matrix_c M;
    fill_level_matrix(b.tables(), 2.0, M);
    check("R-matrix capture: i gamma^2 on the photon channel's diagonal",
          close(M[0][0], 1.0 - 2.0 - 0.25 * b.lo[0] - 1e-4 * b.lo[1] - i * 1e-4, 1e-15));
    check("a width below 1e-12 couples nothing", M[0][1] == complex(0.0, 0.0) && M[1][0] == complex(0.0, 0.0));
    check("a level without widths at E = E_l: unit diagonal", M[1][1] == complex(1.0, 0.0));
    fill_level_matrix(b.tables(), 2.5, M);
    check("  and E_l - E away from it", M[1][1] == complex(-0.5, 0.0));
    b.lo[0] = complex(std::numeric_limits<double>::quiet_NaN(), 0.0);
    fill_level_matrix(b.tables(), 2.5, M);
    check("an element that is not finite is left zero", M[0][0] == complex(0.0, 0.0));
    b.g(0, 0) = 0.0;
    b.g(0, 1) = 0.0;
    b.g(1, 0) = 0.0;
    ChannelTablesAtE ch;
    const complex one[2] = {complex(1.0, 0.0), complex(1.0, 0.0)};
    const double sp[2] = {1.0, 1.0};
    ch.coulomb_phase = one;
    ch.hard_sphere_phase = one;
    ch.sqrt_penetrability = sp;
    b.E[0] = b.E[1] = 2.0;
    // A zero matrix cannot be made by the fill (unit diagonal): one by hand.
    LevelMatrixSolver solver;
    LevelMatrixSolver::matrix_c zero(1, LevelMatrixSolver::vector_c(1, complex(0.0, 0.0)));
    solver.Decompose(zero);
    check("LevelMatrixSolver: a zero matrix is singular", solver.singular());
    const CollisionMatrix cm = collision_matrix(b.tables(), 2.0, ch);
    check("no widths: U = 1, T = 0",
          !cm.singular && cm.u(0, 0) == complex(1.0, 0.0) && cm.u(0, 1) == complex(0.0, 0.0) &&
              cm.t(1, 1) == complex(0.0, 0.0));
  }

  // The element functions.
  {
    const ChannelPhases a{std::exp(i * 0.3), std::exp(-i * 0.7)}, c{std::exp(i * 1.3), std::exp(-i * 0.2)};
    const complex u = pathway_amplitude(0.5, 2.0, complex(0.3, -0.4));
    check("pathway_amplitude = 2 i sqrt(P) sqrt(P') g^T A g'", close(u, 2.0 * i * complex(0.3, -0.4), 1e-15));
    check("t_matrix_element off the diagonal = -U",
          t_matrix_element(false, a, c, u) == -collision_element(false, a, c, u));
    check("rmc_t_matrix_element = -i u", close(rmc_t_matrix_element(u), -i * u, 1e-16));
    check("external_capture_element = gamma sqrt(N_f) conversion amplitude",
          close(external_capture_element(0.2, 0.9, 3.0, complex(1.0, 2.0)), 0.54 * complex(1.0, 2.0), 4e-16));
    check("channel_capture_factor = 2 i sqrt(P) u^T A v",
          close(channel_capture_factor(0.5, complex(0.3, -0.4)), i * complex(0.3, -0.4), 1e-16));
    const double g[3] = {1.0, -2.0, 0.5};
    const complex x[3] = {complex(1.0, 1.0), complex(0.0, 2.0), complex(4.0, 0.0)};
    check("width_bilinear = sum_l g_l x_l", width_bilinear(g, x, 3) == complex(3.0, -3.0));
  }

  std::printf("%s: %d failure(s)\n", failures ? "FAIL" : "ok", failures);
  return failures ? 1 : 0;
}
