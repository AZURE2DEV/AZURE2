// CollisionMatrixAdjoint.h against central differences of the forward functions of CollisionMatrix.h, on a block of
// four levels and three particle channels whose shift functions are cubic polynomials of the energy (so that
// S(E_l), S'(E_l) and S''(E_l) are exact).  No project file, nothing of the engine.
//  - level_matrix_adjoint: for a random complex cotangent Mbar, d/dE_l and d/dgamma_lc of
//    L = sum_rc Re(conj(Mbar_rc) M_rc) with M from fill_level_matrix, for the standard, the Brune and the Park
//    matrix; and Brune with two levels 4e-7 MeV apart (the Taylor forms of LevelPairShift);
//  - inverse_adjoint: d/dRe M_rc and d/dIm M_rc of L = sum Re(conj(Abar) A), A = M^-1 (LevelMatrixSolver);
//  - the chain of one bilinear form, L = Re(conj(beta) g_a^T A g_b): add_explicit_widths, add_rank_one_cotangent and
//    level_matrix_adjoint against differences of fill_level_matrix + LevelMatrixSolver + width_bilinear, in the three
//    formalisms; deferred_bilinear_adjoint (K and N summed, resolved with a zero kernel) gives the same numbers;
//  - level_parameter_adjoint (standard matrix: dK = e_l e_l^T for an energy, the kernel alone for a width) gives the
//    same numbers again, parameter by parameter;
//  - the step from a T-matrix element to its bilinear form (amplitude_cotangent, penetrability_factor,
//    bilinear_cotangent) against differences of t_matrix_element(pathway_amplitude(b)) in Re b and Im b;
//  - the dense adjoint (dense_bilinear_adjoint, then inverse_adjoint and level_matrix_adjoint) against differences of
//    masked_amplitude, without and with a masked level;
//  - external capture: the cotangents of the final level's widths against differences of external_capture_element
//    with the normalisation (1 + sum_d A_d gamma_d^2)^-1/2.
// Differences: Richardson's (4 D(h/2) - D(h)) / 3 of the central difference D, h = 1e-4 (MeV or MeV^1/2); the
// agreement reached is printed with each check.
//   physics_collision_matrix_adjoint_test
#include <cmath>
#include <complex>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "azure2/physics/CollisionMatrixAdjoint.h"

using namespace azure2::physics;
typedef std::complex<double> complex;
typedef LevelMatrixSolver::matrix_c matrix_c;

namespace {
int failures = 0;
// |got - want| <= tol * scale; prints the agreement reached.
void check(const std::string &what, double worst, double tol) {
  const bool ok = worst <= tol;
  std::printf("  %s  %s  (largest difference %.2e of the largest derivative; limit %.0e)\n", ok ? "ok  " : "FAIL",
              what.c_str(), worst, tol);
  if (!ok) failures++;
}

const int n = 4, nc = 3;
const double kEin = 2.2;  // the energy of the matrix (MeV)
const double kB[nc] = {-0.4, 0.3, -1.1};
const double kPen[nc] = {0.8, 0.35, 1.7};
// S_c(E) = s0 + s1 E + s2 E^2 + s3 E^3
const double kS[nc][4] = {{-0.6, 0.21, -0.03, 0.004}, {0.2, -0.15, 0.05, -0.006}, {-1.3, 0.4, 0.02, 0.003}};
double S(int c, double E) { return kS[c][0] + E * (kS[c][1] + E * (kS[c][2] + E * kS[c][3])); }
double S1(int c, double E) { return kS[c][1] + E * (2.0 * kS[c][2] + E * 3.0 * kS[c][3]); }
double S2(int c, double E) { return 2.0 * kS[c][2] + 6.0 * kS[c][3] * E; }

enum Mode { kStandard, kBrune, kPark };

// The parameters: p[l] = E_l, p[n + l * nc + c] = gamma_lc.
typedef std::vector<double> Params;
Params start(bool coincident = false) {
  Params p = {0.9, 1.8, 2.9, 4.6, 0.31, -0.12, 0.44, -0.27, 0.52, 0.18, 0.61, 0.23, -0.35, -0.42, 0.29, 0.5};
  if (coincident) p[2] = p[1] + 4.0e-7;
  return p;
}

// Everything a block's functions read, built from the parameters.
struct Block {
  std::vector<double> E, G, shift, d1, d2, sg, inv, park;
  std::vector<std::vector<double>> sd;
  std::vector<complex> lo;
  std::vector<char> rad, shiftCh;
  LevelMatrixTables forward(Mode m) const {
    LevelMatrixTables t;
    t.n = n;
    t.nc = nc;
    t.E = E.data();
    t.gamma = G.data();
    t.shift = shift.data();
    t.shift_d1 = d1.data();
    t.shift_d2 = d2.data();
    t.park_norm = park.data();
    t.lo = lo.data();
    t.boundary = kB;
    t.rad_type = rad.data();
    t.brune = (m == kBrune);
    t.park = (m == kPark);
    return t;
  }
  LevelMatrixAdjointTables adjoint(Mode m, const complex *kernel, bool withBoundary) const {
    LevelMatrixAdjointTables t;
    t.n = n;
    t.nc = nc;
    t.E = E.data();
    t.gamma = G.data();
    t.shift = shift.data();
    t.shift_gamma = sg.data();
    t.inv_gap = inv.data();
    t.shift_d1 = d1.data();
    t.shift_d2 = d2.data();
    t.park_norm = park.data();
    if (m == kBrune) t.brune_shift_d1 = &sd;
    t.kernel = kernel;
    t.boundary = withBoundary ? kB : nullptr;
    t.shift_channel = shiftCh.data();
    t.brune = (m == kBrune);
    t.park = (m == kPark);
    return t;
  }
};
Block block(const Params &p, Mode m) {
  Block b;
  b.E.assign(p.begin(), p.begin() + n);
  b.G.assign(p.begin() + n, p.end());
  b.shift.resize(n * nc);
  b.d1.resize(n * nc);
  b.d2.resize(n * nc);
  b.sg.resize(n * nc);
  b.park.resize(n);
  b.sd.assign(n, std::vector<double>(nc));
  for (int l = 0; l < n; l++) {
    double J = 1.0;
    for (int c = 0; c < nc; c++) {
      const int i = l * nc + c;
      b.shift[i] = S(c, b.E[l]);
      b.d1[i] = b.sd[l][c] = S1(c, b.E[l]);
      b.d2[i] = S2(c, b.E[l]);
      b.sg[i] = b.shift[i] * b.G[i];
      J -= b.G[i] * b.G[i] * b.d1[i];
    }
    b.park[l] = J;
  }
  b.inv.assign(n * n, 0.0);
  for (int r = 0; r < n; r++)
    for (int c = 0; c < n; c++)
      if (r != c) b.inv[r * n + c] = inverse_level_gap(b.E[r], b.E[c]);
  b.lo.resize(nc);
  for (int c = 0; c < nc; c++) b.lo[c] = complex(S(c, kEin) - kB[c], kPen[c]);
  b.rad.assign(nc, 'P');
  b.shiftCh.assign(nc, m == kStandard ? 0 : 1);
  return b;
}

// A fixed "random" complex number per index.
complex noise(int i) { return complex(std::sin(1.7 * i + 0.3), std::cos(2.9 * i + 1.1)); }

// Richardson's central difference of f in coordinate i.
double difference(const std::function<double(const Params &)> &f, Params p, int i, double h) {
  auto D = [&](double step) {
    Params a = p, b = p;
    a[i] += step;
    b[i] -= step;
    return (f(a) - f(b)) / (2.0 * step);
  };
  return (4.0 * D(0.5 * h) - D(h)) / 3.0;
}

// The largest |adjoint - difference| over the parameters, in units of the largest |difference|.
double worst(const std::function<double(const Params &)> &f, const Params &p, const LevelCotangents &bar, double h,
             bool report = false) {
  double scale = 0.0, diff = 0.0;
  for (int i = 0; i < (int)p.size(); i++) {
    const double fd = difference(f, p, i, h);
    const double an = i < n ? bar.e[i] : bar.gamma[i - n];
    if (report) std::printf("      %2d  adjoint % .12e  difference % .12e\n", i, an, fd);
    scale = std::max(scale, std::fabs(fd));
    diff = std::max(diff, std::fabs(an - fd));
  }
  return diff / scale;
}

matrix_c level_matrix(const Params &p, Mode m) {
  const Block b = block(p, m);
  matrix_c M;
  fill_level_matrix(b.forward(m), kEin, M);
  return M;
}

// g_a^T A g_b of the block.
complex bilinear(const Params &p, Mode m, int a, int bch, std::vector<complex> *ua = nullptr,
                 std::vector<complex> *ub = nullptr) {
  const matrix_c M = level_matrix(p, m);
  LevelMatrixSolver solver;
  solver.Decompose(M);
  std::vector<double> ga(n), gb(n);
  std::vector<complex> xa(n), xb(n);
  for (int l = 0; l < n; l++) {
    ga[l] = p[n + l * nc + a];
    gb[l] = p[n + l * nc + bch];
    xa[l] = ga[l];
    xb[l] = gb[l];
  }
  solver.Solve(xa.data());
  solver.Solve(xb.data());
  if (ua) *ua = xa;
  if (ub) *ub = xb;
  return width_bilinear(ga.data(), xb.data(), n);
}

const char *name(Mode m) { return m == kStandard ? "standard" : m == kBrune ? "Brune" : "Park"; }
}  // namespace

int main() {
  const double h = 1.0e-4;
  matrix_c Mbar(n, std::vector<complex>(n));
  for (int r = 0; r < n; r++)
    for (int c = 0; c < n; c++) Mbar[r][c] = noise(r * n + c);

  std::printf("level_matrix_adjoint against differences of fill_level_matrix\n");
  for (Mode m : {kStandard, kBrune, kPark}) {
    const Params p = start();
    auto L = [&](const Params &q) {
      const matrix_c M = level_matrix(q, m);
      double s = 0.0;
      for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++) s += std::real(std::conj(Mbar[r][c]) * M[r][c]);
      return s;
    };
    const Block b = block(p, m);
    LevelCotangents bar(n, nc);
    LevelMatrixAdjointScratch scratch;
    level_matrix_adjoint(b.adjoint(m, b.lo.data(), true), Mbar, kEin, scratch, bar);
    check(std::string(name(m)) + ": every E_l and gamma_lc", worst(L, p, bar, h), 1.0e-9);
  }
  {
    // Two levels 4e-7 MeV apart: the level matrix takes the Taylor forms of the chord, and so does the adjoint, to
    // first order in the gap times S''' (LevelPairShift.h).  The step keeps the pair coincident.
    const Params p = start(true);
    auto L = [&](const Params &q) {
      const matrix_c M = level_matrix(q, kBrune);
      double s = 0.0;
      for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++) s += std::real(std::conj(Mbar[r][c]) * M[r][c]);
      return s;
    };
    const Block b = block(p, kBrune);
    LevelCotangents bar(n, nc);
    LevelMatrixAdjointScratch scratch;
    level_matrix_adjoint(b.adjoint(kBrune, b.lo.data(), true), Mbar, kEin, scratch, bar);
    check("Brune, two coincident levels (Taylor forms)", worst(L, p, bar, 1.0e-7), 1.0e-6);
  }

  std::printf("inverse_adjoint against differences of LevelMatrixSolver's inverse\n");
  {
    const matrix_c M0 = level_matrix(start(), kStandard);
    matrix_c Abar(n, std::vector<complex>(n));
    for (int r = 0; r < n; r++)
      for (int c = 0; c < n; c++) Abar[r][c] = noise(40 + r * n + c);
    auto L = [&](const matrix_c &M) {
      LevelMatrixSolver solver;
      solver.Decompose(M);
      const matrix_c &A = solver.Inverse();
      double s = 0.0;
      for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++) s += std::real(std::conj(Abar[r][c]) * A[r][c]);
      return s;
    };
    LevelMatrixSolver solver;
    solver.Decompose(M0);
    const matrix_c bar = inverse_adjoint(solver.Inverse(), Abar);
    double scale = 0.0, diff = 0.0;
    for (int r = 0; r < n; r++)
      for (int c = 0; c < n; c++)
        for (int part = 0; part < 2; part++) {
          const complex unit = part ? complex(0.0, 1.0) : complex(1.0, 0.0);
          auto D = [&](double step) {
            matrix_c a = M0, b = M0;
            a[r][c] += step * unit;
            b[r][c] -= step * unit;
            return (L(a) - L(b)) / (2.0 * step);
          };
          const double fd = (4.0 * D(0.5 * h) - D(h)) / 3.0;
          const double an = part ? std::imag(bar[r][c]) : std::real(bar[r][c]);
          scale = std::max(scale, std::fabs(fd));
          diff = std::max(diff, std::fabs(an - fd));
        }
    check("Mbar = -A^H Abar A^H: Re and Im of every element", diff / scale, 1.0e-9);
  }

  std::printf("one bilinear form, L = Re(conj(beta) g_a^T A g_b), against differences of the forward chain\n");
  const complex beta(0.7, -1.3);
  const int cha = 0, chb = 2;
  for (Mode m : {kStandard, kBrune, kPark}) {
    const Params p = start();
    auto L = [&](const Params &q) { return std::real(std::conj(beta) * bilinear(q, m, cha, chb)); };
    std::vector<complex> ua, ub;
    bilinear(p, m, cha, chb, &ua, &ub);
    const Block b = block(p, m);
    const BilinearTerm term{cha + 1, chb + 1, beta};

    // (a) the rank-one cotangent of M and the structural adjoint
    LevelCotangents bar(n, nc);
    matrix_c Mb(n, std::vector<complex>(n, complex(0.0, 0.0)));
    add_explicit_widths(term, ua.data(), ub.data(), n, bar);
    add_rank_one_cotangent(beta, ua.data(), ub.data(), n, Mb);
    LevelMatrixAdjointScratch scratch;
    level_matrix_adjoint(b.adjoint(m, b.lo.data(), true), Mb, kEin, scratch, bar);
    check(std::string(name(m)) + ": explicit widths + rank-one Mbar + level_matrix_adjoint", worst(L, p, bar, h), 1.0e-9);

    // (b) the deferred form: the kernel Lambda (+ B in the shift channels) now, K and N resolved afterwards with a
    // zero kernel (the shifts at the level energies belong to K and N).
    {
      LevelCotangents def(n, nc);
      std::vector<double> K(n * n, 0.0), N(n * n, 0.0);
      std::vector<complex> lam(nc), ga(nc), gb(nc);
      for (int c = 0; c < nc; c++) lam[c] = b.lo[c] + (m == kStandard ? 0.0 : kB[c]);
      add_explicit_widths(term, ua.data(), ub.data(), n, def);
      deferred_bilinear_adjoint(term, ua.data(), ub.data(), b.G.data(), lam.data(), kEin, n, nc, K.data(), N.data(), ga,
                                gb, def);
      matrix_c X1(n, std::vector<complex>(n)), X2(n, std::vector<complex>(n));
      for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++) {
          X1[r][c] = K[r * n + c] - N[r * n + c];
          X2[r][c] = N[r * n + c];
        }
      const std::vector<complex> zero(nc, complex(0.0, 0.0));
      level_matrix_adjoint(b.adjoint(m, zero.data(), false), X1, 0.0, scratch, def);
      level_matrix_adjoint(b.adjoint(m, zero.data(), false), X2, 1.0, scratch, def);
      check(std::string(name(m)) + ": deferred_bilinear_adjoint, K and N resolved with a zero kernel",
            worst(L, p, def, h), 1.0e-9);
    }

    // (c) parameter by parameter (standard: dK = e_l e_l^T for E_l, nothing of K and N for a width)
    if (m == kStandard) {
      std::vector<LevelParameterDerivative> entries;
      for (int l = 0; l < n; l++) {
        LevelParameterDerivative e{l, -1, false, std::vector<double>(n, 0.0), std::vector<double>(n, 0.0)};
        e.k[l] = 0.5;
        entries.push_back(e);
        for (int c = 0; c < nc; c++) entries.push_back(LevelParameterDerivative{l, c, true, {}, {}});
      }
      const std::vector<int> slot = {0, -1, 1};
      const std::vector<const std::vector<complex> *> u = {&ua, &ub};
      std::vector<complex> gu, ku;
      width_solve_products(b.G.data(), u, n, nc, gu);
      LevelCotangents par(n, nc);
      level_parameter_adjoint(entries, std::vector<BilinearTerm>{term}, slot, u, gu, b.lo.data(), kEin, n, nc, ku,
                              [&](int l, int ch, double v) {
                                if (ch < 0)
                                  par.energy(l, v);
                                else
                                  par.width(l, ch, v);
                              });
      check("standard: level_parameter_adjoint, parameter by parameter", worst(L, p, par, h), 1.0e-9);
    }
  }

  std::printf("the T-matrix element of a pathway, L = Re(conj(tbar) T), against differences of t_matrix_element\n");
  {
    const ChannelPhases in{std::polar(1.0, 0.4), std::polar(1.0, -0.9)}, out{std::polar(1.0, -1.2), std::polar(1.0, 0.3)};
    const double spIn = 0.83, spOut = 1.21;
    const complex tbar(-0.6, 0.45), b0(0.37, -0.52);
    for (int diagonal = 0; diagonal < 2; diagonal++) {
      auto L = [&](complex b) {
        return std::real(std::conj(tbar) * t_matrix_element(diagonal, in, diagonal ? in : out,
                                                            pathway_amplitude(spIn, diagonal ? spIn : spOut, b)));
      };
      const complex bb = bilinear_cotangent(penetrability_factor(spIn, diagonal ? spIn : spOut),
                                            amplitude_cotangent(in, diagonal ? in : out, tbar));
      const complex bc = bilinear_cotangent(penetrability_factor(complex(spIn), complex(diagonal ? spIn : spOut)),
                                            amplitude_cotangent(in, diagonal ? in : out, tbar));
      const double dRe = (L(b0 + h) - L(b0 - h)) / (2.0 * h), dIm = (L(b0 + complex(0.0, h)) - L(b0 - complex(0.0, h))) / (2.0 * h);
      check(diagonal ? "diagonal pathway: amplitude_cotangent, penetrability_factor, bilinear_cotangent"
                     : "off-diagonal pathway: the same",
            std::max(std::abs(bb - complex(dRe, dIm)), std::abs(bc - complex(dRe, dIm))) / std::abs(bb), 1.0e-9);
    }
  }

  std::printf("the dense adjoint, L = Re(conj(ubar) u), u = sum factor g_l g_m A_lm, against differences\n");
  for (int masked = 0; masked < 2; masked++) {
    const Params p = start();
    const complex factor(0.0, 1.7), ubar(0.9, -0.4);
    std::vector<unsigned char> zero(n * nc, 0);
    const int maskChannel = masked ? 1 : -1;
    if (masked) zero[2 * nc + 1] = 1;  // level 2 has no input width in channel 1: left out of both indices
    auto L = [&](const Params &q) {
      LevelMatrixSolver solver;
      solver.Decompose(level_matrix(q, kStandard));
      return std::real(std::conj(ubar) * masked_amplitude(factor, &q[n], nc, cha, chb, solver.Inverse(), zero.data(), maskChannel));
    };
    const Block b = block(p, kStandard);
    LevelMatrixSolver solver;
    solver.Decompose(level_matrix(p, kStandard));
    const matrix_c &A = solver.Inverse();
    LevelCotangents bar(n, nc);
    matrix_c Abar(n, std::vector<complex>(n, complex(0.0, 0.0)));
    dense_bilinear_adjoint(factor, b.G.data(), nc, cha, chb, A, zero.data(), maskChannel, ubar, Abar, bar);
    LevelMatrixAdjointScratch scratch;
    level_matrix_adjoint(b.adjoint(kStandard, b.lo.data(), true), inverse_adjoint(A, Abar), kEin, scratch, bar);
    check(masked ? "with a masked level: masked_amplitude, dense_bilinear_adjoint"
                 : "dense_bilinear_adjoint + inverse_adjoint + level_matrix_adjoint",
          worst(L, p, bar, h), 1.0e-9);
    if (!masked) {
      // without a mask the amplitude is factor g_a^T A g_b
      const complex u = masked_amplitude(factor, &p[n], nc, cha, chb, A, zero.data(), -1);
      check("masked_amplitude without a mask is factor g_a^T A g_b",
            std::abs(u - factor * bilinear(p, kStandard, cha, chb)) / std::abs(u), 1.0e-13);
    }
  }

  std::printf("external capture, T = gamma_f s c a with s = (1 + sum_d A_d gamma_d^2)^-1/2, against differences\n");
  {
    const double Ad[nc] = {0.8, 0.0, 2.3}, conv = 1.4;
    const complex amp(0.6, -1.1), tbar(0.35, 0.8);
    const int fc = 2;
    Params g = {0.4, -0.7, 0.55};
    auto norm = [&](const Params &q) {
      double sum = 1.0;
      for (int d = 0; d < nc; d++) sum += Ad[d] * q[d] * q[d];
      return 1.0 / std::sqrt(sum);
    };
    auto L = [&](const Params &q) { return std::real(std::conj(tbar) * external_capture_element(q[fc], norm(q), conv, amp)); };
    const double s = norm(g);
    std::vector<double> bar(nc, 0.0);
    bar[fc] += external_capture_width_cotangent(tbar, s, conv, amp);
    const complex dTds = external_capture_norm_derivative(g[fc], conv, amp);
    for (int d = 0; d < nc; d++) bar[d] += bound_norm_width_cotangent(tbar, dTds, Ad[d], g[d], s * s * s);
    double scale = 0.0, diff = 0.0;
    for (int d = 0; d < nc; d++) {
      const double fd = difference(L, g, d, h);
      scale = std::max(scale, std::fabs(fd));
      diff = std::max(diff, std::fabs(bar[d] - fd));
    }
    check("external_capture_width_cotangent, bound_norm_width_cotangent", diff / scale, 1.0e-9);
  }

  std::printf("%s\n", failures ? "FAILED" : "all checks passed");
  return failures ? 1 : 0;
}
