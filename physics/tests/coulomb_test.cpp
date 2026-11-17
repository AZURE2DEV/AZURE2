// physics/coulomb (docs/upgrade/coulomb.md section 5): every algorithm against mpmath at 40 digits
// (data/coulomb.csv, from pycoulomb's tests/make_reference.py, copied with it), per region of (l, eta, rho):
//   eta = 0;  eta > 0 with rho < rho_tp (classically forbidden);  eta > 0 with rho >= rho_tp (oscillating),
// rho_tp = eta + sqrt(eta^2 + l(l+1)) the outer turning point.  The error of F, F', G, G' is relative where the
// functions do not oscillate, and relative to the amplitude sqrt(F^2 + G^2) (sqrt(F'^2 + G'^2)) where they do, so
// that a zero of F does not make it meaningless.  Also: the Wronskian F'G - FG' = 1; the algorithms against each
// other where all are accurate (identity test); sigma_l against mpmath and against the recurrence
// sigma_l = sigma_{l-1} + atan(eta/l); the argument checks; and the same bits from 8 threads walking the grid in
// different orders (any shared mutable state shows up as a mismatch).
//
// What it detects: an algorithm wired to the wrong source or adaptor, a changed sign or normalisation (the
// Wronskian and the mpmath values), an accuracy regression of more than a decade in any region, state shared
// between calls.  Usage: physics_coulomb_test <path to coulomb.csv> [verbose]
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "azure2/physics/coulomb/Coulomb.h"

using namespace azure2::physics::coulomb;

namespace {

struct Row {
  int l;
  double eta, rho;
  double v[5];  // F, F', G, G', sigma
};

std::vector<Row> Read(const std::string &path) {
  std::vector<Row> rows;
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    Row r{};
    std::string tok;
    ss >> r.l >> tok;
    r.eta = std::strtod(tok.c_str(), nullptr);
    ss >> tok;
    r.rho = std::strtod(tok.c_str(), nullptr);
    for (double &x : r.v) {
      ss >> tok;
      x = std::strtod(tok.c_str(), nullptr);
    }
    rows.push_back(r);
  }
  return rows;
}

enum Region { kNeutral, kForbidden, kAllowed, kNRegions };
const char *kRegionName[kNRegions] = {"eta = 0", "eta > 0, rho < rho_tp", "eta > 0, rho >= rho_tp"};

bool Oscillating(const Row &r) { return r.rho >= r.eta + std::sqrt(r.eta * r.eta + r.l * (r.l + 1.0)); }
Region RegionOf(const Row &r) { return r.eta == 0.0 ? kNeutral : (Oscillating(r) ? kAllowed : kForbidden); }

// Largest error of (F, F') (which = 0) or (G, G') (which = 1) of w against the reference values ref[0..3].
double Error(const CoulombWaves &w, const double *ref, bool osc, int which) {
  const double F = ref[0], Fp = ref[1], G = ref[2], Gp = ref[3];
  const double a = std::hypot(F, G), ap = std::hypot(Fp, Gp);
  const double e0 = which == 0 ? std::fabs(w.F - F) / (osc ? a : std::fabs(F)) : std::fabs(w.G - G) / (osc ? a : std::fabs(G));
  const double e1 = which == 0 ? std::fabs(w.Fp - Fp) / (osc ? ap : std::fabs(Fp)) : std::fabs(w.Gp - Gp) / (osc ? ap : std::fabs(Gp));
  return std::fmax(std::isnan(e0) ? INFINITY : e0, std::isnan(e1) ? INFINITY : e1);
}

// Bounds per algorithm and region: {F and F', G and G', |Wronskian - 1|, failed points}: the values measured on this
// grid (GCC 11.5, x86-64, -O3; pycoulomb README for coulsm, coul, coulfg; this test for coulcc) rounded up to the next
// power of ten and multiplied by ten, so that a regression of more than about a decade fails.  Entries above 1e-2
// mark where an algorithm is known to be inaccurate (physics/coulomb/README.md); they are bounds too, so that a change
// there is noticed.  coulsm's F at eta = 0 came from an unstable upward recurrence in l (wrong for l >= 2 at small
// rho; bound 1e47, Wronskian 1e31) until 2026-11-16, when F and F' became the Riccati-Bessel function (coulomb.md
// section 10; 7.4e-16, Wronskian 6.7e-16 measured; coulomb_neutral_test.cpp covers rho down to 1e-6).
struct Bound {
  double f, g, w;
  int failed;
};
const std::map<Method, std::array<Bound, kNRegions>> kBound = {
    //                 eta = 0                   eta > 0, rho < rho_tp       eta > 0, rho >= rho_tp
    {Method::coulsm, {{{1e-14, 1e-14, 1e-14, 0}, {1e-8, 1e-8, 1e-14, 0}, {1e-8, 1e-8, 1e-9, 0}}}},
    {Method::coul, {{{1e0, 1e-9, 1e-2, 0}, {1e-8, 1e-8, 1e-14, 0}, {1e-8, 1e-8, 1e-9, 0}}}},
    {Method::coulfg, {{{1e1, 1e1, 1e-14, 0}, {1e1, 1e2, 1e-13, 7}, {1e-13, 1e-12, 1e-13, 0}}}},
    {Method::coulcc, {{{1e-12, 1e-11, 1e-14, 0}, {1e-12, 1e-10, 1e-14, 0}, {1e-11, 1e-11, 1e-14, 0}}}},
};

const Method kMethods[4] = {Method::coulsm, Method::coul, Method::coulfg, Method::coulcc};

bool SameBits(const CoulombWaves &a, const CoulombWaves &b) {
  return std::memcmp(&a.F, &b.F, 5 * sizeof(double)) == 0 && a.status == b.status;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("usage: %s coulomb.csv [verbose]\n", argv[0]);
    return 2;
  }
  const bool verbose = argc > 2;
  const std::vector<Row> rows = Read(argv[1]);
  if (rows.size() != 1296) {
    std::printf("FAIL: reference file %s incomplete (%zu rows, 1296 expected)\n", argv[1], rows.size());
    return 1;
  }
  int failures = 0;

  // 1. Accuracy and Wronskian per algorithm and region.
  std::map<Method, std::vector<CoulombWaves>> res;
  std::printf("| Method | Region | F, F' | G, G' | F'G - FG' - 1 | failed |\n|---|---|---|---|---|---|\n");
  for (Method m : kMethods) {
    double ef[kNRegions] = {0, 0, 0}, eg[kNRegions] = {0, 0, 0}, ew[kNRegions] = {0, 0, 0};
    int nfail[kNRegions] = {0, 0, 0}, npts[kNRegions] = {0, 0, 0};
    std::vector<CoulombWaves> &rm = res[m];
    for (const Row &r : rows) {
      const CoulombWaves w = coulomb(r.l, r.eta, r.rho, m);
      rm.push_back(w);
      const Region g = RegionOf(r);
      ++npts[g];
      if (w.status != Status::ok) {
        ++nfail[g];
        if (!std::isnan(w.F) || !std::isnan(w.G) || !std::isnan(w.sigma)) {
          std::printf("FAIL: %s returned numbers with a failed status (l=%d eta=%g rho=%g)\n", method_name(m), r.l, r.eta, r.rho);
          ++failures;
        }
        if (verbose) std::printf("  %s failed at l=%d eta=%g rho=%g\n", method_name(m), r.l, r.eta, r.rho);
        continue;
      }
      ef[g] = std::fmax(ef[g], Error(w, r.v, Oscillating(r), 0));
      eg[g] = std::fmax(eg[g], Error(w, r.v, Oscillating(r), 1));
      ew[g] = std::fmax(ew[g], std::fabs(w.Fp * w.G - w.F * w.Gp - 1.0));
    }
    for (int g = 0; g < kNRegions; ++g) {
      std::printf("| `%s` | %s | %.1e | %.1e | %.1e | %d / %d |\n", method_name(m), kRegionName[g], ef[g], eg[g], ew[g], nfail[g], npts[g]);
      const Bound &b = kBound.at(m)[g];
      if (!(ef[g] <= b.f && eg[g] <= b.g && ew[g] <= b.w && nfail[g] <= b.failed)) {
        std::printf("FAIL: %s, %s above its bound (%.0e, %.0e, %.0e, %d)\n", method_name(m), kRegionName[g], b.f, b.g, b.w, b.failed);
        ++failures;
      }
    }
  }

  // 2. Identity test: in the oscillating region every algorithm is accurate, so every pair must agree.
  for (int i = 0; i < 4; ++i)
    for (int j = i + 1; j < 4; ++j) {
      double d = 0.0;
      for (size_t k = 0; k < rows.size(); ++k) {
        const CoulombWaves &a = res[kMethods[i]][k], &b = res[kMethods[j]][k];
        if (a.status != Status::ok || b.status != Status::ok || RegionOf(rows[k]) != kAllowed) continue;
        const double ref[4] = {a.F, a.Fp, a.G, a.Gp};
        d = std::fmax(d, std::fmax(Error(b, ref, true, 0), Error(b, ref, true, 1)));
      }
      std::printf("identity, oscillating region: %s - %s %.1e\n", method_name(kMethods[i]), method_name(kMethods[j]), d);
      if (!(d <= 1e-8)) {  // 8.3e-10 measured (the 1e-10 precision of cwfcomp)
        std::printf("FAIL: %s and %s differ by %.2e where both are accurate\n", method_name(kMethods[i]), method_name(kMethods[j]), d);
        ++failures;
      }
    }

  // 3. sigma_l: against mpmath, and the recurrence sigma_l - sigma_{l-1} = atan(eta / l); sigma_0(0) = 0.
  double smax = 0.0, rmax = 0.0;
  for (const Row &r : rows) {
    smax = std::fmax(smax, std::fabs(coulomb_phase(r.l, r.eta) - r.v[4]));
    if (r.l > 0) rmax = std::fmax(rmax, std::fabs(coulomb_phase(r.l, r.eta) - coulomb_phase(r.l - 1, r.eta) - std::atan(r.eta / r.l)));
  }
  std::printf("sigma_l: %.1e rad from mpmath, %.1e rad from the recurrence\n", smax, rmax);
  if (!(smax <= 1e-13) || !(rmax <= 1e-13) || coulomb_phase(0, 0.0) != 0.0) {  // 7.1e-15 measured
    std::printf("FAIL: sigma_l\n");
    ++failures;
  }

  // 4. Argument checks: status and NaN, for every algorithm.
  const struct {
    int l;
    double eta, rho;
  } bad[] = {{-1, 1.0, 1.0}, {10001, 1.0, 1.0}, {0, NAN, 1.0}, {0, 1.0, INFINITY}, {0, 1.0, 0.0}, {0, 1.0, -1.0}};
  for (Method m : kMethods)
    for (const auto &b : bad) {
      const CoulombWaves w = coulomb(b.l, b.eta, b.rho, m);
      if (w.status != Status::invalid_argument || !std::isnan(w.F) || !std::isnan(w.Fp) || !std::isnan(w.G) || !std::isnan(w.Gp) || !std::isnan(w.sigma)) {
        std::printf("FAIL: %s accepted l=%d eta=%g rho=%g\n", method_name(m), b.l, b.eta, b.rho);
        ++failures;
      }
    }
  Method m;
  if (!method_from_name("coulcc", m) || m != Method::coulcc || method_from_name("gsl", m)) {
    std::printf("FAIL: method names\n");
    ++failures;
  }

  // 5. Threads: 8 threads, each walking the grid from a different start, against the serial results above.
  std::vector<int> mismatch(8, 0);
  std::vector<std::thread> pool;
  for (int t = 0; t < 8; ++t)
    pool.emplace_back([&, t]() {
      const size_t n = rows.size();
      for (size_t i = 0; i < n; ++i) {
        const size_t k = (i + t * n / 8) % n;
        for (Method mm : kMethods)
          if (!SameBits(coulomb(rows[k].l, rows[k].eta, rows[k].rho, mm), res[mm][k])) ++mismatch[t];
      }
    });
  for (std::thread &th : pool) th.join();
  int nm = 0;
  for (int x : mismatch) nm += x;
  std::printf("threads: %d mismatches over 8 x %zu x 4 evaluations\n", nm, rows.size());
  if (nm != 0) {
    std::printf("FAIL: results depend on the thread or the order of calls\n");
    ++failures;
  }

  if (failures) std::printf("%d failure(s)\n", failures);
  else std::printf("physics/coulomb: all checks passed\n");
  return failures ? 1 : 0;
}
