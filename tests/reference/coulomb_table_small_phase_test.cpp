// CoulombTable::LookupOpen deep below the Coulomb barrier.  The table stores the hard-sphere phase as phi - rho, so
// the looked-up phase is good to the rounding of rho (about 1e-16 rad) and no better; where phi itself is smaller
// (alpha+alpha below about 20 keV: phi = 1e-75 at 1 keV) the lookup returned 0 or +-1e-16 and the elastic amplitude
// exp(-2i phi) - 1 was rounding noise (angle-integrated alpha+alpha elastic at 1 keV: 5e-31 b for 6e-144 b).
// Below CoulombTable::kSmallPhase the lookup now takes the phase from the exact routine.  Checked on a small table
// built here: alpha+alpha at a = 3.5 fm from 1 to 200 keV, l = 0..2.
#include "CoulombTable.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;
static int failures = 0;
static void check(bool ok, const char *what) {
  std::printf("  %s  %s\n", ok ? "ok  " : "FAIL", what);
  if (!ok) failures++;
}

int main() {
  const fs::path dir = fs::temp_directory_path() / ("azure2_cphase_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  fs::remove_all(dir);
  const std::string file = (dir / "table.bin").string();

  CoulombTable::BuildSpec spec;  // a small table over the corner in question: seconds
  spec.lMax = 2;
  spec.etaMin = 1.0, spec.etaMax = 40.0;
  spec.rhoMin = 0.02, spec.rhoMax = 1.0;
  spec.perDecade = 40.0;
  spec.threads = 2;
  std::string error;
  check(CoulombTable::BuildFile(file, spec, error), "the table builds");
  CoulombTable t;
  check(t.Load(file, error), "the table loads");

  // alpha+alpha, a = 3.5 fm
  const double hbarc = 197.3269804, mu = 0.5 * 4.002603 * 931.494, a = 3.5, alpha = 1.0 / 137.035999;
  int small = 0, large = 0, missed = 0;
  double worstSmall = 0.0, worstLarge = 0.0, worstP = 0.0;
  bool positive = true;
  for (double keV : {1.0, 2.0, 5.0, 10.0, 14.0, 16.0, 18.0, 20.0, 30.0, 50.0, 92.0, 150.0, 200.0})
    for (int l = 0; l <= 2; ++l) {
      const double E = keV * 1e-3, rho = std::sqrt(2 * mu * E) / hbarc * a, eta = 4 * alpha * std::sqrt(mu / (2 * E));
      const CoulombTable::Open exact = CoulombTable::ExactOpen(l, eta, rho);
      CoulombTable::Open table;
      if (!t.LookupOpen(l, eta, rho, table, nullptr)) {
        ++missed;
        continue;
      }
      const double rel = std::abs(table.phi / exact.phi - 1.0);
      if (exact.phi < 0.5 * CoulombTable::kSmallPhase) {
        ++small;
        worstSmall = std::max(worstSmall, rel);
        if (!(table.phi > 0.0)) positive = false;
      } else if (exact.phi > 2.0 * CoulombTable::kSmallPhase) {
        ++large;
        worstLarge = std::max(worstLarge, std::abs(table.phi - exact.phi));
      }
      worstP = std::max(worstP, std::abs(table.P / exact.P - 1.0));
      std::printf("    %5.0f keV l=%d  eta %7.3f rho %.5f  phi exact %.6e  lookup %.6e\n", keV, l, eta, rho, exact.phi, table.phi);
    }
  std::printf("  small phases %d (worst relative %.2e), larger %d (worst absolute %.2e), P worst relative %.2e\n", small, worstSmall,
              large, worstLarge, worstP);
  check(missed == 0, "every point is inside the table");
  check(small >= 25 && large >= 3, "points on both sides of the limit");
  check(positive, "a small phase is positive, never 0 or rounding noise");
  check(worstSmall < 1e-12, "below the limit the phase is the exact routine's");
  check(worstLarge < 1e-9, "above the limit the phase is the table's, within its interpolation");
  check(worstP < 1e-5, "the penetrability is the table's throughout");

  fs::remove_all(dir);
  std::printf(failures ? "FAILED: %d\n" : "PASS\n", failures);
  return failures ? 1 : 0;
}
