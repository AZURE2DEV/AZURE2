// CoulFuncCache (classic): the process-wide memo of Coulomb functions, and the
// exact-energy table a key keeps once its near-energy memo is given up.
//
// A key whose energies do not recur is given up after 4096 queries with under
// 5 % hits (the varying-energy-shift case).  It used to release everything,
// so a later session in the same process recomputed all of it; it now keeps
// the values it is handed in a table looked up by the exact bits of E, whose
// hit is what recomputing would give.  Checked here, single-threaded:
//   1. a near energy (1e-13 MeV below a stored one) is a hit while the memo is
//      on, as before;
//   2. after 5000 distinct energies the key is given up (as before), and from
//      then on an energy is a hit only bit for bit: the stored value, never
//      a neighbour's;
//   3. the exact table is capped (32768) and dropped for good when it fills
//      without hits: the memory goes back, nothing is found afterwards;
//   4. another key (another radius) is untouched by all of it.

#include "CoulFuncCache.h"
#include <cmath>
#include <cstdio>

static CoulWaves waves(double e) {
  CoulWaves w;
  w.F = e;
  w.dF = 2.0 * e;
  w.G = 3.0 * e;
  w.dG = 4.0 * e;
  return w;
}

int main() {
  int fails = 0;
  auto check = [&](const char *what, bool ok) {
    std::printf("%s  %s\n", ok ? "ok  " : "FAIL", what);
    if (!ok) fails++;
  };
  CoulFuncCache cache;
  CoulFuncCache::CoulFuncKey key{2, 6, 1.6, 0, 4.15, 0};
  CoulFuncCache::CoulFuncKey other = key;
  other.radius = 5.0;
  CoulWaves out;

  // 1. The near-energy memo, while it is on.
  cache.AddCoulWaves(key, 0.5, waves(0.5));
  cache.AddCoulWaves(key, 1.0, waves(1.0));
  check("1. a stored energy is a hit", cache.TryGetCoulWaves(key, 1.0, out) && out.F == 1.0);
  check("1. 1e-13 MeV below a stored one is a hit while the memo is on (unchanged)",
        cache.TryGetCoulWaves(key, 1.0 - 1.0e-13, out) && out.F == 1.0);
  cache.AddCoulWaves(other, 1.0, waves(7.0));

  // 2. 5000 energies that never recur: the memo is given up.
  for (int i = 1; i <= 5000; i++) {
    const double e = 2.0 + 1.0e-3 * i;
    cache.TryGetCoulWaves(key, e, out);  // a miss, counted
    cache.AddCoulWaves(key, e, waves(e));
  }
  CoulFuncCache::Stats s = cache.GetStats();
  check("2. the key's near-energy memo is given up", s.disabledKeys == 1);
  // The value handed in last before and after the shutoff are found again at their exact energies.
  const double late = 2.0 + 1.0e-3 * 5000;
  check("2. after the shutoff an energy is found at its exact bits",
        cache.TryGetCoulWaves(key, late, out) && out.F == late && out.dG == 4.0 * late);
  check("2. ... and not 1e-13 MeV away (recomputed, as without the memo)",
        !cache.TryGetCoulWaves(key, late + 1.0e-13, out));
  check("2. the memo released at the shutoff is not found (recomputed, as before)",
        !cache.TryGetCoulWaves(key, 1.0, out));

  // 3. Fill the exact table without a hit: capped, then dropped.
  for (int i = 0; i < 40000; i++) {
    const double e = 100.0 + 1.0e-3 * i;
    cache.TryGetCoulWaves(key, e, out);
    cache.AddCoulWaves(key, e, waves(e));
  }
  s = cache.GetStats();
  check("3. the exact table is dropped once full without hits (only the other key's entry left)",
        s.entries == 1);
  check("3. nothing is found afterwards", !cache.TryGetCoulWaves(key, late, out));
  cache.AddCoulWaves(key, 1.5, waves(1.5));
  check("3. and nothing is stored", !cache.TryGetCoulWaves(key, 1.5, out));

  // 4. The other key.
  check("4. another key keeps its memo", cache.TryGetCoulWaves(other, 1.0, out) && out.F == 7.0);

  std::printf(fails ? "FAILED: %d\n" : "all passed\n", fails);
  return fails ? 1 : 0;
}
