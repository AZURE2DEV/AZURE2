// CoulFuncCache::TransientScope on several threads at once (2026-11-05).
//
// A transient scope (CNuc::TransformOut) memoizes Coulomb functions in a
// per-thread scratch map that the end of the scope discards.  The scope depth
// used to be one counter for all threads, and the end of the last scope on ANY
// thread cleared EVERY thread's scratch while the others were inserting into or
// reading theirs: heap corruption in natural-neon fits at four threads, where
// every pooled copy transforms on every evaluation (EData::ReplaceSubPointGrids).
//
// Each thread here opens scopes over and over, memoizes a few hundred energies
// and reads every one back before closing its scope.  With the shared counter
// another thread's scope end wipes the scratch in between (a lookup fails, or
// the map is corrupted); with a per-thread depth every lookup succeeds and
// returns exactly what was stored.  After the scopes, nothing transient is left.
//
// The race needs a narrow timing window, so the second part checks its cause
// deterministically: while thread 0 holds a scope open, thread 1 -- in no scope --
// memoizes an energy.  That entry is permanent and must survive the end of
// thread 0's scope.  With the shared counter it went to thread 1's scratch (the
// counter was > 0) and thread 0's scope end deleted it.

#include "CoulFunc.h"
#include "CoulFuncCache.h"

#include <cstdio>
#include <atomic>
#ifdef _OPENMP
#include <omp.h>
#endif

int main() {
  InitializeCoulFuncCache();
  CoulFuncCache &cache = *g_coulFuncCache;
  std::atomic<long> missing{0}, wrong{0}, checked{0};
  const int rounds = 400, perRound = 200;
#ifdef _OPENMP
  const int threads = omp_get_max_threads() < 4 ? 4 : omp_get_max_threads();
#pragma omp parallel num_threads(threads)
#endif
  {
#ifdef _OPENMP
    const int t = omp_get_thread_num();
#else
    const int t = 0;
#endif
    CoulFuncCache::CoulFuncKey key;
    key.z1 = 1;
    key.z2 = 6 + t;  // one key per thread
    key.redmass = 0.93;
    key.l = 0;
    key.radius = 4.0;
    key.hybridTag = 0;
    for (int r = 0; r < rounds; r++) {
      CoulFuncCache::TransientScope scope;
      for (int k = 0; k < perRound; k++) {
        const double e = 0.1 + 1e-3 * k + 1e-7 * r;
        CoulWaves w = {e, 2 * e, 3 * e, 4 * e};
        cache.AddCoulWaves(key, e, w);
      }
      for (int k = 0; k < perRound; k++) {
        const double e = 0.1 + 1e-3 * k + 1e-7 * r;
        CoulWaves w;
        if (!cache.TryGetCoulWaves(key, e, w)) missing++;
        else if (w.F != e || w.dF != 2 * e || w.G != 3 * e || w.dG != 4 * e) wrong++;
        checked++;
      }
    }
  }
  // Thread 1 outside any scope while thread 0 holds one open.
  bool permanentKept = false;
#ifdef _OPENMP
#pragma omp parallel num_threads(2)
  {
    CoulFuncCache::CoulFuncKey k1;
    k1.z1 = 2;
    k1.z2 = 8;
    k1.redmass = 1.6;
    k1.l = 1;
    k1.radius = 5.0;
    k1.hybridTag = 0;
    const double e1 = 1.2345;
    if (omp_get_thread_num() == 0) {
      {
        CoulFuncCache::TransientScope scope;
#pragma omp barrier
#pragma omp barrier
      }
#pragma omp barrier
    } else {
#pragma omp barrier
      CoulWaves w1 = {1.0, 2.0, 3.0, 4.0};
      cache.AddCoulWaves(k1, e1, w1);
#pragma omp barrier
#pragma omp barrier
      CoulWaves r1;
      permanentKept = cache.TryGetCoulWaves(k1, e1, r1) && r1.F == 1.0 && r1.dG == 4.0;
    }
  }
#else
  permanentKept = true;  // one thread: nothing to check
#endif
  std::printf("entry memoized outside any scope while another thread held one: %s\n",
              permanentKept ? "kept" : "LOST at the end of the other thread's scope");

  // Outside every scope: the transient entries are gone.
  CoulFuncCache::CoulFuncKey key;
  key.z1 = 1;
  key.z2 = 6;
  key.redmass = 0.93;
  key.l = 0;
  key.radius = 4.0;
  key.hybridTag = 0;
  CoulWaves w;
  const bool leftover = cache.TryGetCoulWaves(key, 0.1, w);
  std::printf("%ld lookups inside transient scopes: %ld missing, %ld wrong; transient entry left after the scopes: %s\n",
              checked.load(), missing.load(), wrong.load(), leftover ? "yes" : "no");
  CleanupCoulFuncCache();
  const bool ok = missing == 0 && wrong == 0 && !leftover && permanentKept;
  std::printf("%s\n", ok ? "PASS" : "FAIL");
  return ok ? 0 : 1;
}
