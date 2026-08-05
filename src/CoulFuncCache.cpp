#include "CoulFuncCache.h"
#include <algorithm>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

// Global cache instance
CoulFuncCache* g_coulFuncCache = nullptr;

/* This is an exact-energy memo, so it only pays off when the same energies are
 * queried again: fixed data points evaluated over and over during a fit, or a
 * reaction-rate sweep revisiting grid points.  It is pure overhead when the
 * energies never recur -- which is precisely what a *varying* energy shift does,
 * since it moves every point energy on every likelihood evaluation.
 *
 * Left unbounded that case degrades badly: the grid grew by one entry per point
 * per evaluation, every insertion memmoved the whole vector to stay sorted, and
 * nothing was ever evicted.  An MCMC run with varying energy shifts therefore
 * slowed down quadratically (measured on a 4-segment 3H+d fit: 333 ms/step over
 * 100 steps, 737 ms/step over 200) while leaking memory the whole way.
 *
 * Two guards bound it.  The hard cap stops the vectors growing without limit, so
 * insertion and lookup cost stay constant no matter how long the run goes.  The
 * adaptive shutoff then notices a key whose measured hit rate never justifies the
 * memo and hands its memory back, so the varying-shift case settles into plain
 * recomputation rather than paying to maintain a memo that never hits.
 *
 * Neither guard changes any computed value: a miss just recomputes the Coulomb
 * functions exactly as an uncached build would.
 */
static const std::size_t kMaxEntriesPerKey = 32768;
static const long kMinQueriesBeforeJudging = 4096;
static const double kMinUsefulHitRate = 0.05;

CoulFuncCache::CoulFuncCache() {
#ifdef _OPENMP
    int n = omp_get_max_threads();
#else
    int n = 1;
#endif
    if (n < 1) n = 1;
    threadCaches_.resize(n);
}

/*!
 * Returns the calling thread's private cache map (no locking required).
 */
std::map<CoulFuncCache::CoulFuncKey, CoulFuncCache::CoulFuncData>& CoulFuncCache::localCache() const {
#ifdef _OPENMP
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif
    if (tid < 0 || tid >= (int)threadCaches_.size()) tid = 0;
    return threadCaches_[tid];
}

/*!
 * Initialize the global Coulomb function cache
 */
void InitializeCoulFuncCache() {
    if (!g_coulFuncCache) {
        g_coulFuncCache = new CoulFuncCache();
    }
}

/*!
 * Cleanup the global Coulomb function cache
 */
void CleanupCoulFuncCache() {
    if (g_coulFuncCache) {
        delete g_coulFuncCache;
        g_coulFuncCache = nullptr;
    }
}

/*!
 * Memoize the Coulomb functions computed at one energy, keeping the energy grid
 * sorted.  Subject to the two guards described at the top of this file.
 */
void CoulFuncCache::AddCoulWaves(const CoulFuncKey& key, double energy, const CoulWaves& waves) {
    // Get or create the data structure for this key
    CoulFuncData& data = localCache()[key];

    if (data.disabled) return;

    // Give up on a key whose energies do not recur, and release what it holds.
    // Judged only after enough queries to be past the initial fill, where misses
    // are expected.
    if (data.queries >= kMinQueriesBeforeJudging &&
        double(data.hits) < kMinUsefulHitRate * double(data.queries)) {
        data.disabled = true;
        std::vector<double>().swap(data.energies);
        std::vector<CoulWaves>().swap(data.coulwaves);
        data.minEnergy = 0.0;
        data.maxEnergy = 0.0;
        return;
    }

    // Full: serve what is already memoized, stop taking on more. Keeps insertion
    // and lookup cost constant for the rest of the run.
    if (data.energies.size() >= kMaxEntriesPerKey) return;

    // Find insertion point to keep energies sorted
    auto it = std::lower_bound(data.energies.begin(), data.energies.end(), energy);

    // Check if this exact energy already exists (avoid duplicates)
    if (it != data.energies.end() && fabs(*it - energy) < 1e-12) {
        // Energy already exists, update the value
        size_t idx = std::distance(data.energies.begin(), it);
        data.coulwaves[idx] = waves;
        return;
    }

    // Insert at the correct position to maintain sorted order
    size_t idx = std::distance(data.energies.begin(), it);
    data.energies.insert(it, energy);
    data.coulwaves.insert(data.coulwaves.begin() + idx, waves);

    // Update min/max bounds
    if (data.energies.size() == 1) {
        data.minEnergy = energy;
        data.maxEnergy = energy;
    } else {
        data.minEnergy = data.energies.front();
        data.maxEnergy = data.energies.back();
    }
}

/*!
 * Single-pass exact-energy lookup.  Also records the hit/miss statistics that
 * AddCoulWaves() uses to decide whether this key's memo is worth maintaining.
 */
bool CoulFuncCache::TryGetCoulWaves(const CoulFuncKey& key, double energy, CoulWaves& out) const {
    std::map<CoulFuncKey, CoulFuncData>& cache = localCache();
    auto it = cache.find(key);
    if (it == cache.end()) return false;

    CoulFuncData& data = it->second;
    if (data.disabled) return false;

    data.queries++;

    if (data.energies.empty() || energy < data.minEnergy || energy > data.maxEnergy) return false;

    auto upper_it = std::lower_bound(data.energies.begin(), data.energies.end(), energy);

    // Exact grid point only.  We deliberately do NOT interpolate between nearby
    // energies: for a fixed channel radius the cache is queried while sweeping
    // energy (adaptive reaction-rate / cross-section integrals), and returning an
    // interpolated value for an energy that is merely *close* to a cached one
    // makes sigma(E) inconsistent from one quadrature evaluation to the next.
    // That inconsistency defeats adaptive integrators (gsl reports "divergent")
    // and made the reaction rate depend on cache/thread history.  Exact-match
    // memoization is self-consistent and still fast for repeated evaluations.
    if (upper_it != data.energies.end() && fabs(*upper_it - energy) < 1e-12) {
        out = data.coulwaves[std::distance(data.energies.begin(), upper_it)];
        data.hits++;
        return true;
    }
    return false;
}

/*!
 * Clear all cached data (all per-thread caches)
 */
void CoulFuncCache::Clear() {
    for (auto& c : threadCaches_) c.clear();
}

/*!
 * Aggregate hit statistics over every thread's cache.
 *
 * Summed rather than per-thread: which thread served a lookup is an accident of
 * scheduling, and the quantity of interest -- what fraction of the Coulomb
 * evaluations the memo saved -- is a property of the run.  `disabledKeys`
 * counts the keys whose measured hit rate never justified the memo and which
 * therefore handed their memory back; a run with a free energy shift should
 * show most of its keys there.
 */
CoulFuncCache::Stats CoulFuncCache::GetStats() const {
    Stats s;
    s.threads = static_cast<int>(threadCaches_.size());
    for (const auto& cache : threadCaches_) {
        s.keys += static_cast<long>(cache.size());
        for (const auto& kv : cache) {
            s.queries += kv.second.queries;
            s.hits += kv.second.hits;
            s.entries += static_cast<long>(kv.second.energies.size());
            if (kv.second.disabled) ++s.disabledKeys;
        }
    }
    return s;
}
