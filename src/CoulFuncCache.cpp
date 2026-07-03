#include "CoulFuncCache.h"
#include <algorithm>
#include <iostream>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

// Global cache instance
CoulFuncCache* g_coulFuncCache = nullptr;

// Maximum distance (MeV) between an energy grid point and the requested energy
// for the cached value to be considered usable for interpolation.
static const double kCacheMaxDistance = 0.00001;

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
 * Add pre-computed Coulomb functions for a specific parameter set
 * Inserts in sorted order to maintain energy grid ordering
 */
void CoulFuncCache::AddCoulWaves(const CoulFuncKey& key, double energy, const CoulWaves& waves) {
    // Get or create the data structure for this key
    CoulFuncData& data = localCache()[key];

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
 * Finalize the cache after all data has been added
 * This sorts the energy grids and computes min/max for efficient lookups
 */
void CoulFuncCache::Finalize() {
    for (auto& entry : localCache()) {
        CoulFuncData& data = entry.second;

        // Create paired vector for sorting
        std::vector<std::pair<double, CoulWaves>> paired;
        for (size_t i = 0; i < data.energies.size(); ++i) {
            paired.push_back(std::make_pair(data.energies[i], data.coulwaves[i]));
        }

        // Sort by energy
        std::sort(paired.begin(), paired.end(),
                  [](const std::pair<double, CoulWaves>& a, const std::pair<double, CoulWaves>& b) {
                      return a.first < b.first;
                  });

        // Unpack sorted data
        data.energies.clear();
        data.coulwaves.clear();
        for (const auto& p : paired) {
            data.energies.push_back(p.first);
            data.coulwaves.push_back(p.second);
        }

        // Set min/max for range checking
        if (!data.energies.empty()) {
            data.minEnergy = data.energies.front();
            data.maxEnergy = data.energies.back();
        } else {
            data.minEnergy = 0.0;
            data.maxEnergy = 0.0;
        }
    }
}

/*!
 * Get interpolated Coulomb functions at the specified energy
 */
CoulWaves CoulFuncCache::GetInterpolatedCoulWaves(const CoulFuncKey& key, double energy) const {
    // Default return value
    CoulWaves result = {0.0, 0.0, 0.0, 0.0};

    // Check if we have data for this key
    const std::map<CoulFuncKey, CoulFuncData>& cache = localCache();
    auto it = cache.find(key);
    if (it == cache.end()) {
        return result;
    }

    const CoulFuncData& data = it->second;

    // Check if energy is in range
    if (data.energies.empty() || energy < data.minEnergy || energy > data.maxEnergy) {
        return result;
    }

    // Find the two surrounding energy points using binary search
    // lower_bound returns iterator to first element >= energy
    auto upper_it = std::lower_bound(data.energies.begin(), data.energies.end(), energy);

    // Handle exact match
    if (upper_it != data.energies.end() && fabs(*upper_it - energy) < 1e-12) {
        size_t idx = std::distance(data.energies.begin(), upper_it);
        return data.coulwaves[idx];
    }

    // Handle edge cases
    if (upper_it == data.energies.begin()) {
        return data.coulwaves[0];
    }
    if (upper_it == data.energies.end()) {
        return data.coulwaves.back();
    }

    // Get surrounding points for interpolation
    auto lower_it = upper_it - 1;
    size_t lower_idx = std::distance(data.energies.begin(), lower_it);
    size_t upper_idx = lower_idx + 1;

    double e1 = data.energies[lower_idx];
    double e2 = data.energies[upper_idx];
    const CoulWaves& waves1 = data.coulwaves[lower_idx];
    const CoulWaves& waves2 = data.coulwaves[upper_idx];

    // Perform interpolation
    return InterpolateCoulWaves(energy, e1, e2, waves1, waves2);
}

/*!
 * Perform linear interpolation of CoulWaves structures
 */
CoulWaves CoulFuncCache::InterpolateCoulWaves(double energy, double e1, double e2,
                                               const CoulWaves& waves1, const CoulWaves& waves2) const {
    CoulWaves result;

    // Avoid division by zero
    if (fabs(e2 - e1) < 1e-12) {
        return waves1;
    }

    // Linear interpolation factor
    double t = (energy - e1) / (e2 - e1);

    // Interpolate each component
    result.F = waves1.F + t * (waves2.F - waves1.F);
    result.dF = waves1.dF + t * (waves2.dF - waves1.dF);
    result.G = waves1.G + t * (waves2.G - waves1.G);
    result.dG = waves1.dG + t * (waves2.dG - waves1.dG);

    return result;
}

/*!
 * Single-pass lookup: combines the old IsInRange (acceptance test) and
 * GetInterpolatedCoulWaves (value computation) so the cache map is traversed
 * only once per query.  Preserves the original acceptance rules exactly.
 */
bool CoulFuncCache::TryGetCoulWaves(const CoulFuncKey& key, double energy, CoulWaves& out) const {
    const std::map<CoulFuncKey, CoulFuncData>& cache = localCache();
    auto it = cache.find(key);
    if (it == cache.end()) return false;

    const CoulFuncData& data = it->second;
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
        return true;
    }
    return false;
}

/*!
 * Check if we have cached data for a specific parameter set
 */
bool CoulFuncCache::HasData(const CoulFuncKey& key) const {
    const std::map<CoulFuncKey, CoulFuncData>& cache = localCache();
    return cache.find(key) != cache.end();
}

/*!
 * Check if a specific energy is within the cached energy range
 * and that nearby grid points are close enough for reliable interpolation
 */
bool CoulFuncCache::IsInRange(const CoulFuncKey& key, double energy) const {
    CoulWaves dummy;
    return TryGetCoulWaves(key, energy, dummy);
}

/*!
 * Clear all cached data (all per-thread caches)
 */
void CoulFuncCache::Clear() {
    for (auto& c : threadCaches_) c.clear();
}

/*!
 * Print cache statistics for debugging
 */
void CoulFuncCache::PrintStats() const {
    const std::map<CoulFuncKey, CoulFuncData>& cache = localCache();
    std::cout << "=== CoulFunc Cache Statistics ===" << std::endl;
    std::cout << "Number of cached parameter sets (this thread): " << cache.size() << std::endl;

    for (const auto& entry : cache) {
        const CoulFuncKey& key = entry.first;
        const CoulFuncData& data = entry.second;

        std::cout << "  Key: z1=" << key.z1 << ", z2=" << key.z2
                  << ", redmass=" << key.redmass << ", l=" << key.l
                  << ", radius=" << key.radius << std::endl;
        std::cout << "    Energy points: " << data.energies.size()
                  << " (range: " << data.minEnergy << " - " << data.maxEnergy << " MeV)" << std::endl;
    }

    std::cout << "=================================" << std::endl;
}
