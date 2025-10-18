#include "CoulFuncCache.h"
#include <algorithm>
#include <iostream>
#include <cmath>

// Global cache instance
CoulFuncCache* g_coulFuncCache = nullptr;

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
    std::lock_guard<std::mutex> lock(cache_mutex_);

    // Get or create the data structure for this key
    CoulFuncData& data = cache_[key];

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
    std::lock_guard<std::mutex> lock(cache_mutex_);

    for (auto& entry : cache_) {
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
    std::lock_guard<std::mutex> lock(cache_mutex_);

    // Default return value
    CoulWaves result = {0.0, 0.0, 0.0, 0.0};

    // Check if we have data for this key
    auto it = cache_.find(key);
    if (it == cache_.end()) {
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
 * Check if we have cached data for a specific parameter set
 */
bool CoulFuncCache::HasData(const CoulFuncKey& key) const {
    std::lock_guard<std::mutex> lock(cache_mutex_);
    return cache_.find(key) != cache_.end();
}

/*!
 * Check if a specific energy is within the cached energy range
 * and that nearby grid points are close enough (within 0.001 MeV) for reliable interpolation
 */
bool CoulFuncCache::IsInRange(const CoulFuncKey& key, double energy) const {
    std::lock_guard<std::mutex> lock(cache_mutex_);

    auto it = cache_.find(key);
    if (it == cache_.end()) {
        return false;
    }

    const CoulFuncData& data = it->second;

    // Check if energy is within overall range
    if (energy < data.minEnergy || energy > data.maxEnergy) {
        return false;
    }

    // Check if we have data points close enough for interpolation
    // Find the surrounding energy points
    auto upper_it = std::lower_bound(data.energies.begin(), data.energies.end(), energy);

    // Check for exact match
    if (upper_it != data.energies.end() && fabs(*upper_it - energy) < 1e-12) {
        return true;
    }

    // Check distances to nearest grid points
    const double maxDistance = 0.001; // Maximum distance in MeV for interpolation

    if (upper_it == data.energies.begin()) {
        // Return false
        return false;
    }

    if (upper_it == data.energies.end()) {
        // Energy is after last point, check distance to last point
        return false;
    }

    // Energy is between two points, check distances to both
    auto lower_it = upper_it - 1;
    double distToLower = energy - *lower_it;
    double distToUpper = *upper_it - energy;

    return (distToLower <= maxDistance && distToUpper <= maxDistance);
}

/*!
 * Clear all cached data
 */
void CoulFuncCache::Clear() {
    std::lock_guard<std::mutex> lock(cache_mutex_);
    cache_.clear();
}

/*!
 * Print cache statistics for debugging
 */
void CoulFuncCache::PrintStats() const {
    std::lock_guard<std::mutex> lock(cache_mutex_);

    std::cout << "=== CoulFunc Cache Statistics ===" << std::endl;
    std::cout << "Number of cached parameter sets: " << cache_.size() << std::endl;

    for (const auto& entry : cache_) {
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
