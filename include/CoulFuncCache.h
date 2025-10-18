#ifndef COULFUNCCACHE_H
#define COULFUNCCACHE_H

#include <map>
#include <vector>
#include <mutex>
#include "CoulFunc.h"

/*!
 * Shared cache for Coulomb function calculations across all CoulFunc instances
 * This allows pre-computation on a fine energy grid and interpolation during calculations
 * to significantly improve performance when CoulFunc instances are created many times
 */
class CoulFuncCache {
public:
    /*!
     * Key structure to identify a specific Coulomb function parameter set
     * Based on the physical parameters that determine the Coulomb functions
     */
    struct CoulFuncKey {
        int z1;           // Atomic number of first particle
        int z2;           // Atomic number of second particle
        double redmass;   // Reduced mass of the particle pair
        int l;            // Orbital angular momentum
        double radius;    // Channel radius

        // Comparison operator for use as map key
        bool operator<(const CoulFuncKey& other) const {
            if (z1 != other.z1) return z1 < other.z1;
            if (z2 != other.z2) return z2 < other.z2;
            // Use epsilon comparison for floating point
            if (fabs(redmass - other.redmass) > 1e-10) return redmass < other.redmass;
            if (l != other.l) return l < other.l;
            if (fabs(radius - other.radius) > 1e-10) return radius < other.radius;
            return false;
        }
    };

    /*!
     * Data structure to hold energy-CoulWaves pairs for a specific key
     * Stores pre-computed Coulomb functions on an energy grid
     */
    struct CoulFuncData {
        std::vector<double> energies;       // Energy grid points
        std::vector<CoulWaves> coulwaves;   // Coulomb functions at each energy
        double minEnergy;                   // Minimum energy in grid
        double maxEnergy;                   // Maximum energy in grid
    };

private:
    std::map<CoulFuncKey, CoulFuncData> cache_;
    mutable std::mutex cache_mutex_;  // Protect cache from concurrent access

public:
    /*!
     * Add pre-computed Coulomb functions for a specific parameter set
     * Called during cache initialization to populate the grid
     */
    void AddCoulWaves(const CoulFuncKey& key, double energy, const CoulWaves& waves);

    /*!
     * Get interpolated Coulomb functions at the specified energy
     * Uses linear interpolation between grid points
     * Returns empty CoulWaves {0,0,0,0} if no data available or out of range
     */
    CoulWaves GetInterpolatedCoulWaves(const CoulFuncKey& key, double energy) const;

    /*!
     * Check if we have cached data for a specific parameter set
     */
    bool HasData(const CoulFuncKey& key) const;

    /*!
     * Check if a specific energy is within the cached energy range
     */
    bool IsInRange(const CoulFuncKey& key, double energy) const;

    /*!
     * Clear all cached data
     */
    void Clear();

    /*!
     * Finalize cache data structure after all energies have been added
     * This sorts the energy grid and computes min/max for efficient lookups
     */
    void Finalize();

    /*!
     * Get cache statistics for debugging
     */
    void PrintStats() const;

private:
    /*!
     * Perform linear interpolation of CoulWaves structures
     * Interpolates each component (F, dF, G, dG) separately
     */
    CoulWaves InterpolateCoulWaves(double energy, double e1, double e2,
                                    const CoulWaves& waves1, const CoulWaves& waves2) const;
};

// Global shared cache instance
extern CoulFuncCache* g_coulFuncCache;

/*!
 * Initialize global cache instance
 */
void InitializeCoulFuncCache();

/*!
 * Cleanup global cache instance
 */
void CleanupCoulFuncCache();

#endif // COULFUNCCACHE_H
