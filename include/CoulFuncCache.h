#ifndef COULFUNCCACHE_H
#define COULFUNCCACHE_H

#include <map>
#include <vector>
#include <math.h>
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
        double minEnergy = 0.0;             // Minimum energy in grid
        double maxEnergy = 0.0;             // Maximum energy in grid

        // Measured usefulness of this key's memo, and whether it has been given
        // up on.  See AddCoulWaves() for why an exact-energy memo has to be able
        // to switch itself off.
        long queries = 0;
        long hits = 0;
        bool disabled = false;
    };

private:
    // Per-thread caches.  The cache is a lazily-filled performance
    // approximation, so giving each OpenMP thread its own map removes all
    // locking from the hot path (the previous single mutex-protected map
    // serialized every Coulomb lookup across threads).
    mutable std::vector<std::map<CoulFuncKey, CoulFuncData>> threadCaches_;
    std::map<CoulFuncKey, CoulFuncData>& localCache() const;

public:
    CoulFuncCache();

    /*!
     * Add pre-computed Coulomb functions for a specific parameter set
     * Called during cache initialization to populate the grid
     */
    void AddCoulWaves(const CoulFuncKey& key, double energy, const CoulWaves& waves);

    /*!
     * Look up cached Coulomb functions at the specified energy in a single pass.
     * Returns true and fills \p out only if this exact energy was memoized
     * earlier; nearby energies are deliberately not interpolated (see the
     * implementation).
     */
    bool TryGetCoulWaves(const CoulFuncKey& key, double energy, CoulWaves& out) const;

    /*!
     * Clear all cached data
     */
    void Clear();

    /*!
     * Aggregate hit statistics over every thread's cache.
     *
     * The cache is the main reason a fit is cheaper than a cold calculation, and
     * it is also the thing that quietly switches itself off when a free energy
     * shift makes the memo useless.  Both are invisible from outside unless the
     * counters are exposed, so they are.
     */
    struct Stats {
        long queries = 0;      ///< lookups attempted
        long hits = 0;         ///< lookups that found an exact stored energy
        long entries = 0;      ///< memoized (key, energy) pairs currently held
        long keys = 0;         ///< distinct (Z1,Z2,mu,l,a) keys
        long disabledKeys = 0; ///< keys that gave up and released their memo
        int threads = 0;       ///< per-thread caches in existence
    };
    Stats GetStats() const;
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
