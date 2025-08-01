#ifndef ECINTEGRALCACHE_H
#define ECINTEGRALCACHE_H

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include "Constants.h"

class PPair;
class Config;

///Enhanced caching system for external capture integrals with interpolation

/*!
 * The ECIntegralCache class provides an efficient caching and interpolation
 * system for external capture integrals. It pre-computes integrals on a 
 * fine energy grid and uses interpolation for intermediate energies.
 */

class ECIntegralCache {
public:
    /*!
     * Structure to hold integral calculation parameters as a unique key
     */
    struct IntegralKey {
        int liValue;          // Initial orbital angular momentum
        int lfValue;          // Final orbital angular momentum  
        double siValue;       // Initial channel spin
        double sfValue;       // Final channel spin
        double jInitial;      // Initial total angular momentum
        double jFinal;        // Final total angular momentum
        int multL;            // Gamma ray multipolarity
        char radType;         // Radiation type
        double bindingEnergy; // Final state binding energy
        bool isChannelCapture; // Channel capture flag
        double separationEnergy; // PPair separation energy (GetSepE + GetExE)
        
        // Comparison operators for use as map key
        bool operator<(const IntegralKey& other) const;
        bool operator==(const IntegralKey& other) const;
    };
    
    /*!
     * Structure to hold cached integral data for energy interpolation
     */
    struct CachedIntegrals {
        std::vector<double> energies;      // Energy grid points
        std::vector<complex> integrals;    // Corresponding integral values
        double minEnergy;                  // Minimum energy in grid
        double maxEnergy;                  // Maximum energy in grid
        double deltaEnergy;                // Energy step size
    };

private:
    std::map<IntegralKey, CachedIntegrals> cache_;  // Main cache storage
    std::string cacheFilename_;                     // Cache file name
    bool cacheInitialized_;                         // Cache initialization flag
    double energyTolerance_;                        // Energy matching tolerance
    
public:
    /*!
     * Constructor initializes the cache system
     */
    ECIntegralCache(const std::string& cacheFile = "");
    
    /*!
     * Destructor saves cache if needed
     */
    ~ECIntegralCache();
    
    /*!
     * Get integral value with interpolation if needed
     * If forceAdd is true, calculate and add missing values to cache
     */
    complex GetIntegral(const IntegralKey& key, double energy, PPair* pair, const Config& configure, bool forceAdd = false);
    
    /*!
     * Pre-compute integrals for a given key over energy range
     */
    void PrecomputeIntegrals(const IntegralKey& key, double minEnergy, double maxEnergy, 
                           double deltaEnergy, PPair* pair, const Config& configure);
    
    /*!
     * Load cache from file
     */
    bool LoadCache(const std::string& filename = "");
    
    /*!
     * Save cache to file
     */
    bool SaveCache(const std::string& filename = "") const;
    
    /*!
     * Clear all cached data
     */
    void ClearCache();
    
    /*!
     * Check if integral is cached for given key and energy
     */
    bool IsCached(const IntegralKey& key, double energy) const;
    
    /*!
     * Check if any data exists for given key (regardless of energy)
     */
    bool HasKey(const IntegralKey& key) const;
    
    /*!
     * Set energy tolerance for cache lookup
     */
    void SetEnergyTolerance(double tolerance) { energyTolerance_ = tolerance; }
    
    /*!
     * Get cache statistics
     */
    void PrintCacheStats(std::ostream& out) const;
    
    /*!
     * Pre-populate cache for all relevant integral keys with reasonable energy range
     */
    void PopulateCacheForData(class EData* data, class CNuc* compound, const Config& configure);
    
    /*!
     * Add a single integral point to the cache dynamically
     */
    void AddIntegralPoint(const IntegralKey& key, double energy, complex integral);

private:
    /*!
     * Interpolate integral value at given energy
     */
    complex InterpolateIntegral(const CachedIntegrals& cached, double energy) const;
    
    /*!
     * Find closest energy indices for interpolation
     */
    std::pair<int, int> FindEnergyIndices(const std::vector<double>& energies, double energy) const;
    
    /*!
     * Generate cache key from parameters
     */
    std::string GenerateCacheKey(const IntegralKey& key) const;
    
    /*!
     * Parse cache key from string
     */
    bool ParseCacheKey(const std::string& keyStr, IntegralKey& key) const;
};

// Global cache instance
extern ECIntegralCache* g_ecIntegralCache;

/*!
 * Initialize global cache instance
 */
void InitializeECIntegralCache(const std::string& cacheFile = "");

/*!
 * Cleanup global cache instance
 */
void CleanupECIntegralCache();

#endif // ECINTEGRALCACHE_H