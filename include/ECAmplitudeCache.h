#ifndef ECAMPLITUDECACHE_H
#define ECAMPLITUDECACHE_H

#include <map>
#include <vector>
#include "Constants.h"

/*!
 * Simple shared cache for EC amplitudes and their energies across all EPoints
 * This allows interpolation between different energy points when energy shifts are active
 */
class ECAmplitudeCache {
public:
    /*!
     * Key structure to identify a specific EC amplitude calculation
     */
    struct AmplitudeKey {
        int kGroupNum;
        int ecMGroupNum;
        int entranceKey;
        int exitKey;
        
        bool operator<(const AmplitudeKey& other) const {
            if (kGroupNum != other.kGroupNum) return kGroupNum < other.kGroupNum;
            if (ecMGroupNum != other.ecMGroupNum) return ecMGroupNum < other.ecMGroupNum;
            if (entranceKey != other.entranceKey) return entranceKey < other.entranceKey;
            return exitKey < other.exitKey;
        }
    };
    
    /*!
     * Data structure to hold energy-amplitude pairs for a specific key
     */
    struct AmplitudeData {
        std::vector<double> energies;
        std::vector<complex> amplitudes;
    };
    
private:
    std::map<AmplitudeKey, AmplitudeData> cache_;
    
public:
    /*!
     * Add an amplitude-energy pair to the cache
     */
    void AddAmplitude(const AmplitudeKey& key, double energy, complex amplitude);
    
    /*!
     * Get interpolated amplitude at the specified energy
     * Returns complex(0,0) if no data available for the key
     * If calculateIfMissing is true, will attempt to calculate and cache missing values
     */
    complex GetInterpolatedAmplitude(const AmplitudeKey& key, double energy, bool calculateIfMissing = false, 
                                   class CNuc* theCNuc = nullptr, const class Config* configure = nullptr) const;
    
    /*!
     * Check if we have data for a specific key
     */
    bool HasData(const AmplitudeKey& key) const;
    
    /*!
     * Clear all cached data
     */
    void Clear();
    
    /*!
     * Get cache statistics
     */
    void PrintStats() const;
    
    /*!
     * Calculate and add a missing amplitude value to cache
     * This method computes the full EC amplitude including integrals
     */
    complex CalculateAndAddAmplitude(const AmplitudeKey& key, double energy, 
                                   class CNuc* theCNuc, const class Config& configure) const;
    
private:
    /*!
     * Perform linear interpolation between two points
     */
    complex Interpolate(double energy, double e1, double e2, complex amp1, complex amp2) const;
};

// Global shared cache instance
extern ECAmplitudeCache* g_ecAmplitudeCache;

/*!
 * Initialize global cache instance
 */
void InitializeECAmplitudeCache();

/*!
 * Cleanup global cache instance
 */
void CleanupECAmplitudeCache();

#endif // ECAMPLITUDECACHE_H