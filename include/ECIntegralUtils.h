#ifndef ECINTEGRALUTILS_H
#define ECINTEGRALUTILS_H

#include "Config.h"

class CNuc;
class EData;

///Utility functions for external capture integral caching and management

/*!
 * The ECIntegralUtils namespace provides utility functions for managing
 * external capture integral caching, including pre-computation and 
 * cache management functions.
 */

namespace ECIntegralUtils {
    /*!
     * Pre-compute external capture integrals for all relevant parameter combinations
     * found in the nuclear data and experimental data
     */
    void PrecomputeAllIntegrals(CNuc* compound, EData* data, const Config& configure);
    
    /*!
     * Pre-compute integrals for a specific energy range
     */
    void PrecomputeIntegralsForRange(CNuc* compound, double minEnergy, double maxEnergy, 
                                   double deltaEnergy, const Config& configure);
    
    /*!
     * Initialize the cache system with appropriate file paths
     */
    void InitializeCacheSystem(const Config& configure);
    
    /*!
     * Finalize and save the cache system
     */
    void FinalizeCacheSystem();
    
    /*!
     * Analyze nuclear data to determine appropriate energy ranges for caching
     */
    void AnalyzeEnergyRanges(CNuc* compound, EData* data, 
                           double& minEnergy, double& maxEnergy, double& suggestedDelta);
    
    /*!
     * Print cache utilization statistics
     */
    void PrintCacheStatistics(const Config& configure);
}

#endif // ECINTEGRALUTILS_H