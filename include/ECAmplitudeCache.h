#ifndef ECAMPLITUDECACHE_H
#define ECAMPLITUDECACHE_H

#include <map>
#include <vector>
#include <mutex>
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
    int segmentKey;  // Add segment identifier to prevent cross-segment contamination

    bool operator<(const AmplitudeKey &other) const {
      if (kGroupNum != other.kGroupNum) return kGroupNum < other.kGroupNum;
      if (ecMGroupNum != other.ecMGroupNum) return ecMGroupNum < other.ecMGroupNum;
      if (entranceKey != other.entranceKey) return entranceKey < other.entranceKey;
      if (exitKey != other.exitKey) return exitKey < other.exitKey;
      return segmentKey < other.segmentKey;
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
  mutable std::mutex cache_mutex_;  // Protect cache from concurrent access

 public:
  /*!
   * Add an amplitude-energy pair to the cache
   */
  void AddAmplitude(const AmplitudeKey &key, double energy, complex amplitude);

  /*!
   * Get interpolated amplitude at the specified energy
   * Returns complex(0,0) if no data available for the key
   * Simple interpolation only - no calculation
   */
  complex GetInterpolatedAmplitude(const AmplitudeKey &key, double energy) const;

  /*!
   * Check if we have data for a specific key
   */
  bool HasData(const AmplitudeKey &key) const;

  /*!
   * Clear all cached data
   */
  void Clear();

  /*!
   * Get cache statistics
   */
  void PrintStats() const;

 private:
  /*!
   * Perform linear interpolation between two points
   */
  complex Interpolate(double energy, double e1, double e2, complex amp1, complex amp2) const;
};

// Global shared cache instance
extern ECAmplitudeCache *g_ecAmplitudeCache;

/*!
 * Initialize global cache instance
 */
void InitializeECAmplitudeCache();

/*!
 * Cleanup global cache instance
 */
void CleanupECAmplitudeCache();

#endif  // ECAMPLITUDECACHE_H