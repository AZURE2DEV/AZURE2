#ifndef ECAMPLITUDECACHE_H
#define ECAMPLITUDECACHE_H

#include <map>
#include <vector>
#include <mutex>
#include "Constants.h"

/*!
 * The external-capture amplitudes of a run at the energies of each segment's points, for interpolation at a shifted
 * energy with --no-ec-table (Site::GetECAmplitudeWithShift).  One per run, held by its scratch (azure2::RunScratch::ecAmplitudes):
 * its keys (k-group, EC m-group, entrance, exit, segment key) name positions in one model, so two runs never share one.
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
  std::mutex cache_mutex_;  // the run's threads add and read under it

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
  complex GetInterpolatedAmplitude(const AmplitudeKey &key, double energy);

  /*!
   * Clear all cached data
   */
  void Clear();

 private:
  /*!
   * Perform linear interpolation between two points
   */
  complex Interpolate(double energy, double e1, double e2, complex amp1, complex amp2) const;
};

class Config;

/// Starts the run's cache afresh (azure2::RunScratch::ecAmplitudes): every amplitude of an earlier set-up of the run
/// is gone.  Set-up only (not while an evaluation runs).
void InitializeECAmplitudeCache(const Config &configure);

/// Releases the run's cache: Site::GetECAmplitudeWithShift then returns the point's own amplitudes.
void CleanupECAmplitudeCache(const Config &configure);

/// The run's cache, null when the run has none.
ECAmplitudeCache *RunECAmplitudes(const Config &configure);

#endif  // ECAMPLITUDECACHE_H