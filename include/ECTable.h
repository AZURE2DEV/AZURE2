#ifndef ECTABLE_H
#define ECTABLE_H

#include "Constants.h"
#include <future>
#include <atomic>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

class PPair;
class ChannelTable;
class Config;

/*!
 * Energy tables of the external-capture radial integrals.
 *
 * An external-capture integral (ECIntegral) depends on the energy only through
 * the entrance energy of its final pair, so every data point, sub-point,
 * shifted point and rate energy of one capture term (final pair, l_i -> l_f,
 * multipolarity, spins, final level energy) samples one smooth function.  The
 * direct calculation integrates it from scratch at each energy: ~1 ms each, and
 * resolution or target integration multiplies the energies by 10-30.
 *
 * The table computes the integral exactly at nodes uniform in u = ln(E_out),
 * E_out the c.m. energy above the final pair's threshold, and interpolates with
 * an 8-point Lagrange polynomial.  The steep Coulomb behaviour is divided out
 * first, so the tabulated function is smooth: the integral divided by
 * sqrt(P_li(E_out)) for an ordinary term, and by the hard-sphere phase factor
 * of the final pair for channel capture (whose integral carries it).  Nodes are
 * computed on demand, the first time an energy within reach of them is asked
 * for, and kept for the life of the run (ECTermStore, in the run's scratch:
 * Config::scratch); any number of threads may ask at once (a node being
 * computed by one thread is waited for by the others).
 *
 * A term is keyed by its physical content (charges, masses, spins, channel
 * radius, threshold, quantum numbers, final level energy, the options that
 * change the integrand), never by a position in a model, so the pooled copies
 * of a fit and a session's rebuilds share it and an entry can never be stale.
 * The Coulomb table is fixed within a run (PPair::FunctionSettings), so the key
 * holds whether it is used, not its file: two runs with two tables have two
 * stores.
 *
 * Sub-threshold entrance channels (E_out <= 0, the Whittaker-Whittaker branch)
 * are computed directly.  --no-ec-table (Config::useECTable = false) turns the
 * table off and restores the direct integrals and the interpolated shifts.
 * AZURE2_EC_TABLE_CHECK=1 computes every tabulated integral directly as well
 * and records the largest relative difference (Report).
 */
/// The nodes of one capture term: (reduced) integrals at ln(E_out) = i h, each filled once by the first thread that
/// asks; the statistics of AZURE2_EC_TABLE_CHECK.
struct ECTermNodes {
  /// The tabulated integral at one node; ok = false when it could not be computed (a non-finite value or a GSL
  /// error), which sends every energy that needs it to the direct integral.
  struct Value {
    complex value;
    bool ok;
  };
  std::mutex m;
  std::map<long, std::shared_future<Value>> nodes;
  // AZURE2_EC_TABLE_CHECK: the worst relative difference and where, the largest direct value and the largest
  // absolute difference.
  double maxDev = -1.0, eAtMax = 0.0, absAtMax = 0.0, maxAbs = 0.0, maxAbsDiff = 0.0;
  double eMin = 1e300, eMax = -1e300;
};

/// Every term a run has asked for, keyed by everything the integral depends on besides the energy.  Owned by the
/// run (azure2::RunScratch::ecTerms), shared by its threads behind the lock.
struct ECTermStore {
  std::mutex mutex;
  std::map<std::vector<double>, std::unique_ptr<ECTermNodes>> terms;
  /// What ECTable::Report counts for the run: node values computed (all, and at the last ResetCounters), integrals
  /// answered from the table and directly, and the largest table - direct difference of AZURE2_EC_TABLE_CHECK (-1:
  /// none yet; behind checkMutex).
  std::atomic<long long> nodes{0}, nodesAtReset{0}, lookups{0}, direct{0};
  std::mutex checkMutex;
  double maxDeviation = -1.0;
};

class ECTable {
 public:
  /// The integral ECIntegral(finalPair, configure, radius)(...) returns, from the table; \p radius is the radius
  /// a_c (fm) of the final-state channel.
  static complex Integral(const PPair *finalPair, double radius, const Config &configure,
                          int li, int lf, double si, double sf, double ji, double jf,
                          int multL, char radType, double inEnergy, double levelEnergy,
                          bool isChannelCapture, const ChannelTable *table = nullptr, int finalPairNum = 0);
  /// Node values (exact integrals) the run computed so far, over all terms.
  static long long NodesComputed(const Config &configure);
  /// Integrals the run answered from the table, and directly (sub-threshold or failed nodes).
  static long long Lookups(const Config &configure);
  static long long DirectCalls(const Config &configure);
  /// Largest relative table - direct difference the run saw with AZURE2_EC_TABLE_CHECK=1 (else -1).
  static double MaxCheckDeviation(const Config &configure);
  /// One line for the run log: the run's terms, nodes and lookups since the last ResetCounters.
  static std::string Report(const Config &configure);
  /// Restart the run's counters of Report (the nodes stay).
  static void ResetCounters(const Config &configure);
  /// Forget every node of the run (tests).
  static void Clear(const Config &configure);
};

#endif
