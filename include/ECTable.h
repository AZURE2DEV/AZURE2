#ifndef ECTABLE_H
#define ECTABLE_H

#include "Constants.h"
#include <string>

class PPair;
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
 * for, and kept for the life of the process; any number of threads may ask at
 * once (a node being computed by one thread is waited for by the others).
 *
 * A term is keyed by its physical content (charges, masses, spins, channel
 * radius, threshold, quantum numbers, final level energy, the options that
 * change the integrand), never by a position in a model, so one process can
 * hold several models or sessions and an entry can never be stale.
 *
 * Sub-threshold entrance channels (E_out <= 0, the Whittaker-Whittaker branch)
 * are computed directly.  --no-ec-table (Config::useECTable = false) turns the
 * table off and restores the direct integrals and the interpolated shifts.
 * AZURE2_EC_TABLE_CHECK=1 computes every tabulated integral directly as well
 * and records the largest relative difference (Report).
 */
class ECTable {
 public:
  /// The integral ECIntegral(finalPair, configure)(...) returns, from the table.
  static complex Integral(PPair *finalPair, const Config &configure,
                          int li, int lf, double si, double sf, double ji, double jf,
                          int multL, char radType, double inEnergy, double levelEnergy,
                          bool isChannelCapture);
  /// Node values (exact integrals) computed so far, over all terms.
  static long long NodesComputed();
  /// Integrals answered from the table, and directly (sub-threshold or failed nodes).
  static long long Lookups();
  static long long DirectCalls();
  /// Largest relative table - direct difference seen with AZURE2_EC_TABLE_CHECK=1 (else -1).
  static double MaxCheckDeviation();
  /// One line for the run log: terms, nodes, lookups since the last ResetCounters.
  static std::string Report();
  /// Restart the counters of Report (the nodes stay).
  static void ResetCounters();
  /// Forget every node (tests).
  static void Clear();
};

#endif
