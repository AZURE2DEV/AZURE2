#ifndef COULOMBTABLE_H
#define COULOMBTABLE_H

#include <atomic>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

/// Universal tables of the channel functions of charged two-body channels.

/*!
 * Everything AZURE2 needs from the Coulomb functions of a channel -- the
 * penetrability P_l, the shift function S_l and the hard-sphere phase phi_l of
 * an open channel, and the shift function of a closed one -- depends only on
 * the Sommerfeld parameter eta, rho = k a and l.  One table per l over
 * (log10 eta, log10 rho) therefore serves every pair, charge, mass, radius and
 * energy.  The design is in docs/source/developer/coulomb_table.rst (from
 * 17O+a/10-13-26_natNe_seg1/coulomb_table_spec.md).
 *
 * Tabulated, per l and grid node:
 *   open:  log P - log C0^2(eta), with C0^2 = 2 pi eta / (exp(2 pi eta) - 1);
 *          S;
 *          phi - rho, where phi = arg(G + iF) made continuous from phi(0) = 0
 *          (dphi/drho = P/rho, so phi is monotonic);
 *   bound: S + rho_b of the decaying Whittaker solution, at eta_b and rho_b = kappa a.
 * F and G themselves are never tabulated: above the turning point they
 * oscillate, while the quantities above vary slowly.
 *
 * Interpolation: a bicubic Hermite patch per grid cell, with node derivatives
 * from sixth-order finite differences of the node values (computed when the
 * table is built or loaded, never stored), or tensor-product Lagrange
 * interpolation of order 6 or 8 from the node values alone.  A lookup is a pure function of
 * (eta, rho, l) and of the table file: it never depends on what was evaluated
 * before.  That is the defect that removed interpolation from CoulFuncCache
 * (history-dependent values broke adaptive quadrature; CoulFuncCache.cpp).
 */
namespace azure2 {
namespace support {
class Environment;
}
}  // namespace azure2

/// How often a run's channel functions found their values in the table and how often outside it (the exact routine
/// then answered): the caller's counters (RunDiagnostics::coulombTable), given to each lookup.
struct CoulombTableUse {
  std::atomic<long> lookups{0}, outside{0};
};

class CoulombTable {
 public:
  /// A uniform axis in log10 of the variable.
  struct Axis {
    double x0 = 0.0;  ///< log10 of the first node
    double dx = 0.0;  ///< node spacing in log10
    int n = 0;        ///< number of nodes
    double Last() const { return x0 + dx * (n - 1); }
  };
  /// What to build: ranges in eta and rho (both kinds), l_max, nodes per decade.
  struct BuildSpec {
    int lMax = 10;
    double etaMin = 1e-3, etaMax = 100.0;
    double rhoMin = 1e-3, rhoMax = 100.0;
    double perDecade = 100.0;  // within 1e-7 everywhere (2026-10-18 sweep)
    /// Upper limit of eta * rho = 0.0344 Z1 Z2 mu[u] a[fm], a constant of a channel
    /// (p+12C ~0.8, 12C+a ~7, 12C+12C ~45, 16O+16O ~120); larger values go to the
    /// exact routine.  A channel's Coulomb turning point lies at eta ~ sqrt(eta rho / 2),
    /// so the limit also bounds the eta at which the table meets a turning point.
    double etaRhoMax = 150.0;
    int threads = 1;
    /// Optional progress and cancel hooks (GUI, CLI): rows of BuildRows() done so far,
    /// and a flag that stops the build early (Build then returns false).
    std::atomic<int> *rowsDone = nullptr;
    const std::atomic<bool> *cancel = nullptr;
  };
  /// Open-channel values from the table or from the exact routine.
  struct Open {
    double P, S, phi;
  };

  CoulombTable() = default;

  /// Interpolation order: 8 (the default) or 6 = tensor-product Lagrange through 8 or
  /// 6 nodes per axis, 4 = bicubic Hermite.  The 2026-10-18 sweep: at 60 nodes per decade
  /// order 8 is 5-30x more accurate than 4 or 6 for the same table.
  void SetOrder(int order);
  int Order() const { return order_; }

  /// Build from the exact routines; false when spec.cancel stopped it.
  bool Build(const BuildSpec &spec);
  /// The number of rows Build computes (the unit of spec.rowsDone).
  static int BuildRows(const BuildSpec &spec);
  /// The threads of the command line's build (--build-coulomb-table): every OpenMP thread (1 without OpenMP).
  static int DefaultBuildThreads();
  /// Build and write path (empty: DefaultPath()), creating its directory; the file is
  /// written as path.new and renamed over path, so a process reading the old file keeps
  /// it.  False with a reason in error on cancel or an I/O error; nothing is left behind.
  static bool BuildFile(const std::string &path, const BuildSpec &spec, std::string &error);
  /// Write the table; returns false on an I/O error.
  bool Save(const std::string &path) const;
  /// Read a table; returns false (with a reason in error) on any mismatch.
  bool Load(const std::string &path, std::string &error);

  /// Below this the tabulated phase (stored as phi - rho, so no better than about 1e-16 rad absolute, and within
  /// the interpolation error elsewhere) is replaced by SmallPhase's, which keeps the phase's relative accuracy.
  static constexpr double kSmallPhase = 1e-6;
  /// The hard-sphere phase deep below the barrier from the regular function's power series and the penetrability:
  /// sin(phi) = F sqrt(P / rho), with logP = ln P.  Relative accuracy that of P; for phi < kSmallPhase, where the
  /// series has no cancellation.  NaN when the series does not converge (the caller then uses the exact routine).
  static double SmallPhase(int l, double eta, double rho, double logP);
  /// Table lookup of an open channel; false outside the table.  Counted in `use` (the run's; null: not counted).
  /// P and S are the table's; so is phi, except below kSmallPhase (deep below the barrier), where it is SmallPhase's.
  bool LookupOpen(int l, double eta, double rho, Open &out, CoulombTableUse *use) const;
  /// Table lookup of a closed channel's shift function; false outside the table.  Counted in `use` as above.
  bool LookupBound(int l, double etaB, double rhoB, double &S, CoulombTableUse *use) const;

  /// The table at path (empty: defaultPath, the run's Config::environment.coulombTableDefault), loaded once per
  /// process and shared by all threads; nullptr, after one warning, when it cannot be read.
  static const CoulombTable *Shared(const std::string &path, const std::string &defaultPath);
  /// The same with DefaultPath().
  static const CoulombTable *Shared(const std::string &path) { return Shared(path, DefaultPath()); }
  /// Forget which files Shared() could not read, so the next lookup tries them again
  /// (after a table was built while the process runs, e.g. from the GUI).
  static void ForgetMissing();
  /// $AZURE2_COULOMB_TABLE if set, else ~/.azure2/coulomb_table.bin, in the process's environment record.
  static std::string DefaultPath();
  /// The same in a given record (azure2/support/Environment.h).
  static std::string DefaultPathIn(const azure2::support::Environment &environment);

  /// The exact open-channel values (Michel's routine, AZURE2's default); phi in (-pi, pi].
  static Open ExactOpen(int l, double eta, double rho);
  /// The exact closed-channel shift function from GSL's U: S = z W'(z)/W(z), z = 2 rho_b,
  /// with W = exp(-z/2) z^(l+1) U(l+1+eta_b, 2l+2, z) and U' = -a U(a+1, b+1, z).
  static double ExactBound(int l, double etaB, double rhoB);
  /// log C0^2(eta) = log(2 pi eta / (exp(2 pi eta) - 1)), evaluated stably.
  static double LogGamow(double eta);

  int LMax() const { return lMax_; }
  const Axis &EtaAxis() const { return eta_; }
  const Axis &RhoAxis() const { return rho_; }
  const std::string &Generator() const { return generator_; }
  double EtaRhoMax() const { return etaRhoMax_; }
  /// Is (eta, rho) inside the tabulated domain?
  bool Inside(double eta, double rho) const;
  bool Empty() const { return sheets_.empty(); }

  /// What the table holds, as text: the format, the interpolation order, l_max, the two axes, eta*rho max, the
  /// generator and a 64-bit hash of every tabulated value (the sheets' bytes, 8 at a step).  Two tables give the same
  /// text when they give the same lookups, wherever their files lie and whenever they were written; a table built
  /// with other settings, or rebuilt with another generator, gives another.  The signature of the external-capture
  /// integrals holds it (WorkingCopy::ECSignature).  Computed on the first call (the 88 MB default table: a few
  /// times 1e7 instructions) and fixed after it: the value is a function of the loaded table alone, written once
  /// under std::call_once.  A table that is modified after the first call (Build, Load, SetOrder) is not supported.
  const std::string &ContentIdentity() const;

 private:
  // One tabulated function on (eta, rho) nodes with its Hermite derivative data.
  struct Sheet {
    std::vector<double> f, fx, fy, fxy;  // row-major: index = i * ny + j (i: eta, j: rho)
  };
  enum Kind { LOGP = 0, SHIFT = 1, PHASE = 2, BOUND = 3, NKIND = 4 };

  void Derivatives(Sheet &s) const;
  double Eval(const Sheet &s, double x, double y) const;
  struct Stencil {
    int i0, j0;
    double wx[8], wy[8];
  };
  Stencil MakeStencil(double x, double y) const;
  double Apply(const Sheet &s, const Stencil &st) const;
  bool Cell(double x, double y, int &i, int &j, double &u, double &v) const;
  const Sheet &At(int l, Kind k) const { return sheets_[l * NKIND + k]; }
  Sheet &At(int l, Kind k) { return sheets_[l * NKIND + k]; }

  int lMax_ = -1;
  int order_ = 8;
  double etaRhoMax_ = 0.0;  // 0: no limit
  Axis eta_, rho_;  // the bound sheets share the axes (eta_b, rho_b)
  std::string generator_;
  std::vector<Sheet> sheets_;
  // ContentIdentity(), filled on its first call; held through a pointer so the table stays movable.
  struct Identity {
    std::once_flag once;
    std::string text;
  };
  std::unique_ptr<Identity> identity_ = std::make_unique<Identity>();
};

#endif
