#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <fstream>
#include <memory>
#include <mutex>
#include <set>

#include "ThmOptions.h"

/// A structure holding the reaction rate calculation configuration

/*!
 * The RateParams structure holds the configuration information for a reaction rate calculation.
 */

struct RateParams {
  /// False for looped temperatures, true for temperatures from file.
  bool useFile;
  /// String containing filename with temperatures to use
  std::string temperatureFile;
  /// The entrance pair number for the rate calculation
  int entrancePair;
  /// The exit pair number for the rate calculation
  int exitPair;
  /// The minimum temperature for the rate calculation
  double minTemp;
  /// The maximum temperature for the rate calculation
  double maxTemp;
  /// The temperature step for the rate calculation
  double tempStep;
};

/// A configuration structure for AZURE

/*!
 * The configuration structure is created from the runtime file passed to the
 * AZURE executable, as well as the options specified in the command shell prompt.
 */

class Config {
 public:
  Config(std::ostream &stream);
  /// Restore every option to its default, as a fresh run starts.
  void Reset();
  /*!
   * Bit flags for various options in AZURE2.
   */
  enum ParameterFlags {
    USE_AMATRIX = (1 << 0),
    PERFORM_ERROR_ANALYSIS = (1 << 1),
    PERFORM_FIT = (1 << 2),
    CALCULATE_WITH_DATA = (1 << 3),
    USE_PREVIOUS_PARAMETERS = (1 << 4),
    USE_EXTERNAL_CAPTURE = (1 << 5),
    USE_PREVIOUS_INTEGRALS = (1 << 6),
    CALCULATE_REACTION_RATE = (1 << 7),
    TRANSFORM_PARAMETERS = (1 << 8),
    USE_BRUNE_FORMALISM = (1 << 9),
    IGNORE_ZERO_WIDTHS = (1 << 10),
    USE_RMC_FORMALISM = (1 << 11),
    USE_GSL_COULOMB_FUNC = (1 << 12),
    USE_LONGWAVELENGTH_APPROX = (1 << 13),
    USE_API = (1 << 14),
    PERFORM_MCMC = (1 << 15),
    USE_WIGNER_LIMITS = (1 << 16),
    USE_NLOPT_MINIMIZER = (1 << 17),
    USE_HYBRID_COULOMB = (1 << 18),
    USE_ANALYTIC_GRADIENT = (1 << 19),
    USE_LM_MINIMIZER = (1 << 20),
    CALCULATE_COVARIANCE_BAND = (1 << 21),
    SCALE_COVARIANCE_BY_CHI2 = (1 << 22),
    USE_GSL_LM_MINIMIZER = (1 << 23),
    /// Park (2021) level-dependent boundary conditions: the fit parameters are
    /// the observed reduced width amplitudes.  Always set together with
    /// USE_BRUNE_FORMALISM, whose per-level shift functions it shares.
    USE_PARK_FORMALISM = (1 << 24)
  };
  /*!
   * Bit flags for check file control in AZURE2.
   */
  enum CheckFileFlags {
    CHECK_COMPOUND_NUCLEUS = (1 << 0),
    CHECK_PATHWAYS = (1 << 1),
    CHECK_DATA = (1 << 2),
    CHECK_ENERGY_DEP = (1 << 3),
    CHECK_LEGENDRE = (1 << 4),
    CHECK_BOUNDARY_CONDITIONS = (1 << 5),
    CHECK_ANGULAR_DISTS = (1 << 6),
    CHECK_COUL_AMPLITUDES = (1 << 7)
  };
  /// Output stream
  std::ostream &outStream;
  /// The runtime configuration file name.
  std::string configfile;
  /// A control variable to stop AZURE calculation
  bool stopFlag;
  /// A bitmask for the encoding of configuration flags
  unsigned int paramMask;
  /// A bitmask storing which checks are printed to screen.
  unsigned int screenCheckMask;
  /// A bitmask storing which checks are printed to file.
  unsigned int fileCheckMask;
  /// If performError is true, sets the value of Up (the acceptable variance from the minimum chi-squared.
  double chiVariance;
  /// The path of the output directory
  std::string outputdir;
  /// The path of the check files directory.
  std::string checkdir;
  /// The name of the parameters file from which to read.
  std::string paramfile;
  /// The name of the external capture amplitudes file from which to read.
  std::string integralsfile;
  /// Parameters for calculating reaction rate.
  RateParams rateParams;
  /// NLopt algorithm selection (0=SBPLX, 1=COBYLA, 2=BOBYQA, 3=NEWUOA, 4=PRAXIS, 5=Nelder-Mead)
  int nloptAlgorithm;
  /// Use hybrid Coulomb method with nuclear potential
  bool useHybridMethod;
  /// Extra MIGRAD passes after re-freezing the model of relvector data-covariance terms at the
  /// fitted parameters (--refreeze-covariance[=N]); 0 = none.  The final outputs always use the
  /// model re-frozen at the fitted parameters, so they match a calculation from param.sav.
  int refreezeCovariance;
  /// Function-call limit of a MIGRAD or NLopt fit (--max-calls=N); 0 = none, so MIGRAD's EDM criterion (or NLopt's
  /// tolerances) ends the fit.  Before 2026-10-31 a fixed 50,000 that a fit of ~100+ free parameters reached silently.
  unsigned int maxCalls;
  /// Use adaptive integration grid for target effects (false = uniform grid)
  bool useAdaptiveGrid;
  /// The <thm> block (ThmOptions.h); Config::ThmOptions as before.
  typedef ::ThmOptions ThmOptions;
  ThmOptions thm;
  /// Take P, S and the hard-sphere phase of charged channels from the universal Coulomb
  /// table (CoulombTable) where it covers them; on by default (--no-coulomb-table turns it
  /// off).  Without a table file the run warns once and uses the exact routines.
  bool useCoulombTable;
  /// The table file; empty = $AZURE2_COULOMB_TABLE, else ~/.azure2/coulomb_table.bin.
  std::string coulombTablePath;
  /// External-capture integrals from per-term energy tables (ECTable) instead of a
  /// direct integral at every energy, and shifted points evaluated at their shifted
  /// energies from the same tables; on by default, --no-ec-table restores the direct
  /// integrals and the linear interpolation of shifted amplitudes.
  bool useECTable;
  /// With the long-wavelength approximation off (--no-long-wavelength), external capture also carries the
  /// convection-current (tangential) part of the electric multipole operator, not only its charge (extended
  /// Siegert) part.  On by default; --no-ec-convection-current leaves the charge part alone.
  bool ecConvectionCurrent;
  /// Nodes per decade of E_out of those tables ($AZURE2_EC_TABLE_NODES overrides).
  int ecTableNodesPerDecade;
  /// Lay the adaptive sub-point grids of experimental effects out again on every
  /// evaluation, for the current levels and energy shifts (EData::ReplaceSubPointGrids);
  /// on by default, --no-grid-replacement keeps the grids placed from the input file.
  bool replaceSubPointGrids;
  /// Hard-sphere scattering in every J^pi the particle pairs allow, whether or not a level carries it
  /// (CNuc::AddHardSphereGroups). Read from "includeAllJPi=1" in the <config> block; off when absent.
  bool includeAllJPi;
  /// Orbital momentum limit of that enumeration ("maxOrbitalMomentum=N" in <config>); -1 = the largest l of
  /// any particle channel in the <levels> block.
  int hardSphereMaxL;
  /// Width of the sum-to-one penalty on the isotopic abundances of a project with several compound
  /// nuclei: ((sum w) - 1)^2 / s^2 is added to chi-squared when s > 0.  Read from "abundanceSumSigma=s"
  /// in the <config> block; 0 (off) when absent.
  double abundanceSumSigma;
  /// Relativistic two-body kinematics (Kinematics.h): channel wave numbers, Sommerfeld parameters and
  /// lab <-> c.m. energies from the invariant mass.  Read from "relativistic=1" in the <config> block;
  /// off (non-relativistic, as before) when absent.
  bool relativistic;
  /// Key of the particle pair whose E(s) = (s - m0^2)/(2 m0) relativistic.out uses (the LANL/EDA
  /// level-energy convention); "relativisticReferencePair=<key>" in <config>, 0 = the first pair.
  int relativisticReferencePair;
  /// The Notes a session has printed (see Note()).  Null -- the CLI and the GUI -- prints every Note; a pyazr
  /// session (AZUREAPI) sets one so that the setup Notes of its rebuilds are printed once: a channel-radius
  /// scan rebuilds the model on every likelihood call, and before 2026-11-04 each rebuild printed every Note
  /// again (the 3H+d free-radius MCMC, 2026-11-03).  Not touched by Reset().
  struct NoteLog {
    std::mutex mutex;
    std::set<std::string> printed;
  };
  std::shared_ptr<NoteLog> noteLog;
  /// Prints one informational line ("Note: ...") to outStream, unless noteLog is set and has already seen
  /// exactly this text.  Warnings and errors are not routed through here: they print every time.
  void Note(const std::string &line) const;
  /// Read input files (the .azr, data, covariance and target-effect files) through InputFiles, which keeps
  /// their contents in memory and reads a file again only when it changes on disk.  Set by AZUREAPI: a
  /// rebuild there (SetRadius, a potential change, a mode switch) used to re-read every file, which on
  /// /groups stalled 16 forked MCMC workers (2026-11-03).  False -- read from disk -- for the CLI and GUI.
  /// Not touched by Reset().
  bool cacheInputFiles = false;
  /// A constant indicating the maximum order of the Legendre polynomials to calculate.
  static const int maxLOrder = 20;
  /// Read the <config> block, then the optional <potential> and <thm> blocks.
  /// -1 if the file cannot be read, -2 if an optional block is malformed (already reported).
  int ReadConfigFile();
  /// Reads the <potential> block of the configuration file, if it has one.
  int ReadPotentialBlock();
  /// Reads the optional <thm> block. -1 on an unknown key or an unterminated block.
  int ReadThmBlock();
  /// The R-matrix route (<config> first line false) implements neither IGNORE_ZERO_WIDTHS (the A-matrix route leaves
  /// the levels whose internal width is zero out of the channel-capture term; RMatrixFunc keeps them) nor Park's
  /// parametrization (RMatrixFunc would read the Park amplitudes as Brune's: tests/13N up to 9x off).  When either is
  /// on with the R-matrix route, switch to the A-matrix with a warning and return true.  Called by the CLI, pyazr
  /// (Session) and the GUI's run thread (AZUREMainThread).
  bool ResolveMatrixRoute();
#ifndef NO_STAT
  /// Check that the output and checks directories exist.
  int CheckForInputFiles();
#endif
};

// Global Config instance pointer for use across the codebase
extern Config *g_config;

#endif
