#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <fstream>
#include <map>
#include <memory>
#include <mutex>
#include <set>

#include "azure2/experiment/HybridPotentialSettings.h"
#include "ThmOptions.h"
#include "ProjectOptions.h"

struct ProgramObserver;
namespace AZUREProfile {
class Account;
}
namespace azure2 {
class RunScratch;
struct RunDiagnostics;
namespace experiment {
struct Format3Project;
}
namespace io {
class Document;
}
}  // namespace azure2

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
    // (1 << 12) was USE_GSL_COULOMB_FUNC, the option use_gsl_coul retired on 2026-11-16 (docs/upgrade/coulomb.md 9.1)
    USE_LONGWAVELENGTH_APPROX = (1 << 13),
    USE_API = (1 << 14),
    PERFORM_MCMC = (1 << 15),
    USE_WIGNER_LIMITS = (1 << 16),
    // (1 << 17), (1 << 19), (1 << 20), (1 << 23) were the fitters' bits (USE_NLOPT_MINIMIZER, USE_ANALYTIC_GRADIENT,
    // USE_LM_MINIMIZER, USE_GSL_LM_MINIMIZER), retired on 2026-11-17: the chosen fitter is Config::fitter, a row's
    // name (docs/upgrade/fit.md section 5).  The positions stay unused; a fitter never gets a bit here (check 46).
    USE_HYBRID_COULOMB = (1 << 18),
    CALCULATE_COVARIANCE_BAND = (1 << 21),
    SCALE_COVARIANCE_BY_CHI2 = (1 << 22),
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
  /// The directory a project's relative paths are resolved against (absolute): data, reference, grid, covariance and
  /// table files, the output and checks directories, intEC.dat, parameter, temperature and sample files.  The CLI
  /// sets it to the project file's directory; pyazr to the session's directory (azure2(cwd=...)).  Empty: the
  /// process's working directory (a Config that never named one).  Every path a project or a prompt gives is kept as
  /// written and printed as written; OpenPath gives the file to open.
  std::string baseDirectory;
  /// The file to open for `path` as the project wrote it: `path` itself when it is absolute or no base directory is
  /// set, else baseDirectory/path.  Never printed: messages show `path`.
  std::string OpenPath(const std::string &path) const;
  /// A control variable to stop AZURE calculation
  bool stopFlag;
  /// The front end that runs the program in its own process (include/ProgramObserver.h): its stop request and its
  /// sinks of progress and chi^2.  Null for the command line.  Owned by the caller; not reset by Reset.
  const ProgramObserver *observer = nullptr;
  /// stopFlag, or the observer's stop request.
  bool StopRequested() const;
  /// A bitmask for the encoding of configuration flags
  unsigned int paramMask;
  /// The fitter the program's switches chose: the name of a row of Fitters() (include/Fitters.h), empty for none (the
  /// row without a switch, migrad).  Set only by ChooseFitter (src/Fitters.cpp), which keeps the table's precedence;
  /// read by FitterOfConfig.  A new fitter needs no edit of this header (docs/upgrade/fit.md section 5).
  std::string fitter;
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
  /// Use hybrid Coulomb method with nuclear potential: potentials.AnyEnabled(), set with them
  bool useHybridMethod;
  /// The hybrid nuclear potentials of the project, per pair (its `potential` records; pyazr edits them in a Project)
  azure2::experiment::HybridPotentialSettings potentials;
  /// Extra MIGRAD passes after re-freezing the model of relvector data-covariance terms at the
  /// fitted parameters (--refreeze-covariance[=N]); 0 = none.  The final outputs always use the
  /// model re-frozen at the fitted parameters, so they match a calculation from param.sav.
  int refreezeCovariance;
  /// Function-call limit of a MIGRAD or NLopt fit (--max-calls=N); 0 = none, so MIGRAD's EDM criterion (or NLopt's
  /// tolerances) ends the fit.  Before 2026-10-31 a fixed 50,000 that a fit of ~100+ free parameters reached silently.
  unsigned int maxCalls;
  /// The <thm> block (ThmOptions.h); Config::ThmOptions as before.
  typedef ::ThmOptions ThmOptions;
  ThmOptions thm;
  /// The project as ReadConfigFile read it once (records -> Validate -> <levels> tables -> Schema::Build): the <levels>
  /// rows and the specs the set-up consumes (src/ProjectSetup.cpp).  Set by every successful ReadConfigFile (a file
  /// that is not format 3 is refused before it); null only before it.
  std::shared_ptr<const azure2::experiment::Format3Project> project;
  /// Where the value of each option of the calculation came from (include/ProjectOptions.h): an option absent here has
  /// its default.  Set by the program's flags (parseOptions), a session's keywords and ReadConfigFile's keys.
  std::map<std::string, OptionSource> optionSources;
  /// The options a flag or keyword named itself (not one it implied): an implication does not over-rule them, so that
  /// --use-park with --no-brune, in either order, meets as a contradiction instead of letting the order decide.
  std::set<std::string> optionsStated;
  /// The contradictions of options ApplyProjectOptions refused (their diagnostics, one line each); empty otherwise.
  std::string refusedOptions;
  /// Take P, S and the hard-sphere phase of charged channels from the universal Coulomb
  /// table (CoulombTable) where it covers them; on by default (--no-coulomb-table turns it
  /// off).  Without a table file the run warns once and uses the exact routines.
  bool useCoulombTable;
  /// The table file; empty = environment.coulombTableDefault.
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
  /// Nodes per decade of E_out of those tables (environment.ecTableNodes overrides).
  int ecTableNodesPerDecade;
  /// Lay the adaptive sub-point grids of experimental effects out again on every
  /// evaluation, for the current levels and energy shifts (WorkingCopy::ReplaceSubPointGrids);
  /// on by default, --no-grid-replacement keeps the grids placed from the input file.
  bool replaceSubPointGrids;
  /// Refinement factor k of every sub-point grid (AZURE_REFINE_QUADRATURE=k, a diagnostic of the quadrature gate,
  /// docs/upgrade/energy_grid.md section 7): the grids get k times the points per width and k times the smooth
  /// intervals of the window (the smooth step / k); the tail step per distance follows the points per width (/ k).
  /// 1 (the default): the grids as the project states them.  Set by the CLI's parseOptions and pyazr's session
  /// (RefineQuadratureFromEnvironment, CommandLineValues.h), which refuse a bad value.
  int refineQuadrature;
  /// Panels per Gaussian width of the depth integral of a target with a beam spread or straggling
  /// (experiment::KernelParameters::depthStepsPerSigma, docs/upgrade/energy_grid.md section 3.4): each panel sigma/k
  /// wide, never ending above the next node of the mesh; 0 = the node intervals (the integral before 2026-11-27).
  /// --fold-depth-steps=K, AZURE_FOLD_DEPTH_STEPS, or a pyazr session's fold_depth_steps; set by the CLI's
  /// parseOptions and the session (FoldDepthStepsFromEnvironment, CommandLineValues.h), which refuse a bad value.
  int foldDepthSteps;
  static constexpr int kFoldDepthStepsDefault = 4;
  /// Hard-sphere scattering in every J^pi the particle pairs allow, whether or not a level carries it
  /// (CompoundSystem::AddHardSphereGroups). Read from "includeAllJPi=1" in the <config> block; off when absent.
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
  /// Prior parameter covariance (PriorCovariance.h): "priorCovariance=<path>" in the <config> block, the path relative
  /// to the .azr's directory; empty (the default) = none.  A multivariate Gaussian prior on the listed parameters,
  /// (p - mu)^T C^-1 (p - mu) added to chi-squared (Total-Prior-Covariance-Chi-Squared).
  std::string priorCovariance;
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
  /// The warnings printed once per run (a fit evaluates thousands of times; the first case is enough), one key per
  /// message.  Owned by the run, so that a second run in the same process -- a pyazr session, the web UI's server
  /// -- prints them again (they were process-wide flags before 2026-11-12).  Copies of a Config share it.  Not
  /// touched by Reset().
  std::shared_ptr<NoteLog> printedOnce = std::make_shared<NoteLog>();
  /// True the first time this run asks with `key` (thread-safe): print the warning then, and not again.
  bool FirstTime(const std::string &key) const { return FirstTimeIn(printedOnce, key); }
  /// The same for a holder of the run's log that has no Config at hand (set-up records keep the pointer); a null
  /// log answers true every time.
  static bool FirstTimeIn(const std::shared_ptr<NoteLog> &log, const std::string &key);
  /// The run's per-thread scratch (RunScratch.h): what the external-capture integrals and the THM vertex boundary keep
  /// between calls, keyed by everything it depends on.  Owned by the run, never process-wide; copies of a Config share
  /// it.  Made by the constructor, not touched by Reset().
  std::shared_ptr<azure2::RunScratch> scratch;
  /// What the run's diagnostics count (RunDiagnostics.h: the profile accounting, the Coulomb table's lookups, the
  /// call counters of the printing switches).  Owned by the run, shared by copies of the Config, as scratch is; made
  /// by the constructor, not touched by Reset().
  std::shared_ptr<azure2::RunDiagnostics> diagnostics;
  /// The run's profile account (diagnostics->profile) when AZURE_PROFILE is on, else null: what every Timer of the run
  /// is given (Profile.h).
  AZUREProfile::Account *profile = nullptr;
  /// Read input files (the .azr, data, covariance and target-effect files) through InputFiles, which keeps
  /// their contents in memory and reads a file again only when it changes on disk.  Set by AZUREAPI: a
  /// rebuild there (SetRadius, a potential change, a mode switch) used to re-read every file, which on
  /// /groups stalled 16 forked MCMC workers (2026-11-03).  False -- read from disk -- for the CLI and GUI.
  /// Not touched by Reset().
  bool cacheInputFiles = false;
  /*!
   * The settings the process's environment gives the run: the switches that can move a number
   * (azure2/support/Environment.h; docs/upgrade/evaluation.md section 8.2), copied from the process's record by
   * Reset().  The code reads them here, never the environment, so that a command-line or pyazr option for one of them
   * is one assignment after Reset (the PI's decision; none exists).
   */
  struct EnvironmentSettings {
    std::string coulombTableDefault;   ///< AZURE2_COULOMB_TABLE, else <home>/.azure2/coulomb_table.bin
    std::string srimDataDirectory;     ///< AZURE2_SRIM_DATA_DIR; "" = the stopping library's compiled copy
    int ecTableNodes = 0;              ///< AZURE2_EC_TABLE_NODES; > 0 overrides ecTableNodesPerDecade
    bool boundShiftNumerical = false;  ///< AZURE_BOUND_SHIFT_NUMERICAL: bound S from the numerical derivative
    bool noNeutralClosedForm = false;  ///< AZURE_NO_NEUTRAL_CLOSED_FORM
    bool shiftDerivative1eV = false;   ///< AZURE_SHIFT_DE_1EV: dS/dE by a 1 eV difference of the exact S
    bool levelSpaceOnly = false;       ///< AZURE_LEVEL_SPACE: no channel-space solves
    bool thmByDifferences = false;     ///< AZURE2_THM_FD: THM rows of the gradient by differences
    bool gradientByAdjoint = false;    ///< AZURE_GRAD_ADJOINT: the adjoint gradient, not the Jacobian path
    bool denseAdjoint = false;         ///< AZURE_GRAD_DENSE: the dense level-matrix adjoint
    bool foldWeightsByDifferences = false;  ///< AZURE_TE_WEIGHTS_FD: combiner weights by differences
    bool serialPointPass = false;      ///< AZURE_SERIAL_CHI2: the points one by one (the same bits)
  };
  EnvironmentSettings environment;
  /// A constant indicating the maximum order of the Legendre polynomials to calculate.
  static const int maxLOrder = 20;
  /// What ReadConfigFile returns besides 0.
  static constexpr int kConfigUnreadable = -1;      ///< the file cannot be read or has no <config> section
  static constexpr int kConfigBlockMalformed = -2;  ///< the <thm> section is malformed; reported
  static constexpr int kConfigRefused = -3;         ///< a format-3 project is refused; its diagnostics are printed
  /// the file is not a format-3 project (format 1 or 2, or not a project): refused before anything is read, with
  /// io::NotFormat3Message (docs/upgrade/project_format3.md section 8.4); the command line exits with 3
  static constexpr int kConfigNotFormat3 = -4;
  /// Read the project (format 3 only: the records into `project`), its <config> block, its potentials and <thm>: 0, or
  /// one of the codes above.
  int ReadConfigFile();
  /// The format-3 branch of ReadConfigFile: reads the project into `project`; kConfigRefused (with every diagnostic
  /// printed, the last line from ReportFormat3Refused) when it is refused.
  int ReadFormat3Project(const azure2::io::Document &doc);
  /// The last line of a refused format-3 project, after its diagnostics: the file is refused and --validate lists
  /// the same diagnostics.  One text for every place a format-3 project is refused (reading, its data files).
  void ReportFormat3Refused() const;
  /*!
   * The section <name> of the .azr (configfile) as a stream: its opening tag, its lines and its closing tag (each
   * line ending read as "\n"; to the end of the file if the section is not closed), or null when the file cannot be
   * read or has no such top-level section.  The file is read through InputFiles (cached per cacheInputFiles, so the
   * section is that of the file as it now is) and parsed as an io::Document (docs/upgrade/project.md, Phase 3): a tag
   * is a whole trimmed line.  `openLine`, when given, gets the 1-based line of the section's opening tag.
   */
  std::unique_ptr<std::istream> SectionStream(const std::string &name, int *openLine = nullptr) const;
  /// The hybrid potentials of the project's <potentials> records (project->built.potentials) into `potentials`.
  void ReadPotentials();
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

#endif
