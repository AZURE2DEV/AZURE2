#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <fstream>
#include <map>
#include <memory>
#include <atomic>
#include <vector>

/*!
 * An energy-dependent weight w(E) multiplying the THM (HOES) model cross
 * section of one segment, read from a two-column file named by the
 * `weight[<segment key>]=` key of the <thm> block (Config::ReadThmBlock).
 * It is the hook for a correction AZURE2 does not compute itself, typically
 * Mukhamedzhanov's Coulomb-distortion factor R(E) = |M(E)|^2/|M(E_norm)|^2
 * (docs/source/theory/thm_implementation.rst).  Interpolated log-linearly
 * (linear in E, linear in ln w); constant beyond the ends of the table.
 */
struct ThmWeightTable {
  /// The file as named in the <thm> block, and as resolved against the .azr.
  std::string name, path;
  /// Strictly increasing E_cm (MeV) of the THM entrance pair, w(E) and ln w(E).
  std::vector<double> e, w, lnw;
  /// Set once the first evaluation outside [e.front(), e.back()] has been reported.
  mutable std::atomic<bool> warned{false};
  /// Reads the file; returns an empty string on success, else the reason.
  std::string Read(const std::string &file);
  /// w(E); *outside (if given) is set when E lies beyond the table.
  double operator()(double energy, bool *outside = nullptr) const;
  bool Covers(double energy) const { return energy >= e.front() && energy <= e.back(); }
};

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
    USE_GSL_LM_MINIMIZER = (1 << 23)
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
  /// Use adaptive integration grid for target effects (false = uniform grid)
  bool useAdaptiveGrid;
  /// Choices in the THM (modified R-matrix, HOES) observable, read from the
  /// optional <thm> block.  The defaults reproduce the Typel--Baur / La Cognata
  /// working formula (Tumino et al. 2021, eqs. 50-51) that the code has always
  /// used, except that entrance partial waves are summed incoherently.
  struct ThmOptions {
    /// Boundary in the entrance vertex (B - 1) j_l - rho j_l' (key `vertex=`):
    ///  perlevel  S_c(E_lambda) of each level under Brune, as mrmpy (alias `real`);
    ///  constant  B_c = S_c(E_1) at the lowest-energy level of the J group,
    ///            applied after the level sum.  gamma^T A gamma is
    ///            the same in the Brune and formal representations, so this is
    ///            the formal-R-matrix vertex at that B (La Cognata et al.,
    ///            ApJ 723 (2010) 1512 use B = S(E_1));
    ///  onshell   L_c(E) = S_c(E) + i P_c(E) (Tribble et al. 2014 eq. 2.76),
    ///            independent of any boundary constant.
    /// perlevel and constant coincide for a J group with a single level.
    enum Vertex { PER_LEVEL, CONSTANT, ON_SHELL };
    Vertex vertex = CONSTANT;
    /// Which kinematic factors the model carries, set by how the HOES data
    /// were extracted from the triple cross section (key `kinematics=`):
    ///  lacognata  exit k_f/mu_f, as La Cognata's working formula (Tumino et
    ///             al. 2021 eq. 50) -- the default, what the code always did;
    ///  triple     exit Gamma_f alone: raw d3sigma / |phi|^2 (Mukhamedzhanov
    ///             2017 eq. 34);
    ///  kf3body    exit Gamma_f / (mu_f k_f): data divided by the full
    ///             three-body phase-space KF (Typel & Baur 2003 eq. 16);
    ///  lambda32   exit Gamma_f, entrance 1/k_i on shell: data divided by
    ///             KF = lambda3/lambda2 (Pizzone et al. 2011; Tumino eq. 22).
    enum Kinematics { LA_COGNATA, TRIPLE, KF_THREE_BODY, LAMBDA32 };
    Kinematics kinematics = LA_COGNATA;
    /// Sum entrance partial waves of one channel spin coherently, as mrmpy and
    /// AZURE2 before Sep 2026 did.  For an observable integrated over the exit
    /// direction the l cross terms vanish, so incoherent is the default.
    /// Key `entranceL=coherent|incoherent`.
    bool coherentL = false;
    /// Add the external Coulomb term 2 eta k \int_a^inf O_l(kr)/O_l(ka) j_l(pr) dr to
    /// the vertex (Tribble eq. 2.79; Typel & Baur eq. A.4).  Key `coulombIntegral=0|1`.
    bool coulombIntegral = false;
    /// Mean spectator kinetic energy <p_sx^2>/2mu_sx (MeV) added to E + B in the
    /// half-off-shell momentum (Typel & Baur eq. 11).  Key `spectatorEnergy=` for
    /// every THM pair, `spectatorEnergy[<pair key>]=` for one.
    double spectatorEnergy = 0.0;
    std::map<int, double> spectatorEnergyByPair;
    double SpectatorEnergy(int pairKey) const {
      std::map<int, double>::const_iterator it = spectatorEnergyByPair.find(pairKey);
      return it == spectatorEnergyByPair.end() ? spectatorEnergy : it->second;
    }
    /// Energy-dependent weight of the THM model cross section, per segment:
    /// key `weight[<k>]=<file>` for the k-th line of <segmentsData> (counting
    /// inactive lines, the numbering of segment_<k>_norm and the output
    /// files), `weightTest[<k>]=<file>` for the k-th line of <segmentsTest>.
    /// A relative path is taken from the directory of the .azr file.
    std::map<int, std::shared_ptr<const ThmWeightTable>> weightBySegment, weightByTestSegment;
  };
  ThmOptions thm;
  /// A constant indicating the maximum order of the Legendre polynomials to calculate.
  static const int maxLOrder = 20;
  /// Read the <config> block, then the optional <potential> and <thm> blocks.
  /// -1 if the file cannot be read, -2 if an optional block is malformed (already reported).
  int ReadConfigFile();
  /// Reads the <potential> block of the configuration file, if it has one.
  int ReadPotentialBlock();
  /// Reads the optional <thm> block. -1 on an unknown key or an unterminated block.
  int ReadThmBlock();
#ifndef NO_STAT
  /// Check that the output and checks directories exist.
  int CheckForInputFiles();
#endif
};

// Global Config instance pointer for use across the codebase
extern Config *g_config;

#endif
