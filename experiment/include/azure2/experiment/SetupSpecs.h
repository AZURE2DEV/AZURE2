#ifndef AZURE2_EXPERIMENT_SETUPSPECS_H
#define AZURE2_EXPERIMENT_SETUPSPECS_H

// What the engine's set-up consumes (data_model.md S9, project_format3.md step 3b): one plain struct per thing the
// engine builds -- a data set, an extrapolation, a target / beam / detector effect with the data sets it applies to,
// the abundances of the samples, and a model datum.  A spec is filled once by its producer (today the old-line adaptor of convert/,
// SpecsFromOldLines.h; after step 3b also Schema::Build from format-3 records) and then only read: the engine's
// constructors take it and nothing else.  Plain values, standard library only; energies in MeV, angles in degrees,
// as the members say.  Where the engine needs a number as it was written (a token it checks or quotes in a message),
// the member holds the text.

#include <string>
#include <utility>
#include <vector>

namespace azure2 {
namespace experiment {

/// What a data set or an extrapolation computes at each point (the model's observable, before any apparatus).
enum class Observable {
  CrossSection,          ///< angle-integrated cross section (b), or the transformed observable of a sample block
  CrossSectionE1,        ///< angle-integrated capture, compared with the E1 part only
  CrossSectionE2,        ///< angle-integrated capture, compared with the E2 part only
  TotalCapture,          ///< angle-integrated capture summed over every capture pair of the entrance nucleus
  DifferentialLab,       ///< dsigma/dOmega (b/sr) at a lab angle
  DifferentialCM,        ///< dsigma/dOmega (b/sr) at a c.m. angle
  AnalyzingPower,        ///< vector analyzing power A_y at a c.m. angle
  PhaseShift,            ///< the nuclear phase shift (deg) of one partial wave
  AngularCoefficients,   ///< Legendre coefficients of the c.m. angular distribution as data (data sets only)
  AngularDistribution,   ///< Legendre coefficients up to an order, computed (extrapolations only)
};

/// The Coulomb cross section a ratio-to-Coulomb data set is divided by (data_kinds.md section 4.3): none, Rutherford
/// (a charged, non-identical elastic pair), or Mott (identical particles, the spin from the pair).
enum class CoulombReference { None, Rutherford, Mott };

/*!
 * A factor on the model value of a data set or an extrapolation (format 3: the quantity and the unit of the values;
 * the tree kind Scaled of experiment/Prediction.h): the value compared with the data is
 *
 *   factor x [E exp(2 pi eta(E)) if sFactor] x (the observable's tree),
 *
 * E the c.m. energy of the entrance pair at the point as the model evaluates it (after an energy shift; MeV) and eta
 * its Sommerfeld parameter (physics/AstrophysicalFactor.h), so that a fitted shift moves the factor.  The identity
 * (factor 1, no sFactor) is no node at all, and every old line has it.
 */
struct ValueScale {
  double factor = 1.0;   ///< the constant c: 1e-6 C / (q e) for a yield per uC, 1 / (unit factor) for an extrapolation
  bool sFactor = false;  ///< the astrophysical S factor S = E sigma exp(2 pi eta) of an angle-integrated sigma
  std::string text;      ///< what the factor is, in words, for --describe ("" when factor is 1)
  bool Identity() const { return factor == 1.0 && !sFactor; }
};

/*!
 * How a row of a data file becomes a point (format 3: the units of the column roles, the unit of the values, the
 * frame of the energy column; experiment/DataRows.h, RowToEngine, applied where the points are made).  The identity
 * leaves every number's bits as read, and every old line has it.
 */
struct ColumnConversion {
  double energy = 1.0;     ///< MeV per unit of the energy column (and of window_lo, window_hi)
  bool energyCM = false;   ///< the energy column is the entrance pair's E_cm: taken to E_lab (Frames.h BeamLabFromCM)
  double lightMass = 0.0, heavyMass = 0.0;  ///< u, the entrance pair's masses (energyCM)
  bool relativistic = false;                ///< the run's kinematics (energyCM)
  double angle = 1.0;      ///< deg per unit of the angle column (1 for an order column)
  double value = 1.0;      ///< engine units of the quantity per unit of the value and error columns
  bool Identity() const { return energy == 1.0 && !energyCM && angle == 1.0 && value == 1.0; }
};

/// One term of a composite (SUM / RATIO) data set or extrapolation.
struct ComponentSpec {
  int entranceKey = 0;   ///< pair key of <levels>
  int exitKey = 0;       ///< pair key of <levels>
  bool hasAngle = false; ///< false: the term is evaluated at the data set's own angles
  double angle = 0.0;    ///< fixed lab angle (deg) when hasAngle
  double scaling = 1.0;  ///< the term's factor (1 when none is given)
};

/// The fields an extrapolation and a data set share: the reaction, the observable and the model-side options.
struct ReactionSpec {
  int entranceKey = 0;  ///< pair key of <levels>
  int exitKey = 0;      ///< pair key of <levels>; -1 = every capture pair (TotalCapture)
  Observable observable = Observable::CrossSection;
  bool thm = false;     ///< a THM (half-off-shell) observable of the same kind
  double phaseJ = 0.0;  ///< PhaseShift: total spin J of the partial wave
  int phaseL = 0;       ///< PhaseShift: orbital angular momentum l
  double phaseS = -1.0; ///< PhaseShift: channel spin s, < 0 when not given
  enum class Combination { None, Sum, Ratio };
  Combination combination = Combination::None;  ///< a composite of the leading term and `components`
  std::vector<ComponentSpec> components;
  std::string sample = "0";  ///< "0": no sample; "1": the default sample; else a sample name of the abundances
  bool total = false;        ///< the total cross section of the entrance pair (optical theorem)
  double transmission = 0.0; ///< > 0: transmission exp(-n sigma) of a sample n atoms/b thick
  double yieldThickness = 0.0;  ///< > 0: self-shielded yield of a sample n atoms/b thick
  /// Not None: the value is dsigma/dOmega (c.m.) over this Coulomb cross section at the same c.m. energy and angle,
  /// each folded by the point's apparatus (quantities ratio_to_rutherford, ratio_to_mott; observable DifferentialLab
  /// or DifferentialCM, elastic).  The value is the same in both frames: it is not converted at load.
  CoulombReference coulombReference = CoulombReference::None;
  /// A factor on the model value (an S factor, a yield per unit charge, an extrapolation in another unit).
  ValueScale scale;
};

/// The names of a data set's or an extrapolation's fit parameters, given by the spec's producer
/// (docs/upgrade/parameter_kinds.md: a parameter's name is given where it is made, never read for its kind).  The
/// old-line adaptor writes the old names (segment_<key>_norm, ...); Schema::Build the format-3 names: of the data set
/// (<name>.norm, <name>.shift.offset, <name>.shift.sqrt), of a shared record for the first active member of its group
/// (<norm record>.value, <energy_scale record>.offset, <energy_scale record>.a<j>, <target>.areal_density), and for a
/// later member, whose parameter the engine keeps fixed, <record>.<field>[<data set name>].
struct ParameterNames {
  std::string norm;             ///< used when the norm varies and is not profiled
  std::string energyShift;      ///< every data set has this parameter (fixed at 0 without a shift)
  std::string energyShiftSqrt;  ///< likewise
  std::vector<std::string> shiftFunction;  ///< a0, a1, ... of the function the data set owns
  std::string thickness;        ///< used with a fitted sample thickness
};

/// A data set: points read from a file and compared with the model.
struct DatasetSpec {
  int key = 0;              ///< the data set's number, 1, 2, ... in file order (inactive ones counted)
  bool active = true;
  ParameterNames names;
  ReactionSpec reaction;
  double minE = 0.0, maxE = 0.0;  ///< lab energy window of the points used (MeV)
  double minA = 0.0, maxA = 0.0;  ///< lab angle window (deg), differential observables
  std::string dataFile;
  ColumnConversion columns;       ///< the rows of dataFile into the engine's units and frame
  bool absolute = false;          ///< AngularCoefficients: the data are A_k (b/sr), not a_k = A_k / A_0
  // Normalization n: the model is compared with n x data ... (EData; value, uncertainty in percent, vary).
  double norm = 1.0;
  double normError = 0.0;          ///< prior width (%), 0 = no prior
  bool varyNorm = false;
  bool hasNormCentre = false;      ///< the prior is centred on normCentre, not on `norm`
  double normCentre = 0.0;
  std::string normCentreSource;    ///< the text that set it, for messages ("" when none)
  std::string normGroup;           ///< shared normalization: data sets with one label share the first's norm
  // Energy shift E' = E + shift + b sqrt(E) + f(E) (MeV; E the lab energy in MeV).
  double energyShift = 0.0;
  double energyShiftError = 0.0;   ///< prior width (MeV), 0 = no prior
  bool varyEnergyShift = false;
  bool hasShiftCentre = false;     ///< the prior is centred on shiftCentre, not on 0
  double shiftCentre = 0.0;
  std::string shiftCentreSource;
  bool hasSqrtShift = false;       ///< a sqrt(E) term is declared
  double energyShiftSqrt = 0.0;    ///< b (MeV^1/2)
  double energyShiftSqrtError = 0.0;
  bool varyEnergyShiftSqrt = false;
  std::string shiftFunction;       ///< f(x) with x the lab energy (MeV), "" = none; parameters a0, a1, ...
  std::vector<double> shiftFunctionValues, shiftFunctionErrors;
  std::vector<bool> shiftFunctionVary;
  std::string shiftGroup;          ///< data sets with one label share the first's constant shift
  std::string shiftFunctionGroup;  ///< data sets with one label share the first's shift function
  std::string covarianceFile;      ///< "" = none; data sets naming one file share its correlations
  bool thicknessFit = false;       ///< the sample thickness is a fit parameter
  bool thicknessVary = false;
  double thicknessSigma = 0.0;     ///< prior width (atoms/b), 0 = none
  std::string thicknessGroup;      ///< data sets with one label share the first's thickness
  // Unobserved primary, observed secondary transition (UPOS).
  bool upos = false;
  int secondaryDecayL = 0;
  double secondaryFinalSpin = 0.0;   ///< I_c
  double secondaryMixingRatio = 0.0; ///< delta
};

/// An extrapolation: points on a grid, no data.
struct ExtrapolationSpec {
  int key = 0;  ///< 1, 2, ... in file order (inactive ones counted)
  bool active = true;
  ParameterNames names;
  ReactionSpec reaction;
  double minE = 0.0, maxE = 0.0;  ///< the energy grid's bounds (MeV, lab), as the segment reports them
  double minA = 0.0, maxA = 0.0;  ///< the angle grid's bounds (deg)
  std::vector<double> energies;  ///< lab energies of the grid (MeV), in order
  std::vector<double> angles;    ///< angles of the grid (deg; lab, or c.m. for the c.m. observables), in order; the
                                 ///< energies are repeated at each angle
  int maxAngDistOrder = 0;       ///< AngularDistribution: highest Legendre order
  std::string gridFile;          ///< a `file` energy grid: the path as opened ("" for the other grids)
  /// An adaptive energy grid (project_format3.md section 3.8): `energies` holds its start grid, the uniform grid of
  /// start_step over [minE, maxE] (lab), which the extrapolation run refines at the run's parameters
  /// (src/AdaptiveExtrapolation.cpp with experiment/AdaptiveGrid.h).
  struct Adaptive {
    bool on = false;
    double tolerance = 0.0;  ///< the relative deviation from linear interpolation that splits an interval
    long maxPoints = 0;      ///< energies of the refined grid at most
  } adaptive;
};

/// An apparatus effect (target integration, beam resolution, beam profile, Doppler, UDR, screening) and the data
/// sets it applies to.  Grammar and meaning of every field: docs/source/reference/data_formats.rst, <targetInt>.
struct TargetEffectSpec {
  std::string source;  ///< the text it was read from, quoted in messages
  /// The target's name, which names the effect's fit parameters (<target>.areal_density, <target>.doppler,
  /// <target>.screening); "" for old lines (density_<effect>, doppler_<effect>, screening_<effect>)
  std::string targetName;
  bool active = false;
  std::vector<int> datasetKeys;  ///< DatasetSpec::key / ExtrapolationSpec::key of the data sets it applies to
  int numIntegrationPoints = 0;
  bool convolution = false;      ///< Gaussian beam-energy resolution
  double sigma = 0.0;            ///< its width (MeV; frame: sigmaLab)
  bool targetIntegration = false;
  double density = 0.0;          ///< target thickness (atoms/cm^2); <= 0: an infinitely thick target
  std::string stoppingPowerEquation;
  std::vector<double> stoppingPowerParameters;
  bool qCoefficients = false;
  std::vector<double> qCoefficientValues;
  bool convolutionCoefficients = false;  ///< energy-dependent resolution sigma(E) by an equation
  std::string convolutionEquation;
  std::vector<double> convolutionCoefficientValues;
  bool hasStragglingFlag = false;
  bool straggling = false;
  bool hasStragglingCoefficient = false;
  double stragglingCoefficient = 0.0;
  bool hasResonanceWidthMultiplier = false;
  double resonanceWidthMultiplier = 0.0;
  bool hasPointsPerWidth = false;
  double pointsPerWidth = 0.0;
  bool hasRanges = false;                       ///< the lab-energy windows are given ("" = unrestricted)
  std::vector<std::pair<double, double>> ranges;  ///< lab energy windows (MeV)
  bool autoWindows = false;
  bool hasAutoWindowTolerance = false;
  double autoWindowTolerance = 0.0;
  bool hasAutoWindowMargin = false;
  double autoWindowMargin = 0.0;
  double transitionWidth = 0.0;  ///< MeV, 0 when absent
  double autoTolerance = 0.0;    ///< 0 when absent
  struct ProfileComponent {
    double xi = 0, omega = 0, alpha = 0, weight = 0;
  };
  bool beamProfile = false;
  std::vector<ProfileComponent> beamProfileComponents;
  double beamTpcSigma = 0.0;
  double beamTruncation = 0.0;
  bool beamDetailedBalance = false;
  bool sigmaLab = false;
  bool screening = false;
  bool screeningMalformed = false;  ///< old-line reader: the screening block could not be read
  double screeningUe = 0.0;         ///< electron-screening potential (eV)
  bool screeningVary = false;
  double screeningSigma = 0.0;      ///< prior width (eV)
  bool densityFit = false;          ///< density fit block present
  bool densityVary = false;
  std::string densitySigmaText;     ///< prior width as written: a number (same unit as density) or a percentage
  double densityMax = 0.0;          ///< largest density the grids span (atoms/cm^2), 0 = automatic
  std::string densityReadError;     ///< old-line reader: why the density block could not be read ("" = it could)
  std::string densityGroup;
  bool doppler = false;
  double dopplerTemperature = 0.0;  ///< K
  bool dopplerFit = false;
  std::string dopplerVaryText, dopplerSigmaText, dopplerMaxText;  ///< as written: 0|1, K, K ("" = none)
  std::string dopplerReadError;     ///< old-line reader: why the Doppler block could not be read
  bool udr = false;
  std::string udrFile;
  double udrFlightPath = 0.0;    ///< m
  double udrBurstFwhm = 0.0;     ///< ns
  double udrChannelWidth = 0.0;  ///< ns
  bool udrCentred = true;
};

/// One row of the abundances of the samples, in the order given: rows are applied in that order (a sample exists from
/// the first row that sets one of its abundances).
struct AbundanceRowSpec {
  std::string sample;   ///< "" = the unnamed sample
  /// The isotope as format 3 writes it ("13C"), which names the parameter <sample>.isotope.<isotope>; "" for old
  /// lines (abundance_[<sample>_]<Z>_<A>)
  std::string isotope;
  bool isSum = false;   ///< true: the width of the sample's sum-to-one prior; false: one nucleus's abundance
  double sumSigma = 0.0;
  int z = 0, a = 0;     ///< the compound nucleus
  double value = 0.0;   ///< abundance (fraction)
  bool vary = false;
  double sigma = 0.0;   ///< prior width (fraction)
};

/// The name of an abundance's fit parameter, made here only (CNuc::AbundanceParamName and the build's parameter sources
/// call it): `<sample>.isotope.<isotope>` for a format-3 isotope, else `abundance_[<sample>_]<Z>_<A>` with the Z and A of
/// the compound nucleus (BuildEffects.cpp).
std::string AbundanceParameterName(const std::string &sample, const std::string &isotope, int z, int a);

/// The abundances block: absent (no abundance is a parameter) or present with its rows (possibly none).
struct AbundancesSpec {
  bool present = false;
  std::vector<AbundanceRowSpec> rows;
};

/// What a reader made of one section: the specs in order, each with the reader's messages to print before the engine
/// sets it up, and why the reader stopped after the last item (`failed`; `error` is "" when the old engine stopped
/// without a message).  The engine sets the items up in order and reports `error` after them, so that the messages
/// come in the order of the section.
template <class Spec>
struct SpecList {
  struct Item {
    Spec spec;
    std::string warning;      ///< printed as "WARNING: ..." before the item ("" = none)
    std::string note;         ///< passed to Config::Note before the item ("" = none)
    bool readFailed = false;  ///< the reader could not read all of the item's text (old <targetInt> lines)
  };
  std::vector<Item> items;
  bool failed = false;
  std::string error;  ///< the whole message, "ERROR: ..." ("" = none)
};

/// The data sets of a project and what they link to, as EData::Fill consumes them.
struct DataSpecs {
  SpecList<DatasetSpec> datasets;
  /// A prior-centre row the old-line reader refused; reported where the old engine reported it (after the groups of
  /// the data sets are resolved).  "" = none.
  std::string priorCentreError;
  SpecList<TargetEffectSpec> targetEffects;
};

/// The extrapolations of a project and what they link to, as EData::MakePoints consumes them.
struct ExtrapolationSpecs {
  SpecList<ExtrapolationSpec> extrapolations;
  SpecList<TargetEffectSpec> targetEffects;
};

/// THM data sets (docs/upgrade/THM.md): the members of each experiment[<name>] of <thm> and the weight files of
/// single data sets and extrapolations (<thm> weight[k] / weightTest[k] in old files, `thm_weight` in format 3).  The
/// rest of an experiment is read from <thm> by its own parser.
struct ThmMembershipSpec {
  struct Experiment {
    std::string name;
    std::vector<int> datasetKeys;  ///< DatasetSpec::key, ascending
  };
  std::vector<Experiment> experiments;  ///< in the order the data sets name them first
  std::vector<std::pair<int, std::string>> datasetWeights;        ///< (DatasetSpec::key, weight file as written)
  std::vector<std::pair<int, std::string>> extrapolationWeights;  ///< (ExtrapolationSpec::key, weight file)
};

/// The abundances as CNuc consumes them, and why the reader stopped (`error` printed after the rows are applied).
struct AbundanceSpecs {
  AbundancesSpec abundances;
  bool failed = false;
  std::string error;
};

/// Model data (docs/upgrade/project_format3.md section 3.10, ModelDatum.h): a literature value of a level property or
/// of a scattering length, made by BuildModelData from a `datum` record.  Values are in the engine unit of the
/// quantity's dimension (MeV for energies and widths, 1, fm, fm^-1/2), with the unit and factor as written kept for
/// the output.
/// The quantities a datum can be, each evaluated by the engine (src/ModelData.cpp).  ModelQuantities() maps the
/// table's names (Quantities.cpp, record kind datum) to them; repo_consistency check 27 fails for a datum quantity of
/// the table without a row there and a test named after it.
enum class ModelQuantity {
  ResonanceStrength,  ///< omega gamma = omega Gamma_a Gamma_b / Gamma (data_kinds.md eq. 8)
  ResonanceEnergy,    ///< the observed level energy E_lambda, as E_x, E_cm or E_lab of a pair
  TotalWidth,         ///< Gamma = sum_c |Gamma_c|
  PartialWidth,       ///< Gamma of one pair: the sum over its channels
  BranchingRatio,     ///< Gamma_b / Gamma
  Anc,                ///< |C| of a closed channel (pair, l, s)
  ScatteringLength    ///< Re and Im of the s-wave scattering length a_0 (include/ScatteringLength.h)
};

struct ModelDatum {
  std::string name;
  std::string source;
  int line = 0;  ///< the record's header line
  ModelQuantity quantity = ModelQuantity::ResonanceStrength;
  std::string quantityName;
  double value = 0.0;  ///< in the engine unit of the quantity's dimension
  double error = 0.0;  ///< > 0, same unit
  bool hasImag = false;  ///< scattering length: Im a_0 given
  double imag = 0.0, imagError = 0.0;
  std::string unit;  ///< the unit as written (for the output)
  double unitFactor = 1.0;  ///< the written unit in engine units
  // Level quantities: the level.
  int levelId = 0;
  int group = 0, level = 0;  ///< the engine's indices (LevelInfo)
  std::string jpi;           ///< as written (3-)
  // Pairs, as <levels> keys (0 when the quantity has none).
  int entrancePair = 0;  ///< resonance_strength: the entrance pair a + A
  int exitPair = 0;      ///< resonance_strength, branching_ratio: the exit pair b + B
  int pair = 0;          ///< resonance_energy (E_cm, E_lab), partial_width, anc, scattering_length
  std::string frame;     ///< resonance_energy: E_x, E_cm or E_lab
  int l = -1;            ///< anc
  int twoS = -1;         ///< anc, scattering_length (2s)
  bool coherent = false; ///< scattering_length: the spin-weighted sum over the channel spins
};

}  // namespace experiment
}  // namespace azure2

#endif
