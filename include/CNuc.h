#ifndef CNUC_H
#define CNUC_H

#include <string>
#include <map>
#include "JGroup.h"
#include "PPair.h"

namespace ROOT {
namespace Minuit2 {
class MnUserParameters;
}
}  // namespace ROOT
class Config;

/// An AZURE compound nucleus

/*!
 * The compound nucleus is the fundamental concept of R-Matrix theory.  As such, the CNuc object
 * in AZURE is the top level container object for all structure and reaction objects.  Specifically,
 * the CNuc object is the container object for vectors of PPair and JGroup objects, within which all other
 * nuclear data objects are contained.
 */

/// One compound nucleus of a project: the (Z, A) its pairs form, and the
/// isotopic abundance of the target it stands for in a natural-target
/// segment (see the <abundances> block).
struct CompoundNucleus {
  int Z = 0;
  int A = 0;
  double abundance = 1.0;   ///< current value (a fit parameter when vary is set)
  double nominal = 1.0;     ///< value read from the file, centre of the prior
  double sigma = 0.0;       ///< prior width; 0 = no prior
  bool vary = false;
};

class CNuc {
 public:
  // -- structure ------------------------------------------------------------
  /// Number of compound nuclei the pairs form.  One for every project written
  /// before several nuclei could share a file.
  int NumNuclei() const;
  /// Nucleus \p i, 1-based.
  CompoundNucleus *GetNucleus(int);
  const CompoundNucleus *GetNucleus(int) const;
  /// 1-based nucleus of pair \p pairNum (1-based pair number).
  int NucleusOfPair(int);
  /// "Z=10 A=23" style label of nucleus \p i.
  std::string NucleusLabel(int) const;
  /// Did the project carry an <abundances> block?  Only then are the
  /// abundances Minuit parameters (one per nucleus, after the level block),
  /// so a project without the block keeps its parameter layout.
  bool HasAbundanceBlock() const { return hasAbundanceBlock_; }
  /// Number of abundance parameters FillMnParams emits (NumNuclei() with the block, else 0).
  int NumAbundanceParams() const { return hasAbundanceBlock_ ? (int)nuclei_.size() : 0; }
  /// Sum of the current abundances.
  double AbundanceSum() const;
  /// The abundance penalty added to chi-squared: sum over nuclei with a prior
  /// of ((w - w0)/sigma)^2, plus ((sum w) - 1)^2 / s^2 when \p sumSigma > 0.
  double AbundancePenalty(double sumSigma) const;
  /// Does a pair with this key exist? Keys are what the input files write; the
  /// PPair vector position may differ, so use GetPairNumFromKey to convert.
  bool IsPairKey(int);
  int NumPairs() const;
  int NumJGroups() const;
  /// 1-based position of the pair, or 0 if absent.
  int IsPair(PPair);
  /// 1-based position of the \f$J^\pi\f$ group, or 0 if absent.
  int IsJGroup(JGroup);
  /// Convert an input-file pair key to its 1-based position in the PPair vector.
  int GetPairNumFromKey(int);
  /// Pair \p i, 1-based.
  PPair *GetPair(int);
  /// \f$J^\pi\f$ group \p i, 1-based.
  JGroup *GetJGroup(int);
  void AddPair(PPair);
  void AddJGroup(JGroup);
  /// Largest orbital angular momentum appearing in any channel.
  int GetMaxLValue() const;
  void SetMaxLValue(int);

  // -- building -------------------------------------------------------------
  /// Build the whole object from the nuclear and external-capture input files.
  /// Returns -1 if a file cannot be read. \p radii optionally overrides one
  /// pair's channel radius, as SetRadius does.
  int Fill(const Config &, std::pair<int, double> radii = std::pair<int, double>(0, 0.0));
  /// Adds a level-less J^pi group, with its particle channels, for every J^pi the particle pairs allow up to
  /// \p maxL that no active level provides. Such a group contributes hard-sphere scattering only.
  int AddHardSphereGroups(const Config &, int maxL);
  /// Refuse (return -1) a channel with L+S odd on an identical pair with spin.
  int CheckIdenticalChannels(const Config &);
  /// Read the external-capture file into the ECLevel vector and check each
  /// final state against the nuclear file.
  void ParseExternalCapture(const Config &, std::map<int, int> &);
  /// Boundary conditions, the physical-to-formal transformation, the reaction
  /// pathways and the angular-distribution coefficients -- everything derived
  /// that a calculation needs. Run after Fill, and again after anything that
  /// changes the Coulomb functions.
  void Initialize(const Config &);

  // -- derived quantities ---------------------------------------------------
  /// Boundary condition per channel, evaluated at the energy of the group's
  /// first level.
  void CalcBoundaryConditions(const Config &);
  /// Build the internal (MGroup) and external (ECMGroup) reaction pathways.
  void SortPathways(const Config &);
  /// Build the KLGroup and Interference structure and its \f$Z_1Z_2\f$ coefficients.
  void CalcAngularDists(int);
  /// Build one decay's capture analyzing-power coefficient table (Seyler and
  /// Weller). Called on first use, not from Initialize.
  void CalcCaptureAnalyzingPower(int, int, int);
  /// Recompute the shift functions at the current level energies. Needed every
  /// iteration under the Brune parameterization, where they move with the fit.
  void CalcShiftFunctions(const Config &);
  /// Forget every level's cached shift functions (after a change of the
  /// nuclear potential, which moves them at fixed energy).
  void InvalidateShiftFunctions();
  /// Rescale the level amplitudes of a Minuit parameter set between Brune's
  /// and Park's normalizations (gamma_Park = gamma_Brune sqrt J), in place.
  /// `from`/`to` are AZUREParams::Basis values; false if the conversion is
  /// not level-wise (either side standard) or a level has J <= 0.
  bool ConvertAmplitudeBasis(ROOT::Minuit2::MnUserParameters &, int from, int to, const Config &);
  /// Park formalism: the penalty added to chi-squared for levels whose
  /// overlap J = 1 - sum gamma^2 dS/dE is not positive, sum (J/kParkNormScale)^2
  /// over J < 0.  Zero for every valid parameter set, and in Brune mode.
  double ParkNormPenalty();
  /// Scale of the Park J > 0 penalty: a level at J = -kParkNormScale costs one
  /// unit of chi-squared, so the wall is steep on the scale of J but smooth.
  static constexpr double kParkNormScale = 1.0e-3;
  /// External reduced width amplitude for one channel.
  complex CalcExternalWidth(JGroup *, ALevel *, AChannel *, bool, const Config &);

  // -- parameter transformations --------------------------------------------
  /// Observed widths and energies to formal R-matrix parameters. Returns false
  /// if a level could not be transformed.
  bool TransformIn(const Config &);
  /// Has TransformIn run on the current level scheme?  Until it has, the
  /// level gammas are the input values (partial widths in eV, ANCs, gamma
  /// widths in eV), not reduced width amplitudes -- consumers that need to
  /// estimate a resonance width before initialisation completes must check.
  bool IsTransformedIn() const { return transformedIn_; };
  /// Formal parameters back to observed energies and partial widths, for
  /// reporting. Under Brune the two are the same and this is a near no-op.
  void TransformOut(const Config &);
  /// The transformed physical parameters as a flat vector.
  vector_r GetTransformParams(const Config &configure);
  /// Write the transformed parameters to parameters.out.
  void PrintTransformParams(const Config &);

  // -- fitting --------------------------------------------------------------
  /// Seed the Minuit parameter array from the current level parameters.
  void FillMnParams(ROOT::Minuit2::MnUserParameters &, const Config *config = nullptr);
  /// Write a Minuit parameter vector back into the levels.
  void FillCompoundFromParams(const vector_r &);
  /// As above, for a vector already in physical rather than formal parameters.
  void FillCompoundFromParamsPhysical(const vector_r &);
  void PrintCompoundFromParams();
  /// Warn once if any level has a radiative width not small against its
  /// particle width, where the R-matrix treatment of capture breaks down.
  void CheckRadiativeWidths(const Config &, const vector_r &);

  // -- diagnostics ----------------------------------------------------------
  /// Print pairs, groups, levels and channels as read, before Initialize.
  void PrintNuc(const Config &);
  void PrintPathways(const Config &);
  void PrintBoundaryConditions(const Config &);
  void PrintAngularDists(const Config &);

  /// Deep copy. AZURECalc gives each thread its own, since a calculation
  /// writes the fitted parameters into the object it runs on.
  CNuc *Clone() const;

 private:
  /// 1-based index of the nucleus with this (Z, A), appending it if new.
  int FindOrAddNucleus(int z, int a);
  /// Reads the optional <abundances> block of the project file into nuclei_.
  /// Returns -1 on a malformed line.
  int ReadAbundances(const Config &);
  std::vector<PPair> pairs_;
  std::vector<JGroup> jgroups_;
  std::vector<CompoundNucleus> nuclei_;
  bool hasAbundanceBlock_ = false;
  int maxLValue_;
  bool transformedIn_ = false;
};

extern double DoubleFactorial(int);

#endif
