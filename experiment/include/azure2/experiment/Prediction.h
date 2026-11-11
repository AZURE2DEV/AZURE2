#ifndef AZURE2_EXPERIMENT_PREDICTION_H
#define AZURE2_EXPERIMENT_PREDICTION_H

// What a data set or an extrapolation predicts, as an expression tree (docs/upgrade/data_model.md sections 2.3 and 3,
// step S10): built once at set-up from the data set's reaction and the terms the set-up made, then evaluated at every
// point of every evaluation by the engine's model pass (src/ModelPass.cpp).  The tree holds structure and constants
// only; what a leaf, a factor or a parameter refers to is listed beside it (leaves, factors, params), and the engine
// supplies their values at a point.  One builder per kind of tree, chosen in one table (TreeKinds).

#include <string>
#include <vector>

#include "azure2/experiment/Expression.h"
#include "azure2/experiment/SetupSpecs.h"

namespace azure2 {
namespace experiment {

/// Which value of a term's calculation unit a leaf reads at a point: the observable itself, or the E1 / E2 part of an
/// angle-integrated capture cross section (b) -- the codes 0, 1, 2 are the engine's cross-section component --, or the
/// Rutherford or Mott cross section (b/sr, c.m.) at the term's point, folded as the point's observable is: an analytic
/// leaf, which no level parameter moves (src/CoulombReference.cpp).
enum class LeafPart { Value = 0, E1 = 1, E2 = 2, Rutherford = 3, Mott = 4 };

/// Whether a leaf is the analytic Coulomb reference (no R-matrix calculation, no adjoint to the level parameters).
inline bool IsCoulombReference(LeafPart part) { return part == LeafPart::Rutherford || part == LeafPart::Mott; }

/// One term of a prediction as the set-up resolved it: term 0 is the data set's own points, terms 1, 2, ... its
/// component segments in order (the components of a SUM or RATIO, the further capture pairs of a total capture, the
/// total-cross-section term of a yield last).
struct PredictionTerm {
  int entranceKey = 0;      ///< pair key of <levels>
  int exitKey = 0;          ///< pair key of <levels> (resolved: a total capture's term 0 holds its first capture pair)
  double scaling = 1.0;     ///< the component's factor c (1 for term 0)
  bool hasAngle = false;    ///< the term is evaluated at a fixed angle, not at the data set's
  double angle = 0.0;       ///< that lab angle (deg)
  int nucleus = 0;          ///< the compound nucleus (1-based) of the entrance pair: the index of its abundance
  bool total = false;       ///< the term is the total cross section of its entrance pair (optical theorem)
  bool yieldTotal = false;  ///< the total-cross-section term s_t of a yield
};

/// What the tree of one data set or extrapolation is built from.
struct PredictionSetup {
  ReactionSpec reaction;     ///< the observable and the model options as read (pairs, components: see `terms`)
  bool absolute = false;     ///< AngularCoefficients: the data are A_k, not a_k
  int maxAngDistOrder = 0;   ///< AngularDistribution: highest Legendre order
  bool upos = false;         ///< unobserved primary, observed secondary gamma (UPOS)
  int uposL = 0;
  double uposIc = 0.0, uposDelta = 0.0;
  std::vector<PredictionTerm> terms;  ///< term 0 first (never empty)
  int sample = 0;            ///< the 1-based sample of the abundances weighting the terms, 0 for none
  std::string sampleName;    ///< its name as written ("1" for the default sample), for the text form
  double thickness = 0.0;    ///< transmission / yield: the nominal sample thickness n (atoms/b), for the text form
};

/// What a leaf input of the tree reads: a term's unit, which value, and its text form.
struct LeafSelector {
  int term = 0;
  LeafPart part = LeafPart::Value;
  std::string text;  ///< e.g. "sigma(1->2; mult=E1)" (data_model.md section 3)
};

/// What a per-point factor array of the tree holds at a point.
struct FactorRef {
  enum class Kind {
    TermFrame,   ///< J_0 / J_k of a SUM or transform term k: the term into the frame of term 0 (Segment::TermFrameFactor)
    RatioFrame,  ///< J_k of a RATIO term k, its own lab-to-c.m. factor (Segment::RatioTermFrameFactor)
    SFactor,     ///< E exp(2 pi eta(E) + g E) of term k's entrance pair at its point's c.m. energy E (MeV; S = this x
                 ///< sigma, physics/AstrophysicalFactor.h; g = shape_perMeV, 0 for the S factor)
    Legendre,    ///< LegendreExpansionFactor(legendre, L) of term k's point, L its Legendre order (coefficient data)
  };
  Kind kind = Kind::TermFrame;
  int term = 0;
  double shape_perMeV = 0.0;  ///< SFactor: g of the modified S factor (1/MeV)
  LegendreExpansion legendre = LegendreExpansion::AOverA0;  ///< Legendre: the expansion
};

/// What a parameter of the tree is.
struct ParamRef {
  enum class Kind {
    Abundance,  ///< the abundance of compound nucleus `nucleus` in the setup's sample
    Thickness,  ///< the sample thickness n (atoms/b) of a transmission or yield
  };
  Kind kind = Kind::Abundance;
  int nucleus = 0;
};

/// The kinds of tree, one builder each (TreeKinds).
enum class TreeKind {
  Observable,    ///< one leaf (weighted by its abundance in a sample)
  Sum,           ///< term 0 plus the scaled terms in its frame (SUM, total capture)
  Ratio,         ///< term 0 over term 1, each with its own frame factor
  Transmission,  ///< T = exp(-n s) of the weighted terms, at every sub-point before the kernel
  Yield,         ///< Y0 = (s / s_t)(1 - exp(-n s_t)), s_t the yield's total term
  CoulombRatio,  ///< dsigma/dOmega over the Rutherford or Mott cross section at the same point, both folded alike
  Scaled,        ///< c [x E exp(2 pi eta)] x the tree of another kind (ReactionSpec::scale: S factor, yield per charge)
};

/// A built prediction: the tree and what its inputs, factors and parameters refer to (indices into these lists).
struct Prediction {
  TreeKind kind = TreeKind::Observable;
  Tree tree;
  std::vector<LeafSelector> leaves;
  std::vector<FactorRef> factors;
  std::vector<ParamRef> params;
  std::string sampleName;  ///< for the text form of the abundances
  double thickness = 0.0;  ///< for the text form of a transform
};

/// One row of the table of tree kinds.
struct TreeKindRow {
  TreeKind kind;
  const char *name;
  void (*build)(const PredictionSetup &, Prediction &);
};

/// The table: one row per TreeKind, in the order of the enum.
const std::vector<TreeKindRow> &TreeKinds();

/// The kind of tree a setup needs (the one place that decides it): a factor on the model value (ReactionSpec::scale,
/// or the factor of a Legendre expansion other than a_over_a0 and absolute on coefficient data) first, whose builder wraps the tree of the kind the setup has without it; then a transmission or a yield by its thickness, else a
/// RATIO, else a SUM (a SUM line or a total capture), else a ratio to a Coulomb cross section, else the observable
/// alone.
TreeKind TreeKindOf(const PredictionSetup &setup);

/// The tree of a setup, by its kind's builder.  Throws std::invalid_argument for a setup no builder takes (no terms, a
/// RATIO without a second term, a yield without its total term).
Prediction BuildPrediction(const PredictionSetup &setup);

/// The tree as one line of text (data_model.md section 3): leaves "sigma(1->2)", "dsdo_cm(1->2)", ...; "sum(a, b)",
/// "c*J0/J1*a" (a SUM term: its factor and the frame factor J_0/J_1), "abundance(<sample>, <nucleus>)*a",
/// "ratio(a, b)*J1/J0", "transmission(n=<n>, a)", "yield(n=<n>, a, b)"; the Coulomb leaves "rutherford(1->1)",
/// "mott(1->1)"; a factor on the value "c*S(1)*a" (c omitted when 1, S(k) the S-factor conversion of term k; "S*(k;
/// g=<g>/MeV)" the modified S factor's).
std::string PredictionText(const Prediction &prediction);

}  // namespace experiment
}  // namespace azure2

#endif
