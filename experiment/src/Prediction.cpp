#include "azure2/experiment/Prediction.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>

namespace azure2 {
namespace experiment {

namespace {

/// The shortest decimal form of x that reads back to the same double (without an exponent between 1e-4 and 1e15).
std::string Number(double x) {
  char buf[40];
  const bool plain = std::fabs(x) >= 1e-4 && std::fabs(x) < 1e15;
  for (int p = 1; p <= 17; p++) {
    std::snprintf(buf, sizeof buf, "%.*g", p, x);
    if (std::strtod(buf, nullptr) == x && !(plain && std::strchr(buf, 'e'))) break;
  }
  return buf;
}

/// The leaf of term t: its observable's name, its pairs and its options (data_model.md section 3).
LeafSelector Leaf(const PredictionSetup &s, int t, LeafPart part) {
  const ReactionSpec &r = s.reaction;
  const PredictionTerm &term = s.terms[t];
  std::string kind, options;
  switch (r.observable) {
    case Observable::CrossSection:
    case Observable::CrossSectionE1:
    case Observable::CrossSectionE2:
    case Observable::TotalCapture: kind = term.total ? "sigma_total" : "sigma"; break;
    case Observable::DifferentialLab:  // lab data are compared in the c.m. frame (converted at load)
    case Observable::DifferentialCM: kind = "dsdo_cm"; break;
    case Observable::AnalyzingPower: kind = "ay"; break;
    case Observable::PhaseShift:
      kind = "phase";
      options += "; J=" + Number(r.phaseJ) + "; l=" + std::to_string(r.phaseL);
      if (r.phaseS >= 0.0) options += "; s=" + Number(r.phaseS);
      break;
    case Observable::AngularCoefficients: kind = s.absolute ? "Ak" : "ak"; break;
    case Observable::AngularDistribution:
      kind = "angdist";
      options += "; k=" + std::to_string(s.maxAngDistOrder);
      break;
  }
  if (r.thm) kind = "thm_" + kind;
  if (part == LeafPart::Rutherford) kind = "rutherford";
  if (part == LeafPart::Mott) kind = "mott";
  if (part == LeafPart::E1) options += "; mult=E1";
  if (part == LeafPart::E2) options += "; mult=E2";
  if (s.upos) options += "; upos L=" + std::to_string(s.uposL) + "; Ic=" + Number(s.uposIc) + "; delta=" + Number(s.uposDelta);
  if (term.hasAngle) options += "; angle=" + Number(term.angle);
  LeafSelector leaf;
  leaf.term = t;
  leaf.part = part;
  leaf.text = kind + "(" + std::to_string(term.entranceKey) + "->" + std::to_string(term.exitKey) + options + ")";
  return leaf;
}

/// Term 0's part: the E1 / E2 part for E1 / E2 data, the observable otherwise.  The other terms of a SUM or a RATIO
/// read the observable (ModelPass.h).
LeafPart BasePart(const PredictionSetup &s) {
  if (s.reaction.observable == Observable::CrossSectionE1) return LeafPart::E1;
  if (s.reaction.observable == Observable::CrossSectionE2) return LeafPart::E2;
  return LeafPart::Value;
}

int AddLeaf(const PredictionSetup &s, int t, LeafPart part, Prediction &p) {
  p.leaves.push_back(Leaf(s, t, part));
  return p.tree.Leaf(static_cast<int>(p.leaves.size()) - 1);
}

int AddFactor(FactorRef::Kind kind, int t, Prediction &p) {
  p.factors.push_back({kind, t});
  return static_cast<int>(p.factors.size()) - 1;
}

/// Coefficient data in an expansion whose coefficient is the engine's A_L / A_0 or A_L times a factor of L.
bool LegendreFactorOnData(const ReactionSpec &r) {
  return r.observable == Observable::AngularCoefficients && r.legendre != LegendreExpansion::AOverA0 &&
         r.legendre != LegendreExpansion::Absolute;
}

int AddAbundance(const PredictionSetup &s, int t, Prediction &p) {
  p.params.push_back({ParamRef::Kind::Abundance, s.terms[t].nucleus});
  return static_cast<int>(p.params.size()) - 1;
}

/// v a for a term of a sample (a = the abundance of the term's nucleus), v otherwise; the node's group is the
/// nucleus (0-based) for the decomposition.
int Weigh(const PredictionSetup &s, int t, int node, Prediction &p) {
  if (s.sample <= 0) return node;
  return p.tree.Weighted(AddAbundance(s, t, p), node, s.terms[t].nucleus - 1);
}

// ---------------------------------------------------------------- the builders, one per kind

void BuildObservable(const PredictionSetup &s, Prediction &p) { Weigh(s, 0, AddLeaf(s, 0, BasePart(s), p), p); }

/// sum(term 0, (c_k J_0/J_k) v_k, ...), each term weighted by its abundance in a sample: the order of Segment's SUM
/// (the base first, the components in order; a term is (c J) v, then times a).
void BuildSum(const PredictionSetup &s, Prediction &p) {
  std::vector<int> children{Weigh(s, 0, AddLeaf(s, 0, BasePart(s), p), p)};
  for (int t = 1; t < static_cast<int>(s.terms.size()); t++) {
    const int leaf = AddLeaf(s, t, LeafPart::Value, p);
    const int scaled = p.tree.Scale(s.terms[t].scaling, leaf, AddFactor(FactorRef::Kind::TermFrame, t, p));
    children.push_back(Weigh(s, t, scaled, p));
  }
  p.tree.Sum(children);
}

/// term 0 over term 1, each with its own lab-to-c.m. factor (num/den (J_den/J_num)); further terms are not read.
void BuildRatio(const PredictionSetup &s, Prediction &p) {
  if (s.terms.size() < 2) throw std::invalid_argument("experiment::BuildPrediction: a RATIO needs a second term");
  if (s.sample > 0) throw std::invalid_argument("experiment::BuildPrediction: a RATIO is not weighted by a sample");
  const int num = AddLeaf(s, 0, BasePart(s), p);
  const int den = AddLeaf(s, 1, LeafPart::Value, p);
  const int jNum = AddFactor(FactorRef::Kind::RatioFrame, 0, p);
  const int jDen = AddFactor(FactorRef::Kind::RatioFrame, 1, p);
  p.tree.Ratio(num, den, jNum, jDen);
}

/// term 0 over the Coulomb cross section of term 0's pair at its point: no frame factors (both are c.m. values at the
/// same c.m. angle, and the data are not converted at load).  A sample, components and a THM observable are refused.
void BuildCoulombRatio(const PredictionSetup &s, Prediction &p) {
  const ReactionSpec &r = s.reaction;
  if (s.terms.size() != 1 || s.sample > 0 || r.thm || r.combination != ReactionSpec::Combination::None)
    throw std::invalid_argument("experiment::BuildPrediction: a ratio to a Coulomb cross section is one elastic term");
  if (r.observable != Observable::DifferentialLab && r.observable != Observable::DifferentialCM)
    throw std::invalid_argument("experiment::BuildPrediction: a ratio to a Coulomb cross section divides dsigma/dOmega");
  const int num = AddLeaf(s, 0, LeafPart::Value, p);
  const int den =
      AddLeaf(s, 0, r.coulombReference == CoulombReference::Mott ? LeafPart::Mott : LeafPart::Rutherford, p);
  p.tree.Ratio(num, den);
}

/// The weighted terms of a transmission or a yield (Segment::TransformTerms): term 0 with weight a, a component
/// with (c J_0/J_k) a, the yield's total term with 1.
void BuildTransform(const PredictionSetup &s, Prediction &p, TransformSpec::Kind kind) {
  TransformSpec spec;
  spec.kind = kind;
  p.params.push_back({ParamRef::Kind::Thickness, 0});
  spec.thickness = static_cast<int>(p.params.size()) - 1;
  spec.fold = 0;  // the base unit's kernel, supplied by the engine at the point
  bool hasTotal = false;
  for (int t = 0; t < static_cast<int>(s.terms.size()); t++) {
    TransformTerm term;
    p.leaves.push_back(Leaf(s, t, t == 0 ? BasePart(s) : LeafPart::Value));
    term.input = static_cast<int>(p.leaves.size()) - 1;
    if (s.terms[t].yieldTotal) {
      term.total = true;
      hasTotal = true;
    } else {
      term.c = (t == 0) ? 1.0 : s.terms[t].scaling;
      if (t > 0) term.factor = AddFactor(FactorRef::Kind::TermFrame, t, p);
      if (s.sample > 0) {
        term.abundance = AddAbundance(s, t, p);
        term.group = s.terms[t].nucleus - 1;
      }
    }
    spec.terms.push_back(term);
  }
  if (kind == TransformSpec::Kind::Yield && !hasTotal)
    throw std::invalid_argument("experiment::BuildPrediction: a yield needs its total-cross-section term");
  p.tree.Transform(spec);
}

/// c [x E exp(2 pi eta)] x the tree the setup has without its factor (ReactionSpec::scale): the S factor of a cross
/// section, a yield per unit charge, an extrapolation in another unit than the engine's.
void BuildScaled(const PredictionSetup &s, Prediction &p) {
  PredictionSetup inner = s;
  inner.reaction.scale = ValueScale();
  const bool legendre = LegendreFactorOnData(s.reaction);
  if (legendre)
    inner.reaction.legendre =
        LegendreExpansionNormalised(s.reaction.legendre) ? LegendreExpansion::AOverA0 : LegendreExpansion::Absolute;
  TreeKinds()[static_cast<std::size_t>(TreeKindOf(inner))].build(inner, p);
  const ValueScale &scale = s.reaction.scale;
  int factor = -1;
  if (scale.sFactor) {
    factor = AddFactor(FactorRef::Kind::SFactor, 0, p);
    p.factors[factor].shape_perMeV = scale.sFactorShape_perMeV;
  } else if (legendre) {
    factor = AddFactor(FactorRef::Kind::Legendre, 0, p);
    p.factors[factor].legendre = s.reaction.legendre;
  }
  p.tree.Scale(scale.factor, p.tree.Root(), factor);
}

void BuildTransmission(const PredictionSetup &s, Prediction &p) { BuildTransform(s, p, TransformSpec::Kind::Transmission); }
void BuildYield(const PredictionSetup &s, Prediction &p) { BuildTransform(s, p, TransformSpec::Kind::Yield); }

// ---------------------------------------------------------------- text

std::string FactorText(const Prediction &p, int f) {
  if (f < 0) return "";
  const FactorRef &r = p.factors[f];
  if (r.kind == FactorRef::Kind::Legendre)
    return std::string("c_") + LegendreExpansionName(r.legendre) + "(" + std::to_string(r.term) + ")";
  if (r.kind == FactorRef::Kind::SFactor && r.shape_perMeV != 0.0)
    return "S*(" + std::to_string(r.term) + ";g=" + Number(r.shape_perMeV) + "/MeV)";
  if (r.kind == FactorRef::Kind::SFactor) return "S(" + std::to_string(r.term) + ")";
  return r.kind == FactorRef::Kind::TermFrame ? "J0/J" + std::to_string(r.term) : "J" + std::to_string(r.term);
}

std::string AbundanceText(const Prediction &p, int param) {
  return "abundance(" + p.sampleName + ", " + std::to_string(p.params[param].nucleus) + ")*";
}

/// c*J*: the scaling (omitted when 1) and the frame factor (when there is one).
std::string ScaleText(const Prediction &p, double c, int factor) {
  std::string t = (c == 1.0 && factor >= 0) ? "" : Number(c) + "*";
  if (factor >= 0) t += FactorText(p, factor) + "*";
  return t;
}

std::string NodeText(const Prediction &p, int k) {
  const Node &nd = p.tree.Nodes()[k];
  const int *ch = p.tree.Children().data() + nd.first;
  switch (nd.kind) {
    case NodeKind::Leaf: return p.leaves[nd.input].text;
    case NodeKind::Sum: {
      std::string t = "sum(";
      for (int m = 0; m < nd.count; m++) t += (m ? ", " : "") + NodeText(p, ch[m]);
      return t + ")";
    }
    case NodeKind::Scale: return ScaleText(p, nd.c, nd.factor) + NodeText(p, ch[0]);
    case NodeKind::Weighted: return AbundanceText(p, nd.param) + NodeText(p, ch[0]);
    case NodeKind::Ratio: {
      std::string t = "ratio(" + NodeText(p, ch[0]) + ", " + NodeText(p, ch[1]) + ")";
      if (nd.factor >= 0) t += "*" + FactorText(p, nd.factor2) + "/" + FactorText(p, nd.factor);
      return t;
    }
    case NodeKind::Transform: {
      const TransformSpec &spec = p.tree.Transforms()[nd.input];
      std::string sum, total;
      int n = 0;
      for (const TransformTerm &term : spec.terms) {
        if (term.total) {
          total = p.leaves[term.input].text;
          continue;
        }
        std::string t = (term.abundance >= 0) ? AbundanceText(p, term.abundance) : "";
        t += (term.factor >= 0 || term.c != 1.0) ? ScaleText(p, term.c, term.factor) : "";
        sum += (n++ ? ", " : "") + t + p.leaves[term.input].text;
      }
      if (n > 1) sum = "sum(" + sum + ")";
      const std::string head = (spec.kind == TransformSpec::Kind::Yield ? "yield(n=" : "transmission(n=") + Number(p.thickness);
      return head + ", " + sum + (total.empty() ? "" : ", " + total) + ")";
    }
  }
  return "";
}

}  // namespace

const std::vector<TreeKindRow> &TreeKinds() {
  static const std::vector<TreeKindRow> rows = {
      {TreeKind::Observable, "observable", BuildObservable},
      {TreeKind::Sum, "sum", BuildSum},
      {TreeKind::Ratio, "ratio", BuildRatio},
      {TreeKind::Transmission, "transmission", BuildTransmission},
      {TreeKind::Yield, "yield", BuildYield},
      {TreeKind::CoulombRatio, "coulomb_ratio", BuildCoulombRatio},
      {TreeKind::Scaled, "scaled", BuildScaled},
  };
  return rows;
}

TreeKind TreeKindOf(const PredictionSetup &s) {
  const ReactionSpec &r = s.reaction;
  if (!r.scale.Identity() || LegendreFactorOnData(r)) return TreeKind::Scaled;
  if (r.transmission > 0.0) return TreeKind::Transmission;
  if (r.yieldThickness > 0.0) return TreeKind::Yield;
  if (r.combination == ReactionSpec::Combination::Ratio) return TreeKind::Ratio;
  if (r.combination == ReactionSpec::Combination::Sum || r.observable == Observable::TotalCapture) return TreeKind::Sum;
  if (r.coulombReference != CoulombReference::None) return TreeKind::CoulombRatio;
  return TreeKind::Observable;
}

Prediction BuildPrediction(const PredictionSetup &setup) {
  if (setup.terms.empty()) throw std::invalid_argument("experiment::BuildPrediction: no terms");
  Prediction p;
  p.kind = TreeKindOf(setup);
  p.sampleName = setup.sampleName;
  p.thickness = setup.thickness;
  TreeKinds()[static_cast<std::size_t>(p.kind)].build(setup, p);
  return p;
}

std::string PredictionText(const Prediction &prediction) {
  if (prediction.tree.Nodes().empty()) return "";
  return NodeText(prediction, prediction.tree.Root());
}

}  // namespace experiment
}  // namespace azure2
