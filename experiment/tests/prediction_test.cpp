// The tree builders of experiment/Prediction.h (docs/upgrade/data_model.md section 3, step S10): for every kind of
// tree, the tree chosen, its text form, what its inputs refer to, and its value at random inputs against the engine's
// combination before S10 written out (src/ModelPass.cpp at cdf58d41: ComposeModelValue, TransformedValue), bit for bit.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "azure2/experiment/Prediction.h"

using namespace azure2::experiment;

namespace {

int failures = 0;

void check(const char *what, bool ok, const std::string &detail = "") {
  std::printf("  %s  %s%s%s\n", ok ? "ok  " : "FAIL", what, detail.empty() ? "" : ": ", detail.c_str());
  if (!ok) failures++;
}

std::mt19937_64 rng(20261110);
double rnd(double lo, double hi) { return std::uniform_real_distribution<double>(lo, hi)(rng); }
bool same(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

PredictionTerm Term(int entrance, int exit, int nucleus, double scaling = 1.0) {
  PredictionTerm t;
  t.entranceKey = entrance;
  t.exitKey = exit;
  t.nucleus = nucleus;
  t.scaling = scaling;
  return t;
}

/// The values the engine supplies at one point: leaf k, factor f, parameter q, by what the prediction says they are.
struct PointValues {
  std::vector<double> leaf, factor, param;
  std::vector<PointArray> leafArray, factorArray;
  std::vector<const SubPointFold *> folds{nullptr};  // a transform's kernel: none, the point values themselves
  Inputs in;
  void Fill(const Prediction &p, const std::vector<double> &termValue, const std::vector<double> &termJ,
            const std::vector<double> &abundance, double thickness) {
    leaf.clear();
    factor.clear();
    param.clear();
    for (const LeafSelector &l : p.leaves) leaf.push_back(termValue[l.term]);
    for (const FactorRef &f : p.factors)
      factor.push_back(f.kind == FactorRef::Kind::TermFrame ? TermFrameFactor(termJ[0], termJ[f.term]) : termJ[f.term]);
    for (const ParamRef &q : p.params) param.push_back(q.kind == ParamRef::Kind::Abundance ? abundance[q.nucleus] : thickness);
    leafArray.clear();
    factorArray.clear();
    for (double &v : leaf) leafArray.push_back({&v, 1});
    for (double &v : factor) factorArray.push_back({&v, 1});
    in.leaves = leafArray;
    in.factors = factorArray;
    in.params = param.data();
    in.folds = folds;
  }
};

void TestObservable() {
  std::printf("observable\n");
  PredictionSetup s;
  s.reaction.observable = Observable::CrossSectionE1;
  s.terms = {Term(1, 2, 1)};
  Prediction p = BuildPrediction(s);
  check("kind observable", p.kind == TreeKind::Observable && TreeKindOf(s) == TreeKind::Observable);
  check("text", PredictionText(p) == "sigma(1->2; mult=E1)", PredictionText(p));
  check("E1 leaf", p.leaves.size() == 1 && p.leaves[0].part == LeafPart::E1);

  s.reaction.observable = Observable::DifferentialLab;
  s.reaction.thm = true;
  check("THM lab differential", PredictionText(BuildPrediction(s)) == "thm_dsdo_cm(1->2)", PredictionText(BuildPrediction(s)));
  s.reaction.thm = false;
  s.reaction.observable = Observable::PhaseShift;
  s.reaction.phaseJ = 1.5;
  s.reaction.phaseL = 1;
  s.reaction.phaseS = 0.5;
  s.terms = {Term(1, 1, 1)};
  check("phase shift", PredictionText(BuildPrediction(s)) == "phase(1->1; J=1.5; l=1; s=0.5)", PredictionText(BuildPrediction(s)));
  s.reaction.observable = Observable::AngularCoefficients;
  s.absolute = true;
  check("A_k", PredictionText(BuildPrediction(s)) == "Ak(1->1)");
  s.reaction.observable = Observable::AngularDistribution;
  s.maxAngDistOrder = 4;
  check("angular distribution", PredictionText(BuildPrediction(s)) == "angdist(1->1; k=4)");
  s.reaction.observable = Observable::CrossSection;
  s.terms[0].total = true;
  check("total", PredictionText(BuildPrediction(s)) == "sigma_total(1->1)");

  // a sample: v a
  s = PredictionSetup();
  s.terms = {Term(3, 4, 2)};
  s.sample = 1;
  s.sampleName = "natNe";
  p = BuildPrediction(s);
  check("sample text", PredictionText(p) == "abundance(natNe, 2)*sigma(3->4)", PredictionText(p));
  bool ok = true, okNat = true;
  PointValues pv;
  for (int i = 0; i < 2000; i++) {
    const double v = rnd(-5, 5), a = rnd(0, 1);
    pv.Fill(p, {v}, {1.0}, {0.0, 0.0, a}, 0.0);
    Workspace ws;
    double nat[2] = {0.0, 0.0};
    const double m = p.tree.Value(pv.in, 0, ws, nat);
    double baseValue = v;
    baseValue *= a;  // ComposeModelValue
    ok = ok && same(m, baseValue);
    okNat = okNat && same(nat[1], 0.0 + baseValue) && nat[0] == 0.0;
  }
  check("sample value = v a, bit for bit", ok);
  check("sample decomposition", okNat);
}

void TestSum() {
  std::printf("sum\n");
  PredictionSetup s;
  s.reaction.observable = Observable::DifferentialLab;
  s.reaction.combination = ReactionSpec::Combination::Sum;
  s.terms = {Term(1, 2, 1), Term(1, 3, 1, 0.5), Term(4, 5, 2, 2.0)};
  s.terms[2].hasAngle = true;
  s.terms[2].angle = 30.0;
  Prediction p = BuildPrediction(s);
  check("kind sum", p.kind == TreeKind::Sum);
  check("text", PredictionText(p) == "sum(dsdo_cm(1->2), 0.5*J0/J1*dsdo_cm(1->3), 2*J0/J2*dsdo_cm(4->5; angle=30))",
        PredictionText(p));
  // total capture: a SUM of the capture pairs, without a SUM line
  PredictionSetup tc;
  tc.reaction.observable = Observable::TotalCapture;
  tc.terms = {Term(1, 2, 1), Term(1, 3, 1)};
  check("total capture is a sum", TreeKindOf(tc) == TreeKind::Sum &&
                                      PredictionText(BuildPrediction(tc)) == "sum(sigma(1->2), J0/J1*sigma(1->3))",
        PredictionText(BuildPrediction(tc)));
  tc.terms.resize(1);
  check("total capture of one pair", PredictionText(BuildPrediction(tc)) == "sum(sigma(1->2))");

  for (int natural = 0; natural <= 1; natural++) {
    s.sample = natural;
    s.sampleName = "1";
    p = BuildPrediction(s);
    bool ok = true, okNat = true;
    PointValues pv;
    for (int i = 0; i < 4000; i++) {
      const std::vector<double> v{rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)}, J{rnd(0.5, 2), rnd(0.5, 2), rnd(0.5, 2)},
          a{0.0, rnd(0, 1), rnd(0, 1)};
      pv.Fill(p, v, J, a, 0.0);
      Workspace ws;
      double nat[2] = {0.0, 0.0};
      const double m = p.tree.Value(pv.in, 0, ws, natural ? nat : nullptr);
      // ComposeModelValue's SUM
      double natOld[2] = {0.0, 0.0};
      double baseValue = v[0];
      if (natural) {
        baseValue *= a[1];
        natOld[0] += baseValue;
      }
      double sum = baseValue;
      for (int k = 1; k <= 2; k++) {
        double term = s.terms[k].scaling * TermFrameFactor(J[0], J[k]) * v[k];
        if (natural) {
          term *= a[s.terms[k].nucleus];
          natOld[s.terms[k].nucleus - 1] += term;
        }
        sum += term;
      }
      ok = ok && same(m, sum);
      okNat = okNat && same(nat[0], natOld[0]) && same(nat[1], natOld[1]);
    }
    check(natural ? "sample SUM = Segment's, bit for bit" : "SUM = Segment's, bit for bit", ok);
    if (natural) check("sample SUM decomposition", okNat);
  }
  check("sample SUM text", PredictionText(p) == "sum(abundance(1, 1)*dsdo_cm(1->2), abundance(1, 1)*0.5*J0/J1*dsdo_cm(1->3), "
                                                "abundance(1, 2)*2*J0/J2*dsdo_cm(4->5; angle=30))",
        PredictionText(p));
}

void TestRatio() {
  std::printf("ratio\n");
  PredictionSetup s;
  s.reaction.observable = Observable::CrossSectionE2;
  s.reaction.combination = ReactionSpec::Combination::Ratio;
  s.terms = {Term(1, 2, 1), Term(1, 3, 1), Term(1, 4, 1)};
  Prediction p = BuildPrediction(s);
  check("kind ratio", p.kind == TreeKind::Ratio);
  check("text (a third term is not read)", PredictionText(p) == "ratio(sigma(1->2; mult=E2), sigma(1->3))*J1/J0",
        PredictionText(p));
  check("numerator E2, denominator the observable", p.leaves[0].part == LeafPart::E2 && p.leaves[1].part == LeafPart::Value);
  bool ok = true;
  PointValues pv;
  for (int i = 0; i < 4000; i++) {
    const std::vector<double> v{rnd(-3, 3), (i % 100 == 0) ? 0.0 : rnd(-3, 3), 7.0}, J{rnd(0.5, 2), rnd(0.5, 2), 1.0};
    pv.Fill(p, v, J, {}, 0.0);
    Workspace ws;
    const double m = p.tree.Value(pv.in, 0, ws);
    const double ref = (v[1] == 0.0) ? NAN : v[0] / v[1] * (J[1] / J[0]);  // ComposeModelValue's RATIO
    ok = ok && (std::isnan(ref) ? std::isnan(m) : same(m, ref));
  }
  check("RATIO = Segment's, bit for bit (NaN at a zero denominator)", ok);
  s.terms.resize(1);
  bool threw = false;
  try {
    BuildPrediction(s);
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  check("a RATIO without a second term is refused", threw);
}

void TestTransforms() {
  std::printf("transmission and yield\n");
  // transmission of a sample with a component: the point-value path (no sub-point stage)
  PredictionSetup s;
  s.reaction.transmission = 0.02;
  s.reaction.combination = ReactionSpec::Combination::Sum;
  s.terms = {Term(1, 1, 1), Term(5, 5, 2, 0.5)};
  s.terms[0].total = s.terms[1].total = true;
  s.sample = 2;
  s.sampleName = "mix";
  s.thickness = 0.02;
  Prediction p = BuildPrediction(s);
  check("kind transmission", p.kind == TreeKind::Transmission);
  check("text", PredictionText(p) == "transmission(n=0.02, sum(abundance(mix, 1)*sigma_total(1->1), "
                                     "abundance(mix, 2)*0.5*J0/J1*sigma_total(5->5)))",
        PredictionText(p));
  bool ok = true, okNat = true;
  PointValues pv;
  for (int i = 0; i < 2000; i++) {
    const std::vector<double> v{rnd(0, 30), rnd(0, 30)}, J{rnd(0.5, 2), rnd(0.5, 2)}, a{0.0, rnd(0, 1), rnd(0, 1)};
    const double n = rnd(0, 0.1);
    pv.Fill(p, v, J, a, n);
    Workspace ws;
    double nat[2] = {0.0, 0.0};
    const double m = p.tree.Value(pv.in, 0, ws, nat);
    // Segment::TransformTerms and TransformedValue's `transform`
    const double w0 = a[1];
    const double scale1 = s.terms[1].scaling * TermFrameFactor(J[0], J[1]);
    const double w1 = scale1 * a[2];
    double sum = 0.0;
    sum += w0 * v[0];
    sum += w1 * v[1];
    ok = ok && same(m, node::Transmission(n, sum));
    okNat = okNat && same(nat[0], 0.0 + w0 * v[0]) && same(nat[1], 0.0 + w1 * v[1]);
  }
  check("transmission = TransformedValue, bit for bit", ok);
  check("transmission decomposition", okNat);

  // a yield: the total term last, weight 1
  PredictionSetup y;
  y.reaction.yieldThickness = 0.01;
  y.terms = {Term(1, 2, 1), Term(1, 1, 1)};
  y.terms[1].total = y.terms[1].yieldTotal = true;
  y.thickness = 0.01;
  p = BuildPrediction(y);
  check("kind yield", p.kind == TreeKind::Yield);
  check("text", PredictionText(p) == "yield(n=0.01, sigma(1->2), sigma_total(1->1))", PredictionText(p));
  ok = true;
  for (int i = 0; i < 2000; i++) {
    const std::vector<double> v{rnd(0, 3), rnd(0, 30)};
    const double n = rnd(0, 0.1);
    pv.Fill(p, v, {1.0, 1.0}, {}, n);
    Workspace ws;
    double sum = 0.0, total = 0.0;
    sum += 1.0 * v[0];
    total += 1.0 * v[1];
    ok = ok && same(p.tree.Value(pv.in, 0, ws), node::YieldY0(n, sum, total));
  }
  check("yield = TransformedValue, bit for bit", ok);
  y.terms.resize(1);
  bool threw = false;
  try {
    BuildPrediction(y);
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  check("a yield without its total term is refused", threw);
}

/// A ratio to the Rutherford or Mott cross section: dsdo_cm over the analytic Coulomb leaf of the same point, no frame
/// factors (the value is not converted at load), and its adjoint g/den, -g num/den^2.
void TestCoulombRatio() {
  std::printf("ratio to a Coulomb cross section\n");
  PredictionSetup s;
  s.reaction.observable = Observable::DifferentialLab;
  s.reaction.coulombReference = CoulombReference::Rutherford;
  s.terms = {Term(1, 1, 1)};
  Prediction p = BuildPrediction(s);
  check("kind coulomb_ratio", p.kind == TreeKind::CoulombRatio && TreeKindOf(s) == TreeKind::CoulombRatio);
  check("text", PredictionText(p) == "ratio(dsdo_cm(1->1), rutherford(1->1))", PredictionText(p));
  check("leaves: the observable and the Rutherford leaf of term 0, no factors, no parameters",
        p.leaves.size() == 2 && p.leaves[0].part == LeafPart::Value && p.leaves[1].part == LeafPart::Rutherford &&
            p.leaves[0].term == 0 && p.leaves[1].term == 0 && IsCoulombReference(p.leaves[1].part) &&
            !IsCoulombReference(p.leaves[0].part) && p.factors.empty() && p.params.empty());
  bool ok = true, okBack = true;
  for (int i = 0; i < 2000; i++) {
    double v[2] = {rnd(0, 30), rnd(1e-3, 30)};
    PointArray arrays[2] = {{&v[0], 1}, {&v[1], 1}};
    Inputs in;
    in.leaves = View<PointArray>(arrays, 2);
    Workspace ws;
    ok = ok && same(p.tree.Value(in, 0, ws), v[0] / v[1]);
    const double g = rnd(-2, 2);
    double cn = 0.0, cd = 0.0;
    double *cots[2] = {&cn, &cd};
    Cotangents cot;
    cot.leaves = View<double *>(cots, 2);
    p.tree.Backward(in, 1, &g, cot, ws);
    okBack = okBack && std::fabs(cn - g / v[1]) <= 1e-14 * std::fabs(g / v[1]) &&
             std::fabs(cd + g * v[0] / (v[1] * v[1])) <= 1e-13 * std::fabs(g * v[0] / (v[1] * v[1]));
  }
  check("value = num / den, bit for bit", ok);
  check("backward: g / den and -g num / den^2", okBack);
  s.reaction.observable = Observable::DifferentialCM;
  s.reaction.coulombReference = CoulombReference::Mott;
  check("Mott, c.m. angles", PredictionText(BuildPrediction(s)) == "ratio(dsdo_cm(1->1), mott(1->1))",
        PredictionText(BuildPrediction(s)));
  auto refused = [](PredictionSetup x) {
    try {
      BuildPrediction(x);
    } catch (const std::invalid_argument &) {
      return true;
    }
    return false;
  };
  PredictionSetup x = s;
  x.sample = 1;
  check("refused: a sample", refused(x));
  x = s;
  x.reaction.thm = true;
  check("refused: a THM observable", refused(x));
  x = s;
  x.terms.push_back(Term(1, 2, 1));
  check("refused: a second term", refused(x));
  x = s;
  x.reaction.observable = Observable::CrossSection;
  check("refused: an angle-integrated observable", refused(x));
}

/// A factor on the model value (ReactionSpec::scale): c x [S(0) x] the tree of the setup without it; the S-factor
/// factor is a per-point factor of term 0, so the value is c S sigma and the leaf's cotangent g c S.
void TestScaled() {
  std::printf("scaled (S factor, yield per charge)\n");
  PredictionSetup s;
  s.reaction.observable = Observable::CrossSection;
  s.terms = {Term(1, 2, 1)};
  s.reaction.scale.sFactor = true;
  Prediction p = BuildPrediction(s);
  check("kind scaled", p.kind == TreeKind::Scaled && TreeKindOf(s) == TreeKind::Scaled);
  check("text", PredictionText(p) == "S(0)*sigma(1->2)", PredictionText(p));
  check("one leaf, one S-factor factor of term 0", p.leaves.size() == 1 && p.factors.size() == 1 &&
                                                       p.factors[0].kind == FactorRef::Kind::SFactor && p.factors[0].term == 0);
  bool ok = true, okBack = true;
  for (int i = 0; i < 2000; i++) {
    double v = rnd(0, 30), S = rnd(1e-3, 1e6);
    PointArray leaf = {&v, 1}, factor = {&S, 1};
    Inputs in;
    in.leaves = View<PointArray>(&leaf, 1);
    in.factors = View<PointArray>(&factor, 1);
    Workspace ws;
    ok = ok && same(p.tree.Value(in, 0, ws), 1.0 * S * v);
    const double g = rnd(-2, 2);
    double c = 0.0;
    double *cots[1] = {&c};
    Cotangents cot;
    cot.leaves = View<double *>(cots, 1);
    p.tree.Backward(in, 1, &g, cot, ws);
    okBack = okBack && same(c, g * (1.0 * S));
  }
  check("value = S sigma, bit for bit", ok);
  check("backward: g S", okBack);
  s.reaction.scale.sFactorShape_perMeV = 0.46;
  const Prediction pm = BuildPrediction(s);
  check("modified S factor: text S*(0;g=0.46/MeV)", PredictionText(pm) == "S*(0;g=0.46/MeV)*sigma(1->2)", PredictionText(pm));
  check("modified S factor: the factor of term 0 carries g", pm.factors.size() == 1 &&
                                                                pm.factors[0].kind == FactorRef::Kind::SFactor &&
                                                                pm.factors[0].shape_perMeV == 0.46);
  check("S factor: g = 0 in its factor", p.factors[0].shape_perMeV == 0.0);
  s.reaction.scale.sFactorShape_perMeV = 0.0;
  s.reaction.scale.sFactor = false;
  {
    // Legendre coefficient data: rsl, 2l2, 4pi are a factor of the point's order on the a_k or A_k leaf
    PredictionSetup c;
    c.reaction.observable = Observable::AngularCoefficients;
    c.terms = {Term(1, 2, 1)};
    c.reaction.legendre = LegendreExpansion::Rsl;
    const Prediction pr = BuildPrediction(c);
    check("legendre rsl: c_rsl(0)*ak", pr.kind == TreeKind::Scaled && PredictionText(pr) == "c_rsl(0)*ak(1->2)", PredictionText(pr));
    check("legendre rsl: one Legendre factor", pr.factors.size() == 1 && pr.factors[0].kind == FactorRef::Kind::Legendre &&
                                                    pr.factors[0].legendre == LegendreExpansion::Rsl);
    c.absolute = true;
    c.reaction.legendre = LegendreExpansion::FourPi;
    check("legendre 4pi: c_4pi(0)*Ak", PredictionText(BuildPrediction(c)) == "c_4pi(0)*Ak(1->2)", PredictionText(BuildPrediction(c)));
    c.reaction.legendre = LegendreExpansion::Absolute;
    check("legendre absolute: the leaf alone", TreeKindOf(c) == TreeKind::Observable && PredictionText(BuildPrediction(c)) == "Ak(1->2)",
          PredictionText(BuildPrediction(c)));
    c.reaction.observable = Observable::AngularDistribution;
    c.reaction.legendre = LegendreExpansion::TwoL2;
    check("an extrapolation in 2l2: no factor in the tree (Segment::LegendreCoefficients maps the output)",
          TreeKindOf(c) == TreeKind::Observable);
    check("LegendreExpansionFactor: 1, 1, 1/(2L+1), 2/(2L+1), 4 pi at L = 3",
          LegendreExpansionFactor(LegendreExpansion::AOverA0, 3) == 1.0 && LegendreExpansionFactor(LegendreExpansion::Absolute, 3) == 1.0 &&
              LegendreExpansionFactor(LegendreExpansion::Rsl, 3) == 1.0 / 7.0 &&
              LegendreExpansionFactor(LegendreExpansion::TwoL2, 3) == 2.0 / 7.0 &&
              LegendreExpansionFactor(LegendreExpansion::FourPi, 3) == 4.0 * 3.14159265358979323846);
  }
  s.reaction.scale.factor = 3.1e12;
  check("a constant: c*sigma", PredictionText(BuildPrediction(s)) == "3100000000000*sigma(1->2)", PredictionText(BuildPrediction(s)));
  s.reaction.combination = ReactionSpec::Combination::Sum;
  s.terms.push_back(Term(1, 3, 1));
  check("over a sum: c*sum(...)", PredictionText(BuildPrediction(s)) == "3100000000000*sum(sigma(1->2), J0/J1*sigma(1->3))",
        PredictionText(BuildPrediction(s)));
  s.reaction.scale = ValueScale();
  check("the identity is no node (a sum)", TreeKindOf(s) == TreeKind::Sum);
}

/// A ratio to a reference table: sigma over the reference leaf of the same term, no frame factors.
void TestReferenceRatio() {
  std::printf("reference ratio\n");
  PredictionSetup s;
  s.reaction.observable = Observable::CrossSection;
  s.terms = {Term(1, 2, 1)};
  s.reaction.reference.energyLab_MeV = {0.1, 1.0};
  s.reaction.reference.sigma_b = {1.0, 2.0};
  const Prediction p = BuildPrediction(s);
  check("kind reference_ratio", p.kind == TreeKind::ReferenceRatio && TreeKindOf(s) == TreeKind::ReferenceRatio);
  check("text", PredictionText(p) == "ratio(sigma(1->2), reference(1->2))", PredictionText(p));
  check("two leaves of term 0, the second the reference, no factors",
        p.leaves.size() == 2 && p.leaves[1].part == LeafPart::Reference && p.leaves[1].term == 0 && p.factors.empty() &&
            IsReferenceLeaf(LeafPart::Reference) && !IsCoulombReference(LeafPart::Reference));
}

void TestTable() {
  std::printf("table\n");
  const std::vector<TreeKindRow> &rows = TreeKinds();
  bool ordered = rows.size() == static_cast<std::size_t>(TreeKind::ReferenceRatio) + 1;  // the enum's last kind
  for (std::size_t k = 0; k < rows.size(); k++) ordered = ordered && static_cast<std::size_t>(rows[k].kind) == k && rows[k].build;
  check("one row per kind, in the order of the enum", ordered);
  bool threw = false;
  try {
    BuildPrediction(PredictionSetup());
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  check("no terms is refused", threw);
}

}  // namespace

int main() {
  TestObservable();
  TestSum();
  TestRatio();
  TestTransforms();
  TestCoulombRatio();
  TestScaled();
  TestReferenceRatio();
  TestTable();
  if (failures) {
    std::printf("%d FAILURE(S)\n", failures);
    return 1;
  }
  std::printf("all passed\n");
  return 0;
}
