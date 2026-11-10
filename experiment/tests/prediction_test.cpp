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
    check(natural ? "sample SUM = ESegment's, bit for bit" : "SUM = ESegment's, bit for bit", ok);
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
  check("RATIO = ESegment's, bit for bit (NaN at a zero denominator)", ok);
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
    // ESegment::TransformTerms and TransformedValue's `transform`
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

void TestTable() {
  std::printf("table\n");
  const std::vector<TreeKindRow> &rows = TreeKinds();
  bool ordered = rows.size() == 5;
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
  TestTable();
  if (failures) {
    std::printf("%d FAILURE(S)\n", failures);
    return 1;
  }
  std::printf("all passed\n");
  return 0;
}
