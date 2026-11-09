// Unit test of the model expressions (azure2::experiment::Tree and its nodes), using experiment/ and physics/ alone
// (docs/upgrade/data_model.md step S8).
//
// value: each node, and the compositions the engine forms today (SUM with frame factors and abundances, RATIO,
// transmission and self-shielded yield with and without sub-points and free-gas Doppler), against the engine's
// expressions written out below as they stand in src/ESegment.cpp (CalculateTheoreticalCrossSection,
// TransformedValue, Transmission, YieldY0, TermFrameFactor), bit for bit, on values with full 53-bit mantissas.  A
// reordered sum or product (a + b + c as a + (b + c), (c J) v as c (J v)) changes the last bit at some of the 4000
// random points and fails here.  Edge cases: a zero RATIO denominator (NaN, also -0), an empty Sum (0), a zero sample
// thickness, a vanishing yield total (the series branch), a term with a sub-point grid of another length (NaN).
//
// backward: every cotangent (leaf values, sub-point values, abundances, sample thickness, Doppler temperature)
// against central finite differences of sum_i g_i m_i, to 1e-6 relative; and against the gradient expressions of
// src/AZUREGrad.cpp written out (GradOnePoint's shared-component route, GradTransformedPoint), bit for bit where the
// tree's reverse sweep forms the same products (README: where it does not).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "azure2/experiment/Expression.h"

using namespace azure2::experiment;

namespace {

int failures = 0;

void check(const std::string &name, bool ok, const std::string &detail = "") {
  std::printf("  %s  %s%s%s\n", ok ? "ok  " : "FAIL", name.c_str(), detail.empty() ? "" : "  ", detail.c_str());
  if (!ok) failures++;
}

bool same(double a, double b) {
  if (std::isnan(a) && std::isnan(b)) return true;
  return std::memcmp(&a, &b, sizeof a) == 0;
}

std::mt19937_64 rng(20261109);

/// A double with a random 53-bit mantissa in [lo, hi) (lo > 0), sign optional.
double rnd(double lo, double hi, bool signedValue = false) {
  std::uniform_real_distribution<double> u(std::log(lo), std::log(hi));
  double x = std::exp(u(rng));
  // fill the last bits: exp of a random argument already has a full mantissa; add a random last-bit pattern
  std::uint64_t bits;
  std::memcpy(&bits, &x, sizeof x);
  bits ^= (rng() & 0xFFFu);
  std::memcpy(&x, &bits, sizeof x);
  if (signedValue && (rng() & 1u)) x = -x;
  return x;
}

const double NaN = std::numeric_limits<double>::quiet_NaN();

// ---------------------------------------------------------------- the engine's expressions, written out

// ESegment::TermFrameFactor (jTerm, jBase: GetCrossSectionKinFactor of the term's and the base's point)
double OldTermFrameFactor(double jBase, double jTerm) {
  if (jTerm == 0.0 || jBase == 0.0 || jTerm == jBase) return 1.0;
  return jBase / jTerm;
}

struct OldComponent {
  double scaling, jBase, jTerm, value, abundance;
};

// ESegment::CalculateTheoreticalCrossSection, SUM and natural branches
double OldSum(double baseValue, bool natural, double baseAbundance, const std::vector<OldComponent> &comps,
              std::vector<double> *naturalTerms = nullptr) {
  if (natural) {
    baseValue *= baseAbundance;
    if (naturalTerms) (*naturalTerms)[0] += baseValue;
  }
  if (comps.empty()) return baseValue;
  double sum = baseValue;
  int k = 1;
  for (const OldComponent &c : comps) {
    double compValue = c.value;
    double term = c.scaling * OldTermFrameFactor(c.jBase, c.jTerm) * compValue;
    if (natural) {
      term *= c.abundance;
      if (naturalTerms) (*naturalTerms)[k] += term;
    }
    sum += term;
    k++;
  }
  return sum;
}

// ESegment::CalculateTheoreticalCrossSection, RATIO branch
double OldRatio(double baseValue, double denominatorValue, double baseConversionFactor,
                double denominatorConversionFactor) {
  if (denominatorValue == 0.0) return NaN;
  double ratioValue = baseValue / denominatorValue * (denominatorConversionFactor / baseConversionFactor);
  return ratioValue;
}

// ESegment::Transmission, ESegment::YieldY0
double OldTransmission(double n, double s, double *dTds, double *dTdn) {
  const double t = std::exp(-n * s);
  if (dTds) *dTds = -n * t;
  if (dTdn) *dTdn = -s * t;
  return t;
}
double OldYieldY0(double n, double s, double st, double *dYds, double *dYdst, double *dYdn) {
  if (dYdn) *dYdn = s * std::exp(-n * st);
  const double u = n * st;
  double g, dg;
  if (std::fabs(u) < 1.0e-4) {
    g = n * (1.0 - u / 2.0 + u * u / 6.0);
    dg = n * n * (-0.5 + u / 3.0 - u * u / 8.0);
  } else {
    const double oneMinus = -std::expm1(-u);
    g = oneMinus / st;
    dg = (n * std::exp(-u) * st - oneMinus) / (st * st);
  }
  if (dYds) *dYds = g;
  if (dYdst) *dYdst = s * dg;
  return s * g;
}

// ESegment::TransformTerms + TransformedValue on plain arrays.  terms: base first, then components; a yield's
// total term last.  sub[q][j]: sub-point j of term q.  kernel c: the base point's IntegrateTargetEffect, here the
// linear test kernel sum_j c_j f_j.  doppler: rows of (sub-point, weight) per kernel sub-point (or the one row at
// the point's energy when !resolution).
struct OldTerm {
  double scaling = 1.0, jBase = 1.0, jTerm = 1.0, abundance = 1.0;
  bool isBase = false, total = false;
  double point = 0.0;
  std::vector<double> sub;
};
struct OldTransformSetup {
  bool yield = false, natural = false, subPoints = false, doppler = false, resolution = true;
  double n = 0.0;
  std::vector<OldTerm> terms;
  std::vector<double> kernel;
  std::vector<std::vector<std::pair<int, double>>> dopplerRows;  // 1-based sub-point index, as DopplerWeights
};
struct OldTT {
  double weight, scale;
  bool total;
};
std::vector<OldTT> OldTransformTerms(const OldTransformSetup &s) {
  std::vector<OldTT> terms;
  for (const OldTerm &t : s.terms) {
    if (t.isBase) {
      terms.push_back({s.natural ? t.abundance : 1.0, 1.0, false});
      continue;
    }
    if (t.total) {
      terms.push_back({1.0, 1.0, true});
      continue;
    }
    const double scale = t.scaling * OldTermFrameFactor(t.jBase, t.jTerm);
    const double w = s.natural ? scale * t.abundance : scale;
    terms.push_back({w, scale, false});
  }
  return terms;
}
double OldKernel(const std::vector<double> &c, const std::vector<double> &f) {
  double y = 0.0;
  for (size_t j = 0; j < f.size(); j++) y += c[j] * f[j];
  return y;
}
double OldTransformedValue(const OldTransformSetup &s) {
  const std::vector<OldTT> terms = OldTransformTerms(s);
  const double n = s.n;
  const bool yield = s.yield;
  auto transform = [&](const std::vector<double> &v) {
    double sum = 0.0, total = 0.0;
    for (size_t q = 0; q < terms.size(); q++) (terms[q].total ? total : sum) += terms[q].weight * v[q];
    return yield ? OldYieldY0(n, sum, total, nullptr, nullptr, nullptr) : OldTransmission(n, sum, nullptr, nullptr);
  };
  std::vector<double> v(terms.size());
  if (!s.subPoints) {
    for (size_t q = 0; q < terms.size(); q++) v[q] = s.terms[q].point;
    return transform(v);
  }
  const int nsp = static_cast<int>(s.terms[0].sub.size());
  for (const OldTerm &t : s.terms)
    if (static_cast<int>(t.sub.size()) != nsp) return NaN;
  std::vector<double> f(nsp);
  std::vector<std::vector<double>> raw;
  auto broadenAt = [&](const std::vector<std::pair<int, double>> &dw) {
    for (size_t q = 0; q < terms.size(); q++) {
      double sv = 0.0;
      for (const auto &t : dw) sv += t.second * raw[q][t.first - 1];
      v[q] = sv;
    }
  };
  if (s.doppler) {
    raw.assign(terms.size(), std::vector<double>(nsp));
    for (size_t q = 0; q < terms.size(); q++)
      for (int j = 1; j <= nsp; j++) raw[q][j - 1] = s.terms[q].sub[j - 1];
    if (!s.resolution) {
      broadenAt(s.dopplerRows[0]);
      return transform(v);
    }
  }
  for (int j = 1; j <= nsp; j++) {
    if (s.doppler) {
      broadenAt(s.dopplerRows[j - 1]);
    } else {
      for (size_t q = 0; q < terms.size(); q++) v[q] = s.terms[q].sub[j - 1];
    }
    f[j - 1] = transform(v);
  }
  return OldKernel(s.kernel, f);
}

// AZUREGrad's GradTransformedPoint on plain arrays: the cotangents of every term's point or sub-point values, of the
// abundances (per term, index q), the thickness and the Doppler temperature, for cotangent fitBar of the model.
struct OldTransformGrad {
  std::vector<double> point;               // per term
  std::vector<std::vector<double>> sub;    // per term, per sub-point
  std::vector<double> abundance;           // per term
  double thickness = 0.0, doppler = 0.0;
};
OldTransformGrad OldGradTransformedPoint(const OldTransformSetup &s, double fitBar,
                                         const std::vector<std::vector<std::pair<int, double>>> &dwdTRows) {
  const std::vector<OldTT> terms = OldTransformTerms(s);
  const double n = s.n;
  const bool yield = s.yield;
  const bool natural = s.natural;
  OldTransformGrad out;
  out.point.assign(terms.size(), 0.0);
  out.sub.assign(terms.size(), std::vector<double>(s.terms[0].sub.size(), 0.0));
  out.abundance.assign(terms.size(), 0.0);
  std::vector<double> v(terms.size());
  auto backProp = [&](double c, auto pointOf) {
    double sum = 0.0, total = 0.0, dS = 0.0, dT = 0.0, dN = 0.0;
    for (size_t q = 0; q < terms.size(); q++) (terms[q].total ? total : sum) += terms[q].weight * v[q];
    if (yield)
      OldYieldY0(n, sum, total, &dS, &dT, &dN);
    else
      OldTransmission(n, sum, &dS, &dN);
    out.thickness += fitBar * c * dN;
    for (size_t q = 0; q < terms.size(); q++) {
      const double d = fitBar * c * (terms[q].total ? dT : dS);
      if (natural && !terms[q].total) out.abundance[q] += d * terms[q].scale * v[q];
      const double cot = d * terms[q].weight;
      if (cot != 0.0) pointOf(q, cot);
    }
  };
  if (!s.subPoints) {
    for (size_t q = 0; q < terms.size(); q++) v[q] = s.terms[q].point;
    backProp(1.0, [&](size_t q, double cot) { out.point[q] += cot; });
    return out;
  }
  const int nsp = static_cast<int>(s.terms[0].sub.size());
  if (s.doppler) {
    std::vector<std::vector<double>> raw(terms.size(), std::vector<double>(nsp));
    for (size_t q = 0; q < terms.size(); q++)
      for (int k = 1; k <= nsp; k++) raw[q][k - 1] = s.terms[q].sub[k - 1];
    std::vector<std::pair<int, double>> stage;  // (row, kernel weight)
    if (!s.resolution) {
      stage.emplace_back(0, 1.0);
    } else {
      for (int j = 1; j <= nsp; j++)
        if (s.kernel[j - 1] != 0.0) stage.emplace_back(j - 1, s.kernel[j - 1]);
    }
    std::vector<std::vector<double>> cot(terms.size(), std::vector<double>(nsp, 0.0));
    for (const auto &st : stage) {
      const auto &dw = s.dopplerRows[st.first];
      const auto &dwdT = dwdTRows[st.first];
      for (size_t q = 0; q < terms.size(); q++) {
        double sv = 0.0;
        for (const auto &t : dw) sv += t.second * raw[q][t.first - 1];
        v[q] = sv;
      }
      backProp(st.second, [&](size_t q, double cq) {
        for (const auto &t : dw) cot[q][t.first - 1] += cq * t.second;
        double ds = 0.0;
        for (const auto &t : dwdT) ds += t.second * raw[q][t.first - 1];
        out.doppler += cq * ds;
      });
    }
    for (size_t q = 0; q < terms.size(); q++)
      for (int k = 1; k <= nsp; k++)
        if (cot[q][k - 1] != 0.0) out.sub[q][k - 1] += cot[q][k - 1];
    return out;
  }
  for (int j = 1; j <= nsp; j++) {
    if (s.kernel[j - 1] == 0.0) continue;
    for (size_t q = 0; q < terms.size(); q++) v[q] = s.terms[q].sub[j - 1];
    backProp(s.kernel[j - 1], [&](size_t q, double cot) { out.sub[q][j - 1] += cot; });
  }
  return out;
}

// ---------------------------------------------------------------- test doubles

/// A linear kernel over the stage rows, with an optional premap (Doppler) in CSR form.
class TestFold : public SubPointFold {
 public:
  TransformStage st;
  std::vector<double> c;
  bool linear = true;
  TransformStage stage(std::size_t) const override { return st; }
  double fold(std::size_t, const double *f, std::size_t n) const override {
    double y = 0.0;
    for (std::size_t j = 0; j < n; j++) y += c[j] * f[j];
    return y;
  }
  bool weights(std::size_t, double *w, std::size_t n) const override {
    if (!linear) return false;
    for (std::size_t j = 0; j < n; j++) w[j] = c[j];
    return true;
  }
};

struct Csr {
  std::vector<std::size_t> start{0};
  std::vector<int> col;
  std::vector<double> w;
  void addRow(const std::vector<std::pair<int, double>> &row) {  // 1-based columns, as DopplerWeights
    for (const auto &e : row) {
      col.push_back(e.first - 1);
      w.push_back(e.second);
    }
    start.push_back(col.size());
  }
  SparseRows view() const { return SparseRows{start.size() - 1, start.data(), col.data(), w.data()}; }
};

std::vector<std::pair<int, double>> RandomRow(int nsp) {
  std::vector<std::pair<int, double>> row;
  const int len = 2 + static_cast<int>(rng() % 4);
  int k = 1 + static_cast<int>(rng() % nsp);
  for (int m = 0; m < len; m++) {
    row.emplace_back(k, rnd(0.01, 1.0));
    k = 1 + static_cast<int>((k + rng() % 3) % nsp);
  }
  return row;
}

// ---------------------------------------------------------------- finite differences

/// sum_i g_i m_i over the points with g_i != 0 (a NaN model value elsewhere does not enter), or g_i m_i at one point
double Objective(const Tree &t, const Inputs &in, std::size_t n, const std::vector<double> &g, long only = -1) {
  Workspace ws;
  std::vector<double> m(n);
  t.Evaluate(in, n, m.data(), ws);
  double s = 0.0;
  for (std::size_t i = 0; i < n; i++)
    if (g[i] != 0.0 && (only < 0 || static_cast<long>(i) == only)) s += g[i] * m[i];
  return s;
}

/// Central difference of the objective with respect to *x (a value of point `only`, or a parameter: -1), Richardson
/// extrapolated, (4 D(h/2) - D(h)) / 3: truncation error O(h^4).
double FD(const Tree &t, const Inputs &in, std::size_t n, const std::vector<double> &g, double *x, long only = -1) {
  const double x0 = *x;
  auto D = [&](double h) {
    *x = x0 + h;
    const double fp = Objective(t, in, n, g, only);
    *x = x0 - h;
    const double fm = Objective(t, in, n, g, only);
    *x = x0;
    return (fp - fm) / (2.0 * h);
  };
  const double h = 1.0e-3 * std::max(std::fabs(x0), 1.0);
  return (4.0 * D(h / 2.0) - D(h)) / 3.0;
}

/// Every cotangent against its finite difference: |a - fd| <= 1e-6 max(|a|, |fd|) + 1e-9 max_k |a_k| (the absolute
/// floor for cotangents that vanish, as a sub-point with kernel weight 0).
struct FdCheck {
  std::vector<std::pair<double, double>> pairs;
  int bad = 0;
  double worst = 0.0;
  void one(double analytic, double fd) { pairs.emplace_back(analytic, fd); }
  bool ok() {
    double scale = 0.0;
    for (const auto &p : pairs) scale = std::max(scale, std::fabs(p.first));
    bad = 0;
    worst = 0.0;
    for (const auto &p : pairs) {
      const double diff = std::fabs(p.first - p.second), mag = std::max(std::fabs(p.first), std::fabs(p.second));
      if (!(diff <= 1e-6 * mag + 1e-9 * scale)) bad++;
      worst = std::max(worst, diff / (mag + 1e-9 * scale));
    }
    return bad == 0;
  }
  std::string text() const {
    char b[96];
    std::snprintf(b, sizeof b, "%zu derivatives, worst relative difference %.1e", pairs.size(), worst);
    return b;
  }
};

// ---------------------------------------------------------------- tests

void TestNodesAlone() {
  std::printf("nodes alone, against the engine's expressions (bit for bit)\n");
  bool sumOk = true, scaleOk = true, weightedOk = true, ratioOk = true, tffOk = true;
  for (int k = 0; k < 4000; k++) {
    const double a = rnd(1e-3, 1e3, true), b = rnd(1e-3, 1e3, true), c = rnd(1e-3, 1e3, true);
    const double v3[3] = {a, b, c};
    double ref = a;
    ref += b;
    ref += c;
    sumOk = sumOk && same(node::SumValue(v3, 3), ref);
    const double s = rnd(0.1, 10), J = rnd(0.5, 2), v = rnd(1e-6, 1e3);
    scaleOk = scaleOk && same(node::ScaleValue(s, J, v), s * J * v);
    double base = v;
    base *= J;
    weightedOk = weightedOk && same(node::WeightedValue(J, v), base);
    const double num = rnd(1e-6, 1e3), den = rnd(1e-6, 1e3), jn = rnd(0.5, 2), jd = rnd(0.5, 2);
    ratioOk = ratioOk && same(node::RatioValue(num, den, jn, jd), OldRatio(num, den, jn, jd));
    const double jb = (k % 7 == 0) ? 0.0 : (k % 5 == 0) ? jd : jn;
    tffOk = tffOk && same(TermFrameFactor(jb, jd), OldTermFrameFactor(jb, jd));
  }
  check("Sum: v0 + v1 + v2 from v0", sumOk);
  check("Scale: (c J) v", scaleOk);
  check("Weighted: v a", weightedOk);
  check("Ratio: num / den (J_den / J_num)", ratioOk);
  check("TermFrameFactor: J_base / J_term, 1 for 0 or equal", tffOk);
  check("Ratio: den = 0 gives NaN", std::isnan(node::RatioValue(1.5, 0.0, 1.0, 1.0)));
  check("Ratio: den = -0 gives NaN", std::isnan(node::RatioValue(1.5, -0.0, 1.0, 1.0)));
  check("Ratio: num = 0 gives 0", node::RatioValue(0.0, 2.5, 1.1, 0.9) == 0.0);
  check("Sum: no children gives 0", node::SumValue(nullptr, 0) == 0.0);
  {
    const double s0 = -0.0, s1 = 0.0;
    double v2[2] = {s0, s1};
    double ref = s0;
    ref += s1;
    check("Sum: -0 + 0 as the engine (sum from the first term)", same(node::SumValue(v2, 2), ref));
  }
  bool trOk = true, yOk = true;
  for (int k = 0; k < 4000; k++) {
    const double n = (k % 11 == 0) ? 0.0 : rnd(1e-6, 1.0), s = rnd(1e-3, 1e4);
    double a1, a2, b1, b2;
    trOk = trOk && same(node::Transmission(n, s, &a1, &a2), OldTransmission(n, s, &b1, &b2)) && same(a1, b1) &&
           same(a2, b2);
    // st spans both branches of the series switch (|n st| < 1e-4)
    const double st = (k % 13 == 0) ? 0.0 : (k % 2) ? rnd(1e-9, 1e-2) : rnd(1e-3, 1e4);
    double c1, c2, c3, d1, d2, d3;
    yOk = yOk && same(node::YieldY0(n, s, st, &c1, &c2, &c3), OldYieldY0(n, s, st, &d1, &d2, &d3)) && same(c1, d1) &&
          same(c2, d2) && same(c3, d3);
  }
  check("Transmission exp(-n s) and its derivatives", trOk);
  check("YieldY0 (s/st)(1 - e^{-n st}) and its derivatives, both branches", yOk);
  {
    double dTds, dTdn;
    const double T = node::Transmission(0.0, 7.25, &dTds, &dTdn);
    double dYds, dYdst, dYdn;
    const double Y = node::YieldY0(0.0, 7.25, 3.5, &dYds, &dYdst, &dYdn);
    check("zero thickness: T = 1, dT/ds = 0, dT/dn = -s", T == 1.0 && dTds == 0.0 && dTdn == -7.25);
    check("zero thickness: Y0 = 0, dY0/ds = 0, dY0/dn = s", Y == 0.0 && dYds == 0.0 && dYdn == 7.25);
    const double Y0 = node::YieldY0(0.02, 7.25, 0.0, &dYds, &dYdst, &dYdn);
    check("vanishing total: Y0 = n s (series, no 0/0)", Y0 == 7.25 * 0.02 && std::isfinite(dYdst));
  }
}

/// A SUM of a base and K components with frame factors, natural (abundances) or not, over N points.
void TestSum(bool natural) {
  std::printf("SUM %s, %d points\n", natural ? "with abundances (sample)" : "with frame factors", 1000);
  const std::size_t N = 1000, K = 3;
  std::vector<std::vector<double>> leaf(K + 1, std::vector<double>(N));
  std::vector<std::vector<double>> jTerm(K, std::vector<double>(N)), factor(K, std::vector<double>(N));
  std::vector<double> jBase(N), scaling(K), abundance(K + 1);
  for (auto &l : leaf)
    for (auto &x : l) x = rnd(1e-6, 1e3);
  for (std::size_t i = 0; i < N; i++) {
    jBase[i] = (i % 17 == 0) ? 0.0 : rnd(0.5, 2.0);
    for (std::size_t k = 0; k < K; k++) jTerm[k][i] = (i % 23 == 0) ? jBase[i] : rnd(0.5, 2.0);
  }
  for (auto &s : scaling) s = rnd(0.1, 3.0);
  for (auto &a : abundance) a = rnd(1e-3, 1.0);
  for (std::size_t k = 0; k < K; k++)
    for (std::size_t i = 0; i < N; i++) factor[k][i] = TermFrameFactor(jBase[i], jTerm[k][i]);

  Tree t;
  std::vector<int> terms;
  int b = t.Leaf(0);
  terms.push_back(natural ? t.Weighted(0, b, 0) : b);
  for (std::size_t k = 0; k < K; k++) {
    int s = t.Scale(scaling[k], t.Leaf(static_cast<int>(k + 1)), static_cast<int>(k));
    terms.push_back(natural ? t.Weighted(static_cast<int>(k + 1), s, static_cast<int>(k + 1)) : s);
  }
  t.Sum(terms);
  Inputs in;
  for (auto &l : leaf) in.leaves.push_back({l.data(), 1});
  for (auto &f : factor) in.factors.push_back({f.data(), 1});
  in.params = abundance.data();

  Workspace ws;
  std::vector<double> m(N);
  t.Evaluate(in, N, m.data(), ws);
  bool ok = true, okValue = true, okNatural = true;
  for (std::size_t i = 0; i < N; i++) {
    std::vector<OldComponent> comps;
    for (std::size_t k = 0; k < K; k++) comps.push_back({scaling[k], jBase[i], jTerm[k][i], leaf[k + 1][i], abundance[k + 1]});
    std::vector<double> natOld(K + 1, 0.0), natNew(K + 1, 0.0);
    const double ref = OldSum(leaf[0][i], natural, abundance[0], comps, &natOld);
    ok = ok && same(m[i], ref);
    okValue = okValue && same(t.Value(in, i, ws), ref);
    t.Natural(in, i, natNew.data(), ws);
    for (std::size_t k = 0; k <= K; k++) okNatural = okNatural && same(natNew[k], natOld[k]);
  }
  check("Evaluate = ESegment SUM, bit for bit", ok);
  check("Value(i) = Evaluate", okValue);
  if (natural) check("Natural: per-nucleus terms = ESegment's SetNaturalTerms", okNatural);

  // backward: finite differences on a few points (the leaves and the abundances)
  const std::size_t n = 5;
  std::vector<double> g(n);
  for (auto &x : g) x = rnd(0.1, 2.0, true);
  std::vector<std::vector<double>> leafBar(K + 1, std::vector<double>(N, 0.0));
  std::vector<double> paramBar(K + 1, 0.0);
  Cotangents cot;
  for (auto &l : leafBar) cot.leaves.push_back(l.data());
  cot.params = paramBar.data();
  check("Backward runs", t.Backward(in, n, g.data(), cot, ws));
  FdCheck fd;
  for (std::size_t k = 0; k <= K; k++)
    for (std::size_t i = 0; i < n; i++) fd.one(leafBar[k][i], FD(t, in, n, g, &leaf[k][i], static_cast<long>(i)));
  if (natural)
    for (std::size_t k = 0; k <= K; k++) fd.one(paramBar[k], FD(t, in, n, g, &abundance[k]));
  {
    const bool good = fd.ok();
    check("Backward = finite differences", good, fd.text());
  }

  // AZUREGrad (GradOnePoint, shared-component route): coef = s * w, cotangent fitBar * coef; abundance fitBar * s * v.
  bool agLeaf = true, agAbund = true;
  double worst = 0.0;
  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t q = 0; q <= K; q++) {
      const double s = q > 0 ? scaling[q - 1] * OldTermFrameFactor(jBase[i], jTerm[q - 1][i]) : 1.0;
      const double w = natural ? abundance[q] : 1.0;
      const double coef = s * w;
      const double ag = g[i] * coef;
      // the tree's cotangent of this leaf at this point alone
      std::vector<std::vector<double>> lb(K + 1, std::vector<double>(N, 0.0));
      std::vector<double> pb(K + 1, 0.0);
      Cotangents c1;
      for (auto &l : lb) c1.leaves.push_back(l.data());
      c1.params = pb.data();
      std::vector<double> gi(n, 0.0);
      gi[i] = g[i];
      t.Backward(in, n, gi.data(), c1, ws);
      if (!same(lb[q][i], ag)) agLeaf = false;
      worst = std::max(worst, std::fabs(lb[q][i] - ag) / std::fabs(ag));
      if (natural && !same(pb[q], g[i] * s * leaf[q][i])) agAbund = false;
    }
  }
  char d[64];
  std::snprintf(d, sizeof d, "largest relative difference %.1e", worst);
  if (!natural) {
    check("Backward leaves = AZUREGrad's expressions, bit for bit", agLeaf, d);
  } else {
    // (fitBar a) s here, fitBar (s a) in AZUREGrad: equal to the rounding of one product (README).
    check("Backward leaves = AZUREGrad's expressions to 2 ulp (different product order)", worst <= 4.5e-16, d);
    std::printf("        bitwise equal to AZUREGrad: leaves %s, abundances %s\n", agLeaf ? "yes" : "no",
                agAbund ? "yes" : "no");
  }
}

void TestRatio() {
  std::printf("RATIO with frame factors, 1000 points\n");
  const std::size_t N = 1000;
  std::vector<double> num(N), den(N), jn(N), jd(N);
  for (std::size_t i = 0; i < N; i++) {
    num[i] = (i % 31 == 0) ? 0.0 : rnd(1e-6, 1e3);
    den[i] = (i % 29 == 0) ? 0.0 : (i % 37 == 0) ? -0.0 : rnd(1e-6, 1e3);
    jn[i] = rnd(0.5, 2.0);
    jd[i] = rnd(0.5, 2.0);
  }
  Tree t;
  t.Ratio(t.Leaf(0), t.Leaf(1), 0, 1);
  Inputs in;
  in.leaves = {{num.data(), 1}, {den.data(), 1}};
  in.factors = {{jn.data(), 1}, {jd.data(), 1}};
  Workspace ws;
  std::vector<double> m(N);
  t.Evaluate(in, N, m.data(), ws);
  bool ok = true, nanOk = true;
  for (std::size_t i = 0; i < N; i++) {
    ok = ok && same(m[i], OldRatio(num[i], den[i], jn[i], jd[i]));
    if (den[i] == 0.0) nanOk = nanOk && std::isnan(m[i]);
  }
  check("Evaluate = ESegment RATIO, bit for bit", ok);
  check("zero (and -0) denominator: NaN", nanOk);

  // backward away from den = 0; one point with num = 0
  std::vector<std::size_t> pts;
  for (std::size_t i = 0; i < N && pts.size() < 6; i++)
    if (den[i] != 0.0 && (num[i] == 0.0 || pts.size() < 5)) pts.push_back(i);
  std::vector<double> g(N, 0.0);
  for (std::size_t i : pts) g[i] = rnd(0.1, 2.0, true);
  std::vector<double> nb(N, 0.0), db(N, 0.0);
  Cotangents cot;
  cot.leaves = {nb.data(), db.data()};
  t.Backward(in, N, g.data(), cot, ws);
  FdCheck fd;
  bool ag = true;
  for (std::size_t i : pts) {
    fd.one(nb[i], FD(t, in, N, g, &num[i], static_cast<long>(i)));
    fd.one(db[i], FD(t, in, N, g, &den[i], static_cast<long>(i)));
    // AZUREGrad shared route: coef[0] = model / baseVal, coef[1] = -model / denomVal, cotangent fitBar * coef
    const double model = m[i];
    if (num[i] != 0.0) ag = ag && same(nb[i], g[i] * (model / num[i])) && same(db[i], g[i] * (-model / den[i]));
  }
  {
    const bool good = fd.ok();
    check("Backward = finite differences (incl. a zero numerator)", good, fd.text());
  }
  check("Backward = AZUREGrad's expressions, bit for bit (num != 0)", ag);
  std::vector<double> g1(N, 0.0);
  std::size_t z = 0;
  while (den[z] != 0.0) z++;
  g1[z] = 1.0;
  std::vector<double> nb1(N, 0.0), db1(N, 0.0);
  Cotangents c1;
  c1.leaves = {nb1.data(), db1.data()};
  t.Backward(in, N, g1.data(), c1, ws);
  check("Backward at a zero denominator: NaN", std::isnan(nb1[z]) && std::isnan(db1[z]));
}

/// Transmission or yield, natural or not, without sub-points, with sub-points, with Doppler (with and without a
/// kernel after it).
void TestTransform(bool yield, bool natural, int mode) {
  // mode 0: no sub-point effect; 1: sub-points and kernel; 2: Doppler then kernel; 3: Doppler, no kernel
  static const char *modes[] = {"at the point", "sub-points + kernel", "Doppler + kernel", "Doppler only"};
  std::printf("%s%s, %s\n", yield ? "self-shielded yield" : "transmission", natural ? " of a sample" : "",
              modes[mode]);
  const std::size_t N = 200, K = 2;  // base + K components (+ total for a yield)
  const std::size_t nsp = 9;
  const std::size_t nTerms = 1 + K + (yield ? 1 : 0);
  // leaf inputs: one per term, point values and sub-point values
  std::vector<std::vector<double>> point(nTerms, std::vector<double>(N)), sub(nTerms, std::vector<double>(N * nsp));
  std::vector<std::size_t> offset(N + 1);
  for (std::size_t i = 0; i <= N; i++) offset[i] = i * nsp;
  for (std::size_t q = 0; q < nTerms; q++) {
    for (auto &x : point[q]) x = rnd(1e-3, 50.0);
    for (auto &x : sub[q]) x = rnd(1e-3, 50.0);
  }
  std::vector<std::vector<double>> jTerm(K, std::vector<double>(N)), factor(K, std::vector<double>(N));
  std::vector<double> jBase(N), scaling(K);
  for (std::size_t i = 0; i < N; i++) {
    jBase[i] = rnd(0.5, 2.0);
    for (std::size_t k = 0; k < K; k++) {
      jTerm[k][i] = rnd(0.5, 2.0);
      factor[k][i] = TermFrameFactor(jBase[i], jTerm[k][i]);
    }
  }
  for (auto &s : scaling) s = rnd(0.2, 2.0);
  // params: 0 thickness, 1..K+1 abundances, K+2 Doppler temperature
  std::vector<double> params(K + 3);
  params[0] = (mode == 1 && !yield) ? rnd(1e-3, 0.05) : rnd(1e-3, 0.05);
  for (std::size_t k = 1; k <= K + 1; k++) params[k] = rnd(0.05, 1.0);
  params[K + 2] = 1.0;

  TestFold fold;
  for (std::size_t j = 0; j < nsp; j++) fold.c.push_back((j % 4 == 3) ? 0.0 : rnd(0.01, 1.0));
  Csr premap, premapDT;
  std::vector<std::vector<std::pair<int, double>>> dRows, dTRows;
  if (mode >= 2) {
    const std::size_t rows = (mode == 3) ? 1 : nsp;
    for (std::size_t r = 0; r < rows; r++) {
      dRows.push_back(RandomRow(static_cast<int>(nsp)));
      premap.addRow(dRows.back());
      std::vector<std::pair<int, double>> dt = dRows.back();
      for (auto &e : dt) e.second = rnd(0.01, 1.0, true);
      premapDT.addRow(dt);
      dTRows.push_back(dt);
    }
  }
  fold.st.subPoints = mode >= 1;
  fold.st.nsp = nsp;
  if (mode >= 2) {
    fold.st.premap = premap.view();
    fold.st.premapDT = premapDT.view();
    fold.st.dopplerParam = static_cast<int>(K + 2);
  }
  fold.st.fold = mode != 3;

  TransformSpec spec;
  spec.kind = yield ? TransformSpec::Kind::Yield : TransformSpec::Kind::Transmission;
  spec.thickness = 0;
  spec.fold = 0;
  {
    TransformTerm base;
    base.input = 0;
    if (natural) {
      base.abundance = 1;
      base.group = 0;
    }
    spec.terms.push_back(base);
    for (std::size_t k = 0; k < K; k++) {
      TransformTerm c;
      c.input = static_cast<int>(k + 1);
      c.c = scaling[k];
      c.factor = static_cast<int>(k);
      if (natural) {
        c.abundance = static_cast<int>(k + 2);
        c.group = static_cast<int>(k + 1);
      }
      spec.terms.push_back(c);
    }
    if (yield) {
      TransformTerm tot;
      tot.input = static_cast<int>(K + 1);
      tot.total = true;
      spec.terms.push_back(tot);
    }
  }
  Tree t;
  t.Transform(spec);
  Inputs in;
  for (std::size_t q = 0; q < nTerms; q++) {
    in.leaves.push_back({point[q].data(), 1});
    in.subLeaves.push_back({sub[q].data(), offset.data()});
  }
  for (auto &f : factor) in.factors.push_back({f.data(), 1});
  in.params = params.data();
  in.folds = {&fold};

  Workspace ws;
  std::vector<double> m(N);
  t.Evaluate(in, N, m.data(), ws);
  bool ok = true, okNat = true;
  std::vector<OldTransformSetup> setups(N);
  for (std::size_t i = 0; i < N; i++) {
    OldTransformSetup &s = setups[i];
    s.yield = yield;
    s.natural = natural;
    s.subPoints = mode >= 1;
    s.doppler = mode >= 2;
    s.resolution = mode != 3;
    s.n = params[0];
    s.kernel = fold.c;
    s.dopplerRows = dRows;
    for (std::size_t q = 0; q < nTerms; q++) {
      OldTerm o;
      o.isBase = q == 0;
      o.total = yield && q == nTerms - 1;
      if (q >= 1 && q <= K) {
        o.scaling = scaling[q - 1];
        o.jBase = jBase[i];
        o.jTerm = jTerm[q - 1][i];
      }
      o.abundance = (q <= K) ? params[q + 1] : 1.0;
      o.point = point[q][i];
      o.sub.assign(sub[q].begin() + i * nsp, sub[q].begin() + (i + 1) * nsp);
      s.terms.push_back(o);
    }
    ok = ok && same(m[i], OldTransformedValue(s));
    if (natural) {
      // ESegment: naturalTerms[nucleus] += weight * point value, every non-total term
      std::vector<OldTT> tt = OldTransformTerms(s);
      std::vector<double> natOld(K + 1, 0.0), natNew(K + 1, 0.0);
      for (std::size_t q = 0; q <= K; q++) natOld[q] += tt[q].weight * point[q][i];
      t.Natural(in, i, natNew.data(), ws);
      for (std::size_t q = 0; q <= K; q++) okNat = okNat && same(natOld[q], natNew[q]);
    }
  }
  check("Evaluate = ESegment::TransformedValue, bit for bit", ok);
  if (natural) check("Natural = ESegment's SetNaturalTerms", okNat);

  // backward against finite differences on 3 points
  const std::size_t n = 3;
  std::vector<double> g(n);
  for (auto &x : g) x = rnd(0.1, 2.0, true);
  std::vector<std::vector<double>> pb(nTerms, std::vector<double>(N, 0.0)), sb(nTerms, std::vector<double>(N * nsp, 0.0));
  std::vector<double> parBar(params.size(), 0.0);
  Cotangents cot;
  for (std::size_t q = 0; q < nTerms; q++) {
    cot.leaves.push_back(pb[q].data());
    cot.subLeaves.push_back(sb[q].data());
  }
  cot.params = parBar.data();
  check("Backward runs", t.Backward(in, n, g.data(), cot, ws));
  FdCheck fd;
  for (std::size_t q = 0; q < nTerms; q++)
    for (std::size_t i = 0; i < n; i++) {
      if (mode == 0) {
        fd.one(pb[q][i], FD(t, in, n, g, &point[q][i], static_cast<long>(i)));
      } else {
        for (std::size_t j = 0; j < nsp; j++) fd.one(sb[q][i * nsp + j], FD(t, in, n, g, &sub[q][i * nsp + j], static_cast<long>(i)));
      }
    }
  fd.one(parBar[0], FD(t, in, n, g, &params[0]));
  if (natural)
    for (std::size_t k = 1; k <= K + 1; k++) fd.one(parBar[k], FD(t, in, n, g, &params[k]));
  if (mode >= 2) {
    // the temperature enters through the premap: d/dT of sum_i g_i m_i with D -> D + dT dD/dT
    std::vector<double> saved = premap.w;
    auto D = [&](double h) {
      for (std::size_t k = 0; k < premap.w.size(); k++) premap.w[k] = saved[k] + h * premapDT.w[k];
      fold.st.premap = premap.view();
      const double fp = Objective(t, in, n, g);
      for (std::size_t k = 0; k < premap.w.size(); k++) premap.w[k] = saved[k] - h * premapDT.w[k];
      fold.st.premap = premap.view();
      const double fm = Objective(t, in, n, g);
      premap.w = saved;
      fold.st.premap = premap.view();
      return (fp - fm) / (2 * h);
    };
    fd.one(parBar[K + 2], (4.0 * D(5e-4) - D(1e-3)) / 3.0);
  }
  {
    const bool good = fd.ok();
    check("Backward = finite differences", good, fd.text());
  }

  {
    // AZUREGrad's GradTransformedPoint, point by point: the same products in the same order
    bool same_ = true;
    for (std::size_t i = 0; i < 20; i++) {
      std::vector<double> gi(i + 1, 0.0);
      gi[i] = rnd(0.1, 2.0, true);
      std::vector<std::vector<double>> pb1(nTerms, std::vector<double>(N, 0.0)), sb1(nTerms, std::vector<double>(N * nsp, 0.0));
      std::vector<double> par1(params.size(), 0.0);
      Cotangents c1;
      for (std::size_t q = 0; q < nTerms; q++) {
        c1.leaves.push_back(pb1[q].data());
        c1.subLeaves.push_back(sb1[q].data());
      }
      c1.params = par1.data();
      t.Backward(in, i + 1, gi.data(), c1, ws);
      const OldTransformGrad ag = OldGradTransformedPoint(setups[i], gi[i], dTRows);
      for (std::size_t q = 0; q < nTerms; q++) {
        same_ = same_ && same(pb1[q][i], ag.point[q]);
        for (std::size_t j = 0; j < nsp; j++) same_ = same_ && same(sb1[q][i * nsp + j], ag.sub[q][j]);
        if (natural && q <= K) same_ = same_ && same(par1[q + 1], ag.abundance[q]);
      }
      same_ = same_ && same(par1[0], ag.thickness);
      if (mode >= 2) same_ = same_ && same(par1[K + 2], ag.doppler);
    }
    check("Backward = AZUREGrad's GradTransformedPoint, bit for bit", same_);
  }

  if (mode != 3) {
    // a kernel that is not linear: the backward refuses (AZUREGrad falls back to finite differences)
    TestFold nl = fold;
    nl.linear = false;
    if (mode >= 1) {
      Inputs in2 = in;
      in2.folds = {&nl};
      Cotangents none;
      none.leaves.assign(nTerms, nullptr);
      none.subLeaves.assign(nTerms, nullptr);
      check("non-linear kernel: Backward returns false", !t.Backward(in2, n, g.data(), none, ws));
    }
  }
}

void TestEdges() {
  std::printf("edge cases\n");
  // empty Sum
  {
    Tree t;
    t.Sum({});
    Inputs in;
    Workspace ws;
    check("Sum of nothing = 0", t.Value(in, 0, ws) == 0.0);
  }
  // a term with another sub-point count: NaN (value) and false (backward)
  {
    std::vector<double> p0{1.0}, p1{2.0}, s0{1, 2, 3}, s1{1, 2};
    std::vector<std::size_t> o0{0, 3}, o1{0, 2};
    TestFold fold;
    fold.st.subPoints = true;
    fold.st.nsp = 3;
    fold.c = {0.2, 0.5, 0.3};
    TransformSpec spec;
    spec.thickness = 0;
    spec.fold = 0;
    spec.terms.resize(2);
    spec.terms[0].input = 0;
    spec.terms[1].input = 1;
    Tree t;
    t.Transform(spec);
    double n = 0.01;
    Inputs in;
    in.leaves = {{p0.data(), 1}, {p1.data(), 1}};
    in.subLeaves = {{s0.data(), o0.data()}, {s1.data(), o1.data()}};
    in.params = &n;
    in.folds = {&fold};
    Workspace ws;
    check("sub-point grids of different length: NaN", std::isnan(t.Value(in, 0, ws)));
    double g = 1.0;
    Cotangents cot;
    cot.leaves.assign(2, nullptr);
    cot.subLeaves.assign(2, nullptr);
    check("... and Backward false", !t.Backward(in, 1, &g, cot, ws));
  }
  // zero thickness through the tree: T = 1 at every sub-point, the kernel of 1
  {
    std::vector<double> p0{3.0}, s0{1.1, 2.2, 3.3};
    std::vector<std::size_t> o0{0, 3};
    TestFold fold;
    fold.st.subPoints = true;
    fold.st.nsp = 3;
    fold.c = {0.25, 0.5, 0.25};
    TransformSpec spec;
    spec.thickness = 0;
    spec.fold = 0;
    spec.terms.resize(1);
    spec.terms[0].input = 0;
    Tree t;
    t.Transform(spec);
    double n = 0.0;
    Inputs in;
    in.leaves = {{p0.data(), 1}};
    in.subLeaves = {{s0.data(), o0.data()}};
    in.params = &n;
    in.folds = {&fold};
    Workspace ws;
    check("transmission at zero thickness = kernel of 1", t.Value(in, 0, ws) == 1.0);
  }
  // a refused structure
  {
    Tree t;
    int l = t.Leaf(0);
    t.Scale(2.0, l);
    bool threw = false;
    try {
      t.Scale(3.0, l);
    } catch (const std::invalid_argument &) {
      threw = true;
    }
    check("a node used twice is refused", threw);
  }
  // strided leaves: the slot layout of a unit with several quantities
  {
    std::vector<double> slots{1.5, 9.0, 2.5, 9.0, 3.5, 9.0};
    Tree t;
    t.Scale(2.0, t.Leaf(0));
    Inputs in;
    in.leaves = {{slots.data(), 2}};
    Workspace ws;
    std::vector<double> m(3);
    t.Evaluate(in, 3, m.data(), ws);
    check("strided leaf view", m[0] == 3.0 && m[1] == 5.0 && m[2] == 7.0);
  }
}

}  // namespace

int main() {
  TestNodesAlone();
  TestSum(false);
  TestSum(true);
  TestRatio();
  for (int yield = 0; yield <= 1; yield++)
    for (int natural = 0; natural <= 1 - yield; natural++)  // a sample yield is refused (EData)
      for (int mode = 0; mode <= 3; mode++) TestTransform(yield == 1, natural == 1, mode);
  TestEdges();
  if (failures) {
    std::printf("%d FAILURE(S)\n", failures);
    return 1;
  }
  std::printf("all passed\n");
  return 0;
}
