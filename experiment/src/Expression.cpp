#include "azure2/experiment/Expression.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace azure2 {
namespace experiment {

namespace {
const double kNaN = std::numeric_limits<double>::quiet_NaN();
}  // namespace

namespace node {

double RatioValue(double num, double den, double jNum, double jDen) {
  // No ratio exists where the denominator vanishes (its exit channel closed, say): NaN, never 0 (ESegment, since
  // 2026-11-01).
  if (den == 0.0) return kNaN;
  return num / den * (jDen / jNum);
}

void RatioBackward(double g, double num, double den, double jNum, double jDen, double m, double &cotNum,
                   double &cotDen) {
  if (den == 0.0) {
    cotNum = cotDen = kNaN;
    return;
  }
  // dm/dnum = m / num and dm/dden = -m / den, formed as AZUREGrad's shared-component route forms them (coef[0],
  // coef[1], then fitBar * coef).  m / num is 0/0 at num = 0, where AZUREGrad gives the numerator no cotangent; the
  // derivative there is (J_den / J_num) / den.
  cotNum = (num != 0.0) ? g * (m / num) : g * ((jDen / jNum) / den);
  cotDen = g * (-m / den);
}

double Transmission(double n, double s, double *dTds, double *dTdn) {
  const double t = std::exp(-n * s);
  if (dTds) *dTds = -n * t;
  if (dTdn) *dTdn = -s * t;
  return t;
}

double YieldY0(double n, double s, double st, double *dYds, double *dYdst, double *dYdn) {
  if (dYdn) *dYdn = s * std::exp(-n * st);  // d/dn [(s/st)(1 - e^{-n st})] = s e^{-n st}, also for st -> 0
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

}  // namespace node

// ---------------------------------------------------------------- building

int Tree::Add(Node nd, const std::vector<int> &children) {
  const int id = static_cast<int>(nodes_.size());
  nd.first = static_cast<int>(children_.size());
  nd.count = static_cast<int>(children.size());
  for (int c : children) {
    if (c < 0 || c >= id) throw std::invalid_argument("experiment::Tree: child " + std::to_string(c) + " does not exist");
    if (used_[c]) throw std::invalid_argument("experiment::Tree: node " + std::to_string(c) + " used twice");
    used_[c] = 1;
    children_.push_back(c);
  }
  nodes_.push_back(nd);
  used_.push_back(0);
  return id;
}

int Tree::Leaf(int input) {
  if (input < 0) throw std::invalid_argument("experiment::Tree: negative leaf input");
  Node nd;
  nd.kind = NodeKind::Leaf;
  nd.input = input;
  numInputs_ = std::max(numInputs_, input + 1);
  return Add(nd, {});
}

int Tree::Sum(const std::vector<int> &children) {
  Node nd;
  nd.kind = NodeKind::Sum;
  return Add(nd, children);
}

int Tree::Scale(double c, int child, int factor) {
  Node nd;
  nd.kind = NodeKind::Scale;
  nd.c = c;
  nd.factor = factor;
  numFactors_ = std::max(numFactors_, factor + 1);
  return Add(nd, {child});
}

int Tree::Weighted(int param, int child, int group) {
  if (param < 0) throw std::invalid_argument("experiment::Tree: Weighted needs a parameter");
  Node nd;
  nd.kind = NodeKind::Weighted;
  nd.param = param;
  nd.group = group;
  numParams_ = std::max(numParams_, param + 1);
  return Add(nd, {child});
}

int Tree::Ratio(int num, int den, int jNum, int jDen) {
  if ((jNum < 0) != (jDen < 0)) throw std::invalid_argument("experiment::Tree: Ratio needs both frame factors or neither");
  Node nd;
  nd.kind = NodeKind::Ratio;
  nd.factor = jNum;
  nd.factor2 = jDen;
  numFactors_ = std::max(numFactors_, std::max(jNum, jDen) + 1);
  return Add(nd, {num, den});
}

int Tree::Transform(const TransformSpec &spec) {
  if (spec.thickness < 0) throw std::invalid_argument("experiment::Tree: Transform needs a thickness parameter");
  if (spec.terms.empty()) throw std::invalid_argument("experiment::Tree: Transform without terms");
  numParams_ = std::max(numParams_, spec.thickness + 1);
  for (const TransformTerm &t : spec.terms) {
    if (t.input < 0) throw std::invalid_argument("experiment::Tree: Transform term without a leaf input");
    if (t.total && spec.kind != TransformSpec::Kind::Yield)
      throw std::invalid_argument("experiment::Tree: a total term belongs to a yield");
    numInputs_ = std::max(numInputs_, t.input + 1);
    numFactors_ = std::max(numFactors_, t.factor + 1);
    numParams_ = std::max(numParams_, t.abundance + 1);
  }
  Node nd;
  nd.kind = NodeKind::Transform;
  nd.input = static_cast<int>(transforms_.size());
  transforms_.push_back(spec);
  return Add(nd, {});
}

// ---------------------------------------------------------------- Transform

namespace {

/// The scale c J and the weight (c J) a of term q at point i (ESegment::TransformTerms: the total term 1 and 1).
inline void TermWeight(const TransformTerm &t, const Inputs &in, std::size_t i, double &scale, double &weight) {
  if (t.total) {
    scale = weight = 1.0;
    return;
  }
  scale = (t.factor >= 0) ? t.c * in.factors[t.factor][i] : t.c;
  weight = (t.abundance >= 0) ? scale * in.params[t.abundance] : scale;
}

/// F of the term values v (ESegment::TransformedValue's `transform`).
inline double ApplyF(const TransformSpec &t, const double *w, const double *v, double n) {
  double sum = 0.0, total = 0.0;
  for (std::size_t q = 0; q < t.terms.size(); q++) (t.terms[q].total ? total : sum) += w[q] * v[q];
  return (t.kind == TransformSpec::Kind::Yield) ? node::YieldY0(n, sum, total) : node::Transmission(n, sum);
}

/// Row r of a sparse map applied to x (s = 0; s += w x, in the row's order).
inline double RowDot(const SparseRows &m, std::size_t r, const double *x) {
  double s = 0.0;
  for (std::size_t k = m.rowStart[r]; k < m.rowStart[r + 1]; k++) s += m.w[k] * x[m.col[k]];
  return s;
}

}  // namespace

double Tree::TransformValue(const TransformSpec &t, const Inputs &in, std::size_t i, Workspace &ws) const {
  const std::size_t nt = t.terms.size();
  const double n = in.params[t.thickness];
  std::vector<double> &w = ws.c;
  std::vector<double> &v = ws.v;
  w.resize(nt);
  v.resize(nt);
  double scale;
  for (std::size_t q = 0; q < nt; q++) TermWeight(t.terms[q], in, i, scale, w[q]);

  const SubPointFold *fold = (t.fold >= 0) ? in.folds[t.fold] : nullptr;
  const TransformStage st = fold ? fold->stage(i) : TransformStage{};
  if (!st.subPoints) {
    for (std::size_t q = 0; q < nt; q++) v[q] = in.leaves[t.terms[q].input][i];
    return ApplyF(t, w.data(), v.data(), n);
  }
  const std::size_t nsp = st.nsp;
  for (std::size_t q = 0; q < nt; q++)
    if (in.subLeaves[t.terms[q].input].count(i) != nsp) return kNaN;  // refused at load (EData::CheckTransformedSegments)
  std::vector<double> &f = ws.f;
  if (st.premap.rowStart) {
    // a linear map of every term's sub-point values first (free-gas Doppler), F at each row, then the kernel
    const std::size_t rows = st.premap.rows;
    auto rowValues = [&](std::size_t r) {
      for (std::size_t q = 0; q < nt; q++) v[q] = RowDot(st.premap, r, in.subLeaves[t.terms[q].input].at(i));
    };
    if (!st.fold) {
      rowValues(0);
      return ApplyF(t, w.data(), v.data(), n);
    }
    f.resize(rows);
    for (std::size_t r = 0; r < rows; r++) {
      rowValues(r);
      f[r] = ApplyF(t, w.data(), v.data(), n);
    }
    return fold->fold(i, f.data(), rows);
  }
  f.resize(nsp);
  for (std::size_t j = 0; j < nsp; j++) {
    for (std::size_t q = 0; q < nt; q++) v[q] = in.subLeaves[t.terms[q].input].at(i)[j];
    f[j] = ApplyF(t, w.data(), v.data(), n);
  }
  return fold->fold(i, f.data(), nsp);
}

bool Tree::TransformBackward(const TransformSpec &t, const Inputs &in, std::size_t i, double g, Cotangents &cot,
                             Workspace &) const {
  if (g == 0.0) return true;
  const std::size_t nt = t.terms.size();
  const double n = in.params[t.thickness];
  const bool yield = t.kind == TransformSpec::Kind::Yield;
  std::vector<double> scale(nt), w(nt), v(nt);
  for (std::size_t q = 0; q < nt; q++) TermWeight(t.terms[q], in, i, scale[q], w[q]);

  // AZUREGrad's backProp: at kernel weight c and term values v, every term's cotangent (fitBar c) dF/d(its sum)
  // times its weight, the abundances (fitBar c dF/ds) scale v, the thickness (fitBar c) dF/dn.
  auto backProp = [&](double c, auto deliver) {
    double sum = 0.0, total = 0.0, dS = 0.0, dT = 0.0, dN = 0.0;
    for (std::size_t q = 0; q < nt; q++) (t.terms[q].total ? total : sum) += w[q] * v[q];
    if (yield)
      node::YieldY0(n, sum, total, &dS, &dT, &dN);
    else
      node::Transmission(n, sum, &dS, &dN);
    if (cot.params) cot.params[t.thickness] += g * c * dN;
    for (std::size_t q = 0; q < nt; q++) {
      const TransformTerm &term = t.terms[q];
      const double d = g * c * (term.total ? dT : dS);
      if (term.abundance >= 0 && !term.total && cot.params) cot.params[term.abundance] += d * scale[q] * v[q];
      const double cq = d * w[q];
      if (cq != 0.0) deliver(q, cq);
    }
  };

  const SubPointFold *fold = (t.fold >= 0) ? in.folds[t.fold] : nullptr;
  const TransformStage st = fold ? fold->stage(i) : TransformStage{};
  if (!st.subPoints) {
    for (std::size_t q = 0; q < nt; q++) v[q] = in.leaves[t.terms[q].input][i];
    backProp(1.0, [&](std::size_t q, double cq) {
      if (double *p = cot.leaves[t.terms[q].input]) p[i] += cq;
    });
    return true;
  }
  const std::size_t nsp = st.nsp;
  for (std::size_t q = 0; q < nt; q++)
    if (in.subLeaves[t.terms[q].input].count(i) != nsp) return false;
  auto subCot = [&](std::size_t q) -> double * {
    double *p = cot.subLeaves[t.terms[q].input];
    return p ? p + in.subLeaves[t.terms[q].input].offset[i] : nullptr;
  };

  if (st.premap.rowStart) {
    // Doppler before the transform: Y = sum_r c_r F(D_r sigma), c_r the kernel's weights after the map (or 1 for
    // the one row of a line without a kernel); term q's sub-point k gets sum_r cq_r D_rk.
    std::vector<std::pair<std::size_t, double>> stages;
    if (!st.fold) {
      stages.emplace_back(0, 1.0);
    } else {
      std::vector<double> c(st.premap.rows);
      if (!fold->weights(i, c.data(), c.size())) return false;
      for (std::size_t r = 0; r < c.size(); r++)
        if (c[r] != 0.0) stages.emplace_back(r, c[r]);
    }
    std::vector<double> cr(nt * nsp, 0.0);
    for (const auto &sr : stages) {
      const std::size_t r = sr.first;
      for (std::size_t q = 0; q < nt; q++) v[q] = RowDot(st.premap, r, in.subLeaves[t.terms[q].input].at(i));
      backProp(sr.second, [&](std::size_t q, double cq) {
        for (std::size_t k = st.premap.rowStart[r]; k < st.premap.rowStart[r + 1]; k++)
          cr[q * nsp + st.premap.col[k]] += cq * st.premap.w[k];
        if (st.premapDT.rowStart && st.dopplerParam >= 0 && cot.params)
          cot.params[st.dopplerParam] += cq * RowDot(st.premapDT, r, in.subLeaves[t.terms[q].input].at(i));
      });
    }
    for (std::size_t q = 0; q < nt; q++)
      if (double *p = subCot(q))
        for (std::size_t k = 0; k < nsp; k++)
          if (cr[q * nsp + k] != 0.0) p[k] += cr[q * nsp + k];
    return true;
  }
  std::vector<double> c(nsp);
  if (!fold->weights(i, c.data(), nsp)) return false;  // not a pure linear kernel
  for (std::size_t j = 0; j < nsp; j++) {
    if (c[j] == 0.0) continue;
    for (std::size_t q = 0; q < nt; q++) v[q] = in.subLeaves[t.terms[q].input].at(i)[j];
    backProp(c[j], [&](std::size_t q, double cq) {
      if (double *p = subCot(q)) p[j] += cq;
    });
  }
  return true;
}

// ---------------------------------------------------------------- evaluation

void Tree::Sweep(const Inputs &in, std::size_t i0, std::size_t n, double *rows, Workspace &ws) const {
  for (std::size_t k = 0; k < nodes_.size(); k++) {
    const Node &nd = nodes_[k];
    double *row = rows + k * n;
    const int *ch = children_.data() + nd.first;
    switch (nd.kind) {
      case NodeKind::Leaf: {
        const PointArray &a = in.leaves[nd.input];
        for (std::size_t p = 0; p < n; p++) row[p] = a[i0 + p];
        break;
      }
      case NodeKind::Sum: {
        if (nd.count == 0) {
          for (std::size_t p = 0; p < n; p++) row[p] = 0.0;
          break;
        }
        const double *c0 = rows + ch[0] * n;
        for (std::size_t p = 0; p < n; p++) row[p] = c0[p];
        for (int m = 1; m < nd.count; m++) {
          const double *cm = rows + ch[m] * n;
          for (std::size_t p = 0; p < n; p++) row[p] += cm[p];
        }
        break;
      }
      case NodeKind::Scale: {
        const double *c0 = rows + ch[0] * n;
        if (nd.factor >= 0) {
          const PointArray &J = in.factors[nd.factor];
          for (std::size_t p = 0; p < n; p++) row[p] = node::ScaleValue(nd.c, J[i0 + p], c0[p]);
        } else {
          for (std::size_t p = 0; p < n; p++) row[p] = node::ScaleValue(nd.c, 1.0, c0[p]);
        }
        break;
      }
      case NodeKind::Weighted: {
        const double *c0 = rows + ch[0] * n;
        const double a = in.params[nd.param];
        for (std::size_t p = 0; p < n; p++) row[p] = node::WeightedValue(a, c0[p]);
        break;
      }
      case NodeKind::Ratio: {
        const double *cn = rows + ch[0] * n, *cd = rows + ch[1] * n;
        if (nd.factor >= 0) {
          const PointArray &Jn = in.factors[nd.factor], &Jd = in.factors[nd.factor2];
          for (std::size_t p = 0; p < n; p++) row[p] = node::RatioValue(cn[p], cd[p], Jn[i0 + p], Jd[i0 + p]);
        } else {
          for (std::size_t p = 0; p < n; p++) row[p] = node::RatioValue(cn[p], cd[p], 1.0, 1.0);
        }
        break;
      }
      case NodeKind::Transform: {
        const TransformSpec &t = transforms_[nd.input];
        for (std::size_t p = 0; p < n; p++) row[p] = TransformValue(t, in, i0 + p, ws);
        break;
      }
    }
  }
}

void Tree::Evaluate(const Inputs &in, std::size_t n, double *out, Workspace &ws) const {
  if (nodes_.empty()) throw std::logic_error("experiment::Tree: empty tree");
  ws.rows.resize(nodes_.size() * n);
  Sweep(in, 0, n, ws.rows.data(), ws);
  const double *root = ws.rows.data() + static_cast<std::size_t>(Root()) * n;
  for (std::size_t p = 0; p < n; p++) out[p] = root[p];
}

double Tree::Value(const Inputs &in, std::size_t i, Workspace &ws) const {
  if (nodes_.empty()) throw std::logic_error("experiment::Tree: empty tree");
  ws.rows.resize(nodes_.size());
  Sweep(in, i, 1, ws.rows.data(), ws);
  return ws.rows[Root()];
}

bool Tree::Backward(const Inputs &in, std::size_t n, const double *g, Cotangents &cot, Workspace &ws) const {
  if (nodes_.empty()) throw std::logic_error("experiment::Tree: empty tree");
  ws.rows.resize(nodes_.size() * n);
  Sweep(in, 0, n, ws.rows.data(), ws);
  ws.bar.assign(nodes_.size() * n, 0.0);
  const double *rows = ws.rows.data();
  double *bar = ws.bar.data();
  for (std::size_t p = 0; p < n; p++) bar[static_cast<std::size_t>(Root()) * n + p] = g[p];
  bool ok = true;
  for (std::size_t k = nodes_.size(); k-- > 0;) {
    const Node &nd = nodes_[k];
    const double *b = bar + k * n;
    const int *ch = children_.data() + nd.first;
    switch (nd.kind) {
      case NodeKind::Leaf: {
        double *p = cot.leaves[nd.input];
        if (p)
          for (std::size_t q = 0; q < n; q++)
            if (b[q] != 0.0) p[q] += b[q];
        break;
      }
      case NodeKind::Sum:
        for (int m = 0; m < nd.count; m++) {
          double *cb = bar + ch[m] * n;
          for (std::size_t q = 0; q < n; q++)
            if (b[q] != 0.0) cb[q] += node::SumBackward(b[q]);
        }
        break;
      case NodeKind::Scale: {
        double *cb = bar + ch[0] * n;
        for (std::size_t q = 0; q < n; q++) {
          if (b[q] == 0.0) continue;
          const double J = (nd.factor >= 0) ? in.factors[nd.factor][q] : 1.0;
          cb[q] += node::ScaleBackward(b[q], nd.c, J);
        }
        break;
      }
      case NodeKind::Weighted: {
        double *cb = bar + ch[0] * n;
        const double *cv = rows + ch[0] * n;
        const double a = in.params[nd.param];
        for (std::size_t q = 0; q < n; q++) {
          if (b[q] == 0.0) continue;
          if (cot.params) cot.params[nd.param] += node::WeightedBackwardParam(b[q], cv[q]);
          cb[q] += node::WeightedBackwardChild(b[q], a);
        }
        break;
      }
      case NodeKind::Ratio: {
        double *nb = bar + ch[0] * n, *db = bar + ch[1] * n;
        const double *nv = rows + ch[0] * n, *dv = rows + ch[1] * n, *mv = rows + k * n;
        for (std::size_t q = 0; q < n; q++) {
          if (b[q] == 0.0) continue;
          const double jn = (nd.factor >= 0) ? in.factors[nd.factor][q] : 1.0;
          const double jd = (nd.factor2 >= 0) ? in.factors[nd.factor2][q] : 1.0;
          double cn, cd;
          node::RatioBackward(b[q], nv[q], dv[q], jn, jd, mv[q], cn, cd);
          nb[q] += cn;
          db[q] += cd;
        }
        break;
      }
      case NodeKind::Transform: {
        const TransformSpec &t = transforms_[nd.input];
        for (std::size_t q = 0; q < n; q++) ok = TransformBackward(t, in, q, b[q], cot, ws) && ok;
        break;
      }
    }
  }
  return ok;
}

void Tree::Natural(const Inputs &in, std::size_t i, double *perGroup, Workspace &ws) const {
  if (nodes_.empty()) throw std::logic_error("experiment::Tree: empty tree");
  ws.rows.resize(nodes_.size());
  Sweep(in, i, 1, ws.rows.data(), ws);
  for (std::size_t k = 0; k < nodes_.size(); k++) {
    const Node &nd = nodes_[k];
    if (nd.kind == NodeKind::Weighted && nd.group >= 0) perGroup[nd.group] += ws.rows[k];
    if (nd.kind == NodeKind::Transform) {
      for (const TransformTerm &term : transforms_[nd.input].terms) {
        if (term.total || term.group < 0) continue;
        double scale, weight;
        TermWeight(term, in, i, scale, weight);
        perGroup[term.group] += weight * in.leaves[term.input][i];
      }
    }
  }
}

}  // namespace experiment
}  // namespace azure2
