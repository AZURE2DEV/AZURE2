#ifndef AZURE2_EXPERIMENT_EXPRESSION_H
#define AZURE2_EXPERIMENT_EXPRESSION_H

#include <cstddef>
#include <cstdint>
#include <vector>

namespace azure2 {
namespace experiment {

/*!
 * Model expressions (docs/upgrade/data_model.md section 2.3): what a data set measures, as a tree whose leaves are
 * calculated observables (a cross section in b or b/sr, c.m., or another quantity in its own unit) and whose nodes
 * combine them point by point.  The model value m_i of measured point i is the value of the root at i.
 *
 *   Leaf       v = the leaf's value at point i (read from the caller's arrays; in the engine a calculation slot)
 *   Sum        v = v_0 + v_1 + ... + v_{n-1}, left to right from v_0 (0 for no children)
 *   Scale      v = (c J_i) v_child, c a constant, J_i an optional per-point factor (the frame factor J_base/J_k)
 *   Weighted   v = a v_child, a a parameter (the abundance a_{s,n} of a nucleus n in sample s)
 *   Ratio      v = v_num / v_den (J_den,i / J_num,i), NaN where v_den = 0
 *   Transform  T = exp(-n s) (transmission) or Y0 = (s / s_t)(1 - exp(-n s_t)) (SAMMY's self-shielded yield) of the
 *              weighted term sums s, s_t at every sub-point of the base unit, then the base unit's kernel
 *
 * Every node's arithmetic is the engine's, in the engine's order of operations (src/ModelPass.cpp,
 * ComposeModelValue and TransformedValue), so that the tree reproduces today's model values bit for
 * bit; `backward` is the adjoint of each node in the order src/AZUREGrad.cpp uses (GradOnePoint, shared-component
 * route; GradTransformedPoint).  experiment/README.md says where the two differ and why.
 *
 * Storage: a Tree is a flat vector of Node values in creation order.  A node is created after its children, so the
 * creation order is a post-order and the evaluation is one forward sweep over the vector (no recursion, no pointers,
 * no virtual call); the backward sweep runs it in reverse.  The tree holds only structure and constants; the per-point
 * values (leaves, frame factors), the parameters and the apparatus of a Transform are passed in at each evaluation
 * (Inputs), so one tree serves every evaluation and every working copy.
 */

/// A per-point array: the value at point i is data[i * stride].
struct PointArray {
  const double *data = nullptr;
  std::ptrdiff_t stride = 1;
  double operator[](std::size_t i) const { return data[static_cast<std::ptrdiff_t>(i) * stride]; }
};

/// Per-point runs of sub-point values: point i owns data[offset[i]] .. data[offset[i + 1] - 1] (offset has n + 1
/// entries).  The sub-point order is the base unit's (sub-point 1 first, the highest energy).
struct SubPointArray {
  const double *data = nullptr;
  const std::size_t *offset = nullptr;
  std::size_t count(std::size_t i) const { return offset[i + 1] - offset[i]; }
  const double *at(std::size_t i) const { return data + offset[i]; }
};

/// A sparse linear map in compressed rows: row r is sum_k w[k] x[col[k]] over k in [rowStart[r], rowStart[r + 1]).
struct SparseRows {
  std::size_t rows = 0;
  const std::size_t *rowStart = nullptr;  ///< rows + 1 entries
  const int *col = nullptr;                ///< 0-based columns
  const double *w = nullptr;
};

/*!
 * The apparatus of a Transform at one point, as the engine lays it out (TransformedValue in src/ModelPass.cpp).
 *
 * - subPoints = false: no sub-point effect; F is applied to the terms' point values (their own folded values) and
 *   is the model value.
 * - subPoints = true: every term has nsp sub-point values.  The values at which F is evaluated are those sub-point
 *   values themselves (premap.rowStart == nullptr: rows = nsp), or a linear map of them (premap: free-gas Doppler
 *   broadening, which SAMMY applies before the transform).  With fold = true the values F_r are then folded with the
 *   base unit's kernel (SubPointFold::fold); with fold = false there must be one row and F_0 is the model value (a
 *   Doppler line with no kernel after the Doppler step).
 * - premapDT (optional, same rows): d(premap)/dT for a fitted temperature, parameter dopplerParam (backward only).
 */
struct TransformStage {
  bool subPoints = false;
  std::size_t nsp = 0;
  SparseRows premap;
  bool fold = true;
  SparseRows premapDT;
  int dopplerParam = -1;
};

/*!
 * The base unit's kernel, supplied by the caller (the engine: EPoint::IntegrateKernelsImpl with a reader over the
 * values given here, under EPoint::DopplerAppliedScope).  One call per point; the kernel is the caller's.
 */
class SubPointFold {
 public:
  virtual ~SubPointFold() = default;
  /// The layout at point i.
  virtual TransformStage stage(std::size_t i) const = 0;
  /// The kernel applied to the values f[0..n) at point i (n = the stage's rows).
  virtual double fold(std::size_t i, const double *f, std::size_t n) const = 0;
  /// The exact weights c[0..n) of a linear kernel, fold(f) = sum_r c_r f_r (EPoint::TargetEffectWeights); false if
  /// the kernel is not linear (backward then fails).
  virtual bool weights(std::size_t i, double *c, std::size_t n) const = 0;
};

/// One term of a Transform: weight w_i = (c J_i) a (J_i = 1 without a factor array, a = 1 without an abundance);
/// `total` marks the total term s_t of a yield (weight 1, no scaling: ESegment::TransformTerms).
struct TransformTerm {
  int input = -1;      ///< leaf input (its point values in Inputs::leaves, its sub-point values in Inputs::subLeaves)
  double c = 1.0;      ///< the term's scaling
  int factor = -1;     ///< per-point frame factor J (Inputs::factors), or -1
  int abundance = -1;  ///< abundance parameter (Inputs::params), or -1
  int group = -1;      ///< nucleus for the decomposition (natural), or -1
  bool total = false;
};

struct TransformSpec {
  enum class Kind { Transmission, Yield };
  Kind kind = Kind::Transmission;
  int thickness = -1;  ///< sample thickness n (atoms/b), a parameter (Inputs::params)
  int fold = -1;       ///< the kernel (Inputs::folds); needed only where the stage has sub-points
  std::vector<TransformTerm> terms;
};

/// A view of `size` values the caller owns (an array on its stack, a vector it keeps): Inputs and Cotangents refer to
/// the caller's arrays, so that an evaluation at one point allocates nothing.  A view of a vector is valid while the
/// vector is neither destroyed nor resized; a temporary vector is refused.
template <class T>
struct View {
  const T *data = nullptr;
  std::size_t size = 0;
  View() = default;
  View(const T *values, std::size_t count) : data(values), size(count) {}
  View(const std::vector<T> &values) : data(values.data()), size(values.size()) {}  // NOLINT: implicit by design
  View(std::vector<T> &&) = delete;
  const T &operator[](std::size_t k) const { return data[k]; }
};

/// What an evaluation reads.
struct Inputs {
  View<PointArray> leaves;      ///< leaf input k at point i
  View<SubPointArray> subLeaves;  ///< leaf input k at the sub-points of point i (Transform terms only)
  View<PointArray> factors;     ///< per-point factors (Scale, Ratio, Transform terms)
  const double *params = nullptr;  ///< abundances, sample thicknesses, Doppler temperatures
  View<const SubPointFold *> folds;
};

/// Cotangents accumulated by Tree::backward (+=; null entries are skipped).
struct Cotangents {
  View<double *> leaves;     ///< per leaf input, per point
  View<double *> subLeaves;  ///< per leaf input, laid out as its SubPointArray
  double *params = nullptr;  ///< per parameter
};

/// The arithmetic of each node at one point, value and backward, as used by Tree (and tested alone).
namespace node {

/// v_0 + v_1 + ... from v_0 (as ESegment: sum = base; sum += term).  0 for n = 0.
inline double SumValue(const double *v, std::size_t n) {
  if (n == 0) return 0.0;
  double s = v[0];
  for (std::size_t k = 1; k < n; k++) s += v[k];
  return s;
}
/// d(Sum)/d(child) = 1: every child gets g.
inline double SumBackward(double g) { return g; }

/// (c J) v (ESegment: componentScaling * TermFrameFactor * compValue).
inline double ScaleValue(double c, double J, double v) { return c * J * v; }
/// child cotangent g (c J) (AZUREGrad: fitBar * s, s = scaling * TermFrameFactor).
inline double ScaleBackward(double g, double c, double J) { return g * (c * J); }

/// v a (ESegment: baseValue *= a; term *= a).
inline double WeightedValue(double a, double v) { return v * a; }
/// child cotangent g a; parameter cotangent g v (AZUREGrad: fitBar * weight, fitBar * baseVal).
inline double WeightedBackwardChild(double g, double a) { return g * a; }
inline double WeightedBackwardParam(double g, double v) { return g * v; }

/// num / den (J_den / J_num), NaN where den = 0 (ESegment: baseValue / denominatorValue * (dJ / bJ)).
double RatioValue(double num, double den, double jNum, double jDen);
/// Cotangents of numerator and denominator from the value m: g (m / num) and g (-m / den) (AZUREGrad's shared
/// route).  At num = 0 the numerator's is g (J_den / J_num) / den (AZUREGrad gives 0 there); NaN where den = 0.
void RatioBackward(double g, double num, double den, double jNum, double jDen, double m, double &cotNum,
                   double &cotDen);

/// T = exp(-n s) and its partial derivatives (ESegment::Transmission).
double Transmission(double n, double s, double *dTds = nullptr, double *dTdn = nullptr);
/// Y0 = (s / s_t)(1 - exp(-n s_t)) and its partial derivatives (ESegment::YieldY0); below |n s_t| = 1e-4 the series
/// of g(s_t) = (1 - e^{-n s_t}) / s_t, so s_t -> 0 is no 0/0.
double YieldY0(double n, double s, double st, double *dYds = nullptr, double *dYdst = nullptr, double *dYdn = nullptr);

}  // namespace node

enum class NodeKind : std::uint8_t { Leaf, Sum, Scale, Weighted, Ratio, Transform };

/// One node of a Tree.  Children are listed in Tree::children() [first, first + count).
struct Node {
  NodeKind kind = NodeKind::Leaf;
  int first = 0, count = 0;
  int input = -1;     ///< Leaf: leaf input; Transform: index of its TransformSpec
  int factor = -1;    ///< Scale: J (Inputs::factors) or -1; Ratio: J_num
  int factor2 = -1;   ///< Ratio: J_den
  int param = -1;     ///< Weighted: abundance (Inputs::params)
  int group = -1;     ///< Weighted: nucleus for the decomposition, or -1
  double c = 1.0;     ///< Scale: constant
};

/// Scratch space of an evaluation, reused across calls (one per thread).  The node values and cotangents of a tree of
/// up to kInline nodes at one point live in the fixed arrays, so that Value and a one-point Backward of such a tree
/// allocate nothing; larger trees and runs of points use the vectors.
struct Workspace {
  static constexpr std::size_t kInline = 32;
  double inlineRows[kInline];
  double inlineBar[kInline];
  std::vector<double> rows;  // node x point values (Tree::evaluate, backward)
  std::vector<double> bar;   // node x point cotangents (backward)
  std::vector<double> v, f, c, raw, cot;
};

/*!
 * A model expression.  Built once (the builder functions return node ids; a node's children must exist and may be
 * used once), then evaluated for every point of its measurement at every chi^2 evaluation.
 */
class Tree {
 public:
  int Leaf(int input);
  int Sum(const std::vector<int> &children);
  int Scale(double c, int child, int factor = -1);
  int Weighted(int param, int child, int group = -1);
  int Ratio(int num, int den, int jNum = -1, int jDen = -1);
  int Transform(const TransformSpec &spec);

  /// The root is the last node created.
  int Root() const { return static_cast<int>(nodes_.size()) - 1; }
  const std::vector<Node> &Nodes() const { return nodes_; }
  const std::vector<int> &Children() const { return children_; }
  const std::vector<TransformSpec> &Transforms() const { return transforms_; }
  /// Number of leaf inputs, factor arrays and parameters the tree refers to (one more than the largest index).
  int NumInputs() const { return numInputs_; }
  int NumFactors() const { return numFactors_; }
  int NumParams() const { return numParams_; }

  /// m_i for i in [0, n): out[i].  One pass per node over all points.
  void Evaluate(const Inputs &in, std::size_t n, double *out, Workspace &ws) const;
  /// m_i at one point; with `perGroup`, also the per-nucleus terms of that evaluation added to it (Natural).
  double Value(const Inputs &in, std::size_t i, Workspace &ws, double *perGroup = nullptr) const;
  /// Accumulates g_i d m_i / d(inputs) into `cot` for i in [0, n) (reverse sweep; values recomputed).  Returns false
  /// where a Transform's kernel is not linear (SubPointFold::weights false).
  bool Backward(const Inputs &in, std::size_t n, const double *g, Cotangents &cot, Workspace &ws) const;
  /// The per-nucleus terms at point i (the decomposition of a sample): every Weighted node's value added to
  /// perGroup[group] in creation order, and every Transform term's w_i v_i (point values) likewise.
  void Natural(const Inputs &in, std::size_t i, double *perGroup, Workspace &ws) const { Value(in, i, ws, perGroup); }

 private:
  int Add(Node nd, const std::vector<int> &children);
  double TransformValue(const TransformSpec &t, const Inputs &in, std::size_t i, Workspace &ws) const;
  bool TransformBackward(const TransformSpec &t, const Inputs &in, std::size_t i, double g, Cotangents &cot,
                         Workspace &ws) const;
  void Sweep(const Inputs &in, std::size_t i0, std::size_t n, double *rows, Workspace &ws) const;
  void AddNatural(const Inputs &in, std::size_t i, const double *rows, double *perGroup) const;
  /// node x point storage of the workspace: its fixed array when it fits, else the vector, resized.
  double *Rows(Workspace &ws, std::size_t count) const;
  double *Bar(Workspace &ws, std::size_t count) const;

  std::vector<Node> nodes_;
  std::vector<int> children_;
  std::vector<char> used_;
  std::vector<TransformSpec> transforms_;
  int numInputs_ = 0, numFactors_ = 0, numParams_ = 0;
};

/// The frame factor of a SUM term at a point, J_base / J_term, from the two points' lab-to-c.m. factors; 1 where
/// either is 0 or they are equal (ESegment::TermFrameFactor).
inline double TermFrameFactor(double jBase, double jTerm) {
  if (jTerm == 0.0 || jBase == 0.0 || jTerm == jBase) return 1.0;
  return jBase / jTerm;
}

}  // namespace experiment
}  // namespace azure2

#endif
