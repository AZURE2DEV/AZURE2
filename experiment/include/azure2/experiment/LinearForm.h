// The linear form on which the apparatus kernels (Kernels.h) are evaluated to obtain their weights.
#ifndef AZURE2_EXPERIMENT_LINEARFORM_H
#define AZURE2_EXPERIMENT_LINEARFORM_H

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace azure2 {
namespace experiment {

/*!
 * A linear form sum_k c_k s_{j_k} in the node values s_j of a point's kernel.  Every kernel of Kernels.h is linear in
 * those values, so evaluating it on this type instead of double yields its exact weights in one pass (the gradient's
 * FoldWeights in the engine).  Only the operations the kernels perform are provided; a nonzero constant marks the
 * form invalid.
 *
 * Cost.  The kernels build a form as they build a double: an interpolated value (two terms), scaled, added to a
 * running sum, once per Gauss point.  Two representations keep that at the cost of the double arithmetic:
 *  - up to kInline terms live in the object, so the temporaries of an interpolation never touch the heap;
 *  - a larger form is a dense run of coefficients over the nodes [low, low + size): a running sum over a window of
 *    nodes holds one coefficient per node, and adding a form to it is one pass over that form's own run.
 * (A plain list of terms, appended to at every addition, made the sum of the target kernel with a beam spread or
 * straggling hold O(nodes x window) terms: 3e7 terms, 0.5 GB, per point of 11B+a on the shared mesh.)
 */
class LinearForm {
 public:
  using Term = std::pair<int, double>;
  LinearForm(double v = 0.0) : valid_(v == 0.0) {}
  static LinearForm Unit(int j) {
    LinearForm f;
    f.small_[f.count_++] = Term(j, 1.0);
    return f;
  }
  LinearForm &operator*=(double a) {
    for (int i = 0; i < count_; i++) small_[i].second *= a;
    for (double &c : dense_) c *= a;
    return *this;
  }
  LinearForm &operator/=(double a) {
    for (int i = 0; i < count_; i++) small_[i].second /= a;
    for (double &c : dense_) c /= a;
    return *this;
  }
  LinearForm &operator+=(const LinearForm &o) {
    valid_ = valid_ && o.valid_;
    if (o.dense_.empty()) {
      if (dense_.empty() && count_ + o.count_ <= kInline) {
        for (int i = 0; i < o.count_; i++) small_[count_++] = o.small_[i];
        return *this;
      }
      if (o.count_ == 0) return *this;
      int low = o.small_[0].first, high = low;
      for (int i = 1; i < o.count_; i++) low = std::min(low, o.small_[i].first), high = std::max(high, o.small_[i].first);
      Cover(low, high);
      for (int i = 0; i < o.count_; i++) dense_[(std::size_t)(o.small_[i].first - low_)] += o.small_[i].second;
      return *this;
    }
    Cover(o.low_, o.low_ + (int)o.dense_.size() - 1);
    double *c = dense_.data() + (o.low_ - low_);
    for (std::size_t i = 0; i < o.dense_.size(); i++) c[i] += o.dense_[i];
    return *this;
  }
  /*!
   * num += c (fa a + fb b), the quadrature's step (KernelQuadrature.h, AddInterpolated).  For a and b single terms
   * (the unit forms of two nodes: every kernel that folds the node values themselves) the two products are added
   * straight into num, with the arithmetic and the order of the general expression -- (coefficient fa) c, then num's
   * own += of a two-term form -- so the weights are the same to the last bit, without a temporary form per Gauss point.
   */
  friend void AddInterpolated(LinearForm &num, double c, double fa, const LinearForm &a, double fb, const LinearForm &b) {
    if (a.count_ != 1 || b.count_ != 1 || !a.dense_.empty() || !b.dense_.empty()) {
      LinearForm sigma = fa * a + fb * b;
      num += c * sigma;
      return;
    }
    num.valid_ = num.valid_ && a.valid_ && b.valid_;
    const int ia = a.small_[0].first, ib = b.small_[0].first;
    const double va = (a.small_[0].second * fa) * c, vb = (b.small_[0].second * fb) * c;
    if (num.dense_.empty() && num.count_ == 0) {
      num.small_[0] = Term(ia, va);
      num.small_[1] = Term(ib, vb);
      num.count_ = 2;
      return;
    }
    num.Cover(std::min(ia, ib), std::max(ia, ib));
    num.dense_[(std::size_t)(ia - num.low_)] += va;
    num.dense_[(std::size_t)(ib - num.low_)] += vb;
  }
  friend LinearForm operator*(LinearForm f, double a) { return f *= a; }
  friend LinearForm operator*(double a, LinearForm f) { return f *= a; }
  friend LinearForm operator/(LinearForm f, double a) { return f /= a; }
  friend LinearForm operator+(LinearForm f, const LinearForm &o) { return f += o; }
  bool Valid() const { return valid_; }
  /// The terms (node, coefficient).  A node may appear more than once: its coefficient is the sum.
  std::vector<Term> Terms() const {
    std::vector<Term> terms(small_, small_ + count_);
    for (std::size_t i = 0; i < dense_.size(); i++)
      if (dense_[i] != 0.0) terms.emplace_back(low_ + (int)i, dense_[i]);
    return terms;
  }

 private:
  enum { kInline = 2 };
  /// Makes the form dense over (at least) the nodes [low, high], its inline terms included.
  void Cover(int low, int high) {
    for (int i = 0; i < count_; i++) low = std::min(low, small_[i].first), high = std::max(high, small_[i].first);
    if (dense_.empty()) {
      low_ = low;
      dense_.assign((std::size_t)(high - low + 1), 0.0);
    } else {
      if (low < low_) {
        // grows at the low end by at least its own size, so that a sum built towards lower nodes stays linear
        const int grow = std::max(low_ - low, (int)dense_.size());
        dense_.insert(dense_.begin(), (std::size_t)grow, 0.0);
        low_ -= grow;
      }
      if (high >= low_ + (int)dense_.size()) dense_.resize((std::size_t)(high - low_ + 1), 0.0);
    }
    for (int i = 0; i < count_; i++) dense_[(std::size_t)(small_[i].first - low_)] += small_[i].second;
    count_ = 0;
  }
  Term small_[kInline];
  int count_ = 0;              ///< terms in small_ (the form is small_ plus dense_; dense_ non-empty: count_ = 0)
  int low_ = 0;                ///< node of dense_[0]
  std::vector<double> dense_;  ///< coefficients of the nodes low_, low_ + 1, ...
  bool valid_;
};

}  // namespace experiment
}  // namespace azure2

#endif
