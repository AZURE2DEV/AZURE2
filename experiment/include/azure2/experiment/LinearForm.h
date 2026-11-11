// The linear form on which the apparatus kernels (Kernels.h) are evaluated to obtain their weights.
#ifndef AZURE2_EXPERIMENT_LINEARFORM_H
#define AZURE2_EXPERIMENT_LINEARFORM_H

#include <utility>
#include <vector>

namespace azure2 {
namespace experiment {

/*!
 * A linear form sum_k c_k s_{j_k} in the node values s_j of a point's kernel.  Every kernel of Kernels.h is linear in
 * those values, so evaluating it on this type instead of double yields its exact weights in one pass (the gradient's
 * FoldWeights in the engine).  Only the operations the kernels perform are provided; a nonzero constant marks the
 * form invalid.
 */
class LinearForm {
 public:
  LinearForm(double v = 0.0) : valid_(v == 0.0) {}
  static LinearForm Unit(int j) {
    LinearForm f;
    f.terms_.emplace_back(j, 1.0);
    return f;
  }
  LinearForm &operator*=(double a) {
    for (auto &t : terms_) t.second *= a;
    return *this;
  }
  LinearForm &operator/=(double a) {
    for (auto &t : terms_) t.second /= a;
    return *this;
  }
  LinearForm &operator+=(const LinearForm &o) {
    terms_.insert(terms_.end(), o.terms_.begin(), o.terms_.end());
    valid_ = valid_ && o.valid_;
    return *this;
  }
  friend LinearForm operator*(LinearForm f, double a) { return f *= a; }
  friend LinearForm operator*(double a, LinearForm f) { return f *= a; }
  friend LinearForm operator/(LinearForm f, double a) { return f /= a; }
  friend LinearForm operator+(LinearForm f, const LinearForm &o) { return f += o; }
  bool Valid() const { return valid_; }
  const std::vector<std::pair<int, double>> &Terms() const { return terms_; }

 private:
  std::vector<std::pair<int, double>> terms_;
  bool valid_;
};

}  // namespace experiment
}  // namespace azure2

#endif
