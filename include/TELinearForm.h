// The linear form on which the apparatus folds of a target-effect point are evaluated to obtain their weights.
#ifndef TELINEARFORM_H
#define TELINEARFORM_H

#include <utility>
#include <vector>

/*!
 * A linear form sum_k c_k s_{j_k} in the sub-point values s_j of a
 * target-effect point.  Unit::IntegrateTargetEffect is linear in those
 * values, so evaluating it on this type instead of double yields its exact
 * weights in one pass (Unit::TargetEffectWeights), which the analytic
 * gradient needs for every target-integrated point.  Only the operations the
 * integrator performs are provided; a nonzero constant marks the form invalid.
 */
class TELinearForm {
 public:
  TELinearForm(double v = 0.0) : valid_(v == 0.0) {}
  static TELinearForm Unit(int j) {
    TELinearForm f;
    f.terms_.emplace_back(j, 1.0);
    return f;
  }
  TELinearForm &operator*=(double a) {
    for (auto &t : terms_) t.second *= a;
    return *this;
  }
  TELinearForm &operator/=(double a) {
    for (auto &t : terms_) t.second /= a;
    return *this;
  }
  TELinearForm &operator+=(const TELinearForm &o) {
    terms_.insert(terms_.end(), o.terms_.begin(), o.terms_.end());
    valid_ = valid_ && o.valid_;
    return *this;
  }
  friend TELinearForm operator*(TELinearForm f, double a) { return f *= a; }
  friend TELinearForm operator*(double a, TELinearForm f) { return f *= a; }
  friend TELinearForm operator/(TELinearForm f, double a) { return f /= a; }
  friend TELinearForm operator+(TELinearForm f, const TELinearForm &o) { return f += o; }
  bool Valid() const { return valid_; }
  const std::vector<std::pair<int, double>> &Terms() const { return terms_; }

 private:
  std::vector<std::pair<int, double>> terms_;
  bool valid_;
};

#endif
