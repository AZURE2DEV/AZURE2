#ifndef AZURE2_INFERENCE_TERMLIST_H
#define AZURE2_INFERENCE_TERMLIST_H

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace azure2 {
namespace inference {

/*!
 * One term of an objective chi2 = sum of terms, over the context `Context` its evaluation reads (the engine's
 * ObjectiveContext: the compound system, the level state, the working copy, the caller's results buffer, the vector).
 * A term adds into the running sum, so that a list summed in order gives one floating-point result whatever
 * path asks for it (docs/upgrade/objective.md).
 */
template <class Context>
class Term {
 public:
  virtual ~Term() = default;
  virtual std::string Name() const = 0;
  /// Adds this term to chi2 (some terms add in several steps, in the order the engine always used).
  virtual void Add(const Context &c, double &chi2) const = 0;
  /// A term of the model alone (no data, no parameter vector).
  virtual bool ModelOnly() const { return false; }
  /*!
   * The part of the term the ensemble sampler adds to its -2 ln L, at the model the evaluation filled
   * (docs/upgrade/objective.md section 6).  The default is the whole term, so that a new term reaches the sampler
   * with no edit of the sampler.  A term whose parameters the sampler gives normalised Gaussian priors from their
   * enumeration adds that part here no more, and says so; a term that is a support of the posterior adds nothing and
   * answers ExcludesForSampler.
   */
  virtual void AddForSampler(const Context &c, double &chi2) const { Add(c, chi2); }
  /// True where the sampler's posterior is zero (-infinity in ln P) by this term.
  virtual bool ExcludesForSampler(const Context &) const { return false; }
};

/*!
 * The ordered list of an objective's terms and the rules of its sums.  The order of the list is the order of the
 * floating-point additions; it is part of the result.  The list holds no state of an evaluation: every sum is a
 * const function of the context.
 */
template <class Context>
class TermList {
 public:
  using TermPtr = std::unique_ptr<Term<Context>>;

  /// Appends a term: the list's order is the order of the calls.
  void Append(TermPtr term) { terms_.push_back(std::move(term)); }
  const std::vector<TermPtr> &Terms() const { return terms_; }

  /// Every term added into one running sum, in order.
  double Sum(const Context &c) const {
    double chi2 = 0.0;
    for (const TermPtr &t : terms_) t->Add(c, chi2);
    return chi2;
  }
  /// The first `count` terms added into one running sum, in order.
  double SumOfFirst(const Context &c, std::size_t count) const {
    double chi2 = 0.0;
    for (std::size_t i = 0; i < count && i < terms_.size(); i++) terms_[i]->Add(c, chi2);
    return chi2;
  }
  /// With a breakdown: each term on its own into `parts` (a small prior next to a large data chi^2 keeps its
  /// digits), then summed in order.  The total may differ from Sum in its last bits.
  double SumWithParts(const Context &c, std::vector<double> &parts) const {
    double chi2 = 0.0;
    parts.clear();
    for (const TermPtr &t : terms_) {
      double v = 0.0;
      t->Add(c, v);
      parts.push_back(v);
      chi2 += v;
    }
    return chi2;
  }
  /// Every term's AddForSampler into one running sum, in order.
  double SumForSampler(const Context &c) const {
    double chi2 = 0.0;
    for (const TermPtr &t : terms_) t->AddForSampler(c, chi2);
    return chi2;
  }
  /// Whether a term puts c's point outside the sampler's support.
  bool ExcludesForSampler(const Context &c) const {
    for (const TermPtr &t : terms_)
      if (t->ExcludesForSampler(c)) return true;
    return false;
  }
  /// The terms of the model alone (ModelOnly), added into one running sum, in order.
  double SumOfModelOnly(const Context &c) const {
    double chi2 = 0.0;
    for (const TermPtr &t : terms_)
      if (t->ModelOnly()) t->Add(c, chi2);
    return chi2;
  }

 private:
  std::vector<TermPtr> terms_;
};

}  // namespace inference
}  // namespace azure2

#endif  // AZURE2_INFERENCE_TERMLIST_H
