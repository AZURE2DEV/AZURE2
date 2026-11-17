// inference_term_list: the rules of an objective's sums (azure2/inference/TermList.h) on terms that add known
// numbers.  The order of the list is the order of the floating-point additions: the numbers are chosen so that
// another order, or a breakdown summed as the plain sum, gives another double.
#include "azure2/inference/TermList.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "check.h"

using namespace azure2::inference;

namespace {

// The context of these terms: a scale every term multiplies its number with, and a point to exclude.
struct Context {
  double scale = 1.0;
  bool outside = false;
};

// Adds value * scale; `steps` > 1 adds it in that many equal steps (as the engine's data term adds per segment).
class Number : public Term<Context> {
 public:
  Number(std::string name, double value, int steps = 1, bool modelOnly = false)
      : name_(std::move(name)), value_(value), steps_(steps), modelOnly_(modelOnly) {}
  std::string Name() const override { return name_; }
  void Add(const Context &c, double &chi2) const override {
    for (int i = 0; i < steps_; i++) chi2 += value_ * c.scale / steps_;
  }
  bool ModelOnly() const override { return modelOnly_; }

 private:
  std::string name_;
  double value_;
  int steps_;
  bool modelOnly_;
};

// A prior the sampler takes as a normalised log-prior elsewhere: nothing here.
class SamplerPrior : public Number {
 public:
  using Number::Number;
  void AddForSampler(const Context &, double &) const override {}
};

// A support: adds a wall in the fit, nothing in the sampler, and excludes the point when the context says so.
class Support : public Number {
 public:
  using Number::Number;
  void AddForSampler(const Context &, double &) const override {}
  bool ExcludesForSampler(const Context &c) const override { return c.outside; }
};

}  // namespace

int main() {
  TermList<Context> list;
  check("an empty list sums to 0", list.Sum(Context()) == 0.0 && list.SumForSampler(Context()) == 0.0 &&
                                       list.SumOfModelOnly(Context()) == 0.0 && !list.ExcludesForSampler(Context()));
  // 1e16 + 1 + 1 in this order is 1e16 (each 1 is lost); 1 + 1 + 1e16 would be 1e16 + 2.
  list.Append(std::make_unique<Number>("data", 1.0e16));
  list.Append(std::make_unique<SamplerPrior>("prior", 1.0));
  list.Append(std::make_unique<Support>("wall", 1.0, 1, true));
  list.Append(std::make_unique<Number>("model data", 3.0, 3, true));
  check("the names, in the order of the appends",
        list.Terms().size() == 4 && list.Terms()[0]->Name() == "data" && list.Terms()[1]->Name() == "prior" &&
            list.Terms()[2]->Name() == "wall" && list.Terms()[3]->Name() == "model data");

  const Context c;
  // Running sum: 1e16, +1 -> 1e16, +1 -> 1e16, then three steps of 1 -> 1e16 each (ulp of 1e16 is 2, ties to even).
  check("Sum adds into one running sum, in order", list.Sum(c) == 1.0e16);
  std::vector<double> parts = {99.0};
  const double total = list.SumWithParts(c, parts);
  check("parts: each term on its own", parts.size() == 4 && parts[0] == 1.0e16 && parts[1] == 1.0 && parts[2] == 1.0 && parts[3] == 3.0);
  // 1e16 + 1 -> 1e16; + 1 -> 1e16; + 3 -> 1e16 + 4 (3 above an even neighbour rounds up): not the plain Sum.
  check("SumWithParts sums the parts in order", total == ((1.0e16 + 1.0) + 1.0) + 3.0 && total != list.Sum(c));
  check("SumOfFirst(1) is the first term", list.SumOfFirst(c, 1) == 1.0e16);
  check("SumOfFirst(2) is the first two, SumOfFirst past the end the whole list",
        list.SumOfFirst(c, 2) == 1.0e16 + 1.0 && list.SumOfFirst(c, 9) == list.Sum(c) && list.SumOfFirst(c, 0) == 0.0);
  check("model-only terms: the wall and the model data", list.SumOfModelOnly(c) == 4.0);
  // The sampler: the default is the whole term (data, model data); the prior and the wall add nothing.
  check("the sampler's sum on this list: the data and the model data", list.SumForSampler(c) == 1.0e16);

  // A list of binary fractions, where every sum is exact: 0.5 (data), 0.25 (a prior the sampler normalises itself),
  // 0.125 (a support, model only), 2 (model data).
  TermList<Context> exact;
  exact.Append(std::make_unique<Number>("data", 0.5));
  exact.Append(std::make_unique<SamplerPrior>("prior", 0.25));
  exact.Append(std::make_unique<Support>("wall", 0.125, 1, true));
  exact.Append(std::make_unique<Number>("model data", 2.0, 2, true));
  check("Sum = 2.875", exact.Sum(c) == 2.875);
  check("the sampler's sum leaves out the prior and the wall: 2.5", exact.SumForSampler(c) == 2.5);
  check("model only: 2.125", exact.SumOfModelOnly(c) == 2.125);
  check("the first two: 0.75", exact.SumOfFirst(c, 2) == 0.75);
  Context doubled;
  doubled.scale = 2.0;
  check("the context reaches every term: 5.75", exact.Sum(doubled) == 5.75 && exact.SumWithParts(doubled, parts) == 5.75 &&
                                                    parts.size() == 4 && parts[1] == 0.5);
  Context out;
  out.outside = true;
  check("a term that excludes", list.ExcludesForSampler(out) && !list.ExcludesForSampler(c));
  return result();
}
