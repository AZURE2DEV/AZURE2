#ifndef AZURE2_EXPERIMENT_ADAPTIVEGRID_H
#define AZURE2_EXPERIMENT_ADAPTIVEGRID_H

#include <functional>
#include <vector>

namespace azure2 {
namespace experiment {

/*!
 * The refinement of an adaptive extrapolation grid (docs/upgrade/project_format3.md section 3.8), as a pure function
 * of a model given as a callable: the engine's extrapolation run passes its point pass (src/AdaptiveExtrapolation.cpp),
 * a test a closed form.
 *
 * The grid starts as `start` (sorted, duplicates removed).  Each sweep tests every interval [a, b] made in the sweep
 * before (all intervals in the first): the model is evaluated at the midpoints m = (a + b)/2 of all of them in one
 * batch, and an interval is split -- m joins the grid and [a, m], [m, b] are tested in the next sweep -- when, for any
 * component c of the model's value,
 *
 *     |f_c(m) - (f_c(a) + f_c(b))/2| > tolerance * max(|f_c(m)|, floor_c),   floor_c = kFloorFraction * max |f_c|,
 *
 * the maximum over the values computed so far.  The floor is the absolute scale below which the tolerance is not
 * relative: values smaller than kFloorFraction of the largest are refined to an absolute deviation of
 * tolerance * kFloorFraction * max |f_c|.  The midpoint of an interval that is not split is not kept, so every
 * interval of the final grid has been tested.  A component that is not finite at a, b or m is not tested; an interval
 * shorter than 1e-12 of its energy is not split.  The refinement stops when no interval is split, or when the splits
 * of a sweep would make the grid longer than `maxPoints`: then the splits with the largest deviations are made up to
 * that count (equal deviations: lower energy first), and the largest deviation of the splits left out is reported.
 *
 * Which intervals are split is decided once per sweep from the values of the batch, so the grid does not depend on
 * how the callable distributes a batch over threads, only on the values it returns.
 */
struct AdaptiveGridResult {
  std::vector<double> x;                    ///< the refined grid, ascending
  std::vector<std::vector<double>> values;  ///< the model's values at x (values[i] has the callable's components)
  bool budgetReached = false;               ///< stopped by maxPoints with intervals still to split
  double largestDeviation = 0.0;            ///< when budgetReached: the largest relative deviation left
  double largestDeviationAt = 0.0;          ///< and the midpoint where it was found
  int sweeps = 0;                           ///< batches of midpoints evaluated
  long evaluations = 0;                     ///< model values computed (start grid and midpoints)
};

/// The values of the model at a batch of abscissae: one vector of components per abscissa, all of the same length.
using GridBatchModel = std::function<std::vector<std::vector<double>>(const std::vector<double> &)>;

/// The floor of the relative tolerance, as a fraction of the largest |value| of a component (see above).
constexpr double kFloorFraction = 1e-3;

AdaptiveGridResult RefineGrid(std::vector<double> start, double tolerance, long maxPoints, const GridBatchModel &model);

}  // namespace experiment
}  // namespace azure2

#endif
