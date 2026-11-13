#ifndef AZURE2_EXPERIMENT_LAYOUTBASE_H
#define AZURE2_EXPERIMENT_LAYOUTBASE_H

#include <string>

#include "azure2/experiment/SetupSpecs.h"

namespace azure2 {
namespace experiment {

/*!
 * The range a fitted nuisance parameter whose value lays out a grid may take (docs/upgrade/parameter_kinds.md
 * section 4): the sub-point grids span its highest value, the fit's limits are [lowest, highest].  Both are laid out
 * from one number, the base: lowest = lowestFraction x base; highest = grid_max where written, else
 * max(base + 3 sigma, widenWithPrior x base) with a prior (sigma > 0), widenWithout x base without one.
 */
struct LayoutSpan {
  double lowestFraction = 0.0;
  double widenWithPrior = 0.0;
  double widenWithout = 0.0;
  double Lowest(double base) const { return lowestFraction * base; }
  /// \p sigma the prior's width (0 = none), \p gridMax the written grid_max (0 = none); same unit as the base.
  double Highest(double base, double sigma, double gridMax) const;
};

/// An areal density: [N/100, max(N + 3 sigma, 1.1 N)] (1.5 N without a prior).
constexpr LayoutSpan kDensitySpan{0.01, 1.1, 1.5};
/// A Doppler temperature: [T/4, max(T + 3 sigma, 1.21 T)] (2 T without a prior); the Gauss-Legendre pieces are laid
/// out for T too.
constexpr LayoutSpan kDopplerSpan{0.25, 1.21, 2.0};
/// An apparatus parameter (the stopping scale, a kernel's parameter): [s/10, max(s + 3 sigma, 1.1 s)] (1.5 s).
constexpr LayoutSpan kApparatusSpan{0.1, 1.1, 1.5};

/// What a fitted parameter's grids and limits are laid out from, and the line that says why when it is not the centre.
struct LayoutBase {
  double value = 0.0;
  std::string note;  ///< not empty when a written centre was passed over (the start lies outside its range)
};

/*!
 * The one rule of every kind: the prior's centre where one is written and the start lies inside the range laid out
 * from it, so that a fit started at the centre and a fit continued from its result (the fitted value written as the
 * start, the prior unchanged) are computed on the same grids with the same limits; otherwise (no centre written, a
 * fixed parameter, or a start outside the centre's range) the start, with a note in the last case.  A start equal to
 * the centre returns the start itself, so a project with value = centre, or without a prior, is laid out as before.
 * \p name and \p unit name the parameter in the note.
 */
LayoutBase ChooseLayoutBase(const LayoutSpan &span, double start, const PriorCentre &centre, bool vary, double sigma,
                            double gridMax, const std::string &name, const std::string &unit);

}  // namespace experiment
}  // namespace azure2

#endif
