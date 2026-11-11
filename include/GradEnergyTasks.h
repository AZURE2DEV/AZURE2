#ifndef AZURE2_GRADENERGYTASKS_H
#define AZURE2_GRADENERGYTASKS_H

#include <functional>
#include <vector>

#include "AMatrixFunc.h"
#include "AZUREGrad.h"
#include "CalculationPlan.h"

/*!
 * The analytic gradient over the energy tasks of the point pass (docs/upgrade/energy_grid.md section 8, second
 * branch; docs/upgrade/ARCHITECTURE.md, "The analytic gradient").  A pass has three steps:
 *
 *   1. forward: the point pass of the data segments (WorkingCopy::CalculatePointsInParallel, one level-matrix solve
 *      per compound energy, then the folds) and each point's model from the records it wrote;
 *   2. cotangents: each point's cotangent of its model (the objective's, fitBar, or a row's), taken back through its
 *      expression tree and its fold (ModelPointBackward, GradTargetEffectAdjoint) to cotangents of the observables of
 *      its units -- points or sub-points (UnitCotangent);
 *   3. adjoint: per energy task, one fill and factorisation of the level matrix, one T-matrix per T-matrix group, the
 *      T-matrix cotangents of all units of a group summed (AMatrixFunc::PointTBar), and one back-propagation through
 *      the level matrix per task (AMatrixFunc::BilinearAdjoint of the merged bilinear terms; a decay with external
 *      capture by AMatrixFunc::TBarAdjoint per group).
 *
 * The sums are formed in an order that does not depend on the number of threads: the tasks fall into a fixed number
 * of interleaved chunks, each summed in task order, the chunks then in chunk order; rows are summed task by task in
 * task order.  So a gradient (and a fit) gives the same bits at any thread count.
 */
namespace azure2 {

/// One cotangent of the value a unit's record holds: `value` times the derivative of that value (the observable, or
/// with xsComponent 1 / 2 its E1 / E2 part, -1 the angle-integrated A_0 of an angular distribution; with
/// ayCrossSection the differential cross section an A_y point carries), for data row `row` of the pass.
struct UnitCotangent {
  UnitRef unit;
  int row;
  double value;
  int xsComponent;
  bool ayCrossSection;
};

/// One data point of a gradient pass, `row` its index among the pass's points.
struct GradPoint {
  const Segment *segment;
  int segNum;    ///< 1-based
  int pointIdx;  ///< 0-based
  Unit *point;
};

/// The parts of a point's cotangents that belong to fitted parameters other than the level parameters (abundance,
/// Doppler temperature, sample thickness), in GradAccum's index order (0-based).
struct ParamCotangent {
  char kind;  ///< 'a' abundance, 'd' Doppler temperature, 't' thickness
  int index;
  double value;
};

/*!
 * The cotangents of the units of point `p` for the cotangent `fitBar` of its model (src/AZUREGrad.cpp): its expression
 * tree's (ModelPointBackward) and, for a leaf with a sub-point fold, the fold's (GradTargetEffectAdjoint), appended to
 * `out` with `row`; the parameter parts into `params`.  The records of `results` must hold the forward pass at the
 * parameters filled now.  False where the point is outside the analytic path (finite differences then).
 */
bool PointCotangents(WorkingCopy &data, const GradPoint &p, CNuc *compound, const Config &config,
                     azure2::Results &results, double fitBar, int row, GradAccum &params,
                     std::vector<UnitCotangent> &out);

/*!
 * d(objective)/d(level parameters) over the energy tasks, added to `accum` (and the fold and expression parameters
 * of the points, abundances, Doppler temperatures, thicknesses); fitBarFn gives each point's cotangent from its
 * model, called once per point in data order (its side effects, the norm gradient, see the same order at any thread
 * count).  The THM (HOES) segments are left out (AZURECalc::Gradient adds them).  False, with `accum` untouched and
 * `bailSeg` the segment, where a point is outside the analytic path.
 */
bool EnergyTaskGradient(CNuc *compound, WorkingCopy *data, azure2::Results &results, const Config &config,
                        const vector_matrix_r *shiftDeriv, const FitBarFn &fitBarFn, GradAccum &accum, int &bailSeg);

/*!
 * Rows d(model)/d(theta) times rowFitBar(point, model) over the energy tasks: for every point of the data segments
 * but the THM ones (`points`, in data order), its model (`models`) and its packed row (`rows`, row-major,
 * pmap.NumPacked() columns; the level parameters, abundances, Doppler temperatures, thicknesses).  The level-matrix
 * adjoint is taken by free parameter of `mask` where that is cheaper (AMatrixFunc::LevelBasis).  False, with
 * `bailSeg`, as EnergyTaskGradient; `numTasks` the energy tasks of the pass.
 */
bool EnergyTaskRows(CNuc *compound, WorkingCopy *data, azure2::Results &results, const Config &config,
                    const vector_matrix_r *shiftDeriv, const ParamIndexMap &pmap, const AMatrixFunc::FreeMask *mask,
                    const std::function<double(const GradPoint &, double)> &rowFitBar, std::vector<GradPoint> &points,
                    std::vector<double> &models, std::vector<double> &rows, int &bailSeg, int &numTasks);

}  // namespace azure2

#endif
