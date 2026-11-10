#ifndef AZURE2_MODELPASS_H
#define AZURE2_MODELPASS_H

#include "Results.h"

class CNuc;
class Config;
class EData;
class EPoint;
class ESegment;

/*!
 * The model pass (docs/upgrade/data_model.md section 2.4, step S4): after the point pass has filled the leaf values
 * of an evaluation (Results::slot), what every data or test point predicts goes into Results::model -- the value the
 * residual compares with the datum.  It is the segment's combination of its points' leaf values: the base point's
 * value (or its E1 / E2 part), times its abundance in a sample segment, plus the scaled components of a SUM, over the
 * denominator of a RATIO, or the transmission / self-shielded yield of the weighted terms (at every sub-point, before
 * the base point's kernel).  A leaf value is never overwritten by a model value.  S10 replaces the body of the
 * combination by the expression tree (experiment/, Expression.h); the callers stay.
 *
 * The residual of a datum is written once, here: r = (f - y n)/(dy n), f the model value, y and dy the c.m. datum
 * and its uncertainty, n the segment's normalization (on the data side, decision M5).
 */
namespace azure2 {

/// Whether the leaf values the model reads are in the buffer already (the parallel point pass,
/// EData::CalculatePointsInParallel, ran for these parameters), or each unmapped point is calculated where the model
/// needs it (inside an enclosing parallel region, without point mapping, or serially in data order).
enum class LeafValues { Ready, Calculate };

/// The model value of point `pointIndex` (0-based) of a data or test segment, stored in its model record (with the
/// angular-distribution coefficients of such a point, and the per-nucleus terms of a sample segment in the base
/// point's leaf record) and returned.  A point whose calculation throws, or a RATIO whose denominator is 0, is NaN.
double EvaluateModelPoint(ESegment &segment, int pointIndex, CNuc *compound, const Config &config, Results &results,
                          LeafValues leaves);
/// Every point of one segment, in order.
void EvaluateModelSegment(ESegment &segment, CNuc *compound, const Config &config, Results &results,
                          LeafValues leaves);
/// Every point of every segment of the data set, segment by segment, point by point.
void EvaluateModel(EData &data, CNuc *compound, const Config &config, Results &results, LeafValues leaves);

/// The residual of one datum: f - y n and dy n, as one expression for every path (the objective, the per-segment
/// chi^2, the least-squares rows, the API, the covariance groups, the writers' diagnostics).
struct Residual {
  double difference;  ///< f - y n
  double error;       ///< dy n
  /// (f - y n)/(dy n), 0 where dy n = 0.
  double Pull() const { return error != 0.0 ? difference / error : 0.0; }
  /// The same squared, as (f - y n)^2/(dy n)^2 (the per-segment sums), 0 where dy n = 0.
  double Chi2() const { return error != 0.0 ? (difference * difference) / (error * error) : 0.0; }
};
inline Residual DataResidual(double model, double y, double dy, double norm) { return {model - y * norm, dy * norm}; }
/// The residual of a data point: its model value in `results` against its c.m. datum, with the norm n.
Residual PointResidual(const EPoint &point, const Results &results, double norm);

}  // namespace azure2

#endif
