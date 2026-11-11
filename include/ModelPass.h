#ifndef AZURE2_MODELPASS_H
#define AZURE2_MODELPASS_H

#include <vector>

#include "Results.h"
#include "azure2/experiment/Prediction.h"

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
 * the base point's kernel).  A leaf value is never overwritten by a model value.  The combination is the segment's
 * expression tree (step S10: Predictions.h, experiment/Expression.h), built once at set-up.
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
/// `data` is the data set that holds the segment (its plan, its target effects); below too.
double EvaluateModelPoint(EData &data, ESegment &segment, int pointIndex, CNuc *compound, const Config &config,
                          Results &results, LeafValues leaves);
/// Every point of one segment, in order.
void EvaluateModelSegment(EData &data, ESegment &segment, CNuc *compound, const Config &config, Results &results,
                          LeafValues leaves);
/// Every point of every segment of the data set, segment by segment, point by point.
void EvaluateModel(EData &data, CNuc *compound, const Config &config, Results &results, LeafValues leaves);

/// The cotangents of one point's model value m (step S11): g dm/d(leaf) for every leaf of the segment's tree that the
/// R-matrix computes (a Coulomb reference leaf, which no level parameter moves, is left out), at the
/// term's point or -- a transmission or yield under a sub-point effect -- at each of its sub-points, and g dm/d(param)
/// for the tree's parameters (azure2::experiment::ParamRef: an abundance per weighted term, the thickness) with one
/// slot more, the last, for a fitted Doppler temperature.  The adjoint of the composition, from the tree's backward;
/// the leaves' own adjoint (the R-matrix part) is the caller's.
struct ModelCotangents {
  struct Leaf {
    EPoint *point = nullptr;  ///< the term's point
    int xsComponent = 0;      ///< which value it reads: 0 the observable, 1 / 2 the E1 / E2 part
    double cotangent = 0.0;   ///< at the point (no sub-point stage)
    std::vector<double> sub;  ///< at its sub-points 1..n (subPoints)
  };
  std::vector<Leaf> leaves;
  std::vector<double> params;
  const std::vector<experiment::ParamRef> *refs = nullptr;  ///< what params[0 .. size - 2] are
  bool subPoints = false;     ///< the cotangents are at the sub-points
  bool dopplerOrder = false;  ///< sub-points under a Doppler map (the adjoint takes them term by term)
};
/// The cotangents of point `pointIndex` of a data or test segment for the cotangent g of its model value, from the
/// leaf values in `results` (the point pass, or the caller, calculated them).  False where the composition has no
/// exact adjoint (a missing term point, a kernel that is not linear): the caller falls back to finite differences.
bool ModelPointBackward(const EData &data, ESegment &segment, int pointIndex, CNuc *compound, const Config &config,
                        const Results &results, double g, ModelCotangents &out);

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
/// The residual of row i (0-based) of a segment: its model value in `results` against its c.m. datum, with the norm n.
Residual PointResidual(const ESegment &segment, int i, const Results &results, double norm);

}  // namespace azure2

#endif
