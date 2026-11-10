#ifndef AZURE2_EXPERIMENT_MEASUREMENT_H
#define AZURE2_EXPERIMENT_MEASUREMENT_H

#include <functional>
#include <istream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "azure2/experiment/SetupSpecs.h"

namespace azure2 {
namespace experiment {

struct Prediction;

/*!
 * What was measured (docs/upgrade/data_structure.md section 1.1): the rows of one data set that lie inside its
 * window, as read, and where they came from.  Built once at set-up, held as std::shared_ptr<const Measurement>,
 * shared by every working copy of the data, never written afterwards.  An extrapolation is a measurement too: its
 * rows are the grid (value 0, error 0.1, no file).
 */
struct MeasurementRows {
  bool angleIsLegendreOrder = false;  ///< the angle column holds the Legendre order k (angular coefficients)
  bool absolute = false;              ///< with angleIsLegendreOrder: the values are A_k (b/sr), not a_k = A_k / A_0
  std::vector<double> energy;         ///< MeV: E_lab of the beam, or the lab energy of the light decay particle
  std::vector<double> angle;          ///< deg: lab, or c.m. for the c.m. observables; or the Legendre order k
  std::vector<double> value, error;   ///< the quantity's engine unit (b, b/sr, 1, deg); error after a covariance file
  std::vector<double> binLow, binHigh;  ///< MeV lab, columns 5-6 of the file; NaN without them
  std::vector<int> record;            ///< 1-based data row in the file (comments, blanks not counted); 0 for a grid
  int Size() const { return (int)energy.size(); }
  /// Append one row.
  void Add(double e, double a, double y, double dy, double lo, double hi, int rec);
};

/// Where the rows came from: the file as opened, the data rows it holds, how many lie in the window, and the
/// lab-energy range of those (MeV; NaN when none is kept).
struct MeasurementProvenance {
  std::string file;
  int rowsRead = 0;
  int rowsKept = 0;
  double minEnergy = std::numeric_limits<double>::quiet_NaN();
  double maxEnergy = std::numeric_limits<double>::quiet_NaN();
};

struct Measurement {
  int key = 0;                    ///< the data set's number in the file (1-based), the label of its outputs
  bool extrapolation = false;     ///< the rows are a grid, not data
  std::shared_ptr<const MeasurementRows> rows;
  MeasurementProvenance provenance;
  std::shared_ptr<const Prediction> prediction;  ///< what the model predicts at the rows; null before it is built
  int NumRows() const { return rows ? rows->Size() : 0; }
};

/// The window of a data set (lab energy, and lab angle for the differential observables): the rows kept.
struct RowWindow {
  double minEnergy = 0.0, maxEnergy = 0.0;  ///< MeV lab
  bool angles = false;                      ///< the angle window applies
  double minAngle = 0.0, maxAngle = 0.0;    ///< deg
  bool Contains(double energy, double angle) const {
    if (!(energy >= minEnergy && energy <= maxEnergy)) return false;
    return !angles || (angle >= minAngle && angle <= maxAngle);
  }
};
/// The window of a data set's spec: its energy range, and its angle range for dsigma/dOmega and A_y.
RowWindow WindowOf(const DatasetSpec &spec);

/// One data row of a file as read (DataLine): its four columns, the optional columns 5-6, its 1-based record.
struct DataRow {
  double energy, angle, value, error;
  double binLow, binHigh;  ///< NaN without columns 5-6
  int record;
};

/*!
 * Reads the data rows of `in` (the file `file`, for messages) and keeps those inside `window`, in file order:
 * the rows of a measurement and their provenance.  With `legendreOrders` the angle column must hold an integer
 * order 0-20 (checked on every row, kept or not).  `keep`, when given, is called for each row inside the window
 * before it is kept: it may change the row's error (a covariance file's sqrt(C_ii)) and returns false to stop, having
 * reported why itself.  Returns false on a line that does not parse, a bad order (the message in `error`) or a
 * refusal of `keep` (`error` empty).
 */
bool ReadMeasurementRows(std::istream &in, const std::string &file, const RowWindow &window, bool legendreOrders,
                         const std::function<bool(DataRow &)> &keep, MeasurementRows &rows,
                         MeasurementProvenance &provenance, std::string &error);

}  // namespace experiment
}  // namespace azure2

#endif
