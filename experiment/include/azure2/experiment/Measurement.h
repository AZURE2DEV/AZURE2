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
 * the rows of a measurement and their provenance.  Each row is first taken to the engine's units and frame by
 * `columns` (DataRows.h RowToEngine; the identity leaves every bit), so the window, the order check, `keep` and the
 * rows kept are all in engine units.  With `legendreOrders` the angle column must hold an integer
 * order 0-20 (checked on every row, kept or not).  `keep`, when given, is called for each row inside the window
 * before it is kept: it may change the row's error (a covariance file's sqrt(C_ii)) and returns false to stop, having
 * reported why itself.  Returns false on a line that does not parse, a bad order (the message in `error`) or a
 * refusal of `keep` (`error` empty).
 */
bool ReadMeasurementRows(std::istream &in, const std::string &file, const RowWindow &window,
                         const ColumnConversion &columns, bool legendreOrders,
                         const std::function<bool(DataRow &)> &keep, MeasurementRows &rows,
                         MeasurementProvenance &provenance, std::string &error);

/// Is `k` (the angle column of an angular-coefficient row) a Legendre order: an integer 0-20?
bool IsLegendreOrder(double k);

/// One problem of a row of a data file (CheckDataRows).
struct DataRowProblem {
  enum class Kind {
    NotParsed,       ///< not a comment and not four numbers (the engine's reader stops there)
    TooFewColumns,   ///< fewer numeric columns than the data set's column roles name
    NotAnOrder,      ///< an angular-coefficient row whose second column is no Legendre order 0-20
    RefusedKept      ///< a row inside the window the caller's check refuses (a non-positive error, a variance)
  };
  Kind kind = Kind::NotParsed;
  int line = 0;    ///< 1-based physical line of the file
  int record = 0;  ///< 1-based data row (comments and blanks not counted); 0 when the line did not parse
  int column = 0;  ///< 1-based column the problem is in; 0 for the row as a whole
  std::string found;    ///< the text found (the token, or the line)
  std::string message;  ///< what is wrong, naming the file, line and column
};

/// Every problem of a data file and the rows it would give: what the set-up's reader (ReadMeasurementRows) refuses,
/// found on every line instead of the first.
struct DataRowCheck {
  MeasurementProvenance provenance;  ///< rows read and kept as ReadMeasurementRows counts them
  std::vector<DataRowProblem> problems;
  int problemsNotListed = 0;  ///< problems beyond `maxProblems` (counted, not listed)
};

/*!
 * Reads every line of `in` (the file `file`, for messages) as ReadMeasurementRows does (DataLine, the row taken to
 * engine units by `columns`, the window) and lists, instead of stopping at the first: each line that does not parse,
 * each row with fewer numeric columns than `roleColumns`, with `legendreOrders` each row whose second column is no
 * Legendre order, and each row inside the window that `checkKept` refuses (it returns why, "" to accept; the row's
 * value and error in engine units).  At most `maxProblems` are listed.
 */
DataRowCheck CheckDataRows(std::istream &in, const std::string &file, const RowWindow &window,
                           const ColumnConversion &columns, bool legendreOrders, std::size_t roleColumns,
                           const std::function<std::string(const DataRow &)> &checkKept, int maxProblems = 100);

}  // namespace experiment
}  // namespace azure2

#endif
