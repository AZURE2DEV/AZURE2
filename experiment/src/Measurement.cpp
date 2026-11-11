#include "azure2/experiment/Measurement.h"

#include <cmath>
#include <sstream>

#include "azure2/experiment/DataLine.h"

namespace azure2 {
namespace experiment {

void MeasurementRows::Add(double e, double a, double y, double dy, double lo, double hi, int rec) {
  energy.push_back(e);
  angle.push_back(a);
  value.push_back(y);
  error.push_back(dy);
  binLow.push_back(lo);
  binHigh.push_back(hi);
  record.push_back(rec);
}

RowWindow WindowOf(const DatasetSpec &spec) {
  RowWindow w;
  w.minEnergy = spec.minE;
  w.maxEnergy = spec.maxE;
  const Observable o = spec.reaction.observable;
  w.angles = o == Observable::DifferentialLab || o == Observable::DifferentialCM || o == Observable::AnalyzingPower;
  w.minAngle = spec.minA;
  w.maxAngle = spec.maxA;
  return w;
}

bool ReadMeasurementRows(std::istream &in, const std::string &file, const RowWindow &window,
                         const ColumnConversion &columns, bool legendreOrders,
                         const std::function<bool(DataRow &)> &keep, MeasurementRows &rows,
                         MeasurementProvenance &provenance, std::string &error) {
  constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
  provenance.file = file;
  int lineNumber = 0;
  int record = 0;
  while (true) {
    DataLine line(in);
    lineNumber += line.linesConsumed();
    if (line.atEnd()) break;
    if (!line.valid()) {
      // Not a comment and not a data row: say exactly where; silently skipping a mangled row would hide a corrupt file.
      std::ostringstream m;
      m << "ERROR: Cannot parse line " << lineNumber << " of " << file << ": \"" << line.raw() << "\"\n"
        << "       A data line needs at least four numeric columns (energy, angle, value,"
        << " uncertainty). Lines beginning with '#' and blank lines are ignored.";
      error = m.str();
      return false;
    }
    record++;
    provenance.rowsRead = record;
    line.ToEngineUnits(columns);
    DataRow row{line.energy(), line.angle(), line.crossSection(), line.error(), kNaN, kNaN, record};
    if (line.numExtra() >= 2) {
      row.binLow = line.extra(0);
      row.binHigh = line.extra(1);
    }
    if (legendreOrders) {
      const double k = line.angle();
      if (k < 0.0 || k > 20.0 || std::fabs(k - std::round(k)) > 1.0e-9) {
        std::ostringstream m;
        m << "ERROR: line " << lineNumber << " of " << file << ": an angular-coefficient"
          << " segment's second column is the Legendre order k (an integer, 0-20), got " << k << ".";
        error = m.str();
        return false;
      }
    }
    if (!window.Contains(row.energy, row.angle)) continue;
    if (keep && !keep(row)) {
      error.clear();
      return false;
    }
    rows.Add(row.energy, row.angle, row.value, row.error, row.binLow, row.binHigh, row.record);
    provenance.rowsKept++;
    if (!(row.energy >= provenance.minEnergy)) provenance.minEnergy = row.energy;
    if (!(row.energy <= provenance.maxEnergy)) provenance.maxEnergy = row.energy;
  }
  return true;
}

}  // namespace experiment
}  // namespace azure2
