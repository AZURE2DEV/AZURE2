#include "azure2/experiment/Measurement.h"

#include <cmath>
#include <limits>
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

bool IsLegendreOrder(double k) { return !(k < 0.0 || k > 20.0 || std::fabs(k - std::round(k)) > 1.0e-9); }

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
      if (!IsLegendreOrder(k)) {
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

namespace {

// The 1-based column of the first of the four required tokens of `raw` that is no number, or 0 when there are fewer
// than four tokens (the diagnostic's column; DataLine decided that the line does not parse).
int FirstBadColumn(const std::string &raw, std::string &found, int &tokens) {
  std::istringstream in(raw);
  std::string token;
  tokens = 0;
  while (in >> token) {
    tokens++;
    if (tokens > 4) break;
    std::istringstream number(token);
    double x;
    char rest;
    if (!(number >> x) || (number >> rest)) {
      found = token;
      return tokens;
    }
  }
  return 0;
}

}  // namespace

DataRowCheck CheckDataRows(std::istream &in, const std::string &file, const RowWindow &window,
                           const ColumnConversion &columns, bool legendreOrders, std::size_t roleColumns,
                           const std::function<std::string(const DataRow &)> &checkKept, int maxProblems) {
  constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
  DataRowCheck out;
  out.provenance.file = file;
  auto add = [&](DataRowProblem p) {
    if ((int)out.problems.size() < maxProblems) out.problems.push_back(std::move(p));
    else out.problemsNotListed++;
  };
  auto where = [&](int lineNumber, int column) {
    return file + ", line " + std::to_string(lineNumber) + (column > 0 ? ", column " + std::to_string(column) : "");
  };
  int lineNumber = 0, record = 0;
  while (true) {
    DataLine line(in);
    lineNumber += line.linesConsumed();
    if (line.atEnd()) break;
    if (!line.valid()) {
      DataRowProblem p;
      p.line = lineNumber;
      int tokens = 0;
      p.column = FirstBadColumn(line.raw(), p.found, tokens);
      if (p.column == 0) {
        p.found = line.raw();
        p.message = where(lineNumber, 0) + ": " + std::to_string(tokens) +
                    " column(s); a data row needs four numbers (energy, angle, value, uncertainty)";
      } else {
        p.message = where(lineNumber, p.column) + ": '" + p.found + "' is not a number";
      }
      add(p);
      continue;
    }
    record++;
    out.provenance.rowsRead = record;
    const std::size_t numeric = 4 + (std::size_t)line.numExtra();
    if (numeric < roleColumns) {
      DataRowProblem p;
      p.kind = DataRowProblem::Kind::TooFewColumns;
      p.line = lineNumber;
      p.record = record;
      p.column = (int)numeric + 1;
      p.found = line.raw().empty() ? std::to_string(numeric) + " columns" : line.raw();
      p.message = where(lineNumber, 0) + ": " + std::to_string(numeric) + " numeric column(s), the column roles name " +
                  std::to_string(roleColumns);
      add(p);
    }
    line.ToEngineUnits(columns);
    DataRow row{line.energy(), line.angle(), line.crossSection(), line.error(), kNaN, kNaN, record};
    if (line.numExtra() >= 2) {
      row.binLow = line.extra(0);
      row.binHigh = line.extra(1);
    }
    if (legendreOrders && !IsLegendreOrder(row.angle)) {
      DataRowProblem p;
      p.kind = DataRowProblem::Kind::NotAnOrder;
      p.line = lineNumber;
      p.record = record;
      p.column = 2;
      std::ostringstream k;
      k << row.angle;
      p.found = k.str();
      p.message = where(lineNumber, 2) + ": an angular-coefficient row's second column is the Legendre order k (an "
                                         "integer, 0-20), found " + p.found;
      add(p);
      continue;
    }
    if (!window.Contains(row.energy, row.angle)) continue;
    out.provenance.rowsKept++;
    if (!(row.energy >= out.provenance.minEnergy)) out.provenance.minEnergy = row.energy;
    if (!(row.energy <= out.provenance.maxEnergy)) out.provenance.maxEnergy = row.energy;
    const std::string why = checkKept ? checkKept(row) : std::string();
    if (!why.empty()) {
      DataRowProblem p;
      p.kind = DataRowProblem::Kind::RefusedKept;
      p.line = lineNumber;
      p.record = record;
      p.column = 4;
      std::ostringstream e;
      e << row.error;
      p.found = e.str();
      p.message = where(lineNumber, 4) + " (data row " + std::to_string(record) + "): " + why;
      add(p);
    }
  }
  return out;
}

}  // namespace experiment
}  // namespace azure2
