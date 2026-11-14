// The rows of a data set's file as the engine reads them (experiment/Measurement.h ReadMeasurementRows): every data
// row in the engine's units and frame, its 1-based record number and whether the data set's window keeps it.  The
// covariance builder (pyazr/covariance.py) counts the rows of a covariance file's `file` lines and takes the values
// and errors its parts scale from these; it reads no data file itself.

#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "azure2/experiment/Measurement.h"
#include "azure2/experiment/ProjectPaths.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/experiment/Schema.h"

namespace py = pybind11;

namespace {

namespace ex = azure2::experiment;

py::dict DataRows(const std::string &text, const std::string &projectDir, const std::string &name) {
  ex::Format3ReadOptions read;
  read.projectDir = projectDir;
  read.setUpChecks = false;
  read.allRecordsActive = true;  // an inactive data set's rows are read as when it is switched on
  const ex::Format3Project p = ex::ReadFormat3Text(text, read);
  std::vector<const ex::ValidatedRecord *> sets;
  for (const ex::ValidatedRecord &r : p.validated.records)
    if (r.kind == "dataset") sets.push_back(&r);
  std::size_t index = sets.size();
  for (std::size_t i = 0; i < sets.size(); i++)
    if (sets[i]->name == name) index = i;
  if (index == sets.size()) throw std::invalid_argument("data_rows: no data set named '" + name + "'");
  if (index >= p.built.data.datasets.items.size()) {
    std::string first;
    for (const ex::Diagnostic &d : p.diagnostics)
      if (d.severity == ex::Severity::Error) {
        first = d.message;
        break;
      }
    throw std::invalid_argument("data_rows: the project does not build, so data set '" + name +
                                "' has no rows" + (first.empty() ? std::string() : ": " + first));
  }
  const ex::DatasetSpec &d = p.built.data.datasets.items[index].spec;
  const std::string path = ex::PathToOpen(projectDir, d.dataFile);
  std::ifstream in(path);
  if (!in) throw std::invalid_argument("data_rows: the data file " + path + " of '" + name + "' cannot be opened");
  ex::RowWindow all;  // every row: the window is applied below, row by row, as the engine's reader applies it
  all.minEnergy = -std::numeric_limits<double>::infinity();
  all.maxEnergy = std::numeric_limits<double>::infinity();
  ex::MeasurementRows rows;
  ex::MeasurementProvenance provenance;
  std::string error;
  if (!ex::ReadMeasurementRows(in, d.dataFile, all, d.columns,
                               d.reaction.observable == ex::Observable::AngularCoefficients, nullptr, rows,
                               provenance, error))
    throw std::invalid_argument("data_rows: " + error);
  const ex::RowWindow window = ex::WindowOf(d);
  std::vector<bool> kept(rows.Size());
  for (int i = 0; i < rows.Size(); i++) kept[i] = window.Contains(rows.energy[i], rows.angle[i]);
  py::dict out;
  out["data_file"] = d.dataFile;
  out["path"] = path;
  out["engine_unit"] = sets[index]->quantity ? sets[index]->quantity->engineUnit : std::string();
  out["value_factor"] = d.columns.value;
  out["energy"] = rows.energy;
  out["angle"] = rows.angle;
  out["value"] = rows.value;
  out["error"] = rows.error;
  out["record"] = rows.record;
  out["in_window"] = kept;
  out["covariance"] = d.covarianceFile;
  out["norm_fit"] = d.varyNorm;
  out["norm_prior_percent"] = d.normError;
  out["norm_group"] = d.normGroup;
  return out;
}

}  // namespace

void BindDataRows(py::module_ &m) {
  m.def("data_rows", &DataRows, py::arg("text"), py::arg("project_dir"), py::arg("name"),
        "Every data row of the file of data set `name` of the format-3 project `text` (data files relative to "
        "`project_dir`), read as the engine reads it (Measurement.h ReadMeasurementRows): a dict with data_file (as "
        "written), path (as opened), engine_unit (the quantity's), value_factor (engine units per unit of the value "
        "column), energy (MeV, lab), angle (deg), value and error (engine units), record (1-based data row: blank "
        "and comment lines not counted), in_window (the data set's window keeps the row) and covariance (the "
        "record's covariance file, '' without), norm_fit (the normalization is a fit parameter), norm_prior_percent "
        "(its prior's width in percent, 0 = none) and norm_group (a shared normalization's label, '' = none).  "
        "ValueError when the data set does not exist, the project does not build or the file cannot be read.");
}
