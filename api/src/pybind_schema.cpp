// The table of project format 3 for pyazr and the web UI (docs/upgrade/project_format3.md sections 6.1, 7.3): the
// same JSON document as AZURE2 --schema (experiment/'s Schema::ToJson), returned as a dict.  pyazr spells no field,
// quantity or unit itself; it reads them from here.  And what a project's data sets predict, as AZURE2
// --print-predictions prints it (src/Predictions.cpp).  And the extrapolation grids of section 3.8 (pyazr/grids.py): a
// grid's text read by the table's parser and built by experiment's ReadGridValues, the adaptive refinement RefineGrid
// on a Python model.  And the data files a data covariance file names (its one reader, ReadCovarianceFile).

#include <pybind11/pybind11.h>

#include <pybind11/stl.h>

#include <sstream>

#include "Predictions.h"
#include "azure2/experiment/AdaptiveGrid.h"
#include "azure2/experiment/CovarianceFile.h"
#include "azure2/experiment/Grids.h"
#include "azure2/experiment/Parsers.h"
#include "azure2/experiment/Schema.h"
#include "azure2/io/RecordReader.h"

namespace py = pybind11;

void BindSchema(py::module_ &m) {
  m.def(
      "schema", []() { return py::module_::import("json").attr("loads")(azure2::experiment::Schema::ToJson()); },
      "The table of project format 3 as a dict: record kinds and fields, quantities, units, diagnostic codes "
      "(the JSON printed by AZURE2 --schema).");
  m.def("schema_json", []() { return azure2::experiment::Schema::ToJson(); },
        "The same table as the JSON text, byte for byte as AZURE2 --schema prints it (without the final newline).");
  m.def(
      "predictions",
      [](const std::string &path) {
        std::ostringstream messages;
        std::vector<std::string> lines = azure2::ProjectPredictions(path, messages);
        return py::make_tuple(lines, messages.str());
      },
      py::arg("path"),
      "The expression tree of every data set and extrapolation of the project at path, as the lines AZURE2 "
      "--print-predictions prints, and the set-up's messages: (lines, messages).");
  m.def(
      "extrapolation_grid",
      [](const std::string &text, double m1, double m2, bool relativistic, const std::string &projectDir,
         bool angles) -> py::tuple {
        namespace ex = azure2::experiment;
        std::string comment;
        std::vector<azure2::io::RecordDiagnostic> problems;
        const std::vector<azure2::io::Token> tokens = azure2::io::Tokenize(text, 1, comment, problems);
        if (!problems.empty()) return py::make_tuple(false, std::vector<double>(), problems[0].message);
        if (tokens.empty()) return py::make_tuple(false, std::vector<double>(), std::string("the grid is empty"));
        std::vector<ex::Diagnostic> diagnostics;
        ex::ParseContext context;
        context.out = &diagnostics;
        context.base.recordKind = "extrapolation";
        context.base.field = angles ? "angles" : "energies";
        ex::TokenCursor cursor;
        cursor.tokens = &tokens;
        cursor.line = 1;
        cursor.endCol = tokens.back().col + (int)tokens.back().text.size() + (tokens.back().quoted ? 2 : 0);
        ex::Grid grid;
        if (!ex::ParseGrid(ex::ValueSpec(ex::ValueType::Grid, angles ? "angle" : "energy"), cursor, context, grid))
          return py::make_tuple(false, std::vector<double>(),
                                diagnostics.empty() ? std::string("the grid does not parse") : diagnostics[0].message);
        if (!cursor.AtEnd())
          return py::make_tuple(false, std::vector<double>(), "'" + cursor.Peek().text + "' after the grid");
        ex::PairInfo pair;
        pair.lightMass = m1;
        pair.heavyMass = m2;
        ex::GridSetting setting;
        setting.energies = !angles;
        setting.entrance = m1 > 0.0 && m2 > 0.0 ? &pair : nullptr;
        setting.relativistic = relativistic;
        setting.projectDir = projectDir;
        ex::GridValues values = ex::ReadGridValues(grid, setting);
        return py::make_tuple(values.ok, values.values, values.error);
      },
      py::arg("text"), py::arg("m1"), py::arg("m2"), py::arg("relativistic"), py::arg("project_dir"), py::arg("angles"),
      "The values of an extrapolation grid's text (MeV lab, or deg): (ok, values, message); pyazr.extrapolation_grid.");
  m.def(
      "refine_grid",
      [](const std::vector<double> &start, py::function model, double tolerance, long maxPoints) {
        namespace ex = azure2::experiment;
        const ex::AdaptiveGridResult r = ex::RefineGrid(start, tolerance, maxPoints, [&](const std::vector<double> &x) {
          return model(x).cast<std::vector<std::vector<double>>>();
        });
        py::dict out;
        out["x"] = r.x;
        out["values"] = r.values;
        out["budget_reached"] = r.budgetReached;
        out["largest_deviation"] = r.largestDeviation;
        out["largest_deviation_at"] = r.largestDeviationAt;
        out["sweeps"] = r.sweeps;
        out["evaluations"] = r.evaluations;
        return out;
      },
      py::arg("start"), py::arg("model"), py::arg("tolerance"), py::arg("max_points"),
      "The adaptive refinement of an extrapolation grid (experiment::RefineGrid) on a Python model; pyazr.refine_grid.");
  m.def(
      "covariance_file_names",
      [](const std::string &path, const std::string &directory) {
        namespace ex = azure2::experiment;
        ex::CovarianceFileAccess access;
        access.resolve = [&directory](const std::string &name) {
          return name.empty() || name[0] == '/' || directory.empty() ? name : directory + "/" + name;
        };
        std::string why;
        const std::shared_ptr<ex::CovarianceFile> f = ex::ReadCovarianceFile(path, access, why);
        if (!f) throw py::value_error(why);
        return f->files;
      },
      py::arg("path"), py::arg("directory"),
      "The data files the data covariance file `path` names in its header, as written (experiment::ReadCovarianceFile, "
      "the engine's one reader), with every relative name -- the file's and the data files' -- read from `directory` "
      "(the project's working directory, as the engine reads them); ValueError with the reader's reason.");
}
