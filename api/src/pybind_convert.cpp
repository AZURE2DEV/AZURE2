// Format 3 for pyazr (docs/upgrade/project_format3.md sections 6.2, 8): the validation of a format-3 text, as
// AZURE2 --validate, and the conversion of an old project to format 3, as AZURE2 --convert ... --to 3.  Both are
// convert/'s (CheckFormat3Text, ConvertToFormat3); pyazr only hands texts over and parses the JSON.

#include <map>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "azure2/convert/Convert.h"

namespace py = pybind11;

void BindConvert(py::module_ &m) {
  m.def(
      "validate_format3_text",
      [](const std::string &text, const std::string &fileName) -> py::tuple {
        const std::string refusal = azure2::convert::NotFormat3Message(text, fileName);
        if (!refusal.empty()) return py::make_tuple(refusal, std::string("[]"), 0);  // section 8.4
        const azure2::convert::Format3Check c = azure2::convert::CheckFormat3Text(text, fileName);
        std::string json = "[";
        for (std::size_t i = 0; i < c.diagnostics.size(); i++) json += (i ? ",\n " : "\n ") + c.diagnostics[i].ToJson();
        json += c.diagnostics.empty() ? "]" : "\n]";
        return py::make_tuple(std::string(), json, c.buildChecksSkipped);
      },
      py::arg("text"), py::arg("file_name"),
      "Validate a format-3 project text (AZURE2 --validate): (refusal of a format-1/2 text or '', the diagnostics as a "
      "JSON array, the number of checks left to the format-3 reader because they need <levels> or the data files).");
  m.def(
      "convert_format3_text",
      [](const std::string &text, const std::string &projectDir, const std::map<std::string, std::string> &answers,
         bool exactGrids, bool strict, bool allowQuestions, const std::string &input, const std::string &output) {
        azure2::convert::ConvertOptions o;
        o.projectDir = projectDir;
        o.answers = answers;
        o.exactGrids = exactGrids;
        o.strict = strict;
        o.allowQuestions = allowQuestions;
        const azure2::convert::ConvertResult r = azure2::convert::ConvertToFormat3(text, o);
        py::dict d;
        d["ok"] = r.ok;
        d["refusal"] = r.refusal;
        d["text"] = py::bytes(r.text);
        d["report"] = r.ReportJson(input, output);
        py::list side;
        for (const auto &f : r.sideFiles) side.append(py::make_tuple(f.first, py::bytes(f.second)));
        d["side_files"] = side;
        return d;
      },
      py::arg("text"), py::arg("project_dir"), py::arg("answers"), py::arg("exact_grids"), py::arg("strict"),
      py::arg("allow_questions"), py::arg("input"), py::arg("output"),
      "Convert an old project text to format 3 (AZURE2 --convert ... --to 3): a dict with ok, refusal, text (bytes), "
      "report (the JSON report) and side_files [(path, bytes)].");
}
