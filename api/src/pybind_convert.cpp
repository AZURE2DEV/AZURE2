// Format 3 for pyazr (docs/upgrade/project_format3.md sections 6.2, 6.3, 8): the validation of a format-3 text, as
// AZURE2 --validate, its description, as AZURE2 --describe, and the conversion of an old project to format 3, as AZURE2 --convert ... --to 3.  Both are
// experiment/'s reading path (ReadFormat3Text: records, validation, build) and convert/'s ConvertToFormat3; pyazr
// only hands texts over and parses the JSON.

#include <map>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "azure2/convert/Convert.h"
#include "azure2/experiment/Describe.h"
#include "azure2/experiment/ParameterFileNames.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/io/ProjectFormat.h"

namespace py = pybind11;

void BindConvert(py::module_ &m) {
  m.def(
      "not_format3_message",
      [](const std::string &text, const std::string &fileName) { return azure2::io::NotFormat3Message(text, fileName); },
      py::arg("text"), py::arg("file_name"),
      "The engine's refusal of a text that is not a format-3 project (io::NotFormat3Message: the format, the "
      "first old section, the command AZURE2 --convert OLD NEW), or '' for a format-3 text.");
  m.def(
      "convert_parameter_file",
      [](const std::string &parameterText, const std::string &oldProjectText, const std::string &projectDir) {
        const azure2::convert::ParameterFileConversion c =
            azure2::convert::ConvertParameterFile(parameterText, oldProjectText, projectDir);
        return py::make_tuple(c.refusal, py::bytes(c.text), c.renamed);
      },
      py::arg("parameter_text"), py::arg("old_project_text"), py::arg("project_dir"),
      "A parameter file of an old project renamed for the project's format-3 conversion (AZURE2 --convert-params; "
      "convert::ConvertParameterFile): (refusal or '', the renamed file as bytes, [(old name, new name)]).");
  m.def(
      "match_parameter_file_names",
      [](const std::string &file, const std::vector<std::string> &inFile, const std::vector<std::string> &ofRun,
         const std::vector<std::string> &optionalOfRun, const std::vector<std::string> &ofOtherRuns) {
        azure2::experiment::ParameterFileNames names;
        names.inFile = inFile;
        names.ofRun = ofRun;
        names.optionalOfRun = optionalOfRun;
        names.ofOtherRuns = ofOtherRuns;
        std::vector<std::string> out;
        for (const azure2::experiment::Diagnostic &d : azure2::experiment::MatchParameterFileNames(file, names))
          out.push_back(d.ToJson());
        return out;
      },
      py::arg("file"), py::arg("in_file"), py::arg("of_run"), py::arg("optional_of_run"), py::arg("of_other_runs"),
      "The names of a parameter file against a run's parameters, as AZURE2 matches them when it reads a param.sav "
      "(experiment::MatchParameterFileNames): at most one diagnostic, as JSON -- E_PARAMETER_FILE_NO_MATCH when no "
      "name of the file is a parameter of the run, W_PARAMETER_FILE_NAMES when some are not.");
  m.def(
      "validate_format3_text",
      [](const std::string &text, const std::string &fileName, const std::string &projectDir) -> py::tuple {
        const std::string refusal = azure2::io::NotFormat3Message(text, fileName);
        if (!refusal.empty()) return py::make_tuple(refusal, std::string("[]"), 0);  // section 8.4
        azure2::experiment::Format3ReadOptions read;
        read.fileName = fileName;
        read.projectDir = projectDir;  // "": data files as the engine opens them, relative to the working directory
        const azure2::experiment::Format3Project c = azure2::experiment::ReadFormat3Text(text, read);
        std::string json = "[";
        for (std::size_t i = 0; i < c.diagnostics.size(); i++) json += (i ? ",\n " : "\n ") + c.diagnostics[i].ToJson();
        json += c.diagnostics.empty() ? "]" : "\n]";
        return py::make_tuple(std::string(), json, c.buildChecksSkipped);
      },
      py::arg("text"), py::arg("file_name"), py::arg("project_dir") = "",
      "Validate a format-3 project text (AZURE2 --validate; data files relative to project_dir, '' = the working "
      "directory): (refusal of a format-1/2 text or '', the diagnostics as a "
      "JSON array, the number of the table's build checks not made: all of them when the validation found an error or "
      "<levels> could not be read, else 0).");
  m.def(
      "describe_format3_text",
      [](const std::string &text, const std::string &fileName, bool json, const std::string &projectDir,
         const std::string &record) -> py::tuple {
        const std::string refusal = azure2::io::NotFormat3Message(text, fileName);
        if (!refusal.empty()) return py::make_tuple(refusal, std::string("[]"), std::string());
        azure2::experiment::Format3ReadOptions read;
        read.fileName = fileName;
        read.projectDir = projectDir;
        const azure2::experiment::Format3Project c = azure2::experiment::ReadFormat3Text(text, read);
        std::string diagnostics = "[";
        for (std::size_t i = 0; i < c.diagnostics.size(); i++) diagnostics += (i ? ",\n " : "\n ") + c.diagnostics[i].ToJson();
        diagnostics += c.diagnostics.empty() ? "]" : "\n]";
        if (!c.Ok() || !c.levelsRead) return py::make_tuple(std::string(), diagnostics, std::string());
        return py::make_tuple(std::string(), diagnostics,
                              json ? azure2::experiment::DescribeJson(c, record) : azure2::experiment::DescribeText(c, record));
      },
      py::arg("text"), py::arg("file_name"), py::arg("json"), py::arg("project_dir") = "", py::arg("record") = "",
      "Describe a format-3 project text (AZURE2 --describe [--json]; data files relative to project_dir, '' = the "
      "working directory; only the record named `record`, '' = all): (refusal of a format-1/2 text or '', the "
      "diagnostics as a JSON array, the description: '' when the project has an error).");
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
