// The project file for pyazr (docs/upgrade/project.md section 2, project_format3.md section 4): the lossless document
// io::Document that pyazr's AzrModel and Project read and edit, as the engine reads the file; and the project's
// format-3 records as the engine reads (io::ReadRecords), writes (io::WriteRecords, experiment::NuisanceText) and
// resolves them (the validator's references, the sources of the fit parameters).  pyazr's Project (pyazr/project.py)
// spells no field, unit or grammar: it hands texts over and applies the line edits these functions describe.
// Texts are bytes (UTF-8 as in the file); line numbers and columns are the file's, 1-based, columns in bytes.

#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "azure2/experiment/ProjectReader.h"
#include "azure2/experiment/Schema.h"
#include "azure2/experiment/ValueText.h"
#include "azure2/io/Convert.h"
#include "azure2/io/Document.h"
#include "azure2/io/RecordReader.h"

namespace py = pybind11;

namespace {

namespace io = azure2::io;
namespace ex = azure2::experiment;

// The text of a field's value as it stands in its line: from the first token to the end of the last (with its quotes).
std::string ValueText(const io::Document &doc, const io::Field &f) {
  if (f.tokens.empty() || f.line < 1 || (std::size_t)f.line > doc.NumLines()) return std::string();
  const std::string &line = doc.Line(f.line - 1);
  const io::Token &last = f.tokens.back();
  const std::size_t from = f.tokens.front().col - 1, to = last.col - 1 + last.text.size() + (last.quoted ? 2 : 0);
  return to <= line.size() && from <= to ? line.substr(from, to - from) : std::string();
}

py::list Comments(const std::vector<std::string> &comments) {
  py::list out;
  for (const std::string &c : comments) out.append(py::bytes(c));
  return out;
}

py::dict RecordDict(const io::Document &doc, const io::Record &r) {
  py::dict d;
  d["kind"] = r.kind;
  d["name"] = py::bytes(r.name);
  d["line"] = r.line;
  d["name_col"] = r.nameCol;
  d["comments"] = Comments(r.comments);
  int last = r.line;
  py::list fields;
  for (const io::Field &f : r.fields) {
    py::dict x;
    x["key"] = py::bytes(f.key);
    x["line"] = f.line;
    x["col"] = f.col;
    x["value"] = py::bytes(ValueText(doc, f));
    x["value_col"] = f.tokens.empty() ? 0 : f.tokens.front().col;
    x["comment"] = py::bytes(f.trailingComment);
    x["comments"] = Comments(f.comments);
    fields.append(x);
    last = std::max(last, f.line);
  }
  d["fields"] = fields;
  d["last_line"] = last;
  return d;
}

// One value text as the file's tokens and the comment after them ('#' starts a comment, as in the file); the reader's
// problems (E_UNTERMINATED_QUOTE, E_QUOTE_IN_TOKEN) appended as (code, message).
using Problems = std::vector<std::pair<std::string, std::string>>;
std::vector<io::Token> ValueTokens(const std::string &key, const std::string &text, std::string &comment, Problems &problems) {
  std::vector<io::RecordDiagnostic> found;
  std::vector<io::Token> tokens = io::Tokenize(text, 1, comment, found);
  for (const io::RecordDiagnostic &d : found) problems.push_back({d.code, key + ": " + d.message});
  return tokens;
}

}  // namespace

void BindNuclearPart(py::module_ &m);  // pybind_nuclear_part.cpp: <levels>, <config>, <thm> read and edited

void BindProject(py::module_ &m) {
  BindNuclearPart(m);
  // The project file as a lossless document (io::Document; docs/upgrade/project.md section 2): pyazr's AzrModel
  // reads and edits the file through it, as the engine does.  Lines are str (decoded UTF-8, invalid bytes kept as
  // surrogates on the Python side by the caller); sections are (name, open, close, depth, parent) with close = -1 for
  // an unclosed section and parent = -1 at top level.
  {
    using azure2::io::Document;
    auto sectionTuple = [](const Document &d, std::size_t k) {
      const Document::Section &s = d.Sections()[k];
      const long npos = -1;
      return py::make_tuple(s.name, (long)s.open, s.close == Document::Section::npos ? npos : (long)s.close, s.depth,
                            s.parent == Document::Section::npos ? npos : (long)s.parent);
    };
    py::class_<Document>(m, "Document")
        .def_static("parse", [](const py::bytes &b) { return Document::Parse(std::string(b)); }, py::arg("data"))
        .def("write", [](const Document &d) { return py::bytes(d.Write()); })
        .def("num_lines", &Document::NumLines)
        .def("line", [](const Document &d, std::size_t i) {
          if (i >= d.NumLines()) throw py::index_error("line index out of range");
          return py::bytes(d.Line(i));
        }, py::arg("i"))
        .def("ending", [](const Document &d, std::size_t i) {
          if (i >= d.NumLines()) throw py::index_error("line index out of range");
          return py::bytes(d.Ending(i));
        }, py::arg("i"))
        .def("sections", [sectionTuple](const Document &d) {
          py::list out;
          for (std::size_t k = 0; k < d.Sections().size(); k++) out.append(sectionTuple(d, k));
          return out;
        })
        .def("problems", &Document::Problems)
        .def("replace_lines", [](Document &d, std::size_t first, std::size_t last, const std::vector<py::bytes> &lines) {
          if (first > last || last > d.NumLines()) throw py::index_error("line range out of range");
          std::vector<std::string> v;
          for (const py::bytes &b : lines) v.push_back(std::string(b));
          d.ReplaceLines(first, last, v);
        }, py::arg("first"), py::arg("last"), py::arg("lines"))
        .def("common_ending", [](const Document &d) { return py::bytes(d.CommonEnding()); })
        .def("format", [](const Document &d) { return azure2::io::FormatOf(d); });
  }
  m.def(
      "record_sections",
      [](const io::Document &doc) {
        py::list out;
        for (std::size_t k = 0; k < doc.Sections().size(); k++) {
          const io::Document::Section &s = doc.Sections()[k];
          const std::vector<std::string> kinds = ex::Schema::KindsOfSection(s.name);
          if (s.depth != 0 || kinds.empty()) continue;
          py::dict d;
          d["section"] = s.name;
          d["open"] = (long)s.open + 1;
          d["close"] = s.close == io::Document::Section::npos ? -1L : (long)s.close + 1;
          py::list records, problems;
          if (s.close != io::Document::Section::npos) {
            const io::RecordReadResult r = io::ReadRecords(doc, s, kinds);
            for (const io::Record &rec : r.set.records) records.append(RecordDict(doc, rec));
            for (const io::RecordDiagnostic &p : r.diagnostics) problems.append(py::make_tuple(p.code, p.line, p.message));
          }
          d["records"] = records;
          d["problems"] = problems;
          out.append(d);
        }
        return out;
      },
      py::arg("doc"),
      "The format-3 record sections of a document (top level, in file order) as io::ReadRecords reads them: dicts "
      "with section, open and close (1-based lines of the tags, close -1 when not closed), records and the reader's "
      "problems (code, line, message).  A record: kind, name, line, name_col, last_line, comments (above the header), "
      "fields; a field: key, line, col, value (its text from the first token to the end of the last), value_col, "
      "comment (after the value), comments (above it).  Texts are bytes.");
  m.def(
      "record_lines",
      [](const std::string &kind, const std::string &name, const std::vector<std::pair<std::string, std::string>> &fields) {
        io::RecordSet set;
        io::Record r;
        r.kind = kind;
        r.name = name;
        Problems problems;
        for (const auto &f : fields) {
          io::Field x;
          x.key = f.first;
          x.tokens = ValueTokens(f.first, f.second, x.trailingComment, problems);
          r.fields.push_back(x);
        }
        set.records.push_back(r);
        py::list lines;
        for (const std::string &l : io::WriteRecords(set)) lines.append(py::bytes(l));
        return py::make_tuple(lines, problems);
      },
      py::arg("kind"), py::arg("name"), py::arg("fields"),
      "A record as io::WriteRecords writes it, from (key, value text) pairs, a '#' in a value starting the field's "
      "comment: (lines as bytes, the reader's problems of the values' syntax as (code, message)).");
  m.def(
      "value_text",
      [](const std::string &key, const std::string &text) {
        Problems problems;
        std::string comment;
        const std::vector<io::Token> tokens = ValueTokens(key, text, comment, problems);
        return py::make_tuple(py::bytes(io::WriteTokens(tokens)), py::bytes(comment), problems);
      },
      py::arg("key"), py::arg("text"),
      "A field's value text as the record writer writes it (tokens separated by one blank, quoted where needed) and "
      "the comment a '#' starts: (text, comment, the reader's problems as (code, message)); texts as bytes.");
  m.def(
      "nuisance_text",
      [](const std::string &value, const std::string &unit, const std::string &mode, const std::string &priorCentre,
         const std::string &priorWidth, const std::string &gridMax) {
        ex::NuisanceWords w;
        w.value = value;
        w.unit = unit;
        w.mode = mode;
        w.priorCentre = priorCentre;
        w.priorWidth = priorWidth;
        w.gridMax = gridMax;
        return ex::NuisanceText(w);
      },
      py::arg("value"), py::arg("unit"), py::arg("mode"), py::arg("prior_centre"), py::arg("prior_width"),
      py::arg("grid_max"),
      "The text of a nuisance (experiment::NuisanceText, the converter's writer): '<value> [<unit>] <mode> [prior "
      "<centre> <width>] [grid_max <value>]', or the mode alone without a value; '' leaves a part out.");
  m.def(
      "record_references",
      [](const std::string &text) {
        ex::Format3ReadOptions read;
        read.build = false;
        const ex::Format3Project c = ex::ReadFormat3Text(text, read);
        py::list out;
        for (const ex::RecordReference &r : c.validated.references) {
          py::dict d;
          d["record"] = py::bytes(r.record);
          d["field"] = r.field;
          d["line"] = r.line;
          d["col"] = r.col;
          d["name"] = py::bytes(r.name);
          d["epoch"] = py::bytes(r.epoch);
          out.append(d);
        }
        return out;
      },
      py::arg("text"),
      "Every reference of a record to another in a format-3 project text, as the validator resolves names "
      "(experiment::RecordReference): dicts with record, field, line, col (of the name), name and epoch (bytes).");
  m.def(
      "parameter_names",
      [](const std::string &text, bool allActive) {
        ex::Format3ReadOptions read;
        read.readDataFiles = false;  // the names come from the records and <levels>, not from the data
        read.setUpChecks = false;
        read.allRecordsActive = allActive;
        const ex::Format3Project c = ex::ReadFormat3Text(text, read);
        if (!c.Ok() || !c.levelsRead) throw std::invalid_argument("parameter_names: the project has errors or no <levels>");
        return ex::ParameterNamesOf(c.built);
      },
      py::arg("text"), py::arg("all_active") = false,
      "Every fit-parameter name the format-3 project text gives that is made from a record's name, in record order "
      "(experiment::ParameterNamesOf); with all_active the names of inactive data sets and extrapolations too, as "
      "when they are switched on (BuildOptions::allRecordsActive).  ValueError when the project has an error.  pyazr's "
      "Project.rename maps the names before and after a rename by position.");
  m.def(
      "parameter_value_edits",
      [](const std::string &text, const std::string &projectDir, const std::map<std::string, double> &values) {
        ex::Format3ReadOptions read;
        read.projectDir = projectDir;
        const ex::Format3Project c = ex::ReadFormat3Text(text, read);
        py::list edits, missing;
        for (const auto &v : values) {
          const auto source = c.built.parameterSources.find(v.first);
          std::string why = "no record field holds this parameter's value";
          ex::ValueEdit e;
          if (source != c.built.parameterSources.end() && ex::ParameterValueEdit(c.validated, source->second, v.second, e, why))
            edits.append(py::make_tuple(e.line, e.col, py::bytes(e.old), py::bytes(e.text)));
          else
            missing.append(py::make_tuple(v.first, why));
        }
        return py::make_tuple(edits, missing, c.Ok());
      },
      py::arg("text"), py::arg("project_dir"), py::arg("values"),
      "The edits that write parameter values (by name, in the engine's units) into the nuisances of the format-3 "
      "project text they come from (experiment::ParameterValueEdit): (edits [(line, col, old token, new token)], "
      "names written nowhere [(name, why)], the project read without error).");
}
