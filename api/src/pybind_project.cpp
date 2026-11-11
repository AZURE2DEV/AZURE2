// The project file for pyazr (docs/upgrade/project.md section 2, project_format3.md section 4): the lossless document
// io::Document that pyazr's AzrModel and Project read and edit, as the engine reads the file.

#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "azure2/io/Convert.h"
#include "azure2/io/Document.h"

namespace py = pybind11;

void BindProject(py::module_ &m) {
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
}
