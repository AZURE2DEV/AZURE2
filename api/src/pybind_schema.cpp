// The table of project format 3 for pyazr and the web UI (docs/upgrade/project_format3.md sections 6.1, 7.3): the
// same JSON document as AZURE2 --schema (experiment/'s Schema::ToJson), returned as a dict.  pyazr spells no field,
// quantity or unit itself; it reads them from here.

#include <pybind11/pybind11.h>

#include "azure2/experiment/Schema.h"

namespace py = pybind11;

void BindSchema(py::module_ &m) {
  m.def(
      "schema", []() { return py::module_::import("json").attr("loads")(azure2::experiment::Schema::ToJson()); },
      "The table of project format 3 as a dict: record kinds and fields, quantities, units, diagnostic codes "
      "(the JSON printed by AZURE2 --schema).");
  m.def("schema_json", []() { return azure2::experiment::Schema::ToJson(); },
        "The same table as the JSON text, byte for byte as AZURE2 --schema prints it (without the final newline).");
}
