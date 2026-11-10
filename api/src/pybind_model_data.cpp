// The <modelData> data of a session as Python records (AZUREAPI::ModelDataRWA): the "model_data" entry of
// Session.objective_terms_rwa, one dict per datum, read by pyazr.model_data.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "AZUREAPI.h"

namespace py = pybind11;

void AddModelData(py::dict &terms, const AZUREAPI &api, const vector_r &free) {
  py::list out;
  for (const AZUREAPI::ModelDatumValue &v : api.ModelDataRWA(free)) {
    py::dict r;
    r["name"] = v.name;
    r["quantity"] = v.quantity;
    r["unit"] = v.unit;
    r["source"] = v.source;
    r["model"] = v.model;
    r["value"] = v.value;
    r["error"] = v.error;
    r["residual"] = v.residual;
    if (v.hasImag) {
      r["model_imag"] = v.modelImag;
      r["imag"] = v.imag;
      r["imag_error"] = v.imagError;
      r["residual_imag"] = v.residualImag;
    }
    r["chi2"] = v.chi2;
    out.append(r);
  }
  terms["model_data"] = out;
}
