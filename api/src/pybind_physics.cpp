// Free functions of physics/ in the _azure2 module, with array inputs (docs/upgrade/ARCHITECTURE.md section 8: every
// physics function is plottable).  No session, no project: each call is the physics function applied element-wise.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "azure2/physics/Elements.h"
#include "azure2/physics/coulomb/Whittaker.h"

namespace py = pybind11;
namespace coul = azure2::physics::coulomb;

namespace {

/// The Whittaker function W_{-eta_b, l+1/2}(2 rho) on equal-length arrays l, eta_b, rho (Whittaker.h): W = mantissa
/// 10^log10_scale, rho W'/W, and status (0 ok, 1 invalid_argument, 2 failed); NaN where the status is not ok.
py::dict Whittaker(const py::array_t<int, py::array::forcecast> &l, const py::array_t<double, py::array::forcecast> &eta_b,
                   const py::array_t<double, py::array::forcecast> &rho) {
  const auto L = l.unchecked<1>();
  const auto E = eta_b.unchecked<1>();
  const auto R = rho.unchecked<1>();
  const py::ssize_t n = L.shape(0);
  if (E.shape(0) != n || R.shape(0) != n)
    throw std::invalid_argument("whittaker: l, eta_b and rho must have the same length (got " + std::to_string(n) + ", " +
                                std::to_string(E.shape(0)) + ", " + std::to_string(R.shape(0)) + ")");
  py::array_t<double> mantissa(n), dlog(n);
  py::array_t<long long> scale(n);
  py::array_t<int> status(n);
  auto M = mantissa.mutable_unchecked<1>();
  auto D = dlog.mutable_unchecked<1>();
  auto S = scale.mutable_unchecked<1>();
  auto T = status.mutable_unchecked<1>();
  {
    py::gil_scoped_release release;
    for (py::ssize_t i = 0; i < n; i++) {
      const coul::WhittakerScaled w = coul::whittaker_scaled(L(i), E(i), R(i));
      M(i) = w.mantissa;
      D(i) = w.rho_dlog;
      S(i) = w.log10_scale;
      T(i) = static_cast<int>(w.status);
    }
  }
  py::dict out;
  out["mantissa"] = mantissa;
  out["log10_scale"] = scale;
  out["rho_dlog"] = dlog;
  out["status"] = status;
  return out;
}

}  // namespace

void BindPhysics(py::module_ &m) {
  m.def("whittaker", &Whittaker, py::arg("l"), py::arg("eta_b"), py::arg("rho"),
        "Whittaker function W_{-eta_b, l+1/2}(2 rho) of a bound channel (physics/coulomb, no GSL) on equal-length "
        "arrays: dict of mantissa, log10_scale (W = mantissa * 10**log10_scale), rho_dlog (rho W'/W, the shift "
        "function at rho = kappa a) and status (0 ok, 1 invalid_argument, 2 failed).");
  m.def(
      "element_symbols",
      []() {
        std::vector<std::string> out;
        for (int z = 0; z <= azure2::physics::kMaxElementZ; z++) out.push_back(azure2::physics::ElementSymbol(z));
        return out;
      },
      "Element symbols by charge number, index 0 (the neutron, \"n\") to 118 (physics/Elements.h, the one table).");
}
