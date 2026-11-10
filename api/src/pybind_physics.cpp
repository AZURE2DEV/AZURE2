// Free functions of physics/ in the _azure2 module, with array inputs (docs/upgrade/ARCHITECTURE.md section 8: every
// physics function is plottable).  No session, no project: each call is the physics function applied element-wise.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <initializer_list>
#include <stdexcept>
#include <string>
#include <vector>

#include "azure2/physics/Elements.h"
#include "azure2/physics/ResonanceQuantities.h"
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

using Array = py::array_t<double, py::array::forcecast>;

/// n of equal-length 1-d arrays, or an error naming the function.
py::ssize_t Common(const char *what, std::initializer_list<py::ssize_t> sizes) {
  const py::ssize_t n = *sizes.begin();
  for (py::ssize_t s : sizes)
    if (s != n) throw std::invalid_argument(std::string(what) + ": the arrays must have the same length");
  return n;
}

/// omega gamma = omega Gamma_a Gamma_b / Gamma and its partial derivatives, element-wise (ResonanceQuantities.h).
py::dict ResonanceStrengths(const Array &omega, const Array &a, const Array &b, const Array &g) {
  const auto O = omega.unchecked<1>(), A = a.unchecked<1>(), B = b.unchecked<1>(), G = g.unchecked<1>();
  const py::ssize_t n = Common("resonance_strength", {O.shape(0), A.shape(0), B.shape(0), G.shape(0)});
  py::array_t<double> value(n), dA(n), dB(n), dG(n);
  auto V = value.mutable_unchecked<1>(), DA = dA.mutable_unchecked<1>(), DB = dB.mutable_unchecked<1>(),
       DG = dG.mutable_unchecked<1>();
  for (py::ssize_t i = 0; i < n; i++) {
    const azure2::physics::ResonanceStrength r = azure2::physics::resonance_strength(O(i), A(i), B(i), G(i));
    V(i) = r.value;
    DA(i) = r.d_entrance;
    DB(i) = r.d_exit;
    DG(i) = r.d_total;
  }
  py::dict out;
  out["value"] = value;
  out["d_gamma_a"] = dA;
  out["d_gamma_b"] = dB;
  out["d_gamma"] = dG;
  return out;
}

/// Gamma_b / Gamma and its partial derivatives, element-wise.
py::dict BranchingRatios(const Array &b, const Array &g) {
  const auto B = b.unchecked<1>(), G = g.unchecked<1>();
  const py::ssize_t n = Common("branching_ratio", {B.shape(0), G.shape(0)});
  py::array_t<double> value(n), dB(n), dG(n);
  auto V = value.mutable_unchecked<1>(), DB = dB.mutable_unchecked<1>(), DG = dG.mutable_unchecked<1>();
  for (py::ssize_t i = 0; i < n; i++) {
    const azure2::physics::BranchingRatio r = azure2::physics::branching_ratio(B(i), G(i));
    V(i) = r.value;
    DB(i) = r.d_exit;
    DG(i) = r.d_total;
  }
  py::dict out;
  out["value"] = value;
  out["d_gamma_b"] = dB;
  out["d_gamma"] = dG;
  return out;
}

/// omega of eq. (8), element-wise on 2J of the level and of the two entrance particles.
py::array_t<double> StatisticalFactors(const py::array_t<int, py::array::forcecast> &twoJ,
                                       const py::array_t<int, py::array::forcecast> &twoJa,
                                       const py::array_t<int, py::array::forcecast> &twoJA,
                                       const py::array_t<bool, py::array::forcecast> &identical) {
  const auto J = twoJ.unchecked<1>(), Ja = twoJa.unchecked<1>(), JA = twoJA.unchecked<1>();
  const auto I = identical.unchecked<1>();
  const py::ssize_t n = Common("statistical_factor", {J.shape(0), Ja.shape(0), JA.shape(0), I.shape(0)});
  py::array_t<double> out(n);
  auto O = out.mutable_unchecked<1>();
  for (py::ssize_t i = 0; i < n; i++) O(i) = azure2::physics::statistical_factor(J(i), Ja(i), JA(i), I(i));
  return out;
}

}  // namespace

void BindPhysics(py::module_ &m) {
  m.def("statistical_factor", &StatisticalFactors, py::arg("two_J"), py::arg("two_ja"), py::arg("two_jA"),
        py::arg("identical"),
        "omega = (2J+1)(1+delta_aA)/((2j_a+1)(2j_A+1)) of a resonance strength (physics/ResonanceQuantities.h), "
        "element-wise on equal-length arrays of 2J (level), 2j_a, 2j_A (entrance particles) and identical flags.");
  m.def("resonance_strength", &ResonanceStrengths, py::arg("omega"), py::arg("gamma_a"), py::arg("gamma_b"),
        py::arg("gamma"),
        "omega gamma = omega Gamma_a Gamma_b / Gamma (data_kinds.md eq. 8) from observed widths, element-wise: dict of "
        "value and the partial derivatives d_gamma_a, d_gamma_b, d_gamma (the three widths independent).");
  m.def("branching_ratio", &BranchingRatios, py::arg("gamma_b"), py::arg("gamma"),
        "Gamma_b / Gamma, element-wise: dict of value, d_gamma_b and d_gamma.");
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
