// Free functions of physics/ in the _azure2 module, with array inputs (docs/upgrade/ARCHITECTURE.md section 8: every
// physics function is plottable).  No session, no project: each call is the physics function applied element-wise.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <initializer_list>
#include <stdexcept>
#include <string>
#include <vector>

#include "azure2/physics/Constants.h"
#include "azure2/physics/AstrophysicalFactor.h"
#include "azure2/physics/CoulombScattering.h"
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

/// The kinematics of a pair from its charges and atomic masses (u), as the engine builds it (Kinematics::Channel).
azure2::physics::ChannelKinematics PairKinematics(int z1, int z2, double m1_u, double m2_u, bool relativistic) {
  azure2::physics::ChannelKinematics k;
  k.z1 = z1;
  k.z2 = z2;
  k.mu_u = m1_u * m2_u / (m1_u + m2_u);
  k.m1_MeV = m1_u * azure2::physics::kAtomicMassUnit_MeV;
  k.m2_MeV = m2_u * azure2::physics::kAtomicMassUnit_MeV;
  k.relativistic = relativistic;
  return k;
}

/// The Rutherford (two_s < 0) or Mott (two_s >= 0) cross section and d/dE, element-wise on E_cm and theta_cm.
py::dict CoulombCrossSections(const char *what, const azure2::physics::ChannelKinematics &kin, int two_s,
                              const Array &e_cm, const Array &theta_cm) {
  const auto E = e_cm.unchecked<1>(), T = theta_cm.unchecked<1>();
  const py::ssize_t n = Common(what, {E.shape(0), T.shape(0)});
  py::array_t<double> value(n), dE(n);
  auto V = value.mutable_unchecked<1>(), D = dE.mutable_unchecked<1>();
  for (py::ssize_t i = 0; i < n; i++) {
    const azure2::physics::CoulombCrossSection c = two_s < 0 ? azure2::physics::rutherford_cross_section(kin, E(i), T(i))
                                                             : azure2::physics::mott_cross_section(kin, two_s, E(i), T(i));
    V(i) = c.value;
    D(i) = c.d_dE;
  }
  py::dict out;
  out["value"] = value;
  out["d_dE"] = dE;
  return out;
}

/// E exp(2 pi eta(E)) element-wise on E_cm (physics/AstrophysicalFactor.h).
py::array_t<double> SFactorConversions(const azure2::physics::ChannelKinematics &kin, const Array &e_cm) {
  const auto E = e_cm.unchecked<1>();
  py::array_t<double> out(E.shape(0));
  auto O = out.mutable_unchecked<1>();
  for (py::ssize_t i = 0; i < E.shape(0); i++) O(i) = azure2::physics::sfactor_conversion(kin, E(i));
  return out;
}

}  // namespace

void BindStopping(py::module_ &parent);  // pybind_stopping.cpp: stopping/, the submodule _azure2.stopping

void BindPhysics(py::module_ &m) {
  BindStopping(m);
  m.def(
      "sfactor_conversion",
      [](int z1, int z2, double m1_u, double m2_u, const Array &e_cm) {
        return SFactorConversions(PairKinematics(z1, z2, m1_u, m2_u, false), e_cm);
      },
      py::arg("z1"), py::arg("z2"), py::arg("m1_u"), py::arg("m2_u"), py::arg("e_cm"),
      "E exp(2 pi eta(E)) (MeV) of a pair of charges z1, z2 and atomic masses m1_u, m2_u (u), element-wise on E_cm "
      "(MeV): the factor S = sigma x this (S in MeV b, sigma in b; physics/AstrophysicalFactor.h, the engine's "
      "quantity sfactor).  eta is non-relativistic; NaN for E <= 0 or an uncharged pair.");
  m.def(
      "modified_sfactor_conversion",
      [](int z1, int z2, double m1_u, double m2_u, double g_per_MeV, const Array &e_cm) {
        const azure2::physics::ChannelKinematics kin = PairKinematics(z1, z2, m1_u, m2_u, false);
        const auto E = e_cm.unchecked<1>();
        py::array_t<double> out(E.shape(0));
        auto O = out.mutable_unchecked<1>();
        for (py::ssize_t i = 0; i < E.shape(0); i++)
          O(i) = azure2::physics::modified_sfactor_conversion(kin, E(i), g_per_MeV);
        return out;
      },
      py::arg("z1"), py::arg("z2"), py::arg("m1_u"), py::arg("m2_u"), py::arg("g_per_MeV"), py::arg("e_cm"),
      "E exp(2 pi eta(E) + g E) (MeV), element-wise on E_cm (MeV): the factor S* = sigma x this of the modified S "
      "factor of heavy-ion fusion (physics/AstrophysicalFactor.h, the engine's quantity modified_sfactor; g in 1/MeV, "
      "0.46 for 12C + 12C).  NaN for E <= 0 or an uncharged pair.");
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
  m.def(
      "rutherford_cross_section",
      [](int z1, int z2, double m1_u, double m2_u, const Array &e_cm, const Array &theta_cm, bool relativistic) {
        return CoulombCrossSections("rutherford_cross_section", PairKinematics(z1, z2, m1_u, m2_u, relativistic), -1,
                                    e_cm, theta_cm);
      },
      py::arg("z1"), py::arg("z2"), py::arg("m1_u"), py::arg("m2_u"), py::arg("e_cm"), py::arg("theta_cm"),
      py::arg("relativistic") = false,
      "Rutherford cross section (b/sr, c.m.) of a pair of charges z1, z2 and atomic masses m1_u, m2_u (u), with the "
      "engine's eta and k (physics/CoulombScattering.h), element-wise on equal-length arrays of E_cm (MeV) and "
      "theta_cm (deg): dict of value and d_dE (b/sr/MeV, at fixed angle).  NaN outside E > 0, 0 < theta < 180.");
  m.def(
      "mott_cross_section",
      [](int z, double m_u, int two_s, const Array &e_cm, const Array &theta_cm, bool relativistic) {
        if (two_s < 0) throw std::invalid_argument("mott_cross_section: two_s must be >= 0");
        return CoulombCrossSections("mott_cross_section", PairKinematics(z, z, m_u, m_u, relativistic), two_s, e_cm,
                                    theta_cm);
      },
      py::arg("z"), py::arg("m_u"), py::arg("two_s"), py::arg("e_cm"), py::arg("theta_cm"),
      py::arg("relativistic") = false,
      "Mott cross section (b/sr, c.m.) of two identical particles of charge z, atomic mass m_u (u) and spin two_s/2: "
      "the symmetrised Coulomb amplitude with the interference weight (-1)^(2s)/(2s+1) (physics/CoulombScattering.h), "
      "element-wise on E_cm (MeV) and theta_cm (deg): dict of value and d_dE (b/sr/MeV).");
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
