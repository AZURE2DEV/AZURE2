// The binding's dict of one J^pi block on an energy grid (JpiBlockReport.h, AZUREAPI::JpiBlockAllRWA): pyazr's
// Session.jpi_block (pyazr/session_physics.py) and its row of pyazr/session_functions.py.  Plain copies of the engine's
// arrays; no physics here.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <complex>
#include <vector>

#include "JpiBlockReport.h"

namespace py = pybind11;

namespace {

py::array_t<double> Column(const std::vector<double> &a, std::size_t stride, std::size_t k, std::size_t n) {
  py::array_t<double> out(static_cast<py::ssize_t>(n));
  double *o = out.mutable_data();
  for (std::size_t e = 0; e < n; e++) o[e] = a[e * stride + k];
  return out;
}

template <class F>
py::array_t<double> ComplexColumn(const std::vector<std::complex<double>> &a, std::size_t stride, std::size_t k,
                                  std::size_t n, F part) {
  py::array_t<double> out(static_cast<py::ssize_t>(n));
  double *o = out.mutable_data();
  for (std::size_t e = 0; e < n; e++) o[e] = part(a[e * stride + k]);
  return out;
}

}  // namespace

py::dict JpiBlockDict(const JpiBlockReport &r) {
  const std::size_t ne = r.ex.size(), nc = r.channels.size(), n = r.levels.size();
  const auto re = [](std::complex<double> z) { return z.real(); };
  const auto im = [](std::complex<double> z) { return z.imag(); };
  const auto ab = [](std::complex<double> z) { return std::abs(z); };
  py::dict d;
  d["jgroup"] = r.jgroup;
  d["J"] = r.J;
  d["pi"] = r.pi;
  d["parametrization"] = r.park ? "park" : r.brune ? "brune" : "standard";
  d["E_x"] = py::array_t<double>(static_cast<py::ssize_t>(ne), r.ex.data());
  py::list levels;
  for (std::size_t i = 0; i < n; i++) {
    py::dict l;
    l["level"] = r.levels[i];
    l["E_level"] = r.levelEnergy[i];
    levels.append(l);
  }
  d["levels"] = levels;
  py::list channels, elements, phases;
  for (std::size_t c = 0; c < nc; c++) {
    const JpiBlockReport::Channel &k = r.channels[c];
    py::dict ch;
    ch["channel"] = k.channel;
    ch["pair"] = k.pairKey;
    ch["radiation"] = std::string(1, k.radType);
    ch["l"] = k.l;
    ch["s"] = k.s;
    ch["radius"] = k.radius;
    ch["threshold"] = k.threshold;
    ch["boundary"] = k.boundary;
    ch["E_c"] = Column(r.ec, nc, c, ne);
    ch["P"] = Column(r.P, nc, c, ne);
    ch["S"] = Column(r.S, nc, c, ne);
    ch["L_re"] = ComplexColumn(r.lo, nc, c, ne, re);
    ch["L_im"] = ComplexColumn(r.lo, nc, c, ne, im);
    channels.append(ch);
    if (k.radType == 'P') {
      py::dict p;
      p["channel"] = k.channel;
      p["pair"] = k.pairKey;
      p["l"] = k.l;
      p["s"] = k.s;
      p["delta"] = Column(r.delta, nc, c, ne);
      p["eta"] = Column(r.eta, nc, c, ne);
      phases.append(p);
    }
    for (std::size_t cp = c; cp < nc; cp++) {
      if (k.radType != 'P' || r.channels[cp].radType != 'P') continue;
      py::dict e;
      e["channel"] = k.channel;
      e["channel_out"] = r.channels[cp].channel;
      const std::size_t at = c * nc + cp;
      e["U_re"] = ComplexColumn(r.U, nc * nc, at, ne, re);
      e["U_im"] = ComplexColumn(r.U, nc * nc, at, ne, im);
      e["U_abs"] = ComplexColumn(r.U, nc * nc, at, ne, ab);
      e["T_re"] = ComplexColumn(r.T, nc * nc, at, ne, re);
      e["T_im"] = ComplexColumn(r.T, nc * nc, at, ne, im);
      elements.append(e);
    }
  }
  d["channels"] = channels;
  d["elements"] = elements;
  d["phase_shifts"] = phases;
  py::list pairs;
  for (std::size_t l = 0; l < n; l++)
    for (std::size_t m = l; m < n; m++) {
      py::dict a;
      a["level"] = r.levels[l];
      a["level_out"] = r.levels[m];
      a["A_abs"] = ComplexColumn(r.A, n * n, l * n + m, ne, ab);
      pairs.append(a);
    }
  d["level_pairs"] = pairs;
  d["det_re"] = ComplexColumn(r.detInverse, 1, 0, ne, re);
  d["det_im"] = ComplexColumn(r.detInverse, 1, 0, ne, im);
  d["log10_abs_det"] = ComplexColumn(r.detInverse, 1, 0, ne, [](std::complex<double> z) { return std::log10(std::abs(z)); });
  py::array_t<bool> singular(static_cast<py::ssize_t>(ne));
  for (std::size_t e = 0; e < ne; e++) singular.mutable_data()[e] = r.singular[e] != 0;
  d["singular"] = singular;
  return d;
}
