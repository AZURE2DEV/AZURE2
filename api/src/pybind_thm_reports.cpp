// The THM reports of a pyazr session as Python dicts (Session::thm_experiments, thm_lineshape, thm_vertex and
// thm_distortion of pybind_api.cpp, which ask the engine and refuse a failed call): the experiments as the last
// chi-squared evaluation profiled them, the Coulomb line shape, the window-averaged entrance vertex and the
// distortion factor R(E) (include/ThmReports.h).  Copies of stored numbers; no arithmetic.
#include <string>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "AZUREAPI.h"
#include "ThmExperiment.h"
#include "ThmReports.h"
#include "pybind_arrays.h"

namespace py = pybind11;

// The THM experiments as the last chi-squared evaluation profiled them.
py::list ThmExperimentsList(AZUREAPI &api) {
  py::list out;
  for (const ThmExperimentReport &r : api.GetThmExperiments()) {
    py::dict d;
    d["name"] = r.name;
    d["segments"] = r.segments;
    d["background"] = r.background;
    d["points"] = r.points;
    d["chi2"] = r.chi2;
    d["status"] = r.status;
    d["value"] = std::vector<double>(r.value, r.value + 4);
    d["covariance"] = std::vector<double>(r.covariance, r.covariance + 16);
    if (!r.coherentNames.empty()) {
      d["cbkg_names"] = r.coherentNames;
      d["cbkg_values"] = r.coherentValues;
    }
    out.append(d);
  }
  return out;
}

// The Coulomb line shape of a THM experiment at the current parameters.
py::dict ThmLineshapeDict(const ThmLineshapeReport &r) {
  py::dict d;
  d["experiment"] = r.experiment;
  d["spectator"] = r.spectator;
  d["Zs"] = r.Zs;
  d["ZF"] = r.ZF;
  d["E_aA"] = r.eAA;
  d["B"] = r.bind;
  d["E"] = to_array(r.energy);
  d["E_sF"] = to_array(r.esf);
  d["eta_0"] = to_array(r.eta0);
  py::list exits;
  for (const ThmLineshapeReport::Exit &x : r.exits) {
    py::dict xd;
    xd["pair"] = x.pairKey;
    xd["Zb"] = x.Zb;
    xd["ZB"] = x.ZB;
    xd["mb"] = x.mb;
    xd["mB"] = x.mB;
    xd["zeta"] = to_array(x.zeta);
    xd["eta_sb"] = to_array(x.etaSb);
    py::list levels;
    for (const ThmLineshapeReport::Level &l : x.levels) {
      py::dict ld;
      ld["jgroup"] = l.jgroup;
      ld["level"] = l.level;
      ld["J"] = l.J;
      ld["pi"] = l.pi;
      ld["E_level"] = l.energy;
      ld["Gamma"] = l.width;
      ld["NC2"] = to_array(l.nc2);
      levels.append(ld);
    }
    xd["levels"] = levels;
    exits.append(xd);
  }
  d["exits"] = exits;
  return d;
}

// The window-averaged THM entrance vertex of an experiment at the current parameters.
py::dict ThmVertexDict(const ThmVertexReport &r) {
  py::dict d;
  d["experiment"] = r.experiment;
  d["window"] = r.window;
  d["pair"] = r.pairKey;
  d["mu_sx"] = r.muSx;
  d["B"] = r.bind;
  d["radius"] = r.radius;
  // The nodes at each energy (they follow the accepted directions): lists
  // over E of arrays over the nodes; empty lists with the DW vertex.
  auto rows = [&](const std::vector<std::vector<double>> &v) {
    py::list l;
    for (const std::vector<double> &row : v) l.append(to_array(row));
    return l;
  };
  d["p_s"] = rows(r.p);
  d["weights"] = rows(r.weight);
  d["T_s"] = rows(r.es);
  if (!r.theta.empty()) d["theta_cm"] = rows(r.theta);
  d["E"] = to_array(r.energy);
  {
    // Where the window (or the DW grid) does not reach E, the engine uses
    // the nodes of the nearest data point / grid energy (as for a folding
    // sub-point): flagged False here.
    py::array_t<bool> reached(r.reached.size());
    auto m = reached.mutable_unchecked<1>();
    for (size_t i = 0; i < r.reached.size(); i++) m(i) = r.reached[i] != 0;
    d["reached"] = reached;
  }
  py::list rho;
  for (const std::vector<double> &row : r.rho) rho.append(to_array(row));
  d["rho"] = rho;
  d["model"] = r.model;
  if (r.model == "dw") {
    py::list q, w;
    for (const std::vector<double> &row : r.dwQ) q.append(to_array(row));
    for (const std::vector<double> &row : r.dwWeight) w.append(to_array(row));
    d["dw_q"] = q;
    d["dw_weights"] = w;
    d["dw_q_delta"] = to_array(r.dwQDelta);
    d["dw_p_delta"] = to_array(r.dwPDelta);
    if (!r.dwTheta.empty()) {
      py::list th, aq, aw;
      for (const std::vector<double> &row : r.dwTheta) th.append(to_array(row));
      for (const std::vector<double> &row : r.dwAngleQ) aq.append(to_array(row));
      for (const std::vector<double> &row : r.dwAngleWeight) aw.append(to_array(row));
      d["angle_theta_cm"] = th;
      d["angle_q"] = aq;
      d["angle_weights"] = aw;
    }
  }
  py::list channels;
  for (const ThmVertexReport::Channel &c : r.channels) {
    py::dict cd;
    cd["jgroup"] = c.jgroup;
    cd["channel"] = c.channel;
    cd["J"] = c.J;
    cd["pi"] = c.pi;
    cd["l"] = c.l;
    cd["s"] = c.s;
    py::list levels;
    for (const ThmVertexReport::Level &l : c.levels) {
      py::dict ld;
      ld["level"] = l.level;
      ld["boundary"] = l.boundary;
      ld["M2"] = to_array(l.m2);
      ld["M2_qf"] = to_array(l.m2qf);
      if (!l.m2pw.empty()) ld["M2_pw"] = to_array(l.m2pw);
      levels.append(ld);
    }
    cd["levels"] = levels;
    channels.append(cd);
  }
  d["channels"] = channels;
  return d;
}

// The distortion factor R(E) of a THM experiment (distortion=...).
py::dict ThmDistortionDict(const ThmDistortionReport &r) {
  py::dict d;
  d["experiment"] = r.experiment;
  d["kind"] = r.kind;
  d["description"] = r.description;
  d["E"] = to_array(r.energy);
  d["R_model"] = to_array(r.rModel);
  if (r.kind == "table") return d;
  d["E_ref"] = r.eRef;
  d["E_aA"] = r.eAA;
  d["B"] = r.bind;
  d["k_aA"] = r.kAA;
  d["eta_aA"] = r.etaAA;
  d["kappa"] = r.kappa;
  d["eta_b"] = r.etaB;
  d["beta"] = r.beta;
  d["E_sF"] = to_array(r.esf);
  d["k_sF"] = to_array(r.ksf);
  d["eta_sF"] = to_array(r.etasf);
  d["theta_cm"] = to_array(r.thetaCm);
  d["x"] = to_array(r.x);
  d["q"] = to_array(r.q);
  d["M2"] = to_array(r.m2);
  d["M2_PW"] = to_array(r.mpw2);
  d["R"] = to_array(r.r);
  d["lmax"] = r.lmax;
  return d;
}
