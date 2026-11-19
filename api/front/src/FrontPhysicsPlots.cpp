// What the Physics tab draws from the per-session calls (docs/upgrade/qt_gui.md section 15.2): the calls in words,
// the plots each answer offers (the series names are those FrontPhysics.cpp gives), and the last fit's parameter file.
// No arithmetic: the series' values are the engine's.
#include "azure2/front/Engine.h"

#include <filesystem>

#include "azure2/experiment/ProjectPaths.h"
#include "azure2/physics_table/PhysicsTable.h"
#include "azure2/physics_table/StoppingTable.h"

namespace azure2 {
namespace front {

std::vector<PhysicsCall> PhysicsCalls() {
  std::vector<PhysicsCall> out;
  auto add = [&out](const char *name, const char *title, const char *doc, const char *axis, bool slow) {
    PhysicsCall c;
    c.name = name, c.title = title, c.doc = doc, c.axisTitle = axis, c.slow = slow;
    out.push_back(c);
  };
  add("ChannelFunctions", "Channel functions of a pair",
      "The session's Coulomb functions F, G and their derivatives at the channel radius, the penetrability P, the "
      "shift function S and the hard-sphere phase of one pair and l (the Coulomb table and a pair's nuclear potential "
      "as the calculation uses them).",
      "E_c, c.m. energy of the pair (MeV)", false);
  add("Block", "One J^pi group of the R-matrix",
      "The level matrix A, the collision matrix U and T = exp(2i omega) - U of one J^pi group at the session's "
      "parameters, with the channels' P, S and L - B, the nuclear phase shifts and the inelasticities "
      "(docs/upgrade/physics_tab.md section 4).",
      "E_x, excitation energy of the compound nucleus (MeV)", false);
  add("ThmVertex", "THM vertex",
      "The entrance vertex of a Trojan-horse experiment: <|M_l|^2> over the spectator-momentum window, |M_l|^2 at "
      "p_s = 0, and the plane-wave |M_l(p)|^2 of vertexModel=dw, per entrance channel and level.",
      "E, c.m. energy of x + A (MeV)", true);
  add("ThmLineshape", "THM Coulomb line shape",
      "The Coulomb line shape |N_C|^2 per exit pair and level (lineshape=on), zeta = eta_sB - eta_0 and the "
      "neglected eta_sb.  One chi^2 evaluation runs first.",
      "E, c.m. energy of x + A (MeV)", true);
  add("ThmDistortion", "THM distortion",
      "The distortion factor R(E) of the experiment and its parts (distortion=): a DWBA per energy, about 17 ms "
      "each.",
      "E, c.m. energy of x + A (MeV)", true);
  return out;
}

std::vector<PhysicsPlot> PhysicsPlots(const std::string &call) {
  std::vector<PhysicsPlot> out;
  auto add = [&out](const char *c, const char *name, const char *title, const char *y, const char *yTitle,
                    const char *doc) -> PhysicsPlot & {
    PhysicsPlot p;
    p.call = c, p.name = name, p.title = title, p.y = y, p.yTitle = yTitle, p.doc = doc;
    out.push_back(p);
    return out.back();
  };
  // ChannelFunctions
  add("ChannelFunctions", "P", "Penetrability P_l", "P", "P_l", "P_l = rho / (F_l^2 + G_l^2) at the channel radius.")
      .logY = true;
  add("ChannelFunctions", "S", "Shift function S_l", "S", "S_l",
      "S_l = rho (F F' + G G') / (F^2 + G^2); below threshold the closed channel's S_l.");
  add("ChannelFunctions", "delta_hs", "Hard-sphere phase", "delta_hs", "delta_hs (rad)",
      "The hard-sphere phase shift -phi_l = -arg(G_l + i F_l) (rad).");
  add("ChannelFunctions", "F", "Regular Coulomb function F_l", "F", "F_l", "F_l at the channel radius.");
  add("ChannelFunctions", "G", "Irregular Coulomb function G_l", "G", "G_l", "G_l at the channel radius.").logY = true;
  add("ChannelFunctions", "Fp", "Derivative F'_l", "Fp", "F'_l", "dF_l/drho at the channel radius.");
  add("ChannelFunctions", "Gp", "Derivative G'_l", "Gp", "G'_l", "dG_l/drho at the channel radius.");
  // Block
  add("Block", "delta", "Phase shifts delta_c", "delta", "delta (deg)",
      "The nuclear phase shift arg(U_cc)/2 of each open particle channel, in degrees, as the program's phase_shift "
      "data kind: in (-90, 90], [0, 180) for identical particles.")
      .diagonalOnly = true;
  add("Block", "eta", "Inelasticities eta_c = |U_cc|", "eta", "eta = |U_cc|",
      "1 where the channel's flux stays in it; below 1 where other channels take flux.")
      .diagonalOnly = true;
  add("Block", "U_abs", "|U_cc'|", "U_abs", "|U_cc'|", "The modulus of the collision matrix, every pair c <= c'.");
  {
    PhysicsPlot &a = add("Block", "argand", "Argand plot of U_cc", "U_im", "Im U_cc",
                         "Im U_cc against Re U_cc with the excitation energy as the curve's parameter; a resonance "
                         "is a loop; the markers sit at the levels' energies.");
    a.x = "U_re", a.xTitle = "Re U_cc", a.diagonalOnly = true;
  }
  add("Block", "T_abs", "|T_cc'|", "T_abs", "|T_cc'|", "T = exp(2i omega) - U, every pair c <= c'.").logY = true;
  add("Block", "U_re", "Re U_cc'", "U_re", "Re U_cc'", "The real part of the collision matrix.");
  add("Block", "U_im", "Im U_cc'", "U_im", "Im U_cc'", "The imaginary part of the collision matrix.");
  add("Block", "P", "Penetrability P_c per channel", "P", "P_c", "Each channel's penetrability at its E_c.").logY = true;
  add("Block", "S", "Shift function S_c per channel", "S", "S_c", "Each channel's shift function at its E_c.");
  add("Block", "A_abs", "|A_lambda mu| (level matrix)", "A_abs", "|A| (MeV^-1)",
      "The modulus of the level matrix A, every pair of levels lambda <= mu.")
      .logY = true;
  add("Block", "log10_abs_det", "log10 |det A^-1|", "log10_abs_det", "log10 |det A^-1|",
      "The determinant of the inverse level matrix: its minima on the real axis are where A^-1 is nearly singular, "
      "the resonance energies of narrow levels.");
  // the maps (Pole map page)
  add("Block", "map_det", "log10 |det A^-1|, every J^pi group", "log10_abs_det", "log10 |det A^-1|",
      "One row per J^pi group, one Block per group.  A dark valley is a minimum of |det A^-1|: a resonance.")
      .mapRows = "jgroup";
  add("Block", "map_U", "|U_cc'| of one J^pi group", "U_abs", "|U_cc'|",
      "One row per pair of open particle channels c <= c' of the group.")
      .mapRows = "series";
  {
    PhysicsPlot &m = add("Block", "map_A", "|A_lambda mu| of one J^pi group", "A_abs", "|A| (MeV^-1)",
                         "One row per pair of levels lambda <= mu of the group; a peak of |A| at a level's energy.");
    m.mapRows = "series", m.logY = true;
  }
  // THM
  add("ThmVertex", "M2", "<|M_l|^2> over the window", "M2", "<|M_l|^2>",
      "The vertex averaged over the spectator-momentum window, per entrance channel and level.")
      .logY = true;
  add("ThmVertex", "M2_qf", "|M_l|^2 at p_s = 0", "M2_qf", "|M_l|^2 (p_s = 0)", "The vertex at the quasi-free point.")
      .logY = true;
  add("ThmVertex", "M2_pw", "Plane-wave |M_l(p)|^2", "M2_pw", "|M_l(p)|^2 (plane wave)",
      "Only for vertexModel=dw: the plane-wave vertex for comparison.")
      .logY = true;
  add("ThmVertex", "reached", "Window reaches E", "reached", "reached (1 or 0)",
      "0 where the spectator window does not reach the energy.");
  add("ThmLineshape", "NC2", "Coulomb line shape |N_C|^2", "NC2", "|N_C|^2", "Per exit pair and level.");
  add("ThmLineshape", "zeta", "zeta = eta_sB - eta_0", "zeta", "zeta", "Per exit pair.");
  add("ThmLineshape", "eta_sb", "eta_sb (neglected)", "eta_sb", "eta_sb", "Per exit pair.");
  add("ThmLineshape", "E_sF", "E_sF", "E_sF", "E_sF (MeV)", "The spectator-fragment energy.");
  add("ThmDistortion", "R", "Distortion factor R(E)", "R", "R", "The distortion factor of the experiment.");
  add("ThmDistortion", "R_model", "R of the model", "R_model", "R (model)", "The model's distortion factor.");
  add("ThmDistortion", "M2", "|M|^2 (DWBA)", "M2", "|M|^2 (fm^6)", "The distorted-wave vertex.").logY = true;
  add("ThmDistortion", "M2_PW", "|M|^2 (plane wave)", "M2_PW", "|M|^2 PW (fm^6)", "The plane-wave vertex.").logY = true;
  add("ThmDistortion", "theta_cm", "theta_cm", "theta_cm", "theta_cm (deg)", "The c.m. angle of the kinematics.");
  add("ThmDistortion", "q", "Momentum transfer q", "q", "q (fm^-1)", "The momentum transfer.");
  if (call.empty()) return out;
  std::vector<PhysicsPlot> of;
  for (const PhysicsPlot &p : out)
    if (p.call == call) of.push_back(p);
  return of;
}

std::string StoppingStatusText(const std::string &name, int status) {
  namespace table = azure2::physics_table;
  if (status == table::kOk) return "computed";
  if (status == table::kFailed) return "not converged or failed";
  if (status == table::kJointly) return "the arguments together are outside the domain (a layer at least the range)";
  const table::StoppingRow *row = table::FindStoppingRow(name);
  if (row) {
    const std::vector<table::ArgSpec> spec = table::StoppingNumericArgs(*row);
    if (status >= 1 && status <= (int)spec.size())
      return spec[status - 1].name + " outside its domain (" + table::DomainText(spec[status - 1]) + ")";
  }
  return "status " + std::to_string(status);
}

std::string LastParameterFile(const Project &project) {
  std::string dir = experiment::PathToOpen(project.Directory(), project.ConfigValue("output_dir"));
  if (dir.empty()) return std::string();
  if (dir.back() != '/') dir += '/';
  // the name AZUREParams::WriteUserParameters gives it after a fit (src/AZUREParams.cpp: "<outputdir>param.sav")
  const std::string file = dir + "param.sav";
  std::error_code ec;
  return std::filesystem::is_regular_file(file, ec) ? file : std::string();
}

}  // namespace front
}  // namespace azure2
