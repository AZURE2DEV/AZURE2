// What the THM experiments report: output/thm_experiments.out (WorkingCopy::WriteThmExperiments) and the tables pyazr
// reads (ThmLineshapeTable, ThmDistortionTable, ThmVertexTable).  Members of WorkingCopy; the groups they read are built
// and profiled in ThmNorm.cpp.
#include <gsl/gsl_sf_bessel.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "CNuc.h"
#include "Config.h"
#include "Constants.h"
#include "WorkingCopy.h"
#include "PPair.h"
#include "ThmAngular.h"
#include "ThmDistortion.h"
#include "ThmDwVertex.h"
#include "ThmFunc.h"
#include "ThmLineshape.h"
#include "ThmReports.h"
#include "ThmVertexBoundary.h"

void WorkingCopy::WriteThmExperiments(const Config &configure, const azure2::Results &results) {
  std::string file = configure.outputdir + "thm_experiments.out";
  std::ofstream out(configure.OpenPath(file).c_str());
  if (!out) {
    configure.outStream << "Could not write " << file << "." << std::endl;
    return;
  }
  out << "# THM experiments (<thm> experiment[<name>] lines): segments sharing one profiled norm and\n"
         "# a background b(E) = b0 + b1 E + b2 E^2 added to the folded model (model units, E the c.m.\n"
         "# energy of the entrance pair in MeV).  norm multiplies the data, as the segment norms of\n"
         "# normalizations.out.  Uncertainties and covariance are those of the closed-form profile\n"
         "# at fixed R-matrix parameters (inverse normal matrix of the linear least squares, not\n"
         "# scaled by chi2/nu).\n";
  out.precision(10);
  out << std::scientific;
  for (const ThmExperimentReport &r : ThmExperimentReports(results)) {
    out << "\nexperiment: " << r.name << "\nsegments:";
    for (size_t k = 0; k < r.segments.size(); k++) out << (k ? "," : " ") << r.segments[k];
    out << "\nbackground: " << r.background << "\npoints: " << r.points << "\nchi2: " << r.chi2
        << "\nstatus: " << r.status << "\n";
    const char *names[4] = {"norm", "b0", "b1", "b2"};
    int q = 1 + (r.background == "none" ? 0 : r.background == "const" ? 1 : r.background == "linear" ? 2 : 3);
    for (int i = 0; i < q; i++)
      out << std::left << std::setw(6) << names[i] << std::right << std::setw(18) << r.value[i]
          << std::setw(18) << std::sqrt(std::max(0.0, r.covariance[i * 4 + i])) << "\n";
    out << "covariance:\n";
    for (int i = 0; i < q; i++) {
      out << std::left << std::setw(6) << names[i] << std::right;
      for (int j = 0; j < q; j++) out << std::setw(18) << r.covariance[i * 4 + j];
      out << "\n";
    }
    // Coherent background (cbackground=): its fit parameters (their errors
    // are the fit's, in param.par).
    if (!r.coherentNames.empty()) {
      out << "cbackground: " << r.coherentNames.size() << " parameter(s); c(E) = c0 + c1 E times the vertex M_l "
          << "is added to the resonant amplitude\n";
      for (size_t k = 0; k < r.coherentNames.size(); k++)
        out << "cbkg " << std::left << std::setw(40) << r.coherentNames[k] << std::right << std::setw(18)
            << r.coherentValues[k] << (r.coherentFixed[k] ? "  fixed" : "") << "\n";
    }
    // Spectator-momentum window (ps=...): the nodes.
    for (const ThmGroup &group : set_->thmGroups_) {
      if (group.name != r.name || !group.window) continue;
      const ThmSpectatorWindow &w = *group.window;
      out << "ps: " << w.description << "; mu_sx = " << w.muSx << " MeV\n"
          << "# The model at E is the average of the HOES cross section over the accepted spectator\n"
          << "# directions, weight |phi(p_s)|^2 d cos(theta_cm) (x the acceptance), incoherent; at fixed E\n"
          << "# the direction fixes p_s, and node k adds T_s = p_s^2/2mu_sx to E + B in the vertex.  The\n"
          << "# nodes at the lowest and the highest point; columns: E (MeV), theta_cm (deg), p_s (MeV/c),\n"
          << "# weight, T_s (MeV).\n";
      for (double e : {w.dataE.front(), w.dataE.back()}) {
        std::vector<ThmSpectatorWindow::Node> nodes;
        w.NodesAt(e, nodes);
        for (const ThmSpectatorWindow::Node &k : nodes)
          out << "ps_node" << std::setw(18) << e << std::setw(18) << k.theta << std::setw(18) << k.p
              << std::setw(18) << k.weight << std::setw(18) << k.es << "\n";
        out << std::left << std::setw(16) << "<T_s>" << std::right << std::setw(18) << e << std::setw(18)
            << w.MeanEs(e) << "\n";
      }
      // Memory of the per-point node tables (Unit::ThmPsTable): points
      // including the folding sub-points, the stored entrance channels summed
      // over the points (x nodes = entries), and their bytes.
      size_t tablePoints = 0, tableBytes = 0;
      for (int s : group.segments)
        for (int p = 1; p <= Data().GetSegment(s)->NumPoints(); p++) {
          Unit *point = Data().GetSegment(s)->GetPoint(*this, p);
          tablePoints++;
          tableBytes += point->ThmPsTableBytes();
          for (int q = 1; q <= point->NumSubPoints(); q++) {
            tablePoints++;
            tableBytes += point->GetSubPoint(*this, q)->ThmPsTableBytes();
          }
        }
      out << "# Node tables of the vertex: points (with folding sub-points), bytes.\n"
          << "ps_table" << std::setw(18) << (double)tablePoints << std::setw(18) << (double)tableBytes << "\n";
    }
    // Distortion factor (distortion=...): R at the lowest point, E_ref and the highest point.
    for (const ThmGroup &group : set_->thmGroups_) {
      if (group.name != r.name || !group.distortion) continue;
      const ThmDistortion &d = *group.distortion;
      out << "distortion: " << d.description << "\n";
      if (d.kind == ThmDistortion::TABLE) continue;
      double eLo, eHi;
      ThmGroupEnergyRange(*this, group, eLo, eHi);
      out << "# Zero-range DWBA transfer amplitude M(E) = <chi(-)_sF phi_sx chi(+)_aA(beta r)> (Mukhamedzhanov &\n"
          << "# Pang PRC 99 (2019) 064618 eqs. 20-24; Mukhamedzhanov arXiv:2609.04498 eqs. 22-30), M_PW its\n"
          << "# plane-wave limit; the model is multiplied by R(E) before folding (PWA-extracted S* / R).\n"
          << "# E_aA = " << d.eAA << " MeV, B = " << d.kin.bind << " MeV, k_aA = " << d.aa.k << " fm^-1, eta_aA = "
          << d.aa.eta << ", kappa = " << d.kappa << " fm^-1, eta_b = " << d.etaB << ", beta = " << d.beta << "\n"
          << (d.angWindow ? "# spectatorAngles: |M|^2 and |M_PW|^2 are averages over the accepted directions\n"
                            "# (weight d cos theta_cm x acceptance), theta_cm their weighted mean.\n"
                          : "")
          << "# Columns: E, E_sF (MeV), eta_sF, theta_cm (deg), |M|^2, |M_PW|^2, R, l_max.\n";
      for (double e : {eLo, d.eRef, eHi}) {
        ThmDistortion::Point p = d.Evaluate(e);
        out << "distortion_point" << std::setw(18) << e << std::setw(18) << p.esf << std::setw(18) << p.etasf
            << std::setw(18) << p.thetaCm << std::setw(18) << ThmDistortion::M2(p) << std::setw(18)
            << ThmDistortion::MPW2(p) << std::setw(18) << (p.ok ? d.R(p) : 0.0) << std::setw(6) << p.lmax << "\n";
      }
    }
    // Distorted-wave entrance vertex (vertexModel=dw): the Gram entries at the
    // lowest, the middle and the highest point.
    for (const ThmGroup &group : set_->thmGroups_) {
      if (group.name != r.name || !group.dwVertex) continue;
      const ThmDwVertex &v = *group.dwVertex;
      double eLo, eHi;
      ThmGroupEnergyRange(*this, group, eLo, eHi);
      out << "vertex: " << v.description << "\n"
          << "# Surface term of the prior-form DWBA (Mukhamedzhanov PRC 84 (2011) 044616; Mukhamedzhanov,\n"
          << "# Kadyrov & Pang EPJA 56 (2020) 233 eqs. 28-32) replaces M_l; R(E) is not applied.\n"
          << "# G = (4pi/(2l+1)) sum_m (s_m, d_m)^+ (s_m, d_m), (s_m, d_m) = (S_lm(a), a S_lm'(a))/(4pi phi~(q)),\n"
          << "# |M_l|^2 = c^+ G c with c = (B - 1, -1); plane waves: G11 = j_l(pa)^2, G22 = (pa j_l'(pa))^2.\n"
          << "# k_aA = " << v.dist.aa.k << " fm^-1, eta_aA = " << v.dist.aa.eta << ", kappa = " << v.dist.kappa
          << " fm^-1, alpha = " << v.alpha << ", beta = " << v.beta << ", a = " << v.radius << " fm\n"
          << "# Columns: E (MeV), q (MeV/c), p a, l, G11, G22, Re G12, Im G12 (spectatorAngle direction).\n";
      for (double e : {eLo, 0.5 * (eLo + eHi), eHi}) {
        std::vector<double> w, q, g, gd;
        double qd, pd;
        v.Interpolate(e, w, q, g, gd, qd, pd);
        for (size_t li = 0; li < v.lvals.size(); li++)
          out << "dw_vertex_point" << std::setw(18) << e << std::setw(18) << qd * hbarc << std::setw(18)
              << pd * v.radius << std::setw(4) << v.lvals[li] << std::setw(18) << gd[li * 4] << std::setw(18)
              << gd[li * 4 + 1] << std::setw(18) << gd[li * 4 + 2] << std::setw(18) << gd[li * 4 + 3] << "\n";
      }
    }
    // Coulomb line shape (lineshape=on): the ranges over the experiment's points.
    for (const ThmGroup &group : set_->thmGroups_) {
      if (group.name != r.name || !group.lineshape) continue;
      const ThmLineshape &ls = *group.lineshape;
      double eLo, eHi;
      ThmGroupEnergyRange(*this, group, eLo, eHi);
      out << "lineshape: on (spectator " << ls.spectator << ", Z_s = " << ls.Zs << ", Z_F = " << ls.ZF
          << "; E_aA = " << ls.eAA << " MeV, B = " << ls.bind << " MeV)\n"
          << "# |N_C|^2 = exp[2 zeta arctan(2 (E_lambda - E)/Gamma_lambda)] per level, zeta = eta_sB - eta_0\n"
          << "# (Mukhamedzhanov et al. EPJA 56 (2020) 233 eqs. 56-62, case 2: m_B >> m_s, m_b, eta_sb\n"
          << "# neglected); zeta < 0 moves the peaks up in E.  eta_sb is the neglected s-b term, averaged\n"
          << "# over the b direction: the approximation assumes |eta_sb| << 1.  Rows: the value at the\n"
          << "# lowest and at the highest point energy E (MeV); [k] is the exit pair key.\n";
      auto row = [&](const char *name, double a, double b) {
        out << std::left << std::setw(16) << name << std::right << std::setw(18) << a << std::setw(18) << b << "\n";
      };
      row("E", eLo, eHi);
      row("E_sF", ls.EsF(eLo), ls.EsF(eHi));
      row("eta_0", ls.Eta0(eLo), ls.Eta0(eHi));
      for (const ThmLineshape::Exit &e : ls.exits) {
        std::ostringstream key;
        key << "zeta[" << e.pairKey << "]";
        row(key.str().c_str(), ls.Zeta(eLo, e.ZB, e.mB), ls.Zeta(eHi, e.ZB, e.mB));
        key.str("");
        key << "eta_sb[" << e.pairKey << "]";
        row(key.str().c_str(), ls.EtaSbEstimate(eLo, eLo + e.q, e.Zb, e.mb, e.mB),
            ls.EtaSbEstimate(eHi, eHi + e.q, e.Zb, e.mb, e.mB));
      }
    }
  }
}

bool WorkingCopy::ThmLineshapeTable(const std::string &name, const std::vector<double> &energies, const CNuc *compound, const LevelState &state,
                              const Config &configure, ThmLineshapeReport &out, std::string &why) {
  const ThmGroup *group = nullptr;
  for (const ThmGroup &g : set_->thmGroups_)
    if (g.name == name) group = &g;
  if (!group) {
    why = "no THM experiment '" + name + "' in use";
    return false;
  }
  if (!group->lineshape) {
    why = "THM experiment '" + name + "' has no line shape (lineshape=on)";
    return false;
  }
  const ThmLineshape &ls = *group->lineshape;
  out = ThmLineshapeReport();
  out.experiment = name;
  out.spectator = ls.spectator;
  out.Zs = ls.Zs;
  out.ZF = ls.ZF;
  out.eAA = ls.eAA;
  out.bind = ls.bind;
  out.energy = energies;
  for (double e : energies) {
    out.esf.push_back(ls.EsF(e));
    out.eta0.push_back(ls.Eta0(e));
  }
  int entranceKey = Data().GetSegment(group->segments[0])->GetEntranceKey();
  int aa = compound->GetPairNumFromKey(entranceKey);
  const PPair *entrance = compound->GetPair(aa);
  double threshold = entrance->GetSepE() + entrance->GetExE();
  for (const ThmLineshape::Exit &x : ls.exits) {
    ThmLineshapeReport::Exit ex;
    ex.pairKey = x.pairKey;
    ex.Zb = x.Zb;
    ex.ZB = x.ZB;
    ex.mb = x.mb;
    ex.mB = x.mB;
    for (double e : energies) {
      ex.zeta.push_back(ls.Zeta(e, x.ZB, x.mB));
      ex.etaSb.push_back(ls.EtaSbEstimate(e, e + x.q, x.Zb, x.mb, x.mB));
    }
    int exitNum = compound->GetPairNumFromKey(x.pairKey);
    for (int j = 1; j <= compound->NumJGroups(); j++) {
      const JGroup *jg = compound->GetJGroup(j);
      if (!jg->IsInRMatrix()) continue;
      bool in = false, outCh = false;
      for (int ch = 1; ch <= jg->NumChannels(); ch++) {
        in = in || jg->GetChannel(ch)->GetPairNum() == aa;
        outCh = outCh || jg->GetChannel(ch)->GetPairNum() == exitNum;
      }
      if (!in || !outCh) continue;
      for (int la = 1; la <= jg->NumLevels(); la++) {
        const ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        ThmLineshapeReport::Level lv;
        lv.jgroup = j;
        lv.level = la;
        lv.J = jg->GetJ();
        lv.pi = jg->GetPi();
        lv.energy = state.E(j, la) - threshold;
        lv.width = ThmLevelWidth(compound, state, j, la, configure);
        for (size_t k = 0; k < energies.size(); k++)
          lv.nc2.push_back(ThmLineshapeFactorSq(ex.zeta[k], lv.energy - energies[k], lv.width));
        ex.levels.push_back(lv);
      }
    }
    out.exits.push_back(ex);
  }
  return true;
}

bool WorkingCopy::ThmDistortionTable(const std::string &name, const std::vector<double> &energies, ThmDistortionReport &out,
                               std::string &why) {
  const ThmGroup *group = nullptr;
  for (const ThmGroup &g : set_->thmGroups_)
    if (g.name == name) group = &g;
  if (!group) {
    why = "no THM experiment '" + name + "' in use";
    return false;
  }
  if (!group->distortion) {
    why = "THM experiment '" + name + "' has no distortion (distortion=coulomb|optical|table:<file>)";
    return false;
  }
  const ThmDistortion &d = *group->distortion;
  out = ThmDistortionReport();
  out.experiment = name;
  out.description = d.description;
  out.energy = energies;
  if (d.kind == ThmDistortion::TABLE) {
    out.kind = "table";
    for (double e : energies) out.rModel.push_back(d.Weight(e));
    return true;
  }
  out.kind = d.kind == ThmDistortion::COULOMB ? "coulomb" : "optical";
  out.eRef = d.eRef;
  out.eAA = d.eAA;
  out.bind = d.kin.bind;
  out.kAA = d.aa.k;
  out.etaAA = d.aa.eta;
  out.kappa = d.kappa;
  out.etaB = d.etaB;
  out.beta = d.beta;
  for (double e : energies) {
    ThmDistortion::Point p;
    std::string bad = d.CheckEnergy(e);
    if (bad.empty()) p = d.Evaluate(e);
    if (!bad.empty() || !p.ok) {
      why = "THM experiment '" + name + "': " + (bad.empty() ? p.why : bad);
      return false;
    }
    out.esf.push_back(p.esf);
    out.ksf.push_back(p.ksf);
    out.etasf.push_back(p.etasf);
    out.thetaCm.push_back(p.thetaCm);
    out.x.push_back(p.x);
    out.q.push_back(p.q);
    out.m2.push_back(ThmDistortion::M2(p));
    out.mpw2.push_back(ThmDistortion::MPW2(p));
    out.r.push_back(d.R(p));
    out.rModel.push_back(d.Weight(e));
    out.lmax.push_back(p.lmax);
  }
  return true;
}

bool WorkingCopy::ThmVertexTable(const std::string &name, const std::vector<double> &energies, const CNuc *compound, const LevelState &state,
                           const Config &configure, ThmVertexReport &out, std::string &why) {
  const ThmGroup *group = nullptr;
  for (const ThmGroup &g : set_->thmGroups_)
    if (g.name == name) group = &g;
  if (!group) {
    why = "no THM experiment '" + name + "' in use";
    return false;
  }
  if (group->pairKey == 0 || !compound->IsPairKey(group->pairKey)) {
    why = "THM experiment '" + name + "': its segments have different entrance pairs";
    return false;
  }
  const bool useGSL = !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC);
  const int aa = compound->GetPairNumFromKey(group->pairKey);
  const PPair *pair = compound->GetPair(aa);
  if (pair->GetPType() != 0) {
    why = "THM experiment '" + name + "': the entrance pair is not a particle pair";
    return false;
  }
  out = ThmVertexReport();
  out.experiment = name;
  out.pairKey = group->pairKey;
  out.bind = pair->GetBindingEnergy();
  out.radius = pair->GetChRad();
  const ThmDwVertex *dw = group->dwVertex.get();
  out.energy = energies;
  if (group->window) {
    const ThmSpectatorWindow &w = *group->window;
    out.window = w.description;
    out.muSx = w.muSx;
    out.windowObject = group->window;
    for (double e : energies) {
      std::vector<ThmSpectatorWindow::Node> nodes;
      bool moved = false;
      w.Nodes(e, nodes, &moved);
      out.reached.push_back(moved ? 0 : 1);
      std::vector<double> p, wt, es, th;
      for (const ThmSpectatorWindow::Node &k : nodes) {
        p.push_back(k.p);
        wt.push_back(k.weight);
        es.push_back(k.es);
        th.push_back(k.theta);
      }
      out.p.push_back(p);
      out.weight.push_back(wt);
      out.es.push_back(es);
      out.theta.push_back(th);
    }
  } else if (dw) {
    out.window = dw->angles ? dw->dist.AngleText() : "delta";
    if (dw->angles) out.muSx = dw->dist.muSx;
    for (double e : energies) out.reached.push_back(dw->Reached(e) ? 1 : 0);
  } else {
    out.window = "delta";
    for (size_t i = 0; i < energies.size(); i++) {
      out.p.push_back({0.0});
      out.weight.push_back({1.0});
      out.es.push_back({configure.thm.SpectatorEnergy(group->pairKey)});
      out.reached.push_back(1);
    }
  }
  const double mu = pair->GetRedMass() * uconv;
  // M_l = (B - 1) j_l - rho j_l' + C_l at E with T_s = es added to E + B (Unit::CalcEDependentValues).
  struct Pieces {
    double jl = 0.0, rhoDjl = 0.0;
    complex coul = complex(0.0, 0.0);
  };
  auto pieces = [&](int l, double e, double es) {
    Pieces q;
    double b = out.bind + es;
    if (e + b > 0.0) {
      ThmBesselParts(l, mu, e, b, out.radius, q.jl, q.rhoDjl);
      if (configure.thm.coulombIntegral && pair->GetZ(1) * pair->GetZ(2) != 0)
        q.coul = ThmCoulombTerm(pair, l, e, ThmRho(mu, e, b, 1.0), useGSL);
    }
    return q;
  };
  // vertexModel=dw: the Gram matrices of the DW vertex on the grid.
  std::vector<std::vector<double>> dwG, dwGd;
  if (dw) {
    out.model = "dw";
    for (double e : energies) {
      std::vector<double> w, q, g, gd;
      double qd, pd;
      dw->Interpolate(e, w, q, g, gd, qd, pd);
      if (dw->angles) {
        // The directions (at the nearest grid energy); M2 uses the averaged G (one node).
        std::vector<double> aw, aq, th;
        dw->AngleNodesAt(e, aw, aq, th);
        out.dwTheta.push_back(th);
        out.dwAngleQ.push_back(aq);
        out.dwAngleWeight.push_back(aw);
      }
      const double ks = std::sqrt(2.0 * dw->dist.sf.mu * std::max(dw->dist.EsF(e), 0.0)) / hbarc;
      const double ka = dw->dist.aa.k, kb = dw->beta * ka;
      std::vector<double> row;
      for (double qk : q) {
        double qq = qk / hbarc, x = ks > 0.0 ? (ks * ks + kb * kb - qq * qq) / (2.0 * ks * kb) : 1.0;
        x = std::max(-1.0, std::min(1.0, x));
        row.push_back(std::sqrt(std::max(0.0, ka * ka + dw->alpha * dw->alpha * ks * ks - 2.0 * ka * dw->alpha * ks * x)) *
                      out.radius);
      }
      out.rho.push_back(row);
      out.dwQ.push_back(q);
      out.dwWeight.push_back(w);
      out.dwQDelta.push_back(qd * hbarc);
      out.dwPDelta.push_back(pd);
      dwG.push_back(g);
      dwGd.push_back(gd);
    }
  } else
    for (double e : energies) {
      std::vector<double> row;
      for (double es : out.es[out.rho.size()])
        row.push_back(e + out.bind + es > 0.0 ? ThmRho(mu, e, out.bind + es, out.radius) : 0.0);
      out.rho.push_back(row);
    }
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    const JGroup *jg = compound->GetJGroup(j);
    if (!jg->IsInRMatrix()) continue;
    for (int ch = 1; ch <= jg->NumChannels(); ch++) {
      const AChannel *c = jg->GetChannel(ch);
      if (c->GetPairNum() != aa) continue;
      ThmVertexReport::Channel cr;
      cr.jgroup = j;
      cr.channel = ch;
      cr.J = jg->GetJ();
      cr.pi = jg->GetPi();
      cr.l = c->GetL();
      cr.s = c->GetS();
      // The vertex boundary of the model (THMMatrixFunc::CalculateTHMCrossSection).
      ThmVertexBoundary boundaryOf(configure, pair, jg, state, j, ch);
      // The pieces on the grid: [E][node] and the quasi-free [E].
      std::vector<std::vector<Pieces>> grid(energies.size());
      std::vector<Pieces> qf(energies.size());
      if (!dw)
        for (size_t i = 0; i < energies.size(); i++) {
          for (double es : out.es[i]) grid[i].push_back(pieces(cr.l, energies[i], es));
          qf[i] = pieces(cr.l, energies[i], 0.0);
        }
      const int li = dw ? dw->LIndex(cr.l) : -1;
      for (int la = 1; la <= jg->NumLevels(); la++) {
        const ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        ThmVertexReport::Level lr;
        lr.level = la;
        lr.boundary = boundaryOf.OnShell() ? std::nan("") : boundaryOf.Level(la);
        for (size_t i = 0; i < energies.size(); i++) {
          const complex B = boundaryOf.OnShell() ? boundaryOf.OnShellAt(energies[i]) : complex(lr.boundary, 0.0);
          if (dw) {
            const int nl = (int)dw->lvals.size();
            double avg = 0.0;
            for (size_t k = 0; k < out.dwWeight[i].size(); k++)
              avg += out.dwWeight[i][k] * (li < 0 ? 0.0 : ThmDwVertex::Vertex2Factored(&dwG[i][(k * nl + li) * 4], B));
            lr.m2.push_back(avg);
            lr.m2qf.push_back(li < 0 ? 0.0 : ThmDwVertex::Vertex2Factored(&dwGd[i][li * 4], B));
            // Plane waves at the same kinematics: M_l(p), analytic j_l'.
            double rho = out.dwPDelta[i] * out.radius;
            double jl = gsl_sf_bessel_jl(cr.l, rho);
            double djl = rho > 0.0 ? cr.l * jl / rho - gsl_sf_bessel_jl(cr.l + 1, rho) : (cr.l == 1 ? 1.0 / 3.0 : 0.0);
            lr.m2pw.push_back(std::norm((B - 1.0) * jl - rho * djl));
            continue;
          }
          auto m2 = [&](const Pieces &q) { return std::norm((B - 1.0) * q.jl - q.rhoDjl + q.coul); };
          double avg = 0.0;
          for (size_t k = 0; k < out.es[i].size(); k++) avg += out.weight[i][k] * m2(grid[i][k]);
          lr.m2.push_back(avg);
          lr.m2qf.push_back(m2(qf[i]));
        }
        cr.levels.push_back(lr);
      }
      out.channels.push_back(cr);
    }
  }
  return true;
}
