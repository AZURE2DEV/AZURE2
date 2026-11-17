#include "ThmFormFactor.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <sstream>

#include "ALevel.h"
#include "AChannel.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "Constants.h"
#include "JGroup.h"
#include "PPair.h"
#include "ThmExperiment.h"
#include "ThmGroup.h"
#include "ThmLineshape.h"
#include "ThmVertexBoundary.h"
#include "WorkingCopy.h"


namespace ff = azure2::physics::transfer::form_factor;

const ThmFormFactorTable::Slot *ThmFormFactorTable::Find(int j, int la, int ch) const {
  for (const Slot &s : slots)
    if (s.j == j && s.level == la && s.channel == ch) return &s;
  return nullptr;
}

double ThmFormFactorTable::Delta(const Slot &s, double eps, double Ex, double B, bool *resolved) const {
  if (resolved) *resolved = false;
  if (!volume || !(eps > 0.0) || !s.f.ok) return 0.0;
  const double h2m = well.H2M();
  const double p = std::sqrt(eps / h2m);
  if (Ex != s.Ex || B != s.B) {
    if (resolved) *resolved = true;
    const ff::LevelFunction f = ff::SolveLevel(well, s.l, Ex - threshold, B, s.nodes);
    if (!f.ok) return 0.0;
    return (f.e - eps) / h2m * ff::Overlap(well, f, p) / f.ua;
  }
  const double g = p <= pMax ? s.g(p) : ff::Overlap(well, s.f, p) / s.f.ua;
  return (s.f.e - eps) / h2m * g;
}

double ThmFormFactorTable::Ratio(const Slot &s, double eps) const {
  if (!(eps > 0.0)) return 1.0;
  const double p = std::sqrt(eps / well.H2M());
  const double M = ff::SurfaceM(s.l, p * well.a, s.B);
  return M != 0.0 ? 1.0 + Delta(s, eps, s.Ex, s.B) / M : 1.0;
}

// ---------------------------------------------------------------------------------------------------------------------
// The set-up: the table of an experiment.

namespace {
const int kChebyshevTerms = 48;
}  // namespace

bool ThmBuildFormFactor(WorkingCopy &data, const Config &configure, const CompoundSystem *compound,
                        const ThmExperiment &x, ThmGroup &group, const std::string &where) {
  std::vector<ThmExperiment> one(1, x);
  const std::string refused = CheckThmFormFactorVertex(one, configure.thm.vertex == Config::ThmOptions::PER_LEVEL,
                                                       !!(configure.paramMask & Config::USE_BRUNE_FORMALISM));
  if (!refused.empty()) {
    configure.outStream << "ERROR: <thm> " << refused << "." << std::endl;
    return false;
  }
  if (group.pairKey == 0) {
    configure.outStream << where << "vertexModel=ff needs one entrance pair; the data sets have different ones."
                        << std::endl;
    return false;
  }
  const int aa = compound->GetPairNumFromKey(group.pairKey);
  const PPair *pair = compound->GetPair(aa);
  if (pair->GetPType() != 0) {
    configure.outStream << where << "vertexModel=ff: the entrance pair is not a pair of two particles." << std::endl;
    return false;
  }
  std::shared_ptr<ThmFormFactorTable> t = std::make_shared<ThmFormFactorTable>();
  t->experiment = x.name;
  t->r0 = x.ffR0;
  t->a0 = x.ffA0;
  t->volume = x.ffVolume;
  t->threshold = pair->GetSepE() + pair->GetExE();
  ff::Well &w = t->well;
  w.mu = pair->GetRedMass() * uconv;
  w.a = pair->GetPairRadius();
  const long A1 = std::lround(pair->GetM(1)), A2 = std::lround(pair->GetM(2));
  w.R = t->r0 * std::cbrt((double)std::max(A1, A2));
  w.diff = t->a0;
  w.rc = w.R;
  w.zze2 = pair->GetZ(1) * pair->GetZ(2) * hbarc * fstruc;
  // The range of the Chebyshev series in p: up to the highest E + B + T_s of the points, with 1 MeV to spare for the
  // folding sub-points (beyond it Delta is computed directly).
  double eLo = 0.0, eHi = 0.0;
  WorkingCopy::ThmGroupEnergyRange(data, group, eLo, eHi);
  double tMax = configure.thm.SpectatorEnergy(group.pairKey);
  if (group.window)
    for (double e : group.window->dataE) {
      std::vector<ThmSpectatorWindow::Node> nodes;
      group.window->NodesAt(e, nodes);
      for (const ThmSpectatorWindow::Node &k : nodes) tMax = std::max(tMax, k.es);
    }
  const double epsMax = std::max(eHi, 0.0) + pair->GetBindingEnergy() + tMax + 1.0;
  t->pMax = std::sqrt(epsMax / w.H2M());
  // A nucleon on a core: the Wildermuth rule is defined.
  const int nucleon = A1 == 1 ? 1 : A2 == 1 ? 2 : 0;
  int coreCount = -1;
  std::string nucleonName;
  if (nucleon) {
    const int core = 3 - nucleon;
    const int zx = pair->GetZ(nucleon), zc = pair->GetZ(core);
    const long ac = core == 1 ? A1 : A2;
    if (zx == 1) {
      coreCount = zc;
      nucleonName = "proton";
    } else if (zx == 0) {
      coreCount = (int)ac - zc;
      nucleonName = "neutron";
    }
  }
  std::vector<int> used(x.ffNodes.size(), 0);
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    const JGroup *jg = compound->GetJGroup(j);
    if (!jg->IsInRMatrix()) continue;
    const std::string jpi = ThmSpinText(jg->GetJ()) + (jg->GetPi() > 0 ? "+" : "-");
    for (int ch = 1; ch <= jg->NumChannels(); ch++) {
      const AChannel *c = jg->GetChannel(ch);
      if (c->GetPairNum() != aa) continue;
      for (int la = 1; la <= jg->NumLevels(); la++) {
        const ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        if (!level->EnergyFixed()) {
          configure.outStream << where << "vertexModel=ff: level " << la << " of J = " << jpi << " (E_x = "
                              << level->FileE() << " MeV) has a free energy.  The full form factor solves each "
                              << "level function once, at the level's energy (resonance energies are kept at their "
                              << "literature values); fix the energy of this level, or use vertexModel=pw."
                              << std::endl;
          return false;
        }
        ThmFormFactorTable::Slot s;
        s.j = j;
        s.level = la;
        s.channel = ch;
        s.l = c->GetL();
        const ThmFfNodeTerm *term = MatchThmFfNodes(x.ffNodes, jg->GetJ(), jg->GetPi(), la, s.l);
        if (term) {
          s.nodes = term->nodes;
          used[term - &x.ffNodes[0]] = 1;
        } else if (coreCount >= 0) {
          s.nodes = ff::WildermuthNodes(s.l, coreCount);
          s.wildermuth = true;
        } else {
          configure.outStream << where << "vertexModel=ff: no node count for level " << la << " of J = " << jpi
                              << " (entrance l = " << s.l << "); the Wildermuth rule is defined for a nucleon on a "
                              << "core only, so give ffNodes= (e.g. ffNodes=" << jpi << ":<n>)." << std::endl;
          return false;
        }
        s.Ex = level->FileE();
        s.B = ThmVertexBoundary::ShiftAtLevelEnergy(pair, s.l, pair->SpinOrbitTwoJ(jg->GetJ()), c->GetRadius(), s.Ex,
                                                    nullptr);
        s.f = ff::SolveLevel(w, s.l, s.Ex - t->threshold, s.B, s.nodes);
        if (!s.f.ok) {
          configure.outStream << where << "vertexModel=ff: level " << la << " of J = " << jpi << " (l = " << s.l
                              << ", n_r = " << s.nodes << ", E = " << s.Ex - t->threshold << " MeV, B = " << s.B
                              << "): " << s.f.why << "." << std::endl;
          return false;
        }
        // g(p) = Int u r j_l(p r) dr / u(a) at the Chebyshev nodes of [0, pMax], as a series.
        const std::vector<double> p = ff::Chebyshev::Nodes(0.0, t->pMax, kChebyshevTerms);
        std::vector<double> g(p.size());
        for (size_t k = 0; k < p.size(); k++) g[k] = ff::Overlap(w, s.f, p[k]) / s.f.ua;
        s.g.Fit(0.0, t->pMax, g);
        t->slots.push_back(s);
      }
    }
  }
  for (size_t k = 0; k < x.ffNodes.size(); k++)
    if (!used[k]) {
      std::vector<ThmFfNodeTerm> term(1, x.ffNodes[k]);
      configure.outStream << where << "ffNodes term '" << FormatThmFfNodes(term)
                          << "' matches no level of the R matrix with a channel of entrance pair " << group.pairKey
                          << "; correct or remove it." << std::endl;
      return false;
    }
  if (t->slots.empty()) {
    configure.outStream << where << "vertexModel=ff: the entrance pair has no channel in the R matrix." << std::endl;
    return false;
  }
  std::ostringstream l;
  l.precision(6);
  l << "  Transfer vertex: full form factor (vertexModel=ff; La Cognata et al., ApJ 708 (2010) 796, Eqs. 20-21): "
    << "M_l + Delta_lambda(E), the volume term of a Woods-Saxon level function"
    << (t->volume ? "" : " -- ffVolume=off: Delta = 0, the vertex of vertexModel=pw") << ".\n"
    << "  Well: r0 = " << t->r0 << " fm, a0 = " << t->a0 << " fm (R = " << w.R << " fm, the heavier nucleus), "
    << (w.zze2 != 0.0 ? "uniform-sphere Coulomb R_C = R" : "no Coulomb") << ", no spin-orbit; channel radius " << w.a
    << " fm; depth from a u'(a)/u(a) = B_lambda = S_c(E_lambda) at E_lambda.\n"
    << "  Level functions: J^pi, level, l, n_r, E_lambda (MeV above the entrance threshold), B_lambda, V0 (MeV), "
       "R = 1 + Delta/M_l at E = E_lambda:";
  const double spectator = configure.thm.SpectatorEnergy(group.pairKey);
  for (const ThmFormFactorTable::Slot &s : t->slots) {
    const JGroup *jg = compound->GetJGroup(s.j);
    const double eps = s.f.e + pair->GetBindingEnergy() + spectator;
    l << "\n    " << ThmSpinText(jg->GetJ()) << (jg->GetPi() > 0 ? "+" : "-") << "  " << s.level << "  l = " << s.l
      << "  n_r = " << s.nodes << (s.wildermuth ? " (Wildermuth, " + nucleonName + ")" : " (ffNodes)")
      << "  E = " << s.f.e << "  B = " << s.B << "  V0 = " << s.f.depth << "  R = " << t->Ratio(s, eps);
  }
  configure.outStream << l.str() << std::endl;
  group.formFactor = t;
  for (int seg : group.segments) data.Draft().GetSegment(seg)->SetThmFormFactor(data, t);
  return true;
}
