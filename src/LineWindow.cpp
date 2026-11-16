// The line-shape integration of a capture into a broad unbound final state (LineWindow.h).
#include "LineWindow.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <sstream>
#include <stdexcept>

#include <gsl/gsl_integration.h>

#include "AChannel.h"
#include "ALevel.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "CoulFunc.h"
#include "GenMatrixFunc.h"
#include "JGroup.h"
#include "Kinematics.h"
#include "LevelState.h"
#include "NFIntegral.h"
#include "PPair.h"
#include "PathwayTable.h"
#include "Results.h"
#include "RunScratch.h"
#include "Site.h"

/// The external-capture amplitude of one pathway with its parts by final-state exterior function
/// (UnitExternalCapture.cpp).
complex ExternalCapturePathwayParts(const CompoundSystem *compound, const Config &configure, int entrance,
                                    const azure2::EntrancePathways &P, const azure2::SpinGroupRow &group,
                                    const azure2::ExternalRow &pathway, double levelEnergy, double cmEnergy,
                                    complex parts[2]);

namespace azure2 {

namespace {

constexpr int kLineNodes = 48;

/// The final level of a site's capture: its group, level and its one open particle channel.
struct LineFinal {
  int aa = 0, ir = 0, j = 0, la = 0, ch = 0, exit = 0;
  const PPair *fp = nullptr;
  double thr = 0.0;
  bool ecLinked = true;  ///< false: internal capture alone, the final level found by J^pi and energy
};

/// The final level's energy at set-up (MeV, E_x): the energy external capture takes for it where it has the
/// external-capture link, else its energy in <levels> (the observed energy in the Brune parametrizations, the only
/// ones that serve a final level without the link).
double SetUpEnergy(const CompoundSystem *c, const LineFinal &lf) {
  const ALevel *lv = c->GetJGroup(lf.j)->GetLevel(lf.la);
  return lf.ecLinked ? lv->GetECEnergy() : lv->FileE();
}

std::string FindLineFinal(const CompoundSystem *c, int entranceKey, int exitKey, LineFinal &lf) {
  if (!c->IsPairKey(entranceKey) || !c->IsPairKey(exitKey)) return "its entrance or exit pair is not in the project";
  lf.aa = c->GetPairNumFromKey(entranceKey);
  lf.ir = c->GetPairNumFromKey(exitKey);
  const PPair *xp = c->GetPair(lf.ir);
  if (xp->GetPType() != 10) return "its exit pair (key " + std::to_string(exitKey) + ") is not a capture pair";
  const EntrancePathways &P = c->Pathways().Entrance(lf.aa);
  lf.exit = P.ExitOf(lf.ir);
  if (!lf.exit) return "its exit pair is not a decay of its entrance pair";
  for (int j = 1; j <= c->NumJGroups() && !lf.j; j++)
    for (int la = 1; la <= c->GetJGroup(j)->NumLevels(); la++)
      if (c->GetJGroup(j)->GetLevel(la)->IsECLevel() && c->GetJGroup(j)->GetLevel(la)->GetECPairNum() == lf.ir) {
        lf.j = j;
        lf.la = la;
        break;
      }
  if (!lf.j) {
    // Internal capture alone (no external capture into the pair, or external capture off for the run): the final
    // state is not linked to a level.  It is the level of the pair's J^pi whose energy in <levels> is the pair's
    // excitation energy (within 1 keV, the rule of the external-capture link, experiment::ExternalCaptureFinalStates).
    JGroup probe(xp->GetJ(2), xp->GetPi(2));
    probe.SetNucleus(xp->GetNucleus(), xp->GetNucleusZ(), xp->GetNucleusA());
    const int jn = c->IsJGroup(probe);
    int matches = 0;
    if (jn && c->GetJGroup(jn)->IsInRMatrix())
      for (int la = 1; la <= c->GetJGroup(jn)->NumLevels(); la++) {
        const ALevel *cand = c->GetJGroup(jn)->GetLevel(la);
        if (!cand->IsInRMatrix() || std::fabs(cand->FileE() - xp->GetHeavyExE()) > 1.0e-3) continue;
        matches++;
        lf.j = jn;
        lf.la = la;
      }
    if (matches != 1) {
      lf.j = lf.la = 0;
      std::ostringstream os;
      os << "its final state (J^pi = " << xp->GetJ(2) << (xp->GetPi(2) > 0 ? "+" : "-") << ", E_x = " << xp->GetHeavyExE()
         << " MeV) "
         << (matches ? "is ambiguous: more than one level of that J^pi lies"
                     : "is not a level of the project: the line shape needs that level, with its particle width,")
         << " within 1 keV of the pair's excitation energy";
      return os.str();
    }
    lf.ecLinked = false;
  }
  const JGroup *jg = c->GetJGroup(lf.j);
  const double levelEnergy = SetUpEnergy(c, lf);
  int open = 0;
  for (int ch = 1; ch <= jg->NumChannels(); ch++) {
    const PPair *p = c->GetPair(jg->GetChannel(ch)->GetPairNum());
    if (p->GetPType() != 0) continue;
    if (levelEnergy - p->GetSepE() - p->GetExE() > 0.0) {
      open++;
      lf.ch = ch;
      lf.fp = p;
      lf.thr = p->GetSepE() + p->GetExE();
    }
  }
  if (open == 0) return "its final state is bound (below every particle threshold): it has no line shape";
  if (open > 1) return "its final state lies above more than one particle threshold; the line-shape integration handles one";
  const ExitRow &exit = P.exits[lf.exit - 1];
  for (int k = 1; k <= exit.num_groups; k++) {
    const SpinGroupRow &g = P.Group(exit, k);
    for (int m = 1; m <= g.num_external; m++)
      if (P.External(g, m).final_channel != lf.ch)
        return "its external capture also runs through a closed channel of the final level, which the line-shape "
               "integration does not handle";
  }
  return "";
}

/// 1/N_f - 1 = sum over the particle channels of gamma~_c^2 dS_c/dE at the level energy, per unit gamma~_c^2: dS/dE
/// for an open channel, (2 mu a / hbar^2) int_a^inf W^2/W(a)^2 dr for a closed one (the same quantity,
/// ECFinalStateFactors of ExternalCaptureFinalStates.cpp).  For a final level whose N_f the state does not carry at
/// this energy: one without the external-capture link, or one whose line centre has moved.
std::vector<double> LineNFTerms(const CompoundSystem *c, const Config &configure, const LineFinal &lf, double levelEnergy) {
  const JGroup *jg = c->GetJGroup(lf.j);
  std::vector<double> terms((std::size_t)jg->NumChannels(), 0.0);
  for (int ch = 1; ch <= jg->NumChannels(); ch++) {
    const AChannel *chan = jg->GetChannel(ch);
    const PPair *p = c->GetPair(chan->GetPairNum());
    if (p->GetPType() != 0 || chan->GetRadType() != 'P') continue;
    const double e = levelEnergy - p->GetSepE() - p->GetExE();
    const int twoJ = p->SpinOrbitTwoJ(jg->GetJ());
    if (e > 0.0) {
      CoulFunc coul(p, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC));
      terms[ch - 1] = coul.PEShift_dE(chan->GetL(), twoJ, p->GetChRad(), e);
    } else {
      NFIntegral nf(p);
      terms[ch - 1] = 2.0 * p->GetRedMass() * uconv * p->GetChRad() / (hbarc * hbarc) * nf(chan->GetL(), twoJ, levelEnergy);
    }
  }
  return terms;
}

/// The excitation energy the line shape is centred on: the final level's observed energy from the CURRENT
/// parameters in the Brune parametrizations, where the level's working energy is the observed one, so that a fitted
/// final-level energy moves the line; the energy of the set-up in the formal parametrization.
double LineCentreEnergy(const CompoundSystem *c, const LevelState &state, const Config &configure, const LineFinal &lf) {
  const double setUp = SetUpEnergy(c, lf);
  if (!(configure.paramMask & Config::USE_BRUNE_FORMALISM)) return setUp;
  const double working = state.E(lf.j, lf.la);
  return working != 0.0 ? working : setUp;  // 0: the working values are not filled yet
}

/// The final level's observed reduced width amplitude in its open channel, gamma~ N_f^(1/2), from the current
/// parameters.  With the external-capture link the state carries N_f, at the energy of the set-up; without the link,
/// or once the line centre has moved away from that energy, 1/N_f = 1 + sum_c gamma~_c^2 dS_c/dE (`nfTerms`, at the
/// centre) is formed here.
double LineObservedAmplitude(const CompoundSystem *c, const LevelState &state, const Config &configure, const LineFinal &lf,
                             double centre, const std::vector<double> &nfTerms) {
  const bool brune = !!(configure.paramMask & Config::USE_BRUNE_FORMALISM);
  if (lf.ecLinked && (!brune || centre == SetUpEnergy(c, lf)))
    return state.ECGamma(lf.j, lf.la, lf.ch) * state.SqrtNF(lf.j, lf.la);
  auto amplitude = [&](int ch) { return lf.ecLinked ? state.ECGamma(lf.j, lf.la, ch) : state.Gamma(lf.j, lf.la, ch); };
  double sum = 1.0;
  for (std::size_t ch = 0; ch < nfTerms.size(); ch++)
    if (nfTerms[ch] != 0.0) sum += std::pow(amplitude((int)ch + 1), 2.0) * nfTerms[ch];
  return amplitude(lf.ch) / std::sqrt(sum);
}

std::uint64_t Bits(double x) {
  std::uint64_t b;
  std::memcpy(&b, &x, sizeof b);
  return b;
}

/// The nodes of `site` around the line centre `centre`, made (throws std::runtime_error where they cannot be).
std::shared_ptr<const LineNodeSet> MakeLineNodes(const Site &site, const CompoundSystem *c, const LevelState &state,
                                                 const Config &configure, const LineFinal &lf, double centre) {
  auto set = std::make_shared<LineNodeSet>();
  set->centre = centre;
  set->nfTerms = LineNFTerms(c, configure, lf, centre);
  const JGroup *jg = c->GetJGroup(lf.j);
  const AChannel *fc = jg->GetChannel(lf.ch);
  const int l = fc->GetL(), twoJ = lf.fp->SpinOrbitTwoJ(jg->GetJ());
  const double a = lf.fp->GetChRad();
  CoulFunc coul(lf.fp, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC));
  const double ER = centre - lf.thr;
  if (!(ER > 0.0)) throw std::runtime_error("line_window: the final level's energy has moved below its particle threshold.");
  const double gt = LineObservedAmplitude(c, state, configure, lf, centre, set->nfTerms);
  const double gam = 2.0 * coul.Penetrability(l, twoJ, a, ER) * gt * gt;
  if (!(gam > 0.0)) throw std::runtime_error("line_window: the final level has no width in its open channel.");
  const PPair *entrancePair = c->GetPair(lf.aa);
  const double inEnergy = site.GetCMEnergy() + entrancePair->GetSepE() + entrancePair->GetExE();
  const double W = site.GetLineWindow() * 1e-3;
  double lo = std::max(ER - W, 1e-3);
  const double hi = std::min(ER + W, inEnergy - lf.thr - 1e-3);
  // the outgoing-wave limit of the narrow resonance: no node within about a keV of a charged threshold
  while (lo < hi && Kinematics::Channel(lf.fp).Eta(lo) > 100.0) lo *= 1.05;
  if (!(hi > lo)) return set;
  const std::unique_ptr<gsl_integration_glfixed_table, void (*)(gsl_integration_glfixed_table *)> gl(
      gsl_integration_glfixed_table_alloc(kLineNodes), gsl_integration_glfixed_table_free);
  const double t0 = std::atan(2.0 * (lo - ER) / gam), t1 = std::atan(2.0 * (hi - ER) / gam);
  const double twoMuH2 = 2.0 * lf.fp->GetRedMass() * uconv / (hbarc * hbarc);
  const EntrancePathways &P = c->Pathways().Entrance(lf.aa);
  const ExitRow &exit = P.exits[lf.exit - 1];
  for (std::size_t i = 0; i < gl->n; i++) {
    double t, w;
    gsl_integration_glfixed_point(t0, t1, i, &t, &w, gl.get());
    LineNodeData n;
    n.ef = ER + 0.5 * gam * std::tan(t);
    n.weight = w * 0.5 * gam / (std::cos(t) * std::cos(t));
    n.sqrtDens = std::sqrt(twoMuH2 / (pi * std::sqrt(twoMuH2 * n.ef)));
    const CoulWaves cw = coul(l, twoJ, a, n.ef);
    n.Fa = cw.F;
    n.Ga = cw.G;
    n.P = coul.Penetrability(l, twoJ, a, n.ef);
    n.S = coul.PEShift(l, twoJ, a, n.ef);
    n.parts.resize((std::size_t)exit.num_groups);
    for (int k = 1; k <= exit.num_groups; k++) {
      const SpinGroupRow &g = P.Group(exit, k);
      n.parts[k - 1].resize((std::size_t)g.num_external);
      for (int m = 1; m <= g.num_external; m++) {
        complex pp[2];
        ExternalCapturePathwayParts(c, configure, lf.aa, P, g, P.External(g, m), lf.thr + n.ef, site.GetCMEnergy(), pp);
        n.parts[k - 1][m - 1] = std::make_pair(pp[0], pp[1]);
      }
    }
    set->nodes.push_back(std::move(n));
  }
  return set;
}

/// The nodes of `site` at `centre`: the run's, made and kept when the run has none (the newest two centres stay).
std::shared_ptr<const LineNodeSet> LineNodesOf(const Site &site, const CompoundSystem *c, const LevelState &state,
                                               const Config &configure, const LineFinal &lf, double centre) {
  LineNodeStore &store = configure.scratch->lineNodes;
  const std::vector<double> key = {(double)site.GetEntranceKey(), (double)site.GetExitKey(), site.GetCMEnergy(),
                                   site.GetLineWindow()};
  {
    std::lock_guard<std::mutex> lock(store.mutex);
    const auto found = store.sets.find(key);
    if (found != store.sets.end())
      for (const std::shared_ptr<const LineNodeSet> &s : found->second)
        if (Bits(s->centre) == Bits(centre)) return s;
  }
  std::shared_ptr<const LineNodeSet> made = MakeLineNodes(site, c, state, configure, lf, centre);  // outside the lock
  std::lock_guard<std::mutex> lock(store.mutex);
  std::vector<std::shared_ptr<const LineNodeSet>> &kept = store.sets[key];
  for (const std::shared_ptr<const LineNodeSet> &s : kept)
    if (Bits(s->centre) == Bits(centre)) return s;  // another thread made it meanwhile
  if (kept.size() >= 2) kept.erase(kept.begin());
  kept.push_back(made);
  return made;
}

}  // namespace

std::string LineWindowProblem(int entranceKey, int exitKey, bool otherObservable, const CompoundSystem *compound,
                              const Config &configure) {
  if (configure.paramMask & Config::USE_RMC_FORMALISM) return "the RMC formalism is not supported";
  if (configure.paramMask & Config::USE_ANALYTIC_GRADIENT)
    return "--use-gradient is not supported (its gradient does not run over the line shape); fit without it";
  if (otherObservable)
    return "only a capture cross section (angle-integrated or differential) can be integrated over the line shape";
  LineFinal lf;
  const std::string why = FindLineFinal(compound, entranceKey, exitKey, lf);
  // Without the external-capture link nothing keeps the observed energy and amplitude of a final state whose formal
  // parameters differ from them.
  if (lf.j && !lf.ecLinked && !(configure.paramMask & Config::USE_BRUNE_FORMALISM)) {
    // the boundary conditions are taken at the lowest level of a group: its formal parameters are its observed ones
    const JGroup *jg = compound->GetJGroup(lf.j);
    bool lowest = true;
    for (int la = 1; la <= jg->NumLevels(); la++)
      if (la != lf.la && jg->GetLevel(la)->IsInRMatrix() && jg->GetLevel(la)->FileE() < jg->GetLevel(lf.la)->FileE())
        lowest = false;
    if (!lowest)
      return "without external capture into its pair, the formal parametrization (--no-brune) serves only a final "
             "state that is the lowest level of its J^pi group (the one the boundary conditions are taken at); use the "
             "Brune parametrization, or switch external capture on for the pair";
  }
  return why;
}

std::string LineWindowProblem(const Site &site, const CompoundSystem *compound, const Config &configure) {
  return LineWindowProblem(site.GetEntranceKey(), site.GetExitKey(),
                           site.IsAngularDist() || site.IsPhase() || site.IsAnalyzingPower() || site.IsTHM(), compound,
                           configure);
}

void LineWindowCrossSection(GenMatrixFunc &func, Site *site, Results &results) {
  const CompoundSystem *c = func.compound();
  const LevelState &state = func.state();
  const Config &configure = func.configure();
  const std::string why = LineWindowProblem(*site, c, configure);
  if (!why.empty()) throw std::runtime_error("line_window: " + why + ".");
  LineFinal lf;
  FindLineFinal(c, site->GetEntranceKey(), site->GetExitKey(), lf);
  const double centre = LineCentreEnergy(c, state, configure, lf);
  const std::shared_ptr<const LineNodeSet> set = LineNodesOf(*site, c, state, configure, lf, centre);

  const std::size_t length = site->RecordLength();
  const std::size_t countAt = site->HasAngularPart() ? Results::AngularOffset() : length;  // not a value: kept
  std::vector<double> sum(length, 0.0);
  double count = 0.0;
  if (!set->nodes.empty()) {
    const JGroup *jg = c->GetJGroup(lf.j);
    const AChannel *fc = jg->GetChannel(lf.ch);
    const int l = fc->GetL(), twoJ = lf.fp->SpinOrbitTwoJ(jg->GetJ());
    const double a = lf.fp->GetChRad();
    CoulFunc coul(lf.fp, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC));
    // the current parameters of the final level
    const double ER = centre - lf.thr;
    const double gt = LineObservedAmplitude(c, state, configure, lf, centre, set->nfTerms);
    const double SR = coul.PEShift(l, twoJ, a, ER), dSR = coul.PEShift_dE(l, twoJ, a, ER);
    if (!(1.0 - gt * gt * dSR > 0.0) || gt == 0.0)
      throw std::runtime_error("line_window: the final level's reduced width has no formal (one-level) equivalent.");
    const double g2 = gt * gt / (1.0 - gt * gt * dSR);
    const double s = std::sqrt(2.0 * lf.fp->GetRedMass() * uconv * a) / hbarc;
    const PPair *entrancePair = c->GetPair(lf.aa);
    const double inEnergy = site->GetCMEnergy() + entrancePair->GetSepE() + entrancePair->GetExE();
    LineNode node;
    func.SetActiveLineNode(&node);
    try {
      for (const LineNodeData &n : set->nodes) {
        const double dR = std::atan2(n.P * g2, ER - n.ef - g2 * (n.S - SR));
        const double d = dR - std::atan2(n.Fa, n.Ga);
        const double cd = std::cos(d), sd = std::sin(d);
        node.photonEnergy = inEnergy - (lf.thr + n.ef);
        node.amp = n.sqrtDens * (cd * n.Fa + sd * n.Ga) / (gt * s);
        node.ec.resize(n.parts.size());
        for (std::size_t k = 0; k < n.parts.size(); k++) {
          node.ec[k].resize(n.parts[k].size());
          for (std::size_t m = 0; m < n.parts[k].size(); m++)
            node.ec[k][m] = n.sqrtDens * (cd * n.parts[k][m].first + sd * n.parts[k][m].second);
        }
        func.RecalculateTMatrix(site);
        func.CalculateCrossSection(site, results);
        const double *rec = site->Record(results);
        for (std::size_t i = 0; i < length; i++)
          if (i != countAt) sum[i] += n.weight * rec[i];
        if (countAt < length) count = std::max(count, rec[countAt]);
      }
    } catch (...) {
      func.SetActiveLineNode(nullptr);
      throw;
    }
    func.SetActiveLineNode(nullptr);
  }
  double *rec = site->Record(results);
  for (std::size_t i = 0; i < length; i++) rec[i] = (i == countAt) ? count : sum[i];
}

}  // namespace azure2
