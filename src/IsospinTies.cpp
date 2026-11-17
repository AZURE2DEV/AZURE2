// Isospin-labelled levels across compound nuclei (IsospinTies.h): the records of <isospin> resolved against the
// compound system, the ties as functions of the parameter vector, their penalties and the report.
#include "IsospinTies.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <map>
#include <sstream>

#include "AZUREParams.h"
#include "AngCoeff.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "LevelState.h"
#include "Minuit2/MnUserParameters.h"
#include "azure2/experiment/ProjectReader.h"

namespace ex = azure2::experiment;

namespace {

/// Isospin and projection of a particle of charge z and mass m (the mass number is round(m)): t = |N - Z|/2,
/// t3 = (Z - N)/2.
void ParticleIsospin(int z, double m, double &t, double &t3) {
  const int a = (int)std::lround(m);
  const int n = a - z;
  t3 = 0.5 * (z - n);
  t = std::fabs(t3);
}

std::string Fmt(double x, int prec = 4) {
  std::ostringstream o;
  o << std::setprecision(prec) << x;
  return o.str();
}

bool HoldsID(const ALevel *level, int id) {
  const std::vector<int> &ids = level->IDs();
  return std::find(ids.begin(), ids.end(), id) != ids.end();
}

}  // namespace

const IsospinSet::Shift *IsospinSet::ShiftOf(int nucleus) const {
  for (const Shift &e : nuclei_)
    if (e.nucleus == nucleus) return &e;
  return nullptr;
}

int IsospinSet::ShiftParamOf(int nucleus) const {
  const Shift *e = ShiftOf(nucleus);
  return e ? e->param : -1;
}

double IsospinSet::ShiftValue(int nucleus, const std::vector<double> &values) const {
  const int param = ShiftParamOf(nucleus);
  return param >= 0 ? values[param - firstParam_] : 0.0;
}

double IsospinSet::At(const CompoundSystem &compound, const LevelState &state, int index) const {
  if (index >= firstParam_) return state.IsospinValues()[index - firstParam_];
  const LevelParameterLayout::Slot &slot = compound.ParameterLayout().Slots()[index];
  return slot.channel == 0 ? state.E(slot.group, slot.level) : state.Gamma(slot.group, slot.level, slot.channel);
}

std::vector<double> IsospinSet::FileValues() const {
  std::vector<double> values;
  for (const Shift &e : shifts_) values.push_back(e.spec.shift.value);
  for (const Mix &m : mixes_) values.push_back(m.spec.amplitude.value);
  return values;
}

std::vector<double> IsospinSet::ValuesIn(const vector_r &q) const {
  if (!active_ || (int)q.size() < firstParam_ + NumParams()) return {};
  return std::vector<double>(q.begin() + firstParam_, q.begin() + firstParam_ + NumParams());
}

std::shared_ptr<const IsospinSet> IsospinSet::Read(const CompoundSystem &c, const Config &configure) {
  std::shared_ptr<IsospinSet> set = std::make_shared<IsospinSet>();
  if (!configure.project || configure.project->built.isospin.Empty()) return set;
  const ex::IsospinSpec &spec = configure.project->built.isospin;
  std::ostream &out = configure.outStream;
  auto fail = [&](const std::string &why) {
    out << "ERROR: <isospin>: " << why << std::endl;
    return nullptr;
  };
  set->active_ = true;

  // Two parameter spaces carry the ties, both linear in the fitted parameters:
  //  - formal (--no-brune), Hale's analysis: formal amplitudes and eigenvalues at the fixed boundary values
  //    B_c = -l_c, which CompoundSystem::CalcBoundaryConditions sets on the levels' J^pi groups (Groups());
  //  - Brune (the default): the Brune (on-resonance) amplitudes and the observed energies.  No boundary value enters.
  //    For a level alone in its J^pi group the Brune amplitudes ARE the formal ones (they do not depend on B), so the
  //    amplitude ties are the same statement; with other levels of the same J^pi in the group (background poles) the
  //    two differ through level mixing.  The energy relation differs in any case: observed energies carry each
  //    nucleus's boundary (Thomas-Ehrman) shift, which is level dependent, so one shift per nucleus is a weaker
  //    description here than in the formal space (soft ties with a wide energy width leave the energies to the data).
  // Park's orthogonalized amplitudes are not Clebsch-Gordan projections of one isospin amplitude (--use-park sets the
  // Brune bit too, so it is tested first).
  if (configure.paramMask & Config::USE_PARK_FORMALISM)
    return fail("isospin levels are not available in Park's basis: run in the Brune basis (the default) or the formal "
                "one (--no-brune, Hale's B = -l).");
  if (configure.paramMask & Config::USE_RMC_FORMALISM) return fail("isospin levels are not available with the RMC formalism.");
  set->brune_ = (configure.paramMask & Config::USE_BRUNE_FORMALISM) != 0;
  const bool brune = set->brune_;

  // The nuclei of the reference and shift records, and their thresholds.
  const LevelParameterLayout &layout = c.ParameterLayout();
  set->firstParam_ = layout.Size();
  int next = set->firstParam_;
  int refNucleus = 0;
  for (const ex::IsospinNucleusEntry &entry : spec.nuclei) {
    Shift e;
    e.spec = entry;
    for (int k = 1; k <= c.NumNuclei(); k++)
      if (c.GetNucleus(k)->Z == entry.Z && c.GetNucleus(k)->A == entry.A) e.nucleus = k;
    if (!e.nucleus)
      return fail(entry.record + ": no active particle pair forms the nucleus Z=" + std::to_string(entry.Z) + " A=" +
                  std::to_string(entry.A) + ".");
    const int pn = c.GetPairNumFromKey(entry.pairKey);
    if (!pn || c.GetPair(pn)->GetNucleus() != e.nucleus || c.GetPair(pn)->GetPType() != 0)
      return fail(entry.record + ": pair " + std::to_string(entry.pairKey) + " is not an active particle pair of " +
                  c.NucleusLabel(e.nucleus) + ".");
    e.threshold = c.GetPair(pn)->GetSepE() + c.GetPair(pn)->GetExE();
    if (entry.reference) {
      refNucleus = e.nucleus;
      set->reference_ = e;
    } else {
      e.param = next++;
      set->shifts_.push_back(e);
    }
    set->nuclei_.push_back(e);
  }
  if (!refNucleus) return fail("exactly one isospin_reference record is needed.");
  for (const ex::IsospinMixEntry &entry : spec.mixes) {
    Mix m;
    m.spec = entry;
    m.param = next++;
    set->mixes_.push_back(m);
  }

  auto pairIsospin = [&](int pairNum, double &t1, double &t31, double &t2, double &t32) {
    const PPair *pp = c.GetPair(pairNum);
    double a1, a2;
    ParticleIsospin(pp->GetZ(1), pp->GetM(1), a1, t31);
    ParticleIsospin(pp->GetZ(2), pp->GetM(2), a2, t32);
    t1 = a1;
    t2 = a2;
    for (const ex::IsospinPairEntry &o : spec.pairs)
      if (o.pairKey == pp->GetPairKey()) {
        t1 = 0.5 * o.twoT1;
        t2 = 0.5 * o.twoT2;
      }
  };
  // Which state of its nucleus a particle is in: the rank of its excitation energy among the distinct excitation
  // energies (1 keV) of the particles in the project with the same Z, A, spin and parity.  The ground state (or the
  // lowest state listed) is rank 0.  Mirror partners are ranked in their own nucleus (different Z), so the n-th listed
  // state of a spin and parity pairs with the n-th; a nucleus that is its own mirror (12C for A = 13) pairs each state
  // with itself.  Without the rank, an excited particle with the ground state's spin and parity (12C*(7.65) 0+) fell
  // in the ground-state arrangement and was tied to the partner nucleus's ground-state channel.
  auto stateKey = [&](const PPair *pp, int k) {
    std::ostringstream o;
    o << pp->GetZ(k) << ":" << std::lround(pp->GetM(k)) << ":" << pp->GetJ(k) << ":" << pp->GetPi(k);
    return o.str();
  };
  auto particleEx = [&](const PPair *pp, int k) { return k == 1 ? pp->GetLightExE() : pp->GetHeavyExE(); };
  std::map<std::string, std::set<long>> statesOf;
  for (int p = 1; p <= c.NumPairs(); p++) {
    const PPair *pp = c.GetPair(p);
    if (pp->GetPType() != 0) continue;
    for (int k = 1; k <= 2; k++) statesOf[stateKey(pp, k)].insert(std::lround(1000.0 * particleEx(pp, k)));
  }
  auto stateRank = [&](const PPair *pp, int k) {
    const std::set<long> &s = statesOf[stateKey(pp, k)];
    return (int)std::distance(s.begin(), s.find(std::lround(1000.0 * particleEx(pp, k))));
  };
  // Pairs whose particles match in mass number, spin, parity, isospin and state rank form one arrangement.
  auto arrangement = [&](int pairNum) {
    const PPair *pp = c.GetPair(pairNum);
    double t1, t31, t2, t32;
    pairIsospin(pairNum, t1, t31, t2, t32);
    auto rank = [&](int k) {
      const int r = pp->GetPType() == 0 ? stateRank(pp, k) : 0;
      return r > 0 ? "#" + std::to_string(r + 1) : std::string();
    };
    std::ostringstream o;
    o << "A" << std::lround(pp->GetM(1)) << "+A" << std::lround(pp->GetM(2)) << " j" << pp->GetJ(1) << (pp->GetPi(1) < 0 ? "-" : "+")
      << rank(1) << "," << pp->GetJ(2) << (pp->GetPi(2) < 0 ? "-" : "+") << rank(2) << " t" << t1 << "," << t2;
    return o.str();
  };

  std::map<std::string, double> arrangementRadius;
  for (const ex::IsospinLevelEntry &entry : spec.levels) {
    Level L;
    L.spec = entry;
    L.t = 0.5 * entry.twoT;
    const std::string idText = "level ID " + std::to_string(entry.levelId);
    for (int j = 1; j <= c.NumJGroups(); j++)
      for (int la = 1; la <= c.GetJGroup(j)->NumLevels(); la++) {
        const ALevel *lev = c.GetJGroup(j)->GetLevel(la);
        if (!HoldsID(lev, entry.levelId) || !lev->IsInRMatrix()) continue;
        Member m;
        m.nucleus = c.GetJGroup(j)->GetNucleus();
        m.jGroup = j;
        m.level = la;
        m.energySlot = layout.IndexOf(j, la, 0);
        for (const Member &o : L.members)
          if (o.nucleus == m.nucleus)
            return fail(idText + " has lines in two levels of " + c.NucleusLabel(m.nucleus) +
                        ": the lines of one level share J, parity and (within 1 keV) the energy in a nucleus.");
        const int z = c.GetJGroup(j)->GetNucleusZ(), a = c.GetJGroup(j)->GetNucleusA();
        m.t3 = 0.5 * (z - (a - z));
        if (std::fabs(m.t3) > L.t + 1e-9)
          return fail(idText + " has T=" + Fmt(L.t) + " but lines in " + c.NucleusLabel(m.nucleus) + ", whose T3=" + Fmt(m.t3) + ".");
        L.members.push_back(m);
      }
    if (L.members.empty()) return fail(idText + " has no active level lines in the R-matrix.");
    const JGroup *g0 = c.GetJGroup(L.members[0].jGroup);
    for (const Member &m : L.members) {
      const JGroup *g = c.GetJGroup(m.jGroup);
      if (g->GetJ() != g0->GetJ() || g->GetPi() != g0->GetPi())
        return fail(idText + ": its lines have different J or parity in different nuclei.");
    }
    L.ref = 0;
    for (size_t i = 0; i < L.members.size(); i++)
      if (L.members[i].nucleus == refNucleus) L.ref = (int)i;
    if (L.members.size() > 1)
      for (Member &m : L.members) {
        const Shift *e = set->ShiftOf(m.nucleus);
        if (!e)
          return fail(idText + " has lines in " + c.NucleusLabel(m.nucleus) +
                      ", which has no isospin_reference or isospin_shift record (its energies need a threshold and a shift).");
        m.threshold = e->threshold;
      }

    // Amplitudes: Clebsch-Gordan factors, isospin-channel keys.
    for (Member &m : L.members) {
      const JGroup *g = c.GetJGroup(m.jGroup);
      set->groups_.insert(m.jGroup);
      for (int ch = 1; ch <= g->NumChannels(); ch++) {
        const AChannel *chan = g->GetChannel(ch);
        // Photon and beta channels stay outside the ties: Hale's relations are between the particle reduced-width
        // amplitudes of isospin partners, and an electromagnetic or weak amplitude is not charge independent.  Their
        // slots keep their own parameters in every member (the group's B = -l applies to its particle channels).
        if (chan->GetRadType() != 'P' || c.GetPair(chan->GetPairNum())->GetPType() != 0) {
          m.untied++;
          continue;
        }
        Width w;
        w.channel = ch;
        w.slot = layout.IndexOf(m.jGroup, m.level, ch);
        const int pairKey = c.GetPair(chan->GetPairNum())->GetPairKey();
        double t1, t31, t2, t32;
        pairIsospin(chan->GetPairNum(), t1, t31, t2, t32);
        if (std::fabs(t31 + t32 - m.t3) > 1e-9)
          return fail("pair " + std::to_string(pairKey) + ": the particle isospin projections do not add up to the nucleus's T3.");
        w.cgComputed = (std::fabs(t31) <= t1 + 1e-9 && std::fabs(t32) <= t2 + 1e-9 && L.t <= t1 + t2 + 1e-9 &&
                        L.t >= std::fabs(t1 - t2) - 1e-9)
                           ? AngCoeff::ClebGord(t1, t2, L.t, t31, t32, m.t3)
                           : 0.0;
        if (std::fabs(w.cgComputed) < 1e-12) w.cgComputed = 0.0;
        w.cg = w.cgComputed;
        for (const ex::IsospinFactorEntry &f : spec.factors)
          if (f.levelId == entry.levelId && f.pairKey == pairKey) {
            w.cg = f.value;
            w.overridden = true;
          }
        const std::string arr = arrangement(chan->GetPairNum());
        std::ostringstream key;
        key << arr << " l=" << chan->GetL() << " s=" << chan->GetS();
        w.key = key.str();
        // The channel's own radius (a channel_radius record may set it per l, s, J): partner channels -- one
        // arrangement, the same l, s and J in every nucleus -- must share it.
        const double radius = chan->GetRadius();
        std::ostringstream radKey;
        radKey << w.key << " J=" << g->GetJ();
        auto rad = arrangementRadius.find(radKey.str());
        if (rad == arrangementRadius.end())
          arrangementRadius[radKey.str()] = radius;
        else if (std::fabs(rad->second - radius) > 1e-9)
          return fail("channels " + radKey.str() + " of one arrangement have channel radii " + Fmt(rad->second, 6) +
                      " and " + Fmt(radius, 6) + " fm; Hale's isospin relations need equal radii.");
        m.widths.push_back(w);
      }
    }
    // Owners: the reference member first, then the others in J-group order.
    std::vector<int> order;
    order.push_back(L.ref);
    for (int i = 0; i < (int)L.members.size(); i++)
      if (i != L.ref) order.push_back(i);
    std::map<std::string, std::pair<int, double>> owner;  // key -> (slot, cg)
    for (int i : order)
      for (Width &w : L.members[i].widths)
        if (w.cg != 0.0 && !owner.count(w.key)) {
          owner[w.key] = {w.slot, w.cg};
          w.ownerSlot = -1;
        }
    for (int i : order)
      for (Width &w : L.members[i].widths) {
        if (w.cg == 0.0) {
          w.ownerSlot = -2;
          w.ratio = 0.0;
          continue;
        }
        const auto &o = owner[w.key];
        if (o.first == w.slot) continue;
        w.ownerSlot = o.first;
        w.ratio = w.cg / o.second;
      }
    set->levels_.push_back(L);
  }

  // Mix and factor records name existing channels of declared levels.
  for (size_t i = 0; i < set->mixes_.size(); i++) {
    Mix &mx = set->mixes_[i];
    for (Level &L : set->levels_) {
      if (L.spec.levelId != mx.spec.levelId) continue;
      for (Member &m : L.members)
        for (Width &w : m.widths) {
          const AChannel *chan = c.GetJGroup(m.jGroup)->GetChannel(w.channel);
          if (c.GetPair(chan->GetPairNum())->GetPairKey() != mx.spec.pairKey || chan->GetL() != mx.spec.l ||
              std::fabs(chan->GetS() - 0.5 * mx.spec.twoS) > 1e-9)
            continue;
          if (w.ownerSlot != -2)
            return fail(mx.spec.record + ": the channel of level ID " + std::to_string(mx.spec.levelId) + ", pair " +
                        std::to_string(mx.spec.pairKey) + " is isospin-allowed (factor " + Fmt(w.cg) +
                        "); an isospin_mix is for a forbidden channel.");
          w.mix = (int)i;
          mx.slot = w.slot;
        }
    }
    if (mx.slot < 0)
      return fail(mx.spec.record + ": level ID " + std::to_string(mx.spec.levelId) + " has no channel with pair " +
                  std::to_string(mx.spec.pairKey) + ", l=" + std::to_string(mx.spec.l) + ", s=" + Fmt(0.5 * mx.spec.twoS) + ".");
  }
  for (const ex::IsospinFactorEntry &f : spec.factors) {
    bool used = false;
    for (const Level &L : set->levels_)
      for (const Member &m : L.members)
        for (const Width &w : m.widths)
          if (L.spec.levelId == f.levelId && w.overridden &&
              c.GetPair(c.GetJGroup(m.jGroup)->GetChannel(w.channel)->GetPairNum())->GetPairKey() == f.pairKey)
            used = true;
    if (!used)
      return fail(f.record + ": level ID " + std::to_string(f.levelId) + " has no channel of pair " + std::to_string(f.pairKey) + ".");
  }

  // Derived slots (fixed, overwritten) with their ties, or the penalty terms of a soft level.
  for (const Level &L : set->levels_) {
    const Member &r = L.members[L.ref];
    for (size_t i = 0; i < L.members.size(); i++) {
      const Member &m = L.members[i];
      if ((int)i != L.ref && L.members.size() > 1) {
        std::vector<std::pair<int, double>> sources = {{r.energySlot, 1.0}};
        if (set->ShiftParamOf(m.nucleus) >= 0) sources.push_back({set->ShiftParamOf(m.nucleus), 1.0});
        if (set->ShiftParamOf(r.nucleus) >= 0) sources.push_back({set->ShiftParamOf(r.nucleus), -1.0});
        if (!L.spec.soft) {
          set->derived_.insert(m.energySlot);
          set->ties_.push_back({m.energySlot, sources});
        } else {
          set->softTerms_.push_back({m.energySlot, sources, m.threshold - r.threshold, L.spec.sigmaE});
        }
      }
      for (const Width &w : m.widths) {
        if (w.ownerSlot == -1) continue;
        std::vector<std::pair<int, double>> sources;
        if (w.ownerSlot >= 0) sources.push_back({w.ownerSlot, w.ratio});
        if (!L.spec.soft) {
          if (w.ownerSlot == -2 && w.mix >= 0) sources.push_back({set->mixes_[w.mix].param, 1.0});
          set->derived_.insert(w.slot);
          set->ties_.push_back({w.slot, sources});
        } else {
          set->softTerms_.push_back({w.slot, sources, 0.0, L.spec.sigmaGamma});
        }
      }
    }
  }

  // Report.
  out << "Isospin levels: " << set->levels_.size() << " level(s); reference " << c.NucleusLabel(refNucleus)
      << (brune ? "; Brune basis: the ties act on the observed energies and the Brune (on-resonance) amplitudes."
                : "; B_c = -l_c in their J-groups.")
      << std::endl;
  for (const Level &L : set->levels_) {
    const JGroup *g = c.GetJGroup(L.members[L.ref].jGroup);
    out << "  level ID " << L.spec.levelId << "  T=" << L.t << "  J=" << g->GetJ() << (g->GetPi() < 0 ? "-" : "+")
        << (L.spec.soft ? "  SOFT (sigma E " + Fmt(L.spec.sigmaE) + " MeV, sigma gamma " + Fmt(L.spec.sigmaGamma) + " MeV^1/2)" : "")
        << "  in";
    for (size_t i = 0; i < L.members.size(); i++)
      out << " " << c.NucleusLabel(L.members[i].nucleus) << ((int)i == L.ref ? " (reference)" : "");
    out << std::endl;
    for (const Member &m : L.members)
      for (const Width &w : m.widths) {
        const AChannel *chan = c.GetJGroup(m.jGroup)->GetChannel(w.channel);
        out << "      pair " << std::setw(2) << c.GetPair(chan->GetPairNum())->GetPairKey() << "  l=" << chan->GetL()
            << " s=" << chan->GetS() << "  factor " << std::setw(9) << Fmt(w.cg, 6);
        if (w.overridden) out << " (computed " << Fmt(w.cgComputed, 6) << ")";
        if (w.ownerSlot == -1)
          out << "  free";
        else if (w.ownerSlot == -2)
          out << (w.mix >= 0 ? "  forbidden, mix " + std::to_string(w.mix + 1) : "  forbidden (0)");
        else
          out << "  = " << Fmt(w.ratio, 6) << " x free amplitude";
        out << std::endl;
      }
    for (const Member &m : L.members)
      if (m.untied > 0)
        out << "      " << c.NucleusLabel(m.nucleus) << ": " << m.untied
            << " photon/beta channel(s), not tied (own parameters in each nucleus)" << std::endl;
  }
  // The groups of labelled levels may hold unlabelled levels: say so.
  for (int j : set->groups_)
    for (int la = 1; la <= c.GetJGroup(j)->NumLevels(); la++) {
      const std::vector<int> &ids = c.GetJGroup(j)->GetLevel(la)->IDs();
      bool declared = false;
      for (const int id : ids) declared = declared || spec.Declares(id);
      if (!declared)
        out << "  note: level ID " << (ids.empty() ? 0 : ids[0]) << " in a J-group of " << c.NucleusLabel(c.GetJGroup(j)->GetNucleus())
            << " has no isospin label; it is unconstrained" << (brune ? "." : " (and uses B = -l like its group).") << std::endl;
    }
  // A nucleus the level's isospin allows, with pairs of its arrangements, but no line.
  for (const Level &L : set->levels_) {
    std::set<std::string> arrs;
    std::set<int> present;
    for (const Member &m : L.members) {
      present.insert(m.nucleus);
      for (const Width &w : m.widths) arrs.insert(arrangement(c.GetJGroup(m.jGroup)->GetChannel(w.channel)->GetPairNum()));
    }
    for (int n = 1; n <= c.NumNuclei(); n++) {
      if (present.count(n)) continue;
      const double t3 = 0.5 * (c.GetNucleus(n)->Z - (c.GetNucleus(n)->A - c.GetNucleus(n)->Z));
      if (std::fabs(t3) > L.t + 1e-9) continue;
      for (int p = 1; p <= c.NumPairs(); p++)
        if (c.GetPair(p)->GetNucleus() == n && c.GetPair(p)->GetPType() == 0 && arrs.count(arrangement(p))) {
          out << "  WARNING: level ID " << L.spec.levelId << " (T=" << L.t << ") has no lines in " << c.NucleusLabel(n)
              << ", which its isospin allows; it is absent there." << std::endl;
          break;
        }
    }
  }
  return set;
}

ParameterRecord IsospinSet::RecordOf(int k) const {
  if (k < (int)shifts_.size()) {
    const Shift &e = shifts_[k];
    return ParameterRecord::IsospinShift(e.spec.record, e.nucleus, e.spec.Z, e.spec.A);
  }
  const Mix &m = mixes_[k - (int)shifts_.size()];
  int j = 0, la = 0, ch = 0;
  for (const Level &L : levels_)
    for (const Member &mem : L.members)
      for (const Width &w : mem.widths)
        if (w.slot == m.slot) {
          j = mem.jGroup;
          la = mem.level;
          ch = w.channel;
        }
  return ParameterRecord::IsospinMix(m.spec.record, j, la, ch);
}

void IsospinSet::AddParameters(AZUREParams &params) const {
  if (!active_) return;
  ROOT::Minuit2::MnUserParameters &p = params.GetMinuitParams();
  for (size_t i = 0; i < shifts_.size(); i++) {
    const Shift &e = shifts_[i];
    const std::string name = e.spec.record + ".shift";  // format 3: <record>.<field>
    params.Add(name, e.spec.shift.value, e.spec.shift.hasPrior ? e.spec.shift.priorWidth : 0.1, RecordOf((int)i));
    if (!e.spec.shift.fit) p.Fix(name);
  }
  for (size_t i = 0; i < mixes_.size(); i++) {
    const Mix &m = mixes_[i];
    const std::string name = m.spec.record + ".amplitude";
    params.Add(name, m.spec.amplitude.value,
               m.spec.amplitude.hasPrior ? m.spec.amplitude.priorWidth : 0.1 * std::fabs(m.spec.amplitude.value) + 1e-3,
               RecordOf((int)(shifts_.size() + i)));
    if (!m.spec.amplitude.fit) p.Fix(name);
  }
  this->Sync(p);
  // Soft levels: the derived slots start on their tie targets, free (a width whose input is 0 is otherwise fixed).
  const vector_r q = p.Params();
  for (const SoftTerm &t : softTerms_) {
    if (t.slot >= (int)q.size()) continue;
    double target = t.constant;
    for (const auto &s : t.sources) target += s.second * q[s.first];
    p.SetValue(t.slot, target);
    if (p.Parameter(t.slot).IsFixed()) p.Release(t.slot);
    p.SetError(t.slot, std::max(0.1 * std::fabs(target), t.sigma));
  }
}

void IsospinSet::Apply(vector_r &q, const std::vector<double> &values) const {
  if (!active_ || (int)values.size() != NumParams() || (int)q.size() < firstParam_) return;
  for (const Level &L : levels_) {
    if (L.spec.soft) continue;  // free slots, penalised instead
    const Member &r = L.members[L.ref];
    if (L.members.size() > 1) {
      const double eps = q[r.energySlot] - r.threshold - ShiftValue(r.nucleus, values);
      for (size_t i = 0; i < L.members.size(); i++) {
        if ((int)i == L.ref) continue;
        const Member &m = L.members[i];
        q[m.energySlot] = eps + m.threshold + ShiftValue(m.nucleus, values);
      }
    }
    for (const Member &m : L.members)
      for (const Width &w : m.widths) {
        if (w.ownerSlot >= 0)
          q[w.slot] = w.ratio * q[w.ownerSlot];
        else if (w.ownerSlot == -2)
          q[w.slot] = (w.mix >= 0) ? values[mixes_[w.mix].param - firstParam_] : 0.0;
      }
  }
}

void IsospinSet::Sync(ROOT::Minuit2::MnUserParameters &p) const {
  if (!active_) return;
  vector_r q = p.Params();
  if ((int)q.size() < firstParam_ + NumParams()) return;
  Apply(q);
  for (int i : derived_) {
    if (i >= (int)q.size()) continue;
    if (!p.Parameter(i).IsFixed()) p.Fix(i);
    p.SetValue(i, q[i]);
  }
}

std::vector<IsospinSet::PenaltyRow> IsospinSet::PenaltyRows(const CompoundSystem &compound, const LevelState &state) const {
  std::vector<PenaltyRow> rows;
  if (!active_ || (int)state.IsospinValues().size() != NumParams()) return rows;
  ForEachPrior([&](int index, double centre, double width) {
    rows.push_back({(At(compound, state, index) - centre) / width, {{index, 1.0 / width}}});
  });
  for (const SoftTerm &t : softTerms_) {
    double target = t.constant;
    for (const auto &s : t.sources) target += s.second * At(compound, state, s.first);
    PenaltyRow row;
    row.r = (At(compound, state, t.slot) - target) / t.sigma;
    row.d.push_back({t.slot, 1.0 / t.sigma});
    for (const auto &s : t.sources) row.d.push_back({s.first, -s.second / t.sigma});
    rows.push_back(row);
  }
  return rows;
}

bool IsospinSet::HasPenalty() const {
  if (!active_) return false;
  bool any = !softTerms_.empty();
  ForEachPrior([&](int, double, double) { any = true; });
  return any;
}

double IsospinSet::Penalty(const CompoundSystem &compound, const LevelState &state) const {
  double chi2 = 0.0;
  for (const PenaltyRow &row : PenaltyRows(compound, state)) chi2 += row.r * row.r;
  return chi2;
}

double IsospinSet::SoftPenalty(const CompoundSystem &compound, const LevelState &state) const {
  double chi2 = 0.0;
  if (!active_ || (int)state.IsospinValues().size() != NumParams()) return chi2;
  for (const SoftTerm &t : softTerms_) {
    double target = t.constant;
    for (const auto &s : t.sources) target += s.second * At(compound, state, s.first);
    chi2 += std::pow((At(compound, state, t.slot) - target) / t.sigma, 2.0);
  }
  return chi2;
}

void IsospinSet::AddPenaltyGradient(const CompoundSystem &compound, const LevelState &state,
                                    std::vector<std::pair<int, double>> &grad) const {
  for (const PenaltyRow &row : PenaltyRows(compound, state))
    for (const auto &d : row.d) grad.push_back({d.first, 2.0 * row.r * d.second});
}

std::string IsospinSet::ParamLabel(int index) const {
  if (!active_) return std::string();
  for (const Shift &e : shifts_)
    if (e.param == index)
      return "internal Coulomb shift of isospin levels in Z=" + std::to_string(e.spec.Z) + " A=" + std::to_string(e.spec.A);
  for (size_t k = 0; k < mixes_.size(); k++)
    if (mixes_[k].param == index)
      return "isospin-forbidden amplitude " + std::to_string(k + 1) + ": level ID " + std::to_string(mixes_[k].spec.levelId) +
             ", pair " + std::to_string(mixes_[k].spec.pairKey) + ", l=" + std::to_string(mixes_[k].spec.l) + ", s=" +
             Fmt(0.5 * mixes_[k].spec.twoS);
  return std::string();
}

std::string IsospinSet::MemberLabel(int jGroup, int level) const {
  if (!active_) return std::string();
  for (const Level &L : levels_)
    for (size_t i = 0; i < L.members.size(); i++)
      if (L.members[i].jGroup == jGroup && L.members[i].level == level)
        return "ISOSPIN LEVEL " + std::to_string(L.spec.levelId) + " T=" + Fmt(L.t) +
               ((int)i == L.ref ? " (reference)" : (L.spec.soft ? " (soft)" : " (derived)"));
  return std::string();
}

void IsospinSet::Write(const Config &configure, const CompoundSystem &c, const LevelState &state) const {
  if (!active_) return;
  std::ofstream out(configure.OpenPath(configure.outputdir + "isospin.out").c_str());
  if (!out) return;
  out << (brune_ ? "ISOSPIN LEVELS (Brune parameters: observed energies and on-resonance amplitudes, no boundary values)"
                 : "ISOSPIN LEVELS (formal parameters, B_c = -l_c in their J-groups)")
      << std::endl
      << std::endl;
  const std::vector<double> &values = state.IsospinValues();
  const bool haveVector = (int)values.size() == NumParams();
  if (!softTerms_.empty() && haveVector)
    out << "soft-tie penalty = " << std::setprecision(8) << SoftPenalty(c, state) << std::endl;
  for (const Shift &e : nuclei_) {
    out << (e.spec.reference ? "reference " : "shift     ") << "Z=" << e.spec.Z << " A=" << e.spec.A << "  threshold pair "
        << e.spec.pairKey << " at Ex = " << std::setprecision(6) << e.threshold << " MeV";
    if (!e.spec.reference)
      out << "  shift = " << std::setprecision(6) << (haveVector ? values[e.param - firstParam_] : e.spec.shift.value) << " MeV"
          << (e.spec.shift.fit ? " (varied)" : " (fixed)");
    out << std::endl;
  }
  for (size_t k = 0; k < mixes_.size(); k++)
    out << "mix " << k + 1 << "  level ID " << mixes_[k].spec.levelId << "  pair " << mixes_[k].spec.pairKey << "  l="
        << mixes_[k].spec.l << " s=" << 0.5 * mixes_[k].spec.twoS << "  amplitude = " << std::setprecision(6)
        << (haveVector ? values[mixes_[k].param - firstParam_] : mixes_[k].spec.amplitude.value) << " MeV^1/2"
        << (mixes_[k].spec.amplitude.fit ? " (varied)" : " (fixed)") << std::endl;
  out << std::endl;
  for (const Level &L : levels_) {
    const JGroup *g = c.GetJGroup(L.members[L.ref].jGroup);
    out << "level ID " << L.spec.levelId << "  T=" << L.t << "  J=" << g->GetJ() << (g->GetPi() < 0 ? "-" : "+") << std::endl;
    for (size_t i = 0; i < L.members.size(); i++) {
      const Member &m = L.members[i];
      out << "  " << c.NucleusLabel(m.nucleus) << " (T3=" << m.t3 << ")" << ((int)i == L.ref ? "  reference" : "  derived")
          << (brune_ ? "  observed E = " : "  formal E = ") << std::setprecision(14) << state.E(m.jGroup, m.level) << " MeV"
          << std::endl;
      if (L.spec.soft && (int)i != L.ref) out << "    (soft tie)" << std::endl;
      for (const Width &w : m.widths) {
        const AChannel *chan = c.GetJGroup(m.jGroup)->GetChannel(w.channel);
        out << "    pair " << std::setw(2) << c.GetPair(chan->GetPairNum())->GetPairKey() << "  l=" << chan->GetL()
            << " s=" << chan->GetS();
        if (!brune_) out << "  B=" << chan->GetBoundaryCondition();
        out << "  factor " << std::setw(10) << std::setprecision(6) << w.cg;
        if (w.overridden) out << " (computed " << w.cgComputed << ")";
        out << "  gamma = " << std::setw(22) << std::setprecision(14) << state.Gamma(m.jGroup, m.level, w.channel) << " MeV^1/2";
        if (w.ownerSlot == -1)
          out << "  free";
        else if (L.spec.soft)
          out << (w.ownerSlot == -2 ? "  soft (forbidden, target 0)" : "  soft");
        else if (w.ownerSlot == -2)
          out << (w.mix >= 0 ? "  forbidden, mix " + std::to_string(w.mix + 1) : "  forbidden");
        else
          out << "  derived";
        out << "   [" << w.key << "]" << std::endl;
      }
    }
    out << std::endl;
  }
}
