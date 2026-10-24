#include "IsospinLevels.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>

#include "AngCoeff.h"
#include "CNuc.h"
#include "Config.h"

namespace {

std::string Trim(const std::string &s) {
  size_t a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return std::string();
  size_t b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}

bool IsHalfInteger(double x) { return std::fabs(2.0 * x - std::round(2.0 * x)) < 1e-9 && x >= 0.0; }

/// Isospin and projection of a particle with charge z and mass m (MeV-free:
/// the mass number is round(m)), t = |N - Z|/2, t3 = (Z - N)/2.
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

}  // namespace

int IsospinLevels::Read(const Config &configure) {
  active_ = false;
  shifts_.clear();
  levelIDs_.clear();
  levelT_.clear();
  mixes_.clear();
  factors_.clear();
  pairTs_.clear();
  levels_.clear();
  derived_.clear();
  softTerms_.clear();
  softSpec_.clear();
  std::ifstream in(configure.configfile.c_str());
  if (!in) return -1;
  std::string line;
  bool found = false;
  while (getline(in, line))
    if (Trim(line) == "<isospin>") {
      found = true;
      break;
    }
  if (!found) return 0;
  active_ = true;
  auto fail = [&](const std::string &raw, const std::string &why) {
    configure.outStream << "ERROR: <isospin> line \"" << raw << "\": " << why << std::endl;
    return -1;
  };
  while (getline(in, line)) {
    std::string raw = Trim(line);
    if (raw == "</isospin>") return 0;
    std::string body = raw.substr(0, raw.find('#'));
    if (Trim(body).empty()) continue;
    std::istringstream s(body);
    std::string kind;
    s >> kind;
    if (kind == "reference" || kind == "shift") {
      IsospinShift e;
      s >> e.z >> e.a >> e.pairKey;
      e.reference = (kind == "reference");
      if (!e.reference) {
        int vary = 0;
        s >> e.value >> vary >> e.nominal >> e.sigma;
        e.vary = (vary != 0);
        if (e.sigma < 0.0) return fail(raw, "the prior width must be >= 0");
      }
      if (s.fail()) return fail(raw, kind == "reference" ? "expected: reference <Z> <A> <pair key>"
                                                         : "expected: shift <Z> <A> <pair key> <value> <vary> <prior centre> <prior sigma>");
      shifts_.push_back(e);
    } else if (kind == "level") {
      int id = 0;
      double t = -1.0;
      s >> id >> t;
      if (s.fail()) return fail(raw, "expected: level <level ID> <T>");
      if (!IsHalfInteger(t)) return fail(raw, "T must be a non-negative integer or half-integer");
      if (std::find(levelIDs_.begin(), levelIDs_.end(), id) != levelIDs_.end()) return fail(raw, "level ID declared twice");
      levelIDs_.push_back(id);
      levelT_.push_back(t);
    } else if (kind == "mix") {
      IsospinMix m;
      int vary = 0;
      s >> m.levelID >> m.pairKey >> m.l >> m.s >> m.value >> vary >> m.sigma;
      m.vary = (vary != 0);
      if (s.fail()) return fail(raw, "expected: mix <level ID> <pair key> <l> <s> <value> <vary> <prior sigma>");
      if (m.sigma < 0.0) return fail(raw, "the prior width must be >= 0");
      mixes_.push_back(m);
    } else if (kind == "factor") {
      IsospinFactor f;
      s >> f.levelID >> f.pairKey >> f.value;
      if (s.fail()) return fail(raw, "expected: factor <level ID> <pair key> <value>");
      factors_.push_back(f);
    } else if (kind == "pairT") {
      IsospinPairT p;
      s >> p.pairKey >> p.t1 >> p.t2;
      if (s.fail() || !IsHalfInteger(p.t1) || !IsHalfInteger(p.t2))
        return fail(raw, "expected: pairT <pair key> <t1> <t2> (non-negative integers or half-integers)");
      pairTs_.push_back(p);
    } else if (kind == "soft") {
      int id = 0;
      double sE = 0.0, sG = 0.0;
      s >> id >> sE >> sG;
      if (s.fail() || sE <= 0.0 || sG <= 0.0)
        return fail(raw, "expected: soft <level ID> <sigma E (MeV)> <sigma gamma (MeV^1/2)>, both > 0");
      softSpec_.push_back({id, {sE, sG}});
    } else {
      return fail(raw, "unknown keyword \"" + kind + "\" (reference, shift, level, mix, factor, pairT, soft)");
    }
  }
  configure.outStream << "ERROR: <isospin> block has no closing </isospin>." << std::endl;
  return -1;
}

bool IsospinLevels::IsDeclared(int levelID) const {
  return std::find(levelIDs_.begin(), levelIDs_.end(), levelID) != levelIDs_.end();
}

const IsospinShift *IsospinLevels::ShiftOf(int nucleus) const {
  for (const IsospinShift &e : shifts_)
    if (e.nucleus == nucleus) return &e;
  return nullptr;
}

double IsospinLevels::ShiftValue(int nucleus, const vector_r &q) const {
  const IsospinShift *e = ShiftOf(nucleus);
  if (!e || e->reference) return 0.0;
  return q[e->param];
}

int IsospinLevels::Resolve(const Config &configure, CNuc &c) {
  if (!active_) return 0;
  std::ostream &out = configure.outStream;
  auto fail = [&](const std::string &why) {
    out << "ERROR: <isospin>: " << why << std::endl;
    return -1;
  };

  // Hale's analysis is a formal fit with fixed boundary values; the ties are
  // linear in the formal parameters only.
  if (configure.paramMask & Config::USE_BRUNE_FORMALISM)
    return fail("isospin levels need the standard (formal) parametrization: run with --no-brune (Brune is the default).");
  if (configure.paramMask & Config::USE_PARK_FORMALISM)
    return fail("isospin levels need the standard (formal) parametrization, not the Park formalism.");
  if (configure.paramMask & Config::USE_RMC_FORMALISM)
    return fail("isospin levels are not available with the RMC formalism.");

  // Nuclei of the reference and shift lines, and their thresholds.
  int nRef = 0;
  for (IsospinShift &e : shifts_) {
    for (int k = 1; k <= c.NumNuclei(); k++)
      if (c.GetNucleus(k)->Z == e.z && c.GetNucleus(k)->A == e.a) e.nucleus = k;
    if (!e.nucleus) return fail("no particle pair forms the nucleus Z=" + std::to_string(e.z) + " A=" + std::to_string(e.a) + ".");
    const int pn = c.GetPairNumFromKey(e.pairKey);
    if (!pn || c.GetPair(pn)->GetNucleus() != e.nucleus || c.GetPair(pn)->GetPType() != 0)
      return fail("pair " + std::to_string(e.pairKey) + " is not a particle pair of " + c.NucleusLabel(e.nucleus) + ".");
    e.threshold = c.GetPair(pn)->GetSepE() + c.GetPair(pn)->GetExE();
    if (e.reference) nRef++;
    for (const IsospinShift &o : shifts_)
      if (&o != &e && o.z == e.z && o.a == e.a) return fail(c.NucleusLabel(e.nucleus) + " has two reference/shift lines.");
  }
  if (nRef != 1) return fail("exactly one 'reference <Z> <A> <pair key>' line is needed.");
  int refNucleus = 0;
  for (const IsospinShift &e : shifts_)
    if (e.reference) refNucleus = e.nucleus;

  // Flat parameter slots, in CNuc::FillMnParams order.
  std::map<std::pair<int, int>, int> energySlot;
  int slot = 0;
  for (int j = 1; j <= c.NumJGroups(); j++)
    for (int la = 1; la <= c.GetJGroup(j)->NumLevels(); la++) {
      energySlot[{j, la}] = slot;
      slot += 1 + c.GetJGroup(j)->NumChannels();
    }
  firstParam_ = slot + c.NumAbundanceParams();
  int k = firstParam_;
  for (IsospinShift &e : shifts_)
    if (!e.reference) e.param = k++;
  for (IsospinMix &m : mixes_) m.param = k++;

  auto pairIsospin = [&](int pairNum, double &t1, double &t31, double &t2, double &t32) {
    PPair *pp = c.GetPair(pairNum);
    double a1, a2;
    ParticleIsospin(pp->GetZ(1), pp->GetM(1), a1, t31);
    ParticleIsospin(pp->GetZ(2), pp->GetM(2), a2, t32);
    t1 = a1;
    t2 = a2;
    for (const IsospinPairT &o : pairTs_)
      if (o.pairKey == pp->GetPairKey()) {
        t1 = o.t1;
        t2 = o.t2;
      }
  };
  // Pairs whose particles match in mass number, spin, parity and isospin form one arrangement.
  auto arrangement = [&](int pairNum) {
    PPair *pp = c.GetPair(pairNum);
    double t1, t31, t2, t32;
    pairIsospin(pairNum, t1, t31, t2, t32);
    std::ostringstream o;
    o << "A" << std::lround(pp->GetM(1)) << "+A" << std::lround(pp->GetM(2)) << " j" << pp->GetJ(1) << (pp->GetPi(1) < 0 ? "-" : "+")
      << "," << pp->GetJ(2) << (pp->GetPi(2) < 0 ? "-" : "+") << " t" << t1 << "," << t2;
    return o.str();
  };

  levels_.clear();
  derived_.clear();
  std::map<std::string, double> arrangementRadius;
  for (size_t d = 0; d < levelIDs_.size(); d++) {
    IsospinLevel L;
    L.id = levelIDs_[d];
    L.t = levelT_[d];
    const std::string idText = "level ID " + std::to_string(L.id);
    for (int j = 1; j <= c.NumJGroups(); j++)
      for (int la = 1; la <= c.GetJGroup(j)->NumLevels(); la++) {
        ALevel *lev = c.GetJGroup(j)->GetLevel(la);
        if (lev->GetLevelID() != L.id || !lev->IsInRMatrix()) continue;
        IsospinMember m;
        m.nucleus = c.GetJGroup(j)->GetNucleus();
        m.jGroup = j;
        m.level = la;
        m.energySlot = energySlot[{j, la}];
        for (const IsospinMember &o : L.members)
          if (o.nucleus == m.nucleus)
            return fail(idText + " has lines in two J-groups of " + c.NucleusLabel(m.nucleus) +
                        ": the lines of one level must share J and parity.");
        const int z = c.GetJGroup(j)->GetNucleusZ(), a = c.GetJGroup(j)->GetNucleusA();
        m.t3 = 0.5 * (z - (a - z));
        if (std::fabs(m.t3) > L.t + 1e-9)
          return fail(idText + " has T=" + Fmt(L.t) + " but lines in " + c.NucleusLabel(m.nucleus) + ", whose T3=" + Fmt(m.t3) + ".");
        L.members.push_back(m);
      }
    if (L.members.empty()) return fail(idText + " has no active level lines.");
    const JGroup *g0 = c.GetJGroup(L.members[0].jGroup);
    for (const IsospinMember &m : L.members) {
      JGroup *g = c.GetJGroup(m.jGroup);
      if (g->GetJ() != g0->GetJ() || g->GetPi() != g0->GetPi())
        return fail(idText + ": its lines have different J or parity in different nuclei.");
    }
    L.ref = 0;
    for (size_t i = 0; i < L.members.size(); i++)
      if (L.members[i].nucleus == refNucleus) L.ref = (int)i;
    if (L.members.size() > 1)
      for (IsospinMember &m : L.members) {
        const IsospinShift *e = ShiftOf(m.nucleus);
        if (!e)
          return fail(idText + " has lines in " + c.NucleusLabel(m.nucleus) +
                      ", which has no reference or shift line (its energies need a threshold and a shift).");
        m.threshold = e->threshold;
      }

    // Widths: Clebsch-Gordan factors, isospin-channel keys.
    for (IsospinMember &m : L.members) {
      JGroup *g = c.GetJGroup(m.jGroup);
      g->SetBoundaryMinusL(true);
      for (int ch = 1; ch <= g->NumChannels(); ch++) {
        AChannel *chan = g->GetChannel(ch);
        if (chan->GetRadType() != 'P')
          return fail(idText + ": the J-group in " + c.NucleusLabel(m.nucleus) +
                      " has a photon or beta channel; isospin levels take particle channels only.");
        IsospinWidth w;
        w.channel = ch;
        w.slot = m.energySlot + ch;
        double t1, t31, t2, t32;
        pairIsospin(chan->GetPairNum(), t1, t31, t2, t32);
        if (std::fabs(t31 + t32 - m.t3) > 1e-9)
          return fail("pair " + std::to_string(c.GetPair(chan->GetPairNum())->GetPairKey()) +
                      ": the particle isospin projections do not add up to the nucleus's T3.");
        w.cgComputed = (std::fabs(t31) <= t1 + 1e-9 && std::fabs(t32) <= t2 + 1e-9 && L.t <= t1 + t2 + 1e-9 &&
                        L.t >= std::fabs(t1 - t2) - 1e-9)
                           ? AngCoeff::ClebGord(t1, t2, L.t, t31, t32, m.t3)
                           : 0.0;
        if (std::fabs(w.cgComputed) < 1e-12) w.cgComputed = 0.0;
        w.cg = w.cgComputed;
        for (const IsospinFactor &f : factors_)
          if (f.levelID == L.id && f.pairKey == c.GetPair(chan->GetPairNum())->GetPairKey()) {
            w.cg = f.value;
            w.overridden = true;
          }
        const std::string arr = arrangement(chan->GetPairNum());
        std::ostringstream key;
        key << arr << " l=" << chan->GetL() << " s=" << chan->GetS();
        w.key = key.str();
        const double radius = c.GetPair(chan->GetPairNum())->GetChRad();
        auto rad = arrangementRadius.find(arr);
        if (rad == arrangementRadius.end())
          arrangementRadius[arr] = radius;
        else if (std::fabs(rad->second - radius) > 1e-9)
          return fail("pairs of one arrangement (" + arr + ") have channel radii " + Fmt(rad->second, 6) + " and " +
                      Fmt(radius, 6) + " fm; Hale's isospin relations need equal radii.");
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
      for (IsospinWidth &w : L.members[i].widths)
        if (w.cg != 0.0 && !owner.count(w.key)) {
          owner[w.key] = {w.slot, w.cg};
          w.ownerSlot = -1;
        }
    for (int i : order)
      for (IsospinWidth &w : L.members[i].widths) {
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
    levels_.push_back(L);
  }

  // Soft ties.
  for (const auto &sp : softSpec_) {
    IsospinLevel *L = nullptr;
    for (IsospinLevel &l : levels_)
      if (l.id == sp.first) L = &l;
    if (!L) return fail("soft names level ID " + std::to_string(sp.first) + ", which has no 'level' line.");
    L->soft = true;
    L->sigmaE = sp.second.first;
    L->sigmaGamma = sp.second.second;
    for (const IsospinMix &mx : mixes_)
      if (mx.levelID == L->id)
        return fail("level ID " + std::to_string(L->id) + " has both a soft tie and a mix line; the soft tie already "
                    "frees its forbidden channels.");
  }
  if (!softSpec_.empty() && (configure.paramMask & (Config::USE_LM_MINIMIZER | Config::USE_GSL_LM_MINIMIZER)))
    return fail("the LM / GSL least-squares fitters take per-parameter priors only, and a soft tie couples parameters; "
                "fit with MIGRAD (with or without --use-gradient).");

  // Mix and factor lines must name existing level lines.
  for (size_t i = 0; i < mixes_.size(); i++) {
    IsospinMix &mx = mixes_[i];
    if (!IsDeclared(mx.levelID)) return fail("mix names level ID " + std::to_string(mx.levelID) + ", which has no 'level' line.");
    for (IsospinLevel &L : levels_) {
      if (L.id != mx.levelID) continue;
      for (IsospinMember &m : L.members)
        for (IsospinWidth &w : m.widths) {
          AChannel *chan = c.GetJGroup(m.jGroup)->GetChannel(w.channel);
          if (c.GetPair(chan->GetPairNum())->GetPairKey() != mx.pairKey || chan->GetL() != mx.l ||
              std::fabs(chan->GetS() - mx.s) > 1e-9)
            continue;
          if (w.ownerSlot != -2)
            return fail("mix on level ID " + std::to_string(mx.levelID) + ", pair " + std::to_string(mx.pairKey) +
                        ": this channel is isospin-allowed (factor " + Fmt(w.cg) + "); mix is for forbidden channels.");
          w.mix = (int)i;
          mx.slot = w.slot;
        }
    }
    if (mx.slot < 0)
      return fail("mix on level ID " + std::to_string(mx.levelID) + ": no channel with pair " + std::to_string(mx.pairKey) +
                  ", l=" + std::to_string(mx.l) + ", s=" + Fmt(mx.s) + ".");
  }
  for (const IsospinFactor &f : factors_) {
    bool used = false;
    for (const IsospinLevel &L : levels_)
      for (const IsospinMember &m : L.members)
        for (const IsospinWidth &w : m.widths)
          if (L.id == f.levelID && w.overridden &&
              c.GetPair(c.GetJGroup(m.jGroup)->GetChannel(w.channel)->GetPairNum())->GetPairKey() == f.pairKey)
            used = true;
    if (!used)
      return fail("factor on level ID " + std::to_string(f.levelID) + ", pair " + std::to_string(f.pairKey) +
                  " matches no channel of a declared level.");
  }

  // Derived slots (fixed, overwritten), or soft-tie terms for a soft level.
  auto shiftParamOf = [&](int nucleus) {
    const IsospinShift *e = ShiftOf(nucleus);
    return (e && !e->reference) ? e->param : -1;
  };
  for (const IsospinLevel &L : levels_) {
    const IsospinMember &r = L.members[L.ref];
    for (size_t i = 0; i < L.members.size(); i++) {
      const IsospinMember &m = L.members[i];
      if ((int)i != L.ref && L.members.size() > 1) {
        if (!L.soft) {
          derived_.insert(m.energySlot);
        } else {
          IsospinSoftTerm t;
          t.slot = m.energySlot;
          t.sigma = L.sigmaE;
          t.constant = m.threshold - r.threshold;
          t.sources.push_back({r.energySlot, 1.0});
          if (shiftParamOf(m.nucleus) >= 0) t.sources.push_back({shiftParamOf(m.nucleus), 1.0});
          if (shiftParamOf(r.nucleus) >= 0) t.sources.push_back({shiftParamOf(r.nucleus), -1.0});
          softTerms_.push_back(t);
        }
      }
      for (const IsospinWidth &w : m.widths) {
        if (w.ownerSlot == -1) continue;
        if (!L.soft) {
          derived_.insert(w.slot);
        } else {
          IsospinSoftTerm t;
          t.slot = w.slot;
          t.sigma = L.sigmaGamma;
          if (w.ownerSlot >= 0) t.sources.push_back({w.ownerSlot, w.ratio});
          softTerms_.push_back(t);
        }
      }
    }
  }

  // Report.
  out << "Isospin levels: " << levels_.size() << " level(s); reference " << c.NucleusLabel(refNucleus)
      << "; B_c = -l_c in their J-groups." << std::endl;
  out << "  note: the GUI does not keep the <isospin> block yet; edit this project as text." << std::endl;
  for (const IsospinLevel &L : levels_) {
    JGroup *g = c.GetJGroup(L.members[L.ref].jGroup);
    out << "  level ID " << L.id << "  T=" << L.t << "  J=" << g->GetJ() << (g->GetPi() < 0 ? "-" : "+")
        << (L.soft ? "  SOFT (sigma E " + Fmt(L.sigmaE) + " MeV, sigma gamma " + Fmt(L.sigmaGamma) + " MeV^1/2)" : "") << "  in";
    for (size_t i = 0; i < L.members.size(); i++)
      out << " " << c.NucleusLabel(L.members[i].nucleus) << ((int)i == L.ref ? " (reference)" : "");
    out << std::endl;
    for (const IsospinMember &m : L.members)
      for (const IsospinWidth &w : m.widths) {
        AChannel *chan = c.GetJGroup(m.jGroup)->GetChannel(w.channel);
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
  }
  // Labelled levels' groups may hold unlabelled levels: say so.
  std::set<int> groups;
  for (const IsospinLevel &L : levels_)
    for (const IsospinMember &m : L.members) groups.insert(m.jGroup);
  for (int j : groups)
    for (int la = 1; la <= c.GetJGroup(j)->NumLevels(); la++)
      if (!IsDeclared(c.GetJGroup(j)->GetLevel(la)->GetLevelID()))
        out << "  note: level ID " << c.GetJGroup(j)->GetLevel(la)->GetLevelID() << " in a J-group of "
            << c.NucleusLabel(c.GetJGroup(j)->GetNucleus())
            << " has no isospin label; it is unconstrained (and uses B = -l like its group)." << std::endl;
  // A nucleus the level's isospin allows, with pairs of its arrangements, but no line.
  for (const IsospinLevel &L : levels_) {
    std::set<std::string> arrs;
    std::set<int> present;
    for (const IsospinMember &m : L.members) {
      present.insert(m.nucleus);
      for (const IsospinWidth &w : m.widths) arrs.insert(arrangement(c.GetJGroup(m.jGroup)->GetChannel(w.channel)->GetPairNum()));
    }
    for (int n = 1; n <= c.NumNuclei(); n++) {
      if (present.count(n)) continue;
      const double t3 = 0.5 * (c.GetNucleus(n)->Z - (c.GetNucleus(n)->A - c.GetNucleus(n)->Z));
      if (std::fabs(t3) > L.t + 1e-9) continue;
      for (int p = 1; p <= c.NumPairs(); p++)
        if (c.GetPair(p)->GetNucleus() == n && c.GetPair(p)->GetPType() == 0 && arrs.count(arrangement(p))) {
          out << "  WARNING: level ID " << L.id << " (T=" << L.t << ") has no lines in " << c.NucleusLabel(n)
              << ", which its isospin allows; it is absent there." << std::endl;
          break;
        }
    }
  }
  return 0;
}

int IsospinLevels::NumParams() const {
  if (!active_) return 0;
  int n = (int)mixes_.size();
  for (const IsospinShift &e : shifts_)
    if (!e.reference) n++;
  return n;
}

void IsospinLevels::AddParams(ROOT::Minuit2::MnUserParameters &p) {
  if (!active_) return;
  for (const IsospinShift &e : shifts_) {
    if (e.reference) continue;
    std::string name = "isospin_shift_" + std::to_string(e.z) + "_" + std::to_string(e.a);
    p.Add(name, e.value, e.sigma > 0.0 ? e.sigma : 0.1);
    if (!e.vary) p.Fix(name);
  }
  for (size_t i = 0; i < mixes_.size(); i++) {
    const IsospinMix &m = mixes_[i];
    std::string name = "isospin_mix_" + std::to_string(i + 1);
    p.Add(name, m.value, m.sigma > 0.0 ? m.sigma : 0.1 * std::fabs(m.value) + 1e-3);
    if (!m.vary) p.Fix(name);
  }
}

void IsospinLevels::Apply(vector_r &q) const {
  if (!active_) return;
  for (const IsospinLevel &L : levels_) {
    if (L.soft) continue;  // free slots, penalised instead (SoftPenalty)
    const IsospinMember &r = L.members[L.ref];
    if (L.members.size() > 1) {
      const double eps = q[r.energySlot] - r.threshold - ShiftValue(r.nucleus, q);
      for (size_t i = 0; i < L.members.size(); i++) {
        if ((int)i == L.ref) continue;
        const IsospinMember &m = L.members[i];
        q[m.energySlot] = eps + m.threshold + ShiftValue(m.nucleus, q);
      }
    }
    for (const IsospinMember &m : L.members)
      for (const IsospinWidth &w : m.widths) {
        if (w.ownerSlot >= 0)
          q[w.slot] = w.ratio * q[w.ownerSlot];
        else if (w.ownerSlot == -2)
          q[w.slot] = (w.mix >= 0) ? q[mixes_[w.mix].param] : 0.0;
      }
  }
}

std::vector<IsospinTie> IsospinLevels::Ties() const {
  std::vector<IsospinTie> ties;
  if (!active_) return ties;
  auto shiftParam = [&](int nucleus) {
    const IsospinShift *e = ShiftOf(nucleus);
    return (e && !e->reference) ? e->param : -1;
  };
  for (const IsospinLevel &L : levels_) {
    if (L.soft) continue;
    const IsospinMember &r = L.members[L.ref];
    if (L.members.size() > 1)
      for (size_t i = 0; i < L.members.size(); i++) {
        if ((int)i == L.ref) continue;
        const IsospinMember &m = L.members[i];
        IsospinTie t;
        t.slot = m.energySlot;
        t.sources.push_back({r.energySlot, 1.0});
        if (shiftParam(m.nucleus) >= 0) t.sources.push_back({shiftParam(m.nucleus), 1.0});
        if (shiftParam(r.nucleus) >= 0) t.sources.push_back({shiftParam(r.nucleus), -1.0});
        ties.push_back(t);
      }
    for (const IsospinMember &m : L.members)
      for (const IsospinWidth &w : m.widths) {
        if (w.ownerSlot == -1) continue;
        IsospinTie t;
        t.slot = w.slot;
        if (w.ownerSlot >= 0)
          t.sources.push_back({w.ownerSlot, w.ratio});
        else if (w.mix >= 0)
          t.sources.push_back({mixes_[w.mix].param, 1.0});
        ties.push_back(t);
      }
  }
  return ties;
}

void IsospinLevels::ParamPrior(int k, double &centre, double &sigma) const {
  centre = 0.0;
  sigma = 0.0;
  const int idx = firstParam_ + k;
  for (const IsospinShift &e : shifts_)
    if (!e.reference && e.param == idx) {
      centre = e.nominal;
      sigma = e.sigma;
    }
  for (const IsospinMix &m : mixes_)
    if (m.param == idx) sigma = m.sigma;
}

std::vector<double> IsospinLevels::PenaltyGradient() const {
  std::vector<double> g(NumParams(), 0.0);
  if (!active_) return g;
  for (const IsospinShift &e : shifts_)
    if (!e.reference && e.vary && e.sigma > 0.0) g[e.param - firstParam_] = 2.0 * (e.value - e.nominal) / (e.sigma * e.sigma);
  for (const IsospinMix &m : mixes_)
    if (m.vary && m.sigma > 0.0) g[m.param - firstParam_] = 2.0 * m.value / (m.sigma * m.sigma);
  return g;
}

void IsospinLevels::ReadParams(const vector_r &q) {
  if (!active_) return;
  if (!softTerms_.empty()) lastQ_ = q;
  for (IsospinShift &e : shifts_)
    if (!e.reference && e.param >= 0 && e.param < (int)q.size()) e.value = q[e.param];
  for (IsospinMix &m : mixes_)
    if (m.param >= 0 && m.param < (int)q.size()) m.value = q[m.param];
}

void IsospinLevels::Sync(ROOT::Minuit2::MnUserParameters &p) const {
  if (!active_) return;
  vector_r q = p.Params();
  Apply(q);
  for (int i : derived_) {
    if (i >= (int)q.size()) continue;
    if (!p.Parameter(i).IsFixed()) p.Fix(i);
    p.SetValue(i, q[i]);
  }
}

double IsospinLevels::Penalty() const {
  double chi2 = 0.0;
  if (!active_) return chi2;
  for (const IsospinShift &e : shifts_)
    if (!e.reference && e.vary && e.sigma > 0.0) chi2 += std::pow((e.value - e.nominal) / e.sigma, 2.0);
  for (const IsospinMix &m : mixes_)
    if (m.vary && m.sigma > 0.0) chi2 += std::pow(m.value / m.sigma, 2.0);
  return chi2;
}

bool IsospinLevels::HasPenalty() const {
  if (!active_) return false;
  for (const IsospinShift &e : shifts_)
    if (!e.reference && e.vary && e.sigma > 0.0) return true;
  for (const IsospinMix &m : mixes_)
    if (m.vary && m.sigma > 0.0) return true;
  return false;
}

double IsospinLevels::SoftPenalty(const vector_r &q) const {
  double chi2 = 0.0;
  for (const IsospinSoftTerm &t : softTerms_) {
    if (t.slot >= (int)q.size()) continue;
    double target = t.constant;
    for (const auto &s : t.sources) target += s.second * q[s.first];
    chi2 += std::pow((q[t.slot] - target) / t.sigma, 2.0);
  }
  return chi2;
}

std::vector<std::pair<int, double>> IsospinLevels::SoftPenaltyGradient(const vector_r &q) const {
  std::vector<std::pair<int, double>> g;
  for (const IsospinSoftTerm &t : softTerms_) {
    if (t.slot >= (int)q.size()) continue;
    double target = t.constant;
    for (const auto &s : t.sources) target += s.second * q[s.first];
    const double d = 2.0 * (q[t.slot] - target) / (t.sigma * t.sigma);
    g.push_back({t.slot, d});
    for (const auto &s : t.sources) g.push_back({s.first, -s.second * d});
  }
  return g;
}

void IsospinLevels::AddSoftCurvature(std::vector<double> &g2) const {
  for (const IsospinSoftTerm &t : softTerms_) {
    const double c = 2.0 / (t.sigma * t.sigma);
    if (t.slot < (int)g2.size()) g2[t.slot] += c;
    for (const auto &s : t.sources)
      if (s.first < (int)g2.size()) g2[s.first] += c * s.second * s.second;
  }
}

void IsospinLevels::InitSoftParams(ROOT::Minuit2::MnUserParameters &p) const {
  if (softTerms_.empty()) return;
  const vector_r q = p.Params();
  for (const IsospinSoftTerm &t : softTerms_) {
    if (t.slot >= (int)q.size()) continue;
    double target = t.constant;
    for (const auto &s : t.sources) target += s.second * q[s.first];
    p.SetValue(t.slot, target);
    // A width at exactly 0 is otherwise fixed (CNuc::FillMnParams); a soft slot starts free.
    if (p.Parameter(t.slot).IsFixed()) p.Release(t.slot);
    p.SetError(t.slot, std::max(0.1 * std::fabs(target), t.sigma));
  }
}

std::string IsospinLevels::ParamLabel(int i) const {
  if (!active_) return std::string();
  for (const IsospinShift &e : shifts_)
    if (!e.reference && e.param == i)
      return "internal Coulomb shift of isospin levels in Z=" + std::to_string(e.z) + " A=" + std::to_string(e.a);
  for (size_t k = 0; k < mixes_.size(); k++)
    if (mixes_[k].param == i)
      return "isospin-forbidden amplitude " + std::to_string(k + 1) + ": level ID " + std::to_string(mixes_[k].levelID) +
             ", pair " + std::to_string(mixes_[k].pairKey) + ", l=" + std::to_string(mixes_[k].l) + ", s=" + Fmt(mixes_[k].s);
  return std::string();
}

std::string IsospinLevels::MemberLabel(int jGroup, int level) const {
  if (!active_) return std::string();
  for (const IsospinLevel &L : levels_)
    for (size_t i = 0; i < L.members.size(); i++)
      if (L.members[i].jGroup == jGroup && L.members[i].level == level)
        return "ISOSPIN LEVEL " + std::to_string(L.id) + " T=" + Fmt(L.t) +
               ((int)i == L.ref ? " (reference)" : (L.soft ? " (soft)" : " (derived)"));
  return std::string();
}

void IsospinLevels::Write(const Config &configure, CNuc &c) const {
  if (!active_) return;
  std::ofstream out((configure.outputdir + "isospin.out").c_str());
  if (!out) return;
  out << "ISOSPIN LEVELS (formal parameters, B_c = -l_c in their J-groups)" << std::endl << std::endl;
  if (!softTerms_.empty()) out << "soft-tie penalty = " << std::setprecision(8) << SoftPenalty() << std::endl;
  for (const IsospinShift &e : shifts_) {
    out << (e.reference ? "reference " : "shift     ") << "Z=" << e.z << " A=" << e.a << "  threshold pair " << e.pairKey
        << " at Ex = " << std::setprecision(6) << e.threshold << " MeV";
    if (!e.reference) out << "  shift = " << std::setprecision(6) << e.value << " MeV" << (e.vary ? " (varied)" : " (fixed)");
    out << std::endl;
  }
  for (size_t k = 0; k < mixes_.size(); k++)
    out << "mix " << k + 1 << "  level ID " << mixes_[k].levelID << "  pair " << mixes_[k].pairKey << "  l=" << mixes_[k].l
        << " s=" << mixes_[k].s << "  amplitude = " << std::setprecision(6) << mixes_[k].value << " MeV^1/2"
        << (mixes_[k].vary ? " (varied)" : " (fixed)") << std::endl;
  out << std::endl;
  for (const IsospinLevel &L : levels_) {
    JGroup *g = c.GetJGroup(L.members[L.ref].jGroup);
    out << "level ID " << L.id << "  T=" << L.t << "  J=" << g->GetJ() << (g->GetPi() < 0 ? "-" : "+") << std::endl;
    for (size_t i = 0; i < L.members.size(); i++) {
      const IsospinMember &m = L.members[i];
      ALevel *lev = c.GetJGroup(m.jGroup)->GetLevel(m.level);
      out << "  " << c.NucleusLabel(m.nucleus) << " (T3=" << m.t3 << ")" << ((int)i == L.ref ? "  reference" : "  derived")
          << "  formal E = " << std::setprecision(14) << lev->GetFitE() << " MeV" << std::endl;
      if (L.soft && (int)i != L.ref) out << "    (soft tie)" << std::endl;
      for (const IsospinWidth &w : m.widths) {
        AChannel *chan = c.GetJGroup(m.jGroup)->GetChannel(w.channel);
        out << "    pair " << std::setw(2) << c.GetPair(chan->GetPairNum())->GetPairKey() << "  l=" << chan->GetL()
            << " s=" << chan->GetS() << "  B=" << chan->GetBoundaryCondition() << "  factor " << std::setw(10)
            << std::setprecision(6) << w.cg;
        if (w.overridden) out << " (computed " << w.cgComputed << ")";
        out << "  gamma = " << std::setw(22) << std::setprecision(14) << lev->GetFitGamma(w.channel) << " MeV^1/2";
        if (w.ownerSlot == -1)
          out << "  free";
        else if (L.soft)
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
