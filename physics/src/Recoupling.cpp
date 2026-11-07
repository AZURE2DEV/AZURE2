// Convention conversions of level_parameters (LevelParameters.h; docs/upgrade/level_parameters.md section 5): the
// matrix T with gamma_internal = T gamma_project on the gamma columns of one J^pi block.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

#include "azure2/physics/AngularCoupling.h"
#include "azure2/physics/LevelParameters.h"

namespace azure2 {
namespace physics {

namespace {

constexpr double kZero = 1e-14;  // a coupling coefficient below this is a selection-rule zero

int PhaseOfTwice(int twice) {  // (-1)^(twice/2) for an even argument
  if (twice % 2 != 0) throw std::logic_error("recoupling: half-integer phase exponent");
  return ((twice / 2) % 2 == 0) ? 1 : -1;
}

// <((ja jb) s, l) J | ((ja l) j, jb) J> = (-1)^(4 ja + 3 jb + 3 l - s - j) shat jhat {ja jb s; J l j}.
// Derived from the standard (ab)e,d <-> a,(bd)f recoupling with two order swaps; checked against an m-scheme sum of
// Clebsch-Gordan products for all 2j <= 4, l <= 4 (physics/tests/angular_coupling_test.cpp, 1831 cases here when it
// was derived).
double SpinFirstJJ(int tja, int tjb, int l, int ts, int tJ, int tj) {
  const double sixj = wigner_6j(tja, tjb, ts, tJ, 2 * l, tj);
  if (sixj == 0.0) return 0.0;
  const int phase = PhaseOfTwice(4 * tja + 3 * tjb + 6 * l - ts - tj);
  return phase * std::sqrt((double)((ts + 1) * (tj + 1))) * sixj;
}

// (-1)^floor(l/2): the real sign that represents i^(l - l') between channels of equal l parity (and is common to
// a partition otherwise); level_parameters.md section 5.
int OrbitalPhaseSign(OrbitalPhase p, int l) {
  if (p == OrbitalPhase::None) return 1;
  return ((l / 2) % 2 == 0) ? 1 : -1;
}

struct GroupKey {
  int partition, l;
  bool operator<(const GroupKey &o) const { return std::tie(partition, l) < std::tie(o.partition, o.l); }
};

// The complete set of 2s (internal) or 2x (project) for one (partition, l) at total 2J.
std::vector<int> CompleteInternal(int tj1, int tj2, int l, int tJ) {
  std::vector<int> out;
  for (int ts = std::abs(tj1 - tj2); ts <= tj1 + tj2; ts += 2)
    if (triangle(ts, 2 * l, tJ)) out.push_back(ts);
  return out;
}

std::vector<int> CompleteProject(int tj1, int tj2, int l, int tJ, const Convention &c) {
  if (c.coupling == CouplingScheme::LS) return CompleteInternal(tj1, tj2, l, tJ);
  const int tja = (c.spin_order == SpinOrder::LightFirst) ? tj1 : tj2;
  const int tjb = (c.spin_order == SpinOrder::LightFirst) ? tj2 : tj1;
  std::vector<int> out;
  for (int tj = std::abs(2 * l - tja); tj <= 2 * l + tja; tj += 2)
    if (triangle(tj, tjb, tJ)) out.push_back(tj);
  return out;
}

std::string ChannelText(int partition, int l, int twice, const char *name) {
  std::ostringstream o;
  o << "(partition " << partition << ", l = " << l << ", 2" << name << " = " << twice << ")";
  return o.str();
}

}  // namespace

double recoupling_coefficient(int tj1, int tj2, int l, int ts, int tJ, const Convention &c, int tx) {
  if (!triangle(tj1, tj2, ts) || !triangle(ts, 2 * l, tJ)) return 0.0;
  const int il = OrbitalPhaseSign(c.orbital_phase, l);
  if (c.coupling == CouplingScheme::LS) {
    if (tx != ts) return 0.0;
    int phase = il;
    if (c.spin_order == SpinOrder::HeavyFirst) phase *= PhaseOfTwice(tj1 + tj2 - ts);    // (j2 j1)s = (-1)^(j1+j2-s) (j1 j2)s
    if (c.orbital_order == OrbitalOrder::OrbitalFirst) phase *= PhaseOfTwice(2 * l + ts - tJ);  // (l s)J = (-1)^(l+s-J) (s l)J
    return phase;
  }
  // jj: x is the particle whose spin couples with l (LightFirst: x = 1), y the other.
  const bool light = (c.spin_order == SpinOrder::LightFirst);
  const int tja = light ? tj1 : tj2, tjb = light ? tj2 : tj1;
  if (!triangle(2 * l, tja, tx) || !triangle(tx, tjb, tJ)) return 0.0;
  double value = SpinFirstJJ(tja, tjb, l, ts, tJ, tx);
  if (!light) value *= PhaseOfTwice(tj1 + tj2 - ts);                                    // internal s is (j1 j2)
  if (c.orbital_order == OrbitalOrder::OrbitalFirst) value *= PhaseOfTwice(2 * l + tja - tx);  // (l ja)j = (-1)^(l+ja-j) (ja l)j
  return il * value;
}

std::vector<ProjectChannel> project_channels(const System &system, std::size_t b, const Convention &convention) {
  const JpiBlock &blk = system.block(b);
  std::vector<ProjectChannel> out;
  std::set<GroupKey> done;
  for (std::size_t k = 0; k < blk.num_channels(); k++) {
    const Channel &ch = system.channel(blk.channel(k));
    if (ch.kind != ChannelKind::Particle || convention.coupling == CouplingScheme::LS) {
      out.push_back(ProjectChannel{ch.partition, ch.kind, ch.l, ch.two_s});
      continue;
    }
    const GroupKey key{ch.partition, ch.l};
    if (done.count(key)) continue;
    done.insert(key);
    const Partition &p = system.partition(ch.partition);
    for (int tx : CompleteProject(p.two_j[0], p.two_j[1], ch.l, blk.two_J(), convention))
      out.push_back(ProjectChannel{ch.partition, ChannelKind::Particle, ch.l, tx});
  }
  return out;
}

std::vector<std::string> recoupling_gaps(const System &system, std::size_t b, const Convention &convention,
                                         const std::vector<ProjectChannel> &columns) {
  const JpiBlock &blk = system.block(b);
  const int tJ = blk.two_J();
  std::vector<std::string> gaps;
  const char *xName = (convention.coupling == CouplingScheme::LS) ? "s" : "j";
  std::map<GroupKey, std::set<int> > internal, project;
  std::set<std::tuple<int, int, int, int> > internalOther, projectOther;  // (partition, kind, l, 2s)
  for (std::size_t k = 0; k < blk.num_channels(); k++) {
    const Channel &ch = system.channel(blk.channel(k));
    if (ch.kind == ChannelKind::Particle)
      internal[GroupKey{ch.partition, ch.l}].insert(ch.two_s);
    else
      internalOther.insert(std::make_tuple(ch.partition, (int)ch.kind, ch.l, ch.two_s));
  }
  for (std::size_t m = 0; m < columns.size(); m++) {
    const ProjectChannel &col = columns[m];
    if (col.partition < 0 || col.partition >= (int)system.num_partitions()) {
      gaps.push_back("project column " + std::to_string(m) + ": partition index out of range");
      continue;
    }
    if (col.kind == ChannelKind::Particle) {
      if (!project[GroupKey{col.partition, col.l}].insert(col.two_x).second)
        gaps.push_back("project channel " + ChannelText(col.partition, col.l, col.two_x, xName) + " listed twice");
    } else if (!projectOther.insert(std::make_tuple(col.partition, (int)col.kind, col.l, col.two_x)).second) {
      gaps.push_back("project channel " + ChannelText(col.partition, col.l, col.two_x, "s") + " listed twice");
    }
  }
  for (const auto &t : internalOther)
    if (!projectOther.count(t))
      gaps.push_back("block channel " + ChannelText(std::get<0>(t), std::get<2>(t), std::get<3>(t), "s") +
                     " (photon or beta) has no project column");
  for (const auto &t : projectOther)
    if (!internalOther.count(t))
      gaps.push_back("project channel " + ChannelText(std::get<0>(t), std::get<2>(t), std::get<3>(t), "s") +
                     " (photon or beta) is not a channel of the block");

  std::set<GroupKey> keys;
  for (const auto &kv : internal) keys.insert(kv.first);
  for (const auto &kv : project) keys.insert(kv.first);
  for (const GroupKey &key : keys) {
    const Partition &p = system.partition(key.partition);
    const int tj1 = p.two_j[0], tj2 = p.two_j[1];
    const std::vector<int> allS = CompleteInternal(tj1, tj2, key.l, tJ);
    const std::vector<int> allX = CompleteProject(tj1, tj2, key.l, tJ, convention);
    const std::set<int> &haveS = internal[key];
    const std::set<int> &haveX = project[key];
    for (int tx : haveX)
      if (std::find(allX.begin(), allX.end(), tx) == allX.end())
        gaps.push_back("project channel " + ChannelText(key.partition, key.l, tx, xName) +
                       " is not allowed with these spins and J");
    for (int ts : haveS)
      for (int tx : allX)
        if (!haveX.count(tx) && std::fabs(recoupling_coefficient(tj1, tj2, key.l, ts, tJ, convention, tx)) > kZero)
          gaps.push_back("project lacks channel " + ChannelText(key.partition, key.l, tx, xName) +
                         ", which block channel " + ChannelText(key.partition, key.l, ts, "s") + " needs");
    for (int tx : haveX)
      for (int ts : allS)
        if (!haveS.count(ts) && std::fabs(recoupling_coefficient(tj1, tj2, key.l, ts, tJ, convention, tx)) > kZero)
          gaps.push_back("block lacks channel " + ChannelText(key.partition, key.l, ts, "s") +
                         ", which project channel " + ChannelText(key.partition, key.l, tx, xName) + " feeds");
  }
  return gaps;
}

ConventionTransform convention_transform(const System &system, std::size_t b, const Convention &convention,
                                         const std::vector<ProjectChannel> &columns) {
  const std::vector<std::string> gaps = recoupling_gaps(system, b, convention, columns);
  if (!gaps.empty()) {
    std::ostringstream o;
    o << "level_parameters: the recoupling of block " << b << " is incomplete:";
    for (const std::string &g : gaps) o << "\n    " << g;
    throw std::invalid_argument(o.str());
  }
  const JpiBlock &blk = system.block(b);
  const std::size_t nc = blk.num_channels(), np = columns.size();
  std::vector<double> T(nc * np, 0.0);
  for (std::size_t k = 0; k < nc; k++) {
    const Channel &ch = system.channel(blk.channel(k));
    const Partition &p = system.partition(ch.partition);
    for (std::size_t m = 0; m < np; m++) {
      const ProjectChannel &col = columns[m];
      if (col.partition != ch.partition || col.kind != ch.kind || col.l != ch.l) continue;
      if (ch.kind == ChannelKind::Particle)
        T[k * np + m] = recoupling_coefficient(p.two_j[0], p.two_j[1], ch.l, ch.two_s, blk.two_J(), convention, col.two_x);
      else if (col.two_x == ch.two_s)
        T[k * np + m] = (ch.kind == ChannelKind::Electric) ? convention.electric_sign
            : (ch.kind == ChannelKind::Magnetic)           ? convention.magnetic_sign
                                                           : 1.0;
    }
  }
  return ConventionTransform(nc, np, std::move(T));
}

std::vector<double> ConventionTransform::to_internal(const std::vector<double> &gp, std::size_t nl) const {
  if (gp.size() != nl * n_project_) throw std::invalid_argument("ConventionTransform::to_internal: wrong size");
  std::vector<double> out(nl * n_internal_, 0.0);
  for (std::size_t la = 0; la < nl; la++)
    for (std::size_t k = 0; k < n_internal_; k++) {
      double s = 0.0;
      for (std::size_t m = 0; m < n_project_; m++) {
        const double t = T_[k * n_project_ + m];
        if (t != 0.0) s += t * gp[la * n_project_ + m];
      }
      out[la * n_internal_ + k] = s;
    }
  return out;
}

std::vector<double> ConventionTransform::to_project(const std::vector<double> &gi, std::size_t nl) const {
  if (gi.size() != nl * n_internal_) throw std::invalid_argument("ConventionTransform::to_project: wrong size");
  std::vector<double> out(nl * n_project_, 0.0);
  for (std::size_t la = 0; la < nl; la++)
    for (std::size_t m = 0; m < n_project_; m++) {
      double s = 0.0;
      for (std::size_t k = 0; k < n_internal_; k++) {
        const double t = T_[k * n_project_ + m];
        if (t != 0.0) s += t * gi[la * n_internal_ + k];
      }
      out[la * n_project_ + m] = s;
    }
  return out;
}

}  // namespace physics
}  // namespace azure2
