// The pairs and levels of <levels> as plain tables and the rule that finds the pair a label names (LevelTables.h).
#include "azure2/experiment/LevelTables.h"

#include <cmath>
#include <sstream>

#include "azure2/physics/Elements.h"

namespace azure2 {
namespace experiment {

std::string NuclideLabel(int Z, int A) { return std::to_string(A) + physics::ElementSymbol(Z); }

std::string LightParticleLabel(int kind, int Z, int A) {
  if (kind == 10) return "g";
  if (kind == 20) return Z < 0 ? "b-" : "b+";
  if (Z == 0 && A == 1) return "n";
  if (Z == 1 && A == 1) return "p";
  if (Z == 1 && A == 2) return "d";
  if (Z == 1 && A == 3) return "t";
  if (Z == 2 && A == 3) return "h";
  if (Z == 2 && A == 4) return "a";
  return NuclideLabel(Z, A);
}

namespace {

std::string State(double e, const std::string &text) {
  if (e == 0.0) return "";
  if (!text.empty()) return "[" + text + "]";
  std::ostringstream s;
  s << e;
  return "[" + s.str() + "]";
}

// Digits after the decimal point of a number as written ("4.439" -> 3; "4" -> 0; "4.4e0" -> 1).
int Decimals(const std::string &text) {
  const std::size_t dot = text.find('.');
  if (dot == std::string::npos) return 0;
  const std::size_t e = text.find_first_of("eE", dot);
  return (int)((e == std::string::npos ? text.size() : e) - dot - 1);
}

// A written state [E_x] matches an excitation of <levels> to the digits written: [4.439] matches 4.43891 MeV.
bool StateMatches(bool excited, const Number &written, double e) {
  if (!excited || written.value == 0.0) return e == 0.0;
  const double half = 0.5 * std::pow(10.0, -Decimals(written.text));
  return std::fabs(e - written.value) <= half * (1.0 + 1e-9);
}

bool LightMatches(const Particle &p, const PairInfo &pair) {
  if (p.light == "g") return pair.kind == 10;
  if (p.light == "b-") return pair.kind == 20 && pair.lightZ < 0;
  if (p.light == "b+") return pair.kind == 20 && pair.lightZ > 0;
  return pair.kind == 0 && p.Z == pair.lightZ && p.A == pair.lightA &&
         StateMatches(p.excited, p.excitation, pair.lightExcitation);
}

}  // namespace

std::string PairLabel(const PairInfo &p) {
  return LightParticleLabel(p.kind, p.lightZ, p.lightA) + State(p.lightExcitation, p.lightExcitationText) + "+" +
         NuclideLabel(p.heavyZ, p.heavyA) + State(p.heavyExcitation, p.heavyExcitationText);
}

std::vector<std::string> PairListing(const std::vector<PairInfo> &pairs) {
  std::vector<std::string> labels;
  for (const PairInfo &p : pairs) labels.push_back(std::to_string(p.key) + ": " + PairLabel(p));
  return labels;
}

const PairInfo *PairByKey(int key, const std::vector<PairInfo> &pairs) {
  for (const PairInfo &p : pairs)
    if (p.key == key) return &p;
  return nullptr;
}

std::vector<int> PairsNamed(const Particle &light, const Particle &heavy, const ResidualState &state,
                            const std::vector<PairInfo> &pairs) {
  std::vector<int> keys;
  for (const PairInfo &p : pairs)
    if (LightMatches(light, p) && heavy.Z == p.heavyZ && heavy.A == p.heavyA &&
        StateMatches(state.kind == ResidualState::Excited, state.energy, p.heavyExcitation))
      keys.push_back(p.key);
  return keys;
}

std::string UnknownPairMessage(const std::string &what) {
  return what + " is no pair of <levels> (particles and states); the pairs are listed";
}

std::string AmbiguousPairMessage(const std::string &what, const std::vector<int> &keys, const std::string &remedy) {
  std::string list;
  for (int k : keys) list += (list.empty() ? "" : ", ") + std::to_string(k);
  return what + " names pairs " + list + " of <levels>; " + remedy;
}

}  // namespace experiment
}  // namespace azure2
