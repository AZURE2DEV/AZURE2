// The records of <isospin> through Schema::Build (experiment/src/Isospin.cpp): the six kinds resolved against the
// pairs and level IDs of <levels> into BuiltProject::isospin, and one case per refusal (codes E_ISOSPIN_RECORD,
// E_ISOSPIN_PAIR; check 24).  Beside build_test.cpp, which is at the size the ratchet allows.
//   experiment_isospin_test
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "azure2/experiment/Schema.h"
#include "azure2/io/Document.h"
#include "azure2/io/RecordReader.h"

using namespace azure2::experiment;
using azure2::io::Document;

namespace {

int failures = 0;
int checks = 0;

void check(const std::string &what, bool ok) {
  checks++;
  if (!ok) {
    std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

// p+15N (1, 16O), g+16O (4), p+12C (7, 13N): the pairs the cases name; one level, ID 5, in 13N
LevelTables Tables() {
  auto pair = [](int key, int type, int z1, int a1, int z2, int a2, int twoJ1, int twoJ2) {
    PairInfo p;
    p.key = key;
    p.kind = type;
    p.lightZ = z1;
    p.lightA = a1;
    p.heavyZ = z2;
    p.heavyA = a2;
    p.twoJLight = twoJ1;
    p.twoJHeavy = twoJ2;
    return p;
  };
  LevelTables t;
  t.pairs = {pair(1, 0, 1, 1, 7, 15, 1, 1), pair(4, 10, 0, 0, 8, 16, 2, 0), pair(7, 0, 1, 1, 6, 12, 1, 0)};
  LevelInfo level;
  level.ids = {5};
  level.twoJ = 1;
  level.parity = 1;
  level.group = level.level = 1;
  level.nucleusZ = 7;
  level.nucleusA = 13;
  t.levels = {level};
  return t;
}

/// The records `body` of <isospin>, validated and built against Tables().
BuiltProject build(const std::string &body) {
  const Document doc = Document::Parse("<isospin>\n" + body + "</isospin>\n");
  std::vector<SectionRecords> sections;
  for (const Document::Section &s : doc.Sections())
    sections.push_back({s.name, azure2::io::ReadRecords(doc, s, Schema::KindsOfSection(s.name)).set.records});
  const ValidationResult v = Schema::Validate(sections);
  BuiltProject b;
  if (!v.Ok()) {
    b.diagnostics = v.diagnostics;
    for (Diagnostic &d : b.diagnostics) d.message = "(validation) " + d.message;
    return b;
  }
  return Schema::Build(v, Tables(), BuildOptions());
}

/// `base` with the line of field `key` replaced by `line`.
std::string Set(const std::string &base, const std::string &key, const std::string &line) {
  std::string out;
  std::size_t a = 0;
  bool done = false;
  while (a < base.size()) {
    const std::size_t b = base.find('\n', a);
    const std::string l = base.substr(a, b - a);
    if (!done && l.compare(0, 4 + key.size() + 1, "    " + key + " ") == 0) {
      out += line + "\n";
      done = true;
    } else {
      out += l + "\n";
    }
    a = b + 1;
  }
  return out;
}

void ExpectOne(const std::string &code, const BuiltProject &b) {
  const bool ok = b.diagnostics.size() == 1 && b.diagnostics[0].code == code;
  check(code + " alone", ok);
  if (!ok)
    for (const Diagnostic &d : b.diagnostics) std::printf("        got %s: %s\n", d.code.c_str(), d.message.c_str());
}

}  // namespace

int main() {
  const std::string ref = "isospin_reference ref\n    nucleus 13N\n    threshold p+12C\n";
  const std::string shift = "isospin_shift sh\n    nucleus 16O\n    threshold p+15N\n    shift -500 keV fit prior -500 20%\n";
  const std::string lev = "isospin_level lv\n    level 5\n    T 1/2\n";
  const std::string mix = "isospin_mix mx\n    level 5\n    pair p+12C\n    l 1\n    s 1/2\n    amplitude 0.05 fit prior 0 0.1\n";
  const std::string factor = "isospin_factor fc\n    level 5\n    pair p+12C\n    factor 0.5\n";
  const std::string pairT = "isospin_pair pt\n    pair p+12C\n    t_light 1/2\n    t_heavy 0\n";
  const BuiltProject r = build(ref + shift + lev + mix + factor + pairT);
  const bool all = r.Ok() && r.diagnostics.empty() && r.isospin.nuclei.size() == 2 && r.isospin.levels.size() == 1 &&
                   r.isospin.mixes.size() == 1 && r.isospin.factors.size() == 1 && r.isospin.pairs.size() == 1;
  check("the six kinds of <isospin> build", all);
  if (all) {
    const IsospinSpec &iso = r.isospin;
    check("the reference: nucleus, threshold pair key", iso.nuclei[0].reference && iso.nuclei[0].Z == 7 && iso.nuclei[0].A == 13 &&
                                                            iso.nuclei[0].pairKey == 7);
    check("the shift in MeV, its prior through the unit table and the percentage",
          !iso.nuclei[1].reference && iso.nuclei[1].pairKey == 1 && std::fabs(iso.nuclei[1].shift.value + 0.5) < 1e-12 && iso.nuclei[1].shift.fit &&
              iso.nuclei[1].shift.hasPrior && std::fabs(iso.nuclei[1].shift.priorCentre + 0.5) < 1e-12 &&
              std::fabs(iso.nuclei[1].shift.priorWidth - 0.1) < 1e-12);
    check("the level's ID and 2T; Declares", iso.levels[0].levelId == 5 && iso.levels[0].twoT == 1 && !iso.levels[0].soft &&
                                                 iso.Declares(5) && !iso.Declares(6));
    check("the mix: pair key, l, 2s, amplitude and prior", iso.mixes[0].pairKey == 7 && iso.mixes[0].l == 1 && iso.mixes[0].twoS == 1 &&
                                                               iso.mixes[0].amplitude.value == 0.05 && iso.mixes[0].amplitude.priorWidth == 0.1);
    check("the factor and the pair's isospins", iso.factors[0].value == 0.5 && iso.pairs[0].twoT1 == 1 && iso.pairs[0].twoT2 == 0);
  }
  const BuiltProject soft = build(ref + lev + "    soft_energy 50 keV\n    soft_amplitude 0.02\n");
  check("a soft tie: both widths, the energy width in MeV",
        soft.Ok() && soft.isospin.levels.size() == 1 && soft.isospin.levels[0].soft &&
            std::fabs(soft.isospin.levels[0].sigmaE - 0.05) < 1e-12 && soft.isospin.levels[0].sigmaGamma == 0.02);
  ExpectOne("E_ISOSPIN_RECORD", build(ref + lev + "    soft_energy 50 keV\n"));                       // one soft width
  ExpectOne("E_ISOSPIN_RECORD", build(ref + lev + "isospin_level again\n    level 5\n    T 1/2\n"));  // a level twice
  ExpectOne("E_ISOSPIN_RECORD", build(ref + Set(lev, "level", "    level 6")));                         // no such level
  ExpectOne("E_ISOSPIN_RECORD", build(lev));                                                            // no reference
  ExpectOne("E_ISOSPIN_RECORD", build(ref + "isospin_reference two\n    nucleus 16O\n    threshold p+15N\n"));
  ExpectOne("E_ISOSPIN_RECORD",                                                                         // a nucleus twice
            build(ref + "isospin_shift again\n    nucleus 13N\n    threshold p+12C\n    shift 0 MeV fixed\n"));
  ExpectOne("E_ISOSPIN_RECORD", build(ref + mix));                                                      // a mix without its level
  ExpectOne("E_ISOSPIN_RECORD", build(ref + lev + "    soft_energy 50 keV\n    soft_amplitude 0.02\n" + mix));  // soft and mix
  ExpectOne("E_ISOSPIN_RECORD", build(ref + Set(shift, "shift", "    shift -0.5 MeV fit prior 0 10%")));  // a prior without a width
  ExpectOne("E_ISOSPIN_PAIR", build(Set(ref, "threshold", "    threshold p+15N")));                     // not a pair of 13N
  ExpectOne("E_ISOSPIN_PAIR", build(ref + Set(pairT, "pair", "    pair g+16O")));                       // a photon pair
  ExpectOne("E_UNKNOWN_PAIR", build(Set(ref, "threshold", "    threshold p+13C")));
  std::printf("experiment_isospin: %d checks, %d failed\n", checks, failures);
  return failures == 0 ? 0 : 1;
}
