// The format-3 records of an old project (Build.h): docs/upgrade/project_format3.md section 8.1, construct by
// construct, with the notes of section 8.2 and the questions of section 8.3.  Every value is the one the old engine
// used: numbers are written as the old line wrote them (NumberText), defaults the old engine applied silently are
// written out, and what it ignored is dropped with a note.  This file: reaction labels, phase-shift spins, the order of
// the parts (Builder::Run).
#include <set>

#include "Builder.h"

namespace azure2 {
namespace convert {

// ================================================================ reactions
bool Builder::ReactionTokens(const std::string &record, int entrance, int exit, int isAdvanced, int op,
                             const std::vector<ComponentEntry> &components, const std::string &angleFrame,
                             std::vector<io::Token> &out, std::string &why,
                             std::vector<std::pair<std::string, std::string>> &pairKeys) {
  bool inAmbiguous = false, outAmbiguous = false;
  const std::string base = Reaction(entrance, exit, why, inAmbiguous, outAmbiguous);
  if (base.empty()) return false;
  out.push_back(Tok(base));
  // The pair keys of the terms whose label names two pairs, per side (the record's entrance_pair and exit_pair).
  std::set<int> inKeys, outKeys;
  if (inAmbiguous) inKeys.insert(entrance);
  if (outAmbiguous) outKeys.insert(exit);
  if (isAdvanced == 1) {
    Count(Construct::kComposite);
    if (components.empty()) {
      Note(Construct::kComposite, record, "the composite block lists fewer numbers than its components need; the old "
                                          "engine read the line as a plain segment (a warning), and so does this record");
    } else if (op == 1 && components.size() != 1) {
      why = "a ratio with " + std::to_string(components.size()) + " components (format 3 writes one numerator and one denominator)";
      return false;
    }
    for (const ComponentEntry &c : components) {
      bool a = false, b = false;
      std::string term = Reaction(c.entrance, c.exit, why, a, b);
      if (term.empty()) return false;
      if (a) inKeys.insert(c.entrance);
      if (b) outKeys.insert(c.exit);
      if (c.hasScaling) term = ShortestText(c.scaling) + "*" + term;
      out.push_back(Tok(op == 1 ? "/" : "+"));
      out.push_back(Tok(term));
      if (c.hasAngle) {
        out.push_back(Tok("@"));
        out.push_back(Tok(ShortestText(c.angle)));
        out.push_back(Tok("deg"));
        out.push_back(Tok(angleFrame));
      }
    }
  }
  // A key applies to every term whose label names two pairs on its side (section 3.5): the old line's terms carry
  // their own keys, so a line whose such terms take different pairs of the same particles has no format-3 form.
  for (const auto &side : {std::make_pair("entrance", &inKeys), std::make_pair("exit", &outKeys)}) {
    if (side.second->size() < 2) continue;
    std::string keys;
    for (int k : *side.second) keys += (keys.empty() ? "" : " and ") + std::to_string(k);
    // The old engine computed the line: refused, never dropped (BuildData stops at res_.refusal).
    res_.refusal = record + ": its terms take the " + std::string(side.first) + " pairs " + keys +
          ", which have the same particles and states; format 3 writes one " + side.first +
          "_pair for every term of a sum or ratio whose label names two pairs, so this line has no format-3 form. "
          "Split it into one data set per pair, or give the pairs different particles or states in <levels>";
    return false;
  }
  if (!inKeys.empty()) pairKeys.push_back({"entrance_pair", std::to_string(*inKeys.begin())});
  if (!outKeys.empty()) pairKeys.push_back({"exit_pair", std::to_string(*outKeys.begin())});
  Count(Construct::kReaction);
  return true;
}

// ================================================================ phase shifts
std::string Builder::PhaseSpin(const std::string &record, int entrance, double J, int l) {
  std::vector<double> spins;
  auto add = [&](double sv) {
    for (double v : spins)
      if (std::fabs(v - sv) < 1e-6) return;
    spins.push_back(sv);
  };
  const OldPair *p = Pair(entrance);
  if (p && p->pType == 0) {
    for (const OldChannel &c : old_.channels)
      if (c.active && c.pairKey == entrance && std::fabs(c.J - J) < 1e-6 && c.l == l) add(c.s);
    if (old_.includeAllJPi)  // hard-sphere channels of every (J, l, s) the pair allows
      for (double sv = std::fabs(p->j1 - p->j2); sv <= p->j1 + p->j2 + 1e-6; sv += 1.0)
        if (std::fabs(l - sv) <= J + 1e-6 && J <= l + sv + 1e-6) add(sv);
  }
  if (spins.size() == 1) return SpinText(spins[0]);
  if (spins.size() > 1) {
    res_.refusal = record + ": a phase shift without 'spin s' whose (J, l) has several channel spins, which the old engine refused";
    return "";
  }
  return Ask(Construct::kSpin, record + ".s", "the phase shift's channel spin: the line gives none and the level scheme has no channel of its (J, l)", "");
}

// ================================================================ the whole
std::string Builder::Run() {
  // Names: data sets from their file names (section 8.1), extrapolations x<k>; then the rest.
  for (const OldDataLine &d : old_.data) {
    dataName_[d.key] = Unique(Sanitized(Stem(d.parsed->dataFile())));
    if (d.parsed->isActive() == 1) activeData_[d.key] = &d;
  }
  for (const OldTestLine &t : old_.tests) {
    testName_[t.key] = Unique("x" + std::to_string(t.key));
    if (t.parsed->isActive() == 1) activeTest_[t.key] = &t;
  }
  BuildApparatus();
  BuildSamples();
  MoveDopplerToSamples();
  NameApparatus();
  if (!res_.refusal.empty()) return res_.refusal;
  BuildData();
  if (!res_.refusal.empty()) return res_.refusal;
  BuildExtrapolations();
  if (!res_.refusal.empty()) return res_.refusal;
  BuildScatteringLengths();
  if (!res_.refusal.empty()) return res_.refusal;
  BuildPotentials();
  if (!res_.refusal.empty()) return res_.refusal;
  BuildChannelRadii();
  if (!res_.refusal.empty()) return res_.refusal;
  BuildIsospin();
  if (!res_.refusal.empty()) return res_.refusal;
  // group records before the data sets that name them
  out_.data.records.insert(out_.data.records.begin(), groupRecords_.begin(), groupRecords_.end());
  RewriteThm();
  RewriteSettings();
  RewriteConfigAndPriors();
  if (unanswered_ > 0 && !opt_.allowQuestions)
    return std::to_string(unanswered_) + " question(s) need an answer (--answers), or --allow-questions to write '?'";
  return "";
}

std::string BuildFormat3(const io::Document &doc, const OldProject &old, const ConvertOptions &options, Built &out,
                         ConvertResult &result) {
  Builder b(doc, old, options, out, result);
  return b.Run();
}

}  // namespace convert
}  // namespace azure2
