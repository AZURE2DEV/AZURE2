#ifndef AZURE2_CONVERT_BUILDER_H
#define AZURE2_CONVERT_BUILDER_H

// The builder of format-3 records from an old project (Build.h), shared by its source files: Build.cpp (reactions,
// phase-shift spins, the whole), BuildApparatus.cpp (<targetInt>, <abundances>, names of the shared records),
// BuildData.cpp (data sets, extrapolations, sample thicknesses), BuildSections.cpp (<scatteringLengths>, <thm>,
// <parameterSettings>, <config>, the prior-covariance copy).  Internal to the converter.
#include "Build.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include "azure2/convert/PriorCentreRow.h"
#include "azure2/io/Convert.h"
#include "azure2/physics/Elements.h"

namespace azure2 {
namespace convert {

inline io::Token Tok(const std::string &text) {
  io::Token t;
  t.text = text;
  return t;
}
inline io::Token Quoted(const std::string &text) {
  io::Token t = Tok(text);
  t.quoted = true;
  return t;
}
/// The blank-separated words of `text` as tokens.
inline std::vector<io::Token> Words(const std::string &text) {
  std::vector<io::Token> out;
  std::istringstream s(text);
  std::string w;
  while (s >> w) out.push_back(Tok(w));
  return out;
}
inline void Add(io::Record &r, const std::string &key, std::vector<io::Token> tokens, const std::string &comment = "") {
  io::Field f;
  f.key = key;
  f.tokens = std::move(tokens);
  if (!comment.empty()) f.comments.push_back("# CONVERTED: " + comment);
  r.fields.push_back(f);
}
inline void Add(io::Record &r, const std::string &key, const std::string &words, const std::string &comment = "") {
  Add(r, key, Words(words), comment);
}
/// The record's fields as one string (to find identical records).
inline std::string Signature(const io::Record &r) {
  std::string s = r.kind;
  for (const io::Field &f : r.fields) {
    s += "\n" + f.key;
    for (const io::Token &t : f.tokens) s += (t.quoted ? " \"" : " ") + t.text;
  }
  return s;
}
inline bool HasFreeParameter(const io::Record &r) {
  for (const io::Field &f : r.fields)
    for (const io::Token &t : f.tokens)
      if (t.text == "fit") return true;
  return false;
}

/// A name of format 3 ([A-Za-z][A-Za-z0-9_-]*) made from any text: other characters become '_', and a name that
/// would not start with a letter gets "d" in front.
inline std::string Sanitized(const std::string &text) {
  std::string s;
  for (char c : text) s += (std::isalnum((unsigned char)c) || c == '_' || c == '-') ? c : '_';
  if (s.empty() || !std::isalpha((unsigned char)s[0])) s = "d" + s;
  return s;
}
/// The stem of a data-file path: "data/s01_int_thin.dat" -> "s01_int_thin".
inline std::string Stem(const std::string &path) {
  std::string s = path;
  const std::size_t slash = s.find_last_of("/\\");
  if (slash != std::string::npos) s = s.substr(slash + 1);
  const std::size_t dot = s.find_last_of('.');
  if (dot != std::string::npos && dot > 0) s = s.substr(0, dot);
  return s;
}

inline std::string Abs(const std::string &numberText) {
  return (!numberText.empty() && numberText[0] == '-') ? numberText.substr(1) : numberText;
}

/// Number of whole numeric columns of the first data row of a data file (blank and '#' lines skipped), as DataLine
/// reads it; 0 when the file cannot be read or has no row.
inline int DataColumns(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) return 0;
  std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  std::string line;
  std::size_t a = 0;
  while (a <= all.size()) {
    std::size_t b = all.find_first_of("\r\n", a);
    if (b == std::string::npos) b = all.size();
    line = all.substr(a, b - a);
    a = b + 1;
    if (b < all.size() && all[b] == '\r' && b + 1 < all.size() && all[b + 1] == '\n') a = b + 2;
    const std::size_t k = line.find_first_not_of(" \t");
    if (k == std::string::npos || line[k] == '#') {
      if (b >= all.size()) break;
      continue;
    }
    std::istringstream ls(line);
    double v;
    int n = 0;
    while (ls >> v) n++;
    return n;
  }
  return 0;
}

struct Construct {
  static constexpr const char *kDataLine = "data line";
  static constexpr const char *kInactive = "active flag 0";
  static constexpr const char *kReaction = "entrance, exit keys";
  static constexpr const char *kEnergyWindow = "Emin Emax";
  static constexpr const char *kAngleWindow = "theta min, theta max";
  static constexpr const char *kDataCode = "data isDiff";
  static constexpr const char *kTestCode = "test isDiff";
  static constexpr const char *kThm = "isDiff + 10 (THM)";
  static constexpr const char *kSpin = "spin s";
  static constexpr const char *kNorm = "norm, varyNorm, normError";
  static constexpr const char *kShift = "shift, shiftError, vary";
  static constexpr const char *kFile = "file";
  static constexpr const char *kComposite = "composite";
  static constexpr const char *kUpos = "UPOS";
  static constexpr const char *kSqrtShift = "sqrtshift";
  static constexpr const char *kShiftFn = "shiftfn";
  static constexpr const char *kShiftGroup = "shiftgroup";
  static constexpr const char *kShiftFnGroup = "shiftfngroup";
  static constexpr const char *kNormGroup = "normgroup";
  static constexpr const char *kCovariance = "covariance";
  static constexpr const char *kSample = "sample / natural";
  static constexpr const char *kTotal = "total 1";
  static constexpr const char *kThickness = "transmission / yield n";
  static constexpr const char *kThicknessGroup = "thicknessgroup";
  static constexpr const char *kTestLine = "<segmentsTest> line";
  static constexpr const char *kEffectLine = "<targetInt> line";
  static constexpr const char *kNotASegment = "listed index that is not a segment";
  static constexpr const char *kThick = "density 0 with integration";
  static constexpr const char *kDensityFit = "densityfit";
  static constexpr const char *kDensityGroup = "densitygroup";
  static constexpr const char *kStraggling = "straggling";
  static constexpr const char *kSpread = "Gaussian spread";
  static constexpr const char *kProfile = "beamprofile";
  static constexpr const char *kQ = "Q coefficients";
  static constexpr const char *kUdr = "udr";
  static constexpr const char *kScreening = "screening";
  static constexpr const char *kDoppler = "doppler";
  static constexpr const char *kWindows = "fold windows";
  static constexpr const char *kAbundances = "<abundances>";
  static constexpr const char *kScatteringLengths = "<scatteringLengths>";
  static constexpr const char *kThmSection = "<thm> segments=, weight[k]";
  static constexpr const char *kSettings = "<parameterSettings>";
  static constexpr const char *kMcmc = "<mcmc> priors";
  static constexpr const char *kPriorCovariance = "priorCovariance";
};

inline const char *const kChargedElastic =
    "an angle-integrated cross section of the elastic scattering of a charged pair: the old engine computed the integral "
    "of the nuclear amplitude alone, pi/k^2 sum g |T|^2 without the Coulomb amplitude (GenMatrixFunc); format 3 refuses "
    "quantity sigma there (section 7.2) and has no quantity for this integral";

/// One active <targetInt> line as format-3 records (decision F7: target, beam and detector own their effects).
struct Apparatus {
  const OldEffectLine *line = nullptr;
  bool hasTarget = false, hasBeam = false, hasDetector = false;
  io::Record target, beam, detector;
  std::vector<io::Token> quadrature, foldWindows;
  bool folded = false;
  bool profileWindows = false;  ///< a beam profile: data columns 5-6 are the reconstructed-energy windows
  std::vector<int> dataKeys, testKeys;
  std::string targetName, beamName, detectorName;
};

/// One sample of <abundances>.
struct Sample {
  std::string name;  ///< "" for the unnamed sample
  std::map<std::pair<int, int>, const OldAbundance *> entries;  ///< compound (Z, A) -> its line (the last one)
  bool hasSum = false;
  double sum = 0;
  std::string sumText;
};

class Builder {
 public:
  Builder(const io::Document &doc, const OldProject &old, const ConvertOptions &options, Built &out, ConvertResult &result)
      : doc_(doc), old_(old), opt_(options), out_(out), res_(result) {}

  std::string Run();

 private:
  // ---------------------------------------------------------------- bookkeeping
  void Count(const std::string &construct) { res_.constructs[construct].converted++; }
  void Note(const std::string &construct, const std::string &record, const std::string &text) {
    res_.notes.push_back({construct, record, text});
    res_.constructs[construct].notes++;
  }
  /// An implicit meaning written with a note; with --strict a question whose proposal is `value` (section 8.2, F12).
  std::string Implicit(const std::string &construct, const std::string &record, const std::string &field,
                       const std::string &value, const std::string &text) {
    if (!opt_.strict) {
      Note(construct, record, text);
      pendingComment_ = text;
      return value;
    }
    return Ask(construct, record + "." + field, text, value);
  }
  /// A question (section 8.3): the answer from --answers, or '?' with --allow-questions, or the conversion fails.
  std::string Ask(const std::string &construct, const std::string &key, const std::string &text,
                  const std::string &proposal) {
    res_.questions.push_back({construct, key, text, proposal});
    res_.constructs[construct].questions++;
    const auto a = opt_.answers.find(key);
    if (a != opt_.answers.end()) return a->second;
    unanswered_++;
    return "?";
  }
  std::string Unique(const std::string &base) {
    std::string name = base;
    for (int k = 2; used_.count(name); k++) name = base + "_" + std::to_string(k);
    used_.insert(name);
    return name;
  }
  void Param(const std::string &oldName, const std::string &newName) { res_.parameterNames.push_back({oldName, newName}); }

  // ---------------------------------------------------------------- labels
  const OldPair *Pair(int key) const {
    const auto it = old_.pairs.find(key);
    return it == old_.pairs.end() ? nullptr : &it->second;
  }
  static std::string Nuclide(int Z, int A) { return std::to_string(A) + physics::ElementSymbol(Z); }
  static std::string LightLabel(const OldPair &p) {
    if (p.pType == 10) return "g";
    if (p.pType == 20) return p.z1 < 0 ? "b-" : "b+";
    const int A = p.A1();
    if (p.z1 == 0 && A == 1) return "n";
    if (p.z1 == 1 && A == 1) return "p";
    if (p.z1 == 1 && A == 2) return "d";
    if (p.z1 == 1 && A == 3) return "t";
    if (p.z1 == 2 && A == 3) return "h";
    if (p.z1 == 2 && A == 4) return "a";
    return Nuclide(p.z1, A);
  }
  /// Other pairs with the same particles and states: the label would not identify `p` (section 8.3).
  bool Ambiguous(const OldPair &p) const {
    for (const auto &kv : old_.pairs) {
      const OldPair &o = kv.second;
      if (o.key != p.key && o.pType == p.pType && o.z1 == p.z1 && o.A1() == p.A1() && o.z2 == p.z2 && o.A2() == p.A2() &&
          o.e1 == p.e1 && o.e2 == p.e2)
        return true;
    }
    return false;
  }
  /// The label of one reaction; exit -1 is [total].  "" with `why` set when a pair is unknown.
  std::string Reaction(int entrance, int exit, std::string &why, std::string &ambiguity) const {
    const OldPair *in = Pair(entrance);
    if (!in) {
      why = "entrance pair " + std::to_string(entrance) + " is not a pair of <levels>";
      return "";
    }
    const OldPair *out = exit == -1 ? nullptr : Pair(exit);
    if (exit != -1 && !out) {
      why = "exit pair " + std::to_string(exit) + " is not a pair of <levels>";
      return "";
    }
    if (in->e1 != 0 || in->e2 != 0) ambiguity += "entrance pair " + std::to_string(entrance) + " has an excited particle; ";
    if (Ambiguous(*in)) ambiguity += "pair " + std::to_string(entrance) + " has the same particles and states as another pair; ";
    if (out && Ambiguous(*out)) ambiguity += "pair " + std::to_string(exit) + " has the same particles and states as another pair; ";
    std::string s = Nuclide(in->z2, in->A2()) + "(" + LightLabel(*in) + ",";
    if (!out) {
      s += "g)" + Nuclide(in->z1 + in->z2, in->A1() + in->A2()) + "[total]";
      return s;
    }
    s += LightLabel(*out);
    if (out->e1 != 0) s += "[" + out->e1Text + "]";
    s += ")" + Nuclide(out->z2, out->A2());
    if (out->e2 != 0) s += "[" + out->e2Text + "]";
    return s;
  }
  static std::string PairLabel(const OldPair &p) {
    return LightLabel(p) + "+" + Nuclide(p.z2, p.A2()) + (p.e2 != 0 ? "[" + p.e2Text + "]" : "");
  }

  /// The channel spin of a phase-shift line without "spin s": the one spin of (J, l) of the entrance pair in the
  /// level scheme (EData::Initialize refuses a line with several).  "" with res_.refusal set when there are several.
  std::string PhaseSpin(const std::string &record, int entrance, double J, int l);
  /// True when a term of the line is the elastic scattering of a charged pair, angle-integrated: the old engine
  /// computed the integral of the nuclear amplitude alone (pi/k^2 sum g |T|^2, GenMatrixFunc), which format 3 has no
  /// quantity for (section 7.2 refuses sigma there).
  bool ChargedElasticTerm(int entrance, int exit, const std::vector<ComponentEntry> &components) const {
    auto charged = [&](int in, int out) {
      const OldPair *p = Pair(in);
      return in == out && p && p->pType == 0 && p->z1 > 0 && p->z2 > 0;
    };
    if (charged(entrance, exit)) return true;
    for (const ComponentEntry &c : components)
      if (charged(c.entrance, c.exit)) return true;
    return false;
  }
  /// The CONVERTED comment of the last Implicit call ("" when there was none), taken by the next field written.
  std::string TakeComment() {
    std::string c;
    c.swap(pendingComment_);
    return c;
  }

  /// The target of a transmission / yield sample thickness (section 8.1): kind sample, the areal density in atoms/b.
  std::string SampleThicknessTarget(const std::string &record, int entrance, double n, const std::string &line, bool fit,
                                    int vary, double sigma, const std::string &group, int key, bool isData);

  // ---------------------------------------------------------------- parts
  void BuildApparatus();
  void BuildSamples();
  void BuildData();
  void BuildExtrapolations();
  void NameApparatus();
  void BuildScatteringLengths();
  void RewriteThm();
  void RewriteSettings();
  void RewriteConfigAndPriors();

  /// The reaction field of a line (base pair keys, components, operation, frame of fixed angles).
  bool ReactionTokens(const std::string &record, int entrance, int exit, int isAdvanced, int op,
                      const std::vector<ComponentEntry> &components, const std::string &angleFrame, std::vector<io::Token> &out,
                      std::string &why);
  /// A nuisance: "<v> [unit] fixed|fit [prior c w]".
  static std::vector<io::Token> Nuisance(const std::string &value, const std::string &unit, bool fit, bool prior,
                                         const std::string &centre, const std::string &width) {
    std::string s = value + (unit.empty() ? "" : " " + unit) + (fit ? " fit" : " fixed");
    if (prior) s += " prior " + centre + " " + width;
    return Words(s);
  }

  const io::Document &doc_;
  const OldProject &old_;
  const ConvertOptions &opt_;
  Built &out_;
  ConvertResult &res_;
  std::set<std::string> used_;
  int unanswered_ = 0;
  std::string pendingComment_;

  std::map<int, std::string> dataName_, testName_;  // key -> record name
  std::map<int, const OldDataLine *> activeData_;   // key -> active data line
  std::map<int, const OldTestLine *> activeTest_;
  std::vector<Apparatus> apparatus_;
  std::map<int, int> apparatusOfData_, apparatusOfTest_;  // key -> index into apparatus_
  std::vector<Sample> samples_;
  std::map<int, int> sampleOfData_, sampleOfTest_;  // key -> index into samples_
  std::map<int, std::string> sampleTargetOfData_, sampleTargetOfTest_;  // key -> target name for a sample without <targetInt>
  std::vector<io::Record> sampleTargets_;
  std::vector<io::Record> groupRecords_;  // normalization and energy_scale
  std::map<std::string, std::string> thicknessTarget_;  // thicknessgroup label -> target name
  std::map<std::string, int> thicknessLeader_;          // label -> leader key
  std::map<std::string, std::string> renamed_;  // old parameter name -> new
};

}  // namespace convert
}  // namespace azure2

#endif
