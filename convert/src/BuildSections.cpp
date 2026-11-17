// The sections that are not records: <scatteringLengths> as model data, <thm>, <parameterSettings>, <config> and
// the prior-covariance copy (Builder.h).
#include <cctype>
#include <cmath>
#include "Builder.h"

namespace azure2 {
namespace convert {

// ================================================================ <scatteringLengths>
void Builder::BuildScatteringLengths() {
  for (const OldScatteringLength &l : old_.scatteringLengths) {
    Count(Construct::kScatteringLengths);
    const ScatteringLengthLine &s = l.line;
    if (s.active == 0) {
      Note(Construct::kScatteringLengths, "", "<scatteringLengths> line " + std::to_string(l.number) + " is inactive: dropped");
      continue;
    }
    const OldPair *p = Pair(s.pairKey);
    if (!p) {
      res_.refusal = "<scatteringLengths> line " + std::to_string(l.number) + ": no particle pair with key " + std::to_string(s.pairKey) + " (the old engine refused it)";
      return;
    }
    io::Record r;
    r.kind = "datum";
    r.name = Unique("scattering_length_" + std::to_string(l.number));
    r.comments.push_back("# <scatteringLengths> line " + std::to_string(l.number));
    Add(r, "quantity", "scattering_length");
    Add(r, "pair", PairLabel(*p));
    Add(r, "s", s.s < 0 ? "coherent" : SpinText(s.s));
    Add(r, "value", NumberText(s.value, l.text) + " fm");
    Add(r, "error", (s.sigma > 0 ? NumberText(s.sigma, l.text) : Ask(Construct::kScatteringLengths, r.name + ".error", "Re a is not constrained (sigma 0): format 3 has no such datum", "")) + " fm");
    if (s.hasImaginary) {
      Add(r, "imag", NumberText(s.imValue, l.text) + " fm");
      Add(r, "imag_error", (s.imSigma > 0 ? NumberText(s.imSigma, l.text) : Ask(Construct::kScatteringLengths, r.name + ".imag_error", "Im a is not constrained (sigma 0)", "")) + " fm");
    }
    out_.modelData.records.push_back(r);
  }
}

// ================================================================ the sections that stay
void Builder::RewriteThm() {
  const io::Document::Section *s = doc_.Find("thm");
  if (!s) return;
  out_.hasThm = true;
  std::vector<std::string> body;
  bool changed = false;
  for (const std::string &line : doc_.Body(*s)) {
    const std::size_t hash = line.find('#');
    const std::string code = hash == std::string::npos ? line : line.substr(0, hash);
    const std::size_t a = code.find_first_not_of(" \t");
    if (a != std::string::npos && code.compare(a, 11, "experiment[") == 0) {
      // drop the token segments=...
      std::string kept;
      std::istringstream tokens(code);
      std::string tok;
      bool dropped = false;
      while (tokens >> tok) {
        if (tok.compare(0, 9, "segments=") == 0) {
          dropped = true;
          continue;
        }
        kept += (kept.empty() ? "" : " ") + tok;
      }
      if (dropped) {
        Count(Construct::kThmSection);
        changed = true;
        body.push_back(kept + (hash == std::string::npos ? "" : "  " + line.substr(hash)));
        continue;
      }
    } else if (a != std::string::npos && code.compare(a, 6, "weight") == 0) {
      Count(Construct::kThmSection);
      changed = true;
      continue;  // now thm_weight on the data set or extrapolation
    }
    body.push_back(line);
  }
  if (changed) out_.bodies["thm"] = body;
}

void Builder::RewriteSettings() {
  for (const auto &p : res_.parameterNames) renamed_[p.first] = p.second;
  const io::Document::Section *s = doc_.Find("parameterSettings");
  if (!s) return;
  std::vector<std::string> body;
  int droppedRows = 0, centreRows = 0;
  for (const std::string &line : doc_.Body(*s)) {
    const std::size_t a = line.find_first_not_of(" \t");
    if (a == std::string::npos || line[a] == '#') {
      body.push_back(line);
      continue;
    }
    std::string trimmed = line.substr(a);
    while (!trimmed.empty() && (trimmed.back() == ' ' || trimmed.back() == '\t' || trimmed.back() == '\r')) trimmed.pop_back();
    if (ReadPriorCentreRow(trimmed).kind != PriorCentreRow::NotARow) {
      centreRows++;
      continue;  // folded into the norm and shift fields (F14)
    }
    std::istringstream tokens(trimmed);
    std::string first;
    tokens >> first;
    const auto it = renamed_.find(first);
    if (it != renamed_.end()) {
      body.push_back(line.substr(0, a) + it->second + line.substr(a + first.size()));
      continue;
    }
    if (first.compare(0, 8, "segment_") == 0 || first.compare(0, 8, "density_") == 0 || first.compare(0, 10, "screening_") == 0 ||
        first.compare(0, 8, "doppler_") == 0 || first.compare(0, 10, "abundance_") == 0) {
      droppedRows++;
      continue;
    }
    body.push_back(line);
  }
  Count(Construct::kSettings);
  if (centreRows)
    Note(Construct::kSettings, "",
         std::to_string(centreRows) + " prior_centre row(s) folded into the norm and shift fields and removed: each "
         "prior is centred on the row's value.  An AZURE2 without these rows (the dev branch, and this program before "
         "October 2026) ignored them and centred the prior on the line's start value; if the project was run with "
         "such a program, its priors move here");
  if (droppedRows)
    Note(Construct::kSettings, "", std::to_string(droppedRows) + " row(s) for parameters format 3 does not have (a fixed zero shift, an "
                                   "inactive or dropped line) removed");
  out_.bodies["parameterSettings"] = body;
}

void Builder::RewriteConfigAndPriors() {
  out_.configEdits["format"] = "3";
  if (!old_.abundanceSumSigmaText.empty()) {
    out_.configEdits["abundanceSumSigma"] = "";
    Note(Construct::kAbundances, "", "<config> abundanceSumSigma= removed (format 3 writes isotope_sum on the target)");
  }
  // useAdaptiveGrid: format 3 has one grid for every folded point, the shared mesh (the old value 1).
  if (!old_.adaptiveGridText.empty()) {
    out_.configEdits["useAdaptiveGrid"] = "";
    Count(Construct::kAdaptiveGrid);
    if (old_.adaptiveGridText != "1" && old_.adaptiveGridText != "true")
      Note(Construct::kAdaptiveGrid, "",
           "<config> useAdaptiveGrid=" + old_.adaptiveGridText + " removed: the uniform sub-point grid it chose no longer "
           "exists; every folded point is on the shared mesh, which is never coarser than the uniform grid of the same "
           "number of points, so a folded value changes by the uniform grid's quadrature error (over a narrow resonance "
           "by a large factor, to the converged value)");
  }
  // priorCovariance=f: a copy f.f3 with the parameters renamed; <mcmc> <parameters> likewise (renamed in place).
  const io::Document::Section *c = doc_.Find("config");
  std::string file;
  if (c)
    for (std::string line : doc_.Body(*c)) {
      const std::size_t hash = line.find('#');
      if (hash != std::string::npos) line = line.substr(0, hash);
      const std::size_t eq = line.find('=');
      if (eq == std::string::npos) continue;
      std::string key = line.substr(0, eq), value = line.substr(eq + 1);
      key.erase(0, key.find_first_not_of(" \t"));
      key.erase(key.find_last_not_of(" \t") + 1);
      value.erase(0, value.find_first_not_of(" \t"));
      value.erase(value.find_last_not_of(" \t\r") + 1);
      if (key == "priorCovariance") file = value;
    }
  auto renameTokens = [&](const std::string &line) {
    std::string out;
    std::size_t i = 0;
    while (i < line.size()) {
      if (line[i] == ' ' || line[i] == '\t') {
        out += line[i++];
        continue;
      }
      std::size_t j = i;
      while (j < line.size() && line[j] != ' ' && line[j] != '\t') j++;
      const std::string tok = line.substr(i, j - i);
      const auto it = renamed_.find(tok);
      out += it != renamed_.end() ? it->second : tok;
      i = j;
    }
    return out;
  };
  if (!file.empty()) {
    Count(Construct::kPriorCovariance);
    const std::string path = (file[0] == '/' || opt_.projectDir.empty()) ? file : opt_.projectDir + "/" + file;
    std::ifstream in(path, std::ios::binary);
    if (!in) {
      Note(Construct::kPriorCovariance, "", "the prior-covariance file " + file + " cannot be read: priorCovariance= left as it is");
    } else {
      std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
      std::string converted;
      std::size_t a = 0;
      while (a < all.size()) {
        std::size_t b = all.find('\n', a);
        if (b == std::string::npos) b = all.size();
        converted += renameTokens(all.substr(a, b - a));
        if (b < all.size()) converted += "\n";
        a = b + 1;
      }
      // <stem>.f3<extension>: prior_cov.txt -> prior_cov.f3.txt (the note's section 3.11)
      const std::size_t dot = file.find_last_of('.'), slash = file.find_last_of("/\\");
      const std::string side = (dot != std::string::npos && (slash == std::string::npos || dot > slash))
                                   ? file.substr(0, dot) + ".f3" + file.substr(dot)
                                   : file + ".f3";
      res_.sideFiles.push_back({side, converted});
      out_.configEdits["priorCovariance"] = side;
      Note(Construct::kPriorCovariance, "", "priorCovariance=" + side + ", a copy of " + file + " with the parameters renamed");
    }
  }
  if (const io::Document::Section *m = doc_.Find("mcmc"))
    if (const io::Document::Section *p = doc_.FindChild(*m, "parameters")) {
      std::vector<std::string> body;
      bool changed = false;
      for (const std::string &line : doc_.Body(*p)) {
        const std::string r = renameTokens(line);
        changed = changed || r != line;
        body.push_back(r);
      }
      if (changed) {
        Count(Construct::kMcmc);
        out_.bodies["mcmc/parameters"] = body;
      }
    }
}

// ================================================================ <potential>
// The old block -> one `potential` record per particle pair it gives the model to (PotentialRecords, which pyazr's
// writer uses too), coulomb_radius 3.5 fm and the solver's matching radius and step written out; the engine's warnings
// about the block and every section that writes no record are notes.  An unclosed block was refused by the old engine.
void Builder::BuildPotentials() {
  if (!old_.hasPotential) return;
  if (!old_.potential.closed) {
    res_.refusal = "<potential> is not closed (the old engine refused the project)";
    return;
  }
  for (const std::string &w : old_.potential.warnings) Note(Construct::kPotential, "", "the old engine: " + w);
  std::vector<experiment::PairInfo> pairs;
  for (const auto &p : old_.pairs) pairs.push_back(p.second.Info());
  const PotentialRecordsResult r = PotentialRecords(old_.potential.settings, pairs, used_);
  for (const io::Record &rec : r.records.records) {
    used_.insert(rec.name);
    Count(Construct::kPotential);
  }
  for (const std::string &n : r.notes) Note(Construct::kPotential, "", n);
  out_.potentials = r.records;
}

// ================================================================ <isospin>
// dev's block (2026-10-17) -> isospin_* records in the lines' order (the k of isospin_mix_<k> is the order of the mix
// records): reference -> isospin_reference, shift -> isospin_shift (its prior when the width is > 0), level ->
// isospin_level with the widths of its soft line, mix -> isospin_mix (a prior centred on 0 when the width is > 0),
// factor -> isospin_factor, pairT -> isospin_pair.  Pairs by their labels (pair_key where a label names two pairs);
// what dev's engine refused (an unknown pair, a level without lines) is refused by this build's validation.
void Builder::BuildIsospin() {
  auto integer = [](const std::string &f) { return (int)std::lround(std::stod(f)); };
  auto half = [](const std::string &f) {
    const long twice = std::lround(2.0 * std::stod(f));
    return twice % 2 == 0 ? std::to_string(twice / 2) : std::to_string(twice) + "/2";
  };
  auto refuse = [&](const OldProject::IsospinLine &c, const std::string &why) {
    res_.refusal = "<isospin> line \"" + c.text + "\" " + why;
  };
  auto pairFields = [&](io::Record &r, const std::string &field, const OldProject::IsospinLine &c, const std::string &key) {
    const OldPair *pair = Pair(integer(key));
    if (!pair) {
      refuse(c, "names pair key " + key + ", which <levels> does not have");
      return false;
    }
    Add(r, field, PairLabel(*pair));
    if (Ambiguous(*pair)) Add(r, "pair_key", std::to_string(pair->key));
    return true;
  };
  auto nuisance = [](const std::string &value, const std::string &unit, const std::string &vary, const std::string &centre,
                     const std::string &sigma) {
    std::string text = value + (unit.empty() ? "" : " " + unit) + (std::stod(vary) != 0.0 ? " fit" : " fixed");
    if (std::stod(sigma) > 0.0) text += " prior " + centre + " " + sigma;
    return text;
  };
  int mixes = 0;
  for (const OldProject::IsospinLine &c : old_.isospin) {
    const std::vector<std::string> &f = c.fields;
    io::Record r;
    if (c.kind == "soft") continue;  // with its level, below
    if (c.kind == "reference" || c.kind == "shift") {
      const std::string nuclide = experiment::NuclideLabel(integer(f[0]), integer(f[1]));
      r.kind = c.kind == "reference" ? "isospin_reference" : "isospin_shift";
      r.name = Unique((c.kind == "reference" ? "reference_" : "shift_") + nuclide);
      Add(r, "nucleus", nuclide);
      if (!pairFields(r, "threshold", c, f[2])) return;
      if (c.kind == "shift") {
        Add(r, "shift", nuisance(f[3], "MeV", f[4], f[5], f[6]));
        Param("isospin_shift_" + std::to_string(integer(f[0])) + "_" + std::to_string(integer(f[1])), r.name + ".shift");
      }
    } else if (c.kind == "level") {
      r.kind = "isospin_level";
      r.name = Unique("level_" + std::to_string(integer(f[0])));
      Add(r, "level", std::to_string(integer(f[0])));
      Add(r, "T", half(f[1]));
      for (const OldProject::IsospinLine &s : old_.isospin)
        if (s.kind == "soft" && integer(s.fields[0]) == integer(f[0])) {
          Add(r, "soft_energy", s.fields[1] + " MeV");
          Add(r, "soft_amplitude", s.fields[2]);
        }
    } else if (c.kind == "mix") {
      r.kind = "isospin_mix";
      r.name = Unique("mix_" + std::to_string(++mixes));
      Add(r, "level", std::to_string(integer(f[0])));
      if (!pairFields(r, "pair", c, f[1])) return;
      Add(r, "l", std::to_string(integer(f[2])));
      Add(r, "s", half(f[3]));
      Add(r, "amplitude", nuisance(f[4], "", f[5], "0", f[6]));
      Param("isospin_mix_" + std::to_string(mixes), r.name + ".amplitude");
    } else if (c.kind == "factor") {
      r.kind = "isospin_factor";
      r.name = Unique("factor_" + std::to_string(integer(f[0])) + "_" + std::to_string(integer(f[1])));
      Add(r, "level", std::to_string(integer(f[0])));
      if (!pairFields(r, "pair", c, f[1])) return;
      Add(r, "factor", f[2]);
    } else if (c.kind == "pairT") {
      r.kind = "isospin_pair";
      r.name = Unique("pair_" + std::to_string(integer(f[0])));
      if (!pairFields(r, "pair", c, f[0])) return;
      Add(r, "t_light", half(f[1]));
      Add(r, "t_heavy", half(f[2]));
    }
    out_.isospin.records.push_back(r);
    Count(Construct::kIsospin);
  }
  for (const OldProject::IsospinLine &c : old_.isospin) {
    if (c.kind != "soft") continue;
    bool declared = false;
    for (const OldProject::IsospinLine &l : old_.isospin)
      declared = declared || (l.kind == "level" && integer(l.fields[0]) == integer(c.fields[0]));
    if (!declared) {
      refuse(c, "names level ID " + c.fields[0] + ", which has no 'level' line");
      return;
    }
  }
}

// ================================================================ <channelRadii>
// dev's block (2026-11-10) -> one `channel_radius` record per line, in the lines' order: the pair by its label (with
// pair_key where the label names two pairs), l, the radius in fm as written, and s and J where the line gives them.
// A line on an unknown pair or on a photon or beta pair was refused by dev's engine, and is here.
void Builder::BuildChannelRadii() {
  for (const OldProject::ChannelRadiusLine &c : old_.channelRadii) {
    const OldPair *pair = Pair(c.key);
    if (!pair) {
      res_.refusal = "<channelRadii> line \"" + c.text + "\" names pair key " + std::to_string(c.key) + ", which <levels> does not have";
      return;
    }
    if (pair->pType != 0) {
      res_.refusal = "<channelRadii> line \"" + c.text + "\" names pair " + std::to_string(c.key) + ", a " +
                     (pair->pType == 10 ? "photon" : "beta") + " pair; only a particle channel has a radius";
      return;
    }
    io::Record r;
    r.kind = "channel_radius";
    // "p+12C" -> "p12C" (as the potential records' names): a format-3 name
    std::string base;
    for (char ch : PairLabel(*pair)) {
      if (ch == '+' || ch == ']') continue;
      base += std::isalnum((unsigned char)ch) ? ch : '_';
    }
    if (base.empty() || !std::isalpha((unsigned char)base[0])) base = "p" + base;
    base += "_l" + std::to_string(c.l);
    auto plain = [](std::string v) {
      for (char &ch : v)
        if (ch == '/') ch = '_';
      return v;
    };
    if (!c.s.empty()) base += "_s" + plain(c.s);
    if (!c.J.empty()) base += "_J" + plain(c.J);
    r.name = Unique(base);
    Add(r, "pair", PairLabel(*pair));
    if (Ambiguous(*pair)) Add(r, "pair_key", std::to_string(pair->key));
    Add(r, "l", std::to_string(c.l));
    Add(r, "radius", c.radius + " fm");
    if (!c.s.empty()) Add(r, "s", c.s);
    if (!c.J.empty()) Add(r, "J", c.J);
    out_.channelRadii.records.push_back(r);
    Count(Construct::kChannelRadii);
  }
}

}  // namespace convert
}  // namespace azure2
