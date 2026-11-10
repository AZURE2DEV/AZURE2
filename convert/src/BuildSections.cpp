// The sections that are not records: <scatteringLengths> as model data, <thm>, <parameterSettings>, <config> and
// the prior-covariance copy (Builder.h).
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
  if (centreRows) Note(Construct::kSettings, "", std::to_string(centreRows) + " prior_centre row(s) folded into the norm and shift fields and removed");
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

}  // namespace convert
}  // namespace azure2
