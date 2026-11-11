// One line of the <thm> section as far as the names go (ThmSectionLines.h).
#include "azure2/experiment/ThmSectionLines.h"

namespace azure2 {
namespace experiment {

namespace {

std::string Trim(const std::string &s) {
  const std::size_t a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return std::string();
  return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}

}  // namespace

ThmSectionLine ReadThmSectionLine(const std::string &line) {
  ThmSectionLine out;
  out.text = Trim(line.substr(0, line.find('#')));
  const std::string open = "experiment[";
  if (out.text.compare(0, open.size(), open) == 0) {
    out.experiment = true;
    const std::size_t close = out.text.find(']');
    if (close == std::string::npos) return out;
    out.name = out.text.substr(open.size(), close - open.size());
    for (std::size_t i = close + 1; i < out.text.size();) {
      while (i < out.text.size() && (out.text[i] == ' ' || out.text[i] == '\t')) i++;
      const std::size_t a = i;
      while (i < out.text.size() && out.text[i] != ' ' && out.text[i] != '\t') i++;
      if (i > a) out.keys.push_back(out.text.substr(a, std::min(out.text.find('=', a), i) - a));
    }
    return out;
  }
  const std::size_t eq = out.text.find('=');
  if (eq != std::string::npos) out.keys.push_back(Trim(out.text.substr(0, eq)));
  return out;
}

std::string ThmOldListKey(const ThmSectionLine &line, std::string &field) {
  if (line.experiment) {
    for (const std::string &k : line.keys)
      if (k == "segments") {
        field = "thm " + line.name;
        return k;
      }
    return std::string();
  }
  if (line.keys.empty()) return std::string();
  const std::string &k = line.keys[0];
  if (k.compare(0, 7, "weight[") == 0 || k.compare(0, 11, "weightTest[") == 0) {
    field = "thm_weight <file>";
    return k;
  }
  return std::string();
}

}  // namespace experiment
}  // namespace azure2
