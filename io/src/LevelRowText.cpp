// The text of a <levels> row: its tokens by column, and the row with tokens replaced (LevelRowText.h).
#include "azure2/io/LevelRowText.h"

#include <algorithm>

namespace azure2 {
namespace io {

namespace {

constexpr int kMandatory = 31;
constexpr std::size_t npos = std::string::npos;

bool Blank(char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f' || c == '\v'; }

std::vector<std::pair<std::size_t, std::size_t>> TokenSpans(const std::string &line) {
  std::vector<std::pair<std::size_t, std::size_t>> spans;
  std::size_t i = 0;
  while (i < line.size()) {
    while (i < line.size() && Blank(line[i])) i++;
    if (i >= line.size()) break;
    const std::size_t a = i;
    while (i < line.size() && !Blank(line[i])) i++;
    spans.push_back({a, i});
  }
  return spans;
}

}  // namespace

const std::vector<std::string> &LevelColumnNames() {
  static const std::vector<std::string> names = {
      "levelJ", "levelPi", "levelE", "levelFix", "aa", "ir", "s2", "l2", "levelID", "isActive", "channelFix",
      "gamma", "j1", "pi1", "j2", "pi2", "e2", "m1", "m2", "z1", "z2", "entranceSepE", "sepE", "j3", "pi3", "e3",
      "pType", "chRad", "g1", "g2", "ecMultMask", "bindingE", "gammaIsRWA", "e1"};
  return names;
}

int LevelColumnIndex(const std::string &name) {
  const std::vector<std::string> &n = LevelColumnNames();
  const auto it = std::find(n.begin(), n.end(), name);
  return it == n.end() ? -1 : (int)(it - n.begin());
}

std::string LevelRowText::Token(int column) const {
  if (column < 0 || column >= (int)spans.size() || spans[column].first == npos) return std::string();
  return line.substr(spans[column].first, spans[column].second - spans[column].first);
}

bool SplitLevelRow(const std::string &line, LevelRowText &out, std::string &why) {
  const std::vector<std::pair<std::size_t, std::size_t>> t = TokenSpans(line);
  if (t.size() < (std::size_t)kMandatory || t.size() > (std::size_t)kMandatory + 3) {
    why = "a <levels> row has 31 numbers, then 1 (e1), 2 (bindingE gammaIsRWA) or 3 (both and e1); this line has " +
          std::to_string(t.size());
    return false;
  }
  out = LevelRowText();
  out.line = line;
  const int extra = (int)t.size() - kMandatory;
  out.hasThm = extra >= 2;
  out.hasE1 = extra == 1 || extra == 3;
  out.spans.assign(kMandatory + 3, {npos, npos});
  for (int k = 0; k < kMandatory; k++) out.spans[k] = t[k];
  if (out.hasThm) {
    out.spans[kMandatory] = t[kMandatory];
    out.spans[kMandatory + 1] = t[kMandatory + 1];
  }
  if (out.hasE1) out.spans[kMandatory + 2] = t.back();
  return true;
}

std::string EditLevelRow(const LevelRowText &row, const std::map<int, std::string> &tokens) {
  bool thm = row.hasThm, e1 = row.hasE1;
  for (const auto &kv : tokens) {
    if (kv.first == kMandatory || kv.first == kMandatory + 1) thm = true;
    if (kv.first == kMandatory + 2) e1 = true;
  }
  // The tokens written in place, right to left so that earlier spans stay valid.
  std::string out = row.line;
  std::vector<std::pair<int, std::string>> inPlace;
  for (const auto &kv : tokens)
    if (kv.first >= 0 && kv.first < (int)row.spans.size() && row.spans[kv.first].first != npos)
      inPlace.push_back(kv);
  std::sort(inPlace.begin(), inPlace.end(), [](const auto &a, const auto &b) { return a.first > b.first; });
  if (thm == row.hasThm && e1 == row.hasE1) {
    for (const auto &kv : inPlace) {
      const std::size_t a = row.spans[kv.first].first, b = row.spans[kv.first].second, width = b - a;
      const std::string text = kv.second.size() < width ? std::string(width - kv.second.size(), ' ') + kv.second : kv.second;
      out.replace(a, width, text);
    }
    return out;
  }
  // The optional group changes: the mandatory tokens in place, then the optional ones written after the 31st.
  const std::size_t end31 = row.spans[kMandatory - 1].second;
  std::size_t endAll = end31;
  for (const auto &s : row.spans)
    if (s.first != npos) endAll = std::max(endAll, s.second);
  std::string head = row.line.substr(0, end31);
  for (const auto &kv : inPlace) {
    if (kv.first >= kMandatory) continue;
    const std::size_t a = row.spans[kv.first].first, b = row.spans[kv.first].second, width = b - a;
    const std::string text = kv.second.size() < width ? std::string(width - kv.second.size(), ' ') + kv.second : kv.second;
    head.replace(a, width, text);
  }
  auto value = [&](int column) {
    const auto it = tokens.find(column);
    if (it != tokens.end()) return it->second;
    const std::string t = row.Token(column);
    return t.empty() ? std::string("0") : t;
  };
  std::vector<std::string> tail;
  if (thm) {
    tail.push_back(value(kMandatory));
    tail.push_back(value(kMandatory + 1));
  }
  if (e1) tail.push_back(value(kMandatory + 2));
  const std::string sep = thm ? "  " : " ";
  for (const std::string &t : tail) head += sep + t;
  return head + row.line.substr(endAll);
}

}  // namespace io
}  // namespace azure2
