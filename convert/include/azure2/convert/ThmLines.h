#ifndef AZURE2_CONVERT_THMLINES_H
#define AZURE2_CONVERT_THMLINES_H

#include <algorithm>
#include <cctype>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace azure2 {
namespace convert {

/*!
 * The two parts of the <thm> section that name data lines by their place (src/ThmExperiment.cpp), and that format 3
 * replaces by fields of the data sets (docs/upgrade/project_format3.md section 3.9): an experiment's
 * `segments=1,2,5-7` and the weight tables `weight[<k>]=<file>` / `weightTest[<k>]=<file>`.  The grammar of
 * ParseSegmentList and of the weight branch of ParseThmOptionLine, moved here so that the engine (until the format-3
 * reader, step 3b) and the converter to format 3 read them with one parser.
 */

/// Far above any project's number of segment lines (each is a file).
constexpr int kMaxThmSegment = 100000;

/// "1,2,5-7" -> {1,2,5,6,7}, sorted.  Returns "" or what is wrong.
inline std::string ParseThmSegmentList(const std::string &text, std::vector<int> &out) {
  auto wholeInt = [](const std::string &t, int &x) {
    if (t.empty()) return false;
    for (std::size_t i = 0; i < t.size(); i++)
      if (!std::isdigit((unsigned char)t[i]) && !(i == 0 && (t[i] == '+' || t[i] == '-'))) return false;
    std::istringstream s(t);
    return !!(s >> x);
  };
  std::vector<std::string> items;
  {
    std::string item;
    std::istringstream s(text);
    while (std::getline(s, item, ',')) items.push_back(item);
    if (!text.empty() && text.back() == ',') items.push_back("");
  }
  out.clear();
  std::set<int> seen;
  for (const std::string &item : items) {
    int lo, hi;
    std::size_t dash = item.find('-', 1);
    bool ok;
    if (dash == std::string::npos) {
      ok = wholeInt(item, lo);
      hi = lo;
    } else {
      ok = wholeInt(item.substr(0, dash), lo) && wholeInt(item.substr(dash + 1), hi) && hi >= lo;
    }
    if (!ok || lo < 1) return "segments='" + text + "': expected segment numbers >= 1 like 1,2,4-6";
    // A bound before the range is expanded: 1-100000000 would fill memory.
    if (hi > kMaxThmSegment) {
      std::ostringstream why;
      why << "segments='" << text << "': segment numbers go up to " << kMaxThmSegment;
      return why.str();
    }
    for (int k = lo; k <= hi; k++) {
      if (!seen.insert(k).second) {
        std::ostringstream why;
        why << "segments='" << text << "': segment " << k << " is listed twice";
        return why.str();
      }
      out.push_back(k);
    }
  }
  std::sort(out.begin(), out.end());
  return "";
}

/// A weight-table line of <thm>, key and value already split at '=' and trimmed.
struct ThmWeightKey {
  bool test = false;  ///< weightTest[k] (a <segmentsTest> line), else weight[k] (a <segmentsData> line)
  int segment = 0;    ///< k >= 1
  std::string file;   ///< as written
};

/// True when `key` starts with "weight" (the line is a weight line; `valid` then says whether it is well formed).
inline bool ReadThmWeightKey(const std::string &key, const std::string &value, ThmWeightKey &out, bool &valid) {
  if (key.compare(0, 6, "weight") != 0) return false;
  out.test = key.compare(0, 11, "weightTest[") == 0;
  const std::size_t open = out.test ? 10 : 6;
  valid = false;
  if (!(key.size() > open + 2 && key[open] == '[' && key.back() == ']' && !value.empty())) return true;
  std::istringstream ks(key.substr(open + 1, key.size() - open - 2));
  std::string rest;
  if (!(ks >> out.segment) || out.segment < 1 || (ks >> rest)) return true;
  out.file = value;
  valid = true;
  return true;
}

}  // namespace convert
}  // namespace azure2

#endif
