#ifndef THM_TEXT_H
#define THM_TEXT_H

#include <cctype>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

/*!
 * The small readers of the <thm> grammar's values, shared by the key table (ThmKeys.cpp) and the readers of nuclides
 * and the coherent background (ThmExperiment.cpp).
 */
namespace thm_text {

// A whole token as a number (operator>> that must consume everything).
inline bool ReadWholeDouble(const std::string &text, double &x) {
  std::istringstream s(text);
  std::string rest;
  return !!(s >> x) && !(s >> rest) && std::isfinite(x);
}
inline bool ReadWholeInt(const std::string &text, int &x) {
  if (text.empty()) return false;
  for (size_t i = 0; i < text.size(); i++)
    if (!std::isdigit((unsigned char)text[i]) && !(i == 0 && (text[i] == '+' || text[i] == '-'))) return false;
  std::istringstream s(text);
  return !!(s >> x);
}

inline std::vector<std::string> Split(const std::string &text, char sep) {
  std::vector<std::string> out;
  std::string item;
  std::istringstream s(text);
  while (std::getline(s, item, sep)) out.push_back(item);
  if (!text.empty() && text.back() == sep) out.push_back("");
  return out;
}

// "pmin-pmax" (MeV/c, 0 <= pmin <= pmax): the '-' whose both sides are numbers.
inline bool ReadWindow(const std::string &text, double &lo, double &hi) {
  for (size_t k = 1; k + 1 < text.size(); k++) {
    if (text[k] != '-') continue;
    if (ReadWholeDouble(text.substr(0, k), lo) && ReadWholeDouble(text.substr(k + 1), hi))
      return lo >= 0.0 && hi >= lo;
  }
  return false;
}

}  // namespace thm_text

#endif
