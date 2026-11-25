// A line of an old project that the old grammar cannot read: the diagnosis and the readings proposed (LineRepair.h).
#include "azure2/convert/LineRepair.h"

#include <algorithm>
#include <cctype>
#include <sstream>

#include "azure2/convert/ExtrapLine.h"
#include "azure2/convert/SegLine.h"
#include "azure2/convert/TargetIntLine.h"

namespace azure2 {
namespace convert {

namespace {

std::vector<std::string> Tokens(const std::string &line) {
  std::istringstream in(line);
  std::vector<std::string> out;
  for (std::string t; in >> t;) out.push_back(t);
  return out;
}

std::string Join(const std::vector<std::string> &t) {
  std::string out;
  for (const std::string &s : t) out += (out.empty() ? "" : " ") + s;
  return out;
}

bool Number(const std::string &t) {
  double v = 0.0;
  return WholeNumber(t, v);
}

}  // namespace

bool DecimalCommaNumber(const std::string &t) {
  const std::size_t comma = t.find(',');
  if (comma == std::string::npos || t.find(',', comma + 1) != std::string::npos) return false;
  std::string dotted = t;
  dotted[comma] = '.';
  return comma > 0 && comma + 1 < t.size() && std::isdigit((unsigned char)t[comma - 1]) &&
         std::isdigit((unsigned char)t[comma + 1]) && Number(dotted);
}

bool GluedNumbers(const std::string &t) {
  if (!Number(t)) return false;
  std::size_t i = (t[0] == '-' || t[0] == '+') ? 1 : 0;
  if (i + 1 < t.size() && t[i] == '0' && std::isdigit((unsigned char)t[i + 1])) return true;
  int digits = 0;
  bool leading = true;
  for (; i < t.size() && t[i] != 'e' && t[i] != 'E'; i++) {
    if (!std::isdigit((unsigned char)t[i])) continue;
    if (leading && t[i] == '0') continue;
    leading = false;
    digits++;
  }
  return digits > 17;
}

namespace {

std::vector<std::string> Insert(std::vector<std::string> t, std::size_t at, const std::vector<std::string> &what) {
  t.insert(t.begin() + (long)std::min(at, t.size()), what.begin(), what.end());
  return t;
}

/// The number of tokens before the normalization of a data line (isDiff 2: J and l; dev's isDiff 10: its operator).
std::size_t DataHead(const std::vector<std::string> &t) {
  if (t.size() < 8) return t.size();
  if (t[7] == "2") return 10;
  if (t[7] == "10" && t.size() > 8 && t[8].find(':') != std::string::npos) return 9;
  return 8;
}

}  // namespace

std::string ExpectedFields(const std::string &section) {
  if (section == "segmentsData")
    return "isActive entranceKey exitKey minE maxE minAngle maxAngle isDiff [J l, for isDiff 2] norm varyNorm normError "
           "[energyShift [shiftError varyShift]] dataFile [keyword blocks]";
  if (section == "segmentsTest")
    return "isActive entranceKey exitKey minE maxE eStep minAngle maxAngle angleStep isDiff [J l, for isDiff 2] "
           "[maxOrder, for isDiff 3] [keyword blocks]";
  if (section == "targetInt")
    return "isActive \"segments\" points convolution sigma integration density \"stoppingEquation\" nParameters p1..pn "
           "qCoefficients nQ q1..qn [convolutionCoefficients \"convolutionEquation\" nC c1..cn [straggling coefficient "
           "...]] [keyword blocks]";
  return "";
}

bool OldLineReads(const std::string &section, const std::string &line) {
  std::istringstream stm(line);
  if (section == "segmentsData") {
    SegLine p(stm);
    return !(stm.rdstate() & (std::ios::failbit | std::ios::badbit)) && p.unreadTail().empty() && !p.dataFile().empty() &&
           !DecimalCommaNumber(p.dataFile());
  }
  if (section == "segmentsTest") {
    ExtrapLine p(stm);
    return !(stm.rdstate() & (std::ios::failbit | std::ios::badbit)) && p.unreadTail().empty();
  }
  if (section == "targetInt") {
    TargetIntLine l;
    ReadTargetIntLine(stm, l);
    return !(stm.rdstate() & (std::ios::failbit | std::ios::badbit)) && !l.screeningError && !l.dopplerError &&
           !l.densityFitError && l.unknownKeyword.empty();
  }
  return false;
}

std::vector<LineReading> ProposeReadings(const std::string &section, const std::string &line) {
  std::vector<LineReading> out;
  auto offer = [&](const std::string &label, const std::vector<std::string> &tokens) {
    const std::string proposed = Join(tokens);
    if (proposed == Join(Tokens(line)) || !OldLineReads(section, proposed)) return;
    for (const LineReading &r : out)
      if (r.line == proposed) return;
    out.push_back({label, proposed});
  };
  const std::vector<std::string> t = Tokens(line);

  // Decimal commas: every such token read with a point (a number written in a locale with a decimal comma).
  {
    std::vector<std::string> dotted = t;
    std::string which;
    for (std::string &s : dotted)
      if (DecimalCommaNumber(s)) {
        which += (which.empty() ? "" : ", ") + s;
        s[s.find(',')] = '.';
      }
    if (!which.empty()) offer("decimal commas read as decimal points (" + which + ")", dotted);
  }

  // A number as wide as its 15-character column run into the next field (the setup program of 2026-09-25 .. 10-07):
  // the token split after its first 15 or more characters where both parts are numbers.
  for (std::size_t i = 0; i < t.size(); i++) {
    if (t[i].size() < 16 || !Number(t[i])) continue;
    for (std::size_t k = 15; k < t[i].size(); k++) {
      const std::string a = t[i].substr(0, k), b = t[i].substr(k);
      if (!Number(a) || !Number(b)) continue;
      std::vector<std::string> split = t;
      split[i] = a;
      split.insert(split.begin() + (long)i + 1, b);
      offer("'" + t[i] + "' is two numbers run together by the setup program of 2026-09 (a number as wide as its "
            "column): " + a + " and " + b, split);
    }
  }

  if (section == "segmentsData") {
    const std::size_t head = DataHead(t);
    std::size_t f = head;
    while (f < t.size() && Number(t[f])) f++;
    if (f >= t.size()) return out;  // no data file: nothing to propose, the user names the file
    const std::size_t k = f - head;  // numbers between the head and the data file
    if (k == 1) offer("the normalization's vary flag and error are missing: 0 (fixed) and 0 (no prior)", Insert(t, head + 1, {"0", "0"}));
    if (k == 2) offer("the normalization's error is missing: 0 (no prior)", Insert(t, head + 2, {"0"}));
    if (k == 5) {
      offer("the energy shift's vary flag is missing: 0 (the shift fixed)", Insert(t, head + 5, {"0"}));
      offer("the normalization's error is missing: 0 (no prior); the next three numbers are the energy shift, its "
            "error and its vary flag",
            Insert(t, head + 2, {"0"}));
    }
    // a data file whose name holds blanks, unquoted: the words up to the next number (or keyword-free end) quoted
    std::size_t g = f + 1;
    while (g < t.size() && !Number(t[g])) g++;
    if (g > f + 1) {
      std::vector<std::string> quoted(t.begin(), t.begin() + (long)f);
      std::vector<std::string> name(t.begin() + (long)f, t.begin() + (long)g);
      quoted.push_back("\"" + Join(name) + "\"");
      quoted.insert(quoted.end(), t.begin() + (long)g, t.end());
      offer("the data file's name holds blanks: \"" + Join(name) + "\"", quoted);
    }
  } else if (section == "segmentsTest") {
    if (t.size() == 9) {
      offer("the data type (isDiff) is missing: 0, the angle-integrated cross section", Insert(t, 9, {"0"}));
      offer("the data type (isDiff) is missing: 1, the differential cross section", Insert(t, 9, {"1"}));
    }
  } else if (section == "targetInt") {
    if (t.size() == 8 && t[7] == "\"\"") offer("the stopping-power equation is empty: 0 parameters", Insert(t, 8, {"0"}));
  }
  return out;
}

std::string LineProblem::Text() const {
  std::ostringstream o;
  o << "line " << fileLine << " of the file, in <" << section << ">:\n    " << text << "\n";
  o << "  expected: " << expected << "\n";
  if (readings.empty()) {
    o << "  no reading can be proposed: edit the line in a copy of the project and convert the copy";
  } else {
    o << "  possible readings (none is applied by itself; take one, or edit the line, in a copy of the project and "
         "convert the copy):";
    for (std::size_t i = 0; i < readings.size(); i++)
      o << "\n    " << (i + 1) << ") " << readings[i].label << ":\n       " << readings[i].line;
  }
  return o.str();
}

}  // namespace convert
}  // namespace azure2
