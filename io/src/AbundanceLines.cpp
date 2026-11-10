// The grammar of the <abundances> section (AbundanceLines.h).
#include "azure2/io/AbundanceLines.h"

#include <sstream>

namespace azure2 {
namespace io {

AbundanceSection ReadAbundanceSection(std::istream &in) {
  AbundanceSection section;
  std::string line;
  getline(in, line);  // <abundances>
  line = "";
  while (!in.eof() && line != "</abundances>") {
    getline(in, line);
    std::string entry = line;
    size_t hash = entry.find('#');
    if (hash != std::string::npos) entry = entry.substr(0, hash);
    if (entry.find_first_not_of(" \t\r") == std::string::npos) continue;
    if (line == "</abundances>") break;
    AbundanceLine parsed;
    parsed.text = line;
    // optional "<name>:" prefix (one token ending in a colon)
    std::istringstream probe(entry);
    std::string first;
    probe >> first;
    if (!first.empty() && first.back() == ':') {
      parsed.sample = first.substr(0, first.size() - 1);
      std::string rest;
      getline(probe, rest);
      entry = rest;
    }
    std::istringstream stm(entry);
    std::string head;
    if (!(stm >> head)) continue;
    if (head == "sum") {
      if (!(stm >> parsed.sumSigma)) {
        section.problem = AbundanceSection::Problem::BadSum;
        section.badLine = line;
        return section;
      }
      parsed.isSum = true;
      section.lines.push_back(parsed);
      continue;
    }
    std::istringstream stm2(entry);
    if (!(stm2 >> parsed.Z >> parsed.A >> parsed.value >> parsed.vary >> parsed.sigma)) {
      section.problem = AbundanceSection::Problem::BadEntry;
      section.badLine = line;
      return section;
    }
    section.lines.push_back(parsed);
  }
  section.complete = (line == "</abundances>");
  return section;
}

}  // namespace io
}  // namespace azure2
