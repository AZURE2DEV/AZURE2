// The grammar of a <levels> line and of the section (LevelRows.h).
#include "azure2/io/LevelRows.h"

#include <sstream>

namespace azure2 {
namespace io {

LevelRow::LevelRow(std::istream &stream) {
  stream >> levelJ_ >> levelPi_ >> levelE_ >> levelFix_ >> aa_ >> ir_ >> s_ >> l_ >> levelID_ >> isActive_ >>
      channelFix_ >> gamma_ >> j1_ >> pi1_ >> j2_ >> pi2_ >> e2_ >> m1_ >> m2_ >> z1_ >> z2_ >> entranceSepE_ >> sepE_ >>
      j3_ >> pi3_ >> e3_ >> pType_ >> chRad_ >> g1_ >> g2_ >> ecMultMask_;

  // The optional trailing columns: 1 = e1, 2 = binding energy and flag, 3 = binding energy, flag and e1.  A token that
  // is not a number sets failbit, more than three fields too.
  e1_ = 0.;
  bindingE_ = 0.0;
  gammaIsRWA_ = 0;
  if (!stream.fail() && !stream.eof()) {
    double v[3];
    int n = 0;
    stream >> std::ws;
    while (!stream.eof() && n < 4) {
      double x;
      if (!(stream >> x)) break;  // failbit set: reported by the caller
      if (n < 3) v[n] = x;
      n++;
      if (!stream.eof()) stream >> std::ws;  // not at eof: a newer libstdc++'s ws sets failbit there
    }
    if (!stream.fail()) {
      if (n == 1) {
        e1_ = v[0];
      } else if (n == 2 || n == 3) {
        bindingE_ = v[0];
        gammaIsRWA_ = (int)v[1];
        if (n == 3) e1_ = v[2];
      } else if (n > 3) {
        stream.setstate(std::ios::failbit);
      }
    }
  }
  s_ /= 2.;
  l_ /= 2;
}

LevelSection ReadLevelSection(std::istream &in) {
  LevelSection section;
  std::string line;
  getline(in, line);  // <levels>
  line = "";
  while (!in.eof() && line != "</levels>") {
    getline(in, line);
    bool empty = true;
    for (unsigned int i = 0; i < line.size(); ++i)
      if (line[i] != ' ' && line[i] != '\t') {
        empty = false;
        break;
      }
    if (empty) continue;
    if (!in.eof() && line != "</levels>") {
      std::istringstream stm;
      stm.str(line);
      LevelRow row(stm);
      if (stm.rdstate() & (std::stringstream::failbit | std::stringstream::badbit)) {
        section.malformed = true;
        section.badLine = line;
        return section;
      }
      section.rows.push_back(row);
    }
  }
  section.complete = (line == "</levels>");
  return section;
}

}  // namespace io
}  // namespace azure2
