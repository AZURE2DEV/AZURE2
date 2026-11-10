// The grids of an extrapolation (Grids.h) against values written by hand (project_format3.md section 3.8):
//   1. range step / points / points log / list / file / the start grid of adaptive, in MeV and keV, lab and c.m.
//      (E_lab = E_cm (m1 + m2)/m2, and the relativistic inverse);
//   experiment_grids_test
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "azure2/experiment/Grids.h"
#include "azure2/experiment/Units.h"

using namespace azure2::experiment;

namespace {

int failures = 0, checks = 0;

void check(const std::string &what, bool ok) {
  checks++;
  if (!ok) {
    std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

bool Near(double a, double b, double rel = 1e-14) { return std::fabs(a - b) <= rel * std::max(std::fabs(a), std::fabs(b)); }

const Unit *UnitOf(const std::string &token) {
  for (const Unit &u : Units())
    if (token == u.token) return &u;
  return nullptr;
}

Number Num(double v) {
  Number n;
  n.value = v;
  return n;
}

Grid Range(double lo, double hi, const char *unit, const char *frame) {
  Grid g;
  g.kind = "range";
  g.lo = Num(lo);
  g.hi = Num(hi);
  g.unit = UnitOf(unit);
  g.frame = frame;
  return g;
}

}  // namespace

int main() {
  // ---- 1. static grids
  PairInfo pair;  // p + 12C: E_lab = E_cm (1.00728 + 12)/12
  pair.key = 1;
  pair.lightMass = 1.00728;
  pair.heavyMass = 12.0;
  GridSetting lab;
  GridSetting cm;
  cm.entrance = &pair;
  {
    Grid g = Range(0.1, 0.4, "MeV", "lab");
    g.hasStep = true;
    g.step = Num(0.1);
    const GridValues v = ReadGridValues(g, lab);
    check("range step: 0.1, 0.2, 0.30000000000000004, 0.4 (lo + k step)",
          v.ok && v.values == std::vector<double>({0.1, 0.1 + 0.1, 0.1 + 2 * 0.1, 0.1 + 3 * 0.1}) && v.lo == 0.1 && v.hi == 0.4);
  }
  {
    Grid g = Range(100, 300, "keV", "lab");
    g.hasPoints = true;
    g.points = 3;
    const GridValues v = ReadGridValues(g, lab);
    check("range points 3 in keV: 0.1, 0.2, 0.3 MeV", v.ok && v.values.size() == 3 && Near(v.values[0], 0.1) &&
                                                         Near(v.values[1], 0.2) && Near(v.values[2], 0.3));
    g.log = true;
    const GridValues w = ReadGridValues(g, lab);
    check("range points 3 log: 0.1, 0.17320508075688773, 0.3 MeV",
          w.ok && w.values.size() == 3 && Near(w.values[0], 0.1) && Near(w.values[1], 0.17320508075688773) && Near(w.values[2], 0.3));
    g.lo = Num(0.0);
    check("range points log from 0: refused", !ReadGridValues(g, lab).ok);
  }
  {
    Grid g;
    g.kind = "list";
    g.values = {Num(1.2), Num(0.6)};
    g.unit = UnitOf("MeV");
    g.frame = "cm";
    const GridValues v = ReadGridValues(g, cm);
    check("list in c.m.: 1.2 and 0.6 MeV -> 1.300728, 0.650364 MeV lab, in the order written",
          v.ok && v.values.size() == 2 && Near(v.values[0], 1.300728) && Near(v.values[1], 0.650364) && v.lo == v.values[0]);
    check("list in c.m. without the entrance pair's masses: refused", !ReadGridValues(g, lab).ok);
    GridSetting relativistic = cm;
    relativistic.relativistic = true;
    const GridValues r = ReadGridValues(g, relativistic);
    // E_lab = T (2 (m1 + m2) + T)/(2 m2), masses in MeV (u = 931.49410242 MeV)
    const double m1 = 1.00728 * 931.49410242, m2 = 12.0 * 931.49410242;
    check("list in c.m., relativistic: T (2 (m1 + m2) + T)/(2 m2)",
          r.ok && Near(r.values[0], 1.2 * (2 * (m1 + m2) + 1.2) / (2 * m2), 1e-9));
  }
  {
    std::ofstream("grids_test_grid.txt") << "# E (keV, c.m.)\n100\n\n  250   # a comment\n400\n";
    Grid g;
    g.kind = "file";
    g.path = "grids_test_grid.txt";
    g.unit = UnitOf("keV");
    g.frame = "cm";
    const GridValues v = ReadGridValues(g, cm);
    check("file in keV, c.m.: three values, comments and blank lines skipped, the path as opened",
          v.ok && v.values.size() == 3 && Near(v.values[0], 0.1 * 13.00728 / 12) && Near(v.values[1], 0.25 * 13.00728 / 12) &&
              Near(v.values[2], 0.4 * 13.00728 / 12) && v.file == "grids_test_grid.txt");
    GridSetting elsewhere = cm;
    elsewhere.projectDir = "no_such_directory";
    const GridValues w = ReadGridValues(g, elsewhere);
    check("file resolved in the project's directory (absent there): E_DATA_FILE",
          !w.ok && std::string(w.code) == "E_DATA_FILE" && w.error.find("no_such_directory/grids_test_grid.txt") != std::string::npos);
    std::ofstream("grids_test_bad.txt") << "0.1\n0.2 0.3\n";
    g.path = "grids_test_bad.txt";
    check("file with two numbers on a line: refused", !ReadGridValues(g, cm).ok);
  }
  {
    Grid g = Range(0.2, 0.45, "MeV", "lab");
    g.kind = "adaptive";
    g.startStep = Num(0.1);
    const GridValues v = ReadGridValues(g, lab);
    check("adaptive: the start grid 0.2, 0.3, 0.4 and hi 0.45",
          v.ok && v.values.size() == 4 && v.values[1] == 0.2 + 0.1 && v.values[3] == 0.45 && v.lo == 0.2 && v.hi == 0.45);
  }

  std::printf("%d checks, %d failures\n", checks, failures);
  return failures ? 1 : 0;
}
