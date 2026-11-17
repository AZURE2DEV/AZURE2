// The grids of an extrapolation (Grids.h) against values written by hand, and the refinement of an adaptive grid
// (AdaptiveGrid.h) against the closed form of a Lorentzian line (project_format3.md section 3.8):
//   1. range step / points / points log / list / file / the start grid of adaptive, in MeV and keV, lab and c.m.
//      (E_lab = E_cm (m1 + m2)/m2, and the relativistic inverse);
//   2. a Lorentzian of width 1e-6 (1 eV at 0.3001234 MeV) on adaptive 0.2 .. 0.4, start step 0.01, tolerance 1e-3: the peak
//      to 1e-3, the trapezoid integral to 1e-3 of pi Gamma/2 (less the tails beyond the range), fewer than 400
//      points; the uniform grid of the same count misses the peak by more than 50 %;
//   3. the grid does not depend on how a batch is computed (forward, reversed);
//   4. max_points reached: the flag, the largest deviation left above the tolerance, the count kept;
//   5. two components (two lines at different energies): the union grid refines at both.
//   experiment_grids_test
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "azure2/experiment/AdaptiveGrid.h"
#include "azure2/experiment/Grids.h"
#include "azure2/experiment/ProjectPaths.h"
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

double Lorentzian(double e, double centre, double width) {
  const double h = 0.5 * width;
  return h * h / ((e - centre) * (e - centre) + h * h);
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
    check("file in keV, c.m.: three values, comments and blank lines skipped, the path as written",
          v.ok && v.values.size() == 3 && Near(v.values[0], 0.1 * 13.00728 / 12) && Near(v.values[1], 0.25 * 13.00728 / 12) &&
              Near(v.values[2], 0.4 * 13.00728 / 12) && v.file == "grids_test_grid.txt");
    GridSetting elsewhere = cm;
    elsewhere.projectDir = "no_such_directory";
    const GridValues w = ReadGridValues(g, elsewhere);
    // The file exists in the working directory, so the refusal shows that it was looked for in the project's; the
    // message names it as the project writes it (ProjectPaths.h).
    check("file resolved in the project's directory (absent there): E_DATA_FILE, named as written",
          !w.ok && std::string(w.code) == "E_DATA_FILE" && w.error.find("grid file grids_test_grid.txt") != std::string::npos &&
              w.error.find("no_such_directory") == std::string::npos);
    // One rule for "absolute", the same on every system: a Windows drive path is not relative to the project.
    using azure2::experiment::IsAbsolutePath;
    using azure2::experiment::PathToOpen;
    check("absolute: /x, \\x, C:/x, c:\\x; not x, a/b, C:x, 1:/x, an empty path",
          IsAbsolutePath("/x") && IsAbsolutePath("\\x") && IsAbsolutePath("C:/x") && IsAbsolutePath("c:\\x") &&
              !IsAbsolutePath("x") && !IsAbsolutePath("a/b") && !IsAbsolutePath("C:x") && !IsAbsolutePath("1:/x") &&
              !IsAbsolutePath(""));
    check("the file to open: a drive path as written, a relative one in the project's directory",
          PathToOpen("D:/p", "D:/t/g.dat") == "D:/t/g.dat" && PathToOpen("D:/p", "g.dat") == "D:/p/g.dat" &&
              PathToOpen("/p/", "g.dat") == "/p/g.dat" && PathToOpen("/p", "") == "/p" && PathToOpen("", "g.dat") == "g.dat");
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

  // ---- 2. a narrow line on an adaptive grid
  const double centre = 0.3001234, width = 1e-6, lo = 0.2, hi = 0.4;  // the centre on no uniform grid of [lo, hi]
  auto model = [&](const std::vector<double> &x) {
    std::vector<std::vector<double>> v;
    for (double e : x) v.push_back({Lorentzian(e, centre, width)});
    return v;
  };
  std::vector<double> start;
  for (int k = 0; k <= 20; k++) start.push_back(lo + k * 0.01);
  for (double k : {0.0, -0.5, 0.5, -1.0, 1.0, -2.0, 2.0, -4.0, 4.0}) start.push_back(centre + k * width);
  const AdaptiveGridResult r = RefineGrid(start, 1e-3, 2000, model);
  double peak = 0.0, integral = 0.0;
  for (std::size_t i = 0; i < r.x.size(); i++) {
    peak = std::max(peak, r.values[i][0]);
    if (i) integral += 0.5 * (r.x[i] - r.x[i - 1]) * (r.values[i][0] + r.values[i - 1][0]);
  }
  const double h = 0.5 * width;
  const double exact = h * (std::atan((hi - centre) / h) + std::atan((centre - lo) / h));  // the integral over [lo, hi]
  std::printf("  adaptive: %zu points, %d sweeps, %ld evaluations, peak %.9f, integral %.6e (exact %.6e, %.2e)\n", r.x.size(),
              r.sweeps, r.evaluations, peak, integral, exact, integral / exact - 1);
  check("adaptive: the peak to 1e-3", std::fabs(peak - 1.0) < 1e-3);
  check("adaptive: the integral to 1e-3", std::fabs(integral / exact - 1.0) < 1e-3);
  check("adaptive: fewer than 400 points, not stopped by max_points", r.x.size() < 400 && !r.budgetReached);
  check("adaptive: ascending, inside [lo, hi]",
        std::is_sorted(r.x.begin(), r.x.end()) && r.x.front() == lo && r.x.back() == hi);
  double uniformPeak = 0.0;
  for (std::size_t k = 0; k < r.x.size(); k++)
    uniformPeak = std::max(uniformPeak, Lorentzian(lo + (hi - lo) * (double)k / (double)(r.x.size() - 1), centre, width));
  check("the uniform grid of the same count misses the peak by more than 50 %", uniformPeak < 0.5);

  // ---- 3. the grid does not depend on the order in which a batch is computed
  auto reversed = [&](const std::vector<double> &x) {
    std::vector<std::vector<double>> v(x.size());
    for (std::size_t i = x.size(); i-- > 0;) v[i] = {Lorentzian(x[i], centre, width)};
    return v;
  };
  const AdaptiveGridResult r2 = RefineGrid(std::vector<double>(start.rbegin(), start.rend()), 1e-3, 2000, reversed);
  check("the same grid, bit for bit, from the start grid reversed and the batch computed backwards", r2.x == r.x);

  // ---- 4. max_points
  const AdaptiveGridResult r3 = RefineGrid(start, 1e-3, 100, model);
  check("max_points 100: stopped there, the largest deviation left above the tolerance",
        r3.budgetReached && r3.x.size() == 100 && r3.largestDeviation > 1e-3 && r3.largestDeviationAt > lo &&
            r3.largestDeviationAt < hi);

  // ---- 5. two components
  auto two = [&](const std::vector<double> &x) {
    std::vector<std::vector<double>> v;
    for (double e : x) v.push_back({Lorentzian(e, 0.25, 1e-4), Lorentzian(e, 0.35, 1e-4)});
    return v;
  };
  const AdaptiveGridResult r4 = RefineGrid(std::vector<double>(start.begin(), start.begin() + 21), 1e-3, 2000, two);
  auto near = [&](double c) {
    return std::count_if(r4.x.begin(), r4.x.end(), [&](double e) { return std::fabs(e - c) < 1e-4; });
  };
  check("two components: the union grid is refined at both lines (> 10 points within one width of each)",
        near(0.25) > 10 && near(0.35) > 10 && !r4.budgetReached);

  std::printf("%d checks, %d failures\n", checks, failures);
  return failures ? 1 : 0;
}
