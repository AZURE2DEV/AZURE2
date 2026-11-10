// The values of an extrapolation's grid (Grids.h).
#include "azure2/experiment/Grids.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>

#include "azure2/experiment/Diagnostics.h"
#include "azure2/experiment/Frames.h"

namespace azure2 {
namespace experiment {

namespace {

// The value in the engine's unit of its dimension (BuildContext.h's InBase: bit-exact for factor 1).
double InBase(double value, const Unit *unit) { return unit && unit->factor != 1.0 ? value * unit->factor : value; }

GridValues Fail(const char *code, const std::string &why) {
  GridValues g;
  g.code = code;
  g.error = why;
  return g;
}

// The numbers of a grid file, one per line ("#" starts a comment).
bool ReadGridFile(const std::string &path, std::vector<double> &values, std::string &why) {
  std::ifstream in(path);
  if (!in) {
    why = "cannot read the grid file " + path;
    return false;
  }
  std::string line;
  int n = 0;
  while (std::getline(in, line)) {
    n++;
    const std::size_t hash = line.find('#');
    std::istringstream row(hash == std::string::npos ? line : line.substr(0, hash));
    std::string token, extra;
    if (!(row >> token)) continue;
    char *end = nullptr;
    const double v = std::strtod(token.c_str(), &end);
    if (end == token.c_str() || *end != '\0' || !std::isfinite(v) || (row >> extra)) {
      why = path + ", line " + std::to_string(n) + ": '" + line + "' is not one number";
      return false;
    }
    values.push_back(v);
  }
  if (values.empty()) {
    why = "the grid file " + path + " holds no value";
    return false;
  }
  return true;
}

}  // namespace

GridValues ReadGridValues(const Grid &g, const GridSetting &s) {
  GridValues out;
  const bool cm = s.energies && g.frame == "cm";
  Partition entrance;
  if (cm) {
    if (!s.entrance || s.entrance->kind != 0 || !(s.entrance->lightMass > 0.0) || !(s.entrance->heavyMass > 0.0))
      return Fail(codes::kNotBuilt, "an energy grid in the c.m. frame needs a particle pair as the entrance (a beam on "
                                    "a target at rest, with its masses); give a beta-delayed grid in lab energies");
    entrance.m1 = s.entrance->lightMass;
    entrance.m2 = s.entrance->heavyMass;
  }
  // The engine's unit, then the lab frame of a c.m. energy.
  auto value = [&](double v) {
    const double x = InBase(v, g.unit);
    return cm ? BeamLabFromCM(entrance, x, s.relativistic) : x;
  };
  bool bounds = false;  // lo and hi are the written bounds (a range or adaptive grid), not the first and last value
  if (g.kind == "list") {
    for (const Number &n : g.values) out.values.push_back(value(n.value));
  } else if (g.kind == "file") {
    const bool absolute = !g.path.empty() && (g.path[0] == '/' || g.path[0] == '\\');
    out.file = absolute || s.projectDir.empty() ? g.path : s.projectDir + "/" + g.path;
    std::vector<double> raw;
    std::string why;
    if (!ReadGridFile(out.file, raw, why)) return Fail(codes::kDataFile, why);
    for (double v : raw) out.values.push_back(value(v));
  } else if (g.kind == "range" && g.hasStep) {
    // The expressions of decision F10 in the frame written, then each value converted.
    const double lo = InBase(g.lo.value, g.unit), hi = InBase(g.hi.value, g.unit), step = InBase(g.step.value, g.unit);
    const long n = step > 0 ? (long)std::floor((hi - lo) / step + 1e-9) : 0;
    for (long k = 0; k <= n; k++) out.values.push_back(cm ? value(lo + k * step) : lo + k * step);
    out.lo = value(g.lo.value);
    out.hi = value(g.hi.value);
    bounds = true;
  } else if (g.kind == "range" && g.hasPoints) {
    const double lo = InBase(g.lo.value, g.unit), hi = InBase(g.hi.value, g.unit);
    if (g.log && !(lo > 0.0)) return Fail(codes::kEmptyRange, "a log grid needs lo > 0");
    for (long k = 0; k < g.points; k++) {
      const double t = g.points > 1 ? (double)k / (double)(g.points - 1) : 0.0;
      const double x = g.log ? lo * std::pow(hi / lo, t) : lo + t * (hi - lo);
      out.values.push_back(cm ? value(x) : x);
    }
    out.lo = value(g.lo.value);
    out.hi = value(g.hi.value);
    bounds = true;
  } else if (g.kind == "adaptive") {
    if (!s.energies) return Fail(codes::kNotBuilt, "an adaptive grid is an energy grid; give the angles as a range or list");
    // The uniform start grid of start_step in the frame written, closed by hi.
    const double lo = InBase(g.lo.value, g.unit), hi = InBase(g.hi.value, g.unit), step = InBase(g.startStep.value, g.unit);
    const long n = (long)std::floor((hi - lo) / step + 1e-9);
    for (long k = 0; k <= n; k++) out.values.push_back(value(lo + k * step));
    if (lo + n * step < hi - 1e-9 * step) out.values.push_back(value(hi));
    out.lo = value(g.lo.value);
    out.hi = value(g.hi.value);
    bounds = true;
  } else {
    return Fail(codes::kNotBuilt, "a grid of kind " + g.kind + " without step or points");
  }
  if (!bounds && !out.values.empty()) {
    out.lo = out.values.front();
    out.hi = out.values.back();
  }
  out.ok = true;
  return out;
}

}  // namespace experiment
}  // namespace azure2
