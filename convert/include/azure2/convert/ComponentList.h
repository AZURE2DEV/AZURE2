#ifndef COMPONENTLIST_H
#define COMPONENTLIST_H

#include <cstdlib>
#include <iomanip>
#include <istream>
#include <sstream>
#include <string>
#include <vector>

/*!
 * Reads the component block of a composite (SUM / RATIO) segment line, after the operation type:
 *   <n> e1 x1 angle1 ... en xn anglen              (angle -999 = the segment's own angles)
 *   -1 <n> e1 x1 angle1 scaling1 ...               (scaling -999 = none, i.e. 1)
 *   <n> e1 x1 ... en xn                            (old GUI files: entrance and exit only)
 * componentsList receives "Entrance: e, Exit: x[, Angle: a][, Scaling: s]" entries joined by ';' (parsed by WorkingCopy).
 * Returns "" on success, else what could not be read; the caller fails the line with it.  A non-empty `warning` is
 * a note the caller prints with the segment number.  Shared by SegLine and ExtrapLine.
 *
 * Without the -1 marker the layout is decided as the GUI decides it (SegmentsTab::readSegDataFile): three numbers
 * per component when at least 3n numbers follow the count, else two when at least 2n do.  Too few numbers for n
 * components: no components and a warning (the line is read as a plain segment, as every engine before did).
 *
 * Before 2026-11-01: a scaling that was not a number became 1.0 without a word; the two-number layout was misread
 * (the next component's entrance key taken as an angle, the last component dropped silently: "1 0 2 2 2 2 3" gave the
 * single component 2->2 at a fixed 2 deg instead of 2->2 and 2->3); and the numbers went into the list through
 * std::to_string, i.e. "%f" (six decimals: a scaling of 1.5e-7 was read as 0).  Now the scaling format is strict,
 * the two-number layout is read as the GUI reads it, and the numbers keep 17 significant digits.
 */
/// One component as read: pair keys, and the fixed angle (deg) and scaling when the line gives them (-999 = none).
struct ComponentEntry {
  int entrance = 0, exit = 0;
  bool hasAngle = false;
  double angle = 0.0;
  bool hasScaling = false;
  double scaling = 1.0;
};

/// True when the whole token reads as a number.
inline bool IsNumberToken(const std::string &t) {
  if (t.empty()) return false;
  char *end = nullptr;
  std::strtod(t.c_str(), &end);
  return end != t.c_str() && *end == '\0';
}

/// As ReadComponentList below; `entries` receives the components as values (the converter to format 3 reads them).
inline std::string ReadComponentList(std::istream &in, std::string &componentsList, std::string &warning,
                                     std::vector<ComponentEntry> &entries) {
  auto fmt = [](double v) {
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
  };
  int numComponents;
  if (!(in >> numComponents)) return "";  // no component block (as before)
  bool hasScalingToken = false;
  if (numComponents == -1) {
    hasScalingToken = true;
    if (!(in >> numComponents)) return "the component count after the -1 (scaling) marker is not a whole number";
  }
  std::vector<std::string> components;
  if (hasScalingToken) {
    for (int i = 0; i < numComponents; i++) {
      const std::string which = "component " + std::to_string(i + 1) + " of " + std::to_string(numComponents);
      int entrance, exit;
      double angle, scaling;
      if (!(in >> entrance >> exit >> angle)) return which + ": expected entrance key, exit key, angle, scaling";
      if (!(in >> scaling)) return which + ": the scaling factor is not a number";
      std::string componentStr = "Entrance: " + std::to_string(entrance) + ", Exit: " + std::to_string(exit);
      if (angle > -900.0) componentStr += ", Angle: " + fmt(angle);
      if (scaling > -900.0) componentStr += ", Scaling: " + fmt(scaling);
      components.push_back(componentStr);
      ComponentEntry e;
      e.entrance = entrance;
      e.exit = exit;
      e.hasAngle = angle > -900.0;
      e.angle = angle;
      e.hasScaling = scaling > -900.0;
      e.scaling = scaling;
      entries.push_back(e);
    }
  } else if (numComponents > 0) {
    // Count the numbers that follow (without consuming them) to choose the layout.
    const std::streampos start = in.tellg();
    long available = 0;
    std::string token;
    while (in >> token && IsNumberToken(token)) available++;
    in.clear();
    in.seekg(start);
    const int perComponent = (available >= 3L * numComponents) ? 3 : (available >= 2L * numComponents) ? 2 : 0;
    if (perComponent == 0) {
      warning = "the component list declares " + std::to_string(numComponents) + " component(s) but holds only " +
                std::to_string(available) + " number(s); no components are used (the line is a plain segment)";
      in.setstate(std::ios::failbit);  // as before: nothing after the block is read as positional fields
      return "";
    }
    for (int i = 0; i < numComponents; i++) {
      double e, x, angle = -999.0;
      in >> e >> x;
      if (perComponent == 3) in >> angle;
      std::string componentStr = "Entrance: " + std::to_string((int)e) + ", Exit: " + std::to_string((int)x);
      if (angle > -900.0) componentStr += ", Angle: " + fmt(angle);
      components.push_back(componentStr);
      ComponentEntry c;
      c.entrance = (int)e;
      c.exit = (int)x;
      c.hasAngle = angle > -900.0;
      c.angle = angle;
      entries.push_back(c);
    }
  }
  for (size_t i = 0; i < components.size(); i++) componentsList += (i ? ";" : "") + components[i];
  return "";
}

#endif
