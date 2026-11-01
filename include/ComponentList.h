#ifndef COMPONENTLIST_H
#define COMPONENTLIST_H

#include <iomanip>
#include <istream>
#include <sstream>
#include <string>
#include <vector>

/*!
 * Reads the component block of a composite (SUM / RATIO) segment line, after the operation type:
 *   <n> e1 x1 angle1 ... en xn anglen              (angle -999 = the segment's own angles)
 *   -1 <n> e1 x1 angle1 scaling1 ...               (scaling -999 = none, i.e. 1)
 * and an old two-token form (e x) as a fallback.  componentsList receives "Entrance: e, Exit: x[, Angle: a]
 * [, Scaling: s]" entries joined by ';' (parsed by EData).  Returns "" on success, else what could not be read; the
 * caller fails the line with it.  Shared by SegLine and ExtrapLine.
 *
 * Before 2026-11-01 a scaling that was not a number became 1.0 without a word (and the stream, left failed, then
 * dropped the components after it and the UPOS block), an unreadable component was left out, and the numbers went
 * into the list through std::to_string, i.e. "%f": six decimals, so a scaling of 1.5e-7 was read as 0 and 0.00012345
 * as 0.000123.  Now an unreadable component fails the line, and the numbers keep 17 significant digits.
 */
inline std::string ReadComponentList(std::istream &in, std::string &componentsList) {
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
  for (int i = 0; i < numComponents; i++) {
    const std::string which = "component " + std::to_string(i + 1) + " of " + std::to_string(numComponents);
    int entrance, exit;
    double angle;
    if (in >> entrance >> exit >> angle) {
      double scaling = 1.0;
      if (hasScalingToken && !(in >> scaling)) return which + ": the scaling factor is not a number";
      std::string componentStr = "Entrance: " + std::to_string(entrance) + ", Exit: " + std::to_string(exit);
      if (angle > -900.0) componentStr += ", Angle: " + fmt(angle);
      if (hasScalingToken && scaling > -900.0) componentStr += ", Scaling: " + fmt(scaling);
      components.push_back(componentStr);
    } else {
      // Fallback to the old format (entrance, exit only) for backward compatibility
      in.clear();
      if (!(in >> entrance >> exit))
        return which + ": expected entrance key, exit key, angle" + std::string(hasScalingToken ? ", scaling" : "");
      components.push_back("Entrance: " + std::to_string(entrance) + ", Exit: " + std::to_string(exit));
    }
  }
  for (size_t i = 0; i < components.size(); i++) componentsList += (i ? ";" : "") + components[i];
  return "";
}

#endif
