// The Calculate tab runs the engine on a worker thread that owns its own
// Config, filled field by field from the window's.  The <config> block is
// parsed only by the CLI entry point, so an option the GUI sets on its own
// Config reaches a GUI run only if that copy list carries it.  The
// hard-sphere option (includeAllJPi / maxOrbitalMomentum) was left out when
// it was added, and every GUI calculation silently ran without it.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <iostream>
#include <sstream>
#include "AZUREMainThread.h"
#include "RunTab.h"
#include "Config.h"
struct SegPairs {int firstPair; int secondPair;};

Config* g_config = nullptr;
void exitMessage(const Config&) {}
bool checkExternalCapture(Config&, const std::vector<SegPairs>&) { return true; }
bool readSegmentFile(const Config&, std::vector<SegPairs>&) { return true; }
void startMessage(const Config&) {}

static int fails = 0;
static void ok(const char* what, bool cond) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what << std::endl;
  if(!cond) fails++;
}

int main(int argc, char** argv) {
  QApplication app(argc, argv);
  std::ostringstream sink;
  Config config(sink);
  config.configfile = "project.azr";
  config.includeAllJPi = true;
  config.hardSphereMaxL = 3;
  config.useAdaptiveGrid = true;
  config.chiVariance = 0.25;
  config.relativistic = true;
  config.relativisticReferencePair = 2;
  config.ecConvectionCurrent = false;  // default is true: a copied false shows the field is carried over

  RunTab tab;
  AZUREMainThread thread(&tab, config);
  const Config& copy = thread.configure();
  ok("project file copied", copy.configfile == "project.azr");
  ok("existing fields still copied", copy.useAdaptiveGrid && copy.chiVariance == 0.25);
  ok("hard-sphere switch reaches the worker Config", copy.includeAllJPi);
  ok("hard-sphere L limit reaches the worker Config", copy.hardSphereMaxL == 3);
  ok("relativistic kinematics reaches the worker Config", copy.relativistic && copy.relativisticReferencePair == 2);
  ok("EC convection-current switch reaches the worker Config", !copy.ecConvectionCurrent);

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
