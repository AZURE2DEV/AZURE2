// The base a fitted nuisance parameter's grids and limits are laid out from (azure2/experiment/LayoutBase.h,
// docs/upgrade/parameter_kinds.md section 4), against numbers written by hand:
//   1. the three spans' formulas (with and without a prior, grid_max);
//   2. no centre written, a fixed parameter, a centre equal to the start: the start itself, no note;
//   3. a start inside the centre's range (the fitted value of a fit started at the centre): the centre, no note;
//      at the range's ends inclusive;
//   4. a start outside it (above, below, above a grid_max): the start and a note naming the range and the centre.
//   experiment_layout_base_test
#include <cmath>
#include <cstdio>
#include <string>

#include "azure2/experiment/LayoutBase.h"

using namespace azure2::experiment;

namespace {

int failures = 0;

void Check(const std::string &what, bool ok) {
  std::printf("%s %s\n", ok ? "  ok  " : "  FAIL", what.c_str());
  if (!ok) failures++;
}

PriorCentre Centre(double c) {
  PriorCentre p;
  p.written = true;
  p.value = c;
  return p;
}

}  // namespace

int main() {
  std::printf("1. the spans\n");
  Check("density, a prior: max(N + 3 sigma, 1.1 N)", kDensitySpan.Highest(1e18, 2e16, 0.0) == 1.1 * 1e18 &&
                                                         kDensitySpan.Highest(1e18, 1e17, 0.0) == 1e18 + 3e17);
  Check("density, no prior: 1.5 N; lowest N/100", kDensitySpan.Highest(1e18, 0.0, 0.0) == 1.5e18 &&
                                                      kDensitySpan.Lowest(1e18) == 1e16);
  Check("Doppler: max(T + 3 sigma, 1.21 T), 2 T, T/4", kDopplerSpan.Highest(300.0, 20.0, 0.0) == 1.21 * 300.0 &&
                                                           kDopplerSpan.Highest(300.0, 0.0, 0.0) == 600.0 &&
                                                           kDopplerSpan.Lowest(300.0) == 75.0);
  Check("apparatus: max(s + 3 sigma, 1.1 s), 1.5 s, s/10", kApparatusSpan.Highest(1.0, 0.04, 0.0) == 1.0 + 3.0 * 0.04 &&
                                                               kApparatusSpan.Highest(14.0, 0.0, 0.0) == 21.0 &&
                                                               kApparatusSpan.Lowest(14.0) == 0.1 * 14.0);
  Check("grid_max wins", kDensitySpan.Highest(1e18, 2e16, 3e18) == 3e18 && kDopplerSpan.Highest(300.0, 0.0, 900.0) == 900.0);

  std::printf("2. the start\n");
  const PriorCentre none;
  LayoutBase b = ChooseLayoutBase(kDensitySpan, 1.05e18, none, true, 2e16, 0.0, "t.areal_density", "atoms/cm2");
  Check("no centre written: the start", b.value == 1.05e18 && b.note.empty());
  b = ChooseLayoutBase(kDensitySpan, 1.05e18, Centre(1e18), false, 2e16, 0.0, "t.areal_density", "atoms/cm2");
  Check("a fixed parameter: the start", b.value == 1.05e18 && b.note.empty());
  b = ChooseLayoutBase(kDopplerSpan, 300.0, Centre(300.0), true, 20.0, 0.0, "D.doppler", "K");
  Check("centre = start: the start", b.value == 300.0 && b.note.empty());

  std::printf("3. the centre\n");
  b = ChooseLayoutBase(kDensitySpan, 1.0829246e18, Centre(1e18), true, 2e16, 0.0, "t.areal_density", "atoms/cm2");
  Check("a start inside [N/100, 1.1 N]: the centre", b.value == 1e18 && b.note.empty());
  b = ChooseLayoutBase(kDopplerSpan, 1.21 * 300.0, Centre(300.0), true, 20.0, 0.0, "D.doppler", "K");
  Check("a start at the upper end (363 K): the centre", b.value == 300.0 && b.note.empty());
  b = ChooseLayoutBase(kApparatusSpan, 0.1, Centre(1.0), true, 0.04, 0.0, "stopping_scale.TiN.1H", "1");
  Check("a start at the lower end (s/10): the centre", b.value == 1.0 && b.note.empty());
  b = ChooseLayoutBase(kDopplerSpan, 800.0, Centre(300.0), true, 20.0, 900.0, "D.doppler", "K");
  Check("a start inside a grid_max above the formula: the centre", b.value == 300.0 && b.note.empty());

  std::printf("4. a start outside the centre's range\n");
  b = ChooseLayoutBase(kDopplerSpan, 400.0, Centre(300.0), true, 20.0, 0.0, "D.doppler", "K");
  Check("above (400 K > 363 K): the start and a note", b.value == 400.0 &&
        b.note == "Note: D.doppler starts at 400 K, outside [75, 363] K, the range laid out from its prior's centre "
                  "300 K; its grids and limits are laid out from the start.");
  std::printf("    %s\n", b.note.c_str());
  b = ChooseLayoutBase(kApparatusSpan, 0.05, Centre(1.0), true, 0.04, 0.0, "stopping_scale.TiN.1H", "1");
  Check("below (0.05 < 0.1): the start and a note without a unit", b.value == 0.05 &&
        b.note.find("outside [0.1, 1.12],") != std::string::npos);
  b = ChooseLayoutBase(kDensitySpan, 2.5e18, Centre(1e18), true, 2e16, 2e18, "t.areal_density", "atoms/cm2");
  Check("above a grid_max: the start and a note", b.value == 2.5e18 && !b.note.empty());

  if (failures) {
    std::printf("%d failure(s)\n", failures);
    return 1;
  }
  std::printf("PASS\n");
  return 0;
}
