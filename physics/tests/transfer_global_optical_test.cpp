// transfer_global_optical (docs/upgrade/transfer.md, section "The engine's THM observable on this library"): every row
// of the global optical table against numbers that do not come from its code.
//
//  1. Per row, the ten numbers V, R, a, W, RW, aW, WD, RD, aD, RC (MeV, fm; R = r0 A_t^{1/3}) at two to five
//     (projectile, target, lab energy) points, 1e-9 relative (1e-12 absolute for a zero).  The expected numbers are
//     those of tests/reference/thm_optical_reference.py, which evaluates the RIPL-3 library entries itself
//     (om-parameter-u.dat: An & Cai 6200, Koning & Delaroche 2405/5405, Becchetti & Greenlees 7100/8100, McFadden &
//     Satchler 9100, Avrigeanu 9600) and the FRONT21 front end of TWOFNR for Daehnick and Liang.  The first An & Cai
//     point is the published d + 19F potential of examples/f19_pag_thm.
//  2. Which light ions a row describes, the lookup by name and the list of names.
//
// A new row adds its points to kPoints (and its projectiles to the second part).
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "azure2/physics/transfer/GlobalOptical.h"

using namespace azure2::physics::transfer;

namespace {

int failures = 0;

void Check(bool ok, const std::string &what) {
  if (!ok) {
    std::printf("FAIL: %s\n", what.c_str());
    ++failures;
  }
}

struct Point {
  const char *set;
  int zp, ap, zt, at;
  double e_lab;
  double p[10];  // V, R, a, W, RW, aW, WD, RD, aD, RC
  const char *source;
};

const Point kPoints[] = {
    {"ancai06", 1, 2, 9, 19, 5.83, {92.5676139, 3.066238699, 0.7526218608, 1.466626, 3.581964152, 0.5881598351, 10.651602, 3.711647799, 0.6964409022, 3.476927348}, "RIPL 6200"},
    {"ancai06", 1, 2, 12, 24, 2, {94.02329163, 3.31518301, 0.7553446892, 1.2284, 3.863971379, 0.5665500859, 10.7688, 3.999921854, 0.7098389467, 3.75850238}, "RIPL 6200"},
    {"ancai06", 1, 2, 20, 40, 56, {82.02421478, 3.932024581, 0.7620913939, 4.5872, 4.562737221, 0.5130048107, 9.1164, 4.714215826, 0.7430370174, 4.456197317}, "RIPL 6200"},
    {"ancai06", 1, 2, 82, 208, 100, {76.99507508, 6.817830942, 0.7936549009, 7.324, 7.831814739, 0.2625007863, 7.77, 8.055939511, 0.8983495125, 7.720264754}, "RIPL 6200"},
    {"daehnick80", 1, 2, 20, 40, 56, {79.08627122, 4.001343715, 0.8042, 3.676040966, 4.531436259, 0.729396628, 9.979959034, 4.531436259, 0.729396628, 4.445937461}, "FRONT21"},
    {"daehnick80", 1, 2, 28, 58, 20, {89.66548314, 4.52892567, 0.743, 0.498758334, 5.128911549, 0.7862461872, 12.22124167, 5.128911549, 0.7862461872, 5.032139633}, "FRONT21"},
    {"daehnick80", 1, 2, 82, 208, 80, {79.87891912, 6.9322408, 0.845, 6.750264185, 7.850614581, 0.9047494496, 7.529735815, 7.850614581, 0.9047494496, 7.702489778}, "FRONT21"},
    {"daehnick80", 1, 2, 9, 19, 5.8, {89.96006892, 3.122029929, 0.71886, 0.04147828561, 3.535632185, 0.7020729378, 12.30932171, 3.535632185, 0.7020729378, 3.468922143}, "FRONT21"},
    {"kd03", 0, 1, 10, 20, 2.6, {53.37047263, 3.13392913, 0.674826, 0.3835113907, 3.13392913, 0.674826, 6.86586386, 3.52705041, 0.541288, 0}, "RIPL 2405"},
    {"kd03", 0, 1, 26, 56, 14, {47.45208266, 4.583141938, 0.6694728, 1.118259206, 4.583141938, 0.6694728, 7.168952993, 4.903837657, 0.5353264, 0}, "RIPL 2405"},
    {"kd03", 0, 1, 82, 208, 0.5, {47.7303101, 7.320197247, 0.6468704, 0.07572086866, 7.320197247, 0.6468704, 2.586696471, 7.397286765, 0.5101552, 0}, "RIPL 2405"},
    {"kd03", 0, 1, 82, 208, 150, {13.88561425, 7.320197247, 0.6468704, 11.74180123, 7.320197247, 0.6468704, 0.7592637623, 7.397286765, 0.5101552, 0}, "RIPL 2405"},
    {"kd03", 1, 1, 20, 40, 30, {46.5397075, 4.053875274, 0.671852, 2.946652563, 4.053875274, 0.671852, 6.417163507, 4.405560697, 0.53952, 4.395893121}, "RIPL 5405"},
    {"kd03", 1, 1, 82, 208, 65, {37.50515578, 7.320197247, 0.6468704, 6.33616321, 7.320197247, 0.6468704, 5.302106915, 7.397286765, 0.626964, 7.226321535}, "RIPL 5405"},
    {"kd03", 1, 1, 12, 24, 3, {56.31526791, 3.355698429, 0.6742312, 0.3142593949, 3.355698429, 0.6742312, 6.048115684, 3.740274332, 0.531192, 3.88496505}, "RIPL 5405"},
    {"bg71", 1, 3, 40, 90, 15, {161.7388889, 5.377685696, 0.72, 28.82777778, 6.273966645, 0.84, 0, 0, 0, 5.825826171}, "RIPL 7100"},
    {"bg71", 2, 3, 28, 58, 30, {148.5241379, 4.645051969, 0.72, 33.31724138, 5.419227297, 0.88, 0, 0, 0, 5.032139633}, "RIPL 8100"},
    {"liang09", 2, 3, 8, 16, 10, {118.0948739, 2.977479936, 0.740955184, 0, 0, 0, 21.53873315, 3.032617646, 0.794753626, 3.248076467}, "FRONT21"},
    {"liang09", 2, 3, 50, 120, 100, {108.3687635, 5.78982683, 0.8145389365, 17.4533, 6.958045141, 0.8812200889, 8.919704402, 5.87994698, 0.9532602666, 6.357894728}, "FRONT21"},
    {"mcfadden66", 2, 4, 8, 16, 20, {185, 3.52777894, 0.52, 25, 3.52777894, 0.52, 0, 0, 0, 3.52777894}, "RIPL 9100"},
    {"mcfadden66", 2, 4, 82, 208, 24.7, {185, 8.294988992, 0.52, 25, 8.294988992, 0.52, 0, 0, 0, 8.294988992}, "RIPL 9100"},
    {"avrigeanu94", 2, 4, 20, 40, 20, {131.5323107, 4.257840107, 0.7879304089, 10.80556207, 5.369324473, 0.6236009621, 0, 0, 0, 4.257840107}, "RIPL 9600"},
    {"avrigeanu94", 2, 4, 82, 208, 50, {172.4577483, 7.37661521, 0.7666375668, 12.53196341, 9.302237655, 0.5735001573, 0, 0, 0, 7.37661521}, "RIPL 9600"},
};

void TestRows() {
  static const char *names[10] = {"V", "R", "a", "W", "RW", "aW", "WD", "RD", "aD", "RC"};
  std::vector<int> points_per_row(GlobalOpticalSets().size(), 0);
  for (const Point &pt : kPoints) {
    const int index = GlobalOpticalIndex(pt.set);
    Check(index >= 0, std::string(pt.set) + " is a row");
    if (index < 0) continue;
    ++points_per_row[index];
    const GlobalOpticalSet &set = GlobalOpticalSets()[index];
    Check(Describes(set, pt.zp, pt.ap), std::string(pt.set) + " describes its projectile");
    const OpticalPotential o = set.at(pt.zp, pt.ap, pt.zt, pt.at, pt.e_lab);
    const WoodsSaxonTerm *terms[3] = {&o.real_volume, &o.imaginary_volume, &o.imaginary_surface};
    double got[10];
    for (int t = 0; t < 3; ++t) {
      got[3 * t] = terms[t]->depth_MeV;
      got[3 * t + 1] = RadiusOf(o, terms[t]->r0_fm);
      got[3 * t + 2] = terms[t]->diffuseness_fm;
    }
    got[9] = RadiusOf(o, o.coulomb_r0_fm);
    Check(o.real_spin_orbit.depth_MeV == 0.0 && o.imaginary_spin_orbit.depth_MeV == 0.0,
          std::string(pt.set) + ": central terms only");
    for (int t = 0; t < 10; ++t) {
      const double diff = std::fabs(got[t] - pt.p[t]);
      char what[160];
      std::snprintf(what, sizeof what, "%s (%d,%d) on (%d,%d) at %g MeV: %s = %.12g, %s %.12g", pt.set, pt.zp, pt.ap,
                    pt.zt, pt.at, pt.e_lab, names[t], got[t], pt.source, pt.p[t]);
      Check(diff <= 1e-9 * std::fabs(pt.p[t]) || diff <= 1e-12, what);
    }
  }
  for (std::size_t i = 0; i < points_per_row.size(); ++i)
    Check(points_per_row[i] >= 2, std::string(GlobalOpticalSets()[i].name) + " has at least two test points");
}

void TestLookup() {
  Check(GlobalOpticalNames() == "ancai06, daehnick80, kd03, bg71, liang09, mcfadden66, avrigeanu94", "the names");
  Check(GlobalOpticalIndex("kd04") == -1, "an unknown name");
  const GlobalOpticalSet &kd = GlobalOpticalSets()[GlobalOpticalIndex("kd03")];
  Check(Describes(kd, 0, 1) && Describes(kd, 1, 1) && !Describes(kd, 1, 2), "kd03: n and p, not d");
  const GlobalOpticalSet &ac = GlobalOpticalSets()[GlobalOpticalIndex("ancai06")];
  Check(Describes(ac, 1, 2) && !Describes(ac, 2, 4), "ancai06: d, not alpha");
  const GlobalOpticalSet &bg = GlobalOpticalSets()[GlobalOpticalIndex("bg71")];
  Check(Describes(bg, 1, 3) && Describes(bg, 2, 3) && !Describes(bg, 2, 4), "bg71: t and 3He");
  // Liang's volume absorption is negative below about 23 MeV: set to 0, and the term has no geometry then.
  const OpticalPotential l = GlobalOpticalSets()[GlobalOpticalIndex("liang09")].at(2, 3, 8, 16, 10.0);
  Check(l.imaginary_volume.depth_MeV == 0.0 && l.imaginary_volume.r0_fm == 0.0, "liang09: W < 0 set to 0");
}

}  // namespace

int main() {
  TestRows();
  TestLookup();
  if (failures) {
    std::printf("transfer_global_optical: %d check(s) failed\n", failures);
    return 1;
  }
  std::printf("transfer_global_optical: all checks passed\n");
  return 0;
}
