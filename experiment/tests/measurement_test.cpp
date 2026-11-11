// Unit test of what was measured (azure2/experiment/Measurement.h) and of the one lab -> c.m. conversion of its rows
// (ConvertRow, azure2/experiment/Frames.h), with values written by hand (docs/upgrade/data_structure.md 1.1, 1.2).
//
// 1. ReadMeasurementRows on a small file: comments and blank lines are no rows, the energy and (for differential
//    data) angle windows decide which rows are kept, records count every data row, columns 5-6 are the bins, the
//    hook replaces an error, the provenance counts rows read and kept and the kept energy range; a line that does not
//    parse and an order that is not an integer are refused with their line number.
// 2. ConvertRow against closed forms: elastic scattering of equal masses (theta_cm = 2 theta_lab, J = 1 / (4 cos
//    theta_lab)); an inelastic exit with the entrance's masses and Q = 0 (the reaction map gives the elastic angle);
//    a capture photon at 90 deg lab (cos theta_cm = -beta, J = 1 / (1 - beta^2), beta = p c / E of the c.m.); a
//    beta-delayed row (E_cm = E (m1 + m2) / m2 of the decay pair); relativistic on and off (E_cm from the invariant
//    mass); a ratio, whose value is not converted; a row beyond the kinematic maximum, at load (refused, value as read)
//    and after a shift (the load-time angle kept, the value converted there).
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>

#include "azure2/experiment/Frames.h"
#include "azure2/experiment/Measurement.h"
#include "azure2/physics/Constants.h"

using namespace azure2::experiment;

namespace {

int failures = 0;

void check(const std::string &name, bool ok, const std::string &detail = "") {
  std::printf("  %s  %s%s%s\n", ok ? "ok  " : "FAIL", name.c_str(), detail.empty() ? "" : "  ", detail.c_str());
  if (!ok) failures++;
}

bool near(double a, double b, double rel = 1e-12) { return std::fabs(a - b) <= rel * std::max(1.0, std::fabs(b)); }

std::string num(double x) {
  char b[64];
  std::snprintf(b, sizeof b, "%.17g", x);
  return b;
}

const double kDeg = 3.14159265358979323846 / 180.0;

void TestReader() {
  const std::string file =
      "# E   theta  y   dy\n"
      "1.0  30  10  1\n"
      "\n"
      "2.0  30  20  2  1.9 2.1\n"
      "3.0  95  30  3\n"
      "9.0  30  40  4\n"
      "  # trailing comment\n"
      "2.5  60  50  5\n";
  RowWindow w;
  w.minEnergy = 1.5;
  w.maxEnergy = 5.0;
  w.angles = true;
  w.minAngle = 0.0;
  w.maxAngle = 90.0;
  MeasurementRows rows;
  MeasurementProvenance prov;
  std::string error;
  std::istringstream in(file);
  int calls = 0;
  const bool ok = ReadMeasurementRows(in, "f.dat", w, ColumnConversion(), false,
                                      [&](DataRow &r) {
                                        calls++;
                                        if (r.record == 5) r.error = 0.5;
                                        return true;
                                      },
                                      rows, prov, error);
  check("reader: file read", ok, error);
  check("reader: rows kept 2 (E 2.0 and 2.5; 1.0, 9.0 outside E, 95 deg outside the angle window)",
        rows.Size() == 2 && rows.energy[0] == 2.0 && rows.energy[1] == 2.5);
  check("reader: records count data rows only (2 and 5)", rows.record == std::vector<int>({2, 5}));
  check("reader: bins from columns 5-6, NaN without", rows.binLow[0] == 1.9 && rows.binHigh[0] == 2.1 &&
                                                          std::isnan(rows.binLow[1]) && std::isnan(rows.binHigh[1]));
  check("reader: the hook sees kept rows only and replaces an error", calls == 2 && rows.error[1] == 0.5 &&
                                                                         rows.error[0] == 2.0);
  check("reader: provenance 5 read, 2 kept, E 2.0-2.5",
        prov.file == "f.dat" && prov.rowsRead == 5 && prov.rowsKept == 2 && prov.minEnergy == 2.0 &&
            prov.maxEnergy == 2.5);

  w.angles = false;
  MeasurementRows all;
  MeasurementProvenance p2;
  std::istringstream in2(file);
  ReadMeasurementRows(in2, "f.dat", w, ColumnConversion(), false, nullptr, all, p2, error);
  check("reader: without the angle window the 95 deg row is kept", all.Size() == 3);

  MeasurementRows bad;
  MeasurementProvenance p3;
  std::istringstream in3("1 2 3 4\n# c\n1 2 x 4\n");
  const bool refused = !ReadMeasurementRows(in3, "g.dat", w, ColumnConversion(), false, nullptr, bad, p3, error);
  check("reader: a line that does not parse is refused with its line number",
        refused && error.find("Cannot parse line 3 of g.dat") != std::string::npos, error);
  std::istringstream in4("2 1.5 3 4\n");
  const bool refusedOrder = !ReadMeasurementRows(in4, "h.dat", w, ColumnConversion(), true, nullptr, bad, p3, error);
  check("reader: a Legendre order 1.5 is refused", refusedOrder && error.find("got 1.5") != std::string::npos, error);
}

Partition Equal() {
  Partition p;
  p.m1 = 4.0;
  p.m2 = 4.0;
  return p;
}

FrameRule Rule(ExitKind kind, bool elastic, bool labDifferential = true, bool convertValue = true) {
  FrameRule r;
  r.entrance = Equal();
  r.exit = Equal();
  r.exitKind = kind;
  r.elastic = elastic;
  r.labDifferential = labDifferential;
  r.convertValue = convertValue;
  return r;
}

void TestConvert() {
  const double eLab = 2.0, thetaLab = 30.0;
  // Elastic, equal masses: theta_cm = 2 theta_lab, J = sin^2 tl / sin^2 tc cos(tc - tl) = 1 / (4 cos tl).
  {
    FrameRow row = UnconvertedRow(eLab, thetaLab, 1.9, 2.1);
    AngleLimit lim;
    const bool ok = ConvertRow(Rule(ExitKind::Particle, true), eLab, thetaLab, 1.9, 2.1, true, row, lim);
    const double j = 1.0 / (4.0 * std::cos(thetaLab * kDeg));
    check("elastic: E_cm = E_lab / 2, bins too", ok && near(row.eCM, 1.0) && near(row.binLowCM, 0.95) &&
                                                     near(row.binHighCM, 1.05), num(row.eCM));
    check("elastic: theta_cm = 2 theta_lab", near(row.thetaCM, 60.0), num(row.thetaCM));
    check("elastic: J = 1 / (4 cos theta_lab)", near(row.jacobian, j), num(row.jacobian));
  }
  // Inelastic with the entrance's masses and Q = 0: the reaction map gives the elastic angle and Jacobian.
  {
    FrameRow row = UnconvertedRow(eLab, thetaLab, NAN, NAN);
    AngleLimit lim;
    ConvertRow(Rule(ExitKind::Particle, false), eLab, thetaLab, NAN, NAN, true, row, lim);
    check("reaction, Q = 0, equal masses: theta_cm = 2 theta_lab", near(row.thetaCM, 60.0, 1e-10), num(row.thetaCM));
    check("reaction: c.m. bins NaN without bins", std::isnan(row.binLowCM) && std::isnan(row.binHighCM));
  }
  // Capture photon at 90 deg lab: cos theta_cm = -beta, J = 1 / (1 - beta^2).
  {
    const double u = azure2::physics::kAtomicMassUnit_MeV;
    const double p = std::sqrt(eLab * (eLab + 2.0 * 4.0 * u)), e = eLab + 8.0 * u, beta = p / e;
    FrameRow row = UnconvertedRow(eLab, 90.0, NAN, NAN);
    AngleLimit lim;
    ConvertRow(Rule(ExitKind::Photon, false), eLab, 90.0, NAN, NAN, true, row, lim);
    check("photon at 90 deg: cos theta_cm = -beta", near(std::cos(row.thetaCM * kDeg), -beta, 1e-9),
          num(row.thetaCM));
    check("photon at 90 deg: J = 1 / (1 - beta^2)", near(row.jacobian, 1.0 / (1.0 - beta * beta)), num(row.jacobian));
  }
  // Beta-delayed: E is the light particle's; E_cm = E (m1 + m2) / m2 of the decay (exit) pair, E_x = E_cm + S + E*.
  {
    FrameRule r = Rule(ExitKind::Particle, false, false);
    r.betaDelayed = true;
    r.exit.m1 = 4.0;
    r.exit.m2 = 12.0;
    r.exit.sepE = 7.0;
    r.exit.exE = 0.5;
    FrameRow row = UnconvertedRow(1.2, 0.0, NAN, NAN);
    AngleLimit lim;
    ConvertRow(r, 1.2, 0.0, NAN, NAN, true, row, lim);
    check("beta-delayed: E_cm = E (m1 + m2) / m2, E_x = E_cm + S + E*",
          near(row.eCM, 1.2 * 16.0 / 12.0) && near(row.eX, 1.2 * 16.0 / 12.0 + 7.5), num(row.eCM));
    check("beta-delayed (not lab differential): angle as read, J 1", row.thetaCM == 0.0 && row.jacobian == 1.0);
  }
  // Relativistic: E_cm = sqrt(s) - m1 - m2 with s = (m1 + m2)^2 + 2 m2 E_lab (MeV); off: E_lab m2 / (m1 + m2).
  {
    const double u = azure2::physics::kAtomicMassUnit_MeV, m = 4.0 * u;
    const double eRel = std::sqrt(4.0 * m * m + 2.0 * m * eLab) - 2.0 * m;
    FrameRule r = Rule(ExitKind::Particle, true, false);
    r.relativistic = true;
    FrameRow row = UnconvertedRow(eLab, thetaLab, NAN, NAN);
    AngleLimit lim;
    ConvertRow(r, eLab, thetaLab, NAN, NAN, true, row, lim);
    check("relativistic on: E_cm from the invariant mass", near(row.eCM, eRel, 1e-9), num(row.eCM));
    r.relativistic = false;
    row = UnconvertedRow(eLab, thetaLab, NAN, NAN);
    ConvertRow(r, eLab, thetaLab, NAN, NAN, true, row, lim);
    check("relativistic off: E_cm = E_lab m2 / (m1 + m2)", near(row.eCM, 1.0) && !near(row.eCM, eRel, 1e-9),
          num(row.eCM));
  }
  // A ratio: the angle is converted, the value is not.
  {
    FrameRow row = UnconvertedRow(eLab, thetaLab, NAN, NAN);
    AngleLimit lim;
    ConvertRow(Rule(ExitKind::Particle, true, true, false), eLab, thetaLab, NAN, NAN, true, row, lim);
    check("ratio: theta_cm converted, J 1 (the value is dimensionless)", near(row.thetaCM, 60.0) &&
                                                                              row.jacobian == 1.0);
  }
  // Beyond the kinematic maximum: elastic with m1 = 2 m2, theta_max = asin(m2/m1) = 30 deg.
  {
    FrameRule r = Rule(ExitKind::Particle, true);
    r.entrance.m1 = 8.0;
    r.exit.m1 = 8.0;
    FrameRow load = UnconvertedRow(eLab, 20.0, NAN, NAN);
    AngleLimit lim;
    const bool okLoad = ConvertRow(r, eLab, 20.0, NAN, NAN, true, load, lim);
    check("kinematic maximum: 20 deg converted at load", okLoad && load.thetaCM > 20.0);
    FrameRow beyond = UnconvertedRow(eLab, 40.0, NAN, NAN);
    AngleLimit limLoad;
    const bool refused = !ConvertRow(r, eLab, 40.0, NAN, NAN, true, beyond, limLoad);
    check("kinematic maximum: 40 deg refused at load, J 1, maximum 30 deg reported",
          refused && beyond.jacobian == 1.0 && limLoad.maxSet && near(limLoad.max, 30.0) && !limLoad.why.empty(),
          num(limLoad.max));
    FrameRow shifted = load;  // after a shift a row starts from its load-time conversion
    AngleLimit limShift;
    const bool inRange = ConvertRow(r, 2.5, 40.0, NAN, NAN, false, shifted, limShift);
    check("kinematic maximum after a shift: false, the load-time angle kept, J computed there",
          !inRange && shifted.thetaCM == load.thetaCM && near(shifted.eCM, 2.5 * 4.0 / 12.0) &&
              shifted.jacobian == ParticleJacobian(r.entrance, r.exit, 2.5, 40.0, load.thetaCM));
  }
}

}  // namespace

int main() {
  TestReader();
  TestConvert();
  if (failures) {
    std::printf("%d FAILURE(S)\n", failures);
    return 1;
  }
  std::printf("all passed\n");
  return 0;
}
