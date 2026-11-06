// The EXFOR client's ratio-to-Rutherford conversion (ExforData::parseCsv, quantity ,,RTH) writes a
// cross section in the frame of the angle column: for a lab ANG column the lab angle is taken to the
// centre of mass with the elastic kinematics and the c.m. Rutherford value is divided by the lab-to-c.m.
// factor AZURE2 applies to a lab data value, so the engine reads back ratio x sigma_R,cm(theta_cm).
// Until 2026-10-27 it used the c.m. Rutherford at the lab angle (audit item 8; pyazr.nds, same bug,
// fixed 2026-10-26).  An ANG-CM column is unchanged.  Headless: no network, parseCsv on a canned CSV.
#include <QCoreApplication>
#include <QString>
#include <cmath>
#include <iostream>
#include "ExforData.h"
#include "Config.h"

// stubs the core library expects from the application (as the other GUI tests)
struct SegPairs { int firstPair; int secondPair; };
Config *g_config = nullptr;
void exitMessage(const Config &) {}
bool checkExternalCapture(Config &, const std::vector<SegPairs> &) { return true; }
bool readSegmentFile(const Config &, std::vector<SegPairs> &) { return true; }
void startMessage(const Config &) {}

static int fails = 0;
static void ok(const QString &what, bool cond, const QString &detail = QString()) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what.toStdString();
  if (!cond && !detail.isEmpty()) std::cout << "  -- " << detail.toStdString();
  std::cout << std::endl;
  if (!cond) fails++;
}
static bool close(double a, double b, double tol) { return std::fabs(a - b) <= tol * std::max(1.0, std::fabs(b)); }

int main(int argc, char **argv) {
  QCoreApplication app(argc, argv);
  const double kPi = 3.141592650, kFstruc = 1.0 / 137.0359996790, kHbarc = 197.32696310;
  const double z1 = 1, z2 = 6, m1 = 1.00728, m2 = 12.0;   // p + 12C
  const double eLab = 2.0, eCm = eLab * m2 / (m1 + m2);
  auto rutherford = [&](double thetaDeg) {
    const double s = std::sin(thetaDeg * kPi / 360.0);
    const double amp = z1 * z2 * kFstruc * kHbarc / (4.0 * eCm);
    return amp * amp / (s * s * s * s) * 1.0e-2;
  };
  const double angles[3] = {30.0, 90.0, 150.0}, ratio[3] = {0.9, 1.2, 1.5}, err[3] = {0.02, 0.03, 0.04};

  ExforData d;
  d.setEntranceChannel(z1, z2, m1, m2);
  QString error;
  bool differential = false;

  // --- ANG-CM: ratio x sigma_R,cm(theta)
  QList<ExforPoint> cm;
  QString csvCm = "EN (MEV),ANG-CM (ADEG),DATA (NO-DIM),ERR-S (NO-DIM)\n";
  for (int i = 0; i < 3; i++) csvCm += QString("%1,%2,%3,%4\n").arg(eLab).arg(angles[i]).arg(ratio[i]).arg(err[i]);
  ok("ANG-CM CSV parses", d.parseCsv(csvCm, cm, differential, error) && cm.size() == 3 && differential, error);
  for (int i = 0; i < 3 && cm.size() == 3; i++)
    ok(QString("ANG-CM %1 deg: ratio x sigma_R,cm(theta)").arg(angles[i]),
       close(cm[i].cross, ratio[i] * rutherford(angles[i]), 1e-9) && close(cm[i].angle, angles[i], 1e-12),
       QString("%1 vs %2").arg(cm[i].cross).arg(ratio[i] * rutherford(angles[i])));

  // --- ANG (lab): ratio x sigma_R,cm(theta_cm) / J, angle column written through
  QList<ExforPoint> lab;
  QString csvLab = csvCm;
  csvLab.replace("ANG-CM (ADEG)", "ANG (ADEG)");
  ok("ANG CSV parses", d.parseCsv(csvLab, lab, differential, error) && lab.size() == 3, error);
  const double g = m1 / m2;
  for (int i = 0; i < 3 && lab.size() == 3; i++) {
    const double t = angles[i] * kPi / 180.0;
    const double tc = t + std::asin(g * std::sin(t));
    const double c = std::cos(tc);
    const double jac = std::fabs(1.0 + g * c) / std::pow(1.0 + g * g + 2.0 * g * c, 1.5);
    const double want = ratio[i] * rutherford(tc * 180.0 / kPi) / jac;
    ok(QString("ANG %1 deg: ratio x sigma_R,cm(theta_cm) / J").arg(angles[i]),
       close(lab[i].cross, want, 1e-9) && close(lab[i].angle, angles[i], 1e-12) && close(lab[i].error, err[i] * want / ratio[i], 1e-9),
       QString("%1 vs %2").arg(lab[i].cross).arg(want));
    ok(QString("ANG %1 deg: the old value (c.m. Rutherford at the lab angle) is off").arg(angles[i]),
       std::fabs(ratio[i] * rutherford(angles[i]) / want - 1.0) > 0.04);
  }

  // --- a lab column without masses is refused
  ExforData noMass;
  noMass.setEntranceChannel(z1, z2, 0.0, 0.0);
  QList<ExforPoint> none;
  ok("ANG without the entrance masses is refused",
     !noMass.parseCsv(csvLab, none, differential, error) && error.contains("masses"), error);
  ok("ANG-CM without masses still converts", noMass.parseCsv(csvCm, none, differential, error) && none.size() == 3, error);

  std::cout << (fails ? "FAILED" : "ALL PASSED") << " (" << fails << " failures)" << std::endl;
  return fails ? 1 : 0;
}
