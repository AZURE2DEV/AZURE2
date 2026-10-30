// Headless round-trip checks on the Experimental Effects tab's optional
// trailing tokens.
//
// A targetInt line may end with a quoted lab-energy range list, a blend width
// and an automatic-application tolerance, and after those with a beam-profile
// block introduced by the `beamprofile` keyword.  The GUI has to read them,
// carry them in its model, and write them back -- but only when they are used,
// so a file that never touched either feature stays byte-identical.  Both
// directions are checked here through the tab's own readFile/writeFile.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QTextStream>
#include <QString>
#include <cmath>
#include <iostream>
#include "TargetIntTab.h"
#include "AddTargetIntDialog.h"
#include <QCheckBox>
#include <QLineEdit>
#include <QRadioButton>
#include <QHeaderView>
#include <QTableWidget>
#include <QSpinBox>
#include <cstdlib>
#include "TargetIntModel.h"
#include "ElectronBinding.h"
#include "Config.h"
struct SegPairs {int firstPair; int secondPair;};

// Defined by AZURE2.cpp, which belongs to the executable rather than the GUI
// library, so this test supplies its own. They are never called from here.
Config* g_config = nullptr;
void exitMessage(const Config&) {}
bool checkExternalCapture(Config&, const std::vector<SegPairs>&) { return true; }
bool readSegmentFile(const Config&, std::vector<SegPairs>&) { return true; }
void startMessage(const Config&) {}

static int fails = 0;
static void ok(const char* what, bool cond, const QString& detail = QString()) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what;
  if(!cond && !detail.isEmpty()) std::cout << "  -- " << detail.toStdString();
  std::cout << std::endl;
  if(!cond) fails++;
}

static QString roundTrip(TargetIntTab& tab, const QString& block, bool& readOk) {
  tab.reset();
  QString in(block);
  QTextStream inStream(&in);
  readOk = tab.readFile(inStream);
  QString out;
  QTextStream outStream(&out);
  tab.writeFile(outStream);
  return out;
}

int main(int argc, char** argv) {
  QApplication app(argc, argv);
  TargetIntTab tab;

  // 1. A line with the new tokens survives a read/write round trip.
  bool readOk = false;
  QString withTokens =
      "1 \"1,2\" 50 1 0.030 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 5 50 "
      "\"1.95-2.55,2.8-3.0\" 0.12 0.002\n</targetInt>\n";
  QString out = roundTrip(tab, withTokens, readOk);
  ok("read a line carrying ranges/blend/tolerance", readOk);
  ok("ranges token written back", out.contains("\"1.95-2.55,2.8-3.0\""), out);
  ok("blend width written back", out.contains(" 0.12 "), out);
  ok("tolerance written back", out.contains(" 0.002"), out);

  QList<TargetIntData> lines = tab.getTargetIntModel()->getLines();
  ok("model row present", lines.size() == 1);
  if (lines.size() == 1) {
    ok("model carries ranges", lines.at(0).applyRanges == "1.95-2.55,2.8-3.0",
       lines.at(0).applyRanges);
    ok("model carries blend width", lines.at(0).transitionWidth == 0.12);
    ok("model carries tolerance", lines.at(0).autoTolerance == 0.002);
  }

  // 2. Reading the GUI's own output again reproduces the same fields.
  bool readOk2 = false;
  QString out2 = roundTrip(tab, out + "</targetInt>\n", readOk2);
  ok("re-read the written line", readOk2);
  ok("second round trip stable", out2 == out, out2);

  // 3. A legacy line without the tokens gains nothing on write: files that
  //    never used the feature stay byte-identical in this respect.
  QString legacy = "1 \"1,2\" 50 1 0.030 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 5 50\n</targetInt>\n";
  QString out3 = roundTrip(tab, legacy, readOk);
  ok("read a legacy line", readOk);
  // Exactly the six quotes of the three legacy string fields; a written
  // ranges token would add two more.
  ok("no ranges token invented", out3.count('"') == 6, out3);
  ok("line ends at points-per-width", out3.trimmed().endsWith("50"), out3);

  // 4. Tolerance-only form: empty ranges token, zero width, tolerance set.
  QString tolOnly = "1 \"1\" 50 1 0.030 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 5 50 \"\" 0 0.005\n</targetInt>\n";
  QString out4 = roundTrip(tab, tolOnly, readOk);
  ok("read a tolerance-only line", readOk);
  lines = tab.getTargetIntModel()->getLines();
  ok("empty ranges stay empty", lines.size() == 1 && lines.at(0).applyRanges.isEmpty());
  ok("tolerance-only round trip keeps tolerance", out4.contains(" 0.005"), out4);

  // 5. A beam-profile block with no ranges tokens before it: the keyword is
  //    the first non-numeric thing on the line after the grid parameters.
  //    This is the shape pyazr writes, so it is the one that matters most.
  QString beamOnly =
      "1 \"7\" 150 0 0 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 20 50 "
      "beamprofile 1 3.7802596656384173 0.3760611 -2.06 1 0.07334525 0 1\n</targetInt>\n";
  QString out5 = roundTrip(tab, beamOnly, readOk);
  ok("read a beam-profile line", readOk);
  lines = tab.getTargetIntModel()->getLines();
  ok("beam-profile row present", lines.size() == 1, out5);
  if (lines.size() == 1) {
    const TargetIntData& d = lines.at(0);
    ok("model flags the beam profile", d.isBeamProfile);
    ok("one component, four numbers", d.beamProfile.size() == 4,
       QString::number(d.beamProfile.size()));
    if (d.beamProfile.size() == 4) {
      // 12 significant digits survive the write; the engine's own parser reads
      // the same text, so anything coarser would move the profile on a save.
      ok("location kept to 1e-9", std::fabs(d.beamProfile.at(0) / 3.7802596656384173 - 1.) < 1e-9,
         QString::number(d.beamProfile.at(0), 'g', 17));
      ok("scale kept", d.beamProfile.at(1) == 0.3760611);
      ok("skewness kept", d.beamProfile.at(2) == -2.06);
      ok("weight kept", d.beamProfile.at(3) == 1.);
    }
    ok("resolution sigma kept", d.beamTpcSigma == 0.07334525);
    ok("truncation kept", d.beamTruncation == 0.);
    ok("detailed-balance flag kept", d.beamPhotodissociation);
    // The grid parameters sit just before the keyword and must not be eaten.
    ok("width multiplier still read", d.resonanceWidthMultiplier == 20.);
    ok("points per width still read", d.pointsPerWidth == 50.);
    ok("no ranges invented", d.applyRanges.isEmpty());
  }
  ok("keyword written back", out5.contains("beamprofile 1 "), out5);
  ok("skewness written back", out5.contains("-2.06"), out5);
  ok("no ranges token invented alongside", out5.count('"') == 6, out5);
  bool readOk5 = false;
  QString out5b = roundTrip(tab, out5 + "</targetInt>\n", readOk5);
  ok("re-read the written beam-profile line", readOk5);
  ok("beam-profile round trip stable", out5b == out5, out5b);

  // 6. Several components, and ranges tokens in front of the keyword: both
  //    optional blocks on one line, in the order the engine parses them.
  QString beamMulti =
      "1 \"3-5\" 200 0 0 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 10 25 \"1.9-2.6\" 0.1 0.001 "
      "beamprofile 2 2.5 0.3 -1.5 0.4 2.9 0.31 -1.7 0.6 0.055 2 0\n</targetInt>\n";
  QString out6 = roundTrip(tab, beamMulti, readOk);
  ok("read a two-component line with ranges", readOk);
  lines = tab.getTargetIntModel()->getLines();
  if (lines.size() == 1) {
    const TargetIntData& d = lines.at(0);
    ok("two components, eight numbers", d.beamProfile.size() == 8,
       QString::number(d.beamProfile.size()));
    ok("second component read", d.beamProfile.size() == 8 && d.beamProfile.at(4) == 2.9
                                   && d.beamProfile.at(7) == 0.6);
    ok("truncation read", d.beamTruncation == 2.);
    ok("detailed balance off", !d.beamPhotodissociation);
    ok("ranges still read alongside", d.applyRanges == "1.9-2.6", d.applyRanges);
    ok("blend width still read", d.transitionWidth == 0.1);
  }
  ok("both blocks written back", out6.contains("\"1.9-2.6\"") && out6.contains("beamprofile 2 "), out6);
  bool readOk6 = false;
  QString out6b = roundTrip(tab, out6 + "</targetInt>\n", readOk6);
  ok("two-component round trip stable", readOk6 && out6b == out6, out6b);

  // 7. A block that declares more components than it supplies is dropped
  //    rather than half-applied: the engine would otherwise fold the data
  //    with a profile the tab never showed.
  QString beamShort =
      "1 \"1\" 50 0 0 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 5 50 "
      "beamprofile 2 2.5 0.3 -1.5 0.4\n</targetInt>\n";
  QString out7 = roundTrip(tab, beamShort, readOk);
  lines = tab.getTargetIntModel()->getLines();
  ok("truncated block leaves the profile off",
     lines.size() == 1 && !lines.at(0).isBeamProfile && lines.at(0).beamProfile.isEmpty());
  ok("truncated block is not written back", !out7.contains("beamprofile"), out7);

  // 8. A legacy line still gains no beam-profile token (case 3 checks the
  //    ranges half of the same promise).
  QString out8 = roundTrip(tab, legacy, readOk);
  ok("no beam-profile token invented", !out8.contains("beamprofile"), out8);

  // 9. A user-defined resolution block: quoted path with a space, flight
  //    path, burst, channel width and the centring flag; alone and after a
  //    beam-profile block on the same line.
  QString udrLine =
      "1 \"4\" 200 0 0 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 20 50 "
      "udr \"data/RF EAR1_v2.txt\" 183.5 7 0 1\n</targetInt>\n";
  QString out9 = roundTrip(tab, udrLine, readOk);
  ok("read a line carrying a UDR block", readOk);
  lines = tab.getTargetIntModel()->getLines();
  if (lines.size() == 1) {
    const TargetIntData& d = lines.at(0);
    ok("UDR flag read", d.isUdr);
    ok("UDR path with a space kept", d.udrFile == "data/RF EAR1_v2.txt", d.udrFile);
    ok("flight path read", d.udrFlightPath == 183.5);
    ok("burst read", d.udrBurstFwhm == 7.);
    ok("channel width read", d.udrChannelWidth == 0.);
    ok("centring read", d.udrCentred);
    ok("no beam profile invented", !d.isBeamProfile);
    ok("grid parameters before the keyword still read", d.resonanceWidthMultiplier == 20. && d.pointsPerWidth == 50.);
  }
  ok("UDR block written back", out9.contains("udr \"data/RF EAR1_v2.txt\" 183.5 7 0 1"), out9);
  bool readOk9 = false;
  QString out9b = roundTrip(tab, out9 + "</targetInt>\n", readOk9);
  ok("UDR round trip stable", readOk9 && out9b == out9, out9b);
  ok("no UDR token on a legacy line", !out8.contains(" udr "), out8);

  QString bothBlocks =
      "1 \"3\" 200 0 0 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 20 50 "
      "beamprofile 1 2.5 0.3 -1.5 1 0.07 0 1 udr nTOF.txt 185 0 3.2 0\n</targetInt>\n";
  QString out10 = roundTrip(tab, bothBlocks, readOk);
  ok("read beam-profile and UDR blocks together", readOk);
  lines = tab.getTargetIntModel()->getLines();
  if (lines.size() == 1) {
    const TargetIntData& d = lines.at(0);
    ok("beam profile still read beside UDR", d.isBeamProfile && d.beamProfile.size() == 4);
    ok("unquoted UDR path read", d.isUdr && d.udrFile == "nTOF.txt", d.udrFile);
    ok("uncentred flag read", !d.udrCentred);
    ok("channel width read beside burst zero", d.udrChannelWidth == 3.2 && d.udrBurstFwhm == 0.);
  }
  ok("both blocks written back in order",
     out10.indexOf("beamprofile 1 ") >= 0 && out10.indexOf(" udr \"nTOF.txt\" 185 0 3.2 0") > out10.indexOf("beamprofile 1 "), out10);

  // 11. One Gaussian beam spread: an expression written in the lab frame carries the sigmalab keyword and
  // keeps it through a round trip; a legacy line without it stays without it (and is read as c.m.).
  QString labLine = "1 \"4\" 60 0 0 0 0 \"\" 0 0 0 1 \"a0*x^1.5\" 1 0.0092 0 0.04 20 50 sigmalab\n</targetInt>\n";
  QString out11 = roundTrip(tab, labLine, readOk);
  ok("read a sigmalab line", readOk);
  lines = tab.getTargetIntModel()->getLines();
  if (lines.size() == 1) {
    ok("expression read with the lab flag", lines.at(0).isConvCoefficients && lines.at(0).sigmaLab && lines.at(0).convolutionEq == "a0*x^1.5", lines.at(0).convolutionEq);
    ok("constant-sigma flag off for an expression line", !lines.at(0).isConvolution);
  }
  ok("sigmalab written back, before any keyword block", out11.contains(" 0.04 20 50 sigmalab") || out11.contains(" 50 sigmalab"), out11);
  QString legacyLine = "1 \"4\" 60 0 0 0 0 \"\" 0 0 0 1 \"a0*x^1.5\" 1 0.0092 0 0.04 20 50\n</targetInt>\n";
  QString out12 = roundTrip(tab, legacyLine, readOk);
  lines = tab.getTargetIntModel()->getLines();
  ok("legacy expression line read as c.m.", readOk && lines.size() == 1 && lines.at(0).isConvCoefficients && !lines.at(0).sigmaLab);
  ok("no sigmalab added to a legacy line", !out12.contains("sigmalab"), out12);
  QString ranged = "1 \"4\" 60 0 0 0 0 \"\" 0 0 0 1 \"a0\" 1 0.03 0 0.04 20 50 \"1.5-2.5\" 0.1 0 sigmalab beamprofile 1 2.5 0.3 -1.5 1 0.07 0 1\n</targetInt>\n";
  QString out13 = roundTrip(tab, ranged, readOk);
  lines = tab.getTargetIntModel()->getLines();
  ok("sigmalab read between a ranges block and a beam profile", readOk && lines.size() == 1 && lines.at(0).sigmaLab && lines.at(0).isBeamProfile && lines.at(0).applyRanges == "1.5-2.5",
     lines.size() == 1 ? lines.at(0).applyRanges : QString("(no line)"));
  ok("and written back in that order", out13.indexOf("\"1.5-2.5\"") < out13.indexOf("sigmalab") && out13.indexOf("sigmalab") < out13.indexOf("beamprofile"), out13);

  // 12. The dialog: one checkbox, a constant or an expression; the legacy c.m. choice converts exactly.
  {
    AddTargetIntDialog dialog;
    dialog.isConvolutionCheck->setChecked(true);
    dialog.sigmaText->setText("0.030");
    ok("a number in the width field is a constant", dialog.sigmaIsConstant());
    dialog.sigmaText->setText("a0*x^1.5+exp(-x)");
    ok("an expression in the width field is not", !dialog.sigmaIsConstant());
    dialog.sigmaCmRadio->setChecked(true);
    dialog.setConversionFactor(0.5);
    dialog.convertExpressionToLab();
    ok("c.m. -> lab rewrite substitutes x and divides, leaving exp alone",
       dialog.sigmaText->text() == "(a0*(0.5*x)^1.5+exp(-(0.5*x)))/0.5", dialog.sigmaText->text());
    ok("after the rewrite the line is marked lab", dialog.sigmaLabRadio->isChecked());
    dialog.numConvCoefficientSpin->setValue(5);
    QApplication::processEvents();
    ok("five coefficient rows all visible without scrolling",
       dialog.convCoefficientTable->height() >= dialog.convCoefficientTable->horizontalHeader()->height() + 5 * dialog.convCoefficientTable->rowHeight(0),
       QString::number(dialog.convCoefficientTable->height()));
    if (const char* stem = std::getenv("AZURE2_TARGET_INT_PNG")) {
      dialog.show(); QApplication::processEvents();
      dialog.grab().save(QString("%1_expression.png").arg(stem));
    }
  }

  // Electron screening: the block round-trips wherever it sits (pyazr writes it after udr), is written
  // before the blocks that run to the end of the line, and the adiabatic limit follows NIST binding energies.
  {
    QString scrLine =
        "1 \"6,7\" 200 1 0.002 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 20 50 "
        "udr \"res file.udr\" 185 7 1.6 1 screening 934.7 1 200\n</targetInt>\n";
    bool readS = false;
    QString outS = roundTrip(tab, scrLine, readS);
    ok("read a line with screening after udr", readS);
    QList<TargetIntData> sl = tab.getTargetIntModel()->getLines();
    ok("screening row present", sl.size() == 1, outS);
    if (sl.size() == 1) {
      const TargetIntData &d = sl.at(0);
      ok("screening flagged, U_e 934.7, fitted, prior 200",
         d.isScreening && d.screeningUe == 934.7 && d.screeningVary && d.screeningSigma == 200.);
      ok("udr still read beside it", d.isUdr && d.udrFile == "res file.udr" && d.udrCentred);
      ok("Gaussian width still read", d.isConvolution && d.sigma == 0.002);
    }
    ok("screening written before udr", outS.contains(" screening 934.7 1 200 udr "), outS);
    bool readS2 = false;
    QString outS2 = roundTrip(tab, outS + "</targetInt>\n", readS2);
    ok("second round trip identical", readS2 && outS2 == outS, outS2);
    QString noScr = roundTrip(tab, legacy, readS);
    ok("a line without screening gains no block", !noScr.contains("screening"), noScr);
    ok("adiabatic limit alpha + 13C = 934.73 eV", std::fabs(AdiabaticScreeningEv(2, 6) - 934.73) < 0.01,
       QString::number(AdiabaticScreeningEv(2, 6), 'f', 3));
    ok("adiabatic limit p + 12C = B(N) - B(H) - B(C) = 442.35 eV", std::fabs(AdiabaticScreeningEv(1, 6) - (1486.058 - 13.598 - 1030.108)) < 1e-9);
    ok("no adiabatic limit for a neutron", AdiabaticScreeningEv(0, 6) < 0.0);
  }

  // Fitted target density: densityfit / densitygroup read into the line's fields, written back from them, edited in
  // the dialog (box shown with target integration, the default grid maximum, refusals), summarised in the table.
  {
    QString denLine =
        "1 \"18\" 50 1 0.0009 1 4.09554e+17 \"a1*x+a0\" 2 6.1e-20 1.2e-19 0 0 0 \"\" 0 1 0.7489 20 50 "
        "screening 0 0 0 densityfit 1 2.3% 4.5e17 densitygroup t3b\n"
        "1 \"19\" 50 1 0.0009 1 4.09554e+17 \"a1*x+a0\" 2 6.1e-20 1.2e-19 0 0 0 \"\" 0 1 0.7489 20 50 densitygroup t3b\n"
        "</targetInt>\n";
    bool readD = false;
    QString outD = roundTrip(tab, denLine, readD);
    ok("read a line with densityfit and densitygroup, and a group member", readD);
    QList<TargetIntData> dl = tab.getTargetIntModel()->getLines();
    ok("two rows", dl.size() == 2, outD);
    if (dl.size() == 2) {
      const TargetIntData &d = dl.at(0), &m = dl.at(1);
      ok("leader: parameter, varied, prior 2.3%, maximum 4.5e17, group t3b",
         d.densityFit && d.densityVary && d.densitySigma == "2.3%" && d.densityMax == 4.5e17 && d.densityGroup == "t3b");
      ok("straggling, screening and the density still read beside them",
         d.isStraggling && d.stragglingCoefficient == 0.7489 && d.isScreening && d.density == 4.09554e+17);
      ok("member: the label only", !m.densityFit && m.densityGroup == "t3b" && m.densitySigma == "0" && m.densityMax == 0.);
      const QString s0 = tab.getTargetIntModel()->data(tab.getTargetIntModel()->index(0, 33), Qt::DisplayRole).toString();
      const QString s1 = tab.getTargetIntModel()->data(tab.getTargetIntModel()->index(1, 33), Qt::DisplayRole).toString();
      ok("table summary names the fitted density and the group", s0.contains("density fitted (group t3b)") &&
         s1.contains("density (group t3b)") && !s1.contains("fitted"), s0 + " | " + s1);
    }
    ok("leader written back", outD.contains(" densityfit 1 2.3% 4.5e+17 densitygroup t3b\n") ||
                              outD.contains(" densityfit 1 2.3% 4.5e17 densitygroup t3b\n"), outD);
    ok("member written with the label only", outD.contains("20 50 densitygroup t3b\n") &&
                                             outD.count("densityfit") == 1, outD);
    bool readD2 = false;
    QString outD2 = roundTrip(tab, outD + "</targetInt>\n", readD2);
    ok("second round trip identical", readD2 && outD2 == outD, outD2);
    QString noDen = roundTrip(tab, legacy, readD);
    ok("a line without them gains no block", !noDen.contains("density"), noDen);

    // A line built in the tab (TargetIntTab::addLine(TargetIntData)) is written with its blocks.
    tab.reset();
    TargetIntData n;
    n.isActive = 1; n.segmentsList = "3"; n.numPoints = 50; n.isConvolution = false; n.sigma = 0.;
    n.isTargetIntegration = true; n.density = 1e18; n.stoppingPowerEq = "a0"; n.numParameters = 1; n.parameters << 1e-20;
    n.isQCoefficients = false; n.isConvCoefficients = false; n.isStraggling = false; n.stragglingCoefficient = 0.;
    n.resonanceWidthMultiplier = 20.; n.pointsPerWidth = 50.; n.transitionWidth = 0.; n.autoTolerance = 0.;
    n.densityFit = true; n.densityVary = true; n.densitySigma = "5e16";
    tab.addLine(n);
    QString built;
    QTextStream builtStream(&built);
    tab.writeFile(builtStream);
    ok("a new line: densityfit 1 5e16, no maximum, no group", built.contains(" densityfit 1 5e16\n") &&
                                                               !built.contains("densitygroup"), built);

    // The dialog
    AddTargetIntDialog dlg;
    ok("density box hidden without target integration", !dlg.findChild<QCheckBox*>("densityFitCheck")->isVisibleTo(&dlg));
    dlg.isTargetIntegrationCheck->setChecked(true);
    dlg.densityText->setText("4.09554e17");
    ok("shown with it", dlg.findChild<QCheckBox*>("densityFitCheck")->isVisibleTo(&dlg));
    dlg.setDensityFit(false, false, "0", 0., "");
    ok("prior and maximum disabled until the density is a parameter",
       !dlg.densitySigmaText->isEnabled() && !dlg.densityMaxText->isEnabled() && dlg.densityGroupText->isEnabled());
    dlg.setDensityFit(true, true, "2.3%", 0., "t3b");
    ok("enabled with it", dlg.densitySigmaText->isEnabled() && dlg.densityVaryCheck->isEnabled());
    ok("default maximum shown: max(N + 3 sigma, 1.1 N) = 4.50509e+17",
       dlg.densityMaxText->placeholderText() == "4.50509e+17", dlg.densityMaxText->placeholderText());
    dlg.densitySigmaText->setText("16%");
    ok("16 %: N + 3 sigma = 6.0614e+17", dlg.densityMaxText->placeholderText() == "6.0614e+17", dlg.densityMaxText->placeholderText());
    ok("valid settings pass", dlg.densityFitProblem().isEmpty(), dlg.densityFitProblem());
    dlg.densitySigmaText->setText("2.3 %x");
    ok("a malformed sigma is refused", !dlg.densityFitProblem().isEmpty());
    dlg.densitySigmaText->setText("-1");
    ok("a negative sigma is refused", !dlg.densityFitProblem().isEmpty());
    dlg.densitySigmaText->setText("1e16");
    dlg.densityMaxText->setText("4e17");
    ok("a maximum below the density is refused", !dlg.densityFitProblem().isEmpty());
    dlg.densityMaxText->setText("5e17");
    ok("a maximum above it passes", dlg.densityFitProblem().isEmpty(), dlg.densityFitProblem());
    dlg.densityText->setText("0");
    ok("a thick target (density 0) is refused", !dlg.densityFitProblem().isEmpty());
    dlg.densityText->setText("4.09554e17");
    dlg.isBeamProfileCheck->setChecked(true);
    ok("with a beam profile it is refused", !dlg.densityFitProblem().isEmpty());
    dlg.isBeamProfileCheck->setChecked(false);
    dlg.isDensityFitCheck->setChecked(false);
    ok("a group member (label only) passes", dlg.densityFitProblem().isEmpty(), dlg.densityFitProblem());
    if (const char* stem = std::getenv("AZURE2_TARGET_INT_PNG")) {
      dlg.isDensityFitCheck->setChecked(true);
      dlg.show(); QApplication::processEvents(); dlg.adjustSize(); QApplication::processEvents();
      dlg.grab().save(QString("%1_density.png").arg(stem));
    }
  }

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
