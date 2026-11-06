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
#include <QDoubleSpinBox>
#include <QGroupBox>
#include <QTableView>
#include <QItemSelectionModel>
#include <QTimer>
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

  // 13. A sigma or density whose round-trip text fills the 15-wide column (17
  //     digits, as Python writes 0.1 + 0.2) keeps a space before the next field.
  QString longNumbers =
      "1 \"1\" 50 1 0.30000000000000004 1 0.12345678901234566 \"\" 0 0 0 0 \"\" 0 0 0.04 5 50\n</targetInt>\n";
  // Left-aligned, as AZURESetup::writeConfig leaves the project stream.
  QString outLong;
  {
    tab.reset();
    QString in(longNumbers);
    QTextStream inStream(&in);
    readOk = tab.readFile(inStream);
    QTextStream outStream(&outLong);
    outStream.setFieldAlignment(QTextStream::AlignLeft);
    tab.writeFile(outStream);
  }
  ok("read a line with 17-digit numbers", readOk);
  ok("sigma and the next field stay two tokens", outLong.contains("0.30000000000000004 "), outLong);
  ok("density and the next field stay two tokens", outLong.contains("0.12345678901234566 "), outLong);
  bool readOkLong = false;
  QString outLongB = roundTrip(tab, outLong + "</targetInt>\n", readOkLong);
  lines = tab.getTargetIntModel()->getLines();
  ok("re-read: sigma exact", readOkLong && lines.size() == 1 && lines.at(0).sigma == 0.30000000000000004, outLongB);

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

  // Free-gas Doppler broadening: `doppler <T>` read into the line, written back with the rest of the line
  // unchanged, a line without it gains nothing, and the dialog adds and edits it (Doppler-only lines included).
  {
    const QString dopOnly = "1 \"1,2\" 400 0 0 0 0 \"\" 0 0 0 0 \"\" 0 doppler 300\n";
    const QString dopGauss = "1 \"1,2\" 1000 1 0.0001 0 0 \"\" 0 0 0 0 \"\" 0 doppler 300\n";
    bool readP = false;
    QString outP = roundTrip(tab, dopOnly + dopGauss + "</targetInt>\n", readP);
    ok("read a Doppler-only line and a Doppler + Gaussian line", readP);
    QList<TargetIntData> pl = tab.getTargetIntModel()->getLines();
    ok("two Doppler rows", pl.size() == 2, outP);
    if (pl.size() == 2) {
      ok("Doppler-only: flagged, T = 300 K, no convolution, no target integration",
         pl.at(0).isDoppler && pl.at(0).dopplerTemperature == 300. && !pl.at(0).isConvolution && !pl.at(0).isTargetIntegration);
      ok("Doppler-only: 400 sub-points kept", pl.at(0).numPoints == 400);
      ok("Doppler + Gaussian: flagged, T = 300 K, Gaussian 0.0001",
         pl.at(1).isDoppler && pl.at(1).dopplerTemperature == 300. && pl.at(1).isConvolution && pl.at(1).sigma == 0.0001);
      TargetIntModel *m = tab.getTargetIntModel();
      ok("model columns 44/45 carry the flag and T",
         m->data(m->index(0, 44), Qt::EditRole).toBool() && m->data(m->index(0, 45), Qt::EditRole).toDouble() == 300.);
      ok("Doppler-only line shows its sub-point count (not N/A)",
         m->data(m->index(0, 2), Qt::DisplayRole).toString() == "400", m->data(m->index(0, 2), Qt::DisplayRole).toString());
      const QString sum0 = m->data(m->index(0, 33), Qt::DisplayRole).toString();
      ok("summary names the Doppler broadening", sum0 == "Free-gas Doppler broadening, T = 300 K", sum0);
      const QString sum1 = m->data(m->index(1, 33), Qt::DisplayRole).toString();
      ok("summary of the Gaussian line names both", sum1.contains("Gaussian convolution") && sum1.contains("Doppler"), sum1);
      ok("header names", m->headerData(44, Qt::Horizontal, Qt::DisplayRole).toString() == "Doppler?" &&
                         m->headerData(45, Qt::Horizontal, Qt::DisplayRole).toString() == "Doppler T [K]");
    }
    const QStringList outLines = outP.split('\n', Qt::SkipEmptyParts);
    ok("both lines written with ' doppler 300' at the end", outLines.size() == 2 && outLines.at(0).endsWith(" doppler 300") &&
       outLines.at(1).endsWith(" doppler 300"), outP);
    // The rest of each line is what the same line without the keyword writes.
    QString without = dopOnly;
    without.replace(" doppler 300", "");
    QString without2 = dopGauss;
    without2.replace(" doppler 300", "");
    bool readW = false;
    const QString outW = roundTrip(tab, without + without2 + "</targetInt>\n", readW);
    ok("lines without doppler gain no block", readW && !outW.contains("doppler"), outW);
    QString stripped = outP;
    stripped.replace(" doppler 300", "");
    ok("rest of each line unchanged", stripped == outW, stripped + " | " + outW);
    bool readP2 = false;
    const QString outP2 = roundTrip(tab, outP + "</targetInt>\n", readP2);
    ok("second round trip identical", readP2 && outP2 == outP, outP2);
    // Full precision, and beside the other keyword blocks.
    const QString mixed = "1 \"3\" 200 0 0 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 20 50 screening 934.7 1 200 doppler 293.61234567 "
                          "beamprofile 1 2.5 0.3 -1.5 1 0.07 0 1\n</targetInt>\n";
    const QString outM = roundTrip(tab, mixed, readP);
    pl = tab.getTargetIntModel()->getLines();
    ok("doppler read between screening and a beam profile", readP && pl.size() == 1 && pl.at(0).isDoppler &&
       pl.at(0).dopplerTemperature == 293.61234567 && pl.at(0).isScreening && pl.at(0).isBeamProfile &&
       pl.at(0).resonanceWidthMultiplier == 20. && pl.at(0).pointsPerWidth == 50., outM);
    ok("written before the beam profile, all digits kept",
       outM.contains(" screening 934.7 1 200 doppler 293.61234567 beamprofile 1 "), outM);
    const QString bad = roundTrip(tab, "1 \"3\" 200 0 0 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 20 50 doppler -5\n</targetInt>\n", readP);
    pl = tab.getTargetIntModel()->getLines();
    ok("a non-positive temperature is not taken as a Doppler block", pl.size() == 1 && !pl.at(0).isDoppler &&
       pl.at(0).pointsPerWidth == 50., bad);
    const QString noDop = roundTrip(tab, legacy, readP);
    ok("a legacy line gains no doppler block", !noDop.contains("doppler"), noDop);

    // The dialog on its own: hidden box, default T, sub-points enabled for a Doppler-only line, validation.
    {
      AddTargetIntDialog dlg;
      QGroupBox *box = dlg.findChild<QGroupBox*>("dopplerBox");
      ok("Doppler box hidden until checked", box && !box->isVisibleTo(&dlg));
      ok("sub-points disabled with nothing checked", !dlg.numPointsSpin->isEnabled());
      dlg.isDopplerCheck->setChecked(true);
      ok("Doppler box shown when checked", box && box->isVisibleTo(&dlg));
      ok("T defaults to 293.6 K", dlg.dopplerTemperatureText->text() == "293.6", dlg.dopplerTemperatureText->text());
      ok("Doppler alone enables the sub-points and the adaptive grid", dlg.numPointsSpin->isEnabled() &&
         dlg.resonanceWidthMultiplierSpin->isEnabled() && dlg.pointsPerWidthSpin->isEnabled());
      dlg.isConvolutionCheck->setChecked(true);
      dlg.isConvolutionCheck->setChecked(false);
      dlg.isTargetIntegrationCheck->setChecked(true);
      dlg.isTargetIntegrationCheck->setChecked(false);
      dlg.isBeamProfileCheck->setChecked(true);
      dlg.isBeamProfileCheck->setChecked(false);
      dlg.isUdrCheck->setChecked(true);
      dlg.isUdrCheck->setChecked(false);
      ok("other effects switched off again leave them enabled", dlg.numPointsSpin->isEnabled() &&
         dlg.resonanceWidthMultiplierSpin->isEnabled() && dlg.pointsPerWidthSpin->isEnabled());
      ok("293.6 K passes", dlg.dopplerProblem().isEmpty(), dlg.dopplerProblem());
      dlg.dopplerTemperatureText->setText("0");
      ok("0 K is refused", !dlg.dopplerProblem().isEmpty());
      dlg.dopplerTemperatureText->setText("-10");
      ok("a negative T is refused", !dlg.dopplerProblem().isEmpty());
      dlg.dopplerTemperatureText->setText("warm");
      ok("a non-number is refused", !dlg.dopplerProblem().isEmpty());
      dlg.isDopplerCheck->setChecked(false);
      ok("unchecked: no problem, sub-points disabled again", dlg.dopplerProblem().isEmpty() && !dlg.numPointsSpin->isEnabled() &&
         !dlg.resonanceWidthMultiplierSpin->isEnabled());
      if (const char* stem = std::getenv("AZURE2_TARGET_INT_PNG")) {
        dlg.dopplerTemperatureText->setText("300");
        dlg.isDopplerCheck->setChecked(true);
        dlg.show(); QApplication::processEvents(); dlg.adjustSize(); QApplication::processEvents();
        dlg.grab().save(QString("%1_doppler.png").arg(stem));
      }
    }

    // Through the tab: Add opens the dialog, which is filled and accepted as a user would.
    tab.reset();
    bool sawDialog = false;
    QTimer::singleShot(0, [&]() {
      AddTargetIntDialog *d = qobject_cast<AddTargetIntDialog*>(QApplication::activeModalWidget());
      if (!d) return;
      sawDialog = true;
      d->segmentsListText->setText("5");
      d->isDopplerCheck->setChecked(true);
      d->numPointsSpin->setValue(250);
      d->dopplerTemperatureText->setText("300");
      d->accept();
    });
    tab.addLine();
    ok("Add opened the dialog", sawDialog);
    pl = tab.getTargetIntModel()->getLines();
    ok("dialog-added Doppler-only line in the model", pl.size() == 1 && pl.at(0).isDoppler &&
       pl.at(0).dopplerTemperature == 300. && pl.at(0).numPoints == 250 && !pl.at(0).isConvolution && !pl.at(0).isTargetIntegration);
    QString added;
    QTextStream addedStream(&added);
    tab.writeFile(addedStream);
    ok("dialog-added line written with ' doppler 300'", added.trimmed().endsWith(" doppler 300") && added.contains(" 250"), added);

    // Edit: the dialog comes up with the line's T, keeps the sub-points enabled, and a new T is written.
    QTableView *view = tab.findChild<QTableView*>();
    bool editSaw = false, editT = false, editPoints = false;
    if (view) {
      view->selectionModel()->select(tab.getTargetIntModel()->index(0, 1),
                                     QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
      QTimer::singleShot(0, [&]() {
        AddTargetIntDialog *d = qobject_cast<AddTargetIntDialog*>(QApplication::activeModalWidget());
        if (!d) return;
        editSaw = true;
        editT = d->isDopplerCheck->isChecked() && d->dopplerTemperatureText->text() == "300";
        editPoints = d->numPointsSpin->isEnabled() && d->numPointsSpin->value() == 250;
        d->dopplerTemperatureText->setText("77.5");
        d->accept();
      });
      tab.editLine();
    }
    ok("Edit opened the dialog", editSaw);
    ok("edit dialog shows Doppler on at 300 K", editT);
    ok("edit dialog keeps the 250 sub-points enabled", editPoints);
    pl = tab.getTargetIntModel()->getLines();
    ok("edited T in the model", pl.size() == 1 && pl.at(0).isDoppler && pl.at(0).dopplerTemperature == 77.5);
    QString edited;
    QTextStream editedStream(&edited);
    tab.writeFile(editedStream);
    ok("edited line written with ' doppler 77.5'", edited.trimmed().endsWith(" doppler 77.5"), edited);

    // Switching it off in the edit dialog removes the block.
    QTimer::singleShot(0, [&]() {
      AddTargetIntDialog *d = qobject_cast<AddTargetIntDialog*>(QApplication::activeModalWidget());
      if (!d) return;
      d->isDopplerCheck->setChecked(false);
      d->accept();
    });
    if (view) tab.editLine();
    QString off;
    QTextStream offStream(&off);
    tab.writeFile(offStream);
    ok("unchecked in the edit dialog: no doppler block written", !off.contains("doppler"), off);
  }

  // The fitted temperature, `doppler T vary sigma [Tmax]` (parameter doppler_<n>): read into columns 46-49, written
  // back as read (the plain form stays plain, `0 sigma` stays the fit form), and set from the dialog.
  {
    bool readP = false;
    const QString fitA = "1 \"1,2\" 400 0 0 0 0 \"\" 0 0 0 0 \"\" 0 doppler 300 1 20 900\n";
    const QString fitB = "1 \"3,4\" 1000 1 0.0001 0 0 \"\" 0 0 0 0 \"\" 0 doppler 310 0 0\n";
    const QString plain = "1 \"5\" 400 0 0 0 0 \"\" 0 0 0 0 \"\" 0 doppler 320\n";
    const QString outP = roundTrip(tab, fitA + fitB + plain + "</targetInt>\n", readP);
    QList<TargetIntData> pl = tab.getTargetIntModel()->getLines();
    ok("read two fit-form lines and a plain one", readP && pl.size() == 3, outP);
    if (pl.size() == 3) {
      ok("fit form with Tmax: T 300, fitted, sigma 20, max 900",
         pl.at(0).isDoppler && pl.at(0).dopplerTemperature == 300. && pl.at(0).dopplerFit && pl.at(0).dopplerVary &&
         pl.at(0).dopplerSigma == 20. && pl.at(0).dopplerMax == 900.);
      ok("fit form 'vary 0 sigma 0': a fixed parameter, no Tmax, Gaussian kept",
         pl.at(1).isDoppler && pl.at(1).dopplerTemperature == 310. && pl.at(1).dopplerFit && !pl.at(1).dopplerVary &&
         pl.at(1).dopplerSigma == 0. && pl.at(1).dopplerMax == 0. && pl.at(1).isConvolution && pl.at(1).sigma == 0.0001);
      ok("plain line: no fit form", pl.at(2).isDoppler && pl.at(2).dopplerTemperature == 320. && !pl.at(2).dopplerFit);
      TargetIntModel *m = tab.getTargetIntModel();
      ok("columns 46-49 carry the fit form (edit role)",
         m->data(m->index(0, 46), Qt::EditRole).toBool() && m->data(m->index(0, 47), Qt::EditRole).toBool() &&
         m->data(m->index(0, 48), Qt::EditRole).toDouble() == 20. && m->data(m->index(0, 49), Qt::EditRole).toDouble() == 900.);
      ok("header names", m->headerData(46, Qt::Horizontal, Qt::DisplayRole).toString() == "Doppler T Parameter?" &&
                         m->headerData(49, Qt::Horizontal, Qt::DisplayRole).toString() == "Doppler Grid Maximum [K]");
      const QString sum0 = m->data(m->index(0, 33), Qt::DisplayRole).toString();
      ok("summary says the temperature is fitted", sum0.contains("T = 300 K (fitted)"), sum0);
      ok("defaultDopplerMax: max(T + 3 sigma, 1.21 T) / 2 T", defaultDopplerMax(300, 20) == 363. &&
         defaultDopplerMax(300, 0) == 600. && std::fabs(defaultDopplerMax(300, 5) - 363.) < 1e-9);
    }
    const QStringList outLines = outP.split('\n', Qt::SkipEmptyParts);
    ok("written back as read: ' doppler 300 1 20 900', ' doppler 310 0 0', ' doppler 320'",
       outLines.size() == 3 && outLines.at(0).endsWith(" doppler 300 1 20 900") &&
       outLines.at(1).endsWith(" doppler 310 0 0") && outLines.at(2).endsWith(" doppler 320"), outP);
    const QString mixed = roundTrip(tab, "1 \"3\" 200 0 0 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 20 50 doppler 293.6 1 15 screening 934.7 1 200\n"
                                         "</targetInt>\n", readP);
    pl = tab.getTargetIntModel()->getLines();
    ok("a fit form followed by another keyword: both read, both written",
       readP && pl.size() == 1 && pl.at(0).dopplerFit && pl.at(0).dopplerSigma == 15. && pl.at(0).isScreening &&
       mixed.contains(" screening 934.7 1 200") && mixed.contains(" doppler 293.6 1 15") && !mixed.contains("doppler 293.6 1 15 9"),
       mixed);

    // The dialog: the fields follow the 'Temperature Is a Fit Parameter' box; the default maximum is shown; refusals.
    {
      AddTargetIntDialog dlg;
      dlg.isDopplerCheck->setChecked(true);
      dlg.dopplerTemperatureText->setText("300");
      ok("fit fields disabled until the fit box is checked",
         !dlg.dopplerVaryCheck->isEnabled() && !dlg.dopplerSigmaText->isEnabled() && !dlg.dopplerMaxText->isEnabled());
      dlg.dopplerFitCheck->setChecked(true);
      ok("enabled with it", dlg.dopplerVaryCheck->isEnabled() && dlg.dopplerSigmaText->isEnabled() &&
         dlg.dopplerMaxText->isEnabled());
      dlg.dopplerSigmaText->setText("20");
      ok("blank maximum shows the default 363", dlg.dopplerMaxText->placeholderText() == "363",
         dlg.dopplerMaxText->placeholderText());
      dlg.dopplerSigmaText->setText("0");
      ok("without a prior the default is 2 T = 600", dlg.dopplerMaxText->placeholderText() == "600",
         dlg.dopplerMaxText->placeholderText());
      ok("a valid fit form passes", dlg.dopplerProblem().isEmpty(), dlg.dopplerProblem());
      dlg.dopplerSigmaText->setText("-1");
      ok("a negative prior sigma is refused", !dlg.dopplerProblem().isEmpty());
      dlg.dopplerSigmaText->setText("20");
      dlg.dopplerMaxText->setText("250");
      ok("a grid maximum not above T is refused", !dlg.dopplerProblem().isEmpty());
      dlg.dopplerMaxText->setText("900");
      ok("a grid maximum above T passes", dlg.dopplerProblem().isEmpty(), dlg.dopplerProblem());
      dlg.dopplerFitCheck->setChecked(false);
      dlg.dopplerMaxText->setText("250");
      ok("the fit fields are not checked without the fit box", dlg.dopplerProblem().isEmpty(), dlg.dopplerProblem());
      if (const char* stem = std::getenv("AZURE2_TARGET_INT_PNG")) {
        dlg.dopplerFitCheck->setChecked(true);
        dlg.dopplerMaxText->setText("");
        dlg.show(); QApplication::processEvents(); dlg.adjustSize(); QApplication::processEvents();
        dlg.grab().save(QString("%1_doppler_fit.png").arg(stem));
      }
    }

    // Through the tab: edit the plain line into the fit form, then the dialog comes back with it.
    roundTrip(tab, plain + "</targetInt>\n", readP);
    QTableView *view = tab.findChild<QTableView*>();
    bool sawFit = false, backFit = false;
    if (view) {
      view->selectionModel()->select(tab.getTargetIntModel()->index(0, 1),
                                     QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
      QTimer::singleShot(0, [&]() {
        AddTargetIntDialog *d = qobject_cast<AddTargetIntDialog*>(QApplication::activeModalWidget());
        if (!d) return;
        sawFit = !d->dopplerFitCheck->isChecked();
        d->dopplerFitCheck->setChecked(true);
        d->dopplerVaryCheck->setChecked(true);
        d->dopplerSigmaText->setText("25");
        d->dopplerMaxText->setText("1000");
        d->accept();
      });
      tab.editLine();
      QTimer::singleShot(0, [&]() {
        AddTargetIntDialog *d = qobject_cast<AddTargetIntDialog*>(QApplication::activeModalWidget());
        if (!d) return;
        backFit = d->dopplerFitCheck->isChecked() && d->dopplerVaryCheck->isChecked() &&
                  d->dopplerSigmaText->text() == "25" && d->dopplerMaxText->text() == "1000";
        d->reject();
      });
      tab.editLine();
    }
    ok("the edit dialog opened on the plain line without the fit form", sawFit);
    ok("and comes back with the fit form it was given", backFit);
    QString fitOut;
    QTextStream fitStream(&fitOut);
    tab.writeFile(fitStream);
    ok("edited line written with ' doppler 320 1 25 1000'", fitOut.trimmed().endsWith(" doppler 320 1 25 1000"), fitOut);
  }

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
