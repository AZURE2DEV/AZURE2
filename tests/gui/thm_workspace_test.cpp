// Headless test of the THM workspace (Configure > THM Workspace...).
//
// Checked here, through AZURESetup on projects built from tests/7Li_p_a:
//  1. experiment records: read/compose keep unchanged lines verbatim and keys
//     the page does not show; the nuclide list is the engine's table;
//  2. a project without a THM segment: pages disabled, explanation shown, and
//     the project saves byte for byte;
//  3. a THM project opened and accepted untouched saves byte for byte, and the
//     classic editor no longer shows the binding energy or the RWA flag;
//  4. Model page: an option set there is written into the <thm> block;
//  5. Experiments page: add (segment picker offers only THM segments with a
//     free norm in no other experiment), background, kinematics with the
//     derived energies, save, reopen, edit in place, refusals with the
//     engine's messages, remove;
//  6. Channels page: B and the RWA flag edited there are written on the
//     <levels> lines, exactly, and come back on reopening;
//  7. the engine runs the GUI's experiment block and prints the binding and
//     quasi-free energies the page shows;
//  8. the Coulomb line shape switch: offered only with a complete reaction,
//     written as lineshape=on, read back, refused without Brune; zeta at the
//     ends of the data as the engine prints it (charged spectator);
//  9. the Diagnostics page (builds with the engine API): the snapshot it runs
//     is what Accept + save would write and leaves the project untouched; the
//     engine's curves are finite, the vertex nodes are where the engine's
//     ThmFormFactor changes sign, the line shape and weight panels appear
//     with their options, Compute runs off the GUI thread;
// 10. the spectator-momentum window (ps=, psNodes=): controls offered with a
//     complete reaction, only the chosen distribution's fields shown, write /
//     read back / verbatim for every distribution, a table stored relative to
//     the project, refusals (no reaction, bad table, with spectatorEnergy --
//     the engine refuses with the same words), <T_s> = the engine's;
// 11. Diagnostics with a window: nodes, weights, mu_sx and <T_s> are the
//     engine's, w(p) the Hulthen weight with the nodes on it, <|M_l|^2> =
//     sum_k w_k M_l^2 with the engine's ThmFormFactor, the nodes filled; a
//     relative ps table read from the engine's temporary copy.
//
// 12. the distortion factor (Stage D) on a 12C+12C-like project;
// 13. the exit angle (theta=, Stage E): read/written verbatim, the window
//     row only for a window, 0-0 a single angle, refusals in the engine's
//     words (a reversed window; entranceL=coherent on the Model page, live,
//     and the engine refuses the same file); Diagnostics: the angular
//     distribution card, 4 pi <dsigma/dOmega> over 0-180 = the engine's
//     angle-integrated HOES, the window average = the curve's average.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDir>
#include <QFile>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QPixmap>
#include <QProcess>
#include <QPushButton>
#include <QRadioButton>
#include <QRegExp>
#include <QScrollArea>
#include <QStandardItemModel>
#include <QTabWidget>
#include <QTableView>
#include <QTableWidget>
#include <QTemporaryDir>
#include <iostream>
#include <vector>

#include "AZURESetup.h"
#include "AddPairDialog.h"
#include "ChannelDetails.h"
#include "ChannelsModel.h"
#include "Config.h"
#include "ThmChannelsPage.h"
#include "ThmExperiment.h"
#include "ThmExperimentsPage.h"
#include "ThmModelPage.h"
#include "ThmNumberSpin.h"
#include "ThmWorkspace.h"
#include "Constants.h"
#include "ThmFunc.h"
#include "ThmLineshape.h"
#include <algorithm>
#ifdef AZURE2_THM_DIAGNOSTICS
#include <QEventLoop>
#include <QProgressBar>
#include <QTimer>
#include "ThmDiagnosticsPage.h"
#include "ThmPlotWidget.h"
#endif
#include <cmath>
struct SegPairs {int firstPair; int secondPair;};

// Defined by AZURE2.cpp, which belongs to the executable rather than the GUI
// library, so this test supplies its own. They are never called from here.
Config* g_config = nullptr;
void exitMessage(const Config&) {}
bool checkExternalCapture(Config&, const std::vector<SegPairs>&) { return true; }
bool readSegmentFile(const Config&, std::vector<SegPairs>&) { return true; }
void startMessage(const Config&) {}

#ifndef AZURE2_SOURCE_DIR
#error "AZURE2_SOURCE_DIR must be defined"
#endif
#ifndef AZURE2_BINARY
#error "AZURE2_BINARY must be defined"
#endif

static int fails = 0;
static void ok(const char* what, bool cond, const QString& detail = QString()) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what;
  if(!cond && !detail.isEmpty()) std::cout << "  -- " << detail.toStdString();
  std::cout << std::endl;
  if(!cond) fails++;
}

static QString slurp(const QString& path) {
  QFile f(path);
  if(!f.open(QIODevice::ReadOnly)) return QString();
  return QString::fromUtf8(f.readAll());
}
static void spit(const QString& path, const QString& text) {
  QFile f(path);
  if(f.open(QIODevice::WriteOnly | QIODevice::Truncate)) f.write(text.toUtf8());
}
static QString blockOf(const QString& text) {
  int a = text.indexOf("<thm>\n");
  if(a < 0) return "<none>";
  a += 6;
  int b = text.indexOf("</thm>\n", a);
  return b < 0 ? "<unterminated>" : text.mid(a, b - a);
}
// The <levels> lines of a saved project, split into fields.
static QList<QStringList> levelLines(const QString& text) {
  QList<QStringList> out;
  int a = text.indexOf("<levels>\n"), b = text.indexOf("</levels>");
  for(const QString& line : text.mid(a + 9, b - a - 9).split('\n'))
    if(!line.trimmed().isEmpty()) out << line.split(QRegExp("\\s+"), Qt::SkipEmptyParts);
  return out;
}
static QStringList checkedSegments(QListWidget* list) {
  QStringList s;
  for(int i = 0; i < list->count(); i++)
    if(list->item(i)->checkState() == Qt::Checked) s << list->item(i)->data(Qt::UserRole).toString();
  return s;
}
static QStringList offeredSegments(QListWidget* list) {
  QStringList s;
  for(int i = 0; i < list->count(); i++) s << list->item(i)->data(Qt::UserRole).toString();
  return s;
}
// A number typed into a field of the THM workspace (its text as typed).
static void typeNumber(ThmNumberSpin* e, const QString& text) {
  const double before = e->value();
  e->setValue(text.toDouble());
  if(e->value() == before) emit e->valueChanged(e->value());
}

static QListWidgetItem* segmentItem(QListWidget* list, int key) {
  for(int i = 0; i < list->count(); i++)
    if(list->item(i)->data(Qt::UserRole).toInt() == key) return list->item(i);
  return nullptr;
}

// Runs "Calculate Segments From Data" on dir/name; the console output.
static QString engineRun(const QString& dir, const QString& name, int* exitCode) {
  QDir(dir).mkpath("output");
  QDir(dir).mkpath("checks");
  QProcess p;
  p.setWorkingDirectory(dir);
  p.setProcessChannelMode(QProcess::MergedChannels);
  QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
  if(!env.contains("OMP_NUM_THREADS")) env.insert("OMP_NUM_THREADS", "1");
  p.setProcessEnvironment(env);
  p.start(AZURE2_BINARY, QStringList() << "--no-gui" << "--no-readline" << name);
  *exitCode = -1;
  if(!p.waitForStarted(30000)) return QString();
  p.write("1\n\n\n7\n");
  p.closeWriteChannel();
  if(!p.waitForFinished(600000)) { p.kill(); return QString(); }
  *exitCode = p.exitCode();
  return QString::fromUtf8(p.readAll());
}

int main(int argc, char** argv) {
  QTemporaryDir settingsDir;
  qputenv("XDG_CONFIG_HOME", settingsDir.path().toUtf8());
  QApplication app(argc, argv);
  QCoreApplication::setOrganizationName("AZURE2-tests");
  QCoreApplication::setApplicationName("thm_workspace_test");

  // 1. Records and the nuclide table.
  {
    const QStringList lines = QStringList() << "experiment[A] segments=1-3   # three runs"
                                            << "experiment[B] segments=4 future=3 background=const"
                                            << "experiment[A] beam=7Li target=d spectator=n Ebeam=19";
    QList<ThmExperimentRecord> r = ThmExperimentRecord::read(lines);
    ok("records merged by name", r.size() == 2 && r[0].name == "A" && r[1].name == "B");
    ok("segment ranges expanded", r.size() == 2 && r[0].segments == (QList<int>() << 1 << 2 << 3));
    ok("unshown key kept aside", r.size() == 2 && r[1].extraTokens == QStringList("future=3"));
    QList<ThmExperimentRecord> ls = ThmExperimentRecord::read(QStringList() << "experiment[C] segments=1 lineshape=on"
                                                                            << "experiment[D] segments=2 lineshape=off");
    ok("lineshape read as a field", ls.size() == 2 && ls[0].lineshape && !ls[1].lineshape && ls[0].extraTokens.isEmpty() &&
                                        ls[1].extraTokens.isEmpty());
    QList<ThmExperimentRecord> lsOff = ls;
    lsOff[0].lineshape = false;
    ok("lineshape off is not written (the default)",
       ThmExperimentRecord::compose(QStringList() << "experiment[C] segments=1 lineshape=on", ls.mid(0, 1),
                                    lsOff.mid(0, 1)) == QStringList("experiment[C] segments=1"));
    QList<ThmExperimentRecord> ps = ThmExperimentRecord::read(QStringList() << "experiment[P] psNodes=24 segments=1 ps=gauss:50:0-40.0");
    ok("ps and psNodes read as fields, as written", ps.size() == 1 && ps[0].ps == "gauss:50:0-40.0" &&
                                                        ps[0].psNodes == "24" && ps[0].extraTokens.isEmpty() && ps[0].hasWindow());
    ok("ps=delta is no window", !ThmExperimentRecord::read(QStringList() << "experiment[P] segments=1 ps=delta")[0].hasWindow());
    ok("unchanged records: lines verbatim", ThmExperimentRecord::compose(lines, r, r) == lines);
    QList<ThmExperimentRecord> e = r;
    e[1].background = "linear";
    ok("edited record keeps the key it does not show",
       ThmExperimentRecord::compose(lines, r, e) ==
           (QStringList() << lines[0] << "experiment[B] segments=4 background=linear future=3" << lines[2]),
       ThmExperimentRecord::compose(lines, r, e).join("|"));
    e = r;
    e[0].background = "quadratic";  // two lines of A become one, in place of the first
    ok("edited record: one line at its first place",
       ThmExperimentRecord::compose(lines, r, e) ==
           (QStringList() << "experiment[A] segments=1-3 background=quadratic beam=7Li target=d spectator=n Ebeam=19"
                          << lines[1]),
       ThmExperimentRecord::compose(lines, r, e).join("|"));
    // Stage D: every distortion key, the ten-number lists included, comes
    // back verbatim when the page rewrites the line (another key edited).
    {
      const QStringList stageD = {"distortion=optical", "opticalAA=50.0,4.5,0.60,10,4.5,0.6,0,0,0,4.50",
                                  "opticalSF=12,5.0,0.7,0.0,0,0,8,5.5,0.65,0", "spectatorAngle=cm:8.0",
                                  "distortionRef=2.664", "distortionRatio=dw", "boundState=yukawa:3.0"};
      const QStringList dLines = {"experiment[S] boundState=yukawa:3.0 segments=1 beam=14N target=12C spectator=d "
                                  "Ebeam=30 " + stageD.mid(0, 6).join(' ')};
      QList<ThmExperimentRecord> d = ThmExperimentRecord::read(dLines);
      ok("stage D keys read as fields, as written",
         d.size() == 1 && d[0].distortion == "optical" && d[0].opticalAA == "50.0,4.5,0.60,10,4.5,0.6,0,0,0,4.50" &&
             d[0].opticalSF == "12,5.0,0.7,0.0,0,0,8,5.5,0.65,0" && d[0].spectatorAngle == "cm:8.0" &&
             d[0].distortionRef == "2.664" && d[0].distortionRatio == "dw" && d[0].boundState == "yukawa:3.0" &&
             d[0].extraTokens.isEmpty() && d[0].hasComputedDistortion());
      ok("stage D: unchanged line verbatim", ThmExperimentRecord::compose(dLines, d, d) == dLines);
      QList<ThmExperimentRecord> e2 = d;
      e2[0].background = "const";
      const QStringList rewritten = ThmExperimentRecord::compose(dLines, d, e2);
      bool all = rewritten.size() == 1;
      for(const QString& t : stageD) all = all && rewritten.value(0).split(' ').contains(t);
      ok("stage D: every key verbatim through a rewrite", all, rewritten.join("|"));
      ok("stage D: the rewritten line passes the engine's parser and check",
         ThmSettings::checkExperimentLines(rewritten).isEmpty(), ThmSettings::checkExperimentLines(rewritten));
      ok("table distortion needs no kinematics",
         ThmExperimentRecord::read(QStringList() << "experiment[T] segments=1 distortion=table:w.dat")[0].hasDistortion() &&
             !ThmExperimentRecord::read(QStringList() << "experiment[T] segments=1 distortion=table:w.dat")[0]
                  .hasComputedDistortion() &&
             ThmSettings::checkExperimentLines(QStringList() << "experiment[T] segments=1 distortion=table:w.dat").isEmpty());
    }
    // Stage E: theta= read as a field, verbatim, and through a rewrite.
    {
      const QStringList tLines = {"experiment[T] segments=1 theta=50.0-70 future=1", "experiment[U] segments=2 theta=all"};
      QList<ThmExperimentRecord> t = ThmExperimentRecord::read(tLines);
      ok("theta read as a field, as written", t.size() == 2 && t[0].theta == "50.0-70" && t[0].hasTheta() &&
                                                  t[0].extraTokens == QStringList("future=1") && t[1].theta == "all" &&
                                                  !t[1].hasTheta() && t[1].extraTokens.isEmpty());
      ok("theta: unchanged lines verbatim", ThmExperimentRecord::compose(tLines, t, t) == tLines);
      QList<ThmExperimentRecord> e2 = t;
      e2[0].background = "const";
      e2[1].background = "const";
      const QStringList rewritten = ThmExperimentRecord::compose(tLines, t, e2);
      ok("theta: verbatim through a rewrite",
         rewritten == QStringList({"experiment[T] segments=1 background=const theta=50.0-70 future=1",
                                   "experiment[U] segments=2 background=const theta=all"}),
         rewritten.join("|"));
    }
    ok("segment list text", ThmExperimentRecord::segmentsListText(QList<int>() << 1 << 2 << 3 << 5 << 7 << 8) ==
                                "1-3,5,7,8");

    const QStringList names = ThmExperimentsPage::nuclideNames();
    ThmNuclide n;
    const QString why = QString::fromStdString(ThmNuclide::Parse("nothing", n));
    ok("nuclide list is the engine's table", why.contains("(known: " + names.join(' ') + ";"), why);
    bool all = !names.isEmpty();
    for(const QString& name : names) all = all && ThmNuclide::Parse(name.toStdString(), n).empty();
    ok("every listed nuclide parses", all && names.contains("7Li") && names.contains("d"));
  }

  // Projects.
  QTemporaryDir work;
  const QString src = QString(AZURE2_SOURCE_DIR) + "/tests/7Li_p_a";
  QString plain = slurp(src + "/7Li_p_a.azr");
  ok("found tests/7Li_p_a", !plain.isEmpty());
  if(!plain.endsWith('\n')) plain += '\n';
  QDir(work.path()).mkpath("data");
  QFile::copy(src + "/data/tumino_thm.dat", work.filePath("data/tumino_thm.dat"));
  const QString seg1 = "1  5  4  0  8.2  0  180  10  0.000597139  1  0  0  0  0  data/tumino_thm.dat  0  0";
  ok("7Li_p_a has its THM segment line", plain.contains(seg1));
  spit(work.filePath("plain.azr"), plain);
  QString noThm = plain;
  noThm.replace(seg1, "1  5  4  0  8.2  0  180  0  0.000597139  1  0  0  0  0  data/tumino_thm.dat  0  0");
  spit(work.filePath("nothm.azr"), noThm);
  // Four data segments: 1, 2 THM with a free norm, 3 THM with a fixed norm, 4 not THM.
  QString four = plain;
  four.replace(seg1, seg1 + "\n" + "1  5  4  0  4.1  0  180  10  0.000597139  1  0  0  0  0  data/tumino_thm.dat  0  0\n" +
                         "1  5  4  0  8.2  0  180  10  0.000597139  0  0  0  0  0  data/tumino_thm.dat  0  0\n" +
                         "1  5  4  0  3  0  180  0  0.000597139  1  0  0  0  0  data/tumino_thm.dat  0  0");
  spit(work.filePath("four.azr"), four);

  AZURESetup w;

  // 2. No THM segment.
  {
    const QString path = work.filePath("nothm.azr");
    w.open(path);
    w.saveProject();
    const QString before = slurp(path);
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    ok("no THM: pages disabled", !ws.hasThmSegments() && !ws.pages->isEnabled() && !ws.acceptButton->isEnabled());
    ok("no THM: explanation shown", ws.noThmLabel->isVisibleTo(&ws) && ws.noThmLabel->text().contains("no THM"));
#ifdef AZURE2_THM_DIAGNOSTICS
    ok("no THM: four pages, the diagnostics page disabled with the others",
       ws.pages->count() == 4 && ws.diagnosticsPage && !ws.diagnosticsPage->isEnabled());
#else
    ok("no THM: three pages (no engine API: no diagnostics page)", ws.pages->count() == 3 && !ws.diagnosticsPage);
#endif
    ws.accept();
    w.saveProject();
    ok("no THM: byte-identical save", slurp(path) == before);
  }

  // 3. Untouched THM project; the classic editor.
  {
    const QString block = "# options\nvertex=onshell  # note\nexperiment[E1] segments=1   # the THM run\n";
    const QString path = work.filePath("untouched.azr");
    spit(path, plain + "<thm>\n" + block + "</thm>\n");
    w.open(path);
    w.saveProject();
    const QString before = slurp(path);
    ThmSettings s;
    QString err;
    ok("block with an experiment opens", w.thmSettings(s, &err), err);
    ThmWorkspace ws(&w, s);
    ok("THM project: pages enabled", ws.hasThmSegments() && ws.pages->isEnabled() && !ws.noThmLabel->isVisibleTo(&ws));
    ok("untouched: nothing refused", ws.validate().isEmpty(), ws.validate());
    ws.accept();
    ok("untouched: accepted", ws.result() == QDialog::Accepted);
    w.saveProject();
    ok("untouched: byte-identical save", slurp(path) == before);
    ok("untouched: block verbatim", blockOf(slurp(path)) == block, blockOf(slurp(path)));

    // The classic editor: no binding energy, no RWA radio buttons.
    AddPairDialog pairDialog;
    bool bindingShown = false;
    for(QLabel* l : pairDialog.findChildren<QLabel*>()) bindingShown = bindingShown || l->text().contains("Binding");
    ok("classic pair dialog has no binding energy", !bindingShown);
    ChannelDetails details;
    ok("classic channel details have no RWA radio buttons", details.findChildren<QRadioButton*>().isEmpty());
    bool hidden = false;
    for(QTableView* v : w.getPairsTab()->findChildren<QTableView*>()) hidden = hidden || v->isColumnHidden(15);
    ok("classic pairs table hides the binding energy column", hidden);
  }

  // 4. Model page.
  {
    const QString path = work.filePath("model.azr");
    spit(path, plain);
    w.open(path);
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    ws.modelPage->kinematicsCombo->setCurrentText("triple");
    ws.modelPage->coulombIntegralCheck->setChecked(true);
    ok("model page: accepted", ws.validate().isEmpty());
    ws.accept();
    w.saveProject();
    ok("model page: block written", blockOf(slurp(path)) == "kinematics=triple\ncoulombIntegral=1\n",
       blockOf(slurp(path)));
  }

  // 5. Experiments page.
  const QString fourPath = work.filePath("exp.azr");
  spit(fourPath, four);
  {
    w.open(fourPath);
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    ThmExperimentsPage* p = ws.experimentsPage;
    ok("no experiment yet", p->experimentTable->rowCount() == 0 && !p->editorBox->isEnabled());
    p->addExperiment();
    ok("added E1", p->experimentTable->rowCount() == 1 && p->nameEdit->text() == "E1" && p->editorBox->isEnabled());
    ok("picker offers THM segments with a free norm only", offeredSegments(p->segmentList) == QStringList({"1", "2"}),
       offeredSegments(p->segmentList).join(","));
    segmentItem(p->segmentList, 1)->setCheckState(Qt::Checked);
    segmentItem(p->segmentList, 2)->setCheckState(Qt::Checked);
    p->backgroundCombo->setCurrentText("linear");
    p->kinematicsBox->setChecked(true);
    p->beamCombo->setEditText("7Li");
    p->targetCombo->setEditText("d");
    p->spectatorCombo->setEditText("n");
    typeNumber(p->beamEnergyEdit, "19");
    ok("derived energies shown", p->derivedText().contains("B(x+s) = 2.22") &&
                                     p->derivedText().contains("quasi-free E(x+A)"),
       p->derivedText());
    ok("table row", p->experimentTable->item(0, 1)->text() == "1,2" && p->experimentTable->item(0, 2)->text() == "linear");

    // A second experiment: nothing left to offer; refusals with the engine's messages.
    p->addExperiment();
    ok("E2: segments of E1 not offered", offeredSegments(p->segmentList).isEmpty(),
       offeredSegments(p->segmentList).join(","));
    p->setSegmentsOfCurrent(QList<int>() << 1);
    ok("refused: segment in two experiments", ws.validate().contains("segment 1 is already in experiment[E1]"),
       ws.validate());
    ok("refusal shows the page", ws.pages->currentWidget() == p);
    p->setSegmentsOfCurrent(QList<int>() << 4);
    ok("refused: not a THM segment", ws.validate() == "<thm> experiment[E2]: segment 4 is not a THM segment (isDiff < 10).",
       ws.validate());
    p->setSegmentsOfCurrent(QList<int>() << 3);
    ok("refused: fixed norm", ws.validate().contains("segment 3 has a fixed norm"), ws.validate());
    p->setSegmentsOfCurrent(QList<int>() << 9);
    ok("refused: no such segment", ws.validate().contains("segment 9: <segmentsData> has only 4 line(s)"),
       ws.validate());
    p->nameEdit->setText("E 2");
    emit p->nameEdit->textEdited("E 2");
    ok("refused: bad name (engine's message)", ws.validate().contains("a name is letters, digits"), ws.validate());
    p->removeCurrent();
    ok("E2 removed", p->experimentTable->rowCount() == 1 && p->currentRow() == 0 && p->nameEdit->text() == "E1");

    // Partial kinematics, and a reaction that does not give the entrance pair.
    p->spectatorCombo->setEditText("");
    ok("refused: partial kinematics", ws.validate().contains("all four or none"), ws.validate());
    p->spectatorCombo->setEditText("p");
    ok("refused: reaction not the entrance pair", ws.validate().contains("does not give the entrance pair"),
       ws.validate());
    p->spectatorCombo->setEditText("n");
    ok("accepted", ws.validate().isEmpty(), ws.validate());
    ws.accept();
    w.saveProject();
    ok("experiment written", blockOf(slurp(fourPath)) ==
                                 "experiment[E1] segments=1,2 background=linear beam=7Li target=d spectator=n Ebeam=19\n",
       blockOf(slurp(fourPath)));
  }
  {
    // Reopen: the page shows what was saved.
    w.open(work.filePath("plain.azr"));
    w.open(fourPath);
    ThmSettings s;
    QString err;
    ok("saved experiment block parses", w.thmSettings(s, &err), err);
    ThmWorkspace ws(&w, s);
    ThmExperimentsPage* p = ws.experimentsPage;
    ok("reopened: one experiment, selected", p->experimentTable->rowCount() == 1 && p->currentRow() == 0);
    ok("reopened: editor", p->nameEdit->text() == "E1" && checkedSegments(p->segmentList) == QStringList({"1", "2"}) &&
                               p->backgroundCombo->currentText() == "linear" && p->kinematicsBox->isChecked() &&
                               p->beamCombo->currentText() == "7Li" && p->targetCombo->currentText() == "d" &&
                               p->spectatorCombo->currentText() == "n" && p->beamEnergyEdit->writtenText() == "19");
    // Edit: segment 2 out, background const; the kinematics unticked as a unit.
    segmentItem(p->segmentList, 2)->setCheckState(Qt::Unchecked);
    p->backgroundCombo->setCurrentText("const");
    p->kinematicsBox->setChecked(false);
    ok("edited: accepted", ws.validate().isEmpty(), ws.validate());
    ws.accept();
    w.saveProject();
    ok("edited: rewritten", blockOf(slurp(fourPath)) == "experiment[E1] segments=1 background=const\n",
       blockOf(slurp(fourPath)));
  }
  {
    // Two experiments with comments: editing one leaves the other verbatim.
    const QString block = "# runs\nexperiment[A] segments=1   # first run\nvertex=onshell\n"
                          "experiment[B] segments=2 background=const  # second\n";
    spit(fourPath, four + "<thm>\n" + block + "</thm>\n");
    w.open(fourPath);
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    ThmExperimentsPage* p = ws.experimentsPage;
    p->selectExperiment(1);
    p->backgroundCombo->setCurrentText("quadratic");
    ws.accept();
    w.saveProject();
    ok("edit in place, other experiment verbatim",
       blockOf(slurp(fourPath)) == "# runs\nexperiment[A] segments=1   # first run\nvertex=onshell\n"
                                   "experiment[B] segments=2 background=quadratic\n",
       blockOf(slurp(fourPath)));
    // Remove both: the options stay; comments stay with them.
    ThmSettings s2;
    w.thmSettings(s2);
    ThmWorkspace ws2(&w, s2);
    ws2.experimentsPage->selectExperiment(0);
    ws2.experimentsPage->removeCurrent();
    ws2.experimentsPage->removeCurrent();
    ws2.accept();
    w.saveProject();
    ok("experiments removed", blockOf(slurp(fourPath)) == "# runs\nvertex=onshell\n", blockOf(slurp(fourPath)));
  }
  {
    // lineshape=on as written is shown by the switch and kept through an edit.
    const QString block = "experiment[L] segments=1 lineshape=on beam=7Li target=d spectator=n Ebeam=60\n";
    spit(fourPath, four + "<thm>\n" + block + "</thm>\n");
    w.open(fourPath);
    ThmSettings s;
    QString err;
    ok("lineshape: block opens", w.thmSettings(s, &err), err);
    ThmWorkspace ws(&w, s);
    ok("lineshape: switch on, enabled", ws.experimentsPage->lineshapeCheck->isChecked() &&
                                            ws.experimentsPage->lineshapeCheck->isEnabled());
    ok("lineshape: zeta shown (neutral spectator: 0)",
       ws.experimentsPage->derivedText().contains(QString::fromUtf8("\u03b6 = 0 \u2026 0 into")),
       ws.experimentsPage->derivedText());
    ws.experimentsPage->backgroundCombo->setCurrentText("const");
    ws.accept();
    w.saveProject();
    ok("lineshape kept through an edit",
       blockOf(slurp(fourPath)) ==
           "experiment[L] segments=1 background=const beam=7Li target=d spectator=n Ebeam=60 lineshape=on\n",
       blockOf(slurp(fourPath)));
  }
  {
    // The switch: offered only with a complete reaction, written, read back,
    // dropped with the reaction.
    spit(fourPath, four);
    w.open(fourPath);
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    ThmExperimentsPage* p = ws.experimentsPage;
    p->addExperiment();
    segmentItem(p->segmentList, 1)->setCheckState(Qt::Checked);
    ok("lineshape: disabled without the reaction", !p->lineshapeCheck->isEnabled());
    p->kinematicsBox->setChecked(true);
    p->beamCombo->setEditText("7Li");
    p->targetCombo->setEditText("d");
    p->spectatorCombo->setEditText("n");
    ok("lineshape: disabled with an incomplete reaction", !p->lineshapeCheck->isEnabled());
    typeNumber(p->beamEnergyEdit, "60");
    ok("lineshape: enabled with the complete reaction", p->lineshapeCheck->isEnabled() && !p->lineshapeCheck->isChecked());
    p->lineshapeCheck->setChecked(true);
    ok("lineshape: in the record", p->records().at(0).lineshape);
    ok("lineshape: accepted", ws.validate().isEmpty(), ws.validate());
    ws.accept();
    w.saveProject();
    ok("lineshape: written", blockOf(slurp(fourPath)) ==
                                 "experiment[E1] segments=1 beam=7Li target=d spectator=n Ebeam=60 lineshape=on\n",
       blockOf(slurp(fourPath)));
    ThmSettings s2;
    w.thmSettings(s2);
    ThmWorkspace again(&w, s2);
    ok("lineshape: read back", again.experimentsPage->lineshapeCheck->isChecked());
    again.experimentsPage->kinematicsBox->setChecked(false);
    ok("lineshape: dropped with the reaction", !again.experimentsPage->records().at(0).lineshape &&
                                                   !again.experimentsPage->lineshapeCheck->isChecked());
    again.accept();
    w.saveProject();
    ok("lineshape: removed", blockOf(slurp(fourPath)) == "experiment[E1] segments=1\n", blockOf(slurp(fourPath)));

    // Without Brune the engine refuses the line shape; the page says so.
    ThmExperimentsPage noBrune(QStringList() << "experiment[E1] segments=1 beam=7Li target=d spectator=n Ebeam=60 lineshape=on",
                               w.getSegmentsTab()->getSegmentsDataModel(), w.getPairsTab()->getPairsModel(),
                               w.projectDirectory(), false);
    ok("lineshape: refused without Brune", noBrune.check().contains("needs the Brune parameterization"),
       noBrune.check());
    // E_sF <= 0 at the highest point: refused with the engine's words.
    ThmExperimentsPage low(QStringList() << "experiment[E1] segments=1 beam=3He target=7Li spectator=d Ebeam=8 lineshape=on",
                           w.getSegmentsTab()->getSegmentsDataModel(), w.getPairsTab()->getPairsModel(),
                           w.projectDirectory(), true);
    ok("lineshape: refused when the spectator has no energy left", low.check().contains("has no energy left"),
       low.check());
  }

  // 6. Channels page.
  {
    const QString path = work.filePath("channels.azr");
    spit(path, plain);
    w.open(path);
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    ThmChannelsPage* c = ws.channelsPage;
    ok("channels: the THM entrance pair only", c->pairTable->rowCount() == 1 && c->pairRow(5) == 0);
    QLineEdit* b = qobject_cast<QLineEdit*>(c->pairTable->cellWidget(0, 3));
    ok("channels: B as in the file", b && b->text() == "2.2246", b ? b->text() : QString());
    const QList<ChannelsData> channels = w.getLevelsTab()->getChannelsModel()->getChannels();
    int pair5 = 0, first = -1;
    for(int i = 0; i < channels.size(); i++)
      if(channels[i].pairIndex == 4 && channels[i].radType == 'P') {
        pair5++;
        if(first < 0 && channels[i].reducedWidth != 0.0) first = i;
      }
    ok("channels: every particle channel of the pair listed", c->channelTable->rowCount() == pair5 && pair5 > 0);
    ok("channels: flags as read", c->channelRow(first) >= 0 &&
                                      c->channelTable->item(c->channelRow(first), 5)->checkState() == Qt::Unchecked);
    c->setBindingText(5, "abc");
    ok("refused: B not a number", ws.validate().contains("not a number"), ws.validate());
    c->setBindingText(5, "2.224566");
    c->setReducedWidthFlag(first, true);
    ok("channels: accepted", ws.validate().isEmpty(), ws.validate());
    ws.accept();
    w.saveProject();
    const QList<QStringList> lines = levelLines(slurp(path));
    int withB = 0, flagged = 0, bad = 0;
    const ChannelsData& f = channels[first];
    for(const QStringList& l : lines) {
      if(l.size() < 33) { bad++; continue; }
      if(l[5] == "5") {
        if(l[31] == "2.224566") withB++;
        if(l[32] == "1") {
          flagged++;
          if(l[6].toInt() != int(f.sValue * 2) || l[7].toInt() != f.lValue * 2 ||
             l[11].toDouble() != f.reducedWidth)
            bad++;
        }
      } else if(l[32] != "0") {
        bad++;
      }
    }
    ok("channels: B written exactly on every line of the pair", withB == pair5 && bad == 0,
       QString("%1 of %2, bad %3").arg(withB).arg(pair5).arg(bad));
    ok("channels: one line flagged, the edited channel", flagged == 1);

    w.open(work.filePath("plain.azr"));
    w.open(path);
    ThmSettings s2;
    w.thmSettings(s2);
    ThmWorkspace again(&w, s2);
    QLineEdit* b2 = qobject_cast<QLineEdit*>(again.channelsPage->pairTable->cellWidget(0, 3));
    ok("reopened: B", b2 && b2->text() == "2.224566");
    const QList<ChannelsData> ch2 = w.getLevelsTab()->getChannelsModel()->getChannels();
    int flaggedRows = 0;
    for(int r = 0; r < again.channelsPage->channelTable->rowCount(); r++)
      if(again.channelsPage->channelTable->item(r, 5)->checkState() == Qt::Checked) flaggedRows++;
    int flaggedModel = 0;
    for(const ChannelsData& ch : ch2) flaggedModel += ch.gammaIsRWA == 1;
    ok("reopened: the flag", flaggedRows == 1 && flaggedModel == 1);
    // Untouched page: byte for byte.  (Saved once first: the Fitting tab
    // names a flagged width's unit only when it is refreshed, as on opening.)
    w.saveProject();
    const QString before = slurp(path);
    again.accept();
    w.saveProject();
    ok("channels untouched: byte-identical save", slurp(path) == before);
  }

  // 7. The engine runs the GUI's experiment and prints what the page shows.
  {
    const QString path = work.filePath("engine.azr");
    spit(path, plain);
    w.open(path);
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    ThmExperimentsPage* p = ws.experimentsPage;
    p->addExperiment();
    segmentItem(p->segmentList, 1)->setCheckState(Qt::Checked);
    p->kinematicsBox->setChecked(true);
    p->beamCombo->setEditText("7Li");
    p->targetCombo->setEditText("d");
    p->spectatorCombo->setEditText("n");
    typeNumber(p->beamEnergyEdit, "19");
    const QString shown = p->derivedText();
    ws.accept();
    w.saveProject();
    int code = -1;
    const QString out = engineRun(work.path(), "engine.azr", &code);
    QRegExp rx("B\\(x\\+s\\) = ([-+0-9.eE]+) MeV; quasi-free E\\(x\\+A\\) = ([-+0-9.eE]+) MeV, "
               "E_qf = E\\(x\\+A\\) - B = ([-+0-9.eE]+) MeV");
    ok("engine accepts the GUI's experiment", code == 0 && out.contains("THM experiment 'E1'"), out.right(400));
    const bool found = rx.indexIn(out) >= 0;
    ok("engine prints the energies", found);
    if(found)
      ok("page shows the engine's numbers",
         shown.contains("B(x+s) = " + rx.cap(1) + " MeV") && shown.contains("E(x+A) = " + rx.cap(2) + " MeV") &&
             shown.contains("E_qf = E(x+A) - B = " + rx.cap(3) + " MeV"),
         shown + " | " + rx.cap(0));
    if(found)
      ok("page's B(x+s) and E_qf fields show the engine's numbers",
         p->bindingValue->text() == rx.cap(1) + " MeV" && p->qfValue->text() == rx.cap(3) + " MeV",
         p->bindingValue->text() + " | " + p->qfValue->text());
  }

  // 8. zeta at the ends of the data, charged spectator: 7Li(p,a) via
  //    3He(7Li, a a)d -- made up, the proton carried by 3He = p + d.
  const QString lsBlock = "experiment[E1] segments=1 beam=3He target=7Li spectator=d Ebeam=20 lineshape=on\n";
  const QString lsPath = work.filePath("lineshape.azr");
  spit(lsPath, plain + "<thm>\n" + lsBlock + "</thm>\n");
  {
    w.open(lsPath);
    ThmSettings s;
    QString err;
    ok("charged spectator: block opens", w.thmSettings(s, &err), err);
    ThmWorkspace ws(&w, s);
    ThmExperimentsPage* p = ws.experimentsPage;
    const QString shown = p->derivedText();
    double lo = 0, hi = 0;
    ok("charged spectator: data range read", p->pointRange(QList<int>() << 1, lo, hi) && hi > lo);
    int code = -1;
    const QString out = engineRun(work.path(), "lineshape.azr", &code);
    ok("charged spectator: engine accepts", code == 0, out.right(400));
    const QString report = slurp(work.filePath("output/thm_experiments.out"));
    QRegExp rowE("\\nE\\s+([-+0-9.eE]+)\\s+([-+0-9.eE]+)"), rowZ("zeta\\[\\d+\\]\\s+([-+0-9.eE]+)\\s+([-+0-9.eE]+)");
    const bool haveE = rowE.indexIn(report) >= 0, haveZ = rowZ.indexIn(report) >= 0;
    ok("charged spectator: engine reports E and zeta", haveE && haveZ, report.right(600));
    if(haveE && haveZ) {
      ok("data range = the engine's lowest/highest point",
         std::fabs(lo - rowE.cap(1).toDouble()) < 1e-5 * std::fabs(lo) &&
             std::fabs(hi - rowE.cap(2).toDouble()) < 1e-5 * std::fabs(hi),
         QString("%1 %2 | %3").arg(lo).arg(hi).arg(rowE.cap(0)));
      const QString z = QString::fromUtf8("\u03b6 = %1 \u2026 %2 into")
                            .arg(QString::number(rowZ.cap(1).toDouble(), 'g', 3), QString::number(rowZ.cap(2).toDouble(), 'g', 3));
      ok("page shows the engine's zeta at the ends", shown.contains(z) && rowZ.cap(1).toDouble() < 0.0,
         shown + " | " + rowZ.cap(0));
    }
  }


  // 10. The spectator-momentum window (ps=, psNodes=) on the Experiments page.
  auto type = [](QLineEdit* e, const QString& text) {
    e->setText(text);
    emit e->textEdited(text);
  };
  const QString reaction = "beam=7Li target=d spectator=n Ebeam=19";
  spit(work.filePath("psflat.dat"), "# p_s  w\n20 1\n40 1\n");
  {
    spit(fourPath, four);
    w.open(fourPath);
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    ThmExperimentsPage* p = ws.experimentsPage;
    p->addExperiment();
    segmentItem(p->segmentList, 1)->setCheckState(Qt::Checked);
    ok("ps: disabled without the reaction", !p->psBox->isEnabled());
    p->kinematicsBox->setChecked(true);
    p->beamCombo->setEditText("7Li");
    p->targetCombo->setEditText("d");
    p->spectatorCombo->setEditText("n");
    ok("ps: disabled with an incomplete reaction", !p->psBox->isEnabled());
    typeNumber(p->beamEnergyEdit, "19");
    ok("ps: enabled with the complete reaction", p->psBox->isEnabled());
    ok("ps: point by default, only the distribution shown",
       p->psKindCombo->currentData().toString() == "delta" && !p->psMinEdit->isVisibleTo(p) &&
           !p->psNodesSpin->isVisibleTo(p) && !p->psTableEdit->isVisibleTo(p) && !p->psFwhmEdit->isVisibleTo(p));
    p->psKindCombo->setCurrentIndex(p->psKindCombo->findData("hulthen"));
    ok("ps: Hulthen shows the window, a,b (standard, read-only) and the nodes",
       p->psMinEdit->isVisibleTo(p) && p->psAEdit->isVisibleTo(p) && !p->psAEdit->isEnabled() &&
           p->psAEdit->writtenText() == "0.2317" && p->psBEdit->writtenText() == "1.202" && p->psNodesSpin->isVisibleTo(p) &&
           p->psNodesSpin->value() == 16 && !p->psFwhmEdit->isVisibleTo(p) && !p->psTableEdit->isVisibleTo(p));
    typeNumber(p->psMinEdit, "0");
    typeNumber(p->psMaxEdit, "40");
    ok("ps: hulthen:0-40", p->records().at(0).ps == "hulthen:0-40" && p->records().at(0).psNodes.isEmpty(),
       p->records().at(0).ps);
    ok("ps: <T_s> shown", p->derivedText().contains("<T_s> = ") && p->derivedText().contains("16 Gauss-Legendre"),
       p->derivedText());
    p->psCustomCheck->setChecked(true);
    ok("ps: custom a,b editable", p->psAEdit->isEnabled() && p->psBEdit->isEnabled());
    typeNumber(p->psAEdit, "0.42");
    typeNumber(p->psBEdit, "1.2");
    ok("ps: hulthen:a,b:pmin-pmax", p->records().at(0).ps == "hulthen:0.42,1.2:0-40", p->records().at(0).ps);
    p->psNodesSpin->setValue(24);
    ok("ps: psNodes", p->records().at(0).psNodes == "24");
    p->psKindCombo->setCurrentIndex(p->psKindCombo->findData("gauss"));
    ok("ps: Gaussian shows FWHM, not a,b", p->psFwhmEdit->isVisibleTo(p) && !p->psAEdit->isVisibleTo(p));
    typeNumber(p->psFwhmEdit, "50");
    ok("ps: gauss:FWHM:pmin-pmax", p->records().at(0).ps == "gauss:50:0-40", p->records().at(0).ps);
    p->psKindCombo->setCurrentIndex(p->psKindCombo->findData("table"));
    ok("ps: table shows the file, not the window", p->psTableEdit->isVisibleTo(p) && p->psTableButton->isVisibleTo(p) &&
                                                       !p->psMinEdit->isVisibleTo(p) && p->psNodesSpin->isVisibleTo(p));
    type(p->psTableEdit, "psmissing.dat");
    ok("ps: table:file", p->records().at(0).ps == "table:psmissing.dat");
    ok("refused: unreadable table (engine's reader)", ws.validate().contains("ps: cannot read the ps table"), ws.validate());
    type(p->psTableEdit, "psflat.dat");
    ok("ps: table accepted, <T_s> over its range", ws.validate().isEmpty() && p->derivedText().contains("[20, 40]"),
       ws.validate() + " | " + p->derivedText());
    ok("ps: chosen file stored relative to the project",
       ThmExperimentsPage::projectRelative(work.filePath("psflat.dat"), work.path()) == "psflat.dat" &&
           ThmExperimentsPage::projectRelative("/elsewhere/ps.dat", work.path()) == "/elsewhere/ps.dat");
    p->psKindCombo->setCurrentIndex(p->psKindCombo->findData("delta"));
    ok("ps: point drops ps and psNodes", p->records().at(0).ps.isEmpty() && p->records().at(0).psNodes.isEmpty() &&
                                             p->psNodesSpin->value() == 16);
    p->psKindCombo->setCurrentIndex(p->psKindCombo->findData("hulthen"));
    p->psCustomCheck->setChecked(false);
    ok("ps: unticking custom restores the deuteron's a,b", p->records().at(0).ps == "hulthen:0-40" &&
                                                                p->psAEdit->writtenText() == "0.2317", p->records().at(0).ps);
    p->kinematicsBox->setChecked(false);
    ok("ps: dropped with the reaction", p->records().at(0).ps.isEmpty() && !p->psBox->isEnabled() &&
                                            p->psKindCombo->currentData().toString() == "delta");
  }
  {
    // Write, read back and keep verbatim, for each distribution.
    const QStringList values = QStringList() << "hulthen:0-40" << "hulthen:0.42,1.2:0-40 psNodes=24" << "gauss:50:10-40"
                                             << "table:psflat.dat" << "hulthen:30-30";
    for(const QString& v : values) {
      const QString value = v.section(' ', 0, 0), nodes = v.contains("psNodes=") ? v.section('=', -1) : QString();
      // Through the controls.
      spit(fourPath, four);
      w.open(fourPath);
      ThmSettings s;
      w.thmSettings(s);
      {
        ThmWorkspace ws(&w, s);
        ThmExperimentsPage* p = ws.experimentsPage;
        p->addExperiment();
        segmentItem(p->segmentList, 1)->setCheckState(Qt::Checked);
        p->kinematicsBox->setChecked(true);
        p->beamCombo->setEditText("7Li");
        p->targetCombo->setEditText("d");
        p->spectatorCombo->setEditText("n");
        typeNumber(p->beamEnergyEdit, "19");
        const QStringList f = value.split(':');
        p->psKindCombo->setCurrentIndex(p->psKindCombo->findData(f[0]));
        if(f[0] == "table") {
          type(p->psTableEdit, f[1]);
        } else {
          const QString window = f.last();
          typeNumber(p->psMinEdit, window.section('-', 0, 0));
          typeNumber(p->psMaxEdit, window.section('-', 1));
          if(f[0] == "gauss") typeNumber(p->psFwhmEdit, f[1]);
          if(f[0] == "hulthen" && f.size() == 3) {
            p->psCustomCheck->setChecked(true);
            typeNumber(p->psAEdit, f[1].section(',', 0, 0));
            typeNumber(p->psBEdit, f[1].section(',', 1));
          }
        }
        if(!nodes.isEmpty()) p->psNodesSpin->setValue(nodes.toInt());
        ok(qPrintable("ps " + v + ": accepted"), ws.validate().isEmpty(), ws.validate());
        ws.accept();
      }
      w.saveProject();
      const QString expect = "experiment[E1] segments=1 " + reaction + " ps=" + v + "\n";
      ok(qPrintable("ps " + v + ": written"), blockOf(slurp(fourPath)) == expect, blockOf(slurp(fourPath)));
      // The engine's parser takes it.
      ThmSettings back;
      QString err;
      ok(qPrintable("ps " + v + ": read back by the engine's parser"), w.thmSettings(back, &err), err);
      // Reopened: the controls show it; untouched it stays byte for byte.
      w.open(work.filePath("plain.azr"));
      w.open(fourPath);
      w.saveProject();
      const QString before = slurp(fourPath);
      w.thmSettings(back);
      {
        ThmWorkspace ws(&w, back);
        ThmExperimentsPage* p = ws.experimentsPage;
        ok(qPrintable("ps " + v + ": controls read back"),
           p->psKindCombo->currentData().toString() == value.section(':', 0, 0) && p->psText() == value &&
               p->psNodesSpin->value() == (nodes.isEmpty() ? 16 : nodes.toInt()) && p->psBox->isEnabled(),
           p->psText());
        ws.accept();
      }
      w.saveProject();
      ok(qPrintable("ps " + v + ": untouched, byte-identical"), slurp(fourPath) == before);
    }
    // As written by hand (spacing, key order, ps=delta): verbatim until edited.
    const QString block = "experiment[H] psNodes=8 segments=1   ps=hulthen:0.0-40.00 " + reaction + "  # by hand\n"
                          "experiment[D] segments=2 ps=delta\n";
    spit(fourPath, four + "<thm>\n" + block + "</thm>\n");
    w.open(fourPath);
    w.saveProject();
    ThmSettings s;
    QString err;
    ok("ps by hand: opens", w.thmSettings(s, &err), err);
    {
      ThmWorkspace ws(&w, s);
      ok("ps by hand: controls", ws.experimentsPage->psMinEdit->writtenText() == "0.0" &&
                                     ws.experimentsPage->psMaxEdit->writtenText() == "40.00" &&
                                     ws.experimentsPage->psNodesSpin->value() == 8);
      ws.accept();
    }
    w.saveProject();
    ok("ps by hand: untouched block verbatim", blockOf(slurp(fourPath)) == block, blockOf(slurp(fourPath)));
    w.thmSettings(s);
    {
      ThmWorkspace ws(&w, s);
      typeNumber(ws.experimentsPage->psMaxEdit, "30");
      ws.accept();
    }
    w.saveProject();
    ok("ps by hand: edited window rewritten, the rest as written",
       blockOf(slurp(fourPath)) == "experiment[H] segments=1 " + reaction + " ps=hulthen:0.0-30 psNodes=8\n"
                                   "experiment[D] segments=2 ps=delta\n",
       blockOf(slurp(fourPath)));
  }
  {
    // Refused: a window without the reaction (engine's check), and a window
    // with a spectator energy for its pair (the engine's startup rule).
    ThmExperimentsPage noKin(QStringList() << "experiment[W] segments=1 ps=hulthen:0-40",
                             w.getSegmentsTab()->getSegmentsDataModel(), w.getPairsTab()->getPairsModel(),
                             w.projectDirectory(), true);
    ok("refused: a window without the reaction", noKin.check().contains("needs the kinematics"), noKin.check());
    ok("refused: psNodes without a window",
       ThmSettings::checkExperimentLines(QStringList() << "experiment[W] segments=1 psNodes=8").contains("psNodes= needs a ps window"));
    const QString wBlock = "experiment[W] segments=1 " + reaction + " ps=hulthen:0-40\n";
    const QString path = work.filePath("psrefuse.azr");
    const QStringList energies = {"spectatorEnergy=0.1\n", "spectatorEnergy[5]=0.1\n"};
    for(const QString& energy : energies) {
      spit(path, plain + "<thm>\n" + energy + wBlock + "</thm>\n");
      w.open(path);
      ThmSettings s;
      w.thmSettings(s);
      ThmWorkspace ws(&w, s);
      const QString why = ws.validate();
      ok(qPrintable("refused: window with " + energy.trimmed()),
         why == "<thm> experiment[W]: a ps window and spectatorEnergy both set the spectator motion of entrance pair 5; "
                "use one (ps=delta keeps spectatorEnergy).",
         why);
      ok("refusal shown on the Experiments page", ws.pages->currentWidget() == ws.experimentsPage &&
                                                      ws.experimentsPage->messageLabel->text().contains("spectatorEnergy both"));
    }
    int code = -1;
    const QString out = engineRun(work.path(), "psrefuse.azr", &code);
    ok("the engine refuses it with the same words",
       code != 0 && out.contains("ERROR: <thm> experiment[W]: a ps window and spectatorEnergy both set the spectator "
                                 "motion of entrance pair 5; use one (ps=delta keeps spectatorEnergy)."),
       out.right(400));
    // The Model page's value, changed in the workspace, is what counts.
    spit(path, plain + "<thm>\n" + wBlock + "</thm>\n");
    w.open(path);
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    ok("window alone accepted", ws.validate().isEmpty(), ws.validate());
    ws.modelPage->spectatorEnergySpin->setValue(0.2);
    ok("refused after setting the spectator energy on the Model page",
       ws.validate().contains("a ps window and spectatorEnergy both"), ws.validate());
  }
  {
    // <T_s> on the page = the engine's (startup message), for a Hulthen window
    // and for a table relative to the project.
    const QStringList means = {"hulthen:0-40", "table:psflat.dat"};
    for(const QString& value : means) {
      const QString path = work.filePath("psmean.azr");
      spit(path, plain + "<thm>\nexperiment[W] segments=1 " + reaction + " ps=" + value + "\n</thm>\n");
      w.open(path);
      ThmSettings s;
      w.thmSettings(s);
      ThmWorkspace ws(&w, s);
      const QString shown = ws.experimentsPage->derivedText();
      int code = -1;
      const QString out = engineRun(work.path(), "psmean.azr", &code);
      QRegExp rx("mu_sx = ([-+0-9.eE]+) MeV, T_s = p_s\\^2/2mu_sx from ([-+0-9.eE]+) to ([-+0-9.eE]+) MeV, "
                 "<T_s> = ([-+0-9.eE]+) MeV");
      const bool found = code == 0 && rx.indexIn(out) >= 0;
      ok(qPrintable("ps " + value + ": the engine prints the window"), found, out.right(400));
      if(found)
        ok(qPrintable("ps " + value + ": page shows the engine's mu_sx, T_s range and <T_s>"),
           shown.contains("mu_sx = " + rx.cap(1) + " MeV") && shown.contains("from " + rx.cap(2) + " to " + rx.cap(3)) &&
               shown.contains("<T_s> = " + rx.cap(4) + " MeV"),
           shown + " | " + rx.cap(0));
      if(found)
        ok(qPrintable("ps " + value + ": the <T_s> field shows the engine's"),
           ws.experimentsPage->meanTsValue->text() == rx.cap(4) + " MeV", ws.experimentsPage->meanTsValue->text());
    }
    // THM_EXPERIMENTS_PNG=<file>: keep a rendering of the page, to look at it
    // (charged spectator, line shape and a Hulthen window).
    const QString path = work.filePath("pspng.azr");
    spit(path, plain + "<thm>\nexperiment[E1] segments=1 beam=3He target=7Li spectator=d Ebeam=20 lineshape=on "
                       "ps=hulthen:0-40\n</thm>\n");
    w.open(path);
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    ok("charged spectator with a window: accepted", ws.validate().isEmpty(), ws.validate());
    if(qEnvironmentVariableIsSet("THM_EXPERIMENTS_PNG")) {
      ws.pages->setCurrentWidget(ws.experimentsPage);  // at the default size
      ws.show();  // lays the page out, so that its editor can be scrolled to the window
      QApplication::processEvents();
      if(QScrollArea* a = ws.experimentsPage->findChild<QScrollArea*>()) a->ensureWidgetVisible(ws.experimentsPage->psBox);
      QPixmap shot(ws.size());
      ws.render(&shot);
      shot.save(qEnvironmentVariable("THM_EXPERIMENTS_PNG"));
    }
  }

  {
    // coulombIntegral=1 (Model page) with R(E) on a distorted a + A wave
    // (Experiments page) is refused, as the engine refuses it.
    const QString rPath = work.filePath("coulomb_consistency.azr");
    spit(rPath, plain + "<thm>\nexperiment[E1] segments=1 beam=3He target=7Li spectator=d Ebeam=20 "
                        "distortion=coulomb\n</thm>\n");
    w.open(rPath);
    ThmSettings rs;
    w.thmSettings(rs);
    ThmWorkspace rws(&w, rs);
    ok("R(E) without coulombIntegral: accepted", rws.validate().isEmpty(), rws.validate());
    rws.modelPage->coulombIntegralCheck->setChecked(true);
    ok("R(E) with coulombIntegral=1: refused", rws.validate().contains("would be counted twice"), rws.validate());
    rws.modelPage->coulombIntegralCheck->setChecked(false);
    ok("and accepted again", rws.validate().isEmpty(), rws.validate());
  }

#ifdef AZURE2_THM_DIAGNOSTICS
  // 9. Diagnostics page.
  {
    // A weight ramp on segment 1 as well.
    spit(work.filePath("ramp.dat"), "0.0 1.0\n1.0 4.0\n8.0 40.0\n");
    const QString diagPath = work.filePath("diag.azr");
    spit(diagPath, plain + "<thm>\n" + lsBlock + "weight[1]=ramp.dat\n</thm>\n");
    w.open(diagPath);
    w.saveProject();
    const QString saved = slurp(diagPath);
    ThmSettings s;
    QString err;
    ok("diagnostics: block opens", w.thmSettings(s, &err), err);
    ThmWorkspace ws(&w, s);
    ThmDiagnosticsPage* d = ws.diagnosticsPage;
    ok("diagnostics: page, one THM segment offered", d && d->segmentCombo->count() == 1 &&
                                                        d->segmentCombo->currentText().startsWith("experiment E1: segment 1"),
       d ? d->segmentCombo->currentText() : QString());

    // The snapshot is Accept + save, and leaves the project as it was.
    QLineEdit* b = qobject_cast<QLineEdit*>(ws.channelsPage->pairTable->cellWidget(0, 3));
    b->setText("2.3");
    ws.modelPage->kinematicsCombo->setCurrentText("triple");
    QString snap;
    ok("snapshot written", ws.projectSnapshot(snap, &err), err);
    w.saveProject();
    ok("snapshot leaves the project untouched", slurp(diagPath) == saved);
    ThmSettings s2;
    w.thmSettings(s2);
    {
      ThmWorkspace copy(&w, s2);
      qobject_cast<QLineEdit*>(copy.channelsPage->pairTable->cellWidget(0, 3))->setText("2.3");
      copy.modelPage->kinematicsCombo->setCurrentText("triple");
      copy.accept();
    }
    w.saveProject();
    ok("snapshot = the file Accept + save writes", snap == slurp(diagPath));
    spit(diagPath, saved);  // back to the original project
    w.open(diagPath);
  }
  {
    const QString diagPath = work.filePath("diag.azr");
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    ThmDiagnosticsPage* d = ws.diagnosticsPage;
    ok("diagnostics: computed", d->computeNow(), d->result().error);
    const ThmDiagnosticsResult& r = d->result();
    auto finite = [](const QVector<double>& v, bool positive) {
      if(v.isEmpty()) return false;
      for(double x : v)
        if(!std::isfinite(x) || (positive && x < 0.0)) return false;
      return true;
    };
    ok("diagnostics: grid over the data", r.energy.size() == 201 && std::fabs(r.energy.first() - r.eLo) < 1e-12 &&
                                              std::fabs(r.energy.last() - r.eHi) < 1e-9);
    bool vertexOk = !r.vertex.isEmpty();
    int nodes = 0;
    for(const ThmDiagnosticsResult::VertexGroup& g : r.vertex)
      for(const ThmDiagnosticsCurve& c : g.curves) {
        vertexOk = vertexOk && c.y.size() == r.energy.size() && finite(c.y, true);
        nodes += c.nodes.size();
      }
    ok("diagnostics: vertex curves finite", vertexOk && r.vertexMode == "constant" && !r.vertexComplex);
    ok("diagnostics: HOES and on-shell finite, positive somewhere",
       finite(r.hoes, true) && finite(r.onShell, true) && *std::max_element(r.hoes.begin(), r.hoes.end()) > 0.0 &&
           *std::max_element(r.onShell.begin(), r.onShell.end()) > 0.0);
    bool lsOk = r.lineshape && finite(r.zeta, false) && r.zeta.first() < 0.0 && !r.nc2.isEmpty();
    for(const ThmDiagnosticsCurve& c : r.nc2) lsOk = lsOk && finite(c.y, true);
    ok("diagnostics: line shape finite, zeta < 0 (Z_B < Z_F)", lsOk);
    // |N_C|^2 is the engine's factor at the level's pole and width.
    bool ncOk = !r.nc2.isEmpty();
    for(const ThmDiagnosticsCurve& c : r.nc2)
      for(int i = 0; i < r.energy.size(); i += 50)
        ncOk = ncOk && std::fabs(c.y[i] - ThmLineshapeFactorSq(r.zeta[i], c.pole - r.energy[i], c.width)) <= 1e-12 * c.y[i];
    ok("diagnostics: |N_C|^2 = exp[2 zeta arctan(2 (E_l - E)/G_l)]", ncOk);
    // The weight: the engine's log-linear interpolation of the table.
    bool wOk = r.weight.size() == r.energy.size();
    for(int i = 0; wOk && i < r.energy.size(); i += 40) {
      const double e = r.energy[i];
      const double expect = e <= 1.0 ? std::exp(std::log(4.0) * e) : 4.0 * std::exp(std::log(10.0) * (e - 1.0) / 7.0);
      wOk = std::fabs(r.weight[i] - expect) < 1e-9 * expect;
    }
    ok("diagnostics: weight = the table, log-linear", wOk);

    // Nodes: where the engine's own M_l changes sign.
    const QList<PairsData> pairList = w.getPairsTab()->getPairsModel()->getPairs();
    const PairsData& pr = pairList.at(4);
    const double mu = pr.lightM * pr.heavyM / (pr.lightM + pr.heavyM) * uconv;
    bool nodesOk = true;
    int engineNodes = 0;
    for(const ThmDiagnosticsResult::VertexGroup& g : r.vertex)
      for(const ThmDiagnosticsCurve& c : g.curves) {
        const int l = c.label.mid(4).toInt();
        auto M = [&](double e) { return ThmFormFactor(l, c.boundary, mu, e, pr.bindingEnergy, pr.channelRadius); };
        for(double e : c.nodes) {
          nodesOk = nodesOk && M(e - 1e-6) * M(e + 1e-6) < 0.0;
          std::cout << "        node of " << g.jpi.toStdString() << ", " << c.label.toStdString() << ": E = " << e
                    << " MeV (B_c = " << c.boundary << ")" << std::endl;
        }
        for(int i = 1; i < 4000; i++) {
          const double e0 = r.eLo + (r.eHi - r.eLo) * (i - 1) / 3999.0, e1 = r.eLo + (r.eHi - r.eLo) * i / 3999.0;
          if(M(e0) * M(e1) < 0.0) engineNodes++;
        }
      }
    ok("diagnostics: vertex nodes where the engine's M_l changes sign", nodesOk && nodes == engineNodes,
       QString("%1 nodes, engine %2").arg(nodes).arg(engineNodes));
    std::cout << "        (" << nodes << " vertex node(s) in " << r.vertex.size() << " J^pi group(s))" << std::endl;

    // The panels.
    ok("diagnostics: vertex panel", d->vertexGroupCombo->count() == r.vertex.size() &&
                                        d->vertexPlot->series().size() == r.vertex[0].curves.size() &&
                                        !d->vertexPlot->title().isEmpty());
    int markers = 0;
    for(const ThmDiagnosticsCurve& c : r.vertex[0].curves) markers += c.nodes.size();
    ok("diagnostics: nodes marked", d->vertexPlot->markers().size() == markers);
    ok("diagnostics: HOES panel, log scale, two curves", d->hoesPlot->logY() && d->hoesPlot->series().size() == 2);
    ok("diagnostics: line shape panels shown", d->lineshapePlot->isVisibleTo(d) && d->zetaPlot->isVisibleTo(d) &&
                                                   d->lineshapePlot->series().size() == r.nc2.size());
    ok("diagnostics: weight panel shown", d->weightPlot->isVisibleTo(d) && d->weightPlot->series().size() == 1);
    QPixmap shot(d->size().expandedTo(QSize(900, 800)));
    d->resize(shot.size());
    d->render(&shot);  // paints every panel once
    ok("diagnostics: page paints", !shot.isNull());

    // Off the GUI thread, with the busy bar.
    ThmSettings s3;
    spit(diagPath, plain + "<thm>\nexperiment[E1] segments=1\n</thm>\n");
    w.open(diagPath);
    w.thmSettings(s3);
    ThmWorkspace ws3(&w, s3);
    ThmDiagnosticsPage* d3 = ws3.diagnosticsPage;
    QEventLoop loop;
    QObject::connect(d3, &ThmDiagnosticsPage::computed, &loop, &QEventLoop::quit);
    QTimer::singleShot(300000, &loop, &QEventLoop::quit);
    d3->compute();
    const bool wasBusy = d3->busy() && d3->busyBar->isVisibleTo(d3) && !d3->computeButton->isEnabled();
    loop.exec();
    ok("diagnostics: Compute runs in a thread with the busy bar", wasBusy);
    ok("diagnostics: thread result", !d3->busy() && d3->result().error.isEmpty() && !d3->busyBar->isVisibleTo(d3) &&
                                         d3->computeButton->isEnabled(),
       d3->result().error);
    ok("diagnostics: no line shape, no weight: panels hidden",
       !d3->result().lineshape && d3->result().weight.isEmpty() && !d3->lineshapePlot->isVisibleTo(d3) &&
           !d3->weightPlot->isVisibleTo(d3));
  }
  {
    // 11. Diagnostics with a spectator-momentum window: charged spectator,
    //     line shape, weight and a Hulthen window (every panel).
    const QString line = "experiment[E1] segments=1 beam=3He target=7Li spectator=d Ebeam=20 lineshape=on ps=hulthen:0-40";
    const QString path = work.filePath("diagps.azr");
    spit(path, plain + "<thm>\n" + line + "\nweight[1]=ramp.dat\n</thm>\n");
    w.open(path);
    ThmSettings s;
    QString err;
    ok("window diagnostics: block opens", w.thmSettings(s, &err), err);
    ThmWorkspace ws(&w, s);
    ThmDiagnosticsPage* d = ws.diagnosticsPage;
    ok("window diagnostics: computed", d->computeNow(), d->result().error);
    const ThmDiagnosticsResult& r = d->result();
    double sum = 0.0;
    for(double x : r.nodeWeight) sum += x;
    ok("window: 16 nodes, weights sum to 1", r.window && r.nodeP.size() == 16 && r.nodeWeight.size() == 16 &&
                                                 std::fabs(sum - 1.0) < 1e-12);
    // The nodes and weights are the engine's BuildThmSpectatorWindow; mu_sx of p + d.
    std::vector<ThmExperiment> xs;
    const bool parsed = ParseThmExperimentLine(line.toStdString(), xs).empty() && xs.size() == 1;
    const ThmNuclide *np = ThmNuclide::Find(1, 1), *nd = ThmNuclide::Find(1, 2);
    const double muSx = np->mass * nd->mass / (np->mass + nd->mass) * 931.49410242;
    ThmSpectatorWindow win;
    const bool built = parsed && BuildThmSpectatorWindow(xs[0], muSx, win).empty();
    bool same = built && win.p.size() == (size_t)r.nodeP.size() && std::fabs(r.muSx - muSx) < 1e-12 * muSx;
    for(size_t k = 0; same && k < win.p.size(); k++)
      same = win.p[k] == r.nodeP[k] && win.weight[k] == r.nodeWeight[k] && win.es[k] == r.nodeTs[k];
    ok("window: nodes, weights, T_s and mu_sx are the engine's", same,
       QString("mu_sx %1 vs %2").arg(r.muSx, 0, 'g', 12).arg(muSx, 0, 'g', 12));
    ok("window: <T_s> is the engine's", built && std::fabs(r.meanTs - win.MeanEs()) < 1e-14 * win.MeanEs(),
       QString("%1 vs %2").arg(r.meanTs, 0, 'g', 15).arg(win.MeanEs(), 0, 'g', 15));
    std::cout << "        mu_sx = " << r.muSx << " MeV, <T_s> = " << r.meanTs << " MeV" << std::endl;
    // w(p) = |phi|^2 p^2 of the deuteron Hulthen function, unit area; the
    // engine's normalized weights are omega_k w(p_k) (checked above), so the
    // curve at the nodes is w(p_k) itself.
    auto hulthen = [](double p) {
      const double q2 = (p / hbarc) * (p / hbarc), phi = 1.0 / (0.2317 * 0.2317 + q2) - 1.0 / (1.202 * 1.202 + q2);
      return phi * phi * p * p;
    };
    double area = 0.0, lo = 1e300, hi = -1e300;
    for(int i = 0; i < r.windowP.size(); i++) {
      if(i > 0) area += 0.5 * (r.windowW[i] + r.windowW[i - 1]) * (r.windowP[i] - r.windowP[i - 1]);
      if(r.windowP[i] > 0.0) {
        lo = std::min(lo, r.windowW[i] / hulthen(r.windowP[i]));
        hi = std::max(hi, r.windowW[i] / hulthen(r.windowP[i]));
      }
    }
    bool atNodes = r.nodeW.size() == r.nodeP.size();
    for(int k = 0; atNodes && k < r.nodeP.size(); k++)
      atNodes = std::fabs(r.nodeW[k] / hulthen(r.nodeP[k]) - lo) < 1e-12 * lo;
    ok("window: w(p) is |phi|^2 p^2 over [0, 40] MeV/c, unit area, nodes on the curve",
       r.windowP.size() > 100 && r.windowP.first() == 0.0 && r.windowP.last() == 40.0 && std::fabs(area - 1.0) < 1e-12 &&
           hi - lo < 1e-12 * lo && atNodes,
       QString("area %1, spread %2").arg(area, 0, 'g', 15).arg((hi - lo) / lo));
    // <|M_l|^2> = sum_k w_k M_l(E; B + T_k)^2 with the engine's ThmFormFactor.
    const PairsData& pr = w.getPairsTab()->getPairsModel()->getPairs().at(4);
    const double mu = pr.lightM * pr.heavyM / (pr.lightM + pr.heavyM) * uconv;
    bool avgOk = !r.vertex.isEmpty(), qfOk = avgOk, filled = true;
    double worst = 0.0;
    int curves = 0, nodes = 0;
    for(const ThmDiagnosticsResult::VertexGroup& g : r.vertex)
      for(const ThmDiagnosticsCurve& c : g.curves) {
        curves++;
        const int l = c.label.mid(4).toInt();
        avgOk = avgOk && c.yWindow.size() == r.energy.size();
        if(!avgOk) break;
        const double scale = *std::max_element(c.y.begin(), c.y.end());
        for(int i = 0; i < r.energy.size(); i += 10) {
          double expect = 0.0;
          for(int k = 0; k < r.nodeP.size(); k++) {
            const double m = ThmFormFactor(l, c.boundary, mu, r.energy[i], pr.bindingEnergy + r.nodeTs[k], pr.channelRadius);
            expect += r.nodeWeight[k] * m * m;
          }
          const double qf = ThmFormFactor(l, c.boundary, mu, r.energy[i], pr.bindingEnergy, pr.channelRadius);
          worst = std::max(worst, std::fabs(c.yWindow[i] - expect) / scale);
          qfOk = qfOk && std::fabs(c.y[i] - qf * qf) <= 1e-12 * scale;
        }
        // The window fills the quasi-free nodes.
        for(double e : c.nodes) {
          nodes++;
          int at = 0;
          for(int i = 1; i < r.energy.size(); i++)
            if(std::fabs(r.energy[i] - e) < std::fabs(r.energy[at] - e)) at = i;
          filled = filled && c.yWindow[at] > c.y[at];
          std::cout << "        " << g.jpi.toStdString() << ", " << c.label.toStdString() << ": node at " << e
                    << " MeV, |M|^2 = " << c.y[at] << " -> <|M|^2> = " << c.yWindow[at] << " (max " << scale << ")"
                    << std::endl;
        }
      }
    ok("window: <|M_l|^2> = sum_k w_k M_l(p_xA(p_k))^2 (engine's ThmFormFactor)", avgOk && worst < 1e-10,
       QString("worst %1 of the maximum").arg(worst));
    ok("window: quasi-free curve unchanged (p_s = 0)", qfOk);
    ok("window: the average fills the vertex nodes", filled);
    // The panels: the window with its nodes, and the average next to each quasi-free curve.
    ok("window panel shown: w(p) and the nodes as points",
       d->windowPlot->isVisibleTo(d) && d->windowPlot->series().size() == 2 && !d->windowPlot->series()[0].symbols &&
           d->windowPlot->series()[1].symbols && d->windowPlot->series()[1].x == r.nodeP && !d->windowPlot->title().isEmpty());
    ok("vertex panel: <|M_l|^2> dashed next to each |M_l|^2",
       d->vertexPlot->series().size() == 2 * r.vertex[0].curves.size() &&
           d->vertexPlot->series()[1].style == Qt::DashLine && d->vertexPlot->series()[1].y == r.vertex[0].curves[0].yWindow);
    ok("status names <T_s> (details in its tooltip)", d->statusLabel->toolTip().contains("<T_s> = "), d->statusLabel->toolTip());
    ok("every panel shown", d->lineshapePlot->isVisibleTo(d) && d->zetaPlot->isVisibleTo(d) && d->weightPlot->isVisibleTo(d));
    std::cout << "        (" << curves << " vertex curve(s), " << nodes << " node(s))" << std::endl;
    // Shown: a J^pi group whose vertex has a node, to see it filled.
    for(int g = 0; g < r.vertex.size(); g++)
      if(!r.vertex[g].curves.isEmpty() && !r.vertex[g].curves[0].nodes.isEmpty()) {
        d->vertexGroupCombo->setCurrentIndex(g);
        break;
      }
    ws.resize(1000, 1250);
    ws.pages->setCurrentWidget(d);
    QPixmap shot(ws.size());
    ws.render(&shot);
    ok("window diagnostics: page paints", !shot.isNull());
    // THM_DIAGNOSTICS_PNG=<file>: keep the rendering, to look at it.
    if(qEnvironmentVariableIsSet("THM_DIAGNOSTICS_PNG")) shot.save(qEnvironmentVariable("THM_DIAGNOSTICS_PNG"));
    // THM_PNG_DIR=<dir>: every page at the workspace's default size (model.png,
    // experiments.png, channels.png, diagnostics.png), to look at them.
    if(qEnvironmentVariableIsSet("THM_PNG_DIR")) {
      const QDir dir(qEnvironmentVariable("THM_PNG_DIR"));
      ws.resize(900, 700);
      ws.show();
      const QList<QPair<QWidget*, QString>> shots = {{ws.modelPage, "model.png"},
                                                     {ws.experimentsPage, "experiments.png"},
                                                     {ws.channelsPage, "channels.png"},
                                                     {d, "diagnostics.png"}};
      for(const auto& s : shots) {
        ws.pages->setCurrentWidget(s.first);
        QApplication::processEvents();
        QPixmap page(ws.size());
        ws.render(&page);
        page.save(dir.filePath(s.second));
      }
    }
  }
  {
    // ps=table:<file> relative to the project: made absolute in the engine's
    // temporary copy only.
    const QString path = work.filePath("diagtable.azr");
    spit(path, plain + "<thm>\nexperiment[E1] segments=1 " + reaction + " ps=table:psflat.dat\n</thm>\n");
    w.open(path);
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    QString snap, err;
    ok("table: the snapshot keeps the relative path", ws.projectSnapshot(snap, &err) &&
                                                          snap.contains(" ps=table:psflat.dat\n"), err);
    ThmDiagnosticsPage* d = ws.diagnosticsPage;
    ok("table: computed from a copy elsewhere (path made absolute)", d->computeNow(), d->result().error);
    const ThmDiagnosticsResult& r = d->result();
    bool flat = r.window && !r.windowW.isEmpty() && r.windowP.first() == 20.0 && r.windowP.last() == 40.0;
    for(double v : r.windowW) flat = flat && std::fabs(v - 1.0 / 20.0) < 1e-12;
    for(double p : r.nodeP) flat = flat && p > 20.0 && p < 40.0;
    ok("table: flat w(p) = 1/20 per MeV/c on [20, 40], nodes inside", flat);
  }
#endif

  // 12. Distortion factor (distortion= ...) on a 12C+12C-like project: the
  //     four THM segments of examples/c12c12_tumino2018 (12C(14N,a/p)d at 30
  //     MeV) with free norms and nothing else (no pair 6, the 12C+12C pair
  //     without B that the other data use), no extrapolation, no folding.
  const QString c12Dir = work.filePath("c12");
  QDir(c12Dir).mkpath("data");
  QString c12;
  {
    const QString c12src = QString(AZURE2_SOURCE_DIR) + "/examples/c12c12_tumino2018";
    for(const QString& f : QDir(c12src + "/data").entryList(QDir::Files))
      QFile::copy(c12src + "/data/" + f, c12Dir + "/data/" + f);
    const QStringList in = slurp(c12src + "/c12c12_tumino2018.azr").split('\n');
    QStringList out;
    QString block;
    bool skip = false;
    for(const QString& line : in) {
      if(line.startsWith("<")) block = line;
      if(line == "<thm>") { skip = true; continue; }
      if(line == "</thm>") { skip = false; continue; }
      if(skip) continue;
      QStringList f = line.split(QRegExp("\\s+"), Qt::SkipEmptyParts);
      if(block == "<levels>" && f.size() > 5 && f[5] == "6") continue;  // pair 6 (12C+12C, no B): unused here
      if(block == "<segmentsData>" && f.size() > 10 && !line.startsWith("<")) {
        if(f[7] != "10") continue;  // the THM segments only
        f[9] = "1";                 // with a free norm
        out << f.join("  ");
      } else if((block == "<segmentsTest>" || block == "<targetInt>") && !line.startsWith("<")) {
        continue;  // no extrapolation, no folding
      } else {
        out << line;
      }
    }
    c12 = out.join('\n');
    if(!c12.endsWith('\n')) c12 += '\n';
  }
  ok("c12: project prepared", c12.contains("data/thm_a0.dat") && !c12.contains("<thm>"));
  spit(c12Dir + "/rtable.dat", "# E w\n0.5 1\n3.0 4\n");
  spit(c12Dir + "/rshort.dat", "0.5 1\n1.5 4\n");
  const QString c12Reaction = "beam=14N target=12C spectator=d Ebeam=30";
  const QString c12Path = c12Dir + "/c12.azr";
  auto c12Open = [&](const QString& block) {
    spit(c12Path, c12 + "<thm>\n" + block + "</thm>\n");
    w.open(c12Path);
  };
  // Drives the distortion controls to the given key=value tokens.
  auto setDistortion = [&](ThmExperimentsPage* p, const QString& tokens) {
    for(const QString& t : tokens.split(' ', Qt::SkipEmptyParts)) {
      const QString key = t.section('=', 0, 0), value = t.section('=', 1);
      if(key == "distortion") {
        const QString kind = value.startsWith("table:") ? "table" : value;
        p->distortionCombo->setCurrentIndex(p->distortionCombo->findData(kind));
        if(kind == "table") type(p->distortionTableEdit, value.mid(6));
      } else if(key == "opticalAA" || key == "opticalSF") {
        const int c = key == "opticalAA" ? 0 : 1;
        if(value == "plane" || value == "coulomb") {
          p->opticalCombo[c]->setCurrentIndex(p->opticalCombo[c]->findData(value));
        } else {
          p->opticalCombo[c]->setCurrentIndex(p->opticalCombo[c]->findData("ws"));
          p->setOpticalText(c, value);
        }
      } else if(key == "spectatorAngle") {
        const QString kind = value == "qf" ? "qf" : value.startsWith("cm:") ? "cm" : "lab";
        p->angleKindCombo->setCurrentIndex(p->angleKindCombo->findData(kind));
        if(kind != "qf") typeNumber(p->angleEdit, kind == "cm" ? value.mid(3) : value);
      } else if(key == "distortionRef") {
        typeNumber(p->distortionRefEdit, value);
      } else if(key == "distortionRatio") {
        p->ratioCombo->setCurrentIndex(p->ratioCombo->findData(value));
      } else if(key == "boundState") {
        p->boundCombo->setCurrentIndex(p->boundCombo->findData(value.section(':', 0, 0)));
        if(value.contains(':')) typeNumber(p->rminEdit, value.section(':', 1));
      }
    }
  };
  {
    // Enabling: coulomb and optical need the complete reaction, a table does not.
    c12Open("");
    ThmSettings s;
    QString err;
    ok("c12: opens", w.thmSettings(s, &err), err);
    ThmWorkspace ws(&w, s);
    ThmExperimentsPage* p = ws.experimentsPage;
    p->addExperiment();
    segmentItem(p->segmentList, 3)->setCheckState(Qt::Checked);
    auto itemEnabled = [p](const QString& kind) {
      QStandardItemModel* m = qobject_cast<QStandardItemModel*>(p->distortionCombo->model());
      return m && (m->item(p->distortionCombo->findData(kind))->flags() & Qt::ItemIsEnabled);
    };
    ok("distortion: none by default, only the combo shown",
       p->distortionCombo->currentData().toString() == "none" && p->distortionCombo->isVisibleTo(p) &&
           !p->angleKindCombo->isVisibleTo(p) && !p->opticalCombo[0]->isVisibleTo(p) &&
           !p->distortionTableEdit->isVisibleTo(p) && !p->distortionValue->isVisibleTo(p));
    ok("distortion: without the reaction coulomb and optical disabled, table enabled",
       !itemEnabled("coulomb") && !itemEnabled("optical") && itemEnabled("table") && itemEnabled("none"));
    p->kinematicsBox->setChecked(true);
    p->beamCombo->setEditText("14N");
    p->targetCombo->setEditText("12C");
    p->spectatorCombo->setEditText("d");
    ok("distortion: incomplete reaction, still disabled", !itemEnabled("coulomb") && !itemEnabled("optical"));
    typeNumber(p->beamEnergyEdit, "30");
    ok("distortion: complete reaction, enabled", itemEnabled("coulomb") && itemEnabled("optical"));
    p->distortionCombo->setCurrentIndex(p->distortionCombo->findData("coulomb"));
    ok("distortion: coulomb written, its rows shown, optical row hidden",
       p->records().at(0).distortion == "coulomb" && p->angleKindCombo->isVisibleTo(p) &&
           p->ratioCombo->isVisibleTo(p) && p->boundCombo->isVisibleTo(p) && p->distortionRefEdit->isVisibleTo(p) &&
           !p->opticalCombo[0]->isVisibleTo(p) && !p->distortionTableEdit->isVisibleTo(p) &&
           p->distortionValue->isVisibleTo(p));
    ok("distortion: quasi-free angle, the degrees disabled", !p->angleEdit->isEnabled());
    ok("distortion: R at the ends shown", p->distortionValue->text().contains(QString::fromUtf8(" … ")) &&
                                              p->derivedText().contains("R = "),
       p->distortionValue->text() + " | " + p->derivedText());
    p->distortionCombo->setCurrentIndex(p->distortionCombo->findData("optical"));
    ok("distortion: optical shows both channels, Edit only for Woods-Saxon",
       p->opticalCombo[0]->isVisibleTo(p) && p->opticalCombo[1]->isVisibleTo(p) &&
           p->opticalCombo[0]->currentData().toString() == "coulomb" && !p->opticalButton[0]->isEnabled());
    p->opticalCombo[1]->setCurrentIndex(p->opticalCombo[1]->findData("ws"));
    ok("distortion: Woods-Saxon starts at ten zeros, Edit enabled",
       p->records().at(0).opticalSF == "0,0,0,0,0,0,0,0,0,0" && p->opticalButton[1]->isEnabled());
    typeNumber(p->distortionRefEdit, "2.664");
    p->distortionCombo->setCurrentIndex(p->distortionCombo->findData("coulomb"));
    ok("distortion: coulomb drops opticalAA/SF, keeps E_ref",
       p->records().at(0).opticalSF.isEmpty() && p->records().at(0).distortionRef == "2.664");
    p->distortionCombo->setCurrentIndex(p->distortionCombo->findData("none"));
    ok("distortion: none drops every distortion key",
       p->records().at(0).distortion.isEmpty() && p->records().at(0).distortionRef.isEmpty() &&
           !p->distortionValue->isVisibleTo(p));
    p->distortionCombo->setCurrentIndex(p->distortionCombo->findData("coulomb"));
    p->kinematicsBox->setChecked(false);
    ok("distortion: coulomb dropped with the reaction, now disabled",
       p->records().at(0).distortion.isEmpty() && p->distortionCombo->currentData().toString() == "none" &&
           !itemEnabled("coulomb"));
    p->distortionCombo->setCurrentIndex(p->distortionCombo->findData("table"));
    type(p->distortionTableEdit, "rtable.dat");
    ok("distortion: a table without the reaction", p->records().at(0).distortion == "table:rtable.dat" &&
                                                       ws.validate().isEmpty() && p->distortionValue->isVisibleTo(p),
       ws.validate());
    ok("distortion: table stored relative to the project",
       ThmExperimentsPage::projectRelative(c12Dir + "/rtable.dat", c12Dir) == "rtable.dat");
    // The table's w(E) at the ends: log-linear between (0.5, 1) and (3, 4).
    double lo = 0, hi = 0;
    p->pointRange(QList<int>() << 3, lo, hi);
    auto w = [](double e) { return std::exp(std::log(4.0) * (e - 0.5) / 2.5); };
    ok("distortion: table value field is w at the ends",
       p->distortionValue->text() == QString::fromUtf8("%1 … %2").arg(QString::number(w(lo), 'g', 3),
                                                                            QString::number(w(hi), 'g', 3)),
       p->distortionValue->text());
  }
  {
    // The Woods-Saxon dialog: the ten numbers in their fields, as written until changed.
    ThmOpticalDialog d("test", "50.0,4.5,0.60,10,4.5,0.6,0,0,0,4.50");
    ok("optical dialog: fields as written", d.text() == "50.0,4.5,0.60,10,4.5,0.6,0,0,0,4.50" &&
                                                d.fields[0]->value() == 50.0 && d.fields[9]->value() == 4.5 &&
                                                d.fields[0]->suffix() == " MeV" && d.fields[1]->suffix() == " fm",
       d.text());
    typeNumber(d.fields[3], "12");
    ok("optical dialog: one number changed", d.text() == "50.0,4.5,0.60,12,4.5,0.6,0,0,0,4.50", d.text());
  }
  // Round trip of every form: set through the controls, written, read by the
  // engine's parser, shown again on reopening, untouched byte for byte.
  const QString c12Base = "experiment[E1] segments=3 " + c12Reaction;
  {
    const QStringList forms = {"distortion=coulomb",
                               "distortion=coulomb spectatorAngle=8",
                               "distortion=coulomb spectatorAngle=cm:90",
                               "distortion=coulomb distortionRef=2.664 distortionRatio=dw",
                               "distortion=coulomb boundState=yukawa:3",
                               "distortion=coulomb boundState=whittaker:3",
                               "distortion=coulomb boundState=yukawa",
                               "distortion=optical opticalAA=plane",
                               "distortion=optical opticalAA=50,4.5,0.6,10,4.5,0.6,0,0,0,4.5 opticalSF=plane",
                               "distortion=optical opticalSF=0,0,0,0,0,0,8,5.5,0.65,5 spectatorAngle=cm:8",
                               "distortion=table:rtable.dat"};
    for(const QString& form : forms) {
      c12Open(c12Base + "\n");
      ThmSettings s;
      w.thmSettings(s);
      {
        ThmWorkspace ws(&w, s);
        setDistortion(ws.experimentsPage, form);
        ok(qPrintable("distortion " + form + ": accepted"), ws.validate().isEmpty(), ws.validate());
        if(ws.validate().isEmpty()) ws.accept();
      }
      w.saveProject();
      if(form == "distortion=coulomb")
        ok("c12: a 15-character norm is written apart from the free-norm flag",
           slurp(c12Path).contains("3.659747237e-06 1 "), slurp(c12Path).section("<segmentsData>", 1).left(400));
      const QString expect = c12Base + " " + form + "\n";
      ok(qPrintable("distortion " + form + ": written"), blockOf(slurp(c12Path)) == expect, blockOf(slurp(c12Path)));
      ThmSettings back;
      QString err;
      ok(qPrintable("distortion " + form + ": read back by the engine's parser"), w.thmSettings(back, &err), err);
      w.open(work.filePath("plain.azr"));
      w.open(c12Path);
      w.saveProject();
      const QString before = slurp(c12Path);
      w.thmSettings(back);
      {
        ThmWorkspace ws(&w, back);
        const ThmExperimentRecord& r = ws.experimentsPage->records().at(0);
        const QString kind = r.distortion.startsWith("table:") ? "table" : r.distortion;
        ok(qPrintable("distortion " + form + ": controls read back"),
           ws.experimentsPage->distortionCombo->currentData().toString() == kind &&
               (kind != "table" || ws.experimentsPage->distortionTableEdit->text() == r.distortion.mid(6)) &&
               ws.experimentsPage->experimentLines() == QStringList(c12Base + " " + form),
           ws.experimentsPage->experimentLines().join("|"));
        ok(qPrintable("distortion " + form + ": reopened, accepted"), ws.validate().isEmpty(), ws.validate());
        if(!ws.validate().isEmpty()) continue;
        if(ws.validate().isEmpty()) ws.accept();
      }
      w.saveProject();
      ok(qPrintable("distortion " + form + ": untouched, byte-identical"), slurp(c12Path) == before);
    }
  }
  {
    // Written by hand: verbatim until edited; an edit rewrites only its key.
    const QString hand = c12Base + " distortion=optical   opticalAA=50.0,4.50,0.60,10,4.5,0.6,0,0,0,4.5 "
                                   "spectatorAngle=cm:8.0 boundState=whittaker:3.0 distortionRatio=dwpw  # by hand\n";
    c12Open(hand);
    w.saveProject();
    ThmSettings s;
    w.thmSettings(s);
    {
      ThmWorkspace ws(&w, s);
      ThmExperimentsPage* p = ws.experimentsPage;
      ok("distortion by hand: controls", p->distortionCombo->currentData().toString() == "optical" &&
                                             p->opticalCombo[0]->currentData().toString() == "ws" &&
                                             p->opticalCombo[1]->currentData().toString() == "coulomb" &&
                                             p->angleKindCombo->currentData().toString() == "cm" &&
                                             p->angleEdit->writtenText() == "8.0" && p->rminEdit->writtenText() == "3.0" &&
                                             p->ratioCombo->currentData().toString() == "dwpw");
      if(ws.validate().isEmpty()) ws.accept();
    }
    w.saveProject();
    ok("distortion by hand: untouched block verbatim", blockOf(slurp(c12Path)) == hand, blockOf(slurp(c12Path)));
    w.thmSettings(s);
    {
      ThmWorkspace ws(&w, s);
      ws.experimentsPage->ratioCombo->setCurrentIndex(ws.experimentsPage->ratioCombo->findData("dw"));
      if(ws.validate().isEmpty()) ws.accept();
    }
    w.saveProject();
    ok("distortion by hand: edited ratio rewritten, the rest as written",
       blockOf(slurp(c12Path)) == c12Base + " distortion=optical opticalAA=50.0,4.50,0.60,10,4.5,0.6,0,0,0,4.5 "
                                            "spectatorAngle=cm:8.0 distortionRatio=dw boundState=whittaker:3.0\n",
       blockOf(slurp(c12Path)));
  }
  {
    // Refusals, in the engine's words (the engine refuses the same files).
    struct Refusal {
      QString tokens, words;
      bool run;
    };
    const QList<Refusal> refusals = {
        {"distortion=optical opticalAA=1,2,3", "opticalAA='1,2,3': expected plane, coulomb or ten numbers", true},
        {"distortion=coulomb spectatorAngle=190", "spectatorAngle='190': expected qf, a lab angle", false},
        {"spectatorAngle=8", "spectatorAngle= needs distortion=coulomb or distortion=optical", false},
        {"distortion=coulomb opticalSF=plane", "opticalAA= and opticalSF= need distortion=optical", false},
        {"distortion=table:rshort.dat", "distortion: the points span E_cm = ", true},
        {"distortion=table:nothere.dat", "distortion: cannot read the weight file", false},
        {"distortion=coulomb spectatorAngle=60", "is beyond the reach of the spectator", true},
    };
    for(const Refusal& r : refusals) {
      c12Open(c12Base + " " + r.tokens + "\n");
      ThmSettings s;
      QString err;
      if(!w.thmSettings(s, &err)) {
        // Refused already by the parser: the engine's message, from its own parser.
        ok(qPrintable("refused: " + r.tokens), err.contains(r.words), err);
        if(r.run) {
          int code = -1;
          const QString out = engineRun(c12Dir, "c12.azr", &code);
          ok(qPrintable("the engine refuses " + r.tokens + " with the same words"), code != 0 && out.contains("ERROR: " + err),
             err + " | " + out.right(400));
        }
        continue;
      }
      ThmWorkspace ws(&w, s);
      const QString why = ws.validate();
      ok(qPrintable("refused: " + r.tokens), why.contains(r.words) && why.startsWith("<thm> experiment[E1]: "), why);
      ok(qPrintable("refusal shown on the page: " + r.tokens),
         ws.experimentsPage->messageLabel->isVisibleTo(ws.experimentsPage), ws.experimentsPage->messageLabel->text());
      if(r.run) {
        int code = -1;
        const QString out = engineRun(c12Dir, "c12.azr", &code);
        ok(qPrintable("the engine refuses " + r.tokens + " with the same words"), code != 0 && out.contains("ERROR: " + why),
           why + " | " + out.right(400));
      }
    }
    // On the page: a malformed Woods-Saxon list refused at once, in the parser's words.
    c12Open(c12Base + " distortion=optical\n");
    {
      ThmSettings s;
      w.thmSettings(s);
      ThmWorkspace ws(&w, s);
      ThmExperimentsPage* p = ws.experimentsPage;
      p->opticalCombo[0]->setCurrentIndex(p->opticalCombo[0]->findData("ws"));
      p->setOpticalText(0, "1,2,3");
      ok("optical: a short list refused on the page", p->messageLabel->text().startsWith("opticalAA='1,2,3': expected plane") &&
                                                          p->messageLabel->isVisibleTo(p) &&
                                                          ws.validate().contains("opticalAA='1,2,3'"),
         p->messageLabel->text());
      p->setOpticalText(0, "50,4.5,0.6,10,4.5,0.6,0,0,0,4.5");
      ok("optical: ten numbers accepted", ws.validate().isEmpty() && !p->messageLabel->isVisibleTo(p), ws.validate());
    }
    // Coulomb without the reaction: the engine's check.
    ThmExperimentsPage noKin(QStringList() << "experiment[W] segments=3 distortion=coulomb",
                             w.getSegmentsTab()->getSegmentsDataModel(), w.getPairsTab()->getPairsModel(), c12Dir, true);
    ok("refused: coulomb without the reaction", noKin.check().contains("distortion=coulomb needs the kinematics"),
       noKin.check());
  }
  {
    // R at the ends of the data on the page = the engine's (its startup line),
    // the B(x+s) note hidden (masses and field 32 agree).
    c12Open(c12Base + " distortion=coulomb\n");
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    ThmExperimentsPage* p = ws.experimentsPage;
    const QString shown = p->derivedText();
    int code = -1;
    const QString out = engineRun(c12Dir, "c12.azr", &code);
    QRegExp rx("R = ([-+0-9.eE]+) at E = ([-+0-9.eE]+) MeV \\(E_sF = [^)]*\\), ([-+0-9.eE]+) at E = ([-+0-9.eE]+) "
               "MeV \\(E_sF = [^)]*\\)\\.");
    const bool found = code == 0 && rx.indexIn(out) >= 0;
    ok("distortion: the engine runs it and prints R at the ends", found, out.right(600));
    if(found) {
      ok("distortion: the page shows the engine's R line", shown.contains(rx.cap(0)), shown + " | " + rx.cap(0));
      ok("distortion: the R field shows the engine's numbers",
         p->distortionValue->text() == QString::fromUtf8("%1 … %2")
                                           .arg(QString::number(rx.cap(1).toDouble(), 'g', 3),
                                                QString::number(rx.cap(3).toDouble(), 'g', 3)) &&
             rx.cap(1).toDouble() > 1.0 && rx.cap(3).toDouble() < 1.0,
         p->distortionValue->text());
      std::cout << "        engine: " << rx.cap(0).toStdString() << std::endl;
    }
    ok("B(x+s): masses and field 32 agree, no note", p->bindingWarningLabel->isHidden() &&
                                                         ws.channelsPage->bindingWarning(1).isEmpty() &&
                                                         !out.contains("B(x+s) from the masses"));
    // B edited on the Channels page: the note on both pages, and away again.
    ws.channelsPage->setBindingText(1, "5");
    p->refreshDerived();
    const QString note = QString::fromUtf8("B from masses 10.27 MeV ≠ pair B 5 MeV (vertex uses the pair value)");
    ok("B(x+s): note in the reaction section", !p->bindingWarningLabel->isHidden() && !p->bindingWarningIcon->isHidden() &&
                                                   p->bindingWarningLabel->text() == note,
       p->bindingWarningLabel->text());
    ok("B(x+s): warning icon beside the pair's B", ws.channelsPage->bindingWarning(1) == "experiment[E1]: " + note,
       ws.channelsPage->bindingWarning(1));
    ws.channelsPage->setBindingText(1, "10.272312");
    p->refreshDerived();
    ok("B(x+s): back to the file's B, no note", p->bindingWarningLabel->isHidden() &&
                                                    ws.channelsPage->bindingWarning(1).isEmpty());
    if(qEnvironmentVariableIsSet("THM_PNG_DIR")) {
      ws.resize(900, 700);
      ws.show();
      ws.pages->setCurrentWidget(p);
      QApplication::processEvents();
      if(QScrollArea* a = p->findChild<QScrollArea*>()) a->ensureWidgetVisible(p->distortionBox);
      QApplication::processEvents();
      QPixmap page(ws.size());
      ws.render(&page);
      page.save(QDir(qEnvironmentVariable("THM_PNG_DIR")).filePath("experiments_distortion.png"));
    }
  }
  {
    // The B(x+s) note in the engine's words: 3He = p + d against the pair's
    // 2.2246 MeV (7Li(p,a) through a made-up 3He(7Li,aa)d).
    w.open(lsPath);
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    ThmExperimentsPage* p = ws.experimentsPage;
    int code = -1;
    const QString out = engineRun(work.path(), "lineshape.azr", &code);
    const QString detail = p->bindingWarningLabel->toolTip();
    ok("B(x+s): note shown for 3He = p + d",
       !p->bindingWarningLabel->isHidden() &&
           p->bindingWarningLabel->text() ==
               QString::fromUtf8("B from masses 5.493 MeV ≠ pair B 2.225 MeV (vertex uses the pair value)"),
       p->bindingWarningLabel->text());
    ok("B(x+s): the note's detail is the engine's warning",
       code == 0 && !detail.isEmpty() && out.contains("WARNING: <thm> experiment[E1]: " + detail), detail + " | " + out.right(600));
    ok("B(x+s): Channels page marks pair 5", ws.channelsPage->bindingWarning(5).contains("B from masses 5.493 MeV"),
       ws.channelsPage->bindingWarning(5));
    ok("B(x+s): not a refusal", ws.validate().isEmpty(), ws.validate());
  }
#ifdef AZURE2_THM_DIAGNOSTICS
  {
    // Diagnostics: the distortion card, R the engine's (thm_distortion), the
    // page's R at the ends the same numbers.
    c12Open(c12Base + " distortion=coulomb\n");
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    ThmDiagnosticsPage* d = ws.diagnosticsPage;
    d->segmentCombo->setCurrentIndex(d->segmentCombo->findData(3));
    ok("distortion diagnostics: computed", d->computeNow(), d->result().error);
    const ThmDiagnosticsResult& r = d->result();
    const int n = r.energy.size();
    ok("distortion diagnostics: R, |M|^2, |M_PW|^2 on the grid", r.distortion && r.distortionError.isEmpty() &&
                                                                      r.distortionKind == "coulomb" && r.distortionR.size() == n &&
                                                                      r.distortionDirect.size() == n &&
                                                                      r.distortionM2.size() == n && r.distortionPW2.size() == n,
       r.distortionError);
    double worstModel = 0.0, worstRatio = 0.0;
    for(int i = 0; i < n && r.distortionDirect.size() == n && r.distortionM2.size() == n; i++) {
      worstModel = std::max(worstModel, std::fabs(r.distortionR[i] / r.distortionDirect[i] - 1.0));
      worstRatio = std::max(worstRatio, std::fabs(r.distortionM2[i] / r.distortionPW2[i] / r.distortionDirect[i] - 1.0));
    }
    ok("distortion diagnostics: R as the model uses it = R direct (grid, 1e-5)", n > 0 && worstModel < 1e-5,
       QString::number(worstModel));
    ok("distortion diagnostics: R = (|M|^2/|M_PW|^2) over its value at E_ref (dwpw)", n > 0 && worstRatio < 1e-10,
       QString::number(worstRatio));
    double lo = 0, hi = 0, rLo = 0, rHi = 0;
    QString why;
    ws.experimentsPage->distortionInfo(ws.experimentsPage->records().at(0), &why, &lo, &hi, &rLo, &rHi);
    ok("distortion: the page's R at the ends = the engine's thm_distortion there",
       why.isEmpty() && n > 0 && std::fabs(lo - r.energy.first()) < 1e-12 && std::fabs(hi - r.energy.last()) < 1e-9 &&
           std::fabs(rLo / r.distortionDirect.first() - 1.0) < 1e-9 && std::fabs(rHi / r.distortionDirect.last() - 1.0) < 1e-6,
       QString("%1 %2 | %3 %4").arg(rLo, 0, 'g', 12).arg(rHi, 0, 'g', 12).arg(r.distortionDirect.value(0), 0, 'g', 12)
           .arg(r.distortionDirect.value(n - 1), 0, 'g', 12));
    ok("distortion card: shown, log scale, R + |M|^2 + |M_PW|^2, E_ref marked",
       d->distortionPlot->isVisibleTo(d) && d->distortionPlot->logY() && d->distortionPlot->series().size() == 3 &&
           d->distortionPlot->series()[0].y == r.distortionR && d->distortionPlot->markers().size() == 1 &&
           std::fabs(d->distortionPlot->markers()[0].x - r.distortionRef) < 1e-12 && !d->distortionPlot->title().isEmpty());
    ok("distortion: the status tooltip names it", d->statusLabel->toolTip().contains("Distortion: coulomb"),
       d->statusLabel->toolTip());
    if(qEnvironmentVariableIsSet("THM_PNG_DIR")) {
      ws.resize(900, 700);
      ws.show();
      ws.pages->setCurrentWidget(d);
      QApplication::processEvents();
      QPixmap page(ws.size());
      ws.render(&page);
      page.save(QDir(qEnvironmentVariable("THM_PNG_DIR")).filePath("diagnostics_distortion.png"));
    }
  }
  {
    // distortion=dw: |M|^2 is R, so not drawn twice.
    c12Open(c12Base + " distortion=coulomb distortionRatio=dw\n");
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    ThmDiagnosticsPage* d = ws.diagnosticsPage;
    d->segmentCombo->setCurrentIndex(d->segmentCombo->findData(3));
    ok("distortion dw: computed", d->computeNow(), d->result().error);
    ok("distortion dw: R and |M_PW|^2", d->distortionPlot->series().size() == 2);
  }
  {
    // A table relative to the project: made absolute in the engine's copy only.
    c12Open(c12Base + " distortion=table:rtable.dat\n");
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    QString snap, err;
    ok("distortion table: the snapshot keeps the relative path",
       ws.projectSnapshot(snap, &err) && snap.contains(" distortion=table:rtable.dat\n"), err);
    ThmDiagnosticsPage* d = ws.diagnosticsPage;
    d->segmentCombo->setCurrentIndex(d->segmentCombo->findData(3));
    ok("distortion table: computed from a copy elsewhere (path made absolute)", d->computeNow(), d->result().error);
    const ThmDiagnosticsResult& r = d->result();
    bool same = r.distortionKind == "table" && r.distortionR.size() == r.energy.size();
    for(int i = 0; same && i < r.energy.size(); i++)
      same = std::fabs(r.distortionR[i] / std::exp(std::log(4.0) * (r.energy[i] - 0.5) / 2.5) - 1.0) < 1e-12;
    ok("distortion table: w(E) of the table, one curve", same && d->distortionPlot->series().size() == 1 &&
                                                             d->distortionPlot->isVisibleTo(d));
  }
  {
    // No distortion: the card stays hidden.
    c12Open(c12Base + "\n");
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    ThmDiagnosticsPage* d = ws.diagnosticsPage;
    d->segmentCombo->setCurrentIndex(d->segmentCombo->findData(3));
    ok("no distortion: computed", d->computeNow(), d->result().error);
    ok("no distortion: no card", !d->result().distortion && !d->distortionPlot->isVisibleTo(d));
  }
#endif

  // 13. The exit angle (theta=) on tests/7Li_p_a (Tumino's 50-70 degrees).
  {
    const QString path = work.filePath("theta.azr");
    auto thetaOpen = [&](const QString& block) {
      spit(path, plain + "<thm>\n" + block + "</thm>\n");
      w.open(path);
    };
    const QString hand = "experiment[E1] segments=1 theta=50.0-70   # Tumino 2006\n";
    thetaOpen(hand);
    ThmSettings s;
    QString err;
    ok("theta: block opens", w.thmSettings(s, &err), err);
    {
      ThmWorkspace ws(&w, s);
      ThmExperimentsPage* p = ws.experimentsPage;
      ok("theta: window shown as written", p->thetaCombo->currentData().toString() == "window" &&
                                               p->thetaMinEdit->writtenText() == "50.0" &&
                                               p->thetaMaxEdit->writtenText() == "70" &&
                                               p->thetaMinEdit->isVisibleTo(p) && p->thetaMinEdit->suffix() == QString::fromUtf8("°") &&
                                               p->thetaMaxEdit->maximum() == 180.0 && ws.validate().isEmpty(),
         ws.validate());
      ok("theta: tooltip names the angle and the supplementary window",
         p->thetaCombo->toolTip().contains("exit particle 1") && p->thetaCombo->toolTip().contains("supplementary"));
      ws.accept();
    }
    w.saveProject();
    ok("theta: untouched block verbatim", blockOf(slurp(path)) == hand, blockOf(slurp(path)));
    w.thmSettings(s);
    {
      ThmWorkspace ws(&w, s);
      typeNumber(ws.experimentsPage->thetaMaxEdit, "80");
      ok("theta: edited max accepted", ws.validate().isEmpty(), ws.validate());
      ws.accept();
    }
    w.saveProject();
    ok("theta: edited max rewritten, min as written", blockOf(slurp(path)) == "experiment[E1] segments=1 theta=50.0-80\n",
       blockOf(slurp(path)));
    w.thmSettings(s);
    {
      ThmWorkspace ws(&w, s);
      ThmExperimentsPage* p = ws.experimentsPage;
      p->thetaCombo->setCurrentIndex(p->thetaCombo->findData("all"));
      ok("theta: all hides the window", !p->thetaMinEdit->isVisibleTo(p) && !p->thetaMaxEdit->isVisibleTo(p));
      ws.accept();
    }
    w.saveProject();
    ok("theta: all (the default) not written", blockOf(slurp(path)) == "experiment[E1] segments=1\n", blockOf(slurp(path)));
    w.thmSettings(s);
    {
      ThmWorkspace ws(&w, s);
      ThmExperimentsPage* p = ws.experimentsPage;
      p->thetaCombo->setCurrentIndex(p->thetaCombo->findData("window"));
      ok("theta: a new window starts at 0-180", p->records().at(0).theta == "0-180", p->records().at(0).theta);
      typeNumber(p->thetaMaxEdit, "0");
      ws.accept();
    }
    w.saveProject();
    ok("theta: one angle, 0-0", blockOf(slurp(path)) == "experiment[E1] segments=1 theta=0-0\n", blockOf(slurp(path)));
    ok("theta: 0-0 passes the engine's parser", w.thmSettings(s, &err), err);
    // Refusals in the engine's words: a reversed window (the parser) ...
    thetaOpen("experiment[E1] segments=1\n");
    w.thmSettings(s);
    {
      ThmWorkspace ws(&w, s);
      ThmExperimentsPage* p = ws.experimentsPage;
      p->thetaCombo->setCurrentIndex(p->thetaCombo->findData("window"));
      typeNumber(p->thetaMinEdit, "90");
      typeNumber(p->thetaMaxEdit, "80");
      ok("theta: reversed window refused on the page", p->messageLabel->isVisibleTo(p) &&
                                                          p->messageLabel->text().startsWith("theta='90-80': expected all or thmin-thmax"),
         p->messageLabel->text());
      ok("theta: reversed window refused on Accept", ws.validate().startsWith("<thm> experiment[E1]: theta='90-80'"),
         ws.validate());
      typeNumber(p->thetaMaxEdit, "100");
      ok("theta: 90-100 accepted", ws.validate().isEmpty() && !p->messageLabel->isVisibleTo(p), ws.validate());
    }
    // ... and a window with entranceL=coherent, read live from the Model page.
    thetaOpen("experiment[E1] segments=1 theta=50-70\n");
    w.thmSettings(s);
    {
      ThmWorkspace ws(&w, s);
      ThmExperimentsPage* p = ws.experimentsPage;
      QStandardItemModel* m = qobject_cast<QStandardItemModel*>(p->thetaCombo->model());
      ok("theta: window offered with incoherent", m && (m->item(1)->flags() & Qt::ItemIsEnabled) && ws.validate().isEmpty());
      ws.modelPage->entranceLCombo->setCurrentText("coherent");
      ws.pages->setCurrentWidget(ws.channelsPage);
      ws.pages->setCurrentWidget(p);  // the workspace refreshes the page when it is shown
      const QString words = ThmExperimentsPage::coherentRefusal();
      ok("theta + coherent: refused on the page in the engine's words",
         p->messageLabel->isVisibleTo(p) && p->messageLabel->text() == words, p->messageLabel->text());
      ok("theta + coherent: refused on Accept", ws.validate() == "<thm> experiment[E1]: " + words, ws.validate());
      ok("theta + coherent: the window item stays selectable while current", m && (m->item(1)->flags() & Qt::ItemIsEnabled));
      p->thetaCombo->setCurrentIndex(p->thetaCombo->findData("all"));
      ok("theta + coherent: all accepted, window no longer offered",
         ws.validate().isEmpty() && !p->messageLabel->isVisibleTo(p) && m && !(m->item(1)->flags() & Qt::ItemIsEnabled) &&
             m->item(1)->toolTip() == words,
         ws.validate());
    }
    {
      // The engine refuses the same file with the same words.
      thetaOpen("entranceL=coherent\nexperiment[E1] segments=1 theta=50-70\n");
      w.thmSettings(s);
      ThmWorkspace ws(&w, s);
      const QString why = ws.validate();
      int code = -1;
      const QString out = engineRun(work.path(), "theta.azr", &code);
      ok("theta + coherent: the engine refuses it with the same words",
         why == "<thm> experiment[E1]: " + ThmExperimentsPage::coherentRefusal() && code != 0 && out.contains("ERROR: " + why),
         why + " | " + out.right(400));
    }
    if(qEnvironmentVariableIsSet("THM_PNG_DIR")) {
      thetaOpen("experiment[E1] segments=1 beam=7Li target=d spectator=n Ebeam=19 theta=50-70\n");
      w.thmSettings(s);
      ThmWorkspace ws(&w, s);
      ws.resize(900, 700);
      ws.show();
      ws.pages->setCurrentWidget(ws.experimentsPage);
      QApplication::processEvents();
      QPixmap page(ws.size());
      ws.render(&page);
      page.save(QDir(qEnvironmentVariable("THM_PNG_DIR")).filePath("experiments_theta.png"));
    }
  }
#ifdef AZURE2_THM_DIAGNOSTICS
  {
    // Diagnostics: the angular distribution at one energy.
    const QString path = work.filePath("theta.azr");
    spit(path, plain + "<thm>\nexperiment[E1] segments=1 theta=50-70\n</thm>\n");
    w.open(path);
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    ThmDiagnosticsPage* d = ws.diagnosticsPage;
    ok("angular: computed", d->computeNow(), d->result().error);
    const ThmDiagnosticsResult& r = d->result();
    bool finite = r.angular && r.angularError.isEmpty() && r.angle.size() == 19 && r.dsdo.size() == 19 &&
                  r.angle.first() == 0.0 && r.angle.last() == 180.0 && r.windowMean > 0.0;
    for(double v : r.dsdo) finite = finite && std::isfinite(v) && v >= 0.0;
    ok("angular: 19 angles over 0-180, finite, a window average", finite, r.angularError);
    ok("angular: at the middle of the data by default (rounded)",
       std::fabs(r.angularEnergy - 0.5 * (r.eLo + r.eHi)) <= 0.05 * (r.eHi - r.eLo),
       QString::number(r.angularEnergy, 'g', 12));
    ok("angular card: shown, window shaded, curve and average",
       d->angularPlot->isVisibleTo(d) && d->angularPlot->bands().size() == 1 && d->angularPlot->bands()[0].x0 == 50.0 &&
           d->angularPlot->bands()[0].x1 == 70.0 && d->angularPlot->series().size() == 2 &&
           d->angularPlot->series()[0].y == r.dsdo && d->angularPlot->series()[1].y.value(0) == r.windowMean &&
           !d->angularPlot->title().isEmpty());
    ok("angular: the energy box holds the energy, within the data",
       std::fabs(d->angularEnergyEdit->value() - r.angularEnergy) < 1e-9 &&
           d->angularEnergyEdit->minimum() <= r.eLo + 1e-9 && d->angularEnergyEdit->maximum() >= r.eHi - 1e-9,
       QString("%1 in [%2, %3]").arg(d->angularEnergyEdit->value(), 0, 'g', 12).arg(d->angularEnergyEdit->minimum(), 0, 'g', 12)
           .arg(d->angularEnergyEdit->maximum(), 0, 'g', 12));
    if(qEnvironmentVariableIsSet("THM_PNG_DIR")) {
      ws.resize(900, 700);
      ws.show();
      ws.pages->setCurrentWidget(d);
      QApplication::processEvents();
      if(QScrollArea* a = d->findChild<QScrollArea*>()) a->ensureWidgetVisible(d->angularPlot);
      QApplication::processEvents();
      QPixmap page(ws.size());
      ws.render(&page);
      page.save(QDir(qEnvironmentVariable("THM_PNG_DIR")).filePath("diagnostics_angular.png"));
    }
    // The next Compute uses the energy in the box: a point of the HOES grid
    // (not the middle, where the first one was).
    const int k = r.hoesEnergy.size() / 3;
    const double e = r.hoesEnergy.value(k), hoes = r.hoes.value(k) / r.hoesScale;
    d->angularEnergyEdit->setValue(e);
    ok("angular: an energy in the box marks the result as changed", d->statusLabel->text().contains("Compute again"),
       d->statusLabel->text());
    ok("angular: recomputed at the chosen energy", d->computeNow() && std::fabs(d->result().angularEnergy - e) < 1e-9,
       d->result().error);
    // 4 pi <dsigma/dOmega> over 0-180 is the angle-integrated HOES cross
    // section (docs, "Fixed-angle observable"); the window average is the
    // curve's average over the solid angle of the window.  1-degree steps, Simpson.
    ThmDiagnosticsRequest q;
    QString err;
    ok("angular: snapshot", ws.projectSnapshot(q.projectText, &err), err);
    q.projectDir = work.path();
    q.paramMask = w.GetConfig().paramMask;
    q.segment = 1;
    q.angularEnergy = e;
    q.angularPoints = 181;
    const ThmDiagnosticsResult fine = ComputeThmDiagnostics(q);
    auto simpson = [&](int a, int b) {  // integral of f sin(theta) dtheta over [a, b] degrees
      double sum = 0.0;
      const double h = M_PI / 180.0;
      for(int i = a; i <= b; i++) {
        const double w8 = (i == a || i == b) ? 1.0 : ((i - a) % 2 ? 4.0 : 2.0);
        sum += w8 * fine.dsdo[i] * std::sin(i * h);
      }
      return sum * h / 3.0;
    };
    const bool have = fine.error.isEmpty() && fine.angularError.isEmpty() && fine.dsdo.size() == 181;
    const double total = have ? 2.0 * M_PI * simpson(0, 180) : 0.0;
    const double mean = have ? simpson(50, 70) / (std::cos(50 * M_PI / 180) - std::cos(70 * M_PI / 180)) : 0.0;
    ok("angular: 4 pi <dsigma/dOmega> = the engine's angle-integrated HOES (1e-5)",
       have && hoes > 0.0 && std::fabs(total / hoes - 1.0) < 1e-5,
       QString("%1 vs %2 %3").arg(total, 0, 'g', 12).arg(hoes, 0, 'g', 12).arg(fine.error + fine.angularError));
    ok("angular: the window average = the curve's average over 50-70 (1e-5)",
       have && std::fabs(mean / fine.windowMean - 1.0) < 1e-5,
       QString("%1 vs %2").arg(mean, 0, 'g', 12).arg(fine.windowMean, 0, 'g', 12));
    std::cout << "        E = " << e << " MeV: 4 pi <dsigma/dOmega> = " << total << ", HOES = " << hoes
              << "; <50-70> = " << fine.windowMean << std::endl;
  }
  {
    // No theta: no card.
    w.open(work.filePath("plain.azr"));
    ThmSettings s;
    w.thmSettings(s);
    ThmWorkspace ws(&w, s);
    ThmDiagnosticsPage* d = ws.diagnosticsPage;
    ok("no theta: computed", d->computeNow(), d->result().error);
    ok("no theta: no angular card", !d->result().angular && !d->angularPlot->isVisibleTo(d));
  }
#endif

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
