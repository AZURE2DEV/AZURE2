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
                                            << "experiment[B] segments=4 distortion=3 background=const"
                                            << "experiment[A] beam=7Li target=d spectator=n Ebeam=19";
    QList<ThmExperimentRecord> r = ThmExperimentRecord::read(lines);
    ok("records merged by name", r.size() == 2 && r[0].name == "A" && r[1].name == "B");
    ok("segment ranges expanded", r.size() == 2 && r[0].segments == (QList<int>() << 1 << 2 << 3));
    ok("unshown key kept aside", r.size() == 2 && r[1].extraTokens == QStringList("distortion=3"));
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
           (QStringList() << lines[0] << "experiment[B] segments=4 background=linear distortion=3" << lines[2]),
       ThmExperimentRecord::compose(lines, r, e).join("|"));
    e = r;
    e[0].background = "quadratic";  // two lines of A become one, in place of the first
    ok("edited record: one line at its first place",
       ThmExperimentRecord::compose(lines, r, e) ==
           (QStringList() << "experiment[A] segments=1-3 background=quadratic beam=7Li target=d spectator=n Ebeam=19"
                          << lines[1]),
       ThmExperimentRecord::compose(lines, r, e).join("|"));
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

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
