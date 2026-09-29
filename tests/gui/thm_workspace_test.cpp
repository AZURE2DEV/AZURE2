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
//     quasi-free energies the page shows.
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
#include <QProcess>
#include <QPushButton>
#include <QRadioButton>
#include <QRegExp>
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
#include "ThmWorkspace.h"
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
                                            << "experiment[B] segments=4 ps=3 background=const"
                                            << "experiment[A] beam=7Li target=d spectator=n Ebeam=19";
    QList<ThmExperimentRecord> r = ThmExperimentRecord::read(lines);
    ok("records merged by name", r.size() == 2 && r[0].name == "A" && r[1].name == "B");
    ok("segment ranges expanded", r.size() == 2 && r[0].segments == (QList<int>() << 1 << 2 << 3));
    ok("unshown key kept aside", r.size() == 2 && r[1].extraTokens == QStringList("ps=3"));
    ok("unchanged records: lines verbatim", ThmExperimentRecord::compose(lines, r, r) == lines);
    QList<ThmExperimentRecord> e = r;
    e[1].background = "linear";
    ok("edited record keeps the key it does not show",
       ThmExperimentRecord::compose(lines, r, e) ==
           (QStringList() << lines[0] << "experiment[B] segments=4 background=linear ps=3" << lines[2]),
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
    ok("no THM: three pages, no diagnostics page", ws.pages->count() == 3);
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
    p->beamEnergyEdit->setText("19");
    emit p->beamEnergyEdit->textEdited("19");
    ok("derived energies shown", p->derivedLabel->text().contains("B(x+s) = 2.22") &&
                                     p->derivedLabel->text().contains("quasi-free E(x+A)"),
       p->derivedLabel->text());
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
                               p->spectatorCombo->currentText() == "n" && p->beamEnergyEdit->text() == "19");
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
    // A key the engine knows and the page does not show is kept as written
    // (only where the engine linked here knows one: lineshape, a later stage).
    std::vector<ThmExperiment> x;
    const bool known =
        ParseThmExperimentLine("experiment[L] segments=1 beam=7Li target=d spectator=n Ebeam=19 lineshape=on", x).empty();
    if(known) {
      const QString block = "experiment[L] segments=1 lineshape=on beam=7Li target=d spectator=n Ebeam=19\n";
      spit(fourPath, four + "<thm>\n" + block + "</thm>\n");
      w.open(fourPath);
      ThmSettings s;
      QString err;
      ok("future key: block opens", w.thmSettings(s, &err), err);
      ThmWorkspace ws(&w, s);
      ws.experimentsPage->backgroundCombo->setCurrentText("const");
      ws.accept();
      w.saveProject();
      ok("future key kept through an edit",
         blockOf(slurp(fourPath)) ==
             "experiment[L] segments=1 background=const beam=7Li target=d spectator=n Ebeam=19 lineshape=on\n",
         blockOf(slurp(fourPath)));
    } else {
      std::cout << "  skip  future key through the workspace (the linked engine knows no key the page hides)"
                << std::endl;
    }
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
    p->beamEnergyEdit->setText("19");
    emit p->beamEnergyEdit->textEdited("19");
    const QString shown = p->derivedLabel->text();
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
  }

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
