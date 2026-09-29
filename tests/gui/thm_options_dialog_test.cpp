// Headless test of the THM workspace's Model page (the former THM Options dialog).
//
// The page edits the optional <thm> block through ThmSettings.  Checked here:
//  1. ThmSettings parses what the engine accepts and refuses what it refuses;
//  2. compose() keeps comments and unchanged lines, rewrites changed values,
//     drops defaults and appends new keys;
//  3. through AZURESetup on tests/7Li_p_a: settings set on the dialog's widgets
//     are saved as the expected block, a reopened project shows them again in
//     the dialog, an untouched dialog leaves the block byte for byte, defaults
//     remove the block, validate() refuses bad weight files;
//  4. the engine accepts the GUI's block, weight file included: with a free
//     THM norm a constant weight leaves the chi2 unchanged.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFile>
#include <QLineEdit>
#include <QProcess>
#include <QRegExp>
#include <QSpinBox>
#include <QString>
#include <QStringList>
#include <QTableWidget>
#include <QTemporaryDir>
#include <iostream>
#include "AZURESetup.h"
#include "Config.h"
#include "ThmModelPage.h"
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
// The lines between <thm> and </thm> of a saved project ("<none>" if absent).
static QString blockOf(const QString& text) {
  int a = text.indexOf("<thm>\n");
  if(a < 0) return "<none>";
  a += 6;
  int b = text.indexOf("</thm>\n", a);
  return b < 0 ? "<unterminated>" : text.mid(a, b - a);
}

// Runs "Calculate Segments From Data" on dir/name; the total chi2 printed, or "".
static QString engineChi2(const QString& dir, const QString& name) {
  QDir(dir).mkpath("output");
  QDir(dir).mkpath("checks");
  QProcess p;
  p.setWorkingDirectory(dir);
  p.setProcessChannelMode(QProcess::MergedChannels);
  QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
  if(!env.contains("OMP_NUM_THREADS")) env.insert("OMP_NUM_THREADS", "1");
  p.setProcessEnvironment(env);
  p.start(AZURE2_BINARY, QStringList() << "--no-gui" << "--no-readline" << name);
  if(!p.waitForStarted(30000)) return QString();
  p.write("1\n\n\n7\n");
  p.closeWriteChannel();
  if(!p.waitForFinished(600000)) { p.kill(); return QString(); }
  const QString out = QString::fromUtf8(p.readAll());
  QRegExp rx("Total Chi-Squared: ([-+0-9.eE]+)");
  return (p.exitCode() == 0 && rx.indexIn(out) >= 0) ? rx.cap(1) : QString();
}

int main(int argc, char** argv) {
  QTemporaryDir settingsDir;
  qputenv("XDG_CONFIG_HOME", settingsDir.path().toUtf8());
  QApplication app(argc, argv);
  QCoreApplication::setOrganizationName("AZURE2-tests");
  QCoreApplication::setApplicationName("thm_options_dialog_test");

  // 1. Parsing, with the engine's rules.
  {
    ThmSettings s;
    QString err;
    bool good = ThmSettings::parse(QStringList() << "  # notes" << "vertex = real  # alias"
                                                 << "kinematics=lambda32" << "coulombIntegral=on"
                                                 << "spectatorEnergy=0.4" << "spectatorEnergy[ 5]=0.6"
                                                 << "weight[2]=R.dat" << "weightTest[1]=/abs/R.dat"
                                                 << "entranceL=coherent" << "",
                                   s, &err);
    ok("a full block parses", good, err);
    ok("vertex=real is perlevel", s.vertex == "perlevel");
    ok("values read", s.kinematics == "lambda32" && s.coulombIntegral && s.spectatorEnergy == 0.4 &&
                          s.entranceL == "coherent");
    ok("per-pair spectator energy read", s.spectatorByPair.value(5, -1) == 0.6);
    ok("weights read", s.weight.value(2) == "R.dat" && s.weightTest.value(1) == "/abs/R.dat");
    const char* refused[] = {"vertx=onshell", "kinematics=kf2body", "spectatorEnergy=-0.1",
                             "coulombIntegral=2", "weight[0]=a.dat", "weight[1]=", "weight=a.dat",
                             "weight[1x]=a.dat", "vertex"};
    for(const char* line : refused) {
      ThmSettings r;
      ok(qPrintable(QString("refused like the engine: ") + line),
         !ThmSettings::parse(QStringList() << line, r, &err) && err.contains(line), err);
    }
    ok("defaults are default", ThmSettings().isDefault() && ThmSettings().keyValues().isEmpty());
  }

  // 2. compose().
  {
    const QStringList old = QStringList() << "# my notes" << "kinematics=kf3body   # data / KF3"
                                          << "" << "  spectatorEnergy[5] = 0.5" << "entranceL=incoherent";
    ThmSettings s;
    ThmSettings::parse(old, s);
    s.kinematics = "triple";  // changed: rewritten, inline comment kept
    s.vertex = "onshell";     // new: appended
    // entranceL=incoherent is the default: its line goes; spectatorEnergy[5] unchanged: verbatim.
    const QStringList want = QStringList() << "# my notes" << "kinematics=triple   # data / KF3" << ""
                                           << "  spectatorEnergy[5] = 0.5" << "vertex=onshell";
    QStringList got = s.compose(old);
    ok("compose keeps comments and unchanged lines, drops defaults", got == want, got.join("|"));
    ok("compose of nothing lists the non-default keys",
       s.compose(QStringList()) ==
           (QStringList() << "vertex=onshell" << "kinematics=triple" << "spectatorEnergy[5]=0.5"));
  }

  // 3. Through AZURESetup.
  QTemporaryDir work;
  const QString src = QString(AZURE2_SOURCE_DIR) + "/tests/7Li_p_a";
  QString plain = slurp(src + "/7Li_p_a.azr");
  ok("found tests/7Li_p_a", !plain.isEmpty());
  if(!plain.endsWith('\n')) plain += '\n';
  QDir(work.path()).mkpath("data");
  QFile::copy(src + "/data/tumino_thm.dat", work.filePath("data/tumino_thm.dat"));
  spit(work.filePath("two.dat"), "# constant weight\n0 2\n10 2\n");
  spit(work.filePath("decr.dat"), "0 1\n5 2\n4 3\n");
  spit(work.filePath("plain.azr"), plain);
  const QString blockIn = "# 7Li(p,a) options\nkinematics=kf3body   # data / KF3\n";
  spit(work.filePath("block.azr"), plain + "<thm>\n" + blockIn + "</thm>\n");

  AZURESetup w;

  // An untouched dialog changes nothing, with or without a block.
  for(const char* name : {"plain.azr", "block.azr"}) {
    const QString path = work.filePath(name);
    w.open(path);
    w.saveProject();
    const QString before = slurp(path);
    ThmSettings s;
    ok("the block of the project parses", w.thmSettings(s));
    ThmModelPage d(s, w.projectDirectory());
    w.setThmSettings(d.settings());
    w.saveProject();
    ok(qPrintable(QString("untouched dialog: byte-identical save of ") + name), slurp(path) == before);
  }
  ok("the block read is kept verbatim", blockOf(slurp(work.filePath("block.azr"))) == blockIn,
     blockOf(slurp(work.filePath("block.azr"))));

  // Set options through the dialog's widgets, save, check the block.
  {
    w.open(work.filePath("block.azr"));
    ThmSettings s;
    w.thmSettings(s);
    ThmModelPage d(s, w.projectDirectory());
    ok("dialog shows the block's kinematics", d.kinematicsCombo->currentText() == "kf3body");
    d.vertexCombo->setCurrentText("onshell");
    d.coulombIntegralCheck->setChecked(true);
    d.spectatorEnergySpin->setValue(0.25);
    d.addSpectatorRow(5, 0.5);
    d.addWeightRow(false, 1, "two.dat");
    d.addWeightRow(true, 2, "two.dat");
    ThmSettings edited = d.settings();
    ok("validate accepts readable weight files", edited.validate(w.projectDirectory()).isEmpty(),
       edited.validate(w.projectDirectory()));
    w.setThmSettings(edited);
    w.saveProject();
    const QString want =
        "# 7Li(p,a) options\nkinematics=kf3body   # data / KF3\nvertex=onshell\ncoulombIntegral=1\n"
        "spectatorEnergy=0.25\nspectatorEnergy[5]=0.5\nweight[1]=two.dat\nweightTest[2]=two.dat\n";
    const QString saved = slurp(work.filePath("block.azr"));
    ok("edited block written as expected", blockOf(saved) == want, blockOf(saved));
    ok("written after </targetInt>", saved.contains("</targetInt>\n<thm>\n"));
  }
  // Reopen: the dialog shows what was saved.
  {
    w.open(work.filePath("plain.azr"));  // something else in between
    w.open(work.filePath("block.azr"));
    ThmSettings s;
    ok("saved block parses", w.thmSettings(s));
    ThmModelPage d(s, w.projectDirectory());
    ok("reopened: vertex", d.vertexCombo->currentText() == "onshell");
    ok("reopened: kinematics", d.kinematicsCombo->currentText() == "kf3body");
    ok("reopened: entranceL", d.entranceLCombo->currentText() == "incoherent");
    ok("reopened: coulombIntegral", d.coulombIntegralCheck->isChecked());
    ok("reopened: spectatorEnergy", d.spectatorEnergySpin->value() == 0.25);
    ok("reopened: per-pair table", d.spectatorTable->rowCount() == 1 &&
                                       qobject_cast<QSpinBox*>(d.spectatorTable->cellWidget(0, 0))->value() == 5 &&
                                       qobject_cast<QLineEdit*>(d.spectatorTable->cellWidget(0, 1))->text() == "0.5");
    bool rows = d.weightTable->rowCount() == 2;
    if(rows) {
      rows = qobject_cast<QComboBox*>(d.weightTable->cellWidget(0, 0))->currentIndex() == 0 &&
             qobject_cast<QSpinBox*>(d.weightTable->cellWidget(0, 1))->value() == 1 &&
             qobject_cast<QLineEdit*>(d.weightTable->cellWidget(0, 2))->text() == "two.dat" &&
             qobject_cast<QComboBox*>(d.weightTable->cellWidget(1, 0))->currentIndex() == 1 &&
             qobject_cast<QSpinBox*>(d.weightTable->cellWidget(1, 1))->value() == 2;
    }
    ok("reopened: weight table", rows);
    ok("reopened: dialog gives back the same settings", d.settings() == s);

    // validate refuses what the engine refuses at startup.
    ThmSettings bad = s;
    bad.weight[1] = "nothere.dat";
    ok("validate: missing weight file", bad.validate(w.projectDirectory()).contains("cannot read"),
       bad.validate(w.projectDirectory()));
    bad.weight[1] = "decr.dat";
    ok("validate: non-increasing energies",
       bad.validate(w.projectDirectory()).contains("strictly increasing"));
    bad.weight[1] = work.filePath("two.dat");
    ok("validate: absolute path", bad.validate(QString()).isEmpty());
  }

  // 4. The engine takes the GUI's block: a constant weight on the profiled THM
  // norm leaves chi2 as it is (the same block without the weights).
  {
    const QString withWeights = slurp(work.filePath("block.azr"));
    ThmSettings s;
    w.thmSettings(s);
    s.weight.clear();
    s.weightTest.clear();
    w.setThmSettings(s);
    w.saveProject();
    spit(work.filePath("noweight.azr"), slurp(work.filePath("block.azr")));
    spit(work.filePath("weighted.azr"), withWeights);
    const QString c0 = engineChi2(work.path(), "noweight.azr");
    const QString c1 = engineChi2(work.path(), "weighted.azr");
    ok("engine runs the GUI's block", !c0.isEmpty() && !c1.isEmpty(), c0 + " / " + c1);
    ok("constant weight: same chi2 as without it", c0 == c1, c0 + " / " + c1);
    w.open(work.filePath("weighted.azr"));
  }

  // All defaults: the block goes (comments included).
  {
    w.setThmSettings(ThmSettings());
    w.saveProject();
    const QString saved = slurp(work.filePath("weighted.azr"));
    ok("defaults: no <thm> block written", !saved.contains("<thm>"));
    w.open(work.filePath("plain.azr"));
    w.saveProject();
    ok("defaults: the file is the plain project's save", saved == slurp(work.filePath("plain.azr")));
  }

  // 5. THM experiment lines (experiment[<name>] ...) are not the dialog's to
  // edit: they are accepted, kept verbatim through any edit, and a block that
  // has one is not removed when every option is back to its default.
  {
    ThmSettings s;
    QString err;
    ok("experiment lines parse",
       ThmSettings::parse(QStringList() << "experiment[E1] segments=1 background=linear  # note"
                                        << "vertex=onshell",
                          s, &err),
       err);
    ok("experiment lines kept aside", s.experimentLines.size() == 1 && s.vertex == "onshell");
    ok("a block with only an experiment line is not default",
       ThmSettings::parse(QStringList() << "experiment[E1] segments=1", s) && !s.isDefault());

    const QString expBlock = "# grouped segments\nexperiment[E1] segments=1   # the only THM line\n"
                             "kinematics=kf3body\n";
    spit(work.filePath("exp.azr"), plain + "<thm>\n" + expBlock + "</thm>\n");
    w.open(work.filePath("exp.azr"));
    w.saveProject();
    const QString before = slurp(work.filePath("exp.azr"));
    ok("experiment block kept verbatim on save", blockOf(before) == expBlock, blockOf(before));
    ThmSettings cur;
    ok("experiment block parses in the GUI", w.thmSettings(cur, &err), err);
    {
      ThmModelPage d(cur, w.projectDirectory());
      ok("dialog round trip keeps the experiment lines", d.settings() == cur);
      w.setThmSettings(d.settings());
      w.saveProject();
      ok("untouched dialog: byte-identical save with an experiment", slurp(work.filePath("exp.azr")) == before);
      d.kinematicsCombo->setCurrentText("lacognata");  // every option back to its default
      w.setThmSettings(d.settings());
      w.saveProject();
    }
    const QString saved = slurp(work.filePath("exp.azr"));
    ok("options at default: block kept for its experiment line",
       blockOf(saved) == "# grouped segments\nexperiment[E1] segments=1   # the only THM line\n", blockOf(saved));
    // A one-segment experiment without background is the per-segment profile.
    const QString c0 = engineChi2(work.path(), "plain.azr");
    const QString c1 = engineChi2(work.path(), "exp.azr");
    ok("engine runs the block with the experiment line, same chi2", !c0.isEmpty() && c0 == c1, c0 + " / " + c1);
  }

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
