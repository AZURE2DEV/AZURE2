// The Potential panel of the Particle Pairs tab and the Levels tab's check of the partial waves (docs/upgrade/qt_gui.md
// section 19), headless:  gui_potential_test <source dir> [<screenshot dir>]
// On a copy of tests/potential_background/models/bound_level_a40.azr (3He + 4He, the Gaussian potential VD in seven
// waves, each with `background interior`, the p3/2 one with `subtract 1`; a 3/2- level at the potential's bound state):
//  (a) the pair selected: the panel shows its record and its seven waves; Apply without an edit and a save move no byte;
//  (b) off and on: the record becomes comment lines, the project stays valid, on again gives the text byte for byte;
//  (c) each mode writes the record as the engine reads it: exterior only, exterior + interior (the B kept), none
//      (only bound_state) refused with the engine's words, capture interior from its box, central (the line removed);
//  (d) the pole table of the p3/2 wave equals the facade's (front::PotentialWaveOf, itself equal to the session's
//      AZUREAPI numbers: api_front (h)), the subtracted pole marked; the 3/2- level beside it;
//  (e) tuning: the depth for a0 = 40 fm (scipy's brentq 65.28518406 MeV), used and applied, the read-out then a0 = 40;
//      the depth binding the 3/2- level: the written 83.589554 MeV with the potential's ANC 4.0932 beside the level's;
//  (f) the Levels tab: without the 1/2+ level the panel warns that the s wave is left out; the check lists 1/2+ and its
//      blank level brings it back (widths 0, fixed); with includeAllJPi switched on from the check the s wave is in,
//      and the save writes includeAllJPi=1 and maxOrbitalMomentum;
//  (g) screenshots of the panel and of the check.
// PASS/FAIL per check; exit 1 on a failure.
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include <QAbstractButton>
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QElapsedTimer>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QPushButton>
#include <QTableView>
#include <QTableWidget>
#include <QTimer>

#include "AZURESetup.h"
#include "ChannelsModel.h"
#include "LevelsModel.h"
#include "LevelsTab.h"
#include "PairsTab.h"
#include "PartialWavesDialog.h"
#include "PotentialPanel.h"
#include "ProjectHolder.h"

namespace fs = std::filesystem;
namespace front = azure2::front;

namespace {

int failures = 0;

void Check(bool ok, const std::string &what) {
  std::cout << (ok ? "PASS " : "FAIL ") << what << std::endl;
  if (!ok) failures++;
}

std::string Bytes(const fs::path &p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

fs::path CopyProject(const fs::path &projectFile, const fs::path &to) {
  fs::remove_all(to);
  fs::create_directories(to);
  fs::copy(projectFile.parent_path(), to, fs::copy_options::recursive | fs::copy_options::skip_symlinks);
  for (const char *d : {"output", "checks"}) fs::create_directories(to / d);
  return to / projectFile.filename();
}

void Spin(int ms = 50) {
  QElapsedTimer t;
  t.start();
  while (t.elapsed() < ms) QApplication::processEvents(QEventLoop::AllEvents, 10);
}

/// Until the panel's computations end (60 s at most).
bool Settle(PotentialPanel *p) {
  QElapsedTimer t;
  t.start();
  Spin(30);
  while (p->busy() && t.elapsed() < 60000) Spin(30);
  Spin(30);
  return !p->busy();
}

/// The table row of the wave (l, J as shown); -1 when none.
int WaveRow(PotentialPanel *p, int l, const QString &J) {
  QTableWidget *t = p->wavesTable();
  for (int r = 0; r < t->rowCount(); r++)
    if (t->item(r, 0) && t->item(r, 0)->text() == QString::number(l) && t->item(r, 1)->text() == J) return r;
  return -1;
}

/// The wave fields of the project's record `name`.
std::vector<std::string> Waves(const front::Project &p, const std::string &name) {
  std::vector<std::string> out;
  for (const front::Record &r : p.Records("potentials"))
    if (r.name == name)
      for (const front::Field &f : r.fields)
        if (f.key == "wave") out.push_back(f.value);
  return out;
}

bool Has(const std::vector<std::string> &v, const std::string &s) {
  for (const std::string &x : v)
    if (x == s) return true;
  return false;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: gui_potential_test <source dir> [<screenshot dir>]" << std::endl;
    return 2;
  }
  QApplication app(argc, argv);
  QCoreApplication::setAttribute(Qt::AA_DontUseNativeDialogs, true);
  QCoreApplication::setOrganizationName("jina-test");
  QCoreApplication::setApplicationName("azure2-gui-potential-test");
  const fs::path source = argv[1];
  const QString shots = argc > 2 ? QString::fromLocal8Bit(argv[2]) : QString();
  if (!shots.isEmpty()) fs::create_directories(shots.toStdString());
  const fs::path work = fs::absolute("gui_potential_work");
  const fs::path file = CopyProject(source / "tests/potential_background/models/bound_level_a40.azr", work / "model");
  const std::string original = Bytes(file);
  // message boxes (none expected) are answered so that no modal dialog blocks the test
  QTimer boxes;
  boxes.setInterval(50);
  QObject::connect(&boxes, &QTimer::timeout, []() {
    if (QMessageBox *m = qobject_cast<QMessageBox *>(QApplication::activeModalWidget())) {
      std::cout << "note: a message box: " << m->text().toStdString() << std::endl;
      m->reject();
    }
  });
  boxes.start();

  AZURESetup w;
  w.resize(1500, 1000);
  w.show();
  QString why;
  Check(w.readProjectFile(QString::fromStdString(file.string()), why), "opened bound_level_a40 (" + why.toStdString() + ")");
  w.tabs()->setCurrentWidget(w.getPairsTab());
  PotentialPanel *panel = w.potentialPanel();
  Check(panel != nullptr, "the Pairs tab has the Potential panel");
  w.getPairsTab()->getPairsView()->selectRow(0);
  Spin(100);
  Settle(panel);

  // (a)
  Check(panel->pair() == 1 && panel->hasRecord() && panel->recordOn() && panel->onCheck()->isChecked(),
        "(a) pair 1 selected: its potential he3he4, on");
  int lines = 0;
  for (int r = 0; r < panel->wavesTable()->rowCount(); r++) {
    QComboBox *mode = qobject_cast<QComboBox *>(panel->wavesTable()->cellWidget(r, 2));
    lines += mode && mode->currentIndex() == PotentialPanel::ExteriorInterior;
  }
  Check(panel->wavesTable()->rowCount() == 7 && lines == 7, "(a) seven waves, each exterior + interior (" +
                                                                std::to_string(panel->wavesTable()->rowCount()) + ")");
  Check(panel->apply(why) && w.writeProjectFile(QString::fromStdString(file.string()), why) && Bytes(file) == original,
        "(a) Apply without an edit and a save: byte for byte");

  // (b) off and on
  panel->onCheck()->setChecked(false);
  Spin(50);
  const std::string offText = w.holder()->project().Text();
  Check(!panel->recordOn() && w.holder()->project().Records("potentials").empty() &&
            offText.find("#off potential he3he4") != std::string::npos && w.holder()->project().Validate().ok,
        "(b) off: the record kept as comment lines, the project valid");
  Check(panel->messagesLabel()->text().contains("switched off"), "(b) the panel says the potential is off");
  panel->onCheck()->setChecked(true);
  Spin(50);
  Check(panel->recordOn() && w.holder()->project().Text() == original, "(b) on again: the text byte for byte");
  Settle(panel);

  // (c) the modes
  const int s = WaveRow(panel, 0, "1/2");
  Check(s >= 0, "(c) the s wave's row");
  panel->setMode(s, PotentialPanel::Exterior);
  {
    std::vector<std::string> edited;
    for (const front::Field &f : panel->editedRecord().fields)
      if (f.key == "wave") edited.push_back(f.value);
    Check(Has(edited, "l 0 J 1/2 depth 67.5 MeV width 2.5200031651 fm"), "(c) exterior only: the line without background");
  }
  panel->setMode(WaveRow(panel, 0, "1/2"), PotentialPanel::ExteriorInterior);
  Check(panel->apply(why) && w.holder()->project().Text() == original,
        "(c) exterior + interior again before Apply: B kept, the text as it was");
  panel->setMode(WaveRow(panel, 0, "1/2"), PotentialPanel::Exterior);
  Check(panel->apply(why) && Has(Waves(w.holder()->project(), "he3he4"), "l 0 J 1/2 depth 67.5 MeV width 2.5200031651 fm"),
        "(c) exterior only applied: written so (" + why.toStdString() + ")");
  panel->setMode(WaveRow(panel, 0, "1/2"), PotentialPanel::ExteriorInterior);
  Check(panel->apply(why) && Has(Waves(w.holder()->project(), "he3he4"),
                                 "l 0 J 1/2 depth 67.5 MeV width 2.5200031651 fm background interior boundary 0"),
        "(c) exterior + interior on a line without one: background interior boundary 0 (B to be set)");
  {
    front::Project back = front::Project::FromText(original, file.string());
    for (const front::Record &r : back.Records("potentials")) w.holder()->setRecord("potentials", r, why);
    Spin(50);
  }
  {
    const int r = WaveRow(panel, 2, "5/2");
    QTableWidgetItem *cap = panel->wavesTable()->item(r, 7);
    cap->setCheckState(Qt::Checked);
    Spin(20);
    Check(panel->apply(why) &&
              Has(Waves(w.holder()->project(), "he3he4"),
                  "l 2 J 5/2 depth 69.0 MeV width 2.5200031651 fm background interior boundary 0.0 capture interior"),
          "(c) capture interior from its box: written");
    panel->wavesTable()->item(WaveRow(panel, 2, "5/2"), 7)->setCheckState(Qt::Unchecked);
    Spin(20);
    Check(panel->apply(why) && w.holder()->project().Text() == original, "(c) capture interior off: the text as it was");
  }
  {
    const int r = WaveRow(panel, 3, "7/2");
    panel->wavesTable()->setCurrentCell(r, 2);
    panel->setMode(r, PotentialPanel::PointCoulomb);
    Settle(panel);
    const QString m = panel->messagesLabel()->text();
    Check(m.contains("E_POTENTIAL_WAVE") && m.contains("only bound_state"),
          "(c) none (only bound_state) without a capture interior: the engine's refusal shown: " + m.toStdString().substr(0, 200));
    panel->setMode(WaveRow(panel, 3, "7/2"), PotentialPanel::Central);
    Check(panel->apply(why) && Waves(w.holder()->project(), "he3he4").size() == 6,
          "(c) central: the line removed (" + std::to_string(Waves(w.holder()->project(), "he3he4").size()) + " waves)");
    panel->revertButton()->click();
    Spin(20);
    w.holder()->project();
  }
  // back to the opened text (the central edit was applied: restore the file's record)
  {
    front::Project back = front::Project::FromText(original, file.string());
    for (const front::Record &r : back.Records("potentials")) w.holder()->setRecord("potentials", r, why);
    Spin(50);
    Check(w.holder()->project().Text() == original, "(c) the record restored");
  }
  Settle(panel);

  // (d) the poles of the p3/2 wave
  {
    const int r = WaveRow(panel, 1, "3/2");
    panel->wavesTable()->setCurrentCell(r, 2);
    Settle(panel);
    front::PotentialWaveReport ref;
    const front::PairPotential pp = front::PotentialOfPair(w.holder()->project(), 1);
    std::string refWhy;
    front::PotentialWaveOf(w.holder()->project(), pp.record, 1, 1, 3, -1000.0, 30.0, ref, refWhy);
    QTableWidget *t = panel->polesTable();
    bool same = panel->reportOk() && t->rowCount() == (int)ref.poles.size() && !ref.poles.empty();
    int subtracted = -1;
    for (int i = 0; same && i < t->rowCount(); i++) {
      same = t->item(i, 0)->text().toInt() == ref.poles[i].nodes &&
             t->item(i, 1)->text() == QString::number(ref.poles[i].energy, 'g', 7) &&
             t->item(i, 3)->text() == QString::number(ref.poles[i].gamma2, 'g', 5);
      if (t->item(i, 5)->text().startsWith("subtracted")) subtracted = t->item(i, 0)->text().toInt();
    }
    Check(same && subtracted == 1, "(d) p3/2: " + std::to_string(t->rowCount()) + " poles as the facade gives them, node 1 "
                                   "marked subtracted");
    Check(panel->levelsTable()->rowCount() == 1 && panel->levelsTable()->item(0, 3)->text().contains("4.093154"),
          "(d) the 3/2- level beside the poles with its ANC");
    if (!shots.isEmpty()) panel->window()->grab().save(shots + "/potential_panel_p32.png");

    // (e) binding the 3/2- level
    Check(panel->bindButton()->isEnabled() && panel->bindCombo()->count() == 1, "(e) the 3/2- level offered to bind");
    panel->bindButton()->click();
    Settle(panel);
    const front::PotentialDepth &d = panel->depthFound();
    Check(d.found && std::fabs(d.depth / 83.589554 - 1) < 1e-6 && std::fabs(d.anc - 4.0931604450) < 2e-4 &&
              panel->tuningLabel()->text().contains("4.09315"),
          "(e) the depth binding 3/2-: " + std::to_string(d.depth) + " MeV, the potential's ANC " + std::to_string(d.anc) +
              " beside the level's");
  }
  {
    const int r = WaveRow(panel, 0, "1/2");
    panel->wavesTable()->setCurrentCell(r, 2);
    Settle(panel);
    Check(panel->reportOk() && std::fabs(panel->report().a - 30.953078) < 1e-3 && panel->readOutLabel()->text().contains("30.95"),
          "(e) the s wave's own a0 read out: " + std::to_string(panel->report().a) + " fm (scipy 30.953078)");
    panel->lengthEdit()->setText("40");
    panel->lengthButton()->click();
    Settle(panel);
    const front::PotentialDepth d = panel->depthFound();
    Check(d.found && std::fabs(d.depth / 65.28518406 - 1) < 1e-6 && d.boundStates == 2,
          "(e) the depth for a0 = 40 fm: " + std::to_string(d.depth) + " MeV (brentq 65.28518406), two bound states kept");
    if (!shots.isEmpty()) panel->window()->grab().save(shots + "/potential_panel_tuning.png");
    panel->useDepthButton()->click();
    Settle(panel);
    Check(std::fabs(panel->report().a - 40.0) < 1e-4, "(e) used: the read-out gives a0 = " + std::to_string(panel->report().a));
    Check(panel->apply(why) && Waves(w.holder()->project(), "he3he4")[0].find("depth 65.2851") != std::string::npos,
          "(e) applied: " + Waves(w.holder()->project(), "he3he4")[0]);
    front::Project back = front::Project::FromText(original, file.string());
    for (const front::Record &rec : back.Records("potentials")) w.holder()->setRecord("potentials", rec, why);
    Spin(50);
  }

  // (f) the Levels tab: the 1/2+ level removed
  LevelsTab *levels = w.getLevelsTab();
  {
    const QList<LevelsData> all = levels->getLevelsModel()->getLevels();
    int row = -1;
    for (int i = 0; i < all.size(); i++)
      if (all[i].jValue == 0.5 && all[i].piValue == 1) row = i;
    levels->removeLevels({row}, false);
    Spin(50);
    Settle(panel);
    const QString m = panel->messagesLabel()->text();
    Check(m.contains("Left out of the calculation") && m.contains("l = 0 (1/2+)"),
          "(f) the panel warns: the s wave is left out (no 1/2+ level, includeAllJPi off)");
    if (!shots.isEmpty()) {
      w.tabs()->setCurrentWidget(w.getPairsTab());
      Spin(100);
      w.grab().save(shots + "/potential_panel_left_out.png");
    }
    levels->showPartialWaves();
    Spin(50);
    PartialWavesDialog *d = levels->partialWavesDialog();
    QListWidget *list = d->list();
    Check(list->count() >= 1 && list->item(0)->text().startsWith("1/2") && d->explanation()->text().contains("includeAllJPi"),
          "(f) the check lists 1/2+ and names includeAllJPi: " + (list->count() ? list->item(0)->text().toStdString() : ""));
    if (!shots.isEmpty()) d->grab().save(shots + "/levels_partial_waves.png");
    for (int i = 0; i < list->count(); i++) list->item(i)->setCheckState(i == 0 ? Qt::Checked : Qt::Unchecked);
    const int before = list->count();
    d->addBlankButton()->click();
    Spin(50);
    int blank = -1;
    const QList<LevelsData> now = levels->getLevelsModel()->getLevels();
    for (int i = 0; i < now.size(); i++)
      if (now[i].jValue == 0.5 && now[i].piValue == 1) blank = i;
    bool zero = blank >= 0 && now[blank].isFixed == 1 && now[blank].energy == 50.0;
    int channels = 0;
    for (const ChannelsData &c : levels->getChannelsModel()->getChannels())
      if (c.levelIndex == blank) {
        channels++;
        zero = zero && c.reducedWidth == 0.0 && c.isFixed == 1;
      }
    Check(zero && channels == 1 && d->list()->count() == before - 1,
          "(f) a blank 1/2+ level at 50 MeV: its width 0, energy and width fixed; the list one shorter");
    Settle(panel);
    Check(!panel->messagesLabel()->text().contains("Left out"), "(f) the panel's warning gone");
    // the other remedy: without the level, includeAllJPi switched on from the check
    levels->removeLevels({blank}, false);
    Spin(50);
    levels->showPartialWaves();
    d->switchOnButton()->click();
    Spin(50);
    Settle(panel);
    Check(levels->includeAllJPi() && !panel->messagesLabel()->text().contains("Left out") &&
              panel->wavesTable()->item(WaveRow(panel, 0, "1/2"), 8)->text() == "yes",
          "(f) includeAllJPi on from the check: the s wave is in the calculation");
    Check(w.writeProjectFile(QString::fromStdString(file.string()), why) &&
              w.holder()->configValue("includeAllJPi") == "1" &&
              w.holder()->configValue("maxOrbitalMomentum") == QString::number(levels->maxL()),
          "(f) saved: includeAllJPi=1 and maxOrbitalMomentum=" + std::to_string(levels->maxL()));
    front::Project saved;
    std::string w2;
    Check(front::Project::Open(file.string(), saved, w2) && saved.Validate().ok, "(f) the saved project is valid");
    if (!shots.isEmpty()) {
      w.tabs()->setCurrentWidget(levels);
      Spin(100);
      w.grab().save(shots + "/levels_tab.png");
      w.tabs()->setCurrentWidget(w.getPairsTab());
      Spin(100);
      w.grab().save(shots + "/pairs_tab_potential.png");
    }
  }
  std::cout << (failures ? "FAILED " : "all passed ") << failures << std::endl;
  return failures ? 1 : 0;
}
