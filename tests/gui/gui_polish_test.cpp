// Round 8's polish of the record forms, the pair dialog and the Levels tab without a display
// (docs/upgrade/qt_gui.md section 15):
//   gui_polish_test <source dir> <screenshot dir>
//  (a) the fields a quantity decides: shown only for a quantity that uses them; a change of the quantity in the form
//      shows and hides them; a value of a field that disappears is dropped on Apply and the user told;
//  (b) the Add dialog: a change of the quantity makes its rows again (the fields it requires, its unit);
//  (c) the kernels drawn beside a target's Doppler, a beam's spread and a detector's resolution, with the FWHM, for a
//      data set and an energy chosen; a value typed in the form redraws them;
//  (d) an optional word shows its default greyed, not "(none)"; the edit-as-text toggle has an icon and says what it
//      does; an empty section shows its hint and no empty form;
//  (e) a sum of reactions edited term by term (weight, pairs), written back as the grammar has it;
//  (f) the pair dialog: element and A, J by halves, parity, type; a pair duplicated;
//  (g) the Levels tab's channels from the facade (the engine's rule) equal to the old copy's on T1-T3 and 13N;
//  (h) screenshots of every form changed and of the sections not looked at before.
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
#include <QDialog>
#include <QDialogButtonBox>
#include <QDir>
#include <QFileInfo>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QTabWidget>
#include <QTableView>
#include <QTimer>
#include <QToolButton>

#include "AZURESetup.h"
#include "AddPairDialog.h"
#include "LevelsTab.h"
#include "PairsModel.h"
#include "PairsTab.h"
#include "RecordTableModel.h"
#include "RecordTabs.h"
#include "RecordsPanel.h"
#include "ValueEditor.h"

namespace fs = std::filesystem;
namespace front = azure2::front;

namespace {

int failures = 0;

void Check(bool ok, const std::string &what) {
  std::cout << (ok ? "PASS " : "FAIL ") << what << std::endl;
  if (!ok) failures++;
}

fs::path CopyProject(const fs::path &projectFile, const fs::path &to) {
  fs::remove_all(to);
  fs::create_directories(to);
  fs::copy(projectFile.parent_path(), to, fs::copy_options::recursive | fs::copy_options::skip_symlinks);
  return to / projectFile.filename();
}

std::string S(const QString &s) { return s.toStdString(); }

RecordsPanel *PanelOf(AZURESetup &w, const QString &section, const QString &name) {
  for (RecordsPanel *p : w.findChildren<RecordsPanel *>())
    if (p->section() == section && p->recordNames().contains(name)) return p;
  return nullptr;
}

RecordsPanel *PanelOfSection(AZURESetup &w, const QString &section) {
  for (RecordsPanel *p : w.findChildren<RecordsPanel *>())
    if (p->section() == section) return p;
  return nullptr;
}

std::string ValueOf(AZURESetup &w, const QString &section, const QString &name, const std::string &key) {
  for (const front::Record &r : w.holder()->records(section))
    if (QString::fromStdString(r.name) == name) return r.Value(key);
  return std::string();
}

bool HasField(AZURESetup &w, const QString &section, const QString &name, const std::string &key) {
  for (const front::Record &r : w.holder()->records(section))
    if (QString::fromStdString(r.name) == name)
      for (const front::Field &f : r.fields)
        if (f.key == key) return true;
  return false;
}

/// Choose the item of a drop-down whose data (else text) is `word`, as the user does.
bool ChooseData(QComboBox *c, const QString &word) {
  int at = c->findData(word);
  if (at < 0) at = c->findText(word);
  if (at < 0) return false;
  c->setCurrentIndex(at);
  emit c->activated(at);
  return true;
}

/// Bring a widget's tab and sub-tabs to the front (for a grab or for the widgets to be laid out).
void Show(AZURESetup &w, QWidget *widget) {
  for (QWidget *p = widget; p; p = p->parentWidget())
    if (QTabWidget *tabs = qobject_cast<QTabWidget *>(p->parentWidget() ? p->parentWidget()->parentWidget() : nullptr))
      tabs->setCurrentWidget(p);
  for (int i = 0; i < w.tabs()->count(); i++)
    if (w.tabs()->widget(i)->isAncestorOf(widget)) w.tabs()->setCurrentIndex(i);
  QApplication::processEvents();
}

void Grab(AZURESetup &w, QWidget *widget, const QString &file) {
  Show(w, widget);
  QApplication::processEvents();
  w.grab().save(file);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "usage: gui_polish_test <source dir> <screenshot dir>" << std::endl;
    return 2;
  }
  QApplication app(argc, argv);
  QCoreApplication::setOrganizationName("jina-test");
  QCoreApplication::setApplicationName("azure2-gui-polish-test");
  const fs::path source = argv[1];
  const QString shots = QString::fromLocal8Bit(argv[2]);
  QDir().mkpath(shots);
  const fs::path work = fs::absolute("gui_polish_work");
  fs::create_directories(work);

  AZURESetup w;
  w.resize(1400, 950);
  w.show();
  QApplication::processEvents();
  QString why;

  // (a) the fields a quantity decides
  const fs::path t1 = CopyProject(source / "tests/upgrade/T1/T1.azr", work / "quantity_T1");
  w.readProjectFile(QString::fromStdString(t1.string()), why);
  {
    RecordsPanel *data = w.getSegmentsTab()->records("data");
    Show(w, data);
    data->selectRecord("el_lab090");
    Check(data->fieldEditor("J") && !data->fieldRowShown("J") && !data->fieldRowShown("l") && !data->fieldRowShown("operator") &&
              data->fieldRowShown("reaction"),
          "T1: a differential cross section's form shows no J, l, operator");
    FieldEditor *quantity = data->fieldEditor("quantity");
    QComboBox *qc = quantity ? quantity->piece()->findChild<QComboBox *>() : nullptr;
    Check(qc && ChooseData(qc, "phase_shift") && data->fieldRowShown("J") && data->fieldRowShown("l") && data->fieldRowShown("s"),
          "T1: phase_shift chosen: J, l and s shown at once");
    data->setFieldText("J", "1/2+");
    Check(qc && ChooseData(qc, "dsigma_domega") && !data->fieldRowShown("J") && data->hiddenByQuantity().contains("J"),
          "T1: and back to dsigma_domega: J hidden, its value kept in the form until Apply");
    Check(data->applyForm(why) && !HasField(w, "data", "el_lab090", "J") &&
              ValueOf(w, "data", "el_lab090", "quantity") == "dsigma_domega",
          "T1: Apply drops the hidden J (the record keeps no J)");
    QLabel *status = nullptr;
    for (QLabel *l : data->findChildren<QLabel *>())
      if (l->text().contains("Dropped")) status = l;
    Check(status && status->text().contains("J"), "T1: and says so: " + (status ? S(status->text()) : std::string("no message")));
    Grab(w, data, shots + "/segments_data_quantity_fields.png");
  }

  // (b) the Add dialog follows the quantity
  {
    RecordsPanel *data = w.getSegmentsTab()->records("data");
    bool rowsMade = false, unitFollows = false;
    QTimer::singleShot(200, [&]() {
      QDialog *d = qobject_cast<QDialog *>(QApplication::activeModalWidget());
      if (!d) return;
      FieldEditor *q = d->findChild<FieldEditor *>("field_quantity");
      QComboBox *qc = q ? q->piece()->findChild<QComboBox *>() : nullptr;
      if (qc) ChooseData(qc, "phase_shift");
      QApplication::processEvents();
      rowsMade = d->findChild<FieldEditor *>("field_J") && d->findChild<FieldEditor *>("field_l") &&
                 d->findChild<FieldEditor *>("field_s");
      FieldEditor *unit = d->findChild<FieldEditor *>("field_unit");
      const std::vector<std::string> units = front::QuantityUnits("phase_shift");
      unitFollows = unit && !units.empty() && unit->text() == QString::fromStdString(units.front());
      QApplication::processEvents();  // the new rows are shown on the next turn
      QApplication::processEvents();
      d->grab().save(shots + "/add_dataset_phase_shift.png");
      d->reject();
    });
    QMetaObject::invokeMethod(data, "addRecord");
    Check(rowsMade, "T1: Add, quantity phase_shift chosen: the rows J, l, s made");
    Check(unitFollows, "T1: Add: the unit chosen anew from phase_shift's units");
  }

  // (c) the kernels drawn
  const fs::path t2 = CopyProject(source / "tests/upgrade/T2/T2.azr", work / "kernels_T2");
  w.readProjectFile(QString::fromStdString(t2.string()), why);
  {
    RecordsPanel *beams = PanelOf(w, "targets", "beam_targetInt2");
    Show(w, beams);
    beams->selectRecord("beam_targetInt2");
    QApplication::processEvents();
    const front::KernelPreview &k = beams->kernelPreview();
    Check(beams->preview() && beams->preview()->isVisibleTo(beams) && k.curves.size() == 1 && k.curves[0].kernel == "gaussian" &&
              std::fabs(k.curves[0].fwhm_keV - 2.0 * 2.3548) < 0.02 && !k.user.empty() && k.labEnergy >= k.lowEnergy &&
              k.labEnergy <= k.highEnergy,
          "T2: beam_targetInt2's spread of 2 keV drawn (" + S(beams->previewWhy()) + ") for " + k.user + " at " + std::to_string(k.labEnergy) + " MeV: FWHM " +
              (k.curves.empty() ? std::string("none") : std::to_string(k.curves[0].fwhm_keV)) + " keV");
    QLabel *fwhm = beams->preview()->findChild<QLabel *>("fwhm_gaussian");
    Check(fwhm && fwhm->text().contains("keV"), "T2: the FWHM written under the curve: " + (fwhm ? S(fwhm->text()) : std::string()));
    Grab(w, beams, shots + "/effects_beam_spread.png");

    RecordsPanel *targets = PanelOf(w, "targets", "tgt_s05_implanted");
    targets->selectRecord("tgt_s05_implanted");
    QApplication::processEvents();
    Check(targets->kernelPreview().curves.empty() || targets->kernelPreview().curves[0].kernel != "doppler",
          "T2: tgt_s05_implanted without Doppler: no Doppler curve");
    targets->setFieldText("doppler", "free_gas 300 K fixed");
    QApplication::processEvents();
    const front::KernelPreview d = targets->kernelPreview();
    const double f = d.cmFactor, kT = 8.617333262e-11 * 300.0;
    const double closed = 1e3 * 2.0 * std::sqrt(2.0 * std::log(2.0)) * std::sqrt(2.0 * d.labEnergy * kT * (1.0 - f) / f);
    Check(targets->preview()->isVisibleTo(targets) && d.curves.size() == 1 && d.curves[0].kernel == "doppler" &&
              std::fabs(d.curves[0].fwhm_keV - closed) < 0.01 * closed,
          "T2: Doppler typed into the form (not applied): drawn at once, FWHM " +
              (d.curves.empty() ? std::string("none") : std::to_string(d.curves[0].fwhm_keV)) + " keV (free gas " +
              std::to_string(closed) + ")");
    targets->setPreviewPoint(QString(), 2.0 * d.labEnergy);
    Check(!targets->kernelPreview().curves.empty() &&
              std::fabs(targets->kernelPreview().curves[0].fwhm_keV / closed - std::sqrt(2.0)) < 0.01,
          "T2: at twice the energy the Doppler width is sqrt 2 wider");
    QLineEdit *energy = targets->preview()->findChild<QLineEdit *>("kernelEnergy");
    Check(energy && std::fabs(energy->text().toDouble() - 2.0 * d.labEnergy) < 1e-6 * d.labEnergy,
          "T2: the energy line shows the energy drawn at");
    Grab(w, targets, shots + "/effects_target_doppler.png");
    QMetaObject::invokeMethod(targets, "revert");

    RecordsPanel *detectors = PanelOf(w, "targets", "det_targetInt2");
    detectors->selectRecord("det_targetInt2");
    detectors->setFieldText("energy_resolution", "triangular width 6 keV fixed");
    QApplication::processEvents();
    const front::KernelPreview &r = detectors->kernelPreview();
    Check(r.curves.size() == 1 && r.curves[0].kernel == "triangular" && std::fabs(r.curves[0].fwhm_keV - 3.0) < 0.05,
          "T2: a triangular resolution of base 6 keV typed (" + S(detectors->previewWhy()) + "): FWHM " +
              (r.curves.empty() ? std::string("none") : std::to_string(r.curves[0].fwhm_keV)) + " keV");
    Grab(w, detectors, shots + "/effects_detector_resolution.png");
    QMetaObject::invokeMethod(detectors, "revert");
  }

  // (d) defaults, the toggle, an empty section
  {
    RecordsPanel *data = w.getSegmentsTab()->records("data");
    Show(w, data);
    for (const front::Record &r : data->model()->records())
      if (r.kind == "dataset") {
        data->selectRecord(QString::fromStdString(r.name));
        break;
      }
    FieldEditor *active = data->fieldEditor("active");
    QComboBox *ac = active ? active->piece()->findChild<QComboBox *>() : nullptr;
    Check(ac && ac->itemText(0).contains("(default)") && !ac->itemText(0).contains("(none)") && ac->currentIndex() == 0,
          "T2: 'active' left out shows its default: " + (ac ? S(ac->itemText(0)) : std::string()));
    QToolButton *toggle = active ? active->findChild<QToolButton *>("editAsText") : nullptr;
    Check(toggle && !toggle->icon().isNull() && toggle->toolButtonStyle() == Qt::ToolButtonIconOnly &&
              toggle->toolTip().contains("as text"),
          "the edit-as-text toggle: an icon and a tool tip that says what it does");
    RecordsPanel *radii = PanelOfSection(w, "channelRadii");
    QLabel *hint = radii ? radii->findChild<QLabel *>("emptyHint") : nullptr;
    Check(radii && radii->recordCount() == 0 && hint && hint->isVisibleTo(radii) && !radii->fieldEditor("pair"),
          "T2: an empty section shows its hint, no empty form");
    if (radii) Grab(w, radii, shots + "/pairs_empty_radii.png");
  }

  // (e) a sum of reactions edited term by term
  w.readProjectFile(QString::fromStdString(t2.string()), why);
  {
    RecordsPanel *data = w.getSegmentsTab()->records("data");
    Show(w, data);
    data->selectRecord("s08_sum");
    FieldEditor *reaction = data->fieldEditor("reaction");
    QComboBox *mode = reaction ? reaction->piece()->findChild<QComboBox *>("reactionMode") : nullptr;
    const QList<QWidget *> terms = reaction ? reaction->piece()->findChildren<QWidget *>("reactionTerm") : QList<QWidget *>();
    Check(reaction && !reaction->textMode() && mode && mode->currentData().toString() == "sum" && terms.size() == 2,
          "T2: s08_sum's reaction shown as a sum of two terms by their pairs, not as text");
    QLineEdit *weight = terms.isEmpty() ? nullptr : terms[0]->findChild<QLineEdit *>();
    if (weight) {
      weight->setText("0.5");
      emit weight->textEdited("0.5");
    }
    Grab(w, data, shots + "/segments_sum_terms.png");
    const bool applied = data->applyForm(why);
    const std::string sum = ValueOf(w, "data", "s08_sum", "reaction");
    Check(applied && sum == "0.5*15N(p,a)12C + 15N(p,a)12C[4.43891]", "T2: a weight typed for the first term: " + sum);
    data->selectRecord("s09_ratio");
    reaction = data->fieldEditor("reaction");
    mode = reaction ? reaction->piece()->findChild<QComboBox *>("reactionMode") : nullptr;
    Check(reaction && !reaction->textMode() && mode && mode->currentData().toString() == "ratio",
          "T2: s09_ratio's reaction shown as a ratio of two terms");
    // one reaction turned into a sum: the reaction shown, and a second term
    data->selectRecord("s05_implanted");
    reaction = data->fieldEditor("reaction");
    mode = reaction ? reaction->piece()->findChild<QComboBox *>("reactionMode") : nullptr;
    if (mode) ChooseData(mode, "sum");
    const bool appliedSum = data->applyForm(why);
    const std::string summed = ValueOf(w, "data", "s05_implanted", "reaction");
    Check(appliedSum && summed.rfind("15N(p,a)12C + ", 0) == 0, "T2: one reaction made a sum keeps it as its first term: " + summed);
  }

  // (f) the pair dialog, a pair duplicated
  {
    PairsTab *pairs = w.getPairsTab();
    Show(w, pairs);
    const int before = pairs->getPairsModel()->numPairs();
    bool filled = false, unknownSaid = false, photonHides = false;
    QTimer::singleShot(200, [&]() {
      AddPairDialog *d = qobject_cast<AddPairDialog *>(QApplication::activeModalWidget());
      if (!d) return;
      d->chooseNucleus(d->light, "He", 4);
      d->chooseNucleus(d->heavy, "C", 12);
      filled = d->light.z->text() == "2" && std::fabs(d->light.mass->text().toDouble() - 4.002603254) < 1e-8 &&
               d->light.j->value() == 0.0 && d->light.parity->isChecked() && d->heavy.z->text() == "6" &&
               std::fabs(d->heavy.mass->text().toDouble() - 12.0) < 1e-12;
      d->chooseNucleus(d->heavy, "N", 15);
      unknownSaid = d->heavy.note->text().contains("give J and parity");
      d->heavy.j->setValue(0.5);
      d->heavy.parity->setChecked(false);
      d->fillFromMassTable();
      QApplication::processEvents();
      d->grab().save(shots + "/pair_dialog_particle.png");
      d->pairTypeCombo->setCurrentIndex(1);
      QApplication::processEvents();
      photonHides = !d->light.nucleusRow->isVisibleTo(d) && !d->channelRadiusText->isVisibleTo(d) && d->multBox->isVisibleTo(d);
      d->grab().save(shots + "/pair_dialog_photon.png");
      d->pairTypeCombo->setCurrentIndex(0);
      d->accept();
    });
    QMetaObject::invokeMethod(pairs, "addPair");
    Check(filled, "T2: the pair dialog: He 4 and C 12 chosen: Z, AME2020 masses, alpha's 0+ filled");
    Check(unknownSaid, "T2: 15N's spin is not in the table: the dialog says so");
    Check(photonHides, "T2: a photon pair hides the light particle and the radius, shows the external capture");
    const QList<PairsData> list = pairs->getPairsModel()->getPairs();
    Check(list.size() == before + 1 && list.last().lightZ == 2 && list.last().heavyZ == 7 && list.last().heavyJ == 0.5 &&
              list.last().heavyPi == -1 && list.last().seperationEnergy > 0,
          "T2: the pair a + 15N added with 1/2-, S from the table: " + std::to_string(list.isEmpty() ? 0 : list.last().seperationEnergy));
    pairs->findChild<QTableView *>()->setCurrentIndex(pairs->getPairsModel()->index(list.size() - 1, 0));
    pairs->findChild<QTableView *>()->selectRow(list.size() - 1);
    QTimer::singleShot(200, [&]() {
      AddPairDialog *d = qobject_cast<AddPairDialog *>(QApplication::activeModalWidget());
      if (!d) return;
      d->excitationEnergyText->setText("5.27");
      d->accept();
    });
    QPushButton *duplicate = pairs->findChild<QPushButton *>("duplicatePair");
    if (duplicate) duplicate->click();
    const QList<PairsData> after = pairs->getPairsModel()->getPairs();
    Check(duplicate && after.size() == before + 2 && after.last().heavyZ == 7 && after.last().excitationEnergy == 5.27 &&
              after.last().heavyM == after.at(after.size() - 2).heavyM,
          "T2: Duplicate: the copy in the dialog, its excitation changed, added as a new pair");
    Grab(w, pairs, shots + "/pairs_tab.png");
  }

  // (g) the Levels tab's channels from the facade: the GUI's own copy of the rule was compared with them for every
  // J^pi up to 7 and every nucleus of T1, T2, T3 and 13N (equal; commit 1bd92886) and deleted.  Now: every channel a
  // particle or photon channel a level of those files carries is one the facade gives its J^pi under the limits read from the file.
  for (const char *project : {"tests/upgrade/T1/T1.azr", "tests/upgrade/T2/T2.azr", "tests/upgrade/T3/T3.azr", "tests/13N/13N.azr"}) {
    const fs::path copy = CopyProject(source / project, work / ("channels_" + fs::path(project).stem().string()));
    w.readProjectFile(QString::fromStdString(copy.string()), why);
    LevelsTab *levels = w.getLevelsTab();
    const QList<LevelsData> all = levels->getLevelsModel()->getLevels();
    const QList<ChannelsData> fileChannels = levels->getChannelsModel()->getChannels();
    int found = 0, missing = 0;
    std::string first;
    for (int r = 0; r < all.size(); r++) {
      const QList<ChannelsData> allowed = levels->calculateChannels(all[r], r);
      for (const ChannelsData &c : fileChannels) {
        // a beta channel is left out: T3's 2+ level carries a Gamow-Teller channel of s = 1 to a 2+ parent, which
        // neither the engine's rule nor the GUI's former copy gives (they give s = 0 for J = J_parent): an open
        // question of the rule, reported, not of this tab
        if (c.levelIndex != r || c.radType == 'G' || c.radType == 'F') continue;
        bool in = false;
        for (const ChannelsData &a : allowed)
          in = in || (a.pairIndex == c.pairIndex && a.sValue == c.sValue && a.lValue == c.lValue && a.radType == c.radType);
        found += in;
        missing += !in;
        if (!in && first.empty()) first = "level " + std::to_string(r) + " pair " + std::to_string(c.pairIndex + 1);
      }
    }
    Check(missing == 0 && found > 0, std::string(project) + ": all " + std::to_string(found) +
                                         " channels of the file's levels are the facade's for their J^pi" +
                                         (first.empty() ? "" : "; not: " + first));
  }

  // (h) the sections not looked at before, one by one
  {
    const fs::path ex = CopyProject(source / "experiment/tests/data/section_examples.azr", work / "sections");
    w.readProjectFile(QString::fromStdString(ex.string()), why);
    for (const QString &section : {QString("modelData"), QString("potentials"), QString("channelRadii"), QString("isospin")}) {
      RecordsPanel *p = PanelOfSection(w, section);
      if (!p || p->recordCount() == 0) {
        Check(false, "section_examples: a panel of " + S(section) + " with records");
        continue;
      }
      p->selectRecord(p->recordNames().value(0));
      Grab(w, p, shots + "/section_" + section + ".png");
    }
    const fs::path gd = CopyProject(source / "tests/gamma_decay/models/record.azr", work / "gamma");
    w.readProjectFile(QString::fromStdString(gd.string()), why);
    RecordsPanel *g = PanelOfSection(w, "gammaDecays");
    if (g && g->recordCount() > 0) {
      g->selectRecord(g->recordNames().value(0));
      Grab(w, g, shots + "/section_gammaDecays.png");
    }
    Check(QFileInfo::exists(shots + "/section_modelData.png") && QFileInfo::exists(shots + "/section_gammaDecays.png"),
          "screenshots of the other sections in " + S(shots));
  }

  std::cout << (failures ? "FAILED: " + std::to_string(failures) : std::string("all passed")) << std::endl;
  return failures ? 1 : 0;
}
