// The forms of the record tabs without a display (docs/upgrade/qt_gui.md section 11):
//   gui_forms_test <source dir> <screenshot dir>
//  (a) every record of T1, T2, T3, 13N and section_examples.azr shown in its form: every field by its widgets (none
//      falls back to text), the widgets' value written from its pieces equal to the field as written, and Apply of the
//      unedited record moving no byte;
//  (b) values set through the widgets read back from the project text: a nuisance's "Fit" and its prior, a beam
//      profile's omega, the stopping's alternative, a reaction by its exit pair, a unit of the quantity;
//  (c) the table: its columns from the schema's summary fields, a norm and its "fit" edited in place, a data set left
//      out and used again, a record moved down and up (the same bytes), several selected;
//  (d) the engine's word beside the field as the user types, and Apply refusing a record the engine refuses;
//  (e) "used by" and "assign to" of a target, the Add dialog pre-filled from the project;
//  (f) screenshots of the Segments and Experimental Effects tabs into the screenshot directory.
// PASS/FAIL per check; exit 1 on a failure.
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
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QScrollArea>
#include <QTabWidget>
#include <QTableView>
#include <QTimer>

#include "AZURESetup.h"
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
  return to / projectFile.filename();
}

std::string S(const QString &s) { return s.toStdString(); }

/// The panel of the window that shows `name` of `section`.
RecordsPanel *PanelOf(AZURESetup &w, const QString &section, const QString &name) {
  for (RecordsPanel *p : w.findChildren<RecordsPanel *>())
    if (p->section() == section && p->recordNames().contains(name)) return p;
  return nullptr;
}

/// Show a panel's tab and sub-tab (for a grab or for the widgets to be laid out).
void ShowPanel(AZURESetup &w, QWidget *panel) {
  for (QWidget *p = panel; p; p = p->parentWidget())
    if (QTabWidget *tabs = qobject_cast<QTabWidget *>(p->parentWidget() ? p->parentWidget()->parentWidget() : nullptr))
      tabs->setCurrentWidget(p);
  QApplication::processEvents();
  (void)w;
}

/// The value of a field of a record of the project.
std::string ValueOf(AZURESetup &w, const QString &section, const QString &name, const std::string &key) {
  for (const front::Record &r : w.holder()->records(section))
    if (QString::fromStdString(r.name) == name) return r.Value(key);
  return std::string();
}

int DifferingLines(const std::string &a, const std::string &b) {
  std::istringstream x(a), y(b);
  std::string la, lb;
  int n = 0;
  while (true) {
    const bool ga = (bool)std::getline(x, la), gb = (bool)std::getline(y, lb);
    if (!ga && !gb) break;
    if (!ga || !gb || la != lb) n++;
  }
  return n;
}

/// Type into a line edit as the user does.
void Type(QLineEdit *e, const QString &text) {
  e->setText(text);
  emit e->textEdited(text);
}

QLineEdit *LineWithText(QWidget *in, const QString &text) {
  for (QLineEdit *e : in->findChildren<QLineEdit *>())
    if (e->text() == text && e->isVisibleTo(in)) return e;
  return nullptr;
}

QCheckBox *CheckBoxNamed(QWidget *in, const QString &text) {
  for (QCheckBox *b : in->findChildren<QCheckBox *>())
    if (b->text() == text) return b;
  return nullptr;
}

/// Choose the item of a drop-down whose text contains `part`, as the user does.
bool Choose(QComboBox *c, const QString &part) {
  for (int i = 0; i < c->count(); i++)
    if (c->itemText(i).contains(part)) {
      c->setCurrentIndex(i);
      emit c->activated(i);
      return true;
    }
  return false;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "usage: gui_forms_test <source dir> <screenshot dir>" << std::endl;
    return 2;
  }
  QApplication app(argc, argv);
  QCoreApplication::setOrganizationName("jina-test");
  QCoreApplication::setApplicationName("azure2-gui-forms-test");
  const fs::path source = argv[1];
  const QString shots = QString::fromLocal8Bit(argv[2]);
  QDir().mkpath(shots);
  const fs::path work = fs::absolute("gui_forms_work");
  fs::create_directories(work);

  AZURESetup w;
  w.resize(1400, 950);
  w.show();
  QApplication::processEvents();

  // (a) every record through its form
  const std::vector<fs::path> projects = {source / "tests/upgrade/T1/T1.azr", source / "tests/upgrade/T2/T2.azr",
                                          source / "tests/upgrade/T3/T3.azr", source / "tests/13N/13N.azr",
                                          source / "experiment/tests/data/section_examples.azr"};
  for (const fs::path &f : projects) {
    const std::string name = f.stem().string();
    const fs::path copy = CopyProject(f, work / ("rt_" + name));
    QString why;
    if (!w.readProjectFile(QString::fromStdString(copy.string()), why)) {
      Check(false, name + ": opened (" + S(why) + ")");
      continue;
    }
    const std::string original = Bytes(f);
    int records = 0, fields = 0, asText = 0, same = 0, applied = 0;
    std::vector<std::string> notes;
    for (const auto &section : front::RecordSections()) {
      const QString sec = QString::fromStdString(section.first);
      for (const front::Record &r : w.holder()->records(sec)) {
        const QString rn = QString::fromStdString(r.name);
        RecordsPanel *panel = PanelOf(w, sec, rn);
        if (!panel || !panel->selectRecord(rn)) {
          notes.push_back(r.name + ": no panel shows it");
          continue;
        }
        records++;
        int occurrence = 0;
        std::string last;
        for (const front::Field &x : r.fields) {
          occurrence = x.key == last ? occurrence + 1 : 0;
          last = x.key;
          FieldEditor *e = panel->fieldEditor(QString::fromStdString(x.key), occurrence);
          if (!e) {
            notes.push_back(r.name + "." + x.key + ": no editor");
            continue;
          }
          fields++;
          const bool reaction = e->field().value.type == "reaction_expr";
          int in = 0, out = 0;
          const bool pairs = reaction && front::ReactionPairs(w.holder()->project(), x.value, in, out) && in > 0 && out != 0;
          if (e->textMode() && reaction && !pairs) {
            same++;  // a sum, a ratio, a beta decay or a project without the pairs: a reaction the pair drop-downs
                     // cannot show stays as text
            continue;
          }
          if (e->textMode()) {
            asText++;
            notes.push_back(r.name + "." + x.key + " shown as text: " + x.value);
            continue;
          }
          const QString fromWidgets = QString::fromStdString(front::ValueText(e->field().value, e->piece()->get()));
          int win = 0, wout = 0;
          const bool samePairs = reaction && front::ReactionPairs(w.holder()->project(), S(fromWidgets), win, wout) &&
                                 win == in && wout == out;  // the label may spell an excitation otherwise ("[1]")
          if ((fromWidgets == OneBlank(QString::fromStdString(x.value)) || samePairs) && e->text() == QString::fromStdString(x.value))
            same++;
          else
            notes.push_back(r.name + "." + x.key + ": '" + x.value + "' -> widgets '" + S(fromWidgets) + "'");
        }
        applied += panel->applyForm(why);
      }
    }
    for (std::size_t i = 0; i < notes.size() && i < 10; i++) std::cout << "    " << notes[i] << std::endl;
    Check(records > 0 && asText == 0 && same == fields,
          name + ": " + std::to_string(records) + " records, " + std::to_string(same) + " of " + std::to_string(fields) +
              " fields shown by their widgets and written back from them unchanged (" + std::to_string(asText) + " as text)");
    Check(applied == records && w.holder()->project().Text() == original,
          name + ": " + std::to_string(applied) + " records applied unedited: not a byte moved");
  }

  // (b) values set through the widgets
  const fs::path t2 = CopyProject(source / "tests/upgrade/T2/T2.azr", work / "edit_T2");
  QString why;
  w.readProjectFile(QString::fromStdString(t2.string()), why);
  {
    const std::string before = w.holder()->project().Text();
    RecordsPanel *targets = PanelOf(w, "targets", "tgt_s05_implanted");
    Check(targets && targets->selectRecord("tgt_s05_implanted") && targets->kinds() == QStringList{"target"},
          "T2: the targets have a list of their own");
    FieldEditor *areal = targets ? targets->fieldEditor("areal_density") : nullptr;
    QCheckBox *fit = areal ? CheckBoxNamed(areal->piece(), QObject::tr("Fit")) : nullptr;
    QCheckBox *prior = areal ? CheckBoxNamed(areal->piece(), QObject::tr("Gaussian prior: centre")) : nullptr;
    Check(fit && fit->isChecked() && prior && prior->isChecked(), "T2: the areal density's Fit and prior shown checked");
    if (fit) fit->setChecked(false);
    if (prior) prior->setChecked(false);
    Check(targets && targets->applyForm(why) &&
              ValueOf(w, "targets", "tgt_s05_implanted", "areal_density") == "2e+17 atoms/cm2 fixed grid_max 3e+17" &&
              DifferingLines(before, w.holder()->project().Text()) == 1,
          "T2: Fit and the prior unchecked: '" + ValueOf(w, "targets", "tgt_s05_implanted", "areal_density") +
              "', one line changed");
    // the stopping's alternative and its formula's parameter
    targets->selectRecord("tgt_s05_implanted");
    FieldEditor *stopping = targets->fieldEditor("stopping");
    QLineEdit *a1 = stopping ? LineWithText(stopping->piece(), "-0.30") : nullptr;
    if (a1) Type(a1, "-0.35");
    targets->applyForm(why);
    Check(ValueOf(w, "targets", "tgt_s05_implanted", "stopping") == "formula \"a0*x^a1\" a0=6.0e-20 a1=-0.35 MeV*cm2/atom",
          "T2: a formula's parameter typed in its own field: " + ValueOf(w, "targets", "tgt_s05_implanted", "stopping"));
    targets->selectRecord("tgt_s05_implanted");
    stopping = targets->fieldEditor("stopping");
    QComboBox *alternative = stopping ? stopping->piece()->findChild<QComboBox *>() : nullptr;
    Check(alternative && alternative->count() >= 3 && alternative->currentText() == "formula",
          "T2: the stopping's alternatives in a drop-down (" + std::to_string(alternative ? alternative->count() : 0) + ")");
    // the beam profile's omega
    RecordsPanel *beams = PanelOf(w, "targets", "beam_s07_profile");
    Check(beams && beams != targets && beams->selectRecord("beam_s07_profile"), "T2: the beams have a list of their own");
    FieldEditor *profile = beams ? beams->fieldEditor("profile") : nullptr;
    QLineEdit *omega = profile ? LineWithText(profile->piece(), "0.035") : nullptr;
    if (omega) Type(omega, "0.04");
    if (beams) beams->applyForm(why);
    Check(ValueOf(w, "targets", "beam_s07_profile", "profile") == "skewed_gaussian xi 0.34 MeV omega 0.04 MeV alpha 1.5 weight 1.0",
          "T2: the beam profile's omega typed in its own field: " + ValueOf(w, "targets", "beam_s07_profile", "profile"));
    beams->selectRecord("beam_s07_profile");
    Check(beams->preview() && beams->preview()->isVisibleTo(beams), "T2: the beam profile drawn beside its form");
    QWidget *profileBox = beams ? beams->box("profile") : nullptr;
    Check(profileBox && profileBox->findChild<FieldEditor *>("field_truncation"),
          "T2: the profile's truncation in the profile's box");
  }
  {
    // a data set: its unit from the quantity's units, its reaction by the exit pair
    RecordsPanel *data = w.getSegmentsTab()->records("data");
    QString first;
    for (const front::Record &r : data->model()->records())
      if (r.kind == "dataset" && first.isEmpty() && r.Value("reaction").find('+') == std::string::npos &&
          r.Value("reaction").find('/') == std::string::npos)
        first = QString::fromStdString(r.name);
    data->selectRecord(first);
    FieldEditor *unit = data->fieldEditor("unit");
    QComboBox *units = unit ? unit->piece()->findChild<QComboBox *>() : nullptr;
    const std::vector<std::string> allowed = front::QuantityUnits(ValueOf(w, "data", first, "quantity"));
    Check(units && !allowed.empty() && units->findData(QString::fromStdString(allowed.back())) >= 0 && units->findData("MeV") < 0,
          "T2: the unit drop-down offers the units of " + ValueOf(w, "data", first, "quantity"));
    FieldEditor *reaction = data->fieldEditor("reaction");
    Check(reaction && !reaction->textMode(), "T2: the reaction shown by its pair drop-downs");
    const std::vector<front::PairChoice> pairs = front::Pairs(w.holder()->project());
    QList<QComboBox *> combos = reaction ? reaction->piece()->findChildren<QComboBox *>() : QList<QComboBox *>();
    std::string other;
    if (combos.size() == 2)
      for (int i = 0; i < combos[1]->count(); i++)
        if (combos[1]->itemData(i).toInt() > 0 && combos[1]->itemData(i).toInt() != combos[1]->currentData().toInt()) {
          combos[1]->setCurrentIndex(i);
          emit combos[1]->activated(i);
          other = front::ReactionLabel(w.holder()->project(), combos[0]->currentData().toInt(), combos[1]->currentData().toInt());
          break;
        }
    data->applyForm(why);
    Check(!other.empty() && ValueOf(w, "data", first, "reaction") == other,
          "T2: another exit pair chosen: the reaction " + ValueOf(w, "data", first, "reaction"));
  }

  // (c) the table
  const fs::path t1 = CopyProject(source / "tests/upgrade/T1/T1.azr", work / "table_T1");
  w.readProjectFile(QString::fromStdString(t1.string()), why);
  {
    const std::string original = w.holder()->project().Text();
    RecordsPanel *data = w.getSegmentsTab()->records("data");
    RecordTableModel *m = data->model();
    QStringList titles;
    for (int c = 0; c < m->columnCount(); c++) titles << m->headerData(c, Qt::Horizontal).toString();
    const int norm = m->columnOf(RecordTableModel::Column::Field, "norm"), fit = m->columnOf(RecordTableModel::Column::Fit, "norm");
    const int use = m->columnOf(RecordTableModel::Column::Include);
    Check(use == 0 && m->columnOf(RecordTableModel::Column::Field, "reaction") >= 0 && norm > 0 && fit == norm + 1 &&
              m->columnOf(RecordTableModel::Column::Field, "columns") < 0,
          "T1: the data sets' columns: " + S(titles.join(" | ")));
    Check(m->data(m->index(0, norm)).toString() == "1" && m->data(m->index(0, fit), Qt::CheckStateRole).toInt() == Qt::Unchecked &&
              m->data(m->index(0, 1), Qt::ToolTipRole).toString().contains("columns"),
          "T1: a row's norm, its fit box and the record as its tooltip");
    const QString first = QString::fromStdString(m->records()[0].name);
    m->setData(m->index(0, norm), "1.05", Qt::EditRole);
    m->setData(m->index(0, fit), Qt::Checked, Qt::CheckStateRole);
    Check(ValueOf(w, "data", first, "norm") == "1.05 fit" && DifferingLines(original, w.holder()->project().Text()) == 1,
          "T1: the norm edited in its cell and its fit box checked: '" + ValueOf(w, "data", first, "norm") + "'");
    m->setData(m->index(0, norm), "1", Qt::EditRole);
    m->setData(m->index(0, fit), Qt::Unchecked, Qt::CheckStateRole);
    Check(w.holder()->project().Text() == original, "T1: and back: the same bytes");
    const int quantity = m->columnOf(RecordTableModel::Column::Field, "quantity");
    Check(quantity >= 0 && m->choices(m->index(0, quantity)).contains("sigma"), "T1: the quantity cell offers the quantities");
    Check(!m->setData(m->index(0, norm), "abc", Qt::EditRole) && w.holder()->project().Text() == original,
          "T1: a norm the parser refuses is not written");
    // left out and used again
    data->selectRecords({first});
    data->includeSelected(false, why);
    Check(ValueOf(w, "data", first, "active") == "no" && m->data(m->index(0, use), Qt::CheckStateRole).toInt() == Qt::Unchecked,
          "T1: the data set left out: active no");
    data->selectRecords({first});
    data->includeSelected(true, why);
    Check(ValueOf(w, "data", first, "active").empty() && w.holder()->project().Text() == original,
          "T1: and used again: the field gone, the same bytes");
    // moved down and up
    const QString second = data->recordNames().value(1);
    data->selectRecord(first);
    Check(data->moveShown(+1, why) && data->recordNames().value(1) == first && data->recordNames().value(0) == second,
          "T1: a data set moved down");
    Check(data->moveShown(-1, why) && w.holder()->project().Text() == original, "T1: and up again: the same bytes");
    data->selectRecords({first, second});
    Check(data->selectedNames().size() == 2, "T1: two data sets selected");
    data->setFilter(second);
    Check(data->table()->model()->rowCount() == 1, "T1: the filter line keeps one row");
    data->setFilter(QString());
  }

  // (d) the engine's words as the user types; Apply refuses what the engine refuses
  {
    const std::string original = w.holder()->project().Text();
    RecordsPanel *data = w.getSegmentsTab()->records("data");
    const QString first = data->recordNames().value(0);
    data->selectRecord(first);
    data->setFieldText("norm", "1 fitt");
    FieldEditor *norm = data->fieldEditor("norm");
    Check(norm && !norm->message().isEmpty(), "T1: a mistyped norm: '" + (norm ? S(norm->message()) : std::string()) + "'");
    Check(!data->applyForm(why, true) && why.contains("norm") && w.holder()->project().Text() == original,
          "T1: Apply refuses it, naming the field: " + S(why));
    data->selectRecord(first);
    data->setFieldText("unit", "barn_per_furlong");
    Check(data->fieldEditor("unit") && !data->fieldEditor("unit")->message().isEmpty(),
          "T1: a wrong unit: " + S(data->fieldEditor("unit")->message()));
    data->selectRecord(first);  // the form keeps the record shown while the engine refuses it
    Check(data->shownName() == first, "T1: the record with a refused edit stays shown");
    QMetaObject::invokeMethod(data, "revert");
    Check(w.holder()->project().Text() == original, "T1: reverted: nothing written");
  }

  // (e) "used by", "assign to", the Add dialog
  w.readProjectFile(QString::fromStdString(t2.string()), why);
  {
    RecordsPanel *targets = PanelOf(w, "targets", "tgt_s05_implanted");
    targets->selectRecord("tgt_s05_implanted");
    const QStringList users = targets->usedBy();
    Check(!users.isEmpty(), "T2: tgt_s05_implanted used by " + S(users.join(", ")));
    QStringList names;
    for (const QString &u : users) names << u.section('/', 1);
    QString more;
    for (const front::Record &r : w.getSegmentsTab()->records("data")->model()->records())
      if (r.kind == "dataset" && r.Value("target").empty() && more.isEmpty()) more = QString::fromStdString(r.name);
    names << more;
    Check(targets->assignShownTo(names, why) && ValueOf(w, "data", more, "target") == "tgt_s05_implanted",
          "T2: assigned to " + S(more) + " too");
    names.removeAll(more);
    Check(targets->assignShownTo(names, why) && targets->usedBy() == users, "T2: and taken back");

    RecordsPanel *data = w.getSegmentsTab()->records("data");
    const int before = data->recordCount();
    QTimer::singleShot(200, [&]() {
      if (QDialog *d = qobject_cast<QDialog *>(QApplication::activeModalWidget())) {
        d->grab().save(shots + "/add_dataset_dialog.png");
        d->accept();
      }
    });
    QMetaObject::invokeMethod(data, "addRecord");
    const QString added = data->shownName();
    Check(data->recordCount() == before + 1 && ValueOf(w, "data", added, "reaction") != "?" &&
              ValueOf(w, "data", added, "quantity") != "?" && ValueOf(w, "data", added, "unit") != "?" &&
              ValueOf(w, "data", added, "norm") == "1 fixed" && ValueOf(w, "data", added, "columns") != "?" &&
              ValueOf(w, "data", added, "quadrature").empty(),
          "T2: Add pre-fills a data set: " + ValueOf(w, "data", added, "reaction") + ", " + ValueOf(w, "data", added, "quantity") +
              ", " + ValueOf(w, "data", added, "unit") + ", norm " + ValueOf(w, "data", added, "norm") +
              ", columns " + ValueOf(w, "data", added, "columns"));
    data->removeShown(why);
  }

  // (f) screenshots
  {
    w.readProjectFile(QString::fromStdString(t2.string()), why);
    auto grab = [&](RecordsPanel *panel, const QString &record, const QString &file, const QStringList &open) {
      if (!panel) return;
      panel->selectRecord(record);
      for (const QString &key : open)
        if (QWidget *b = panel->box(key))
          for (QAbstractButton *h : b->findChildren<QAbstractButton *>())
            if (h->isCheckable() && h->parentWidget() == b) h->setChecked(true);
      w.tabs()->setCurrentWidget(panel->window() == &w ? w.tabs()->widget(0) : nullptr);
      for (int i = 0; i < w.tabs()->count(); i++)
        if (w.tabs()->widget(i)->isAncestorOf(panel)) w.tabs()->setCurrentIndex(i);
      ShowPanel(w, panel);
      w.grab().save(shots + "/" + file);
    };
    grab(w.getSegmentsTab()->records("data"), "s02_diff090", "segments_data.png", {});
    grab(w.getSegmentsTab()->records("extrapolations"), w.getSegmentsTab()->records("extrapolations")->recordNames().value(0),
         "segments_extrapolations.png", {});
    grab(PanelOf(w, "targets", "tgt_s05_implanted"), "tgt_s05_implanted", "effects_target.png", {});
    grab(PanelOf(w, "targets", "beam_s07_profile"), "beam_s07_profile", "effects_beam_profile.png", {});
    grab(PanelOf(w, "targets", "det_s07_profile"), "det_s07_profile", "effects_detector.png", {});
    Check(QFileInfo::exists(shots + "/effects_target.png") && QFileInfo::exists(shots + "/segments_data.png"),
          "screenshots in " + S(shots));
  }

  std::cout << (failures ? "FAILED: " + std::to_string(failures) : std::string("all passed")) << std::endl;
  return failures ? 1 : 0;
}
