// The Physics tab without a display (docs/upgrade/qt_gui.md section 15):
//   gui_physqt_test <source dir> <screenshot dir>
// Every page driven through its widgets, the plotted arrays compared with a direct call of the facade:
//  (a) no project: Functions (a row drawn, an argument changed and a curve added, points outside the domain left as
//      gaps and counted in the status line), Stopping (a material's stopping against EvaluateStopping, the citation),
//      the project pages saying what they need;
//  (b) T1: Channels (P of pair 1, l = 1, against ChannelFunctions), R-matrix (the 1/2+ group's phase shift against
//      Block, the pole table against Poles(), a click on a pole centres the plot on it, the Argand plot's arrays),
//      Pole map (the map's rows against Block's log10_abs_det of every group, the poles on it, the cursor read-out of a
//      cell, a click on a row gives its cut);
//  (c) T2: THM (LC09's vertex against ThmVertex);
//  (d) the tab in the window; screenshots of every page in the light and the dark look.
// PASS/FAIL per check; exit 1 on a failure.
#include <clocale>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>

#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDir>
#include <QElapsedTimer>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMouseEvent>
#include <QPushButton>
#include <QSettings>
#include <QSpinBox>
#include <QStackedWidget>
#include <QTabWidget>
#include <QTableWidget>

#include "AZURESetup.h"
#include "AzureTheme.h"
#include "HeatmapWidget.h"
#include "PhysicsTab.h"
#include "ProjectHolder.h"

namespace fs = std::filesystem;
namespace front = azure2::front;
using namespace physicsqt;

namespace {

int failures = 0;

void Check(bool ok, const std::string &what) {
  std::cout << (ok ? "PASS " : "FAIL ") << what << std::endl;
  if (!ok) failures++;
}

void Settle() {
  for (int i = 0; i < 5; i++) QApplication::processEvents();
}

/// Turn the event loop until the tab's workers are idle (and a little after, for the queued results).
bool WaitIdle(PhysicsTab &tab, int ms = 180000) {
  QElapsedTimer t;
  t.start();
  do {
    QApplication::processEvents(QEventLoop::AllEvents, 30);
    if (!tab.busy()) {
      Settle();
      if (!tab.busy()) return true;
    }
  } while (t.elapsed() < ms);
  return false;
}

fs::path CopyProject(const fs::path &projectFile, const fs::path &to) {
  fs::remove_all(to);
  fs::create_directories(to);
  fs::copy(projectFile.parent_path(), to, fs::copy_options::recursive | fs::copy_options::skip_symlinks);
  for (const char *d : {"output", "checks"}) {
    fs::remove_all(to / d);
    fs::create_directories(to / d);
  }
  return to / projectFile.filename();
}

/// The largest relative difference of two arrays (NaN where both are NaN counts as equal; a NaN on one side as 1).
double Differ(const QVector<double> &a, const std::vector<double> &b) {
  if ((std::size_t)a.size() != b.size()) return INFINITY;
  double d = 0;
  for (int i = 0; i < a.size(); i++) {
    if (std::isnan(a[i]) && std::isnan(b[i])) continue;
    if (std::isnan(a[i]) != std::isnan(b[i])) return 1;
    d = std::max(d, std::fabs(a[i] - b[i]) / std::max(1e-300, std::fabs(b[i])));
  }
  return d;
}

int SelectData(QComboBox *box, const QVariant &data) {
  const int i = box->findData(data);
  if (i >= 0) box->setCurrentIndex(i);
  return i;
}

void Shot(QWidget &w, const QString &file) { w.grab().save(file); }

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "usage: gui_physqt_test <source dir> <screenshot dir>" << std::endl;
    return 2;
  }
  QApplication app(argc, argv);
  // as the program does after making its application: this test builds the Physics pages without the main window
  std::setlocale(LC_NUMERIC, "C");
  QCoreApplication::setAttribute(Qt::AA_DontUseNativeDialogs, true);
  const fs::path source = argv[1];
  const QString shots = QString::fromLocal8Bit(argv[2]);
  const fs::path work = fs::absolute("gui_physqt_work");
  fs::remove_all(work);
  fs::create_directories(work);
  QSettings::setDefaultFormat(QSettings::IniFormat);
  QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, QString::fromStdString((work / "settings").string()));
  QCoreApplication::setOrganizationName("jina-physqt-test");
  QCoreApplication::setApplicationName("azure2-physqt-test");
  QDir().mkpath(shots);
  qputenv("AZURE2_THEME", "");
  applyAzureTheme(app);
  setAzureTheme(ThemeChoice::Light);

  ProjectHolder holder;
  PhysicsTab tab(&holder);
  tab.resize(1350, 860);
  tab.show();
  Settle();

  // ---------------------------------------------------------------- (a) no project
  FunctionsPage *fp = tab.functionsPage();
  Check(WaitIdle(tab), "functions: the first row drawn without waiting on the window");
  {
    const std::vector<front::PhysicsFunction> rows = front::PhysicsFunctions();
    Check(fp->functionChoice()->count() == (int)rows.size(), "functions: every row of PhysicsFunctions in the drop-down (" +
                                                                 std::to_string(rows.size()) + ")");
    // the penetrability over e_cm from -1 to 2 MeV: the points <= 0 outside the domain, gaps with their words
    int row = -1;
    for (int i = 0; i < (int)rows.size(); i++)
      if (rows[i].name == "penetrability") row = i;
    fp->functionChoice()->setCurrentIndex(row);
    WaitIdle(tab);
    const front::PhysicsFunction &f = rows[row];
    int eArg = -1, lArg = -1;
    for (int k = 0; k < (int)f.arguments.size(); k++) {
      if (f.arguments[k].name == "e_cm") eArg = k;
      if (f.arguments[k].name == "l") lArg = k;
    }
    Check(fp->form()->swept() == eArg, "functions: the row's axis e_cm is the swept argument");
    fp->form()->sweep()->from->setText("-1");
    fp->form()->sweep()->to->setText("2");
    fp->form()->sweep()->points->setValue(31);
    fp->form()->sweep()->log->setChecked(false);
    QPushButton *drawButton = fp->findChild<QPushButton *>("functionsDraw");
    drawButton->click();
    WaitIdle(tab);
    std::vector<std::vector<double>> columns;
    QVector<double> grid;
    QString why;
    fp->form()->columns(columns, grid, why, false);
    front::PhysicsValues direct;
    std::string w;
    front::EvaluatePhysics(f.name, columns, direct, w);
    const QList<PlotSeries> &drawn = fp->panel()->plot()->series();
    Check(drawn.size() == 1 && Differ(drawn[0].modelY, direct.values[0]) == 0 && Differ(drawn[0].modelX, columns[eArg]) == 0,
          "functions: the drawn P_l(e_cm) is EvaluatePhysics's, element by element");
    int gaps = 0;
    for (double v : drawn.value(0).modelY) gaps += std::isnan(v);
    const QString status = fp->panel()->status();
    Check(gaps == 11 && status.contains("11 of 31") && status.contains("e_cm outside its domain"),
          "functions: 11 points at e_cm <= 0 are gaps, counted with their words: " + status.toStdString());
    // l = 2 added over it
    fp->form()->setValue(lArg, 2);
    fp->findChild<QPushButton *>("functionsAdd")->click();
    WaitIdle(tab);
    fp->form()->columns(columns, grid, why, false);
    front::EvaluatePhysics(f.name, columns, direct, w);
    const QList<PlotSeries> &two = fp->panel()->plot()->series();
    Check(two.size() == 2 && Differ(two[1].modelY, direct.values[0]) == 0 && two[1].label.contains("l = 2"),
          "functions: a curve of l = 2 added, its values the facade's (" + (two.size() == 2 ? two[1].label.toStdString() : "") + ")");
    Shot(tab, shots + "/physics_functions_light.png");
  }

  // stopping
  tab.showPage(PhysicsTab::Stopping);
  Settle();
  {
    FunctionsPage *sp = tab.stoppingPage();
    WaitIdle(tab);
    const std::vector<front::PhysicsFunction> rows = front::StoppingFunctions();
    sp->form()->setText(0, "Ta2O5");
    sp->findChild<QPushButton *>("stoppingDraw")->click();
    WaitIdle(tab);
    std::vector<std::vector<double>> columns;
    QVector<double> grid;
    QString why;
    sp->form()->columns(columns, grid, why, true);
    front::PhysicsValues direct;
    std::string w;
    const int row = sp->functionChoice()->currentIndex();
    const bool ok = front::EvaluateStopping(rows[row].name, "Ta2O5", columns, direct, w);
    const QList<PlotSeries> &drawn = sp->panel()->plot()->series();
    Check(ok && drawn.size() == 1 && Differ(drawn[0].modelY, direct.values[0]) == 0,
          "stopping: " + rows[row].name + " of Ta2O5 drawn as EvaluateStopping gives it (" + w + ")");
    QLabel *citation = sp->findChild<QLabel *>("stoppingPanelCitation");
    Check(citation && citation->isVisible() && citation->text().contains(QString::fromStdString(front::StoppingCitation())),
          "stopping: SRIM's citation under the plot");
    Shot(tab, shots + "/physics_stopping_light.png");
  }

  // the project pages without a project
  for (int page : {PhysicsTab::Channels, PhysicsTab::RMatrix, PhysicsTab::PoleMap, PhysicsTab::Thm}) {
    tab.showPage(page);
    Settle();
  }
  Check(tab.channelsPage()->needsText().contains("open a project") && !tab.rmatrixPage()->contentShown() &&
            !tab.poleMapPage()->contentShown() && !tab.thmPage()->contentShown(),
        "no project: the project pages say what they need: " + tab.channelsPage()->needsText().left(60).toStdString());
  tab.showPage(PhysicsTab::RMatrix);
  Settle();
  Shot(tab, shots + "/physics_noproject_light.png");

  // ---------------------------------------------------------------- (b) T1
  const fs::path t1 = CopyProject(source / "tests/upgrade/T1/T1.azr", work / "T1");
  {
    front::Project p;
    std::string why;
    Check(front::Project::Open(t1.string(), p, why), "T1 opened");
    holder.setProject(std::move(p));
  }
  tab.showPage(PhysicsTab::Channels);
  Check(WaitIdle(tab) && tab.context()->info().ok, "T1: the session set up in the worker (" +
                                                       tab.context()->info().why.toStdString() + ")");
  front::PhysicsSession direct;
  {
    std::string why;
    Check(front::PhysicsSession::Open(t1.string(), "", false, direct, why), "T1: a direct session for the comparisons");
  }
  {
    ChannelsPage *cp = tab.channelsPage();
    Check(cp->contentShown() && SelectData(cp->pair(), 1) >= 0, "channels: pair 1 offered (" +
                                                                   cp->pair()->currentText().toStdString() + ")");
    cp->l()->setValue(1);
    SelectData(cp->quantity(), "P");
    cp->energies()->set(0.1, 3, 120, true);
    cp->findChild<QPushButton *>("channelsDraw")->click();
    WaitIdle(tab);
    QVector<double> grid;
    QString why;
    cp->energies()->grid(grid, why);
    front::PhysicsCurves c;
    std::string w;
    direct.ChannelFunctions(1, 1, 0, -1, std::vector<double>(grid.begin(), grid.end()), c, w);
    const QList<PlotSeries> &drawn = cp->panel()->plot()->series();
    Check(drawn.size() == 1 && c.Find("P") && Differ(drawn[0].modelY, c.Find("P")->values) == 0,
          "channels: P of pair 1, l = 1 as ChannelFunctions gives it");
    Check(cp->panel()->logY()->isChecked(), "channels: P on a log axis, as the facade's plot suggests");
    Shot(tab, shots + "/physics_channels_light.png");
  }
  {
    tab.showPage(PhysicsTab::RMatrix);
    WaitIdle(tab);
    RMatrixPage *rp = tab.rmatrixPage();
    int half = -1;
    for (const front::PhysicsGroup &g : direct.Groups())
      if (g.J == 0.5 && g.parity == 1) half = g.jgroup;
    Check(half > 0 && SelectData(rp->group(), half) >= 0, "rmatrix: the 1/2+ group offered");
    WaitIdle(tab);
    rp->energies()->set(2.0, 2.8, 200, false);
    emit rp->energies()->changed();
    WaitIdle(tab);
    SelectData(rp->quantity(), "delta");
    Settle();
    QVector<double> grid;
    QString why;
    rp->energies()->grid(grid, why);
    front::PhysicsCurves b;
    std::string w;
    direct.Block(half, std::vector<double>(grid.begin(), grid.end()), b, w);
    const QList<PlotSeries> &drawn = rp->panel()->plot()->series();
    const front::PhysicsSeries *delta = b.Find("delta", {{"channel", 1}, {"channel_out", 1}});
    Check(!drawn.isEmpty() && delta && Differ(drawn[0].modelY, delta->values) == 0 &&
              Differ(drawn[0].modelX, b.axisValues) == 0,
          "rmatrix: the 1/2+ group's delta of channel 1 as Block gives it");
    // the pole table against Poles()
    std::vector<front::PhysicsPole> poles;
    for (const front::PhysicsPole &p : direct.Poles())
      if (p.jgroup == half) poles.push_back(p);
    QTableWidget *t = rp->poleTable();
    bool same = t->rowCount() == (int)poles.size() && !poles.empty();
    for (int r = 0; same && r < t->rowCount(); r++)
      same = t->item(r, 2)->data(Qt::UserRole).toDouble() == poles[r].energy &&
             t->item(r, 3)->data(Qt::UserRole).toDouble() == poles[r].totalWidth &&
             t->columnCount() == 4 + (int)direct.Groups()[half - 1].channels.size();
    Check(same, "rmatrix: the pole table is Poles() of the group (E_x, Gamma, a column per channel)");
    const QVector<double> &marks = rp->panel()->plot()->markers();
    Check(!poles.empty() && marks.size() == (int)poles.size() && marks[0] == poles[0].energy,
          "rmatrix: a dashed marker at each pole");
    // a click on the pole centres the plot on it
    QRect cell = t->visualItemRect(t->item(0, 2));
    QMouseEvent press(QEvent::MouseButtonPress, cell.center(), t->viewport()->mapToGlobal(cell.center()), Qt::LeftButton,
                      Qt::LeftButton, Qt::NoModifier);
    QMouseEvent release(QEvent::MouseButtonRelease, cell.center(), t->viewport()->mapToGlobal(cell.center()), Qt::LeftButton,
                        Qt::NoButton, Qt::NoModifier);
    QApplication::sendEvent(t->viewport(), &press);
    QApplication::sendEvent(t->viewport(), &release);
    Settle();
    const QRectF v = rp->panel()->plot()->view();
    Check(!poles.empty() && std::fabs(v.center().x() - poles[0].energy) < 1e-9 && v.width() < 0.8,
          "rmatrix: a click on the pole centres the plot on E_x = " + std::to_string(poles.empty() ? 0 : poles[0].energy));
    Shot(tab, shots + "/physics_rmatrix_light.png");
    // the Argand plot: Im U against Re U of the diagonal, markers at the levels
    SelectData(rp->quantity(), "argand");
    Settle();
    const QList<PlotSeries> &argand = rp->panel()->plot()->series();
    const front::PhysicsSeries *ure = b.Find("U_re", {{"channel", 1}, {"channel_out", 1}}),
                               *uim = b.Find("U_im", {{"channel", 1}, {"channel_out", 1}});
    Check(!argand.isEmpty() && ure && uim && Differ(argand[0].modelX, ure->values) == 0 &&
              Differ(argand[0].modelY, uim->values) == 0 && argand[0].dataX.size() == 1,
          "rmatrix: the Argand plot's x = Re U_11, y = Im U_11, one marker at the level in range");
    Shot(tab, shots + "/physics_argand_light.png");
    SelectData(rp->quantity(), "delta");
  }
  {
    tab.showPage(PhysicsTab::PoleMap);
    WaitIdle(tab);
    PoleMapPage *pm = tab.poleMapPage();
    SelectData(pm->kind(), "map_det");
    pm->energies()->set(1.5, 4.5, 300, false);
    pm->findChild<QPushButton *>("poleMapDraw")->click();
    WaitIdle(tab);
    QVector<double> grid;
    QString why;
    pm->energies()->grid(grid, why);
    const std::vector<front::PhysicsGroup> groups = direct.Groups();
    HeatmapWidget *map = pm->map();
    bool rows = map->rows().size() == (int)groups.size() && Differ(map->x(), std::vector<double>(grid.begin(), grid.end())) == 0;
    for (std::size_t g = 0; rows && g < groups.size(); g++) {
      front::PhysicsCurves b;
      std::string w;
      direct.Block(groups[g].jgroup, std::vector<double>(grid.begin(), grid.end()), b, w);
      rows = rows && b.Find("log10_abs_det") && Differ(map->values()[g], b.Find("log10_abs_det")->values) == 0 &&
             map->rows()[g] == JPi(groups[g].J, groups[g].parity);
    }
    Check(rows, "pole map: one row per J group, each Block's log10 |det A^-1| on the grid (" +
                    std::to_string(map->rows().size()) + " rows)");
    Check(map->points().size() == (int)direct.Poles().size(), "pole map: every pole drawn as a point on its group's row");
    // the cursor's read-out of one cell: column 150, row 1
    const QRect r = map->mapRect();
    const int col = 150, row = 1;
    const QRectF view = map->view();
    const QPointF pixel(r.left() + (map->x()[col] - view.left()) / view.width() * r.width(),
                        r.top() + (row - view.top()) / view.height() * r.height());
    const QString readout = map->readoutAt(pixel);
    Check(readout.contains(QString::number(map->x()[col], 'g', 6)) &&
              readout.contains(QString::number(map->values()[row][col], 'g', 5)) && readout.contains(map->rows()[row]),
          "pole map: the read-out at a cell is its (x, row, value): " + readout.toStdString());
    // a click on row 1: its cut below
    QMouseEvent press(QEvent::MouseButtonPress, pixel, map->mapToGlobal(pixel), Qt::LeftButton, Qt::LeftButton, Qt::NoModifier);
    QMouseEvent release(QEvent::MouseButtonRelease, pixel, map->mapToGlobal(pixel), Qt::LeftButton, Qt::NoButton, Qt::NoModifier);
    QApplication::sendEvent(map, &press);
    QApplication::sendEvent(map, &release);
    Settle();
    const QList<PlotSeries> &cut = pm->cut()->plot()->series();
    Check(cut.size() == 1 && cut[0].modelY == map->values()[row] && cut[0].label == map->rows()[row],
          "pole map: a click on a row draws its cut below (" + (cut.isEmpty() ? std::string() : cut[0].label.toStdString()) + ")");
    Check(pm->explanation()->text().contains("real-energy map") && pm->explanation()->text().contains("not available"),
          "pole map: the page says it is a real-energy map and that the complex plane is not available");
    Check(map->colourAt(0).lightness() < map->colourAt(0.5).lightness() &&
              map->colourAt(0.5).lightness() < map->colourAt(1).lightness(),
          "pole map: the colour map's lightness rises from low to high values");
    Shot(tab, shots + "/physics_polemap_light.png");
    // |U| of one group
    SelectData(pm->kind(), "map_U");
    WaitIdle(tab);
    Check(!map->rows().isEmpty() && map->points().isEmpty(), "pole map: |U_cc'| of one group, one row per channel pair (" +
                                                                 std::to_string(map->rows().size()) + ")");
  }
  // the last fit: none yet, so the start values stay
  {
    QComboBox *params = tab.parameterChoice();
    params->setCurrentIndex(1);
    emit params->activated(1);
    Settle();
    Check(params->currentIndex() == 0 && tab.context()->parameterFile().isEmpty(),
          "parameters: no param.sav in the output directory: the start values kept");
  }

  // ---------------------------------------------------------------- (c) T2: THM
  const fs::path t2 = CopyProject(source / "tests/upgrade/T2/T2.azr", work / "T2");
  {
    front::Project p;
    std::string why;
    front::Project::Open(t2.string(), p, why);
    holder.setProject(std::move(p));
    tab.showPage(PhysicsTab::Thm);
    Check(WaitIdle(tab) && tab.thmPage()->contentShown(), "T2: the THM page shows the experiments");
    ThmPage *tp = tab.thmPage();
    Check(tp->experiment()->currentText() == "LC09", "thm: the experiment LC09 offered");
    SelectData(tp->part(), "ThmVertex");
    SelectData(tp->quantity(), "M2");
    tp->energies()->set(0.1, 1.0, 40, false);
    tp->findChild<QPushButton *>("thmDraw")->click();
    WaitIdle(tab);
    front::PhysicsSession s2;
    front::PhysicsSession::Open(t2.string(), "", false, s2, why);
    QVector<double> grid;
    QString gwhy;
    tp->energies()->grid(grid, gwhy);
    front::PhysicsCurves v;
    s2.ThmVertex("LC09", std::vector<double>(grid.begin(), grid.end()), v, why);
    const QList<PlotSeries> &drawn = tp->panel()->plot()->series();
    int k = 0;
    bool same = !drawn.isEmpty();
    for (const front::PhysicsSeries &s : v.series)
      if (s.name == "M2") same = same && k < drawn.size() && Differ(drawn[k++].modelY, s.values) == 0;
    Check(same && k == drawn.size(), "thm: every <|M_l|^2> line of LC09 as ThmVertex gives it (" + std::to_string(k) + " lines)");
    Shot(tab, shots + "/physics_thm_light.png");
  }

  // ---------------------------------------------------------------- (d) dark look, the window
  setAzureTheme(ThemeChoice::Dark);
  Settle();
  {
    front::Project p;
    std::string why;
    front::Project::Open(t1.string(), p, why);
    holder.setProject(std::move(p));
  }
  for (auto page : std::vector<std::pair<int, QString>>{{PhysicsTab::Functions, "functions"},
                                                        {PhysicsTab::Channels, "channels"},
                                                        {PhysicsTab::RMatrix, "rmatrix"},
                                                        {PhysicsTab::PoleMap, "polemap"},
                                                        {PhysicsTab::Stopping, "stopping"}}) {
    tab.showPage(page.first);
    WaitIdle(tab);
    if (page.first == PhysicsTab::Channels) tab.channelsPage()->draw(false), WaitIdle(tab);
    if (page.first == PhysicsTab::PoleMap) SelectData(tab.poleMapPage()->kind(), "map_det"), WaitIdle(tab);
    Shot(tab, shots + "/physics_" + page.second + "_dark.png");
  }
  setAzureTheme(ThemeChoice::Light);
  {
    AZURESetup window;
    QTabWidget *tabs = window.findChild<QTabWidget *>();
    int at = -1;
    for (int i = 0; tabs && i < tabs->count(); i++)
      if (tabs->tabText(i).remove('&') == "Physics") at = i;
    Check(at > 0 && tabs->tabText(at - 1).remove('&') == "Plot", "the window: a Physics tab after Plot");
  }

  std::cout << (failures ? "FAILED: " + std::to_string(failures) : std::string("all passed")) << std::endl;
  return failures ? 1 : 0;
}
