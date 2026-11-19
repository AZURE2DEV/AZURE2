#ifndef OVERVIEWPAGE_H
#define OVERVIEWPAGE_H

#include <QList>
#include <QMap>
#include <QPointer>
#include <QStringList>
#include <QWidget>
#include <memory>

#include "azure2/front/Engine.h"

QT_BEGIN_NAMESPACE
class QGridLayout;
class QLabel;
class QFrame;
class QPushButton;
class QScrollArea;
class QStackedWidget;
class QTableWidget;
class QThread;
class QTimer;
class QVBoxLayout;
QT_END_NAMESPACE

class BarChart;
class ChiSquaredTrace;
class ProjectHolder;

/*!
 * The Overview page (docs/upgrade/qt_gui.md section 15): the first tab of an open project.  What the evaluation is and
 * where its fit stands, at a glance: what to look at next (links to the tab and record), the project (pairs and
 * thresholds, levels per J^pi, data sets, effects, free parameters by kind, the basis and the options set, the
 * validity), the last run read from the output directory (command, time, stale mark, the fit's status, chi^2, N,
 * free parameters, chi^2 per degree of freedom, a sortable table per data set with a bar of chi^2/N, the parameters
 * at a limit or without an error) and two charts.  Everything comes from the facade (SummarizeProject, ReadLastRun,
 * ReadRunParameters); the page counts and divides, nothing else.  refresh() reads files only; the parameters' set-up
 * runs in a worker thread.  Cards in a scroll area, in three, two or one columns by the page's width.
 */
class OverviewPage : public QWidget {
  Q_OBJECT

 public:
  explicit OverviewPage(ProjectHolder *holder, QWidget *parent = nullptr);
  ~OverviewPage() override;

  /// One line of "What to look at next": its words and its link ("diagnostics", "run:calculate", "run:fit",
  /// "record:<section>:<name>", "plot:<name>", "tab:fitting").
  struct Step {
    QString text;
    QString link;
  };

  /// Read the project and the output directory again (files only); the parameters' worker after it when needed.
  void refresh();
  /// refresh() after a short quiet time (several edits in a row refresh once).
  void refreshLater();
  /// A run of the window is in progress: the output files are not read until it ends.
  void setRunning(bool running);
  /// The last validation (DiagnosticsPanel::validated).
  void setValidity(bool ok, int errors, int warnings);
  /// The fitter chosen in Fitting Settings.
  void setFitterChoice(const QString &fitter);
  /// A progress record of the window's run (RunProgress): kept for the convergence chart while a fit runs.
  void addProgress(const QString &stage, qlonglong evaluation, double best);
  void clearProgress();

  // What the cards show, for the tests.
  const azure2::front::ProjectSummary &projectSummary() const { return project_; }
  const azure2::front::RunSummary &runSummary() const { return run_; }
  /// The parameters of the last worker job (valid when haveParameters()).
  const azure2::front::RunParameters &runParameters() const { return parameters_; }
  bool haveParameters() const { return haveParameters_; }
  bool parametersBusy() const { return worker_ != nullptr; }
  const QList<Step> &steps() const { return steps_; }
  QTableWidget *dataSetTable() const { return table_; }
  /// The run card's numbers as shown ("chi2", "points", "free", "dof") and the stale mark.
  QString tileText(const QString &name) const;
  bool staleShown() const;
  /// The run card shows "No run yet" with its two buttons.
  bool noRunShown() const;
  /// Follow a step's link, as a click does.
  void followLink(const QString &link);

 signals:
  void commandRequested(const QString &command);
  void validateRequested();
  void diagnosticsRequested();
  void plotRequested(const QString &record);
  void recordRequested(const QString &section, const QString &record);
  void tabRequested(const QString &name);

 protected:
  void resizeEvent(QResizeEvent *) override;
  void showEvent(QShowEvent *) override;
  void changeEvent(QEvent *) override;

 private:
  void build();
  void reflow();
  void fillProject();
  void fillRun();
  void fillCharts();
  void fillSteps();
  void startParameters();
  void parametersDone();

  ProjectHolder *holder_;
  azure2::front::ProjectSummary project_;
  azure2::front::RunSummary run_;
  azure2::front::RunParameters parameters_;
  bool haveParameters_ = false;
  QString parametersKey_, parametersWhy_;
  QThread *worker_ = nullptr;
  bool pending_ = false;
  struct Job;
  std::shared_ptr<Job> job_;
  bool running_ = false;
  bool validated_ = false, valid_ = true;
  int errors_ = 0, warnings_ = 0;
  QString fitter_;
  QList<Step> steps_;
  int columns_ = 0;

  QTimer *refreshTimer_;
  QScrollArea *scroll_;
  QGridLayout *grid_;
  QFrame *stepsCard_, *projectCard_, *runCard_, *chartsCard_;
  QLabel *stepsText_, *projectHead_, *projectText_, *levelsTitle_, *projectMore_;
  QPushButton *validity_;
  BarChart *levelBars_;
  QStackedWidget *runStack_;
  QLabel *runHead_, *stale_, *runStatus_, *totalsText_, *flagsText_, *notesText_;
  QMap<QString, QLabel *> tiles_;
  QTableWidget *table_;
  BarChart *chiBars_;
  ChiSquaredTrace *trace_;
  QLabel *traceNote_, *traceTitle_;
};

#endif
