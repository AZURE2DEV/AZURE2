#ifndef AZURESETUP_H
#define AZURESETUP_H

#include <QMainWindow>
#include <QThread>

#include "ProjectHolder.h"
#include "RunWorker.h"

class PairsTab;
class LevelsTab;
class SectionsTab;
class FittingTab;
class RunTab;
class PlotTab;
class DiagnosticsPanel;

QT_BEGIN_NAMESPACE
class QAction;
class QActionGroup;
class QDockWidget;
class QMenu;
class QTabWidget;
class QToolButton;
QT_END_NAMESPACE

/*!
 * The setup program's main window (docs/upgrade/qt_gui.md).  It owns the open project (ProjectHolder: the file's
 * lossless document, edited with the engine's grammars through azure2::front) and the tabs the old program had, in
 * its order: Particle Pairs (with the per-channel radii), Levels and Channels (with gamma decays and isospin),
 * Segments (from data, without data, model data), Experimental Effects, Nuclear Potential, THM, Fitting Settings and
 * Calculate.  A run is the program's run (azure2::front::Run) in a worker thread; the window stays responsive and a
 * run can be stopped.  Diagnostics come from the engine (Validate) and are listed in the panel below the tabs.
 */
class AZURESetup : public QMainWindow {
  Q_OBJECT

 public:
  AZURESetup();
  ~AZURESetup() override;

  /// Open a project: a format-3 file, or an old one after asking to convert it (the converter of --convert).
  void open(QString filename);
  /// Read / write a project without any dialog (the headless tests); false with the reason.
  bool readProjectFile(const QString &filename, QString &why);
  bool writeProjectFile(const QString &filename, QString &why);
  /// Write what the tabs hold into the project (the levels when they were edited, the text sections).
  bool commitTabs(QString &why);
  /// Grab the window once per tab into dir/tab<n>_<name>.png (AZURE2_QT_SHOTS, the tests).
  void saveTabScreenshots(const QString &dir);
  /// Start a run of the Calculate tab's command (after saving); false with the reason when it cannot start.
  bool startRun(QString &why);
  bool running() const { return running_; }

  ProjectHolder *holder() const { return holder_; }
  PairsTab *getPairsTab() const { return pairsTab; }
  LevelsTab *getLevelsTab() const { return levelsTab; }
  SectionsTab *getSegmentsTab() const { return segmentsTab; }
  SectionsTab *getEffectsTab() const { return effectsTab; }
  FittingTab *getFittingTab() const { return fittingTab; }
  RunTab *getRunTab() const { return runTab; }
  PlotTab *getPlotTab() const { return plotTab; }
  DiagnosticsPanel *getDiagnostics() const { return diagnostics; }
  QTabWidget *tabs() const { return tabWidget; }

 signals:
  void runFinished(const RunOutcome &outcome);
  void startWorker(const QString &command, const QStringList &options, const QString &project, int threads);

 public slots:
  void SaveAndRun();
  void stopRun();

 private slots:
  void newProject();
  void open();
  void openRecent();
  void clearRecent();
  void save();
  void saveAs();
  void convertProject();
  void validate();
  void matrixChanged(QAction *action);
  void editChecks();
  void editDirs();
  void editOptions();
  void showAbout();
  void showTabInfo();
  void openWebsite();
  void projectOpened();
  void levelsEdited();
  void runEnded(const RunOutcome &outcome);
  void showRecord(const QString &section, const QString &record);

 private:
  void createActions();
  void createMenus();
  void updateRecent();
  void addRecent(const QString &file);
  void updateTitle();
  /// Fill the Pairs and Levels tabs from the project's <levels> rows.
  bool loadLevels(QString &why);
  /// The logo with the polarization arrow while a data set or extrapolation is a polarization observable.
  void updateLogo();
  bool maybeSave();
  void closeEvent(QCloseEvent *event) override;

  ProjectHolder *holder_;
  bool levelsDirty_ = false;  ///< the Pairs or Levels tab changed since the rows were read
  bool loadingLevels_ = false;
  bool loadedIncludeAllJPi_ = false;  ///< <config> includeAllJPi as read (the Levels tab's box)
  int loadedMaxL_ = 0;                ///< the Levels tab's maximum l as read
  bool running_ = false;
  QThread runThread_;
  RunWorker *worker_ = nullptr;

  QToolButton *logoButton_ = nullptr;
  bool polarizedLogo_ = false;
  QAction *aMatrixAction, *rMatrixAction;
  QActionGroup *matrixActionGroup;
  enum { numRecent = 5 };
  QAction *recentFileActions[numRecent];
  QAction *recentSeparator, *clearRecentAction;
  QMenu *recentFileMenu;
  QTabWidget *tabWidget;
  PairsTab *pairsTab;
  LevelsTab *levelsTab;
  SectionsTab *segmentsTab;
  SectionsTab *effectsTab;
  SectionsTab *potentialTab;
  SectionsTab *thmTab;
  FittingTab *fittingTab;
  RunTab *runTab;
  PlotTab *plotTab;
  DiagnosticsPanel *diagnostics;
  QDockWidget *diagnosticsDock;
};

#endif
