#ifndef AZURESETUP_H
#define AZURESETUP_H

#include <QMainWindow>
#include <QMap>
#include <QUrl>
#include <QThread>

#include "ProjectHolder.h"
#include "RunWorker.h"

class PairsTab;
class LevelsTab;
class SectionsTab;
class FittingTab;
class RunTab;
class PlotTab;
class ThmTab;
class PythonTools;
class CoulombTableControl;
class MCMCTab;
class DiagnosticsPanel;
class StartPage;
class NotesBar;

QT_BEGIN_NAMESPACE
class QAction;
class QActionGroup;
class QDockWidget;
class QLabel;
class QMenu;
class QProgressBar;
class QStackedWidget;
class QTabWidget;
class QToolBar;
class QToolButton;
QT_END_NAMESPACE

/*!
 * The setup program's main window (docs/upgrade/qt_gui.md).  It owns the open project (ProjectHolder: the file's
 * lossless document, edited with the engine's grammars through azure2::front) and the tabs the old program had, in
 * its order: Particle Pairs (with the per-channel radii), Levels and Channels (with gamma decays and isospin),
 * Segments (from data, without data, model data), Experimental Effects, Nuclear Potential, THM, Fitting Settings and
 * Calculate.  A run is the program's run (azure2::front::Run) in a worker thread; the window stays responsive and a
 * run can be stopped.  Diagnostics come from the engine (Validate) and are listed in a dock below the tabs, hidden
 * while there is nothing to show and opened by a problem.
 *
 * The shell (docs/upgrade/qt_gui.md section 11): a start page while no project is open, a tool bar of the main actions
 * (New, Open, Save, Validate, Calculate, Fit, Stop, Plot), a status bar (the file and its modified mark, the
 * validity, the running job with its progress and a Stop button), the View menu's look (System, Light, Dark), the
 * window's state and the recent files kept in QSettings, an .azr dropped on the window opened.
 */
class PotentialPanel;

class AZURESetup : public QMainWindow {
  Q_OBJECT

 public:
  AZURESetup();
  ~AZURESetup() override;

  /// Open a project: a format-3 file, or an old one after asking to convert it (the converter of --convert).
  void open(QString filename);
  /// Convert the old project `from` to `to` (AZURE2 --convert) and open the result.  A line the converter cannot read
  /// is shown with its number, the fields expected and the readings proposed (ConvertRepairDialog), to be edited and
  /// converted again from a copy (`from` is kept): qt_gui.md section 19.5.  False with the converter's report.
  bool convertOldProject(const QString &from, const QString &to, QString &report);
  /// Read / write a project without any dialog (the headless tests); false with the reason.
  bool readProjectFile(const QString &filename, QString &why);
  bool writeProjectFile(const QString &filename, QString &why);
  /// Write what the tabs hold into the project (the levels when they were edited, the text sections).
  bool commitTabs(QString &why);
  /// The Pairs and Levels tabs' part of commitTabs: the <levels> rows when they were edited, includeAllJPi and
  /// maxOrbitalMomentum when changed (the Potential panel asks for it before the engine reads the pairs).
  bool commitLevels(QString &why);
  /// The l up to which includeAllJPi's groups go as the engine will read the project after commitLevels:
  /// maxOrbitalMomentum when it is (or will be) written, else -1 (the largest l of the <levels> rows).
  int keyMaxOrbitalMomentum() const;
  /// The Potential panel of the Pairs tab.
  PotentialPanel *potentialPanel() const { return potentialPanel_; }
  /// Grab the window once per tab into dir/tab<n>_<name>.png (AZURE2_QT_SHOTS, the tests).
  void saveTabScreenshots(const QString &dir);
  /// Start a run of the Calculate tab's command (after saving); false with the reason when it cannot start.
  bool startRun(QString &why);
  /// The MCMC tab's run: the command mcmc with the tab's options (and its priors written into the project).
  bool startMCMC(QString &why);
  bool running() const { return running_; }

  ProjectHolder *holder() const { return holder_; }
  PairsTab *getPairsTab() const { return pairsTab; }
  LevelsTab *getLevelsTab() const { return levelsTab; }
  SectionsTab *getSegmentsTab() const { return segmentsTab; }
  SectionsTab *getEffectsTab() const { return effectsTab; }
  ThmTab *getThmTab() const { return thmTab; }
  FittingTab *getFittingTab() const { return fittingTab; }
  RunTab *getRunTab() const { return runTab; }
  PlotTab *getPlotTab() const { return plotTab; }
  /// The MCMC tab, null when the engine has no mcmc command.
  MCMCTab *getMCMCTab() const { return mcmcTab; }
  PythonTools *tools() const { return tools_; }
  /// The Coulomb table of this computer: its notices (start page, Overview) and its dialog (qt_gui.md section 19.12).
  CoulombTableControl *coulombTable() const { return coulomb_; }
  DiagnosticsPanel *getDiagnostics() const { return diagnostics; }
  QTabWidget *tabs() const { return tabWidget; }
  StartPage *startPage() const { return startPage_; }
  /// A project is open (the tabs are shown; otherwise the start page).
  bool hasProject() const { return hasProject_; }
  /// The tool bar's actions by name: new, open, save, validate, calculate, fit, stop, plot.
  QAction *action(const QString &name) const;
  QToolBar *mainToolBar() const { return toolBar_; }
  /// The status bar's parts: the file and modified mark, the validity (a button: opens the Diagnostics), the job.
  QLabel *statusFile() const { return statusFile_; }
  QToolButton *statusValidity() const { return statusValidity_; }
  QLabel *statusJob() const { return jobLabel_; }
  QProgressBar *statusProgress() const { return jobBar_; }
  QDockWidget *diagnosticsDockWidget() const { return diagnosticsDock; }
  /// Standard or expert view (qt_gui.md section 19.4): View > Show Expert Options, the start page's box and the status
  /// bar's switch; kept in QSettings.  It changes what is shown, never the project.
  void setExpertMode(bool on);
  bool expertMode() const;
  QAction *expertAction() const { return expertAction_; }
  /// The data files the active data sets name that this computer does not have at that path (as written in the
  /// records), "record: file" each (qt_gui.md section 19.10).
  QStringList missingDataFiles() const;
  /// Look for each missing data file by its name under `folder` (four levels deep) and point its data set to the file
  /// found (a path relative to the project when it is below the project's folder's parent, else absolute); returns how
  /// many were pointed, `still` gets those not found.  The project is modified, not saved.
  int locateDataFiles(const QString &folder, QStringList *still = nullptr);
  /// The notes about the project just opened (one bar, never a modal box: qt_gui.md section 19.2).
  NotesBar *notesBar() const { return notes_; }
  /// Close the project (after asking to save) and show the start page.
  bool closeProject();
  /// Run `command` (calculate, fit, ...) from the tool bar: the Calculate tab's command set, then startRun.
  bool runCommand(const QString &command, QString &why);
  /// Whether a dropped file would be opened (an .azr).
  static bool acceptsDrop(const QList<QUrl> &urls);
  /// Open the first .azr of dropped files (after asking to save); false when there is none.
  bool openDropped(const QList<QUrl> &urls);

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
  void showJobProgress(const QString &stage, qlonglong evaluation, qlonglong calls, double objective, double best);
  void diagnosticsDone(bool ok, int errors, int warnings);

 protected:
  void dragEnterEvent(QDragEnterEvent *event) override;
  void dropEvent(QDropEvent *event) override;

 private:
  void createToolBar();
  void createStatusBar();
  void updateActions();
  void showProject(bool open);
  void restoreWindowState();
  void saveWindowState();
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
  bool includeAllJPiSwitchedOn_ = false;  ///< includeAllJPi switched on in this window since the project was read
  bool running_ = false;
  CoulombTableControl *coulomb_ = nullptr;
  QThread runThread_;
  RunWorker *worker_ = nullptr;

  QToolButton *logoButton_ = nullptr;
  bool polarizedLogo_ = false;
  QAction *aMatrixAction, *rMatrixAction;
  QActionGroup *matrixActionGroup;
  enum { numRecent = 8 };
  QAction *recentFileActions[numRecent];
  QAction *recentSeparator, *clearRecentAction;
  QMenu *recentFileMenu;
  QTabWidget *tabWidget;
  PairsTab *pairsTab;
  LevelsTab *levelsTab;
  SectionsTab *segmentsTab;
  SectionsTab *effectsTab;
  SectionsTab *potentialTab;
  PotentialPanel *potentialPanel_ = nullptr;
  ThmTab *thmTab;
  FittingTab *fittingTab;
  RunTab *runTab;
  PlotTab *plotTab;
  MCMCTab *mcmcTab = nullptr;
  PythonTools *tools_ = nullptr;
  DiagnosticsPanel *diagnostics;
  QDockWidget *diagnosticsDock;
  NotesBar *notes_ = nullptr;
  QAction *expertAction_ = nullptr;
  QToolButton *statusExpert_ = nullptr;
  QWidget *physicsTab_ = nullptr;
  QMenu *toolsMenu_ = nullptr;
  QMap<QWidget *, QString> tabTips_;  ///< what each tab is for (its tool tip)
  /// The expert tabs (Nuclear Potential, THM, MCMC) shown in standard view only when the project uses them.
  void updateExpertTabs();
  bool openValidation_ = false;      ///< the next validation is the one of a project just opened (its notes)
  bool validateRequested_ = false;   ///< the next validation was asked for by the user (Validate)

  QStackedWidget *central_ = nullptr;
  StartPage *startPage_ = nullptr;
  bool hasProject_ = false;
  QToolBar *toolBar_ = nullptr;
  QMap<QString, QAction *> actions_;
  QLabel *statusFile_ = nullptr;
  QToolButton *statusValidity_ = nullptr;
  QLabel *jobLabel_ = nullptr;
  QProgressBar *jobBar_ = nullptr;
  QToolButton *jobStop_ = nullptr;
  QWidget *jobBox_ = nullptr;
};

#endif
