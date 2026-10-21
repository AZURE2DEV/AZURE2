#ifndef AZURESETUP_H
#define AZURESETUP_H

#include <QMainWindow>
#include <QMap>
#include <QPair>

#include "PairsTab.h"
#include "LevelsTab.h"
#include "SegmentsTab.h"
#include "TargetIntTab.h"
#include "Config.h"

class RunTab;
class FittingTab;
class NuclearPotentialTab;
#ifdef USE_QWT
class PlotTab;
#endif
#ifdef USE_MCMC
class MCMCTab;
#endif
class AZUREMainThread;

QT_BEGIN_NAMESPACE

class QTabWidget;
class QMenu;
class QAction;
class QActionGroup;
class QTextEdit;

QT_END_NAMESPACE

class Directories {
 public:
  Directories() :
    outputDir(QString("")),
    checksDir(QString("")) {};
  QString outputDir;
  QString checksDir;
};

/*!
 * The setup program's main window: it owns every tab, and reads and writes the .azr project file.
 *
 * The tabs are the sections of that file -- particle pairs, levels, segments, target integration -- plus the ones that drive a run: Run, Fitting, MCMC, Plot and Nuclear Potential. Saving walks them in file order, so a tab's read and write methods are the authority on its section's format.
 */
class AZURESetup : public QMainWindow {
  Q_OBJECT

 public:
  AZURESetup();
  Config &GetConfig();
  void open(QString filename);
  void saveProject() { save(); }  // Public wrapper for MCMCTab to use
  /*! Reads / writes a project file without any dialog (the headless GUI tests). */
  bool readProjectFile(const QString &filename) { return readFile(filename); }
  bool writeProjectFile(const QString &filename) { return writeFile(filename); }
  SegmentsTab *getSegmentsTab() const { return segmentsTab; }
  PairsTab *getPairsTab() const { return pairsTab; }
  /*! Grabs the window once per tab into dir/tab<n>_<name>.png (see AZURE2_SCREENSHOT_DIR in main.cpp). */
  void saveTabScreenshots(const QString &dir);

  // Convert RWA parameter to physical value using proper R-Matrix transformation
  double ConvertRWAToPhysical(const QString &paramName, double rwaValue);

  // Batch convert multiple RWA parameters to physical values (much more efficient)
  std::vector<double> BatchConvertRWAToPhysical(const QStringList &paramNames, const std::vector<double> &rwaValues);

  // Batch convert RWA to physical using OLD compound structure (for loading old param files)
  // oldLevelChannelCounts maps energyIndex -> number of channels in that level from the old structure
  std::vector<double> BatchConvertRWAToPhysicalWithOldStructure(
      const QStringList &paramNames, const std::vector<double> &rwaValues, const QMap<int, int> &oldLevelChannelCounts);

  /*!
   * The engine's parameter names (energy_N, width_N_M) of the project file on disk, keyed by the Levels tab's
   * model rows (level row, channel row or -1 for the energy).  Matched by identity -- J, parity, nucleus and
   * energy of a level; pair, l and s of a channel -- so a .sav written by a run of that file lands on the right
   * levels whatever the order of the file and of the GUI's rows.  Empty if the file cannot be read.
   */
  QMap<QPair<int, int>, QString> EngineParameterNamesOnDisk();

  // Getter for FittingTab (for MCMCTab access)
  FittingTab *getFittingTab() const { return fittingTab; }
  TargetIntTab *getTargetIntTab() const { return targetIntTab; }
  LevelsTab *getLevelsTab() const { return levelsTab; }

 public slots:
  void SaveAndRun();
#ifdef USE_MCMC
  void SaveAndRunMCMC();
  void DeleteMCMCThread();
#endif
  void DeleteThread();

 private slots:
  void reset();
  void open();
  void openRecent();
  void clearRecent();
  void save();
  void saveAs();
  void matrixChanged(QAction *action);
  void editChecks();
  void editDirs();
  void editOptions();
  void buildCoulombTable();
  void showAbout();
  void showTabInfo();
  void openWebsite();

 private:
  bool readFile(QString filename);
  bool readConfig(QTextStream &inStream);
  bool writeFile(QString filename);
  bool writeConfig(QTextStream &outStream, QString directory);
  bool readLastRun(QTextStream &inStream);
  bool writeLastRun(QTextStream &outStream);
  void createActions();
  void createMenus();
  void updateRecent();
  void updateNuclearPotentialTabVisibility();  // Show/hide Nuclear Potential tab based on config

  Config config;

  QAction *aboutAction;
  QAction *resetAction;
  QAction *quitAction;
  QAction *openAction;
  QAction *saveAction;
  QAction *saveAsAction;
  QAction *editChecksAction;
  QAction *editDirsAction;
  QAction *copyAction;
  QAction *aMatrixAction;
  QAction *rMatrixAction;
  QAction *editOptionsAction;
  QAction *buildCoulombTableAction;
  QAction *recentSeparator;
  QAction *clearRecentAction;
  enum { numRecent = 5 };
  QAction *recentFileActions[numRecent];
  QAction *showTabInfoAction;
  QAction *openAZURESiteAction;

  QActionGroup *matrixActionGroup;

  QMenu *fileMenu;
  QMenu *editMenu;
  QMenu *configMenu;
  QMenu *formalismMenu;
  QMenu *recentFileMenu;
  QMenu *helpMenu;
  QTabWidget *tabWidget;
  bool fileIncludeAllJPi = false;
  int fileMaxL = -1;
  double fileAbundanceSumSigma = 0.0;
  PairsTab *pairsTab;
  LevelsTab *levelsTab;
  SegmentsTab *segmentsTab;
  TargetIntTab *targetIntTab;
  NuclearPotentialTab *nuclearPotentialTab;
  int nuclearPotentialTabIndex;  // Store tab index for show/hide
  FittingTab *fittingTab;
  RunTab *runTab;
  AZUREMainThread *azureMain;
#ifdef USE_QWT
  PlotTab *plotTab;
#endif
#ifdef USE_MCMC
  MCMCTab *mcmcTab;
  class AZUREMCMCThread *azureMCMC;
#endif
};

#endif
