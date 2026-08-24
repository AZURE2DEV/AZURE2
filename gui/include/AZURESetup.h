#ifndef AZURESETUP_H
#define AZURESETUP_H

#include <QMainWindow>
#include <QTemporaryFile>

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

  // Convert RWA parameter to physical value using proper R-Matrix transformation
  /*! A Config whose configfile is a snapshot of the current GUI state.

      CNuc::Fill reads the .azr from disk, so a conversion built straight on the
      live Config sees the last *saved* file.  That matters most for the
      per-channel reduced-width-amplitude flag, which decides whether a width is
      converted at all: toggling it without saving would otherwise change the
      value by orders of magnitude.  Falls back to the unsnapshotted Config if
      the temporary file cannot be written, so conversion degrades rather than
      fails. */
  struct GuiStateSnapshot {
    explicit GuiStateSnapshot(AZURESetup *setup);
    QTemporaryFile file;  // removed when this goes out of scope
    Config config;
    bool ok;
  };

  double ConvertRWAToPhysical(const QString &paramName, double rwaValue);

  // Batch convert multiple RWA parameters to physical values (much more efficient)
  std::vector<double> BatchConvertRWAToPhysical(const QStringList &paramNames, const std::vector<double> &rwaValues);

  // Batch convert RWA to physical using OLD compound structure (for loading old param files)
  // oldLevelChannelCounts maps energyIndex -> number of channels in that level from the old structure
  std::vector<double> BatchConvertRWAToPhysicalWithOldStructure(
      const QStringList &paramNames, const std::vector<double> &rwaValues, const QMap<int, int> &oldLevelChannelCounts);

  // Getter for FittingTab (for MCMCTab access)
  FittingTab *getFittingTab() const { return fittingTab; }

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
