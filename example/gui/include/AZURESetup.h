#ifndef AZURESETUP_H
#define AZURESETUP_H

#include "PairsTab.h"
#include "LevelsTab.h"
#include "SegmentsTab.h"
#include "TargetIntTab.h"
#include "Config.h"

class RunTab;
#ifdef USE_QWT
class PlotTab;
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
 Directories() : outputDir(QString("")), checksDir(QString("")) {};
  QString outputDir;
  QString checksDir;
};

class AZURESetup : public QMainWindow {
  
 Q_OBJECT

 public:
  AZURESetup();
  Config& GetConfig();
  void open(QString filename);

 public slots:
  void SaveAndRun();
  void DeleteThread();

 private slots:
  void reset();
  void open();
  void openRecent();
  void clearRecent();
  void save();
  void saveAs();
  void matrixChanged(QAction* action);
  void editChecks();
  void editDirs(); 
  void editOptions();
  void showAbout();
  void showTabInfo();
  void openWebsite();

 private:
  bool readFile(QString filename);
  bool readConfig(QTextStream& inStream);
  bool writeFile(QString filename);
  bool writeConfig(QTextStream& outStream, QString directory);
  bool readLastRun(QTextStream& inStream);
  bool writeLastRun(QTextStream& outStream);
  void createActions();
  void createMenus();
  void updateRecent();

  Config config;

  QAction* aboutAction;
  QAction* resetAction;
  QAction* quitAction;
  QAction* openAction;
  QAction* saveAction;
  QAction* saveAsAction;
  QAction* editChecksAction;
  QAction* editDirsAction;
  QAction* copyAction;
  QAction* aMatrixAction;
  QAction* rMatrixAction;
  QAction* editOptionsAction;
  QAction* recentSeparator;
  QAction* clearRecentAction;
  enum { numRecent = 5 };
  QAction* recentFileActions[numRecent];
  QAction* showTabInfoAction; 
  QAction* openAZURESiteAction;

  QActionGroup* matrixActionGroup;

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
  RunTab *runTab;
  AZUREMainThread *azureMain;
#ifdef USE_QWT
  PlotTab* plotTab;
#endif
};

#endif 
