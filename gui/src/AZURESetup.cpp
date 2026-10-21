#include <QTextEdit>
#include <QMessageBox>
#include <QScrollBar>
#include <QFileDialog>
#include <QMenu>
#include <QMenuBar>
#include <QRadioButton>
#include <QAction>
#include <QActionGroup>
#include <QSettings>
#include <QTextStream>
#include <QDesktopServices>
#include <QStatusBar>
#include <algorithm>
#include <QProgressDialog>
#include <algorithm>
#include <QCoreApplication>
#include <QFileInfo>
#include <QElapsedTimer>
#include <atomic>
#include <thread>
#include <chrono>
#include "CoulombTable.h"

#include "AZURESetup.h"
#include "FittingTab.h"
#include "EditChecksDialog.h"
#include "EditDirsDialog.h"
#include "RunTab.h"
#include "Config.h"
#include "EditOptionsDialog.h"
#include "AZUREMainThread.h"
#include "AboutAZURE2Dialog.h"
#include "NuclearPotentialTab.h"
#include "NuclearPotentialManager.h"
#include "CNuc.h"
#include <cmath>
#include "PPair.h"
#include "AChannel.h"
#include "ALevel.h"
#include "JGroup.h"
#include "EData.h"
#include "AZUREParams.h"

#ifdef USE_QWT
#include "PlotTab.h"
#endif
#ifdef USE_MCMC
#include "MCMCTab.h"
#include "AZUREMCMCThread.h"
#include <QMetaType>
#include <vector>
#endif
#include <iostream>

struct SegPairs {
  int firstPair;
  int secondPair;
};
extern bool readSegmentFile(const Config &configure, std::vector<SegPairs> &segPairs);
extern bool checkExternalCapture(Config &configure, const std::vector<SegPairs> &segPairs);
extern void startMessage(const Config &configure);
extern void exitMessage(const Config &configure);


AZURESetup::AZURESetup() :
  config(std::cout) {
  // Set global Config pointer for use across codebase
  extern Config *g_config;
  g_config = &config;

#ifdef USE_MCMC
  // Register meta types for Qt signal/slot system
  qRegisterMetaType<std::vector<std::vector<double>>>("std::vector<std::vector<double>>");
#endif

  setMinimumSize(1000, 640);

  tabWidget = new QTabWidget();
  // The logo at the right end of the tab bar; clicking it opens the About box.
  QToolButton *logoButton = new QToolButton;
  logoButton->setIcon(QIcon(":/azure-icon.png"));
  logoButton->setIconSize(QSize(26, 26));
  logoButton->setText("AZURE2");
  logoButton->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
  logoButton->setAutoRaise(true);
  logoButton->setCursor(Qt::PointingHandCursor);
  logoButton->setToolTip(tr("AZURE2 -- about this program"));
  logoButton->setStyleSheet("QToolButton { border: none; padding: 0 8px; color: #1f4f8f; font-weight: bold; }"
                            "QToolButton:hover { background: #e3edf8; border-radius: 4px; }");
  connect(logoButton, &QToolButton::clicked, this, [this]() { AboutAZURE2Dialog about(this); about.exec(); });
  tabWidget->setCornerWidget(logoButton, Qt::TopRightCorner);

  pairsTab = new PairsTab;

  levelsTab = new LevelsTab;
  levelsTab->setPairsModel(pairsTab->getPairsModel());
  connect(pairsTab, SIGNAL(pairAdded(int)), levelsTab, SLOT(updateChannelsPairAddedEdited()));
  connect(pairsTab, SIGNAL(pairEdited(int)), levelsTab, SLOT(updateChannelsPairAddedEdited()));
  connect(pairsTab, SIGNAL(pairRemoved(int)), levelsTab, SLOT(updateChannelsPairRemoved(int)));
  connect(levelsTab, SIGNAL(readNewPair(PairsData, int, bool)), pairsTab, SLOT(addPair(PairsData, int, bool)));
  connect(levelsTab, SIGNAL(readExistingPair(PairsData, int, bool)), pairsTab, SLOT(editPair(PairsData, int, bool)));

  segmentsTab = new SegmentsTab;
  segmentsTab->setPairsModel(pairsTab->getPairsModel());
  segmentsTab->setAbundancesModel(pairsTab->getAbundancesModel());

  targetIntTab = new TargetIntTab;
  targetIntTab->setModels(segmentsTab->getSegmentsDataModel(), pairsTab->getPairsModel());

  nuclearPotentialTab = new NuclearPotentialTab();
  // The tab edits one particle pair at a time, so it needs the pair list.  The
  // pairs do not exist yet -- they are built while <levels> is read -- so the
  // tab follows the model from here rather than taking a snapshot.
  nuclearPotentialTab->setPairsModel(pairsTab->getPairsModel());
  // Pair keys are positional: deleting one shifts the rest down, and the
  // potentials have to move with them or they end up on the wrong channels.
  connect(pairsTab, SIGNAL(pairRemoved(int)), nuclearPotentialTab, SLOT(onPairRemoved(int)));

  fittingTab = new FittingTab();
  fittingTab->setTabReferences(levelsTab, segmentsTab);
  fittingTab->setAbundancesModel(pairsTab->getAbundancesModel());
  fittingTab->setConfig(&GetConfig());
  // Segment lines are numbered by position: what names them by number follows a move or a delete.
  connect(segmentsTab, SIGNAL(dataSegmentsRenumbered(QVector<int>)), this, SLOT(followDataSegments(QVector<int>)));
  connect(segmentsTab, SIGNAL(testSegmentsRenumbered(QVector<int>)), this, SLOT(followTestSegments(QVector<int>)));

  runTab = new RunTab();
  runTab->setPairsModel(pairsTab->getPairsModel());
  connect(runTab->calcButton, SIGNAL(clicked()), this, SLOT(SaveAndRun()));

#ifdef USE_QWT
  plotTab = new PlotTab(config, segmentsTab->getSegmentsDataModel(), segmentsTab->getSegmentsTestModel(), levelsTab->getLevelsModel());
  plotTab->setPairsModel(pairsTab->getPairsModel());
#endif

#ifdef USE_MCMC
  mcmcTab = new MCMCTab();
  mcmcTab->setTabReferences(levelsTab, segmentsTab);
  connect(mcmcTab->runButton, SIGNAL(clicked()), this, SLOT(SaveAndRunMCMC()));
#endif

  tabWidget->addTab(pairsTab, tr("&Particle Pairs"));
  tabWidget->addTab(levelsTab, tr("&Levels and Channels"));
  tabWidget->addTab(segmentsTab, tr("&Segments"));
  tabWidget->addTab(targetIntTab, tr("&Experimental Effects"));
  nuclearPotentialTabIndex = tabWidget->addTab(nuclearPotentialTab, tr("&Nuclear Potential"));
  tabWidget->addTab(fittingTab, tr("&Fitting Settings"));
  tabWidget->addTab(runTab, tr("&Calculate"));
#ifdef USE_MCMC
  tabWidget->addTab(mcmcTab, tr("&MCMC"));
#endif
#ifdef USE_QWT
  tabWidget->addTab(plotTab, tr("Pl&ot"));
#endif

  // Initially hide Nuclear Potential tab if hybrid method is disabled
  updateNuclearPotentialTabVisibility();

  setCentralWidget(tabWidget);

  createActions();
  createMenus();
  setThmEnabled(false);


  setWindowTitle(tr("AZURE2 -- untitled"));
}

Config &AZURESetup::GetConfig() {
  return config;
}

void AZURESetup::createActions() {
  aboutAction = new QAction(tr("&About AZURE2..."), this);
  connect(aboutAction, SIGNAL(triggered()), this, SLOT(showAbout()));

  resetAction = new QAction(tr("&New Project"), this);
  resetAction->setShortcuts(QKeySequence::New);
  connect(resetAction, SIGNAL(triggered()), this, SLOT(reset()));

  quitAction = new QAction(tr("&Quit"), this);
  quitAction->setShortcuts(QKeySequence::Quit);
  connect(quitAction, SIGNAL(triggered()), this, SLOT(close()));

  openAction = new QAction(tr("&Open..."), this);
  openAction->setShortcuts(QKeySequence::Open);
  connect(openAction, SIGNAL(triggered()), this, SLOT(open()));

  saveAction = new QAction(tr("&Save"), this);
  saveAction->setShortcuts(QKeySequence::Save);
  connect(saveAction, SIGNAL(triggered()), this, SLOT(save()));

  saveAsAction = new QAction(tr("Save &As..."), this);
  saveAsAction->setShortcuts(QKeySequence::SaveAs);
  connect(saveAsAction, SIGNAL(triggered()), this, SLOT(saveAs()));

  for (int i = 0; i < numRecent; i++) {
    recentFileActions[i] = new QAction(this);
    recentFileActions[i]->setVisible(false);
    connect(recentFileActions[i], SIGNAL(triggered()), this, SLOT(openRecent()));
  }
  clearRecentAction = new QAction(tr("&Clear"), this);
  clearRecentAction->setVisible(false);
  connect(clearRecentAction, SIGNAL(triggered()), this, SLOT(clearRecent()));

  copyAction = new QAction(tr("&Copy"), this);
  copyAction->setShortcuts(QKeySequence::Copy);

  matrixActionGroup = new QActionGroup(this);
  aMatrixAction = new QAction(tr("&A-Matrix"), this);
  aMatrixAction->setCheckable(true);
  matrixActionGroup->addAction(aMatrixAction);
  rMatrixAction = new QAction(tr("&R-Matrix"), this);
  rMatrixAction->setCheckable(true);
  matrixActionGroup->addAction(rMatrixAction);
  aMatrixAction->setChecked(true);
  connect(matrixActionGroup, SIGNAL(triggered(QAction *)), this, SLOT(matrixChanged(QAction *)));

  editChecksAction = new QAction(tr("&Checks..."), this);
  connect(editChecksAction, SIGNAL(triggered()), this, SLOT(editChecks()));

  editDirsAction = new QAction(tr("&Directories..."), this);
  connect(editDirsAction, SIGNAL(triggered()), this, SLOT(editDirs()));

  editOptionsAction = new QAction(tr("&Runtime Options..."), this);
  connect(editOptionsAction, SIGNAL(triggered()), this, SLOT(editOptions()));

  thmWorkspaceAction = new QAction(tr("&THM Workspace..."), this);
  thmWorkspaceAction->setToolTip(
      tr("The Trojan Horse (HOES) observable: model options, experiments, binding energies"));
  connect(thmWorkspaceAction, SIGNAL(triggered()), this, SLOT(editThmWorkspace()));

  buildCoulombTableAction = new QAction(tr("Build &Coulomb Table..."), this);
  buildCoulombTableAction->setStatusTip(tr("Build the table AZURE2 takes charged-particle Coulomb functions from"));
  connect(buildCoulombTableAction, SIGNAL(triggered()), this, SLOT(buildCoulombTable()));

  showTabInfoAction = new QAction(tr("Show Documentation For Current Tab"), this);
  showTabInfoAction->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_D));
  connect(showTabInfoAction, SIGNAL(triggered()), this, SLOT(showTabInfo()));

  openAZURESiteAction = new QAction(tr("Open AZURE Website..."), this);
  connect(openAZURESiteAction, SIGNAL(triggered()), this, SLOT(openWebsite()));
}

void AZURESetup::createMenus() {
  fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(aboutAction);
  fileMenu->addSeparator();
  fileMenu->addAction(resetAction);
  fileMenu->addAction(openAction);
  recentFileMenu = fileMenu->addMenu(tr("Open &Recent..."));
  for (int i = 0; i < numRecent; i++) recentFileMenu->addAction(recentFileActions[i]);
  recentSeparator = recentFileMenu->addSeparator();
  recentFileMenu->addAction(clearRecentAction);
  updateRecent();
  fileMenu->addAction(saveAction);
  fileMenu->addAction(saveAsAction);
  fileMenu->addSeparator();
  fileMenu->addAction(quitAction);

  configMenu = menuBar()->addMenu(tr("Co&nfigure"));
  formalismMenu = configMenu->addMenu(tr("&Formalism"));
  formalismMenu->addAction(aMatrixAction);
  formalismMenu->addAction(rMatrixAction);
  configMenu->addAction(editChecksAction);
  configMenu->addAction(editDirsAction);
  configMenu->addAction(editOptionsAction);
  configMenu->addAction(thmWorkspaceAction);
  configMenu->addAction(buildCoulombTableAction);

  helpMenu = menuBar()->addMenu(tr("&Documentation"));
  helpMenu->addAction(showTabInfoAction);
  helpMenu->addAction(openAZURESiteAction);
}

void AZURESetup::updateRecent() {
  QSettings settings;
  QStringList files = settings.value("recentFileList").toStringList();

  int numFiles = qMin(files.size(), (int)numRecent);

  for (int i = 0; i < numFiles; i++) {
    recentFileActions[i]->setText(tr("&%1 %2").arg(i + 1).arg(QFileInfo(files[i]).fileName()));
    recentFileActions[i]->setData(files[i]);
    recentFileActions[i]->setVisible(true);
  }

  for (int i = numFiles; i < numRecent; i++) recentFileActions[i]->setVisible(false);
  recentSeparator->setVisible(numFiles > 0);
  clearRecentAction->setVisible(numFiles > 0);
}

void AZURESetup::open() {
  QString filename = QFileDialog::getOpenFileName(this);
  if (!filename.isEmpty()) {
    if (!this->readFile(filename)) {
      reset();
      QMessageBox::information(this,
                               tr("Can't Access File"),
                               tr("An error was encountered while reading the file."));
    }
  }
}

void AZURESetup::open(QString filename) {
  if (!filename.isEmpty()) {
    if (!this->readFile(filename)) {
      reset();
      QMessageBox::information(this,
                               tr("Can't Access File"),
                               tr("An error was encountered while reading the file."));
    }
  }
}

void AZURESetup::openRecent() {
  QString filename = qobject_cast<QAction *>(sender())->data().toString();
  open(filename);
}

void AZURESetup::clearRecent() {
  QSettings settings;
  QStringList files = settings.value("recentFileList").toStringList();
  files.clear();
  settings.setValue("recentFileList", files);
  updateRecent();
}

bool AZURESetup::readFile(QString filename) {
  QFile file(filename);
  if (!file.open(QIODevice::ReadOnly)) return false;
  QFileInfo info(file);
  QString directory = info.absolutePath();

  reset();

  QTextStream in(&file);
  QString line("");

  while (line.trimmed() != QString("<config>") && !in.atEnd()) line = in.readLine();
  if (in.atEnd()) return false;
  if (!this->readConfig(in)) return false;

  // Read potential section if available (after config, before levels)
  line = QString("");
  while (line.trimmed() != QString("<potential>") && line.trimmed() != QString("<levels>") && !in.atEnd()) line = in.readLine();
  if (!in.atEnd() && line.trimmed() == QString("<potential>")) {
    if (!nuclearPotentialTab->readPotentialSettings(in, GetConfig())) {
      // Continue even if reading fails (old format without potential section)
    }
    // Look for levels section
    line = QString("");
    while (line.trimmed() != QString("<levels>") && !in.atEnd()) line = in.readLine();
  }
  if (in.atEnd()) return false;
  if (!levelsTab->readNuclearFile(in)) return false;
  levelsTab->setIncludeAllJPi(fileIncludeAllJPi);
  if (fileMaxL > levelsTab->maxL()) levelsTab->setMaxL(fileMaxL);
  GetConfig().includeAllJPi = fileIncludeAllJPi;
  GetConfig().hardSphereMaxL = levelsTab->maxL();
  pairsTab->setAbundanceSumSigma(fileAbundanceSumSigma);
  GetConfig().abundanceSumSigma = fileAbundanceSumSigma;

  // The optional <abundances> block sits between </levels> and <segmentsData>.
  line = QString("");
  while (line.trimmed() != QString("<abundances>") && line.trimmed() != QString("<segmentsData>") && !in.atEnd()) line = in.readLine();
  if (!in.atEnd() && line.trimmed() == QString("<abundances>")) {
    if (!pairsTab->readAbundances(in)) return false;
    line = QString("");
    while (line.trimmed() != QString("<segmentsData>") && !in.atEnd()) line = in.readLine();
  }
  if (in.atEnd()) return false;
  if (!segmentsTab->readSegDataFile(in)) return false;

  line = QString("");
  while (line.trimmed() != QString("<segmentsTest>") && !in.atEnd()) line = in.readLine();
  if (in.atEnd()) return false;
  if (!segmentsTab->readSegTestFile(in)) return false;

  line = QString("");
  while (line.trimmed() != QString("<targetInt>") && !in.atEnd()) line = in.readLine();
  if (in.atEnd()) return false;
  if (!targetIntTab->readFile(in)) return false;

  line = QString("");
  while (line.trimmed() != QString("<parameterSettings>") && !in.atEnd()) line = in.readLine();
  if (!in.atEnd() && line.trimmed() == QString("<parameterSettings>")) {
    if (!fittingTab->readParameterSettings(in)) return false;
  }

  line = QString("");
  while (line.trimmed() != QString("<lastRun>") && !in.atEnd()) line = in.readLine();
  if (!in.atEnd())
    if (!this->readLastRun(in)) return false;

#ifdef USE_MCMC
  // Try to read MCMC settings if available
  line = QString("");
  while (line.trimmed() != QString("<mcmc>") && !in.atEnd()) line = in.readLine();
  if (!in.atEnd() && line.trimmed() == QString("<mcmc>")) {
    if (!mcmcTab->readMCMCSettings(in)) return false;
  }
#endif

  // Update tab visibility based on useHybridMethod setting (was read in readConfig)
  updateNuclearPotentialTabVisibility();

  file.close();

  QFile file2(filename);
  if (!file2.open(QIODevice::ReadOnly)) return false;
  QTextStream in2(&file2);
  line = QString("");
  while (line.trimmed() != QString("<externalCapture>") && !in2.atEnd()) line = in2.readLine();
  if (!in2.atEnd()) {
    if (!pairsTab->parseOldECSection(in2)) return false;
  }
  file2.close();

  if (!readThmContent(filename)) return false;

  GetConfig().configfile = QDir::fromNativeSeparators(info.absoluteFilePath()).toStdString();
  setWindowTitle(QString("AZURE2 -- %1").arg(QString::fromStdString(GetConfig().configfile)));
  QDir::setCurrent(directory);

  // Populate FittingTab with current GUI state after file is loaded
  fittingTab->populateFromCurrentGUIState();

  QSettings settings;
  QStringList files = settings.value("recentFileList").toStringList();
  QString fullFileName = QDir::fromNativeSeparators(info.absoluteFilePath());
  files.removeAll(fullFileName);
  files.prepend(fullFileName);
  while (files.size() > numRecent) files.removeLast();

  settings.setValue("recentFileList", files);
  updateRecent();

  return true;
}

QString AZURESetup::projectDirectory() {
  QString file = QString::fromStdString(GetConfig().configfile);
  return file.isEmpty() ? QDir::currentPath() : QFileInfo(file).absolutePath();
}

namespace {
// An Experimental Effects segment list ("1,3-5", TargetEffect::GetSegmentsList)
// after a renumbering; the text unchanged when its set is, "" when no segment
// is left, and `ok` false (text unchanged) for a list this does not read.
QString renumberedSegmentList(const QString &text, const QVector<int> &newNumber, bool &ok) {
  ok = false;
  QList<int> before;
  for (const QString &part : text.split(',', Qt::SkipEmptyParts)) {
    const QStringList ends = part.split('-');
    bool a = false, b = true;
    const int lo = ends[0].trimmed().toInt(&a);
    const int hi = ends.size() == 2 ? ends[1].trimmed().toInt(&b) : lo;
    if (!a || !b || ends.size() > 2 || lo < 1 || hi < lo || hi - lo > 100000) return text;
    for (int k = lo; k <= hi; k++) before << k;
  }
  ok = true;
  QList<int> after;
  for (int k : before) {
    const int n = k <= newNumber.size() ? newNumber[k - 1] : k;
    if (n > 0 && !after.contains(n)) after << n;
  }
  std::sort(after.begin(), after.end());
  std::sort(before.begin(), before.end());
  if (after == before) return text;
  return ThmExperimentRecord::segmentsListText(after);
}
}  // namespace

void AZURESetup::followDataSegments(const QVector<int> &newNumber) {
  QStringList notes = followThmSegments(newNumber, false);
  TargetIntModel *effects = targetIntTab->getTargetIntModel();
  const QList<TargetIntData> lines = effects->getLines();
  int dropped = 0;
  for (int i = lines.size() - 1; i >= 0; i--) {
    bool ok;
    const QString list = renumberedSegmentList(lines[i].segmentsList, newNumber, ok);
    if (!ok || list == lines[i].segmentsList) continue;
    if (list.isEmpty()) {
      effects->removeRows(i, 1);
      dropped++;
    } else {
      effects->setData(effects->index(i, 1), list, Qt::EditRole);
    }
  }
  if (dropped) notes << tr("%n experimental effect line(s) removed: no segments left.", "", dropped);
  fittingTab->followSegments(newNumber);
  if (!notes.isEmpty()) statusBar()->showMessage(notes.join(' '));
}

bool AZURESetup::readLastRun(QTextStream &inStream) {
  unsigned int paramMask;
  unsigned int useTempFile;
  unsigned int rateEntrancePair;
  unsigned int rateExitPair;
  QString paramFile;
  QString integralsFile;
  QString temperatureFile;
  QString dummyString;
  double minTemp;
  double maxTemp;
  double tempStep;

  inStream >> paramMask;
  dummyString = inStream.readLine();
  dummyString = inStream.readLine();
  paramFile = dummyString.trimmed();
  dummyString = inStream.readLine();
  integralsFile = dummyString.trimmed();
  inStream >> rateEntrancePair >> rateExitPair;
  dummyString = inStream.readLine();
  inStream >> useTempFile;
  dummyString = inStream.readLine();
  temperatureFile = dummyString.trimmed();
  inStream >> minTemp >> maxTemp >> tempStep;

  QString line("");
  while (line.trimmed() != QString("</lastRun>") && !inStream.atEnd())
    line = inStream.readLine();
  if (line.trimmed() != QString("</lastRun>")) return false;

  if (paramMask & Config::CALCULATE_WITH_DATA) {
    if (paramMask & Config::PERFORM_FIT)
      runTab->calcType->setCurrentIndex(1);
    else if (paramMask & Config::PERFORM_ERROR_ANALYSIS)
      runTab->calcType->setCurrentIndex(3);
    else
      runTab->calcType->setCurrentIndex(0);
  } else {
    if (paramMask & Config::CALCULATE_REACTION_RATE)
      runTab->calcType->setCurrentIndex(4);
    else
      runTab->calcType->setCurrentIndex(2);
  }

  // Restore the uncertainty-band checkboxes (after calcType, so their enabled
  // state is already set).  Band first, then scaling (which depends on it).
  runTab->uncertaintyBandCheck->setChecked(paramMask & Config::CALCULATE_COVARIANCE_BAND);
  runTab->scaleCovarianceCheck->setChecked(paramMask & Config::SCALE_COVARIANCE_BY_CHI2);
  runTab->wignerLimitsCheck->setChecked(paramMask & Config::USE_WIGNER_LIMITS);

  // Set minimizer selection (0 Minuit2, 1 Minuit2+analytic grad, 2 Levenberg-
  // Marquardt, 3 GSL trust-region, 4+ NLopt).
#ifdef USE_NLOPT
  if (paramMask & Config::USE_NLOPT_MINIMIZER) {
    runTab->minimizerType->setCurrentIndex(GetConfig().nloptAlgorithm + 4);
  } else
#endif
      if (paramMask & Config::USE_GSL_LM_MINIMIZER) {
    runTab->minimizerType->setCurrentIndex(3);
  } else if (paramMask & Config::USE_LM_MINIMIZER) {
    runTab->minimizerType->setCurrentIndex(2);
  } else if (paramMask & Config::USE_ANALYTIC_GRADIENT) {
    runTab->minimizerType->setCurrentIndex(1);
  } else {
    runTab->minimizerType->setCurrentIndex(0);  // Minuit2
  }

  if (paramMask & Config::USE_GSL_COULOMB_FUNC)
    GetConfig().paramMask |= Config::USE_GSL_COULOMB_FUNC;
  else
    GetConfig().paramMask &= ~Config::USE_GSL_COULOMB_FUNC;

  if (paramMask & Config::USE_BRUNE_FORMALISM)
    GetConfig().paramMask |= Config::USE_BRUNE_FORMALISM;
  else
    GetConfig().paramMask &= ~Config::USE_BRUNE_FORMALISM;

  if (paramMask & Config::USE_PARK_FORMALISM)
    GetConfig().paramMask |= (Config::USE_PARK_FORMALISM | Config::USE_BRUNE_FORMALISM);
  else
    GetConfig().paramMask &= ~Config::USE_PARK_FORMALISM;

  if (paramMask & Config::IGNORE_ZERO_WIDTHS)
    GetConfig().paramMask |= Config::IGNORE_ZERO_WIDTHS;
  else
    GetConfig().paramMask &= ~Config::IGNORE_ZERO_WIDTHS;

  if (paramMask & Config::USE_RMC_FORMALISM)
    GetConfig().paramMask |= Config::USE_RMC_FORMALISM;
  else
    GetConfig().paramMask &= ~Config::USE_RMC_FORMALISM;

  if (paramMask & Config::TRANSFORM_PARAMETERS)
    GetConfig().paramMask |= Config::TRANSFORM_PARAMETERS;
  else
    GetConfig().paramMask &= ~Config::TRANSFORM_PARAMETERS;

  if (paramMask & Config::USE_LONGWAVELENGTH_APPROX)
    GetConfig().paramMask |= Config::USE_LONGWAVELENGTH_APPROX;
  else
    GetConfig().paramMask &= ~Config::USE_LONGWAVELENGTH_APPROX;

  if (paramMask & Config::USE_WIGNER_LIMITS)
    GetConfig().paramMask |= Config::USE_WIGNER_LIMITS;
  else
    GetConfig().paramMask &= ~Config::USE_WIGNER_LIMITS;

  if (paramMask & Config::SCALE_COVARIANCE_BY_CHI2)
    GetConfig().paramMask |= Config::SCALE_COVARIANCE_BY_CHI2;
  else
    GetConfig().paramMask &= ~Config::SCALE_COVARIANCE_BY_CHI2;

  runTab->rateEntranceKey->setKey(rateEntrancePair);
  runTab->rateExitKey->setKey(rateExitPair);

  if (minTemp != -1.) runTab->minTempText->setText(QString("%1").arg(minTemp));
  if (maxTemp != -1.) runTab->maxTempText->setText(QString("%1").arg(maxTemp));
  if (tempStep != -1.) runTab->tempStepText->setText(QString("%1").arg(tempStep));

  if (useTempFile == 1)
    runTab->fileTempButton->setChecked(true);
  else
    runTab->gridTempButton->setChecked(true);
  if (temperatureFile[0] == QChar('"')) temperatureFile.remove(0, 1);
  if (temperatureFile[temperatureFile.length() - 1] == QChar('"')) temperatureFile.remove(temperatureFile.length() - 1, 1);
  if (!temperatureFile.trimmed().isEmpty()) runTab->fileTempText->setText(temperatureFile.trimmed());

  if (paramMask & Config::USE_PREVIOUS_PARAMETERS)
    runTab->oldParamFileButton->setChecked(true);
  else
    runTab->newParamFileButton->setChecked(true);
  if (paramFile[0] == QChar('"')) paramFile.remove(0, 1);
  if (paramFile[paramFile.length() - 1] == QChar('"')) paramFile.remove(paramFile.length() - 1, 1);
  if (!paramFile.trimmed().isEmpty()) runTab->paramFileText->setText(paramFile.trimmed());

  if (paramMask & Config::USE_PREVIOUS_INTEGRALS)
    runTab->oldIntegralsFileButton->setChecked(true);
  else
    runTab->newIntegralsFileButton->setChecked(true);
  if (integralsFile[0] == QChar('"')) integralsFile.remove(0, 1);
  if (integralsFile[integralsFile.length() - 1] == QChar('"')) integralsFile.remove(integralsFile.length() - 1, 1);
  if (!integralsFile.trimmed().isEmpty()) runTab->integralsFileText->setText(integralsFile.trimmed());

  return true;
}

bool AZURESetup::readConfig(QTextStream &inStream) {
  QString isAMatrix;
  QString outputDirectory;
  QString checksDirectory;
  QString compoundCheck;
  QString boundaryCheck;
  QString dataCheck;
  QString lMatrixCheck;
  QString legendreCheck;
  QString coulAmpCheck;
  QString pathwaysCheck;
  QString angDistsCheck;
  QString dummyString;

  inStream >> isAMatrix;
  dummyString = inStream.readLine();
  dummyString = inStream.readLine();
  int poundSignPos = dummyString.lastIndexOf('#');
  if (poundSignPos == -1)
    outputDirectory = dummyString.trimmed();
  else
    outputDirectory = dummyString.left(poundSignPos).trimmed();
  dummyString = inStream.readLine();
  poundSignPos = dummyString.lastIndexOf('#');
  if (poundSignPos == -1)
    checksDirectory = dummyString.trimmed();
  else
    checksDirectory = dummyString.left(poundSignPos).trimmed();
  inStream >> compoundCheck;
  dummyString = inStream.readLine();
  inStream >> boundaryCheck;
  dummyString = inStream.readLine();
  inStream >> dataCheck;
  dummyString = inStream.readLine();
  inStream >> lMatrixCheck;
  dummyString = inStream.readLine();
  inStream >> legendreCheck;
  dummyString = inStream.readLine();
  inStream >> coulAmpCheck;
  dummyString = inStream.readLine();
  inStream >> pathwaysCheck;
  dummyString = inStream.readLine();
  inStream >> angDistsCheck;
  dummyString = inStream.readLine();

  // Optional key=value lines after the fixed ones; absent keys keep the defaults of an older file.
  fileIncludeAllJPi = false;
  fileMaxL = -1;
  fileAbundanceSumSigma = 0.0;
  QString line("");
  while (line.trimmed() != QString("</config>") && !inStream.atEnd()) {
    line = inStream.readLine();
    QString entry = line.section('#', 0, 0);
    int eq = entry.indexOf('=');
    if (eq < 0) continue;
    QString key = entry.left(eq).trimmed(), value = entry.mid(eq + 1).trimmed();
    if (key == "includeAllJPi")
      fileIncludeAllJPi = (value == "1" || value == "true");
    else if (key == "maxOrbitalMomentum")
      fileMaxL = value.toInt();
    else if (key == "abundanceSumSigma")
      fileAbundanceSumSigma = value.toDouble();
  }
  if (line.trimmed() != QString("</config>")) return false;

  if (isAMatrix == "false")
    rMatrixAction->activate(QAction::Trigger);
  else
    aMatrixAction->activate(QAction::Trigger);
  GetConfig().outputdir = outputDirectory.toStdString();
  GetConfig().checkdir = checksDirectory.toStdString();
  if (compoundCheck == "file")
    GetConfig().fileCheckMask |= Config::CHECK_COMPOUND_NUCLEUS;
  else if (compoundCheck == "screen")
    GetConfig().screenCheckMask |= Config::CHECK_COMPOUND_NUCLEUS;
  if (boundaryCheck == "file")
    GetConfig().fileCheckMask |= Config::CHECK_BOUNDARY_CONDITIONS;
  else if (boundaryCheck == "screen")
    GetConfig().screenCheckMask |= Config::CHECK_BOUNDARY_CONDITIONS;
  if (dataCheck == "file")
    GetConfig().fileCheckMask |= Config::CHECK_DATA;
  else if (dataCheck == "screen")
    GetConfig().screenCheckMask |= Config::CHECK_DATA;
  if (lMatrixCheck == "file")
    GetConfig().fileCheckMask |= Config::CHECK_ENERGY_DEP;
  else if (lMatrixCheck == "screen")
    GetConfig().screenCheckMask |= Config::CHECK_ENERGY_DEP;
  if (legendreCheck == "file")
    GetConfig().fileCheckMask |= Config::CHECK_LEGENDRE;
  else if (legendreCheck == "screen")
    GetConfig().screenCheckMask |= Config::CHECK_LEGENDRE;
  if (coulAmpCheck == "file")
    GetConfig().fileCheckMask |= Config::CHECK_COUL_AMPLITUDES;
  else if (coulAmpCheck == "screen")
    GetConfig().screenCheckMask |= Config::CHECK_COUL_AMPLITUDES;
  if (pathwaysCheck == "file")
    GetConfig().fileCheckMask |= Config::CHECK_PATHWAYS;
  else if (pathwaysCheck == "screen")
    GetConfig().screenCheckMask |= Config::CHECK_PATHWAYS;
  if (angDistsCheck == "file")
    GetConfig().fileCheckMask |= Config::CHECK_ANGULAR_DISTS;
  else if (angDistsCheck == "screen")
    GetConfig().screenCheckMask |= Config::CHECK_ANGULAR_DISTS;

  return true;
}

void AZURESetup::save() {
  if (!GetConfig().configfile.empty()) {
    if (!this->writeFile(QString::fromStdString(GetConfig().configfile)))
      QMessageBox::information(this,
                               tr("Can't Access File"),
                               tr("An error occured while writing the file."));
  } else
    saveAs();
}

void AZURESetup::saveAs() {
  QString filename = QFileDialog::getSaveFileName(this);
  if (!filename.isEmpty()) {
    if (!this->writeFile(filename))
      QMessageBox::information(this,
                               tr("Can't Access File"),
                               tr("An error occured while writing the file."));
    else {
      QSettings settings;
      QStringList files = settings.value("recentFileList").toStringList();
      QFile file(filename);
      QFileInfo info(file);
      QString fullFileName = QDir::fromNativeSeparators(info.absoluteFilePath());
      files.removeAll(fullFileName);
      files.prepend(fullFileName);
      while (files.size() > numRecent) files.removeLast();
      settings.setValue("recentFileList", files);
      updateRecent();
    }
  }
}

bool AZURESetup::writeFile(QString filename) {
  QFile file(filename);
  if (!file.open(QIODevice::WriteOnly)) return false;
  QFileInfo info(file);
  QString directory = info.absolutePath();

  QTextStream out(&file);
  if (!this->writeProject(out, directory)) return false;

  GetConfig().configfile = QDir::fromNativeSeparators(info.absoluteFilePath()).toStdString();
  setWindowTitle(QString("AZURE2 -- %1").arg(QString::fromStdString(GetConfig().configfile)));
  QDir::setCurrent(directory);

  out.flush();
  file.close();
  return true;
}

bool AZURESetup::writeProject(QTextStream &out, const QString &directory) {
  out << "<config>" << Qt::endl;
  if (!this->writeConfig(out, directory)) return false;
  // The hard-sphere option and its limit live in the file so that --no-gui and pyazr runs agree with the GUI.
  GetConfig().includeAllJPi = levelsTab->includeAllJPi();
  GetConfig().hardSphereMaxL = levelsTab->maxL();
  out << "includeAllJPi=" << (levelsTab->includeAllJPi() ? 1 : 0) << " #Hard-sphere scattering in every J^pi up to maxOrbitalMomentum" << Qt::endl;
  out << "maxOrbitalMomentum=" << levelsTab->maxL() << Qt::endl;
  // Only a project with a sum-to-one prior on its abundances carries the key, so other files are unchanged.
  GetConfig().abundanceSumSigma = pairsTab->abundanceSumSigma();
  if (pairsTab->abundanceSumSigma() > 0.0)
    out << "abundanceSumSigma=" << pairsTab->abundanceSumSigma() << " #Sum-to-one prior on the isotopic abundances" << Qt::endl;
  out << "</config>" << Qt::endl;

  // Write potential section after config.  A pair may be switched on while the
  // default is off, so the master switch has to follow the pairs; otherwise the
  // file would come back with the whole hybrid model disabled.
  if (NuclearPotentialManager::instance().isAnyEnabled()) GetConfig().useHybridMethod = true;
  out << "<potential>" << Qt::endl;
  out << "useHybridPotential="
      << (NuclearPotentialManager::instance().getDefaultEnabled() ? "1" : "0") << Qt::endl;
  out << "useAdaptiveGrid=" << (GetConfig().useAdaptiveGrid ? "1" : "0") << Qt::endl;
  nuclearPotentialTab->writePotentialSettings(out);
  out << "</potential>" << Qt::endl;

  out << "<levels>" << Qt::endl;
  if (!levelsTab->writeNuclearFile(out)) return false;
  out << "</levels>" << Qt::endl;

  // Isotopic abundances, one line per compound nucleus, only for a project
  // that needs them (several nuclei, or a non-default entry).
  if (pairsTab->hasAbundanceBlock()) {
    out << "<abundances>" << Qt::endl;
    if (!pairsTab->writeAbundances(out)) return false;
    out << "</abundances>" << Qt::endl;
  }

  out << "<segmentsData>" << Qt::endl;
  if (!segmentsTab->writeSegDataFile(out)) return false;
  out << "</segmentsData>" << Qt::endl;

  out << "<segmentsTest>" << Qt::endl;
  if (!segmentsTab->writeSegTestFile(out)) return false;
  out << "</segmentsTest>" << Qt::endl;

  out << "<targetInt>" << Qt::endl;
  if (!targetIntTab->writeFile(out)) return false;
  out << "</targetInt>" << Qt::endl;

  // Only if the file that was read had one: a project without the block must
  // come back without it.
  if (hasThmBlock) writeThmBlock(out, thmBlockLines);

  out << "<parameterSettings>" << Qt::endl;
  if (!fittingTab->writeParameterSettings(out)) return false;
  out << "</parameterSettings>" << Qt::endl;

  out << "<lastRun>" << Qt::endl;
  if (!writeLastRun(out)) return false;
  out << "</lastRun>" << Qt::endl;

#ifdef USE_MCMC
  if (!mcmcTab->writeMCMCSettings(out)) return false;
#endif
  return true;
}

bool AZURESetup::writeConfig(QTextStream &outStream, QString directory) {
  QString isAMatrix;
  QString outputDirectory;
  QString checksDirectory;
  QString compoundCheck;
  QString boundaryCheck;
  QString dataCheck;
  QString lMatrixCheck;
  QString legendreCheck;
  QString coulAmpCheck;
  QString pathwaysCheck;
  QString angDistsCheck;

  if (GetConfig().paramMask & Config::USE_AMATRIX)
    isAMatrix = "true";
  else
    isAMatrix = "false";
  bool emptyCheckDir = false;
  bool emptyOutputDir = false;
  if (!GetConfig().outputdir.empty())
    outputDirectory = QString::fromStdString(GetConfig().outputdir);
  else {
    outputDirectory = QDir::fromNativeSeparators(directory) + '/';
    GetConfig().outputdir = outputDirectory.toStdString();
    emptyOutputDir = true;
  }
  if (!GetConfig().checkdir.empty())
    checksDirectory = QString::fromStdString(GetConfig().checkdir);
  else {
    checksDirectory = QDir::fromNativeSeparators(directory) + '/';
    GetConfig().checkdir = checksDirectory.toStdString();
    emptyCheckDir = true;
  }
  if (emptyCheckDir && emptyOutputDir) {
    QMessageBox::information(this, tr("Unspecified Directories"),
                             QString("The output and checks directories are unspecified. "
                                     "They will be set to %1.")
                                 .arg(outputDirectory.trimmed()));
  } else if (emptyCheckDir) {
    QMessageBox::information(this, tr("Unspecified Directory"),
                             QString("The checks directory is unspecified. "
                                     "It will be set to %1.")
                                 .arg(checksDirectory.trimmed()));
  } else if (emptyOutputDir) {
    QMessageBox::information(this, tr("Unspecified Directory"),
                             QString("The output directory is unspecified. "
                                     "It will be set to %1.")
                                 .arg(outputDirectory.trimmed()));
  }
  if (GetConfig().fileCheckMask & Config::CHECK_COMPOUND_NUCLEUS)
    compoundCheck = "file";
  else if (GetConfig().screenCheckMask & Config::CHECK_COMPOUND_NUCLEUS)
    compoundCheck = "screen";
  else
    compoundCheck = "none";
  if (GetConfig().fileCheckMask & Config::CHECK_BOUNDARY_CONDITIONS)
    boundaryCheck = "file";
  else if (GetConfig().screenCheckMask & Config::CHECK_BOUNDARY_CONDITIONS)
    boundaryCheck = "screen";
  else
    boundaryCheck = "none";
  if (GetConfig().fileCheckMask & Config::CHECK_DATA)
    dataCheck = "file";
  else if (GetConfig().screenCheckMask & Config::CHECK_DATA)
    dataCheck = "screen";
  else
    dataCheck = "none";
  if (GetConfig().fileCheckMask & Config::CHECK_ENERGY_DEP)
    lMatrixCheck = "file";
  else if (GetConfig().screenCheckMask & Config::CHECK_ENERGY_DEP)
    lMatrixCheck = "screen";
  else
    lMatrixCheck = "none";
  if (GetConfig().fileCheckMask & Config::CHECK_LEGENDRE)
    legendreCheck = "file";
  else if (GetConfig().screenCheckMask & Config::CHECK_LEGENDRE)
    legendreCheck = "screen";
  else
    legendreCheck = "none";
  if (GetConfig().fileCheckMask & Config::CHECK_COUL_AMPLITUDES)
    coulAmpCheck = "file";
  else if (GetConfig().screenCheckMask & Config::CHECK_COUL_AMPLITUDES)
    coulAmpCheck = "screen";
  else
    coulAmpCheck = "none";
  if (GetConfig().fileCheckMask & Config::CHECK_PATHWAYS)
    pathwaysCheck = "file";
  else if (GetConfig().screenCheckMask & Config::CHECK_PATHWAYS)
    pathwaysCheck = "screen";
  else
    pathwaysCheck = "none";
  if (GetConfig().fileCheckMask & Config::CHECK_ANGULAR_DISTS)
    angDistsCheck = "file";
  else if (GetConfig().screenCheckMask & Config::CHECK_ANGULAR_DISTS)
    angDistsCheck = "screen";
  else
    angDistsCheck = "none";

  outStream.setFieldAlignment(QTextStream::AlignLeft);
  outStream << qSetFieldWidth(100) << isAMatrix << qSetFieldWidth(0) << "#Perform A-Matrix Calculation" << Qt::endl;
  outStream << qSetFieldWidth(100) << outputDirectory << qSetFieldWidth(0) << "#Full Path to Output Directory" << Qt::endl;
  outStream << qSetFieldWidth(100) << checksDirectory << qSetFieldWidth(0) << "#Full Path to Checks Directory" << Qt::endl;
  outStream << qSetFieldWidth(100) << compoundCheck << qSetFieldWidth(0) << "#Compond Nucleus Check" << Qt::endl;
  outStream << qSetFieldWidth(100) << boundaryCheck << qSetFieldWidth(0) << "#Boundary Condition Check" << Qt::endl;
  outStream << qSetFieldWidth(100) << dataCheck << qSetFieldWidth(0) << "#Data Check" << Qt::endl;
  outStream << qSetFieldWidth(100) << lMatrixCheck << qSetFieldWidth(0) << "#Lo-Matrix and Penetrability Check" << Qt::endl;
  outStream << qSetFieldWidth(100) << legendreCheck << qSetFieldWidth(0) << "#Legendre Polynomial Check" << Qt::endl;
  outStream << qSetFieldWidth(100) << coulAmpCheck << qSetFieldWidth(0) << "#Coulomb Amplitudes Check" << Qt::endl;
  outStream << qSetFieldWidth(100) << pathwaysCheck << qSetFieldWidth(0) << "#Reaction Pathway Check" << Qt::endl;
  outStream << qSetFieldWidth(100) << angDistsCheck << qSetFieldWidth(0) << "#Angular Distributions Check" << Qt::endl;

  return true;
}

bool AZURESetup::writeLastRun(QTextStream &outStream) {
  unsigned int paramMask = GetConfig().paramMask;

  if (runTab->calcType->currentIndex() == 1 ||
      runTab->calcType->currentIndex() == 3)
    paramMask |= Config::PERFORM_FIT;
  else
    paramMask &= ~Config::PERFORM_FIT;
  if (runTab->calcType->currentIndex() == 2 ||
      runTab->calcType->currentIndex() == 4)
    paramMask &= ~Config::CALCULATE_WITH_DATA;
  else
    paramMask |= Config::CALCULATE_WITH_DATA;
  if (runTab->calcType->currentIndex() == 3)
    paramMask |= Config::PERFORM_ERROR_ANALYSIS;
  else
    paramMask &= ~Config::PERFORM_ERROR_ANALYSIS;
  if (runTab->calcType->currentIndex() == 4)
    paramMask |= Config::CALCULATE_REACTION_RATE;
  else
    paramMask &= ~Config::CALCULATE_REACTION_RATE;
  // Analytic cross-section uncertainty band: user-selected via the Run-tab
  // checkboxes (enabled only for fit / extrapolation / MINOS modes).
  if (runTab->uncertaintyBandCheck->isChecked())
    paramMask |= Config::CALCULATE_COVARIANCE_BAND;
  else
    paramMask &= ~Config::CALCULATE_COVARIANCE_BAND;
  if (runTab->scaleCovarianceCheck->isChecked())
    paramMask |= Config::SCALE_COVARIANCE_BY_CHI2;
  else
    paramMask &= ~Config::SCALE_COVARIANCE_BY_CHI2;
  if (runTab->wignerLimitsCheck->isChecked())
    paramMask |= Config::USE_WIGNER_LIMITS;
  else
    paramMask &= ~Config::USE_WIGNER_LIMITS;

  if (runTab->oldParamFileButton->isChecked())
    paramMask |= Config::USE_PREVIOUS_PARAMETERS;
  else
    paramMask &= ~Config::USE_PREVIOUS_PARAMETERS;
  if (runTab->oldIntegralsFileButton->isChecked())
    paramMask |= Config::USE_PREVIOUS_INTEGRALS;
  else
    paramMask &= ~Config::USE_PREVIOUS_INTEGRALS;

  outStream << paramMask << Qt::endl;
  outStream << '"' << runTab->paramFileText->text() << '"' << Qt::endl;
  outStream << '"' << runTab->integralsFileText->text() << '"' << Qt::endl;
  outStream << runTab->rateEntranceKey->key() << ' ' << runTab->rateExitKey->key();
  outStream << Qt::endl;
  if (runTab->fileTempButton->isChecked())
    outStream << "1 ";
  else
    outStream << "0 ";
  outStream << '"' << runTab->fileTempText->text() << '"' << Qt::endl;
  if (!runTab->minTempText->text().isEmpty())
    outStream << runTab->minTempText->text() << ' ';
  else
    outStream << "-1. ";
  if (!runTab->maxTempText->text().isEmpty())
    outStream << runTab->maxTempText->text() << ' ';
  else
    outStream << "-1. ";
  if (!runTab->tempStepText->text().isEmpty())
    outStream << runTab->tempStepText->text();
  else
    outStream << "-1.";
  outStream << Qt::endl;

  return true;
}

void AZURESetup::matrixChanged(QAction *action) {
  if (action == aMatrixAction)
    GetConfig().paramMask |= Config::USE_AMATRIX;
  else {
    if (GetConfig().paramMask & Config::IGNORE_ZERO_WIDTHS) {
      QMessageBox::information(this, tr("Incompatible Option"),
                               tr("The option to ignore external widths is not possible for R-Matrix formalism. Remove option to use R-Matrix formalism."));
      aMatrixAction->activate(QAction::Trigger);
    } else
      GetConfig().paramMask &= ~Config::USE_AMATRIX;
  }
}

void AZURESetup::editChecks() {
  EditChecksDialog aDialog;
  if (GetConfig().fileCheckMask & Config::CHECK_COMPOUND_NUCLEUS)
    aDialog.compoundCheckCombo->setCurrentIndex(2);
  else if (GetConfig().screenCheckMask & Config::CHECK_COMPOUND_NUCLEUS)
    aDialog.compoundCheckCombo->setCurrentIndex(1);
  else
    aDialog.compoundCheckCombo->setCurrentIndex(0);

  if (GetConfig().fileCheckMask & Config::CHECK_BOUNDARY_CONDITIONS)
    aDialog.boundaryCheckCombo->setCurrentIndex(2);
  else if (GetConfig().screenCheckMask & Config::CHECK_BOUNDARY_CONDITIONS)
    aDialog.boundaryCheckCombo->setCurrentIndex(1);
  else
    aDialog.boundaryCheckCombo->setCurrentIndex(0);

  if (GetConfig().fileCheckMask & Config::CHECK_DATA)
    aDialog.dataCheckCombo->setCurrentIndex(2);
  else if (GetConfig().screenCheckMask & Config::CHECK_DATA)
    aDialog.dataCheckCombo->setCurrentIndex(1);
  else
    aDialog.dataCheckCombo->setCurrentIndex(0);

  if (GetConfig().fileCheckMask & Config::CHECK_ENERGY_DEP)
    aDialog.lMatrixCheckCombo->setCurrentIndex(2);
  else if (GetConfig().screenCheckMask & Config::CHECK_ENERGY_DEP)
    aDialog.lMatrixCheckCombo->setCurrentIndex(1);
  else
    aDialog.lMatrixCheckCombo->setCurrentIndex(0);

  if (GetConfig().fileCheckMask & Config::CHECK_LEGENDRE)
    aDialog.legendreCheckCombo->setCurrentIndex(2);
  else if (GetConfig().screenCheckMask & Config::CHECK_LEGENDRE)
    aDialog.legendreCheckCombo->setCurrentIndex(1);
  else
    aDialog.legendreCheckCombo->setCurrentIndex(0);

  if (GetConfig().fileCheckMask & Config::CHECK_COUL_AMPLITUDES)
    aDialog.coulAmpCheckCombo->setCurrentIndex(2);
  else if (GetConfig().screenCheckMask & Config::CHECK_COUL_AMPLITUDES)
    aDialog.coulAmpCheckCombo->setCurrentIndex(1);
  else
    aDialog.coulAmpCheckCombo->setCurrentIndex(0);

  if (GetConfig().fileCheckMask & Config::CHECK_PATHWAYS)
    aDialog.pathwaysCheckCombo->setCurrentIndex(2);
  else if (GetConfig().screenCheckMask & Config::CHECK_PATHWAYS)
    aDialog.pathwaysCheckCombo->setCurrentIndex(1);
  else
    aDialog.pathwaysCheckCombo->setCurrentIndex(0);

  if (GetConfig().fileCheckMask & Config::CHECK_ANGULAR_DISTS)
    aDialog.angDistsCheckCombo->setCurrentIndex(2);
  else if (GetConfig().screenCheckMask & Config::CHECK_ANGULAR_DISTS)
    aDialog.angDistsCheckCombo->setCurrentIndex(1);
  else
    aDialog.angDistsCheckCombo->setCurrentIndex(0);

  if (aDialog.exec()) {
    GetConfig().fileCheckMask = 0;
    GetConfig().screenCheckMask = 0;
    if (aDialog.compoundCheckCombo->currentIndex() == 1)
      GetConfig().screenCheckMask |= Config::CHECK_COMPOUND_NUCLEUS;
    else if (aDialog.compoundCheckCombo->currentIndex() == 2)
      GetConfig().fileCheckMask |= Config::CHECK_COMPOUND_NUCLEUS;
    if (aDialog.boundaryCheckCombo->currentIndex() == 1)
      GetConfig().screenCheckMask |= Config::CHECK_BOUNDARY_CONDITIONS;
    else if (aDialog.boundaryCheckCombo->currentIndex() == 2)
      GetConfig().fileCheckMask |= Config::CHECK_BOUNDARY_CONDITIONS;
    if (aDialog.dataCheckCombo->currentIndex() == 1)
      GetConfig().screenCheckMask |= Config::CHECK_DATA;
    else if (aDialog.dataCheckCombo->currentIndex() == 2)
      GetConfig().fileCheckMask |= Config::CHECK_DATA;
    if (aDialog.lMatrixCheckCombo->currentIndex() == 1)
      GetConfig().screenCheckMask |= Config::CHECK_ENERGY_DEP;
    else if (aDialog.lMatrixCheckCombo->currentIndex() == 2)
      GetConfig().fileCheckMask |= Config::CHECK_ENERGY_DEP;
    if (aDialog.legendreCheckCombo->currentIndex() == 1)
      GetConfig().screenCheckMask |= Config::CHECK_LEGENDRE;
    else if (aDialog.legendreCheckCombo->currentIndex() == 2)
      GetConfig().fileCheckMask |= Config::CHECK_LEGENDRE;
    if (aDialog.coulAmpCheckCombo->currentIndex() == 1)
      GetConfig().screenCheckMask |= Config::CHECK_COUL_AMPLITUDES;
    else if (aDialog.coulAmpCheckCombo->currentIndex() == 2)
      GetConfig().fileCheckMask |= Config::CHECK_COUL_AMPLITUDES;
    if (aDialog.pathwaysCheckCombo->currentIndex() == 1)
      GetConfig().screenCheckMask |= Config::CHECK_PATHWAYS;
    else if (aDialog.pathwaysCheckCombo->currentIndex() == 2)
      GetConfig().fileCheckMask |= Config::CHECK_PATHWAYS;
    if (aDialog.angDistsCheckCombo->currentIndex() == 1)
      GetConfig().screenCheckMask |= Config::CHECK_ANGULAR_DISTS;
    else if (aDialog.angDistsCheckCombo->currentIndex() == 2)
      GetConfig().fileCheckMask |= Config::CHECK_ANGULAR_DISTS;
  }
}

void AZURESetup::editDirs() {
  EditDirsDialog aDialog;

  aDialog.outputDirectoryText->setText(QString::fromStdString(GetConfig().outputdir));
  aDialog.checksDirectoryText->setText(QString::fromStdString(GetConfig().checkdir));

  if (aDialog.exec()) {
    GetConfig().outputdir = aDialog.outputDirectoryText->text().toStdString();
    GetConfig().checkdir = aDialog.checksDirectoryText->text().toStdString();
  }
}

// Configure > Build Coulomb Table...: the table is not distributed with AZURE2, so a user
// of a ready-built AZURE2 makes it here (the same build as AZURE2 --build-coulomb-table).
// The build runs on a worker thread; this loop keeps the window responsive and offers
// Cancel.  Calculations started afterwards use the new table (CoulombTable::ForgetMissing).
void AZURESetup::buildCoulombTable() {
  const QString path = QString::fromStdString(CoulombTable::DefaultPath());
  const QString intro =
      tr("AZURE2 takes the Coulomb functions of charged-particle channels from a table when it "
         "finds one. The results are the same and calculations are several times faster.");
  const QString question = QFileInfo::exists(path)
      ? tr("%1\n\nA table already exists at\n%2\n\nBuild it again?").arg(intro, path)
      : tr("%1\n\nThere is no table at\n%2\n\nBuild it now? It takes about 5 minutes on four "
           "threads and needs 88 MB of disk space.").arg(intro, path);
  if (QMessageBox::question(this, tr("Build Coulomb Table"), question, QMessageBox::Yes | QMessageBox::No,
                            QMessageBox::Yes) != QMessageBox::Yes)
    return;

  CoulombTable::BuildSpec spec;
  spec.threads = std::max(1u, std::thread::hardware_concurrency());
  std::atomic<int> done{0};
  std::atomic<bool> cancel{false}, finished{false};
  spec.rowsDone = &done;
  spec.cancel = &cancel;
  const int rows = CoulombTable::BuildRows(spec);

  QProgressDialog progress(tr("Building the Coulomb table on %1 thread(s)...").arg(spec.threads), tr("Cancel"), 0,
                           rows, this);
  progress.setWindowTitle(tr("Build Coulomb Table"));
  progress.setWindowModality(Qt::WindowModal);
  progress.setAutoClose(false);
  progress.setAutoReset(false);
  progress.setMinimumDuration(0);
  progress.setValue(0);

  bool ok = false;
  std::string error;
  const std::string file = path.toStdString();
  QElapsedTimer timer;
  timer.start();
  std::thread worker([&] {
    ok = CoulombTable::BuildFile(file, spec, error);
    finished = true;
  });
  while (!finished) {
    QCoreApplication::processEvents(QEventLoop::AllEvents, 100);
    progress.setValue(std::min(done.load(), rows - 1));
    if (progress.wasCanceled()) {
      cancel = true;
      progress.setLabelText(tr("Cancelling..."));
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
  }
  worker.join();
  progress.close();

  if (ok) {
    CoulombTable::ForgetMissing();
    QMessageBox::information(this, tr("Build Coulomb Table"),
                             tr("The Coulomb table was written to\n%1\nin %2 s. Calculations started from now on "
                                "use it.").arg(path).arg(timer.elapsed() / 1000));
  } else if (cancel) {
    QMessageBox::information(this, tr("Build Coulomb Table"), tr("Cancelled; no table was written."));
  } else {
    QMessageBox::warning(this, tr("Build Coulomb Table"),
                         tr("The Coulomb table was not built: %1").arg(QString::fromStdString(error)));
  }
}

void AZURESetup::editOptions() {
  EditOptionsDialog aDialog;

  if (GetConfig().paramMask & Config::USE_GSL_COULOMB_FUNC)
    aDialog.useGSLCoulCheck->setChecked(true);
  else
    aDialog.useGSLCoulCheck->setChecked(false);

  if (GetConfig().paramMask & Config::USE_BRUNE_FORMALISM)
    aDialog.useBruneCheck->setChecked(true);
  else
    aDialog.useBruneCheck->setChecked(false);

  if (GetConfig().paramMask & Config::USE_PARK_FORMALISM)
    aDialog.useParkCheck->setChecked(true);
  else
    aDialog.useParkCheck->setChecked(false);

  if (GetConfig().paramMask & Config::IGNORE_ZERO_WIDTHS)
    aDialog.ignoreExternalsCheck->setChecked(true);
  else
    aDialog.ignoreExternalsCheck->setChecked(false);

  if (GetConfig().paramMask & Config::USE_RMC_FORMALISM)
    aDialog.useRMCCheck->setChecked(true);
  else
    aDialog.useRMCCheck->setChecked(false);

  if (!(GetConfig().paramMask & Config::TRANSFORM_PARAMETERS))
    aDialog.noTransformCheck->setChecked(true);
  else
    aDialog.noTransformCheck->setChecked(false);

  if (GetConfig().useHybridMethod)
    aDialog.useHybridMethodCheck->setChecked(true);
  else
    aDialog.useHybridMethodCheck->setChecked(false);

  aDialog.useAdaptiveGridCheck->setChecked(GetConfig().useAdaptiveGrid);
  aDialog.useThmCheck->setChecked(thmEnabled());

  // if(!(GetConfig().paramMask & Config::USE_LONGWAVELENGTH_APPROX)) aDialog.noLongWavelengthCheck->setChecked(true);
  // else aDialog.noLongWavelengthCheck->setChecked(false);

  if (aDialog.exec()) {
    if (aDialog.useGSLCoulCheck->isChecked())
      GetConfig().paramMask |= Config::USE_GSL_COULOMB_FUNC;
    else
      GetConfig().paramMask &= ~Config::USE_GSL_COULOMB_FUNC;

    if (aDialog.useBruneCheck->isChecked())
      GetConfig().paramMask |= Config::USE_BRUNE_FORMALISM;
    else
      GetConfig().paramMask &= ~Config::USE_BRUNE_FORMALISM;

    if (aDialog.useParkCheck->isChecked())
      GetConfig().paramMask |= (Config::USE_PARK_FORMALISM | Config::USE_BRUNE_FORMALISM);
    else
      GetConfig().paramMask &= ~Config::USE_PARK_FORMALISM;

    if (aDialog.ignoreExternalsCheck->isChecked()) {
      GetConfig().paramMask |= Config::IGNORE_ZERO_WIDTHS;
      if (!(GetConfig().paramMask & Config::USE_AMATRIX)) {
        QMessageBox::information(this, tr("Incompatible Option"),
                                 tr("The option to ignore external widths is not possible for R-Matrix formalism.  The formalism will be changed to A-Matrix. "));
        aMatrixAction->activate(QAction::Trigger);
      }
    } else
      GetConfig().paramMask &= ~Config::IGNORE_ZERO_WIDTHS;

    if (aDialog.useRMCCheck->isChecked())
      GetConfig().paramMask |= Config::USE_RMC_FORMALISM;
    else
      GetConfig().paramMask &= ~Config::USE_RMC_FORMALISM;

    if (aDialog.noTransformCheck->isChecked())
      GetConfig().paramMask &= ~Config::TRANSFORM_PARAMETERS;
    else
      GetConfig().paramMask |= Config::TRANSFORM_PARAMETERS;

    if (aDialog.useHybridMethodCheck->isChecked()) {
      GetConfig().useHybridMethod = true;
    } else {
      GetConfig().useHybridMethod = false;
    }

    GetConfig().useAdaptiveGrid = aDialog.useAdaptiveGridCheck->isChecked();

    // Update tab visibility based on hybrid method setting
    updateNuclearPotentialTabVisibility();
    applyThmOption(aDialog.useThmCheck->isChecked());

    // if(aDialog.noLongWavelengthCheck->isChecked()) GetConfig().paramMask &= ~Config::USE_LONGWAVELENGTH_APPROX;
    // else GetConfig().paramMask |= Config::USE_LONGWAVELENGTH_APPROX;
  }
}

void AZURESetup::updateNuclearPotentialTabVisibility() {
  if (GetConfig().useHybridMethod) {
    // Make sure the tab is visible
    if (tabWidget->indexOf(nuclearPotentialTab) == -1) {
      // Tab was removed, insert it back at its original position
      // Insert before Fitting Settings tab
      int fittingTabIndex = tabWidget->indexOf(fittingTab);
      if (fittingTabIndex != -1) {
        tabWidget->insertTab(fittingTabIndex, nuclearPotentialTab, tr("&Nuclear Potential"));
        nuclearPotentialTabIndex = tabWidget->indexOf(nuclearPotentialTab);
      }
    }
  } else {
    // Hide the tab by removing it from the tab widget
    int index = tabWidget->indexOf(nuclearPotentialTab);
    if (index != -1) {
      tabWidget->removeTab(index);
    }
  }
}

void AZURESetup::SaveAndRun() {
  save();
  if (GetConfig().configfile.empty()) return;
  runTab->runtimeText->clear();
  QFile file(QString::fromStdString(GetConfig().configfile));
  QFileInfo info(file);
  QString directory = info.absolutePath();
  if (GetConfig().outputdir.empty()) GetConfig().outputdir = QDir::fromNativeSeparators(directory).toStdString() + '/';
  if (GetConfig().checkdir.empty()) GetConfig().checkdir = QDir::fromNativeSeparators(directory).toStdString() + '/';

  GetConfig().chiVariance = runTab->chiVarianceText->text().toDouble();

  if (runTab->calcType->currentIndex() == 1 ||
      runTab->calcType->currentIndex() == 3)
    GetConfig().paramMask |= Config::PERFORM_FIT;
  else
    GetConfig().paramMask &= ~Config::PERFORM_FIT;
  if (runTab->calcType->currentIndex() == 2 ||
      runTab->calcType->currentIndex() == 4)
    GetConfig().paramMask &= ~Config::CALCULATE_WITH_DATA;
  else
    GetConfig().paramMask |= Config::CALCULATE_WITH_DATA;
  if (runTab->calcType->currentIndex() == 3)
    GetConfig().paramMask |= Config::PERFORM_ERROR_ANALYSIS;
  else
    GetConfig().paramMask &= ~Config::PERFORM_ERROR_ANALYSIS;
  if (runTab->calcType->currentIndex() == 4)
    GetConfig().paramMask |= Config::CALCULATE_REACTION_RATE;
  else
    GetConfig().paramMask &= ~Config::CALCULATE_REACTION_RATE;
  if (runTab->uncertaintyBandCheck->isChecked())
    GetConfig().paramMask |= Config::CALCULATE_COVARIANCE_BAND;
  else
    GetConfig().paramMask &= ~Config::CALCULATE_COVARIANCE_BAND;
  if (runTab->scaleCovarianceCheck->isChecked())
    GetConfig().paramMask |= Config::SCALE_COVARIANCE_BY_CHI2;
  else
    GetConfig().paramMask &= ~Config::SCALE_COVARIANCE_BY_CHI2;
  if (runTab->wignerLimitsCheck->isChecked())
    GetConfig().paramMask |= Config::USE_WIGNER_LIMITS;
  else
    GetConfig().paramMask &= ~Config::USE_WIGNER_LIMITS;

  // Handle minimizer selection (only for fitting operations).  Combo layout:
  //   0 Minuit2, 1 Minuit2 + analytic gradient, 2 Levenberg-Marquardt,
  //   3 GSL trust-region (geodesic), 4+ NLopt.
  GetConfig().paramMask &= ~(Config::USE_NLOPT_MINIMIZER |
                             Config::USE_ANALYTIC_GRADIENT |
                             Config::USE_LM_MINIMIZER |
                             Config::USE_GSL_LM_MINIMIZER);
  if (runTab->calcType->currentIndex() == 1 || runTab->calcType->currentIndex() == 3) {
    int mIdx = runTab->minimizerType->currentIndex();
    if (mIdx == 1) {
      GetConfig().paramMask |= Config::USE_ANALYTIC_GRADIENT;
    } else if (mIdx == 2) {
      GetConfig().paramMask |= Config::USE_LM_MINIMIZER;
    } else if (mIdx == 3) {
      GetConfig().paramMask |= Config::USE_GSL_LM_MINIMIZER;
    }
#ifdef USE_NLOPT
    else if (mIdx >= 4) {
      GetConfig().paramMask |= Config::USE_NLOPT_MINIMIZER;
      GetConfig().nloptAlgorithm = mIdx - 4;
    }
#endif
  }

  if (runTab->oldParamFileButton->isChecked()) {
    GetConfig().paramMask |= Config::USE_PREVIOUS_PARAMETERS;
    GetConfig().paramfile = runTab->paramFileText->text().toStdString();
  } else
    GetConfig().paramMask &= ~Config::USE_PREVIOUS_PARAMETERS;

  std::vector<SegPairs> segPairs;
  if (!(GetConfig().paramMask & Config::CALCULATE_REACTION_RATE)) {
    if (!readSegmentFile(GetConfig(), segPairs)) return;
  } else {
    GetConfig().rateParams.entrancePair = runTab->rateEntranceKey->key();
    GetConfig().rateParams.exitPair = runTab->rateExitKey->key();
    if (GetConfig().rateParams.entrancePair == GetConfig().rateParams.exitPair) {
      QMessageBox::information(this, tr("No Scattering Rates"),
                               tr("Reaction rates cannot be calculated for elastic scattering."));
      return;
    }
    if (!runTab->fileTempButton->isChecked()) {
      GetConfig().rateParams.useFile = false;
      GetConfig().rateParams.minTemp = runTab->minTempText->text().toDouble();
      GetConfig().rateParams.maxTemp = runTab->maxTempText->text().toDouble();
      GetConfig().rateParams.tempStep = runTab->tempStepText->text().toDouble();
    } else {
      GetConfig().rateParams.useFile = true;
      GetConfig().rateParams.temperatureFile = runTab->fileTempText->text().toStdString();
    }
    SegPairs tempPair = {runTab->rateEntranceKey->key(),
                         runTab->rateExitKey->key()};
    segPairs.push_back(tempPair);
  }
  if (segPairs.size() == 0) {
    QMessageBox::information(this, tr("Empty Segments"), tr("No active segments have been found."));
    return;
  }
  int maxPairs = pairsTab->getPairsModel()->getPairs().size();
  for (std::vector<SegPairs>::const_iterator it = segPairs.begin(); it < segPairs.end(); it++) {
    if (it->secondPair == -1) {
      QList<PairsData> pairsList = pairsTab->getPairsModel()->getPairs();
      int i;
      for (i = 0; i < pairsList.size(); i++)
        if (pairsList[i].pairType == 10) break;
      if (i == pairsList.size()) {
        QMessageBox::information(this, tr("No Capture Pairs"),
                                 tr("Total capture is specified, but no capture pairs exist."));
        return;
      }
    } else if (it->firstPair > maxPairs || it->secondPair > maxPairs || it->firstPair < 1 || it->secondPair < 1) {
      QMessageBox::information(this, tr("Undefined Key"), tr("An undefined pair key is specified."));
      return;
    }
  }

  GetConfig().paramMask &= ~Config::USE_EXTERNAL_CAPTURE;
  if (!checkExternalCapture(GetConfig(), segPairs)) return;
  if (GetConfig().paramMask & Config::USE_EXTERNAL_CAPTURE) {
    if (runTab->oldIntegralsFileButton->isChecked() &&
        !(GetConfig().paramMask & Config::CALCULATE_REACTION_RATE)) {
      GetConfig().paramMask |= Config::USE_PREVIOUS_INTEGRALS;
      GetConfig().integralsfile = runTab->integralsFileText->text().toStdString();
    } else
      GetConfig().paramMask &= ~Config::USE_PREVIOUS_INTEGRALS;
  }

  if (!QDir(QString::fromStdString(GetConfig().outputdir)).exists()) {
    QMessageBox::information(this, tr("Directory Doesn't Exist"),
                             tr("The specified output directory doesn't exist."));
    return;
  }
  if (!QDir(QString::fromStdString(GetConfig().checkdir)).exists()) {
    QMessageBox::information(this, tr("Directory Doesn't Exist"),
                             tr("The specified checks directory doesn't exist."));
    return;
  }
  if ((GetConfig().paramMask & Config::USE_PREVIOUS_PARAMETERS) &&
      !QFile(QString::fromStdString(GetConfig().paramfile)).exists()) {
    QMessageBox::information(this, tr("File Doesn't Exist"),
                             tr("The specified parameter file doesn't exist."));
    return;
  }
  if (((GetConfig().paramMask & Config::USE_PREVIOUS_INTEGRALS) &&
       (GetConfig().paramMask & Config::USE_EXTERNAL_CAPTURE)) &&
      !QFile(QString::fromStdString(GetConfig().integralsfile)).exists()) {
    QMessageBox::information(this, tr("File Doesn't Exist"),
                             tr("The specified integrals file doesn't exist."));
    return;
  }
  if ((GetConfig().paramMask & Config::CALCULATE_REACTION_RATE &&
       GetConfig().rateParams.useFile) &&
      !QFile(QString::fromStdString(GetConfig().rateParams.temperatureFile)).exists()) {
    QMessageBox::information(this, tr("File Doesn't Exist"),
                             tr("The specified rate temperature file doesn't exist."));
    return;
  }

  if (!(GetConfig().paramMask & Config::CALCULATE_WITH_DATA) &&
      !(GetConfig().paramMask & Config::CALCULATE_REACTION_RATE)) {
    QList<TargetIntData> targetIntData = targetIntTab->getTargetIntModel()->getLines();
    QList<SegmentsTestData> segmentsTestData = segmentsTab->getSegmentsTestModel()->getLines();
    for (unsigned int i = 0; i < targetIntData.size(); i++) {
      if (targetIntData.at(i).isActive == 1 &&
          (targetIntData.at(i).isTargetIntegration || targetIntData.at(i).isConvolution)) {
        unsigned int j = 0;
        unsigned int lastSegNum = 0;
        bool inclusive = false;
        QList<unsigned int> tempList;
        QString segmentsList = targetIntData.at(i).segmentsList;
        while (j < segmentsList.length()) {
          if (segmentsList[j] >= '0' && segmentsList[j] <= '9') {
            QString tempString;
            while (segmentsList[j] != ',' && segmentsList[j] != '-' &&
                   j < segmentsList.length()) {
              tempString += segmentsList[j];
              j++;
            }
            QTextStream stm(&tempString);
            unsigned int tempSegNum;
            stm >> tempSegNum;
            if (inclusive == true)
              for (int k = lastSegNum + 1; k <= tempSegNum; k++)
                tempList.push_back(k);
            else
              tempList.push_back(tempSegNum);
            lastSegNum = tempSegNum;
          }
          if (segmentsList[j] == '-')
            inclusive = true;
          else
            inclusive = false;
          j++;
        }
        bool isAngularDistribution = false;
        for (j = 0; j < tempList.size(); j++) {
          if (tempList.at(j) <= segmentsTestData.size()) {
            for (int k = 0; k < segmentsTestData.size(); k++) {
              if (segmentsTestData.at(k).isActive == 1 &&
                  tempList.at(j) - 1 == k &&
                  segmentsTestData.at(k).dataType == 3) {
                isAngularDistribution = true;
                break;
              }
            }
          }
          if (isAngularDistribution) break;
        }
        if (isAngularDistribution) {
          QMessageBox::information(this, tr("Incompatable Options"),
                                   tr("Angular distribution coefficients cannot be used with convolution or target integration."));
          return;
        }
      }
    }
  }

  GetConfig().includeAllJPi = levelsTab->includeAllJPi();
  GetConfig().hardSphereMaxL = levelsTab->maxL();
  azureMain = new AZUREMainThread(runTab, GetConfig());
  connect(azureMain, SIGNAL(finished()), this, SLOT(DeleteThread()));
  setWindowTitle(QString("AZURE2 -- %1 -- Running").arg(QString::fromStdString(GetConfig().configfile)));
  runTab->calcButton->setEnabled(false);
  runTab->stopAZUREButton->setEnabled(true);
  runTab->runtimeText->SetMouseFiltered(true);
  startMessage(azureMain->configure());
  azureMain->start();
}

void AZURESetup::DeleteThread() {
  exitMessage(azureMain->configure());
  QScrollBar *sb = runTab->runtimeText->verticalScrollBar();
  sb->setValue(sb->maximum());

  setWindowTitle(QString("AZURE2 -- %1").arg(QString::fromStdString(GetConfig().configfile)));
  runTab->calcButton->setEnabled(true);
  runTab->stopAZUREButton->setEnabled(false);
  runTab->runtimeText->SetMouseFiltered(false);
  delete azureMain;
}

#ifdef USE_MCMC
void AZURESetup::SaveAndRunMCMC() {
  save();
  if (GetConfig().configfile.empty()) return;
  mcmcTab->logTextEdit->clear();
  QFile file(QString::fromStdString(GetConfig().configfile));
  QFileInfo info(file);
  QString directory = info.absolutePath();
  if (GetConfig().outputdir.empty()) GetConfig().outputdir = QDir::fromNativeSeparators(directory).toStdString() + '/';
  if (GetConfig().checkdir.empty()) GetConfig().checkdir = QDir::fromNativeSeparators(directory).toStdString() + '/';

  // Set MCMC-specific config flags (MCMC always fits with data)
  GetConfig().paramMask |= Config::PERFORM_FIT;
  GetConfig().paramMask |= Config::CALCULATE_WITH_DATA;
  GetConfig().paramMask &= ~Config::PERFORM_ERROR_ANALYSIS;
  GetConfig().paramMask &= ~Config::CALCULATE_REACTION_RATE;

  // Read segments same way as SaveAndRun
  std::vector<SegPairs> segPairs;
  if (!readSegmentFile(GetConfig(), segPairs)) return;

  if (segPairs.size() == 0) {
    QMessageBox::information(this, tr("Empty Segments"), tr("No active segments have been found."));
    return;
  }

  // Check segments validity (same as SaveAndRun)
  int maxPairs = pairsTab->getPairsModel()->getPairs().size();
  for (std::vector<SegPairs>::const_iterator it = segPairs.begin(); it < segPairs.end(); it++) {
    if (it->secondPair == -1) {
      QList<PairsData> pairsList = pairsTab->getPairsModel()->getPairs();
      int i;
      for (i = 0; i < pairsList.size(); i++)
        if (pairsList[i].pairType == 10) break;
      if (i == pairsList.size()) {
        QMessageBox::information(this, tr("No Capture Pairs"),
                                 tr("Total capture is specified, but no capture pairs exist."));
        return;
      }
    } else if (it->firstPair > maxPairs || it->secondPair > maxPairs || it->firstPair < 1 || it->secondPair < 1) {
      QMessageBox::information(this, tr("Undefined Key"), tr("An undefined pair key is specified."));
      return;
    }
  }

  // Check external capture (same as SaveAndRun)
  GetConfig().paramMask &= ~Config::USE_EXTERNAL_CAPTURE;
  if (!checkExternalCapture(GetConfig(), segPairs)) return;

  // Check directories exist
  if (!QDir(QString::fromStdString(GetConfig().outputdir)).exists()) {
    QMessageBox::information(this, tr("Directory Doesn't Exist"),
                             tr("The specified output directory doesn't exist."));
    return;
  }
  if (!QDir(QString::fromStdString(GetConfig().checkdir)).exists()) {
    QMessageBox::information(this, tr("Directory Doesn't Exist"),
                             tr("The specified checks directory doesn't exist."));
    return;
  }

  // Create dedicated MCMC thread (simplified like AZUREMainThread)
  azureMCMC = new AZUREMCMCThread(mcmcTab, GetConfig(), this);
  connect(azureMCMC, &AZUREMCMCThread::finished, this, &AZURESetup::DeleteMCMCThread);
  // Note: logMessage, samplingError, and samplingComplete are handled internally by the thread

  setWindowTitle(QString("AZURE2 -- %1 -- Running MCMC").arg(QString::fromStdString(GetConfig().configfile)));
  mcmcTab->runButton->setEnabled(false);
  mcmcTab->stopButton->setEnabled(true);

  // Connect stop button to thread stop method
  connect(mcmcTab->stopButton, &QPushButton::clicked, this, [this]() {
    if (azureMCMC) {
      azureMCMC->stop();
      mcmcTab->logTextEdit->append("Stopping MCMC sampling...");
    }
  });

  startMessage(azureMCMC->configure());
  azureMCMC->start();
}

void AZURESetup::DeleteMCMCThread() {
  exitMessage(azureMCMC->configure());
  QScrollBar *sb = mcmcTab->logTextEdit->verticalScrollBar();
  sb->setValue(sb->maximum());

  setWindowTitle(QString("AZURE2 -- %1").arg(QString::fromStdString(GetConfig().configfile)));
  mcmcTab->runButton->setEnabled(true);
  mcmcTab->stopButton->setEnabled(false);

  // Call mcmcFinished to reset MCMCTab state
  mcmcTab->mcmcFinished();

  delete azureMCMC;
  azureMCMC = nullptr;
}
#endif

void AZURESetup::showAbout() {
  AboutAZURE2Dialog aboutDialog;
  aboutDialog.exec();
}

void AZURESetup::reset() {
  GetConfig().Reset();
  aMatrixAction->activate(QAction::Trigger);
  levelsTab->reset();
  pairsTab->reset();
  segmentsTab->reset();
  targetIntTab->reset();
  fittingTab->reset();
  runTab->reset();
#ifdef USE_QWT
  plotTab->reset();
#endif
  hasThmBlock = false;
  thmBlockLines.clear();
  setThmEnabled(false);
  setWindowTitle(tr("AZURE2 -- untitled"));
  GetConfig().configfile = "";
}

void AZURESetup::showTabInfo() {
  QString tabTitle = tabWidget->tabText(tabWidget->currentIndex()).remove(QChar('&'));
  if (tabWidget->currentIndex() == 0) pairsTab->showInfo(0, tabTitle);
  if (tabWidget->currentIndex() == 1) levelsTab->showInfo(0, tabTitle);
  if (tabWidget->currentIndex() == 2) segmentsTab->showInfo(0, tabTitle);
  if (tabWidget->currentIndex() == 3) targetIntTab->showInfo(0, tabTitle);
  if (tabWidget->currentIndex() == 4) runTab->showInfo(0, tabTitle);
#ifdef USE_QWT
  if (tabWidget->currentIndex() == 5) plotTab->showInfo(0, tabTitle);
#endif
}

void AZURESetup::openWebsite() {
  if (!QDesktopServices::openUrl(QUrl("https://azure.nd.edu")))
    QMessageBox::information(this,
                             tr("Can't Open Browser"),
                             tr("AZURE2 could not access your web browser.  "
                                "Please navitgate to https://azure.nd.edu/ "
                                "to visit the website."));
}

/*!
 * A Config whose configfile is a snapshot of the *current* GUI state.
 *
 * CNuc::Fill reads the .azr from disk, so any conversion built on it sees the
 * last saved file rather than what the user is looking at.  That is wrong for
 * anything the conversion depends on -- channel radii, spins -- and badly wrong
 * for the per-channel "Reduced Width Amplitude" flag, since CNuc::GetTransformParams
 * uses it to decide whether to convert an amplitude at all.  Toggling that
 * button and loading a param.sav without saving first would otherwise convert a
 * value that should have passed through, or pass through one that should have
 * been converted, changing it by orders of magnitude.
 *
 * The snapshot is written to a temporary file that is removed when the returned
 * guard goes out of scope; the open project (its file name, the window title,
 * the working directory) is left as it was.  If writing fails the original Config is returned
 * unchanged, so the conversion degrades to the old behaviour rather than
 * failing outright.
 */
AZURESetup::GuiStateSnapshot::GuiStateSnapshot(AZURESetup *setup) :
  config(setup->GetConfig()),
  ok(false) {
  file.setFileTemplate(QDir::tempPath() + "/azure2_gui_state_XXXXXX.azr");
  if (!file.open()) return;
  // writeProject, not writeFile: writeFile makes the file the open project
  // (configfile, window title, working directory), and this one is deleted
  // with the guard.  Relative data paths stay relative to the project.
  QTextStream out(&file);
  const bool written = setup->writeProject(out, setup->projectDirectory());
  out.flush();
  file.close();
  if (!written) return;
  config.configfile = file.fileName().toStdString();
  ok = true;
}

double AZURESetup::ConvertRWAToPhysical(const QString &paramName, double rwaValue) {
  // Convert RWA parameter to physical value using proper R-Matrix transformation
  // This is the inverse of ParameterLimitsManager::ConvertPhysicalLimitToReduced

  try {
    // Convert against a snapshot of the current GUI state, not the last saved
    // file -- see GuiStateSnapshot.
    GuiStateSnapshot snap(this);
    const Config &cfg = snap.ok ? snap.config : config;

    // Create compound nucleus and data objects
    CNuc *compound = new CNuc();
    EData *data = new EData();

    if (compound->Fill(cfg) == -1) {
      delete compound;
      delete data;
      return rwaValue;  // Return original value if compound creation fails
    }

    // Fill data object if needed for parameter context
    if (config.paramMask & Config::CALCULATE_WITH_DATA) {
      if (data->Fill(cfg, compound) == -1) {
        delete compound;
        delete data;
        return rwaValue;
      }
    } else {
      if (data->MakePoints(cfg, compound) == -1) {
        delete compound;
        delete data;
        return rwaValue;
      }
    }

    // Initialize compound nucleus
    compound->Initialize(cfg);

    // Create parameter objects
    AZUREParams params;
    compound->FillMnParams(params.GetMinuitParams(), &cfg);
    data->FillMnParams(params.GetMinuitParams());

    // Find the parameter index
    int paramIndex = -1;
    std::string paramNameStd = paramName.toStdString();
    for (int i = 0; i < params.GetMinuitParams().Params().size(); ++i) {
      if (params.GetMinuitParams().Parameter(i).GetName() == paramNameStd) {
        paramIndex = i;
        break;
      }
    }

    if (paramIndex == -1) {
      delete compound;
      delete data;
      return rwaValue;  // Parameter not found
    }

    // Get current RWA parameters and set our specific value
    vector_r rwaParams = params.GetMinuitParams().Params();
    rwaParams[paramIndex] = rwaValue;

    // Fill compound with RWA parameters and transform to physical
    compound->FillCompoundFromParams(rwaParams);
    compound->CalcShiftFunctions(cfg);
    compound->TransformOut(cfg);

    // Get the transformed (physical) parameters
    vector_r physicalParams = compound->GetTransformParams(cfg);

    double result = rwaValue;  // Default fallback
    if (paramIndex < physicalParams.size()) {
      result = physicalParams[paramIndex];
    }

    delete compound;
    delete data;

    return result;

  } catch (...) {
    return rwaValue;  // Return original value on any error
  }
}

std::vector<double> AZURESetup::BatchConvertRWAToPhysical(const QStringList &paramNames, const std::vector<double> &rwaValues) {
  std::vector<double> results;

  if (paramNames.size() != rwaValues.size()) {
    // Return original values if sizes don't match
    for (double val : rwaValues) {
      results.push_back(val);
    }
    return results;
  }

  try {
    // Convert against a snapshot of the current GUI state -- see GuiStateSnapshot.
    GuiStateSnapshot snap(this);
    const Config &cfg = snap.ok ? snap.config : config;

    // Create compound nucleus and data objects ONCE
    CNuc *compound = new CNuc();
    EData *data = new EData();

    if (compound->Fill(cfg) == -1) {
      delete compound;
      delete data;
      // Return original values if compound creation fails
      for (double val : rwaValues) {
        results.push_back(val);
      }
      return results;
    }

    // Fill data object if needed for parameter context
    if (config.paramMask & Config::CALCULATE_WITH_DATA) {
      if (data->Fill(cfg, compound) == -1) {
        delete compound;
        delete data;
        for (double val : rwaValues) {
          results.push_back(val);
        }
        return results;
      }
    } else {
      if (data->MakePoints(cfg, compound) == -1) {
        delete compound;
        delete data;
        for (double val : rwaValues) {
          results.push_back(val);
        }
        return results;
      }
    }

    // Initialize compound nucleus
    compound->Initialize(cfg);

    // Create parameter objects
    AZUREParams params;
    compound->FillMnParams(params.GetMinuitParams(), &cfg);
    data->FillMnParams(params.GetMinuitParams());

    // Get current RWA parameters as base
    vector_r rwaParams = params.GetMinuitParams().Params();

    // Find parameter indices and set values
    std::vector<int> paramIndices;
    for (int i = 0; i < paramNames.size(); i++) {
      const QString &paramName = paramNames[i];
      std::string paramNameStd = paramName.toStdString();

      int paramIndex = -1;
      for (int j = 0; j < params.GetMinuitParams().Params().size(); ++j) {
        if (params.GetMinuitParams().Parameter(j).GetName() == paramNameStd) {
          paramIndex = j;
          break;
        }
      }

      paramIndices.push_back(paramIndex);
      if (paramIndex >= 0 && paramIndex < rwaParams.size()) {
        rwaParams[paramIndex] = rwaValues[i];
      }
    }

    // Fill compound with ALL RWA parameters and transform to physical ONCE
    compound->FillCompoundFromParams(rwaParams);
    compound->CalcShiftFunctions(cfg);
    compound->TransformOut(cfg);

    // Get the transformed (physical) parameters
    vector_r physicalParams = compound->GetTransformParams(cfg);

    // Extract results for each requested parameter
    for (int i = 0; i < paramIndices.size(); i++) {
      int paramIndex = paramIndices[i];
      if (paramIndex >= 0 && paramIndex < physicalParams.size()) {
        results.push_back(physicalParams[paramIndex]);
      } else {
        results.push_back(rwaValues[i]);  // Fallback to original value
      }
    }

    delete compound;
    delete data;

    return results;

  } catch (...) {
    // Return original values on any error
    for (double val : rwaValues) {
      results.push_back(val);
    }
    return results;
  }
}

/*!
 * Batch convert RWA parameters to physical using the OLD compound structure
 * This is crucial when loading old parameter files where the level structure may have changed
 *
 * @param paramNames List of RWA parameter names (e.g., "width_1_1", "width_1_2", "width_2_1")
 * @param rwaValues The RWA parameter values
 * @param oldLevelChannelCounts Map of energyIndex -> number of channels in that level (from the .sav file)
 * @return Vector of physical parameters in the same order as paramNames
 */
QMap<QPair<int, int>, QString> AZURESetup::EngineParameterNamesOnDisk() {
  QMap<QPair<int, int>, QString> names;
  QList<LevelsData> levels = levelsTab->getLevelsModel()->getLevels();
  QList<ChannelsData> channels = levelsTab->getChannelsModel()->getChannels();
  try {
    CNuc compound;
    if (compound.Fill(config) == -1) return names;
    // CNuc::FillMnParams numbering: energy_N over the levels group by group, width_N_M over the group's channels.
    int energyIndex = 1;
    for (int j = 1; j <= compound.NumJGroups(); j++) {
      JGroup *group = compound.GetJGroup(j);
      for (int la = 1; la <= group->NumLevels(); la++) {
        ALevel *level = group->GetLevel(la);
        int row = -1;
        for (int r = 0; r < levels.size() && row < 0; r++) {
          const LevelsData &l = levels.at(r);
          if (l.isActive == 1 && l.jValue == group->GetJ() && l.piValue == group->GetPi() &&
              std::fabs(l.energy - level->GetE()) <= 1.0e-9 * std::max(1.0, std::fabs(l.energy)))
            row = r;
        }
        if (row >= 0) {
          names.insert(qMakePair(row, -1), QString("energy_%1").arg(energyIndex));
          for (int ch = 1; ch <= group->NumChannels(); ch++) {
            AChannel *channel = group->GetChannel(ch);
            const int pairKey = compound.GetPair(channel->GetPairNum())->GetPairKey();
            for (int c = 0; c < channels.size(); c++) {
              const ChannelsData &d = channels.at(c);
              if (d.levelIndex == row && d.pairIndex + 1 == pairKey && d.lValue == channel->GetL() &&
                  std::fabs(d.sValue - channel->GetS()) < 1.0e-6) {
                names.insert(qMakePair(row, c), QString("width_%1_%2").arg(energyIndex).arg(ch));
                break;
              }
            }
          }
        }
        energyIndex++;
      }
    }
  } catch (...) {
    names.clear();
  }
  return names;
}

std::vector<double> AZURESetup::BatchConvertRWAToPhysicalWithOldStructure(
    const QStringList &paramNames, const std::vector<double> &rwaValues, const QMap<int, int> &oldLevelChannelCounts) {
  std::vector<double> results;

  if (paramNames.size() != rwaValues.size()) {
    // Return original values if sizes don't match
    for (double val : rwaValues) {
      results.push_back(val);
    }
    return results;
  }

  try {
    // Convert against a snapshot of the current GUI state -- see GuiStateSnapshot.
    GuiStateSnapshot snap(this);
    const Config &cfg = snap.ok ? snap.config : config;

    // Create compound nucleus and data objects ONCE
    CNuc *compound = new CNuc();
    EData *data = new EData();

    if (compound->Fill(cfg) == -1) {
      delete compound;
      delete data;
      // Return original values if compound creation fails
      for (double val : rwaValues) {
        results.push_back(val);
      }
      return results;
    }

    // Initialize compound nucleus
    compound->Initialize(cfg);

    // Create parameter objects using the CURRENT structure
    AZUREParams params;
    compound->FillMnParams(params.GetMinuitParams(), &cfg);
    data->FillMnParams(params.GetMinuitParams());

    // Get current RWA parameters as base
    vector_r rwaParams = params.GetMinuitParams().Params();

    // CRITICAL: Map old RWA parameter names to new parameter indices
    // We need to handle the case where level ordering/structure has changed
    std::vector<int> paramIndices;

    for (int i = 0; i < paramNames.size(); i++) {
      const QString &paramName = paramNames[i];
      std::string paramNameStd = paramName.toStdString();

      // Extract energyIndex and widthIndex from old parameter name (e.g., "width_1_2" -> energyIndex=1, widthIndex=2)
      int oldEnergyIndex = -1;
      int oldWidthIndex = -1;

      QRegExp widthRegex("^width_(\\d+)_(\\d+)$");
      if (widthRegex.indexIn(paramName) != -1) {
        oldEnergyIndex = widthRegex.cap(1).toInt();
        oldWidthIndex = widthRegex.cap(2).toInt();
      } else if (paramName.contains("energy")) {
        QRegExp energyRegex("^energy_(\\d+)$");
        if (energyRegex.indexIn(paramName) != -1) {
          oldEnergyIndex = energyRegex.cap(1).toInt();
          oldWidthIndex = -1;  // Energy parameter, not width
        }
      }

      // Find the corresponding parameter in the NEW structure
      // For now, try direct name matching first (works if structure didn't change much)
      int paramIndex = -1;
      for (int j = 0; j < params.GetMinuitParams().Params().size(); ++j) {
        if (params.GetMinuitParams().Parameter(j).GetName() == paramNameStd) {
          paramIndex = j;
          break;
        }
      }

      paramIndices.push_back(paramIndex);
      if (paramIndex >= 0 && paramIndex < rwaParams.size()) {
        rwaParams[paramIndex] = rwaValues[i];
      }
    }

    // Fill compound with ALL RWA parameters and transform to physical ONCE
    compound->FillCompoundFromParams(rwaParams);
    compound->CalcShiftFunctions(cfg);
    compound->TransformOut(cfg);

    // Get the transformed (physical) parameters
    vector_r physicalParams = compound->GetTransformParams(cfg);

    // Extract results for each requested parameter
    for (int i = 0; i < paramIndices.size(); i++) {
      int paramIndex = paramIndices[i];
      if (paramIndex >= 0 && paramIndex < physicalParams.size()) {
        results.push_back(physicalParams[paramIndex]);
      } else {
        results.push_back(rwaValues[i]);  // Fallback to original value
      }
    }

    delete compound;
    delete data;

    return results;

  } catch (...) {
    // Return original values on any error
    for (double val : rwaValues) {
      results.push_back(val);
    }
    return results;
  }
}

void AZURESetup::saveTabScreenshots(const QString &dir) {
  QDir().mkpath(dir);
  int current = tabWidget->currentIndex();
  for (int i = 0; i < tabWidget->count(); i++) {
    tabWidget->setCurrentIndex(i);
    QApplication::processEvents();
    QString name = tabWidget->tabText(i).remove('&').toLower().replace(' ', '_');
    this->grab().save(QString("%1/tab%2_%3.png").arg(dir).arg(i).arg(name));
  }
  tabWidget->setCurrentIndex(current);
}
