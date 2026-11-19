// The main window: tabs, menus, the run's worker thread (docs/upgrade/qt_gui.md).  Files: AZURESetupFile.cpp.
#include "AZURESetup.h"

#include <QAction>
#include <QActionGroup>
#include <QApplication>
#include <QCloseEvent>
#include <QDesktopServices>
#include <QDir>
#include <QDockWidget>
#include <QFileInfo>
#include <QLineEdit>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QRegularExpression>
#include <QSettings>
#include <QTabWidget>
#include <QToolButton>
#include <QUrl>

#include "AboutAZURE2Dialog.h"
#include "ConfigDialogs.h"
#include "DiagnosticsPanel.h"
#include "FittingTab.h"
#include "LevelsTab.h"
#include "PairsTab.h"
#include "RecordTabs.h"
#include "RecordsPanel.h"
#include "MCMCTab.h"
#include "PlotTab.h"
#include "PythonTools.h"
#include "RunTab.h"
#include "ThmTab.h"

AZURESetup::AZURESetup() {
  qRegisterMetaType<RunOutcome>("RunOutcome");
  holder_ = new ProjectHolder(this);
  setMinimumSize(1000, 640);

  tabWidget = new QTabWidget();
  // The logo at the right end of the tab bar; clicking it opens the About box.
  logoButton_ = new QToolButton;
  logoButton_->setObjectName("logoButton");
  logoButton_->setIcon(QIcon(":/azure-icon.png"));
  logoButton_->setIconSize(QSize(26, 26));
  logoButton_->setText("AZURE2");
  logoButton_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
  logoButton_->setAutoRaise(true);
  logoButton_->setCursor(Qt::PointingHandCursor);
  logoButton_->setToolTip(tr("AZURE2 -- about this program"));
  logoButton_->setStyleSheet("QToolButton { border: none; padding: 0 8px; color: #1f4f8f; font-weight: bold; }"
                             "QToolButton:hover { background: #e3edf8; border-radius: 4px; }");
  connect(logoButton_, &QToolButton::clicked, this, [this]() {
    AboutAZURE2Dialog about(this, polarizedLogo_);
    about.exec();
  });
  tabWidget->setCornerWidget(logoButton_, Qt::TopRightCorner);

  pairsTab = new PairsTab;
  levelsTab = new LevelsTab;
  levelsTab->setPairsModel(pairsTab->getPairsModel());
  connect(pairsTab, SIGNAL(pairAdded(int)), levelsTab, SLOT(updateChannelsPairAddedEdited()));
  connect(pairsTab, SIGNAL(pairEdited(int)), levelsTab, SLOT(updateChannelsPairAddedEdited()));
  connect(pairsTab, SIGNAL(pairRemoved(int)), levelsTab, SLOT(updateChannelsPairRemoved(int)));
  connect(levelsTab, SIGNAL(readNewPair(PairsData, int, bool)), pairsTab, SLOT(addPair(PairsData, int, bool)));
  connect(levelsTab, SIGNAL(readExistingPair(PairsData, int, bool)), pairsTab, SLOT(editPair(PairsData, int, bool)));
  // the rows of <levels> are written back only when the pairs, levels or channels changed (levelsDirty_)
  for (QAbstractItemModel *m : {static_cast<QAbstractItemModel *>(pairsTab->getPairsModel()),
                                static_cast<QAbstractItemModel *>(levelsTab->getLevelsModel()),
                                static_cast<QAbstractItemModel *>(levelsTab->getChannelsModel())}) {
    connect(m, &QAbstractItemModel::dataChanged, this, &AZURESetup::levelsEdited);
    connect(m, &QAbstractItemModel::rowsInserted, this, &AZURESetup::levelsEdited);
    connect(m, &QAbstractItemModel::rowsRemoved, this, &AZURESetup::levelsEdited);
  }
  connect(levelsTab->getChannelsModel(), &ChannelsModel::widthEdited, this, &AZURESetup::levelsEdited);
  pairsTab->addExtraPanel(
      tr("Per-channel radii"),
      new RecordsPanel(holder_, "channelRadii",
                       tr("A channel's own radius (<channelRadii>): a record per pair and, optionally, l, s and J; every "
                          "channel it does not name keeps its pair's radius (the column a of the pairs).")));
  levelsTab->addExtraPanel(
      tr("Gamma Decays"),
      new RecordsPanel(holder_, "gammaDecays",
                       tr("A level's gamma decay written as its total radiative width and branchings (<gammaDecays>) "
                          "in place of its photon amplitudes; fitted as such.")));
  levelsTab->addExtraPanel(
      tr("Isospin"),
      new RecordsPanel(holder_, "isospin",
                       tr("Isospin-labelled levels (<isospin>): a level of one ID over the nuclei of an isospin "
                          "multiplet, its isospin, the Coulomb shifts and mixing, the reference and factors.")));

  segmentsTab = new SectionsTab(holder_, SegmentsParts());
  effectsTab = new SectionsTab(holder_, EffectsParts());
  potentialTab = new SectionsTab(holder_, PotentialParts());
  thmTab = new ThmTab(holder_);
  fittingTab = new FittingTab(holder_);
  runTab = new RunTab;
  plotTab = new PlotTab(holder_);
  connect(runTab, &RunTab::runRequested, this, &AZURESetup::SaveAndRun);
  connect(fittingTab, &FittingTab::startFileChosen, this,
          [this](const QString &file) { runTab->presetOption("--parameter-file=", file); });
  connect(fittingTab, &FittingTab::parametersRequested, this, [this]() {
    QString why;
    if (running_) why = tr("a run is in progress");
    else if (writeProjectFile(holder_->path(), why)) fittingTab->readParameters(holder_->path(), why);
    if (!why.isEmpty()) QMessageBox::warning(this, tr("Parameters"), why);
  });
  connect(runTab, &RunTab::stopRequested, this, &AZURESetup::stopRun);
  // the rate's pair drop-downs follow the Pairs tab (key = place in the list, as <levels> numbers the pairs)
  auto pairsToRun = [this]() {
    QList<QPair<int, QString>> pairs;
    const QList<PairsData> list = pairsTab->getPairsModel()->getPairs();
    for (int i = 0; i < list.size(); i++)
      pairs << qMakePair(i + 1, QString(pairsTab->getPairsModel()->getPairLabel(list[i])).remove(QRegularExpression("<[^>]*>")));
    runTab->setPairs(pairs);
  };
  for (auto signal : {&QAbstractItemModel::rowsInserted, &QAbstractItemModel::rowsRemoved})
    connect(pairsTab->getPairsModel(), signal, this, pairsToRun);
  connect(pairsTab->getPairsModel(), &QAbstractItemModel::dataChanged, this, pairsToRun);
  connect(pairsTab->getPairsModel(), &QAbstractItemModel::modelReset, this, pairsToRun);

  tabWidget->addTab(pairsTab, tr("&Particle Pairs"));
  tabWidget->addTab(levelsTab, tr("&Levels and Channels"));
  tabWidget->addTab(segmentsTab, tr("&Segments"));
  tabWidget->addTab(effectsTab, tr("&Experimental Effects"));
  tabWidget->addTab(potentialTab, tr("&Nuclear Potential"));
  tabWidget->addTab(thmTab, tr("&THM"));
  tabWidget->addTab(fittingTab, tr("&Fitting Settings"));
  tabWidget->addTab(runTab, tr("&Calculate"));
  // The MCMC tab only where the engine has the command (front::HasCommand): a sampler removed from the engine takes
  // the tab with it, without a change here.
  if (azure2::front::HasCommand("mcmc")) {
    mcmcTab = new MCMCTab(holder_);
    tabWidget->addTab(mcmcTab, tr("&MCMC"));
    connect(mcmcTab, &MCMCTab::parametersRequested, this, [this]() {
      QString why;
      if (running_) why = tr("a run is in progress");
      else if (writeProjectFile(holder_->path(), why)) mcmcTab->readParameters(holder_->path(), why);
      if (!why.isEmpty()) QMessageBox::warning(this, tr("Parameters"), why);
    });
    connect(mcmcTab, &MCMCTab::runRequested, this, [this]() {
      QString why;
      if (!startMCMC(why)) QMessageBox::warning(this, tr("Not run"), why);
    });
    connect(mcmcTab, &MCMCTab::stopRequested, this, &AZURESetup::stopRun);
  }
  tabWidget->addTab(plotTab, tr("Pl&ot"));
  setCentralWidget(tabWidget);

  diagnostics = new DiagnosticsPanel(holder_);
  connect(diagnostics, &DiagnosticsPanel::showRecord, this, &AZURESetup::showRecord);
  diagnosticsDock = new QDockWidget(tr("Diagnostics"), this);
  diagnosticsDock->setObjectName("diagnosticsDock");
  diagnosticsDock->setWidget(diagnostics);
  diagnosticsDock->setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
  addDockWidget(Qt::BottomDockWidgetArea, diagnosticsDock);
  resizeDocks({diagnosticsDock}, {140}, Qt::Vertical);

  // The run's worker in its own thread; queued connections only (RunWorker.h).
  worker_ = new RunWorker;
  worker_->moveToThread(&runThread_);
  connect(&runThread_, &QThread::finished, worker_, &QObject::deleteLater);
  connect(this, &AZURESetup::startWorker, worker_, &RunWorker::run, Qt::QueuedConnection);
  connect(worker_, &RunWorker::message, runTab, &RunTab::appendLine, Qt::QueuedConnection);
  connect(worker_, &RunWorker::progress, runTab, &RunTab::showProgress, Qt::QueuedConnection);
  if (mcmcTab) {
    connect(worker_, &RunWorker::progress, mcmcTab, &MCMCTab::showProgress, Qt::QueuedConnection);
    connect(worker_, &RunWorker::message, mcmcTab, &MCMCTab::appendLine, Qt::QueuedConnection);
  }
  connect(worker_, &RunWorker::finished, this, &AZURESetup::runEnded, Qt::QueuedConnection);
  runThread_.start();

  connect(holder_, &ProjectHolder::opened, this, &AZURESetup::projectOpened);
  connect(holder_, &ProjectHolder::sectionEdited, this, [this]() {
    updateTitle();
    updateLogo();
  });

  tools_ = new PythonTools(this);
  connect(tools_, &PythonTools::line, runTab, &RunTab::appendLine);
  connect(tools_, &PythonTools::webPageReady, this, [](const QString &url) { QDesktopServices::openUrl(QUrl(url)); });
  createActions();
  createMenus();
  projectOpened();
}

AZURESetup::~AZURESetup() {
  if (worker_) worker_->requestStop();
  runThread_.quit();
  runThread_.wait();
}

void AZURESetup::createActions() {
  for (int i = 0; i < numRecent; i++) {
    recentFileActions[i] = new QAction(this);
    recentFileActions[i]->setVisible(false);
    connect(recentFileActions[i], &QAction::triggered, this, &AZURESetup::openRecent);
  }
  clearRecentAction = new QAction(tr("&Clear"), this);
  clearRecentAction->setVisible(false);
  connect(clearRecentAction, &QAction::triggered, this, &AZURESetup::clearRecent);
  matrixActionGroup = new QActionGroup(this);
  aMatrixAction = new QAction(tr("&A-Matrix"), this);
  aMatrixAction->setCheckable(true);
  matrixActionGroup->addAction(aMatrixAction);
  rMatrixAction = new QAction(tr("&R-Matrix"), this);
  rMatrixAction->setCheckable(true);
  matrixActionGroup->addAction(rMatrixAction);
  aMatrixAction->setChecked(true);
  connect(matrixActionGroup, &QActionGroup::triggered, this, &AZURESetup::matrixChanged);
}

void AZURESetup::createMenus() {
  QMenu *fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(tr("&About AZURE2..."), this, &AZURESetup::showAbout);
  fileMenu->addSeparator();
  fileMenu->addAction(tr("&New Project"), QKeySequence::New, this, &AZURESetup::newProject);
  fileMenu->addAction(tr("&Open..."), QKeySequence::Open, this, qOverload<>(&AZURESetup::open));
  recentFileMenu = fileMenu->addMenu(tr("Open &Recent..."));
  for (int i = 0; i < numRecent; i++) recentFileMenu->addAction(recentFileActions[i]);
  recentSeparator = recentFileMenu->addSeparator();
  recentFileMenu->addAction(clearRecentAction);
  updateRecent();
  fileMenu->addAction(tr("&Save"), QKeySequence::Save, this, &AZURESetup::save);
  fileMenu->addAction(tr("Save &As..."), QKeySequence::SaveAs, this, &AZURESetup::saveAs);
  fileMenu->addAction(tr("Con&vert an Old Project..."), this, &AZURESetup::convertProject);
  fileMenu->addSeparator();
  fileMenu->addAction(tr("&Quit"), QKeySequence::Quit, this, &QWidget::close);

  QMenu *configMenu = menuBar()->addMenu(tr("Co&nfigure"));
  QMenu *formalismMenu = configMenu->addMenu(tr("&Formalism"));
  formalismMenu->addAction(aMatrixAction);
  formalismMenu->addAction(rMatrixAction);
  configMenu->addAction(tr("&Checks..."), this, &AZURESetup::editChecks);
  configMenu->addAction(tr("&Directories..."), this, &AZURESetup::editDirs);
  configMenu->addAction(tr("&Runtime Options..."), this, &AZURESetup::editOptions);
  configMenu->addSeparator();
  configMenu->addAction(tr("&Validate Project"), QKeySequence(Qt::Key_F7), this, &AZURESetup::validate);
  configMenu->addAction(diagnosticsDock->toggleViewAction());

  // Tools: pyazr's commands as child processes (PythonTools.h), filled when the menu opens (pyazr importable or not)
  QMenu *toolsMenu = menuBar()->addMenu(tr("&Tools"));
  toolsMenu->setObjectName("toolsMenu");
  connect(toolsMenu, &QMenu::aboutToShow, this, [this, toolsMenu]() { tools_->fillMenu(toolsMenu, holder_->path()); });

  QMenu *helpMenu = menuBar()->addMenu(tr("&Documentation"));
  helpMenu->addAction(tr("Show Documentation For Current Tab"), QKeySequence(Qt::CTRL | Qt::Key_D), this,
                      &AZURESetup::showTabInfo);
  helpMenu->addAction(tr("Open AZURE Website..."), this, &AZURESetup::openWebsite);
}

void AZURESetup::updateRecent() {
  const QStringList files = QSettings().value("recentFileList").toStringList();
  const int numFiles = qMin(files.size(), (int)numRecent);
  for (int i = 0; i < numFiles; i++) {
    recentFileActions[i]->setText(tr("&%1 %2").arg(i + 1).arg(QFileInfo(files[i]).fileName()));
    recentFileActions[i]->setData(files[i]);
    recentFileActions[i]->setVisible(true);
  }
  for (int i = numFiles; i < numRecent; i++) recentFileActions[i]->setVisible(false);
  recentSeparator->setVisible(numFiles > 0);
  clearRecentAction->setVisible(numFiles > 0);
}

void AZURESetup::addRecent(const QString &file) {
  QSettings settings;
  QStringList files = settings.value("recentFileList").toStringList();
  files.removeAll(file);
  files.prepend(file);
  while (files.size() > numRecent) files.removeLast();
  settings.setValue("recentFileList", files);
  updateRecent();
}

void AZURESetup::openRecent() {
  if (QAction *a = qobject_cast<QAction *>(sender())) open(a->data().toString());
}

void AZURESetup::clearRecent() {
  QSettings().setValue("recentFileList", QStringList());
  updateRecent();
}

void AZURESetup::updateTitle() {
  const bool modified = holder_->project().Modified() || levelsDirty_;
  setWindowTitle(tr("AZURE2 -- %1%2").arg(QFileInfo(holder_->path()).fileName(), modified ? " *" : ""));
}

void AZURESetup::updateLogo() {
  bool polarized = false;
  for (const char *section : {"data", "extrapolations"})
    for (const azure2::front::Record &r : holder_->records(section)) {
      const std::string q = r.Value("quantity");
      polarized = polarized || q == "analyzing_power" || q == "ay" || q == "spin_observable";
    }
  if (polarized == polarizedLogo_) return;
  polarizedLogo_ = polarized;
  const QIcon icon(polarized ? ":/azure-icon-polarized.png" : ":/azure-icon.png");
  logoButton_->setIcon(icon);
  setWindowIcon(icon);
  qApp->setWindowIcon(icon);
}

void AZURESetup::matrixChanged(QAction *action) {
  QString why;
  if (!holder_->setConfig("amatrix", action == aMatrixAction ? "true" : "false", why))
    QMessageBox::warning(this, tr("Formalism"), why);
}

void AZURESetup::editChecks() {
  ConfigDialog dialog(holder_, tr("Checks"), {"Checks"}, {}, this);
  dialog.exec();
}

void AZURESetup::editDirs() {
  ConfigDialog dialog(holder_, tr("Directories"), {"Directories"}, {}, this);
  QLineEdit *python = new QLineEdit(PythonTools::interpreter());
  python->setObjectName("pythonInterpreter");
  python->setToolTip(tr("The Python interpreter of the Tools menu (python -m pyazr ...): a setting of this computer, "
                        "not of the project."));
  dialog.addMachineRow(tr("Python interpreter (Tools):"), python);
  if (dialog.exec() == QDialog::Accepted) PythonTools::setInterpreter(python->text());
}

void AZURESetup::editOptions() {
  ConfigDialog dialog(holder_, tr("Runtime Options"), {}, {"Checks", "Directories", "File"}, this);
  dialog.exec();
}

void AZURESetup::validate() {
  QString why;
  commitTabs(why);
  diagnosticsDock->show();
  diagnostics->validate();
}

void AZURESetup::showRecord(const QString &section, const QString &record) {
  if (section == "thm") {
    tabWidget->setCurrentWidget(thmTab);
    return;
  }
  for (SectionsTab *tab : {segmentsTab, effectsTab, potentialTab})
    if (tab->show(section, record)) {
      tabWidget->setCurrentWidget(tab);
      return;
    }
  if (section == "channelRadii") tabWidget->setCurrentWidget(pairsTab);
  if (section == "gammaDecays" || section == "isospin" || section == "levels") tabWidget->setCurrentWidget(levelsTab);
  if (section == "parameterSettings" || section == "mcmc") tabWidget->setCurrentWidget(fittingTab);
}

void AZURESetup::showAbout() {
  AboutAZURE2Dialog about(this, polarizedLogo_);
  about.exec();
}

void AZURESetup::showTabInfo() {
  const QString title = tabWidget->tabText(tabWidget->currentIndex()).remove(QChar('&'));
  if (tabWidget->currentWidget() == pairsTab) pairsTab->showInfo(0, title);
  else if (tabWidget->currentWidget() == levelsTab) levelsTab->showInfo(0, title);
  else
    QMessageBox::information(this, title,
                             tr("Each field of a record says what it is in its tool tip (the engine's table of "
                                "project format 3); the user guide's pages follow the tabs."));
}

void AZURESetup::openWebsite() {
  if (!QDesktopServices::openUrl(QUrl("https://azure.nd.edu")))
    QMessageBox::information(this, tr("Can't Open Browser"),
                             tr("AZURE2 could not access your web browser.  Please navigate to https://azure.nd.edu/ "
                                "to visit the website."));
}

void AZURESetup::levelsEdited() {
  if (loadingLevels_) return;
  levelsDirty_ = true;
  updateTitle();
}

void AZURESetup::saveTabScreenshots(const QString &dir) {
  QDir().mkpath(dir);
  const int current = tabWidget->currentIndex();
  for (int i = 0; i < tabWidget->count(); i++) {
    tabWidget->setCurrentIndex(i);
    QApplication::processEvents();
    const QString name = tabWidget->tabText(i).remove('&').toLower().replace(' ', '_');
    grab().save(QString("%1/tab%2_%3.png").arg(dir).arg(i).arg(name));
  }
  tabWidget->setCurrentIndex(current);
}

void AZURESetup::closeEvent(QCloseEvent *event) {
  if (running_) {
    if (QMessageBox::question(this, tr("A run is in progress"), tr("Stop the run and quit?")) != QMessageBox::Yes) {
      event->ignore();
      return;
    }
    worker_->requestStop();
  }
  if (!maybeSave()) {
    event->ignore();
    return;
  }
  event->accept();
}
