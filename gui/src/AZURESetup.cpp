// The main window: tabs, menus, the run's worker thread (docs/upgrade/qt_gui.md).  Files: AZURESetupFile.cpp.
#include "AZURESetup.h"

#include <QAction>
#include <QActionGroup>
#include <QApplication>
#include <QComboBox>
#include <QLocale>
#include <QDialog>
#include <QFontDatabase>
#include <QPlainTextEdit>
#include <QVBoxLayout>
#include <QCloseEvent>
#include <QDesktopServices>
#include <QDir>
#include <QDockWidget>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QFileInfo>
#include <QLabel>
#include <QMimeData>
#include <QProgressBar>
#include <QStackedWidget>
#include <QStatusBar>
#include <QToolBar>
#include <QLineEdit>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QRegularExpression>
#include <QSettings>
#include <QTabWidget>
#include <QToolButton>
#include <QUrl>
#include <QHBoxLayout>
#include <QStyle>
#include <clocale>
#include <cmath>

#include "AboutAZURE2Dialog.h"
#include "AzureTheme.h"
#include "StartPage.h"
#include "ConfigDialogs.h"
#include "DiagnosticsPanel.h"
#include "FittingTab.h"
#include "LevelsTab.h"
#include "PairsTab.h"
#include "RecordTabs.h"
#include "RecordsPanel.h"
#include "MCMCTab.h"
#include "NotesBar.h"
#include "ExpertMode.h"
#include "OverviewPage.h"
#include "PlotWidget.h"
#include "ProjectWatcher.h"
#include "PlotTab.h"
#include "PhysicsTab.h"
#include "PythonTools.h"
#include "RunTab.h"
#include "ThmTab.h"

// the resources of the library (a static library's resources are linked only when named)
static void InitShellResources() {
  Q_INIT_RESOURCE(AZURESetup);
  Q_INIT_RESOURCE(ShellIcons);
}

AZURESetup::AZURESetup() {
  // Numbers are read and written with a decimal point, whatever the user's language: QApplication takes the C
  // library's locale from the environment, and with a decimal-comma locale (de_DE, it_IT, pl_PL, ...) strtod then
  // stops at the point, so the engine's readers refused every "5.441" of a project (E_NOT_A_NUMBER, E_BAD_NAME,
  // E_REACTION_SYNTAX).  Qt's documentation asks for exactly this reset after the application is made.
  std::setlocale(LC_NUMERIC, "C");
  QLocale::setDefault(QLocale::c());
  InitShellResources();
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
  logoButton_->setStyleSheet("QToolButton { border: none; padding: 0 0.6em; color: palette(link); font-weight: bold; }"
                             "QToolButton:hover { background: palette(midlight); border-radius: 4px; }");
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
                       tr("A channel's own radius: a record per pair and, optionally, l, s and J; every "
                          "channel it does not name keeps its pair's radius (the column a of the pairs).")));
  levelsTab->addExtraPanel(
      tr("Gamma Decays"),
      new RecordsPanel(holder_, "gammaDecays",
                       tr("A level's gamma decay written as its total radiative width and branchings "
                          "in place of its photon amplitudes; fitted as such.")));
  levelsTab->addExtraPanel(
      tr("Isospin"),
      new RecordsPanel(holder_, "isospin",
                       tr("Isospin-labelled levels: a level of one ID over the nuclei of an isospin "
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
  physicsTab_ = new PhysicsTab(holder_);
  tabWidget->addTab(physicsTab_, tr("P&hysics"));
  tabWidget->setDocumentMode(false);
  startPage_ = new StartPage;
  connect(startPage_, &StartPage::newRequested, this, &AZURESetup::newProject);
  connect(startPage_, &StartPage::openRequested, this, qOverload<>(&AZURESetup::open));
  connect(startPage_, &StartPage::openFile, this, [this](const QString &file) { open(file); });
  connect(startPage_, &StartPage::expertToggled, this, &AZURESetup::setExpertMode);
  central_ = new QStackedWidget;
  central_->addWidget(startPage_);
  central_->addWidget(tabWidget);
  setCentralWidget(central_);
  setAcceptDrops(true);

  diagnostics = new DiagnosticsPanel(holder_);
  connect(diagnostics, &DiagnosticsPanel::showRecord, this, &AZURESetup::showRecord);
  diagnosticsDock = new QDockWidget(tr("Diagnostics"), this);
  diagnosticsDock->setObjectName("diagnosticsDock");
  diagnosticsDock->setWidget(diagnostics);
  diagnosticsDock->setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable |
                               QDockWidget::DockWidgetClosable);
  addDockWidget(Qt::BottomDockWidgetArea, diagnosticsDock);
  resizeDocks({diagnosticsDock}, {140}, Qt::Vertical);
  // hidden while there is nothing to show; a validation that finds a problem opens it (diagnosticsDone)
  diagnosticsDock->hide();
  connect(diagnostics, &DiagnosticsPanel::validated, this, &AZURESetup::diagnosticsDone);

  // The run's worker in its own thread; queued connections only (RunWorker.h).
  worker_ = new RunWorker;
  worker_->moveToThread(&runThread_);
  connect(&runThread_, &QThread::finished, worker_, &QObject::deleteLater);
  connect(this, &AZURESetup::startWorker, worker_, &RunWorker::run, Qt::QueuedConnection);
  connect(worker_, &RunWorker::message, runTab, &RunTab::appendLine, Qt::QueuedConnection);
  connect(worker_, &RunWorker::progress, runTab, &RunTab::showProgress, Qt::QueuedConnection);
  connect(worker_, &RunWorker::progress, this, &AZURESetup::showJobProgress, Qt::QueuedConnection);
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
    updateExpertTabs();
  });

  tools_ = new PythonTools(this);
  connect(tools_, &PythonTools::line, runTab, &RunTab::appendLine);
  connect(tools_, &PythonTools::webPageReady, this, [](const QString &url) { QDesktopServices::openUrl(QUrl(url)); });
  createActions();
  createToolBar();
  createMenus();
  createStatusBar();
  connect(runTab, &RunTab::plotRequested, this, [this]() { tabWidget->setCurrentWidget(plotTab); });
  // The Overview page (docs/upgrade/qt_gui.md section 15): the first tab, selected when a project opens; it reads the
  // project and the output directory again after a run, a save, an edit, a validation, and when it is shown.
  {
    OverviewPage *overview = new OverviewPage(holder_);
    tabWidget->insertTab(0, overview, tr("O&verview"));
    connect(holder_, &ProjectHolder::opened, overview, [this, overview]() {
      overview->clearProgress();
      overview->refresh();
      tabWidget->setCurrentWidget(overview);
    });
    connect(holder_, &ProjectHolder::sectionEdited, overview, &OverviewPage::refreshLater);
    connect(diagnostics, &DiagnosticsPanel::validated, overview, &OverviewPage::setValidity);
    connect(this, &AZURESetup::startWorker, overview, [overview]() { overview->setRunning(true); });
    connect(this, &AZURESetup::runFinished, overview, [overview]() { overview->setRunning(false); });
    connect(worker_, &RunWorker::progress, overview,
            [overview](const QString &stage, qlonglong evaluation, qlonglong, double, double best) {
              overview->addProgress(stage, evaluation, best);
            });
    for (const char *name : {"save", "saveAs"}) connect(actions_[name], &QAction::triggered, overview, &OverviewPage::refresh);
    connect(tabWidget, &QTabWidget::currentChanged, overview, [this, overview]() {
      if (tabWidget->currentWidget() != overview || !hasProject_ || running_) return;
      QString why;
      commitTabs(why);  // an edit of the Levels tab counted
      overview->refresh();
    });
    connect(overview, &OverviewPage::commandRequested, this, [this](const QString &command) {
      QString why;
      if (!runCommand(command, why)) QMessageBox::warning(this, tr("Not run"), why);
    });
    connect(overview, &OverviewPage::validateRequested, this, &AZURESetup::validate);
    connect(overview, &OverviewPage::diagnosticsRequested, this, [this]() {
      diagnosticsDock->show();
      diagnosticsDock->raise();
    });
    connect(overview, &OverviewPage::plotRequested, this, [this](const QString &record) {
      if (!record.isEmpty() && plotTab->availableRecords().contains(record)) plotTab->selectRecords({record});
      tabWidget->setCurrentWidget(plotTab);
    });
    connect(overview, &OverviewPage::recordRequested, this, &AZURESetup::showRecord);
    connect(overview, &OverviewPage::tabRequested, this, [this](const QString &name) {
      QWidget *to = name == "fitting"    ? static_cast<QWidget *>(fittingTab)
                    : name == "pairs"    ? static_cast<QWidget *>(pairsTab)
                    : name == "levels"   ? static_cast<QWidget *>(levelsTab)
                    : name == "segments" ? static_cast<QWidget *>(segmentsTab)
                                         : static_cast<QWidget *>(runTab);
      tabWidget->setCurrentWidget(to);
    });
    if (QComboBox *fitter = fittingTab->findChild<QComboBox *>("fitterCombo")) {
      connect(fitter, &QComboBox::currentTextChanged, overview, &OverviewPage::setFitterChoice);
      overview->setFitterChoice(fitter->currentText());
    }

    // Changes on disk (docs/upgrade/qt_gui.md section 16; ProjectWatcher.h): the project file changed outside is
    // re-read when the window holds no unsaved edit, else a bar asks; outputs of an outside run and changed data
    // files refresh the Overview and the Plot tab (its zoom kept).  View > Follow changes on disk.
    ProjectWatcher *watcher = new ProjectWatcher(holder_, this);
    watcher->setObjectName("projectWatcher");
    QToolBar *bar = new QToolBar(tr("Changed on disk"), this);
    bar->setObjectName("diskChangeBar");
    bar->setMovable(false);
    QLabel *barText = new QLabel;
    barText->setObjectName("diskChangeText");
    barText->setContentsMargins(8, 0, 12, 0);
    barText->setStyleSheet("font-weight: bold;");
    bar->addWidget(barText);
    QAction *reload = bar->addAction(tr("Reload (discard mine)"));
    reload->setObjectName("diskReload");
    QAction *keep = bar->addAction(tr("Keep mine"));
    keep->setObjectName("diskKeep");
    QAction *differences = bar->addAction(tr("Show differences"));
    differences->setObjectName("diskDifferences");
    addToolBarBreak();
    addToolBar(bar);
    bar->hide();
    // The notes about the project just opened (qt_gui.md section 19.2): one bar, not modal, shown once per open.
    addToolBarBreak();
    notes_ = new NotesBar(this);
    addToolBar(notes_);
    connect(notes_, &NotesBar::diagnosticsRequested, this, [this]() {
      diagnosticsDock->show();
      diagnosticsDock->raise();
    });
    auto reloadProject = [this, bar]() {
      const int tab = tabWidget->currentIndex();
      QMap<RecordsPanel *, QString> shown;
      for (RecordsPanel *p : findChildren<RecordsPanel *>()) shown[p] = p->shownName();
      QString why;
      if (!readProjectFile(holder_->path(), why)) {
        statusBar()->showMessage(tr("The project changed on disk and could not be read: %1").arg(why), 30000);
        return;
      }
      tabWidget->setCurrentIndex(tab);
      for (auto it = shown.begin(); it != shown.end(); ++it)
        if (!it.value().isEmpty()) it.key()->selectRecord(it.value());
      bar->hide();
      statusBar()->showMessage(tr("Project reloaded: changed on disk."), 15000);
    };
    connect(watcher, &ProjectWatcher::projectChanged, this, [this, bar, barText, reloadProject]() {
      if (!holder_->project().Modified() && !levelsDirty_) {
        reloadProject();
        return;
      }
      barText->setText(tr("The project changed on disk, and this window has unsaved changes."));
      bar->show();
    });
    connect(reload, &QAction::triggered, this, reloadProject);
    connect(keep, &QAction::triggered, this, [watcher, bar]() {
      watcher->acceptDisk();
      bar->hide();
    });
    connect(differences, &QAction::triggered, this, [this, watcher]() {
      QString why;
      commitTabs(why);
      QDialog *d = new QDialog(this);
      d->setObjectName("diskDifferencesDialog");
      d->setAttribute(Qt::WA_DeleteOnClose);
      d->setWindowTitle(tr("The project on disk and in this window"));
      QPlainTextEdit *text = new QPlainTextEdit(ProjectWatcher::lineDiff(
          QString::fromUtf8(watcher->diskText()), QString::fromStdString(holder_->project().Text()), tr("on disk"),
          tr("in this window")));
      text->setReadOnly(true);
      text->setLineWrapMode(QPlainTextEdit::NoWrap);
      text->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
      (new QVBoxLayout(d))->addWidget(text);
      d->resize(900, 600);
      d->show();
    });
    auto refreshFromDisk = [this, overview]() {
      if (running_) return;
      // the plot keeps a zoom the user made (a view other than the whole range)
      PlotWidget *plot = plotTab->plot();
      const QRectF before = plot->view();
      plot->resetZoom();
      const bool zoomed = before != plot->view();
      plotTab->refresh();
      if (zoomed) plot->zoomTo(before);
      runTab->setDataSetNames(plotTab->dataRecords());
      overview->refresh();
    };
    connect(watcher, &ProjectWatcher::outputsChanged, this, [this, refreshFromDisk]() {
      refreshFromDisk();
      statusBar()->showMessage(tr("The outputs changed on disk: the Overview and the plot read again."), 15000);
    });
    connect(watcher, &ProjectWatcher::inputsChanged, this, [this, refreshFromDisk](const QStringList &files) {
      refreshFromDisk();
      statusBar()->showMessage(tr("Data changed on disk: %1").arg(files.join(", ")), 15000);
    });
    connect(watcher, &ProjectWatcher::parameterFileChanged, overview, &OverviewPage::refresh);
    connect(holder_, &ProjectHolder::opened, watcher, &ProjectWatcher::follow);
    connect(holder_, &ProjectHolder::sectionEdited, watcher, &ProjectWatcher::follow);
    for (const char *name : {"save", "saveAs"}) connect(actions_[name], &QAction::triggered, watcher, &ProjectWatcher::follow);
    connect(this, &AZURESetup::startWorker, watcher, [watcher](const QString &, const QStringList &options) {
      QString file;
      for (const QString &o : options)
        if (o.startsWith("--parameter-file=")) file = o.mid(17);
      watcher->setParameterFile(file);
      watcher->setRunning(true);
    });
    connect(this, &AZURESetup::runFinished, watcher, [watcher]() { watcher->setRunning(false); });
    QAction *followDisk = new QAction(tr("&Follow Changes on Disk"), this);
    followDisk->setObjectName("followDisk");
    followDisk->setCheckable(true);
    followDisk->setToolTip(tr("Re-read the project, its data and the run's outputs when another program changes them"));
    followDisk->setChecked(QSettings().value("view/followDisk", true).toBool());
    watcher->setEnabled(followDisk->isChecked());
    connect(followDisk, &QAction::toggled, watcher, [watcher](bool on) {
      QSettings().setValue("view/followDisk", on);
      watcher->setEnabled(on);
    });
    for (QAction *menu : menuBar()->actions())
      if (menu->menu() && menu->text() == tr("&View")) menu->menu()->addAction(followDisk);
  }
  projectOpened();
  showProject(false);
  restoreWindowState();
  setExpertMode(expert::Enabled());  // the view of the last session (standard the first time)
}

AZURESetup::~AZURESetup() {
  if (worker_) worker_->requestStop();
  runThread_.quit();
  runThread_.wait();
}

void AZURESetup::createActions() {
  expertAction_ = new QAction(tr("Show &Expert Options"), this);
  expertAction_->setObjectName("expertAction");
  expertAction_->setCheckable(true);
  expertAction_->setChecked(expert::Enabled());
  expertAction_->setToolTip(tr("Show the options, columns, tabs and sections an occasional user does not need.  It "
                               "changes what is shown, never the project: an expert option the project sets is shown "
                               "either way."));
  connect(expertAction_, &QAction::toggled, this, &AZURESetup::setExpertMode);
  // the main actions: the tool bar's and the menus' (one QAction each, so a shortcut has one owner)
  auto make = [this](const QString &name, const QString &icon, const QString &text, const QKeySequence &key,
                     const QString &tip) {
    QAction *a = new QAction(QIcon(":/actions/" + icon + ".svg"), text, this);
    a->setObjectName("action_" + name);
    a->setShortcut(key);
    a->setToolTip(QString("%1 (%2)").arg(tip, key.toString(QKeySequence::NativeText)));
    a->setStatusTip(tip);
    actions_.insert(name, a);
    return a;
  };
  connect(make("new", "new", tr("&New Project"), QKeySequence::New, tr("Start an empty project")), &QAction::triggered,
          this, &AZURESetup::newProject);
  connect(make("open", "open", tr("&Open..."), QKeySequence::Open, tr("Open a project (.azr)")), &QAction::triggered,
          this, qOverload<>(&AZURESetup::open));
  connect(make("save", "save", tr("&Save"), QKeySequence::Save, tr("Save the project")), &QAction::triggered, this,
          &AZURESetup::save);
  connect(make("validate", "validate", tr("&Validate"), QKeySequence(Qt::Key_F7),
               tr("Check the project with the engine (AZURE2 --validate)")),
          &QAction::triggered, this, &AZURESetup::validate);
  auto runFromToolBar = [this](const QString &command) {
    QString why;
    if (!runCommand(command, why)) QMessageBox::warning(this, tr("Not run"), why);
  };
  connect(make("calculate", "calculate", tr("&Calculate"), QKeySequence(Qt::Key_F5),
               tr("Save, then calculate the data sets from the parameters (no fit)")),
          &QAction::triggered, this, [runFromToolBar]() { runFromToolBar("calculate"); });
  connect(make("fit", "fit", tr("&Fit"), QKeySequence(Qt::Key_F6), tr("Save, then fit the parameters to the data")),
          &QAction::triggered, this, [runFromToolBar]() { runFromToolBar("fit"); });
  connect(make("stop", "stop", tr("S&top"), QKeySequence(Qt::SHIFT | Qt::Key_F5),
               tr("Stop the run (a fit ends at its best point and writes its outputs there)")),
          &QAction::triggered, this, &AZURESetup::stopRun);
  connect(make("plot", "plot", tr("&Plot"), QKeySequence(Qt::Key_F8), tr("Show the last run's curves")),
          &QAction::triggered, this, [this]() { tabWidget->setCurrentWidget(plotTab); });
  connect(make("close", "close", tr("&Close Project"), QKeySequence::Close, tr("Close the project")), &QAction::triggered,
          this, [this]() { closeProject(); });
  for (const QString &name : {"calculate", "fit", "stop"}) actions_[name]->setShortcutContext(Qt::WindowShortcut);
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
  fileMenu->addAction(actions_["new"]);
  fileMenu->addAction(actions_["open"]);
  recentFileMenu = fileMenu->addMenu(tr("Open &Recent..."));
  for (int i = 0; i < numRecent; i++) recentFileMenu->addAction(recentFileActions[i]);
  recentSeparator = recentFileMenu->addSeparator();
  recentFileMenu->addAction(clearRecentAction);
  updateRecent();
  fileMenu->addAction(actions_["save"]);
  QAction *saveAsAction = fileMenu->addAction(tr("Save &As..."), QKeySequence::SaveAs, this, &AZURESetup::saveAs);
  actions_.insert("saveAs", saveAsAction);
  fileMenu->addAction(tr("Con&vert an Old Project..."), this, &AZURESetup::convertProject);
  fileMenu->addAction(actions_["close"]);
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
  configMenu->addAction(actions_["validate"]);

  QMenu *runMenu = menuBar()->addMenu(tr("&Run"));
  for (const QString &name : {"calculate", "fit", "stop"}) runMenu->addAction(actions_[name]);
  runMenu->addSeparator();
  runMenu->addAction(actions_["plot"]);

  // View: the look (AzureTheme.h), the tool bar and the Diagnostics
  QMenu *viewMenu = menuBar()->addMenu(tr("&View"));
  QMenu *lookMenu = viewMenu->addMenu(tr("&Appearance"));
  QActionGroup *looks = new QActionGroup(this);
  for (const ThemeChoice choice : {ThemeChoice::System, ThemeChoice::Light, ThemeChoice::Dark}) {
    QAction *a = lookMenu->addAction(choice == ThemeChoice::System  ? tr("Follow the &System")
                                     : choice == ThemeChoice::Light ? tr("&Light")
                                                                    : tr("&Dark"));
    a->setObjectName("theme_" + themeChoiceName(choice));
    a->setCheckable(true);
    a->setChecked(azureThemeChoice() == choice);
    looks->addAction(a);
    connect(a, &QAction::triggered, this, [choice]() { setAzureTheme(choice); });
  }
  viewMenu->addAction(toolBar_->toggleViewAction());
  QAction *dockAction = diagnosticsDock->toggleViewAction();
  dockAction->setShortcut(QKeySequence(Qt::CTRL | Qt::SHIFT | Qt::Key_D));
  viewMenu->addAction(dockAction);

  // Tools: pyazr's commands as child processes (PythonTools.h), filled when the menu opens (pyazr importable or not)
  // standard or expert view (qt_gui.md section 19.4)
  viewMenu->addSeparator();
  viewMenu->addAction(expertAction_);

  QMenu *toolsMenu = menuBar()->addMenu(tr("&Tools"));
  toolsMenu->setObjectName("toolsMenu");
  toolsMenu_ = toolsMenu;
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
  if (startPage_) startPage_->setRecent(files);
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
  // [*] is Qt's place of the modified mark in the title
  setWindowTitle(hasProject_ ? tr("%1[*] -- AZURE2").arg(QFileInfo(holder_->path()).fileName()) : tr("AZURE2"));
  setWindowModified(hasProject_ && modified);
  if (statusFile_) {
    statusFile_->setText(hasProject_ ? QString("%1%2").arg(QDir::toNativeSeparators(holder_->path()),
                                                           modified ? tr("  (modified)") : QString())
                                     : tr("No project open"));
    statusFile_->setToolTip(hasProject_ ? (modified ? tr("Changed since it was saved") : tr("Saved")) : QString());
  }
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
  validateRequested_ = true;
  diagnostics->validate();
}

void AZURESetup::diagnosticsDone(bool ok, int errors, int warnings) {
  if (statusValidity_) {
    QString text;
    if (!ok && errors > 0) text = tr("%n error(s)", nullptr, errors);
    else text = tr("Valid");
    if (warnings > 0) text += tr(", %n warning(s)", nullptr, warnings);
    statusValidity_->setText(text);
    statusValidity_->setIcon(style()->standardIcon(!ok ? QStyle::SP_MessageBoxCritical
                                                   : warnings > 0 ? QStyle::SP_MessageBoxWarning
                                                                  : QStyle::SP_DialogApplyButton));
    statusValidity_->setToolTip(diagnostics->summary() + "\n" + tr("Click to show the Diagnostics."));
    statusValidity_->setVisible(hasProject_);
  }
  if (!hasProject_) return;
  const bool asked = validateRequested_;
  validateRequested_ = false;
  // The check of a project just opened is one note of the notes bar, never a box (qt_gui.md section 19.2).
  if (openValidation_) {
    openValidation_ = false;
    if (!ok || errors > 0 || warnings > 0)
      notes_->addNote(tr("The engine's check: %1").arg(statusValidity_ ? statusValidity_->text() : diagnostics->summary()),
                      diagnostics->summary() + "\n\n" + diagnostics->listing(30), true);
  }
  // An error (the project would not run) opens the Diagnostics; warnings open them only when the user asked for the
  // check (Validate): otherwise they are counted in the status bar, whose click opens them.  Nothing to show: hidden.
  if (!ok || errors > 0 || (asked && warnings > 0)) diagnosticsDock->setVisible(true);
  else if (warnings == 0) diagnosticsDock->setVisible(false);
}

void AZURESetup::showRecord(const QString &section, const QString &record) {
  showProject(true);
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

QAction *AZURESetup::action(const QString &name) const { return actions_.value(name, nullptr); }

void AZURESetup::createToolBar() {
  toolBar_ = addToolBar(tr("Main"));
  toolBar_->setObjectName("mainToolBar");
  toolBar_->setMovable(false);
  toolBar_->setIconSize(QSize(22, 22));
  toolBar_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
  for (const QString &name : {"new", "open", "save"}) toolBar_->addAction(actions_[name]);
  toolBar_->addSeparator();
  toolBar_->addAction(actions_["validate"]);
  toolBar_->addSeparator();
  for (const QString &name : {"calculate", "fit", "stop"}) toolBar_->addAction(actions_[name]);
  toolBar_->addSeparator();
  toolBar_->addAction(actions_["plot"]);
  // the menu's words without the mnemonic
  for (QAction *a : toolBar_->actions())
    if (QWidget *w = toolBar_->widgetForAction(a)) w->setObjectName("tool_" + a->objectName().mid(7));
}

void AZURESetup::createStatusBar() {
  statusFile_ = new QLabel;
  statusFile_->setObjectName("statusFile");
  statusFile_->setTextInteractionFlags(Qt::TextSelectableByMouse);
  statusValidity_ = new QToolButton;
  statusValidity_->setObjectName("statusValidity");
  statusValidity_->setAutoRaise(true);
  statusValidity_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
  statusValidity_->setIconSize(QSize(16, 16));
  statusValidity_->setText(tr("Not validated"));
  connect(statusValidity_, &QToolButton::clicked, this, [this]() {
    diagnosticsDock->show();
    diagnosticsDock->raise();
  });
  jobLabel_ = new QLabel;
  jobLabel_->setObjectName("statusJob");
  jobBar_ = new QProgressBar;
  jobBar_->setObjectName("statusProgress");
  jobBar_->setMaximumWidth(fontMetrics().horizontalAdvance("M") * 14);
  jobBar_->setTextVisible(false);
  jobStop_ = new QToolButton;
  jobStop_->setObjectName("statusStop");
  jobStop_->setDefaultAction(actions_["stop"]);
  jobStop_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
  statusExpert_ = new QToolButton;
  statusExpert_->setObjectName("statusExpert");
  statusExpert_->setAutoRaise(true);
  statusExpert_->setToolTip(expertAction_->toolTip());
  connect(statusExpert_, &QToolButton::clicked, expertAction_, &QAction::toggle);
  jobBox_ = new QWidget;
  QHBoxLayout *job = new QHBoxLayout(jobBox_);
  job->setContentsMargins(0, 0, 0, 0);
  job->addWidget(jobLabel_);
  job->addWidget(jobBar_);
  job->addWidget(jobStop_);
  jobBox_->hide();
  statusBar()->addWidget(statusFile_, 1);
  statusBar()->addPermanentWidget(jobBox_);
  statusBar()->addPermanentWidget(statusValidity_);
  statusBar()->addPermanentWidget(statusExpert_);
  statusValidity_->hide();
  updateTitle();
}

bool AZURESetup::expertMode() const { return expertAction_->isChecked(); }

void AZURESetup::setExpertMode(bool on) {
  if (expertAction_->isChecked() != on) {
    expertAction_->setChecked(on);  // toggled() comes back here
    return;
  }
  expert::SetEnabled(on);
  statusExpert_->setText(on ? tr("Expert view") : tr("Standard view"));
  if (startPage_) startPage_->setExpertChecked(on);
  pairsTab->setExpertMode(on);
  levelsTab->setExpertMode(on);
  for (SectionsTab *t : {segmentsTab, effectsTab, potentialTab}) t->setExpertMode(on);
  fittingTab->setExpertMode(on);
  runTab->setExpertMode(on);
  thmTab->setExpertMode(on);
  static_cast<PhysicsTab *>(physicsTab_)->setExpertMode(on);
  toolsMenu_->menuAction()->setVisible(on);
  updateExpertTabs();
}

void AZURESetup::updateExpertTabs() {
  // a section other than comments is in use: its tab is shown in either view (nothing hidden is in effect)
  auto inUse = [this](const QString &section) {
    for (const QString &line : holder_->sectionBody(section)) {
      const QString t = line.trimmed();
      if (!t.isEmpty() && !t.startsWith('#')) return true;
    }
    return false;
  };
  const bool on = expertMode();
  const QList<QPair<QWidget *, bool>> tabs = {{potentialTab, potentialTab->hasRecords()},
                                              {thmTab, inUse("thm")},
                                              {mcmcTab, mcmcTab && inUse("mcmc")}};
  for (const auto &t : tabs) {
    if (!t.first) continue;
    const int at = tabWidget->indexOf(t.first);
    tabWidget->setTabVisible(at, on || t.second);
    tabWidget->setTabToolTip(at, tabTips_.value(t.first) + (!on && t.second ? "\n" + expert::ShownBecauseSet() : QString()));
  }
  for (auto it = tabTips_.constBegin(); it != tabTips_.constEnd(); ++it) {
    bool expertTab = false;
    for (const auto &t : tabs) expertTab = expertTab || t.first == it.key();
    if (!expertTab) tabWidget->setTabToolTip(tabWidget->indexOf(it.key()), it.value());
  }
  static_cast<PhysicsTab *>(physicsTab_)->setExpertMode(on);
}

void AZURESetup::updateActions() {
  const bool idle = !running_;
  actions_["new"]->setEnabled(idle);
  actions_["open"]->setEnabled(idle);
  for (const QString &name : {"save", "saveAs", "validate", "close"})
    if (QAction *a = actions_.value(name)) a->setEnabled(hasProject_ && idle);
  actions_["calculate"]->setEnabled(hasProject_ && idle);
  actions_["fit"]->setEnabled(hasProject_ && idle);
  actions_["stop"]->setEnabled(running_);
  actions_["plot"]->setEnabled(hasProject_);
}

void AZURESetup::showProject(bool open) {
  hasProject_ = open;
  central_->setCurrentWidget(open ? static_cast<QWidget *>(tabWidget) : static_cast<QWidget *>(startPage_));
  if (!open) {
    diagnosticsDock->hide();
    QList<QPair<QString, QString>> notes;
    const QMap<QWidget *, QString> what = {
        {pairsTab, tr("The pairs of nuclei: masses, spins, parities, charges, separation energies, channel radii.")},
        {levelsTab, tr("The levels of the compound nucleus (J, parity, energy) and the widths of their channels.")},
        {segmentsTab, tr("The data sets to fit and the extrapolations to calculate, one record each.")},
        {effectsTab, tr("Targets, beams and detectors: energy loss, straggling, resolution, the folds of the data.")},
        {potentialTab, tr("Potentials and the other format-3 sections of the nuclear part.")},
        {thmTab, tr("The Trojan Horse experiments and the options of their vertex.")},
        {fittingTab, tr("The minimizer, the priors, and every fit parameter with its limits.")},
        {runTab, tr("Run calculate, fit, extrapolate, MINOS or the rate; the log and the chi-squared of each data set.")},
        {plotTab, tr("The data and the calculation of the last run, the residuals; export as PNG, SVG, PDF or CSV.")}};
    for (int i = 0; i < tabWidget->count(); i++) {
      QWidget *tab = tabWidget->widget(i);
      const QString text = what.contains(tab) ? what[tab]
                           : tab == mcmcTab   ? tr("Sample the posterior of the parameters (the ensemble sampler).")
                           : tab->objectName() == "physicsTab"
                               ? tr("Physics functions and stopping without a project; the project's channels, R-matrix, "
                                    "pole map and THM parts.")
                               : QString();
      if (!text.isEmpty()) notes << qMakePair(tabWidget->tabText(i).remove('&'), text);
      if (!text.isEmpty()) tabTips_[tab] = text;  // what the tab is for: its tool tip
    }
    startPage_->setTabNotes(notes);
    updateExpertTabs();
  }
  updateTitle();
  updateActions();
  if (statusValidity_) statusValidity_->setVisible(open);
  // the last validation (projectOpened's), now that a project is shown: the dock opens on a problem
  if (open) diagnosticsDone(diagnostics->errorCount() == 0, diagnostics->errorCount(), diagnostics->warningCount());
}

bool AZURESetup::closeProject() {
  if (running_ || !maybeSave()) return false;
  notes_->clearNotes();
  holder_->setProject(azure2::front::Project::Empty(QDir::current().filePath("untitled.azr").toStdString()));
  showProject(false);
  return true;
}

bool AZURESetup::runCommand(const QString &command, QString &why) {
  if (!hasProject_) {
    why = tr("no project is open");
    return false;
  }
  if (!runTab->setCommand(command)) {
    why = tr("this build's engine has no command %1").arg(command);
    return false;
  }
  tabWidget->setCurrentWidget(runTab);
  return startRun(why);
}

void AZURESetup::showJobProgress(const QString &stage, qlonglong evaluation, qlonglong calls, double objective,
                                 double best) {
  if (!running_) return;
  if (stage == "mcmc") {
    jobBar_->setRange(0, (int)std::min<qlonglong>(calls, 1000000000));
    jobBar_->setValue((int)std::min<qlonglong>(evaluation, 1000000000));
    jobLabel_->setText(tr("mcmc: step %1 of %2").arg(evaluation).arg(calls));
    return;
  }
  QString text = tr("%1: evaluation %2, chi-squared %3").arg(stage).arg(evaluation).arg(objective, 0, 'g', 7);
  if (std::isfinite(best)) text += tr(" (best %1)").arg(best, 0, 'g', 7);
  jobLabel_->setText(text);
}

bool AZURESetup::acceptsDrop(const QList<QUrl> &urls) {
  for (const QUrl &u : urls)
    if (u.isLocalFile() && u.toLocalFile().endsWith(".azr", Qt::CaseInsensitive)) return true;
  return false;
}

void AZURESetup::dragEnterEvent(QDragEnterEvent *event) {
  if (!running_ && event->mimeData()->hasUrls() && acceptsDrop(event->mimeData()->urls())) event->acceptProposedAction();
}

void AZURESetup::dropEvent(QDropEvent *event) {
  if (openDropped(event->mimeData()->urls())) event->acceptProposedAction();
}

bool AZURESetup::openDropped(const QList<QUrl> &urls) {
  if (running_) return false;
  for (const QUrl &u : urls)
    if (u.isLocalFile() && u.toLocalFile().endsWith(".azr", Qt::CaseInsensitive)) {
      if (maybeSave()) open(u.toLocalFile());
      return true;
    }
  return false;
}

void AZURESetup::restoreWindowState() {
  QSettings settings;
  if (settings.contains("window/geometry")) restoreGeometry(settings.value("window/geometry").toByteArray());
  else resize(1280, 820);
  // the dock's visibility is the validation's (diagnosticsDone), not the last session's
  if (settings.contains("window/state")) restoreState(settings.value("window/state").toByteArray());
  diagnosticsDock->hide();
  if (notes_) notes_->clearNotes();  // the notes are those of a project opened in this session
}

void AZURESetup::saveWindowState() {
  QSettings settings;
  settings.setValue("window/geometry", saveGeometry());
  settings.setValue("window/state", saveState());
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
  saveWindowState();
  event->accept();
}
