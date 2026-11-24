// The window's project file and runs (docs/upgrade/qt_gui.md sections 2 and 3): open (a format-3 file, or an old one
// through the converter), the tabs' state written into the project, save, and the run in the worker thread.
#include <QDir>
#include <QFileDialog>
#include <QFileInfo>
#include <QMessageBox>
#include <QProgressBar>
#include <QLabel>
#include <QStatusBar>

#include "AZURESetup.h"
#include "ConvertRepairDialog.h"
#include "DiagnosticsPanel.h"
#include "FittingTab.h"
#include "LevelsTab.h"
#include "PairsTab.h"
#include "RecordTabs.h"
#include "MCMCTab.h"
#include "NotesBar.h"
#include "PlotTab.h"
#include "RunTab.h"
#include "ThmTab.h"

void AZURESetup::newProject() {
  if (running_ || !maybeSave()) return;
  holder_->setProject(azure2::front::Project::Empty(QDir::current().filePath("untitled.azr").toStdString()));
  showProject(true);
}

void AZURESetup::open() {
  if (running_ || !maybeSave()) return;
  const QString filename = QFileDialog::getOpenFileName(this, tr("Open a project"), QString(),
                                                        tr("AZURE2 projects (*.azr);;All files (*)"));
  if (!filename.isEmpty()) open(filename);
}

void AZURESetup::open(QString filename) {
  if (filename.isEmpty() || running_) return;
  // the start page's lists may ask twice for the file just opened (a click that also activates)
  if (hasProject_ && QFileInfo(filename).absoluteFilePath() == QFileInfo(holder_->path()).absoluteFilePath() &&
      !holder_->project().Modified() && !levelsDirty_)
    return;
  QString why;
  if (readProjectFile(filename, why)) return;
  if (why.contains("--convert")) {
    // an old project: the converter of AZURE2 --convert writes a format-3 copy beside it (nothing is overwritten)
    const QFileInfo info(filename);
    const QString suggested = info.dir().filePath(info.completeBaseName() + ".f3.azr");
    if (QMessageBox::question(this, tr("An older project"),
                              tr("%1\n\nConvert it to project format 3 now (the old file is kept)?").arg(why)) !=
        QMessageBox::Yes)
      return;
    const QString to = QFileDialog::getSaveFileName(this, tr("The converted project"), suggested,
                                                    tr("AZURE2 projects (*.azr)"));
    if (to.isEmpty()) return;
    QString report;
    if (!convertOldProject(filename, to, report) && !report.isEmpty()) QMessageBox::warning(this, tr("Not converted"), report);
    return;
  }
  QMessageBox::warning(this, tr("Can't open the project"), why);
}

void AZURESetup::convertProject() {
  const QString from = QFileDialog::getOpenFileName(this, tr("An old project (format 1 or 2)"), QString(),
                                                    tr("AZURE2 projects (*.azr);;All files (*)"));
  if (from.isEmpty()) return;
  const QFileInfo info(from);
  const QString to = QFileDialog::getSaveFileName(this, tr("The converted project"),
                                                  info.dir().filePath(info.completeBaseName() + ".f3.azr"),
                                                  tr("AZURE2 projects (*.azr)"));
  if (to.isEmpty()) return;
  if (!maybeSave()) return;
  QString report;
  if (!convertOldProject(from, to, report) && !report.isEmpty()) QMessageBox::warning(this, tr("Not converted"), report);
}

bool AZURESetup::convertOldProject(const QString &from, const QString &to, QString &report) {
  report.clear();
  QString source = from;  // the old file, then the copy with the lines the user repaired
  const QFileInfo info(from);
  QString copy = info.dir().filePath(info.completeBaseName() + ".repaired." + info.suffix());
  for (int n = 2; QFileInfo::exists(copy); n++)  // a copy of an earlier repair is kept too
    copy = info.dir().filePath(info.completeBaseName() + QString(".repaired%1.").arg(n) + info.suffix());
  QStringList repaired;
  for (int round = 0; round < 50; round++) {
    std::string text;
    azure2::front::ConvertLineProblem problem;
    if (azure2::front::ConvertProject(source.toStdString(), to.toStdString(), text, problem)) {
      report.clear();
      QString why;
      if (!readProjectFile(to, why)) {
        report = why;
        return false;
      }
      // what the converter assumed and what the user repaired: notes of the notes bar, not boxes
      if (!repaired.isEmpty())
        notes_->addNote(tr("%n line(s) of %1 repaired before the conversion (in the copy %2)", nullptr, repaired.size())
                            .arg(info.fileName(), QFileInfo(copy).fileName()),
                        repaired.join("\n"));
      if (!text.empty())
        notes_->addNote(tr("Converted from %1 to project format 3").arg(info.fileName()), QString::fromStdString(text));
      return true;
    }
    report = QString::fromStdString(text);
    if (problem.line <= 0) return false;  // a refusal about no single line: the report says why
    ConvertRepairDialog dialog(source, problem, this);
    if (dialog.exec() != QDialog::Accepted) {
      report.clear();  // cancelled: nothing to say
      return false;
    }
    // the copy holds every repair so far; the old file is never written
    const QString next = copy + ".new";
    QFile::remove(next);
    std::string why;
    if (!azure2::front::WriteRepairedCopy(source.toStdString(), problem.line, dialog.lineText().toStdString(),
                                          next.toStdString(), why)) {
      report = QString::fromStdString(why);
      return false;
    }
    if (source == copy) QFile::remove(copy);
    if (QFileInfo::exists(copy) || !QFile::rename(next, copy)) {
      report = tr("%1 exists: nothing is overwritten").arg(copy);
      QFile::remove(next);
      return false;
    }
    repaired << tr("line %1: \"%2\" -> \"%3\"").arg(problem.line).arg(QString::fromStdString(problem.text).simplified(),
                                                       dialog.lineText());
    source = copy;
  }
  report = tr("too many repairs");
  return false;
}

bool AZURESetup::readProjectFile(const QString &filename, QString &why) {
  azure2::front::Project project;
  std::string reason;
  if (!azure2::front::Project::Open(QFileInfo(filename).absoluteFilePath().toStdString(), project, reason)) {
    why = QString::fromStdString(reason);
    return false;
  }
  holder_->setProject(std::move(project));  // projectOpened() fills the tabs
  addRecent(QFileInfo(filename).absoluteFilePath());
  showProject(true);
  return true;
}

bool AZURESetup::loadLevels(QString &why) {
  loadingLevels_ = true;
  levelsTab->reset();
  pairsTab->reset();
  const bool ok = levelsTab->readLevels(holder_->project().Levels(), why);
  loadedIncludeAllJPi_ = holder_->configValue("includeAllJPi", "0") == "1" || holder_->configValue("includeAllJPi") == "true";
  levelsTab->setIncludeAllJPi(loadedIncludeAllJPi_);
  const QString maxL = holder_->configValue("maxOrbitalMomentum");
  if (!maxL.isEmpty()) levelsTab->setMaxL(maxL.toInt());
  loadedMaxL_ = levelsTab->maxL();
  loadingLevels_ = false;
  levelsDirty_ = false;
  return ok;
}

void AZURESetup::projectOpened() {
  // the notes of this project only, shown once (the bar is cleared when another project opens)
  notes_->clearNotes();
  QString why;
  if (!loadLevels(why)) notes_->addNote(tr("The Levels tab could not read <levels>: %1").arg(why.section('\n', 0, 0)), why);
  openValidation_ = true;
  updateExpertTabs();
  levelsTab->setExpertMode(expertMode());  // the photon multipolarity rows follow the pairs just read
  const QString amatrix = holder_->configValue("amatrix", "true");
  (amatrix == "false" ? rMatrixAction : aMatrixAction)->setChecked(true);
  runTab->clearLog();
  updateTitle();
  updateLogo();
  diagnostics->validate();
}

bool AZURESetup::commitTabs(QString &why) {
  if (levelsDirty_) {
    std::string reason;
    if (!holder_->project().SetLevels(levelsTab->levelRows(), reason)) {
      why = QString::fromStdString(reason);
      return false;
    }
    levelsDirty_ = false;
  }
  // the Levels tab's hard-sphere settings: written only when changed
  if (levelsTab->includeAllJPi() != loadedIncludeAllJPi_) {
    if (!holder_->setConfig("includeAllJPi", levelsTab->includeAllJPi() ? "1" : "0", why)) return false;
    loadedIncludeAllJPi_ = levelsTab->includeAllJPi();
  }
  if (levelsTab->maxL() != loadedMaxL_) {
    if (!holder_->setConfig("maxOrbitalMomentum", QString::number(levelsTab->maxL()), why)) return false;
    loadedMaxL_ = levelsTab->maxL();
  }
  thmTab->applyText();
  if (!fittingTab->apply(why)) return false;
  if (mcmcTab) mcmcTab->apply();
  updateTitle();
  return true;
}

bool AZURESetup::writeProjectFile(const QString &filename, QString &why) {
  if (!commitTabs(why)) return false;
  std::string reason;
  if (!holder_->project().Save(filename.toStdString(), reason)) {
    why = QString::fromStdString(reason);
    return false;
  }
  addRecent(holder_->path());
  updateTitle();
  return true;
}

void AZURESetup::save() {
  if (QFileInfo(holder_->path()).fileName() == "untitled.azr" || !QFileInfo::exists(holder_->path())) {
    saveAs();
    return;
  }
  QString why;
  if (!writeProjectFile(holder_->path(), why)) QMessageBox::warning(this, tr("Not saved"), why);
}

void AZURESetup::saveAs() {
  const QString filename = QFileDialog::getSaveFileName(this, tr("Save the project"), holder_->path(),
                                                        tr("AZURE2 projects (*.azr)"));
  if (filename.isEmpty()) return;
  QString why;
  if (!writeProjectFile(filename, why)) QMessageBox::warning(this, tr("Not saved"), why);
}

bool AZURESetup::maybeSave() {
  if (!hasProject_) return true;  // the start page's empty project holds nothing of the user's
  if (!holder_->project().Modified() && !levelsDirty_) return true;
  const QMessageBox::StandardButton answer =
      QMessageBox::question(this, tr("Unsaved changes"), tr("Save the changes to %1?").arg(holder_->path()),
                            QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
  if (answer == QMessageBox::Cancel) return false;
  if (answer == QMessageBox::Save) {
    save();
    return !holder_->project().Modified() && !levelsDirty_;
  }
  return true;
}

void AZURESetup::SaveAndRun() {
  QString why;
  if (!startRun(why)) QMessageBox::warning(this, tr("Not run"), why);
}

bool AZURESetup::startRun(QString &why) {
  if (running_) {
    why = tr("a run is in progress");
    return false;
  }
  // the run reads the file, as the program does: the project is saved first
  if (QFileInfo(holder_->path()).fileName() == "untitled.azr") {
    saveAs();
    if (QFileInfo(holder_->path()).fileName() == "untitled.azr") {
      why = tr("the project has no file yet");
      return false;
    }
  }
  if (!writeProjectFile(holder_->path(), why)) return false;
  // the output and checks directories the run needs (the program refuses a run without them)
  for (const char *key : {"output_dir", "checks_dir"}) {
    const QString dir = holder_->configValue(key);
    if (!dir.isEmpty()) QDir(holder_->directory()).mkpath(dir);
  }
  const QString command = runTab->command();
  const QStringList options = runTab->options() + fittingTab->programOptions(command);
  running_ = true;
  updateActions();
  jobLabel_->setText(tr("%1: starting").arg(command));
  jobBar_->setRange(0, 0);  // busy until the run tells its steps (the sampler)
  jobBox_->show();
  runTab->setRunning(true);
  if (mcmcTab) mcmcTab->setRunning(true);
  runTab->clearLog();
  runTab->appendLine(QString("AZURE2 %1 %2 %3").arg(command, options.join(' '), holder_->path()), false);
  worker_->clearStop();
  emit startWorker(command, options, holder_->path(), runTab->threads());
  return true;
}

void AZURESetup::stopRun() {
  if (running_) worker_->requestStop();
}

bool AZURESetup::startMCMC(QString &why) {
  if (!mcmcTab || !runTab->setCommand("mcmc")) {
    why = tr("this build's engine has no mcmc command");
    return false;
  }
  for (const QString &o : mcmcTab->options()) {
    const int eq = o.indexOf('=');
    if (eq < 0) runTab->setOption(o, "1");
    else runTab->setOption(o.left(eq + 1), o.mid(eq + 1));
  }
  if (!mcmcTab->options().contains("--overwrite-samples")) runTab->setOption("--overwrite-samples", "0");
  return startRun(why);
}

void AZURESetup::runEnded(const RunOutcome &outcome) {
  running_ = false;
  updateActions();
  jobBox_->hide();
  QString end = outcome.exitCode == 0 ? tr("The run ended") : tr("The run failed (exit %1)").arg(outcome.exitCode);
  if (outcome.stopped) end += tr(", stopped on request");
  if (outcome.haveChiSquared) end += tr(": total chi-squared %1").arg(outcome.totalChiSquared, 0, 'g', 10);
  statusBar()->showMessage(end + ".", 30000);
  runTab->setRunning(false);
  if (mcmcTab) mcmcTab->setRunning(false);
  plotTab->refresh();
  runTab->setDataSetNames(plotTab->dataRecords());
  runTab->showOutcome(outcome);
  emit runFinished(outcome);
}
