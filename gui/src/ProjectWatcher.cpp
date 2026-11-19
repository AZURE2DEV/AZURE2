// Changes on disk to the open project's files (ProjectWatcher.h; docs/upgrade/qt_gui.md section 16).
#include "ProjectWatcher.h"

#include <QCryptographicHash>
#include <algorithm>
#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QFileSystemWatcher>
#include <QTimer>
#include <vector>

#include "ProjectHolder.h"

namespace {

QByteArray Hash(const QByteArray &bytes) { return QCryptographicHash::hash(bytes, QCryptographicHash::Sha1); }

/// Size and time of a file ("-" when it is not there): what the quiet time compares.
QByteArray Stamp(const QString &path) {
  const QFileInfo f(path);
  if (!f.exists()) return QByteArray("-");
  return QByteArray::number(f.size()) + ":" + QByteArray::number(f.lastModified().toMSecsSinceEpoch());
}

/// The output files of a run the window shows (the Overview, the Plot tab, the sampler's page).
bool IsRunOutput(const QString &name) {
  return name == "chiSquared.out" || name == "fitStatus.out" || name == "param.sav" || name == "samples.mcmc" ||
         name.startsWith("AZUREOut_");
}

}  // namespace

ProjectWatcher::ProjectWatcher(ProjectHolder *holder, QObject *parent) : QObject(parent), holder_(holder) {
  watcher_ = new QFileSystemWatcher(this);
  quiet_ = new QTimer(this);
  quiet_->setSingleShot(true);
  connect(quiet_, &QTimer::timeout, this, &ProjectWatcher::settle);
  connect(watcher_, &QFileSystemWatcher::fileChanged, this, &ProjectWatcher::touched);
  connect(watcher_, &QFileSystemWatcher::directoryChanged, this, &ProjectWatcher::touched);
}

void ProjectWatcher::follow() {
  const QString previousOutputs = outputDir_, previousProject = projectFile_;
  projectFile_ = QFileInfo(holder_->path()).absoluteFilePath();
  if (!QFileInfo::exists(projectFile_)) projectFile_.clear();
  outputDir_.clear();
  inputs_.clear();
  inputSeen_.clear();
  if (!projectFile_.isEmpty()) {
    const azure2::front::RunSummary run = azure2::front::ReadLastRun(holder_->project());
    outputDir_ = QDir::cleanPath(QString::fromStdString(run.outputDirectory));
    for (const azure2::front::InputFile &f : azure2::front::InputFiles(holder_->project())) {
      const QString path = QFileInfo(QString::fromStdString(f.path)).absoluteFilePath();
      if (!inputs_.contains(path)) inputs_ << path;
    }
  }
  projectHash_ = Hash(diskText());
  // an outside run's outputs not yet reported stay so when only the records changed
  if (outputDir_ != previousOutputs || projectFile_ != previousProject) outputs_ = outputSignature();
  for (const QString &f : inputs_) inputSeen_[f] = inputSignature(f);
  parameterSeen_ = inputSignature(parameterFile_);
  lastProbe_.clear();
  rewatch();
}

void ProjectWatcher::setParameterFile(const QString &file) {
  const QString path = file.isEmpty() ? QString() : QFileInfo(QDir(holder_->directory()), file).absoluteFilePath();
  if (path == parameterFile_) return;
  parameterFile_ = path;
  parameterSeen_ = inputSignature(parameterFile_);
  rewatch();
}

void ProjectWatcher::setEnabled(bool on) {
  enabled_ = on;
  if (on) follow();
  else {
    quiet_->stop();
    rewatch();
  }
}

void ProjectWatcher::setRunning(bool running) {
  running_ = running;
  // the window refreshes itself at the run's end: the outputs as they then stand are seen
  if (!running) {
    outputs_ = outputSignature();
    rewatch();
  }
}

void ProjectWatcher::acceptDisk() { projectHash_ = Hash(diskText()); }

QByteArray ProjectWatcher::diskText() const {
  QFile f(projectFile_);
  if (projectFile_.isEmpty() || !f.open(QIODevice::ReadOnly)) return QByteArray();
  return f.readAll();
}

QStringList ProjectWatcher::watched() const { return watcher_->files() + watcher_->directories(); }

void ProjectWatcher::rewatch() {
  const QStringList old = watched();
  if (!old.isEmpty()) watcher_->removePaths(old);
  if (!enabled_ || projectFile_.isEmpty()) return;
  QStringList paths = {projectFile_, QFileInfo(projectFile_).absolutePath()};
  if (!outputDir_.isEmpty() && QFileInfo(outputDir_).isDir()) {
    paths << outputDir_;
    for (const QFileInfo &f : QDir(outputDir_).entryInfoList(QDir::Files))
      if (IsRunOutput(f.fileName())) paths << f.absoluteFilePath();
  }
  paths << inputs_;
  if (!parameterFile_.isEmpty()) paths << parameterFile_;
  QStringList existing;
  for (const QString &p : paths)
    if (QFileInfo::exists(p) && !existing.contains(p)) existing << p;
  watcher_->addPaths(existing);
}

void ProjectWatcher::touched(const QString &path) {
  if (!enabled_) return;
  // the window's own save: the bytes on disk are its text at once (the notification comes before another edit)
  if (path == projectFile_) {
    const QByteArray disk = Hash(diskText());
    if (disk == Hash(QByteArray::fromStdString(holder_->project().Text()))) projectHash_ = disk;
  }
  pending_ = true;
  quiet_->start(quietMs);
}

QByteArray ProjectWatcher::outputSignature() const {
  QByteArray s;
  if (outputDir_.isEmpty()) return s;
  for (const QFileInfo &f : QDir(outputDir_).entryInfoList(QDir::Files, QDir::Name))
    if (IsRunOutput(f.fileName()))
      s += f.fileName().toUtf8() + "=" + QByteArray::number(f.size()) + ":" +
           QByteArray::number(f.lastModified().toMSecsSinceEpoch()) + ";";
  return s;
}

QByteArray ProjectWatcher::inputSignature(const QString &file) const { return file.isEmpty() ? QByteArray() : Stamp(file); }

void ProjectWatcher::settle() {
  if (!enabled_ || !pending_) return;
  // still being written: wait until every watched file has kept its size and time for the quiet time
  QByteArray probe = Stamp(projectFile_) + outputSignature() + Stamp(parameterFile_);
  for (const QString &f : inputs_) probe += Stamp(f);
  if (probe != lastProbe_) {
    lastProbe_ = probe;
    quiet_->start(quietMs);
    return;
  }
  pending_ = false;
  lastProbe_.clear();

  const QByteArray disk = Hash(diskText());
  const bool exists = QFileInfo::exists(projectFile_);
  if (exists && disk != projectHash_) {
    projectHash_ = disk;
    if (disk != Hash(QByteArray::fromStdString(holder_->project().Text()))) emit projectChanged();
  }
  if (!running_) {
    const QByteArray outputs = outputSignature();
    if (outputs != outputs_) {
      outputs_ = outputs;
      emit outputsChanged();
    }
  }
  QStringList changed;
  for (const QString &f : inputs_) {
    const QByteArray s = inputSignature(f);
    if (s != inputSeen_.value(f)) {
      inputSeen_[f] = s;
      changed << f;
    }
  }
  if (!changed.isEmpty()) emit inputsChanged(changed);
  const QByteArray parameter = inputSignature(parameterFile_);
  if (parameter != parameterSeen_) {
    parameterSeen_ = parameter;
    emit parameterFileChanged();
  }
  rewatch();  // a file replaced (written and renamed) has lost its watch
}

QString ProjectWatcher::lineDiff(const QString &a, const QString &b, const QString &nameA, const QString &nameB) {
  const QStringList x = a.split('\n'), y = b.split('\n');
  QStringList out = {"--- " + nameA, "+++ " + nameB};
  const int n = x.size(), m = y.size();
  if ((qint64)n * m > 25000000) {  // too large for the table: the two texts one after the other
    for (const QString &l : x) out << "- " + l;
    for (const QString &l : y) out << "+ " + l;
    return out.join('\n');
  }
  // the longest common subsequence of lines, from the end
  std::vector<std::vector<int>> lcs(n + 1, std::vector<int>(m + 1, 0));
  for (int i = n - 1; i >= 0; i--)
    for (int j = m - 1; j >= 0; j--)
      lcs[i][j] = x[i] == y[j] ? lcs[i + 1][j + 1] + 1 : std::max(lcs[i + 1][j], lcs[i][j + 1]);
  QStringList all;
  int i = 0, j = 0;
  while (i < n || j < m) {
    if (i < n && j < m && x[i] == y[j]) all << "  " + x[i++], j++;
    else if (i < n && (j == m || lcs[i + 1][j] >= lcs[i][j + 1])) all << "- " + x[i++];  // the old line first
    else all << "+ " + y[j++];
  }
  // the changed lines with three lines of context each; "..." where lines are left out
  const int context = 3;
  std::vector<bool> keep(all.size(), false);
  for (int k = 0; k < all.size(); k++)
    if (!all[k].startsWith("  "))
      for (int c = std::max(0, k - context); c <= std::min((int)all.size() - 1, k + context); c++) keep[c] = true;
  bool gap = false;
  for (int k = 0; k < all.size(); k++) {
    if (keep[k]) {
      if (gap) out << "...";
      out << all[k];
      gap = false;
    } else {
      gap = true;
    }
  }
  if (gap) out << "...";
  return out.join('\n');
}
