#ifndef PROJECTWATCHER_H
#define PROJECTWATCHER_H

#include <QByteArray>
#include <QMap>
#include <QObject>
#include <QStringList>

QT_BEGIN_NAMESPACE
class QFileSystemWatcher;
class QTimer;
QT_END_NAMESPACE

class ProjectHolder;

/*!
 * Changes on disk to the open project's files (docs/upgrade/qt_gui.md section 16): the project file, the input files
 * its records name (the facade's InputFiles: data files, covariance files), the parameter file of the next run, and
 * the output files of a run (chiSquared.out, fitStatus.out, AZUREOut_*, param.sav, samples.mcmc, through the output
 * directory).  QFileSystemWatcher with each path added again after an atomic replace (a file written to a temporary
 * name and renamed loses its watch), and a quiet time: a change is reported once its files have kept their size and
 * time for `quietMs` (a file being written is never read half-way).
 *
 * The project file is compared by content, not by time: a change whose bytes are the window's own text (its own save)
 * or the bytes last seen is no change.  The watcher decides nothing about the window: it reports, and the window
 * reloads, asks or refreshes.  While the window's own run is in progress (setRunning) the outputs are not reported; the
 * run's end refreshes the window, and the watcher takes the outputs as they then stand as seen.
 */
class ProjectWatcher : public QObject {
  Q_OBJECT

 public:
  explicit ProjectWatcher(ProjectHolder *holder, QObject *parent = nullptr);

  /// Watch the holder's project as it now stands (its file, inputs and output directory); the project file's bytes
  /// on disk become the ones seen.  Called after an open, a save and an edit of the records.
  void follow();
  /// The parameter file a run will read ("" none).
  void setParameterFile(const QString &file);
  /// Following on or off (View > Follow changes on disk); off: nothing is watched.
  void setEnabled(bool on);
  bool enabled() const { return enabled_; }
  void setRunning(bool running);
  /// The project file's bytes on disk are what the window decided to keep (Keep mine): not reported again until they
  /// change again.
  void acceptDisk();
  /// The project file as it is on disk now.
  QByteArray diskText() const;
  /// The paths watched (for the tests).
  QStringList watched() const;

  int quietMs = 400;

  /// A plain line diff of two texts: "  " same, "- " only in `a`, "+ " only in `b`, with headers.
  static QString lineDiff(const QString &a, const QString &b, const QString &nameA, const QString &nameB);

 signals:
  /// The project file changed on disk to bytes that are neither the window's text nor the bytes last seen.
  void projectChanged();
  /// Output files of a run changed (an outside run), after the quiet time.
  void outputsChanged();
  /// Input files of the project changed (the data of a record): their paths.
  void inputsChanged(const QStringList &files);
  /// The parameter file of the next run changed.
  void parameterFileChanged();

 private:
  void rewatch();
  void touched(const QString &path);
  void settle();
  QByteArray outputSignature() const;
  QByteArray inputSignature(const QString &file) const;

  ProjectHolder *holder_;
  QFileSystemWatcher *watcher_;
  QTimer *quiet_;
  bool enabled_ = true;
  bool running_ = false;
  QString projectFile_, outputDir_, parameterFile_;
  QStringList inputs_;
  QByteArray projectHash_;    ///< the project file's bytes last seen (sha1)
  QByteArray outputs_;        ///< the outputs' signature last reported or seen
  QMap<QString, QByteArray> inputSeen_;
  QByteArray parameterSeen_;
  bool pending_ = false;
  QByteArray lastProbe_;      ///< the signature of every watched file at the last probe of the quiet time
};

#endif
