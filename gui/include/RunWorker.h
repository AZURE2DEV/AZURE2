#ifndef RUNWORKER_H
#define RUNWORKER_H

#include <atomic>

#include <QList>
#include <QMetaType>
#include <QObject>
#include <QPair>
#include <QString>
#include <QStringList>

/// What a run gave back (azure2::front::RunResult in Qt types), carried by a queued signal to the window.
struct RunOutcome {
  int exitCode = 0;
  bool stopped = false;
  bool haveChiSquared = false;
  double totalChiSquared = 0;
  QList<QPair<int, double>> dataSetChiSquared;  ///< (segment, chi^2)
  QList<int> dataSetPoints;
  QString outputDirectory;
  QStringList outputFiles;
};
Q_DECLARE_METATYPE(RunOutcome)

/*!
 * The worker of a run (docs/upgrade/qt_gui.md section 3): lives in its own QThread and runs one command at a time
 * through azure2::front::Run, the program's run.  Its signals carry the program's lines, the fit's progress and the
 * outcome to the window by queued connections; the window touches the worker only through run() (a queued call) and
 * requestStop(), which sets the run's atomic stop request.  The project is saved before the run starts and the run
 * reads the file: nothing of the window is shared with the run's threads.
 */
class RunWorker : public QObject {
  Q_OBJECT

 public:
  RunWorker() = default;
  /// Thread-safe: the run in progress ends at its next check (a fit at its best point).
  void requestStop() { stop_.store(true); }
  /// Before a run is queued (the window's thread): a stop asked before the run starts is not lost.
  void clearStop() { stop_.store(false); }

 public slots:
  void run(const QString &command, const QStringList &options, const QString &project, int threads);

 signals:
  void message(const QString &line, bool transient);
  void progress(const QString &stage, qlonglong evaluation, double objective, double best);
  void finished(const RunOutcome &outcome);

 private:
  std::atomic<bool> stop_{false};
};

#endif
