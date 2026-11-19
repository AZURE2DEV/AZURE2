#ifndef PHYSICSWORKER_H
#define PHYSICSWORKER_H

#include <atomic>
#include <functional>
#include <memory>
#include <string>

#include <QHash>
#include <QObject>
#include <QString>
#include <QThread>

#include "azure2/front/Engine.h"

/// The engine's session of a saved project as the worker thread holds it (only that thread touches it).
struct PhysicsSessionState {
  azure2::front::PhysicsSession session;
  std::string key;            ///< what it was opened on (file, its time stamp, the parameter file)
  bool extrapolations = false;
  /// Open the session on `file` at `parameterFile` ("" start values) unless it is open on `key` already; a project
  /// the engine cannot set up with its data sets is tried with its extrapolations (a project without data).  False
  /// with the engine's words in `why`.
  bool Ensure(const std::string &key, const std::string &file, const std::string &parameterFile, std::string &why);
};

/*!
 * A worker thread of the Physics tab (docs/upgrade/qt_gui.md section 15.3).  A job is a function run in the worker
 * thread, one at a time, in the order submitted; when it ends, its `done` runs in the GUI thread.  A job submitted on
 * a lane supersedes the lane's earlier job: that one is told to stop (its cancel flag) and its `done` is dropped.
 * cancel() stops every job; a job reads its flag between chunks of its work.  The tab has two workers: one holds the
 * project's PhysicsSession (the facade's rule: one thread at a time per session), one the session-free functions.
 */
class PhysicsWorker : public QObject {
  Q_OBJECT

 public:
  using Job = std::function<void(PhysicsSessionState &state, const std::atomic<bool> &cancel)>;
  explicit PhysicsWorker(QObject *parent = nullptr);
  ~PhysicsWorker() override;

  /// Run `job` in the worker thread, then `done` in the GUI thread (unless cancelled or superseded on `lane`).
  void submit(const QString &lane, Job job, std::function<void()> done);
  /// Stop every job and drop their results.
  void cancel();
  /// Jobs submitted and not yet ended.
  bool busy() const { return pending_ > 0; }

 signals:
  void busyChanged(bool busy);

 private:
  QThread thread_;
  QObject *runner_ = nullptr;  ///< lives in thread_; jobs are queued calls on it
  std::shared_ptr<PhysicsSessionState> state_;
  QHash<QString, std::shared_ptr<std::atomic<bool>>> lanes_;  ///< the cancel flag of each lane's last job
  int pending_ = 0;
};

#endif
