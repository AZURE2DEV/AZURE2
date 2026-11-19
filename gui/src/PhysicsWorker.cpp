// The Physics tab's worker threads (PhysicsWorker.h, docs/upgrade/qt_gui.md section 15.3).
#include "PhysicsWorker.h"

#include <QMetaObject>

bool PhysicsSessionState::Ensure(const std::string &wanted, const std::string &file, const std::string &parameterFile,
                                 std::string &why) {
  if (session.IsOpen() && key == wanted) return true;
  key.clear();
  azure2::front::PhysicsSession fresh;
  std::string first;
  extrapolations = false;
  if (!azure2::front::PhysicsSession::Open(file, parameterFile, false, fresh, first)) {
    std::string second;
    if (!azure2::front::PhysicsSession::Open(file, parameterFile, true, fresh, second)) {
      why = first.empty() ? second : first;
      session = azure2::front::PhysicsSession();
      return false;
    }
    extrapolations = true;
  }
  session = std::move(fresh);
  key = wanted;
  return true;
}

PhysicsWorker::PhysicsWorker(QObject *parent) : QObject(parent), state_(std::make_shared<PhysicsSessionState>()) {
  runner_ = new QObject;
  runner_->moveToThread(&thread_);
  connect(&thread_, &QThread::finished, runner_, &QObject::deleteLater);
  thread_.start();
}

PhysicsWorker::~PhysicsWorker() {
  cancel();
  thread_.quit();
  thread_.wait();
}

void PhysicsWorker::submit(const QString &lane, Job job, std::function<void()> done) {
  if (lanes_.contains(lane)) *lanes_[lane] = true;  // the lane's earlier job: stop, and its result is dropped
  auto flag = std::make_shared<std::atomic<bool>>(false);
  lanes_[lane] = flag;
  if (pending_++ == 0) emit busyChanged(true);
  std::shared_ptr<PhysicsSessionState> state = state_;
  QMetaObject::invokeMethod(
      runner_,
      [this, state, flag, job, done]() {
        if (!*flag) job(*state, *flag);
        // back in the GUI thread; dropped by Qt if this worker is gone
        QMetaObject::invokeMethod(
            this,
            [this, flag, done]() {
              if (!*flag && done) done();
              if (--pending_ == 0) emit busyChanged(false);
            },
            Qt::QueuedConnection);
      },
      Qt::QueuedConnection);
}

void PhysicsWorker::cancel() {
  for (auto &flag : lanes_) *flag = true;
}
