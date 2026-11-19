#include "RunWorker.h"

#include <chrono>
#include <cmath>
#include <mutex>

#include "azure2/front/Engine.h"

namespace {

/// The facade's observer: each line and progress record becomes a signal of the worker (queued to the window, as the
/// worker lives in another thread than the window's widgets).  Progress records come from several threads at once:
/// they are thinned to one every 50 ms, under a lock, plus every new best.
class SignalObserver : public azure2::front::RunObserver {
 public:
  explicit SignalObserver(RunWorker *worker) : worker_(worker) {}
  void Message(const std::string &line, bool transient) override {
    emit worker_->message(QString::fromStdString(line), transient);
  }
  void Progress(const azure2::front::RunProgress &p) override {
    std::lock_guard<std::mutex> lock(mutex_);
    const auto now = std::chrono::steady_clock::now();
    const bool sampler = p.stage == "mcmc";  // a sampler's step: its last always arrives
    const bool better = !sampler && !std::isnan(p.best) && (std::isnan(best_) || p.best < best_);
    const bool last = sampler && p.evaluation == p.calls;
    if (!better && !last && now - last_ < std::chrono::milliseconds(50)) return;
    last_ = now;
    if (better) best_ = p.best;
    emit worker_->progress(QString::fromStdString(p.stage), p.evaluation, p.calls, p.objective, p.best);
  }

 private:
  RunWorker *worker_;
  std::mutex mutex_;
  std::chrono::steady_clock::time_point last_{};
  double best_ = NAN;
};

}  // namespace

void RunWorker::run(const QString &command, const QStringList &options, const QString &project, int threads) {
  azure2::front::RunRequest request;
  request.command = command.toStdString();
  for (const QString &o : options) request.options.push_back(o.toStdString());
  request.project = project.toStdString();
  request.threads = threads;
  SignalObserver observer(this);
  const azure2::front::RunResult r = azure2::front::Run(request, observer, &stop_);
  RunOutcome out;
  out.exitCode = r.exitCode;
  out.stopped = r.stopped;
  out.haveChiSquared = r.haveChiSquared;
  out.totalChiSquared = r.totalChiSquared;
  for (const azure2::front::DataSetChiSquared &d : r.dataSets) {
    out.dataSetChiSquared << qMakePair(d.segment, d.chiSquared);
    out.dataSetPoints << d.points;
  }
  out.outputDirectory = QString::fromStdString(r.outputDirectory);
  for (const std::string &f : r.outputFiles) out.outputFiles << QString::fromStdString(f);
  emit finished(out);
}
