// The Coulomb table of this computer in the window (CoulombTableUi.h, docs/upgrade/qt_gui.md section 19.12).
#include "CoulombTableUi.h"

#include <QDialogButtonBox>
#include <QDir>
#include <QHBoxLayout>
#include <QLabel>
#include <QLocale>
#include <QProgressBar>
#include <QPushButton>
#include <QThread>
#include <QTimer>
#include <QVBoxLayout>
#include <algorithm>
#include <atomic>
#include <string>

namespace front = azure2::front;

namespace {

QString Native(const std::string &path) { return QDir::toNativeSeparators(QString::fromStdString(path)); }

/// The line's words for a state the runs cannot use.
QString NoticeText(const front::CoulombTableState &s) {
  if (s.exists)
    return QObject::tr("Coulomb table not used: %1. Calculations run slower; build it again for this computer "
                       "(about 5 min, 88 MB).")
        .arg(QString::fromStdString(s.problem));
  return QObject::tr("Coulomb table not built: calculations run slower. Build it once for this computer "
                     "(about 5 min, 88 MB).");
}

}  // namespace

// ---- the notice ----

CoulombTableNotice::CoulombTableNotice(QWidget *parent) : QWidget(parent) {
  setObjectName("coulombTableNotice");
  QHBoxLayout *row = new QHBoxLayout(this);
  row->setContentsMargins(0, 0, 0, 0);
  text_ = new QLabel;
  text_->setObjectName("coulombTableNoticeText");
  text_->setProperty("role", "hint");  // quiet: the placeholder colour, never red
  text_->setWordWrap(true);
  build_ = new QPushButton(tr("Build..."));
  build_->setObjectName("coulombTableNoticeBuild");
  connect(build_, &QPushButton::clicked, this, &CoulombTableNotice::buildRequested);
  row->addWidget(text_, 1);
  row->addWidget(build_, 0, Qt::AlignVCenter);
  hide();
}

void CoulombTableNotice::setState(const front::CoulombTableState &state, bool buildAllowed, const QString &why) {
  text_->setText(NoticeText(state));
  text_->setToolTip(tr("Runs look for the table at %1 (or at the file AZURE2_COULOMB_TABLE names). Without it they "
                       "compute the exact Coulomb functions: the same results, about 4 times slower.")
                        .arg(Native(state.path)));
  build_->setEnabled(buildAllowed);
  build_->setToolTip(buildAllowed ? tr("Open the Coulomb table's window: build it there, with a progress bar and Cancel.")
                                  : why);
  setVisible(!state.usable);
}

// ---- the dialog ----

struct CoulombTableDialog::Job {
  front::CoulombTableBuild spec;
  std::atomic<int> done{0};
  std::atomic<bool> cancel{false}, finished{false};
  bool ok = false;
  std::string error;
  int rows = 1;
};

CoulombTableDialog::CoulombTableDialog(QWidget *parent) : QDialog(parent) {
  setObjectName("coulombTableDialog");
  setWindowTitle(tr("The Coulomb Table"));
  setModal(false);
  spec_ = front::DefaultCoulombTableBuild();
  QVBoxLayout *box = new QVBoxLayout(this);
  QLabel *what = new QLabel(tr("The Coulomb table holds the penetrability, shift function and hard-sphere phase of "
                               "every charged channel over (eta, rho, l). Runs that find it compute the same results "
                               "(within 1e-7) about 4 times faster. It is built once per computer: about 5 minutes on "
                               "4 threads, 88 MB."));
  what->setWordWrap(true);
  box->addWidget(what);
  text_ = new QLabel;
  text_->setObjectName("coulombTableState");
  text_->setWordWrap(true);
  text_->setTextInteractionFlags(Qt::TextSelectableByMouse);
  box->addWidget(text_);
  bar_ = new QProgressBar;
  bar_->setObjectName("coulombTableProgress");
  bar_->setRange(0, 100);
  bar_->setValue(0);
  bar_->hide();
  box->addWidget(bar_);
  result_ = new QLabel;
  result_->setObjectName("coulombTableResult");
  result_->setWordWrap(true);
  result_->hide();
  box->addWidget(result_);
  QDialogButtonBox *buttons = new QDialogButtonBox;
  build_ = buttons->addButton(tr("Build"), QDialogButtonBox::ActionRole);
  build_->setObjectName("coulombTableBuildButton");
  cancel_ = buttons->addButton(tr("Cancel Build"), QDialogButtonBox::ActionRole);
  cancel_->setObjectName("coulombTableCancelButton");
  close_ = buttons->addButton(QDialogButtonBox::Close);
  box->addWidget(buttons);
  connect(build_, &QPushButton::clicked, this, [this]() { startBuild(); });
  connect(cancel_, &QPushButton::clicked, this, &CoulombTableDialog::cancelBuild);
  connect(close_, &QPushButton::clicked, this, &QDialog::hide);
  timer_ = new QTimer(this);
  timer_->setInterval(200);
  connect(timer_, &QTimer::timeout, this, &CoulombTableDialog::poll);
  resize(560, sizeHint().height());
  refresh();
}

CoulombTableDialog::~CoulombTableDialog() {
  // a build still running when the window closes: stopped, and nothing left at the path
  if (worker_) {
    job_->cancel = true;
    worker_->wait();
    delete worker_;
  }
}

int CoulombTableDialog::buildThreads() const {
  const int wanted = spec_.threads > 0 ? spec_.threads : front::DefaultCoulombTableBuild().threads;
  return std::max(1, std::min(wanted, 4));
}

void CoulombTableDialog::refresh() {
  state_ = front::CoulombTableStatus();
  QString text = tr("Runs look for it at: %1").arg(Native(state_.path));
  if (state_.usable)
    text += "\n" + tr("Found and used: l up to %1, %2 nodes per decade of eta and rho, eta*rho up to %3 (%4).")
                       .arg(state_.lMax)
                       .arg(state_.nodesPerDecade, 0, 'g', 4)
                       .arg(state_.etaRhoMax, 0, 'g', 4)
                       .arg(QLocale::c().formattedDataSize(state_.bytes));
  else if (state_.exists)
    text += "\n" + tr("A file is there but runs do not use it: %1").arg(QString::fromStdString(state_.problem));
  else
    text += "\n" + tr("None there yet: runs compute the exact Coulomb functions.");
  text_->setText(text);
  updateButtons();
  emit stateChanged(state_);
}

void CoulombTableDialog::setCalculationRunning(bool running) {
  calculationRunning_ = running;
  updateButtons();
}

void CoulombTableDialog::updateButtons() {
  build_->setEnabled(!building() && !calculationRunning_);
  build_->setText(state_.usable || state_.exists ? tr("Build Again") : tr("Build"));
  build_->setToolTip(calculationRunning_ ? tr("A calculation is running: build the table when it has ended.")
                                         : tr("Build the table at the path above on %1 thread(s).").arg(buildThreads()));
  cancel_->setEnabled(building());
  cancel_->setVisible(building());
}

bool CoulombTableDialog::startBuild(QString *why) {
  QString reason;
  if (building()) reason = tr("A build of the Coulomb table is already running.");
  else if (calculationRunning_) reason = tr("A calculation is running: build the Coulomb table when it has ended.");
  if (!reason.isEmpty()) {
    if (why) *why = reason;
    return false;
  }
  job_ = std::make_shared<Job>();
  job_->spec = spec_;
  if (job_->spec.path.empty()) job_->spec.path = state_.path;
  job_->spec.threads = buildThreads();
  job_->rows = std::max(1, front::CoulombTableBuildRows(job_->spec));
  std::shared_ptr<Job> job = job_;
  worker_ = QThread::create([job]() {
    job->ok = front::BuildCoulombTable(job->spec, &job->done, &job->cancel, job->error);
    job->finished = true;
  });
  bar_->setValue(0);
  bar_->show();
  result_->setText(tr("Building %1 on %2 thread(s)...").arg(Native(job_->spec.path)).arg(job_->spec.threads));
  result_->show();
  worker_->start();
  timer_->start();
  updateButtons();
  emit buildStarted();
  return true;
}

void CoulombTableDialog::cancelBuild() {
  if (job_ && worker_) job_->cancel = true;
}

void CoulombTableDialog::poll() {
  if (!job_) return;
  bar_->setValue(int(100.0 * job_->done.load() / job_->rows));
  if (job_->finished) finish();
}

void CoulombTableDialog::finish() {
  timer_->stop();
  worker_->wait();
  delete worker_;
  worker_ = nullptr;
  const bool ok = job_->ok;
  const QString error = QString::fromStdString(job_->error);
  if (ok) {
    bar_->setValue(100);
    result_->setText(tr("Built: the next run uses it."));
  } else {
    bar_->hide();
    result_->setText(error == "cancelled" ? tr("Cancelled: nothing was written.") : tr("Not built: %1").arg(error));
  }
  job_.reset();
  refresh();
  emit buildFinished(ok, error);
}

// ---- the window's place for the table ----

CoulombTableControl::CoulombTableControl(QWidget *window) : QObject(window), window_(window) {}

CoulombTableNotice *CoulombTableControl::makeNotice() {
  CoulombTableNotice *n = new CoulombTableNotice;
  connect(n, &CoulombTableNotice::buildRequested, this, [this]() { showDialog(); });
  notices_ << n;
  updateNotices();
  return n;
}

void CoulombTableControl::refresh() {
  if (dialog_) {
    dialog_->refresh();  // stateChanged comes back to updateNotices
    return;
  }
  state_ = front::CoulombTableStatus();
  updateNotices();
}

void CoulombTableControl::setCalculationRunning(bool running) {
  running_ = running;
  if (dialog_) dialog_->setCalculationRunning(running);
  updateNotices();
}

CoulombTableDialog *CoulombTableControl::showDialog() {
  if (!dialog_) {
    dialog_ = new CoulombTableDialog(window_);
    dialog_->setCalculationRunning(running_);
    connect(dialog_, &CoulombTableDialog::stateChanged, this, [this](const front::CoulombTableState &s) {
      state_ = s;
      updateNotices();
    });
    connect(dialog_, &CoulombTableDialog::buildStarted, this, &CoulombTableControl::updateNotices);
    state_ = dialog_->state();
  } else {
    dialog_->refresh();
  }
  dialog_->show();
  dialog_->raise();
  dialog_->activateWindow();
  return dialog_;
}

void CoulombTableControl::updateNotices() {
  const bool building = dialog_ && dialog_->building();
  const QString why = building ? tr("The Coulomb table is being built (its window shows the progress).")
                               : tr("A calculation is running: build the Coulomb table when it has ended.");
  notices_.removeAll(nullptr);
  for (CoulombTableNotice *n : notices_) n->setState(state_, !running_ && !building, why);
}
