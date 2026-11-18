#include "RunTab.h"

#include <cmath>

#include <QCheckBox>
#include <QComboBox>
#include <QDesktopServices>
#include <QDir>
#include <QFileDialog>
#include <QFormLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QSplitter>
#include <QTableWidget>
#include <QThread>
#include <QToolButton>
#include <QUrl>
#include <QVBoxLayout>

#include "azure2/front/Engine.h"

namespace {

/// The old program's words for the commands (its Calculation Type list), the engine's name otherwise.
QString Title(const QString &command) {
  if (command == "calculate") return QObject::tr("Calculate Segments From Data");
  if (command == "fit") return QObject::tr("Fit Segments From Data");
  if (command == "extrapolate") return QObject::tr("Calculate Segments Without Data");
  if (command == "minos") return QObject::tr("Perform MINOS Error Analysis");
  if (command == "rate") return QObject::tr("Calculate Reaction Rate");
  if (command == "mcmc") return QObject::tr("MCMC Bayesian Sampling");
  return command;
}

}  // namespace

RunTab::RunTab(QWidget *parent) : QWidget(parent) {
  command_ = new QComboBox;
  command_->setObjectName("calcType");
  for (const azure2::front::Command &c : azure2::front::Commands()) {
    if (!c.built) continue;
    command_->addItem(Title(QString::fromStdString(c.name)), QString::fromStdString(c.name));
    command_->setItemData(command_->count() - 1, QString::fromStdString(c.help), Qt::ToolTipRole);
  }
  connect(command_, &QComboBox::currentIndexChanged, this, &RunTab::commandChanged);
  threads_ = new QSpinBox;
  threads_->setObjectName("threads");
  threads_->setRange(1, std::max(1, QThread::idealThreadCount()));
  threads_->setValue(std::min(4, threads_->maximum()));
  threads_->setToolTip(tr("The OpenMP threads of the run (--threads); results do not depend on it."));
  runButton_ = new QPushButton(tr("Save and &Run"));
  runButton_->setObjectName("runButton");
  stopButton_ = new QPushButton(tr("Stop AZURE2"));
  stopButton_->setObjectName("stopButton");
  stopButton_->setEnabled(false);
  stopButton_->setToolTip(tr("A fit ends at its best point and writes its outputs there."));
  connect(runButton_, &QPushButton::clicked, this, &RunTab::runRequested);
  connect(stopButton_, &QPushButton::clicked, this, &RunTab::stopRequested);

  QGridLayout *calcLayout = new QGridLayout;
  calcLayout->addWidget(new QLabel(tr("Calculation Type:")), 0, 0, Qt::AlignRight);
  calcLayout->addWidget(command_, 0, 1);
  calcLayout->addWidget(new QLabel(tr("Threads:")), 0, 2, Qt::AlignRight);
  calcLayout->addWidget(threads_, 0, 3);
  calcLayout->setColumnStretch(4, 1);
  calcLayout->addWidget(runButton_, 0, 5);
  calcLayout->addWidget(stopButton_, 0, 6);

  optionsBox_ = new QGroupBox(tr("Options of the calculation"));
  optionsForm_ = new QFormLayout(optionsBox_);

  progress_ = new QLabel;
  progress_->setObjectName("runProgress");
  log_ = new QPlainTextEdit;
  log_->setObjectName("runLog");
  log_->setReadOnly(true);
  log_->setMaximumBlockCount(20000);
  log_->setFont(QFont("Monospace"));
  chi_ = new QTableWidget(0, 3);
  chi_->setObjectName("chiTable");
  chi_->setHorizontalHeaderLabels({tr("Segment"), tr("Chi-squared"), tr("Chi-squared/N")});
  chi_->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
  chi_->verticalHeader()->setVisible(false);
  chi_->setEditTriggers(QAbstractItemView::NoEditTriggers);
  files_ = new QListWidget;
  files_->setObjectName("outputFiles");
  files_->setToolTip(tr("The files the run wrote in the output directory; double-click opens one."));
  connect(files_, &QListWidget::itemDoubleClicked, this, [this](QListWidgetItem *item) {
    QDesktopServices::openUrl(QUrl::fromLocalFile(QDir(outputDirectory_).filePath(item->text())));
  });
  QSplitter *results = new QSplitter(Qt::Vertical);
  results->addWidget(chi_);
  results->addWidget(files_);
  QSplitter *bottom = new QSplitter(Qt::Horizontal);
  bottom->addWidget(log_);
  bottom->addWidget(results);
  bottom->setStretchFactor(0, 3);
  bottom->setStretchFactor(1, 1);

  QVBoxLayout *layout = new QVBoxLayout(this);
  layout->addLayout(calcLayout);
  layout->addWidget(optionsBox_);
  layout->addWidget(progress_);
  layout->addWidget(bottom, 1);
  commandChanged();
}

QString RunTab::command() const { return command_->currentData().toString(); }

bool RunTab::setCommand(const QString &command) {
  const int i = command_->findData(command);
  if (i < 0) return false;
  command_->setCurrentIndex(i);
  return true;
}

int RunTab::threads() const { return threads_->value(); }

void RunTab::commandChanged() {
  // keep what was typed, by option
  for (auto it = editors_.constBegin(); it != editors_.constEnd(); ++it) {
    if (QCheckBox *b = qobject_cast<QCheckBox *>(it.value())) kept_[it.key()] = b->isChecked() ? "1" : "0";
    else if (QLineEdit *e = it.value()->findChild<QLineEdit *>(); e) kept_[it.key()] = e->text();
    else if (QLineEdit *e2 = qobject_cast<QLineEdit *>(it.value())) kept_[it.key()] = e2->text();
  }
  while (optionsForm_->rowCount() > 0) optionsForm_->removeRow(0);
  editors_.clear();
  const QString command = this->command();
  for (const azure2::front::Option &o : azure2::front::Options()) {
    bool applies = o.commands.empty();
    for (const std::string &c : o.commands) applies = applies || QString::fromStdString(c) == command;
    if (!applies || o.option == "--max-calls=") continue;  // the call limit is the Fitting Settings tab's
    const QString option = QString::fromStdString(o.option);
    const QString help = QString::fromStdString(o.help);
    QString label = option;
    label.remove('=');
    QWidget *editor = nullptr;
    if (o.kind == "flag") {
      QCheckBox *box = new QCheckBox(help.left(90));
      box->setChecked(kept_.value(option) == "1");
      editor = box;
    } else {
      QLineEdit *edit = new QLineEdit(kept_.value(option));
      edit->setPlaceholderText(o.required ? tr("required") : tr("default %1").arg(o.defaultValue));
      if (o.kind == "file") {
        QWidget *host = new QWidget;
        QHBoxLayout *h = new QHBoxLayout(host);
        h->setContentsMargins(0, 0, 0, 0);
        h->addWidget(edit, 1);
        QToolButton *browse = new QToolButton;
        browse->setText("...");
        connect(browse, &QToolButton::clicked, this, [this, edit]() {
          const QString f = QFileDialog::getOpenFileName(this, tr("Choose a file"));
          if (!f.isEmpty()) edit->setText(f);
        });
        h->addWidget(browse);
        editor = host;
      } else {
        editor = edit;
      }
    }
    editor->setObjectName("option_" + label.mid(2));
    editor->setToolTip(help);
    optionsForm_->addRow(label + ":", editor);
    editors_.insert(option, editor);
  }
  optionsBox_->setVisible(!editors_.isEmpty());
}

QStringList RunTab::options() const {
  QStringList out;
  for (auto it = editors_.constBegin(); it != editors_.constEnd(); ++it) {
    if (QCheckBox *b = qobject_cast<QCheckBox *>(it.value())) {
      if (b->isChecked()) out << it.key();
      continue;
    }
    QLineEdit *e = qobject_cast<QLineEdit *>(it.value());
    if (!e) e = it.value()->findChild<QLineEdit *>();
    if (e && !e->text().trimmed().isEmpty()) out << it.key() + e->text().trimmed();
  }
  return out;
}

bool RunTab::setOption(const QString &option, const QString &value) {
  QWidget *w = editors_.value(option, nullptr);
  if (!w) return false;
  if (QCheckBox *b = qobject_cast<QCheckBox *>(w)) {
    b->setChecked(value == "1");
    return true;
  }
  QLineEdit *e = qobject_cast<QLineEdit *>(w);
  if (!e) e = w->findChild<QLineEdit *>();
  if (!e) return false;
  e->setText(value);
  return true;
}

void RunTab::setRunning(bool running) {
  runButton_->setEnabled(!running);
  stopButton_->setEnabled(running);
  command_->setEnabled(!running);
  optionsBox_->setEnabled(!running);
  threads_->setEnabled(!running);
}

QString RunTab::log() const { return log_->toPlainText(); }

void RunTab::clearLog() {
  log_->clear();
  progress_->clear();
  chi_->setRowCount(0);
  files_->clear();
}

void RunTab::appendLine(const QString &line, bool transient) {
  if (transient) {
    progress_->setText(line.trimmed());
    return;
  }
  log_->appendPlainText(line);
}

void RunTab::showProgress(const QString &stage, qlonglong evaluation, double objective, double best) {
  QString text = tr("%1: evaluation %2, chi-squared %3").arg(stage).arg(evaluation).arg(objective, 0, 'g', 10);
  if (!std::isnan(best)) text += tr(" (best %1)").arg(best, 0, 'g', 10);
  progress_->setText(text);
}

void RunTab::showOutcome(const RunOutcome &outcome) {
  outputDirectory_ = outcome.outputDirectory;
  chi_->setRowCount(0);
  for (int i = 0; i < outcome.dataSetChiSquared.size(); i++) {
    const int row = chi_->rowCount();
    chi_->insertRow(row);
    const int n = outcome.dataSetPoints.value(i);
    chi_->setItem(row, 0, new QTableWidgetItem(QString::number(outcome.dataSetChiSquared[i].first)));
    chi_->setItem(row, 1, new QTableWidgetItem(QString::number(outcome.dataSetChiSquared[i].second, 'g', 10)));
    chi_->setItem(row, 2, new QTableWidgetItem(n > 0 ? QString::number(outcome.dataSetChiSquared[i].second / n, 'g', 6)
                                                     : QString()));
  }
  if (outcome.haveChiSquared) {
    const int row = chi_->rowCount();
    chi_->insertRow(row);
    chi_->setItem(row, 0, new QTableWidgetItem(tr("Total")));
    chi_->setItem(row, 1, new QTableWidgetItem(QString::number(outcome.totalChiSquared, 'g', 12)));
  }
  files_->clear();
  files_->addItems(outcome.outputFiles);
  QString end = outcome.exitCode == 0 ? tr("The run ended.") : tr("The run failed (exit %1).").arg(outcome.exitCode);
  if (outcome.stopped) end += tr("  It was stopped on request.");
  progress_->setText(end);
}
