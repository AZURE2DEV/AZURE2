#ifndef RUNTAB_H
#define RUNTAB_H

#include <QMap>
#include <QWidget>

#include "RunWorker.h"

QT_BEGIN_NAMESPACE
class QComboBox;
class QFormLayout;
class QGroupBox;
class QLabel;
class QListWidget;
class QPlainTextEdit;
class QPushButton;
class QSpinBox;
class QTableWidget;
QT_END_NAMESPACE

/*!
 * The Calculate tab: which command to run (the engine's commands: calculate, fit, extrapolate, minos, rate, mcmc where
 * built), its options drawn from the engine's table of options (azure2::front::Options: the uncertainty band, the
 * parameter and capture files, the rate's pairs and temperatures, the sampler's walkers and steps), the threads, the
 * run's lines and progress, and what it returned: the chi^2 of each data set and the files it wrote.
 */
class RunTab : public QWidget {
  Q_OBJECT

 public:
  explicit RunTab(QWidget *parent = nullptr);
  QString command() const;
  bool setCommand(const QString &command);
  /// The command's options as the command line writes them.
  QStringList options() const;
  /// Set an option's editor (a flag: "1" or "0"), for the tests.
  bool setOption(const QString &option, const QString &value);
  int threads() const;
  void setRunning(bool running);
  QString log() const;
  QPushButton *runButton() const { return runButton_; }
  QPushButton *stopButton() const { return stopButton_; }
  QTableWidget *chiTable() const { return chi_; }

 public slots:
  void appendLine(const QString &line, bool transient);
  void showProgress(const QString &stage, qlonglong evaluation, double objective, double best);
  void showOutcome(const RunOutcome &outcome);
  void clearLog();

 signals:
  void runRequested();
  void stopRequested();

 private slots:
  void commandChanged();

 private:
  QComboBox *command_;
  QSpinBox *threads_;
  QGroupBox *optionsBox_;
  QFormLayout *optionsForm_ = nullptr;
  QMap<QString, QWidget *> editors_;  ///< option ("--covariance-band", "--steps=") -> editor
  QMap<QString, QString> kept_;       ///< the values typed, kept across command changes
  QPushButton *runButton_, *stopButton_;
  QLabel *progress_;
  QPlainTextEdit *log_;
  QTableWidget *chi_;
  QListWidget *files_;
  QString outputDirectory_;
};

#endif
