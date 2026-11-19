#ifndef RUNTAB_H
#define RUNTAB_H

#include <QMap>
#include <QVector>
#include <QWidget>

#include "RunWorker.h"

QT_BEGIN_NAMESPACE
class QComboBox;
class QGridLayout;
class QGroupBox;
class QLabel;
class QListWidget;
class QPlainTextEdit;
class QPushButton;
class QSpinBox;
class QTableWidget;
QT_END_NAMESPACE

/*!
 * A small curve of the best chi^2 against the evaluation while a fit runs (log y), from the run's progress records.
 */
class ChiSquaredTrace : public QWidget {
  Q_OBJECT

 public:
  explicit ChiSquaredTrace(QWidget *parent = nullptr);
  void clear();
  void add(double evaluation, double chiSquared);
  int points() const { return (int)x_.size(); }
  QSize sizeHint() const override;

 protected:
  void paintEvent(QPaintEvent *) override;

 private:
  QVector<double> x_, y_;
};

/*!
 * The Calculate tab (docs/upgrade/qt_gui.md sections 3 and 11): which command to run (the engine's commands:
 * calculate, fit, extrapolate, minos, rate, mcmc where built) with its help in one line; every option of the command
 * an editor of its kind made from the engine's table of options (azure2::front::Options): a flag a check box, a number
 * a spin box with the row's range, default and unit, a word a drop-down, a file a path with a chooser, a pair key a
 * drop-down of the project's pairs, each with its help in one line under it, an option that applies only with another
 * enabled only then; the threads; Run and Stop; while a fit runs its evaluation, chi^2 and best so far and a curve of
 * the best chi^2; the run's lines with warnings and errors highlighted; and at the end a summary: the total chi^2 and,
 * per data set, its name, N, chi^2 and chi^2/N in a sortable table, the files the run wrote and "Open Plot".
 */
class RunTab : public QWidget {
  Q_OBJECT

 public:
  explicit RunTab(QWidget *parent = nullptr);
  QString command() const;
  bool setCommand(const QString &command);
  /// The command's options as the command line writes them.
  QStringList options() const;
  /// Set an option's editor (a flag: "1" or "0"; "" leaves a value unset), for the tests.
  bool setOption(const QString &option, const QString &value);
  /// Set an option for every command that takes it (kept across command changes), e.g. --parameter-file=.
  void presetOption(const QString &option, const QString &value);
  /// The project's pairs (key, label) for the rate's drop-downs of entrance and exit pair.
  void setPairs(const QList<QPair<int, QString>> &pairs);
  /// The names of the active data sets in the order of the program's "Segment #" (the summary's names).
  void setDataSetNames(const QStringList &names);
  int threads() const;
  void setRunning(bool running);
  QString log() const;
  QPushButton *runButton() const { return runButton_; }
  QPushButton *stopButton() const { return stopButton_; }
  QPushButton *plotButton() const { return plotButton_; }
  QTableWidget *chiTable() const { return chi_; }
  QLabel *totalLabel() const { return total_; }
  ChiSquaredTrace *trace() const { return trace_; }
  /// The editor of an option ("--steps="), nullptr when the command has none.
  QWidget *editor(const QString &option) const { return editors_.value(option, nullptr); }
  /// Lines of the log shown as a warning or an error (the highlighter's rule), for the tests.
  static int lineSeverity(const QString &line);

 public slots:
  void appendLine(const QString &line, bool transient);
  void showProgress(const QString &stage, qlonglong evaluation, qlonglong calls, double objective, double best);
  void showOutcome(const RunOutcome &outcome);
  void clearLog();

 signals:
  void runRequested();
  void stopRequested();
  void plotRequested();

 private slots:
  void commandChanged();
  void updateDependents();

 private:
  QString valueOf(const QString &option) const;

  QComboBox *command_;
  QLabel *commandHelp_;
  QSpinBox *threads_;
  QGroupBox *optionsBox_;
  QGridLayout *optionsGrid_ = nullptr;
  QMap<QString, QWidget *> editors_;  ///< option ("--covariance-band", "--steps=") -> editor
  QMap<QString, QString> needs_;      ///< option -> the flag it applies with
  QList<QWidget *> rowWidgets_;       ///< the grid's widgets, removed on a command change
  QMap<QString, QString> kept_;       ///< the values typed, kept across command changes
  QPushButton *runButton_, *stopButton_, *plotButton_;
  QLabel *progress_;
  ChiSquaredTrace *trace_;
  QPlainTextEdit *log_;
  QTableWidget *chi_;
  QLabel *total_;
  QListWidget *files_;
  QString outputDirectory_;
  QList<QPair<int, QString>> pairs_;
  QStringList dataSetNames_;
};

#endif
