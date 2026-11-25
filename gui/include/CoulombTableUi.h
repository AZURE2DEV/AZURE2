#ifndef COULOMBTABLEUI_H
#define COULOMBTABLEUI_H

#include <QDialog>
#include <QList>
#include <QPointer>
#include <QWidget>
#include <memory>

#include "azure2/front/Engine.h"

QT_BEGIN_NAMESPACE
class QLabel;
class QProgressBar;
class QPushButton;
class QThread;
class QTimer;
QT_END_NAMESPACE

/*!
 * The Coulomb table of this computer in the window (docs/upgrade/qt_gui.md section 19.12).  A run without the table
 * computes the exact Coulomb functions: the same results, about 4 times slower.  The table is not distributed; it is
 * built once per computer (about 5 min, 88 MB) with the command line's routine (front::BuildCoulombTable).
 *
 * CoulombTableNotice: one quiet line with a Build... button, on the start page and on the Overview, shown only while
 * the runs have no usable table (none at the path, or a file there that cannot be read: then its reason).
 */
class CoulombTableNotice : public QWidget {
  Q_OBJECT

 public:
  explicit CoulombTableNotice(QWidget *parent = nullptr);
  /// Show the line for `state` (hidden when the table is usable); `buildAllowed` false with `why` while a calculation
  /// or a build runs.
  void setState(const azure2::front::CoulombTableState &state, bool buildAllowed, const QString &why);
  QLabel *label() const { return text_; }
  QPushButton *buildButton() const { return build_; }

 signals:
  void buildRequested();

 private:
  QLabel *text_;
  QPushButton *build_;
};

/*!
 * The table's dialog (Run > Build the Coulomb Table... and Check the Coulomb Table; not modal): where runs look for the
 * table and what they find there, a Build button, a progress bar fed by the build's row counter and Cancel.  The build
 * runs in a worker thread on the command line's thread count, at most 4 (front::BuildCoulombTable); it writes
 * <path>.new and renames it at the end, so a cancelled build leaves no file.  It does not start while a calculation
 * of the window runs; a calculation started during a build runs without the table.
 */
class CoulombTableDialog : public QDialog {
  Q_OBJECT

 public:
  explicit CoulombTableDialog(QWidget *parent = nullptr);
  ~CoulombTableDialog() override;
  /// Read the state again (front::CoulombTableStatus) and show it.
  void refresh();
  const azure2::front::CoulombTableState &state() const { return state_; }
  /// A calculation of the window runs (no build may start) or not.
  void setCalculationRunning(bool running);
  bool calculationRunning() const { return calculationRunning_; }
  /// What a build builds (default: the command line's table at the state's path); the tests give a small one.
  void setSpec(const azure2::front::CoulombTableBuild &spec) { spec_ = spec; }
  const azure2::front::CoulombTableBuild &spec() const { return spec_; }
  /// Start the build; false with `why` when it cannot start (a calculation or a build runs).
  bool startBuild(QString *why = nullptr);
  /// Ask the running build to stop (it ends at the next row; nothing is left at the path).
  void cancelBuild();
  bool building() const { return worker_ != nullptr; }
  /// The threads the build uses (the command line's default, at most 4).
  int buildThreads() const;
  QProgressBar *progressBar() const { return bar_; }
  QPushButton *buildButton() const { return build_; }
  QPushButton *cancelButton() const { return cancel_; }
  QLabel *stateLabel() const { return text_; }

 signals:
  /// The state was read again (after refresh() or at the end of a build).
  void stateChanged(const azure2::front::CoulombTableState &state);
  /// A build began / ended (ok, or the reason: "cancelled", an I/O error).
  void buildStarted();
  void buildFinished(bool ok, const QString &error);

 private:
  struct Job;
  void poll();
  void finish();
  void updateButtons();

  azure2::front::CoulombTableState state_;
  azure2::front::CoulombTableBuild spec_;
  bool calculationRunning_ = false;
  std::shared_ptr<Job> job_;
  QThread *worker_ = nullptr;
  QTimer *timer_;
  QLabel *text_, *result_;
  QProgressBar *bar_;
  QPushButton *build_, *cancel_, *close_;
};

/*!
 * The window's one place for the table: its dialog (made when first asked), the notices it keeps in step, and
 * whether a calculation runs.  AZURESetup owns one.
 */
class CoulombTableControl : public QObject {
  Q_OBJECT

 public:
  explicit CoulombTableControl(QWidget *window);
  /// A new notice kept in step with the state (its parent is the caller's).
  CoulombTableNotice *makeNotice();
  /// Read the state again and update the notices.
  void refresh();
  void setCalculationRunning(bool running);
  /// Show the dialog (made on first use).
  CoulombTableDialog *showDialog();
  /// The dialog (null until first shown).
  CoulombTableDialog *dialog() const { return dialog_; }
  const azure2::front::CoulombTableState &state() const { return state_; }

 private:
  void updateNotices();
  QWidget *window_;
  CoulombTableDialog *dialog_ = nullptr;
  QList<QPointer<CoulombTableNotice>> notices_;
  azure2::front::CoulombTableState state_;
  bool running_ = false;
};

#endif
