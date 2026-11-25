#ifndef DIAGNOSTICSPANEL_H
#define DIAGNOSTICSPANEL_H

#include <QWidget>

#include "ProjectHolder.h"

QT_BEGIN_NAMESPACE
class QLabel;
class QTreeWidget;
QT_END_NAMESPACE

/*!
 * The engine's diagnostics of the project as it stands (AZURE2 --validate through azure2::front::Project::Validate):
 * one row each, severity, code, line, record and field, message.  A double click asks the window to show the record.
 */
class DiagnosticsPanel : public QWidget {
  Q_OBJECT

 public:
  explicit DiagnosticsPanel(ProjectHolder *holder, QWidget *parent = nullptr);
  int errorCount() const { return errors_; }
  int warningCount() const { return warnings_; }
  /// The summary line ("Valid: no diagnostic", "Not valid: 2 error(s) ...").
  QString summary() const;
  /// Double-click on row `row` (the tests).
  void activateRow(int row);
  int rowCount() const;
  /// The codes shown, for the tests.
  QStringList codes() const;
  /// The rows as text, "severity code (line N): message", the first `maxRows` of them (the notes bar's Details).
  QString listing(int maxRows) const;

 public slots:
  /// Validate the project and show the result; returns whether it is valid.
  bool validate();

 signals:
  void showRecord(const QString &section, const QString &record);
  /// After every validation: valid or not, the errors and the other diagnostics.
  void validated(bool ok, int errors, int warnings);

 private:
  ProjectHolder *holder_;
  QTreeWidget *tree_;
  QLabel *summary_;
  int errors_ = 0;
  int warnings_ = 0;
};

#endif
