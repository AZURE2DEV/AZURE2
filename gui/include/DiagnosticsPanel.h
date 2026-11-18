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
  int rowCount() const;
  /// The codes shown, for the tests.
  QStringList codes() const;

 public slots:
  /// Validate the project and show the result; returns whether it is valid.
  bool validate();

 signals:
  void showRecord(const QString &section, const QString &record);

 private:
  ProjectHolder *holder_;
  QTreeWidget *tree_;
  QLabel *summary_;
  int errors_ = 0;
};

#endif
