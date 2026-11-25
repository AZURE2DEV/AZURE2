#ifndef CONVERTREPAIRDIALOG_H
#define CONVERTREPAIRDIALOG_H

#include <QDialog>
#include <QList>

#include "azure2/front/Engine.h"

QT_BEGIN_NAMESPACE
class QLineEdit;
class QPushButton;
QT_END_NAMESPACE

/*!
 * A line of an old project the converter could not read (docs/upgrade/qt_gui.md section 19.5): the file and the line
 * number, the line itself in an editable field, the refusal, the fields the old grammar expects, and one button per
 * reading the converter proposes (a button puts its line into the field; nothing is taken without the user).  "Convert
 * again" converts a copy of the old project with the line as the field holds it; the original file is kept.
 */
class ConvertRepairDialog : public QDialog {
  Q_OBJECT

 public:
  ConvertRepairDialog(const QString &file, const azure2::front::ConvertLineProblem &problem, QWidget *parent = nullptr);
  /// The line as the field holds it (the user's edit, or a reading taken).
  QString lineText() const;
  /// The buttons of the readings, in the converter's order (for the tests).
  QList<QPushButton *> readingButtons() const { return readings_; }
  QPushButton *retryButton() const { return retry_; }

 private:
  QLineEdit *edit_ = nullptr;
  QList<QPushButton *> readings_;
  QPushButton *retry_ = nullptr;
};

#endif
