#ifndef EDITCHECKSDIALOG_H
#define EDITCHECKSDIALOG_H

#include <QDialog>

QT_BEGIN_NAMESPACE

class QPushButton;
class QGroupBox;
class QComboBox;

QT_END_NAMESPACE

/*!
 * Dialog choosing which check files a run writes, and whether to screen or file.
 */
class EditChecksDialog : public QDialog {
  Q_OBJECT

 public:
  EditChecksDialog(QWidget *parent = 0);
  QComboBox *compoundCheckCombo;
  QComboBox *boundaryCheckCombo;
  QComboBox *dataCheckCombo;
  QComboBox *lMatrixCheckCombo;
  QComboBox *legendreCheckCombo;
  QComboBox *coulAmpCheckCombo;
  QComboBox *pathwaysCheckCombo;
  QComboBox *angDistsCheckCombo;

 private:
  QPushButton *okButton;
  QPushButton *cancelButton;
  QGroupBox *stoppingPowerBox;
};

#endif
