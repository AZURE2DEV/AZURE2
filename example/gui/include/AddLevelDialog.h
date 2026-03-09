#ifndef ADDLEVELDIALOG_H
#define ADDLEVELDIALOG_H

#include <QDialog>
#include <QtGui>

QT_BEGIN_NAMESPACE

class QLineEdit;
class QLabel;

QT_END_NAMESPACE

class AddLevelDialog : public QDialog {
  Q_OBJECT

 public:
  AddLevelDialog(QWidget *parent=0);
  QLineEdit *jValueText;
  QComboBox *piValueCombo;
  QLineEdit *energyText;
 
 private:
  QLabel *jValueLabel;
  QLabel *energyLabel;

  QPushButton *okButton;
  QPushButton *cancelButton;
};

#endif
