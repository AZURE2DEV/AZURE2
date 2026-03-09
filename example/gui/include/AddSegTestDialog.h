#ifndef ADDSEGTESTDIALOG_H
#define ADDSEGTESTDIALOG_H

#include <QDialog>
#include <QLineEdit>
#include <QSpinBox>
#include <QComboBox>
#include <QLabel>
#include <QCheckBox>

QT_BEGIN_NAMESPACE

class QLabel;

QT_END_NAMESPACE

class AddSegTestDialog : public QDialog {
  Q_OBJECT

 public:
  AddSegTestDialog(QWidget *parent=0);
  QSpinBox *entrancePairIndexSpin;
  QSpinBox *exitPairIndexSpin;
  QLineEdit *lowEnergyText;
  QLineEdit *highEnergyText;
  QLineEdit *energyStepText;
  QLineEdit *lowAngleText;
  QLineEdit *highAngleText;
  QLineEdit *angleStepText;
  QComboBox *dataTypeCombo;
  QLineEdit *phaseJValueText;
  QLineEdit *phaseLValueText;
  QLabel *phaseJValueLabel;
  QLabel *phaseLValueLabel;
  QLabel *angDistLabel;
  QLabel *totalCaptureLabel;
  QSpinBox *angDistSpin;
  QCheckBox *uposCheck;
  QLineEdit *secondaryLText;
  QLineEdit *finalJText;
  QLineEdit *deltaText;
  QLabel *secondaryLLabel;
  QLabel *finalJLabel;
  QLabel *deltaLabel;
 
 public slots:
  void dataTypeChanged(int);
 
 private:
  QPushButton *okButton;
  QPushButton *cancelButton;
};

#endif
