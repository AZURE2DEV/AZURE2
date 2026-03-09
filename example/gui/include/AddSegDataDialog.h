#ifndef ADDSEGDATADIALOG_H
#define ADDSEGDATADIALOG_H

#include <QDialog>
#include <QLineEdit>
#include <QSpinBox>
#include <QComboBox>
#include <QLabel>
#include <QCheckBox>

QT_BEGIN_NAMESPACE

class QLabel;

QT_END_NAMESPACE

class AddSegDataDialog : public QDialog {
  Q_OBJECT

 public:
  AddSegDataDialog(QWidget *parent=0);
  QSpinBox *entrancePairIndexSpin;
  QSpinBox *exitPairIndexSpin;
  QLineEdit *lowEnergyText;
  QLineEdit *highEnergyText;
  QLineEdit *lowAngleText;
  QLineEdit *highAngleText;
  QComboBox *dataTypeCombo;
  QLineEdit *dataFileText;
  QLineEdit *dataNormText;
  QLineEdit *dataNormErrorText;
  QLabel *dataNormErrorLabel;
  QCheckBox *varyNormCheck;
  QLineEdit *phaseJValueText;
  QLineEdit *phaseLValueText;
  QLabel* phaseLValueLabel;
  QLabel* phaseJValueLabel;
  QLabel* totalCaptureLabel;
  QCheckBox *uposCheck;
  QLineEdit *secondaryLText;
  QLabel* secondaryLLabel;
  QLineEdit *finalJText;
  QLabel* finalJLabel;
  QLineEdit* deltaText;
  QLabel* deltaLabel;

 public slots:
  void setChooseFile();
  void dataTypeChanged(int);
  void varyNormChanged(int);
  
 private:
  QPushButton *okButton;
  QPushButton *cancelButton;
};

#endif
