#ifndef ADDSEGTESTDIALOG_H
#define ADDSEGTESTDIALOG_H

#include <QDialog>
#include <QLineEdit>
#include <QSpinBox>
#include <QComboBox>
#include <QLabel>
#include <QCheckBox>
#include <QListWidget>
#include <QGroupBox>
#include <QPushButton>

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
  QCheckBox *thmCheck;
  QLineEdit *phaseJValueText;
  QLineEdit *phaseLValueText;
  QLabel *phaseJValueLabel;
  QLabel *phaseLValueLabel;
  QLabel *angDistLabel;
  QLabel *totalCaptureLabel;
  QSpinBox *angDistSpin;

  QCheckBox *advancedModeCheck;
  QGroupBox *advancedModeBox;
  QComboBox *operationCombo;
  QListWidget *componentsList;
  QSpinBox *componentEntranceSpin;
  QSpinBox *componentExitSpin;
  QCheckBox *useFixedAngleCheck;
  QLineEdit *fixedAngleText;
  QLabel *fixedAngleLabel;
  QLabel *componentScalingLabel;
  QLineEdit *componentScalingText;
  QPushButton *addComponentButton;
  QPushButton *removeComponentButton;

 public slots:
  void dataTypeChanged(int);
  void advancedModeChanged(int);
  void operationTypeChanged(int);
  void useFixedAngleChanged(int);
  void addComponent();
  void removeComponent();
 
 private:
  QPushButton *okButton;
  QPushButton *cancelButton;
};

#endif
