#include "AddSegTestDialog.h"

#include <QPushButton>
#include <QGroupBox>
#include <QGridLayout>
#include <QHBoxLayout>


AddSegTestDialog::AddSegTestDialog(QWidget *parent) : QDialog(parent) {

  //  this->setMaximumSize(370,420);
  //  this->setMinimumSize(370,420);

  entrancePairIndexSpin = new QSpinBox;
  entrancePairIndexSpin->setMinimum(1);
  entrancePairIndexSpin->setMaximum(100);
  entrancePairIndexSpin->setSingleStep(1);
  exitPairIndexSpin = new QSpinBox;
  exitPairIndexSpin->setMinimum(1);
  exitPairIndexSpin->setMaximum(100);
  exitPairIndexSpin->setSingleStep(1);
  lowEnergyText = new QLineEdit;
  highEnergyText = new QLineEdit;
  energyStepText = new QLineEdit;
  lowAngleText = new QLineEdit;
  lowAngleText->setText("0");
  lowAngleText->setEnabled(false);
  highAngleText = new QLineEdit;
  highAngleText->setText("0");
  highAngleText->setEnabled(false);
  angleStepText = new QLineEdit;
  angleStepText->setText("0");
  angleStepText->setEnabled(false);
  dataTypeCombo = new QComboBox;
  dataTypeCombo->addItem(tr("Angle Integrated"));
  dataTypeCombo->addItem(tr("Differential"));
  dataTypeCombo->addItem(tr("Phase Shift"));
  dataTypeCombo->addItem(tr("Angular Distribution Coefficients"));
  dataTypeCombo->addItem(tr("Angle Integrated Total Capture"));
  dataTypeCombo->addItem(tr("C.M. Differential"));
  connect(dataTypeCombo,SIGNAL(currentIndexChanged(int)),this,SLOT(dataTypeChanged(int)));
  QRegExp spinRX("^\\d{0,2}(\\.[05]{0,1})?$");
  QValidator *spinValidator = new QRegExpValidator(spinRX, this);
  phaseJValueText = new QLineEdit;
  phaseJValueText->setValidator(spinValidator);
  phaseJValueText->setVisible(false);
  phaseJValueText->setMaximumWidth(50);
  QRegExp intRX("^[0-6]$");
  QValidator *intValidator = new QRegExpValidator(intRX, this);
  phaseLValueText = new QLineEdit;
  phaseLValueText->setValidator(intValidator);
  phaseLValueText->setVisible(false);
  phaseLValueText->setMaximumWidth(50);

  advancedModeCheck = new QCheckBox(tr("Advanced Mode (Sum/Ratio)"));
  connect(advancedModeCheck,SIGNAL(stateChanged(int)),this,SLOT(advancedModeChanged(int)));

  advancedModeBox = new QGroupBox(tr("Advanced Segment Definition"));
  advancedModeBox->setVisible(false);
  
  operationCombo = new QComboBox;
  operationCombo->addItem(tr("Sum"));
  operationCombo->addItem(tr("Ratio"));

  componentsList = new QListWidget;
  componentsList->setMaximumHeight(100);

  componentEntranceSpin = new QSpinBox;
  componentEntranceSpin->setMinimum(1);
  componentEntranceSpin->setMaximum(100);
  componentEntranceSpin->setSingleStep(1);
  componentEntranceSpin->setValue(1);

  componentExitSpin = new QSpinBox;
  componentExitSpin->setMinimum(1);
  componentExitSpin->setMaximum(100);
  componentExitSpin->setSingleStep(1);
  componentExitSpin->setValue(1);

  addComponentButton = new QPushButton(tr("Add Component"));
  connect(addComponentButton,SIGNAL(clicked()),this,SLOT(addComponent()));
  
  removeComponentButton = new QPushButton(tr("Remove Component"));
  connect(removeComponentButton,SIGNAL(clicked()),this,SLOT(removeComponent()));

  cancelButton = new QPushButton(tr("Cancel"));
  okButton = new QPushButton(tr("Accept"));
  okButton->setDefault(true);

  QGroupBox *valueBox = new QGroupBox;
  QGridLayout *valueLayout = new QGridLayout;
  QGridLayout *pairLayout = new QGridLayout;
  pairLayout->addWidget(new QLabel(tr("Entrance Pair Key:")),0,0,Qt::AlignRight);
  pairLayout->addWidget(entrancePairIndexSpin,0,1);
  pairLayout->addWidget(new QLabel(tr("Exit Pair Key:")),0,2,Qt::AlignRight);
  pairLayout->addWidget(exitPairIndexSpin,0,3);
  totalCaptureLabel = new QLabel(tr("Total Capture"));
  totalCaptureLabel->setVisible(false);
  pairLayout->addWidget(totalCaptureLabel,0,4);
  valueLayout->addLayout(pairLayout,0,0,1,2);
  QGroupBox *energyBox = new QGroupBox(tr("Lab Energy [MeV]"));
  QGridLayout *energyLayout = new QGridLayout;
  energyLayout->addWidget(new QLabel(tr("Low Energy:")),0,0,Qt::AlignRight);
  energyLayout->addWidget(lowEnergyText,0,1);
  energyLayout->addWidget(new QLabel(tr("High Energy:")),1,0,Qt::AlignRight);
  energyLayout->addWidget(highEnergyText,1,1);
  energyLayout->addWidget(new QLabel(tr("Energy Step:")),2,0,Qt::AlignRight);
  energyLayout->addWidget(energyStepText,2,1);
  energyBox->setLayout(energyLayout);
  valueLayout->addWidget(energyBox,1,0);
  QGroupBox *angleBox = new QGroupBox(tr("Lab Angle [degrees]"));
  QGridLayout *angleLayout = new QGridLayout;
  angleLayout->addWidget(new QLabel(tr("Low Angle:")),0,0,Qt::AlignRight);
  angleLayout->addWidget(lowAngleText,0,1);
  angleLayout->addWidget(new QLabel(tr("High Angle:")),1,0,Qt::AlignRight);
  angleLayout->addWidget(highAngleText,1,1);
  angleLayout->addWidget(new QLabel(tr("Angle Step:")),2,0,Qt::AlignRight);
  angleLayout->addWidget(angleStepText,2,1);
  angleBox->setLayout(angleLayout);
  valueLayout->addWidget(angleBox,1,1);

  QGridLayout* lowerLayout = new QGridLayout;
  lowerLayout->addWidget(new QLabel(tr("Data Type:")),0,0,Qt::AlignRight);
  lowerLayout->addWidget(dataTypeCombo,0,1);
  lowerLayout->addItem(new QSpacerItem(1,25),0,2);
  lowerLayout->setColumnStretch(2,1);

  QHBoxLayout* phaseLayout = new QHBoxLayout;
  phaseJValueLabel = new QLabel(tr("J:"));
  phaseJValueLabel->setVisible(false);
  phaseLayout->addWidget(phaseJValueLabel);
  phaseLayout->addWidget(phaseJValueText);
  phaseLValueLabel = new QLabel(tr("l:"));
  phaseLValueLabel->setVisible(false);
  phaseLayout->addWidget(phaseLValueLabel);
  phaseLayout->addWidget(phaseLValueText);
  angDistLabel = new QLabel(tr("Maximum Order"));
  angDistSpin = new QSpinBox;
  angDistSpin->setMinimum(0);
  angDistSpin->setMaximum(10);
  angDistSpin->setSingleStep(1);
  angDistLabel->setVisible(false);
  angDistSpin->setVisible(false);
  phaseLayout->addWidget(angDistLabel);
  phaseLayout->addWidget(angDistSpin);
  lowerLayout->addLayout(phaseLayout,0,3);

  lowerLayout->addWidget(advancedModeCheck,1,0,1,4);
  
  QGridLayout *advancedLayout = new QGridLayout;
  advancedLayout->addWidget(new QLabel(tr("Operation:")),0,0);
  advancedLayout->addWidget(operationCombo,0,1);
  advancedLayout->addWidget(new QLabel(tr("Components:")),1,0);
  advancedLayout->addWidget(componentsList,1,1,3,1);

  // Add component entrance/exit spinboxes above buttons
  advancedLayout->addWidget(new QLabel(tr("Entrance Pair Key:")),1,2);
  advancedLayout->addWidget(componentEntranceSpin,1,3);
  advancedLayout->addWidget(new QLabel(tr("Exit Pair Key:")),2,2);
  advancedLayout->addWidget(componentExitSpin,2,3);

  QGridLayout *buttonLayout = new QGridLayout;
  buttonLayout->addWidget(addComponentButton,0,0);
  buttonLayout->addWidget(removeComponentButton,0,1);
  advancedLayout->addLayout(buttonLayout,3,2,1,2);
  advancedModeBox->setLayout(advancedLayout);

  lowerLayout->addWidget(advancedModeBox,2,0,1,4);

  valueLayout->addLayout(lowerLayout,2,0,1,2);
  valueBox->setLayout(valueLayout);

  QHBoxLayout *buttonBox = new QHBoxLayout;
  buttonBox->addWidget(cancelButton);
  buttonBox->addWidget(okButton);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(valueBox);
  mainLayout->addLayout(buttonBox);

  setLayout(mainLayout);

  connect(okButton, SIGNAL(clicked()),this,SLOT(accept()));
  connect(cancelButton,SIGNAL(clicked()),this,SLOT(reject()));

  setWindowTitle(tr("Add a Segment Without Data"));
}

void AddSegTestDialog::dataTypeChanged(int index) {
  if(index==2) {
    angDistLabel->setVisible(false);
    angDistSpin->setVisible(false);
    phaseJValueLabel->setVisible(true);
    phaseLValueLabel->setVisible(true);
    phaseJValueText->setVisible(true);
    phaseLValueText->setVisible(true);
  } else if(index==3) {
    phaseJValueLabel->setVisible(false);
    phaseLValueLabel->setVisible(false);
    phaseJValueText->setVisible(false);
    phaseLValueText->setVisible(false);
    angDistLabel->setVisible(true);
    angDistSpin->setVisible(true);
  } else {
    phaseJValueLabel->setVisible(false);
    phaseLValueLabel->setVisible(false);
    phaseJValueText->setVisible(false);
    phaseLValueText->setVisible(false);
    angDistLabel->setVisible(false);
    angDistSpin->setVisible(false);
  }
  if(index==1) {
    lowAngleText->setEnabled(true);
    highAngleText->setEnabled(true);
    angleStepText->setEnabled(true);
  } else {
    lowAngleText->setEnabled(false);
    highAngleText->setEnabled(false);
    angleStepText->setEnabled(false);
  }
  if(index==4) {
    exitPairIndexSpin->setVisible(false);
    totalCaptureLabel->setVisible(true);
  } else {
    totalCaptureLabel->setVisible(false);
    exitPairIndexSpin->setVisible(true);
  }
  if(index==5) {
    lowAngleText->setEnabled(true);
    highAngleText->setEnabled(true);
    angleStepText->setEnabled(true);
  }
}

void AddSegTestDialog::advancedModeChanged(int state) {
  if(state==Qt::Checked) {
    advancedModeBox->setVisible(true);
    // Keep the spinboxes enabled so users can select pairs to add as components
    entrancePairIndexSpin->setEnabled(true);
    exitPairIndexSpin->setEnabled(true);
  } else {
    advancedModeBox->setVisible(false);
    entrancePairIndexSpin->setEnabled(true);
    exitPairIndexSpin->setEnabled(true);
  }
  adjustSize();
}

void AddSegTestDialog::addComponent() {
  int entrance = componentEntranceSpin->value();
  int exit = componentExitSpin->value();
  QString component = QString("Entrance: %1, Exit: %2").arg(entrance).arg(exit);
  componentsList->addItem(component);
}

void AddSegTestDialog::removeComponent() {
  int currentRow = componentsList->currentRow();
  if(currentRow >= 0) {
    delete componentsList->takeItem(currentRow);
  }
}
