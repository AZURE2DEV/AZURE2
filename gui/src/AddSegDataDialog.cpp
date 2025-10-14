#include "AddSegDataDialog.h"

#include <QPushButton>
#include <QGroupBox>
#include <QGridLayout>
#include <QFileDialog>

AddSegDataDialog::AddSegDataDialog(QWidget *parent) : QDialog(parent) {

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
  lowAngleText = new QLineEdit;
  lowAngleText->setText("0");
  lowAngleText->setEnabled(false);
  highAngleText = new QLineEdit;
  highAngleText->setText("180");
  highAngleText->setEnabled(false);
  dataTypeCombo = new QComboBox;
  dataTypeCombo->addItem(tr("Angle Integrated"));
  dataTypeCombo->addItem(tr("Differential"));
  dataTypeCombo->addItem(tr("Phase Shift"));
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
  dataFileText = new QLineEdit;
  QPushButton *chooseFileButton = new QPushButton(tr("Choose..."));
  connect(chooseFileButton,SIGNAL(clicked()),this,SLOT(setChooseFile()));
  dataNormText = new QLineEdit;
  dataNormText->setText("1.0");
  dataNormErrorLabel = new QLabel(tr("Norm. Error [%]:"));
  //  dataNormErrorLabel->setVisible(false);
  dataNormErrorText = new QLineEdit(this);
  //  dataNormErrorText->setVisible(false);
  dataNormErrorText->setText("0.0");
  dataNormErrorText->setMaximumWidth(50);
  varyNormCheck = new QCheckBox(tr("Vary Norm?"));
  //connect(varyNormCheck,SIGNAL(stateChanged(int)),this,SLOT(varyNormChanged(int)));
  
  energyShiftText = new QLineEdit;
  energyShiftText->setText("0.0");
  energyShiftLabel = new QLabel(tr("Energy Shift [MeV]:"));
  energyShiftErrorText = new QLineEdit;
  energyShiftErrorText->setText("0.0");
  energyShiftErrorLabel = new QLabel(tr("Energy Shift Error [MeV]:"));
  varyEnergyShiftCheck = new QCheckBox(tr("Vary Energy Shift?"));

  uposCheck = new QCheckBox(tr("Unobserved Primary, Observed Secondary?"));
  secondaryLText = new QLineEdit;
  secondaryLText->setText("0");
  secondaryLText->setMaximumWidth(50);
  secondaryLText->setVisible(false);
  secondaryLLabel = new QLabel(tr("Secondary L:"));
  secondaryLLabel->setVisible(false);
  finalJText = new QLineEdit;
  finalJText->setText("0.0");
  finalJText->setMaximumWidth(50);
  finalJText->setVisible(false);
  finalJLabel = new QLabel(tr("Final J:"));
  finalJLabel->setVisible(false);
  deltaText = new QLineEdit;
  deltaText->setText("0.0");
  deltaText->setMaximumWidth(50);
  deltaText->setVisible(false);
  deltaLabel = new QLabel(tr("Mixing Ratio:"));
  deltaLabel->setVisible(false);
  connect(uposCheck,SIGNAL(stateChanged(int)),this,SLOT(uposChanged(int)));

  advancedModeCheck = new QCheckBox(tr("Advanced Settings (Sum or Ratio of Components)"));
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
  QGroupBox* energyBox = new QGroupBox(tr("Lab Energy [MeV]"));
  QGridLayout *energyLayout = new QGridLayout;
  energyLayout->addWidget(new QLabel(tr("Low Energy:")),0,0,Qt::AlignRight);
  energyLayout->addWidget(lowEnergyText,0,1);
  energyLayout->addWidget(new QLabel(tr("High Energy:")),1,0,Qt::AlignRight);
  energyLayout->addWidget(highEnergyText,1,1);
  energyBox->setLayout(energyLayout);
  valueLayout->addWidget(energyBox,1,0);
  QGroupBox* angleBox = new QGroupBox(tr("Lab Angle [degrees]"));
  QGridLayout *angleLayout = new QGridLayout;
  angleLayout->addWidget(new QLabel(tr("Low Angle:")),0,0,Qt::AlignRight);
  angleLayout->addWidget(lowAngleText,0,1);
  angleLayout->addWidget(new QLabel(tr("High Angle:")),1,0,Qt::AlignRight);
  angleLayout->addWidget(highAngleText,1,1);
  angleBox->setLayout(angleLayout);
  valueLayout->addWidget(angleBox,1,1);

  QGridLayout* lowerLayout = new QGridLayout;
  lowerLayout->addWidget(new QLabel(tr("Data Type:")),0,0,Qt::AlignRight);
  lowerLayout->addWidget(dataTypeCombo,0,1);

  QGridLayout* phaseLayout = new QGridLayout;
  phaseLayout->addItem(new QSpacerItem(1,25),0,0);
  phaseLayout->setColumnStretch(0,1);
  phaseJValueLabel = new QLabel(tr("J:"));
  phaseJValueLabel->setVisible(false);
  phaseLayout->addWidget(phaseJValueLabel,0,1);
  phaseLayout->addWidget(phaseJValueText,0,2);
  phaseLValueLabel = new QLabel(tr("l:"));
  phaseLValueLabel->setVisible(false);
  phaseLayout->addWidget(phaseLValueLabel,0,3);
  phaseLayout->addWidget(phaseLValueText,0,4);
  lowerLayout->addLayout(phaseLayout,0,2);

  lowerLayout->addWidget(new QLabel(tr("Data Norm.:")),1,0,Qt::AlignRight);
  QHBoxLayout *normLayout = new QHBoxLayout;
  normLayout->addWidget(dataNormText);
  normLayout->addWidget(varyNormCheck);
  lowerLayout->addLayout(normLayout,1,1);
  QGridLayout *normErrorLayout = new QGridLayout;
  normErrorLayout->addItem(new QSpacerItem(1,25),0,0);
  normErrorLayout->setColumnStretch(0,1);
  normErrorLayout->addWidget(dataNormErrorLabel,0,1,Qt::AlignRight);
  normErrorLayout->addWidget(dataNormErrorText,0,2);
  lowerLayout->addLayout(normErrorLayout,1,2);  
  
  lowerLayout->addWidget(new QLabel(tr("Data File:")),2,0,Qt::AlignRight);
  QGridLayout *fileLayout = new QGridLayout;
  fileLayout->addWidget(dataFileText,0,0);
  fileLayout->addWidget(chooseFileButton,0,1);
  fileLayout->setColumnStretch(0,1);
  lowerLayout->addLayout(fileLayout,2,1,1,2);
  
  lowerLayout->addWidget(energyShiftLabel,3,0,Qt::AlignRight);
  QHBoxLayout *energyShiftLayout = new QHBoxLayout;
  energyShiftLayout->addWidget(energyShiftText);
  energyShiftLayout->addWidget(varyEnergyShiftCheck);
  lowerLayout->addLayout(energyShiftLayout,3,1);
  QGridLayout *energyShiftErrorLayout = new QGridLayout;
  energyShiftErrorLayout->addItem(new QSpacerItem(1,25),0,0);
  energyShiftErrorLayout->setColumnStretch(0,1);
  energyShiftErrorLayout->addWidget(energyShiftErrorLabel,0,1,Qt::AlignRight);
  energyShiftErrorLayout->addWidget(energyShiftErrorText,0,2);
  lowerLayout->addLayout(energyShiftErrorLayout,3,2);

  lowerLayout->addWidget(uposCheck,4,0,1,3);
  QHBoxLayout *uposParamsLayout = new QHBoxLayout;
  uposParamsLayout->addWidget(secondaryLLabel);
  uposParamsLayout->addWidget(secondaryLText);
  uposParamsLayout->addWidget(finalJLabel);
  uposParamsLayout->addWidget(finalJText);
  uposParamsLayout->addWidget(deltaLabel);
  uposParamsLayout->addWidget(deltaText);
  lowerLayout->addLayout(uposParamsLayout,5,0,1,3);

  lowerLayout->addWidget(advancedModeCheck,6,0,1,3);
  
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

  lowerLayout->addWidget(advancedModeBox,7,0,1,3);
  
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

  setWindowTitle(tr("Add a Segment From Data"));
}

void AddSegDataDialog::setChooseFile() {
  QString filename = QFileDialog::getOpenFileName(this);
  if(!filename.isEmpty()) {
    dataFileText->setText(QDir::fromNativeSeparators(filename));
  }
}

void AddSegDataDialog::dataTypeChanged(int index) {
  if(index==2) {
    phaseJValueLabel->setVisible(true);
    phaseLValueLabel->setVisible(true);
    phaseJValueText->setVisible(true);
    phaseLValueText->setVisible(true);
  } else {
    phaseJValueLabel->setVisible(false);
    phaseLValueLabel->setVisible(false);
    phaseJValueText->setVisible(false);
    phaseLValueText->setVisible(false);
  }
  if(index==1) {
    lowAngleText->setEnabled(true);
    highAngleText->setEnabled(true);
  } else {
    lowAngleText->setEnabled(false);
    highAngleText->setEnabled(false);
  }
  if(index==3) {
    exitPairIndexSpin->setVisible(false);
    totalCaptureLabel->setVisible(true);
  } else {
    totalCaptureLabel->setVisible(false);
    exitPairIndexSpin->setVisible(true);
  }
  if(index==4) {
    lowAngleText->setEnabled(true);
    highAngleText->setEnabled(true);
  }
}

void AddSegDataDialog::varyNormChanged(int state) {
  if(state==Qt::Checked) {
    dataNormErrorLabel->setVisible(true);
    dataNormErrorText->setVisible(true);
  } else {
    dataNormErrorLabel->setVisible(false);
    dataNormErrorText->setVisible(false);
  }
}

void AddSegDataDialog::advancedModeChanged(int state) {
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

void AddSegDataDialog::addComponent() {
  int entrance = componentEntranceSpin->value();
  int exit = componentExitSpin->value();
  QString component = QString("Entrance: %1, Exit: %2").arg(entrance).arg(exit);
  componentsList->addItem(component);
}

void AddSegDataDialog::removeComponent() {
  int currentRow = componentsList->currentRow();
  if(currentRow >= 0) {
    delete componentsList->takeItem(currentRow);
  }
}

void AddSegDataDialog::uposChanged(int state) {
  if(state==Qt::Checked) {
    secondaryLLabel->setVisible(true);
    secondaryLText->setVisible(true);
    finalJLabel->setVisible(true);
    finalJText->setVisible(true);
    deltaLabel->setVisible(true);
    deltaText->setVisible(true);
  } else {
    secondaryLLabel->setVisible(false);
    secondaryLText->setVisible(false);
    finalJLabel->setVisible(false);
    finalJText->setVisible(false);
    deltaLabel->setVisible(false);
    deltaText->setVisible(false);
  }
  adjustSize();
}
