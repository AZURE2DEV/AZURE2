#include "AddSegDataDialog.h"

#include <QPushButton>
#include <QGroupBox>
#include <QGridLayout>
#include "DialogScrollArea.h"
#include "PairsModel.h"
#include <QFileDialog>

AddSegDataDialog::AddSegDataDialog(QWidget *parent) :
  QDialog(parent) {
  //  this->setMaximumSize(370,420);
  //  this->setMinimumSize(370,420);

  entrancePairCombo = new PairKeyCombo(PairKeyCombo::Entrance);
  exitPairCombo = new PairKeyCombo(PairKeyCombo::Exit);
  // The exit list reads as reactions with the chosen entrance pair.
  connect(entrancePairCombo, SIGNAL(keyChanged(int)), exitPairCombo, SLOT(setEntranceKey(int)));
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
  dataTypeCombo->addItem(tr("Angle Integrated (E1 only)"));
  dataTypeCombo->addItem(tr("Angle Integrated (E2 only)"));
  dataTypeCombo->addItem(tr("Analyzing Power"));
  // Codes 0-6 happen to equal their position; the analyzing power is 7.
  for (int i = 0; i < dataTypeCombo->count(); i++) dataTypeCombo->setItemData(i, i);
  dataTypeCombo->setItemData(dataTypeCombo->count() - 1, 7);
  connect(dataTypeCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(dataTypeChanged(int)));
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
  connect(chooseFileButton, SIGNAL(clicked()), this, SLOT(setChooseFile()));
  dataNormText = new QLineEdit;
  dataNormText->setText("1.0");
  dataNormErrorLabel = new QLabel(tr("Norm. Error [%]:"));
  //  dataNormErrorLabel->setVisible(false);
  dataNormErrorText = new QLineEdit(this);
  //  dataNormErrorText->setVisible(false);
  dataNormErrorText->setText("0.0");
  dataNormErrorText->setMaximumWidth(50);
  varyNormCheck = new QCheckBox(tr("Vary Norm?"));
  // connect(varyNormCheck,SIGNAL(stateChanged(int)),this,SLOT(varyNormChanged(int)));

  energyShiftText = new QLineEdit;
  energyShiftText->setText("0.0");
  energyShiftLabel = new QLabel(tr("Energy Shift [MeV]:"));
  energyShiftErrorText = new QLineEdit;
  energyShiftErrorText->setText("0.0");
  energyShiftErrorLabel = new QLabel(tr("Energy Shift Error [MeV]:"));
  varyEnergyShiftCheck = new QCheckBox(tr("Vary Energy Shift?"));

  // UPOS (Unobserved Primary, Observed Secondary) widgets
  uposCheck = new QCheckBox(tr("Unobserved Primary, Observed Secondary (UPOS)"));
  connect(uposCheck, SIGNAL(stateChanged(int)), this, SLOT(uposChanged(int)));
  secondaryLLabel = new QLabel(tr("Secondary Decay L:"));
  secondaryLLabel->setVisible(false);
  secondaryLSpin = new QSpinBox;
  secondaryLSpin->setMinimum(0);
  secondaryLSpin->setMaximum(10);
  secondaryLSpin->setValue(1);
  secondaryLSpin->setVisible(false);
  finalJLabel = new QLabel(tr("Final State J (Ic):"));
  finalJLabel->setVisible(false);
  finalJText = new QLineEdit;
  finalJText->setText("0.0");
  finalJText->setVisible(false);
  finalJText->setMaximumWidth(60);
  deltaLabel = new QLabel(tr("Mixing Ratio (δ):"));
  deltaLabel->setVisible(false);
  deltaText = new QLineEdit;
  deltaText->setText("0.0");
  deltaText->setVisible(false);
  deltaText->setMaximumWidth(60);

  advancedModeCheck = new QCheckBox(tr("Advanced Settings (Sum or Ratio of Components)"));
  connect(advancedModeCheck, SIGNAL(stateChanged(int)), this, SLOT(advancedModeChanged(int)));

  naturalLabel = new QLabel(tr("Natural target (sample):"));
  naturalCombo = new QComboBox;
  naturalCombo->setObjectName("naturalCombo");
  naturalCombo->setToolTip(tr("Weight this segment's term and each of its components by the isotopic abundance\n"
                              "of the compound nucleus its entrance pair forms, in the chosen sample (the Particle\n"
                              "Pairs tab sets the samples and abundances). Not available for ratios, phase shifts\n"
                              "or analyzing powers."));
  naturalLabel->setVisible(false);
  naturalCombo->setVisible(false);
  setSamples(QStringList() << QString());
  connect(advancedModeCheck, SIGNAL(stateChanged(int)), this, SLOT(updateNaturalCheck()));

  advancedModeBox = new QGroupBox(tr("Advanced Segment Definition"));
  advancedModeBox->setVisible(false);

  operationCombo = new QComboBox;
  operationCombo->addItem(tr("Sum"));
  operationCombo->addItem(tr("Ratio"));

  componentsList = new QListWidget;
  componentsList->setMaximumHeight(100);

  componentEntranceCombo = new PairKeyCombo(PairKeyCombo::Entrance);
  componentExitCombo = new PairKeyCombo(PairKeyCombo::Exit);
  connect(componentEntranceCombo, SIGNAL(keyChanged(int)), componentExitCombo, SLOT(setEntranceKey(int)));

  useFixedAngleCheck = new QCheckBox(tr("Use fixed angle for denominator"));
  useFixedAngleCheck->setVisible(false);
  connect(useFixedAngleCheck, SIGNAL(stateChanged(int)), this, SLOT(useFixedAngleChanged(int)));

  fixedAngleLabel = new QLabel(tr("Angle [degrees]:"));
  fixedAngleLabel->setVisible(false);
  fixedAngleText = new QLineEdit;
  fixedAngleText->setText("90.0");
  fixedAngleText->setVisible(false);
  fixedAngleText->setMaximumWidth(80);

  componentScalingLabel = new QLabel(tr("Scaling:"));
  componentScalingText = new QLineEdit;
  componentScalingText->setText("1.0");
  componentScalingText->setMaximumWidth(80);

  addComponentButton = new QPushButton(tr("Add Component"));
  connect(addComponentButton, SIGNAL(clicked()), this, SLOT(addComponent()));

  removeComponentButton = new QPushButton(tr("Remove Component"));
  connect(removeComponentButton, SIGNAL(clicked()), this, SLOT(removeComponent()));

  cancelButton = new QPushButton(tr("Cancel"));
  okButton = new QPushButton(tr("Accept"));
  okButton->setDefault(true);

  QGroupBox *valueBox = new QGroupBox;
  QGridLayout *valueLayout = new QGridLayout;
  QGridLayout *pairLayout = new QGridLayout;
  pairLayout->addWidget(new QLabel(tr("Entrance Pair:")), 0, 0, Qt::AlignRight);
  pairLayout->addWidget(entrancePairCombo, 0, 1);
  pairLayout->addWidget(new QLabel(tr("Exit Pair:")), 0, 2, Qt::AlignRight);
  pairLayout->addWidget(exitPairCombo, 0, 3);
  totalCaptureLabel = new QLabel(tr("Total Capture"));
  totalCaptureLabel->setVisible(false);
  pairLayout->addWidget(totalCaptureLabel, 0, 4);
  valueLayout->addLayout(pairLayout, 0, 0, 1, 2);
  QGroupBox *energyBox = new QGroupBox(tr("Lab Energy [MeV]"));
  QGridLayout *energyLayout = new QGridLayout;
  energyLayout->addWidget(new QLabel(tr("Low Energy:")), 0, 0, Qt::AlignRight);
  energyLayout->addWidget(lowEnergyText, 0, 1);
  energyLayout->addWidget(new QLabel(tr("High Energy:")), 1, 0, Qt::AlignRight);
  energyLayout->addWidget(highEnergyText, 1, 1);
  energyBox->setLayout(energyLayout);
  valueLayout->addWidget(energyBox, 1, 0);
  QGroupBox *angleBox = new QGroupBox(tr("Lab Angle [degrees]"));
  QGridLayout *angleLayout = new QGridLayout;
  angleLayout->addWidget(new QLabel(tr("Low Angle:")), 0, 0, Qt::AlignRight);
  angleLayout->addWidget(lowAngleText, 0, 1);
  angleLayout->addWidget(new QLabel(tr("High Angle:")), 1, 0, Qt::AlignRight);
  angleLayout->addWidget(highAngleText, 1, 1);
  angleBox->setLayout(angleLayout);
  valueLayout->addWidget(angleBox, 1, 1);

  QGridLayout *lowerLayout = new QGridLayout;
  lowerLayout->addWidget(new QLabel(tr("Data Type:")), 0, 0, Qt::AlignRight);
  lowerLayout->addWidget(dataTypeCombo, 0, 1);

  QGridLayout *phaseLayout = new QGridLayout;
  phaseLayout->addItem(new QSpacerItem(1, 25), 0, 0);
  phaseLayout->setColumnStretch(0, 1);
  phaseJValueLabel = new QLabel(tr("J:"));
  phaseJValueLabel->setVisible(false);
  phaseLayout->addWidget(phaseJValueLabel, 0, 1);
  phaseLayout->addWidget(phaseJValueText, 0, 2);
  phaseLValueLabel = new QLabel(tr("l:"));
  phaseLValueLabel->setVisible(false);
  phaseLayout->addWidget(phaseLValueLabel, 0, 3);
  phaseLayout->addWidget(phaseLValueText, 0, 4);
  lowerLayout->addLayout(phaseLayout, 0, 2);

  lowerLayout->addWidget(new QLabel(tr("Data Norm.:")), 1, 0, Qt::AlignRight);
  QHBoxLayout *normLayout = new QHBoxLayout;
  normLayout->addWidget(dataNormText);
  normLayout->addWidget(varyNormCheck);
  lowerLayout->addLayout(normLayout, 1, 1);
  QGridLayout *normErrorLayout = new QGridLayout;
  normErrorLayout->addItem(new QSpacerItem(1, 25), 0, 0);
  normErrorLayout->setColumnStretch(0, 1);
  normErrorLayout->addWidget(dataNormErrorLabel, 0, 1, Qt::AlignRight);
  normErrorLayout->addWidget(dataNormErrorText, 0, 2);
  lowerLayout->addLayout(normErrorLayout, 1, 2);

  lowerLayout->addWidget(new QLabel(tr("Data File:")), 2, 0, Qt::AlignRight);
  QGridLayout *fileLayout = new QGridLayout;
  fileLayout->addWidget(dataFileText, 0, 0);
  fileLayout->addWidget(chooseFileButton, 0, 1);
  fileLayout->setColumnStretch(0, 1);
  lowerLayout->addLayout(fileLayout, 2, 1, 1, 2);

  lowerLayout->addWidget(energyShiftLabel, 3, 0, Qt::AlignRight);
  QHBoxLayout *energyShiftLayout = new QHBoxLayout;
  energyShiftLayout->addWidget(energyShiftText);
  energyShiftLayout->addWidget(varyEnergyShiftCheck);
  lowerLayout->addLayout(energyShiftLayout, 3, 1);
  QGridLayout *energyShiftErrorLayout = new QGridLayout;
  energyShiftErrorLayout->addItem(new QSpacerItem(1, 25), 0, 0);
  energyShiftErrorLayout->setColumnStretch(0, 1);
  energyShiftErrorLayout->addWidget(energyShiftErrorLabel, 0, 1, Qt::AlignRight);
  energyShiftErrorLayout->addWidget(energyShiftErrorText, 0, 2);
  lowerLayout->addLayout(energyShiftErrorLayout, 3, 2);

  // UPOS row
  lowerLayout->addWidget(uposCheck, 4, 0, 1, 3);
  QGridLayout *uposParamLayout = new QGridLayout;
  uposParamLayout->addItem(new QSpacerItem(1, 5), 0, 0);
  uposParamLayout->setColumnStretch(0, 1);
  uposParamLayout->addWidget(secondaryLLabel, 0, 1, Qt::AlignRight);
  uposParamLayout->addWidget(secondaryLSpin, 0, 2);
  uposParamLayout->addWidget(finalJLabel, 0, 3, Qt::AlignRight);
  uposParamLayout->addWidget(finalJText, 0, 4);
  uposParamLayout->addWidget(deltaLabel, 0, 5, Qt::AlignRight);
  uposParamLayout->addWidget(deltaText, 0, 6);
  lowerLayout->addLayout(uposParamLayout, 5, 0, 1, 3);

  lowerLayout->addWidget(advancedModeCheck, 6, 0, 1, 3);

  QGridLayout *advancedLayout = new QGridLayout;
  advancedLayout->addWidget(new QLabel(tr("Operation:")), 0, 0);
  advancedLayout->addWidget(operationCombo, 0, 1);
  connect(operationCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(operationTypeChanged(int)));

  advancedLayout->addWidget(new QLabel(tr("Components:")), 1, 0);
  advancedLayout->addWidget(componentsList, 1, 1, 4, 1);

  // Add component entrance/exit spinboxes above buttons
  advancedLayout->addWidget(new QLabel(tr("Entrance Pair:")), 1, 2);
  advancedLayout->addWidget(componentEntranceCombo, 1, 3);
  advancedLayout->addWidget(new QLabel(tr("Exit Pair:")), 2, 2);
  advancedLayout->addWidget(componentExitCombo, 2, 3);

  // Add fixed angle checkbox and input
  advancedLayout->addWidget(useFixedAngleCheck, 3, 2, 1, 2);
  advancedLayout->addWidget(fixedAngleLabel, 4, 2);
  advancedLayout->addWidget(fixedAngleText, 4, 3);

  // Per-component scaling input (multiplies this component before sum)
  advancedLayout->addWidget(componentScalingLabel, 5, 2);
  advancedLayout->addWidget(componentScalingText, 5, 3);

  QGridLayout *buttonLayout = new QGridLayout;
  buttonLayout->addWidget(addComponentButton, 0, 0);
  buttonLayout->addWidget(removeComponentButton, 0, 1);
  advancedLayout->addLayout(buttonLayout, 6, 2, 1, 2);
  advancedModeBox->setLayout(advancedLayout);

  lowerLayout->addWidget(advancedModeBox, 7, 0, 1, 3);
  lowerLayout->addWidget(naturalLabel, 8, 0, Qt::AlignRight);
  lowerLayout->addWidget(naturalCombo, 8, 1, 1, 2);

  valueLayout->addLayout(lowerLayout, 2, 0, 1, 2);
  valueBox->setLayout(valueLayout);

  QHBoxLayout *buttonBox = new QHBoxLayout;
  buttonBox->addWidget(cancelButton);
  buttonBox->addWidget(okButton);

  // The options outgrew small screens (UPOS, composites, mixture): the body
  // scrolls, the buttons stay put, and the dialog never exceeds the screen.
  DialogScrollArea *scroll = new DialogScrollArea(valueBox);
  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(scroll, 1);
  mainLayout->addLayout(buttonBox);

  setLayout(mainLayout);
  capDialogToScreen(this);

  connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(reject()));

  setWindowTitle(tr("Add a Segment From Data"));
}

void AddSegDataDialog::setChooseFile() {
  QString filename = QFileDialog::getOpenFileName(this);
  if (!filename.isEmpty()) {
    dataFileText->setText(QDir::fromNativeSeparators(filename));
  }
}

int AddSegDataDialog::dataTypeCode() const {
  const QVariant code = dataTypeCombo->currentData();
  return code.isValid() ? code.toInt() : dataTypeCombo->currentIndex();
}

void AddSegDataDialog::setDataTypeCode(int code) {
  const int i = dataTypeCombo->findData(code);
  dataTypeCombo->setCurrentIndex(i >= 0 ? i : code);
}

void AddSegDataDialog::dataTypeChanged(int index) {
  if (index == 2) {
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
  if (index == 1) {
    lowAngleText->setEnabled(true);
    highAngleText->setEnabled(true);
  } else {
    lowAngleText->setEnabled(false);
    highAngleText->setEnabled(false);
  }
  if (index == 3) {
    exitPairCombo->setVisible(false);
    totalCaptureLabel->setVisible(true);
  } else {
    totalCaptureLabel->setVisible(false);
    exitPairCombo->setVisible(true);
  }
  if (index == 4 || index == 7) {  // C.M. differential or A_y
    lowAngleText->setEnabled(true);
    highAngleText->setEnabled(true);
  }
  // A_y is a ratio: a normalization factor has no meaning for it.
  varyNormCheck->setEnabled(index != 7);
  if (index == 7) varyNormCheck->setChecked(false);

  // Update fixed angle visibility when data type changes
  operationTypeChanged(operationCombo->currentIndex());
  updateNaturalCheck();
}

void AddSegDataDialog::varyNormChanged(int state) {
  if (state == Qt::Checked) {
    dataNormErrorLabel->setVisible(true);
    dataNormErrorText->setVisible(true);
  } else {
    dataNormErrorLabel->setVisible(false);
    dataNormErrorText->setVisible(false);
  }
}

void AddSegDataDialog::advancedModeChanged(int state) {
  if (state == Qt::Checked) {
    advancedModeBox->setVisible(true);
    // Keep the spinboxes enabled so users can select pairs to add as components
    entrancePairCombo->setEnabled(true);
    exitPairCombo->setEnabled(true);
  } else {
    advancedModeBox->setVisible(false);
    entrancePairCombo->setEnabled(true);
    exitPairCombo->setEnabled(true);
  }
  fitDialogToContent(this);
}

void AddSegDataDialog::operationTypeChanged(int index) {
  // Show fixed angle option only when Ratio is selected AND segment is differential
  bool isRatio = (index == 1);
  bool isDifferential = (dataTypeCombo->currentIndex() == 1) || (dataTypeCombo->currentIndex() == 4);  // Differential or CM Differential

  if (isRatio && isDifferential) {
    useFixedAngleCheck->setVisible(true);
  } else {
    useFixedAngleCheck->setVisible(false);
    useFixedAngleCheck->setChecked(false);  // Uncheck if hidden
  }

  // Per-component scaling factor only applies to Sum mode
  componentScalingLabel->setVisible(!isRatio);
  componentScalingText->setVisible(!isRatio);
  updateNaturalCheck();
}

void AddSegDataDialog::useFixedAngleChanged(int state) {
  if (state == Qt::Checked) {
    fixedAngleLabel->setVisible(true);
    fixedAngleText->setVisible(true);
  } else {
    fixedAngleLabel->setVisible(false);
    fixedAngleText->setVisible(false);
  }
}

void AddSegDataDialog::addComponent() {
  int entrance = componentEntranceCombo->key();
  int exit = componentExitCombo->key();
  if (entrance < 1 || exit < 1) return;
  QString component = QString("Entrance: %1, Exit: %2").arg(entrance).arg(exit);

  // For ratio, only allow one component (the denominator)
  bool isRatio = (operationCombo->currentIndex() == 1);
  bool useFixedAngle = useFixedAngleCheck->isChecked();

  if (isRatio) {
    // For ratio, this is the denominator (only one component needed)
    if (useFixedAngle) {
      double angle = fixedAngleText->text().toDouble();
      component += QString(", Angle: %1").arg(angle);
    }

    // Only allow one component for ratio
    if (componentsList->count() >= 1) {
      return;  // Don't add more than one component for ratio
    }
  } else {
    // Only sums benefit from a scaling factor (applied to the component before
    // it is added to the base segment). Skip when the value is exactly 1.0 so
    // legacy files keep their existing component-string format.
    bool scalingOk = false;
    double scaling = componentScalingText->text().toDouble(&scalingOk);
    if (scalingOk && scaling != 1.0) {
      component += QString(", Scaling: %1").arg(scaling);
    }
  }

  componentsList->addItem(component);
}

void AddSegDataDialog::removeComponent() {
  int currentRow = componentsList->currentRow();
  if (currentRow >= 0) {
    delete componentsList->takeItem(currentRow);
  }
}

void AddSegDataDialog::uposChanged(int state) {
  bool visible = (state == Qt::Checked);
  secondaryLLabel->setVisible(visible);
  secondaryLSpin->setVisible(visible);
  finalJLabel->setVisible(visible);
  finalJText->setVisible(visible);
  deltaLabel->setVisible(visible);
  deltaText->setVisible(visible);
  fitDialogToContent(this);
}

bool AddSegDataDialog::naturalAllowed() const {
  // The engine's NaturalAllowed: a weighted sum is only meaningful for additive observables.
  int code = dataTypeCode();
  if (code == 2 || code == 7) return false;
  if (advancedModeCheck->isChecked() && operationCombo->currentIndex() == 1) return false;  // ratio
  return true;
}

void AddSegDataDialog::updateNaturalCheck() {
  bool allowed = naturalAllowed();
  naturalCombo->setEnabled(allowed);
  if (!allowed) naturalCombo->setCurrentIndex(0);
}

void AddSegDataDialog::setSamples(const QStringList &samples) {
  QString keep = naturalToken();
  naturalCombo->blockSignals(true);
  naturalCombo->clear();
  naturalCombo->addItem(tr("No"), QString());
  for (int s = 0; s < samples.size(); s++) {
    if (samples.at(s).isEmpty())
      naturalCombo->addItem(tr("Yes: default sample"), QString("1"));
    else
      naturalCombo->addItem(tr("Yes: sample %1").arg(samples.at(s)), samples.at(s));
  }
  if (!samples.contains(QString()) && naturalCombo->findData(QString("1")) < 0 && samples.size() == 1)
    naturalCombo->setItemData(1, QString("1"));   // a single named sample is also what "natural 1" means
  naturalCombo->blockSignals(false);
  setNaturalToken(keep);
}

QString AddSegDataDialog::naturalToken() const {
  return naturalCombo->currentData().toString();
}

void AddSegDataDialog::setNaturalToken(const QString &token) {
  int at = naturalCombo->findData(token);
  if (at < 0 && !token.isEmpty()) {
    // a sample the table does not list (yet): keep it selectable so an edit does not lose it
    naturalCombo->addItem(tr("Yes: sample %1 (not in the Pairs tab)").arg(token), token);
    at = naturalCombo->count() - 1;
  }
  naturalCombo->setCurrentIndex(std::max(0, at));
}

void AddSegDataDialog::setPairsModel(PairsModel *model) {
  // The check box only means something when the pairs form more than one compound nucleus.
  naturalLabel->setVisible(model && model->nuclei().size() > 1);
  naturalCombo->setVisible(model && model->nuclei().size() > 1);
  updateNaturalCheck();
  for (PairKeyCombo *combo : {entrancePairCombo, exitPairCombo, componentEntranceCombo, componentExitCombo}) combo->setPairsModel(model);
  // A fresh dialog starts on the first pair, as the key boxes started on 1.
  if (entrancePairCombo->key() == 0) entrancePairCombo->setKey(1);
  if (exitPairCombo->key() == 0) exitPairCombo->setKey(1);
  if (componentEntranceCombo->key() == 0) componentEntranceCombo->setKey(1);
  if (componentExitCombo->key() == 0) componentExitCombo->setKey(1);
}

void AddSegDataDialog::showEvent(QShowEvent *event) {
  QDialog::showEvent(event);
  fitDialogToContent(this);
}
