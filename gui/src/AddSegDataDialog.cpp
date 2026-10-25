#include "AddSegDataDialog.h"

#include <QPushButton>
#include <QGroupBox>
#include <QGridLayout>
#include "DialogScrollArea.h"
#include "PairsModel.h"
#include <QFileDialog>
#include <QStandardItemModel>

AddSegDataDialog::AddSegDataDialog(QWidget *parent) :
  QDialog(parent) {
  //  this->setMaximumSize(370,420);
  //  this->setMinimumSize(370,420);

  entrancePairCombo = new PairKeyCombo(PairKeyCombo::Entrance);
  exitPairCombo = new PairKeyCombo(PairKeyCombo::Exit);
  // The exit list reads as reactions with the chosen entrance pair.
  connect(entrancePairCombo, SIGNAL(keyChanged(int)), exitPairCombo, SLOT(setEntranceKey(int)));
  // Total capture (isDiff 3) is the last exit item when the project has a gamma pair.
  exitPairCombo->setTotalCaptureOffered(true);
  connect(exitPairCombo, SIGNAL(keyChanged(int)), this, SLOT(exitPairChanged(int)));
  lowEnergyText = new QLineEdit;
  highEnergyText = new QLineEdit;
  lowAngleText = new QLineEdit;
  lowAngleText->setText("0");
  lowAngleText->setEnabled(false);
  highAngleText = new QLineEdit;
  highAngleText->setText("180");
  highAngleText->setEnabled(false);
  // The observable. Each item carries its base isDiff code as item data, and items are only ever looked up
  // by that code (findData), never by position. The variants are separate controls: the frame of a
  // differential cross section (Lab 1, C.M. 4), the multipolarity of an angle-integrated capture cross
  // section (all 0, E1 only 5, E2 only 6), and total capture (3), which is an item of the exit-pair list.
  dataTypeCombo = new QComboBox;
  dataTypeCombo->setObjectName("dataTypeCombo");
  dataTypeCombo->addItem(tr("Angle Integrated"), 0);
  dataTypeCombo->addItem(tr("Differential"), 1);
  dataTypeCombo->addItem(tr("Phase Shift"), 2);
  dataTypeCombo->addItem(tr("Analyzing Power"), 7);
  connect(dataTypeCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(dataTypeChanged(int)));
  frameLabel = new QLabel(tr("Frame:"));
  frameLabel->setVisible(false);
  frameCombo = new QComboBox;
  frameCombo->setObjectName("frameCombo");
  frameCombo->addItem(tr("Lab"), 1);
  frameCombo->addItem(tr("C.M."), 4);
  frameCombo->setToolTip(tr("Frame of the data file's angles and cross sections"));
  frameCombo->setVisible(false);
  connect(frameCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(dataTypeChanged(int)));
  multipolarityLabel = new QLabel(tr("Multipolarity:"));
  multipolarityLabel->setVisible(false);
  multipolarityCombo = new QComboBox;
  multipolarityCombo->setObjectName("multipolarityCombo");
  multipolarityCombo->addItem(tr("All"), 0);
  multipolarityCombo->addItem(tr("E1 only"), 5);
  multipolarityCombo->addItem(tr("E2 only"), 6);
  multipolarityCombo->setToolTip(tr("Fit the capture data to the E1 or E2 part of the calculated cross section only"));
  multipolarityCombo->setVisible(false);
  connect(multipolarityCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(dataTypeChanged(int)));
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
  covarianceFileText = new QLineEdit;
  covarianceFileText->setObjectName("covarianceFileText");
  covarianceFileText->setPlaceholderText(tr("none (uncorrelated errors)"));
  covarianceFileText->setToolTip(tr("Optional data covariance file (\"covariance <file>\" on the segment line).\n"
                                    "Every data segment naming the same file shares one covariance matrix, and its\n"
                                    "points enter chi-squared as r^T M^-1 r.  The file lists its data files and rows\n"
                                    "(format in the documentation, Data Formats).  The path is resolved like the data\n"
                                    "file's.  Leave empty for uncorrelated errors."));
  QPushButton *chooseCovarianceButton = new QPushButton(tr("Choose..."));
  chooseCovarianceButton->setObjectName("chooseCovarianceButton");
  connect(chooseCovarianceButton, SIGNAL(clicked()), this, SLOT(setChooseCovarianceFile()));
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

  naturalLabel = new QLabel(tr("Sample (abundance weights):"));
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
  valueLayout->addLayout(pairLayout, 0, 0, 1, 2);
  QGroupBox *energyBox = new QGroupBox(tr("Lab Energy [MeV]"));
  QGridLayout *energyLayout = new QGridLayout;
  energyLayout->addWidget(new QLabel(tr("Low Energy:")), 0, 0, Qt::AlignRight);
  energyLayout->addWidget(lowEnergyText, 0, 1);
  energyLayout->addWidget(new QLabel(tr("High Energy:")), 1, 0, Qt::AlignRight);
  energyLayout->addWidget(highEnergyText, 1, 1);
  energyBox->setLayout(energyLayout);
  valueLayout->addWidget(energyBox, 1, 0);
  angleBox = new QGroupBox(tr("Lab Angle [degrees]"));
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
  // the observable's variant (frame or multipolarity) sits next to it, the phase-shift J and l on the right
  phaseLayout->addWidget(frameLabel, 0, 0, Qt::AlignRight);
  phaseLayout->addWidget(frameCombo, 0, 1);
  phaseLayout->addWidget(multipolarityLabel, 0, 2, Qt::AlignRight);
  phaseLayout->addWidget(multipolarityCombo, 0, 3);
  phaseLayout->addItem(new QSpacerItem(1, 25), 0, 4);
  phaseLayout->setColumnStretch(4, 1);
  phaseJValueLabel = new QLabel(tr("J:"));
  phaseJValueLabel->setVisible(false);
  phaseLayout->addWidget(phaseJValueLabel, 0, 5);
  phaseLayout->addWidget(phaseJValueText, 0, 6);
  phaseLValueLabel = new QLabel(tr("l:"));
  phaseLValueLabel->setVisible(false);
  phaseLayout->addWidget(phaseLValueLabel, 0, 7);
  phaseLayout->addWidget(phaseLValueText, 0, 8);
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

  lowerLayout->addWidget(new QLabel(tr("Covariance File:")), 3, 0, Qt::AlignRight);
  QGridLayout *covarianceLayout = new QGridLayout;
  covarianceLayout->addWidget(covarianceFileText, 0, 0);
  covarianceLayout->addWidget(chooseCovarianceButton, 0, 1);
  covarianceLayout->setColumnStretch(0, 1);
  lowerLayout->addLayout(covarianceLayout, 3, 1, 1, 2);

  lowerLayout->addWidget(energyShiftLabel, 4, 0, Qt::AlignRight);
  QHBoxLayout *energyShiftLayout = new QHBoxLayout;
  energyShiftLayout->addWidget(energyShiftText);
  energyShiftLayout->addWidget(varyEnergyShiftCheck);
  lowerLayout->addLayout(energyShiftLayout, 4, 1);
  QGridLayout *energyShiftErrorLayout = new QGridLayout;
  energyShiftErrorLayout->addItem(new QSpacerItem(1, 25), 0, 0);
  energyShiftErrorLayout->setColumnStretch(0, 1);
  energyShiftErrorLayout->addWidget(energyShiftErrorLabel, 0, 1, Qt::AlignRight);
  energyShiftErrorLayout->addWidget(energyShiftErrorText, 0, 2);
  lowerLayout->addLayout(energyShiftErrorLayout, 4, 2);

  // UPOS row
  lowerLayout->addWidget(uposCheck, 5, 0, 1, 3);
  QGridLayout *uposParamLayout = new QGridLayout;
  uposParamLayout->addItem(new QSpacerItem(1, 5), 0, 0);
  uposParamLayout->setColumnStretch(0, 1);
  uposParamLayout->addWidget(secondaryLLabel, 0, 1, Qt::AlignRight);
  uposParamLayout->addWidget(secondaryLSpin, 0, 2);
  uposParamLayout->addWidget(finalJLabel, 0, 3, Qt::AlignRight);
  uposParamLayout->addWidget(finalJText, 0, 4);
  uposParamLayout->addWidget(deltaLabel, 0, 5, Qt::AlignRight);
  uposParamLayout->addWidget(deltaText, 0, 6);
  lowerLayout->addLayout(uposParamLayout, 6, 0, 1, 3);

  lowerLayout->addWidget(advancedModeCheck, 7, 0, 1, 3);

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

  lowerLayout->addWidget(advancedModeBox, 8, 0, 1, 3);
  lowerLayout->addWidget(naturalLabel, 9, 0, Qt::AlignRight);
  lowerLayout->addWidget(naturalCombo, 9, 1, 1, 2);

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
  dataTypeChanged(0);
}

void AddSegDataDialog::setChooseFile() {
  QString filename = QFileDialog::getOpenFileName(this);
  if (!filename.isEmpty()) {
    dataFileText->setText(QDir::fromNativeSeparators(filename));
  }
}

void AddSegDataDialog::setChooseCovarianceFile() {
  QString filename = QFileDialog::getOpenFileName(this, tr("Data Covariance File"));
  if (!filename.isEmpty()) {
    covarianceFileText->setText(QDir::fromNativeSeparators(filename));
  }
}

int AddSegDataDialog::observableCode() const {
  const QVariant code = dataTypeCombo->currentData();
  return code.isValid() ? code.toInt() : 0;
}

bool AddSegDataDialog::multipolarityApplies() const {
  return observableCode() == 0 && exitPairCombo->isGammaPair(exitPairCombo->key());
}

int AddSegDataDialog::dataTypeCode() const {
  const int observable = observableCode();
  if (observable == 0) {
    if (exitPairCombo->key() == PairKeyCombo::TotalCapture) return 3;
    if (multipolarityApplies()) return multipolarityCombo->currentData().toInt();
    return 0;
  }
  if (observable == 1) return frameCombo->currentData().toInt();
  return observable;
}

void AddSegDataDialog::setDataTypeCode(int code) {
  int observable = code;
  if (code == 3 || code == 5 || code == 6) observable = 0;
  if (code == 4) observable = 1;
  if (code == 3)
    exitPairCombo->setKey(PairKeyCombo::TotalCapture);
  else if (exitPairCombo->key() == PairKeyCombo::TotalCapture)
    exitPairCombo->setKey(0);  // only a total-capture line has no exit pair
  frameCombo->setCurrentIndex(frameCombo->findData(code == 4 ? 4 : 1));
  multipolarityCombo->setCurrentIndex(multipolarityCombo->findData(code == 5 || code == 6 ? code : 0));
  int i = dataTypeCombo->findData(observable);
  if (i < 0) {
    // a code this dialog does not offer: keep it selectable so an edit does not change it
    dataTypeCombo->addItem(tr("Observable %1").arg(code), code);
    i = dataTypeCombo->count() - 1;
  }
  dataTypeCombo->setCurrentIndex(i);
  dataTypeChanged(i);
}

void AddSegDataDialog::exitPairChanged(int key) {
  // Total capture is an angle-integrated observable only.
  if (key == PairKeyCombo::TotalCapture && observableCode() != 0) dataTypeCombo->setCurrentIndex(dataTypeCombo->findData(0));
  dataTypeChanged(dataTypeCombo->currentIndex());
}

void AddSegDataDialog::dataTypeChanged(int) {
  const int code = dataTypeCode();
  const int observable = observableCode();
  const bool phase = (code == 2);
  phaseJValueLabel->setVisible(phase);
  phaseLValueLabel->setVisible(phase);
  phaseJValueText->setVisible(phase);
  phaseLValueText->setVisible(phase);
  // lab or C.M. differential, or A_y
  const bool angles = (code == 1 || code == 4 || code == 7);
  lowAngleText->setEnabled(angles);
  highAngleText->setEnabled(angles);
  // C.M. differential and A_y take centre-of-mass angles
  angleBox->setTitle((code == 4 || code == 7) ? tr("C.M. Angle [degrees]") : tr("Lab Angle [degrees]"));
  frameLabel->setVisible(observable == 1);
  frameCombo->setVisible(observable == 1);
  const bool multipolarity = multipolarityApplies();
  multipolarityLabel->setVisible(multipolarity);
  multipolarityCombo->setVisible(multipolarity);
  // With "Total capture" as the exit pair only the angle-integrated observable is offered.
  const bool totalCapture = (exitPairCombo->key() == PairKeyCombo::TotalCapture);
  if (QStandardItemModel *items = qobject_cast<QStandardItemModel *>(dataTypeCombo->model())) {
    for (int i = 0; i < items->rowCount(); i++) {
      QStandardItem *item = items->item(i);
      const bool allowed = !totalCapture || item->data(Qt::UserRole).toInt() == 0;
      item->setEnabled(allowed);
      item->setToolTip(allowed ? QString() : tr("Total capture is angle integrated: choose a single exit pair for this observable"));
    }
  }
  // A_y is a ratio: a normalization factor has no meaning for it.
  varyNormCheck->setEnabled(code != 7);
  if (code == 7) varyNormCheck->setChecked(false);

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
  const int code = dataTypeCode();
  bool isDifferential = (code == 1) || (code == 4);  // Differential, lab or C.M.

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
  naturalCombo->addItem(tr("none (single nucleus, unweighted)"), QString());
  for (int s = 0; s < samples.size(); s++) {
    if (samples.at(s).isEmpty())
      naturalCombo->addItem(tr("default sample"), QString("1"));
    else
      naturalCombo->addItem(samples.at(s), samples.at(s));
  }
  if (!samples.contains(QString()) && naturalCombo->findData(QString("1")) < 0 && samples.size() == 1)
    naturalCombo->setItemData(1, QString("1"));   // a single named sample is also what "sample 1" means
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
    naturalCombo->addItem(tr("%1 (not in the Pairs tab)").arg(token), token);
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
  dataTypeChanged(dataTypeCombo->currentIndex());
}

void AddSegDataDialog::showEvent(QShowEvent *event) {
  QDialog::showEvent(event);
  fitDialogToContent(this);
}
