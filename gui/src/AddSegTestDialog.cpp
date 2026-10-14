#include "AddSegTestDialog.h"

#include <QPushButton>
#include <QGroupBox>
#include <QGridLayout>
#include "DialogScrollArea.h"
#include "PairsModel.h"
#include <QHBoxLayout>


AddSegTestDialog::AddSegTestDialog(QWidget *parent) :
  QDialog(parent) {
  //  this->setMaximumSize(370,420);
  //  this->setMinimumSize(370,420);

  entrancePairCombo = new PairKeyCombo(PairKeyCombo::Entrance);
  exitPairCombo = new PairKeyCombo(PairKeyCombo::Exit);
  // The exit list reads as reactions with the chosen entrance pair.
  connect(entrancePairCombo, SIGNAL(keyChanged(int)), exitPairCombo, SLOT(setEntranceKey(int)));
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
  dataTypeCombo->addItem(tr("Analyzing Power"));
  // Codes 0-5 happen to equal their position; the analyzing power is 7.
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

  advancedModeCheck = new QCheckBox(tr("Advanced Mode (Sum/Ratio)"));
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
  energyLayout->addWidget(new QLabel(tr("Energy Step:")), 2, 0, Qt::AlignRight);
  energyLayout->addWidget(energyStepText, 2, 1);
  energyBox->setLayout(energyLayout);
  valueLayout->addWidget(energyBox, 1, 0);
  QGroupBox *angleBox = new QGroupBox(tr("Lab Angle [degrees]"));
  QGridLayout *angleLayout = new QGridLayout;
  angleLayout->addWidget(new QLabel(tr("Low Angle:")), 0, 0, Qt::AlignRight);
  angleLayout->addWidget(lowAngleText, 0, 1);
  angleLayout->addWidget(new QLabel(tr("High Angle:")), 1, 0, Qt::AlignRight);
  angleLayout->addWidget(highAngleText, 1, 1);
  angleLayout->addWidget(new QLabel(tr("Angle Step:")), 2, 0, Qt::AlignRight);
  angleLayout->addWidget(angleStepText, 2, 1);
  angleBox->setLayout(angleLayout);
  valueLayout->addWidget(angleBox, 1, 1);

  QGridLayout *lowerLayout = new QGridLayout;
  lowerLayout->addWidget(new QLabel(tr("Data Type:")), 0, 0, Qt::AlignRight);
  lowerLayout->addWidget(dataTypeCombo, 0, 1);
  lowerLayout->addItem(new QSpacerItem(1, 25), 0, 2);
  lowerLayout->setColumnStretch(2, 1);

  QHBoxLayout *phaseLayout = new QHBoxLayout;
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
  lowerLayout->addLayout(phaseLayout, 0, 3);

  lowerLayout->addWidget(advancedModeCheck, 1, 0, 1, 4);

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

  lowerLayout->addWidget(advancedModeBox, 2, 0, 1, 4);
  lowerLayout->addWidget(naturalLabel, 3, 0, Qt::AlignRight);
  lowerLayout->addWidget(naturalCombo, 3, 1, 1, 3);

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

  setWindowTitle(tr("Add a Segment Without Data"));
}

int AddSegTestDialog::dataTypeCode() const {
  const QVariant code = dataTypeCombo->currentData();
  return code.isValid() ? code.toInt() : dataTypeCombo->currentIndex();
}

void AddSegTestDialog::setDataTypeCode(int code) {
  const int i = dataTypeCombo->findData(code);
  dataTypeCombo->setCurrentIndex(i >= 0 ? i : code);
}

void AddSegTestDialog::dataTypeChanged(int index) {
  if (index == 2) {
    angDistLabel->setVisible(false);
    angDistSpin->setVisible(false);
    phaseJValueLabel->setVisible(true);
    phaseLValueLabel->setVisible(true);
    phaseJValueText->setVisible(true);
    phaseLValueText->setVisible(true);
  } else if (index == 3) {
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
  if (index == 1) {
    lowAngleText->setEnabled(true);
    highAngleText->setEnabled(true);
    angleStepText->setEnabled(true);
  } else {
    lowAngleText->setEnabled(false);
    highAngleText->setEnabled(false);
    angleStepText->setEnabled(false);
  }
  if (index == 4) {
    exitPairCombo->setVisible(false);
    totalCaptureLabel->setVisible(true);
  } else {
    totalCaptureLabel->setVisible(false);
    exitPairCombo->setVisible(true);
  }
  if (index == 5 || index == 6) {  // C.M. differential, or the analyzing power
    lowAngleText->setEnabled(true);
    highAngleText->setEnabled(true);
    angleStepText->setEnabled(true);
  }

  // Update fixed angle visibility when data type changes
  operationTypeChanged(operationCombo->currentIndex());
  updateNaturalCheck();
}

void AddSegTestDialog::advancedModeChanged(int state) {
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

void AddSegTestDialog::operationTypeChanged(int index) {
  // Show fixed angle option only when Ratio is selected AND segment is differential
  bool isRatio = (index == 1);
  bool isDifferential = (dataTypeCombo->currentIndex() == 1) || (dataTypeCombo->currentIndex() == 5);  // Differential or CM Differential

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

void AddSegTestDialog::useFixedAngleChanged(int state) {
  if (state == Qt::Checked) {
    fixedAngleLabel->setVisible(true);
    fixedAngleText->setVisible(true);
  } else {
    fixedAngleLabel->setVisible(false);
    fixedAngleText->setVisible(false);
  }
}

void AddSegTestDialog::addComponent() {
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
    // Sum mode: append optional scaling so this component is multiplied
    // before it is added to the base segment's cross section.
    bool scalingOk = false;
    double scaling = componentScalingText->text().toDouble(&scalingOk);
    if (scalingOk && scaling != 1.0) {
      component += QString(", Scaling: %1").arg(scaling);
    }
  }

  componentsList->addItem(component);
}

void AddSegTestDialog::removeComponent() {
  int currentRow = componentsList->currentRow();
  if (currentRow >= 0) {
    delete componentsList->takeItem(currentRow);
  }
}

bool AddSegTestDialog::naturalAllowed() const {
  // The engine's NaturalAllowed: a weighted sum is only meaningful for additive observables.
  int code = dataTypeCode();
  if (code == 2 || code == 3 || code == 7) return false;
  if (advancedModeCheck->isChecked() && operationCombo->currentIndex() == 1) return false;  // ratio
  return true;
}

void AddSegTestDialog::updateNaturalCheck() {
  bool allowed = naturalAllowed();
  naturalCombo->setEnabled(allowed);
  if (!allowed) naturalCombo->setCurrentIndex(0);
}

void AddSegTestDialog::setSamples(const QStringList &samples) {
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

QString AddSegTestDialog::naturalToken() const {
  return naturalCombo->currentData().toString();
}

void AddSegTestDialog::setNaturalToken(const QString &token) {
  int at = naturalCombo->findData(token);
  if (at < 0 && !token.isEmpty()) {
    // a sample the table does not list (yet): keep it selectable so an edit does not lose it
    naturalCombo->addItem(tr("%1 (not in the Pairs tab)").arg(token), token);
    at = naturalCombo->count() - 1;
  }
  naturalCombo->setCurrentIndex(std::max(0, at));
}

void AddSegTestDialog::setPairsModel(PairsModel *model) {
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

void AddSegTestDialog::showEvent(QShowEvent *event) {
  QDialog::showEvent(event);
  fitDialogToContent(this);
}
