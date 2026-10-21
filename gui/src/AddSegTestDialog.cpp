#include "AddSegTestDialog.h"

#include <QPushButton>
#include <QGroupBox>
#include <QGridLayout>
#include "DialogScrollArea.h"
#include "PairsModel.h"
#include <QHBoxLayout>
#include <QStandardItemModel>


AddSegTestDialog::AddSegTestDialog(QWidget *parent) :
  QDialog(parent) {
  //  this->setMaximumSize(370,420);
  //  this->setMinimumSize(370,420);

  entrancePairCombo = new PairKeyCombo(PairKeyCombo::Entrance);
  exitPairCombo = new PairKeyCombo(PairKeyCombo::Exit);
  // The exit list reads as reactions with the chosen entrance pair.
  connect(entrancePairCombo, SIGNAL(keyChanged(int)), exitPairCombo, SLOT(setEntranceKey(int)));
  // Total capture (isDiff 4) is the last exit item when the project has a gamma pair.
  exitPairCombo->setTotalCaptureOffered(true);
  connect(exitPairCombo, SIGNAL(keyChanged(int)), this, SLOT(exitPairChanged(int)));
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
  // The observable. Each item carries its base isDiff code as item data, and items are only ever looked up
  // by that code (findData), never by position. The frame of a differential cross section (Lab 1, C.M. 5)
  // is a separate drop-down, and total capture (4) is an item of the exit-pair list.
  dataTypeCombo = new QComboBox;
  dataTypeCombo->setObjectName("dataTypeCombo");
  dataTypeCombo->addItem(tr("Angle Integrated"), 0);
  dataTypeCombo->addItem(tr("Differential"), 1);
  dataTypeCombo->addItem(tr("Phase Shift"), 2);
  dataTypeCombo->addItem(tr("Angular Distribution Coefficients"), 3);
  dataTypeCombo->addItem(tr("Analyzing Power"), 7);
  connect(dataTypeCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(dataTypeChanged(int)));
  frameLabel = new QLabel(tr("Frame:"));
  frameLabel->setVisible(false);
  frameCombo = new QComboBox;
  frameCombo->setObjectName("frameCombo");
  frameCombo->addItem(tr("Lab"), 1);
  frameCombo->addItem(tr("C.M."), 5);
  frameCombo->setToolTip(tr("Frame of the calculated angles and cross sections"));
  frameCombo->setVisible(false);
  connect(frameCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(dataTypeChanged(int)));
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
  angleBox = new QGroupBox(tr("Lab Angle [degrees]"));
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
  thmCheck = new QCheckBox(tr("THM (modified R-matrix)"));
  thmCheck->setToolTip(tr("Compute this segment as a Trojan Horse Method "
                          "(half-off-shell) cross section. Requires the "
                          "entrance pair's THM binding energy (Configure > "
                          "THM Workspace, Channels)."));
  lowerLayout->addWidget(thmCheck, 0, 2);
  lowerLayout->addItem(new QSpacerItem(1, 25), 0, 2);
  lowerLayout->setColumnStretch(2, 1);

  QHBoxLayout *phaseLayout = new QHBoxLayout;
  phaseLayout->addWidget(frameLabel);
  phaseLayout->addWidget(frameCombo);
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
  dataTypeChanged(0);
}

int AddSegTestDialog::observableCode() const {
  const QVariant code = dataTypeCombo->currentData();
  return code.isValid() ? code.toInt() : 0;
}

int AddSegTestDialog::dataTypeCode() const {
  const int observable = observableCode();
  if (observable == 0 && exitPairCombo->key() == PairKeyCombo::TotalCapture) return 4;
  if (observable == 1) return frameCombo->currentData().toInt();
  return observable;
}

void AddSegTestDialog::setDataTypeCode(int code) {
  int observable = code;
  if (code == 4) observable = 0;
  if (code == 5) observable = 1;
  if (code == 4)
    exitPairCombo->setKey(PairKeyCombo::TotalCapture);
  else if (exitPairCombo->key() == PairKeyCombo::TotalCapture)
    exitPairCombo->setKey(0);  // only a total-capture line has no exit pair
  frameCombo->setCurrentIndex(frameCombo->findData(code == 5 ? 5 : 1));
  int i = dataTypeCombo->findData(observable);
  if (i < 0) {
    // a code this dialog does not offer: keep it selectable so an edit does not change it
    dataTypeCombo->addItem(tr("Observable %1").arg(code), code);
    i = dataTypeCombo->count() - 1;
  }
  dataTypeCombo->setCurrentIndex(i);
  dataTypeChanged(i);
}

void AddSegTestDialog::exitPairChanged(int key) {
  // Total capture is an angle-integrated observable only.
  if (key == PairKeyCombo::TotalCapture && observableCode() != 0) dataTypeCombo->setCurrentIndex(dataTypeCombo->findData(0));
  dataTypeChanged(dataTypeCombo->currentIndex());
}

void AddSegTestDialog::dataTypeChanged(int) {
  const int code = dataTypeCode();
  const int observable = observableCode();
  const bool phase = (code == 2), angDist = (code == 3);
  // Phase shifts and angular-distribution coefficients have no THM form (the
  // engine reads their extra fields only for an isDiff of exactly 2 or 3).
  thmCheck->setEnabled(!phase && !angDist);
  if (phase || angDist) thmCheck->setChecked(false);
  phaseJValueLabel->setVisible(phase);
  phaseLValueLabel->setVisible(phase);
  phaseJValueText->setVisible(phase);
  phaseLValueText->setVisible(phase);
  angDistLabel->setVisible(angDist);
  angDistSpin->setVisible(angDist);
  // lab or C.M. differential, or A_y
  const bool angles = (code == 1 || code == 5 || code == 7);
  lowAngleText->setEnabled(angles);
  highAngleText->setEnabled(angles);
  angleStepText->setEnabled(angles);
  // C.M. differential and A_y take centre-of-mass angles
  angleBox->setTitle((code == 5 || code == 7) ? tr("C.M. Angle [degrees]") : tr("Lab Angle [degrees]"));
  frameLabel->setVisible(observable == 1);
  frameCombo->setVisible(observable == 1);
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
  const int code = dataTypeCode();
  bool isDifferential = (code == 1) || (code == 5);  // Differential, lab or C.M.

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
  dataTypeChanged(dataTypeCombo->currentIndex());
}

void AddSegTestDialog::showEvent(QShowEvent *event) {
  QDialog::showEvent(event);
  fitDialogToContent(this);
}
