#include "NuclearPotentialTab.h"
#include "NuclearPotentialManager.h"
#include "Config.h"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QLabel>
#include <QGroupBox>
#include <QComboBox>
#include <QLineEdit>
#include <QPushButton>
#include <QTextStream>
#include <QMessageBox>
#include <QDoubleValidator>
#include <iostream>

NuclearPotentialTab::NuclearPotentialTab(QWidget* parent)
    : QWidget(parent) {
  createUI();
  loadCurrentSettings();
  setWindowTitle(tr("Nuclear Potential"));
}

void NuclearPotentialTab::createUI() {
  QVBoxLayout* mainLayout = new QVBoxLayout(this);

  // Potential Type Section
  QHBoxLayout* typeLayout = new QHBoxLayout;
  potentialTypeLabel_ = new QLabel(tr("Potential Type:"));
  potentialTypeCombo_ = new QComboBox;
  potentialTypeCombo_->addItem(tr("Woods-Saxon"));
  potentialTypeCombo_->addItem(tr("Gaussian"));
  typeLayout->addWidget(potentialTypeLabel_);
  typeLayout->addWidget(potentialTypeCombo_);
  typeLayout->addStretch();

  connect(potentialTypeCombo_, QOverload<int>::of(&QComboBox::currentIndexChanged),
          this, &NuclearPotentialTab::onPotentialTypeChanged);

  mainLayout->addLayout(typeLayout);
  mainLayout->addSpacing(15);

  // Woods-Saxon Group
  woodsSaxonGroup_ = new QGroupBox(tr("Woods-Saxon Potential Parameters"));
  QGridLayout* wsLayout = new QGridLayout(woodsSaxonGroup_);

  v0Label_ = new QLabel(tr("Depth (V₀):"));
  v0Input_ = new QLineEdit;
  v0Input_->setValidator(new QDoubleValidator(0.0, 1000.0, 3, this));
  v0UnitLabel_ = new QLabel(tr("MeV"));
  wsLayout->addWidget(v0Label_, 0, 0);
  wsLayout->addWidget(v0Input_, 0, 1);
  wsLayout->addWidget(v0UnitLabel_, 0, 2);

  rLabel_ = new QLabel(tr("Radius (R):"));
  rInput_ = new QLineEdit;
  rInput_->setValidator(new QDoubleValidator(0.1, 100.0, 3, this));
  rUnitLabel_ = new QLabel(tr("fm"));
  wsLayout->addWidget(rLabel_, 1, 0);
  wsLayout->addWidget(rInput_, 1, 1);
  wsLayout->addWidget(rUnitLabel_, 1, 2);

  aLabel_ = new QLabel(tr("Diffuseness (a):"));
  aInput_ = new QLineEdit;
  aInput_->setValidator(new QDoubleValidator(0.01, 10.0, 3, this));
  aUnitLabel_ = new QLabel(tr("fm"));
  wsLayout->addWidget(aLabel_, 2, 0);
  wsLayout->addWidget(aInput_, 2, 1);
  wsLayout->addWidget(aUnitLabel_, 2, 2);

  wsLayout->setColumnStretch(1, 1);

  connect(v0Input_, &QLineEdit::textChanged, this, &NuclearPotentialTab::onParameterChanged);
  connect(rInput_, &QLineEdit::textChanged, this, &NuclearPotentialTab::onParameterChanged);
  connect(aInput_, &QLineEdit::textChanged, this, &NuclearPotentialTab::onParameterChanged);

  mainLayout->addWidget(woodsSaxonGroup_);

  // Gaussian Group
  gaussianGroup_ = new QGroupBox(tr("Gaussian Potential Parameters"));
  gaussianGroup_->setVisible(false);
  QGridLayout* gaussLayout = new QGridLayout(gaussianGroup_);

  gaussV0Label_ = new QLabel(tr("Depth (V₀):"));
  gaussV0Input_ = new QLineEdit;
  gaussV0Input_->setValidator(new QDoubleValidator(0.0, 1000.0, 3, this));
  gaussV0UnitLabel_ = new QLabel(tr("MeV"));
  gaussLayout->addWidget(gaussV0Label_, 0, 0);
  gaussLayout->addWidget(gaussV0Input_, 0, 1);
  gaussLayout->addWidget(gaussV0UnitLabel_, 0, 2);

  gaussR0Label_ = new QLabel(tr("Width (r₀):"));
  gaussR0Input_ = new QLineEdit;
  gaussR0Input_->setValidator(new QDoubleValidator(0.1, 100.0, 3, this));
  gaussR0UnitLabel_ = new QLabel(tr("fm"));
  gaussLayout->addWidget(gaussR0Label_, 1, 0);
  gaussLayout->addWidget(gaussR0Input_, 1, 1);
  gaussLayout->addWidget(gaussR0UnitLabel_, 1, 2);

  gaussLayout->setColumnStretch(1, 1);

  connect(gaussV0Input_, &QLineEdit::textChanged, this, &NuclearPotentialTab::onParameterChanged);
  connect(gaussR0Input_, &QLineEdit::textChanged, this, &NuclearPotentialTab::onParameterChanged);

  mainLayout->addWidget(gaussianGroup_);

  mainLayout->addSpacing(15);

  // Control Buttons
  QHBoxLayout* buttonLayout = new QHBoxLayout;
  applyButton_ = new QPushButton(tr("Apply Settings"));
  resetButton_ = new QPushButton(tr("Reset to Default"));

  connect(applyButton_, &QPushButton::clicked, this, &NuclearPotentialTab::onApplySettings);
  connect(resetButton_, &QPushButton::clicked, this, &NuclearPotentialTab::onResetToDefault);

  buttonLayout->addStretch();
  buttonLayout->addWidget(applyButton_);
  buttonLayout->addWidget(resetButton_);

  mainLayout->addLayout(buttonLayout);
  mainLayout->addStretch();
}

void NuclearPotentialTab::onPotentialTypeChanged(int index) {
  if(index == 0) { // Woods-Saxon
    woodsSaxonGroup_->setVisible(true);
    gaussianGroup_->setVisible(false);
  } else { // Gaussian
    woodsSaxonGroup_->setVisible(false);
    gaussianGroup_->setVisible(true);
  }
}

void NuclearPotentialTab::updateParameterLabels() {
  // Labels are already set in createUI
}

void NuclearPotentialTab::loadCurrentSettings() {
  auto& manager = NuclearPotentialManager::instance();
  std::string potentialType = manager.getCurrentPotentialType();

  if(potentialType == "WoodsSaxon") {
    potentialTypeCombo_->setCurrentIndex(0);
    double V0, R, a;
    if(manager.getWoodsSaxonParameters(V0, R, a)) {
      v0Input_->setText(QString::number(V0, 'f', 3));
      rInput_->setText(QString::number(R, 'f', 3));
      aInput_->setText(QString::number(a, 'f', 3));
    }
  } else if(potentialType == "Gaussian") {
    potentialTypeCombo_->setCurrentIndex(1);
    double V0, r0;
    if(manager.getGaussianParameters(V0, r0)) {
      gaussV0Input_->setText(QString::number(V0, 'f', 3));
      gaussR0Input_->setText(QString::number(r0, 'f', 3));
    }
  }
}

bool NuclearPotentialTab::validateParameters() {
  if(potentialTypeCombo_->currentIndex() == 0) { // Woods-Saxon
    bool ok1, ok2, ok3;
    double v0 = v0Input_->text().toDouble(&ok1);
    double r = rInput_->text().toDouble(&ok2);
    double a = aInput_->text().toDouble(&ok3);

    if(!ok1 || !ok2 || !ok3) {
      QMessageBox::warning(this, tr("Invalid Input"),
                          tr("Please enter valid numbers for all parameters."));
      return false;
    }

    if(v0 <= 0.0 || r <= 0.0 || a <= 0.0) {
      QMessageBox::warning(this, tr("Invalid Parameters"),
                          tr("All parameters must be positive values."));
      return false;
    }
  } else { // Gaussian
    bool ok1, ok2;
    double v0 = gaussV0Input_->text().toDouble(&ok1);
    double r0 = gaussR0Input_->text().toDouble(&ok2);

    if(!ok1 || !ok2) {
      QMessageBox::warning(this, tr("Invalid Input"),
                          tr("Please enter valid numbers for all parameters."));
      return false;
    }

    if(v0 <= 0.0 || r0 <= 0.0) {
      QMessageBox::warning(this, tr("Invalid Parameters"),
                          tr("All parameters must be positive values."));
      return false;
    }
  }

  return true;
}

void NuclearPotentialTab::onApplySettings() {
  if(!validateParameters()) {
    return;
  }

  try {
    auto& manager = NuclearPotentialManager::instance();

    if(potentialTypeCombo_->currentIndex() == 0) { // Woods-Saxon
      double v0 = v0Input_->text().toDouble();
      double r = rInput_->text().toDouble();
      double a = aInput_->text().toDouble();
      manager.setWoodsSaxonPotential(v0, r, a);
    } else { // Gaussian
      double v0 = gaussV0Input_->text().toDouble();
      double r0 = gaussR0Input_->text().toDouble();
      manager.setGaussianPotential(v0, r0);
    }

    applyButton_->setEnabled(true);
  } catch(const std::exception& e) {
    QMessageBox::critical(this, tr("Error"),
                         tr("Failed to apply settings: %1").arg(e.what()));
  }
}

void NuclearPotentialTab::onResetToDefault() {
  try {
    auto& manager = NuclearPotentialManager::instance();
    manager.resetToDefault();
    loadCurrentSettings();
  } catch(const std::exception& e) {
    QMessageBox::critical(this, tr("Error"),
                         tr("Failed to reset: %1").arg(e.what()));
  }
}

void NuclearPotentialTab::onParameterChanged() {
}

bool NuclearPotentialTab::readPotentialSettings(QTextStream& inStream, Config& config) {
  // Read potential configuration from file
  // Format (inside <potential> section):
  // useHybridPotential=1
  // potentialType=0
  // V0=150.0
  // R=3.6
  // a=0.6
  // or
  // useHybridPotential=1
  // potentialType=1
  // V0=50.0
  // r0=5.0

  int typeCode = 0;
  double v0 = 150.0, r = 3.6, a = 0.6, r0 = 5.0;
  bool hasType = false;

  while(!inStream.atEnd()) {
    QString line = inStream.readLine();
    QString trimmedLine = line.trimmed();

    // Stop at end of potential section
    if(trimmedLine == "</potential>") {
      break;
    }

    if(trimmedLine.startsWith("useHybridPotential=")) {
      int value = trimmedLine.mid(19).trimmed().toInt();
      config.useHybridMethod = (value == 1);
    } else if(trimmedLine.startsWith("useAdaptiveGrid=")) {
      int value = trimmedLine.mid(16).trimmed().toInt();
      config.useAdaptiveGrid = (value == 1);
    } else if(trimmedLine.startsWith("potentialType=")) {
      typeCode = trimmedLine.mid(14).trimmed().toInt();
      hasType = true;
    } else if(trimmedLine.startsWith("V0=")) {
      v0 = trimmedLine.mid(3).trimmed().toDouble();
    } else if(trimmedLine.startsWith("R=")) {
      r = trimmedLine.mid(2).trimmed().toDouble();
    } else if(trimmedLine.startsWith("a=")) {
      a = trimmedLine.mid(2).trimmed().toDouble();
    } else if(trimmedLine.startsWith("r0=")) {
      r0 = trimmedLine.mid(3).trimmed().toDouble();
    }
  }

  if(!hasType) {
    return false;
  }

  try {
    auto& manager = NuclearPotentialManager::instance();

    if(typeCode == 0) {  // WoodsSaxon
      manager.setWoodsSaxonPotential(v0, r, a);
      potentialTypeCombo_->setCurrentIndex(0);
      v0Input_->setText(QString::number(v0, 'f', 3));
      rInput_->setText(QString::number(r, 'f', 3));
      aInput_->setText(QString::number(a, 'f', 3));
    } else if(typeCode == 1) {  // Gaussian
      manager.setGaussianPotential(v0, r0);
      potentialTypeCombo_->setCurrentIndex(1);
      gaussV0Input_->setText(QString::number(v0, 'f', 3));
      gaussR0Input_->setText(QString::number(r0, 'f', 3));
    } else {
      return false;
    }

    return true;
  } catch(...) {
    return false;
  }
}

bool NuclearPotentialTab::writePotentialSettings(QTextStream& outStream) {
  auto& manager = NuclearPotentialManager::instance();
  std::string potentialType = manager.getCurrentPotentialType();

  // Write potential type as integer: 0=WoodsSaxon, 1=Gaussian
  int typeCode = 0;
  if(potentialType == "Gaussian") typeCode = 1;
  outStream << "potentialType=" << typeCode << "\n";

  if(potentialType == "WoodsSaxon") {
    double V0, R, a;
    if(manager.getWoodsSaxonParameters(V0, R, a)) {
      outStream << "V0=" << V0 << "\n";
      outStream << "R=" << R << "\n";
      outStream << "a=" << a << "\n";
    }
  } else if(potentialType == "Gaussian") {
    double V0, r0;
    if(manager.getGaussianParameters(V0, r0)) {
      outStream << "V0=" << V0 << "\n";
      outStream << "r0=" << r0 << "\n";
    }
  }

  return true;
}
