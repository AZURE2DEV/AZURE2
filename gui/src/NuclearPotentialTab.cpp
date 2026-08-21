#include "NuclearPotentialTab.h"
#include "NuclearPotentialManager.h"
#include "Config.h"
#include "PairsModel.h"
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
#include <QCheckBox>
#include <iostream>

NuclearPotentialTab::NuclearPotentialTab(QWidget* parent)
    : QWidget(parent) {
  createUI();
  loadCurrentSettings();
  setWindowTitle(tr("Nuclear Potential"));
}

void NuclearPotentialTab::createUI() {
  QVBoxLayout* mainLayout = new QVBoxLayout(this);

  // Which pair is being edited.  A potential belongs to a pair, so the tab
  // shows one pair at a time; "Default" is what a pair inherits when it has no
  // setting of its own.
  QHBoxLayout* pairLayout = new QHBoxLayout;
  pairLabel_ = new QLabel(tr("Applies to:"));
  pairCombo_ = new QComboBox;
  enabledCheck_ = new QCheckBox(tr("Use hybrid potential here"));
  useDefaultButton_ = new QPushButton(tr("Follow Default"));
  useDefaultButton_->setToolTip(tr("Drop this pair's own setting so it follows "
                                   "the default again."));
  pairLayout->addWidget(pairLabel_);
  pairLayout->addWidget(pairCombo_, 1);
  pairLayout->addWidget(enabledCheck_);
  pairLayout->addWidget(useDefaultButton_);

  connect(pairCombo_, QOverload<int>::of(&QComboBox::currentIndexChanged),
          this, &NuclearPotentialTab::onPairSelectionChanged);
  connect(enabledCheck_, &QCheckBox::toggled,
          this, &NuclearPotentialTab::onEnabledToggled);
  connect(useDefaultButton_, &QPushButton::clicked,
          this, &NuclearPotentialTab::onUseDefaultForPair);

  mainLayout->addLayout(pairLayout);
  mainLayout->addSpacing(10);

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

  summaryLabel_ = new QLabel;
  summaryLabel_->setWordWrap(true);
  mainLayout->addWidget(summaryLabel_);

  mainLayout->addStretch();
}

void NuclearPotentialTab::setPairsModel(PairsModel* model) {
  pairsModel_ = model;
  refreshPairCombo();
}

void NuclearPotentialTab::refreshPairCombo() {
  loading_ = true;
  int keep = currentPairKey_;
  pairCombo_->clear();
  pairCombo_->addItem(tr("Default (every pair without its own)"), 0);
  if(pairsModel_) {
    QList<PairsData> pairs = pairsModel_->getPairs();
    for(int i = 0; i < pairs.size(); i++) {
      // A pair's key is its 1-based position, matching the .azr and
      // PPair::GetPairKey().
      pairCombo_->addItem(tr("Pair %1: %2").arg(i + 1)
                            .arg(pairsModel_->getParticleLabel(pairs.at(i))),
                          i + 1);
    }
  }
  int index = pairCombo_->findData(keep);
  pairCombo_->setCurrentIndex(index >= 0 ? index : 0);
  currentPairKey_ = pairCombo_->currentData().toInt();
  loading_ = false;
  loadSettingsFor(currentPairKey_);
}

void NuclearPotentialTab::onPairSelectionChanged(int index) {
  if(loading_ || index < 0) return;
  int next = pairCombo_->itemData(index).toInt();
  if(next == currentPairKey_) return;
  // Commit what is on screen to the pair being left, so switching the selector
  // never silently throws typing away.  If it does not validate, stay put.
  if(!commitCurrent()) {
    loading_ = true;
    pairCombo_->setCurrentIndex(pairCombo_->findData(currentPairKey_));
    loading_ = false;
    return;
  }
  currentPairKey_ = next;
  loadSettingsFor(currentPairKey_);
}

void NuclearPotentialTab::onEnabledToggled(bool) {
  if(loading_) return;
  commitCurrent();
  refreshSummary();
}

void NuclearPotentialTab::onUseDefaultForPair() {
  if(currentPairKey_ == 0) {
    QMessageBox::information(this, tr("Nuclear Potential"),
                             tr("This is the default; the pairs follow it."));
    return;
  }
  NuclearPotentialManager::instance().clearPairSetting(currentPairKey_);
  loadSettingsFor(currentPairKey_);
  refreshSummary();
}

void NuclearPotentialTab::loadSettingsFor(int pairKey) {
  loading_ = true;
  NuclearPotentialManager& manager = NuclearPotentialManager::instance();
  NuclearPotentialSetting s = pairKey ? manager.getSetting(pairKey)
                                      : manager.getDefaultSetting();
  enabledCheck_->setChecked(s.enabled);
  if(s.type == "Gaussian") {
    potentialTypeCombo_->setCurrentIndex(1);
    gaussV0Input_->setText(QString::number(s.V0, 'f', 3));
    gaussR0Input_->setText(QString::number(s.r0, 'f', 3));
  } else {
    potentialTypeCombo_->setCurrentIndex(0);
    v0Input_->setText(QString::number(s.V0, 'f', 3));
    rInput_->setText(QString::number(s.R, 'f', 3));
    aInput_->setText(QString::number(s.a, 'f', 3));
  }
  onPotentialTypeChanged(potentialTypeCombo_->currentIndex());
  bool own = pairKey && manager.hasPairSetting(pairKey);
  useDefaultButton_->setEnabled(own);
  loading_ = false;
  refreshSummary();
}

NuclearPotentialSetting NuclearPotentialTab::settingFromWidgets() const {
  NuclearPotentialSetting s;
  s.enabled = enabledCheck_->isChecked();
  if(potentialTypeCombo_->currentIndex() == 1) {
    s.type = "Gaussian";
    s.V0 = gaussV0Input_->text().toDouble();
    s.r0 = gaussR0Input_->text().toDouble();
  } else {
    s.type = "WoodsSaxon";
    s.V0 = v0Input_->text().toDouble();
    s.R = rInput_->text().toDouble();
    s.a = aInput_->text().toDouble();
  }
  return s;
}

bool NuclearPotentialTab::commitCurrent() {
  if(!validateParameters()) return false;
  try {
    NuclearPotentialManager& manager = NuclearPotentialManager::instance();
    NuclearPotentialSetting s = settingFromWidgets();
    if(currentPairKey_) manager.setSetting(currentPairKey_, s);
    else manager.setDefaultSetting(s);
    useDefaultButton_->setEnabled(currentPairKey_ != 0);
    return true;
  } catch(const std::exception& e) {
    QMessageBox::critical(this, tr("Error"),
                          tr("Failed to apply settings: %1").arg(e.what()));
    return false;
  }
}

void NuclearPotentialTab::refreshSummary() {
  NuclearPotentialManager& manager = NuclearPotentialManager::instance();
  QStringList on;
  if(pairsModel_) {
    QList<PairsData> pairs = pairsModel_->getPairs();
    for(int i = 0; i < pairs.size(); i++)
      if(manager.isPairEnabled(i + 1)) on << QString::number(i + 1);
  }
  QString text;
  if(!pairsModel_ || pairsModel_->getPairs().isEmpty())
    text = tr("Default is %1.").arg(manager.getDefaultEnabled() ? tr("on") : tr("off"));
  else if(on.isEmpty())
    text = tr("The hybrid potential is off for every pair.");
  else
    text = tr("Hybrid potential active for pair(s): %1.").arg(on.join(", "));
  summaryLabel_->setText(text);
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
  if(commitCurrent()) refreshSummary();
}

void NuclearPotentialTab::onResetToDefault() {
  try {
    auto& manager = NuclearPotentialManager::instance();
    manager.resetToDefault();          // also drops every per-pair setting
    currentPairKey_ = 0;
    refreshPairCombo();
  } catch(const std::exception& e) {
    QMessageBox::critical(this, tr("Error"),
                         tr("Failed to reset: %1").arg(e.what()));
  }
}

void NuclearPotentialTab::onParameterChanged() {
}

bool NuclearPotentialTab::readPotentialSettings(QTextStream& inStream, Config& config) {
  // Read the body of the <potential> section.  Keys before the first pair= are
  // the default that unnamed pairs inherit; a pair= line opens a section that
  // starts from that default and only states what differs.  A file with no
  // pair= line therefore reads exactly as it did when the model was global.
  // Config::ReadPotentialBlock parses the same format for --no-gui and pyazr.
  NuclearPotentialManager& manager = NuclearPotentialManager::instance();
  manager.resetToDefault();

  NuclearPotentialSetting current = manager.getDefaultSetting();
  int currentPair = 0;
  bool currentHasType = false, defaultHasType = false, sawAny = false;

  auto flush = [&]() {
    if(current.enabled && !currentHasType) current.enabled = false;
    try {
      if(currentPair) manager.setSetting(currentPair, current);
      else manager.setDefaultSetting(current);
    } catch(...) {
      if(currentPair) manager.clearPairSetting(currentPair);
      else manager.setDefaultEnabled(false);
    }
  };

  while(!inStream.atEnd()) {
    QString trimmedLine = inStream.readLine().trimmed();
    if(trimmedLine == "</potential>") break;

    int eq = trimmedLine.indexOf('=');
    if(eq < 0) continue;
    QString key = trimmedLine.left(eq);
    QString value = trimmedLine.mid(eq + 1).trimmed();

    if(key == "pair") {
      int next = value.toInt();
      if(next <= 0) continue;
      flush();
      currentPair = next;
      current = manager.getDefaultSetting();
      currentHasType = defaultHasType;
      continue;
    }

    if(key == "useHybridPotential") {
      current.enabled = (value.toInt() == 1);
      // The default section's flag is also the master switch the Runtime
      // Options checkbox shows.
      if(!currentPair) config.useHybridMethod = current.enabled;
    } else if(key == "useAdaptiveGrid") {
      config.useAdaptiveGrid = (value.toInt() == 1);
    } else if(key == "potentialType") {
      int typeCode = value.toInt();
      if(typeCode == 0) { current.type = "WoodsSaxon"; currentHasType = true; }
      else if(typeCode == 1) { current.type = "Gaussian"; currentHasType = true; }
      else { current.enabled = false; currentHasType = true; }
      if(!currentPair) defaultHasType = currentHasType;
      sawAny = true;
    } else if(key == "V0") current.V0 = value.toDouble();
    else if(key == "R") current.R = value.toDouble();
    else if(key == "a") current.a = value.toDouble();
    else if(key == "r0") current.r0 = value.toDouble();
  }
  flush();

  // A pair may be on while the default is off, so the master switch has to
  // follow the pairs as well as its own flag.
  if(manager.isAnyEnabled()) config.useHybridMethod = true;

  currentPairKey_ = 0;
  refreshPairCombo();
  return sawAny;
}

bool NuclearPotentialTab::writePotentialSettings(QTextStream& outStream) {
  NuclearPotentialManager& manager = NuclearPotentialManager::instance();

  // The default section first -- AZURESetup has already written this block's
  // useHybridPotential= and useAdaptiveGrid= lines -- then one section per pair
  // that carries a setting of its own.  A project where nobody touched a
  // single pair writes exactly the block it always did.
  auto writeShape = [&outStream](const NuclearPotentialSetting& s) {
    outStream << "potentialType=" << (s.type == "Gaussian" ? 1 : 0) << "\n";
    outStream << "V0=" << s.V0 << "\n";
    if(s.type == "Gaussian") {
      outStream << "r0=" << s.r0 << "\n";
    } else {
      outStream << "R=" << s.R << "\n";
      outStream << "a=" << s.a << "\n";
    }
  };

  writeShape(manager.getDefaultSetting());

  std::vector<int> pairs = manager.configuredPairs();
  for(size_t i = 0; i < pairs.size(); i++) {
    NuclearPotentialSetting s = manager.getSetting(pairs[i]);
    outStream << "pair=" << pairs[i] << "\n";
    outStream << "useHybridPotential=" << (s.enabled ? 1 : 0) << "\n";
    writeShape(s);
  }

  return true;
}
