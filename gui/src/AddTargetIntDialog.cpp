#include <QLineEdit>
#include <QSpinBox>
#include <QCheckBox>
#include <QTableWidget>
#include <QHBoxLayout>
#include <QLabel>
#include <QGroupBox>
#include <QHeaderView>
#include <QPushButton>
#include <QComboBox>
#include <QMessageBox>
#include <vector>
#include <cmath>
#include <iomanip>

#include "AddTargetIntDialog.h"
#include "ElementMap.h"

#ifdef USE_ERYA
// ERYA SRIM integration - include utilities
#include "../../erya/include/SRIMUtilities.h"
#endif

#include <iostream>
#include <cctype>
#include <stdexcept>

#ifdef USE_ERYA
// Global SRIM utilities instance
static SRIMUtilities* srimUtils = nullptr;

// Initialize SRIM utilities on first use
static SRIMUtilities* getSRIMUtilities() {
    if (!srimUtils) {
        srimUtils = new SRIMUtilities();
    }
    return srimUtils;
}
#endif

AddTargetIntDialog::AddTargetIntDialog(QWidget *parent) : QDialog(parent), selectedElement_(13) {
  segmentsListText = new QLineEdit;
  numPointsSpin = new QSpinBox;
  numPointsSpin->setEnabled(false);
  numPointsSpin -> setMinimum(1);
  numPointsSpin -> setMaximum(10000);
  numPointsSpin -> setSingleStep(10);
  numPointsSpin -> setValue(10);

  isConvolutionCheck = new QCheckBox(tr("Include Gaussian Convolution"));
  isConvolutionCheck->setChecked(false);
  connect(isConvolutionCheck,SIGNAL(toggled(bool)),this,SLOT(convolutionCheckChanged(bool)));
  sigmaText = new QLineEdit;
  sigmaText->setEnabled(false);
  
  isTargetIntegrationCheck = new QCheckBox(tr("Include Target Integration"));
  isTargetIntegrationCheck->setChecked(false);
  connect(isTargetIntegrationCheck,SIGNAL(toggled(bool)),this,SLOT(targetIntCheckChanged(bool)));
  densityText = new QLineEdit;
  densityText->setEnabled(false);
  stoppingPowerEqText = new QLineEdit;
  
  // Energy and deltaE calculation widgets
  energyText = new QLineEdit;
  energyText->setPlaceholderText("1000.0");
  energyText->setEnabled(false);
  deltaEText = new QLineEdit;
  deltaEText->setReadOnly(true);
  deltaEText->setEnabled(false);
  calculateDeltaEButton = new QPushButton(tr("Calculate ΔE"));
  calculateDeltaEButton->setEnabled(false);
  connect(calculateDeltaEButton, SIGNAL(clicked()), this, SLOT(calculateDeltaE()));
  
  // Element selection for ERYA integration
  elementComboBox = new QComboBox;
  populateElementComboBox();
  connect(elementComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(elementSelectionChanged(int)));
  
  // Compound formula input (e.g., CH4, SiO2)
  compoundText = new QLineEdit;
  compoundText->setPlaceholderText("Enter compound (e.g., CH4, SiO2)");
  compoundText->setToolTip("Enter chemical formula for compounds. Leave empty for single elements.");
  
  fetchStoppingPowerButton = new QPushButton(tr("Fetch from ERYA"));
  connect(fetchStoppingPowerButton, SIGNAL(clicked()), this, SLOT(fetchStoppingPowerParameters()));
  
  isStraggling = new QCheckBox(tr("Include Energy Straggling"));
  isStraggling->setChecked(false);
  numParametersSpin = new QSpinBox;
  numParametersSpin -> setMinimum(0);
  numParametersSpin -> setMaximum(20);
  numParametersSpin -> setSingleStep(1);
  numParametersSpin -> setValue(0);
  connect(numParametersSpin,SIGNAL(valueChanged(int)),this,SLOT(parameterSpinChanged(int)));
  parametersTable = new QTableWidget(this);
  parametersTable->setColumnCount(2);
  parametersTable->setRowCount(0);
  parametersTable->verticalHeader()->hide();
  parametersTable->verticalHeader()->setHighlightSections(false);
  parametersTable->horizontalHeader()->setHighlightSections(false);
  parametersTable->horizontalHeader()->setSectionResizeMode(0,QHeaderView::Stretch);
  parametersTable->horizontalHeader()->setSectionResizeMode(1,QHeaderView::Stretch);
  parametersTable->setShowGrid(false);
  connect(parametersTable,SIGNAL(cellChanged(int,int)),this,SLOT(parameterChanged(int,int)));
  QStringList labelList;
  labelList.append(QString(tr("Parameter")));
  labelList.append(QString(tr("Value")));
  parametersTable->setHorizontalHeaderLabels(labelList);
  
  isQCoefficientCheck = new QCheckBox(tr("Include Attenuation Coefficients"));
  isQCoefficientCheck->setChecked(false);
  connect(isQCoefficientCheck,SIGNAL(toggled(bool)),this,SLOT(qCoefficientCheckChanged(bool)));

  numQCoefficientSpin= new QSpinBox;
  numQCoefficientSpin->setMinimum(0);
  numQCoefficientSpin->setMaximum(6);
  numQCoefficientSpin->setSingleStep(1);
  numQCoefficientSpin->setValue(0);
  numQCoefficientSpin->setEnabled(false);
  connect(numQCoefficientSpin,SIGNAL(valueChanged(int)),this,SLOT(qCoefficientSpinChanged(int)));
  
  qCoefficientTable = new QTableWidget(this);
  qCoefficientTable->setColumnCount(2);
  qCoefficientTable->setRowCount(0);
  qCoefficientTable->verticalHeader()->hide();
  qCoefficientTable->verticalHeader()->setHighlightSections(false);
  qCoefficientTable->horizontalHeader()->setHighlightSections(false);
  qCoefficientTable->horizontalHeader()->setSectionResizeMode(0,QHeaderView::Stretch);
  qCoefficientTable->horizontalHeader()->setSectionResizeMode(1,QHeaderView::Stretch);
  qCoefficientTable->setShowGrid(false);
  connect(qCoefficientTable,SIGNAL(cellChanged(int,int)),this,SLOT(qCoefficientChanged(int,int)));

  QStringList qlabelList;
  qlabelList.append(QString(tr("Coefficent")));
  qlabelList.append(QString(tr("Value")));
  qCoefficientTable->setHorizontalHeaderLabels(qlabelList);

  isConvolutionDependentCheck = new QCheckBox(tr("Include Energy Dependent Gaussian Convolution"));
  isConvolutionDependentCheck->setChecked(false);
  connect(isConvolutionDependentCheck,SIGNAL(toggled(bool)),this,SLOT(convCoefficientCheckChanged(bool)));

  convolutionEqText = new QLineEdit;
  numConvCoefficientSpin= new QSpinBox;
  numConvCoefficientSpin->setMinimum(0);
  numConvCoefficientSpin->setMaximum(6);
  numConvCoefficientSpin->setSingleStep(1);
  numConvCoefficientSpin->setValue(0);
  numConvCoefficientSpin->setEnabled(false);
  connect(numConvCoefficientSpin,SIGNAL(valueChanged(int)),this,SLOT(convCoefficientSpinChanged(int)));

  convCoefficientTable = new QTableWidget(this);
  convCoefficientTable->setColumnCount(2);
  convCoefficientTable->setRowCount(0);
  convCoefficientTable->verticalHeader()->hide();
  convCoefficientTable->verticalHeader()->setHighlightSections(false);
  convCoefficientTable->horizontalHeader()->setHighlightSections(false);
  convCoefficientTable->horizontalHeader()->setSectionResizeMode(0,QHeaderView::Stretch);
  convCoefficientTable->horizontalHeader()->setSectionResizeMode(1,QHeaderView::Stretch);
  convCoefficientTable->setShowGrid(false);
  connect(convCoefficientTable,SIGNAL(cellChanged(int,int)),this,SLOT(convCoefficientChanged(int,int)));

  QStringList convlabelList;
  convlabelList.append(QString(tr("Coefficent")));
  convlabelList.append(QString(tr("Value")));
  convCoefficientTable->setHorizontalHeaderLabels(convlabelList);

  cancelButton = new QPushButton(tr("Cancel"));
  okButton = new QPushButton(tr("Accept"));
  okButton->setDefault(true);

  QHBoxLayout *segListLayout = new QHBoxLayout;
  segListLayout->addWidget(new QLabel(tr("Segments List:")));
  segListLayout->addWidget(segmentsListText);

  QHBoxLayout *numPointsLayout = new QHBoxLayout;
  numPointsLayout->addWidget(new QLabel(tr("Number of Integration Points:")));
  numPointsLayout->addWidget(numPointsSpin);
  
  QHBoxLayout *topLayout = new QHBoxLayout;
  topLayout->addLayout(segListLayout);
  topLayout->addLayout(numPointsLayout);

  QGridLayout *checkBoxLayout = new QGridLayout;
  checkBoxLayout->addWidget(isConvolutionCheck,0,0);
  checkBoxLayout->addWidget(new QLabel(tr("Sigma [MeV]:")),0,1,Qt::AlignRight);
  checkBoxLayout->addWidget(sigmaText,0,2);

  checkBoxLayout->addWidget(isTargetIntegrationCheck,1,0);
  checkBoxLayout->addWidget(new QLabel(tr("Active Density [atoms/cm^2]:")),1,1,Qt::AlignRight);
  checkBoxLayout->addWidget(densityText,1,2);
  
  // Add energy and deltaE calculation on the same row
  QHBoxLayout *energyDeltaELayout = new QHBoxLayout;
  energyDeltaELayout->addWidget(new QLabel(tr("Energy [keV]:")));
  energyDeltaELayout->addWidget(energyText);
  energyDeltaELayout->addWidget(new QLabel(tr("ΔE [keV]:")));
  energyDeltaELayout->addWidget(deltaEText);
  energyDeltaELayout->addWidget(calculateDeltaEButton);
  
  QWidget *energyDeltaEWidget = new QWidget;
  energyDeltaEWidget->setLayout(energyDeltaELayout);
  checkBoxLayout->addWidget(energyDeltaEWidget,2,1,1,2); // Span 2 columns

  stoppingPowerBox = new QGroupBox(tr("Effective Stopping Cross Section [MeV cm^2/atoms]"));
  QVBoxLayout *stoppingPowerLayout = new QVBoxLayout;
  
  // Element selection layout
  QHBoxLayout *elementLayout = new QHBoxLayout;
  elementLayout->addWidget(new QLabel(tr("Target Element:")));
  elementLayout->addWidget(elementComboBox);
  stoppingPowerLayout->addLayout(elementLayout);
  
  // Compound input layout
  QHBoxLayout *compoundLayout = new QHBoxLayout;
  compoundLayout->addWidget(new QLabel(tr("Or Compound Formula:")));
  compoundLayout->addWidget(compoundText);
  compoundLayout->addWidget(fetchStoppingPowerButton);
  stoppingPowerLayout->addLayout(compoundLayout);
  
  // Straggling checkbox
  stoppingPowerLayout->addWidget(isStraggling);
  
  QHBoxLayout *stoppingPowerTopLayout = new QHBoxLayout;
  QHBoxLayout *equationLayout = new QHBoxLayout;
  equationLayout->addWidget(new QLabel(tr("y=")));
  equationLayout->addWidget(stoppingPowerEqText);
  QHBoxLayout *numParamLayout = new QHBoxLayout;
  numParamLayout->addWidget(new QLabel(tr("Number of Parameters:")));
  numParamLayout->addWidget(numParametersSpin);
  stoppingPowerTopLayout->addLayout(equationLayout);
  stoppingPowerTopLayout->addLayout(numParamLayout);
  stoppingPowerLayout->addLayout(stoppingPowerTopLayout);
  stoppingPowerLayout->addWidget(parametersTable);
  stoppingPowerBox->setLayout(stoppingPowerLayout);
  stoppingPowerBox->hide();

  QGridLayout *qCoefficientCheckBoxLayout = new QGridLayout;
  qCoefficientCheckBoxLayout->addWidget(isQCoefficientCheck,0,0);
  qCoefficientCheckBoxLayout->addItem(new QSpacerItem(1,20),0,1);
  qCoefficientCheckBoxLayout->addWidget(new QLabel(tr("Number of Coefficients:")),0,2,Qt::AlignRight);
  qCoefficientCheckBoxLayout->addWidget(numQCoefficientSpin,0,3);
  qCoefficientCheckBoxLayout->setColumnStretch(1,1);

  qCoefficientBox = new QGroupBox(tr("Attenuation Coefficients"));
  QVBoxLayout *qCoefficientLayout = new QVBoxLayout;
  qCoefficientLayout->addWidget(qCoefficientTable);
  qCoefficientBox->setLayout(qCoefficientLayout);
  qCoefficientBox->hide();

  QGridLayout *energyConvolutionCheckBoxLayout = new QGridLayout;
  energyConvolutionCheckBoxLayout->addWidget(isConvolutionDependentCheck,0,0);
  energyConvolutionCheckBoxLayout->addItem(new QSpacerItem(1,20),0,1);
  energyConvolutionCheckBoxLayout->addWidget(new QLabel(tr("Number of Coefficients:")),0,2,Qt::AlignRight);
  energyConvolutionCheckBoxLayout->addWidget(numConvCoefficientSpin,0,3);
  energyConvolutionCheckBoxLayout->setColumnStretch(1,1);

  convCoefficientBox = new QGroupBox(tr("Convolution Coefficients"));
  QVBoxLayout *convCoefficientLayout = new QVBoxLayout;
  QHBoxLayout *equationConvLayout = new QHBoxLayout;
  equationConvLayout->addWidget(new QLabel(tr("y=")));
  equationConvLayout->addWidget(convolutionEqText);
  convCoefficientLayout->addLayout(equationConvLayout);
  convCoefficientLayout->addWidget(convCoefficientTable);
  convCoefficientBox->setLayout(convCoefficientLayout);
  convCoefficientBox->hide();

  QHBoxLayout *buttonBox = new QHBoxLayout;
  buttonBox->addWidget(cancelButton);
  buttonBox->addWidget(okButton);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addLayout(topLayout);
  mainLayout->addLayout(checkBoxLayout);
  mainLayout->addWidget(stoppingPowerBox);
  mainLayout->addLayout(qCoefficientCheckBoxLayout);
  mainLayout->addWidget(qCoefficientBox);
  mainLayout->addLayout(energyConvolutionCheckBoxLayout);
  mainLayout->addWidget(convCoefficientBox);
  mainLayout->addLayout(buttonBox);

  setLayout(mainLayout);

  connect(okButton, SIGNAL(clicked()),this,SLOT(accept()));
  connect(cancelButton,SIGNAL(clicked()),this,SLOT(reject()));

  setWindowTitle(tr("Add an Experimental Effect Line"));
}

void AddTargetIntDialog::createParameterItem(int row, double value) {
  QTableWidgetItem *labelItem = new QTableWidgetItem(QString("a%1").arg(row));
  
  // Use adaptive formatting for SRIM parameters
  QString valueStr;
  if (std::abs(value) < 1e-6 || std::abs(value) > 1e6) {
    // Use scientific notation for very small or very large numbers
    valueStr = QString::number(value, 'e', 8);
  } else if (std::abs(value) < 0.001) {
    // Use scientific notation for small numbers
    valueStr = QString::number(value, 'e', 6);
  } else if (std::abs(value) < 1.0) {
    // Use fixed notation with more decimal places for numbers < 1
    valueStr = QString::number(value, 'f', 8);
  } else {
    // Use fixed notation for normal numbers
    valueStr = QString::number(value, 'f', 6);
  }
  
  QTableWidgetItem *valueItem = new QTableWidgetItem(valueStr);
  labelItem->setTextAlignment(Qt::AlignCenter);
  labelItem->setFlags(Qt::ItemIsEditable);
  valueItem->setTextAlignment(Qt::AlignCenter);
  parametersTable->setItem(row,0,labelItem);
  parametersTable->setItem(row,1,valueItem);
  parametersTable->resizeRowsToContents();
}

void AddTargetIntDialog::createQCoefficientItem(int row, double value) {
  QTableWidgetItem *labelItem = new QTableWidgetItem(QString("q%1").arg(row));
  QTableWidgetItem *valueItem = new QTableWidgetItem(QString("%1").arg(value));
  labelItem->setTextAlignment(Qt::AlignCenter);
  labelItem->setFlags(Qt::ItemIsEditable);
  valueItem->setTextAlignment(Qt::AlignCenter);
  qCoefficientTable->setItem(row,0,labelItem);
  qCoefficientTable->setItem(row,1,valueItem);
  qCoefficientTable->resizeRowsToContents();
}

void AddTargetIntDialog::createConvCoefficientItem(int row, double value) {
  QTableWidgetItem *labelItem = new QTableWidgetItem(QString("a%1").arg(row));
  QTableWidgetItem *valueItem = new QTableWidgetItem(QString("%1").arg(value));
  labelItem->setTextAlignment(Qt::AlignCenter);
  labelItem->setFlags(Qt::ItemIsEditable);
  valueItem->setTextAlignment(Qt::AlignCenter);
  convCoefficientTable->setItem(row,0,labelItem);
  convCoefficientTable->setItem(row,1,valueItem);
  convCoefficientTable->resizeRowsToContents();
}

void AddTargetIntDialog::convolutionCheckChanged(bool checked) {
  if(checked) {
    sigmaText->setEnabled(true);
    numPointsSpin->setEnabled(true);
  } else {
    sigmaText->setEnabled(false);
    if(!isTargetIntegrationCheck->isChecked()) numPointsSpin->setEnabled(false);
  }
}

void AddTargetIntDialog::targetIntCheckChanged(bool checked) {
  if(checked) {
    stoppingPowerBox->show();
    densityText->setEnabled(true);
    energyText->setEnabled(true);
    calculateDeltaEButton->setEnabled(true);
    numPointsSpin->setEnabled(true);
  } else {
    stoppingPowerBox->hide();
    densityText->setEnabled(false);
    energyText->setEnabled(false);
    deltaEText->clear();
    calculateDeltaEButton->setEnabled(false);
    if(!isConvolutionCheck->isChecked()) numPointsSpin->setEnabled(false);
    this->adjustSize();
  }
}

void AddTargetIntDialog::parameterSpinChanged(int newNumber) {
  parametersTable->clearContents();
  parametersTable->setRowCount(newNumber);
  for(int i=0;i<newNumber;i++) {
    if(i<tempParameters.size())
      createParameterItem(i,tempParameters.at(i));
    else 
      createParameterItem(i);
  }
}

void AddTargetIntDialog::parameterChanged(int row, int column) {
  if(column==1) {
    while(row>=tempParameters.size()) {
      double tempDouble=0.0;
      tempParameters.append(tempDouble);
    }
    tempParameters[row]=parametersTable->item(row,column)->text().toDouble();
  }
}

void AddTargetIntDialog::qCoefficientCheckChanged(bool checked) {
  if(checked) {
    qCoefficientBox->show();
    numQCoefficientSpin->setEnabled(true);
  } else {
    qCoefficientBox->hide();
    numQCoefficientSpin->setEnabled(false);
    this->adjustSize();
  }
}

void AddTargetIntDialog::convCoefficientCheckChanged(bool checked) {
  if(checked) {
    convCoefficientBox->show();
    numConvCoefficientSpin->setEnabled(true);
    numPointsSpin->setEnabled(true);
  } else {
    convCoefficientBox->hide();
    numConvCoefficientSpin->setEnabled(false);
    numPointsSpin->setEnabled(false);
    this->adjustSize();
  }
}

void AddTargetIntDialog::qCoefficientSpinChanged(int newNumber){
  qCoefficientTable->clearContents();
  qCoefficientTable->setRowCount(newNumber);
  for(int i=0;i<newNumber;i++) {
    if(i<tempQCoefficients.size())
      createQCoefficientItem(i,tempQCoefficients.at(i));
    else 
      createQCoefficientItem(i);
  }
}

void AddTargetIntDialog::convCoefficientSpinChanged(int newNumber){
  convCoefficientTable->clearContents();
  convCoefficientTable->setRowCount(newNumber);
  for(int i=0;i<newNumber;i++) {
    if(i<tempConvCoefficients.size())
      createConvCoefficientItem(i,tempConvCoefficients.at(i));
    else 
      createConvCoefficientItem(i);
  }
}

void AddTargetIntDialog::qCoefficientChanged(int row, int column) {
  if(column==1) {
    while(row>=tempQCoefficients.size()) {
      double tempDouble=0.0;
      tempQCoefficients.append(tempDouble);
    }
    tempQCoefficients[row]=qCoefficientTable->item(row,column)->text().toDouble();
  }
}

void AddTargetIntDialog::convCoefficientChanged(int row, int column) {
  if(column==1) {
    while(row>=tempConvCoefficients.size()) {
      double tempDouble=0.0;
      tempConvCoefficients.append(tempDouble);
    }
    tempConvCoefficients[row]=convCoefficientTable->item(row,column)->text().toDouble();
  }
}

void AddTargetIntDialog::populateElementComboBox() {
  // Populate combo box with elements from ElementMap, organized by common usage
  elementComboBox->clear();
  
  // Most common target elements in nuclear physics, grouped logically
  struct ElementGroup {
    QString groupName;
    std::vector<int> elements;
  };
  
  std::vector<ElementGroup> elementGroups = {
    {"Light Elements", {1, 2, 3, 4}},                    // H, He, Li, Be
    {"Common Targets", {6, 7, 8, 9, 10}},               // C, N, O, F, Ne  
    {"Metals (Light)", {11, 12, 13, 14, 15, 16}},       // Na, Mg, Al, Si, P, S
    {"Metals (Medium)", {17, 18, 19, 20, 21, 22, 23, 24, 25, 26}}, // Cl-Fe
    {"Heavy Metals", {27, 28, 29, 30, 47, 48, 49, 50}}, // Co, Ni, Cu, Zn, Ag, Cd, In, Sn
    {"Actinides", {90, 92}}                             // Th, U
  };  
  // Add elements organized by groups for better user experience
  for (const auto& group : elementGroups) {
    // Add separator for visual grouping (except for first group)
    if (&group != &elementGroups[0]) {
      elementComboBox->insertSeparator(elementComboBox->count());
    }
    
    for (int Z : group.elements) {
      if (elementMap.find(Z) != elementMap.end()) {
        QString elementText = QString("%1 (%2) - %3").arg(elementMap.at(Z)).arg(Z).arg(group.groupName);
        elementComboBox->addItem(elementText, Z);
      }
    }
  }
  
  // Set default to Aluminum (Z=13) as it's very common in nuclear physics
  int alIndex = elementComboBox->findData(13);
  if (alIndex >= 0) {
    elementComboBox->setCurrentIndex(alIndex);
  }
}

void AddTargetIntDialog::elementSelectionChanged(int index) {
  if (index >= 0) {
    selectedElement_ = elementComboBox->itemData(index).toInt();
  }
}

void AddTargetIntDialog::fetchStoppingPowerParameters() {
#ifdef USE_ERYA
  // Check if compound formula is provided
  QString compoundFormula = compoundText->text().trimmed();
  
  if (!compoundFormula.isEmpty()) {
    // Handle compound
    updateStoppingPowerFromCompound(compoundFormula.toStdString());
  } else {
    // Handle single element
    updateStoppingPowerFromElement(selectedElement_);
  }
#else
  QMessageBox::warning(this, tr("ERYA Not Available"), 
                      tr("ERYA SRIM utilities are not enabled in this build.\n\n"
                         "To use ERYA stopping power data, please rebuild with USE_ERYA=ON."));
#endif
}

void AddTargetIntDialog::updateStoppingPowerFromElement(int element) {
#ifdef USE_ERYA
  SRIMUtilities* srimUtils = getSRIMUtilities();
  
  if (!srimUtils->isDataLoaded()) {
    QMessageBox::warning(this, tr("ERYA SRIM Error"), 
                        tr("Failed to load SRIM data.\n\n"
                           "Please check that the SRIM2013.xml file is available."));
    return;
  }
  
  SRIMElementData srimData = srimUtils->readSRIMDataForElement(element);
  
  if (!srimData.isValid) {
    QMessageBox::warning(this, tr("ERYA SRIM Error"), 
                        tr("Failed to load SRIM data for element Z=%1.\n\n"
                           "Please check that the element is available in the database.").arg(element));
    return;
  }
  
  QString elementName = QString::fromStdString(srimData.elementName);
  
  // Show detailed information from actual SRIM data
  QMessageBox::information(this, tr("ERYA SRIM Integration"), 
                          tr("Loaded SRIM data for %1 (Z=%2):\n\n"
                             "Atomic mass: %3 u\n"
                             "Atomic density: %4 g/cm³\n"
                             "Bloch parameter: %5 eV\n\n"
                             "Ziegler parameters loaded from SRIM 2013 database.")
                             .arg(elementName).arg(element)
                             .arg(srimData.atomicMass, 0, 'f', 4)
                             .arg(srimData.atomicDensity, 0, 'e', 3)
                             .arg(srimData.blochParameter, 0, 'f', 1));
  
  // Generate AZURE2-compatible stopping power equation
  std::vector<double> parameters;
  std::string equation = srimUtils->generateAZUREEquation(element, parameters);
  
  if (equation.empty()) {
    QMessageBox::warning(this, tr("ERYA SRIM Error"), 
                        tr("Failed to generate stopping power equation for element Z=%1.").arg(element));
    return;
  }
  
  updateStoppingPowerGUI(QString::fromStdString(equation), parameters, elementName);
  
  // Auto-enable straggling for elements heavier than boron
  if(element >= 6 && !isStraggling->isChecked()) {
    isStraggling->setChecked(true);
  }
#endif
}

void AddTargetIntDialog::updateStoppingPowerFromCompound(const std::string& formula) {
#ifdef USE_ERYA
  SRIMUtilities* srimUtils = getSRIMUtilities();
  
  if (!srimUtils->isDataLoaded()) {
    QMessageBox::warning(this, tr("ERYA SRIM Error"), 
                        tr("Failed to load SRIM data.\n\n"
                           "Please check that the SRIM2013.xml file is available."));
    return;
  }
  
  // Parse compound formula
  std::vector<CompoundElement> elements = srimUtils->parseCompoundFormula(formula);
  
  if (elements.empty()) {
    QMessageBox::warning(this, tr("Compound Formula Error"), 
                        tr("Failed to parse compound formula '%1'.\n\n"
                           "Please use standard chemical notation (e.g., CH4, SiO2).").arg(QString::fromStdString(formula)));
    return;
  }
  
  // Show compound information
  QString compoundInfo = tr("Parsed compound '%1':\n\n").arg(QString::fromStdString(formula));
  for (const auto& element : elements) {
    SRIMElementData data = srimUtils->readSRIMDataForElement(element.elementNumber);
    if (data.isValid) {
      compoundInfo += tr("- %1 (Z=%2): stoichiometry = %3\n")
                      .arg(QString::fromStdString(data.elementName))
                      .arg(element.elementNumber)
                      .arg(element.stoichiometry);
    }
  }
  compoundInfo += tr("\nStopping power will be calculated as weighted average.");
  
  QMessageBox::information(this, tr("ERYA Compound Integration"), compoundInfo);
  
  // Generate AZURE2-compatible stopping power equation for compound
  std::vector<double> parameters;
  std::string equation = srimUtils->generateCompoundAZUREEquation(elements, parameters);
  
  if (equation.empty()) {
    QMessageBox::warning(this, tr("ERYA SRIM Error"), 
                        tr("Failed to generate stopping power equation for compound '%1'.").arg(QString::fromStdString(formula)));
    return;
  }
  
  updateStoppingPowerGUI(QString::fromStdString(equation), parameters, QString::fromStdString(formula));
  
  // Auto-enable straggling for compounds
  if (!isStraggling->isChecked()) {
    isStraggling->setChecked(true);
  }
#endif
}

void AddTargetIntDialog::updateStoppingPowerGUI(const QString& equation, const std::vector<double>& parameters, const QString& materialName) {
  // Update GUI with the SRIM-based stopping power equation
  stoppingPowerEqText->setText(equation);
  
  // Set number of parameters
  int numParams = parameters.size();
  numParametersSpin->setValue(numParams);
  
  // Update parameter table with actual SRIM values
  parametersTable->clearContents();
  parametersTable->setRowCount(numParams);
  
  for(int i = 0; i < numParams; i++) {
    createParameterItem(i, parameters[i]);
  }
  
  // Update temporary parameters
  tempParameters.clear();
  for(int i = 0; i < numParams; i++) {
    tempParameters.append(parameters[i]);
  }
  
  // Auto-enable target integration checkbox since we have stopping power data
  if(!isTargetIntegrationCheck->isChecked()) {
    isTargetIntegrationCheck->setChecked(true);
    targetIntCheckChanged(true);
  }
}

void AddTargetIntDialog::calculateDeltaE() {
  // Calculate energy loss (deltaE) using stopping power parameters
  
  // Get the energy from input field
  bool ok;
  double energy_keV = energyText->text().toDouble(&ok);
  if (!ok || energy_keV <= 0) {
    QMessageBox::warning(this, tr("Invalid Energy"), 
                        tr("Please enter a valid energy value in keV (e.g., 1000.0)"));
    return;
  }
  
  // Get the target density
  double density_atoms_cm2 = densityText->text().toDouble(&ok);
  if (!ok || density_atoms_cm2 <= 0) {
    QMessageBox::warning(this, tr("Invalid Density"), 
                        tr("Please enter a valid target density in atoms/cm²"));
    return;
  }
  
  // Get parameters from the table
  if (tempParameters.size() < 1) {
    QMessageBox::warning(this, tr("Insufficient Parameters"), 
                        tr("Please fetch stopping power parameters first using 'Fetch from ERYA'.\n\nCurrent parameters: %1").arg(tempParameters.size()));
    return;
  }
  
  try {
#ifdef USE_ERYA
    // Use ERYA utilities for calculation
    SRIMUtilities* srimUtils = getSRIMUtilities();
    
    // Check if compound formula is provided
    QString compoundFormula = compoundText->text().trimmed();
    
    double stoppingPower = 0.0;
    
    if (!compoundFormula.isEmpty()) {
      // Calculate stopping power for compound
      std::vector<CompoundElement> elements = srimUtils->parseCompoundFormula(compoundFormula.toStdString());
      if (!elements.empty()) {
        stoppingPower = srimUtils->calculateCompoundStoppingPower(energy_keV, elements);
      }
    } else {
      // Calculate stopping power for single element
      SRIMElementData data = srimUtils->readSRIMDataForElement(selectedElement_);
      if (data.isValid) {
        stoppingPower = srimUtils->calculateZieglerStoppingPower(energy_keV, data);
      }
    }
    
    if (stoppingPower <= 0.0) {
      QMessageBox::warning(this, tr("Calculation Error"), 
                          tr("Failed to calculate stopping power. Please check your parameters."));
      return;
    }
    
    // Calculate energy loss: ΔE = StoppingPower × Density
    // SRIM stopping power is now converted to MeV⋅cm²/atom in SRIMUtilities
    double deltaE_keV = stoppingPower * density_atoms_cm2 * 1000.0;  // Convert MeV to keV
    
    // Display the result
    deltaEText->setText(QString::number(deltaE_keV, 'f', 3));

#else
    // Fallback calculation without ERYA utilities
    // Use basic stopping power calculation based on parameters
    if (tempParameters.size() >= 8) {
      // Use high-energy compound formula for all energies to match equation
      double stoppingLow = tempParameters[0] * std::pow(energy_keV, tempParameters[1]) +
                          tempParameters[2] * std::pow(energy_keV, tempParameters[3]);
      double stoppingHigh = (tempParameters[4] / std::pow(energy_keV, tempParameters[5])) *
                            std::log((tempParameters[6] / energy_keV) + (tempParameters[7] * energy_keV));
      double stoppingPower = (stoppingLow * stoppingHigh) / (stoppingHigh + stoppingLow);
      
      // Apply unit conversion and calculate energy loss
      double deltaE_keV = (stoppingPower / 1e+21) * density_atoms_cm2 * 1000.0;  // Convert to keV
      deltaEText->setText(QString::number(deltaE_keV, 'f', 3));
    } else {
      QMessageBox::warning(this, tr("Insufficient Parameters"), 
                          tr("Need at least 8 parameters for stopping power calculation."));
      return;
    }
#endif
                            
  } catch (const std::exception& e) {
    QMessageBox::critical(this, tr("Calculation Error"), 
                         tr("Error during ΔE calculation: %1").arg(e.what()));
  }
}

