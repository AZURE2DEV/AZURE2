#include <QtGui>
#include "AddTargetIntDialog.h"

AddTargetIntDialog::AddTargetIntDialog(QWidget *parent) : QDialog(parent) {
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
  numParametersSpin = new QSpinBox;
  numParametersSpin -> setMinimum(0);
  numParametersSpin -> setMaximum(10);
  numParametersSpin -> setSingleStep(1);
  numParametersSpin -> setValue(0);
  connect(numParametersSpin,SIGNAL(valueChanged(int)),this,SLOT(parameterSpinChanged(int)));
  parametersTable = new QTableWidget(this);
  parametersTable->setColumnCount(2);
  parametersTable->setRowCount(0);
  parametersTable->verticalHeader()->hide();
  parametersTable->verticalHeader()->setHighlightSections(false);
  parametersTable->horizontalHeader()->setHighlightSections(false);
  parametersTable->horizontalHeader()->setResizeMode(0,QHeaderView::Stretch);
  parametersTable->horizontalHeader()->setResizeMode(1,QHeaderView::Stretch);
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
  qCoefficientTable->horizontalHeader()->setResizeMode(0,QHeaderView::Stretch);
  qCoefficientTable->horizontalHeader()->setResizeMode(1,QHeaderView::Stretch);
  qCoefficientTable->setShowGrid(false);
  connect(qCoefficientTable,SIGNAL(cellChanged(int,int)),this,SLOT(qCoefficientChanged(int,int)));
  QStringList qlabelList;
  qlabelList.append(QString(tr("Coefficent")));
  qlabelList.append(QString(tr("Value")));
  qCoefficientTable->setHorizontalHeaderLabels(qlabelList);


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

  stoppingPowerBox = new QGroupBox(tr("Effective Stopping Cross Section [MeV cm^2/atoms]"));
  QVBoxLayout *stoppingPowerLayout = new QVBoxLayout;
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

  QHBoxLayout *buttonBox = new QHBoxLayout;
  buttonBox->addWidget(cancelButton);
  buttonBox->addWidget(okButton);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addLayout(topLayout);
  mainLayout->addLayout(checkBoxLayout);
  mainLayout->addWidget(stoppingPowerBox);
  mainLayout->addLayout(qCoefficientCheckBoxLayout);
  mainLayout->addWidget(qCoefficientBox);
  mainLayout->addLayout(buttonBox);

  setLayout(mainLayout);

  connect(okButton, SIGNAL(clicked()),this,SLOT(accept()));
  connect(cancelButton,SIGNAL(clicked()),this,SLOT(reject()));

  setWindowTitle(tr("Add an Experimental Effect Line"));
}

void AddTargetIntDialog::createParameterItem(int row, double value) {
  QTableWidgetItem *labelItem = new QTableWidgetItem(QString("a%1").arg(row));
  QTableWidgetItem *valueItem = new QTableWidgetItem(QString("%1").arg(value));
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
    numPointsSpin->setEnabled(true);
  } else {
    stoppingPowerBox->hide();
    densityText->setEnabled(false);
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

void AddTargetIntDialog::qCoefficientChanged(int row, int column) {
  if(column==1) {
    while(row>=tempQCoefficients.size()) {
      double tempDouble=0.0;
      tempQCoefficients.append(tempDouble);
    }
    tempQCoefficients[row]=qCoefficientTable->item(row,column)->text().toDouble();
  }
}
