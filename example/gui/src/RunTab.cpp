#include "RunTab.h"
#include "ChooseFileButton.h"
#include "FilteredTextEdit.h"
#include "InfoDialog.h"
#include <QtGui>

RunTab::RunTab(QWidget* parent) : QWidget(parent) { 
  calcType = new QComboBox;
  calcType->addItem(tr("Calculate Segments From Data"));
  calcType->addItem(tr("Fit Segments From Data"));
  calcType->addItem(tr("Calculate Segments Without Data"));
  calcType->addItem(tr("Perform MINOS Error Analysis"));
  calcType->addItem(tr("Calculate Reaction Rate"));
  connect(calcType,SIGNAL(currentIndexChanged(int)),this,SLOT(calculationTypeChanged(int)));

  chiVarianceText = new QLineEdit;
  chiVarianceText->setText("1.0");
  chiVarianceText->setEnabled(false);
  chiVarianceText->setMinimumWidth(50);
  chiVarianceText->setMaximumWidth(50);
 
  calcButton = new QPushButton(tr("Save and &Run"));
  stopAZUREButton = new QPushButton(tr("Stop AZURE2"));
  stopAZUREButton->setEnabled(false);

  QGridLayout* calcLayout = new QGridLayout;
  calcLayout->addWidget(new QLabel(tr("Calculation Type:")),0,0,Qt::AlignRight);
  calcLayout->addWidget(calcType,0,1);
  calcLayout->setColumnStretch(1,0);
  calcLayout->addItem(new QSpacerItem(28,28),0,2);
  calcLayout->setColumnStretch(2,1);
  calcLayout->addWidget(new QLabel(tr("Chi-Squared Variance:")),0,3,Qt::AlignRight);
  calcLayout->addWidget(chiVarianceText,0,4);
  calcLayout->setColumnStretch(4,0);
  calcLayout->addWidget(calcButton,0,5);
  calcLayout->setColumnStretch(5,0);
  calcLayout->addWidget(stopAZUREButton,0,6);
  calcLayout->setColumnStretch(6,0);

  paramFileText = new QLineEdit;
  paramFileText->setEnabled(false);
  newParamFileButton = new QRadioButton(tr("Create New Parameters File"));
  newParamFileButton->setMinimumWidth(200);
  newParamFileButton->setMaximumWidth(200);
  oldParamFileButton = new QRadioButton(tr("Use: "));
  connect(oldParamFileButton,SIGNAL(toggled(bool)),this,SLOT(paramFileButtonChanged(bool)));
  newParamFileButton->setChecked(true);
  
  QGroupBox* paramFileGroup = new QGroupBox(tr("Parameters File"));
  QHBoxLayout* paramFileLayout = new QHBoxLayout;
  paramFileLayout->setContentsMargins(5,5,5,5);
  paramFileLayout->addWidget(newParamFileButton);
  paramFileLayout->addWidget(oldParamFileButton);
  paramFileLayout->addWidget(paramFileText);
  paramFileChoose=new ChooseFileButton(tr("Choose..."));
  paramFileChoose->setEnabled(false);
  paramFileChoose->setLineEdit(paramFileText);
  connect(paramFileChoose,SIGNAL(clicked(QLineEdit*)),this,SLOT(setChooseFile(QLineEdit*)));
  paramFileLayout->addWidget(paramFileChoose);
  paramFileGroup->setLayout(paramFileLayout);

  integralsFileText = new QLineEdit;
  integralsFileText->setEnabled(false);
  newIntegralsFileButton = new QRadioButton(tr("Create New Integrals File"));
  newIntegralsFileButton->setMinimumWidth(200);
  newIntegralsFileButton->setMaximumWidth(200);
  oldIntegralsFileButton = new QRadioButton(tr("Use: "));
  connect(oldIntegralsFileButton,SIGNAL(toggled(bool)),this,SLOT(integralsFileButtonChanged(bool)));
  newIntegralsFileButton->setChecked(true);
 
  integralsFileGroup = new QGroupBox(tr("External Capture Integrals File"));
  QHBoxLayout* integralsFileLayout = new QHBoxLayout;
  integralsFileLayout->setContentsMargins(5,5,5,5);
  integralsFileLayout->addWidget(newIntegralsFileButton);
  integralsFileLayout->addWidget(oldIntegralsFileButton);
  integralsFileLayout->addWidget(integralsFileText);
  integralsFileChoose=new ChooseFileButton(tr("Choose..."));
  integralsFileChoose->setEnabled(false);
  integralsFileChoose->setLineEdit(integralsFileText);
  connect(integralsFileChoose,SIGNAL(clicked(QLineEdit*)),this,SLOT(setChooseFile(QLineEdit*)));
  integralsFileLayout->addWidget(integralsFileChoose);
  integralsFileGroup->setLayout(integralsFileLayout);

  gridTempButton = new QRadioButton(tr("Create Temperatures"));
  rateEntranceKey = new QLineEdit;
  rateEntranceKey->setMinimumWidth(50);
  rateEntranceKey->setMaximumWidth(50);
  rateExitKey = new QLineEdit;
  rateExitKey->setMinimumWidth(50);
  rateExitKey->setMaximumWidth(50);
  minTempText = new QLineEdit;
  minTempText->setMinimumWidth(50);
  minTempText->setMaximumWidth(50);
  maxTempText = new QLineEdit;
  maxTempText->setMinimumWidth(50);
  maxTempText->setMaximumWidth(50);
  tempStepText = new QLineEdit;
  tempStepText->setMinimumWidth(50);
  tempStepText->setMaximumWidth(50);
  fileTempButton = new QRadioButton(tr("Use Temperature File: "));
  connect(fileTempButton,SIGNAL(toggled(bool)),this,SLOT(fileTempButtonChanged(bool)));
  fileTempText = new QLineEdit;
  fileTempText->setEnabled(false);
  gridTempButton->setChecked(true);

  rateParamsGroup = new QGroupBox(tr("Reaction Rate Parameters"));
  rateParamsGroup->hide();
  QGridLayout* keyLayout = new QGridLayout;
  keyLayout->addWidget(new QLabel(tr("Entrance Key:")),0,0,Qt::AlignRight);
  keyLayout->addWidget(rateEntranceKey,0,1);
  keyLayout->addWidget(new QLabel(tr("Exit Key:")),0,2,Qt::AlignRight);
  keyLayout->addWidget(rateExitKey,0,3);
  keyLayout->addItem(new QSpacerItem(28,28),0,4);
  keyLayout->setColumnStretch(4,1);
  QGridLayout* gridTempLayout = new QGridLayout;
  gridTempLayout->addWidget(gridTempButton,0,0);
  gridTempLayout->setColumnStretch(0,0);
  gridTempLayout->addItem(new QSpacerItem(28,28),0,1);
  gridTempLayout->setColumnStretch(1,1);
  gridTempLayout->addWidget(new QLabel(tr("Minimum Temperature (GK):")),0,2,Qt::AlignRight);
  gridTempLayout->addWidget(minTempText,0,3);
  gridTempLayout->setColumnStretch(3,0);
  gridTempLayout->addWidget(new QLabel(tr("Maximum Temperature (GK):")),0,4,Qt::AlignRight);
  gridTempLayout->addWidget(maxTempText,0,5);
  gridTempLayout->setColumnStretch(5,0);
  gridTempLayout->addWidget(new QLabel(tr("Temperature Step (GK):")),0,6,Qt::AlignRight);
  gridTempLayout->addWidget(tempStepText,0,7);
  gridTempLayout->setColumnStretch(7,0);
  QHBoxLayout *fileTempLayout = new QHBoxLayout;
  fileTempLayout->addWidget(fileTempButton);
  fileTempLayout->addWidget(fileTempText);
  rateParamsChoose=new ChooseFileButton(tr("Choose..."));
  rateParamsChoose->setEnabled(false);
  rateParamsChoose->setLineEdit(fileTempText);
  connect(rateParamsChoose,SIGNAL(clicked(QLineEdit*)),this,SLOT(setChooseFile(QLineEdit*)));
  fileTempLayout->addWidget(rateParamsChoose);

  QVBoxLayout* rateParamsLayout = new QVBoxLayout;
  rateParamsLayout->setContentsMargins(5,5,5,5);
  rateParamsLayout->addLayout(keyLayout);
  rateParamsLayout->addLayout(gridTempLayout);
  rateParamsLayout->addLayout(fileTempLayout);
  rateParamsGroup->setLayout(rateParamsLayout);

  runtimeText = new FilteredTextEdit;
  runtimeText->setReadOnly(true);
  runtimeText->setAcceptRichText(false);

  QVBoxLayout* mainLayout = new QVBoxLayout;
  mainLayout->addLayout(calcLayout);
  mainLayout->addWidget(paramFileGroup);
  mainLayout->addWidget(integralsFileGroup);
  mainLayout->addWidget(rateParamsGroup);
  mainLayout->addWidget(runtimeText);
  setLayout(mainLayout);
}


void RunTab::calculationTypeChanged(int index) {
  if(index==4) {
    integralsFileGroup->hide();
    rateParamsGroup->show();  
  } else {
    rateParamsGroup->hide();
    integralsFileGroup->show();
  }
  if(index==3) chiVarianceText->setEnabled(true);
  else {
    chiVarianceText->setText("1.0");
    chiVarianceText->setEnabled(false);
  }
}

void RunTab::paramFileButtonChanged(bool checked) {
  if(checked) {
    paramFileText->setEnabled(true);
    paramFileChoose->setEnabled(true);
  } else {
    paramFileText->setEnabled(false);
    paramFileChoose->setEnabled(false);
  }
}

void RunTab::integralsFileButtonChanged(bool checked) {
  if(checked) {
    integralsFileText->setEnabled(true);
    integralsFileChoose->setEnabled(true);
  } else {
    integralsFileText->setEnabled(false);
    integralsFileChoose->setEnabled(false);
  }
}

void RunTab::fileTempButtonChanged(bool checked) {
  if(checked) {
    fileTempText->setEnabled(true);
    rateParamsChoose->setEnabled(true);
    minTempText->setEnabled(false);
    maxTempText->setEnabled(false);
    tempStepText->setEnabled(false);
  } else {
    fileTempText->setEnabled(false);
    rateParamsChoose->setEnabled(false);
    minTempText->setEnabled(true);
    maxTempText->setEnabled(true);
    tempStepText->setEnabled(true);
  }
}

void RunTab::setChooseFile(QLineEdit *lineEdit) {
  QString filename = QFileDialog::getOpenFileName(this);
  if(!filename.isEmpty()) {
    lineEdit->setText(QDir::fromNativeSeparators(filename));
  }
}

void RunTab::reset() {
  calcType->setCurrentIndex(0);
  newParamFileButton->setChecked(true);
  paramFileText->setText("");
  newIntegralsFileButton->setChecked(true);
  integralsFileText->setText("");
  rateEntranceKey->setText("");
  rateExitKey->setText("");
  gridTempButton->setChecked(true);
  minTempText->setText("");
  maxTempText->setText("");
  tempStepText->setText("");
  fileTempText->setText("");
  runtimeText->clear();
}

void RunTab::showInfo(int which, QString title) {
  if(which<infoText.size()) {
    if(!infoDialog[which]) {
      infoDialog[which] = new InfoDialog(infoText[which],this,title);
      infoDialog[which]->setAttribute(Qt:: WA_DeleteOnClose);
      infoDialog[which]->show();
    } else infoDialog[which]->raise();
  }
}
