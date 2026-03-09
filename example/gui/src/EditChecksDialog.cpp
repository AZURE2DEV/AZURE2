#include <QtGui>
#include "EditChecksDialog.h"

EditChecksDialog::EditChecksDialog(QWidget *parent) : QDialog(parent) {
  compoundCheckCombo=new QComboBox;
  compoundCheckCombo->addItem(tr("None"));
  compoundCheckCombo->addItem(tr("Screen"));
  compoundCheckCombo->addItem(tr("File"));
  boundaryCheckCombo=new QComboBox;
  boundaryCheckCombo->addItem(tr("None"));
  boundaryCheckCombo->addItem(tr("Screen"));
  boundaryCheckCombo->addItem(tr("File"));
  dataCheckCombo=new QComboBox;
  dataCheckCombo->addItem(tr("None"));
  dataCheckCombo->addItem(tr("Screen"));
  dataCheckCombo->addItem(tr("File"));
  lMatrixCheckCombo=new QComboBox;
  lMatrixCheckCombo->addItem(tr("None"));
  lMatrixCheckCombo->addItem(tr("Screen"));
  lMatrixCheckCombo->addItem(tr("File"));
  legendreCheckCombo=new QComboBox;
  legendreCheckCombo->addItem(tr("None"));
  legendreCheckCombo->addItem(tr("Screen"));
  legendreCheckCombo->addItem(tr("File"));
  coulAmpCheckCombo=new QComboBox;
  coulAmpCheckCombo->addItem(tr("None"));
  coulAmpCheckCombo->addItem(tr("Screen"));
  coulAmpCheckCombo->addItem(tr("File"));
  pathwaysCheckCombo=new QComboBox;
  pathwaysCheckCombo->addItem(tr("None"));
  pathwaysCheckCombo->addItem(tr("Screen"));
  pathwaysCheckCombo->addItem(tr("File"));
  angDistsCheckCombo=new QComboBox;
  angDistsCheckCombo->addItem(tr("None"));
  angDistsCheckCombo->addItem(tr("Screen"));
  angDistsCheckCombo->addItem(tr("File"));
  
  QGroupBox *checkFilesBox=new QGroupBox(tr("Check Configuration"));
  QGridLayout *checkFilesLayout = new QGridLayout;
  checkFilesLayout->addWidget(new QLabel(tr("Coumpound Nucleus:")),0,0,Qt::AlignRight);
  checkFilesLayout->addWidget(compoundCheckCombo,0,1);
  checkFilesLayout->addWidget(new QLabel(tr("Boundary Conditions:")),1,0,Qt::AlignRight);
  checkFilesLayout->addWidget(boundaryCheckCombo,1,1);
  checkFilesLayout->addWidget(new QLabel(tr("Data:")),2,0,Qt::AlignRight);
  checkFilesLayout->addWidget(dataCheckCombo,2,1);
  checkFilesLayout->addWidget(new QLabel(tr("L-Matrix,Phases,Penetrabilities:")),3,0,Qt::AlignRight);
  checkFilesLayout->addWidget(lMatrixCheckCombo,3,1);
  checkFilesLayout->addWidget(new QLabel(tr("Legendre Polynomials:")),4,0,Qt::AlignRight);
  checkFilesLayout->addWidget(legendreCheckCombo,4,1);
  checkFilesLayout->addWidget(new QLabel(tr("Coulomb Amplitudes:")),5,0,Qt::AlignRight);
  checkFilesLayout->addWidget(coulAmpCheckCombo,5,1);
  checkFilesLayout->addWidget(new QLabel(tr("Reaction Pathways:")),6,0,Qt::AlignRight);
  checkFilesLayout->addWidget(pathwaysCheckCombo,6,1);
  checkFilesLayout->addWidget(new QLabel(tr("Angular Distributions:")),7,0,Qt::AlignRight);
  checkFilesLayout->addWidget(angDistsCheckCombo,7,1);
  checkFilesLayout->setColumnStretch(0,0);
  checkFilesLayout->setColumnStretch(1,1);
  checkFilesBox->setLayout(checkFilesLayout);

  cancelButton = new QPushButton(tr("Cancel"));
  okButton = new QPushButton(tr("Accept"));
  okButton->setDefault(true);
  connect(okButton, SIGNAL(clicked()),this,SLOT(accept()));
  connect(cancelButton,SIGNAL(clicked()),this,SLOT(reject()));

  QHBoxLayout *buttonBox = new QHBoxLayout;
  buttonBox->addWidget(cancelButton);
  buttonBox->addWidget(okButton);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(checkFilesBox);
  mainLayout->addLayout(buttonBox);

  setWindowTitle(tr("Edit Check Configuration"));

  setLayout(mainLayout);
}
