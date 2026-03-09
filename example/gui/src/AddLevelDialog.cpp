#include "AddLevelDialog.h"

AddLevelDialog::AddLevelDialog(QWidget *parent) : QDialog(parent) {

  this->setMaximumSize(250,150);
  //this->setMinimumSize(220,110);

  QRegExp rx("^\\d{0,2}(\\.[05]{0,1})?$");
  QValidator *validator = new QRegExpValidator(rx, this);
  jValueLabel = new QLabel(tr("Level Spin:"));
  jValueText = new QLineEdit;
  jValueText->setValidator(validator);
  piValueCombo = new QComboBox;
  piValueCombo->addItem("-");
  piValueCombo->addItem("+");

  energyLabel = new QLabel(tr("Excitation Energy [MeV]:"));
  energyText = new QLineEdit;

  cancelButton = new QPushButton(tr("Cancel"));
  okButton = new QPushButton(tr("Accept"));
  okButton->setDefault(true);

  QHBoxLayout *energyLayout = new QHBoxLayout;
  energyLayout->addWidget(energyLabel);
  energyLayout->addWidget(energyText);
  
  QGridLayout *spinLayout = new QGridLayout;
  spinLayout->addWidget(jValueLabel,0,0);
  spinLayout->addWidget(jValueText,0,1);
  spinLayout->addWidget(piValueCombo,0,2);
  spinLayout->setColumnStretch(1,1);

  QHBoxLayout *buttonBox = new QHBoxLayout;
  buttonBox->addWidget(cancelButton);
  buttonBox->addWidget(okButton);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addLayout(energyLayout);
  mainLayout->addLayout(spinLayout);
  mainLayout->addLayout(buttonBox);

  setLayout(mainLayout);

  connect(okButton, SIGNAL(clicked()),this,SLOT(accept()));
  connect(cancelButton,SIGNAL(clicked()),this,SLOT(reject()));

  setWindowTitle(tr("Add a Level"));
}
