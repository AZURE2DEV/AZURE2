#include <QComboBox>
#include <QGridLayout>
#include <QSpacerItem>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QGroupBox>
#include <QCheckBox>


#include "AddPairDialog.h"

AddPairDialog::AddPairDialog(QWidget *parent) :
  QDialog(parent) {
  // this->setMaximumSize(370,440);
  // this->setMinimumSize(370,440);
  this->setMaximumWidth(370);
  this->setMinimumWidth(370);

  excitationEnergyLabel = new QLabel(tr("Excitation Energy [MeV]:"));
  excitationEnergyText = new QLineEdit;
  seperationEnergyLabel = new QLabel(tr("Separation Energy [MeV]:"));
  seperationEnergyText = new QLineEdit;
  channelRadiusLabel = new QLabel(tr("Channel Radius [fm]:"));
  channelRadiusText = new QLineEdit;
  bindingEnergyLabel = new QLabel(tr("THM Binding Energy [MeV]:"));
  bindingEnergyText = new QLineEdit;
  bindingEnergyText->setText("0");
  pairTypeLabel = new QLabel(tr("Particle Pair Type:"));
  pairTypeCombo = new QComboBox;
  pairTypeCombo->addItem(tr("Particle, Particle"));
  pairTypeCombo->addItem(tr("Particle, Gamma"));
  pairTypeCombo->addItem(tr("Beta Decay"));
  connect(pairTypeCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateLightParticle(int)));

  QRegExp rx("^\\d{0,2}(\\.[05]{0,1})?$");
  QValidator *validator = new QRegExpValidator(rx, this);
  lightJLabel = new QLabel(tr("J:"));
  lightJText = new QLineEdit;
  lightJText->setValidator(validator);
  lightPiLabel = new QLabel(tr("Pi:"));
  lightPiCombo = new QComboBox;
  lightPiCombo->addItem(tr("-"));
  lightPiCombo->addItem(tr("+"));
  lightZLabel = new QLabel(tr("Z:"));
  lightZText = new QLineEdit;
  lightMLabel = new QLabel(tr("M:"));
  lightMText = new QLineEdit;
  // lightGLabel = new QLabel(tr("g:"));
  // lightGText = new QLineEdit;

  heavyJLabel = new QLabel(tr("J:"));
  heavyJText = new QLineEdit;
  heavyJText->setValidator(validator);
  heavyPiLabel = new QLabel(tr("Pi:"));
  heavyPiCombo = new QComboBox;
  heavyPiCombo->addItem(tr("-"));
  heavyPiCombo->addItem(tr("+"));
  heavyZLabel = new QLabel(tr("Z:"));
  heavyZText = new QLineEdit;
  heavyMLabel = new QLabel(tr("M:"));
  heavyMText = new QLineEdit;
  // heavyGLabel = new QLabel(tr("g:"));
  // heavyGText = new QLineEdit;

  e1Check = new QCheckBox(tr("E1"));
  e1Check->setChecked(false);
  // m1Check = new QCheckBox(tr("M1"));
  // m1Check->setChecked(false);
  e2Check = new QCheckBox(tr("E2"));
  e2Check->setChecked(false);

  cancelButton = new QPushButton(tr("Cancel"));
  okButton = new QPushButton(tr("Accept"));
  okButton->setDefault(true);


  QGridLayout *pairTypeLayout = new QGridLayout;
  pairTypeLayout->addWidget(pairTypeLabel, 0, 0);
  pairTypeLayout->addWidget(pairTypeCombo, 0, 1);
  pairTypeLayout->setColumnStretch(1, 1);

  QGroupBox *channelGroup = new QGroupBox(tr("Channel Properties"));
  QGridLayout *channelLayout = new QGridLayout;
  channelLayout->addWidget(excitationEnergyLabel, 0, 0, Qt::AlignRight);
  channelLayout->addWidget(excitationEnergyText, 0, 1);
  channelLayout->addWidget(seperationEnergyLabel, 1, 0, Qt::AlignRight);
  channelLayout->addWidget(seperationEnergyText, 1, 1);
  channelLayout->addWidget(channelRadiusLabel, 2, 0, Qt::AlignRight);
  channelLayout->addWidget(channelRadiusText, 2, 1);
  channelLayout->addWidget(bindingEnergyLabel, 3, 0, Qt::AlignRight);
  channelLayout->addWidget(bindingEnergyText, 3, 1);
  channelGroup->setLayout(channelLayout);

  QGroupBox *lightGroup = new QGroupBox(tr("Light Particle"));
  QGridLayout *lightLayout = new QGridLayout;
  lightLayout->addWidget(lightJLabel, 0, 0, Qt::AlignRight);
  QHBoxLayout *lightSpinLayout = new QHBoxLayout;
  lightSpinLayout->addWidget(lightJText);
  lightSpinLayout->addWidget(lightPiCombo);
  lightLayout->addLayout(lightSpinLayout, 0, 1);
  lightLayout->addWidget(lightZLabel, 2, 0, Qt::AlignRight);
  lightLayout->addWidget(lightZText, 2, 1);
  lightLayout->addWidget(lightMLabel, 3, 0, Qt::AlignRight);
  lightLayout->addWidget(lightMText, 3, 1);
  // lightLayout->addWidget(lightGLabel,4,0,Qt::AlignRight);
  // lightLayout->addWidget(lightGText,4,1);
  lightGroup->setLayout(lightLayout);

  QGroupBox *heavyGroup = new QGroupBox(tr("Heavy Particle"));
  QGridLayout *heavyLayout = new QGridLayout;
  heavyLayout->addWidget(heavyJLabel, 0, 0, Qt::AlignRight);
  QHBoxLayout *heavySpinLayout = new QHBoxLayout;
  heavySpinLayout->addWidget(heavyJText);
  heavySpinLayout->addWidget(heavyPiCombo);
  heavyLayout->addLayout(heavySpinLayout, 0, 1);
  heavyLayout->addWidget(heavyZLabel, 2, 0, Qt::AlignRight);
  heavyLayout->addWidget(heavyZText, 2, 1);
  heavyLayout->addWidget(heavyMLabel, 3, 0, Qt::AlignRight);
  heavyLayout->addWidget(heavyMText, 3, 1);
  // heavyLayout->addWidget(heavyGLabel,4,0,Qt::AlignRight);
  // heavyLayout->addWidget(heavyGText,4,1);
  heavyGroup->setLayout(heavyLayout);

  QHBoxLayout *entryLayout = new QHBoxLayout;
  entryLayout->addWidget(lightGroup);
  entryLayout->addWidget(heavyGroup);

  QHBoxLayout *buttonBox = new QHBoxLayout;
  buttonBox->addWidget(cancelButton);
  buttonBox->addWidget(okButton);

  multBox = new QGroupBox(tr("External Capture Multipolarities"));
  multBox->hide();
  QHBoxLayout *multLayout = new QHBoxLayout;
  multLayout->addWidget(e1Check);
  // multLayout->addWidget(m1Check);
  multLayout->addWidget(e2Check);
  multBox->setLayout(multLayout);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addLayout(pairTypeLayout);
  mainLayout->addLayout(entryLayout);
  mainLayout->addWidget(channelGroup);
  mainLayout->addWidget(multBox);
  mainLayout->addLayout(buttonBox);

  setLayout(mainLayout);

  connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(reject()));

  setWindowTitle(tr("Add a Particle Pair"));
}

void AddPairDialog::updateLightParticle(int index) {
  if (index == 1) {
    lightJText->setText("1.0");
    lightJText->setEnabled(false);
    lightPiCombo->setCurrentIndex(1);
    lightPiCombo->setEnabled(false);
    lightZText->setText("0");
    lightZText->setEnabled(false);
    lightMText->setText("0.0");
    lightMText->setEnabled(false);
    // lightGText->setText("0.0");
    // lightGText->setEnabled(false);
    seperationEnergyText->setText("0.0");
    seperationEnergyText->setEnabled(false);
    excitationEnergyText->setEnabled(true);
    channelRadiusText->setText("0");
    channelRadiusText->setEnabled(false);
    multBox->show();
  } else if (index == 2) {
    lightJText->setEnabled(false);
    lightJText->setText("0.5");
    lightPiCombo->setEnabled(false);
    lightPiCombo->setCurrentIndex(1);
    lightZText->setEnabled(true);
    lightMText->setEnabled(false);
    lightMText->setText("0.0005");
    // lightGText->setEnabled(false);
    // lightGText->setText("2.0023");
    seperationEnergyText->setEnabled(true);
    excitationEnergyText->setEnabled(false);
    excitationEnergyText->setText("0.000");
    channelRadiusText->setEnabled(true);
    multBox->hide();
    this->adjustSize();
  } else {
    lightJText->setEnabled(true);
    lightPiCombo->setEnabled(true);
    lightZText->setEnabled(true);
    lightMText->setEnabled(true);
    // lightGText->setEnabled(true);
    seperationEnergyText->setEnabled(true);
    excitationEnergyText->setEnabled(true);
    channelRadiusText->setEnabled(true);
    multBox->hide();
    this->adjustSize();
  }
}
