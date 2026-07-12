#include <QLabel>
#include <QLineEdit>
#include <QGridLayout>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QRadioButton>
#include <QButtonGroup>

#include "ChannelDetails.h"

ChannelDetails::ChannelDetails(QWidget *parent) : QWidget(parent) {
  details=new QLabel;
  QFont font("Monospace");
  font.setStyleHint(QFont::TypeWriter);
  details->setFont(font);
  reducedWidthText = new QLineEdit;
  reducedWidthText->setMaximumWidth(100);
  normParam=new QLabel;
  normUnits=new QLabel;
  normParamWhich_=0;

  // Input-convention selector (particle channels only): 
  // choose if the input width is a partial/ANC
  // or a reduced one.
  physicalButton=new QRadioButton(tr("Physical Width / ANC"));
  rwaButton=new QRadioButton(tr("Reduced Width Amplitude"));
  physicalButton->setChecked(true);
  conventionGroup=new QButtonGroup(this);
  conventionGroup->addButton(physicalButton,0);
  conventionGroup->addButton(rwaButton,1);

  QHBoxLayout *conventionLayout=new QHBoxLayout;
  conventionLayout->addWidget(physicalButton);
  conventionLayout->addWidget(rwaButton);
  conventionLayout->addStretch(1);

  QGridLayout *reducedWidthLayout=new QGridLayout;
  reducedWidthLayout->addWidget(normParam,0,0);
  reducedWidthLayout->addWidget(reducedWidthText,0,1);
  reducedWidthLayout->addWidget(normUnits,0,2);
  reducedWidthLayout->addItem(new QSpacerItem(20,20),0,3);
  reducedWidthLayout->setColumnStretch(3,1);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(details);
  mainLayout->addLayout(conventionLayout);
  mainLayout->addLayout(reducedWidthLayout);
  setLayout(mainLayout);
}

void ChannelDetails::setNormParam(int which) {
  normParamWhich_=which;
  if(which==1) {
    normParam->setText("ANC:");
    normUnits->setText("fm^(-1/2)");
  } else if(which==2) {
    normParam->setText("Mu:");
    normUnits->setText("nm");
  } else if(which==3) {
    normParam->setText("Q:");
    normUnits->setText("b");
  } else if(which==4) {
    normParam->setText("B:");
    normUnits->setText("");
  } else {
    normParam->setText("Partial Width:");
    normUnits->setText("eV");
  }
}

void ChannelDetails::setConventionChoice(bool isParticle, bool isRWA) {
  physicalButton->setVisible(isParticle);
  rwaButton->setVisible(isParticle);
  // block signals so seeding the radios from the model does not bounce a
  // toggled() edit back into the model
  physicalButton->blockSignals(true);
  rwaButton->blockSignals(true);
  if(isParticle&&isRWA) rwaButton->setChecked(true);
  else physicalButton->setChecked(true);
  physicalButton->blockSignals(false);
  rwaButton->blockSignals(false);
  if(isParticle&&isRWA) {
    normParam->setText("Reduced Width Amplitude:");
    normUnits->setText("MeV^(1/2)");
  } else setNormParam(normParamWhich_);
}
