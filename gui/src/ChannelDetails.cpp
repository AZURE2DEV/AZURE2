#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
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

  // Stacked rather than side by side: the two labels together are wider than
  // the dialog wants to be, and stacking costs one row instead of forcing the
  // whole channel panel wider.
  QVBoxLayout *conventionLayout=new QVBoxLayout;
  conventionLayout->setSpacing(2);
  conventionLayout->addWidget(physicalButton);
  conventionLayout->addWidget(rwaButton);

  // Wigner-limit calculator: a button with its read-only, selectable result
  // field immediately to its right (so the user can copy the value).
  wignerButton = new QPushButton(tr("Wigner Limit"));
  wignerButton->setToolTip(tr("Wigner limit (single-particle width) for this channel, quoted in "
                              "whichever convention is selected above: the reduced width "
                              "amplitude limit in MeV^1/2, or the partial width limit in eV. "
                              "Available for particle channels only, not radiative capture."));
  wignerLimitText = new QLineEdit;
  wignerLimitText->setReadOnly(true);
  wignerLimitText->setMaximumWidth(170);

  QGridLayout *reducedWidthLayout=new QGridLayout;
  reducedWidthLayout->addWidget(normParam,0,0);
  reducedWidthLayout->addWidget(reducedWidthText,0,1);
  reducedWidthLayout->addWidget(normUnits,0,2);
  reducedWidthLayout->addItem(new QSpacerItem(20,20),0,3);
  reducedWidthLayout->setColumnStretch(3,1);

  QHBoxLayout *wignerButtonRow = new QHBoxLayout;
  wignerButtonRow->addWidget(wignerButton);
  wignerButtonRow->addWidget(wignerLimitText);
  wignerButtonRow->addStretch();

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(details);
  mainLayout->addLayout(conventionLayout);
  mainLayout->addLayout(wignerButtonRow);
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
    normParam->setText("Width:");
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
    // Just "Width": the convention is already stated by the radio buttons above
    // and by the units beside the field, and the full name is wide enough to
    // stretch the whole panel.
    normParam->setText("Width:");
    normUnits->setText("MeV^(1/2)");
  } else setNormParam(normParamWhich_);
}
