#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QGridLayout>
#include <QVBoxLayout>
#include <QHBoxLayout>

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

  // Wigner-limit calculator: a button with its read-only, selectable result
  // field immediately to its right (so the user can copy the value).
  wignerButton = new QPushButton(tr("Wigner Limit"));
  wignerButton->setToolTip(tr("Physical Wigner limit (single-particle width) for this channel. "
                              "Available for particle channels only, not radiative capture."));
  wignerLimitText = new QLineEdit;
  wignerLimitText->setReadOnly(true);
  wignerLimitText->setMaximumWidth(140);

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
  mainLayout->addLayout(wignerButtonRow);
  mainLayout->addLayout(reducedWidthLayout);
  setLayout(mainLayout);
}

void ChannelDetails::setNormParam(int which) {
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
