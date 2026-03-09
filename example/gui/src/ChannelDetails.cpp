#include "ChannelDetails.h"
#include <QtGui>

ChannelDetails::ChannelDetails(QWidget *parent) : QWidget(parent) {
  details=new QLabel;
  QFont font("Monospace");
  font.setStyleHint(QFont::TypeWriter);  
  details->setFont(font);
  reducedWidthText = new QLineEdit;
  reducedWidthText->setMaximumWidth(100);
  normParam=new QLabel;
  normUnits=new QLabel;

  QGridLayout *reducedWidthLayout=new QGridLayout;
  reducedWidthLayout->addWidget(normParam,0,0);
  reducedWidthLayout->addWidget(reducedWidthText,0,1);
  reducedWidthLayout->addWidget(normUnits,0,2);
  reducedWidthLayout->addItem(new QSpacerItem(20,20),0,3);
  reducedWidthLayout->setColumnStretch(3,1);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(details);
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
