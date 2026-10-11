#include <QLabel>
#include <QPushButton>
#include <QVBoxLayout>

#include "AboutAZURE2Dialog.h"

AboutAZURE2Dialog::AboutAZURE2Dialog(QWidget *parent) :
  QDialog(parent) {
  setWindowTitle(tr("About AZURE2"));

#ifndef AZURE2_VERSION_STRING
#define AZURE2_VERSION_STRING "unknown"
#endif
#ifndef AZURE2_GIT_VERSION
#define AZURE2_GIT_VERSION ""
#endif
  QString version = QString("v%1").arg(AZURE2_VERSION_STRING);
  if (QString(AZURE2_GIT_VERSION).length() > 0) version += QString(" (%1)").arg(AZURE2_GIT_VERSION);
  QLabel *label = new QLabel(QString("<center><img src=\":/azure-icon.png\" width=\"128\" height=\"128\" /><br/>"
                                     "<b>AZURE2</b><br/>%1<br/><br/>"
                                     "J. Skowronski, E. Uberseder, R.J. deBoer, R.E. Azuma<br/>"
                                     "Joint Institute For Nuclear Astrophysics (JINA)</center>").arg(version),
                             this);
  QPushButton *okButton = new QPushButton(tr("OK"), this);
  okButton->setMaximumSize(80, 30);

  QVBoxLayout *layout = new QVBoxLayout(this);

  layout->addWidget(label);
  layout->addWidget(okButton);
  layout->setAlignment(okButton, Qt::AlignHCenter);

  connect(okButton, SIGNAL(clicked()), this, SLOT(close()));

  setLayout(layout);
}
