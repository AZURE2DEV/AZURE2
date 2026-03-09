#include <QtGui>
#include "EditDirsDialog.h"
#include "ChooseFileButton.h"

EditDirsDialog::EditDirsDialog(QWidget *parent) : QDialog(parent) {
  this->setMinimumWidth(500);
  QGroupBox *directoryBox=new QGroupBox(tr("Directory Configuration"));
  QGridLayout *directoryLayout = new QGridLayout;  
  directoryLayout->addWidget(new QLabel(tr("Output Directory:")),0,0,Qt::AlignRight);
  outputDirectoryText = new QLineEdit;
  directoryLayout->addWidget(outputDirectoryText,0,1);
  ChooseFileButton *chooseButton=new ChooseFileButton(tr("Choose..."));
  chooseButton->setLineEdit(outputDirectoryText);
  connect(chooseButton,SIGNAL(clicked(QLineEdit*)),this,SLOT(setChooseDirectory(QLineEdit*)));
  directoryLayout->addWidget(chooseButton,0,2);  
  directoryLayout->addWidget(new QLabel(tr("Checks Directory:")),1,0,Qt::AlignRight);
  checksDirectoryText = new QLineEdit;
  directoryLayout->addWidget(checksDirectoryText,1,1);
  chooseButton=new ChooseFileButton(tr("Choose..."));
  chooseButton->setLineEdit(checksDirectoryText);
  connect(chooseButton,SIGNAL(clicked(QLineEdit*)),this,SLOT(setChooseDirectory(QLineEdit*)));
  directoryLayout->addWidget(chooseButton,1,2);  
  directoryBox->setLayout(directoryLayout);

  cancelButton = new QPushButton(tr("Cancel"));
  okButton = new QPushButton(tr("Accept"));
  okButton->setDefault(true);
  connect(okButton, SIGNAL(clicked()),this,SLOT(accept()));
  connect(cancelButton,SIGNAL(clicked()),this,SLOT(reject()));

  QHBoxLayout *buttonBox = new QHBoxLayout;
  buttonBox->addWidget(cancelButton);
  buttonBox->addWidget(okButton);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(directoryBox);
  mainLayout->addLayout(buttonBox);

  setWindowTitle(tr("Edit Directory Configuration"));

  setLayout(mainLayout);
}

void EditDirsDialog::setChooseDirectory(QLineEdit *lineEdit) {
  QString filename = QFileDialog::getExistingDirectory(this);
  if(!filename.isEmpty()) {
    lineEdit->setText(QDir::fromNativeSeparators(filename)+'/');
  }
}
