#include "EditOptionsDialog.h"
#include <QtGui>

EditOptionsDialog::EditOptionsDialog(QWidget* parent) : QDialog(parent) {
  
  useGSLCoulCheck = new QCheckBox(tr("Use GSL Coulomb functions"));
  useBruneCheck = new QCheckBox(tr("Use Brune formalism"));
  ignoreExternalsCheck = new QCheckBox(tr("Ignore external width\nif internal width is zeroed"));
  useRMCCheck = new QCheckBox(tr("Use RMC capture formalism\n(neutron capture only)"));
  noTransformCheck = new QCheckBox(tr("Do not perform parameter\ntransformations"));
  noLongWavelengthCheck = new QCheckBox(tr("Do not use long wavelength\n"
					   "approximation for EL external capture"));

  connect(useBruneCheck,SIGNAL(stateChanged(int)),this,SLOT(useBruneCheckChanged(int)));
  connect(useRMCCheck,SIGNAL(stateChanged(int)),this,SLOT(useRMCCheckChanged(int)));

  QGroupBox* optionsBox = new QGroupBox(tr("AZURE2 Options"));
  QVBoxLayout* optionsLayout = new QVBoxLayout;
  optionsLayout->addWidget(useGSLCoulCheck);
  optionsLayout->addWidget(useBruneCheck);
  optionsLayout->addWidget(ignoreExternalsCheck);
  optionsLayout->addWidget(useRMCCheck);
  optionsLayout->addWidget(noTransformCheck);
  optionsLayout->addWidget(noLongWavelengthCheck);
  optionsBox->setLayout(optionsLayout);

  cancelButton = new QPushButton(tr("Cancel"));
  okButton = new QPushButton(tr("Accept"));
  okButton->setDefault(true);
  connect(okButton, SIGNAL(clicked()),this,SLOT(accept()));
  connect(cancelButton,SIGNAL(clicked()),this,SLOT(reject()));

  QHBoxLayout *buttonBox = new QHBoxLayout;
  buttonBox->addWidget(cancelButton);
  buttonBox->addWidget(okButton);
  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(optionsBox);
  mainLayout->addLayout(buttonBox);
  setWindowTitle(tr("Edit Options"));
  setLayout(mainLayout);
}

void EditOptionsDialog::useBruneCheckChanged(int state) {
  if(state==Qt::Checked) {
    useRMCCheck->setChecked(false);
    useRMCCheck->setEnabled(false);
  } else useRMCCheck->setEnabled(true);
}

void EditOptionsDialog::useRMCCheckChanged(int state) {
  if(state==Qt::Checked) {
    useBruneCheck->setChecked(false);
    useBruneCheck->setEnabled(false);
  } else useBruneCheck->setEnabled(true);
}
