#include "EditOptionsDialog.h"

#include <QCheckBox>
#include <QVBoxLayout>
#include <QPushButton>
#include <QGroupBox>

EditOptionsDialog::EditOptionsDialog(QWidget *parent) :
  QDialog(parent) {
  useGSLCoulCheck = new QCheckBox(tr("Use GSL Coulomb functions"));
  useBruneCheck = new QCheckBox(tr("Use Brune formalism"));
  useParkCheck = new QCheckBox(tr("Use Park parametrization\n(observed reduced widths; implies Brune)"));
  ignoreExternalsCheck = new QCheckBox(tr("Ignore external width\nif internal width is zeroed"));
  useRMCCheck = new QCheckBox(tr("Use RMC capture formalism\n(neutron capture only)"));
  noTransformCheck = new QCheckBox(tr("Do not perform parameter\ntransformations"));
  useHybridMethodCheck = new QCheckBox(tr("Use Hybrid Coulomb method"));
  useAdaptiveGridCheck = new QCheckBox(tr("Use adaptive integration grid\nfor target effects"));
  useThmCheck = new QCheckBox(tr("Use Trojan Horse Method (THM)"));
  useThmCheck->setToolTip(tr("Shows the THM Workspace and the THM controls of the tabs"));
  // noLongWavelengthCheck = new QCheckBox(tr("Do not use long wavelength\n"
  //					   "approximation for EL external capture"));

  connect(useBruneCheck, SIGNAL(stateChanged(int)), this, SLOT(useBruneCheckChanged(int)));
  connect(useParkCheck, SIGNAL(stateChanged(int)), this, SLOT(useParkCheckChanged(int)));
  connect(useRMCCheck, SIGNAL(stateChanged(int)), this, SLOT(useRMCCheckChanged(int)));

  QGroupBox *optionsBox = new QGroupBox(tr("AZURE2 Options"));
  QVBoxLayout *optionsLayout = new QVBoxLayout;
  optionsLayout->addWidget(useGSLCoulCheck);
  optionsLayout->addWidget(useBruneCheck);
  optionsLayout->addWidget(useParkCheck);
  optionsLayout->addWidget(ignoreExternalsCheck);
  optionsLayout->addWidget(useRMCCheck);
  optionsLayout->addWidget(noTransformCheck);
  optionsLayout->addWidget(useHybridMethodCheck);
  optionsLayout->addWidget(useAdaptiveGridCheck);
  optionsLayout->addWidget(useThmCheck);
  // optionsLayout->addWidget(noLongWavelengthCheck);
  optionsBox->setLayout(optionsLayout);

  cancelButton = new QPushButton(tr("Cancel"));
  okButton = new QPushButton(tr("Accept"));
  okButton->setDefault(true);
  connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(reject()));

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
  if (state == Qt::Checked) {
    useRMCCheck->setChecked(false);
    useRMCCheck->setEnabled(false);
  } else {
    // Park's parametrization is Brune's level matrix with rescaled amplitudes;
    // it cannot be on without Brune.
    useParkCheck->setChecked(false);
    useRMCCheck->setEnabled(true);
  }
}

void EditOptionsDialog::useParkCheckChanged(int state) {
  if (state == Qt::Checked) {
    // Park is Brune's level matrix with rescaled amplitudes: Brune is implied
    // and shown as such (checked, greyed out) while Park is on.
    useBruneCheck->setChecked(true);
    useBruneCheck->setEnabled(false);
    useRMCCheck->setChecked(false);
    useRMCCheck->setEnabled(false);
  } else {
    useBruneCheck->setEnabled(!useRMCCheck->isChecked());
  }
}

void EditOptionsDialog::useRMCCheckChanged(int state) {
  if (state == Qt::Checked) {
    useBruneCheck->setChecked(false);
    useBruneCheck->setEnabled(false);
    useParkCheck->setChecked(false);
    useParkCheck->setEnabled(false);
  } else {
    useBruneCheck->setEnabled(true);
    useParkCheck->setEnabled(true);
  }
}
