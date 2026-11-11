#include "EditOptionsDialog.h"
#include "Config.h"

#include <QCheckBox>
#include <QComboBox>
#include <QVBoxLayout>
#include <QPushButton>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QHBoxLayout>
#include <QFileDialog>
#include <QDir>

EditOptionsDialog::EditOptionsDialog(QWidget *parent) :
  QDialog(parent) {
  useGSLCoulCheck = new QCheckBox(tr("Use GSL Coulomb functions"));
  parametrizationCombo = new QComboBox;
  parametrizationCombo->setObjectName("parametrizationCombo");
  parametrizationCombo->addItem(tr("Formal"), Formal);
  parametrizationCombo->setItemData(Formal, tr("Formal R-matrix parameters (Lane and Thomas): pole energies and reduced-width "
                                               "amplitudes that depend on the boundary conditions; the physical parameters of "
                                               "the project are transformed to them per J^pi group."), Qt::ToolTipRole);
  parametrizationCombo->addItem(tr("Brune"), Brune);
  parametrizationCombo->setItemData(Brune, tr("Brune's alternative parametrization (Phys. Rev. C 66, 044611, 2002): level "
                                              "energies and observed reduced-width amplitudes, no boundary-condition "
                                              "dependence.  Recommended."), Qt::ToolTipRole);
  parametrizationCombo->addItem(tr("Brune-Park"), BrunePark);
  parametrizationCombo->setItemData(BrunePark, tr("Park's rescaled amplitudes on Brune's level matrix (gamma_Park = gamma_Brune "
                                                  "sqrt(J)); the command line's --use-park.  Not with the R-matrix route: a "
                                                  "calculation runs the A-matrix with a warning."), Qt::ToolTipRole);
  parametrizationCombo->setToolTip(tr("R-matrix parametrization.  Like the other options in this box it is kept with the "
                                      "GUI's last-run settings (<lastRun>); command-line runs take it from --use-brune / "
                                      "--no-brune / --use-park instead."));
  parametrizationCombo->setCurrentIndex(Brune);
  ignoreExternalsCheck = new QCheckBox(tr("Ignore external width\nif internal width is zeroed"));
  useRMCCheck = new QCheckBox(tr("Use RMC capture formalism\n(neutron capture only)"));
  noTransformCheck = new QCheckBox(tr("Do not perform parameter\ntransformations"));
  useHybridMethodCheck = new QCheckBox(tr("Use Hybrid Coulomb method"));
  useAdaptiveGridCheck = new QCheckBox(tr("Use adaptive integration grid\nfor target effects"));
  relativisticCheck = new QCheckBox(tr("Relativistic kinematics"));
  relativisticCheck->setObjectName("relativisticCheck");
  relativisticCheck->setToolTip(tr("Channel wave numbers, Sommerfeld parameters and lab <-> c.m. energies from the invariant mass "
                                   "(the kinematics of the LANL code EDA, ENDF KRL=1).  Level energies stay excitation energies; "
                                   "relativistic.out lists them in EDA's E(s) convention.  Saved as relativistic=1 in <config>, "
                                   "so command-line and pyazr runs of the project use it too."));
  // noLongWavelengthCheck = new QCheckBox(tr("Do not use long wavelength\n"
  //					   "approximation for EL external capture"));

  connect(parametrizationCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(parametrizationChanged(int)));
  connect(useRMCCheck, SIGNAL(stateChanged(int)), this, SLOT(useRMCCheckChanged(int)));

  QGroupBox *optionsBox = new QGroupBox(tr("AZURE2 Options"));
  QVBoxLayout *optionsLayout = new QVBoxLayout;
  QHBoxLayout *parametrizationRow = new QHBoxLayout;
  QLabel *parametrizationLabel = new QLabel(tr("Parametrization:"));
  parametrizationLabel->setBuddy(parametrizationCombo);
  parametrizationRow->addWidget(parametrizationLabel);
  parametrizationRow->addWidget(parametrizationCombo, 1);
  optionsLayout->addLayout(parametrizationRow);
  optionsLayout->addWidget(useGSLCoulCheck);
  optionsLayout->addWidget(ignoreExternalsCheck);
  optionsLayout->addWidget(useRMCCheck);
  optionsLayout->addWidget(noTransformCheck);
  optionsLayout->addWidget(useHybridMethodCheck);
  optionsLayout->addWidget(useAdaptiveGridCheck);
  // optionsLayout->addWidget(noLongWavelengthCheck);
  optionsBox->setLayout(optionsLayout);

  // Kinematics is part of the project (<config>), not a GUI-only switch: its own box says so.
  QGroupBox *kinematicsBox = new QGroupBox(tr("Kinematics"));
  kinematicsBox->setObjectName("kinematicsBox");
  QLabel *kinematicsNote = new QLabel(tr("Saved in the project: command-line and\npyazr runs of it use this setting too."));
  kinematicsNote->setObjectName("kinematicsNote");
  QVBoxLayout *kinematicsLayout = new QVBoxLayout;
  kinematicsLayout->addWidget(relativisticCheck);
  kinematicsLayout->addWidget(kinematicsNote);
  kinematicsBox->setLayout(kinematicsLayout);

  // Prior parameter covariance (a fit starts from an earlier evaluation's parameters and covariance): a file, no
  // editor -- pyazr.prior_from_fit writes one from a fit's param.sav and covariance.dat.
  QGroupBox *priorBox = new QGroupBox(tr("Prior Parameter Covariance"));
  priorBox->setObjectName("priorCovarianceBox");
  priorCovarianceEdit = new QLineEdit;
  priorCovarianceEdit->setObjectName("priorCovarianceEdit");
  priorCovarianceEdit->setPlaceholderText(tr("none"));
  priorCovarianceEdit->setToolTip(tr("Multivariate Gaussian prior (p - mu)^T C^-1 (p - mu) on the parameters the file "
                                     "lists (names, means, covariance; see the user guide, \"Prior parameter "
                                     "covariance\").  Saved as priorCovariance=<path> in <config>, relative to the "
                                     "project's directory; empty = none."));
  priorCovarianceButton = new QPushButton(tr("Choose..."));
  priorCovarianceButton->setObjectName("priorCovarianceChoose");
  connect(priorCovarianceButton, SIGNAL(clicked()), this, SLOT(choosePriorCovariance()));
  QHBoxLayout *priorRow = new QHBoxLayout;
  priorRow->addWidget(priorCovarianceEdit);
  priorRow->addWidget(priorCovarianceButton);
  QLabel *priorNote = new QLabel(tr("Saved in the project: command-line and\npyazr runs of it use this file too."));
  QVBoxLayout *priorLayout = new QVBoxLayout;
  priorLayout->addLayout(priorRow);
  priorLayout->addWidget(priorNote);
  priorBox->setLayout(priorLayout);

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
  mainLayout->addWidget(kinematicsBox);
  mainLayout->addWidget(priorBox);
  mainLayout->addLayout(buttonBox);
  setWindowTitle(tr("Edit Options"));
  setLayout(mainLayout);
}

int EditOptionsDialog::parametrization() const { return parametrizationCombo->currentIndex(); }

void EditOptionsDialog::setParametrization(int choice) {
  parametrizationCombo->setCurrentIndex(choice >= Formal && choice <= BrunePark ? choice : Brune);
  parametrizationChanged(parametrizationCombo->currentIndex());  // no signal when the index does not change
}

int EditOptionsDialog::parametrizationFromMask(unsigned int mask) {
  if (mask & Config::USE_PARK_FORMALISM) return BrunePark;
  if (mask & Config::USE_BRUNE_FORMALISM) return Brune;
  return Formal;
}

unsigned int EditOptionsDialog::maskWithParametrization(unsigned int mask, int choice) {
  mask &= ~(unsigned int)(Config::USE_BRUNE_FORMALISM | Config::USE_PARK_FORMALISM);
  if (choice == Brune) mask |= Config::USE_BRUNE_FORMALISM;
  else if (choice == BrunePark) mask |= (Config::USE_BRUNE_FORMALISM | Config::USE_PARK_FORMALISM);
  return mask;
}

void EditOptionsDialog::parametrizationChanged(int choice) {
  // RMC capture is a formal-parametrization option: it cannot be on with Brune or Brune-Park.
  if (choice != Formal) {
    useRMCCheck->setChecked(false);
    useRMCCheck->setEnabled(false);
  } else {
    useRMCCheck->setEnabled(true);
  }
}

void EditOptionsDialog::useRMCCheckChanged(int state) {
  if (state == Qt::Checked) {
    parametrizationCombo->setCurrentIndex(Formal);
    parametrizationCombo->setEnabled(false);
  } else {
    parametrizationCombo->setEnabled(true);
  }
}

void EditOptionsDialog::choosePriorCovariance() {
  QString file = QFileDialog::getOpenFileName(this, tr("Prior Parameter Covariance File"),
                                              projectDirectory_.isEmpty() ? QString() : projectDirectory_);
  if (file.isEmpty()) return;
  if (!projectDirectory_.isEmpty()) {
    const QString rel = QDir(projectDirectory_).relativeFilePath(file);
    if (!rel.startsWith("..")) file = rel;  // inside the project's directory: relative, as the data files
  }
  priorCovarianceEdit->setText(file);
}
