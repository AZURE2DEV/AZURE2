#include <QComboBox>
#include <QGridLayout>
#include <QSpacerItem>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QGroupBox>
#include <QCheckBox>


#include "AddPairDialog.h"
#include "MassTable.h"
#include "PairsModel.h"

#include <cmath>

AddPairDialog::AddPairDialog(QWidget *parent) :
  QDialog(parent) {
  // this->setMaximumSize(370,440);
  // this->setMinimumSize(370,440);
  this->setMaximumWidth(420);
  this->setMinimumWidth(420);

  excitationEnergyLabel = new QLabel(tr("Excitation Energy [MeV]:"));
  excitationEnergyText = new QLineEdit;
  seperationEnergyLabel = new QLabel(tr("Separation Energy [MeV]:"));
  seperationEnergyText = new QLineEdit;
  channelRadiusLabel = new QLabel(tr("Channel Radius [fm]:"));
  channelRadiusText = new QLineEdit;
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

  lightIsotopeLabel = new QLabel(tr("Isotope:"));
  lightIsotopeText = new QLineEdit;
  lightIsotopeText->setObjectName("lightIsotopeText");
  lightIsotopeText->setPlaceholderText(tr("p, n, d, a, 3He ..."));
  lightIsotopeText->setToolTip(tr("Optional: an isotope (\"p\", \"a\", \"7Li\" or \"Li7\") fills Z and the atomic mass from %1.")
                                   .arg(MassTable::version()));
  heavyIsotopeLabel = new QLabel(tr("Isotope:"));
  heavyIsotopeText = new QLineEdit;
  heavyIsotopeText->setObjectName("heavyIsotopeText");
  heavyIsotopeText->setPlaceholderText(tr("12C, 16O ..."));
  heavyIsotopeText->setToolTip(lightIsotopeText->toolTip());
  connect(lightIsotopeText, SIGNAL(editingFinished()), this, SLOT(lightIsotopeEdited()));
  connect(heavyIsotopeText, SIGNAL(editingFinished()), this, SLOT(heavyIsotopeEdited()));

  fillButton = new QPushButton(tr("Fill masses and S from %1").arg(MassTable::version()));
  fillButton->setObjectName("fillMassesButton");
  fillButton->setToolTip(tr("Atomic masses of both particles, identified by Z and the mass rounded to a mass number, and the "
                            "separation energy S = m1 + m2 - M of the compound nucleus, all from %1.  A gamma or beta pair gets "
                            "the heavy particle's mass only.").arg(MassTable::version()));
  connect(fillButton, SIGNAL(clicked()), this, SLOT(fillFromMassTable()));
  massCheckLabel = new QLabel;
  massCheckLabel->setObjectName("massCheckLabel");
  massCheckLabel->setWordWrap(true);
  massCheckLabel->setTextFormat(Qt::RichText);

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
  channelLayout->addWidget(fillButton, 3, 0, 1, 2);
  channelLayout->addWidget(massCheckLabel, 4, 0, 1, 2);
  channelGroup->setLayout(channelLayout);

  QGroupBox *lightGroup = new QGroupBox(tr("Light Particle"));
  QGridLayout *lightLayout = new QGridLayout;
  lightLayout->addWidget(lightJLabel, 0, 0, Qt::AlignRight);
  QHBoxLayout *lightSpinLayout = new QHBoxLayout;
  lightSpinLayout->addWidget(lightJText);
  lightSpinLayout->addWidget(lightPiCombo);
  lightLayout->addLayout(lightSpinLayout, 0, 1);
  lightLayout->addWidget(lightIsotopeLabel, 1, 0, Qt::AlignRight);
  lightLayout->addWidget(lightIsotopeText, 1, 1);
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
  heavyLayout->addWidget(heavyIsotopeLabel, 1, 0, Qt::AlignRight);
  heavyLayout->addWidget(heavyIsotopeText, 1, 1);
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

  for (QLineEdit *edit : {lightZText, lightMText, heavyZText, heavyMText, seperationEnergyText})
    connect(edit, SIGNAL(textChanged(QString)), this, SLOT(updateMassCheck()));
  connect(pairTypeCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateMassCheck()));

  setWindowTitle(tr("Add a Particle Pair"));
  updateMassCheck();
}

namespace {

// Full precision, as the .azr writer keeps it (LevelsTab::writeNuclearFile).
QString full(double v) { return QString::number(v, 'g', 15); }

}  // namespace

bool AddPairDialog::isotopeEdited(QLineEdit *isotope, QLineEdit *zText, QLineEdit *mText) {
  // editingFinished also fires when focus merely passes through: act only on a new entry, so tabbing past an
  // isotope field never overwrites a mass typed by hand afterwards
  const QString text = isotope->text().trimmed();
  if (text.isEmpty() || !isotope->isEnabled() || !isotope->isModified()) return false;
  isotope->setModified(false);
  const MassTable &table = MassTable::instance();
  int Z = 0, A = 0;
  if (!table.parseNuclide(text, &Z, &A)) {
    massCheckLabel->setText(tr("<span style=\"color:#d65000\">\"%1\" is not an isotope name (e.g. p, a, 12C, C12).</span>")
                                .arg(text.toHtmlEscaped()));
    return false;
  }
  const MassEntry *e = table.lookup(Z, A);
  if (!e) {
    massCheckLabel->setText(tr("<span style=\"color:#d65000\">%1 has no nuclide Z = %2, A = %3.</span>")
                                .arg(MassTable::version()).arg(Z).arg(A));
    return false;
  }
  zText->setText(QString::number(e->Z));
  mText->setText(full(e->atomicMass));
  return true;
}

void AddPairDialog::lightIsotopeEdited() {
  isotopeEdited(lightIsotopeText, lightZText, lightMText);
}

void AddPairDialog::heavyIsotopeEdited() {
  isotopeEdited(heavyIsotopeText, heavyZText, heavyMText);
}

void AddPairDialog::fillFromMassTable() {
  const MassTable &table = MassTable::instance();
  const int type = pairTypeCombo->currentIndex();
  const MassEntry *heavy = table.identify(heavyZText->text().toInt(), heavyMText->text().toDouble());
  const MassEntry *light = type == 0 ? table.identify(lightZText->text().toInt(), lightMText->text().toDouble()) : nullptr;
  if (!heavy || (type == 0 && !light)) {
    massCheckLabel->setText(tr("<span style=\"color:#d65000\">Enter Z and a mass (or an isotope) for %1 first: %2 has no "
                               "nuclide with that charge and mass number.</span>")
                                .arg(type == 0 ? tr("both particles") : tr("the heavy particle"))
                                .arg(MassTable::version()));
    return;
  }
  heavyMText->setText(full(heavy->atomicMass));
  if (type == 0) {
    lightMText->setText(full(light->atomicMass));
    double S = 0.0;
    if (table.separationEnergy(light->Z, light->A, heavy->Z, heavy->A, &S))
      seperationEnergyText->setText(full(S));
  }
  updateMassCheck();
}

void AddPairDialog::updateMassCheck() {
  const MassTable &table = MassTable::instance();
  if (!table.isLoaded()) {
    massCheckLabel->setText(tr("The mass table could not be read."));
    fillButton->setEnabled(false);
    return;
  }
  PairsData pair{};
  const int type = pairTypeCombo->currentIndex();
  pair.pairType = type == 1 ? 10 : (type == 2 ? 20 : 0);
  pair.lightZ = lightZText->text().toInt();
  pair.lightM = lightMText->text().toDouble();
  pair.heavyZ = heavyZText->text().toInt();
  pair.heavyM = heavyMText->text().toDouble();
  pair.seperationEnergy = seperationEnergyText->text().toDouble();
  const PairMassCheck c = table.check(pair);
  QStringList lines;
  if (c.hasS) {
    const bool off = std::fabs(c.dSkeV) > MassTable::kToleranceKeV;
    QString line = tr("%1: S(%2 + %3) = %4 MeV").arg(MassTable::version()).arg(c.light->label()).arg(c.heavy->label())
                       .arg(c.sTable, 0, 'f', 6);
    if (!seperationEnergyText->text().trimmed().isEmpty()) {
      if (off)
        line = QString("<span style=\"color:#d65000\">%1; the S entered differs by %2 keV.</span>").arg(line).arg(c.dSkeV, 0, 'f', 2);
      else
        line += tr(", consistent.");
    }
    lines << line;
  } else if (pair.pairType != 0 && c.heavy) {
    lines << tr("%1: M(%2) = %3 u").arg(MassTable::version()).arg(c.heavy->label()).arg(full(c.heavy->atomicMass));
  }
  if (c.massOffset > 5e-7)
    lines << tr("The masses entered differ from %1 by up to %2 u.").arg(MassTable::version()).arg(c.massOffset, 0, 'g', 3);
  massCheckLabel->setText(lines.join("<br/>"));
  fillButton->setEnabled(c.fillable);
}

void AddPairDialog::updateLightParticle(int index) {
  lightIsotopeText->setEnabled(index == 0);
  if (index != 0) lightIsotopeText->clear();
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
