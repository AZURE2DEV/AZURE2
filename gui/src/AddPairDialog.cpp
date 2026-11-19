#include <QCheckBox>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QDoubleValidator>
#include <QFormLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QIntValidator>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QRegularExpression>
#include <QSpinBox>
#include <QToolButton>
#include <QVBoxLayout>

#include "AddPairDialog.h"
#include "MassTable.h"
#include "PairsModel.h"

#include <cmath>

namespace {

// Full precision, as the .azr writer keeps it (LevelsTab::writeNuclearFile).
QString full(double v) { return QString::number(v, 'g', 15); }

/// The elements of the mass table, by Z: (Z, symbol), the neutron first as "n".
QList<QPair<int, QString>> Elements() {
  QList<QPair<int, QString>> out;
  QSet<int> seen;
  for (const QString &label : MassTable::instance().labels()) {
    int Z = 0, A = 0;
    if (!MassTable::instance().parseNuclide(label, &Z, &A) || seen.contains(Z)) continue;
    seen.insert(Z);
    QString symbol = label;
    symbol.remove(QRegularExpression("^\\d+"));
    out << qMakePair(Z, symbol);
  }
  std::sort(out.begin(), out.end(), [](const QPair<int, QString> &a, const QPair<int, QString> &b) { return a.first < b.first; });
  return out;
}

/// The ground state's J and parity where they are known without a spin table: the light particles (n, p, d, t, 3He,
/// alpha) and the 0+ of every even-even nucleus.  False otherwise (the user gives them).
bool KnownGroundState(int Z, int A, double &J, int &parity) {
  struct Light {
    int Z, A;
    double J;
    int parity;
  };
  static const Light light[] = {{0, 1, 0.5, 1}, {1, 1, 0.5, 1}, {1, 2, 1.0, 1}, {1, 3, 0.5, 1}, {2, 3, 0.5, 1}, {2, 4, 0.0, 1}};
  for (const Light &l : light)
    if (l.Z == Z && l.A == A) {
      J = l.J;
      parity = l.parity;
      return true;
    }
  if (Z % 2 == 0 && (A - Z) % 2 == 0 && A > 0) {
    J = 0.0;
    parity = 1;
    return true;
  }
  return false;
}

void SetParity(QToolButton *b, int parity) {
  b->setChecked(parity >= 0);
  b->setText(parity >= 0 ? "+" : QString::fromUtf8("−"));
}

}  // namespace

void AddPairDialog::buildParticle(Particle &p, const QString &title, QGroupBox *&group) {
  group = new QGroupBox(title);
  QFormLayout *form = new QFormLayout(group);
  p.element = new QComboBox;
  p.element->setEditable(true);
  p.element->setInsertPolicy(QComboBox::NoInsert);
  p.element->setToolTip(tr("The element, by its symbol (or type it)"));
  for (const auto &e : Elements()) p.element->addItem(QString("%1  (Z = %2)").arg(e.second).arg(e.first), e.first);
  p.element->setCurrentIndex(-1);
  p.massNumber = new QSpinBox;
  p.massNumber->setRange(1, 320);
  p.massNumber->setToolTip(tr("The mass number A"));
  p.nucleusRow = new QWidget;
  QHBoxLayout *n = new QHBoxLayout(p.nucleusRow);
  n->setContentsMargins(0, 0, 0, 0);
  n->addWidget(p.element, 1);
  n->addWidget(new QLabel(tr("A")));
  n->addWidget(p.massNumber);
  form->addRow(tr("Nucleus:"), p.nucleusRow);
  p.z = new QLineEdit;
  p.z->setValidator(new QIntValidator(-1, 130, p.z));
  p.mass = new QLineEdit;
  p.mass->setValidator(new QDoubleValidator(0.0, 1e4, 15, p.mass));
  p.mass->setToolTip(tr("The atomic mass in u (%1 when the nucleus is chosen)").arg(MassTable::version()));
  form->addRow(tr("Z:"), p.z);
  form->addRow(tr("M [u]:"), p.mass);
  p.j = new QDoubleSpinBox;
  p.j->setDecimals(1);
  p.j->setSingleStep(0.5);
  p.j->setRange(0.0, 30.0);
  p.j->setToolTip(tr("The spin J (in halves)"));
  p.parity = new QToolButton;
  p.parity->setCheckable(true);
  p.parity->setToolTip(tr("The parity: click to switch between + and -"));
  SetParity(p.parity, 1);
  connect(p.parity, &QToolButton::toggled, this, [b = p.parity](bool on) { SetParity(b, on ? 1 : -1); });
  QWidget *spin = new QWidget;
  QHBoxLayout *s = new QHBoxLayout(spin);
  s->setContentsMargins(0, 0, 0, 0);
  s->addWidget(p.j);
  s->addWidget(p.parity);
  s->addStretch(1);
  form->addRow(tr("J and parity:"), spin);
  p.note = new QLabel;
  p.note->setWordWrap(false);
  p.note->setStyleSheet("color: palette(mid);");
  form->addRow(p.note);
  Particle *which = &p;
  connect(p.element, QOverload<int>::of(&QComboBox::activated), this, [this, which]() { nucleusChosen(*which); });
  connect(p.element->lineEdit(), &QLineEdit::editingFinished, this, [this, which]() {
    // a symbol typed: the item of that symbol
    const QString typed = which->element->currentText().trimmed().section(' ', 0, 0);
    for (int i = 0; i < which->element->count(); i++)
      if (which->element->itemText(i).section(' ', 0, 0).compare(typed, Qt::CaseInsensitive) == 0) {
        if (which->element->currentIndex() != i || which->z->text() != QString::number(which->element->itemData(i).toInt())) {
          which->element->setCurrentIndex(i);
          nucleusChosen(*which);
        }
        return;
      }
  });
  connect(p.massNumber, QOverload<int>::of(&QSpinBox::valueChanged), this, [this, which]() {
    if (which->element->currentIndex() >= 0) nucleusChosen(*which);
  });
  connect(p.z, &QLineEdit::textEdited, this, [which]() {
    const int at = which->element->findData(which->z->text().toInt());
    which->element->blockSignals(true);
    which->element->setCurrentIndex(at);
    which->element->blockSignals(false);
  });
}

bool AddPairDialog::chooseNucleus(Particle &p, const QString &symbol, int A) {
  for (int i = 0; i < p.element->count(); i++)
    if (p.element->itemText(i).section(' ', 0, 0) == symbol) {
      p.massNumber->blockSignals(true);
      p.massNumber->setValue(A);
      p.massNumber->blockSignals(false);
      p.element->setCurrentIndex(i);
      nucleusChosen(p);
      return true;
    }
  return false;
}

void AddPairDialog::nucleusChosen(Particle &p) {
  if (p.element->currentIndex() < 0) return;
  const int Z = p.element->currentData().toInt(), A = p.massNumber->value();
  const MassEntry *e = MassTable::instance().lookup(Z, A);
  p.z->setText(QString::number(Z));
  if (!e) {
    p.note->setText(tr("%1 has no nuclide Z = %2, A = %3: give the mass by hand.").arg(MassTable::version()).arg(Z).arg(A));
    updateMassCheck();
    return;
  }
  p.mass->setText(full(e->atomicMass));
  double J = 0.0;
  int parity = 1;
  if (KnownGroundState(Z, A, J, parity)) {
    p.j->setValue(J);
    SetParity(p.parity, parity);
    p.note->setText(tr("%1: M, J and parity filled").arg(e->label()));
  } else {
    p.note->setText(tr("%1: M filled; give J and parity").arg(e->label()));
    p.note->setToolTip(tr("The mass table has no spins: give the J and parity of the state of the pair (the ground "
                          "state's from a level table)."));
  }
  updateMassCheck();
}

void AddPairDialog::showNucleus(Particle &p) {
  const MassEntry *e = MassTable::instance().identify(p.z->text().toInt(), p.mass->text().toDouble());
  p.element->blockSignals(true);
  p.massNumber->blockSignals(true);
  p.element->setCurrentIndex(p.element->findData(p.z->text().toInt()));
  if (e) p.massNumber->setValue(e->A);
  else if (p.mass->text().toDouble() >= 1.0) p.massNumber->setValue((int)std::lround(p.mass->text().toDouble()));
  p.element->blockSignals(false);
  p.massNumber->blockSignals(false);
  p.note->clear();
}

AddPairDialog::AddPairDialog(QWidget *parent) : QDialog(parent) {
  setObjectName("addPairDialog");
  setMinimumWidth(560);
  pairTypeCombo = new QComboBox;
  pairTypeCombo->setObjectName("pairType");
  pairTypeCombo->addItem(tr("Particle, Particle"));
  pairTypeCombo->addItem(tr("Particle, Gamma (capture to a level)"));
  pairTypeCombo->addItem(tr("Beta Decay"));
  pairTypeCombo->setToolTip(tr("A particle pair (two nuclei), a photon pair (the compound nucleus in a final level) or "
                               "a beta decay"));

  buildParticle(light, tr("Light Particle"), lightGroup);
  buildParticle(heavy, tr("Heavy Particle"), heavyGroup);

  excitationEnergyText = new QLineEdit("0");
  excitationEnergyText->setToolTip(tr("Excitation energy of the heavy particle (for a gamma pair: of the final level)."));
  lightExcitationText = new QLineEdit("0");
  lightExcitationText->setObjectName("lightExcitationText");
  lightExcitationText->setToolTip(tr("Excitation energy of the light particle (particle pairs only): a scattering or "
                                     "reaction of a particle in an excited state.  Give the light particle the J and "
                                     "parity of that state.  The channel opens at S + heavy + light excitation; two "
                                     "particles with the same Z, M, J, parity and excitation are identical."));
  seperationEnergyText = new QLineEdit;
  channelRadiusText = new QLineEdit;
  for (QLineEdit *e : {excitationEnergyText, lightExcitationText, seperationEnergyText, channelRadiusText})
    e->setValidator(new QDoubleValidator(-1e6, 1e6, 15, e));

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
  e2Check = new QCheckBox(tr("E2"));

  QFormLayout *top = new QFormLayout;
  top->addRow(tr("Pair type:"), pairTypeCombo);

  QGroupBox *channelGroup = new QGroupBox(tr("Channel Properties"));
  QFormLayout *channel = new QFormLayout(channelGroup);
  auto row = [channel](const QString &label, QWidget *field, QWidget *&host) {
    host = field;
    channel->addRow(label, field);
  };
  row(tr("Heavy excitation [MeV]:"), excitationEnergyText, excitationRow);
  row(tr("Light excitation [MeV]:"), lightExcitationText, lightExcitationRow);
  row(tr("Separation energy [MeV]:"), seperationEnergyText, separationRow);
  row(tr("Channel radius [fm]:"), channelRadiusText, radiusRow);
  channel->addRow(fillButton);
  channel->addRow(massCheckLabel);

  QHBoxLayout *particles = new QHBoxLayout;
  particles->addWidget(lightGroup);
  particles->addWidget(heavyGroup);

  multBox = new QGroupBox(tr("External Capture Multipolarities"));
  multBox->hide();
  QHBoxLayout *multLayout = new QHBoxLayout(multBox);
  multLayout->addWidget(e1Check);
  multLayout->addWidget(e2Check);
  multLayout->addStretch(1);

  QDialogButtonBox *buttons = new QDialogButtonBox;
  okButton = buttons->addButton(tr("Accept"), QDialogButtonBox::AcceptRole);
  cancelButton = buttons->addButton(QDialogButtonBox::Cancel);
  okButton->setDefault(true);
  connect(buttons, &QDialogButtonBox::accepted, this, &QDialog::accept);
  connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);

  QVBoxLayout *mainLayout = new QVBoxLayout(this);
  mainLayout->addLayout(top);
  mainLayout->addLayout(particles);
  mainLayout->addWidget(channelGroup);
  mainLayout->addWidget(multBox);
  mainLayout->addWidget(buttons);

  connect(pairTypeCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateLightParticle(int)));
  for (QLineEdit *edit : {light.z, light.mass, heavy.z, heavy.mass, seperationEnergyText})
    connect(edit, SIGNAL(textChanged(QString)), this, SLOT(updateMassCheck()));
  connect(pairTypeCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateMassCheck()));

  setWindowTitle(tr("Add a Particle Pair"));
  updateLightParticle(0);
  updateMassCheck();
}

PairsData AddPairDialog::pair() const {
  PairsData p{};
  const int type = pairTypeCombo->currentIndex();
  p.pairType = type == 1 ? 10 : type == 2 ? 20 : 0;
  p.lightJ = light.j->value();
  p.lightPi = light.parity->isChecked() ? 1 : -1;
  p.lightZ = light.z->text().toInt();
  p.lightM = light.mass->text().toDouble();
  p.lightG = lightG_;
  p.heavyJ = heavy.j->value();
  p.heavyPi = heavy.parity->isChecked() ? 1 : -1;
  p.heavyZ = heavy.z->text().toInt();
  p.heavyM = heavy.mass->text().toDouble();
  p.heavyG = heavyG_;
  p.excitationEnergy = excitationEnergyText->text().toDouble();
  p.lightExcitation = type == 0 ? lightExcitationText->text().toDouble() : 0.;
  p.seperationEnergy = seperationEnergyText->text().toDouble();
  p.channelRadius = channelRadiusText->text().toDouble();
  unsigned int mask = 0;
  if (e1Check->isChecked()) mask |= (1 << 0);
  if (e2Check->isChecked()) mask |= (1 << 2);
  p.ecMultMask = (int)mask;
  return p;
}

void AddPairDialog::setPair(const PairsData &p) {
  pairTypeCombo->setCurrentIndex(p.pairType == 10 ? 1 : p.pairType == 20 ? 2 : 0);
  light.j->setValue(p.lightJ);
  SetParity(light.parity, p.lightPi);
  light.z->setText(QString::number(p.lightZ));
  light.mass->setText(full(p.lightM));
  lightG_ = p.lightG;
  heavy.j->setValue(p.heavyJ);
  SetParity(heavy.parity, p.heavyPi);
  heavy.z->setText(QString::number(p.heavyZ));
  heavy.mass->setText(full(p.heavyM));
  heavyG_ = p.heavyG;
  excitationEnergyText->setText(full(p.excitationEnergy));
  lightExcitationText->setText(full(p.lightExcitation));
  seperationEnergyText->setText(full(p.seperationEnergy));
  channelRadiusText->setText(full(p.channelRadius));
  e1Check->setChecked(p.ecMultMask & (1 << 0));
  e2Check->setChecked(p.ecMultMask & (1 << 2));
  showNucleus(light);
  showNucleus(heavy);
  updateMassCheck();
}

void AddPairDialog::fillFromMassTable() {
  const MassTable &table = MassTable::instance();
  const int type = pairTypeCombo->currentIndex();
  const MassEntry *h = table.identify(heavy.z->text().toInt(), heavy.mass->text().toDouble());
  const MassEntry *l = type == 0 ? table.identify(light.z->text().toInt(), light.mass->text().toDouble()) : nullptr;
  if (!h || (type == 0 && !l)) {
    massCheckLabel->setText(tr("<span style=\"color:#d65000\">Choose the nucleus (or give Z and a mass) of %1 first: %2 has "
                               "no nuclide with that charge and mass number.</span>")
                                .arg(type == 0 ? tr("both particles") : tr("the heavy particle"))
                                .arg(MassTable::version()));
    return;
  }
  heavy.mass->setText(full(h->atomicMass));
  if (type == 0) {
    light.mass->setText(full(l->atomicMass));
    double S = 0.0;
    if (table.separationEnergy(l->Z, l->A, h->Z, h->A, &S)) seperationEnergyText->setText(full(S));
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
  const PairsData pair = this->pair();
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
  if (c.massOffset > MassTable::kMassFlagU)
    lines << QString("<span style=\"color:#d65000\">%1</span>")
                 .arg(tr("The masses entered differ from %1 by up to %2 u.").arg(MassTable::version()).arg(c.massOffset, 0, 'g', 3));
  else if (c.massOffset > 5e-7)
    lines << tr("The masses entered differ from %1 by up to %2 u.").arg(MassTable::version()).arg(c.massOffset, 0, 'g', 3);
  massCheckLabel->setText(lines.join("<br/>"));
  fillButton->setEnabled(c.fillable);
}

void AddPairDialog::updateLightParticle(int index) {
  // what the type does not use is hidden; the values a type fixes are set (the GUI's conventions of <levels>)
  auto setRow = [](QWidget *field, bool on) {
    if (QFormLayout *f = qobject_cast<QFormLayout *>(field->parentWidget() ? field->parentWidget()->layout() : nullptr))
      f->setRowVisible(field, on);
    else
      field->setVisible(on);
  };
  const bool particle = index == 0, photon = index == 1, beta = index == 2;
  // leaving the particle type keeps what it held, and coming back gives it back (a look at another type undoes nothing)
  if (!particle && lastType_ == 0) {
    held_ = pair();
    haveHeld_ = true;
  }
  if (particle && lastType_ != 0 && haveHeld_) {
    light.j->setValue(held_.lightJ);
    SetParity(light.parity, held_.lightPi);
    light.z->setText(QString::number(held_.lightZ));
    light.mass->setText(full(held_.lightM));
    seperationEnergyText->setText(full(held_.seperationEnergy));
    channelRadiusText->setText(full(held_.channelRadius));
    lightExcitationText->setText(full(held_.lightExcitation));
    excitationEnergyText->setText(full(held_.excitationEnergy));
    showNucleus(light);
  }
  lastType_ = index;
  lightGroup->setVisible(!photon);
  setRow(lightExcitationRow, particle);
  setRow(separationRow, !photon);
  setRow(radiusRow, !photon);
  setRow(excitationRow, !beta);
  if (QFormLayout *f = qobject_cast<QFormLayout *>(lightGroup->layout())) {
    f->setRowVisible(light.nucleusRow, particle);
    f->setRowVisible(light.mass, particle);
    f->setRowVisible(light.j->parentWidget(), particle);
    f->setRowVisible(light.note, particle);
  }
  multBox->setVisible(photon);
  if (photon) {
    light.j->setValue(1.0);
    SetParity(light.parity, 1);
    light.z->setText("0");
    light.mass->setText("0.0");
    seperationEnergyText->setText("0.0");
    lightExcitationText->setText("0");
    channelRadiusText->setText("0");
  } else if (beta) {
    light.j->setValue(0.5);
    SetParity(light.parity, 1);
    light.mass->setText("0.0005");
    light.z->setToolTip(tr("The beta's charge: -1 for beta-, +1 for beta+"));
    excitationEnergyText->setText("0.000");
    lightExcitationText->setText("0");
  }
  adjustSize();
}
