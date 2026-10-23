#include <QHeaderView>
#include <QGridLayout>
#include <QComboBox>
#include <QLineEdit>
#include <QCheckBox>
#include <QTextStream>
#include <QMessageBox>

#include <QInputDialog>
#include "PairsTab.h"
#include "RichTextDelegate.h"
#include "InfoDialog.h"
#include <QLabel>
#include <QRegularExpression>
#include <iostream>

PairsTab::PairsTab(QWidget *parent) :
  QWidget(parent) {
  pairsModel = new PairsModel(this);

  pairsView = new QTableView;
  pairsView->setModel(pairsModel);
  pairsView->verticalHeader()->setHighlightSections(false);
  pairsView->horizontalHeader()->setHighlightSections(false);
  pairsView->setColumnHidden(1, true);
  pairsView->setColumnHidden(3, true);
  pairsView->setColumnHidden(4, true);
  pairsView->setColumnHidden(6, true);
  pairsView->setColumnHidden(8, true);
  pairsView->setColumnHidden(9, true);
  pairsView->setColumnHidden(13, true);
  pairsView->setColumnHidden(14, true);
  RichTextDelegate *rt = new RichTextDelegate();
  pairsView->setItemDelegateForColumn(0, rt);
  pairsView->setItemDelegateForColumn(5, rt);
  pairsView->setItemDelegateForColumn(PairsModel::NucleusColumn, rt);
  pairsView->setColumnHidden(PairsModel::NucleusColumn, true);   // shown once the pairs form more than one nucleus
  connect(pairsModel, SIGNAL(rowsInserted(QModelIndex, int, int)), this, SLOT(updateNucleusColumn()));
  connect(pairsModel, SIGNAL(rowsRemoved(QModelIndex, int, int)), this, SLOT(updateNucleusColumn()));
  connect(pairsModel, SIGNAL(dataChanged(QModelIndex, QModelIndex)), this, SLOT(updateNucleusColumn()));
  connect(pairsModel, SIGNAL(modelReset()), this, SLOT(updateNucleusColumn()));
  pairsView->setSelectionBehavior(QAbstractItemView::SelectRows);
  pairsView->setSelectionMode(QAbstractItemView::SingleSelection);
  pairsView->setEditTriggers(QAbstractItemView::NoEditTriggers);
  pairsView->setShowGrid(false);

  connect(pairsView->selectionModel(), SIGNAL(selectionChanged(QItemSelection, QItemSelection)), this, SLOT(updateButtons(QItemSelection)));
  connect(pairsView, SIGNAL(doubleClicked(QModelIndex)), this, SLOT(editPair()));

  for (int i = 0; i < pairsModel->columnCount(QModelIndex()); i++) pairsView->horizontalHeader()->setSectionResizeMode(i, QHeaderView::Stretch);

  addButton = new QPushButton(tr("+"));
  addButton->setMaximumSize(28, 28);
  connect(addButton, SIGNAL(clicked()), this, SLOT(addPair()));
  deleteButton = new QPushButton(tr("-"));
  deleteButton->setMaximumSize(28, 28);
  deleteButton->setEnabled(false);
  connect(deleteButton, SIGNAL(clicked()), this, SLOT(removePair()));

  /*
  mapper = new QSignalMapper(this);
  connect(mapper,SIGNAL(mapped(int)),this,SLOT(showInfo(int)));

  infoButton[0] = new QPushButton(this);
  infoButton[0]->setMaximumSize(28,28);
  infoButton[0]->setIcon(style()->standardIcon(QStyle::SP_MessageBoxInformation));
  mapper->setMapping(infoButton[0],0);
  connect(infoButton[0],SIGNAL(clicked()),mapper,SLOT(map()));
  */

  QGridLayout *buttonBox = new QGridLayout;
  buttonBox->addWidget(addButton, 0, 0);
  buttonBox->addWidget(deleteButton, 0, 1);
  buttonBox->addItem(new QSpacerItem(28, 28), 0, 2);
  // buttonBox->addWidget(infoButton[0],0,3);
  buttonBox->setColumnStretch(0, 0);
  buttonBox->setColumnStretch(1, 0);
  buttonBox->setColumnStretch(2, 1);
  buttonBox->setColumnStretch(3, 0);
#ifdef MACX_SPACING
  buttonBox->setHorizontalSpacing(11);
#else
  buttonBox->setHorizontalSpacing(0);
#endif

  // Isotopic abundances: one row per compound nucleus, shown once the pairs
  // form more than one (a natural-abundance target), hidden otherwise.
  abundancesModel = new AbundancesModel(this);
  abundancesModel->setPairsModel(pairsModel);
  abundancesView = new QTableView;
  abundancesView->setModel(abundancesModel);
  abundancesView->verticalHeader()->setVisible(false);
  abundancesView->horizontalHeader()->setHighlightSections(false);
  abundancesView->setSelectionBehavior(QAbstractItemView::SelectRows);
  abundancesView->setSelectionMode(QAbstractItemView::SingleSelection);
  abundancesView->setShowGrid(false);
  for (int i = 0; i < AbundancesModel::NumCols; i++) abundancesView->horizontalHeader()->setSectionResizeMode(i, QHeaderView::Stretch);
  abundancesView->setMaximumHeight(140);
  sumSigmaSpin = new QDoubleSpinBox;
  sumSigmaSpin->setDecimals(4);
  sumSigmaSpin->setRange(0.0, 10.0);
  sumSigmaSpin->setSingleStep(0.001);
  sumSigmaSpin->setValue(0.0);
  sumSigmaSpin->setToolTip(tr("Width of the sum-to-one prior on the fitted abundances, ((sum w) - 1)^2 / sigma^2; 0 switches it off."));
  QHBoxLayout *sumRow = new QHBoxLayout;
  sumRow->addWidget(new QLabel(tr("Default sum-to-one prior sigma (0 = off):")));
  sumRow->addWidget(sumSigmaSpin);
  sumRow->addStretch(1);
  // Samples: a project may hold several abundance sets (a natural target and an
  // enriched one, say); the table shows one at a time.
  sampleCombo = new QComboBox;
  sampleCombo->setObjectName("sampleCombo");
  sampleCombo->setToolTip(tr("The sample whose abundances the table shows.  A segment chooses its sample in the segment dialog "
                             "(\"sample 1\" is the default sample, \"sample <name>\" a named one)."));
  addSampleButton = new QPushButton(tr("+ sample"));
  addSampleButton->setToolTip(tr("Add a named sample (one word), e.g. an enriched target."));
  removeSampleButton = new QPushButton(tr("- sample"));
  removeSampleButton->setEnabled(false);
  sampleSumSigmaSpin = new QDoubleSpinBox;
  sampleSumSigmaSpin->setDecimals(4);
  sampleSumSigmaSpin->setRange(-1.0, 10.0);
  sampleSumSigmaSpin->setSingleStep(0.001);
  sampleSumSigmaSpin->setSpecialValueText(tr("default"));
  sampleSumSigmaSpin->setValue(-1.0);
  sampleSumSigmaSpin->setToolTip(tr("This sample's own sum-to-one prior sigma: ((sum w) - 1)^2 / sigma^2 over its abundances; "
                                    "0 switches it off for this sample, \"default\" uses the default sigma above."));
  QHBoxLayout *sampleRow = new QHBoxLayout;
  sampleRow->addWidget(new QLabel(tr("Sample:")));
  sampleRow->addWidget(sampleCombo, 1);
  sampleRow->addWidget(addSampleButton);
  sampleRow->addWidget(removeSampleButton);
  sampleRow->addWidget(new QLabel(tr("  its sum-to-one sigma:")));
  sampleRow->addWidget(sampleSumSigmaSpin);
  connect(sampleCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(sampleChosen(int)));
  connect(addSampleButton, SIGNAL(clicked()), this, SLOT(addSample()));
  connect(removeSampleButton, SIGNAL(clicked()), this, SLOT(removeSample()));
  connect(sampleSumSigmaSpin, SIGNAL(valueChanged(double)), this, SLOT(sampleSumSigmaEdited(double)));
  connect(abundancesModel, SIGNAL(samplesChanged()), this, SLOT(updateSampleCombo()));
  QVBoxLayout *abundancesLayout = new QVBoxLayout;
  abundancesLayout->addWidget(new QLabel(tr("Each nucleus's term in a sample segment (\"Sample\" in the segment dialog) is weighted by the "
                                            "abundance of its target isotope in the segment's sample.  A fitted abundance becomes the parameter "
                                            "abundance[_<sample>]_<Z>_<A>.")));
  abundancesLayout->addLayout(sampleRow);
  abundancesLayout->addWidget(abundancesView);
  abundancesLayout->addLayout(sumRow);
  updateSampleCombo();
  abundancesBox = new QGroupBox(tr("Isotopic abundances (natural-abundance target: one row per compound nucleus)"));
  abundancesBox->setLayout(abundancesLayout);
  abundancesBox->setVisible(false);
  connect(abundancesModel, SIGNAL(modelReset()), this, SLOT(updateAbundancesBox()));

  QVBoxLayout *layout = new QVBoxLayout();
  layout->addWidget(pairsView);
  layout->addLayout(buttonBox);
  relativisticCheck = new QCheckBox(tr("Relativistic kinematics"));
  relativisticCheck->setObjectName("relativisticCheck");
  relativisticCheck->setToolTip(tr("Channel wave numbers, Sommerfeld parameters and lab <-> c.m. energies from the invariant mass "
                                   "(the kinematics of the LANL code EDA, ENDF KRL=1).  Level energies stay excitation energies; "
                                   "relativistic.out lists them in EDA's E(s) convention.  Saved as relativistic=1 in <config>."));
  layout->addWidget(relativisticCheck);
  layout->addWidget(abundancesBox);
  setLayout(layout);
}

AbundancesModel *PairsTab::getAbundancesModel() {
  return abundancesModel;
}

void PairsTab::updateAbundancesBox() {
  abundancesBox->setVisible(pairsModel->nuclei().size() > 1);
}

bool PairsTab::readAbundances(QTextStream &inStream) {
  QString line("");
  while (!inStream.atEnd()) {
    line = inStream.readLine();
    if (line.trimmed() == QString("</abundances>")) return true;
    QString entry = line.section('#', 0, 0).trimmed();
    if (entry.isEmpty()) continue;
    QStringList t = entry.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
    // optional "<sample>:" prefix
    QString sample;
    if (!t.isEmpty() && t.first().endsWith(':')) {
      sample = t.first().chopped(1);
      t.removeFirst();
    }
    if (t.size() >= 2 && t.at(0) == "sum") {
      abundancesModel->setSampleSumSigma(sample, t.at(1).toDouble());
      continue;
    }
    if (t.size() < 5) return false;
    Nucleus nucleus;
    nucleus.Z = t.at(0).toInt();
    nucleus.A = t.at(1).toInt();
    if (!abundancesModel->setEntry(sample, nucleus, t.at(2).toDouble(), t.at(3).toInt() == 1, t.at(4).toDouble())) return false;
  }
  return false;
}

bool PairsTab::writeAbundances(QTextStream &outStream) {
  QStringList samples = abundancesModel->samples();
  for (int s = 0; s < samples.size(); s++) {
    const QString prefix = samples.at(s).isEmpty() ? QString() : samples.at(s) + ": ";
    QList<AbundanceEntry> entries = abundancesModel->entries(samples.at(s));
    for (int i = 0; i < entries.size(); i++) {
      const AbundanceEntry &e = entries.at(i);
      outStream << prefix << e.nucleus.Z << " " << e.nucleus.A << " " << e.value << " " << (e.vary ? 1 : 0) << " " << e.sigma
                << "   # " << PairsModel::nucleusLabel(e.nucleus, false) << Qt::endl;
    }
    double sigma = abundancesModel->sampleSumSigma(samples.at(s));
    if (sigma >= 0.0) outStream << prefix << "sum " << sigma << "   # this sample's sum-to-one prior sigma" << Qt::endl;
  }
  return true;
}

void PairsTab::updateSampleCombo() {
  QStringList samples = abundancesModel->samples();
  QString current = abundancesModel->currentSample();
  sampleCombo->blockSignals(true);
  sampleCombo->clear();
  for (int s = 0; s < samples.size(); s++) {
    sampleCombo->addItem(samples.at(s).isEmpty() ? tr("(default sample)") : samples.at(s));
    sampleCombo->setItemData(s, samples.at(s), Qt::UserRole);
  }
  sampleCombo->setCurrentIndex(std::max(0, samples.indexOf(current)));
  sampleCombo->blockSignals(false);
  removeSampleButton->setEnabled(!current.isEmpty());
  sampleSumSigmaSpin->blockSignals(true);
  sampleSumSigmaSpin->setValue(abundancesModel->sampleSumSigma(current));
  sampleSumSigmaSpin->blockSignals(false);
}

void PairsTab::sampleChosen(int index) {
  if (index < 0) return;
  abundancesModel->setCurrentSample(sampleCombo->itemData(index, Qt::UserRole).toString());
  updateSampleCombo();
}

void PairsTab::addSample() {
  bool ok = false;
  QString name = QInputDialog::getText(this, tr("New sample"), tr("Name of the sample (one word, e.g. enr22Ne):"),
                                       QLineEdit::Normal, QString(), &ok).trimmed();
  if (!ok || name.isEmpty()) return;
  if (!abundancesModel->addSample(name)) {
    QMessageBox::warning(this, tr("New sample"), tr("\"%1\" is not a new one-word name without ':' or '#'.").arg(name));
    return;
  }
  abundancesModel->setCurrentSample(name);
  updateSampleCombo();
}

void PairsTab::removeSample() {
  QString current = abundancesModel->currentSample();
  if (current.isEmpty()) return;
  if (QMessageBox::question(this, tr("Remove sample"), tr("Remove sample \"%1\" and its abundances?  Segments that name it "
                                                          "will be refused by the engine until they are changed.").arg(current)) != QMessageBox::Yes)
    return;
  abundancesModel->removeSample(current);
  updateSampleCombo();
}

void PairsTab::sampleSumSigmaEdited(double sigma) {
  abundancesModel->setSampleSumSigma(abundancesModel->currentSample(), sigma);
}

bool PairsTab::hasAbundanceBlock() const {
  return abundancesModel->hasBlock();
}

double PairsTab::abundanceSumSigma() const {
  return sumSigmaSpin->value();
}

void PairsTab::setAbundanceSumSigma(double sigma) {
  sumSigmaSpin->setValue(sigma);
}

bool PairsTab::relativistic() const {
  return relativisticCheck->isChecked();
}

void PairsTab::setRelativistic(bool on) {
  relativisticCheck->setChecked(on);
}

void PairsTab::reset() {
  abundancesModel->clear();
  sumSigmaSpin->setValue(0.0);
  relativisticCheck->setChecked(false);
}

PairsModel *PairsTab::getPairsModel() {
  return pairsModel;
}

/*!
 * The compound-nucleus column only matters when the pairs form more than one
 * nucleus; a single-nucleus project looks as it always did.
 */
void PairsTab::updateNucleusColumn() {
  pairsView->setColumnHidden(PairsModel::NucleusColumn, pairsModel->nuclei().size() < 2);
}

void PairsTab::addPair() {
  AddPairDialog aDialog;
  if (aDialog.exec()) {
    if (!((pairsModel->getPairs()).size() == 0 && aDialog.pairTypeCombo->currentIndex() != 0)) {
      PairsData newPair;
      newPair.lightJ = (aDialog.lightJText->text()).toDouble();
      if (aDialog.lightPiCombo->currentIndex() == 0)
        newPair.lightPi = -1;
      else
        newPair.lightPi = 1;
      newPair.lightZ = (aDialog.lightZText->text()).toInt();
      newPair.lightM = (aDialog.lightMText->text()).toDouble();
      newPair.lightG = 0.;  //(aDialog.lightGText->text()).toDouble();
      newPair.heavyJ = (aDialog.heavyJText->text()).toDouble();
      if (aDialog.heavyPiCombo->currentIndex() == 0)
        newPair.heavyPi = -1;
      else
        newPair.heavyPi = 1;
      newPair.heavyZ = (aDialog.heavyZText->text()).toInt();
      newPair.heavyM = (aDialog.heavyMText->text()).toDouble();
      newPair.heavyG = 0.;  //(aDialog.heavyGText->text()).toDouble();
      newPair.seperationEnergy = (aDialog.seperationEnergyText->text()).toDouble();
      newPair.excitationEnergy = (aDialog.excitationEnergyText->text()).toDouble();
      newPair.channelRadius = (aDialog.channelRadiusText->text()).toDouble();
      if (aDialog.pairTypeCombo->currentIndex() == 1)
        newPair.pairType = 10;
      else if (aDialog.pairTypeCombo->currentIndex() == 2)
        newPair.pairType = 20;
      else
        newPair.pairType = 0;
      unsigned int newMask = 0;
      if (aDialog.e1Check->isChecked()) newMask |= (1 << 0);
      // if(aDialog.m1Check->isChecked()) newMask |= (1<<1);
      if (aDialog.e2Check->isChecked()) newMask |= (1 << 2);
      newPair.ecMultMask = newMask;
      addPair(newPair, pairsModel->numPairs(), false);
    } else {
      QMessageBox::information(this, tr("Pair Type Error"),
                               tr("The first pair must be a particle,particle pair."));
    }
  }
}

void PairsTab::addPair(PairsData pair, int pairIndex, bool fromFile) {
  if (pairsModel->isPair(pair) == -1) {
    pairsModel->insertRows(pairsModel->numPairs(), pairIndex + 1 - pairsModel->numPairs(), QModelIndex());

    QModelIndex index = pairsModel->index(pairIndex, 0, QModelIndex());
    pairsModel->setData(index, pair.lightJ, Qt::EditRole);
    index = pairsModel->index(pairIndex, 1, QModelIndex());
    pairsModel->setData(index, pair.lightPi, Qt::EditRole);
    index = pairsModel->index(pairIndex, 2, QModelIndex());
    pairsModel->setData(index, pair.lightZ, Qt::EditRole);
    index = pairsModel->index(pairIndex, 3, QModelIndex());
    pairsModel->setData(index, pair.lightM, Qt::EditRole);
    index = pairsModel->index(pairIndex, 4, QModelIndex());
    pairsModel->setData(index, pair.lightG, Qt::EditRole);
    index = pairsModel->index(pairIndex, 5, QModelIndex());
    pairsModel->setData(index, pair.heavyJ, Qt::EditRole);
    index = pairsModel->index(pairIndex, 6, QModelIndex());
    pairsModel->setData(index, pair.heavyPi, Qt::EditRole);
    index = pairsModel->index(pairIndex, 7, QModelIndex());
    pairsModel->setData(index, pair.heavyZ, Qt::EditRole);
    index = pairsModel->index(pairIndex, 8, QModelIndex());
    pairsModel->setData(index, pair.heavyM, Qt::EditRole);
    index = pairsModel->index(pairIndex, 9, QModelIndex());
    pairsModel->setData(index, pair.heavyG, Qt::EditRole);
    index = pairsModel->index(pairIndex, 10, QModelIndex());
    pairsModel->setData(index, pair.excitationEnergy, Qt::EditRole);
    index = pairsModel->index(pairIndex, 11, QModelIndex());
    pairsModel->setData(index, pair.seperationEnergy, Qt::EditRole);
    index = pairsModel->index(pairIndex, 12, QModelIndex());
    pairsModel->setData(index, pair.channelRadius, Qt::EditRole);
    index = pairsModel->index(pairIndex, 13, QModelIndex());
    pairsModel->setData(index, pair.pairType, Qt::EditRole);
    index = pairsModel->index(pairIndex, 14, QModelIndex());
    pairsModel->setData(index, pair.ecMultMask, Qt::EditRole);

    pairsView->resizeRowsToContents();
    if (!fromFile) emit(pairAdded(pairIndex));
  } else {
    QMessageBox::information(this, tr("Duplicate Pair"),
                             tr("This pair already exists."));
  }
}

void PairsTab::removePair() {
  QItemSelectionModel *selectionModel = pairsView->selectionModel();
  QModelIndexList indexes = selectionModel->selectedRows();
  QModelIndex index = indexes[0];
  bool previousIsGamma = false;
  if (index.row() == 0 && (pairsModel->getPairs()).size() != 1) {
    QModelIndex previousIndex = pairsModel->index(index.row() + 1, 13, QModelIndex());
    QVariant previousPairType = pairsModel->data(previousIndex, Qt::EditRole);
    if (previousPairType.toInt() == 10) previousIsGamma = true;
  }
  if (!previousIsGamma) {
    pairsModel->removeRows(index.row(), 1, QModelIndex());
    emit(pairRemoved(index.row()));
  } else {
    QMessageBox::information(this, tr("Entrance Channel Error"),
                             tr("This delete will result in a forbidden particle,gamma entrance pair."));
  }
}

void PairsTab::editPair() {
  QItemSelectionModel *selectionModel = pairsView->selectionModel();
  QModelIndexList indexes = selectionModel->selectedRows();

  QModelIndex index = indexes[0];
  QModelIndex i = pairsModel->index(index.row(), 0, QModelIndex());
  QVariant var = pairsModel->data(i, Qt::EditRole);
  QString lightJ = var.toString();
  i = pairsModel->index(index.row(), 1, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  int lightPi = var.toInt();
  i = pairsModel->index(index.row(), 2, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString lightZ = var.toString();
  i = pairsModel->index(index.row(), 3, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString lightM = var.toString();
  i = pairsModel->index(index.row(), 4, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString lightG = var.toString();
  i = pairsModel->index(index.row(), 5, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString heavyJ = var.toString();
  i = pairsModel->index(index.row(), 6, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  int heavyPi = var.toInt();
  i = pairsModel->index(index.row(), 7, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString heavyZ = var.toString();
  i = pairsModel->index(index.row(), 8, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString heavyM = var.toString();
  i = pairsModel->index(index.row(), 9, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString heavyG = var.toString();
  i = pairsModel->index(index.row(), 10, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString excitationEnergy = var.toString();
  i = pairsModel->index(index.row(), 11, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString seperationEnergy = var.toString();
  i = pairsModel->index(index.row(), 12, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString channelRadius = var.toString();
  i = pairsModel->index(index.row(), 13, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  int pairType = var.toInt();
  i = pairsModel->index(index.row(), 14, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  int ecMultMask = var.toInt();


  AddPairDialog aDialog;
  aDialog.setWindowTitle(tr("Edit a Particle Pair"));
  aDialog.lightJText->setText(lightJ);
  if (lightPi == -1)
    aDialog.lightPiCombo->setCurrentIndex(0);
  else
    aDialog.lightPiCombo->setCurrentIndex(1);
  aDialog.lightZText->setText(lightZ);
  aDialog.lightMText->setText(lightM);
  // aDialog.lightGText->setText(lightG);
  aDialog.heavyJText->setText(heavyJ);
  if (heavyPi == -1)
    aDialog.heavyPiCombo->setCurrentIndex(0);
  else
    aDialog.heavyPiCombo->setCurrentIndex(1);
  aDialog.heavyZText->setText(heavyZ);
  aDialog.heavyMText->setText(heavyM);
  // aDialog.heavyGText->setText(heavyG);
  aDialog.excitationEnergyText->setText(excitationEnergy);
  aDialog.seperationEnergyText->setText(seperationEnergy);
  aDialog.channelRadiusText->setText(channelRadius);
  if (pairType == 10)
    aDialog.pairTypeCombo->setCurrentIndex(1);
  else if (pairType == 20)
    aDialog.pairTypeCombo->setCurrentIndex(2);
  else
    aDialog.pairTypeCombo->setCurrentIndex(0);
  if (ecMultMask & (1 << 0))
    aDialog.e1Check->setChecked(true);
  else
    aDialog.e1Check->setChecked(false);
  // if(ecMultMask&(1<<1)) aDialog.m1Check->setChecked(true);
  // else aDialog.m1Check->setChecked(false);
  if (ecMultMask & (1 << 2))
    aDialog.e2Check->setChecked(true);
  else
    aDialog.e2Check->setChecked(false);

  if (aDialog.exec()) {
    if (!(index.row() == 0 && aDialog.pairTypeCombo->currentIndex() != 0)) {
      PairsData pair;
      pair.lightJ = aDialog.lightJText->text().toDouble();
      if (aDialog.lightPiCombo->currentIndex() == 0)
        pair.lightPi = -1;
      else
        pair.lightPi = 1;
      pair.lightZ = aDialog.lightZText->text().toInt();
      pair.lightM = aDialog.lightMText->text().toDouble();
      pair.lightG = 0.;  // aDialog.lightGText->text().toDouble();
      pair.heavyJ = aDialog.heavyJText->text().toDouble();
      if (aDialog.heavyPiCombo->currentIndex() == 0)
        pair.heavyPi = -1;
      else
        pair.heavyPi = 1;
      pair.heavyZ = aDialog.heavyZText->text().toInt();
      pair.heavyM = aDialog.heavyMText->text().toDouble();
      pair.heavyG = 0.;  // aDialog.heavyGText->text().toDouble();
      pair.excitationEnergy = aDialog.excitationEnergyText->text().toDouble();
      pair.seperationEnergy = aDialog.seperationEnergyText->text().toDouble();
      pair.channelRadius = aDialog.channelRadiusText->text().toDouble();
      if (aDialog.pairTypeCombo->currentIndex() == 1)
        pair.pairType = 10;
      else if (aDialog.pairTypeCombo->currentIndex() == 2)
        pair.pairType = 20;
      else
        pair.pairType = 0;
      unsigned char newECMultMask = 0;
      if (aDialog.e1Check->isChecked()) newECMultMask |= (1 << 0);
      // if(aDialog.m1Check->isChecked()) newECMultMask |= (1<<1);
      if (aDialog.e2Check->isChecked()) newECMultMask |= (1 << 2);
      pair.ecMultMask = newECMultMask;
      editPair(pair, index.row(), false);
    } else {
      QMessageBox::information(this, tr("Pair Type Error"),
                               tr("The first pair must be a particle,particle pair."));
    }
  }
}

void PairsTab::editPair(PairsData pair, int pairIndex, bool fromFile) {
  QModelIndex i = pairsModel->index(pairIndex, 0, QModelIndex());
  QVariant var = pairsModel->data(i, Qt::EditRole);
  if (pair.lightJ != var.toDouble()) pairsModel->setData(i, pair.lightJ, Qt::EditRole);
  i = pairsModel->index(pairIndex, 1, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.lightPi != var.toInt()) pairsModel->setData(i, pair.lightPi, Qt::EditRole);
  i = pairsModel->index(pairIndex, 2, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.lightZ != var.toInt()) pairsModel->setData(i, pair.lightZ, Qt::EditRole);
  i = pairsModel->index(pairIndex, 3, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.lightM != var.toDouble()) pairsModel->setData(i, pair.lightM, Qt::EditRole);
  i = pairsModel->index(pairIndex, 4, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.lightG != var.toDouble()) pairsModel->setData(i, pair.lightG, Qt::EditRole);
  i = pairsModel->index(pairIndex, 5, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.heavyJ != var.toDouble()) pairsModel->setData(i, pair.heavyJ, Qt::EditRole);
  i = pairsModel->index(pairIndex, 6, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.heavyPi != var.toInt()) pairsModel->setData(i, pair.heavyPi, Qt::EditRole);
  i = pairsModel->index(pairIndex, 7, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.heavyZ != var.toInt()) pairsModel->setData(i, pair.heavyZ, Qt::EditRole);
  i = pairsModel->index(pairIndex, 8, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.heavyM != var.toDouble()) pairsModel->setData(i, pair.heavyM, Qt::EditRole);
  i = pairsModel->index(pairIndex, 9, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.heavyG != var.toDouble()) pairsModel->setData(i, pair.heavyG, Qt::EditRole);
  i = pairsModel->index(pairIndex, 10, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.excitationEnergy != var.toDouble()) pairsModel->setData(i, pair.excitationEnergy, Qt::EditRole);
  i = pairsModel->index(pairIndex, 11, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.seperationEnergy != var.toDouble()) pairsModel->setData(i, pair.seperationEnergy, Qt::EditRole);
  i = pairsModel->index(pairIndex, 12, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.channelRadius != var.toDouble()) pairsModel->setData(i, pair.channelRadius, Qt::EditRole);
  i = pairsModel->index(pairIndex, 13, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.pairType != var.toInt()) pairsModel->setData(i, pair.pairType, Qt::EditRole);
  i = pairsModel->index(pairIndex, 14, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.pairType != var.toInt()) pairsModel->setData(i, pair.ecMultMask, Qt::EditRole);

  if (!fromFile) emit(pairEdited(pairIndex));
}

void PairsTab::updateButtons(const QItemSelection &selection) {
  QModelIndexList indexes = selection.indexes();

  if (!indexes.isEmpty()) {
    deleteButton->setEnabled(true);
  } else {
    deleteButton->setEnabled(false);
  }
}

bool PairsTab::parseOldECSection(QTextStream &inStream) {
  int isActive;
  int exitPairIndex;
  double minJ;
  double maxJ;
  int multMask;
  QString line("");
  while (line.trimmed() != QString("</externalCapture>") && !inStream.atEnd()) {
    line = inStream.readLine();
    if (line.trimmed().isEmpty()) continue;
    if (line.trimmed() != QString("</externalCapture>") && !inStream.atEnd()) {
      QTextStream in(&line);
      in >> isActive >> exitPairIndex >> minJ >> maxJ >> multMask;
      if (in.status() != QTextStream::Ok) return false;
      if (exitPairIndex - 1 < pairsModel->getPairs().size()) {
        QModelIndex i = pairsModel->index(exitPairIndex - 1, 14, QModelIndex());
        QVariant var = pairsModel->data(i, Qt::EditRole);
        pairsModel->setData(i, multMask, Qt::EditRole);
      }
    }
  }
  if (line.trimmed() != QString("</externalCapture>")) return false;
  return true;
}

void PairsTab::showInfo(int which, QString title) {
  if (which < infoText.size()) {
    if (!infoDialog[which]) {
      infoDialog[which] = new InfoDialog(infoText[which], this, title);
      infoDialog[which]->setAttribute(Qt::WA_DeleteOnClose);
      infoDialog[which]->show();
    } else
      infoDialog[which]->raise();
  }
}
