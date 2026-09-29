#include "ThmExperimentsPage.h"

#include <QCheckBox>
#include <QComboBox>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QItemSelectionModel>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QPushButton>
#include <QSet>
#include <QTableWidget>
#include <QTextDocument>
#include <QVBoxLayout>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include "PairsModel.h"
#include "SegmentsDataModel.h"
#include "DataLine.h"
#include "ThmExperiment.h"
#include "ThmLineshape.h"

namespace {

// MeV/u, the value EData::SetupThmExperiments converts the binding with.
const double kAmu = 931.49410242;

// A whole token as a number, as the engine reads Ebeam.
bool readWholeDouble(const QString &text, double &x) {
  std::istringstream s(text.toStdString());
  std::string rest;
  return !!(s >> x) && !(s >> rest) && std::isfinite(x);
}

QString plain(const QString &html) {
  QTextDocument d;
  d.setHtml(html);
  return d.toPlainText().simplified();
}

const char *kDocs =
    "docs/source/theory/thm_implementation.rst, section \"THM experiments\"; online: "
    "<a href=\"https://rdeboer1.github.io/AZURE2/\">rdeboer1.github.io/AZURE2</a> "
    "(Theory &gt; Trojan Horse (HOES) Observable).";

}  // namespace

QStringList ThmExperimentsPage::nuclideNames() {
  // The engine's table, looked up rather than copied (ThmNuclide::Find).
  QStringList names;
  for (int Z = 0; Z <= 118; Z++)
    for (int A = std::max(Z, 1); A <= 300; A++)
      if (const ThmNuclide *n = ThmNuclide::Find(Z, A)) names << QString::fromStdString(n->name);
  return names;
}

ThmExperimentsPage::ThmExperimentsPage(const QStringList &experimentLines, SegmentsDataModel *segments,
                                       PairsModel *pairs, const QString &projectDir, bool brune, QWidget *parent) :
  QWidget(parent),
  segments_(segments),
  pairs_(pairs),
  projectDir_(projectDir),
  brune_(brune),
  oldLines_(experimentLines) {
  oldRecords_ = ThmExperimentRecord::read(experimentLines);
  records_ = oldRecords_;

  experimentTable = new QTableWidget(0, 4);
  experimentTable->setHorizontalHeaderLabels(QStringList() << tr("Name") << tr("Segments") << tr("Background")
                                                           << tr("Three-body reaction"));
  experimentTable->horizontalHeader()->setStretchLastSection(true);
  experimentTable->verticalHeader()->hide();
  experimentTable->setSelectionBehavior(QAbstractItemView::SelectRows);
  experimentTable->setSelectionMode(QAbstractItemView::SingleSelection);
  experimentTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
  experimentTable->setToolTip(
      tr("experiment[<name>]: THM data segments measured together (exit channels, angular bins or runs "
         "of one three-body reaction) share one profiled normalization and, optionally, a background."));
  addButton = new QPushButton(tr("Add"));
  removeButton = new QPushButton(tr("Remove"));
  connect(addButton, &QPushButton::clicked, this, [this]() { addExperiment(); });
  connect(removeButton, &QPushButton::clicked, this, [this]() { removeCurrent(); });
  connect(experimentTable->selectionModel(), SIGNAL(selectionChanged(const QItemSelection &, const QItemSelection &)),
          this, SLOT(tableSelectionChanged()));

  QGroupBox *listBox = new QGroupBox(tr("THM experiments"));
  QGridLayout *ll = new QGridLayout;
  ll->addWidget(experimentTable, 0, 0, 1, 3);
  ll->addWidget(addButton, 1, 1);
  ll->addWidget(removeButton, 1, 2);
  ll->setColumnStretch(0, 1);
  listBox->setLayout(ll);

  nameEdit = new QLineEdit;
  nameEdit->setToolTip(tr("The experiment's name: letters, digits and _ - . + only."));
  connect(nameEdit, SIGNAL(textEdited(const QString &)), this, SLOT(nameEdited(const QString &)));
  segmentList = new QListWidget;
  segmentList->setToolTip(
      tr("segments=: the data segments of the experiment (numbers count every line of the Data segments "
         "table, active or not). Offered: THM segments (isDiff >= 10) with a free norm that are in no other "
         "experiment; all of them share one norm."));
  connect(segmentList, SIGNAL(itemChanged(QListWidgetItem *)), this, SLOT(segmentItemChanged(QListWidgetItem *)));
  backgroundCombo = new QComboBox;
  backgroundCombo->addItems(QStringList() << "none" << "const" << "linear" << "quadratic");
  backgroundCombo->setToolTip(
      tr("background=: a smooth b0 + b1 E + b2 E^2 (E the c.m. energy of the THM entrance pair) added to the "
         "folded model, profiled with the norm by linear least squares; it may come out negative."));
  connect(backgroundCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(backgroundChanged(int)));

  const QStringList nuclides = nuclideNames();
  auto nuclideCombo = [&](const QString &tip) {
    QComboBox *c = new QComboBox;
    c->setEditable(true);
    c->addItems(nuclides);
    c->setCurrentIndex(-1);
    c->setInsertPolicy(QComboBox::NoInsert);
    c->setToolTip(tip + tr(" A nuclide of the built-in table (AME2020 nuclear masses) or Z,A,mass "
                           "(nuclear mass in u)."));
    connect(c, SIGNAL(editTextChanged(const QString &)), this, SLOT(kinematicsEdited()));
    return c;
  };
  beamCombo = nuclideCombo(tr("beam=: the projectile of the three-body reaction."));
  targetCombo = nuclideCombo(tr("target=: the target of the three-body reaction."));
  spectatorCombo = nuclideCombo(tr("spectator=: the spectator s of the Trojan horse a = x + s."));
  beamEnergyEdit = new QLineEdit;
  beamEnergyEdit->setToolTip(tr("Ebeam=: the lab beam energy, MeV (> 0)."));
  connect(beamEnergyEdit, SIGNAL(textEdited(const QString &)), this, SLOT(kinematicsEdited()));
  lineshapeCheck = new QCheckBox(tr("Coulomb line shape of the spectator"));
  lineshapeCheck->setToolTip(
      tr("lineshape=on: the final-state Coulomb interaction of the charged spectator with the resonance and its "
         "decay products skews and shifts each resonance (the factor N_C per level, inside the coherent level sum; "
         "Mukhamedzhanov, Kadyrov & Pang, EPJA 56 (2020) 233; Mukhamedzhanov, EPJA 58 (2022) 71). Needs the "
         "three-body reaction and the Brune parameterization. zeta is shown at the ends of the data."));
  connect(lineshapeCheck, SIGNAL(toggled(bool)), this, SLOT(lineshapeToggled(bool)));
  derivedLabel = new QLabel;
  derivedLabel->setWordWrap(true);
  derivedLabel->setTextInteractionFlags(Qt::TextSelectableByMouse);
  derivedLabel->setMinimumHeight(3 * derivedLabel->fontMetrics().lineSpacing());  // three lines, always

  kinematicsBox = new QGroupBox(tr("Three-body reaction"));
  kinematicsBox->setCheckable(true);
  kinematicsBox->setChecked(false);
  kinematicsBox->setToolTip(
      tr("beam, target, spectator and Ebeam go together (all four or none). One of beam and target must be a "
         "nucleus of the segments' entrance pair and the other the second nucleus plus the spectator. AZURE2 "
         "reports B(x+s) and the quasi-free energy; the Coulomb line shape uses them."));
  connect(kinematicsBox, SIGNAL(toggled(bool)), this, SLOT(kinematicsEdited()));
  QGridLayout *kl = new QGridLayout;
  kl->addWidget(new QLabel(tr("Beam:")), 0, 0, Qt::AlignRight);
  kl->addWidget(beamCombo, 0, 1);
  kl->addWidget(new QLabel(tr("Target:")), 0, 2, Qt::AlignRight);
  kl->addWidget(targetCombo, 0, 3);
  kl->addWidget(new QLabel(tr("Spectator:")), 1, 0, Qt::AlignRight);
  kl->addWidget(spectatorCombo, 1, 1);
  kl->addWidget(new QLabel(tr("Beam energy (lab, MeV):")), 1, 2, Qt::AlignRight);
  kl->addWidget(beamEnergyEdit, 1, 3);
  kl->addWidget(lineshapeCheck, 2, 0, 1, 4);
  kl->addWidget(derivedLabel, 3, 0, 1, 4);
  kl->setColumnStretch(1, 1);
  kl->setColumnStretch(3, 1);
  kinematicsBox->setLayout(kl);

  editorBox = new QGroupBox(tr("Selected experiment"));
  QGridLayout *el = new QGridLayout;
  el->addWidget(new QLabel(tr("Name:")), 0, 0, Qt::AlignRight);
  el->addWidget(nameEdit, 0, 1);
  el->addWidget(new QLabel(tr("Segments:")), 1, 0, Qt::AlignRight | Qt::AlignTop);
  el->addWidget(segmentList, 1, 1);
  el->addWidget(new QLabel(tr("Background:")), 2, 0, Qt::AlignRight);
  el->addWidget(backgroundCombo, 2, 1);
  el->addWidget(kinematicsBox, 3, 0, 1, 2);
  el->setColumnStretch(1, 1);
  editorBox->setLayout(el);

  QLabel *docs = new QLabel(tr("The segments of an experiment share one profiled norm; segments in no "
                               "experiment keep their own. Documentation: ") +
                            QString(kDocs));
  docs->setWordWrap(true);
  docs->setOpenExternalLinks(true);
  docs->setTextFormat(Qt::RichText);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(listBox);
  mainLayout->addWidget(editorBox, 1);
  mainLayout->addWidget(docs);
  setLayout(mainLayout);

  refreshTable();
  selectExperiment(records_.isEmpty() ? -1 : 0);
}

QStringList ThmExperimentsPage::experimentLines() const {
  return ThmExperimentRecord::compose(oldLines_, oldRecords_, records_);
}

// ---------------------------------------------------------------------------
// Table and selection

void ThmExperimentsPage::refreshRow(int row) {
  const ThmExperimentRecord &r = records_.at(row);
  QString reaction = QString::fromUtf8("—");
  if (r.hasKinematics())
    reaction = tr("%1 + %2 at %3 MeV, spectator %4")
                   .arg(r.beam.isEmpty() ? "?" : r.beam, r.target.isEmpty() ? "?" : r.target,
                        r.beamEnergy.isEmpty() ? "?" : r.beamEnergy, r.spectator.isEmpty() ? "?" : r.spectator) +
               (r.lineshape ? tr(", line shape") : QString());
  const QString cells[4] = {r.name, ThmExperimentRecord::segmentsListText(r.segments), r.background, reaction};
  for (int c = 0; c < 4; c++) {
    QTableWidgetItem *item = experimentTable->item(row, c);
    if (!item) {
      item = new QTableWidgetItem;
      experimentTable->setItem(row, c, item);
    }
    item->setText(cells[c]);
    item->setToolTip(r.extraTokens.isEmpty() ? QString()
                                             : tr("Also kept as written: %1").arg(r.extraTokens.join(' ')));
  }
}

void ThmExperimentsPage::refreshTable() {
  experimentTable->blockSignals(true);
  experimentTable->setRowCount(records_.size());
  for (int row = 0; row < records_.size(); row++) refreshRow(row);
  experimentTable->blockSignals(false);
}

void ThmExperimentsPage::selectExperiment(int row) {
  current_ = (row >= 0 && row < records_.size()) ? row : -1;
  loading_ = true;
  if (current_ >= 0)
    experimentTable->selectRow(current_);
  else
    experimentTable->clearSelection();
  loading_ = false;
  loadEditor();
}

void ThmExperimentsPage::tableSelectionChanged() {
  if (loading_) return;
  QList<QModelIndex> rows = experimentTable->selectionModel()->selectedRows();
  current_ = rows.isEmpty() ? -1 : rows.first().row();
  loadEditor();
}

void ThmExperimentsPage::addExperiment() {
  QSet<QString> names;
  for (const ThmExperimentRecord &r : records_) names.insert(r.name);
  int n = 1;
  while (names.contains(QString("E%1").arg(n))) n++;
  ThmExperimentRecord r;
  r.name = QString("E%1").arg(n);
  records_ << r;
  refreshTable();
  selectExperiment(records_.size() - 1);
}

void ThmExperimentsPage::removeCurrent() {
  if (current_ < 0) return;
  int row = current_;
  records_.removeAt(row);
  refreshTable();
  selectExperiment(std::min(row, (int)records_.size() - 1));
}

// ---------------------------------------------------------------------------
// Editor

QString ThmExperimentsPage::segmentLabel(int key) const {
  const QList<SegmentsDataData> lines = segments_->getLines();
  if (key < 1 || key > lines.size()) return tr("%1: no such line in the Data segments table").arg(key);
  const SegmentsDataData &s = lines.at(key - 1);
  const QList<PairsData> pairs = pairs_->getPairs();
  auto pairText = [&](int pairKey) {
    if (pairKey < 1 || pairKey > pairs.size()) return QString("?");
    const PairsData &p = pairs.at(pairKey - 1);
    return plain(pairs_->getParticleLabel(p, 0)) + "+" + plain(pairs_->getParticleLabel(p, 1));
  };
  QString text = QString("%1: %2 -> %3, %4-%5 MeV, %6")
                     .arg(key)
                     .arg(pairText(s.entrancePairIndex), pairText(s.exitPairIndex))
                     .arg(s.lowEnergy)
                     .arg(s.highEnergy)
                     .arg(s.dataFile);
  QStringList notes;
  if (!s.isActive) notes << tr("inactive");
  if (!s.isTHM) notes << tr("not THM");
  if (!s.varyNorm) notes << tr("fixed norm");
  if (!notes.isEmpty()) text += " (" + notes.join(", ") + ")";
  return text;
}

void ThmExperimentsPage::fillSegmentList() {
  segmentList->blockSignals(true);
  segmentList->clear();
  if (current_ >= 0) {
    const ThmExperimentRecord &r = records_.at(current_);
    QSet<int> others;
    for (int i = 0; i < records_.size(); i++)
      if (i != current_)
        for (int k : records_.at(i).segments) others.insert(k);
    const QList<SegmentsDataData> lines = segments_->getLines();
    QList<int> keys;
    for (int k = 1; k <= lines.size(); k++)
      if ((lines.at(k - 1).isTHM && lines.at(k - 1).varyNorm && !others.contains(k)) || r.segments.contains(k))
        keys << k;
    for (int k : r.segments)
      if (!keys.contains(k)) keys << k;  // beyond the table: shown so that it can be unticked
    for (int k : keys) {
      QListWidgetItem *item = new QListWidgetItem(segmentLabel(k), segmentList);
      item->setData(Qt::UserRole, k);
      item->setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
      item->setCheckState(r.segments.contains(k) ? Qt::Checked : Qt::Unchecked);
    }
  }
  segmentList->blockSignals(false);
}

void ThmExperimentsPage::loadEditor() {
  loading_ = true;
  editorBox->setEnabled(current_ >= 0);
  removeButton->setEnabled(current_ >= 0);
  ThmExperimentRecord r = current_ >= 0 ? records_.at(current_) : ThmExperimentRecord();
  nameEdit->setText(r.name);
  fillSegmentList();
  int bg = backgroundCombo->findText(r.background);
  backgroundCombo->setCurrentIndex(bg >= 0 ? bg : 0);
  kinematicsBox->setChecked(r.hasKinematics());
  beamCombo->setEditText(r.beam);
  targetCombo->setEditText(r.target);
  spectatorCombo->setEditText(r.spectator);
  beamEnergyEdit->setText(r.beamEnergy);
  lineshapeCheck->setChecked(r.lineshape);
  lineshapeCheck->setEnabled(!r.beam.isEmpty() && !r.target.isEmpty() && !r.spectator.isEmpty() &&
                             !r.beamEnergy.isEmpty());
  showDerived(r);
  loading_ = false;
}

void ThmExperimentsPage::showDerived(const ThmExperimentRecord &r) {
  QString why, info = derivedInfo(r, &why);
  derivedLabel->setText(info.isEmpty() ? why : why.isEmpty() ? info : info + "\n" + why);
}

void ThmExperimentsPage::nameEdited(const QString &text) {
  if (loading_ || current_ < 0) return;
  records_[current_].name = text.trimmed();
  refreshRow(current_);
}

void ThmExperimentsPage::setSegmentsOfCurrent(const QList<int> &segments) {
  if (current_ < 0) return;
  storeSegments(segments);
  loadEditor();
}

void ThmExperimentsPage::storeSegments(const QList<int> &segments) {
  ThmExperimentRecord &r = records_[current_];
  QList<int> sorted = segments;
  std::sort(sorted.begin(), sorted.end());
  sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
  QList<int> written;
  if (!ThmExperimentRecord::expandSegments(r.segmentsText, written) || written != sorted) r.segmentsText.clear();
  r.segments = sorted;
  refreshRow(current_);
  showDerived(r);
}

void ThmExperimentsPage::segmentItemChanged(QListWidgetItem *) {
  if (loading_ || current_ < 0) return;
  QList<int> chosen;
  for (int i = 0; i < segmentList->count(); i++)
    if (segmentList->item(i)->checkState() == Qt::Checked) chosen << segmentList->item(i)->data(Qt::UserRole).toInt();
  storeSegments(chosen);  // the list stays as it is while one of its items signals
}

void ThmExperimentsPage::backgroundChanged(int) {
  if (loading_ || current_ < 0) return;
  records_[current_].background = backgroundCombo->currentText();
  refreshRow(current_);
}

void ThmExperimentsPage::kinematicsEdited() {
  if (loading_ || current_ < 0) return;
  ThmExperimentRecord &r = records_[current_];
  const bool on = kinematicsBox->isChecked();
  r.beam = on ? beamCombo->currentText().trimmed() : QString();
  r.target = on ? targetCombo->currentText().trimmed() : QString();
  r.spectator = on ? spectatorCombo->currentText().trimmed() : QString();
  r.beamEnergy = on ? beamEnergyEdit->text().trimmed() : QString();
  // The line shape needs the reaction: offered once all four keys are given,
  // and dropped with them.
  lineshapeCheck->setEnabled(on && !r.beam.isEmpty() && !r.target.isEmpty() && !r.spectator.isEmpty() &&
                             !r.beamEnergy.isEmpty());
  if (!on && lineshapeCheck->isChecked()) {
    lineshapeCheck->blockSignals(true);
    lineshapeCheck->setChecked(false);
    lineshapeCheck->blockSignals(false);
  }
  r.lineshape = on && lineshapeCheck->isChecked();
  refreshRow(current_);
  showDerived(r);
}

void ThmExperimentsPage::lineshapeToggled(bool on) {
  if (loading_ || current_ < 0) return;
  records_[current_].lineshape = on && kinematicsBox->isChecked();
  refreshRow(current_);
  showDerived(records_.at(current_));
}

// ---------------------------------------------------------------------------
// The engine's rules

bool ThmExperimentsPage::reaction(const ThmExperimentRecord &x, Reaction &out, QString *error) const {
  if (error) error->clear();
  if (x.beam.isEmpty() || x.target.isEmpty() || x.spectator.isEmpty() || x.beamEnergy.isEmpty()) return false;
  ThmNuclide b, t, sp;
  std::string why;
  if (!(why = ThmNuclide::Parse(x.beam.toStdString(), b)).empty() ||
      !(why = ThmNuclide::Parse(x.target.toStdString(), t)).empty() ||
      !(why = ThmNuclide::Parse(x.spectator.toStdString(), sp)).empty()) {
    if (error) *error = QString::fromStdString(why);
    return false;
  }
  double ebeam;
  if (!readWholeDouble(x.beamEnergy, ebeam) || !(ebeam > 0.0)) {
    if (error) *error = tr("Ebeam='%1': expected the lab beam energy in MeV, > 0").arg(x.beamEnergy);
    return false;
  }
  // The entrance pair x + A of the segments in use (EData::SetupThmExperiments).
  const QList<SegmentsDataData> lines = segments_->getLines();
  const QList<PairsData> pairs = pairs_->getPairs();
  int pairKey = -1;
  for (int k : x.segments) {
    if (k < 1 || k > lines.size() || !lines.at(k - 1).isActive) continue;
    int e = lines.at(k - 1).entrancePairIndex;
    if (pairKey < 0)
      pairKey = e;
    else if (e != pairKey) {
      if (error)
        *error = tr("beam/target/spectator describe one reaction, but its segments have different entrance pairs.");
      return false;
    }
  }
  if (pairKey < 1 || pairKey > pairs.size()) return false;  // nothing to compare with yet
  const PairsData &pair = pairs.at(pairKey - 1);
  const int Z[2] = {pair.lightZ, pair.heavyZ};
  const int A[2] = {(int)std::lround(pair.lightM), (int)std::lround(pair.heavyM)};
  const double M[2] = {pair.lightM, pair.heavyM};
  int horse = -1, other = -1;
  for (int h = 0; h < 2 && horse < 0; h++) {
    const ThmNuclide &th = h == 0 ? b : t, &tg = h == 0 ? t : b;
    for (int k = 0; k < 2; k++)
      if (tg.Z == Z[k] && tg.A == A[k] && th.Z - sp.Z == Z[1 - k] && th.A - sp.A == A[1 - k]) {
        horse = h;
        other = k;
        break;
      }
  }
  if (horse < 0) {
    if (error)
      *error = tr("beam %1 + target %2 with spectator %3 does not give the entrance pair of its segments (Z,A) = "
                  "(%4,%5) + (%6,%7): one of beam/target must be a nucleus of the pair and the other the second "
                  "nucleus plus the spectator.")
                   .arg(QString::fromStdString(b.name), QString::fromStdString(t.name),
                        QString::fromStdString(sp.name))
                   .arg(Z[0])
                   .arg(A[0])
                   .arg(Z[1])
                   .arg(A[1]);
    return false;
  }
  const ThmNuclide &th = horse == 0 ? b : t, &nA = horse == 0 ? t : b;
  const ThmNuclide *tabX = ThmNuclide::Find(th.Z - sp.Z, th.A - sp.A);
  const double mX = tabX ? tabX->mass : M[1 - other];
  const double mA = nA.mass;
  out.beam = b;
  out.target = t;
  out.spectator = sp;
  out.horse = th;
  out.pairKey = pairKey;
  out.beamEnergy = ebeam;
  out.bind = (mX + sp.mass - th.mass) * kAmu;
  out.exa = horse == 0 ? ebeam * mX / th.mass * mA / (mX + mA) : ebeam * mX / (mA + mX);
  return true;
}

QString ThmExperimentsPage::derivedInfo(const ThmExperimentRecord &x, QString *error) const {
  Reaction r;
  if (!reaction(x, r, error)) return QString();
  QString text =
      tr("B(x+s) = %3 MeV (Trojan horse %1 = x + %2)\nquasi-free E(x+A) = %4 MeV, E_qf = E(x+A) - B = %5 MeV")
          .arg(QString::fromStdString(r.horse.name), QString::fromStdString(r.spectator.name))
          .arg(QString::number(r.bind, 'g', 6))
          .arg(QString::number(r.exa, 'g', 6))
          .arg(QString::number(r.exa - r.bind, 'g', 6));
  if (x.lineshape) {
    const QString shape = lineshapeInfo(x, error);
    if (!shape.isEmpty()) text += "\n" + shape;
  }
  return text;
}

bool ThmExperimentsPage::pointRange(const QList<int> &segments, double &lo, double &hi) const {
  // The c.m. energies of the points AZURE2 reads (ESegment::FillData): rows
  // of the data file inside the segment's lab energy range (and angle range
  // if differential), converted with the entrance pair's masses.
  const QList<SegmentsDataData> lines = segments_->getLines();
  const QList<PairsData> pairs = pairs_->getPairs();
  lo = 1.0e300;
  hi = -1.0e300;
  for (int k : segments) {
    if (k < 1 || k > lines.size() || !lines.at(k - 1).isActive) continue;
    const SegmentsDataData &s = lines.at(k - 1);
    if (s.entrancePairIndex < 1 || s.entrancePairIndex > pairs.size()) return false;
    const PairsData &pair = pairs.at(s.entrancePairIndex - 1);
    const double factor = pair.heavyM / (pair.lightM + pair.heavyM);
    const bool differential = s.dataType == 1 || s.dataType == 4 || s.dataType == 7 || s.dataType == 8;
    QString path = s.dataFile;
    if (QFileInfo(path).isRelative()) path = QDir(projectDir_).filePath(path);
    std::ifstream in(QFile::encodeName(path).constData());
    if (!in) return false;
    while (true) {
      DataLine line(in);
      if (line.atEnd()) break;
      if (!line.valid()) return false;
      if (line.energy() < s.lowEnergy || line.energy() > s.highEnergy) continue;
      if (differential && (line.angle() < s.lowAngle || line.angle() > s.highAngle)) continue;
      lo = std::min(lo, line.energy() * factor);
      hi = std::max(hi, line.energy() * factor);
    }
  }
  return lo <= hi;
}

QString ThmExperimentsPage::lineshapeInfo(const ThmExperimentRecord &x, QString *error) const {
  Reaction r;
  if (!reaction(x, r, error)) return QString();
  if (!brune_) {
    if (error)
      *error = tr("lineshape=on uses the observed level energies and widths as the resonance poles; it needs the "
                  "Brune parameterization.");
    return QString();
  }
  const QList<PairsData> pairs = pairs_->getPairs();
  const PairsData &entrance = pairs.at(r.pairKey - 1);
  ThmLineshape shape;  // as EData::BuildThmGroups sets it up
  shape.Zs = r.spectator.Z;
  shape.ms = r.spectator.mass;
  shape.ZF = entrance.lightZ + entrance.heavyZ;
  shape.mF = entrance.lightM + entrance.heavyM;
  shape.eAA = r.beamEnergy * r.target.mass / (r.beam.mass + r.target.mass);
  shape.bind = r.bind;
  double lo, hi;
  if (!pointRange(x.segments, lo, hi)) return QString();  // the engine reports an unreadable file itself
  if (!(shape.EsF(hi) > 0.0)) {
    if (error)
      *error = tr("lineshape=on: at E = %1 MeV the spectator has no energy left (E_sF = E_aA - B - E = %2 - %3 - "
                  "%1 MeV <= 0); check Ebeam.")
                   .arg(QString::number(hi, 'g', 6), QString::number(shape.eAA, 'g', 6),
                        QString::number(shape.bind, 'g', 6));
    return QString();
  }
  // zeta per exit pair of the segments, at the lowest and highest point energy.
  const QList<SegmentsDataData> lines = segments_->getLines();
  QList<int> exits;
  QStringList parts;
  for (int k : x.segments) {
    if (k < 1 || k > lines.size() || !lines.at(k - 1).isActive) continue;
    const int key = lines.at(k - 1).exitPairIndex;
    if (exits.contains(key) || key < 1 || key > pairs.size()) continue;
    exits << key;
    const PairsData &p = pairs.at(key - 1);
    const bool lightFirst = p.lightM <= p.heavyM;
    const int ZB = lightFirst ? p.heavyZ : p.lightZ;
    const double mB = lightFirst ? p.heavyM : p.lightM;
    parts << tr("%1 into %2")
                 .arg(QString::fromUtf8("%1 \u2026 %2")
                          .arg(QString::number(shape.Zeta(lo, ZB, mB), 'g', 3),
                               QString::number(shape.Zeta(hi, ZB, mB), 'g', 3)),
                      plain(pairs_->getParticleLabel(p, 0)) + "+" + plain(pairs_->getParticleLabel(p, 1)));
  }
  return QString::fromUtf8("\u03b6 = ") + parts.join("; ") +
         tr(" at E = %1 \u2026 %2 MeV (the ends of the data)")
             .arg(QString::number(lo, 'g', 4), QString::number(hi, 'g', 4));
}

QString ThmExperimentsPage::check() const {
  QString why = ThmSettings::checkExperimentLines(experimentLines());
  if (!why.isEmpty()) return why;
  // The startup checks of EData::SetupThmExperiments that need the project.
  const QList<SegmentsDataData> lines = segments_->getLines();
  for (const ThmExperimentRecord &x : records_) {
    const QString where = QString("<thm> experiment[%1]: ").arg(x.name);
    for (int k : x.segments) {
      if (k > lines.size())
        return where + tr("segment %1: <segmentsData> has only %2 line(s).").arg(k).arg(lines.size());
      const SegmentsDataData &s = lines.at(k - 1);
      if (!s.isActive) continue;  // left out with a warning
      if (!s.isTHM) return where + tr("segment %1 is not a THM segment (isDiff < 10).").arg(k);
      if (!s.varyNorm)
        return where + tr("segment %1 has a fixed norm; the segments of an experiment share one free "
                          "(profiled) norm, so free it.")
                           .arg(k);
    }
    QString error;
    derivedInfo(x, &error);
    if (!error.isEmpty()) return where + error;
  }
  return QString();
}
