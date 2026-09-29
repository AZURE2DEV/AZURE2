#include "ThmExperimentsPage.h"

#include <QCheckBox>
#include <QComboBox>
#include <QDir>
#include <QFile>
#include <QFileDialog>
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
#include <QScrollArea>
#include <QSet>
#include <QSignalBlocker>
#include <QSpinBox>
#include <QStyle>
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
#include "ThmNumberSpin.h"

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
  // A compact list: name and segments; the rest is in the row's tooltip.
  experimentTable->setColumnHidden(2, true);
  experimentTable->setColumnHidden(3, true);
  experimentTable->horizontalHeader()->setStretchLastSection(false);
  experimentTable->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
  experimentTable->horizontalHeader()->setSectionResizeMode(1, QHeaderView::ResizeToContents);
  experimentTable->verticalHeader()->hide();
  experimentTable->setSelectionBehavior(QAbstractItemView::SelectRows);
  experimentTable->setSelectionMode(QAbstractItemView::SingleSelection);
  experimentTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
  experimentTable->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  experimentTable->setToolTip(
      tr("experiment[<name>]: THM data segments measured together share one profiled norm and, optionally, a "
         "background. Segments in no experiment keep their own norm."));
  addButton = new QPushButton(tr("Add"));
  removeButton = new QPushButton(tr("Remove"));
  connect(addButton, &QPushButton::clicked, this, [this]() { addExperiment(); });
  connect(removeButton, &QPushButton::clicked, this, [this]() { removeCurrent(); });
  connect(experimentTable->selectionModel(), SIGNAL(selectionChanged(const QItemSelection &, const QItemSelection &)),
          this, SLOT(tableSelectionChanged()));

  QGroupBox *listBox = new QGroupBox(tr("Experiments"));
  listBox->setMaximumWidth(fontMetrics().horizontalAdvance("M") * 12);
  QGridLayout *ll = new QGridLayout;
  ll->addWidget(experimentTable, 0, 0, 1, 2);
  ll->addWidget(addButton, 1, 0);
  ll->addWidget(removeButton, 1, 1);
  listBox->setLayout(ll);

  nameEdit = new QLineEdit;
  nameEdit->setToolTip(tr("Letters, digits and _ - . + only."));
  connect(nameEdit, SIGNAL(textEdited(const QString &)), this, SLOT(nameEdited(const QString &)));
  segmentList = new QListWidget;
  segmentList->setToolTip(
      tr("segments=: THM data segments (isDiff >= 10) with a free norm that are in no other experiment. Numbers "
         "count every line of the Data segments table, active or not."));
  // About three rows: the sections below keep their room.
  segmentList->setMinimumHeight(2 * segmentList->fontMetrics().lineSpacing() + 10);
  segmentList->setMaximumHeight(4 * segmentList->fontMetrics().lineSpacing() + 10);
  segmentList->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Minimum);
  connect(segmentList, SIGNAL(itemChanged(QListWidgetItem *)), this, SLOT(segmentItemChanged(QListWidgetItem *)));
  backgroundCombo = new QComboBox;
  backgroundCombo->addItems(QStringList() << "none" << "const" << "linear" << "quadratic");
  backgroundCombo->setToolTip(
      tr("background=: b0 + b1 E + b2 E^2 (E: c.m. energy of the THM entrance pair) added to the folded model, "
         "profiled with the norm."));
  connect(backgroundCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(backgroundChanged(int)));

  const QStringList nuclides = nuclideNames();
  auto nuclideCombo = [&](const QString &tip) {
    QComboBox *c = new QComboBox;
    c->setEditable(true);
    c->addItems(nuclides);
    c->setCurrentIndex(-1);
    c->setInsertPolicy(QComboBox::NoInsert);
    c->setToolTip(tip + tr(" A nuclide of the built-in table (AME2020) or Z,A,mass (nuclear mass in u)."));
    connect(c, SIGNAL(editTextChanged(const QString &)), this, SLOT(kinematicsEdited()));
    return c;
  };
  beamCombo = nuclideCombo(tr("beam=: the projectile."));
  targetCombo = nuclideCombo(tr("target=: the target."));
  spectatorCombo = nuclideCombo(tr("spectator=: the spectator s of the Trojan horse a = x + s."));
  beamEnergyEdit = new ThmNumberSpin(" MeV", 0.0, 1.0e5, 1.0);
  beamEnergyEdit->setSpecialValueText(QString::fromUtf8("—"));  // 0: not given
  beamEnergyEdit->setToolTip(tr("Ebeam=: the lab beam energy (> 0)."));
  connect(beamEnergyEdit, SIGNAL(valueChanged(double)), this, SLOT(kinematicsEdited()));
  lineshapeCheck = new QCheckBox(tr("Coulomb (lineshape=on)"));
  lineshapeCheck->setToolTip(
      tr("lineshape=on: the Coulomb field of the charged spectator skews and shifts each resonance (factor N_C per "
         "level; Mukhamedzhanov, Kadyrov & Pang, EPJA 56 (2020) 233). Needs the three-body reaction and the Brune "
         "parameterization."));
  connect(lineshapeCheck, SIGNAL(toggled(bool)), this, SLOT(lineshapeToggled(bool)));

  // Spectator-momentum window (ps=, psNodes=).
  const QString psPhysics =
      tr("The vertex is averaged over the accepted spectator momenta p_s with the weight |phi(p_s)|^2 p_s^2 "
         "(p_xA^2/2mu_xA = E + B + p_s^2/2mu_sx). Not together with a non-zero spectator energy (Model page) for "
         "the same entrance pair.");
  psKindCombo = new QComboBox;
  psKindCombo->addItem(tr("point (quasi-free)"), "delta");
  psKindCombo->addItem(QString::fromUtf8("Hulthén"), "hulthen");
  psKindCombo->addItem(tr("Gaussian"), "gauss");
  psKindCombo->addItem(tr("table"), "table");
  psKindCombo->setToolTip(tr("ps=: the momentum distribution |phi(p_s)|^2 of the spectator. Point: p_s = 0, "
                             "nothing written. ") +
                          psPhysics);
  connect(psKindCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(psEdited()));
  auto psSpin = [&](const QString &suffix, double hi, double step, const QString &tip) {
    ThmNumberSpin *e = new ThmNumberSpin(suffix, 0.0, hi, step);
    e->setToolTip(tip);
    connect(e, SIGNAL(valueChanged(double)), this, SLOT(psEdited()));
    return e;
  };
  const QString perFm = QString::fromUtf8(" fm⁻¹");
  psMinEdit = psSpin(" MeV/c", 1.0e4, 1.0, tr("p_min of the accepted window (>= 0)."));
  psMaxEdit = psSpin(" MeV/c", 1.0e4, 1.0, tr("p_max of the accepted window (>= p_min; = p_min: one point)."));
  psCustomCheck = new QCheckBox(tr("custom a, b"));
  psCustomCheck->setToolTip(QString::fromUtf8("Hulthén phi(p) ~ 1/(a^2+q^2) - 1/(b^2+q^2). Unticked: the "
                                              "deuteron, a = 0.2317, b = 1.202 fm^-1 (Tribble 2014 eq. 4.4)."));
  connect(psCustomCheck, SIGNAL(toggled(bool)), this, SLOT(psEdited()));
  psAEdit = psSpin(perFm, 100.0, 0.01, tr("Hulthen a (> 0)."));
  psBEdit = psSpin(perFm, 100.0, 0.01, tr("Hulthen b (> a)."));
  psAEdit->setWrittenText("0.2317");
  psBEdit->setWrittenText("1.202");
  psFwhmEdit = psSpin(" MeV/c", 1.0e4, 1.0, tr("FWHM of |phi(p_s)|^2 = exp(-4 ln2 p^2/FWHM^2) (> 0)."));
  psTableEdit = new QLineEdit;
  psTableEdit->setPlaceholderText(tr("file"));
  psTableEdit->setToolTip(tr("ps=table:<file>: two columns, p_s (MeV/c, increasing) and the weight per unit p_s; "
                             "the window is its range. Relative to the project directory."));
  connect(psTableEdit, SIGNAL(textEdited(const QString &)), this, SLOT(psEdited()));
  psTableButton = new QPushButton("...");
  psTableButton->setToolTip(tr("Choose the table; a file inside the project directory is stored relative to it."));
  connect(psTableButton, SIGNAL(clicked()), this, SLOT(chooseTable()));
  psNodesSpin = new QSpinBox;
  psNodesSpin->setRange(1, 64);
  psNodesSpin->setValue(16);
  psNodesSpin->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
  psNodesSpin->setToolTip(tr("psNodes=: Gauss-Legendre nodes on the window (1-64, default 16)."));
  connect(psNodesSpin, SIGNAL(valueChanged(int)), this, SLOT(psNodesChanged(int)));

  // Derived values, as AZURE2 prints them (the whole text: derivedText, in the tooltips).
  auto value = [this]() {
    QLabel *l = new QLabel(QString::fromUtf8("—"));
    l->setTextInteractionFlags(Qt::TextSelectableByMouse);
    return l;
  };
  bindingValue = value();
  qfValue = value();
  zetaValue = value();
  meanTsValue = value();
  messageIcon = new QLabel;
  const int icon = style()->pixelMetric(QStyle::PM_SmallIconSize);
  messageIcon->setPixmap(style()->standardIcon(QStyle::SP_MessageBoxWarning).pixmap(icon, icon));
  messageIcon->setAlignment(Qt::AlignTop);
  messageLabel = new QLabel;
  messageLabel->setWordWrap(true);
  messageLabel->setTextInteractionFlags(Qt::TextSelectableByMouse);

  // Sections: four-column forms (label, field, label, field) whose label
  // columns have one width in every section, so that the fields line up.
  QList<QGridLayout *> forms;
  QList<QLabel *> leftLabels, rightLabels;
  auto label = [&](const QString &text, bool left, const QString &tip = QString()) {
    QLabel *l = new QLabel(text);
    l->setToolTip(tip);
    (left ? leftLabels : rightLabels) << l;
    return l;
  };
  auto form = [&]() {
    QGridLayout *g = new QGridLayout;
    g->setHorizontalSpacing(8);
    g->setVerticalSpacing(6);
    g->setColumnStretch(1, 1);
    g->setColumnStretch(3, 1);
    forms << g;
    return g;
  };
  const Qt::Alignment right = Qt::AlignRight | Qt::AlignVCenter;

  QGroupBox *experimentBox = new QGroupBox(tr("Experiment"));
  QGridLayout *xl = form();
  xl->addWidget(label(tr("Name:"), true), 0, 0, right);
  xl->addWidget(nameEdit, 0, 1);
  xl->addWidget(label(tr("Background:"), false), 0, 2, right);
  xl->addWidget(backgroundCombo, 0, 3);
  xl->addWidget(label(tr("Segments:"), true), 1, 0, Qt::AlignRight | Qt::AlignTop);
  xl->addWidget(segmentList, 1, 1, 1, 3);
  experimentBox->setLayout(xl);

  kinematicsBox = new QGroupBox(tr("Three-body reaction"));
  kinematicsBox->setCheckable(true);
  kinematicsBox->setChecked(false);
  kinematicsBox->setToolTip(
      tr("beam, target, spectator and Ebeam: all four or none. One of beam and target is a nucleus of the "
         "segments' entrance pair, the other the second nucleus plus the spectator."));
  connect(kinematicsBox, SIGNAL(toggled(bool)), this, SLOT(kinematicsEdited()));
  QGridLayout *kl = form();
  kl->addWidget(label(tr("Beam:"), true), 0, 0, right);
  kl->addWidget(beamCombo, 0, 1);
  kl->addWidget(label(tr("Target:"), false), 0, 2, right);
  kl->addWidget(targetCombo, 0, 3);
  kl->addWidget(label(tr("Spectator:"), true), 1, 0, right);
  kl->addWidget(spectatorCombo, 1, 1);
  kl->addWidget(label(tr("Beam energy:"), false, tr("Lab frame")), 1, 2, right);
  kl->addWidget(beamEnergyEdit, 1, 3);
  kl->addWidget(label("B(x+s):", true, tr("Binding energy of the Trojan horse a = x + s")), 2, 0, right);
  kl->addWidget(bindingValue, 2, 1);
  kl->addWidget(label(QString::fromUtf8("E<sub>qf</sub>:"), false, tr("Quasi-free energy E(x+A) - B")), 2, 2, right);
  kl->addWidget(qfValue, 2, 3);
  // The line shape needs the reaction: in its section.
  kl->addWidget(label(tr("Line shape:"), true), 3, 0, right);
  kl->addWidget(lineshapeCheck, 3, 1);
  kl->addWidget(label(QString::fromUtf8("\u03b6:"), false, tr("zeta at the lowest and highest data point")), 3, 2,
                right);
  kl->addWidget(zetaValue, 3, 3);
  kinematicsBox->setLayout(kl);

  psBox = new QGroupBox(tr("Spectator momentum window"));
  psBox->setToolTip(psPhysics);
  // Only the fields of the chosen distribution are shown (showPsRows).
  QGridLayout *pl = form();
  pl->addWidget(label(tr("Distribution:"), true), 0, 0, right);
  pl->addWidget(psKindCombo, 0, 1);
  pl->addWidget(psCustomCheck, 0, 3);
  QLabel *minLabel = label(QString::fromUtf8("p<sub>min</sub>:"), true);
  QLabel *maxLabel = label(QString::fromUtf8("p<sub>max</sub>:"), false);
  pl->addWidget(minLabel, 1, 0, right);
  pl->addWidget(psMinEdit, 1, 1);
  pl->addWidget(maxLabel, 1, 2, right);
  pl->addWidget(psMaxEdit, 1, 3);
  QLabel *tableLabel = label(tr("Table:"), true);
  QHBoxLayout *tl = new QHBoxLayout;
  tl->setContentsMargins(0, 0, 0, 0);
  tl->addWidget(psTableEdit, 1);
  tl->addWidget(psTableButton);
  pl->addWidget(tableLabel, 1, 0, right);
  pl->addLayout(tl, 1, 1, 1, 3);
  QLabel *aLabel = label("a:", true), *bLabel = label("b:", false);
  pl->addWidget(aLabel, 2, 0, right);
  pl->addWidget(psAEdit, 2, 1);
  pl->addWidget(bLabel, 2, 2, right);
  pl->addWidget(psBEdit, 2, 3);
  QLabel *fwhmLabel = label(tr("FWHM:"), true);
  pl->addWidget(fwhmLabel, 2, 0, right);
  pl->addWidget(psFwhmEdit, 2, 1);
  QLabel *meanLabel = label(QString::fromUtf8("\u27e8T<sub>s</sub>\u27e9:"), true, tr("Mean spectator energy"));
  QLabel *nodesLabel = label(tr("Nodes:"), false);
  pl->addWidget(meanLabel, 3, 0, right);
  pl->addWidget(meanTsValue, 3, 1);
  pl->addWidget(nodesLabel, 3, 2, right);
  pl->addWidget(psNodesSpin, 3, 3);
  psWindowRow_ = {minLabel, psMinEdit, maxLabel, psMaxEdit};
  psHulthenRow_ = {psCustomCheck, aLabel, psAEdit, bLabel, psBEdit};
  psGaussRow_ = {fwhmLabel, psFwhmEdit};
  psTableRow_ = {tableLabel, psTableEdit, psTableButton};
  psNodesRow_ = {nodesLabel, psNodesSpin, meanLabel, meanTsValue};
  psBox->setLayout(pl);

  int leftWidth = 0, rightWidth = 0;
  for (QLabel *l : leftLabels) leftWidth = std::max(leftWidth, l->sizeHint().width());
  for (QLabel *l : rightLabels) rightWidth = std::max(rightWidth, l->sizeHint().width());
  for (QGridLayout *g : forms) {
    g->setColumnMinimumWidth(0, leftWidth);
    g->setColumnMinimumWidth(2, rightWidth);
  }

  QHBoxLayout *message = new QHBoxLayout;
  message->setContentsMargins(0, 0, 0, 0);
  message->addWidget(messageIcon);
  message->addWidget(messageLabel, 1);

  editorBox = new QWidget;
  QVBoxLayout *el = new QVBoxLayout;
  el->setContentsMargins(0, 0, 0, 0);
  el->setSpacing(8);
  el->addWidget(experimentBox);
  el->addWidget(kinematicsBox);
  el->addWidget(psBox);
  el->addLayout(message);
  el->addStretch(1);
  editorBox->setLayout(el);

  // The editor scrolls rather than squeezing its sections when space is short.
  QScrollArea *editorScroll = new QScrollArea;
  editorScroll->setWidget(editorBox);
  editorScroll->setWidgetResizable(true);
  editorScroll->setFrameShape(QFrame::NoFrame);
  QHBoxLayout *mainLayout = new QHBoxLayout;
  mainLayout->setSpacing(12);
  mainLayout->addWidget(listBox, 0);
  mainLayout->addWidget(editorScroll, 1);
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
               (r.lineshape ? tr(", line shape") : QString()) + (r.hasWindow() ? tr(", p_s window") : QString());
  const QString cells[4] = {r.name, ThmExperimentRecord::segmentsListText(r.segments), r.background, reaction};
  for (int c = 0; c < 4; c++) {
    QTableWidgetItem *item = experimentTable->item(row, c);
    if (!item) {
      item = new QTableWidgetItem;
      experimentTable->setItem(row, c, item);
    }
    item->setText(cells[c]);
    QString tip = tr("Background: %1\nReaction: %2").arg(r.background, reaction);
    if (!r.extraTokens.isEmpty()) tip += "\n" + tr("Also kept as written: %1").arg(r.extraTokens.join(' '));
    item->setToolTip(tip);
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
  beamEnergyEdit->setWrittenText(r.beamEnergy);
  lineshapeCheck->setChecked(r.lineshape);
  const bool complete = !r.beam.isEmpty() && !r.target.isEmpty() && !r.spectator.isEmpty() && !r.beamEnergy.isEmpty();
  lineshapeCheck->setEnabled(complete);
  loadPs(r);
  psBox->setEnabled(complete);
  showDerived(r);
  loading_ = false;
}

void ThmExperimentsPage::showDerived(const ThmExperimentRecord &x) {
  QString why, info = derivedInfo(x, &why);
  derivedText_ = info.isEmpty() ? why : why.isEmpty() ? info : info + "\n" + why;
  // The values in the sections, compact; the whole text in their tooltips.
  const QString none = QString::fromUtf8("\u2014");
  auto mev = [](double v) { return QString::number(v, 'g', 6) + " MeV"; };
  Reaction r;
  QString ignored;
  const bool haveReaction = reaction(x, r, &ignored);
  bindingValue->setText(haveReaction ? mev(r.bind) : none);
  qfValue->setText(haveReaction ? mev(r.exa - r.bind) : none);
  qfValue->setToolTip(haveReaction ? tr("E(x+A) = %1, E_qf = E(x+A) - B").arg(mev(r.exa)) : QString());
  QString zeta;
  if (haveReaction && x.lineshape) lineshapeInfo(x, &ignored, &zeta);
  zetaValue->setText(zeta.isEmpty() ? none : zeta);
  ThmSpectatorWindow window;
  const bool haveWindow = haveReaction && x.hasWindow() && !windowInfo(x, &ignored, &window).isEmpty();
  meanTsValue->setText(haveWindow ? mev(window.MeanEs()) : none);
  for (QLabel *l : {bindingValue, zetaValue, meanTsValue}) l->setToolTip(info);
  messageLabel->setText(why);
  messageLabel->setVisible(!why.isEmpty());
  messageIcon->setVisible(!why.isEmpty());
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
  r.beamEnergy = on ? beamEnergyEdit->writtenText() : QString();
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
  // So does the spectator-momentum window.
  psBox->setEnabled(lineshapeCheck->isEnabled());
  if (!on && (!r.ps.isEmpty() || !r.psNodes.isEmpty())) {
    r.ps.clear();
    r.psNodes.clear();
    loading_ = true;
    loadPs(r);
    loading_ = false;
  }
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
  out.mX = mX;
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
  QString shapeError, windowError;
  if (x.lineshape) {
    const QString shape = lineshapeInfo(x, &shapeError);
    if (!shape.isEmpty()) text += "\n" + shape;
  }
  if (x.hasWindow()) {
    const QString window = windowInfo(x, &windowError);
    if (!window.isEmpty()) text += "\n" + window;
  }
  if (error) *error = !shapeError.isEmpty() ? shapeError : windowError;
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

QString ThmExperimentsPage::lineshapeInfo(const ThmExperimentRecord &x, QString *error, QString *compact) const {
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
  if (compact) *compact = parts.join("; ").replace(" into ", " (") .replace("; ", "); ") + ")";
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

// ---------------------------------------------------------------------------
// Spectator-momentum window

namespace {

// "pmin-pmax" as written: split at the first '-' with a number on both sides
// (the engine's ReadWindow); otherwise at the first '-' after the start.
void splitWindow(const QString &text, QString &lo, QString &hi) {
  for (int k = 1; k + 1 < text.size(); k++) {
    if (text[k] != '-') continue;
    double a, b;
    if (readWholeDouble(text.left(k), a) && readWholeDouble(text.mid(k + 1), b)) {
      lo = text.left(k);
      hi = text.mid(k + 1);
      return;
    }
  }
  const int dash = text.indexOf('-', 1);
  lo = dash < 0 ? text : text.left(dash);
  hi = dash < 0 ? QString() : text.mid(dash + 1);
}

}  // namespace

QString ThmExperimentsPage::projectRelative(const QString &file, const QString &projectDir) {
  if (projectDir.isEmpty() || file.isEmpty()) return file;
  const QString rel = QDir(projectDir).relativeFilePath(file);
  return rel.startsWith("..") || QDir::isAbsolutePath(rel) ? file : rel;
}

void ThmExperimentsPage::loadPs(const ThmExperimentRecord &r) {
  // The controls show the value as written, so that composing them again
  // (psText) gives it back unchanged.
  QString kind = "delta", lo, hi, a = "0.2317", b = "1.202", fwhm, table;
  bool custom = false;
  if (r.ps.startsWith("table:")) {
    kind = "table";
    table = r.ps.mid(6);
  } else if (!r.ps.isEmpty() && r.ps != "delta") {
    const QStringList f = r.ps.split(':');
    kind = f[0];
    if (kind == "hulthen" && f.size() == 3) {
      custom = true;
      const QStringList ab = f[1].split(',');
      a = ab.value(0);
      b = ab.value(1);
    } else if (kind == "gauss" && f.size() == 3) {
      fwhm = f[1];
    }
    splitWindow(f.last(), lo, hi);
    if (f.size() < 2) lo = hi = QString();
  }
  const bool was = loading_;
  loading_ = true;
  int at = psKindCombo->findData(kind);
  psKindCombo->setCurrentIndex(at >= 0 ? at : 0);
  psMinEdit->setWrittenText(lo);
  psMaxEdit->setWrittenText(hi);
  psCustomCheck->setChecked(custom);
  psAEdit->setWrittenText(a);
  psBEdit->setWrittenText(b);
  psFwhmEdit->setWrittenText(fwhm);
  psTableEdit->setText(table);
  int nodes = 16;
  if (!r.psNodes.isEmpty()) {
    bool ok = false;
    const int n = r.psNodes.trimmed().toInt(&ok);
    if (ok) nodes = n;
  }
  psNodesSpin->setValue(std::min(64, std::max(1, nodes)));
  loading_ = was;
  showPsRows();
}

void ThmExperimentsPage::showPsRows() {
  // Only the fields of the chosen distribution are shown.
  const QString kind = psKindCombo->currentData().toString();
  for (QWidget *w : psWindowRow_) w->setVisible(kind == "hulthen" || kind == "gauss");
  for (QWidget *w : psHulthenRow_) w->setVisible(kind == "hulthen");
  for (QWidget *w : psGaussRow_) w->setVisible(kind == "gauss");
  for (QWidget *w : psTableRow_) w->setVisible(kind == "table");
  for (QWidget *w : psNodesRow_) w->setVisible(kind != "delta");
  psAEdit->setEnabled(psCustomCheck->isChecked());
  psBEdit->setEnabled(psCustomCheck->isChecked());
}

QString ThmExperimentsPage::psText() const {
  const QString kind = psKindCombo->currentData().toString();
  const QString window = psMinEdit->writtenText() + "-" + psMaxEdit->writtenText();
  if (kind == "hulthen")
    return psCustomCheck->isChecked()
               ? "hulthen:" + psAEdit->writtenText() + "," + psBEdit->writtenText() + ":" + window
               : "hulthen:" + window;
  if (kind == "gauss") return "gauss:" + psFwhmEdit->writtenText() + ":" + window;
  if (kind == "table") return "table:" + psTableEdit->text().trimmed();
  return QString();
}

void ThmExperimentsPage::psEdited() {
  showPsRows();
  if (loading_ || current_ < 0) return;
  ThmExperimentRecord &r = records_[current_];
  if (!psCustomCheck->isChecked()) {
    // Back to the deuteron's values, shown as the engine's defaults.
    const QSignalBlocker blockA(psAEdit), blockB(psBEdit);
    psAEdit->setWrittenText("0.2317");
    psBEdit->setWrittenText("1.202");
  }
  r.ps = psText();
  if (r.ps.isEmpty() && !r.psNodes.isEmpty()) {
    // psNodes needs a window: dropped with it.
    r.psNodes.clear();
    loading_ = true;
    psNodesSpin->setValue(16);
    loading_ = false;
  }
  refreshRow(current_);
  showDerived(r);
}

void ThmExperimentsPage::psNodesChanged(int n) {
  if (loading_ || current_ < 0) return;
  ThmExperimentRecord &r = records_[current_];
  r.psNodes = n == 16 ? QString() : QString::number(n);  // 16 is the engine's default
  showDerived(r);
}

void ThmExperimentsPage::chooseTable() {
  const QString start = projectDir_.isEmpty() ? QDir::currentPath() : projectDir_;
  const QString file = QFileDialog::getOpenFileName(this, tr("Spectator-momentum table"), start,
                                                    tr("Tables (*.dat *.txt);;All files (*)"));
  if (file.isEmpty()) return;
  psTableEdit->setText(projectRelative(file, projectDir_));
  psEdited();
}

void ThmExperimentsPage::refreshDerived() {
  if (current_ >= 0) showDerived(records_.at(current_));
}

QString ThmExperimentsPage::windowInfo(const ThmExperimentRecord &x, QString *error, ThmSpectatorWindow *out) const {
  Reaction r;
  if (!reaction(x, r, error)) return QString();
  // The engine's parse of the line, and its table reader (Config::ReadThmBlock).
  std::vector<ThmExperiment> parsed;
  std::string why = ParseThmExperimentLine(x.line().toStdString(), parsed);
  if (!why.empty() || parsed.empty()) {
    if (error) *error = "<thm> " + QString::fromStdString(why);
    return QString();
  }
  ThmExperiment &e = parsed.front();
  if (e.psKind == ThmExperiment::PS_DELTA) return QString();
  if (spectatorEnergy_(r.pairKey) != 0.0) {
    if (error)
      *error = tr("a ps window and spectatorEnergy both set the spectator motion of entrance pair %1; use one "
                  "(ps=delta keeps spectatorEnergy).")
                   .arg(r.pairKey);
    return QString();
  }
  if (e.psKind == ThmExperiment::PS_TABLE) {
    QString path = QString::fromStdString(e.psTable);
    if (QFileInfo(path).isRelative() && !projectDir_.isEmpty()) path = QDir(projectDir_).filePath(path);
    why = ReadThmPsTable(QFile::encodeName(path).toStdString(), e.psTableP, e.psTableW);
    if (!why.empty()) {
      if (error) *error = "ps: " + QString::fromStdString(why);
      return QString();
    }
    e.psMin = e.psTableP.front();
    e.psMax = e.psTableP.back();
  }
  // mu_sx as EData::BuildThmGroups takes it: x and the spectator.
  const double muSx = r.mX * r.spectator.mass / (r.mX + r.spectator.mass) * kAmu;
  ThmSpectatorWindow window;
  why = BuildThmSpectatorWindow(e, muSx, window);
  if (!why.empty()) {
    if (error) *error = "ps: " + QString::fromStdString(why);
    return QString();
  }
  if (out) *out = window;
  return tr("p_s window: %1; mu_sx = %2 MeV, T_s = p_s^2/2mu_sx from %3 to %4 MeV, <T_s> = %5 MeV")
      .arg(QString::fromStdString(window.description))
      .arg(QString::number(muSx, 'g', 6))
      .arg(QString::number(window.es.front(), 'g', 6))
      .arg(QString::number(window.es.back(), 'g', 6))
      .arg(QString::number(window.MeanEs(), 'g', 6));
}
