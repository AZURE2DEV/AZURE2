// The Potential panel of the Particle Pairs tab (PotentialPanel.h, docs/upgrade/qt_gui.md section 19).
#include "PotentialPanel.h"

#include <cmath>
#include <map>

#include <QCheckBox>
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#include <QRegularExpression>
#include <QSpinBox>
#include <QSplitter>
#include <QTableWidget>
#include <QTimer>
#include <QVBoxLayout>

#include "InfoDialog.h"
#include "PhysicsWorker.h"
#include "ProjectHolder.h"
#include "RecordDialog.h"

namespace front = azure2::front;

namespace {

const char *const kSection = "potentials";

QString Q(const std::string &s) { return QString::fromStdString(s); }

QString Num(double x, int digits = 6) { return QString::number(x, 'g', digits); }

/// The columns of the waves' table.
enum Column { kL = 0, kJ, kMode, kDepth, kWidth, kBoundary, kSubtract, kCapture, kIncluded, kColumns };

QString Spin(int twoJ) { return twoJ % 2 ? QString("%1/2").arg(twoJ) : QString::number(twoJ / 2); }

}  // namespace

PotentialPanel::PotentialPanel(ProjectHolder *holder, Hooks hooks, QWidget *parent)
    : QWidget(parent), holder_(holder), hooks_(std::move(hooks)), worker_(new PhysicsWorker(this)), debounce_(new QTimer(this)) {
  setObjectName("potentialPanel");
  debounce_->setSingleShot(true);
  debounce_->setInterval(250);
  connect(debounce_, &QTimer::timeout, this, &PotentialPanel::refresh);

  // the head: the pair, on/off, the record's buttons, one line of what this is
  title_ = new QLabel;
  title_->setTextFormat(Qt::RichText);
  onCheck_ = new QCheckBox(tr("Potential on"));
  onCheck_->setObjectName("potentialOn");
  onCheck_->setToolTip(tr("Off: the record stays in the file as comment lines (\"#off ...\") and does nothing; on again "
                          "gives it back exactly as it was."));
  connect(onCheck_, &QCheckBox::toggled, this, &PotentialPanel::switchOn);
  create_ = new QPushButton(tr("Add a potential..."));
  create_->setObjectName("potentialCreate");
  create_->setToolTip(tr("A potential record for this pair (<potentials>): its Coulomb radius and central term"));
  connect(create_, &QPushButton::clicked, this, &PotentialPanel::createPotential);
  edit_ = new QPushButton(tr("Edit record..."));
  edit_->setToolTip(tr("The whole record in the record pop-up: Coulomb radius, central term, spin-orbit term, matching "
                       "radius, step"));
  connect(edit_, &QPushButton::clicked, this, &PotentialPanel::editRecord);
  QPushButton *helpButton = new QPushButton(tr("?"));
  helpButton->setMaximumWidth(28);
  helpButton->setToolTip(tr("What the modes, the poles and the tuning are"));
  connect(helpButton, &QPushButton::clicked, this, [this]() {
    InfoDialog d(helpText, this, tr("The potential of a pair"));
    d.exec();
  });
  help_ = new QLabel(tr("The nuclear potential of the selected pair, wave by wave: the mode of each partial wave, the "
                        "poles the potential creates inside the channel radius, and its own scattering length."));
  help_->setWordWrap(true);
  help_->setObjectName("potentialHelp");
  busy_ = new QLabel;
  QHBoxLayout *head = new QHBoxLayout;
  head->addWidget(title_);
  head->addWidget(onCheck_);
  head->addStretch(1);
  head->addWidget(busy_);
  head->addWidget(create_);
  head->addWidget(edit_);
  head->addWidget(helpButton);

  // the waves
  waves_ = new QTableWidget(0, kColumns);
  waves_->setObjectName("potentialWaves");
  waves_->setHorizontalHeaderLabels({tr("l"), tr("J"), tr("Mode"), tr("Depth (MeV)"), tr("Width (fm)"), tr("B_bg"),
                                     tr("Subtract n"), tr("Capture interior"), tr("In the calculation")});
  waves_->horizontalHeaderItem(kMode)->setToolTip(
      tr("Central: no wave line, the central potential acts beyond a_c.  Exterior only: the wave's own depth acts beyond "
         "a_c (hybrid channel functions).  Exterior + interior: also the potential's R-function inside a_c as the "
         "background (background interior boundary B).  None: point-Coulomb channel functions; the depth gives the "
         "bound-state shape of interior capture only (only bound_state)."));
  waves_->horizontalHeaderItem(kBoundary)->setToolTip(tr("B_bg, the boundary condition the background is defined with "
                                                          "(dimensionless): the split level/background depends on it"));
  waves_->horizontalHeaderItem(kSubtract)->setToolTip(tr("The node number of a box state whose pole is subtracted from "
                                                          "the background: a state a level of the project represents"));
  waves_->horizontalHeaderItem(kIncluded)->setToolTip(tr("A J^pi is in the calculation only when an active level carries "
                                                          "it, or with \"hard-sphere scattering in every J^pi\" "
                                                          "(includeAllJPi) on"));
  waves_->verticalHeader()->setVisible(false);
  waves_->setSelectionBehavior(QAbstractItemView::SelectRows);
  waves_->setSelectionMode(QAbstractItemView::SingleSelection);
  waves_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
  waves_->horizontalHeader()->setStretchLastSection(true);
  connect(waves_, &QTableWidget::itemChanged, this, [this](QTableWidgetItem *item) {
    if (!filling_) cellEdited(item->row(), item->column());
  });
  connect(waves_, &QTableWidget::currentCellChanged, this, [this](int row, int, int previous, int) {
    if (!filling_ && row != previous) schedule();
  });
  addWave_ = new QPushButton(tr("Add wave..."));
  addWave_->setToolTip(tr("A wave line of its own for a partial wave (l, and J when the potential depends on J)"));
  connect(addWave_, &QPushButton::clicked, this, &PotentialPanel::addWaveAsked);
  removeWave_ = new QPushButton(tr("Remove wave line"));
  removeWave_->setToolTip(tr("The selected wave back to the central potential (its line removed)"));
  connect(removeWave_, &QPushButton::clicked, this, &PotentialPanel::removeWaveLine);
  apply_ = new QPushButton(tr("Apply"));
  apply_->setObjectName("potentialApply");
  apply_->setProperty("primary", true);
  connect(apply_, &QPushButton::clicked, this, [this]() {
    QString why;
    if (!apply(why)) QMessageBox::warning(this, tr("Potential"), why);
  });
  revert_ = new QPushButton(tr("Revert"));
  connect(revert_, &QPushButton::clicked, this, &PotentialPanel::revert);
  QHBoxLayout *waveButtons = new QHBoxLayout;
  waveButtons->addWidget(addWave_);
  waveButtons->addWidget(removeWave_);
  waveButtons->addStretch(1);
  waveButtons->addWidget(revert_);
  waveButtons->addWidget(apply_);
  messages_ = new QLabel;
  messages_->setObjectName("potentialMessages");
  messages_->setWordWrap(true);
  messages_->setTextFormat(Qt::RichText);
  messages_->setTextInteractionFlags(Qt::TextSelectableByMouse);

  // the poles of the selected wave and the levels of its J^pi
  waveTitle_ = new QLabel;
  waveTitle_->setWordWrap(true);
  waveTitle_->setTextFormat(Qt::RichText);
  eMax_ = new QDoubleSpinBox;
  eMax_->setRange(0.0, 500.0);
  eMax_->setValue(30.0);
  eMax_->setSuffix(tr(" MeV"));
  eMax_->setToolTip(tr("The highest pole listed (E_c, c.m. of the pair)"));
  connect(eMax_, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &PotentialPanel::schedule);
  poles_ = new QTableWidget(0, 6);
  poles_->setObjectName("potentialPoles");
  poles_->setHorizontalHeaderLabels({tr("Nodes"), tr("E_c (MeV)"), tr("E_x (MeV)"), tr("γ² (MeV)"), tr("θ²"), tr("")});
  poles_->horizontalHeaderItem(4)->setToolTip(tr("γ_p² / γ_W², γ_W² = 3ħ²/(2μa²) (the engine's Wigner limit)"));
  poles_->verticalHeader()->setVisible(false);
  poles_->setEditTriggers(QAbstractItemView::NoEditTriggers);
  poles_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
  poles_->horizontalHeader()->setStretchLastSection(true);
  levels_ = new QTableWidget(0, 5);
  levels_->setObjectName("potentialLevels");
  levels_->setHorizontalHeaderLabels({tr("Level"), tr("E_x (MeV)"), tr("E_c (MeV)"), tr("Channel value"), tr("Active")});
  levels_->verticalHeader()->setVisible(false);
  levels_->setEditTriggers(QAbstractItemView::NoEditTriggers);
  levels_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
  levels_->horizontalHeader()->setStretchLastSection(true);
  QGroupBox *polesBox = new QGroupBox(tr("Poles of the potential inside a_c (box states)"));
  QVBoxLayout *pl = new QVBoxLayout(polesBox);
  QHBoxLayout *polesHead = new QHBoxLayout;
  polesHead->addWidget(new QLabel(tr("Up to E_c =")));
  polesHead->addWidget(eMax_);
  polesHead->addStretch(1);
  pl->addLayout(polesHead);
  pl->addWidget(poles_);
  QGroupBox *levelsBox = new QGroupBox(tr("The project's levels of this J\u03c0"));
  QVBoxLayout *ll = new QVBoxLayout(levelsBox);
  ll->addWidget(levels_);

  // tuning
  readOut_ = new QLabel;
  readOut_->setObjectName("potentialReadOut");
  readOut_->setWordWrap(true);
  readOut_->setTextFormat(Qt::RichText);
  lengthEdit_ = new QLineEdit;
  lengthEdit_->setPlaceholderText(tr("a (fm)"));
  lengthEdit_->setMaximumWidth(110);
  lengthButton_ = new QPushButton(tr("Find depth"));
  lengthButton_->setToolTip(tr("The depth of this wave that gives the scattering length a, with the number of bound "
                               "states it has now (its branch)"));
  connect(lengthButton_, &QPushButton::clicked, this, &PotentialPanel::findLength);
  bindCombo_ = new QComboBox;
  bindCombo_->setMinimumContentsLength(24);
  bindCombo_->setToolTip(tr("A level of this J\u03c0 below the pair's threshold"));
  bindButton_ = new QPushButton(tr("Find depth"));
  bindButton_->setToolTip(tr("The depth (within ±30 % of the present one) that binds this wave at the level's energy, "
                             "with the potential's own ANC"));
  connect(bindButton_, &QPushButton::clicked, this, &PotentialPanel::findBinding);
  tuning_ = new QLabel;
  tuning_->setObjectName("potentialTuning");
  tuning_->setWordWrap(true);
  tuning_->setTextFormat(Qt::RichText);
  useDepth_ = new QPushButton(tr("Use this depth"));
  useDepth_->setEnabled(false);
  connect(useDepth_, &QPushButton::clicked, this, &PotentialPanel::useDepth);
  QGroupBox *tuneBox = new QGroupBox(tr("Tuning the wave's depth"));
  QFormLayout *tl = new QFormLayout(tuneBox);
  tl->addRow(tr("The potential alone:"), readOut_);
  QHBoxLayout *lengthRow = new QHBoxLayout;
  lengthRow->addWidget(lengthEdit_);
  lengthRow->addWidget(new QLabel(tr("fm")));
  lengthRow->addWidget(lengthButton_);
  lengthRow->addStretch(1);
  tl->addRow(tr("Depth for a ="), lengthRow);
  QHBoxLayout *bindRow = new QHBoxLayout;
  bindRow->addWidget(bindCombo_, 1);
  bindRow->addWidget(bindButton_);
  tl->addRow(tr("Depth to bind:"), bindRow);
  QHBoxLayout *resultRow = new QHBoxLayout;
  resultRow->addWidget(tuning_, 1);
  resultRow->addWidget(useDepth_);
  tl->addRow(resultRow);

  QWidget *lower = new QWidget;
  QGridLayout *lg = new QGridLayout(lower);
  lg->setContentsMargins(0, 0, 0, 0);
  lg->addWidget(waveTitle_, 0, 0, 1, 2);
  lg->addWidget(polesBox, 1, 0);
  lg->addWidget(levelsBox, 1, 1);
  lg->addWidget(tuneBox, 2, 0, 1, 2);
  lg->setColumnStretch(0, 3);
  lg->setColumnStretch(1, 2);

  QWidget *upper = new QWidget;
  QVBoxLayout *ul = new QVBoxLayout(upper);
  ul->setContentsMargins(0, 0, 0, 0);
  ul->addWidget(waves_);
  ul->addLayout(waveButtons);
  ul->addWidget(messages_);
  QSplitter *split = new QSplitter(Qt::Vertical);
  split->addWidget(upper);
  split->addWidget(lower);
  split->setStretchFactor(0, 1);
  split->setStretchFactor(1, 1);
  split->setChildrenCollapsible(false);
  waves_->setMinimumHeight(8 * waves_->verticalHeader()->defaultSectionSize());

  QVBoxLayout *layout = new QVBoxLayout(this);
  layout->addLayout(head);
  layout->addWidget(help_);
  layout->addWidget(split, 1);

  connect(worker_, &PhysicsWorker::busyChanged, this, [this](bool on) { busy_->setText(on ? tr("computing...") : QString()); });
  connect(holder_, &ProjectHolder::opened, this, [this]() { setPair(0, QString()); });
  connect(holder_, &ProjectHolder::sectionEdited, this, [this](const QString &section) {
    if (section == kSection && !filling_) reload();
  });
  setPair(0, QString());
}

PotentialPanel::~PotentialPanel() = default;

bool PotentialPanel::busy() const { return worker_->busy() || debounce_->isActive(); }

void PotentialPanel::setPair(int pairKey, const QString &label) {
  pair_ = pairKey;
  pairLabel_ = label;
  reload();
}

void PotentialPanel::reload() {
  hasRecord_ = false;
  on_ = true;
  record_ = front::Record();
  if (pair_ > 0) {
    QString why;
    if (hooks_.commitLevels) hooks_.commitLevels(why);
    const front::PairPotential p = front::PotentialOfPair(holder_->project(), pair_);
    hasRecord_ = p.found;
    on_ = p.on;
    record_ = p.record;
  }
  title_->setText(pair_ <= 0 ? tr("<b>No pair selected</b>")
                             : tr("<b>Pair %1 (%2)</b>: %3").arg(pair_).arg(pairLabel_,
                                   hasRecord_ ? tr("potential <i>%1</i>").arg(Q(record_.name)) : tr("no potential")));
  filling_ = true;
  onCheck_->setChecked(on_);
  filling_ = false;
  onCheck_->setEnabled(hasRecord_);
  create_->setVisible(pair_ > 0 && !hasRecord_);
  edit_->setEnabled(hasRecord_ && on_);
  reportOk_ = false;
  reportRow_ = -1;
  depth_ = front::PotentialDepth();
  depthRow_ = -1;
  recordMessages_.clear();
  buildRows();
  fillTable();
  showReport();
  showMessages();
  if (hasRecord_) schedule();
}

void PotentialPanel::buildRows() {
  rows_.clear();
  if (!hasRecord_) return;
  // the record's wave lines by l, and one row of the central potential for each l without a line, up to the larger of
  // the Levels tab's maximum l and the lines' l
  std::map<int, QList<Row>> byL;
  int maxL = hooks_.maxL ? hooks_.maxL() : 2;
  for (std::size_t i = 0; i < record_.fields.size(); i++) {
    const front::Field &f = record_.fields[i];
    if (f.key != "wave") continue;
    Row r;
    r.line = true;
    r.field = (int)i;
    r.original = Q(f.value);
    std::string why;
    front::ParsePotentialWave(f.value, r.wave, why);
    maxL = std::max(maxL, r.wave.l);
    byL[r.wave.l] << r;
  }
  for (int l = 0; l <= maxL; l++) {
    if (byL.count(l)) {
      rows_ << byL[l];
      continue;
    }
    Row r;
    r.wave.l = l;
    rows_ << r;
  }
}

QString PotentialPanel::jpiText(int twoJ, int parity) const {
  return QString("%1%2").arg(Spin(twoJ), parity > 0 ? "+" : "-");
}

std::vector<front::PartialWaveState> PotentialPanel::waveStates(int row) const {
  if (row < 0 || row >= rows_.size() || !hooks_.levelRows) return {};
  const Row &r = rows_[row];
  return front::PartialWaveStates(hooks_.levelRows(), hooks_.includeAllJPi ? hooks_.includeAllJPi() : false,
                                  hooks_.maxOrbitalMomentum ? hooks_.maxOrbitalMomentum() : -1, pair_, r.wave.l, r.wave.twoJ);
}

void PotentialPanel::fillTable() {
  filling_ = true;
  const int keep = std::max(0, waves_->currentRow());
  waves_->clearContents();
  waves_->setRowCount(rows_.size());
  for (int i = 0; i < rows_.size(); i++) fillRow(i);
  if (!rows_.isEmpty()) waves_->setCurrentCell(std::min(keep, (int)rows_.size() - 1), kMode);
  const bool editable = hasRecord_ && on_;
  waves_->setEnabled(hasRecord_);
  addWave_->setEnabled(editable);
  removeWave_->setEnabled(editable);
  apply_->setEnabled(editable);
  revert_->setEnabled(editable);
  filling_ = false;
  updateStates();
}

void PotentialPanel::fillRow(int i) {
  const bool was = filling_;
  filling_ = true;
  const Row &r = rows_[i];
  const front::PotentialWave &w = r.wave;
  const bool editable = hasRecord_ && on_;
  const Mode mode = !r.line ? Central : w.background ? ExteriorInterior : w.boundStateOnly ? PointCoulomb : Exterior;
  auto item = [&](int column, const QString &text, bool canEdit, const QString &tip = QString()) {
    QTableWidgetItem *it = new QTableWidgetItem(text);
    Qt::ItemFlags flags = Qt::ItemIsSelectable | Qt::ItemIsEnabled;
    if (canEdit && editable) flags |= Qt::ItemIsEditable;
    it->setFlags(flags);
    if (!tip.isEmpty()) it->setToolTip(tip);
    if (!canEdit) it->setForeground(palette().color(QPalette::Disabled, QPalette::Text));
    waves_->setItem(i, column, it);
    return it;
  };
  QTableWidgetItem *lItem = item(kL, QString::number(w.l), false);
  lItem->setForeground(palette().color(QPalette::Text));
  QTableWidgetItem *jItem = item(kJ, w.J.empty() ? tr("every J") : Q(w.J), false);
  jItem->setForeground(palette().color(QPalette::Text));
  const QString central = tr("central");
  const front::PotentialWave units;  // the units the columns' headers name
  item(kDepth, w.depth.empty() ? central : Q(w.depth) + (w.depthUnit == units.depthUnit ? QString() : " " + Q(w.depthUnit)), true,
       tr("This wave's own V0 (empty: the central depth)"));
  item(kWidth, w.width.empty() ? central : Q(w.width) + (w.widthUnit == units.widthUnit ? QString() : " " + Q(w.widthUnit)), true,
       tr("This wave's own Gaussian r0 (empty: the central width)"));
  if (w.depth.empty()) waves_->item(i, kDepth)->setForeground(palette().color(QPalette::Disabled, QPalette::Text));
  if (w.width.empty()) waves_->item(i, kWidth)->setForeground(palette().color(QPalette::Disabled, QPalette::Text));
  const bool bg = mode == ExteriorInterior;
  item(kBoundary, bg ? Q(w.boundary) : QString(), bg);
  item(kSubtract, bg ? Q(w.subtract) : QString(), bg, tr("A node number of the poles below (empty: none)"));
  QTableWidgetItem *cap = item(kCapture, QString(), false);
  cap->setFlags(bg && editable ? (Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsUserCheckable) : Qt::ItemIsSelectable);
  cap->setCheckState(bg && w.capture ? Qt::Checked : Qt::Unchecked);
  cap->setToolTip(tr("capture interior: the interior direct capture of the background (needs external capture)"));
  item(kIncluded, QString(), false);
  QComboBox *combo = new QComboBox;
  combo->addItems({tr("Central (no wave line)"), tr("Exterior only"), tr("Exterior + interior"),
                   tr("None (point Coulomb)")});
  combo->setItemData(Central, tr("No wave line: the central potential acts beyond a_c in this wave"), Qt::ToolTipRole);
  combo->setItemData(Exterior, tr("A wave line without background: its depth acts beyond a_c"), Qt::ToolTipRole);
  combo->setItemData(ExteriorInterior, tr("background interior boundary B: also the potential inside a_c, as the "
                                          "background of the wave's R-matrix"), Qt::ToolTipRole);
  combo->setItemData(PointCoulomb, tr("only bound_state: point-Coulomb channel functions; the depth gives the bound-state "
                                      "shape of interior capture (needs a wave with capture interior)"), Qt::ToolTipRole);
  combo->setCurrentIndex(mode);
  combo->setEnabled(editable);
  connect(combo, QOverload<int>::of(&QComboBox::activated), this, [this, combo](int index) {
    for (int row = 0; row < waves_->rowCount(); row++)
      if (waves_->cellWidget(row, kMode) == combo) modeChosen(row, index);
  });
  waves_->setCellWidget(i, kMode, combo);
  filling_ = was;
}

void PotentialPanel::updateStates() {
  // whether each wave is in the calculation (the Levels tab's levels as they stand, the key)
  const bool was = filling_;
  filling_ = true;
  for (int i = 0; i < rows_.size() && i < waves_->rowCount(); i++) {
    QTableWidgetItem *it = waves_->item(i, kIncluded);
    if (!it) continue;
    QStringList in, out;
    for (const front::PartialWaveState &s : waveStates(i))
      (s.included ? in : out) << jpiText(s.twoJ, s.parity);
    if (out.isEmpty()) {
      it->setText(in.isEmpty() ? QString() : tr("yes"));
      it->setForeground(palette().color(QPalette::Text));
      it->setToolTip(tr("%1: an active level, or hard-sphere scattering in every J\u03c0").arg(in.join(", ")));
    } else {
      it->setText(tr("left out: %1").arg(out.join(", ")));
      it->setForeground(QColor(0xb0, 0x00, 0x20));
      it->setToolTip(tr("No active level carries %1 and hard-sphere scattering in every J\u03c0 (includeAllJPi) is off "
                        "(or does not reach this l): the J\u03c0 is not in the calculation, and this wave's potential does "
                        "nothing.").arg(out.join(", ")));
    }
  }
  filling_ = was;
}

int PotentialPanel::currentRow() const {
  const int r = waves_->currentRow();
  return r >= 0 && r < rows_.size() ? r : (rows_.isEmpty() ? -1 : 0);
}

int PotentialPanel::requestTwoJ(const Row &row) const { return row.wave.twoJ; }

QString PotentialPanel::rowText(const Row &row) const {
  return row.edited ? Q(front::PotentialWaveText(row.wave)) : row.original;
}

front::Record PotentialPanel::editedRecord() const {
  bool changed = false;
  for (const Row &r : rows_) changed = changed || r.edited || (r.line != (r.field >= 0));
  if (!changed) return record_;
  front::Record out = record_;
  out.fields.clear();
  std::map<int, const Row *> byField;
  QList<const Row *> fresh;
  for (const Row &r : rows_) {
    if (r.line && r.field >= 0) byField[r.field] = &r;
    else if (r.line) fresh << &r;
  }
  int lastWave = -1, firstTail = -1;
  for (std::size_t i = 0; i < record_.fields.size(); i++) {
    if (record_.fields[i].key == "wave") lastWave = (int)i;
    if (firstTail < 0 && (record_.fields[i].key == "matching_radius" || record_.fields[i].key == "step")) firstTail = (int)i;
  }
  // new lines after the last wave line, or before matching_radius and step, or at the end
  const int insertAfter = lastWave >= 0 ? lastWave : (firstTail >= 0 ? firstTail - 1 : (int)record_.fields.size() - 1);
  auto addFresh = [&]() {
    for (const Row *r : fresh) out.fields.push_back({"wave", rowText(*r).toStdString(), ""});
  };
  if (insertAfter < 0) addFresh();
  for (std::size_t i = 0; i < record_.fields.size(); i++) {
    const front::Field &f = record_.fields[i];
    if (f.key == "wave") {
      const auto it = byField.find((int)i);
      if (it != byField.end()) out.fields.push_back({"wave", rowText(*it->second).toStdString(), f.comment});
    } else {
      out.fields.push_back(f);
    }
    if ((int)i == insertAfter) addFresh();
  }
  return out;
}

void PotentialPanel::cellEdited(int row, int column) {
  if (row < 0 || row >= rows_.size()) return;
  Row &r = rows_[row];
  front::PotentialWave &w = r.wave;
  const QString text = waves_->item(row, column) ? waves_->item(row, column)->text().trimmed() : QString();
  const QString central = tr("central");
  if (column == kDepth || column == kWidth) {
    const QString value = text == central ? QString() : text.section(' ', 0, 0);
    if (column == kDepth) w.depth = value.toStdString();
    else w.width = value.toStdString();
    if (!r.line && !value.isEmpty()) r.line = true;  // a depth of its own makes a wave line
  } else if (column == kBoundary) {
    w.boundary = text.toStdString();
  } else if (column == kSubtract) {
    w.subtract = text.toStdString();
  } else if (column == kCapture) {
    w.capture = waves_->item(row, column)->checkState() == Qt::Checked;
  } else {
    return;
  }
  r.edited = true;
  fillRow(row);
  updateStates();
  schedule();
}

void PotentialPanel::modeChosen(int row, int index) { setMode(row, (Mode)index); }

void PotentialPanel::setMode(int row, Mode mode) {
  if (row < 0 || row >= rows_.size()) return;
  Row &r = rows_[row];
  front::PotentialWave &w = r.wave;
  r.line = mode != Central;
  w.background = mode == ExteriorInterior;
  w.boundStateOnly = mode == PointCoulomb;
  if (mode == ExteriorInterior && w.boundary.empty()) w.boundary = "0";  // the background needs its B (start: 0)
  if (mode != ExteriorInterior) {
    w.capture = false;
    w.subtract.clear();
  }
  r.edited = true;
  fillRow(row);
  waves_->setCurrentCell(row, kMode);
  updateStates();
  schedule();
}

void PotentialPanel::captureToggled(int row, bool on) {
  if (row < 0 || row >= rows_.size()) return;
  rows_[row].wave.capture = on;
  rows_[row].edited = true;
  schedule();
}

int PotentialPanel::addWave(int l, const QString &J) {
  Row r;
  r.line = true;
  r.edited = true;
  r.wave.l = l;
  r.wave.J = J.toStdString();
  if (!J.isEmpty()) {
    const QStringList p = J.split('/');
    r.wave.twoJ = p.size() == 2 ? p[0].toInt() * 2 / p[1].toInt() : 2 * J.toInt();
  }
  // after the rows of its l; a central row of that l (every J) gives way to it when it names every J
  int at = rows_.size();
  for (int i = 0; i < rows_.size(); i++) {
    if (rows_[i].wave.l == l && !rows_[i].line && J.isEmpty()) {
      rows_[i] = r;
      fillTable();
      waves_->setCurrentCell(i, kMode);
      schedule();
      return i;
    }
    if (rows_[i].wave.l <= l) at = i + 1;
  }
  rows_.insert(at, r);
  fillTable();
  waves_->setCurrentCell(at, kMode);
  schedule();
  return at;
}

void PotentialPanel::addWaveAsked() {
  QDialog d(this);
  d.setWindowTitle(tr("Add a wave line"));
  QFormLayout *f = new QFormLayout(&d);
  QSpinBox *l = new QSpinBox;
  l->setRange(0, 20);
  QComboBox *j = new QComboBox;
  auto fillJ = [this, l, j]() {
    j->clear();
    j->addItem(tr("every J of l"), QString());
    // the J the pair's channel spins make with l (the <levels> row's spins)
    const int c1 = front::LevelColumnIndex("j1"), c2 = front::LevelColumnIndex("j2"), ck = front::LevelColumnIndex("ir");
    for (const front::LevelRow &row : holder_->project().Levels()) {
      if ((int)row.tokens.size() <= std::max(c1, std::max(c2, ck)) || row.tokens[ck] != std::to_string(pair_)) continue;
      const int a = (int)std::lround(2 * std::stod(row.tokens[c1])), b = (int)std::lround(2 * std::stod(row.tokens[c2]));
      QList<int> seen;
      for (int s = std::abs(a - b); s <= a + b; s += 2)
        for (int tj = std::abs(2 * l->value() - s); tj <= 2 * l->value() + s; tj += 2)
          if (!seen.contains(tj)) {
            seen << tj;
            j->addItem(Spin(tj), Spin(tj));
          }
      break;
    }
  };
  fillJ();
  connect(l, QOverload<int>::of(&QSpinBox::valueChanged), &d, fillJ);
  f->addRow(tr("l"), l);
  f->addRow(tr("J"), j);
  QDialogButtonBox *b = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
  connect(b, &QDialogButtonBox::accepted, &d, &QDialog::accept);
  connect(b, &QDialogButtonBox::rejected, &d, &QDialog::reject);
  f->addRow(b);
  if (d.exec() == QDialog::Accepted) addWave(l->value(), j->currentData().toString());
}

void PotentialPanel::removeWaveLine() {
  const int row = currentRow();
  if (row >= 0) setMode(row, Central);
}

void PotentialPanel::revert() { reload(); }

bool PotentialPanel::apply(QString &why) {
  if (!hasRecord_ || !on_) {
    why = tr("no potential switched on to write");
    return false;
  }
  filling_ = true;
  const bool ok = holder_->setRecord(kSection, editedRecord(), why);
  filling_ = false;
  if (ok) reload();
  return ok;
}

void PotentialPanel::switchOn(bool on) {
  if (filling_ || !hasRecord_ || on == on_) return;
  QString why;
  filling_ = true;
  const bool ok = holder_->setRecordOn(kSection, Q(record_.name), on, why);
  filling_ = false;
  if (!ok) QMessageBox::warning(this, tr("Potential"), why);
  reload();
}

void PotentialPanel::createPotential() {
  if (pair_ <= 0) return;
  QString label;
  for (const front::PairChoice &c : front::Pairs(holder_->project()))
    if (c.key == pair_) label = Q(c.label);
  front::Record r;
  r.kind = "potential";
  QString name = "potential_" + QString(label).replace(QRegularExpression("[^A-Za-z0-9]"), "");
  r.name = name.toStdString();
  r.fields.push_back({"pair", label.toStdString(), ""});
  RecordDialog dialog(holder_, kSection, {"potential"},
                      [this](const front::Record &rec, QString &why) { return holder_->setRecord(kSection, rec, why); }, this);
  dialog.setWindowTitle(tr("A potential for pair %1").arg(label));
  dialog.setRecord(r);
  dialog.exec();
  reload();
}

void PotentialPanel::editRecord() {
  if (!hasRecord_ || !on_) return;
  const front::Record shown = editedRecord();
  RecordDialog dialog(holder_, kSection, {"potential"},
                      [this, shown](const front::Record &rec, QString &why) {
                        front::Record same = rec;
                        same.name = shown.name;
                        return holder_->setRecord(kSection, same, why);
                      },
                      this);
  dialog.setWindowTitle(tr("Potential %1").arg(Q(shown.name)));
  dialog.setRecord(shown);
  dialog.exec();
  reload();
}

void PotentialPanel::levelsChanged() {
  if (!hasRecord_) return;
  updateStates();
  showMessages();
  schedule();
}

void PotentialPanel::schedule() {
  if (hasRecord_) debounce_->start();
}

void PotentialPanel::refresh() {
  const int row = currentRow();
  if (!hasRecord_ || row < 0) return;
  QString commitWhy;
  if (hooks_.commitLevels) hooks_.commitLevels(commitWhy);
  updateStates();
  const std::string text = holder_->project().Text(), path = holder_->project().Path();
  const front::Record record = editedRecord();
  const int key = pair_, l = rows_[row].wave.l, twoJ = requestTwoJ(rows_[row]);
  const double eMax = eMax_->value();
  const bool validate = on_;
  const int ticket = ++ticket_;
  auto report = std::make_shared<front::PotentialWaveReport>();
  auto ok = std::make_shared<bool>(false);
  auto why = std::make_shared<std::string>();
  auto messages = std::make_shared<QStringList>();
  worker_->submit(
      "wave",
      [=](PhysicsSessionState &, const std::atomic<bool> &) {
        const front::Project copy = front::Project::FromText(text, path);
        *ok = front::PotentialWaveOf(copy, record, key, l, twoJ, -1000.0, eMax, *report, *why);
        if (!validate) return;
        // the engine's diagnostics of the record as edited (Validate of a copy with it written)
        front::Project v = front::Project::FromText(text, path);
        std::string w;
        if (!v.SetRecord(kSection, record, w)) {
          *messages << QString("<font color='#b00020'>%1</font>").arg(Q(w).toHtmlEscaped());
          return;
        }
        for (const front::Diagnostic &d : v.Validate().diagnostics)
          if (d.record == record.name)
            *messages << QString(d.severity == "error" ? "<font color='#b00020'>%1</font>" : "%1")
                             .arg(Q(d.code + ": " + d.message).toHtmlEscaped());
      },
      [=]() {
        if (ticket != ticket_) return;
        reportOk_ = *ok;
        report_ = *report;
        reportWhy_ = Q(*why);
        reportRow_ = row;
        recordMessages_ = *messages;
        showReport();
        showMessages();
        emit computed();
      });
}

void PotentialPanel::showMessages() {
  QStringList lines;
  if (!hasRecord_) {
    if (pair_ > 0)
      lines << tr("This pair has no potential: its channels take the point-Coulomb (or free) functions.  \"Add a "
                  "potential...\" makes one.");
  } else if (!on_) {
    lines << tr("The potential is switched off: its record is kept in the file as comments and does nothing.  Switch it "
                "on to edit it.");
  }
  lines << recordMessages_;
  if (reportOk_) {
    if (!report_.radiusError.empty()) lines << QString("<font color='#b00020'>%1</font>").arg(Q(report_.radiusError).toHtmlEscaped());
    if (!report_.subtractError.empty())
      lines << QString("<font color='#b00020'>%1</font>").arg(Q(report_.subtractError).toHtmlEscaped());
  } else if (!reportWhy_.isEmpty()) {
    lines << QString("<font color='#b00020'>%1</font>").arg(reportWhy_.toHtmlEscaped());
  }
  // the trap: a wave whose J^pi is not in the calculation does nothing
  if (hasRecord_ && on_) {
    QStringList left;
    for (int i = 0; i < rows_.size(); i++) {
      if (!rows_[i].line) continue;
      for (const front::PartialWaveState &s : waveStates(i))
        if (!s.included) left << tr("l = %1 (%2)").arg(rows_[i].wave.l).arg(jpiText(s.twoJ, s.parity));
    }
    if (!left.isEmpty())
      lines << tr("<font color='#b00020'><b>Left out of the calculation:</b> %1.</font> No active level carries these "
                  "J<sup>π</sup> and hard-sphere scattering in every J<sup>π</sup> (includeAllJPi) is off or does not "
                  "reach their l, so the engine makes no channel there and the potential of these waves does nothing.  "
                  "Add a level of each (the Levels tab's \"Partial waves...\" adds blank ones) or switch the option on.")
                   .arg(left.join(", "));
  }
  messages_->setText(lines.join("<br>"));
  messages_->setVisible(!lines.isEmpty());
}

void PotentialPanel::showReport() {
  poles_->setRowCount(0);
  levels_->setRowCount(0);
  bindCombo_->clear();
  const bool ok = reportOk_ && reportRow_ >= 0 && reportRow_ < rows_.size();
  lengthButton_->setEnabled(ok && report_.effectiveRange);
  bindButton_->setEnabled(false);
  if (!ok) {
    waveTitle_->setText(hasRecord_ ? tr("<i>Select a wave.</i>") : QString());
    readOut_->setText(QString());
    tuning_->setText(QString());
    useDepth_->setEnabled(false);
    return;
  }
  const front::PotentialWaveReport &r = report_;
  const QString jpi = jpiText(r.twoJ >= 0 ? r.twoJ : 2 * r.l + 1, r.parity);
  waveTitle_->setText(tr("<b>Wave l = %1, J<sup>π</sup> = %2</b>%3: V0 = %4 MeV%5, a_c = %6 fm%7%8")
                          .arg(r.l)
                          .arg(r.twoJ >= 0 ? jpi : tr("%1 (parity)").arg(r.parity > 0 ? "+" : "-"))
                          .arg(rows_[reportRow_].wave.J.empty() && r.twoJ >= 0 ? tr(" (shown for J = %1)").arg(Spin(r.twoJ)) : QString())
                          .arg(Num(r.depth))
                          .arg(r.width > 0 ? tr(", r0 = %1 fm").arg(Num(r.width)) : QString())
                          .arg(Num(r.radius))
                          .arg(r.background ? tr(", B_bg = %1").arg(Num(r.boundary)) : tr(" (no background: the poles are "
                                                                                          "those B = 0 would give)"))
                          .arg(r.radiusNote.empty() ? QString() : " — " + Q(r.radiusNote)));
  poles_->setRowCount((int)r.poles.size());
  for (int i = 0; i < (int)r.poles.size(); i++) {
    const front::PotentialPole &p = r.poles[i];
    poles_->setItem(i, 0, new QTableWidgetItem(QString::number(p.nodes)));
    poles_->setItem(i, 1, new QTableWidgetItem(Num(p.energy, 7)));
    poles_->setItem(i, 2, new QTableWidgetItem(Num(p.excitation, 7)));
    poles_->setItem(i, 3, new QTableWidgetItem(Num(p.gamma2, 5)));
    poles_->setItem(i, 4, new QTableWidgetItem(Num(p.wignerRatio, 4)));
    poles_->setItem(i, 5, new QTableWidgetItem(p.subtracted ? tr("subtracted (a level replaces it)")
                                                            : (p.energy < 0 ? tr("below threshold") : QString())));
    if (p.subtracted)
      for (int c = 0; c < 6; c++) {
        QFont f = poles_->item(i, c)->font();
        f.setBold(true);
        poles_->item(i, c)->setFont(f);
      }
  }
  levels_->setRowCount((int)r.levels.size());
  for (int i = 0; i < (int)r.levels.size(); i++) {
    const front::PotentialLevel &lv = r.levels[i];
    const QString unit = lv.kind == "closed" ? tr("ANC %1 fm^-1/2") : lv.kind == "open" ? tr("Γ %1 eV")
                         : lv.kind == "amplitude" ? tr("γ %1 MeV^1/2") : tr("%1");
    levels_->setItem(i, 0, new QTableWidgetItem(jpiText((int)std::lround(2 * lv.J), lv.parity) + tr(" s = %1").arg(lv.s)));
    levels_->setItem(i, 1, new QTableWidgetItem(Num(lv.excitation, 7)));
    levels_->setItem(i, 2, new QTableWidgetItem(Num(lv.energy, 7)));
    levels_->setItem(i, 3, new QTableWidgetItem(unit.arg(Num(lv.value, 7))));
    levels_->setItem(i, 4, new QTableWidgetItem(lv.active ? tr("yes") : tr("no")));
    if (lv.energy < 0 && lv.active)
      bindCombo_->addItem(tr("%1 at E_x %2 MeV (E_c %3 MeV), %4").arg(jpiText((int)std::lround(2 * lv.J), lv.parity), Num(lv.excitation),
                                                                      Num(lv.energy), unit.arg(Num(lv.value))),
                          QVariant::fromValue(i));
  }
  bindButton_->setEnabled(bindCombo_->count() > 0);
  if (bindCombo_->count() == 0) bindCombo_->addItem(tr("(no active level of this J\u03c0 below threshold)"));
  // the potential alone
  QString a;
  if (r.effectiveRange) {
    const QString unitA = r.l == 0 ? tr("fm") : tr("fm^%1").arg(2 * r.l + 1), unitR = r.l == 0 ? tr("fm") : tr("fm^%1").arg(1 - 2 * r.l);
    a = tr("a = <b>%1</b> ± %2 %3, r = %4 ± %5 %6 (the effective-range expansion of its own phase shifts, %7-%8 MeV)")
            .arg(Num(r.a, 7), Num(r.aSpread, 2), unitA, Num(r.r, 5), Num(r.rSpread, 2), unitR, Num(r.eLow, 3), Num(r.eHigh, 3));
  } else {
    a = tr("a, r: %1").arg(Q(r.effectiveRangeWhy).toHtmlEscaped());
  }
  readOut_->setText(tr("%1<br>%2 bound state(s) in this wave (forbidden ones included); the zero-energy wave has %3 "
                       "node(s) inside the matching radius.").arg(a).arg(r.boundStates).arg(r.zeroEnergyNodes));
  if (depthRow_ != reportRow_) {
    tuning_->setText(QString());
    useDepth_->setEnabled(false);
  }
}

void PotentialPanel::findLength() {
  const int row = currentRow();
  bool okNumber = false;
  const double a = lengthEdit_->text().trimmed().toDouble(&okNumber);
  if (row < 0 || !okNumber) {
    tuning_->setText(tr("<font color='#b00020'>Type the scattering length (a number, fm).</font>"));
    return;
  }
  const std::string text = holder_->project().Text(), path = holder_->project().Path();
  const front::Record record = editedRecord();
  const int key = pair_, l = rows_[row].wave.l, twoJ = requestTwoJ(rows_[row]);
  auto out = std::make_shared<front::PotentialDepth>();
  auto ok = std::make_shared<bool>(false);
  auto why = std::make_shared<std::string>();
  tuning_->setText(tr("searching..."));
  worker_->submit(
      "depth",
      [=](PhysicsSessionState &, const std::atomic<bool> &) {
        *ok = front::PotentialDepthForScatteringLength(front::Project::FromText(text, path), record, key, l, twoJ, a, *out, *why);
      },
      [=]() {
        depth_ = *out;
        depthRow_ = row;
        if (!*ok || !out->found) {
          tuning_->setText(QString("<font color='#b00020'>%1</font>").arg(Q(*ok ? out->why : *why).toHtmlEscaped()));
          useDepth_->setEnabled(false);
        } else {
          tuning_->setText(tr("V0 = <b>%1 MeV</b> gives a = %2 (r = %3).  Branch kept: %4 bound state(s), as now; the "
                              "zero-energy wave has %5 node(s).")
                               .arg(Num(out->depth, 10), Num(out->a, 7), Num(out->r, 5))
                               .arg(out->boundStates)
                               .arg(out->nodes));
          useDepth_->setEnabled(true);
        }
        emit computed();
      });
}

void PotentialPanel::findBinding() {
  const int row = currentRow();
  if (row < 0 || !reportOk_ || bindCombo_->currentData().isNull()) return;
  const int index = bindCombo_->currentData().toInt();
  if (index < 0 || index >= (int)report_.levels.size()) return;
  const front::PotentialLevel level = report_.levels[index];
  const std::string text = holder_->project().Text(), path = holder_->project().Path();
  const front::Record record = editedRecord();
  const int key = pair_, l = rows_[row].wave.l, twoJ = requestTwoJ(rows_[row]);
  auto out = std::make_shared<front::PotentialDepth>();
  auto ok = std::make_shared<bool>(false);
  auto why = std::make_shared<std::string>();
  tuning_->setText(tr("searching..."));
  worker_->submit(
      "depth",
      [=](PhysicsSessionState &, const std::atomic<bool> &) {
        *ok = front::PotentialDepthForBinding(front::Project::FromText(text, path), record, key, l,
                                              twoJ >= 0 ? twoJ : (int)std::lround(2 * level.J), level.energy, *out, *why);
      },
      [=]() {
        depth_ = *out;
        depthRow_ = row;
        if (!*ok || !out->found) {
          tuning_->setText(QString("<font color='#b00020'>%1</font>").arg(Q(*ok ? out->why : *why).toHtmlEscaped()));
          useDepth_->setEnabled(false);
        } else {
          const QString levelAnc = level.kind == "closed" ? tr("%1 fm^-1/2 (as written in &lt;levels&gt;)").arg(Num(level.value, 6))
                                                          : tr("not an ANC (%1)").arg(Q(level.kind));
          tuning_->setText(tr("V0 = <b>%1 MeV</b> binds the wave at E_c = %2 MeV (a state with %3 node(s); %4 bound "
                              "state(s) in all).  The potential's own ANC: <b>%5 fm^-1/2</b>; the level's: %6.")
                               .arg(Num(out->depth, 10), Num(out->energy, 9))
                               .arg(out->nodes)
                               .arg(out->boundStates)
                               .arg(Num(out->anc, 6), levelAnc));
          useDepth_->setEnabled(true);
        }
        emit computed();
      });
}

void PotentialPanel::useDepth() {
  if (depthRow_ < 0 || depthRow_ >= rows_.size() || !depth_.found) return;
  Row &r = rows_[depthRow_];
  r.line = true;
  r.wave.depth = QString::number(depth_.depth, 'g', 10).toStdString();
  r.wave.depthUnit = front::PotentialWave().depthUnit;  // the facade's unit of the depth found
  r.edited = true;
  fillRow(depthRow_);
  waves_->setCurrentCell(depthRow_, kDepth);
  useDepth_->setEnabled(false);
  updateStates();
  schedule();
}
