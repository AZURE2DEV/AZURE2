#include "MixtureEditor.h"
#include "DialogScrollArea.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QHeaderView>
#include <QFileDialog>
#include <QDir>

MixtureEditor::MixtureEditor(QWidget *parent) : QWidget(parent) {
  enableCheck = new QCheckBox(tr("Natural-abundance or compound target (mixture of nuclei)"));
  connect(enableCheck, SIGNAL(stateChanged(int)), this, SLOT(enableChanged(int)));

  box_ = new QGroupBox(tr("Mixture: other nuclei frozen as cross-section tables"));
  box_->setVisible(false);

  selfWeightText = new QLineEdit;
  selfWeightText->setText("1.0");
  selfWeightText->setMaximumWidth(100);
  connect(selfWeightText, SIGNAL(textChanged(QString)), this, SLOT(updateSum()));

  tablesTable = new QTableWidget(0, 2);
  tablesTable->setHorizontalHeaderLabels(QStringList() << tr("Table file  (E_lab [MeV], \xcf\x83 [b])") << tr("Weight"));
  tablesTable->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
  tablesTable->horizontalHeader()->setSectionResizeMode(1, QHeaderView::ResizeToContents);
  tablesTable->verticalHeader()->setVisible(false);
  tablesTable->setSelectionBehavior(QAbstractItemView::SelectRows);
  tablesTable->setMinimumHeight(110);   // three rows and the header
  tablesTable->setMaximumHeight(160);
  connect(tablesTable, SIGNAL(cellChanged(int, int)), this, SLOT(updateSum()));

  addTableButton = new QPushButton(tr("Add Table..."));
  connect(addTableButton, SIGNAL(clicked()), this, SLOT(addTable()));
  removeTableButton = new QPushButton(tr("Remove Table"));
  connect(removeTableButton, SIGNAL(clicked()), this, SLOT(removeTable()));

  sumLabel = new QLabel;
  rawLabel = new QLabel;
  rawLabel->setWordWrap(true);
  rawLabel->setVisible(false);

  QLabel *hint = new QLabel(tr("Model = w_self \xc3\x97 this segment + \xce\xa3 weight \xc3\x97 table (each table folded with the "
                               "segment's experimental effects); abundances sum to 1, a compound counts atoms per "
                               "atom of this nucleus; files are relative to the .azr."));
  hint->setWordWrap(true);
  hint->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Minimum);

  QGridLayout *grid = new QGridLayout;
  grid->addWidget(new QLabel(tr("Abundance of this nucleus (w_self):")), 0, 0, Qt::AlignRight);
  grid->addWidget(selfWeightText, 0, 1, Qt::AlignLeft);
  grid->addWidget(tablesTable, 1, 0, 1, 3);
  QHBoxLayout *buttons = new QHBoxLayout;
  buttons->addWidget(addTableButton);
  buttons->addWidget(removeTableButton);
  buttons->addStretch(1);
  buttons->addWidget(sumLabel);
  grid->addLayout(buttons, 2, 0, 1, 3);
  grid->addWidget(hint, 3, 0, 1, 3);
  grid->setColumnStretch(2, 1);
  box_->setLayout(grid);

  QVBoxLayout *layout = new QVBoxLayout;
  layout->setContentsMargins(0, 0, 0, 0);
  layout->addWidget(enableCheck);
  // Outside the box: it must show while the (unreadable) block keeps the box closed.
  layout->addWidget(rawLabel);
  layout->addWidget(box_);
  setLayout(layout);
  updateSum();
}

// ---------------------------------------------------------------------------
// token helpers

QStringList MixtureEditor::splitTokens(const QString &text) {
  QStringList out;
  int i = 0;
  const int n = text.size();
  while (i < n) {
    while (i < n && text[i].isSpace()) i++;
    if (i >= n) break;
    if (text[i] == '"') {
      // As the engine's ReadPathToken: a token that opens with a quote runs to
      // the closing quote, spaces included.
      int j = text.indexOf('"', i + 1);
      if (j < 0) j = n;
      out << text.mid(i + 1, j - i - 1);
      i = j + 1;
    } else {
      int j = i;
      while (j < n && !text[j].isSpace()) j++;
      out << text.mid(i, j - i);
      i = j;
    }
  }
  return out;
}

static QString quoteIfNeeded(const QString &token) {
  if (token.contains(' ') || token.contains('\t')) return "\"" + token + "\"";
  return token;
}

static QString formatNumber(double v) {
  return QString::number(v, 'g', 10);
}

bool MixtureEditor::parseBlock(const QString &extraTokens, bool &present, double &selfWeight,
                               QList<Table> &tables, QString &other, QString *error) {
  present = false;
  selfWeight = 1.0;
  tables.clear();
  other = extraTokens.trimmed();
  const QStringList tok = splitTokens(extraTokens);
  const int start = tok.indexOf("mixture");
  if (start < 0) return true;
  present = true;
  if (start + 2 >= tok.size()) {
    if (error) *error = QObject::tr("the mixture block ends before its table count");
    return false;
  }
  bool ok = false;
  selfWeight = tok[start + 1].toDouble(&ok);
  if (!ok) {
    if (error) *error = QObject::tr("'%1' is not a number (w_self)").arg(tok[start + 1]);
    return false;
  }
  const int count = tok[start + 2].toInt(&ok);
  if (!ok || count < 0) {
    if (error) *error = QObject::tr("'%1' is not a table count").arg(tok[start + 2]);
    return false;
  }
  int i = start + 3;
  for (int k = 0; k < count; k++, i += 2) {
    if (i + 1 >= tok.size()) {
      if (error) *error = QObject::tr("the mixture block names %1 tables but ends after %2").arg(count).arg(k);
      return false;
    }
    Table t;
    t.file = tok[i];
    t.weight = tok[i + 1].toDouble(&ok);
    if (!ok) {
      if (error) *error = QObject::tr("'%1' is not a number (weight of %2)").arg(tok[i + 1], t.file);
      return false;
    }
    tables << t;
  }
  QStringList rest;
  for (int k = 0; k < start; k++) rest << quoteIfNeeded(tok[k]);
  for (int k = i; k < tok.size(); k++) rest << quoteIfNeeded(tok[k]);
  other = rest.join(" ");
  return true;
}

QString MixtureEditor::formatBlock(double selfWeight, const QList<Table> &tables) {
  QStringList out;
  out << "mixture" << formatNumber(selfWeight) << QString::number(tables.size());
  for (const Table &t : tables) out << "\"" + t.file + "\"" << formatNumber(t.weight);
  return out.join(" ");
}

// ---------------------------------------------------------------------------
// widget state

void MixtureEditor::setTokens(const QString &extraTokens) {
  bool present = false;
  double w = 1.0;
  QList<Table> tables;
  QString other, error;
  const bool ok = parseBlock(extraTokens, present, w, tables, other, &error);
  tablesTable->blockSignals(true);
  tablesTable->setRowCount(0);
  raw_.clear();
  rawLabel->setVisible(false);
  if (!ok) {
    // Keep the whole text verbatim: a block the engine would refuse is still
    // the user's text, and dropping it on an unrelated edit would be worse.
    other_ = QString();
    raw_ = extraTokens.trimmed();
    rawLabel->setText(tr("This row carries a mixture block that could not be read (%1). It is kept as written:\n%2")
                          .arg(error, raw_));
    rawLabel->setVisible(true);
    selfWeightText->setText("1.0");
    tablesTable->blockSignals(false);
    enableCheck->setChecked(false);
    enableCheck->setEnabled(false);
    updateSum();
    return;
  }
  enableCheck->setEnabled(true);
  other_ = other;
  selfWeightText->setText(formatNumber(w));
  for (const Table &t : tables) {
    const int row = tablesTable->rowCount();
    tablesTable->insertRow(row);
    tablesTable->setItem(row, 0, new QTableWidgetItem(t.file));
    tablesTable->setItem(row, 1, new QTableWidgetItem(formatNumber(t.weight)));
  }
  tablesTable->blockSignals(false);
  enableCheck->setChecked(present);
  updateSum();
}

QList<MixtureEditor::Table> MixtureEditor::currentTables(bool *ok) const {
  QList<Table> tables;
  bool good = true;
  for (int row = 0; row < tablesTable->rowCount(); row++) {
    Table t;
    const QTableWidgetItem *f = tablesTable->item(row, 0);
    const QTableWidgetItem *w = tablesTable->item(row, 1);
    t.file = f ? f->text().trimmed() : QString();
    bool numOk = false;
    t.weight = w ? w->text().trimmed().toDouble(&numOk) : 0.0;
    if (t.file.isEmpty() || !numOk) good = false;
    tables << t;
  }
  if (ok) *ok = good;
  return tables;
}

QString MixtureEditor::block() const {
  if (!isMixtureEnabled()) return QString();
  bool ok = false;
  const double w = selfWeightText->text().trimmed().toDouble(&ok);
  return formatBlock(ok ? w : 1.0, currentTables());
}

QString MixtureEditor::tokens() const {
  if (!raw_.isEmpty()) return raw_;
  const QString b = block();
  if (b.isEmpty()) return other_;
  if (other_.isEmpty()) return b;
  return other_ + " " + b;
}

bool MixtureEditor::isMixtureEnabled() const {
  return enableCheck->isEnabled() && enableCheck->isChecked();
}

bool MixtureEditor::validate(QString *error) const {
  if (!isMixtureEnabled()) return true;
  bool ok = false;
  const double w = selfWeightText->text().trimmed().toDouble(&ok);
  if (!ok || w < 0.0) {
    if (error) *error = tr("The abundance of this nucleus (w_self) must be a number \xe2\x89\xa5 0.");
    return false;
  }
  for (int row = 0; row < tablesTable->rowCount(); row++) {
    const QTableWidgetItem *f = tablesTable->item(row, 0);
    const QTableWidgetItem *wi = tablesTable->item(row, 1);
    if (!f || f->text().trimmed().isEmpty()) {
      if (error) *error = tr("Table %1 has no file.").arg(row + 1);
      return false;
    }
    bool numOk = false;
    const double tw = wi ? wi->text().trimmed().toDouble(&numOk) : 0.0;
    if (!numOk || tw < 0.0) {
      if (error) *error = tr("The weight of table %1 must be a number \xe2\x89\xa5 0.").arg(row + 1);
      return false;
    }
  }
  return true;
}

void MixtureEditor::enableChanged(int state) {
  box_->setVisible(state == Qt::Checked);
  if (QDialog *d = qobject_cast<QDialog *>(window()))
    fitDialogToContent(d);
  else if (QWidget *w = window())
    w->adjustSize();
}

void MixtureEditor::addTable() {
  QString filename = QFileDialog::getOpenFileName(this, tr("Cross-section table (E_lab [MeV], sigma [b])"));
  const int row = tablesTable->rowCount();
  tablesTable->insertRow(row);
  tablesTable->setItem(row, 0, new QTableWidgetItem(filename.isEmpty() ? QString() : QDir::fromNativeSeparators(filename)));
  tablesTable->setItem(row, 1, new QTableWidgetItem("0.0"));
  tablesTable->setCurrentCell(row, filename.isEmpty() ? 0 : 1);
  updateSum();
}

void MixtureEditor::removeTable() {
  const int row = tablesTable->currentRow();
  if (row >= 0) tablesTable->removeRow(row);
  updateSum();
}

void MixtureEditor::updateSum() {
  bool ok = false;
  double sum = selfWeightText->text().trimmed().toDouble(&ok);
  if (!ok) sum = 0.0;
  bool tablesOk = true;
  for (const Table &t : currentTables(&tablesOk)) sum += t.weight;
  sumLabel->setText(tr("Sum of weights: %1").arg(formatNumber(sum)));
}
