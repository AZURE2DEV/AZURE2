#include "MassCheckDialog.h"

#include <QCheckBox>
#include <QDialogButtonBox>
#include <QHeaderView>
#include <QLabel>
#include <QPushButton>
#include <QTableWidget>
#include <QVBoxLayout>
#include <QBrush>
#include <QColor>
#include <algorithm>
#include <cmath>

#include "MassTable.h"

namespace {

QString full(double v) { return QString::number(v, 'g', 15); }

// "1.00728 -> 1.007825031898" when the table differs, the value alone otherwise
QString massCell(double now, const MassEntry *e) {
  if (!e) return full(now);
  if (now == e->atomicMass) return full(now);
  return QString::fromUtf8("%1 → %2").arg(full(now)).arg(full(e->atomicMass));
}

enum Column { Fill = 0, Pair, Light, Heavy, S, STable, DS, NumColumns };

}  // namespace

MassCheckDialog::MassCheckDialog(const PairsModel *model, QWidget *parent) : QDialog(parent) {
  setWindowTitle(tr("Masses and separation energies from %1").arg(MassTable::version()));
  pairs_ = model->getPairs();
  const MassTable &mt = MassTable::instance();

  QLabel *intro = new QLabel(tr("Each pair against the %1 atomic mass table.  A particle is identified by its charge and its mass "
                                "rounded to a mass number.  Ticked pairs get the table's atomic masses and, with the box below, "
                                "the table's separation energy S = m1 + m2 - M.  Moving S moves that pair's threshold: compare "
                                "the fit before and after.  A gamma or beta pair only gets its heavy particle's mass.")
                                 .arg(MassTable::version()));
  intro->setWordWrap(true);

  table_ = new QTableWidget(pairs_.size(), NumColumns, this);
  table_->setObjectName("massCheckTable");
  table_->setHorizontalHeaderLabels({tr("Fill"), tr("Pair"), tr("Light mass [u]"), tr("Heavy mass [u]"), tr("S [MeV]"),
                                     tr("S from table [MeV]"), tr("S - table [keV]")});
  table_->verticalHeader()->setVisible(false);
  table_->setEditTriggers(QAbstractItemView::NoEditTriggers);
  table_->setSelectionMode(QAbstractItemView::NoSelection);

  for (int row = 0; row < pairs_.size(); row++) {
    const PairsData &p = pairs_.at(row);
    const PairMassCheck c = mt.check(p);
    QTableWidgetItem *fill = new QTableWidgetItem;
    fill->setFlags(c.fillable ? (Qt::ItemIsUserCheckable | Qt::ItemIsEnabled) : Qt::NoItemFlags);
    const bool differs = c.massOffset > 0.0 || (c.hasS && c.dSkeV != 0.0);
    fill->setCheckState(c.fillable && differs ? Qt::Checked : Qt::Unchecked);
    if (!c.fillable) fill->setToolTip(tr("%1 has no nuclide for this pair's charge and mass number.").arg(MassTable::version()));
    table_->setItem(row, Fill, fill);
    table_->setItem(row, Pair, new QTableWidgetItem(QString("%1: %2").arg(row + 1).arg(model->getPairLabel(p))));
    table_->setItem(row, Light, new QTableWidgetItem(p.pairType == 0 ? massCell(p.lightM, c.light) : QString("-")));
    table_->setItem(row, Heavy, new QTableWidgetItem(massCell(p.heavyM, c.heavy)));
    table_->setItem(row, S, new QTableWidgetItem(full(p.seperationEnergy)));
    table_->setItem(row, STable, new QTableWidgetItem(c.hasS ? full(c.sTable) : QString("-")));
    QTableWidgetItem *ds = new QTableWidgetItem(c.hasS ? QString::number(c.dSkeV, 'f', 3) : QString("-"));
    if (c.hasS && std::fabs(c.dSkeV) > MassTable::kToleranceKeV) ds->setForeground(QBrush(QColor(214, 80, 0)));
    table_->setItem(row, DS, ds);
  }
  table_->resizeColumnsToContents();
  table_->horizontalHeader()->setStretchLastSection(true);
  // show every column and up to twelve rows without scrolling
  const int frame = 2 * table_->frameWidth() + 4;
  const int width = table_->horizontalHeader()->length() + frame;
  const int rows = std::min(12, static_cast<int>(pairs_.size()));
  const int height = table_->horizontalHeader()->sizeHint().height() + rows * table_->verticalHeader()->defaultSectionSize() + frame;
  table_->setMinimumSize(std::min(width, 1300), height);

  separationCheck_ = new QCheckBox(tr("Also set the separation energies of particle pairs from the table"));
  separationCheck_->setObjectName("massCheckSeparation");
  separationCheck_->setChecked(true);

  summary_ = new QLabel;
  summary_->setWordWrap(true);

  QDialogButtonBox *buttons = new QDialogButtonBox(QDialogButtonBox::Cancel);
  QPushButton *apply = buttons->addButton(tr("Fill ticked pairs"), QDialogButtonBox::AcceptRole);
  apply->setObjectName("massCheckApply");
  connect(buttons, SIGNAL(accepted()), this, SLOT(accept()));
  connect(buttons, SIGNAL(rejected()), this, SLOT(reject()));

  auto updateSummary = [this, apply]() {
    const int n = filledPairs().size();
    summary_->setText(n == 0 ? tr("Nothing to change for the ticked pairs.")
                             : tr("%n pair(s) will change.", "", n));
    apply->setEnabled(n > 0);
  };
  connect(table_, &QTableWidget::itemChanged, this, updateSummary);
  connect(separationCheck_, &QCheckBox::toggled, this, updateSummary);
  updateSummary();

  QVBoxLayout *layout = new QVBoxLayout;
  layout->addWidget(intro);
  layout->addWidget(table_);
  layout->addWidget(separationCheck_);
  layout->addWidget(summary_);
  layout->addWidget(buttons);
  setLayout(layout);
}

PairsData MassCheckDialog::filled(int row, bool separation) const {
  PairsData p = pairs_.at(row);
  const MassTable &mt = MassTable::instance();
  const PairMassCheck c = mt.check(p);
  if (!c.fillable) return p;
  if (c.heavy) p.heavyM = c.heavy->atomicMass;
  if (p.pairType == 0) {
    p.lightM = c.light->atomicMass;
    if (separation && c.hasS) p.seperationEnergy = c.sTable;
  }
  return p;
}

QList<QPair<int, PairsData>> MassCheckDialog::filledPairs() const {
  QList<QPair<int, PairsData>> out;
  for (int row = 0; row < pairs_.size(); row++) {
    if (table_->item(row, Fill)->checkState() != Qt::Checked) continue;
    const PairsData now = pairs_.at(row);
    const PairsData p = filled(row, separationCheck_->isChecked());
    if (p.lightM != now.lightM || p.heavyM != now.heavyM || p.seperationEnergy != now.seperationEnergy) out.append({row, p});
  }
  return out;
}
