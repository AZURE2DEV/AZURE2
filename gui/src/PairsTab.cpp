#include <cmath>
#include <QSet>
#include <QShortcut>
#include <QUndoStack>
#include <functional>
#include <QHeaderView>
#include <QGridLayout>
#include <QComboBox>
#include <QLineEdit>
#include <QCheckBox>
#include <QTextStream>
#include <QMessageBox>

#include <QInputDialog>
#include <QTabWidget>
#include "PairsTab.h"
#include "LevelsTab.h"
#include "LevelsViews.h"
#include "ProjectHolder.h"
#include <QLocale>
#include <QToolTip>
#include "RichTextDelegate.h"
#include "InfoDialog.h"
#include "MassCheckDialog.h"
#include "MassTable.h"
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
  // the pair type in words (particle, gamma, beta) is shown
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
  // a double click edits a number in place (excitation, separation energy, radius) and opens the dialog elsewhere
  connect(pairsView, &QAbstractItemView::doubleClicked, this, [this](const QModelIndex &i) {
    if (!PairsModel::editableInPlace(i.column())) editPair();
  });
  pairsView->setEditTriggers(QAbstractItemView::DoubleClicked | QAbstractItemView::EditKeyPressed |
                             QAbstractItemView::AnyKeyPressed);
  {
    auto *numbers = new CheckedEditDelegate(
        CheckedEditDelegate::Kind::Number,
        [this](const QModelIndex &index, const QVariant &value, QString &why) {
          return setPairNumber(index.row(), index.column(), value.toString(), why);
        },
        this);
    connect(numbers, &CheckedEditDelegate::refused, this, [this](const QModelIndex &index, const QString &why) {
      const QRect r = pairsView->visualRect(index);
      QToolTip::showText(pairsView->viewport()->mapToGlobal(r.bottomLeft()), why, pairsView, r, 6000);
    });
    for (int c : {10, 11, 12, (int)PairsModel::LightExcitationColumn}) pairsView->setItemDelegateForColumn(c, numbers);
  }

  for (int i = 0; i < pairsModel->columnCount(QModelIndex()); i++) pairsView->horizontalHeader()->setSectionResizeMode(i, QHeaderView::Stretch);

  addButton = new QPushButton(tr("+"));
  addButton->setMaximumSize(28, 28);
  connect(addButton, SIGNAL(clicked()), this, SLOT(addPair()));
  deleteButton = new QPushButton(tr("-"));
  deleteButton->setMaximumSize(28, 28);
  deleteButton->setEnabled(false);
  connect(deleteButton, SIGNAL(clicked()), this, SLOT(removePair()));
  duplicateButton = new QPushButton(tr("Duplicate"));
  duplicateButton->setObjectName("duplicatePair");
  duplicateButton->setToolTip(tr("A copy of the selected pair in the pair dialog, to change into another pair (Ctrl+D)"));
  duplicateButton->setEnabled(false);
  connect(duplicateButton, &QPushButton::clicked, this, &PairsTab::duplicatePair);
  {
    auto *s = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_D), pairsView);
    s->setContext(Qt::WidgetShortcut);
    connect(s, &QShortcut::activated, this, &PairsTab::duplicatePair);
  }

  /*
  mapper = new QSignalMapper(this);
  connect(mapper,SIGNAL(mapped(int)),this,SLOT(showInfo(int)));

  infoButton[0] = new QPushButton(this);
  infoButton[0]->setMaximumSize(28,28);
  infoButton[0]->setIcon(style()->standardIcon(QStyle::SP_MessageBoxInformation));
  mapper->setMapping(infoButton[0],0);
  connect(infoButton[0],SIGNAL(clicked()),mapper,SLOT(map()));
  */

  massButton = new QPushButton(tr("Masses from %1...").arg(MassTable::version()));
  massButton->setObjectName("massesButton");
  massButton->setToolTip(tr("Compare every pair's masses and separation energy with the %1 atomic mass table, and fill the "
                            "pairs you tick from it.").arg(MassTable::version()));
  connect(massButton, SIGNAL(clicked()), this, SLOT(checkMasses()));

  QGridLayout *buttonBox = new QGridLayout;
  buttonBox->addWidget(addButton, 0, 0);
  buttonBox->addWidget(deleteButton, 0, 1);
  buttonBox->addWidget(duplicateButton, 0, 3);
  buttonBox->addItem(new QSpacerItem(28, 28), 0, 2);
  buttonBox->addWidget(massButton, 0, 4);
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

  QVBoxLayout *layout = new QVBoxLayout();
  layout->addWidget(pairsView);
  layout->addLayout(buttonBox);
  extras = new QTabWidget;
  extras->setVisible(false);
  layout->addWidget(extras);
  layout->setStretch(0, 3);
  setLayout(layout);
}

void PairsTab::addExtraPanel(const QString &title, QWidget *panel) {
  extras->addTab(panel, title);
  extras->setVisible(true);
}

void PairsTab::reset() {
  pairsModel->removeRows(0, pairsModel->numPairs(), QModelIndex());
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
  AddPairDialog aDialog(this);
  if (aDialog.exec()) {
    if (!((pairsModel->getPairs()).size() == 0 && aDialog.pairTypeCombo->currentIndex() != 0)) {
      addPair(aDialog.pair(), pairsModel->numPairs(), false);
    } else {
      QMessageBox::information(this, tr("Pair Type Error"),
                               tr("The first pair must be a particle,particle pair."));
    }
  }
}

void PairsTab::duplicatePair() {
  const QModelIndex index = pairsView->currentIndex();
  if (!index.isValid() || index.row() >= pairsModel->numPairs()) return;
  AddPairDialog aDialog(this);
  aDialog.setWindowTitle(tr("Duplicate a Particle Pair"));
  aDialog.setPair(pairsModel->getPairs().at(index.row()));
  if (!aDialog.exec()) return;
  const PairsData pair = aDialog.pair();
  if (pairsModel->isPair(pair) != -1) {
    QMessageBox::information(this, tr("Pair not added"),
                             tr("The copy is the same pair as pair %1: change what makes it another pair (a state, a "
                                "radius, a spin) before accepting it.")
                                 .arg(pairsModel->isPair(pair) + 1));
    return;
  }
  addPair(pair, pairsModel->numPairs(), false);
  pairsView->setCurrentIndex(pairsModel->index(pairsModel->numPairs() - 1, index.column()));
}

void PairsTab::addPair(PairsData pair, int pairIndex, bool fromFile) {
  // A pair typed by the user that repeats another is refused; one read from the project is taken at its key (the
  // file may hold the same pair under two keys, LevelsTab::readLevels).
  if (fromFile || pairsModel->isPair(pair) == -1) {
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
    index = pairsModel->index(pairIndex, PairsModel::LightExcitationColumn, QModelIndex());
    pairsModel->setData(index, pair.lightExcitation, Qt::EditRole);

    pairsView->resizeRowsToContents();
    if (!fromFile) emit(pairAdded(pairIndex));
  } else {
    QMessageBox::information(this, tr("Duplicate Pair"),
                             tr("This pair already exists."));
  }
}

void PairsTab::removePairAt(int row) {
  const bool confirm = confirmRemoval_;
  confirmRemoval_ = false;  // the tests' removal: no question
  pairsView->selectionModel()->select(pairsModel->index(row, 0),
                                      QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
  removePair();
  confirmRemoval_ = confirm;
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
    const QStringList goes = pairConsequences(index.row());
    if (confirmRemoval_ && !goes.isEmpty() &&
        QMessageBox::question(this, tr("Remove Pair"),
                              tr("Remove pair %1 (%2)?  The pairs after it move up one key.").arg(index.row() + 1)
                                      .arg(PairsModel::plainText(pairsModel->getPairLabel(pairsModel->getPairs().at(index.row())))) +
                                  "\n\n" + goes.join("\n") + "\n\n" + tr("Undo with Ctrl+Z."),
                              QMessageBox::Yes | QMessageBox::Cancel, QMessageBox::Yes) != QMessageBox::Yes)
      return;
    pairsModel->removeRows(index.row(), 1, QModelIndex());
    emit(pairRemoved(index.row()));
  } else {
    QMessageBox::information(this, tr("Entrance Channel Error"),
                             tr("This delete will result in a forbidden particle,gamma entrance pair."));
  }
}

void PairsTab::editPair() {
  const QModelIndexList indexes = pairsView->selectionModel()->selectedRows();
  if (indexes.isEmpty()) return;
  const QModelIndex index = indexes[0];
  AddPairDialog aDialog(this);
  aDialog.setWindowTitle(tr("Edit a Particle Pair"));
  aDialog.setPair(pairsModel->getPairs().at(index.row()));
  if (aDialog.exec()) {
    if (!(index.row() == 0 && aDialog.pairTypeCombo->currentIndex() != 0)) {
      editPair(aDialog.pair(), index.row(), false);
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
  // compared the pair type with the stored mask before 2026-11-05, so an edited mask was dropped whenever the
  // pair type equalled the old mask (a particle pair, type 0, with mask 0)
  if (pair.ecMultMask != var.toInt()) pairsModel->setData(i, pair.ecMultMask, Qt::EditRole);
  i = pairsModel->index(pairIndex, PairsModel::LightExcitationColumn, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.lightExcitation != var.toDouble()) pairsModel->setData(i, pair.lightExcitation, Qt::EditRole);

  if (!fromFile) emit(pairEdited(pairIndex));
}

void PairsTab::checkMasses() {
  if (pairsModel->numPairs() == 0) {
    QMessageBox::information(this, tr("Masses"), tr("There are no pairs yet."));
    return;
  }
  if (!MassTable::instance().isLoaded()) {
    QMessageBox::warning(this, tr("Masses"), tr("The %1 mass table could not be read.").arg(MassTable::version()));
    return;
  }
  MassCheckDialog dialog(pairsModel, this);
  if (dialog.exec() != QDialog::Accepted) return;
  const QList<QPair<int, PairsData>> changes = dialog.filledPairs();
  for (const auto &change : changes) editPair(change.second, change.first, false);
  pairsView->viewport()->update();
}

void PairsTab::updateButtons(const QItemSelection &selection) {
  QModelIndexList indexes = selection.indexes();

  if (!indexes.isEmpty()) {
    deleteButton->setEnabled(true);
    duplicateButton->setEnabled(true);
  } else {
    deleteButton->setEnabled(false);
    duplicateButton->setEnabled(false);
  }
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

void PairsTab::setUndoStack(QUndoStack *stack) {
  if (!stack) return;
  auto shortcut = [this](const QKeySequence &key, std::function<void()> f) {
    auto *s = new QShortcut(key, this);
    s->setContext(Qt::WidgetWithChildrenShortcut);
    connect(s, &QShortcut::activated, this, f);
  };
  shortcut(QKeySequence(Qt::CTRL | Qt::Key_Z), [stack]() { stack->undo(); });
  shortcut(QKeySequence(Qt::CTRL | Qt::SHIFT | Qt::Key_Z), [stack]() { stack->redo(); });
  shortcut(QKeySequence(Qt::CTRL | Qt::Key_Y), [stack]() { stack->redo(); });
}

bool PairsTab::setPairNumber(int row, int column, const QString &text, QString &why) {
  if (row < 0 || row >= pairsModel->numPairs() || !PairsModel::editableInPlace(column)) return false;
  bool ok = false;
  double v = QLocale::c().toDouble(text.trimmed(), &ok);
  if (!ok) v = QLocale().toDouble(text.trimmed(), &ok);
  const QString what = PairsModel::plainText(pairsModel->headerData(column, Qt::Horizontal, Qt::DisplayRole).toString()).replace('\n', ' ');
  if (!ok || !std::isfinite(v)) {
    why = tr("'%1' is not a number: %2").arg(text, what);
    return false;
  }
  if (column == 12 && pairsModel->getPairs().at(row).pairType == 0 && v <= 0) {
    why = tr("a particle pair's channel radius must be positive (fm)");
    return false;
  }
  if ((column == 10 || column == PairsModel::LightExcitationColumn) && v < 0) {
    why = tr("an excitation energy is not negative (MeV)");
    return false;
  }
  if (pairsModel->data(pairsModel->index(row, column), Qt::EditRole).toDouble() == v) return true;
  pairsModel->setData(pairsModel->index(row, column), v, Qt::EditRole);
  emit pairEdited(row);
  return true;
}

QStringList PairsTab::pairConsequences(int row) const {
  QStringList out;
  QWidget *w = window();
  if (LevelsTab *levels = w ? w->findChild<LevelsTab *>() : nullptr) {
    QSet<int> named;
    int channels = 0;
    for (const ChannelsData &c : levels->getChannelsModel()->getChannels())
      if (c.pairIndex == row) {
        channels++;
        named.insert(c.levelIndex);
      }
    if (channels > 0) out << tr("%1 channel(s) of %2 level(s) go with it.").arg(channels).arg(named.size());
  }
  if (ProjectHolder *h = w ? w->findChild<ProjectHolder *>() : nullptr) {
    const std::vector<std::string> uses = azure2::front::PairUses(h->project(), row + 1);
    if (!uses.empty()) {
      QStringList names;
      for (const std::string &u : uses) names << QString::fromStdString(u);
      out << tr("Records that name it (as the project was last saved or applied; the engine refuses them once the pair "
                "is gone, and those naming a later pair by key point one pair off): %1").arg(names.join(", "));
    }
  }
  return out;
}
