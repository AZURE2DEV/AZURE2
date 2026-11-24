// The Levels tab's navigation (docs/upgrade/qt_gui.md section 11.3): the filter bar, the add bar, the shortcuts, the
// context menus and what follows the current level.
#include <QAction>
#include <QApplication>
#include <QClipboard>
#include <QHBoxLayout>
#include <QKeyEvent>
#include <QMenu>
#include <QShortcut>
#include <QUndoStack>
#include <cmath>

#include "LevelsTab.h"
#include "PairsTab.h"

QLayout *LevelsTab::buildFilterBar() {
  search_ = new QLineEdit;
  search_->setObjectName("levelSearch");
  search_->setPlaceholderText(tr("3/2- or 7.1"));  // the field's purpose in its tip; the hint whole
  search_->setClearButtonEnabled(true);
  search_->setToolTip(tr("Show the levels whose J^pi starts with this (3/2-, 2+) or whose energy contains it"));
  connect(search_, &QLineEdit::textChanged, this, [this](const QString &t) {
    levelsModelProxy->setSearch(t);
    updateCount();
  });
  jpiFilter_ = new QComboBox;
  jpiFilter_->setObjectName("jpiFilter");
  jpiFilter_->setToolTip(tr("Show the levels of one J^pi; the number of levels of each"));
  jpiFilter_->setSizeAdjustPolicy(QComboBox::AdjustToContents);  // "All Jπ (34)" whole, not "All Jπ (34"
  connect(jpiFilter_, qOverload<int>(&QComboBox::currentIndexChanged), this, &LevelsTab::jpiFilterChosen);
  onlyIncluded_ = new QCheckBox(tr("Only included"));
  onlyIncluded_->setObjectName("onlyIncluded");
  onlyIncluded_->setToolTip(tr("Hide the levels left out of the calculation (Include? off)"));
  connect(onlyIncluded_, &QCheckBox::toggled, this, [this](bool on) {
    levelsModelProxy->setOnlyIncluded(on);
    updateCount();
  });
  count_ = new QLabel;
  count_->setObjectName("levelCount");
  status_ = new QLabel;
  status_->setObjectName("levelsStatus");
  status_->setWordWrap(true);
  status_->setTextFormat(Qt::RichText);
  auto *bar = new QHBoxLayout;
  bar->setContentsMargins(0, 0, 0, 0);
  bar->addWidget(search_, 1);
  bar->addWidget(jpiFilter_, 0);
  bar->addWidget(onlyIncluded_);
  bar->addWidget(count_);
  return bar;
}

QLayout *LevelsTab::buildAddBar() {
  addJ_ = new HalfSpinBox;
  addJ_->setObjectName("addLevelJ");
  addJ_->setValue(0.5);
  addPi_ = new ParityToggle;
  addPi_->setObjectName("addLevelParity");
  addE_ = new QLineEdit;
  addE_->setObjectName("addLevelEnergy");
  addE_->setPlaceholderText(tr("MeV"));  // after the label "E": a hint that fits the field (was cut to "E_x [...")
  addE_->setMinimumWidth(addE_->fontMetrics().horizontalAdvance("00.0000") + 16);
  addE_->setToolTip(tr("The level's excitation energy E_x (MeV); Enter adds the level"));
  addNucleus_ = new QComboBox;
  addNucleus_->setObjectName("addLevelNucleus");
  addNucleus_->setVisible(false);
  addBarButton_ = new QPushButton(tr("Add"));
  addBarButton_->setObjectName("addLevelButton");
  addBarButton_->setToolTip(tr("Add a level of this J^pi at this energy, with the channels its J^pi gives (Enter)"));
  connect(addBarButton_, &QPushButton::clicked, this, &LevelsTab::addLevelFromBarClicked);
  connect(addE_, &QLineEdit::returnPressed, this, &LevelsTab::addLevelFromBarClicked);
  duplicateButton_ = new QPushButton(tr("Duplicate"));
  duplicateButton_->setObjectName("duplicateLevelButton");
  duplicateButton_->setToolTip(tr("A copy of the selected levels 10 keV above, with their widths (Ctrl+D)"));
  duplicateButton_->setEnabled(false);
  connect(duplicateButton_, &QPushButton::clicked, this, &LevelsTab::duplicateSelectedLevels);
  auto *bar = new QHBoxLayout;
  bar->setContentsMargins(0, 0, 0, 0);
  bar->addWidget(new QLabel(tr("New level: J")));
  bar->addWidget(addJ_);
  bar->addWidget(addPi_);
  bar->addWidget(new QLabel(tr("E")));
  bar->addWidget(addE_, 1);
  bar->addWidget(addNucleus_);
  bar->addWidget(addBarButton_);
  bar->addWidget(duplicateButton_);
  return bar;
}

void LevelsTab::setUpShortcuts() {
  undo_ = new QUndoStack(this);
  auto shortcut = [this](const QKeySequence &key, QWidget *on, Qt::ShortcutContext context, auto slot) {
    auto *s = new QShortcut(key, on);
    s->setContext(context);
    connect(s, &QShortcut::activated, this, slot);
    return s;
  };
  shortcut(QKeySequence::Delete, levelsView, Qt::WidgetShortcut, &LevelsTab::removeSelectedLevels);
  shortcut(QKeySequence(Qt::Key_Insert), this, Qt::WidgetWithChildrenShortcut, &LevelsTab::focusAddBar);
  // Ctrl+D in the levels table (the window's Help menu has Ctrl+D too: the table takes it first, ShortcutOverride)
  levelsView->installEventFilter(this);
  shortcut(QKeySequence(Qt::CTRL | Qt::Key_Z), this, Qt::WidgetWithChildrenShortcut, [this]() { undo_->undo(); });
  shortcut(QKeySequence(Qt::CTRL | Qt::SHIFT | Qt::Key_Z), this, Qt::WidgetWithChildrenShortcut, [this]() { undo_->redo(); });
  shortcut(QKeySequence(Qt::CTRL | Qt::Key_Y), this, Qt::WidgetWithChildrenShortcut, [this]() { undo_->redo(); });
  // between the levels and the channels of the current level
  shortcut(QKeySequence(Qt::CTRL | Qt::Key_Right), levelsView, Qt::WidgetShortcut, [this]() {
    channelsView->setFocus();
    if (!channelsView->currentIndex().isValid() && proxyModel->rowCount() > 0)
      channelsView->setCurrentIndex(proxyModel->index(0, ChannelsModel::ValueColumn));
  });
  shortcut(QKeySequence(Qt::CTRL | Qt::Key_Left), channelsView, Qt::WidgetShortcut, [this]() { levelsView->setFocus(); });
  setTabOrder(levelsView, channelsView);
  // every change of the three models and the limits is one undo entry once the event loop turns (LevelsUndo.h)
  for (QAbstractItemModel *m : {static_cast<QAbstractItemModel *>(levelsModel), static_cast<QAbstractItemModel *>(channelsModel)}) {
    connect(m, &QAbstractItemModel::dataChanged, this, &LevelsTab::undoLater);
    connect(m, &QAbstractItemModel::rowsInserted, this, &LevelsTab::undoLater);
    connect(m, &QAbstractItemModel::rowsRemoved, this, &LevelsTab::undoLater);
  }
  connect(channelsModel, &ChannelsModel::widthEdited, this, &LevelsTab::undoLater);
  for (QSpinBox *b : {maxLSpin, maxMultSpin, maxNumMultSpin})
    connect(b, qOverload<int>(&QSpinBox::valueChanged), this, &LevelsTab::undoLater);
  connect(includeAllJPiCheck, &QCheckBox::toggled, this, &LevelsTab::undoLater);
  // the J^pi list and the count follow the levels
  for (auto signal : {&QAbstractItemModel::rowsInserted, &QAbstractItemModel::rowsRemoved})
    connect(levelsModel, signal, this, [this]() {
      if (!loading_) updateJpiFilter();
    });
  connect(levelsModel, &QAbstractItemModel::dataChanged, this, [this]() {
    if (!loading_) updateJpiFilter();
  });
  connect(channelsModel, &QAbstractItemModel::rowsInserted, this, [this]() {
    if (!loading_ && !restoring_) updateFilter(QItemSelection());
  });
}

void LevelsTab::addLevelFromBarClicked() {
  QString why;
  if (!addLevelFromBar(why)) {
    setStatus(tr("Not added: %1").arg(why), true);
    return;
  }
  // ready for the next energy of the same J^pi
  addE_->clear();
  addE_->setFocus();
}

void LevelsTab::focusAddBar() {
  const int row = currentLevelRow();
  if (row >= 0) {
    const LevelsData &l = levelsModel->getLevels().at(row);
    addJ_->setValue(l.jValue);
    addPi_->setParity(l.piValue);
  }
  addE_->setFocus();
  addE_->selectAll();
}

void LevelsTab::fillAddBar(const QString &link) {
  // "jpi:<2J>:<pi>:<Z>:<A>"
  const QStringList p = link.split(':');
  if (p.size() < 3 || p.first() != "jpi") return;
  addJ_->setValue(p.at(1).toInt() / 2.0);
  addPi_->setParity(p.at(2).toInt());
  if (p.size() >= 5 && pairsModel) {
    const int n = pairsModel->nucleusIndex(Nucleus{p.at(3).toInt(), p.at(4).toInt()});
    if (n >= 0 && n < addNucleus_->count()) addNucleus_->setCurrentIndex(n);
  }
  addE_->setFocus();
  setStatus(tr("Type the energy of the new %1 level and press Enter.")
                .arg(AddLevelDialog::spinParityLabel(addJ_->value(), addPi_->parity())));
}

void LevelsTab::currentLevelChanged() {
  updateFilter(QItemSelection());
  const int row = currentLevelRow();
  if (row >= 0 && row < levelsModel->getLevels().size()) {
    const LevelsData &l = levelsModel->getLevels().at(row);
    levelsModelProxy->setHighlight((int)std::lround(2 * l.jValue), l.piValue);
  } else {
    levelsModelProxy->setHighlight(-1, 0);
  }
  const bool any = !levelsView->selectionModel()->selectedRows().isEmpty();
  removeLevelButton->setEnabled(any);
  duplicateButton_->setEnabled(any);
}

void LevelsTab::updateJpiFilter() {
  // every J^pi of the levels with its count, the choice kept
  const QString kept = jpiFilter_->currentData().toString();
  QMap<QPair<int, int>, int> counts;
  for (const LevelsData &l : levelsModel->getLevels()) counts[{(int)std::lround(2 * l.jValue), -l.piValue}]++;
  jpiFilter_->blockSignals(true);
  jpiFilter_->clear();
  jpiFilter_->addItem(tr("All Jπ (%1)").arg(levelsModel->getLevels().size()), QString());
  for (auto it = counts.begin(); it != counts.end(); ++it) {
    const int twoJ = it.key().first, pi = -it.key().second;
    jpiFilter_->addItem(QString("%1 (%2)").arg(AddLevelDialog::spinParityLabel(twoJ / 2.0, pi)).arg(it.value()),
                        QString("%1:%2").arg(twoJ).arg(pi));
  }
  const int k = jpiFilter_->findData(kept);
  jpiFilter_->setCurrentIndex(k >= 0 ? k : 0);
  jpiFilter_->blockSignals(false);
  if (k < 0 && !kept.isEmpty()) jpiFilterChosen(0);
  updateCount();
}

void LevelsTab::jpiFilterChosen(int index) {
  const QString d = jpiFilter_->itemData(index).toString();
  if (d.isEmpty()) levelsModelProxy->setJpi(-1, 0);
  else levelsModelProxy->setJpi(d.section(':', 0, 0).toInt(), d.section(':', 1, 1).toInt());
  updateCount();
}

void LevelsTab::updateCount() {
  if (!count_) return;
  const int shown = levelsModelProxy->rowCount(), all = levelsModel->rowCount(QModelIndex());
  count_->setText(shown == all ? tr("%n level(s)", "", all) : tr("%1 of %2 levels").arg(shown).arg(all));
}

void LevelsTab::showLevelsMenu(const QPoint &pos) {
  QMenu menu(this);
  const QList<int> rows = selectedLevelRows();
  const bool any = !rows.isEmpty();
  const int row = currentLevelRow();
  QString jpi;
  if (row >= 0) {
    const LevelsData &l = levelsModel->getLevels().at(row);
    jpi = AddLevelDialog::spinParityLabel(l.jValue, l.piValue);
  }
  menu.addAction(tr("Add a level of %1...").arg(jpi.isEmpty() ? tr("a J^pi") : jpi), this, &LevelsTab::focusAddBar)->setShortcut(Qt::Key_Insert);
  menu.addAction(tr("Add a level in the dialog..."), this, qOverload<>(&LevelsTab::addLevel));
  QAction *a = menu.addAction(tr("Duplicate"), this, &LevelsTab::duplicateSelectedLevels);
  a->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_D));
  a->setEnabled(any);
  a = menu.addAction(tr("Edit in the dialog..."), this, &LevelsTab::editLevel);
  a->setEnabled(any);
  a = menu.addAction(tr("Remove %n level(s)", "", rows.size()), this, &LevelsTab::removeSelectedLevels);
  a->setShortcut(QKeySequence::Delete);
  a->setEnabled(any);
  menu.addSeparator();
  auto bulk = [this, &menu, any, rows](const QString &text, std::function<void()> f) {
    QAction *b = menu.addAction(text, this, f);
    b->setEnabled(any);
  };
  bulk(tr("Include"), [this, rows]() {
    undoText_ = tr("include levels");
    for (int r : rows) levelsModel->setData(levelsModel->index(r, 0), 1, Qt::EditRole);
    updateMissingSpinParities();
  });
  bulk(tr("Leave out of the calculation"), [this, rows]() {
    undoText_ = tr("leave out levels");
    for (int r : rows) levelsModel->setData(levelsModel->index(r, 0), 0, Qt::EditRole);
    updateMissingSpinParities();
  });
  bulk(tr("Fix the energies"), [this, rows]() { setFixed(rows, true, true); });
  bulk(tr("Free the energies"), [this, rows]() { setFixed(rows, true, false); });
  bulk(tr("Fix the widths"), [this, rows]() { setFixed(rows, false, true); });
  bulk(tr("Free the widths"), [this, rows]() { setFixed(rows, false, false); });
  menu.addSeparator();
  a = menu.addAction(tr("Show only %1").arg(jpi), this, [this, row]() {
    const LevelsData &l = levelsModel->getLevels().at(row);
    const int k = jpiFilter_->findData(QString("%1:%2").arg((int)std::lround(2 * l.jValue)).arg(l.piValue));
    if (k >= 0) jpiFilter_->setCurrentIndex(k);
  });
  a->setEnabled(row >= 0);
  a = menu.addAction(tr("Add a background pole of %1").arg(jpi), this, &LevelsTab::addBackgroundPole);
  a->setEnabled(row >= 0);
  a->setToolTip(tr("A level 10 MeV above the J^pi's highest, energy fixed, widths 0: edit it"));
  menu.addSeparator();
  menu.addAction(undo_->createUndoAction(&menu, tr("Undo")));
  menu.addAction(undo_->createRedoAction(&menu, tr("Redo")));
  menu.exec(levelsView->viewport()->mapToGlobal(pos));
}

void LevelsTab::showChannelsMenu(const QPoint &pos) {
  const QModelIndex at = channelsView->indexAt(pos);
  if (!at.isValid()) return;
  const int row = proxyModel->mapToSource(at).row();
  QMenu menu(this);
  const ChannelsData c = channelsModel->getChannels().at(row);
  menu.addAction(c.isFixed ? tr("Free") : tr("Fix"), this, [this, row, c]() {
    undoText_ = c.isFixed ? tr("free a width") : tr("fix a width");
    channelsModel->setData(channelsModel->index(row, 0), c.isFixed ? 0 : 1, Qt::EditRole);
  });
  menu.addAction(tr("Edit the value"), this, [this, at]() { channelsView->edit(at.sibling(at.row(), ChannelsModel::ValueColumn)); });
  menu.addAction(tr("Set to 0"), this, [this, row]() {
    QString why;
    setChannelValue(row, "0", why);
  });
  menu.addAction(tr("Copy the value"), this, [c]() { QApplication::clipboard()->setText(QString::number(c.reducedWidth, 'g', 17)); });
  menu.exec(channelsView->viewport()->mapToGlobal(pos));
}

bool LevelsTab::eventFilter(QObject *object, QEvent *event) {
  if (object == levelsView && (event->type() == QEvent::ShortcutOverride || event->type() == QEvent::KeyPress)) {
    auto *key = static_cast<QKeyEvent *>(event);
    if (key->key() == Qt::Key_D && key->modifiers() == Qt::ControlModifier) {
      if (event->type() == QEvent::ShortcutOverride) event->accept();
      else duplicateSelectedLevels();
      return true;
    }
  }
  return QWidget::eventFilter(object, event);
}
