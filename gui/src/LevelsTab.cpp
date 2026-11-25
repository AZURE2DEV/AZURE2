#include <QHeaderView>
#include <QMessageBox>
#include <QInputDialog>
#include <QFile>
#include <QRegularExpression>
#include <QTimer>
#include <algorithm>

#include "LevelsTab.h"
#include "ButtonSizes.h"
#include "ExpertMode.h"
#include "RecordExtras.h"
#include "RecordsPanel.h"
#include "LevelsHeaderView.h"
#include "AddLevelDialog.h"
#include "IAEALevelDialog.h"
#include "RichTextDelegate.h"
#include "InfoDialog.h"
#include "ElementMap.h"
#include "azure2/front/Engine.h"
#include "PairsTab.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <cmath>

LevelsTab::LevelsTab(QWidget *parent) :
  QWidget(parent) {
  pairsModel = nullptr;
  levelsModel = new LevelsModel(this);
  levelsModelProxy = new LevelsFilterProxy(this);
  levelsModelProxy->setSourceModel(levelsModel);
  levelsView = new QTableView;
  levelsView->setHorizontalHeader(new LevelsHeaderView(Qt::Horizontal, levelsView));
  levelsView->setModel(levelsModelProxy);
  levelsView->horizontalHeader()->setSortIndicator(4, Qt::AscendingOrder);
  levelsView->setSortingEnabled(true);
  levelsView->verticalHeader()->hide();
  levelsView->horizontalHeader()->setHighlightSections(false);
  levelsView->setColumnWidth(0, 60);
  levelsView->setColumnWidth(1, 40);
  levelsView->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Fixed);
  levelsView->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Fixed);
  levelsView->horizontalHeader()->setSectionResizeMode(2, QHeaderView::Stretch);
  levelsView->horizontalHeader()->setSectionResizeMode(4, QHeaderView::Stretch);
  levelsView->setColumnHidden(3, true);
  // The compound nucleus of each level: only shown when the pairs form more than one (updateNucleusColumn).
  levelsView->setColumnHidden(LevelsModel::NucleusColumn, true);
  levelsView->setItemDelegateForColumn(LevelsModel::NucleusColumn, new RichTextDelegate());
  levelsView->setSelectionBehavior(QAbstractItemView::SelectRows);
  // several levels at once (Del, bulk fix); J^pi and E edited in place (docs/upgrade/qt_gui.md section 11.2)
  levelsView->setSelectionMode(QAbstractItemView::ExtendedSelection);
  levelsView->setEditTriggers(QAbstractItemView::DoubleClicked | QAbstractItemView::EditKeyPressed |
                              QAbstractItemView::AnyKeyPressed);
  levelsView->setShowGrid(false);
  levelsView->setObjectName("levelsView");
  levelsView->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(levelsView, &QWidget::customContextMenuRequested, this, &LevelsTab::showLevelsMenu);
  connect(levelsView->selectionModel(), SIGNAL(selectionChanged(QItemSelection, QItemSelection)), this, SLOT(updateButtons(QItemSelection)));
  connect(levelsView->selectionModel(), &QItemSelectionModel::currentRowChanged, this, &LevelsTab::currentLevelChanged);
  connect(levelsView->selectionModel(), &QItemSelectionModel::selectionChanged, this, &LevelsTab::currentLevelChanged);
  {
    auto *jpiDelegate = new CheckedEditDelegate(
        CheckedEditDelegate::Kind::Jpi,
        [this](const QModelIndex &index, const QVariant &value, QString &why) {
          const QVariantList v = value.toList();
          return setLevelJpi(levelsModelProxy->mapToSource(index).row(), v.value(0).toDouble(), v.value(1).toInt(), why);
        },
        this);
    auto *energyDelegate = new CheckedEditDelegate(
        CheckedEditDelegate::Kind::Number,
        [this](const QModelIndex &index, const QVariant &value, QString &why) {
          return setLevelEnergy(levelsModelProxy->mapToSource(index).row(), value.toString(), why);
        },
        this);
    for (CheckedEditDelegate *d : {jpiDelegate, energyDelegate}) {
      connect(d, &CheckedEditDelegate::refused, this, &LevelsTab::showRefusal);
      connect(d, &CheckedEditDelegate::moveDown, this, [this](const QModelIndex &) {
        const QModelIndex c = levelsView->currentIndex();
        const QModelIndex next = c.sibling(c.row() + 1, c.column());
        if (next.isValid()) levelsView->setCurrentIndex(next);
      });
    }
    levelsView->setItemDelegateForColumn(2, jpiDelegate);
    levelsView->setItemDelegateForColumn(4, energyDelegate);
  }

  maxLSpin = new QSpinBox;
  maxLSpin->setMinimum(0);
  maxLSpin->setMaximum(10);
  maxLSpin->setSingleStep(1);
  maxLSpin->setValue(2);
  QLabel *maxLLabel = new QLabel(tr("Maximum Orbital Momentum"));
  maxMultSpin = new QSpinBox;
  maxMultSpin->setObjectName("maxMultSpin");
  maxMultSpin->setMinimum(1);
  maxMultSpin->setMaximum(10);
  maxMultSpin->setSingleStep(1);
  maxMultSpin->setValue(2);
  QLabel *maxMultLabel = new QLabel(tr("Maximum Gamma Multipolarity"));
  maxMultLabelW_ = maxMultLabel;
  maxNumMultSpin = new QSpinBox;
  maxNumMultSpin->setMinimum(1);
  maxNumMultSpin->setMaximum(10);
  maxNumMultSpin->setSingleStep(1);
  maxNumMultSpin->setValue(2);
  QLabel *maxNumMultLabel = new QLabel(tr("Maximum Gamma Multipolarities\nPer Decay"));
  maxNumMultLabelW_ = maxNumMultLabel;
  // two lines: the configuration column narrow, the channel table wide (qt_gui.md section 19.8)
  includeAllJPiCheck = new QCheckBox(tr("Hard-sphere scattering in every J%1\n(levels or not)").arg(QChar(0x03C0)));
  includeAllJPiCheck->setToolTip(tr("AZURE2 only includes a J\u03c0 in the calculation when a level carries it. Without this, every J\u03c0 "
                                    "that has no level is missing from elastic scattering (its hard-sphere phase shift), unless a fixed "
                                    "zero-width level is added by hand. With it, every J\u03c0 the particle pairs allow up to the maximum "
                                    "orbital momentum scatters, levels or not; switching a level off then removes its resonance and keeps "
                                    "the partial wave. Switch it off to recover the old convention, where a J\u03c0 disappears with its levels."));
  includeAllJPiCheck->setChecked(true);
  missingJPiLabel = new QLabel;
  missingJPiLabel->setWordWrap(true);
  missingJPiLabel->setTextFormat(Qt::RichText);
  missingJPiLabel->setObjectName("missingJPiLabel");
  connect(missingJPiLabel, &QLabel::linkActivated, this, &LevelsTab::fillAddBar);
  connect(includeAllJPiCheck, SIGNAL(toggled(bool)), this, SLOT(updateMissingSpinParities()));
  connect(maxLSpin, SIGNAL(valueChanged(int)), this, SLOT(updateChannelsPairAddedEdited()));
  connect(maxMultSpin, SIGNAL(valueChanged(int)), this, SLOT(updateChannelsPairAddedEdited()));
  connect(maxNumMultSpin, SIGNAL(valueChanged(int)), this, SLOT(updateChannelsPairAddedEdited()));

  channelsModel = new ChannelsModel(this);
  proxyModel = new QSortFilterProxyModel(this);
  proxyModel->setSourceModel(channelsModel);
  proxyModel->setDynamicSortFilter(true);
  proxyModel->setFilterKeyColumn(1);
  proxyModel->setFilterRegularExpression("-1");
  proxyModel->sort(1, Qt::AscendingOrder);
  channelsView = new QTableView;
  channelsView->setModel(proxyModel);
  channelsView->verticalHeader()->hide();
  channelsView->horizontalHeader()->setHighlightSections(false);
  channelsView->setColumnWidth(0, 40);
  channelsView->setColumnWidth(2, 150);
  channelsView->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Fixed);
  channelsView->horizontalHeader()->setSectionResizeMode(2, QHeaderView::Fixed);
  // s and l: wide enough for "3/2" and "E2" in one line
  channelsView->horizontalHeader()->setSectionResizeMode(3, QHeaderView::Fixed);
  channelsView->horizontalHeader()->setSectionResizeMode(4, QHeaderView::Fixed);
  channelsView->setColumnWidth(3, 38);
  channelsView->setColumnWidth(4, 34);
  channelsView->setColumnHidden(1, true);
  channelsView->setColumnHidden(5, true);
  // the value in place, what it is and its unit, theta^2 (docs/upgrade/qt_gui.md section 11.4)
  channelsView->setObjectName("channelsView");
  channelsView->horizontalHeader()->setSectionResizeMode(ChannelsModel::ValueColumn, QHeaderView::Stretch);
  channelsView->horizontalHeader()->setSectionResizeMode(ChannelsModel::MeaningColumn, QHeaderView::Fixed);
  channelsView->horizontalHeader()->setSectionResizeMode(ChannelsModel::RatioColumn, QHeaderView::Fixed);
  channelsView->setColumnWidth(ChannelsModel::MeaningColumn, 78);
  channelsView->setColumnWidth(ChannelsModel::RatioColumn, channelsView->fontMetrics().horizontalAdvance("0.000000") + 12);
  channelsView->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(channelsView, &QWidget::customContextMenuRequested, this, &LevelsTab::showChannelsMenu);
  {
    auto *valueDelegate = new CheckedEditDelegate(
        CheckedEditDelegate::Kind::Number,
        [this](const QModelIndex &index, const QVariant &value, QString &why) {
          return setChannelValue(proxyModel->mapToSource(index).row(), value.toString(), why);
        },
        this);
    connect(valueDelegate, &CheckedEditDelegate::refused, this, &LevelsTab::showRefusal);
    connect(valueDelegate, &CheckedEditDelegate::moveDown, this, [this](const QModelIndex &) {
      const QModelIndex c = channelsView->currentIndex();
      const QModelIndex next = c.sibling(c.row() + 1, c.column());
      if (next.isValid()) channelsView->setCurrentIndex(next);
    });
    channelsView->setItemDelegateForColumn(ChannelsModel::ValueColumn, valueDelegate);
  }
  channelsView->setItemDelegateForColumn(2, new RichTextDelegate());
  channelsView->setItemDelegateForColumn(3, new RichTextDelegate());
  channelsView->setItemDelegateForColumn(4, new RichTextDelegate());
  channelsView->setSelectionBehavior(QAbstractItemView::SelectRows);
  channelsView->setSelectionMode(QAbstractItemView::SingleSelection);
  channelsView->setEditTriggers(QAbstractItemView::DoubleClicked | QAbstractItemView::EditKeyPressed |
                                QAbstractItemView::AnyKeyPressed);
  channelsView->setShowGrid(false);
  connect(channelsView->selectionModel(), SIGNAL(selectionChanged(QItemSelection, QItemSelection)), this, SLOT(updateDetails(QItemSelection)));

  channelDetails = new ChannelDetails(this);
  channelDetails->hide();
  connect(channelDetails->reducedWidthText, SIGNAL(textEdited(const QString &)), this, SLOT(updateReducedWidth(const QString &)));
  connect(channelDetails->wignerButton, SIGNAL(clicked()), this, SLOT(calculateWignerLimit()));

  addLevelButton = new QPushButton(tr("+"));
  connect(addLevelButton, SIGNAL(clicked()), this, SLOT(addLevel()));
  removeLevelButton = new QPushButton(tr("-"));
  removeLevelButton->setEnabled(false);
  connect(removeLevelButton, SIGNAL(clicked()), this, SLOT(removeLevel()));
  iaeaQueryButton = new QPushButton(tr("NDS"));
  iaeaQueryButton->setToolTip(tr("Query NDS database for level schemes"));
  connect(iaeaQueryButton, SIGNAL(clicked()), this, SLOT(addLevelsFromIAEA()));

  fixAllWidthsButton = new QPushButton(tr("Fix All Widths"));
  connect(fixAllWidthsButton, SIGNAL(clicked()), this, SLOT(fixAllWidths()));

  fixAllEnergiesButton = new QPushButton(tr("Fix All Energies"));
  connect(fixAllEnergiesButton, SIGNAL(clicked()), this, SLOT(fixAllEnergies()));

  nucleusFilterCombo = new QComboBox;
  nucleusFilterCombo->setObjectName("nucleusFilterCombo");
  nucleusFilterCombo->setToolTip(tr("Show the levels of one compound nucleus only."));
  nucleusFilterCombo->setVisible(false);
  connect(nucleusFilterCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(filterByNucleus(int)));
  exportLatexButton = new QPushButton(tr("Export LaTeX"));
  exportLatexButton->setToolTip(tr("Export levels and widths to LaTeX table"));
  connect(exportLatexButton, SIGNAL(clicked()), this, SLOT(exportLatexTable()));
  // the row's buttons sized from the font, not in pixels: "NDS" was cut under a larger font (qt_gui.md 19.11)
  {
    const int row = buttons::RowHeight(fixAllWidthsButton);
    for (QPushButton *b : {addLevelButton, removeLevelButton}) buttons::FitGlyph(b, row);
    for (QPushButton *b : {iaeaQueryButton, fixAllWidthsButton, fixAllEnergiesButton, exportLatexButton}) buttons::FitText(b, row);
  }

  /*
  mapper = new QSignalMapper(this);
  connect(mapper,SIGNAL(mapped(int)),this,SLOT(showInfo(int)));

  infoButton[0] = new QPushButton(this);
  infoButton[0]->setMaximumSize(28,28);
  infoButton[0]->setIcon(style()->standardIcon(QStyle::SP_MessageBoxInformation));
  mapper->setMapping(infoButton[0],1);
  connect(infoButton[0],SIGNAL(clicked()),mapper,SLOT(map()));
  */

  QGridLayout *buttonBox = new QGridLayout;
  buttonBox->addWidget(addLevelButton, 0, 0);
  buttonBox->addWidget(removeLevelButton, 0, 1);
  buttonBox->addWidget(iaeaQueryButton, 0, 2);
  buttonBox->addItem(new QSpacerItem(28, 28), 0, 3);
  // buttonBox->addWidget(infoButton[0],0,4);
  buttonBox->setColumnStretch(0, 0);
  buttonBox->setColumnStretch(1, 0);
  buttonBox->setColumnStretch(2, 0);
  buttonBox->setColumnStretch(3, 1);
  buttonBox->setColumnStretch(4, 0);
#ifdef MACX_SPACING
  buttonBox->setHorizontalSpacing(11);
#else
  buttonBox->setHorizontalSpacing(0);
#endif

  QGridLayout *fixButtonBox = new QGridLayout;
  fixButtonBox->addWidget(fixAllWidthsButton, 0, 0);
  fixButtonBox->addItem(new QSpacerItem(10, 28), 0, 1);
  fixButtonBox->addWidget(fixAllEnergiesButton, 0, 2);
  fixButtonBox->addItem(new QSpacerItem(10, 28), 0, 3);
  fixButtonBox->addWidget(exportLatexButton, 0, 4);
  fixButtonBox->addItem(new QSpacerItem(28, 28), 0, 5);
  fixButtonBox->setColumnStretch(0, 0);
  fixButtonBox->setColumnStretch(1, 0);
  fixButtonBox->setColumnStretch(2, 0);
  fixButtonBox->setColumnStretch(3, 0);
  fixButtonBox->setColumnStretch(4, 0);
  fixButtonBox->setColumnStretch(5, 1);
#ifdef MACX_SPACING
  fixButtonBox->setHorizontalSpacing(11);
#else
  fixButtonBox->setHorizontalSpacing(5);
#endif

  QGroupBox *levelsBox = new QGroupBox(tr("Compound Nucleus Levels"));
  QGridLayout *levelsLayout = new QGridLayout;
  levelsLayout->setContentsMargins(6, 6, 6, 6);
  levelsLayout->addWidget(nucleusFilterCombo, 0, 0);
  levelsLayout->addLayout(buildFilterBar(), 1, 0);
  levelsLayout->addWidget(levelsView, 2, 0);
  levelsLayout->addLayout(buildAddBar(), 3, 0);
  levelsLayout->addLayout(buttonBox, 4, 0);
  levelsLayout->addLayout(fixButtonBox, 5, 0);
  levelsLayout->addWidget(status_, 6, 0);
#ifdef MACX_SPACING
  levelsLayout->setVerticalSpacing(0);
#endif
  levelsBox->setLayout(levelsLayout);

  QGroupBox *configBox = new QGroupBox(tr("Channel Configuration"));
  QGridLayout *configLayout = new QGridLayout;
  configLayout->addWidget(maxLSpin, 0, 0);
  configLayout->addWidget(maxLLabel, 0, 1);
  configLayout->addWidget(maxMultSpin, 1, 0);
  configLayout->addWidget(maxMultLabel, 1, 1);
  configLayout->addWidget(maxNumMultSpin, 2, 0);
  configLayout->addWidget(maxNumMultLabel, 2, 1);
  configLayout->addWidget(includeAllJPiCheck, 3, 0, 1, 2);
  configLayout->addWidget(missingJPiLabel, 4, 0, 1, 2);
  configLayout->setRowStretch(5, 1);  // the rows at the top, the spare height below them (rows hidden in standard view)
  configLayout->setColumnStretch(0, 0);
  configLayout->setColumnStretch(1, 1);
  configBox->setLayout(configLayout);

  QGroupBox *channelsBox = new QGroupBox(tr("Channels In Selected Level"));
  QGridLayout *channelsLayout = new QGridLayout;
  channelsLayout->setContentsMargins(6, 6, 6, 6);
  channelsLayout->addWidget(channelsView, 0, 0);
#ifdef MACX_SPACING
  channelsLayout->addItem(new QSpacerItem(40, 40), 1, 0);
  channelsLayout->setVerticalSpacing(0);
#else
  channelsLayout->addItem(new QSpacerItem(34, 34), 1, 0);
#endif
  channelsBox->setLayout(channelsLayout);

  QGridLayout *detailsBox = new QGridLayout;
  QLabel *detailsLabel = new QLabel(tr("Channel Details (select from list to view):"));
  detailsLabel->setAlignment(Qt::AlignHCenter);
  detailsBox->addWidget(detailsLabel, 0, 0);
  detailsBox->addWidget(channelDetails, 1, 0);
  detailsBox->setRowStretch(0, 0);
  detailsBox->setRowStretch(1, 1);

  QGridLayout *mainLayout = new QGridLayout;
  mainLayout->addWidget(levelsBox, 0, 0, 2, 1);
  mainLayout->addWidget(channelsBox, 0, 1, 2, 1);
  mainLayout->addWidget(configBox, 0, 2, 1, 1);
  mainLayout->addLayout(detailsBox, 1, 2, 1, 1);
  // the channel table holds nine columns, the configuration a few lines: the middle column wider, so that the Value
  // column (the width or ANC) is not cut to "2..." (qt_gui.md section 19.8)
  mainLayout->setColumnStretch(0, 4);
  mainLayout->setColumnStretch(1, 6);
  mainLayout->setColumnStretch(2, 2);  // the configuration at its natural width, the rest to the channels
  mainLayout->setRowStretch(0, 3);
  mainLayout->setRowStretch(1, 4);
  // the records that belong to the level scheme (gamma decays, isospin): sub-tabs below the levels in a splitter, a
  // panel without records hidden in standard view (RecordExtras; qt_gui.md sections 19.3, 19.4)
  QWidget *content = new QWidget;
  content->setLayout(mainLayout);
  extras = new RecordExtras(content, "levelsExtras");
  QVBoxLayout *outer = new QVBoxLayout;
  outer->setContentsMargins(0, 0, 0, 0);
  outer->addWidget(extras);
  setLayout(outer);
  expert_ = expert::Enabled();
  setUpShortcuts();
}

void LevelsTab::setPairsModel(PairsModel *model) {
  pairsModel = model;
  channelsModel->setPairsModel(model);
  if (model) {
    connect(model, SIGNAL(rowsInserted(QModelIndex, int, int)), this, SLOT(updateNucleusColumn()));
    connect(model, SIGNAL(rowsRemoved(QModelIndex, int, int)), this, SLOT(updateNucleusColumn()));
    connect(model, SIGNAL(dataChanged(QModelIndex, QModelIndex)), this, SLOT(updateNucleusColumn()));
    connect(model, SIGNAL(modelReset()), this, SLOT(updateNucleusColumn()));
    // a pair's edit is an undo entry of this tab's stack (the Particle Pairs tab shares it)
    connect(model, &QAbstractItemModel::dataChanged, this, &LevelsTab::undoLater);
    connect(model, &QAbstractItemModel::rowsInserted, this, &LevelsTab::undoLater);
    connect(model, &QAbstractItemModel::rowsRemoved, this, &LevelsTab::undoLater);
    if (PairsTab *pairs = qobject_cast<PairsTab *>(model->parent())) pairs->setUndoStack(undo_);
  }
  updateNucleusColumn();
}

void LevelsTab::updateNucleusColumn() {
  QList<Nucleus> nuclei = pairsModel ? pairsModel->nuclei() : QList<Nucleus>();
  bool several = nuclei.size() > 1;
  levelsView->setColumnHidden(LevelsModel::NucleusColumn, !several);
  // The filter drop-down follows the nuclei, keeping the selection when its nucleus is still there.
  QString selected = nucleusFilterCombo->currentIndex() > 0 ? nucleusFilterCombo->currentText() : QString();
  nucleusFilterCombo->blockSignals(true);
  nucleusFilterCombo->clear();
  nucleusFilterCombo->addItem(tr("All compound nuclei"));
  for (const Nucleus &nucleus : nuclei) {
    nucleusFilterCombo->addItem(pairsModel->nucleusLongLabel(nucleus));
    nucleusFilterCombo->setItemData(nucleusFilterCombo->count() - 1, PairsModel::nucleusLabel(nucleus), Qt::UserRole);
  }
  int keep = selected.isEmpty() ? -1 : nucleusFilterCombo->findText(selected);
  nucleusFilterCombo->setCurrentIndex(keep > 0 ? keep : 0);
  nucleusFilterCombo->blockSignals(false);
  nucleusFilterCombo->setVisible(several);
  filterByNucleus(nucleusFilterCombo->currentIndex());
  // the add bar's nucleus: shown with several
  if (addNucleus_) {
    const int keepAdd = addNucleus_->currentIndex();
    addNucleus_->clear();
    for (const Nucleus &nucleus : nuclei) addNucleus_->addItem(PairsModel::nucleusLabel(nucleus, false));
    addNucleus_->setCurrentIndex(keepAdd >= 0 && keepAdd < addNucleus_->count() ? keepAdd : 0);
    addNucleus_->setVisible(several);
  }
  if (!several || restoring_) return;
  // A level added before a second nucleus appeared has no nucleus: it would take channels from every pair, which
  // the engine refuses (a level of two compound nuclei), so it is given the first nucleus.
  bool assigned = false;
  QList<LevelsData> levels = levelsModel->getLevels();
  for (int i = 0; i < levels.size(); i++) {
    if (levels.at(i).hasNucleus()) continue;
    levelsModel->setNucleus(i, nuclei.first());
    assigned = true;
  }
  if (assigned) updateChannelsPairAddedEdited();
}

void LevelsTab::filterByNucleus(int index) {
  if (index <= 0 || index >= nucleusFilterCombo->count()) levelsModelProxy->setNucleus(QString());
  else levelsModelProxy->setNucleus(nucleusFilterCombo->itemData(index, Qt::UserRole).toString());
  updateCount();
}

void LevelsTab::addLevel() {
  AddLevelDialog aDialog;
  QList<Nucleus> nuclei = pairsModel ? pairsModel->nuclei() : QList<Nucleus>();
  QStringList labels;
  QList<QList<SpinParityOption> > perNucleus;
  for (const Nucleus &nucleus : nuclei) {
    labels << pairsModel->nucleusLongLabel(nucleus);
    perNucleus << spinParityOptions(nucleus);
  }
  aDialog.setSpinParityOptions(spinParityOptions());
  aDialog.setNucleusOptions(labels, perNucleus, 0, true);
  if (aDialog.exec()) {
    LevelsData newLevel;
    newLevel.isActive = 1;
    newLevel.isFixed = 0;
    newLevel.jValue = aDialog.jValue();
    newLevel.piValue = aDialog.piValue();
    newLevel.energy = (aDialog.energyText->text()).toDouble();
    if (!nuclei.isEmpty()) {
      Nucleus nucleus = nuclei.at(aDialog.nucleusIndex());
      newLevel.nucZ = nucleus.Z;
      newLevel.nucA = nucleus.A;
    }
    addLevel(newLevel, false);
  }
}

void LevelsTab::addLevel(LevelsData level, bool fromFile) {
  QList<LevelsData> levels = levelsModel->getLevels();
  // A level of a single-nucleus project belongs to that nucleus; with several, callers say which.
  if (!level.hasNucleus() && pairsModel && pairsModel->nuclei().size() == 1) {
    level.nucZ = pairsModel->nuclei().first().Z;
    level.nucA = pairsModel->nuclei().first().A;
  }
  if (levelsModel->isLevel(level) == -1) {
    levelsModel->insertRows(levels.size(), 1, QModelIndex());
    QModelIndex index = levelsModel->index(levels.size(), 0, QModelIndex());
    levelsModel->setData(index, level.isActive, Qt::EditRole);
    index = levelsModel->index(levels.size(), 1, QModelIndex());
    levelsModel->setData(index, level.isFixed, Qt::EditRole);
    index = levelsModel->index(levels.size(), 2, QModelIndex());
    levelsModel->setData(index, level.jValue, Qt::EditRole);
    index = levelsModel->index(levels.size(), 3, QModelIndex());
    levelsModel->setData(index, level.piValue, Qt::EditRole);
    index = levelsModel->index(levels.size(), 4, QModelIndex());
    levelsModel->setData(index, level.energy, Qt::EditRole);
    levelsModel->setNucleus(levels.size(), level.nucleus());
    levelsModel->setFileId(levels.size(), level.fileId);
    levelsView->resizeRowsToContents();
    if (!fromFile) updateChannelsLevelAdded(levels.size());
    updateMissingSpinParities();
  } else {
    QMessageBox::information(this, tr("Duplicate Level"), tr("This level already exists."));
  }
}

void LevelsTab::removeLevelAt(int row) {
  levelsView->selectionModel()->select(levelsModelProxy->mapFromSource(levelsModel->index(row, 0)),
                                       QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
  removeLevel();
}

void LevelsTab::removeLevel() {
  // the button "-" and removeLevelAt: the selected levels, without a question for one (as dev)
  const QList<int> rows = selectedLevelRows();
  if (!rows.isEmpty()) removeLevels(rows, rows.size() > 1);
}

void LevelsTab::editLevel() {
  QItemSelectionModel *selectionModel = levelsView->selectionModel();
  QModelIndexList indexes = selectionModel->selectedRows();
  if (indexes.isEmpty()) return;
  QModelIndex index = levelsModelProxy->mapToSource(levelsView->currentIndex().isValid() ? levelsView->currentIndex() : indexes[0]);
  undoText_ = tr("edit level");

  QModelIndex i = levelsModel->index(index.row(), 2, QModelIndex());
  QVariant var = levelsModel->data(i, Qt::EditRole);
  double jValue = var.toDouble();
  i = levelsModel->index(index.row(), 3, QModelIndex());
  var = levelsModel->data(i, Qt::EditRole);
  int piValue = var.toInt();
  i = levelsModel->index(index.row(), 4, QModelIndex());
  var = levelsModel->data(i, Qt::EditRole);
  QString energy = var.toString();

  Nucleus nucleus = levelsModel->getLevels().at(index.row()).nucleus();

  AddLevelDialog aDialog;
  aDialog.setWindowTitle(tr("Edit a Level"));
  aDialog.setSpinParityOptions(spinParityOptions(nucleus), jValue, piValue, true);
  aDialog.energyText->setText(energy);
  // The nucleus is shown but fixed: the level's channels belong to it.
  QList<Nucleus> nuclei = pairsModel ? pairsModel->nuclei() : QList<Nucleus>();
  QStringList labels;
  QList<QList<SpinParityOption> > perNucleus;
  for (const Nucleus &n : nuclei) {
    labels << pairsModel->nucleusLongLabel(n);
    perNucleus << spinParityOptions(n);
  }
  aDialog.setNucleusOptions(labels, perNucleus, std::max(0, pairsModel ? pairsModel->nucleusIndex(nucleus) : 0), false);

  if (aDialog.exec()) {
    double newJValue = aDialog.jValue();
    if (newJValue != jValue) {
      i = levelsModel->index(index.row(), 2, QModelIndex());
      levelsModel->setData(i, newJValue, Qt::EditRole);
    }
    int newPiValue = aDialog.piValue();
    if (newPiValue != piValue) {
      i = levelsModel->index(index.row(), 3, QModelIndex());
      levelsModel->setData(i, newPiValue, Qt::EditRole);
    }
    QString newEnergy = aDialog.energyText->text();
    if (newEnergy != energy) {
      i = levelsModel->index(index.row(), 4, QModelIndex());
      levelsModel->setData(i, newEnergy, Qt::EditRole);
    }
    updateChannelsLevelEdited(index.row());
    updateMissingSpinParities();
  }
}

void LevelsTab::updateButtons(const QItemSelection &selection) {
  QModelIndexList indexes = selection.indexes();

  if (indexes.isEmpty()) {
    removeLevelButton->setEnabled(false);
  } else {
    removeLevelButton->setEnabled(true);
  }
}

void LevelsTab::updateChannelsLevelAdded(int levelIndex) {
  QList<ChannelsData> newChannels = calculateChannels(levelIndex);
  channelsModel->insertRows(0, newChannels.size(), QModelIndex());
  for (int i = 0; i < newChannels.size(); i++) {
    QModelIndex index = channelsModel->index(i, 0, QModelIndex());
    channelsModel->setData(index, newChannels.at(i).isFixed, Qt::EditRole);
    index = channelsModel->index(i, 1, QModelIndex());
    channelsModel->setData(index, newChannels.at(i).levelIndex, Qt::EditRole);
    index = channelsModel->index(i, 2, QModelIndex());
    channelsModel->setData(index, newChannels.at(i).pairIndex, Qt::EditRole);
    index = channelsModel->index(i, 3, QModelIndex());
    channelsModel->setData(index, newChannels.at(i).sValue, Qt::EditRole);
    index = channelsModel->index(i, 4, QModelIndex());
    channelsModel->setData(index, newChannels.at(i).lValue, Qt::EditRole);
    index = channelsModel->index(i, 5, QModelIndex());
    channelsModel->setData(index, newChannels.at(i).radType, Qt::EditRole);
    index = channelsModel->index(i, 6, QModelIndex());
    channelsModel->setData(index, newChannels.at(i).reducedWidth, Qt::EditRole);
  }
  channelsView->resizeRowsToContents();
}

void LevelsTab::updateChannelsLevelDeleted(int levelIndex) {
  QList<ChannelsData> channels = channelsModel->getChannels();
  int deleted = 0;
  for (int i = 0; i < channels.size(); i++) {
    if (channels.at(i).levelIndex == levelIndex) {
      channelsModel->removeRows(i - deleted, 1, QModelIndex());
      deleted++;
    }
  }
  channels = channelsModel->getChannels();
  for (int i = 0; i < channels.size(); i++) {
    if (channels.at(i).levelIndex > levelIndex) {
      QModelIndex index = channelsModel->index(i, 1, QModelIndex());
      channelsModel->setData(index, channels.at(i).levelIndex - 1, Qt::EditRole);
    }
  }
}

void LevelsTab::updateChannelsLevelEdited(int levelIndex) {
  QList<ChannelsData> channels = channelsModel->getChannels();
  int deleted = 0;
  for (int i = 0; i < channels.size(); i++) {
    if (channels.at(i).levelIndex == levelIndex) {
      channelsModel->removeRows(i - deleted, 1, QModelIndex());
      deleted++;
    }
  }
  QList<ChannelsData> newChannels = calculateChannels(levelIndex);
  channelsModel->insertRows(0, newChannels.size(), QModelIndex());
  for (int i = 0; i < newChannels.size(); i++) {
    QModelIndex index = channelsModel->index(i, 1, QModelIndex());
    channelsModel->setData(index, newChannels.at(i).levelIndex, Qt::EditRole);
    index = channelsModel->index(i, 2, QModelIndex());
    channelsModel->setData(index, newChannels.at(i).pairIndex, Qt::EditRole);
    index = channelsModel->index(i, 3, QModelIndex());
    channelsModel->setData(index, newChannels.at(i).sValue, Qt::EditRole);
    index = channelsModel->index(i, 4, QModelIndex());
    channelsModel->setData(index, newChannels.at(i).lValue, Qt::EditRole);
    index = channelsModel->index(i, 5, QModelIndex());
    channelsModel->setData(index, newChannels.at(i).radType, Qt::EditRole);
    for (int ii = 0; ii < channels.size(); ii++) {
      if (channels.at(ii).levelIndex == newChannels.at(i).levelIndex &&
          channels.at(ii).pairIndex == newChannels.at(i).pairIndex &&
          channels.at(ii).sValue == newChannels.at(i).sValue &&
          channels.at(ii).lValue == newChannels.at(i).lValue &&
          channels.at(ii).radType == newChannels.at(i).radType) {
        newChannels[i].isFixed = channels.at(ii).isFixed;
        newChannels[i].reducedWidth = channels.at(ii).reducedWidth;
        channelsModel->setFileTokens(i, channels.at(ii).fileTokens);  // a value kept keeps its token
        break;
      }
    }
    index = channelsModel->index(i, 6, QModelIndex());
    channelsModel->setData(index, newChannels.at(i).reducedWidth, Qt::EditRole);
    index = channelsModel->index(i, 0, QModelIndex());
    channelsModel->setData(index, newChannels.at(i).isFixed, Qt::EditRole);
  }
  channelsView->resizeRowsToContents();
}

QList<ChannelsData> LevelsTab::calculateChannels(int levelIndex) {
  return calculateChannels((levelsModel->getLevels()).at(levelIndex), levelIndex);
}

QList<SpinParityOption> LevelsTab::spinParityOptions() {
  return spinParityOptions(Nucleus{0, 0});
}

QList<SpinParityOption> LevelsTab::spinParityOptions(const Nucleus &nucleus) {
  // Every J^pi that gives at least one channel under the current pairs and limits, by trying each candidate through
  // calculateChannels itself, so the list cannot drift from the channels the tab will actually build.
  QList<SpinParityOption> options;
  QList<PairsData> pairs = pairsModel->getPairs();
  bool filter = nucleus.Z > 0 || nucleus.A > 0;
  double maxJ = 0.;
  for (const PairsData &pair : pairs) {
    if (filter && PairsModel::nucleusOf(pair) != nucleus) continue;
    double pairMaxJ;
    if (pair.pairType == 0)
      pairMaxJ = pair.heavyJ + pair.lightJ + maxLSpin->value();
    else if (pair.pairType == 20)
      pairMaxJ = pair.heavyJ + 1.;
    else
      pairMaxJ = pair.heavyJ + maxMultSpin->value();
    maxJ = std::max(maxJ, pairMaxJ);
  }
  for (int twoJ = 0; twoJ <= int(2. * maxJ + 0.5); twoJ++) {
    for (int pi : {1, -1}) {
      LevelsData level = {1, 0, twoJ / 2.0, pi, 0., nucleus.Z, nucleus.A};
      int numChannels = calculateChannels(level, -1).size();
      if (numChannels > 0) options.append({level.jValue, pi, numChannels});
    }
  }
  return options;
}

QList<ChannelsData> LevelsTab::calculateChannels(const LevelsData &level, int levelIndex) {
  QList<ChannelsData> channels;
  const QList<PairsData> pairs = pairsModel->getPairs();
  // one pair row per pair, as the tab writes a pair (channelRow), the key its place + 1; the level's nucleus's keys
  std::vector<azure2::front::LevelRow> rows;
  std::vector<int> keys;
  const LevelsData pairLevel = {0, 0, 0., 1, 0.};
  for (int i = 0; i < pairs.size(); i++) {
    const ChannelsData ch = {0, -1, i, 0., 0, 'P', 0.};
    rows.push_back(channelRow(ch, pairLevel, pairs.at(i), pairs.at(0), 0));
    if (!level.hasNucleus() || PairsModel::nucleusOf(pairs.at(i)) == level.nucleus()) keys.push_back(i + 1);
  }
  if (rows.empty() || keys.empty()) return channels;
  std::vector<azure2::front::AllowedChannelInfo> allowed;
  std::string why;
  const int twoJ = int(std::lround(2.0 * level.jValue));
  if (!azure2::front::AllowedChannels(rows, twoJ, level.piValue, maxLSpin->value(), maxMultSpin->value(),
                                      maxNumMultSpin->value(), keys, allowed, why))
    return channels;
  for (const azure2::front::AllowedChannelInfo &c : allowed) {
    const int i = c.pairKey - 1;
    if (i < 0 || i >= pairs.size()) continue;
    const PairsData &pair = pairs.at(i);
    QChar radType = 'P';
    if (pair.pairType == 20) {
      radType = c.l == 1 ? 'G' : 'F';  // Gamow-Teller, Fermi
    } else if (pair.pairType == 10) {
      // E or M from the parities: an odd L with a change of parity is electric, an even L without one too
      const int parityChange = pair.heavyPi * level.piValue;
      radType = (c.l % 2 != 0) == (parityChange == -1) ? 'E' : 'M';
    }
    channels.push_back({0, levelIndex, i, c.twoS / 2.0, c.l, radType, 0.0});
  }
  return channels;
}

void LevelsTab::updateFilter(const QItemSelection &selection) {
  Q_UNUSED(selection);
  // the channels of the current level (of a multi-selection: the one the cursor is on)
  const int row = currentLevelRow();
  proxyModel->setFilterRegularExpression(row < 0 ? QString("^-1$") : QString("^%1$").arg(row));
  channelsView->resizeRowsToContents();
  refreshChannelInfo();
}

QStringList LevelsTab::missingSpinParities(QStringList *links) {
  QStringList missing;
  QList<LevelsData> levels = levelsModel->getLevels();
  // One pass per compound nucleus (the J^pi of one nucleus's pairs against its own levels); a single nucleus is
  // reported as before, without a prefix.
  QList<Nucleus> nuclei = pairsModel ? pairsModel->nuclei() : QList<Nucleus>();
  bool several = nuclei.size() > 1;
  if (!several) nuclei = QList<Nucleus>() << Nucleus{0, 0};
  for (const Nucleus &nucleus : nuclei) {
    QString prefix = several ? PairsModel::nucleusLabel(nucleus) + " " : QString();
    for (const SpinParityOption &option : spinParityOptions(nucleus)) {
      // only J^pi with a particle channel scatter; a capture-only J^pi has nothing to add without a level
      LevelsData probe = {1, 0, option.jValue, option.piValue, 0., nucleus.Z, nucleus.A};
      bool particle = false;
      for (const ChannelsData &c : calculateChannels(probe, -1)) particle |= (c.radType == 'P');
      if (!particle) continue;
      bool carried = false;
      for (const LevelsData &l : levels)
        carried |= (l.isActive == 1 && l.jValue == option.jValue && l.piValue == option.piValue &&
                    (!several || !l.hasNucleus() || l.nucleus() == nucleus));
      if (!carried) {
        missing << prefix + AddLevelDialog::spinParityLabel(option.jValue, option.piValue);
        if (links)
          *links << QString("jpi:%1:%2:%3:%4").arg((int)std::lround(2 * option.jValue)).arg(option.piValue).arg(nucleus.Z).arg(nucleus.A);
      }
    }
  }
  return missing;
}

void LevelsTab::updateMissingSpinParities() {
  // What is missing and what to do about it (docs/upgrade/qt_gui.md section 11.3): each J^pi a link that fills the
  // add bar.
  QStringList links;
  const QStringList missing = missingSpinParities(&links);
  QStringList linked;
  for (int i = 0; i < missing.size(); i++) linked << QString("<a href='%1'>%2</a>").arg(links.value(i), missing.at(i));
  QString jpi = QString("J<sup>%1</sup>").arg(QChar(0x03C0));
  if (pairsModel == nullptr || pairsModel->getPairs().isEmpty())
    missingJPiLabel->setText(tr("Add a particle pair first (tab Particle Pairs): the levels' channels come from the pairs."));
  else if (missing.isEmpty())
    missingJPiLabel->setText(tr("Every %1 up to l = %2 has an active level.").arg(jpi).arg(maxLSpin->value()));
  else if (includeAllJPiCheck->isChecked())
    missingJPiLabel->setText(tr("Hard sphere only (no level): %1").arg(linked.join(", ")));
  else
    missingJPiLabel->setText(tr("<font color='#b00020'>Not calculated (no level): %1</font>").arg(linked.join(", ")));
  // the explanation in the tip, not on the page (qt_gui.md section 19.6)
  missingJPiLabel->setToolTip(
      includeAllJPiCheck->isChecked()
          ? tr("No active level carries these J^pi: they scatter as hard spheres only (the option above).  Click one to "
               "add a level of it.")
          : tr("No active level carries these J^pi and the option above is off, so they have no hard-sphere phase "
               "either.  To include them, switch the option on, or click one to add a level of it (a background pole, "
               "far from the data with its energy fixed, is the usual choice)."));
}

void LevelsTab::updateChannelsPairAddedEdited() {
  QList<LevelsData> levels = levelsModel->getLevels();
  QList<ChannelsData> channels = channelsModel->getChannels();
  channelsModel->removeRows(0, channels.size(), QModelIndex());
  for (int levelIndex = 0; levelIndex < levels.size(); levelIndex++) {
    QList<ChannelsData> newChannels = calculateChannels(levelIndex);
    channelsModel->insertRows(0, newChannels.size(), QModelIndex());
    for (int i = 0; i < newChannels.size(); i++) {
      QModelIndex index = channelsModel->index(i, 1, QModelIndex());
      channelsModel->setData(index, newChannels.at(i).levelIndex, Qt::EditRole);
      index = channelsModel->index(i, 2, QModelIndex());
      channelsModel->setData(index, newChannels.at(i).pairIndex, Qt::EditRole);
      index = channelsModel->index(i, 3, QModelIndex());
      channelsModel->setData(index, newChannels.at(i).sValue, Qt::EditRole);
      index = channelsModel->index(i, 4, QModelIndex());
      channelsModel->setData(index, newChannels.at(i).lValue, Qt::EditRole);
      index = channelsModel->index(i, 5, QModelIndex());
      channelsModel->setData(index, newChannels.at(i).radType, Qt::EditRole);
      for (int ii = 0; ii < channels.size(); ii++) {
        if (channels.at(ii).levelIndex == newChannels.at(i).levelIndex &&
            channels.at(ii).pairIndex == newChannels.at(i).pairIndex &&
            channels.at(ii).sValue == newChannels.at(i).sValue &&
            channels.at(ii).lValue == newChannels.at(i).lValue &&
            channels.at(ii).radType == newChannels.at(i).radType) {
          newChannels[i].reducedWidth = channels.at(ii).reducedWidth;
          newChannels[i].isFixed = channels.at(ii).isFixed;
          channelsModel->setFileTokens(i, channels.at(ii).fileTokens);  // a value kept keeps its token
          break;
        }
      }
      index = channelsModel->index(i, 6, QModelIndex());
      channelsModel->setData(index, newChannels.at(i).reducedWidth, Qt::EditRole);
      index = channelsModel->index(i, 0, QModelIndex());
      channelsModel->setData(index, newChannels.at(i).isFixed, Qt::EditRole);
    }
  }
  channelsView->resizeRowsToContents();
  updateMissingSpinParities();
  updateExpertRows();  // a photon pair added in standard view: its multipolarity limits shown at once
}

void LevelsTab::updateChannelsPairRemoved(int pairIndex) {
  QList<ChannelsData> channels = channelsModel->getChannels();
  int deleted = 0;
  for (int i = 0; i < channels.size(); i++) {
    if (channels.at(i).pairIndex == pairIndex) {
      channelsModel->removeRows(i - deleted, 1, QModelIndex());
      deleted++;
    }
  }
  channels = channelsModel->getChannels();
  for (int i = 0; i < channels.size(); i++) {
    if (channels.at(i).pairIndex > pairIndex) {
      QModelIndex index = channelsModel->index(i, 2, QModelIndex());
      channelsModel->setData(index, channels.at(i).pairIndex - 1, Qt::EditRole);
    }
  }
  channelsView->resizeRowsToContents();
  updateExpertRows();
}

void LevelsTab::updateDetails(const QItemSelection &selection) {
  QModelIndexList indexes = selection.indexes();

  if (indexes.isEmpty()) {
    channelDetails->hide();
  } else {
    QModelIndex index = proxyModel->mapToSource(indexes.at(0));

    QModelIndex i = channelsModel->index(index.row(), 1, QModelIndex());
    QVariant var = channelsModel->data(i, Qt::EditRole);
    int levelIndex = var.toInt();
    i = channelsModel->index(index.row(), 2, QModelIndex());
    var = channelsModel->data(i, Qt::EditRole);
    int pairIndex = var.toInt();

    ChannelsData channel = channelsModel->getChannels().at(index.row());
    LevelsData level = levelsModel->getLevels().at(levelIndex);
    PairsData pair = pairsModel->getPairs().at(pairIndex);

    QString details = "";
    QTextStream stm(&details, QIODevice::Append);
    stm << QString("%1 MeV level with spin %2\n   transitioning via pair key #%3").arg(level.energy).arg(levelsModel->getSpinLabel(level)).arg(pairIndex + 1)
        << Qt::endl;
    if (channel.radType == 'P') {
      stm << QString("Channel configuration is\n   s = %1, l = %2").arg(channelsModel->getSpinLabel(channel)).arg(channel.lValue)
          << Qt::endl
          << Qt::endl;
      stm << qSetFieldWidth(21) << Qt::right << "Light Particle Spin: "
          << qSetFieldWidth(0) << Qt::left << QString("%1").arg(pairsModel->getSpinLabel(pair, 0)) << Qt::endl;
      stm << qSetFieldWidth(21) << Qt::right << "Light Particle Z: "
          << qSetFieldWidth(0) << Qt::left << QString("%1").arg(pair.lightZ) << Qt::endl;
      stm << qSetFieldWidth(21) << Qt::right << "Light Particle M: "
          << qSetFieldWidth(0) << Qt::left << QString("%1").arg(pair.lightM) << Qt::endl;
      stm << qSetFieldWidth(21) << Qt::right << "Light Particle G: "
          << qSetFieldWidth(0) << Qt::left << QString("%1").arg(pair.lightG) << Qt::endl;
    } else if (channel.radType == 'G' || channel.radType == 'F') {
      if (channel.radType == 'G')
        stm << QString("Channel is Gamow-Teller beta decay") << Qt::endl
            << Qt::endl;
      else
        stm << QString("Channel is Fermi beta decay") << Qt::endl
            << Qt::endl;
      stm << qSetFieldWidth(21) << Qt::right << "Fermion Charge: "
          << qSetFieldWidth(0) << Qt::left << QString("%1").arg(pair.lightZ) << Qt::endl;
    } else {
      stm << QString("Capture gamma is %1%2 radiation").arg(channel.radType).arg(channel.lValue) << Qt::endl;
      if (((channel.radType == 'E' && channel.lValue == 1) &&
           (pair.ecMultMask & (1 << 0))) ||
          ((channel.radType == 'M' && channel.lValue == 1) &&
           (pair.ecMultMask & (1 << 1))) ||
          ((channel.radType == 'E' && channel.lValue == 2) &&
           (pair.ecMultMask & (1 << 2))))
        stm << "Capture is internal and external" << Qt::endl;
      else
        stm << "Capture is internal only" << Qt::endl;
      stm << Qt::endl;
    }
    stm << qSetFieldWidth(21) << Qt::right << "Heavy Particle Spin: "
        << qSetFieldWidth(0) << Qt::left << QString("%1").arg(pairsModel->getSpinLabel(pair, 1)) << Qt::endl;
    stm << qSetFieldWidth(21) << Qt::right << "Heavy Particle Z: "
        << qSetFieldWidth(0) << Qt::left << QString("%1").arg(pair.heavyZ) << Qt::endl;
    stm << qSetFieldWidth(21) << Qt::right << "Heavy Particle M: "
        << qSetFieldWidth(0) << Qt::left << QString("%1").arg(pair.heavyM) << Qt::endl;
    stm << qSetFieldWidth(21) << Qt::right << "Heavy Particle G: "
        << qSetFieldWidth(0) << Qt::left << QString("%1").arg(pair.heavyG) << Qt::endl;
    if (channel.radType != 'G' && channel.radType != 'F')
      stm << qSetFieldWidth(21) << Qt::right << "Excitation Energy: "
          << qSetFieldWidth(0) << Qt::left << QString("%1").arg(pair.excitationEnergy) << Qt::endl;
    if (pair.lightExcitation != 0.)
      stm << qSetFieldWidth(21) << Qt::right << "Light Excitation: "
          << qSetFieldWidth(0) << Qt::left << QString("%1").arg(pair.lightExcitation) << Qt::endl;
    // the channel at its level's energy, from the facade (azure2::front::ChannelAtLevel: the engine's channel
    // functions; the channel's own radius from <channelRadii>): docs/upgrade/qt_gui.md section 11.4
    const azure2::front::ChannelQuantities q = channelQuantities(index.row());
    const double channelRadius = q.radius;
    auto line = [&stm](const QString &label, const QString &value) {
      stm << qSetFieldWidth(21) << Qt::right << label << qSetFieldWidth(0) << Qt::left << value << Qt::endl;
    };
    auto num = [](double v) { return std::isfinite(v) ? QString::number(v, 'g', 6) : QString("-"); };
    if (channel.radType != 'M' && channel.radType != 'E') {
      line("Separation Energy: ", QString("%1").arg(pair.seperationEnergy));
      line("Channel Radius: ", QString("%1 fm").arg(channelRadius) +
                                   (channelRadius != pair.channelRadius ? QString(" (own; pair %1)").arg(pair.channelRadius)
                                                                        : QString()));
    }
    if (q.error.empty()) {
      stm << Qt::endl;
      line("Value is: ", QString::fromStdString(q.symbol + " [" + q.unit + "], " + q.meaning));
      if (q.particle) {
        line("E_c = E_x - S - e: ", num(q.channelEnergy) + " MeV");
        line(q.channelEnergy > 0 ? "Penetrability P: " : "Bound factor P: ",
             num(q.P) + (q.channelEnergy > 0 ? QString() : QString(" MeV^-1 fm^-1")));
        line("Shift S: ", num(q.S));
        line("dS/dE: ", num(q.dSdE) + " MeV^-1");
        line("Wigner limit gW^2: ", num(q.wignerLimit) + " MeV");
        if (std::isfinite(q.singleParticleWidth)) line("Single-particle G: ", num(q.singleParticleWidth * 1e6) + " eV");
        line("Reduced width g^2: ", num(q.reducedWidth) + " MeV");
        line("theta^2 = g^2/gW^2: ", num(q.wignerRatio));
      }
    } else {
      line("Not read: ", QString::fromStdString(q.error));
    }
    stm.flush();
    channelDetails->details->setText(details);
    const QString kind = QString::fromStdString(q.kind);
    if (kind == "closed") channelDetails->setNormParam(1);
    else if (kind == "magnetic_moment") channelDetails->setNormParam(2);
    else if (kind == "quadrupole_moment") channelDetails->setNormParam(3);
    else if (kind == "beta") channelDetails->setNormParam(4);
    else channelDetails->setNormParam(0);
    channelDetails->reducedWidthText->setText(QString("%1").arg(channel.reducedWidth));
    // the Wigner limit button: the single-particle width of an open particle channel (eV)
    wignerApplicable_ = q.particle && std::isfinite(q.singleParticleWidth);
    wignerEcm_ = q.singleParticleWidth;  // MeV, shown by calculateWignerLimit
    channelDetails->wignerButton->setEnabled(wignerApplicable_);
    channelDetails->wignerLimitText->clear();
    channelDetails->wignerLimitText->setPlaceholderText(wignerApplicable_ ? QString() : tr("N/A"));
    if (wignerApplicable_) calculateWignerLimit();
    channelDetails->show();
  }
}

void LevelsTab::updateReducedWidth(const QString &string) {
  QItemSelectionModel *selectionModel = channelsView->selectionModel();
  QModelIndex index = proxyModel->mapToSource(selectionModel->selectedRows().at(0));

  if (index.isValid()) {
    double reducedWidth = string.toDouble();
    QModelIndex i = channelsModel->index(index.row(), 6, QModelIndex());
    channelsModel->setData(i, reducedWidth, Qt::EditRole);
  }
}

// Physical Wigner limit (single-particle partial width) of the selected channel:
// Gamma_W = 2 P(E) gamma^2_W, with the reduced limit gamma^2_W = 3 hbar^2/(2 mu a^2)
// (Lane & Thomas 1958) and the penetrability P computed as in CoulFunc.  Result
// in eV, so it matches the partial-width box.  Only meaningful for unbound
// particle channels.
void LevelsTab::calculateWignerLimit() {
  if (!wignerApplicable_) {
    channelDetails->wignerLimitText->setText(tr("N/A"));
    return;
  }
  // the facade's single-particle width of the channel at its level (azure2::front::ChannelAtLevel), in eV
  const double gammaW_eV = wignerEcm_ * 1.0e6;
  if (!std::isfinite(gammaW_eV)) {
    channelDetails->wignerLimitText->setText(tr("error"));
    return;
  }
  double av = std::fabs(gammaW_eV);
  QString text = (av != 0. && (av >= 1.e4 || av < 1.e-4)) ? QString::number(gammaW_eV, 'e', 4) : QString::number(gammaW_eV, 'f', 4);
  channelDetails->wignerLimitText->setText(text);
}

QList<int> LevelsTab::writtenLevelOrder() {
  QList<LevelsData> levels = levelsModel->getLevels();
  double lowJ = 0;
  double highJ = 0;
  for (int la = 0; la < levels.size(); la++) {
    double tempJ = levels.at(la).jValue;
    if (la == 0) {
      lowJ = tempJ;
      highJ = tempJ;
    } else if (tempJ <= lowJ)
      lowJ = tempJ;
    else if (tempJ >= highJ)
      highJ = tempJ;
  }
  QList<QList<LevelsData>> sortedLevels;
  QList<QList<int>> levelsMap;
  // One block per compound nucleus, in order of first appearance, each sorted by J, parity and energy as before
  // (a single-nucleus project is written exactly as it always was).
  QList<Nucleus> nucleusOrder;
  for (const LevelsData &level : levels)
    if (!nucleusOrder.contains(level.nucleus())) nucleusOrder.append(level.nucleus());
  for (const Nucleus &nucleus : nucleusOrder)
  for (double j = lowJ; j <= highJ; j += 0.5) {
    for (int pi = -1; pi <= 1; pi += 2) {
      for (int la = 0; la < levels.size(); la++) {
        if (levels.at(la).jValue == j && levels.at(la).piValue == pi && levels.at(la).nucleus() == nucleus) {
          if (sortedLevels.size() == 0 || (sortedLevels.at(sortedLevels.size() - 1).at(0).jValue != j || sortedLevels.at(sortedLevels.size() - 1).at(0).piValue != pi || sortedLevels.at(sortedLevels.size() - 1).at(0).nucleus() != nucleus)) {
            QList<LevelsData> tempLevelList;
            sortedLevels.append(tempLevelList);
            QList<int> tempKeyList;
            levelsMap.append(tempKeyList);
          }
          if (sortedLevels.at(sortedLevels.size() - 1).size() == 0) {
            sortedLevels[sortedLevels.size() - 1].append(levels.at(la));
            levelsMap[levelsMap.size() - 1].append(la);
          } else
            for (int mu = 0; mu < sortedLevels.at(sortedLevels.size() - 1).size(); mu++) {
              if (levels.at(la).energy <= sortedLevels.at(sortedLevels.size() - 1).at(mu).energy) {
                sortedLevels[sortedLevels.size() - 1].insert(mu, levels.at(la));
                levelsMap[levelsMap.size() - 1].insert(mu, la);
                break;
              } else if (mu == sortedLevels.at(sortedLevels.size() - 1).size() - 1) {
                sortedLevels[sortedLevels.size() - 1].append(levels.at(la));
                levelsMap[levelsMap.size() - 1].append(la);
                break;
              }
            }
        }
      }
    }
  }
  QList<int> order;
  for (const QList<int> &group : levelsMap) order.append(group);
  // An isospin level (one level ID whose lines span several compound nuclei; the <isospin> records name it): the
  // GUI's levels that were read under one ID are written together, consecutively, under that ID (levelRows).
  QMap<int, int> membersOfId;
  for (int row : order)
    if (levels.at(row).fileId > 0) membersOfId[levels.at(row).fileId]++;
  QList<int> grouped;
  QSet<int> written;
  for (int row : order) {
    const int fileId = levels.at(row).fileId;
    if (fileId <= 0 || membersOfId.value(fileId) < 2) {
      grouped.append(row);
      continue;
    }
    if (written.contains(fileId)) continue;
    written.insert(fileId);
    for (int r : order)
      if (levels.at(r).fileId == fileId) grouped.append(r);
  }
  return grouped;
}

LevelsTab::EngineNumbering LevelsTab::engineNumbering(const QList<int> &order) {
  QList<LevelsData> levels = levelsModel->getLevels();
  QList<ChannelsData> channels = channelsModel->getChannels();
  EngineNumbering out;
  out.level.fill(0, levels.size());
  out.channel.fill(0, channels.size());
  // J-groups in order of first appearance among the active levels, each with its levels in order.
  struct Group {
    Nucleus nucleus;
    double j;
    int pi;
    QList<int> rows;
  };
  QList<Group> groups;
  for (int row : order) {
    if (row < 0 || row >= levels.size() || levels.at(row).isActive != 1) continue;
    const LevelsData &l = levels.at(row);
    int g = 0;
    while (g < groups.size() && !(groups[g].nucleus == l.nucleus() && groups[g].j == l.jValue && groups[g].pi == l.piValue)) g++;
    if (g == groups.size()) groups.append({l.nucleus(), l.jValue, l.piValue, {}});
    groups[g].rows.append(row);
  }
  int n = 0;
  for (const Group &group : groups) {
    // The group's channels (pair, l, s) in order of first appearance over its levels' channel lines.
    QList<ChannelsData> keys;
    for (int row : group.rows) {
      out.level[row] = ++n;
      for (int ch = 0; ch < channels.size(); ch++) {
        const ChannelsData &c = channels.at(ch);
        if (c.levelIndex != row) continue;
        int m = 0;
        while (m < keys.size() && !(keys[m].pairIndex == c.pairIndex && keys[m].lValue == c.lValue && keys[m].sValue == c.sValue)) m++;
        if (m == keys.size()) keys.append(c);
        out.channel[ch] = m + 1;
      }
    }
  }
  return out;
}

void LevelsTab::reset() {
  // a project is being read: no undo entries until the event loop turns (finishLoading)
  if (!loading_) {
    loading_ = true;
    QTimer::singleShot(0, this, &LevelsTab::finishLoading);
  }
  pairsModel->removeRows(0, pairsModel->getPairs().size(), QModelIndex());
  levelsModel->removeRows(0, levelsModel->getLevels().size(), QModelIndex());
  channelsModel->removeRows(0, channelsModel->getChannels().size(), QModelIndex());
  includeAllJPiCheck->setChecked(true);
  missingJPiLabel->setText("");
  maxLSpin->blockSignals(true);
  maxMultSpin->blockSignals(true);
  maxNumMultSpin->blockSignals(true);
  maxLSpin->setValue(2);
  maxMultSpin->setValue(2);
  maxNumMultSpin->setValue(2);
  maxLSpin->blockSignals(false);
  maxMultSpin->blockSignals(false);
  maxNumMultSpin->blockSignals(false);
}

void LevelsTab::showInfo(int which, QString title) {
  if (which < infoText.size()) {
    if (!infoDialog[which]) {
      infoDialog[which] = new InfoDialog(infoText[which], this, title);
      infoDialog[which]->setAttribute(Qt::WA_DeleteOnClose);
      infoDialog[which]->show();
    } else
      infoDialog[which]->raise();
  }
}

void LevelsTab::fixAllWidths() {
  // several levels selected: theirs; else every level (as dev)
  if (selectedLevelRows().size() > 1) {
    setFixed(selectedLevelRows(), false, true);
    return;
  }
  undoText_ = tr("fix all widths");
  QList<ChannelsData> channels = channelsModel->getChannels();
  for (int i = 0; i < channels.size(); i++) {
    QModelIndex index = channelsModel->index(i, 0, QModelIndex());
    channelsModel->setData(index, 1, Qt::EditRole);
  }
}

void LevelsTab::fixAllEnergies() {
  if (selectedLevelRows().size() > 1) {
    setFixed(selectedLevelRows(), true, true);
    return;
  }
  undoText_ = tr("fix all energies");
  QList<LevelsData> levels = levelsModel->getLevels();
  for (int i = 0; i < levels.size(); i++) {
    QModelIndex index = levelsModel->index(i, 1, QModelIndex());
    levelsModel->setData(index, 1, Qt::EditRole);
  }
}

void LevelsTab::addLevelsFromIAEA() {
  // Check if pairs are defined
  if (pairsModel->getPairs().isEmpty()) {
    QMessageBox::warning(this, tr("No Pairs Defined"),
                         tr("Please define at least one particle pair before querying IAEA database.\n"
                            "The compound nucleus is calculated from the entrance channel."));
    return;
  }

  // The compound nucleus from the pairs; with several, the user picks which one to import.
  QList<Nucleus> nuclei = pairsModel->nuclei();
  if (nuclei.isEmpty()) {
    QMessageBox::warning(this, tr("Error"), tr("Could not determine compound nucleus from pairs."));
    return;
  }
  Nucleus chosen = nuclei.first();
  if (nuclei.size() > 1) {
    QStringList labels;
    for (const Nucleus &nucleus : nuclei) labels << pairsModel->nucleusLongLabel(nucleus);
    bool picked = false;
    QString pick = QInputDialog::getItem(this, tr("Compound Nucleus"), tr("Import the levels of which compound nucleus?"),
                                         labels, 0, false, &picked);
    if (!picked) return;
    chosen = nuclei.at(std::max<qsizetype>(0, labels.indexOf(pick)));
  }
  int massNumber = chosen.A;
  int atomicNumber = chosen.Z;

  // Create and show the IAEA dialog
  IAEALevelDialog dialog(this);
  dialog.setCompoundNucleus(massNumber, atomicNumber);

  if (dialog.exec() == QDialog::Accepted) {
    QList<LevelsData> newLevels = dialog.getSelectedLevels();

    int addedCount = 0;
    int duplicateCount = 0;

    for (LevelsData level : newLevels) {
      level.nucZ = chosen.Z;
      level.nucA = chosen.A;
      if (levelsModel->isLevel(level) == -1) {
        addLevel(level, false);
        addedCount++;
      } else {
        duplicateCount++;
      }
    }

    QString message;
    if (addedCount > 0) {
      message = tr("Added %1 level(s) from IAEA database.").arg(addedCount);
      if (duplicateCount > 0) {
        message += tr("\n%1 duplicate level(s) were skipped.").arg(duplicateCount);
      }
      QMessageBox::information(this, tr("Levels Added"), message);
    } else if (duplicateCount > 0) {
      QMessageBox::information(this, tr("No Levels Added"),
                               tr("All selected levels already exist."));
    }
  }
}

void LevelsTab::exportLatexTable() {
  QList<LevelsData> levels = levelsModel->getLevels();
  QList<ChannelsData> channels = channelsModel->getChannels();
  QList<PairsData> pairs = pairsModel->getPairs();

  if (levels.isEmpty()) {
    QMessageBox::warning(this, tr("No Levels"),
                         tr("There are no levels to export."));
    return;
  }

  // Helper lambda to get LaTeX particle label (no HTML)
  auto getLatexParticleLabel = [](const PairsData &pair, int which) -> QString {
    if (which != -1) {
      if (pair.pairType == 10 && which == 0)
        return "$\\gamma$";
      else if (pair.pairType == 20 && which == 0) {
        if (pair.lightZ < 0)
          return "$\\beta^-$";
        else
          return "$\\beta^+$";
      } else {
        int tempZ = (which == 0) ? pair.lightZ : pair.heavyZ;
        int tempM = (which == 0) ? (int)round(pair.lightM) : (int)round(pair.heavyM);
        std::map<int, QString>::const_iterator it = elementMap.find(tempZ);
        if (it != elementMap.end()) {
          if (tempM == 1) {
            if (tempZ == 1)
              return "$p$";
            else
              return QString("$%1$").arg(it->second);
          } else if (tempZ == 2 && tempM == 4)
            return "$\\alpha$";
          else
            return QString("$^{%1}$%2").arg(tempM).arg(it->second);
        } else
          return "?";
      }
    }
    return "?";
  };

  // Helper lambda to format spin label for LaTeX
  auto formatSpinLabel = [](double spin) -> QString {
    if (spin == (int)spin) {
      return QString::number((int)spin);
    } else {
      return QString("%1/2").arg((int)(spin * 2));
    }
  };

  // Build unique channel types (pair+s+l+radType combinations)
  struct ChannelType {
    int pairIndex;
    double sValue;
    int lValue;
    QChar radType;
    bool operator==(const ChannelType &other) const {
      return pairIndex == other.pairIndex && sValue == other.sValue &&
          lValue == other.lValue && radType == other.radType;
    }
  };
  QList<ChannelType> channelTypes;

  for (const ChannelsData &ch : channels) {
    ChannelType ct = {ch.pairIndex, ch.sValue, ch.lValue, ch.radType};
    bool found = false;
    for (const ChannelType &existing : channelTypes) {
      if (existing == ct) {
        found = true;
        break;
      }
    }
    if (!found) {
      channelTypes.append(ct);
    }
  }

  // Sort levels by energy - create index mapping
  QList<int> sortedIndices;
  for (int i = 0; i < levels.size(); i++) {
    sortedIndices.append(i);
  }
  std::sort(sortedIndices.begin(), sortedIndices.end(), [&levels](int a, int b) {
    return levels.at(a).energy < levels.at(b).energy;
  });

  // Generate LaTeX
  QString latex;
  QTextStream out(&latex);

  // Build column specification
  QString colSpec = "c|c|";
  for (int i = 0; i < channelTypes.size(); i++) {
    colSpec += "c";
  }

  out << "\\begin{table}[htbp]\n";
  out << "\\centering\n";
  out << "\\caption{R-Matrix levels and reduced width amplitudes}\n";
  out << "\\begin{tabular}{" << colSpec << "}\n";
  out << "\\hline\\hline\n";

  // First header row: J^pi, E (MeV), then channel names
  out << "$J^\\pi$ & $E$ (MeV)";
  for (const ChannelType &ct : channelTypes) {
    QString channelLabel;
    if (ct.radType == 'P') {
      // Particle channel: show light+heavy particle
      const PairsData &pair = pairs.at(ct.pairIndex);
      QString lightLabel = getLatexParticleLabel(pair, 0);
      QString heavyLabel = getLatexParticleLabel(pair, 1);
      channelLabel = QString("%1+%2").arg(lightLabel).arg(heavyLabel);
    } else if (ct.radType == 'E' || ct.radType == 'M') {
      // Gamma channel: show multipolarity and final state
      const PairsData &pair = pairs.at(ct.pairIndex);
      QString heavyLabel = getLatexParticleLabel(pair, 1);
      channelLabel = QString("$%1%2$ $\\to$ %3").arg(ct.radType).arg(ct.lValue).arg(heavyLabel);
    } else if (ct.radType == 'G') {
      channelLabel = "GT";
    } else if (ct.radType == 'F') {
      channelLabel = "Fermi";
    } else {
      channelLabel = QString("%1").arg(ct.radType);
    }
    out << " & " << channelLabel;
  }
  out << " \\\\\n";

  // Second header row: empty for J^pi and E, then s,l or excitation energy
  out << " & ";
  for (const ChannelType &ct : channelTypes) {
    QString subLabel;
    if (ct.radType == 'P') {
      // Particle channel: show s, l
      subLabel = QString("$s=%1$, $\\ell=%2$").arg(formatSpinLabel(ct.sValue)).arg(ct.lValue);
    } else if (ct.radType == 'E' || ct.radType == 'M') {
      // Gamma channel: show excitation energy (final state)
      const PairsData &pair = pairs.at(ct.pairIndex);
      subLabel = QString("$E_x=%1$ MeV").arg(pair.excitationEnergy, 0, 'f', 3);
    } else if (ct.radType == 'G' || ct.radType == 'F') {
      const PairsData &pair = pairs.at(ct.pairIndex);
      subLabel = QString("$E_x=%1$ MeV").arg(pair.excitationEnergy, 0, 'f', 3);
    } else {
      subLabel = "";
    }
    out << " & " << subLabel;
  }
  out << " \\\\\n";
  out << "\\hline\n";

  // Data rows: one per level, sorted by energy
  for (int sortedIdx = 0; sortedIdx < sortedIndices.size(); sortedIdx++) {
    int li = sortedIndices.at(sortedIdx);
    const LevelsData &level = levels.at(li);

    // Format J^pi
    QString jLabel = formatSpinLabel(level.jValue);
    QString piLabel = (level.piValue == -1) ? "-" : "+";

    out << "$" << jLabel << "^{" << piLabel << "}$ & "
        << QString::number(level.energy, 'f', 4);

    // Find width for each channel type
    for (const ChannelType &ct : channelTypes) {
      QString widthStr = "---";
      for (const ChannelsData &ch : channels) {
        if (ch.levelIndex == li && ch.pairIndex == ct.pairIndex &&
            ch.sValue == ct.sValue && ch.lValue == ct.lValue &&
            ch.radType == ct.radType) {
          widthStr = QString::number(ch.reducedWidth, 'g', 5);
          break;
        }
      }
      out << " & " << widthStr;
    }
    out << " \\\\\n";
  }

  out << "\\hline\\hline\n";
  out << "\\end{tabular}\n";
  out << "\\end{table}\n";

  // Ask user what to do with the LaTeX
  QMessageBox msgBox(this);
  msgBox.setWindowTitle(tr("Export LaTeX Table"));
  msgBox.setText(tr("LaTeX table generated successfully."));
  msgBox.setInformativeText(tr("Choose an action:"));
  QPushButton *copyButton = msgBox.addButton(tr("Copy to Clipboard"), QMessageBox::ActionRole);
  QPushButton *saveButton = msgBox.addButton(tr("Save to File"), QMessageBox::ActionRole);
  msgBox.addButton(QMessageBox::Cancel);
  msgBox.exec();

  if (msgBox.clickedButton() == copyButton) {
    QApplication::clipboard()->setText(latex);
    QMessageBox::information(this, tr("Copied"),
                             tr("LaTeX table copied to clipboard."));
  } else if (msgBox.clickedButton() == saveButton) {
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save LaTeX Table"),
                                                    "levels_table.tex",
                                                    tr("LaTeX files (*.tex);;All files (*)"));
    if (!fileName.isEmpty()) {
      QFile file(fileName);
      if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream fileOut(&file);
        fileOut << latex;
        file.close();
        QMessageBox::information(this, tr("Saved"),
                                 tr("LaTeX table saved to %1").arg(fileName));
      } else {
        QMessageBox::warning(this, tr("Error"),
                             tr("Could not save file: %1").arg(file.errorString()));
      }
    }
  }
}
