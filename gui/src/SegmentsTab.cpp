#include <QGridLayout>
#include <QHBoxLayout>
#include <QMessageBox>
#include <QGroupBox>
#include <QHeaderView>
#include <QTextStream>
#include <QLabel>
#include <QSet>
#include <QTextDocument>
#include <algorithm>

#include "SegmentsTab.h"
#include "RichTextDelegate.h"
#include "InfoDialog.h"
#include "PairsModel.h"
#include "ExforDialog.h"

SegmentsTab::SegmentsTab(QWidget *parent) :
  QWidget(parent) {
  pairsModel = nullptr;
  segmentsDataModel = new SegmentsDataModel;

  // Connect model signals to update filters when data changes
  connect(segmentsDataModel, SIGNAL(dataChanged(QModelIndex, QModelIndex)), this, SLOT(updateFilterComboboxes()));
  connect(segmentsDataModel, SIGNAL(rowsInserted(QModelIndex, int, int)), this, SLOT(updateFilterComboboxes()));
  connect(segmentsDataModel, SIGNAL(rowsRemoved(QModelIndex, int, int)), this, SLOT(updateFilterComboboxes()));

  segmentsDataView = new QTableView;
  segmentsDataView->setModel(segmentsDataModel);
  RichTextDelegate *rt = new RichTextDelegate();
  segmentsDataView->setItemDelegateForColumn(1, rt);
  segmentsDataView->setColumnHidden(2, true);
  segmentsDataView->setColumnHidden(4, true);
  segmentsDataView->setColumnHidden(6, true);
  segmentsDataView->horizontalHeader()->setStretchLastSection(true);
  segmentsDataView->horizontalHeader()->setHighlightSections(false);
  segmentsDataView->setSelectionBehavior(QAbstractItemView::SelectRows);
  segmentsDataView->setSelectionMode(QAbstractItemView::SingleSelection);
  segmentsDataView->setEditTriggers(QAbstractItemView::NoEditTriggers);
  segmentsDataView->setShowGrid(false);
  segmentsDataView->setColumnWidth(0, 27);
  segmentsDataView->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Fixed);
  segmentsDataView->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Fixed);
  segmentsDataView->horizontalHeader()->setSectionResizeMode(2, QHeaderView::Fixed);
  for (int i = 1; i < 3; i++) segmentsDataView->setColumnWidth(i, 200);
  for (int i = 3; i < 7; i++) segmentsDataView->setColumnWidth(i, 120);
  segmentsDataView->setColumnWidth(7, 140);
  segmentsDataView->setItemDelegateForColumn(7, rt);
  segmentsDataView->setColumnWidth(8, 220);
  segmentsDataView->setItemDelegateForColumn(9, rt);
  segmentsDataView->setColumnHidden(10, true);
  segmentsDataView->setColumnHidden(11, true);
  segmentsDataView->setColumnHidden(12, true);
  segmentsDataView->setColumnHidden(13, true);
  segmentsDataView->setColumnWidth(14, 120);
  segmentsDataView->setItemDelegateForColumn(14, rt);
  segmentsDataView->setColumnHidden(15, true);
  segmentsDataView->setColumnHidden(16, true);
  connect(segmentsDataView->selectionModel(), SIGNAL(selectionChanged(QItemSelection, QItemSelection)), this, SLOT(updateSegDataButtons(QItemSelection)));
  connect(segmentsDataView, SIGNAL(doubleClicked(QModelIndex)), this, SLOT(editSegDataLine()));

  segmentsTestModel = new SegmentsTestModel;
  segmentsTestView = new QTableView;
  segmentsTestView->setModel(segmentsTestModel);
  segmentsTestView->setItemDelegateForColumn(1, rt);
  segmentsTestView->setColumnHidden(2, true);
  segmentsTestView->setColumnHidden(4, true);
  segmentsTestView->setColumnHidden(7, true);
  segmentsTestView->horizontalHeader()->setStretchLastSection(true);
  segmentsTestView->horizontalHeader()->setHighlightSections(false);
  segmentsTestView->setSelectionBehavior(QAbstractItemView::SelectRows);
  segmentsTestView->setSelectionMode(QAbstractItemView::SingleSelection);
  segmentsTestView->setEditTriggers(QAbstractItemView::NoEditTriggers);
  segmentsTestView->setShowGrid(false);
  segmentsTestView->setColumnWidth(0, 27);
  segmentsTestView->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Fixed);
  segmentsTestView->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Fixed);
  segmentsTestView->horizontalHeader()->setSectionResizeMode(2, QHeaderView::Fixed);
  for (int i = 1; i < 3; i++) segmentsTestView->setColumnWidth(i, 200);
  segmentsTestView->setColumnWidth(3, 120);
  segmentsTestView->setColumnWidth(5, 90);
  segmentsTestView->setColumnWidth(6, 120);
  for (int i = 8; i < SegmentsTestData::SIZE; i++) segmentsTestView->setColumnWidth(i, 90);
  segmentsTestView->setItemDelegateForColumn(9, rt);
  segmentsTestView->setColumnHidden(10, true);
  segmentsTestView->setColumnHidden(11, true);
  segmentsTestView->setColumnHidden(12, true);
  connect(segmentsTestView->selectionModel(), SIGNAL(selectionChanged(QItemSelection, QItemSelection)), this, SLOT(updateSegTestButtons(QItemSelection)));
  connect(segmentsTestView, SIGNAL(doubleClicked(QModelIndex)), this, SLOT(editSegTestLine()));

  // segDataAddButton=new QPushButton(tr("Add Line"));
  // segDataDeleteButton = new QPushButton(tr("Delete Line"));
  segDataAddButton = new QPushButton(tr("+"));
  segDataAddButton->setMaximumSize(28, 28);
  segDataDeleteButton = new QPushButton(tr("-"));
  segDataDeleteButton->setEnabled(false);
  segDataDeleteButton->setMaximumSize(28, 28);
  segDataUpButton = new QPushButton;
  segDataUpButton->setIcon(style()->standardIcon(QStyle::SP_ArrowUp));
  segDataUpButton->setEnabled(false);
  segDataUpButton->setMaximumSize(28, 28);
  segDataDownButton = new QPushButton;
  segDataDownButton->setIcon(style()->standardIcon(QStyle::SP_ArrowDown));
  segDataDownButton->setEnabled(false);
  segDataDownButton->setMaximumSize(28, 28);
  connect(segDataAddButton, SIGNAL(clicked()), this, SLOT(addSegDataLine()));
  connect(segDataDeleteButton, SIGNAL(clicked()), this, SLOT(deleteSegDataLine()));
  connect(segDataUpButton, SIGNAL(clicked()), this, SLOT(moveSegDataLineUp()));
  connect(segDataDownButton, SIGNAL(clicked()), this, SLOT(moveSegDataLineDown()));

  segDataCheckAllButton = new QPushButton(tr("Check All"));
  segDataCheckAllButton->setMaximumHeight(28);
  connect(segDataCheckAllButton, SIGNAL(clicked()), this, SLOT(checkAllSegData()));
  segDataUncheckAllButton = new QPushButton(tr("Uncheck All"));
  segDataUncheckAllButton->setMaximumHeight(28);
  connect(segDataUncheckAllButton, SIGNAL(clicked()), this, SLOT(uncheckAllSegData()));

  segDataExforButton = new QPushButton(tr("EXFOR"));
  segDataExforButton->setMaximumHeight(28);
  segDataExforButton->setToolTip(tr("Retrieve experimental data from the EXFOR database"));
  connect(segDataExforButton, SIGNAL(clicked()), this, SLOT(getDataFromExfor()));

  // Add filter comboboxes
  segDataEntranceFilter = new QComboBox;
  segDataEntranceFilter->addItem(tr("All Entrance Pairs"));
  segDataEntranceFilter->setMaximumHeight(28);
  connect(segDataEntranceFilter, SIGNAL(currentIndexChanged(int)), this, SLOT(filterSegDataByPairs()));

  segDataExitFilter = new QComboBox;
  segDataExitFilter->addItem(tr("All Exit Pairs"));
  segDataExitFilter->setMaximumHeight(28);
  connect(segDataExitFilter, SIGNAL(currentIndexChanged(int)), this, SLOT(filterSegDataByPairs()));

  segTestAddButton = new QPushButton(tr("+"));
  segTestAddButton->setMaximumSize(28, 28);
  segTestDeleteButton = new QPushButton(tr("-"));
  segTestDeleteButton->setMaximumSize(28, 28);
  segTestDeleteButton->setEnabled(false);
  segTestUpButton = new QPushButton;
  segTestUpButton->setIcon(style()->standardIcon(QStyle::SP_ArrowUp));
  segTestUpButton->setEnabled(false);
  segTestUpButton->setMaximumSize(28, 28);
  segTestDownButton = new QPushButton;
  segTestDownButton->setIcon(style()->standardIcon(QStyle::SP_ArrowDown));
  segTestDownButton->setEnabled(false);
  segTestDownButton->setMaximumSize(28, 28);
  connect(segTestAddButton, SIGNAL(clicked()), this, SLOT(addSegTestLine()));
  connect(segTestDeleteButton, SIGNAL(clicked()), this, SLOT(deleteSegTestLine()));
  connect(segTestUpButton, SIGNAL(clicked()), this, SLOT(moveSegTestLineUp()));
  connect(segTestDownButton, SIGNAL(clicked()), this, SLOT(moveSegTestLineDown()));

  QGroupBox *segDataBox = new QGroupBox(tr("Segments From Data"));
  QGridLayout *segDataLayout = new QGridLayout;
  segDataLayout->addWidget(segmentsDataView, 0, 0);

  // Button row: [+] [-]   [Check All] [Uncheck All]   [up] [down]   <stretch>   Filter: [Entrance] [Exit]   [EXFOR]
  // A horizontal box vertically centers every widget so the buttons, filter
  // comboboxes and label all line up.  On macOS native push-button bezels are
  // drawn slightly outside the widget geometry, so adjacent buttons overlap
  // unless the gap between them is widened (matches setHorizontalSpacing(11)
  // used by the other tabs' button rows).
#ifdef MACX_SPACING
  const int kBtnGap = 11;    // gap between neighbouring buttons
  const int kGroupGap = 20;  // gap between logical button groups
#else
  const int kBtnGap = 6;
  const int kGroupGap = 14;
#endif
  QHBoxLayout *segDataButtonBox = new QHBoxLayout;
  segDataButtonBox->setSpacing(0);  // all gaps set explicitly below
  segDataButtonBox->addWidget(segDataAddButton);
  segDataButtonBox->addSpacing(kBtnGap);
  segDataButtonBox->addWidget(segDataDeleteButton);
  segDataButtonBox->addSpacing(kGroupGap);
  segDataButtonBox->addWidget(segDataCheckAllButton);
  segDataButtonBox->addSpacing(kBtnGap);
  segDataButtonBox->addWidget(segDataUncheckAllButton);
  segDataButtonBox->addSpacing(kGroupGap);
  segDataButtonBox->addWidget(segDataUpButton);
  segDataButtonBox->addSpacing(kBtnGap);
  segDataButtonBox->addWidget(segDataDownButton);
  segDataButtonBox->addStretch(1);
  segDataButtonBox->addWidget(new QLabel(tr("Filter:")));
  segDataButtonBox->addSpacing(kBtnGap);
  segDataButtonBox->addWidget(segDataEntranceFilter);
  segDataButtonBox->addSpacing(kBtnGap);
  segDataButtonBox->addWidget(segDataExitFilter);
  segDataButtonBox->addSpacing(kGroupGap);
  segDataButtonBox->addWidget(segDataExforButton);
  segDataLayout->addLayout(segDataButtonBox, 1, 0);

  segDataBox->setLayout(segDataLayout);

  QGroupBox *segTestBox = new QGroupBox(tr("Segments Without Data"));
  QGridLayout *segTestLayout = new QGridLayout;
  segTestLayout->addWidget(segmentsTestView, 0, 0);
  QGridLayout *segTestButtonBox = new QGridLayout;
  segTestButtonBox->addWidget(segTestAddButton, 0, 0);
  segTestButtonBox->addWidget(segTestDeleteButton, 0, 1);
  segTestButtonBox->addItem(new QSpacerItem(28, 28), 0, 2);
  segTestButtonBox->addWidget(segTestUpButton, 0, 3);
  segTestButtonBox->addWidget(segTestDownButton, 0, 4);
  segTestButtonBox->setColumnStretch(0, 0);
  segTestButtonBox->setColumnStretch(1, 0);
  segTestButtonBox->setColumnStretch(2, 1);
  segTestButtonBox->setColumnStretch(3, 0);
  segTestButtonBox->setColumnStretch(4, 0);
#ifdef MACX_SPACING
  segTestButtonBox->setHorizontalSpacing(11);
#else
  segTestButtonBox->setHorizontalSpacing(0);
#endif
  segTestLayout->addLayout(segTestButtonBox, 1, 0);
  segTestBox->setLayout(segTestLayout);

  QGridLayout *mainLayout = new QGridLayout;
  mainLayout->addWidget(segDataBox, 0, 0);
  mainLayout->addWidget(segTestBox, 1, 0);

  setLayout(mainLayout);
}

SegmentsTestModel *SegmentsTab::getSegmentsTestModel() {
  return segmentsTestModel;
}

SegmentsDataModel *SegmentsTab::getSegmentsDataModel() {
  return segmentsDataModel;
}

void SegmentsTab::deleteSegDataLine() {
  QItemSelectionModel *selectionModel = segmentsDataView->selectionModel();
  QModelIndexList indexes = selectionModel->selectedRows();
  QModelIndex index = indexes.at(0);

  segmentsDataModel->removeRows(index.row(), 1, QModelIndex());
  updateSegDataButtons(selectionModel->selection());
}

void SegmentsTab::deleteSegTestLine() {
  QItemSelectionModel *selectionModel = segmentsTestView->selectionModel();
  QModelIndexList indexes = selectionModel->selectedRows();
  QModelIndex index = indexes.at(0);

  segmentsTestModel->removeRows(index.row(), 1, QModelIndex());
  updateSegTestButtons(selectionModel->selection());
}

void SegmentsTab::addSegDataLine() {
  AddSegDataDialog aDialog;
  if (aDialog.exec()) {
    SegmentsDataData newLine;
    newLine.isActive = 1;
    newLine.entrancePairIndex = aDialog.entrancePairIndexSpin->value();
    if (aDialog.dataTypeCode() == 3)
      newLine.exitPairIndex = -1;
    else
      newLine.exitPairIndex = aDialog.exitPairIndexSpin->value();
    newLine.lowEnergy = aDialog.lowEnergyText->text().toDouble();
    newLine.highEnergy = aDialog.highEnergyText->text().toDouble();
    newLine.lowAngle = aDialog.lowAngleText->text().toDouble();
    newLine.highAngle = aDialog.highAngleText->text().toDouble();
    newLine.dataType = aDialog.dataTypeCode();
    newLine.dataFile = aDialog.dataFileText->text();
    newLine.dataNorm = aDialog.dataNormText->text().toDouble();
    newLine.dataNormError = aDialog.dataNormErrorText->text().toDouble();
    if (aDialog.varyNormCheck->isChecked())
      newLine.varyNorm = 1;
    else
      newLine.varyNorm = 0;
    newLine.phaseJ = aDialog.phaseJValueText->text().toDouble();
    newLine.phaseL = aDialog.phaseLValueText->text().toInt();
    newLine.energyShift = aDialog.energyShiftText->text().toDouble();
    newLine.energyShiftError = aDialog.energyShiftErrorText->text().toDouble();
    if (aDialog.varyEnergyShiftCheck->isChecked())
      newLine.varyEnergyShift = 1;
    else
      newLine.varyEnergyShift = 0;

    if (aDialog.advancedModeCheck->isChecked()) {
      newLine.isAdvanced = 1;
      newLine.operationType = aDialog.operationCombo->currentIndex();
      QStringList components;
      for (int i = 0; i < aDialog.componentsList->count(); i++) {
        components.append(aDialog.componentsList->item(i)->text());
      }
      newLine.componentsList = components.join(";");
    } else {
      newLine.isAdvanced = 0;
      newLine.operationType = 0;
      newLine.componentsList = "";
    }
    if (aDialog.uposCheck->isChecked()) {
      newLine.isUPOS = 1;
      newLine.secondaryDecayL = aDialog.secondaryLSpin->value();
      newLine.finalJ = aDialog.finalJText->text().toDouble();
      newLine.delta = aDialog.deltaText->text().toDouble();
    } else {
      newLine.isUPOS = 0;
      newLine.secondaryDecayL = 0;
      newLine.finalJ = 0.0;
      newLine.delta = 0.0;
    }
    addSegDataLine(newLine);
  }
}

void SegmentsTab::addSegDataLine(SegmentsDataData line) {
  QList<SegmentsDataData> lines = segmentsDataModel->getLines();
  if (segmentsDataModel->isSegDataLine(line) == -1) {
    segmentsDataModel->insertRows(lines.size(), 1, QModelIndex());
    QModelIndex index = segmentsDataModel->index(lines.size(), 0, QModelIndex());
    segmentsDataModel->setData(index, line.isActive, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 1, QModelIndex());
    segmentsDataModel->setData(index, line.entrancePairIndex, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 2, QModelIndex());
    segmentsDataModel->setData(index, line.exitPairIndex, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 3, QModelIndex());
    segmentsDataModel->setData(index, line.lowEnergy, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 4, QModelIndex());
    segmentsDataModel->setData(index, line.highEnergy, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 5, QModelIndex());
    segmentsDataModel->setData(index, line.lowAngle, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 6, QModelIndex());
    segmentsDataModel->setData(index, line.highAngle, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 7, QModelIndex());
    segmentsDataModel->setData(index, line.dataType, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 8, QModelIndex());
    segmentsDataModel->setData(index, line.dataFile, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 9, QModelIndex());
    segmentsDataModel->setData(index, line.dataNorm, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 10, QModelIndex());
    segmentsDataModel->setData(index, line.dataNormError, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 11, QModelIndex());
    segmentsDataModel->setData(index, line.varyNorm, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 12, QModelIndex());
    segmentsDataModel->setData(index, line.phaseJ, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 13, QModelIndex());
    segmentsDataModel->setData(index, line.phaseL, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 14, QModelIndex());
    segmentsDataModel->setData(index, line.energyShift, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 15, QModelIndex());
    segmentsDataModel->setData(index, line.energyShiftError, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 16, QModelIndex());
    segmentsDataModel->setData(index, line.varyEnergyShift, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 17, QModelIndex());
    segmentsDataModel->setData(index, line.isAdvanced, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 18, QModelIndex());
    segmentsDataModel->setData(index, line.operationType, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 19, QModelIndex());
    segmentsDataModel->setData(index, line.componentsList, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 20, QModelIndex());
    segmentsDataModel->setData(index, line.isUPOS, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 21, QModelIndex());
    segmentsDataModel->setData(index, line.secondaryDecayL, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 22, QModelIndex());
    segmentsDataModel->setData(index, line.finalJ, Qt::EditRole);
    index = segmentsDataModel->index(lines.size(), 23, QModelIndex());
    segmentsDataModel->setData(index, line.delta, Qt::EditRole);
    segmentsDataView->resizeRowToContents(lines.size());
    updateSegDataButtons(segmentsDataView->selectionModel()->selection());
  } else {
    QMessageBox::information(this, tr("Duplicate Line"), tr("This line already exists."));
  }
}

void SegmentsTab::addSegTestLine() {
  AddSegTestDialog aDialog;
  if (aDialog.exec()) {
    SegmentsTestData newLine;
    newLine.isActive = 1;
    newLine.entrancePairIndex = aDialog.entrancePairIndexSpin->value();
    if (aDialog.dataTypeCode() == 4)
      newLine.exitPairIndex = -1;
    else
      newLine.exitPairIndex = aDialog.exitPairIndexSpin->value();
    newLine.lowEnergy = aDialog.lowEnergyText->text().toDouble();
    newLine.highEnergy = aDialog.highEnergyText->text().toDouble();
    newLine.energyStep = aDialog.energyStepText->text().toDouble();
    newLine.lowAngle = aDialog.lowAngleText->text().toDouble();
    newLine.highAngle = aDialog.highAngleText->text().toDouble();
    newLine.angleStep = aDialog.angleStepText->text().toDouble();
    newLine.dataType = aDialog.dataTypeCode();
    newLine.phaseJ = aDialog.phaseJValueText->text().toDouble();
    newLine.phaseL = aDialog.phaseLValueText->text().toInt();
    newLine.maxAngDistOrder = aDialog.angDistSpin->value();

    if (aDialog.advancedModeCheck->isChecked()) {
      newLine.isAdvanced = 1;
      newLine.operationType = aDialog.operationCombo->currentIndex();
      QStringList components;
      for (int i = 0; i < aDialog.componentsList->count(); i++) {
        components.append(aDialog.componentsList->item(i)->text());
      }
      newLine.componentsList = components.join(";");
    } else {
      newLine.isAdvanced = 0;
      newLine.operationType = 0;
      newLine.componentsList = "";
    }
    addSegTestLine(newLine);
  }
}

void SegmentsTab::addSegTestLine(SegmentsTestData line) {
  QList<SegmentsTestData> lines = segmentsTestModel->getLines();
  if (segmentsTestModel->isSegTestLine(line) == -1) {
    segmentsTestModel->insertRows(lines.size(), 1, QModelIndex());
    QModelIndex index = segmentsTestModel->index(lines.size(), 0, QModelIndex());
    segmentsTestModel->setData(index, line.isActive, Qt::EditRole);
    index = segmentsTestModel->index(lines.size(), 1, QModelIndex());
    segmentsTestModel->setData(index, line.entrancePairIndex, Qt::EditRole);
    index = segmentsTestModel->index(lines.size(), 2, QModelIndex());
    segmentsTestModel->setData(index, line.exitPairIndex, Qt::EditRole);
    index = segmentsTestModel->index(lines.size(), 3, QModelIndex());
    segmentsTestModel->setData(index, line.lowEnergy, Qt::EditRole);
    index = segmentsTestModel->index(lines.size(), 4, QModelIndex());
    segmentsTestModel->setData(index, line.highEnergy, Qt::EditRole);
    index = segmentsTestModel->index(lines.size(), 5, QModelIndex());
    segmentsTestModel->setData(index, line.energyStep, Qt::EditRole);
    index = segmentsTestModel->index(lines.size(), 6, QModelIndex());
    segmentsTestModel->setData(index, line.lowAngle, Qt::EditRole);
    index = segmentsTestModel->index(lines.size(), 7, QModelIndex());
    segmentsTestModel->setData(index, line.highAngle, Qt::EditRole);
    index = segmentsTestModel->index(lines.size(), 8, QModelIndex());
    segmentsTestModel->setData(index, line.angleStep, Qt::EditRole);
    index = segmentsTestModel->index(lines.size(), 9, QModelIndex());
    segmentsTestModel->setData(index, line.dataType, Qt::EditRole);
    index = segmentsTestModel->index(lines.size(), 10, QModelIndex());
    segmentsTestModel->setData(index, line.phaseJ, Qt::EditRole);
    index = segmentsTestModel->index(lines.size(), 11, QModelIndex());
    segmentsTestModel->setData(index, line.phaseL, Qt::EditRole);
    index = segmentsTestModel->index(lines.size(), 12, QModelIndex());
    segmentsTestModel->setData(index, line.maxAngDistOrder, Qt::EditRole);
    index = segmentsTestModel->index(lines.size(), 13, QModelIndex());
    segmentsTestModel->setData(index, line.isAdvanced, Qt::EditRole);
    index = segmentsTestModel->index(lines.size(), 14, QModelIndex());
    segmentsTestModel->setData(index, line.operationType, Qt::EditRole);
    index = segmentsTestModel->index(lines.size(), 15, QModelIndex());
    segmentsTestModel->setData(index, line.componentsList, Qt::EditRole);
    segmentsTestView->resizeRowToContents(lines.size());
    updateSegTestButtons(segmentsTestView->selectionModel()->selection());
  } else {
    QMessageBox::information(this, tr("Duplicate Line"), tr("This line already exists."));
  }
}

void SegmentsTab::editSegDataLine() {
  QItemSelectionModel *selectionModel = segmentsDataView->selectionModel();
  QModelIndexList indexes = selectionModel->selectedRows();
  QModelIndex index = indexes[0];

  QModelIndex i = segmentsDataModel->index(index.row(), 1, QModelIndex());
  QVariant var = segmentsDataModel->data(i, Qt::EditRole);
  int entrancePairIndex = var.toInt();
  i = segmentsDataModel->index(index.row(), 2, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  int exitPairIndex = var.toInt();
  i = segmentsDataModel->index(index.row(), 3, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  QString lowEnergy = var.toString();
  i = segmentsDataModel->index(index.row(), 4, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  QString highEnergy = var.toString();
  i = segmentsDataModel->index(index.row(), 5, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  QString lowAngle = var.toString();
  i = segmentsDataModel->index(index.row(), 6, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  QString highAngle = var.toString();
  i = segmentsDataModel->index(index.row(), 7, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  int dataType = var.toInt();
  i = segmentsDataModel->index(index.row(), 8, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  QString dataFile = var.toString();
  i = segmentsDataModel->index(index.row(), 9, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  QString dataNorm = var.toString();
  i = segmentsDataModel->index(index.row(), 10, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  QString dataNormError = var.toString();
  i = segmentsDataModel->index(index.row(), 11, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  int varyNorm = var.toInt();
  i = segmentsDataModel->index(index.row(), 12, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  QString phaseJ = var.toString();
  i = segmentsDataModel->index(index.row(), 13, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  QString phaseL = var.toString();
  i = segmentsDataModel->index(index.row(), 14, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  QString energyShift = var.toString();
  i = segmentsDataModel->index(index.row(), 15, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  QString energyShiftError = var.toString();
  i = segmentsDataModel->index(index.row(), 16, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  int varyEnergyShift = var.toInt();
  i = segmentsDataModel->index(index.row(), 17, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  int isAdvanced = var.toInt();
  i = segmentsDataModel->index(index.row(), 18, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  int operationType = var.toInt();
  i = segmentsDataModel->index(index.row(), 19, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  QString componentsList = var.toString();
  i = segmentsDataModel->index(index.row(), 20, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  int isUPOS = var.toInt();
  i = segmentsDataModel->index(index.row(), 21, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  int secondaryDecayL = var.toInt();
  i = segmentsDataModel->index(index.row(), 22, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  QString finalJ = var.toString();
  i = segmentsDataModel->index(index.row(), 23, QModelIndex());
  var = segmentsDataModel->data(i, Qt::EditRole);
  QString delta = var.toString();

  AddSegDataDialog aDialog;
  aDialog.setWindowTitle(tr("Edit a Segment From Data"));
  aDialog.entrancePairIndexSpin->setValue(entrancePairIndex);
  aDialog.exitPairIndexSpin->setValue(exitPairIndex);
  aDialog.lowEnergyText->setText(lowEnergy);
  aDialog.highEnergyText->setText(highEnergy);
  aDialog.lowAngleText->setText(lowAngle);
  aDialog.highAngleText->setText(highAngle);
  aDialog.setDataTypeCode(dataType);
  aDialog.dataFileText->setText(dataFile);
  aDialog.dataNormText->setText(dataNorm);
  aDialog.dataNormErrorText->setText(dataNormError);
  if (varyNorm == 1)
    aDialog.varyNormCheck->setChecked(true);
  else
    aDialog.varyNormCheck->setChecked(false);
  aDialog.phaseJValueText->setText(phaseJ);
  aDialog.phaseLValueText->setText(phaseL);
  aDialog.energyShiftText->setText(energyShift);
  aDialog.energyShiftErrorText->setText(energyShiftError);
  if (varyEnergyShift == 1)
    aDialog.varyEnergyShiftCheck->setChecked(true);
  else
    aDialog.varyEnergyShiftCheck->setChecked(false);

  // Set advanced mode data
  if (isAdvanced == 1)
    aDialog.advancedModeCheck->setChecked(true);
  else
    aDialog.advancedModeCheck->setChecked(false);
  aDialog.operationCombo->setCurrentIndex(operationType);

  // Load components list
  if (!componentsList.isEmpty()) {
    QStringList components = componentsList.split(";", Qt::SkipEmptyParts);
    for (const QString &component : components) {
      aDialog.componentsList->addItem(component);
    }
  }

  // Set UPOS data
  if (isUPOS == 1)
    aDialog.uposCheck->setChecked(true);
  else
    aDialog.uposCheck->setChecked(false);
  aDialog.secondaryLSpin->setValue(secondaryDecayL);
  aDialog.finalJText->setText(finalJ);
  aDialog.deltaText->setText(delta);

  if (aDialog.exec()) {
    int newEntrancePairIndex = aDialog.entrancePairIndexSpin->value();
    if (newEntrancePairIndex != entrancePairIndex) {
      i = segmentsDataModel->index(index.row(), 1, QModelIndex());
      segmentsDataModel->setData(i, newEntrancePairIndex, Qt::EditRole);
    }
    int newExitPairIndex = (aDialog.dataTypeCode() == 3) ? -1 : aDialog.exitPairIndexSpin->value();
    if (newExitPairIndex != exitPairIndex) {
      i = segmentsDataModel->index(index.row(), 2, QModelIndex());
      segmentsDataModel->setData(i, newExitPairIndex, Qt::EditRole);
    }
    QString newLowEnergy = aDialog.lowEnergyText->text();
    if (newLowEnergy != lowEnergy) {
      i = segmentsDataModel->index(index.row(), 3, QModelIndex());
      segmentsDataModel->setData(i, newLowEnergy, Qt::EditRole);
    }
    QString newHighEnergy = aDialog.highEnergyText->text();
    if (newHighEnergy != highEnergy) {
      i = segmentsDataModel->index(index.row(), 4, QModelIndex());
      segmentsDataModel->setData(i, newHighEnergy, Qt::EditRole);
    }
    QString newLowAngle = aDialog.lowAngleText->text();
    if (newLowAngle != lowAngle) {
      i = segmentsDataModel->index(index.row(), 5, QModelIndex());
      segmentsDataModel->setData(i, newLowAngle, Qt::EditRole);
    }
    QString newHighAngle = aDialog.highAngleText->text();
    if (newHighAngle != highAngle) {
      i = segmentsDataModel->index(index.row(), 6, QModelIndex());
      segmentsDataModel->setData(i, newHighAngle, Qt::EditRole);
    }
    int newDataType = aDialog.dataTypeCode();
    if (newDataType != dataType) {
      i = segmentsDataModel->index(index.row(), 7, QModelIndex());
      segmentsDataModel->setData(i, newDataType, Qt::EditRole);
    }
    QString newDataFile = aDialog.dataFileText->text();
    if (newDataFile != dataFile) {
      i = segmentsDataModel->index(index.row(), 8, QModelIndex());
      segmentsDataModel->setData(i, newDataFile, Qt::EditRole);
    }
    QString newDataNorm = aDialog.dataNormText->text();
    if (newDataNorm != dataNorm) {
      i = segmentsDataModel->index(index.row(), 9, QModelIndex());
      segmentsDataModel->setData(i, newDataNorm, Qt::EditRole);
    }
    QString newDataNormError = aDialog.dataNormErrorText->text();
    if (newDataNormError != dataNormError) {
      i = segmentsDataModel->index(index.row(), 10, QModelIndex());
      segmentsDataModel->setData(i, newDataNormError, Qt::EditRole);
    }
    int newVaryNorm = 0;
    if (aDialog.varyNormCheck->isChecked()) newVaryNorm = 1;
    if (newVaryNorm != varyNorm) {
      i = segmentsDataModel->index(index.row(), 11, QModelIndex());
      segmentsDataModel->setData(i, newVaryNorm, Qt::EditRole);
    }
    QString newPhaseJ = aDialog.phaseJValueText->text();
    if (newPhaseJ != phaseJ) {
      i = segmentsDataModel->index(index.row(), 12, QModelIndex());
      segmentsDataModel->setData(i, newPhaseJ, Qt::EditRole);
    }
    QString newPhaseL = aDialog.phaseLValueText->text();
    if (newPhaseL != phaseL) {
      i = segmentsDataModel->index(index.row(), 13, QModelIndex());
      segmentsDataModel->setData(i, newPhaseL, Qt::EditRole);
    }
    QString newEnergyShift = aDialog.energyShiftText->text();
    if (newEnergyShift != energyShift) {
      i = segmentsDataModel->index(index.row(), 14, QModelIndex());
      segmentsDataModel->setData(i, newEnergyShift, Qt::EditRole);
    }
    QString newEnergyShiftError = aDialog.energyShiftErrorText->text();
    if (newEnergyShiftError != energyShiftError) {
      i = segmentsDataModel->index(index.row(), 15, QModelIndex());
      segmentsDataModel->setData(i, newEnergyShiftError, Qt::EditRole);
    }
    int newVaryEnergyShift = 0;
    if (aDialog.varyEnergyShiftCheck->isChecked()) newVaryEnergyShift = 1;
    if (newVaryEnergyShift != varyEnergyShift) {
      i = segmentsDataModel->index(index.row(), 16, QModelIndex());
      segmentsDataModel->setData(i, newVaryEnergyShift, Qt::EditRole);
    }

    // Save advanced segment data
    int newIsAdvanced = 0;
    if (aDialog.advancedModeCheck->isChecked()) newIsAdvanced = 1;
    if (newIsAdvanced != isAdvanced) {
      i = segmentsDataModel->index(index.row(), 17, QModelIndex());
      segmentsDataModel->setData(i, newIsAdvanced, Qt::EditRole);
    }

    int newOperationType = aDialog.operationCombo->currentIndex();
    if (newOperationType != operationType) {
      i = segmentsDataModel->index(index.row(), 18, QModelIndex());
      segmentsDataModel->setData(i, newOperationType, Qt::EditRole);
    }

    QStringList newComponents;
    for (int j = 0; j < aDialog.componentsList->count(); j++) {
      newComponents.append(aDialog.componentsList->item(j)->text());
    }
    QString newComponentsList = newComponents.join(";");
    if (newComponentsList != componentsList) {
      i = segmentsDataModel->index(index.row(), 19, QModelIndex());
      segmentsDataModel->setData(i, newComponentsList, Qt::EditRole);
    }

    int newIsUPOS = 0;
    if (aDialog.uposCheck->isChecked()) newIsUPOS = 1;
    if (newIsUPOS != isUPOS) {
      i = segmentsDataModel->index(index.row(), 20, QModelIndex());
      segmentsDataModel->setData(i, newIsUPOS, Qt::EditRole);
    }
    int newSecondaryDecayL = aDialog.secondaryLSpin->value();
    if (newSecondaryDecayL != secondaryDecayL) {
      i = segmentsDataModel->index(index.row(), 21, QModelIndex());
      segmentsDataModel->setData(i, newSecondaryDecayL, Qt::EditRole);
    }
    QString newFinalJ = aDialog.finalJText->text();
    if (newFinalJ != finalJ) {
      i = segmentsDataModel->index(index.row(), 22, QModelIndex());
      segmentsDataModel->setData(i, newFinalJ, Qt::EditRole);
    }
    QString newDelta = aDialog.deltaText->text();
    if (newDelta != delta) {
      i = segmentsDataModel->index(index.row(), 23, QModelIndex());
      segmentsDataModel->setData(i, newDelta, Qt::EditRole);
    }
  }
}

void SegmentsTab::editSegTestLine() {
  QItemSelectionModel *selectionModel = segmentsTestView->selectionModel();
  QModelIndexList indexes = selectionModel->selectedRows();
  QModelIndex index = indexes[0];

  QModelIndex i = segmentsTestModel->index(index.row(), 1, QModelIndex());
  QVariant var = segmentsTestModel->data(i, Qt::EditRole);
  int entrancePairIndex = var.toInt();
  i = segmentsTestModel->index(index.row(), 2, QModelIndex());
  var = segmentsTestModel->data(i, Qt::EditRole);
  int exitPairIndex = var.toInt();
  i = segmentsTestModel->index(index.row(), 3, QModelIndex());
  var = segmentsTestModel->data(i, Qt::EditRole);
  QString lowEnergy = var.toString();
  i = segmentsTestModel->index(index.row(), 4, QModelIndex());
  var = segmentsTestModel->data(i, Qt::EditRole);
  QString highEnergy = var.toString();
  i = segmentsTestModel->index(index.row(), 5, QModelIndex());
  var = segmentsTestModel->data(i, Qt::EditRole);
  QString energyStep = var.toString();
  i = segmentsTestModel->index(index.row(), 6, QModelIndex());
  var = segmentsTestModel->data(i, Qt::EditRole);
  QString lowAngle = var.toString();
  i = segmentsTestModel->index(index.row(), 7, QModelIndex());
  var = segmentsTestModel->data(i, Qt::EditRole);
  QString highAngle = var.toString();
  i = segmentsTestModel->index(index.row(), 8, QModelIndex());
  var = segmentsTestModel->data(i, Qt::EditRole);
  QString angleStep = var.toString();
  i = segmentsTestModel->index(index.row(), 9, QModelIndex());
  var = segmentsTestModel->data(i, Qt::EditRole);
  int dataType = var.toInt();
  i = segmentsTestModel->index(index.row(), 10, QModelIndex());
  var = segmentsTestModel->data(i, Qt::EditRole);
  QString phaseJ = var.toString();
  i = segmentsTestModel->index(index.row(), 11, QModelIndex());
  var = segmentsTestModel->data(i, Qt::EditRole);
  QString phaseL = var.toString();
  i = segmentsTestModel->index(index.row(), 12, QModelIndex());
  var = segmentsTestModel->data(i, Qt::EditRole);
  int maxAngDistOrder = var.toInt();
  i = segmentsTestModel->index(index.row(), 13, QModelIndex());
  var = segmentsTestModel->data(i, Qt::EditRole);
  int isAdvanced = var.toInt();
  i = segmentsTestModel->index(index.row(), 14, QModelIndex());
  var = segmentsTestModel->data(i, Qt::EditRole);
  int operationType = var.toInt();
  i = segmentsTestModel->index(index.row(), 15, QModelIndex());
  var = segmentsTestModel->data(i, Qt::EditRole);
  QString componentsList = var.toString();

  AddSegTestDialog aDialog;
  aDialog.setWindowTitle(tr("Edit a Segment Without Data"));
  aDialog.entrancePairIndexSpin->setValue(entrancePairIndex);
  aDialog.exitPairIndexSpin->setValue(exitPairIndex);
  aDialog.lowEnergyText->setText(lowEnergy);
  aDialog.highEnergyText->setText(highEnergy);
  aDialog.energyStepText->setText(energyStep);
  aDialog.lowAngleText->setText(lowAngle);
  aDialog.highAngleText->setText(highAngle);
  aDialog.angleStepText->setText(angleStep);
  aDialog.setDataTypeCode(dataType);
  aDialog.phaseJValueText->setText(phaseJ);
  aDialog.phaseLValueText->setText(phaseL);
  aDialog.angDistSpin->setValue(maxAngDistOrder);

  // Set advanced mode data
  if (isAdvanced == 1)
    aDialog.advancedModeCheck->setChecked(true);
  else
    aDialog.advancedModeCheck->setChecked(false);
  aDialog.operationCombo->setCurrentIndex(operationType);

  // Load components list
  if (!componentsList.isEmpty()) {
    QStringList components = componentsList.split(";", Qt::SkipEmptyParts);
    for (const QString &component : components) {
      aDialog.componentsList->addItem(component);
    }
  }

  if (aDialog.exec()) {
    int newEntrancePairIndex = aDialog.entrancePairIndexSpin->value();
    if (newEntrancePairIndex != entrancePairIndex) {
      i = segmentsTestModel->index(index.row(), 1, QModelIndex());
      segmentsTestModel->setData(i, newEntrancePairIndex, Qt::EditRole);
    }
    int newExitPairIndex = (aDialog.dataTypeCode() == 4) ? -1 : aDialog.exitPairIndexSpin->value();
    if (newExitPairIndex != exitPairIndex) {
      i = segmentsTestModel->index(index.row(), 2, QModelIndex());
      segmentsTestModel->setData(i, newExitPairIndex, Qt::EditRole);
    }
    QString newLowEnergy = aDialog.lowEnergyText->text();
    if (newLowEnergy != lowEnergy) {
      i = segmentsTestModel->index(index.row(), 3, QModelIndex());
      segmentsTestModel->setData(i, newLowEnergy, Qt::EditRole);
    }
    QString newHighEnergy = aDialog.highEnergyText->text();
    if (newHighEnergy != highEnergy) {
      i = segmentsTestModel->index(index.row(), 4, QModelIndex());
      segmentsTestModel->setData(i, newHighEnergy, Qt::EditRole);
    }
    QString newEnergyStep = aDialog.energyStepText->text();
    if (newEnergyStep != energyStep) {
      i = segmentsTestModel->index(index.row(), 5, QModelIndex());
      segmentsTestModel->setData(i, newEnergyStep, Qt::EditRole);
    }
    QString newLowAngle = aDialog.lowAngleText->text();
    if (newLowAngle != lowAngle) {
      i = segmentsTestModel->index(index.row(), 6, QModelIndex());
      segmentsTestModel->setData(i, newLowAngle, Qt::EditRole);
    }
    QString newHighAngle = aDialog.highAngleText->text();
    if (newHighAngle != highAngle) {
      i = segmentsTestModel->index(index.row(), 7, QModelIndex());
      segmentsTestModel->setData(i, newHighAngle, Qt::EditRole);
    }
    QString newAngleStep = aDialog.angleStepText->text();
    if (newAngleStep != angleStep) {
      i = segmentsTestModel->index(index.row(), 8, QModelIndex());
      segmentsTestModel->setData(i, newAngleStep, Qt::EditRole);
    }
    int newDataType = aDialog.dataTypeCode();
    if (newDataType != dataType) {
      i = segmentsTestModel->index(index.row(), 9, QModelIndex());
      segmentsTestModel->setData(i, newDataType, Qt::EditRole);
    }
    QString newPhaseJ = aDialog.phaseJValueText->text();
    if (newPhaseJ != phaseJ) {
      i = segmentsTestModel->index(index.row(), 10, QModelIndex());
      segmentsTestModel->setData(i, newPhaseJ, Qt::EditRole);
    }
    QString newPhaseL = aDialog.phaseLValueText->text();
    if (newPhaseL != phaseL) {
      i = segmentsTestModel->index(index.row(), 11, QModelIndex());
      segmentsTestModel->setData(i, newPhaseL, Qt::EditRole);
    }
    int newMaxAngDistOrder = aDialog.angDistSpin->value();
    if (newMaxAngDistOrder != maxAngDistOrder) {
      i = segmentsTestModel->index(index.row(), 12, QModelIndex());
      segmentsTestModel->setData(i, newMaxAngDistOrder, Qt::EditRole);
    }

    // Save advanced segment data
    int newIsAdvanced = 0;
    if (aDialog.advancedModeCheck->isChecked()) newIsAdvanced = 1;
    if (newIsAdvanced != isAdvanced) {
      i = segmentsTestModel->index(index.row(), 13, QModelIndex());
      segmentsTestModel->setData(i, newIsAdvanced, Qt::EditRole);
    }

    int newOperationType = aDialog.operationCombo->currentIndex();
    if (newOperationType != operationType) {
      i = segmentsTestModel->index(index.row(), 14, QModelIndex());
      segmentsTestModel->setData(i, newOperationType, Qt::EditRole);
    }

    QStringList newComponents;
    for (int j = 0; j < aDialog.componentsList->count(); j++) {
      newComponents.append(aDialog.componentsList->item(j)->text());
    }
    QString newComponentsList = newComponents.join(";");
    if (newComponentsList != componentsList) {
      i = segmentsTestModel->index(index.row(), 15, QModelIndex());
      segmentsTestModel->setData(i, newComponentsList, Qt::EditRole);
    }
  }
}

void SegmentsTab::moveSegDataLineUp() {
  moveSegDataLine(1);
}

void SegmentsTab::moveSegDataLineDown() {
  moveSegDataLine(0);
}

void SegmentsTab::moveSegDataLine(unsigned int upDown) {
  QItemSelectionModel *selectionModel = segmentsDataView->selectionModel();
  QModelIndexList selectionList = selectionModel->selectedRows();
  QModelIndex selectionIndex = selectionList.at(0);

  int previous = selectionIndex.row();
  int future;
  if (upDown == 0)
    future = previous + 1;
  else
    future = previous - 1;
  if (!segmentsDataModel->moveLine(previous, future)) return;
  segmentsDataView->resizeRowToContents(future);

  selectionModel->select(segmentsDataModel->index(future, 0, QModelIndex()),
                         QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
}

void SegmentsTab::moveSegTestLineUp() {
  moveSegTestLine(1);
}

void SegmentsTab::moveSegTestLineDown() {
  moveSegTestLine(0);
}

void SegmentsTab::moveSegTestLine(unsigned int upDown) {
  QItemSelectionModel *selectionModel = segmentsTestView->selectionModel();
  QModelIndexList selectionList = selectionModel->selectedRows();
  QModelIndex selectionIndex = selectionList.at(0);

  int previous = selectionIndex.row();
  int future;
  if (upDown == 0)
    future = previous + 1;
  else
    future = previous - 1;

  if (!segmentsTestModel->moveLine(previous, future)) return;
  segmentsTestView->resizeRowToContents(future);

  selectionModel->select(segmentsTestModel->index(future, 0, QModelIndex()),
                         QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
}

void SegmentsTab::updateSegDataButtons(const QItemSelection &selection) {
  QModelIndexList indexes = selection.indexes();

  if (indexes.isEmpty()) {
    segDataDeleteButton->setEnabled(false);
    segDataUpButton->setEnabled(false);
    segDataDownButton->setEnabled(false);
  } else {
    segDataDeleteButton->setEnabled(true);
    if (indexes.at(0).row() == 0)
      segDataUpButton->setEnabled(false);
    else
      segDataUpButton->setEnabled(true);
    if (indexes.at(0).row() == segmentsDataModel->rowCount(QModelIndex()) - 1)
      segDataDownButton->setEnabled(false);
    else
      segDataDownButton->setEnabled(true);
  }
}

void SegmentsTab::updateSegTestButtons(const QItemSelection &selection) {
  QModelIndexList indexes = selection.indexes();

  if (indexes.isEmpty()) {
    segTestDeleteButton->setEnabled(false);
    segTestUpButton->setEnabled(false);
    segTestDownButton->setEnabled(false);
  } else {
    segTestDeleteButton->setEnabled(true);
    if (indexes.at(0).row() == 0)
      segTestUpButton->setEnabled(false);
    else
      segTestUpButton->setEnabled(true);
    if (indexes.at(0).row() == segmentsTestModel->rowCount(QModelIndex()) - 1)
      segTestDownButton->setEnabled(false);
    else
      segTestDownButton->setEnabled(true);
  }
}

void SegmentsTab::checkAllSegData() {
  QList<SegmentsDataData> lines = segmentsDataModel->getLines();
  for (int i = 0; i < lines.size(); i++) {
    // Only check visible (non-hidden) rows
    if (!segmentsDataView->isRowHidden(i)) {
      QModelIndex index = segmentsDataModel->index(i, 0, QModelIndex());
      segmentsDataModel->setData(index, 1, Qt::EditRole);
    }
  }
}

void SegmentsTab::uncheckAllSegData() {
  QList<SegmentsDataData> lines = segmentsDataModel->getLines();
  for (int i = 0; i < lines.size(); i++) {
    // Only uncheck visible (non-hidden) rows
    if (!segmentsDataView->isRowHidden(i)) {
      QModelIndex index = segmentsDataModel->index(i, 0, QModelIndex());
      segmentsDataModel->setData(index, 0, Qt::EditRole);
    }
  }
}

void SegmentsTab::getDataFromExfor() {
  if (!pairsModel || pairsModel->numPairs() == 0) {
    QMessageBox::warning(this, tr("No Pairs Defined"),
                         tr("Please define at least one particle pair before "
                            "retrieving data from EXFOR."));
    return;
  }
  ExforDialog dialog(pairsModel, this);
  dialog.exec();
}

void SegmentsTab::filterSegDataByPairs() {
  int entranceFilterIndex = segDataEntranceFilter->currentIndex();
  int exitFilterIndex = segDataExitFilter->currentIndex();

  // Get actual pair indices from item data (0 = "All", >0 = specific pair)
  int entranceFilter = -1;
  int exitFilter = -1;
  if (entranceFilterIndex > 0) {
    entranceFilter = segDataEntranceFilter->itemData(entranceFilterIndex).toInt();
  }
  if (exitFilterIndex > 0) {
    exitFilter = segDataExitFilter->itemData(exitFilterIndex).toInt();
  }

  QList<SegmentsDataData> lines = segmentsDataModel->getLines();
  for (int i = 0; i < lines.size(); i++) {
    bool showRow = true;

    // Check entrance pair filter
    if (entranceFilter > 0 && lines[i].entrancePairIndex != entranceFilter) {
      showRow = false;
    }

    // Check exit pair filter
    if (exitFilter > 0 && lines[i].exitPairIndex != exitFilter) {
      showRow = false;
    }

    // Show or hide the row
    segmentsDataView->setRowHidden(i, !showRow);
  }
}


bool SegmentsTab::readSegDataFile(QTextStream &inStream) {
  int isActive;
  int entrancePairIndex;
  int exitPairIndex;
  double lowEnergy;
  double highEnergy;
  double lowAngle;
  double highAngle;
  int dataType;
  QString dataFile;
  double dataNorm;
  double dataNormError;
  int varyNorm;
  double phaseJ;
  int phaseL;
  double energyShift;
  double energyShiftError;
  int varyEnergyShift;

  QString line("");
  while (!inStream.atEnd() && line.trimmed() != QString("</segmentsData>")) {
    line = inStream.readLine();
    if (line.trimmed().isEmpty()) continue;
    if (!inStream.atEnd() && line.trimmed() != QString("</segmentsData>")) {
      QTextStream in(&line);
      in >> isActive >> entrancePairIndex >> exitPairIndex >> lowEnergy >> highEnergy >> lowAngle >> highAngle >> dataType;
      if (dataType == 2)
        in >> phaseJ >> phaseL;
      else {
        phaseJ = 0.;
        phaseL = 0;
      }
      in >> dataNorm >> varyNorm;
      QString restOfLine = in.readLine();
      QTextStream stm(&restOfLine);
      stm >> dataNormError;
      if (stm.status() != QTextStream::Ok) {
        dataNormError = 0.;
        dataFile = restOfLine.trimmed();
        energyShift = 0.0;
        energyShiftError = 0.0;
        varyEnergyShift = 0;
      } else {
        QString restOfLine2 = stm.readLine();
        QTextStream stm2(&restOfLine2);
        stm2 >> energyShift;
        if (stm2.status() != QTextStream::Ok) {
          energyShift = 0.0;
          energyShiftError = 0.0;
          varyEnergyShift = 0;
          dataFile = restOfLine2.trimmed();
        } else {
          stm2 >> energyShiftError >> varyEnergyShift;
          if (stm2.status() != QTextStream::Ok) {
            energyShiftError = 0.0;
            varyEnergyShift = 0;
          }
          dataFile = stm2.readLine().trimmed();
        }
      }
      if (in.status() != QTextStream::Ok) return false;
      // Initialize advanced segment data to defaults for backwards compatibility
      int isAdvanced = 0;
      int operationType = 0;
      QString componentsList = "";
      // Initialize UPOS segment data to defaults
      int isUPOS = 0;
      int secondaryDecayL = 0;
      double finalJ = 0.0;
      double delta = 0.0;
      int uposIdx = -1;  // tracks where UPOS data starts in dataParts

      // Check if there's additional advanced segment data after the file path
      // The data is appended as: " isAdvanced operationType numComponents entrance1 exit1 entrance2 exit2 ..."
      QString remainingData = dataFile;
      QStringList dataParts = remainingData.split(" ", Qt::SkipEmptyParts);

      if (dataParts.size() > 1) {
        // Extract the actual file path (first part)
        dataFile = dataParts[0];

        // Try to parse advanced segment data from remaining parts
        if (dataParts.size() >= 2) {
          bool ok;
          int advancedFlag = dataParts[1].toInt(&ok);
          if (ok && advancedFlag == 1) {
            isAdvanced = 1;
            if (dataParts.size() >= 3) {
              operationType = dataParts[2].toInt(&ok);
              if (ok && dataParts.size() >= 4) {
                int numComponents = dataParts[3].toInt(&ok);

                // A leading -1 marker means the newer format: components carry an extra
                // scaling token (4 per component instead of 3). After the marker the
                // actual component count follows at index 4.
                bool hasScalingToken = false;
                int componentsBase = 4;  // index of first component token
                if (ok && numComponents == -1 && dataParts.size() >= 5) {
                  hasScalingToken = true;
                  numComponents = dataParts[4].toInt(&ok);
                  componentsBase = 5;
                }
                int tokensPerComponent = hasScalingToken ? 4 : 3;
                bool tried4 = false;
                if (hasScalingToken && ok && dataParts.size() >= componentsBase + numComponents * tokensPerComponent) {
                  QStringList components;
                  bool fourTokenOk = true;
                  for (int i = 0; i < numComponents; i++) {
                    bool okE = true, okX = true, okA = true, okS = true;
                    int entrance = dataParts[componentsBase + i * 4].toInt(&okE);
                    int exit = dataParts[componentsBase + i * 4 + 1].toInt(&okX);
                    double angle = dataParts[componentsBase + i * 4 + 2].toDouble(&okA);
                    double scaling = dataParts[componentsBase + i * 4 + 3].toDouble(&okS);
                    if (!okE || !okX || !okA || !okS) {
                      fourTokenOk = false;
                      break;
                    }
                    QString component = QString("Entrance: %1, Exit: %2").arg(entrance).arg(exit);
                    if (angle > -900.0) component += QString(", Angle: %1").arg(angle);
                    if (scaling > -900.0) component += QString(", Scaling: %1").arg(scaling);
                    components.append(component);
                  }
                  if (fourTokenOk) {
                    componentsList = components.join(";");
                    uposIdx = componentsBase + numComponents * 4;
                    tried4 = true;
                  }
                }
                // Legacy 3-tokens-per-component format
                if (!tried4 && ok && !hasScalingToken && dataParts.size() >= 4 + numComponents * 3) {
                  QStringList components;
                  bool newFormatSuccess = true;

                  for (int i = 0; i < numComponents; i++) {
                    int entranceIdx = 4 + i * 3;
                    int exitIdx = 5 + i * 3;
                    int angleIdx = 6 + i * 3;
                    if (entranceIdx < dataParts.size() && exitIdx < dataParts.size() && angleIdx < dataParts.size()) {
                      int entrance = dataParts[entranceIdx].toInt(&ok);
                      if (ok) {
                        int exit = dataParts[exitIdx].toInt(&ok);
                        if (ok) {
                          double angle = dataParts[angleIdx].toDouble(&ok);
                          if (ok) {
                            QString component;
                            if (angle > -900.0) {  // Check if angle is not the sentinel value
                              component = QString("Entrance: %1, Exit: %2, Angle: %3").arg(entrance).arg(exit).arg(angle);
                            } else {
                              component = QString("Entrance: %1, Exit: %2").arg(entrance).arg(exit);
                            }
                            components.append(component);
                          } else {
                            newFormatSuccess = false;
                            break;
                          }
                        } else {
                          newFormatSuccess = false;
                          break;
                        }
                      } else {
                        newFormatSuccess = false;
                        break;
                      }
                    } else {
                      newFormatSuccess = false;
                      break;
                    }
                  }

                  if (newFormatSuccess) {
                    componentsList = components.join(";");
                    uposIdx = 4 + numComponents * 3;
                  }
                }
                // Fallback to old format (2 values per component: entrance, exit)
                else if (!tried4 && !hasScalingToken && ok && dataParts.size() >= 4 + numComponents * 2) {
                  QStringList components;
                  for (int i = 0; i < numComponents; i++) {
                    int entranceIdx = 4 + i * 2;
                    int exitIdx = 5 + i * 2;
                    if (entranceIdx < dataParts.size() && exitIdx < dataParts.size()) {
                      int entrance = dataParts[entranceIdx].toInt(&ok);
                      if (ok) {
                        int exit = dataParts[exitIdx].toInt(&ok);
                        if (ok) {
                          components.append(QString("Entrance: %1, Exit: %2").arg(entrance).arg(exit));
                        }
                      }
                    }
                  }
                  componentsList = components.join(";");
                  uposIdx = 4 + numComponents * 2;
                }
              }
            }
          } else {
            // Non-advanced: UPOS data starts right after the "0" flag
            uposIdx = 2;
          }
        }
      }

      // Parse UPOS data if available
      if (uposIdx >= 0 && dataParts.size() > uposIdx) {
        bool ok;
        int uposFlag = dataParts[uposIdx].toInt(&ok);
        if (ok && uposFlag == 1) {
          isUPOS = 1;
          if (dataParts.size() >= uposIdx + 4) {
            secondaryDecayL = dataParts[uposIdx + 1].toInt(&ok);
            if (ok) finalJ = dataParts[uposIdx + 2].toDouble(&ok);
            if (ok) delta = dataParts[uposIdx + 3].toDouble(&ok);
          }
        }
      }

      SegmentsDataData newLine = {isActive, entrancePairIndex, exitPairIndex, lowEnergy, highEnergy, lowAngle,
                                  highAngle, dataType, dataFile, dataNorm, dataNormError, varyNorm, phaseJ, phaseL, energyShift, energyShiftError, varyEnergyShift, isAdvanced, operationType, componentsList, isUPOS, secondaryDecayL, finalJ, delta};
      addSegDataLine(newLine);
    }
  }
  if (line.trimmed() != QString("</segmentsData>")) return false;
  return true;
}

bool SegmentsTab::writeSegDataFile(QTextStream &outStream) {
  QList<SegmentsDataData> lines = segmentsDataModel->getLines();

  for (int i = 0; i < lines.size(); i++) {
    outStream << qSetFieldWidth(15) << lines.at(i).isActive
              << qSetFieldWidth(15) << lines.at(i).entrancePairIndex
              << qSetFieldWidth(15) << lines.at(i).exitPairIndex
              << qSetFieldWidth(15) << lines.at(i).lowEnergy
              << qSetFieldWidth(15) << lines.at(i).highEnergy
              << qSetFieldWidth(15) << lines.at(i).lowAngle
              << qSetFieldWidth(15) << lines.at(i).highAngle
              << qSetFieldWidth(15) << lines.at(i).dataType;
    if (lines.at(i).dataType == 2) outStream << qSetFieldWidth(15) << lines.at(i).phaseJ
                                             << qSetFieldWidth(15) << lines.at(i).phaseL;
    outStream << qSetFieldWidth(15) << lines.at(i).dataNorm
              << qSetFieldWidth(15) << lines.at(i).varyNorm
              << qSetFieldWidth(15) << lines.at(i).dataNormError
              << qSetFieldWidth(15) << lines.at(i).energyShift
              << qSetFieldWidth(15) << lines.at(i).energyShiftError
              << qSetFieldWidth(15) << lines.at(i).varyEnergyShift
              << qSetFieldWidth(0) << " " << lines.at(i).dataFile;

    // Add advanced segment data after the file path for backwards compatibility
    if (lines.at(i).isAdvanced == 1) {
      outStream << " " << lines.at(i).isAdvanced << " " << lines.at(i).operationType;

      // Parse components and write as numbers (with optional angle and scaling)
      QStringList components = lines.at(i).componentsList.split(";", Qt::SkipEmptyParts);
      // If any component carries a scaling factor, upgrade the on-disk format
      // by writing a -1 marker before numComponents and 4 tokens per component
      // (entrance, exit, angle, scaling). Otherwise keep the legacy 3-token
      // format so older AZURE2 builds can still parse the file.
      bool anyScaling = false;
      for (const QString &component : components) {
        if (component.contains("Scaling:")) {
          anyScaling = true;
          break;
        }
      }
      auto extractAngleScaling = [&](const QString &component, double &angle, double &scaling) {
        angle = -999.0;
        scaling = -999.0;
        QStringList parts = component.split(", ");
        for (int p = 2; p < parts.size(); p++) {
          QString tok = parts[p].trimmed();
          if (tok.startsWith("Angle:"))
            angle = tok.section(":", 1).trimmed().toDouble();
          else if (tok.startsWith("Scaling:"))
            scaling = tok.section(":", 1).trimmed().toDouble();
        }
      };
      if (anyScaling) outStream << " -1";
      if (lines.at(i).operationType == 0) {  // Sum - can have multiple components
        outStream << " " << components.size();
        for (const QString &component : components) {
          QStringList parts = component.split(", ");
          if (parts.size() >= 2) {
            int entrance = parts[0].split(": ")[1].trimmed().toInt();
            int exit = parts[1].split(": ")[1].trimmed().toInt();
            double angle, scaling;
            extractAngleScaling(component, angle, scaling);
            outStream << " " << entrance << " " << exit << " " << angle;
            if (anyScaling) outStream << " " << scaling;
          }
        }
      } else {  // Ratio - only 1 component (denominator, main segment is numerator)
        // For ratio, always write exactly 1 component (only the first one)
        int numComponentsToWrite = qMin(1, components.size());
        outStream << " " << numComponentsToWrite;

        if (numComponentsToWrite > 0) {
          QStringList parts = components[0].split(", ");  // Only use first component
          if (parts.size() >= 2) {
            int entrance = parts[0].split(": ")[1].trimmed().toInt();
            int exit = parts[1].split(": ")[1].trimmed().toInt();
            double angle, scaling;
            extractAngleScaling(components[0], angle, scaling);
            outStream << " " << entrance << " " << exit << " " << angle;
            if (anyScaling) outStream << " " << scaling;
          }
        }
      }
    } else {
      // Write 0 for non-advanced segments for backwards compatibility
      outStream << " 0";
    }
    // Write UPOS data after advanced segment data
    if (lines.at(i).isUPOS == 1) {
      outStream << " 1 " << lines.at(i).secondaryDecayL
                << " " << lines.at(i).finalJ
                << " " << lines.at(i).delta;
    } else {
      outStream << " 0";
    }
    outStream << Qt::endl;
  }

  return true;
}

bool SegmentsTab::readSegTestFile(QTextStream &inStream) {
  int isActive;
  int entrancePairIndex;
  int exitPairIndex;
  double lowEnergy;
  double highEnergy;
  double energyStep;
  double lowAngle;
  double highAngle;
  double angleStep;
  int dataType;
  double phaseJ;
  int phaseL;
  int maxAngDistOrder;

  QString line("");
  while (!inStream.atEnd() && line.trimmed() != QString("</segmentsTest>")) {
    line = inStream.readLine();
    if (line.trimmed().isEmpty()) continue;
    if (!inStream.atEnd() && line.trimmed() != QString("</segmentsTest>")) {
      QTextStream in(&line);
      in >> isActive >> entrancePairIndex >> exitPairIndex >> lowEnergy >> highEnergy >> energyStep >> lowAngle >> highAngle >> angleStep >> dataType;
      if (dataType == 2) {
        in >> phaseJ >> phaseL;
      } else {
        phaseJ = 0.;
        phaseL = 0;
      }
      if (dataType == 3) {
        in >> maxAngDistOrder;
      } else {
        maxAngDistOrder = 0;
      }
      if (in.status() != QTextStream::Ok) return false;

      // Initialize advanced segment data to defaults for backwards compatibility
      int isAdvanced = 0;
      int operationType = 0;
      QString componentsList = "";

      // Try to read advanced segment parameters (if they exist)
      // For test segments, the advanced data is appended directly to the line
      in >> isAdvanced;
      if (in.status() == QTextStream::Ok && isAdvanced == 1) {
        in >> operationType;
        if (in.status() == QTextStream::Ok) {
          int numComponents = 0;
          in >> numComponents;
          // A leading -1 marker means the newer format adds a scaling token
          bool hasScalingToken = false;
          if (in.status() == QTextStream::Ok && numComponents == -1) {
            hasScalingToken = true;
            in >> numComponents;
          }
          if (in.status() == QTextStream::Ok) {
            QStringList components;

            // Try reading new format (entrance, exit, angle [, scaling])
            bool newFormatSuccess = true;
            qint64 initialPos = in.pos();

            for (int i = 0; i < numComponents; i++) {
              int entrance, exit;
              double angle;
              in >> entrance >> exit >> angle;
              if (in.status() == QTextStream::Ok) {
                double scaling = 1.0;
                bool scalingOk = true;
                if (hasScalingToken) {
                  in >> scaling;
                  scalingOk = (in.status() == QTextStream::Ok);
                }
                if (!scalingOk) {
                  newFormatSuccess = false;
                  break;
                }
                QString component;
                if (angle > -900.0) {  // Check if angle is not the sentinel value
                  component = QString("Entrance: %1, Exit: %2, Angle: %3").arg(entrance).arg(exit).arg(angle);
                } else {
                  component = QString("Entrance: %1, Exit: %2").arg(entrance).arg(exit);
                }
                if (hasScalingToken && scaling > -900.0) {
                  component += QString(", Scaling: %1").arg(scaling);
                }
                components.append(component);
              } else {
                newFormatSuccess = false;
                break;
              }
            }

            // If new format failed, try old format (entrance, exit only)
            if (!newFormatSuccess && !hasScalingToken) {
              // Reset stream position to before we tried to read components
              in.seek(initialPos);
              in.resetStatus();
              components.clear();

              for (int i = 0; i < numComponents; i++) {
                int entrance, exit;
                in >> entrance >> exit;
                if (in.status() == QTextStream::Ok) {
                  components.append(QString("Entrance: %1, Exit: %2").arg(entrance).arg(exit));
                } else {
                  break;
                }
              }
            }

            componentsList = components.join(";");
          }
        }
      } else {
        // Reset status if isAdvanced read failed (old file format)
        isAdvanced = 0;
      }

      SegmentsTestData newLine = {isActive, entrancePairIndex, exitPairIndex, lowEnergy, highEnergy, energyStep, lowAngle,
                                  highAngle, angleStep, dataType, phaseJ, phaseL, maxAngDistOrder, isAdvanced, operationType, componentsList};
      addSegTestLine(newLine);
    }
  }
  if (line.trimmed() != QString("</segmentsTest>")) return false;
  return true;
}

bool SegmentsTab::writeSegTestFile(QTextStream &outStream) {
  QList<SegmentsTestData> lines = segmentsTestModel->getLines();

  for (int i = 0; i < lines.size(); i++) {
    outStream << qSetFieldWidth(15) << lines.at(i).isActive
              << qSetFieldWidth(15) << lines.at(i).entrancePairIndex
              << qSetFieldWidth(15) << lines.at(i).exitPairIndex
              << qSetFieldWidth(15) << lines.at(i).lowEnergy
              << qSetFieldWidth(15) << lines.at(i).highEnergy
              << qSetFieldWidth(15) << lines.at(i).energyStep
              << qSetFieldWidth(15) << lines.at(i).lowAngle
              << qSetFieldWidth(15) << lines.at(i).highAngle
              << qSetFieldWidth(15) << lines.at(i).angleStep;
    if (lines.at(i).dataType == 2) {
      outStream << qSetFieldWidth(15) << lines.at(i).dataType
                << qSetFieldWidth(15) << lines.at(i).phaseJ
                << qSetFieldWidth(0) << lines.at(i).phaseL;
    } else if (lines.at(i).dataType == 3) {
      outStream << qSetFieldWidth(15) << lines.at(i).dataType
                << qSetFieldWidth(0) << lines.at(i).maxAngDistOrder;
    } else
      outStream << qSetFieldWidth(0) << lines.at(i).dataType;

    // Add advanced segment data after the existing data for backwards compatibility
    if (lines.at(i).isAdvanced == 1) {
      outStream << " " << lines.at(i).isAdvanced << " " << lines.at(i).operationType;

      // Parse components and write as numbers (with optional angle and scaling)
      QStringList components = lines.at(i).componentsList.split(";", Qt::SkipEmptyParts);
      bool anyScaling = false;
      for (const QString &component : components) {
        if (component.contains("Scaling:")) {
          anyScaling = true;
          break;
        }
      }
      auto extractAngleScaling = [&](const QString &component, double &angle, double &scaling) {
        angle = -999.0;
        scaling = -999.0;
        QStringList parts = component.split(", ");
        for (int p = 2; p < parts.size(); p++) {
          QString tok = parts[p].trimmed();
          if (tok.startsWith("Angle:"))
            angle = tok.section(":", 1).trimmed().toDouble();
          else if (tok.startsWith("Scaling:"))
            scaling = tok.section(":", 1).trimmed().toDouble();
        }
      };
      // -1 marker selects the newer "with scaling" format (see SegLine.h)
      if (anyScaling) outStream << " -1";
      if (lines.at(i).operationType == 0) {  // Sum - can have multiple components
        outStream << " " << components.size();
        for (const QString &component : components) {
          QStringList parts = component.split(", ");
          if (parts.size() >= 2) {
            int entrance = parts[0].split(": ")[1].trimmed().toInt();
            int exit = parts[1].split(": ")[1].trimmed().toInt();
            double angle, scaling;
            extractAngleScaling(component, angle, scaling);
            outStream << " " << entrance << " " << exit << " " << angle;
            if (anyScaling) outStream << " " << scaling;
          }
        }
      } else {  // Ratio - only 1 component (denominator, main segment is numerator)
        // For ratio, always write exactly 1 component (only the first one)
        int numComponentsToWrite = qMin(1, components.size());
        outStream << " " << numComponentsToWrite;

        if (numComponentsToWrite > 0) {
          QStringList parts = components[0].split(", ");  // Only use first component
          if (parts.size() >= 2) {
            int entrance = parts[0].split(": ")[1].trimmed().toInt();
            int exit = parts[1].split(": ")[1].trimmed().toInt();
            double angle, scaling;
            extractAngleScaling(components[0], angle, scaling);
            outStream << " " << entrance << " " << exit << " " << angle;
            if (anyScaling) outStream << " " << scaling;
          }
        }
      }
    } else {
      // Write 0 for non-advanced segments for backwards compatibility
      outStream << " 0";
    }
    outStream << Qt::endl;
  }

  return true;
}

void SegmentsTab::setPairsModel(PairsModel *model) {
  pairsModel = model;
  segmentsDataModel->setPairsModel(model);
  segmentsTestModel->setPairsModel(model);
  updateFilterComboboxes(model);
}

void SegmentsTab::updateFilterComboboxes(PairsModel *model) {
  pairsModel = model;
  updateFilterComboboxes();
}

void SegmentsTab::updateFilterComboboxes() {
  // Preserve the current filter selection across the rebuild, and block signals
  // so removing/re-adding items does not fire filterSegDataByPairs() with a
  // transient "All" selection (which would clear the filter on any checkbox toggle).
  int savedEntrance = segDataEntranceFilter->currentIndex() > 0
      ? segDataEntranceFilter->currentData().toInt()
      : -1;
  int savedExit = segDataExitFilter->currentIndex() > 0
      ? segDataExitFilter->currentData().toInt()
      : -1;
  segDataEntranceFilter->blockSignals(true);
  segDataExitFilter->blockSignals(true);

  // Clear existing items except "All"
  while (segDataEntranceFilter->count() > 1) {
    segDataEntranceFilter->removeItem(1);
  }
  while (segDataExitFilter->count() > 1) {
    segDataExitFilter->removeItem(1);
  }

  // Extract unique entrance and exit pairs from actual segment data
  QList<SegmentsDataData> lines = segmentsDataModel->getLines();
  QSet<int> entrancePairs;
  QSet<int> exitPairs;

  for (int i = 0; i < lines.size(); i++) {
    // Only add valid pair indices (skip -1 and other negative values)
    if (lines[i].entrancePairIndex > 0) {
      entrancePairs.insert(lines[i].entrancePairIndex);
    }
    if (lines[i].exitPairIndex > 0) {
      exitPairs.insert(lines[i].exitPairIndex);
    }
  }

  // Convert to sorted lists
  QList<int> sortedEntrancePairs = entrancePairs.values();
  QList<int> sortedExitPairs = exitPairs.values();
  std::sort(sortedEntrancePairs.begin(), sortedEntrancePairs.end());
  std::sort(sortedExitPairs.begin(), sortedExitPairs.end());

  // Add entrance pairs to filter
  for (int i = 0; i < sortedEntrancePairs.size(); i++) {
    int pairIndex = sortedEntrancePairs[i];
    QString label;
    if (pairsModel && pairIndex > 0 && pairIndex <= pairsModel->numPairs()) {
      QList<PairsData> pairs = pairsModel->getPairs();
      QString rawLabel = pairsModel->getParticleLabel(pairs[pairIndex - 1]);
      // Strip HTML tags for display
      QTextDocument doc;
      doc.setHtml(rawLabel);
      label = doc.toPlainText();
      segDataEntranceFilter->addItem(QString("Pair %1: %2").arg(pairIndex).arg(label), pairIndex);
    } else {
      segDataEntranceFilter->addItem(QString("Pair %1").arg(pairIndex), pairIndex);
    }
  }

  // Add exit pairs to filter
  for (int i = 0; i < sortedExitPairs.size(); i++) {
    int pairIndex = sortedExitPairs[i];
    QString label;
    if (pairsModel && pairIndex > 0 && pairIndex <= pairsModel->numPairs()) {
      QList<PairsData> pairs = pairsModel->getPairs();
      QString rawLabel = pairsModel->getParticleLabel(pairs[pairIndex - 1]);
      // Strip HTML tags for display
      QTextDocument doc;
      doc.setHtml(rawLabel);
      label = doc.toPlainText();
      segDataExitFilter->addItem(QString("Pair %1: %2").arg(pairIndex).arg(label), pairIndex);
    } else {
      segDataExitFilter->addItem(QString("Pair %1").arg(pairIndex), pairIndex);
    }
  }

  // Restore the previous selection (fall back to "All" if that pair is gone),
  // then re-apply the filter once to the current rows.
  int ei = savedEntrance > 0 ? segDataEntranceFilter->findData(savedEntrance) : 0;
  segDataEntranceFilter->setCurrentIndex(ei >= 0 ? ei : 0);
  int xi = savedExit > 0 ? segDataExitFilter->findData(savedExit) : 0;
  segDataExitFilter->setCurrentIndex(xi >= 0 ? xi : 0);
  segDataEntranceFilter->blockSignals(false);
  segDataExitFilter->blockSignals(false);
  filterSegDataByPairs();
}

void SegmentsTab::reset() {
  segmentsDataModel->removeRows(0, segmentsDataModel->getLines().size(), QModelIndex());
  segmentsTestModel->removeRows(0, segmentsTestModel->getLines().size(), QModelIndex());
}

void SegmentsTab::showInfo(int which, QString title) {
  if (which < infoText.size()) {
    if (!infoDialog[which]) {
      infoDialog[which] = new InfoDialog(infoText[which], this, title);
      infoDialog[which]->setAttribute(Qt::WA_DeleteOnClose);
      infoDialog[which]->show();
    } else
      infoDialog[which]->raise();
  }
}
