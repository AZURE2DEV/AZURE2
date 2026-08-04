#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QFormLayout>
#include <QGroupBox>
#include <QComboBox>
#include <QRadioButton>
#include <QCheckBox>
#include <QPushButton>
#include <QListView>
#include <QListWidget>
#include <QListWidgetItem>
#include <QSplitter>
#include <QLabel>
#include <QLineEdit>
#include <QSpinBox>
#include <QDialogButtonBox>
#include <QColorDialog>
#include <QPixmap>
#include <QIcon>
#include <QSet>

#include "PlotTab.h"
#include "Config.h"
#include "AZUREPlot.h"
#include "SegmentsDataModel.h"
#include "SegmentsTestModel.h"
#include "LevelsModel.h"
#include "RichTextDelegate.h"
#include "InfoDialog.h"
#include <iostream>

void ChannelFilterProxy::setEntranceFilter(int pairIndex) {
  entranceFilter_ = pairIndex;
  invalidateFilter();
}

void ChannelFilterProxy::setExitFilter(int pairIndex) {
  exitFilter_ = pairIndex;
  invalidateFilter();
}

bool ChannelFilterProxy::entranceExitAccepts(int source_row, const QModelIndex& source_parent) const {
  if(entranceFilter_<0 && exitFilter_<0) return true;
  QAbstractItemModel* src = sourceModel();
  if(!src) return true;
  if(entranceFilter_>=0) {
    QModelIndex idx = src->index(source_row, 1, source_parent);
    if(src->data(idx, Qt::EditRole).toInt() != entranceFilter_) return false;
  }
  if(exitFilter_>=0) {
    QModelIndex idx = src->index(source_row, 2, source_parent);
    if(src->data(idx, Qt::EditRole).toInt() != exitFilter_) return false;
  }
  return true;
}

QVariant SegTestProxyModel::data(const QModelIndex& index, int role) const {
  if (index.isValid() && role == Qt::DisplayRole) {
    QModelIndex sourceIndex = mapToSource(index);
    return QString("#%1: %2").arg(sourceIndex.row()+1).arg(static_cast<SegmentsTestModel*>(sourceModel())->getReactionLabel(sourceIndex));
  }
  return QVariant();
}


bool SegTestProxyModel::filterAcceptsRow(int source_row, const QModelIndex &source_parent) const {
  if(!QSortFilterProxyModel::filterAcceptsRow(source_row,source_parent)) return false;
  if(!entranceExitAccepts(source_row, source_parent)) return false;
  SegmentsTestModel* model = static_cast<SegmentsTestModel*>(sourceModel());
  QModelIndex source_index = model->index(source_row, 9, source_parent);
  int dataType = model->data(source_index,Qt::EditRole).toInt();
  return dataType != 3;
}

QVariant SegDataProxyModel::data(const QModelIndex& index, int role) const {
  if (index.isValid() && role == Qt::DisplayRole) {
    QModelIndex sourceIndex = mapToSource(index);
    return QString("#%1: %2").arg(sourceIndex.row()+1).arg(static_cast<SegmentsDataModel*>(sourceModel())->getReactionLabel(sourceIndex));
  }
  return QVariant();
}

bool SegDataProxyModel::filterAcceptsRow(int source_row, const QModelIndex &source_parent) const {
  if(!QSortFilterProxyModel::filterAcceptsRow(source_row, source_parent)) return false;
  return entranceExitAccepts(source_row, source_parent);
}


// ---------------- CurveStyleDialog ----------------

static QIcon colorSwatch(const QColor& c, int size = 16) {
  QPixmap pm(size, size);
  pm.fill(c);
  return QIcon(pm);
}

CurveStyleDialog::CurveStyleDialog(PlotEntry* entry, QWidget* parent) :
  QDialog(parent), entry_(entry),
  color_(entry ? entry->color() : QColor(Qt::black)),
  fitColor_(entry && entry->fitColor().isValid() ? entry->fitColor() : color_) {
  setWindowTitle(tr("Curve Style"));

  QFormLayout* form = new QFormLayout;

  labelEdit_ = new QLineEdit(entry_ ? entry_->label() : QString());
  form->addRow(tr("Label:"), labelEdit_);

  // For data entries (type==0) the first color controls the data points and
  // the second controls the calculation/fit line. For calculation-only
  // entries (type==1) only a single color is shown.
  const bool hasFitColor = entry_ && entry_->type()==0;

  colorButton_ = new QPushButton;
  colorButton_->setMinimumWidth(80);
  connect(colorButton_, SIGNAL(clicked()), this, SLOT(chooseColor()));
  updateColorButton();
  form->addRow(hasFitColor ? tr("Data color:") : tr("Color:"), colorButton_);

  fitColorButton_ = new QPushButton;
  fitColorButton_->setMinimumWidth(80);
  connect(fitColorButton_, SIGNAL(clicked()), this, SLOT(chooseFitColor()));
  updateFitColorButton();
  if(hasFitColor) form->addRow(tr("Calculation color:"), fitColorButton_);
  else fitColorButton_->hide();

  symbolCombo_ = new QComboBox;
  symbolCombo_->addItem(tr("Circle"), QwtSymbol::Ellipse);
  symbolCombo_->addItem(tr("Square"), QwtSymbol::Rect);
  symbolCombo_->addItem(tr("Diamond"), QwtSymbol::Diamond);
  symbolCombo_->addItem(tr("Triangle"), QwtSymbol::Triangle);
  symbolCombo_->addItem(tr("Triangle (down)"), QwtSymbol::DTriangle);
  symbolCombo_->addItem(tr("Triangle (up)"), QwtSymbol::UTriangle);
  symbolCombo_->addItem(tr("Cross (+)"), QwtSymbol::Cross);
  symbolCombo_->addItem(tr("X"), QwtSymbol::XCross);
  symbolCombo_->addItem(tr("Hexagon"), QwtSymbol::Hexagon);
  symbolCombo_->addItem(tr("Star"), QwtSymbol::Star1);
  symbolCombo_->addItem(tr("None"), QwtSymbol::NoSymbol);
  if(entry_) {
    for(int i=0;i<symbolCombo_->count();i++) {
      if(symbolCombo_->itemData(i).toInt() == (int)entry_->symbolStyle()) {
        symbolCombo_->setCurrentIndex(i);
        break;
      }
    }
  }
  form->addRow(tr("Marker:"), symbolCombo_);

  symbolSizeSpin_ = new QSpinBox;
  symbolSizeSpin_->setRange(2, 30);
  symbolSizeSpin_->setValue(entry_ ? entry_->symbolSize() : 6);
  form->addRow(tr("Marker size:"), symbolSizeSpin_);

  lineWidthSpin_ = new QSpinBox;
  lineWidthSpin_->setRange(1, 10);
  lineWidthSpin_->setValue(entry_ ? entry_->lineWidth() : 2);
  form->addRow(tr("Line width:"), lineWidthSpin_);

  // For test entries (type 1), there are no symbols.
  if(entry_ && entry_->type()==1) {
    symbolCombo_->setEnabled(false);
    symbolSizeSpin_->setEnabled(false);
  }

  QDialogButtonBox* buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
  connect(buttons, SIGNAL(accepted()), this, SLOT(accept()));
  connect(buttons, SIGNAL(rejected()), this, SLOT(reject()));

  QVBoxLayout* layout = new QVBoxLayout;
  layout->addLayout(form);
  layout->addWidget(buttons);
  setLayout(layout);
}

void CurveStyleDialog::updateColorButton() {
  colorButton_->setIcon(colorSwatch(color_, 20));
  colorButton_->setText(color_.name());
}

void CurveStyleDialog::updateFitColorButton() {
  fitColorButton_->setIcon(colorSwatch(fitColor_, 20));
  fitColorButton_->setText(fitColor_.name());
}

void CurveStyleDialog::chooseColor() {
  QColor chosen = QColorDialog::getColor(color_, this, tr("Select Curve Color"));
  if(chosen.isValid()) {
    color_ = chosen;
    updateColorButton();
  }
}

void CurveStyleDialog::chooseFitColor() {
  QColor chosen = QColorDialog::getColor(fitColor_, this, tr("Select Calculation Color"));
  if(chosen.isValid()) {
    fitColor_ = chosen;
    updateFitColorButton();
  }
}

void CurveStyleDialog::accept() {
  if(entry_) {
    QString trimmed = labelEdit_->text().trimmed();
    if(!trimmed.isEmpty()) entry_->setLabel(trimmed);
    entry_->setColor(color_);
    // For calculation-only entries (type==1) the fit follows color_; keep
    // fitColor_ in sync. For type==0 store the user-selected fit color.
    if(entry_->type()==0) entry_->setFitColor(fitColor_);
    else entry_->setFitColor(color_);
    entry_->setSymbolStyle(static_cast<QwtSymbol::Style>(symbolCombo_->currentData().toInt()));
    entry_->setSymbolSize(symbolSizeSpin_->value());
    entry_->setLineWidth(lineWidthSpin_->value());
  }
  QDialog::accept();
}


// ---------------- PlotTab ----------------

PlotTab::PlotTab(Config& config, SegmentsDataModel* dataModel, SegmentsTestModel* testModel, LevelsModel* levelsModel, QWidget* parent) :
  QWidget(parent), configure(config) {
  azurePlot = new AZUREPlot(this,this);
  azurePlot->setLevelsModel(levelsModel);

  segDataProxyModel = new SegDataProxyModel(this);
  segDataProxyModel->setSourceModel(dataModel);
  segDataProxyModel->setDynamicSortFilter(true);
  segDataProxyModel->setFilterKeyColumn(0);
  segDataProxyModel->setFilterRole(Qt::CheckStateRole);
  segDataProxyModel->setFilterRegExp(QString("%1").arg(Qt::Checked));
  segTestProxyModel = new SegTestProxyModel(this);
  segTestProxyModel->setSourceModel(testModel);
  segTestProxyModel->setDynamicSortFilter(true);
  segTestProxyModel->setFilterKeyColumn(0);
  segTestProxyModel->setFilterRole(Qt::CheckStateRole);
  segTestProxyModel->setFilterRegExp(QString("%1").arg(Qt::Checked));

  // ---- Right: plot area + controls ----
  QVBoxLayout* rightLayout = new QVBoxLayout;

  QGridLayout *topLayout = new QGridLayout;

  QGroupBox *xAxisBox = new QGroupBox(tr("X Axis"));
  xAxisTypeCombo=new QComboBox;
  xAxisTypeCombo->addItem(tr("CoM Energy"));
  xAxisTypeCombo->addItem(tr("Excitation Energy"));
  xAxisTypeCombo->addItem(tr("CoM Angle"));
  connect(xAxisTypeCombo,SIGNAL(activated(int)),this,SLOT(xAxisTypeChanged()));
  xAxisTypeCombo->setCurrentIndex(0);
  azurePlot->setXAxisType(0);
  xAxisIsLogCheck = new QCheckBox(tr("Log scale"));
  connect(xAxisIsLogCheck,SIGNAL(toggled(bool)),this,SLOT(xAxisLogScaleChanged(bool)));
  QHBoxLayout* xAxisLayout = new QHBoxLayout;
  xAxisLayout->setContentsMargins(8,4,8,4);
  xAxisLayout->addWidget(xAxisTypeCombo,1);
  xAxisLayout->addWidget(xAxisIsLogCheck);
  xAxisBox->setLayout(xAxisLayout);

  QGroupBox *yAxisBox = new QGroupBox(tr("Y Axis"));
  yAxisXSButton = new QRadioButton(tr("Cross Section"));
  connect(yAxisXSButton,SIGNAL(toggled(bool)),this,SLOT(yAxisTypeChanged()));
  yAxisXSButton->setChecked(true);
  yAxisSFButton = new QRadioButton(tr("S-Factor"));
  connect(yAxisSFButton,SIGNAL(toggled(bool)),this,SLOT(yAxisTypeChanged()));
  yAxisIsLogCheck = new QCheckBox(tr("Log scale"));
  connect(yAxisIsLogCheck,SIGNAL(toggled(bool)),this,SLOT(yAxisLogScaleChanged(bool)));
  yAxisIsLogCheck->setChecked(true);
  QHBoxLayout* yAxisLayout = new QHBoxLayout;
  yAxisLayout->setContentsMargins(8,4,8,4);
  yAxisLayout->addWidget(yAxisXSButton);
  yAxisLayout->addWidget(yAxisSFButton);
  yAxisLayout->addWidget(yAxisIsLogCheck);
  yAxisBox->setLayout(yAxisLayout);

  QGroupBox *displayBox = new QGroupBox(tr("Display"));
  gridCheck = new QCheckBox(tr("Grid"));
  connect(gridCheck, SIGNAL(toggled(bool)), this, SLOT(gridToggled(bool)));
  legendCheck = new QCheckBox(tr("Legend"));
  legendCheck->setChecked(true);
  connect(legendCheck, SIGNAL(toggled(bool)), this, SLOT(legendToggled(bool)));
  levelsCheck = new QCheckBox(tr("Levels"));
  connect(levelsCheck, SIGNAL(toggled(bool)), this, SLOT(levelsToggled(bool)));
  bandCheck = new QCheckBox(tr("Uncertainty"));
  bandCheck->setToolTip(tr("Shade the 1-sigma analytic uncertainty band around the "
                           "calculation. Requires that the run was performed with the "
                           "Run-tab \"Uncertainty band\" option enabled (which writes the "
                           ".band files this reads)."));
  connect(bandCheck, SIGNAL(toggled(bool)), this, SLOT(bandToggled(bool)));
  QHBoxLayout* displayLayout = new QHBoxLayout;
  displayLayout->setContentsMargins(8,4,8,4);
  displayLayout->addWidget(gridCheck);
  displayLayout->addWidget(legendCheck);
  displayLayout->addWidget(levelsCheck);
  displayLayout->addWidget(bandCheck);
  displayBox->setLayout(displayLayout);
  displayBox->setSizePolicy(QSizePolicy::Maximum, QSizePolicy::Preferred);

  // Axis boxes on the first row; Display on its own row below, left-aligned.
  topLayout->addWidget(xAxisBox,0,0);
  topLayout->addWidget(yAxisBox,0,1);
  topLayout->addWidget(displayBox,1,0,1,2,Qt::AlignLeft);
  topLayout->setColumnStretch(0,1);
  topLayout->setColumnStretch(1,1);

  rightLayout->addLayout(topLayout);
  rightLayout->addWidget(azurePlot, 1);

  // Plotted curves customization panel
  QGroupBox* curvesBox = new QGroupBox(tr("Plotted Curves (double-click to customize)"));
  curveList = new QListWidget;
  curveList->setMaximumHeight(110);
  curveList->setSelectionMode(QAbstractItemView::SingleSelection);
  connect(curveList, SIGNAL(itemDoubleClicked(QListWidgetItem*)), this, SLOT(editSelectedCurveStyle()));
  QHBoxLayout* curvesLayout = new QHBoxLayout;
  curvesLayout->setContentsMargins(8,4,8,4);
  curvesLayout->addWidget(curveList);
  QPushButton* customizeBtn = new QPushButton(tr("Customize..."));
  connect(customizeBtn, SIGNAL(clicked()), this, SLOT(editSelectedCurveStyle()));
  QVBoxLayout* curveBtnLayout = new QVBoxLayout;
  curveBtnLayout->addWidget(customizeBtn);
  curveBtnLayout->addStretch();
  curvesLayout->addLayout(curveBtnLayout);
  curvesBox->setLayout(curvesLayout);
  rightLayout->addWidget(curvesBox);

  // ---- Left: source selection ----
  dataSegmentSelectorList = new QListView;
  testSegmentSelectorList = new QListView;
  dataSegmentSelectorList->setAttribute(Qt::WA_MacShowFocusRect, 0);
  testSegmentSelectorList->setAttribute(Qt::WA_MacShowFocusRect, 0);
  dataSegmentSelectorList->setModel(segDataProxyModel);
  testSegmentSelectorList->setModel(segTestProxyModel);
  dataSegmentSelectorList->setItemDelegate(new RichTextDelegate());
  testSegmentSelectorList->setItemDelegate(new RichTextDelegate());
  dataSegmentSelectorList->setSelectionMode(QAbstractItemView::MultiSelection);
  testSegmentSelectorList->setSelectionMode(QAbstractItemView::MultiSelection);
  dataSegmentSelectorList->setResizeMode(QListView::Adjust);
  testSegmentSelectorList->setResizeMode(QListView::Adjust);

  dataInChannelCombo = new QComboBox;
  dataOutChannelCombo = new QComboBox;
  testInChannelCombo = new QComboBox;
  testOutChannelCombo = new QComboBox;

  QGroupBox *dataSegmentSelectorBox = new QGroupBox(tr("Segments From Data"));
  QGridLayout *dataSegmentSelectorLayout = new QGridLayout;
  dataSegmentSelectorLayout->setContentsMargins(8,6,8,6);
  dataSegmentSelectorLayout->setHorizontalSpacing(6);
  dataSegmentSelectorLayout->addWidget(new QLabel(tr("In:")),  0,0);
  dataSegmentSelectorLayout->addWidget(dataInChannelCombo,     0,1);
  dataSegmentSelectorLayout->addWidget(new QLabel(tr("Out:")), 0,2);
  dataSegmentSelectorLayout->addWidget(dataOutChannelCombo,    0,3);
  dataSegmentSelectorLayout->addWidget(dataSegmentSelectorList,1,0,1,4);
  dataSegmentSelectorLayout->setColumnStretch(1,1);
  dataSegmentSelectorLayout->setColumnStretch(3,1);
  dataSegmentSelectorBox->setLayout(dataSegmentSelectorLayout);

  QGroupBox *testSegmentSelectorBox = new QGroupBox(tr("Segments Without Data"));
  QGridLayout *testSegmentSelectorLayout = new QGridLayout;
  testSegmentSelectorLayout->setContentsMargins(8,6,8,6);
  testSegmentSelectorLayout->setHorizontalSpacing(6);
  testSegmentSelectorLayout->addWidget(new QLabel(tr("In:")),  0,0);
  testSegmentSelectorLayout->addWidget(testInChannelCombo,     0,1);
  testSegmentSelectorLayout->addWidget(new QLabel(tr("Out:")), 0,2);
  testSegmentSelectorLayout->addWidget(testOutChannelCombo,    0,3);
  testSegmentSelectorLayout->addWidget(testSegmentSelectorList,1,0,1,4);
  testSegmentSelectorLayout->setColumnStretch(1,1);
  testSegmentSelectorLayout->setColumnStretch(3,1);
  testSegmentSelectorBox->setLayout(testSegmentSelectorLayout);

  refreshButton = new QPushButton(tr("&Draw"));
  refreshButton->setDefault(true);
  connect(refreshButton,SIGNAL(clicked()),this,SLOT(draw()));
  clearButton = new QPushButton(tr("Clear"));
  connect(clearButton,SIGNAL(clicked()),azurePlot,SLOT(clearEntries()));
  connect(clearButton,SIGNAL(clicked()),this,SLOT(rebuildCurveList()));
  exportButton = new QPushButton(tr("Export..."));
  connect(exportButton,SIGNAL(clicked()),azurePlot,SLOT(exportPlot()));
  printButton = new QPushButton(tr("Print..."));
  connect(printButton,SIGNAL(clicked()),azurePlot,SLOT(print()));
  QGridLayout *buttonLayout = new QGridLayout;
  buttonLayout->addWidget(refreshButton, 0, 0);
  buttonLayout->addWidget(clearButton,   0, 1);
  buttonLayout->addWidget(exportButton,  1, 0);
  buttonLayout->addWidget(printButton,   1, 1);

  QWidget* leftWidget = new QWidget;
  QVBoxLayout* leftLayout = new QVBoxLayout(leftWidget);
  leftLayout->setContentsMargins(0,0,0,0);
  leftLayout->addWidget(dataSegmentSelectorBox);
  leftLayout->addWidget(testSegmentSelectorBox);
  leftLayout->addLayout(buttonLayout);

  QWidget* rightWidget = new QWidget;
  rightWidget->setLayout(rightLayout);

  QSplitter* splitter = new QSplitter(Qt::Horizontal);
  splitter->addWidget(leftWidget);
  splitter->addWidget(rightWidget);
  splitter->setStretchFactor(0, 0);
  splitter->setStretchFactor(1, 1);
  splitter->setChildrenCollapsible(false);
  splitter->setSizes(QList<int>() << 280 << 700);

  QHBoxLayout* mainLayout = new QHBoxLayout;
  mainLayout->setContentsMargins(4,4,4,4);
  mainLayout->addWidget(splitter);
  setLayout(mainLayout);

  // Populate channel filters initially and keep them in sync.
  rebuildChannelFilterOptions();
  connect(dataModel, SIGNAL(rowsInserted(const QModelIndex&, int, int)),
          this, SLOT(rebuildChannelFilterOptions()));
  connect(dataModel, SIGNAL(rowsRemoved(const QModelIndex&, int, int)),
          this, SLOT(rebuildChannelFilterOptions()));
  connect(dataModel, SIGNAL(dataChanged(const QModelIndex&, const QModelIndex&, const QVector<int>&)),
          this, SLOT(rebuildChannelFilterOptions()));
  connect(testModel, SIGNAL(rowsInserted(const QModelIndex&, int, int)),
          this, SLOT(rebuildChannelFilterOptions()));
  connect(testModel, SIGNAL(rowsRemoved(const QModelIndex&, int, int)),
          this, SLOT(rebuildChannelFilterOptions()));
  connect(testModel, SIGNAL(dataChanged(const QModelIndex&, const QModelIndex&, const QVector<int>&)),
          this, SLOT(rebuildChannelFilterOptions()));

  connect(dataInChannelCombo,  SIGNAL(activated(int)), this, SLOT(dataChannelFilterChanged()));
  connect(dataOutChannelCombo, SIGNAL(activated(int)), this, SLOT(dataChannelFilterChanged()));
  connect(testInChannelCombo,  SIGNAL(activated(int)), this, SLOT(testChannelFilterChanged()));
  connect(testOutChannelCombo, SIGNAL(activated(int)), this, SLOT(testChannelFilterChanged()));
}

void PlotTab::populateChannelFilter(QComboBox* combo, QAbstractItemModel* model, int column) {
  int previous = combo->currentData().isValid() ? combo->currentData().toInt() : -1;
  combo->blockSignals(true);
  combo->clear();
  combo->addItem(tr("All"), -1);
  QList<int> values;
  for(int r=0; r<model->rowCount(); r++) {
    int v = model->data(model->index(r, column), Qt::EditRole).toInt();
    if(!values.contains(v)) values.append(v);
  }
  std::sort(values.begin(), values.end());
  for(int v : values) combo->addItem(QString::number(v), v);
  // Restore previous selection if still present.
  int idx = combo->findData(previous);
  combo->setCurrentIndex(idx>=0 ? idx : 0);
  combo->blockSignals(false);
}

void PlotTab::rebuildChannelFilterOptions() {
  populateChannelFilter(dataInChannelCombo,  segDataProxyModel->sourceModel(), 1);
  populateChannelFilter(dataOutChannelCombo, segDataProxyModel->sourceModel(), 2);
  populateChannelFilter(testInChannelCombo,  segTestProxyModel->sourceModel(), 1);
  populateChannelFilter(testOutChannelCombo, segTestProxyModel->sourceModel(), 2);
}

void PlotTab::dataChannelFilterChanged() {
  segDataProxyModel->setEntranceFilter(dataInChannelCombo->currentData().toInt());
  segDataProxyModel->setExitFilter(dataOutChannelCombo->currentData().toInt());
}

void PlotTab::testChannelFilterChanged() {
  segTestProxyModel->setEntranceFilter(testInChannelCombo->currentData().toInt());
  segTestProxyModel->setExitFilter(testOutChannelCombo->currentData().toInt());
}

QList<PlotEntry*> PlotTab::getDataSegments() {
  QList<PlotEntry*> dataSegmentPlotEntries;
  QModelIndexList indexes = dataSegmentSelectorList->selectionModel()->selectedIndexes();
  for(int i = 0; i< indexes.size(); i++) {
    QModelIndex sourceIndex =
      segDataProxyModel->mapToSource(segDataProxyModel->index(indexes[i].row(),1,QModelIndex()));
    int entranceKey = segDataProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toInt();
    sourceIndex = segDataProxyModel->mapToSource(segDataProxyModel->index(indexes[i].row(),2,QModelIndex()));
    int exitKey = segDataProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toInt();
    sourceIndex = segDataProxyModel->mapToSource(segDataProxyModel->index(indexes[i].row(),7,QModelIndex()));
    int dataType = segDataProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toInt();
    if(dataType==7) selectionHasAnalyzingPower_ = true;
    QString filename = (dataType==3) ?
      QString::fromStdString(configure.outputdir)+QString("AZUREOut_aa=%1_TOTAL_CAPTURE.out").arg(entranceKey) :
      QString::fromStdString(configure.outputdir)+QString("AZUREOut_aa=%1_R=%2.out").arg(entranceKey).arg(exitKey);
    sourceIndex = segDataProxyModel->mapToSource(segDataProxyModel->index(indexes[i].row(),8,QModelIndex()));
    QString segmentDataFile = segDataProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toString();
    int numPreviousInBlock = 0;
    for(int j =0; j<indexes[i].row(); j++) {
      sourceIndex = segDataProxyModel->mapToSource(segDataProxyModel->index(j,1,QModelIndex()));
      int previousEntranceKey = segDataProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toInt();
      sourceIndex = segDataProxyModel->mapToSource(segDataProxyModel->index(j,2,QModelIndex()));
      int previousExitKey = segDataProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toInt();
      if(previousEntranceKey==entranceKey&&previousExitKey==exitKey) numPreviousInBlock++;
    }
    PlotEntry* newPlotEntry = new PlotEntry(0,entranceKey,exitKey,numPreviousInBlock,filename);
    if(!segmentDataFile.isEmpty()) {
      newPlotEntry->setLabel(PlotEntry::labelFromFilename(segmentDataFile));
    }
    dataSegmentPlotEntries.push_back(newPlotEntry);
  }
  return dataSegmentPlotEntries;
}

QList<PlotEntry*> PlotTab::getTestSegments() {
  QList<PlotEntry*> testSegmentPlotEntries;
  QModelIndexList indexes = testSegmentSelectorList->selectionModel()->selectedIndexes();
  for(int i = 0; i< indexes.size(); i++) {
    QModelIndex sourceIndex =
      segTestProxyModel->mapToSource(segTestProxyModel->index(indexes[i].row(),1,QModelIndex()));
    int entranceKey = segTestProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toInt();
    sourceIndex = segTestProxyModel->mapToSource(segTestProxyModel->index(indexes[i].row(),2,QModelIndex()));
    int exitKey = segTestProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toInt();
    sourceIndex = segTestProxyModel->mapToSource(segTestProxyModel->index(indexes[i].row(),9,QModelIndex()));
    int dataType = segTestProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toInt();
    if(dataType==7) selectionHasAnalyzingPower_ = true;
    QString filename = (dataType==4) ?
      QString::fromStdString(configure.outputdir)+QString("AZUREOut_aa=%1_TOTAL_CAPTURE.extrap").arg(entranceKey) :
      QString::fromStdString(configure.outputdir)+QString("AZUREOut_aa=%1_R=%2.extrap").arg(entranceKey).arg(exitKey);
    int numPreviousInBlock = 0;
    for(int j =0; j<indexes[i].row(); j++) {
      sourceIndex = segTestProxyModel->mapToSource(segTestProxyModel->index(j,1,QModelIndex()));
      int previousEntranceKey = segTestProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toInt();
      sourceIndex = segTestProxyModel->mapToSource(segTestProxyModel->index(j,2,QModelIndex()));
      int previousExitKey = segTestProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toInt();
      if(previousEntranceKey==entranceKey&&previousExitKey==exitKey) numPreviousInBlock++;
    }
    PlotEntry* newPlotEntry = new PlotEntry(1,entranceKey,exitKey,numPreviousInBlock,filename);
    testSegmentPlotEntries.push_back(newPlotEntry);
  }
  return testSegmentPlotEntries;
}

void PlotTab::draw() {
  selectionHasAnalyzingPower_ = false;
  QList<PlotEntry*> entries = getDataSegments();
  entries.append(getTestSegments());
  // An analyzing power is a ratio lying in [-1,1] and negative over much of its
  // range. A logarithmic axis -- the default here -- simply cannot show it, and
  // an S-factor conversion has no meaning for it. Switch both off rather than
  // leave the user with a plot that looks empty.
  if(selectionHasAnalyzingPower_) {
    if(yAxisIsLogCheck->isChecked()) yAxisIsLogCheck->setChecked(false);
    if(yAxisSFButton->isChecked()) yAxisXSButton->setChecked(true);
  }
  azurePlot->draw(entries);
  rebuildCurveList();
}

void PlotTab::rebuildCurveList() {
  curveList->clear();
  const QList<PlotEntry*>& entries = azurePlot->getEntries();
  for(int i=0; i<entries.size(); i++) {
    PlotEntry* e = entries[i];
    QString prefix = (e->type()==0) ? tr("data") : tr("fit");
    QString text = QString("[%1] %2").arg(prefix).arg(e->label());
    QListWidgetItem* item = new QListWidgetItem(colorSwatch(e->color(), 14), text);
    item->setData(Qt::UserRole, QVariant::fromValue(reinterpret_cast<quintptr>(e)));
    curveList->addItem(item);
  }
}

void PlotTab::editSelectedCurveStyle() {
  QListWidgetItem* item = curveList->currentItem();
  if(!item) return;
  PlotEntry* entry = reinterpret_cast<PlotEntry*>(item->data(Qt::UserRole).value<quintptr>());
  if(!entry) return;
  CurveStyleDialog dlg(entry, this);
  if(dlg.exec()==QDialog::Accepted) {
    azurePlot->redrawEntries();
    rebuildCurveList();
  }
}

void PlotTab::xAxisTypeChanged() {
  azurePlot->setXAxisType(xAxisTypeCombo->currentIndex());
}

void PlotTab::yAxisTypeChanged() {
  if(yAxisXSButton->isChecked())
    azurePlot->setYAxisType(0);
  else if(yAxisSFButton->isChecked())
    azurePlot->setYAxisType(1);
}

void PlotTab::xAxisLogScaleChanged(bool checked) {
  azurePlot->setXAxisLog(checked);
}

void PlotTab::yAxisLogScaleChanged(bool checked) {
  azurePlot->setYAxisLog(checked);
}

void PlotTab::gridToggled(bool checked) {
  azurePlot->setGridVisible(checked);
}

void PlotTab::legendToggled(bool checked) {
  azurePlot->setLegendVisible(checked);
}

void PlotTab::levelsToggled(bool checked) {
  azurePlot->setLevelsVisible(checked);
}

void PlotTab::bandToggled(bool checked) {
  azurePlot->setBandVisible(checked);
}

void PlotTab::reset() {
  azurePlot->clearEntries();
  rebuildCurveList();
  xAxisTypeCombo->setCurrentIndex(0);
  xAxisIsLogCheck->setChecked(false);
  yAxisXSButton->setChecked(true);
  yAxisIsLogCheck->setChecked(true);
  gridCheck->setChecked(false);
  legendCheck->setChecked(true);
  levelsCheck->setChecked(false);
  bandCheck->setChecked(false);
}

void PlotTab::showInfo(int which,QString title) {
  if(which<(int)infoText.size()) {
    if(!infoDialog[which]) {
      infoDialog[which] = new InfoDialog(infoText[which],this,title);
      infoDialog[which]->setAttribute(Qt:: WA_DeleteOnClose);
      infoDialog[which]->show();
    } else infoDialog[which]->raise();
  }
}
