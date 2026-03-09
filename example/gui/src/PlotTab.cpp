#include "PlotTab.h"
#include "Config.h"
#include "AZUREPlot.h"
#include "SegmentsDataModel.h"
#include "SegmentsTestModel.h"
#include "RichTextDelegate.h"
#include "InfoDialog.h"
#include <QtGui>
#include <iostream>

QVariant SegTestProxyModel::data(const QModelIndex& index, int role) const {
  if (index.isValid() && role == Qt::DisplayRole) {
    QModelIndex sourceIndex = mapToSource(index);
    return QString("#%1: %2").arg(sourceIndex.row()+1).arg(static_cast<SegmentsTestModel*>(sourceModel())->getReactionLabel(sourceIndex));
  }
  return QVariant();
}


bool SegTestProxyModel::filterAcceptsRow(int source_row, const QModelIndex &source_parent) const {
  if(QSortFilterProxyModel::filterAcceptsRow(source_row,source_parent)) {
    SegmentsTestModel* model = static_cast<SegmentsTestModel*>(sourceModel());
    QModelIndex source_index = model->index(source_row, 9, source_parent);
    int dataType = model->data(source_index,Qt::EditRole).toInt();
    if(dataType!=3) return true;
    return false;
  }
  return false;
}

QVariant SegDataProxyModel::data(const QModelIndex& index, int role) const {
  if (index.isValid() && role == Qt::DisplayRole) {
    QModelIndex sourceIndex = mapToSource(index);
    return QString("#%1: %2").arg(sourceIndex.row()+1).arg(static_cast<SegmentsDataModel*>(sourceModel())->getReactionLabel(sourceIndex));
  } 
  return QVariant();
}



PlotTab::PlotTab(Config& config, SegmentsDataModel* dataModel, SegmentsTestModel* testModel, QWidget* parent) :  
  configure(config), QWidget(parent)  {
  azurePlot = new AZUREPlot(this,this);

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

  QVBoxLayout* rightLayout = new QVBoxLayout;

  QGridLayout *topLayout = new QGridLayout;

  QGroupBox *xAxisBox = new QGroupBox(tr("X-Axis Configuration"));
  xAxisTypeCombo=new QComboBox;
  xAxisTypeCombo->addItem(tr("CoM Energy"));
  xAxisTypeCombo->addItem(tr("Excitation Energy"));
  xAxisTypeCombo->addItem(tr("CoM Angle"));
  connect(xAxisTypeCombo,SIGNAL(activated(int)),this,SLOT(xAxisTypeChanged()));
  xAxisTypeCombo->setCurrentIndex(0);
  azurePlot->setXAxisType(0);
#ifdef MACX_SPACING
  xAxisIsLogCheck = new QCheckBox(tr("Use Log Scale"));
#else
  xAxisIsLogCheck = new QCheckBox(tr("Use Logarithmic Scale"));
#endif
  connect(xAxisIsLogCheck,SIGNAL(toggled(bool)),this,SLOT(xAxisLogScaleChanged(bool)));
  QHBoxLayout* xAxisLayout = new QHBoxLayout;
  xAxisLayout->setContentsMargins(5,5,5,5);
  xAxisLayout->addWidget(xAxisTypeCombo);
  xAxisLayout->addWidget(xAxisIsLogCheck);
  xAxisBox->setLayout(xAxisLayout);
  QGroupBox *yAxisBox = new QGroupBox(tr("Y-Axis Configuration"));
  yAxisXSButton = new QRadioButton(tr("Cross Section"));
  connect(yAxisXSButton,SIGNAL(toggled(bool)),this,SLOT(yAxisTypeChanged()));  
  yAxisXSButton->setChecked(true);
  yAxisSFButton = new QRadioButton(tr("S-Factor"));
  connect(yAxisSFButton,SIGNAL(toggled(bool)),this,SLOT(yAxisTypeChanged()));
#ifdef MACX_SPACING
  yAxisIsLogCheck = new QCheckBox(tr("Use Log Scale"));
#else
  yAxisIsLogCheck = new QCheckBox(tr("Use Logarithmic Scale"));
#endif
  connect(yAxisIsLogCheck,SIGNAL(toggled(bool)),this,SLOT(yAxisLogScaleChanged(bool)));
  yAxisIsLogCheck->setChecked(true);
  QHBoxLayout* yAxisLayout = new QHBoxLayout;
  yAxisLayout->setContentsMargins(5,5,5,5);
  yAxisLayout->addWidget(yAxisXSButton);
  yAxisLayout->addWidget(yAxisSFButton);
  yAxisLayout->addWidget(yAxisIsLogCheck);
  yAxisBox->setLayout(yAxisLayout);

  topLayout->addWidget(xAxisBox,0,0);
  topLayout->addWidget(yAxisBox,0,1);

  rightLayout->addLayout(topLayout);
  rightLayout->addWidget(azurePlot);
  
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

  QGroupBox *dataSegmentSelectorBox = new QGroupBox(tr("Segments From Data"));
  QGridLayout *dataSegmentSelectorLayout = new QGridLayout;
  dataSegmentSelectorLayout->setContentsMargins(5,5,5,5);
  QGroupBox *testSegmentSelectorBox = new QGroupBox(tr("Segments Without Data"));
  QGridLayout *testSegmentSelectorLayout = new QGridLayout;
  testSegmentSelectorLayout->setContentsMargins(5,5,5,5);

  dataSegmentSelectorLayout->addWidget(dataSegmentSelectorList,0,0);
  testSegmentSelectorLayout->addWidget(testSegmentSelectorList,0,0);
  dataSegmentSelectorBox->setLayout(dataSegmentSelectorLayout);
  testSegmentSelectorBox->setLayout(testSegmentSelectorLayout);

  refreshButton = new QPushButton(tr("&Draw"));
  connect(refreshButton,SIGNAL(clicked()),this,SLOT(draw()));
  exportButton = new QPushButton(tr("Export..."));
  connect(exportButton,SIGNAL(clicked()),azurePlot,SLOT(exportPlot()));
  printButton = new QPushButton(tr("Print..."));
  connect(printButton,SIGNAL(clicked()),azurePlot,SLOT(print()));
  QHBoxLayout *buttonLayout = new QHBoxLayout; 
  buttonLayout->addWidget(refreshButton);
  buttonLayout->addWidget(exportButton);
  buttonLayout->addWidget(printButton);

  QVBoxLayout* leftLayout = new QVBoxLayout;
  leftLayout->addWidget(dataSegmentSelectorBox);
  leftLayout->addWidget(testSegmentSelectorBox);
  leftLayout->addLayout(buttonLayout);

  QGridLayout* mainLayout = new QGridLayout;
  mainLayout->addLayout(leftLayout,0,0);
  mainLayout->addLayout(rightLayout,0,1);
  mainLayout->setColumnStretch(1,1);

  setLayout(mainLayout);
}

QList<PlotEntry*> PlotTab::getDataSegments() {
  QList<PlotEntry*> dataSegmentPlotEntries;
  QModelIndexList indexes = dataSegmentSelectorList->selectionModel()->selectedIndexes();
  for(int i = 0; i< indexes.size(); i++) {
    int sourceRow = segDataProxyModel->mapToSource(indexes[i]).row();
    QModelIndex sourceIndex = 
      segDataProxyModel->mapToSource(segDataProxyModel->index(indexes[i].row(),1,QModelIndex()));
    int entranceKey = segDataProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toInt();
    sourceIndex = segDataProxyModel->mapToSource(segDataProxyModel->index(indexes[i].row(),2,QModelIndex()));
    int exitKey = segDataProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toInt();
    sourceIndex = segDataProxyModel->mapToSource(segDataProxyModel->index(indexes[i].row(),7,QModelIndex()));
    int dataType = segDataProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toInt();
    QString filename = (dataType==3) ? 
      QString::fromStdString(configure.outputdir)+QString("AZUREOut_aa=%1_TOTAL_CAPTURE.out").arg(entranceKey) :
      QString::fromStdString(configure.outputdir)+QString("AZUREOut_aa=%1_R=%2.out").arg(entranceKey).arg(exitKey);
    int numPreviousInBlock = 0;
    for(int j =0; j<indexes[i].row(); j++) {
      sourceIndex = segDataProxyModel->mapToSource(segDataProxyModel->index(j,1,QModelIndex()));
      int previousEntranceKey = segDataProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toInt();
      sourceIndex = segDataProxyModel->mapToSource(segDataProxyModel->index(j,2,QModelIndex()));
      int previousExitKey = segDataProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toInt();
      if(previousEntranceKey==entranceKey&&previousExitKey==exitKey) numPreviousInBlock++;
    }
    PlotEntry* newPlotEntry = new PlotEntry(0,entranceKey,exitKey,numPreviousInBlock,filename);
    dataSegmentPlotEntries.push_back(newPlotEntry);
  }
  return dataSegmentPlotEntries;
}

QList<PlotEntry*> PlotTab::getTestSegments() {
  QList<PlotEntry*> testSegmentPlotEntries;
  QModelIndexList indexes = testSegmentSelectorList->selectionModel()->selectedIndexes();
  for(int i = 0; i< indexes.size(); i++) {
    int sourceRow = segTestProxyModel->mapToSource(indexes[i]).row();
    QModelIndex sourceIndex = 
      segTestProxyModel->mapToSource(segTestProxyModel->index(indexes[i].row(),1,QModelIndex()));
    int entranceKey = segTestProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toInt();
    sourceIndex = segTestProxyModel->mapToSource(segTestProxyModel->index(indexes[i].row(),2,QModelIndex()));
    int exitKey = segTestProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toInt();
    sourceIndex = segTestProxyModel->mapToSource(segTestProxyModel->index(indexes[i].row(),9,QModelIndex()));
    int dataType = segTestProxyModel->sourceModel()->data(sourceIndex,Qt::EditRole).toInt();
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
  QList<PlotEntry*> entries = getDataSegments();
  entries.append(getTestSegments());
  azurePlot->draw(entries);
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

void PlotTab::reset() {
  azurePlot->clearEntries();
  xAxisTypeCombo->setCurrentIndex(0);
  xAxisIsLogCheck->setChecked(false);
  yAxisXSButton->setChecked(true);
  yAxisIsLogCheck->setChecked(true);
}

void PlotTab::showInfo(int which,QString title) {
  if(which<infoText.size()) {
    if(!infoDialog[which]) {
      infoDialog[which] = new InfoDialog(infoText[which],this,title);
      infoDialog[which]->setAttribute(Qt:: WA_DeleteOnClose);
      infoDialog[which]->show();
    } else infoDialog[which]->raise();
  }
}
