#include <QtGui>
#include "LevelsTab.h"
#include "LevelsHeaderView.h"
#include "AddLevelDialog.h"
#include "RichTextDelegate.h"
#include "InfoDialog.h"
#include <iostream>

LevelsTab::LevelsTab(QWidget *parent) : QWidget(parent) {
  levelsModel=new LevelsModel(this);
  levelsModelProxy = new QSortFilterProxyModel(this);
  levelsModelProxy->setSourceModel(levelsModel);
  levelsModelProxy->setDynamicSortFilter(true);
  levelsView=new QTableView;
  levelsView->setHorizontalHeader(new LevelsHeaderView(Qt::Horizontal, levelsView));
  levelsView->setModel(levelsModelProxy);
  levelsView->horizontalHeader()->setSortIndicator(4,Qt::AscendingOrder);
  levelsView->setSortingEnabled(true);
  levelsView->verticalHeader()->hide();
  levelsView->horizontalHeader()->setHighlightSections(false);
  levelsView->setColumnWidth(0,60);
  levelsView->setColumnWidth(1,40);
  levelsView->horizontalHeader()->setResizeMode(0,QHeaderView::Fixed);
  levelsView->horizontalHeader()->setResizeMode(1,QHeaderView::Fixed);
  levelsView->horizontalHeader()->setResizeMode(2,QHeaderView::Stretch);
  levelsView->horizontalHeader()->setResizeMode(4,QHeaderView::Stretch);
  levelsView->setColumnHidden(3,true);
  levelsView->setSelectionBehavior(QAbstractItemView::SelectRows);
  levelsView->setSelectionMode(QAbstractItemView::SingleSelection);
  levelsView->setEditTriggers(QAbstractItemView::NoEditTriggers);
  levelsView->setShowGrid(false);
  connect(levelsView->selectionModel(),SIGNAL(selectionChanged(QItemSelection,QItemSelection)),this,SLOT(updateButtons(QItemSelection)));
  connect(levelsView->selectionModel(),SIGNAL(selectionChanged(QItemSelection,QItemSelection)),this,SLOT(updateFilter(QItemSelection)));
  connect(levelsView,SIGNAL(doubleClicked(QModelIndex)),this,SLOT(editLevel()));

  maxLSpin = new QSpinBox;
  maxLSpin->setMinimum(0);
  maxLSpin->setMaximum(10);
  maxLSpin->setSingleStep(1);
  maxLSpin->setValue(2);
  QLabel *maxLLabel = new QLabel(tr("Maximum Orbital Momentum"));
  maxMultSpin = new QSpinBox;
  maxMultSpin->setMinimum(1);
  maxMultSpin->setMaximum(10);
  maxMultSpin->setSingleStep(1);
  maxMultSpin->setValue(2);
  QLabel *maxMultLabel = new QLabel(tr("Maximum Gamma Multipolarity"));  
  maxNumMultSpin = new QSpinBox;
  maxNumMultSpin->setMinimum(1);
  maxNumMultSpin->setMaximum(10);
  maxNumMultSpin->setSingleStep(1);
  maxNumMultSpin->setValue(2);  
  QLabel *maxNumMultLabel = new QLabel(tr("Maximum Gamma Multipolarities\nPer Decay"));
  connect(maxLSpin,SIGNAL(valueChanged(int)),this,SLOT(updateChannelsPairAddedEdited()));
  connect(maxMultSpin,SIGNAL(valueChanged(int)),this,SLOT(updateChannelsPairAddedEdited()));
  connect(maxNumMultSpin,SIGNAL(valueChanged(int)),this,SLOT(updateChannelsPairAddedEdited()));

  channelsModel = new ChannelsModel(this);
  proxyModel = new QSortFilterProxyModel(this);
  proxyModel->setSourceModel(channelsModel);
  proxyModel->setDynamicSortFilter(true);
  proxyModel->setFilterKeyColumn(1);
  proxyModel->setFilterRegExp("-1");
  proxyModel->sort(1,Qt::AscendingOrder);
  channelsView = new QTableView;
  channelsView->setModel(proxyModel);
  channelsView->verticalHeader()->hide();
  channelsView->horizontalHeader()->setHighlightSections(false);
  channelsView->setColumnWidth(0,40);
  channelsView->setColumnWidth(2,160);
  channelsView->horizontalHeader()->setResizeMode(0,QHeaderView::Fixed);
  channelsView->horizontalHeader()->setResizeMode(2,QHeaderView::Fixed);
  channelsView->horizontalHeader()->setResizeMode(3,QHeaderView::Stretch);
  channelsView->horizontalHeader()->setResizeMode(4,QHeaderView::Stretch);
  channelsView->setColumnHidden(1,true);
  channelsView->setColumnHidden(5,true);
  channelsView->setColumnHidden(6,true);
  channelsView->setItemDelegateForColumn(2,new RichTextDelegate());
  channelsView->setItemDelegateForColumn(3,new RichTextDelegate());
  channelsView->setItemDelegateForColumn(4,new RichTextDelegate());
  channelsView->setSelectionBehavior(QAbstractItemView::SelectRows);
  channelsView->setSelectionMode(QAbstractItemView::SingleSelection);
  channelsView->setEditTriggers(QAbstractItemView::NoEditTriggers);
  channelsView->setShowGrid(false);
  connect(channelsView->selectionModel(),SIGNAL(selectionChanged(QItemSelection,QItemSelection)),this,SLOT(updateDetails(QItemSelection)));

  channelDetails=new ChannelDetails(this);
  channelDetails->hide();
  connect(channelDetails->reducedWidthText,SIGNAL(textEdited(const QString&)),this,SLOT(updateReducedWidth(const QString&)));

  addLevelButton = new QPushButton(tr("+"));
  addLevelButton->setMaximumSize(28,28);
  connect(addLevelButton,SIGNAL(clicked()),this,SLOT(addLevel()));
  removeLevelButton = new QPushButton(tr("-"));
  removeLevelButton->setMaximumSize(28,28);
  removeLevelButton->setEnabled(false);
  connect(removeLevelButton,SIGNAL(clicked()),this,SLOT(removeLevel()));

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
  buttonBox->addWidget(addLevelButton,0,0);
  buttonBox->addWidget(removeLevelButton,0,1);
  buttonBox->addItem(new QSpacerItem(28,28),0,2);
  //buttonBox->addWidget(infoButton[0],0,3);
  buttonBox->setColumnStretch(0,0);
  buttonBox->setColumnStretch(1,0);
  buttonBox->setColumnStretch(2,1);
  buttonBox->setColumnStretch(3,0);
#ifdef MACX_SPACING
  buttonBox->setHorizontalSpacing(11);
#else 
  buttonBox->setHorizontalSpacing(0);
#endif

  QGroupBox *levelsBox=new QGroupBox(tr("Compound Nucleus Levels"));
  QGridLayout *levelsLayout=new QGridLayout;
  levelsLayout->setContentsMargins(6,6,6,6);
  levelsLayout->addWidget(levelsView,0,0);
  levelsLayout->addLayout(buttonBox,1,0);
#ifdef MACX_SPACING
  levelsLayout->setVerticalSpacing(0);
#endif
  levelsBox->setLayout(levelsLayout);  

  QGroupBox *configBox = new QGroupBox(tr("Channel Configuration"));
  QGridLayout *configLayout = new QGridLayout;
  configLayout->addWidget(maxLSpin,0,0);
  configLayout->addWidget(maxLLabel,0,1);
  configLayout->addWidget(maxMultSpin,1,0);
  configLayout->addWidget(maxMultLabel,1,1);
  configLayout->addWidget(maxNumMultSpin,2,0);
  configLayout->addWidget(maxNumMultLabel,2,1);
  configLayout->setColumnStretch(0,0);
  configLayout->setColumnStretch(1,1);
  configBox->setLayout(configLayout);
 
  QGroupBox *channelsBox=new QGroupBox(tr("Channels In Selected Level"));
  QGridLayout *channelsLayout=new QGridLayout;
  channelsLayout->setContentsMargins(6,6,6,6);
  channelsLayout->addWidget(channelsView,0,0);
#ifdef MACX_SPACING
  channelsLayout->addItem(new QSpacerItem(40,40),1,0);
  channelsLayout->setVerticalSpacing(0);
#else
  channelsLayout->addItem(new QSpacerItem(34,34),1,0);
#endif
  channelsBox->setLayout(channelsLayout);

  QGridLayout *detailsBox = new QGridLayout;
  QLabel *detailsLabel = new QLabel(tr("Channel Details (select from list to view):"));
  detailsLabel->setAlignment(Qt::AlignHCenter);
  detailsBox->addWidget(detailsLabel,0,0);
  detailsBox->addWidget(channelDetails,1,0);
  detailsBox->setRowStretch(0,0);
  detailsBox->setRowStretch(1,1);  

  QGridLayout *mainLayout = new QGridLayout;
  mainLayout->addWidget(levelsBox,0,0,2,1);
  mainLayout->addWidget(channelsBox,0,1,2,1);
  mainLayout->addWidget(configBox,0,2,1,1);
  mainLayout->addLayout(detailsBox,1,2,1,1);
  mainLayout->setColumnStretch(0,1);
  mainLayout->setColumnStretch(1,1);
  mainLayout->setColumnStretch(2,1);  

  setLayout(mainLayout);
}

void LevelsTab::setPairsModel(PairsModel *model) {
  pairsModel=model;
  channelsModel->setPairsModel(model);
}

void LevelsTab::addLevel() {
  AddLevelDialog aDialog;
  if(aDialog.exec()) {
    LevelsData newLevel;
    newLevel.isActive=1;
    newLevel.isFixed=0;
    newLevel.jValue=(aDialog.jValueText->text()).toDouble();
    if(aDialog.piValueCombo->currentIndex()==0) newLevel.piValue=-1;
    else newLevel.piValue=1;
    newLevel.energy=(aDialog.energyText->text()).toDouble();
    addLevel(newLevel,false);
  }  
}

void LevelsTab::addLevel(LevelsData level, bool fromFile) {
  QList<LevelsData> levels = levelsModel->getLevels();
  if(levelsModel->isLevel(level)==-1) {
    levelsModel->insertRows(levels.size(),1,QModelIndex());
    QModelIndex index = levelsModel->index(levels.size(),0,QModelIndex());
    levelsModel->setData(index,level.isActive,Qt::EditRole);
    index = levelsModel->index(levels.size(),1,QModelIndex());
    levelsModel->setData(index,level.isFixed,Qt::EditRole);
    index = levelsModel->index(levels.size(),2,QModelIndex());
    levelsModel->setData(index,level.jValue,Qt::EditRole);
    index = levelsModel->index(levels.size(),3,QModelIndex());
    levelsModel->setData(index,level.piValue,Qt::EditRole);
    index = levelsModel->index(levels.size(),4,QModelIndex());
    levelsModel->setData(index,level.energy,Qt::EditRole);
    levelsView->resizeRowsToContents();
    if(!fromFile) updateChannelsLevelAdded(levels.size());
  } else {
    QMessageBox::information(this,tr("Duplicate Level"),tr("This level already exists."));
  }
}

void LevelsTab::removeLevel() {
  QItemSelectionModel *selectionModel = levelsView->selectionModel();
  QModelIndexList indexes = selectionModel->selectedRows();
  QModelIndex index=levelsModelProxy->mapToSource(indexes[0]);
  
  levelsModel->removeRows(index.row(),1,QModelIndex());
  updateChannelsLevelDeleted(index.row());
  selectionModel->clearSelection();
}

void LevelsTab::editLevel() {
  QItemSelectionModel *selectionModel = levelsView->selectionModel();
  QModelIndexList indexes = selectionModel->selectedRows();
  QModelIndex index=levelsModelProxy->mapToSource(indexes[0]);
  
  QModelIndex i=levelsModel->index(index.row(),2,QModelIndex());
  QVariant var=levelsModel->data(i,Qt::EditRole);
  QString jValue=var.toString();
  i=levelsModel->index(index.row(),3,QModelIndex());
  var=levelsModel->data(i,Qt::EditRole);
  int piValue=var.toInt();
  i=levelsModel->index(index.row(),4,QModelIndex());
  var=levelsModel->data(i,Qt::EditRole);
  QString energy=var.toString();

  AddLevelDialog aDialog;
  aDialog.setWindowTitle(tr("Edit a Level"));
  aDialog.jValueText->setText(jValue);
  if(piValue==-1) aDialog.piValueCombo->setCurrentIndex(0);
  else aDialog.piValueCombo->setCurrentIndex(1);
  aDialog.energyText->setText(energy);
  
  if(aDialog.exec()) {
    QString newJValue=aDialog.jValueText->text();
    if(newJValue!=jValue) {
      i=levelsModel->index(index.row(),2,QModelIndex());
      levelsModel->setData(i,newJValue,Qt::EditRole);
    }
    int newPiValue;
    if(aDialog.piValueCombo->currentIndex()==0) newPiValue=-1;
    else newPiValue =1;
    if(newPiValue!=piValue) {
      i=levelsModel->index(index.row(),3,QModelIndex());
      levelsModel->setData(i,newPiValue,Qt::EditRole);
    }
    QString newEnergy=aDialog.energyText->text();
    if(newEnergy!=energy) {
      i=levelsModel->index(index.row(),4,QModelIndex());
      levelsModel->setData(i,newEnergy,Qt::EditRole);
    }
    updateChannelsLevelEdited(index.row());
  }
}

void LevelsTab::updateButtons(const QItemSelection &selection) {
  QModelIndexList indexes=selection.indexes();
  
  if(indexes.isEmpty()) {
    removeLevelButton->setEnabled(false);
  } else {
    removeLevelButton->setEnabled(true);
  }
}

void LevelsTab::updateChannelsLevelAdded(int levelIndex) {
  QList<ChannelsData> newChannels = calculateChannels(levelIndex);
  channelsModel->insertRows(0,newChannels.size(),QModelIndex());
  for(int i=0;i<newChannels.size();i++) {
    QModelIndex index = channelsModel->index(i,0,QModelIndex());
    channelsModel->setData(index,newChannels.at(i).isFixed,Qt::EditRole);
    index = channelsModel->index(i,1,QModelIndex());
    channelsModel->setData(index,newChannels.at(i).levelIndex,Qt::EditRole);
    index = channelsModel->index(i,2,QModelIndex());
    channelsModel->setData(index,newChannels.at(i).pairIndex,Qt::EditRole);
    index = channelsModel->index(i,3,QModelIndex());
    channelsModel->setData(index,newChannels.at(i).sValue,Qt::EditRole);
    index = channelsModel->index(i,4,QModelIndex());
    channelsModel->setData(index,newChannels.at(i).lValue,Qt::EditRole);
    index = channelsModel->index(i,5,QModelIndex());
    channelsModel->setData(index,newChannels.at(i).radType,Qt::EditRole);
    index = channelsModel->index(i,6,QModelIndex());
    channelsModel->setData(index,newChannels.at(i).reducedWidth,Qt::EditRole);
  }
  channelsView->resizeRowsToContents();
}

void LevelsTab::updateChannelsLevelDeleted(int levelIndex) {
  QList<ChannelsData> channels=channelsModel->getChannels();
  int deleted=0;
  for(int i=0;i<channels.size();i++) {
    if(channels.at(i).levelIndex==levelIndex) {
      channelsModel->removeRows(i-deleted,1,QModelIndex());
      deleted++;
    }
  }
  channels=channelsModel->getChannels();
  for(int i=0;i<channels.size();i++) {
    if(channels.at(i).levelIndex>levelIndex) {
      QModelIndex index=channelsModel->index(i,1,QModelIndex());
      channelsModel->setData(index,channels.at(i).levelIndex-1,Qt::EditRole);
    }
  }
}

void LevelsTab::updateChannelsLevelEdited(int levelIndex) {
  QList<ChannelsData> channels=channelsModel->getChannels();
  int deleted=0;
  for(int i=0;i<channels.size();i++) {
    if(channels.at(i).levelIndex==levelIndex) {
      channelsModel->removeRows(i-deleted,1,QModelIndex());
      deleted++;
    }
  }
  QList<ChannelsData> newChannels = calculateChannels(levelIndex);
  channelsModel->insertRows(0,newChannels.size(),QModelIndex());
  for(int i=0;i<newChannels.size();i++) {
    QModelIndex index = channelsModel->index(i,1,QModelIndex());
    channelsModel->setData(index,newChannels.at(i).levelIndex,Qt::EditRole);
    index = channelsModel->index(i,2,QModelIndex());
    channelsModel->setData(index,newChannels.at(i).pairIndex,Qt::EditRole);
    index = channelsModel->index(i,3,QModelIndex());
    channelsModel->setData(index,newChannels.at(i).sValue,Qt::EditRole);
    index = channelsModel->index(i,4,QModelIndex());
    channelsModel->setData(index,newChannels.at(i).lValue,Qt::EditRole);
    index = channelsModel->index(i,5,QModelIndex());
    channelsModel->setData(index,newChannels.at(i).radType,Qt::EditRole);
    for(int ii=0;ii<channels.size();ii++) {
      if(channels.at(ii).levelIndex==newChannels.at(i).levelIndex&&
	 channels.at(ii).pairIndex==newChannels.at(i).pairIndex&&
	 channels.at(ii).sValue==newChannels.at(i).sValue&&
	 channels.at(ii).lValue==newChannels.at(i).lValue&&
	 channels.at(ii).radType==newChannels.at(i).radType) {
	newChannels[i].isFixed=channels.at(ii).isFixed;
	newChannels[i].reducedWidth=channels.at(ii).reducedWidth;
	break;
      }
    }
    index = channelsModel->index(i,6,QModelIndex());
    channelsModel->setData(index,newChannels.at(i).reducedWidth,Qt::EditRole);
    index = channelsModel->index(i,0,QModelIndex());
    channelsModel->setData(index,newChannels.at(i).isFixed,Qt::EditRole);
  }
  channelsView->resizeRowsToContents();
}

QList<ChannelsData> LevelsTab::calculateChannels(int levelIndex) {
  QList<ChannelsData> channels;
  QList<PairsData> pairs = pairsModel->getPairs();
  LevelsData level = (levelsModel->getLevels()).at(levelIndex);

  int maxL=maxLSpin->value();
  int maxMult=maxMultSpin->value();
  int maxNumMult=maxNumMultSpin->value();

  for(int i=0;i<pairs.size();i++) {
    PairsData pair=pairs.at(i);
    if(pair.pairType==0) {
      for(double s=fabs(pair.heavyJ-pair.lightJ);s<=pair.heavyJ+pair.lightJ;s+=1.0) {
	for(double l=fabs(s-level.jValue);l<=s+level.jValue;l+=1.0) {
	  if(int(l*2.0)%2==0&&pair.lightPi*pair.heavyPi*pow(-1,int(l))==level.piValue&&int(l)<=maxL) {
	    ChannelsData channel={0,levelIndex,i,s,int(l),'P',0.0};
	    channels.push_back(channel);
	  }
	}
      }
    } else if(pair.pairType==20) {
      if(fabs(pair.heavyJ-level.jValue)==0.&&pair.heavyPi==level.piValue) {
	ChannelsData gtChannel = {0,levelIndex,i,0.,1,'G',0.0};
	ChannelsData fChannel = {0,levelIndex,i,0.,0,'F',0.0}; 
	channels.push_back(gtChannel);
	channels.push_back(fChannel);
      } else if(fabs(pair.heavyJ-level.jValue)==1.&&pair.heavyPi==level.piValue) {
	ChannelsData gtChannel = {0,levelIndex,i,1.,1,'G',0.0};
	channels.push_back(gtChannel);
      }
    } else {
      int numMult=1;
      for(int l=1;l<=maxMult;l++) {
	if(fabs(l-pair.heavyJ)<=level.jValue&&level.jValue<=l+pair.heavyJ&&numMult<=maxNumMult) {
	  QChar radType;
	  int parityChange=pair.heavyPi*level.piValue;
	  if(l%2!=0) {
	    if(parityChange==-1) radType='E';
	    else radType='M';
	  } else {
	    if(parityChange==-1) radType='M';
	    else radType='E';
	  }
	  ChannelsData channel={0,levelIndex,i,pair.heavyJ,l,radType,0.0};
	  channels.push_back(channel);
	  numMult++;
	}
      }
    }
  }
  return channels;
}

void LevelsTab::updateFilter(const QItemSelection &selection) {
  QModelIndexList indexes = selection.indexes();
  if(indexes.isEmpty())
    proxyModel->setFilterRegExp("-1");
  else {
    QModelIndex index = levelsModelProxy->mapToSource(indexes.at(0));
    int row=index.row();
    proxyModel->setFilterRegExp(QString("\\b%1\\b").arg(row));
  } 
  channelsView->resizeRowsToContents();
}

void LevelsTab::updateChannelsPairAddedEdited() {
  QList<LevelsData> levels=levelsModel->getLevels();
  QList<ChannelsData> channels=channelsModel->getChannels();
  channelsModel->removeRows(0,channels.size(),QModelIndex());
  for(int levelIndex=0;levelIndex<levels.size();levelIndex++) {
    QList<ChannelsData> newChannels=calculateChannels(levelIndex);
    channelsModel->insertRows(0,newChannels.size(),QModelIndex());
    for(int i=0;i<newChannels.size();i++) {
      QModelIndex index = channelsModel->index(i,1,QModelIndex());
      channelsModel->setData(index,newChannels.at(i).levelIndex,Qt::EditRole);
      index = channelsModel->index(i,2,QModelIndex());
      channelsModel->setData(index,newChannels.at(i).pairIndex,Qt::EditRole);
      index = channelsModel->index(i,3,QModelIndex());
      channelsModel->setData(index,newChannels.at(i).sValue,Qt::EditRole);
      index = channelsModel->index(i,4,QModelIndex());
      channelsModel->setData(index,newChannels.at(i).lValue,Qt::EditRole);
      index = channelsModel->index(i,5,QModelIndex());
      channelsModel->setData(index,newChannels.at(i).radType,Qt::EditRole);
      for(int ii=0;ii<channels.size();ii++) {
	if(channels.at(ii).levelIndex==newChannels.at(i).levelIndex&&
	   channels.at(ii).pairIndex==newChannels.at(i).pairIndex&&
	   channels.at(ii).sValue==newChannels.at(i).sValue&&
	   channels.at(ii).lValue==newChannels.at(i).lValue&&
	   channels.at(ii).radType==newChannels.at(i).radType) {
	  newChannels[i].reducedWidth=channels.at(ii).reducedWidth;
	  newChannels[i].isFixed=channels.at(ii).isFixed;
	  break;
	}
      }
      index = channelsModel->index(i,6,QModelIndex());
      channelsModel->setData(index,newChannels.at(i).reducedWidth,Qt::EditRole);
      index = channelsModel->index(i,0,QModelIndex());
      channelsModel->setData(index,newChannels.at(i).isFixed,Qt::EditRole);
    }
  }
  channelsView->resizeRowsToContents();
}

void LevelsTab::updateChannelsPairRemoved(int pairIndex) {
  QList<ChannelsData> channels=channelsModel->getChannels();
  int deleted=0;
  for(int i=0;i<channels.size();i++) {
    if(channels.at(i).pairIndex==pairIndex) {
      channelsModel->removeRows(i-deleted,1,QModelIndex());
      deleted++;
    }
  }
  channels=channelsModel->getChannels();
  for(int i=0;i<channels.size();i++) {
    if(channels.at(i).pairIndex>pairIndex)  {
      QModelIndex index=channelsModel->index(i,2,QModelIndex());
      channelsModel->setData(index,channels.at(i).pairIndex-1,Qt::EditRole);
    }
  }
  channelsView->resizeRowsToContents();
}

void LevelsTab::updateDetails(const QItemSelection &selection) {
  QModelIndexList indexes=selection.indexes();
  
  if(indexes.isEmpty()) {
    channelDetails->hide();
  } else {
    QModelIndex index=proxyModel->mapToSource(indexes.at(0));

    QModelIndex i=channelsModel->index(index.row(),1,QModelIndex());
    QVariant var=channelsModel->data(i,Qt::EditRole);
    int levelIndex = var.toInt();
    i=channelsModel->index(index.row(),2,QModelIndex());
    var=channelsModel->data(i,Qt::EditRole);
    int pairIndex = var.toInt();

    ChannelsData channel=channelsModel->getChannels().at(index.row());
    LevelsData level=levelsModel->getLevels().at(levelIndex);
    PairsData pair=pairsModel->getPairs().at(pairIndex);

    QString details="";
    QTextStream stm(&details,QIODevice::Append);
    stm << QString("%1 MeV level with spin %2\n   transitioning via pair key #%3").arg(level.energy).arg(levelsModel->getSpinLabel(level)).arg(pairIndex+1) 
	<< endl;
    if(channel.radType=='P') {
      stm << QString("Channel configuration is\n   s = %1, l = %2").arg(channelsModel->getSpinLabel(channel)).arg(channel.lValue) 
	  << endl << endl;
      stm << qSetFieldWidth(21) << right << "Light Particle Spin: " 
	  << qSetFieldWidth(0) << left << QString("%1").arg(pairsModel->getSpinLabel(pair,0)) << endl;
      stm << qSetFieldWidth(21) << right << "Light Particle Z: " 
	  << qSetFieldWidth(0) << left << QString("%1").arg(pair.lightZ) << endl;
      stm << qSetFieldWidth(21) << right << "Light Particle M: " 
	  << qSetFieldWidth(0) << left << QString("%1").arg(pair.lightM) << endl;
      stm << qSetFieldWidth(21) << right << "Light Particle G: " 
	  << qSetFieldWidth(0) << left << QString("%1").arg(pair.lightG) << endl;
    } else if(channel.radType=='G' || channel.radType=='F') {
      if(channel.radType=='G') 
	stm << QString("Channel is Gamow-Teller beta decay") << endl << endl;
      else 
	stm << QString("Channel is Fermi beta decay") << endl << endl;
      stm << qSetFieldWidth(21) << right << "Fermion Charge: " 
	  << qSetFieldWidth(0) << left << QString("%1").arg(pair.lightZ) << endl;
    } else {
      stm << QString("Capture gamma is %1%2 radiation").arg(channel.radType).arg(channel.lValue) << endl;
      if(((channel.radType=='E'&&channel.lValue==1)&&
	  (pair.ecMultMask & (1<<0)))||
	 ((channel.radType=='M'&&channel.lValue==1)&&
	  (pair.ecMultMask & (1<<1)))||
	 ((channel.radType=='E'&&channel.lValue==2)&&
	  (pair.ecMultMask & (1<<2)))) 
	stm << "Capture is internal and external" << endl;
      else stm << "Capture is internal only" << endl;
      stm << endl;
    }
    stm << qSetFieldWidth(21) << right << "Heavy Particle Spin: "
	<< qSetFieldWidth(0) << left <<QString("%1").arg(pairsModel->getSpinLabel(pair,1)) << endl;
    stm << qSetFieldWidth(21) << right << "Heavy Particle Z: "
	<< qSetFieldWidth(0) << left <<QString("%1").arg(pair.heavyZ) << endl;
    stm << qSetFieldWidth(21) << right << "Heavy Particle M: "
	<< qSetFieldWidth(0) << left <<QString("%1").arg(pair.heavyM) << endl;
    stm << qSetFieldWidth(21) << right << "Heavy Particle G: "
	<< qSetFieldWidth(0) << left <<QString("%1").arg(pair.heavyG) << endl;
    if(channel.radType!='G'&&channel.radType!='F')
      stm << qSetFieldWidth(21) << right << "Excitation Energy: "
	  << qSetFieldWidth(0) << left <<QString("%1").arg(pair.excitationEnergy) << endl;
    if(channel.radType!='M'&&channel.radType!='E') {
      stm << qSetFieldWidth(21) << right << "Separation Energy: "
	  << qSetFieldWidth(0) << left <<QString("%1").arg(pair.seperationEnergy) << endl;
      stm << qSetFieldWidth(21) << right << "Channel Radius: "
	  << qSetFieldWidth(0) << left <<QString("%1").arg(pair.channelRadius) << endl;
    }
    stm.flush();
    channelDetails->details->setText(details);
    if(level.energy<(pair.seperationEnergy+pair.excitationEnergy)&&pair.pairType==0) channelDetails->setNormParam(1);
    else if (pair.pairType==10&&level.energy==pair.excitationEnergy&&level.jValue==pair.heavyJ&&
	     level.piValue==pair.heavyPi&&channel.radType=='M'&&channel.lValue==1) 
      channelDetails->setNormParam(2);
    else if (pair.pairType==10&&level.energy==pair.excitationEnergy&&level.jValue==pair.heavyJ&&
	     level.piValue==pair.heavyPi&&channel.radType=='E'&&channel.lValue==2) 
      channelDetails->setNormParam(3);    
    else if(pair.pairType==20) channelDetails->setNormParam(4); 
    else channelDetails->setNormParam(0);
    channelDetails->reducedWidthText->setText(QString("%1").arg(channel.reducedWidth));
    channelDetails->show();
  }
}

void LevelsTab::updateReducedWidth(const QString &string) {
  QItemSelectionModel *selectionModel=channelsView->selectionModel();
  QModelIndex index=proxyModel->mapToSource(selectionModel->selectedRows().at(0));
  
  if(index.isValid()) {
    double reducedWidth=string.toDouble();
    QModelIndex i = channelsModel->index(index.row(),6,QModelIndex());
    channelsModel->setData(i,reducedWidth,Qt::EditRole);
  }
}

bool LevelsTab::writeNuclearFile(QTextStream& outStream) {
  QList<PairsData> pairs=pairsModel->getPairs();
  QList<LevelsData> levels=levelsModel->getLevels();
  QList<ChannelsData> channels=channelsModel->getChannels();

  outStream.setFieldAlignment(QTextStream::AlignRight);

  double lowJ=0;
  double highJ=0;
  for(int la=0;la<levels.size();la++) {
    double tempJ=levels.at(la).jValue;
    if(la==0) {
      lowJ=tempJ;
      highJ=tempJ;
    } else if(tempJ<=lowJ) lowJ=tempJ;
    else if(tempJ>=highJ) highJ=tempJ;
  }
  QList<QList<LevelsData> > sortedLevels;
  QList<QList<int> > levelsMap;
  for(double j=lowJ;j<=highJ;j+=0.5) {
    for(int pi=-1;pi<=1;pi+=2) {
      for(int la=0;la<levels.size();la++) {
	if(levels.at(la).jValue==j&&levels.at(la).piValue==pi) {
	  if(sortedLevels.size()==0||(sortedLevels.at(sortedLevels.size()-1).at(0).jValue!=j||
				      sortedLevels.at(sortedLevels.size()-1).at(0).piValue!=pi)) {
	    QList<LevelsData> tempLevelList;
	    sortedLevels.append(tempLevelList);
	    QList<int> tempKeyList;
	    levelsMap.append(tempKeyList);
	  }
	  if(sortedLevels.at(sortedLevels.size()-1).size()==0) {
	    sortedLevels[sortedLevels.size()-1].append(levels.at(la));
	    levelsMap[levelsMap.size()-1].append(la);
	  } else for(int mu=0;mu<sortedLevels.at(sortedLevels.size()-1).size();mu++) {
	    if(levels.at(la).energy<=sortedLevels.at(sortedLevels.size()-1).at(mu).energy) {
	      sortedLevels[sortedLevels.size()-1].insert(mu,levels.at(la));
	      levelsMap[levelsMap.size()-1].insert(mu,la);
	      break;
	    }
	    else if(mu==sortedLevels.at(sortedLevels.size()-1).size()-1) {
	      sortedLevels[sortedLevels.size()-1].append(levels.at(la)); 
	      levelsMap[levelsMap.size()-1].append(la);
	      break;
	    }
	  }
	}
      }
    }
  }
  int levelId=1;
  for(int i=0;i<sortedLevels.size();i++){
    for(int ii=0;ii<sortedLevels.at(i).size();ii++) {
      for(int ch=0;ch<channels.size();ch++) {
	if(channels.at(ch).levelIndex==levelsMap.at(i).at(ii)) {
	  outStream << qSetFieldWidth(4)  << sortedLevels.at(i).at(ii).jValue  
		    << qSetFieldWidth(5)  << sortedLevels.at(i).at(ii).piValue 
		    << qSetFieldWidth(13) << qSetRealNumberPrecision(10) << sortedLevels.at(i).at(ii).energy 
		    << qSetFieldWidth(5)  << sortedLevels.at(i).at(ii).isFixed
		    << qSetFieldWidth(5)  << "1" 
		    << qSetFieldWidth(5)  << channels.at(ch).pairIndex+1
		    << qSetFieldWidth(5)  << int(channels.at(ch).sValue*2)
		    << qSetFieldWidth(5)  << int(channels.at(ch).lValue*2)
		    << qSetFieldWidth(5)  << levelId
		    << qSetFieldWidth(5)  << sortedLevels.at(i).at(ii).isActive
		    << qSetFieldWidth(5)  << channels.at(ch).isFixed
		    << qSetFieldWidth(20) << qSetRealNumberPrecision(15) << channels.at(ch).reducedWidth
		    << qSetFieldWidth(5)  << pairs.at(channels.at(ch).pairIndex).lightJ
		    << qSetFieldWidth(5)  << pairs.at(channels.at(ch).pairIndex).lightPi
		    << qSetFieldWidth(5)  << pairs.at(channels.at(ch).pairIndex).heavyJ
		    << qSetFieldWidth(5)  << pairs.at(channels.at(ch).pairIndex).heavyPi
		    << qSetFieldWidth(13) << qSetRealNumberPrecision(10) << pairs.at(channels.at(ch).pairIndex).excitationEnergy
		    << qSetFieldWidth(13) << qSetRealNumberPrecision(10) << pairs.at(channels.at(ch).pairIndex).lightM
		    << qSetFieldWidth(13) << qSetRealNumberPrecision(10) << pairs.at(channels.at(ch).pairIndex).heavyM
		    << qSetFieldWidth(5)  << pairs.at(channels.at(ch).pairIndex).lightZ
		    << qSetFieldWidth(5)  << pairs.at(channels.at(ch).pairIndex).heavyZ
		    << qSetFieldWidth(13) << qSetRealNumberPrecision(10) << pairs.at(0).seperationEnergy
		    << qSetFieldWidth(13) << qSetRealNumberPrecision(10) << pairs.at(channels.at(ch).pairIndex).seperationEnergy
		    << qSetFieldWidth(5)  << pairs.at(channels.at(ch).pairIndex).isUPOS
                    << qSetFieldWidth(5)  << "0" //pairs.at(channels.at(ch).pairIndex).secondaryDecayL
                    << qSetFieldWidth(13) << "0" //qSetRealNumberPrecision(10) << pairs.at(channels.at(ch).pairIndex).Ic
		    << qSetFieldWidth(6)  << pairs.at(channels.at(ch).pairIndex).pairType
		    << qSetFieldWidth(13) << qSetRealNumberPrecision(10) << pairs.at(channels.at(ch).pairIndex).channelRadius 
		    << qSetFieldWidth(13) << qSetRealNumberPrecision(10) << pairs.at(channels.at(ch).pairIndex).lightG 
		    << qSetFieldWidth(13) << qSetRealNumberPrecision(10) << pairs.at(channels.at(ch).pairIndex).heavyG 
		    << qSetFieldWidth(8)  << pairs.at(channels.at(ch).pairIndex).ecMultMask
		    << qSetFieldWidth(0)  << endl;
	}
      }  
      outStream << endl;
      levelId++;
    }
  }

  outStream.setFieldAlignment(QTextStream::AlignLeft);

  return true;
}

bool LevelsTab::readNuclearFile(QTextStream &inStream) {

  double levelJ;
  int levelPi;
  double levelEnergy;
  int levelFix;
  int aa;
  int ir;
  double channelS;
  int channelL;
  int levelId;
  int levelYN;
  int channelFix;
  double channelReducedWidth;
  double lightJ;
  int lightPi;
  double heavyJ;
  int heavyPi;
  double excitationEnergy;
  double lightM;
  double heavyM;
  int lightZ;
  int heavyZ;
  double seperationEnergyIn;
  double seperationEnergyOut;
  int isupos;
  int secondaryDecayL;
  double Ic;
  int pairType;
  double channelRadius;
  double lightG;
  double heavyG;
  int ecMultMask;

  maxLSpin->blockSignals(true);
  maxMultSpin->blockSignals(true);
  maxNumMultSpin->blockSignals(true);
  maxLSpin->setValue(0);
  maxMultSpin->setValue(0);
  maxNumMultSpin->setValue(0);

  int maxLValue=0;
  int maxMultValue=0;
  int maxNumMultValue=0;
  bool firstLine=true;;
  int lastPair;
  int currentPair=0;
  int thisNumMult=0;

  QString line("");
  while(!inStream.atEnd()&&line.trimmed()!=QString("</levels>")) {
    line=inStream.readLine();
    if(line.trimmed().isEmpty()) continue;
    if(!inStream.atEnd()&&line.trimmed()!=QString("</levels>")) {
      QTextStream in(&line);
      in >> levelJ >> levelPi >> levelEnergy >> levelFix >> aa >> ir >> channelS >> channelL >> levelId >> levelYN >> channelFix
	 >> channelReducedWidth >> lightJ >> lightPi >> heavyJ >> heavyPi >> excitationEnergy 
	 >> lightM >> heavyM >> lightZ >> heavyZ >> seperationEnergyIn >> seperationEnergyOut
	 >> isupos >> secondaryDecayL >> Ic >> pairType >> channelRadius >> lightG >> heavyG;
      if(in.status()!=QTextStream::Ok) return false;
      in >> ecMultMask;
      if(in.status()!=QTextStream::Ok) ecMultMask=0;
      if(firstLine) {
	lastPair=ir;
	firstLine=false;
      } else lastPair=currentPair;
      currentPair=ir;
      if(lastPair==currentPair) thisNumMult++;
      else thisNumMult=1;
      if(pairType==0) {
	if(channelL/2>maxLValue) maxLValue=channelL/2;
      } else if (pairType==10) {
	if(channelL/2>maxMultValue) maxMultValue=channelL/2;
	if(thisNumMult>maxNumMultValue) maxNumMultValue=thisNumMult;
      }
      
      PairsData newPair={lightJ,lightPi,lightZ,lightM,lightG,heavyJ,heavyPi,heavyZ,heavyM,
			 heavyG,excitationEnergy,seperationEnergyOut,channelRadius,pairType,ecMultMask,isupos};
      int pairIndex=ir-1;
      if(pairsModel->numPairs()<ir) {
	emit(readNewPair(newPair,pairIndex,true));
      } else if(pairsModel->isPair(newPair) == -1) {
	emit(readExistingPair(newPair,pairIndex,true));
      }
      LevelsData newLevel = {levelYN,levelFix,levelJ,levelPi,levelEnergy};
      int levelIndex = levelsModel->isLevel(newLevel);
      if(levelIndex == -1) {
	addLevel(newLevel,true);
	levelIndex = levelsModel->isLevel(newLevel);
      }
      QChar radType;
      if(pairType==0) radType='P';
      else if(pairType==20) {
	if(channelL==0) radType = 'F';
	else radType = 'G';
      }
      else {
	int parityChange=heavyPi*levelPi;
	if((channelL/2)%2!=0) {
	  if(parityChange==-1) radType='E';
	  else radType='M';
	} else {
	  if(parityChange==-1) radType='M';
	  else radType='E';
	}
      }
      
      int channelIndex=channelsModel->getChannels().size();
      channelsModel->insertRows(channelIndex,1,QModelIndex());
      QModelIndex index = channelsModel->index(channelIndex,0,QModelIndex());
      channelsModel->setData(index,channelFix,Qt::EditRole);
      index = channelsModel->index(channelIndex,1,QModelIndex());
      channelsModel->setData(index,levelIndex,Qt::EditRole);
      index = channelsModel->index(channelIndex,2,QModelIndex());
      channelsModel->setData(index,pairIndex,Qt::EditRole);
      index = channelsModel->index(channelIndex,3,QModelIndex());
      channelsModel->setData(index,channelS/2.0,Qt::EditRole);
      index = channelsModel->index(channelIndex,4,QModelIndex());
      channelsModel->setData(index,channelL/2,Qt::EditRole);
      index = channelsModel->index(channelIndex,5,QModelIndex());
      channelsModel->setData(index,radType,Qt::EditRole);
      index = channelsModel->index(channelIndex,6,QModelIndex());
      channelsModel->setData(index,channelReducedWidth,Qt::EditRole);
      channelsView->resizeRowsToContents();
    }
  }
  if(line.trimmed()!=QString("</levels>")) return false;

  if(maxLValue>maxLSpin->value()) maxLSpin->setValue(maxLValue);
  if(maxMultValue>maxMultSpin->value()) maxMultSpin->setValue(maxMultValue);
  if(maxNumMultValue>maxNumMultSpin->value()) maxNumMultSpin->setValue(maxNumMultValue);
  maxLSpin->blockSignals(false);
  maxMultSpin->blockSignals(false);
  maxNumMultSpin->blockSignals(false);
  levelsView->resizeRowsToContents();

  return true;
}

void LevelsTab::reset() {
  pairsModel->removeRows(0,pairsModel->getPairs().size(),QModelIndex());
  levelsModel->removeRows(0,levelsModel->getLevels().size(),QModelIndex());
  channelsModel->removeRows(0,channelsModel->getChannels().size(),QModelIndex());
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

void LevelsTab::showInfo(int which,QString title) {
  if(which<infoText.size()) {
    if(!infoDialog[which]) {
      infoDialog[which] = new InfoDialog(infoText[which],this,title);
      infoDialog[which]->setAttribute(Qt:: WA_DeleteOnClose);
      infoDialog[which]->show();
    } else infoDialog[which]->raise();
  }
}
