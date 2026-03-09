#include <QtGui>
#include "PairsTab.h"
#include "RichTextDelegate.h"
#include "InfoDialog.h"
#include <iostream>

PairsTab::PairsTab(QWidget *parent) : QWidget(parent) {
  pairsModel = new PairsModel(this);

  pairsView = new QTableView;
  pairsView->setModel(pairsModel);
  pairsView->verticalHeader()->setHighlightSections(false);
  pairsView->horizontalHeader()->setHighlightSections(false);
  pairsView->setColumnHidden(1,true);
  pairsView->setColumnHidden(3,true);
  pairsView->setColumnHidden(6,true);
  pairsView->setColumnHidden(8,true);
  pairsView->setColumnHidden(13,true);
  pairsView->setColumnHidden(14,true);
  pairsView->setColumnHidden(15,true);
//  pairsView->setColumnHidden(16,true);
//  pairsView->setColumnHidden(17,true);
  RichTextDelegate *rt = new RichTextDelegate();
  pairsView->setItemDelegateForColumn(0,rt);
  pairsView->setItemDelegateForColumn(5,rt);
  pairsView->setSelectionBehavior(QAbstractItemView::SelectRows);
  pairsView->setSelectionMode(QAbstractItemView::SingleSelection);
  pairsView->setEditTriggers(QAbstractItemView::NoEditTriggers);
  pairsView->setShowGrid(false);

  connect(pairsView->selectionModel(),SIGNAL(selectionChanged(QItemSelection,QItemSelection)),this,SLOT(updateButtons(QItemSelection)));
  connect(pairsView,SIGNAL(doubleClicked(QModelIndex)),this,SLOT(editPair()));

  for(int i = 0; i<PairsData::SIZE;i++) pairsView->horizontalHeader()->setResizeMode(i,QHeaderView::Stretch);

  addButton=new QPushButton(tr("+"));
  addButton->setMaximumSize(28,28);
  connect(addButton,SIGNAL(clicked()),this,SLOT(addPair()));
  deleteButton = new QPushButton(tr("-"));
  deleteButton->setMaximumSize(28,28);
  deleteButton->setEnabled(false);
  connect(deleteButton,SIGNAL(clicked()),this,SLOT(removePair()));
  
  /*
  mapper = new QSignalMapper(this);
  connect(mapper,SIGNAL(mapped(int)),this,SLOT(showInfo(int)));
  
  infoButton[0] = new QPushButton(this);
  infoButton[0]->setMaximumSize(28,28);
  infoButton[0]->setIcon(style()->standardIcon(QStyle::SP_MessageBoxInformation));
  mapper->setMapping(infoButton[0],0);
  connect(infoButton[0],SIGNAL(clicked()),mapper,SLOT(map()));
  */

  QGridLayout *buttonBox = new QGridLayout;
  buttonBox->addWidget(addButton,0,0);
  buttonBox->addWidget(deleteButton,0,1);
  buttonBox->addItem(new QSpacerItem(28,28),0,2);
  // buttonBox->addWidget(infoButton[0],0,3);
  buttonBox->setColumnStretch(0,0);
  buttonBox->setColumnStretch(1,0);
  buttonBox->setColumnStretch(2,1);
  buttonBox->setColumnStretch(3,0);
#ifdef MACX_SPACING
  buttonBox->setHorizontalSpacing(11);
#else 
  buttonBox->setHorizontalSpacing(0);
#endif

  QVBoxLayout *layout = new QVBoxLayout();
  layout->addWidget(pairsView);
  layout->addLayout(buttonBox);
  setLayout(layout);
}

PairsModel *PairsTab::getPairsModel() {
  return pairsModel;
}

void PairsTab::addPair() {
  AddPairDialog aDialog;
  if(aDialog.exec()) {
    if(!((pairsModel->getPairs()).size()==0&&aDialog.pairTypeCombo->currentIndex()!=0)) {
      PairsData newPair;
      newPair.lightJ=(aDialog.lightJText->text()).toDouble();
      if(aDialog.lightPiCombo->currentIndex() == 0) newPair.lightPi=-1;
      else newPair.lightPi=1;
      newPair.lightZ=(aDialog.lightZText->text()).toInt();
      newPair.lightM=(aDialog.lightMText->text()).toDouble();
      newPair.lightG=(aDialog.lightGText->text()).toDouble();
      newPair.heavyJ=(aDialog.heavyJText->text()).toDouble();
      if(aDialog.heavyPiCombo->currentIndex() == 0) newPair.heavyPi=-1;
      else newPair.heavyPi=1;
      newPair.heavyZ=(aDialog.heavyZText->text()).toInt();
      newPair.heavyM=(aDialog.heavyMText->text()).toDouble();
      newPair.heavyG=(aDialog.heavyGText->text()).toDouble();
      newPair.seperationEnergy=(aDialog.seperationEnergyText->text()).toDouble();
      newPair.excitationEnergy=(aDialog.excitationEnergyText->text()).toDouble();
      newPair.channelRadius=(aDialog.channelRadiusText->text()).toDouble();
      if(aDialog.pairTypeCombo->currentIndex() == 1) newPair.pairType=10;
      else if(aDialog.pairTypeCombo->currentIndex() == 2) newPair.pairType=20;
      else newPair.pairType=0;
      unsigned int newMask=0;
      if(aDialog.e1Check->isChecked()) newMask |= (1<<0);
      if(aDialog.m1Check->isChecked()) newMask |= (1<<1);
      if(aDialog.e2Check->isChecked()) newMask |= (1<<2);
      newPair.ecMultMask=newMask;
//      if(aDialog.uposCheck->isChecked()) newPair.isUPOS=1;
//      else newPair.isUPOS=0;
      newPair.isUPOS=1;
//      newPair.secondaryDecayL=(aDialog.secondaryDecayLText->text()).toInt();
//      newPair.Ic=(aDialog.icText->text()).toDouble();
      addPair(newPair,pairsModel->numPairs(),false);
    } else {
      QMessageBox::information(this, tr("Pair Type Error"),
			       tr("The first pair must be a particle,particle pair."));
    }
  }
}

void PairsTab::addPair(PairsData pair,int pairIndex,bool fromFile) {
  if(pairsModel->isPair(pair)==-1) {
    pairsModel->insertRows(pairsModel->numPairs(),pairIndex+1-pairsModel->numPairs(),QModelIndex());

    QModelIndex index = pairsModel->index(pairIndex,0,QModelIndex());
    pairsModel->setData(index,pair.lightJ,Qt::EditRole);
    index = pairsModel->index(pairIndex,1,QModelIndex());
    pairsModel->setData(index,pair.lightPi,Qt::EditRole);
    index = pairsModel->index(pairIndex,2,QModelIndex());
    pairsModel->setData(index,pair.lightZ,Qt::EditRole);
    index = pairsModel->index(pairIndex,3,QModelIndex());
    pairsModel->setData(index,pair.lightM,Qt::EditRole);
    index = pairsModel->index(pairIndex,4,QModelIndex());
    pairsModel->setData(index,pair.lightG,Qt::EditRole);
    index = pairsModel->index(pairIndex,5,QModelIndex());
    pairsModel->setData(index,pair.heavyJ,Qt::EditRole);
    index = pairsModel->index(pairIndex,6,QModelIndex());
    pairsModel->setData(index,pair.heavyPi,Qt::EditRole);
    index = pairsModel->index(pairIndex,7,QModelIndex());
    pairsModel->setData(index,pair.heavyZ,Qt::EditRole);
    index = pairsModel->index(pairIndex,8,QModelIndex());
    pairsModel->setData(index,pair.heavyM,Qt::EditRole);
    index = pairsModel->index(pairIndex,9,QModelIndex());
    pairsModel->setData(index,pair.heavyG,Qt::EditRole);
    index = pairsModel->index(pairIndex,10,QModelIndex());
    pairsModel->setData(index,pair.excitationEnergy,Qt::EditRole);
    index = pairsModel->index(pairIndex,11,QModelIndex());
    pairsModel->setData(index,pair.seperationEnergy,Qt::EditRole);
    index = pairsModel->index(pairIndex,12,QModelIndex());
    pairsModel->setData(index,pair.channelRadius,Qt::EditRole);
    index = pairsModel->index(pairIndex,13,QModelIndex());
    pairsModel->setData(index,pair.pairType,Qt::EditRole);
    index = pairsModel->index(pairIndex,14,QModelIndex());
    pairsModel->setData(index,pair.ecMultMask,Qt::EditRole);
    index = pairsModel->index(pairIndex,15,QModelIndex());
    pairsModel->setData(index,pair.isUPOS,Qt::EditRole);
//    index = pairsModel->index(pairIndex,16,QModelIndex());
//    pairsModel->setData(index,pair.secondaryDecayL,Qt::EditRole);
//    index = pairsModel->index(pairIndex,17,QModelIndex());
//    pairsModel->setData(index,pair.Ic,Qt::EditRole);

    pairsView->resizeRowsToContents();
    if(!fromFile) emit(pairAdded(pairIndex));
  } else {
    QMessageBox::information(this,tr("Duplicate Pair"),
			     tr("This pair already exists."));
  }
}

void PairsTab::removePair() {
  QItemSelectionModel *selectionModel = pairsView->selectionModel();
  QModelIndexList indexes = selectionModel->selectedRows();
  QModelIndex index=indexes[0];
  bool previousIsGamma=false;
  if(index.row()==0&&(pairsModel->getPairs()).size()!=1) {
    QModelIndex previousIndex=pairsModel->index(index.row()+1,16,QModelIndex()); //was 13, changing to 16?
    QVariant previousPairType=pairsModel->data(previousIndex,Qt::EditRole);
    if (previousPairType.toInt()==10) previousIsGamma=true;
  }
  if(!previousIsGamma) {
    pairsModel->removeRows(index.row(),1,QModelIndex());
    emit(pairRemoved(index.row()));
  } else {
    QMessageBox::information(this,tr("Entrance Channel Error"),
			     tr("This delete will result in a forbidden particle,gamma entrance pair."));
  }
}

void PairsTab::editPair() {
  QItemSelectionModel *selectionModel = pairsView->selectionModel();
  QModelIndexList indexes = selectionModel->selectedRows(); 
  
  QModelIndex index=indexes[0];
  QModelIndex i = pairsModel->index(index.row(), 0, QModelIndex());
  QVariant var = pairsModel->data(i, Qt::EditRole);
  QString lightJ = var.toString();
  i = pairsModel->index(index.row(), 1, QModelIndex());
  var =pairsModel->data(i, Qt::EditRole);
  int lightPi = var.toInt();
  i = pairsModel->index(index.row(), 2, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString lightZ = var.toString();
  i = pairsModel->index(index.row(), 3, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString lightM = var.toString();
  i = pairsModel->index(index.row(), 4, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString lightG = var.toString();
  i = pairsModel->index(index.row(), 5, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString  heavyJ = var.toString();
  i = pairsModel->index(index.row(), 6, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  int  heavyPi = var.toInt();
  i = pairsModel->index(index.row(), 7, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString heavyZ = var.toString();
  i = pairsModel->index(index.row(), 8, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString heavyM = var.toString();
  i = pairsModel->index(index.row(), 9, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString heavyG = var.toString();
  i = pairsModel->index(index.row(), 10, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString excitationEnergy = var.toString();
  i = pairsModel->index(index.row(), 11, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString seperationEnergy = var.toString();
  i = pairsModel->index(index.row(), 12, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  QString channelRadius = var.toString();
  i = pairsModel->index(index.row(), 13, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  int pairType = var.toInt();
  i = pairsModel->index(index.row(), 14, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  int ecMultMask = var.toInt();
  i = pairsModel->index(index.row(), 15, QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  int isUPOS = var.toInt();
//  std::cout<<isUPOS<<","<<var.toString().toUtf8().constData()<<std::endl; 
//  i = pairsModel->index(index.row(), 16, QModelIndex());
//  var = pairsModel->data(i, Qt::EditRole);
//  QString secondaryDecayL = var.toString();
//  i = pairsModel->index(index.row(), 17, QModelIndex());
//  var = pairsModel->data(i, Qt::EditRole);
//  QString Ic = var.toString();

  
  AddPairDialog aDialog;
  aDialog.setWindowTitle(tr("Edit a Particle Pair"));
  aDialog.lightJText->setText(lightJ);
  if(lightPi==-1) aDialog.lightPiCombo->setCurrentIndex(0);
  else aDialog.lightPiCombo->setCurrentIndex(1);
  aDialog.lightZText->setText(lightZ);
  aDialog.lightMText->setText(lightM);
  aDialog.lightGText->setText(lightG);
  aDialog.heavyJText->setText(heavyJ);
  if(heavyPi==-1) aDialog.heavyPiCombo->setCurrentIndex(0);
  else aDialog.heavyPiCombo->setCurrentIndex(1);
  aDialog.heavyZText->setText(heavyZ);
  aDialog.heavyMText->setText(heavyM);
  aDialog.heavyGText->setText(heavyG);
  aDialog.excitationEnergyText->setText(excitationEnergy);
  aDialog.seperationEnergyText->setText(seperationEnergy);
  aDialog.channelRadiusText->setText(channelRadius);
  if(pairType == 10) aDialog.pairTypeCombo->setCurrentIndex(1);
  else if(pairType == 20) aDialog.pairTypeCombo->setCurrentIndex(2);
  else aDialog.pairTypeCombo->setCurrentIndex(0);
  if(ecMultMask&(1<<0)) aDialog.e1Check->setChecked(true);
  else aDialog.e1Check->setChecked(false);
  if(ecMultMask&(1<<1)) aDialog.m1Check->setChecked(true);
  else aDialog.m1Check->setChecked(false);
  if(ecMultMask&(1<<2)) aDialog.e2Check->setChecked(true);
  else aDialog.e2Check->setChecked(false);
//  if(isUPOS==1) aDialog.uposCheck->setChecked(true);
//  else aDialog.uposCheck->setChecked(false);
//  aDialog.secondaryDecayLText->setText(secondaryDecayL);
//  aDialog.icText->setText(Ic);

  if (aDialog.exec()) {
    if(!(index.row()==0&&aDialog.pairTypeCombo->currentIndex()!=0)) {
      PairsData pair;
      pair.lightJ = aDialog.lightJText->text().toDouble();
      if(aDialog.lightPiCombo->currentIndex()==0) pair.lightPi=-1;
      else pair.lightPi=1;
      pair.lightZ = aDialog.lightZText->text().toInt();
      pair.lightM = aDialog.lightMText->text().toDouble();
      pair.lightG = aDialog.lightGText->text().toDouble();
      pair.heavyJ = aDialog.heavyJText->text().toDouble();
      if(aDialog.heavyPiCombo->currentIndex()==0) pair.heavyPi=-1;
      else pair.heavyPi=1;
      pair.heavyZ = aDialog.heavyZText->text().toInt();
      pair.heavyM = aDialog.heavyMText->text().toDouble();
      pair.heavyG = aDialog.heavyGText->text().toDouble();
      pair.excitationEnergy = aDialog.excitationEnergyText->text().toDouble();
      pair.seperationEnergy = aDialog.seperationEnergyText->text().toDouble();
      pair.channelRadius = aDialog.channelRadiusText->text().toDouble();
      if(aDialog.pairTypeCombo->currentIndex()==1) pair.pairType=10;
      else if(aDialog.pairTypeCombo->currentIndex()==2) pair.pairType=20;
      else pair.pairType=0;
      unsigned char newECMultMask=0;
      if(aDialog.e1Check->isChecked()) newECMultMask |= (1<<0);
      if(aDialog.m1Check->isChecked()) newECMultMask |= (1<<1);
      if(aDialog.e2Check->isChecked()) newECMultMask |= (1<<2);
      pair.ecMultMask=newECMultMask;
//      if(aDialog.uposCheck->isChecked()) pair.isUPOS = 1;
//      else pair.isUPOS = 0;
      pair.isUPOS = 0;
//      std::cout<<pair.isUPOS<<std::endl;
//      pair.secondaryDecayL = aDialog.secondaryDecayLText->text().toInt();
//      pair.Ic = aDialog.icText->text().toDouble(); 
      editPair(pair,index.row(),false);
    } else {
      QMessageBox::information(this,tr("Pair Type Error"),
			       tr("The first pair must be a particle,particle pair."));
    }
  }
}

void PairsTab::editPair(PairsData pair,int pairIndex,bool fromFile) {
  QModelIndex i = pairsModel->index(pairIndex,0,QModelIndex());
  QVariant var = pairsModel->data(i, Qt::EditRole);
  if (pair.lightJ != var.toDouble()) pairsModel->setData(i,pair.lightJ, Qt::EditRole);
  i = pairsModel->index(pairIndex,1,QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.lightPi != var.toInt()) pairsModel->setData(i,pair.lightPi, Qt::EditRole);
  i = pairsModel->index(pairIndex,2,QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.lightZ != var.toInt()) pairsModel->setData(i,pair.lightZ, Qt::EditRole);
  i = pairsModel->index(pairIndex,3,QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.lightM != var.toDouble()) pairsModel->setData(i,pair.lightM, Qt::EditRole);
  i = pairsModel->index(pairIndex,4,QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.lightG != var.toDouble()) pairsModel->setData(i,pair.lightG, Qt::EditRole);
  i = pairsModel->index(pairIndex,5,QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.heavyJ != var.toDouble()) pairsModel->setData(i,pair.heavyJ, Qt::EditRole);
  i = pairsModel->index(pairIndex,6,QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.heavyPi != var.toInt()) pairsModel->setData(i,pair.heavyPi, Qt::EditRole);
  i = pairsModel->index(pairIndex,7,QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.heavyZ != var.toInt()) pairsModel->setData(i,pair.heavyZ, Qt::EditRole);
  i = pairsModel->index(pairIndex,8,QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.heavyM != var.toDouble()) pairsModel->setData(i,pair.heavyM, Qt::EditRole);
  i = pairsModel->index(pairIndex,9,QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.heavyG != var.toDouble()) pairsModel->setData(i,pair.heavyG, Qt::EditRole);
  i = pairsModel->index(pairIndex,10,QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.excitationEnergy != var.toDouble()) pairsModel->setData(i,pair.excitationEnergy, Qt::EditRole);
  i = pairsModel->index(pairIndex,11,QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.seperationEnergy != var.toDouble()) pairsModel->setData(i,pair.seperationEnergy, Qt::EditRole);
  i = pairsModel->index(pairIndex,12,QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.channelRadius != var.toDouble()) pairsModel->setData(i,pair.channelRadius, Qt::EditRole);
  i = pairsModel->index(pairIndex,13,QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.pairType != var.toInt()) pairsModel->setData(i,pair.pairType, Qt::EditRole);
  i = pairsModel->index(pairIndex,14,QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.ecMultMask != var.toInt()) pairsModel->setData(i,pair.ecMultMask, Qt::EditRole);
  i = pairsModel->index(pairIndex,15,QModelIndex());
  var = pairsModel->data(i, Qt::EditRole);
  if (pair.isUPOS != var.toInt()) pairsModel->setData(i,pair.isUPOS, Qt::EditRole);
//  i = pairsModel->index(pairIndex,16,QModelIndex());
//  var = pairsModel->data(i, Qt::EditRole);
//  if (pair.secondaryDecayL != var.toInt()) pairsModel->setData(i,pair.secondaryDecayL, Qt::EditRole);
//  i = pairsModel->index(pairIndex,17,QModelIndex());
//  var = pairsModel->data(i, Qt::EditRole); 
//  if (pair.Ic != var.toDouble()) pairsModel->setData(i,pair.Ic, Qt::EditRole);

  if(!fromFile) emit(pairEdited(pairIndex));
}

void PairsTab::updateButtons(const QItemSelection &selection) {
  QModelIndexList indexes = selection.indexes();
  
  if (!indexes.isEmpty()) {
    deleteButton->setEnabled(true);
  } else {
    deleteButton->setEnabled(false);
  }
}

bool PairsTab::parseOldECSection(QTextStream& inStream) {
  int isActive;
  int exitPairIndex;
  double minJ;
  double maxJ;
  int multMask;
  QString line("");
  while(line.trimmed()!=QString("</externalCapture>")&&!inStream.atEnd()) {
    line = inStream.readLine();
    if(line.trimmed().isEmpty()) continue;
    if(line.trimmed()!=QString("</externalCapture>")&&!inStream.atEnd()) {
      QTextStream in(&line);
      in >> isActive >> exitPairIndex >> minJ >> maxJ >> multMask;
      if(in.status()!=QTextStream::Ok) return false;
      if(exitPairIndex-1<pairsModel->getPairs().size()) {
	QModelIndex i = pairsModel->index(exitPairIndex-1,14,QModelIndex());
	QVariant var = pairsModel->data(i, Qt::EditRole);
	pairsModel->setData(i,multMask, Qt::EditRole);
      }
    }
  }
  if(line.trimmed()!=QString("</externalCapture>")) return false;
  return true;
}

void PairsTab::showInfo(int which, QString title) {
  if(which<infoText.size()) {
    if(!infoDialog[which]) {
      infoDialog[which] = new InfoDialog(infoText[which],this,title);
      infoDialog[which]->setAttribute(Qt:: WA_DeleteOnClose);
      infoDialog[which]->show();
    } else infoDialog[which]->raise();
  }
}
