#include "PairsModel.h"
#include "ElementMap.h"
#include <math.h>
#include <iostream>

PairsModel::PairsModel(QObject *parent) : QAbstractTableModel(parent) {
}

int PairsModel::rowCount(const QModelIndex &parent) const {
  Q_UNUSED(parent);
  return pairsList.size();
}

int PairsModel::columnCount(const QModelIndex &parent) const {
  Q_UNUSED(parent);
  return PairsData::SIZE;
}

QVariant PairsModel::data(const QModelIndex &index, int role) const {
  if(!index.isValid()) return QVariant();

  if(index.row() >= pairsList.size() || index.row() < 0)  return QVariant();

  if (role == Qt::DisplayRole) {
    PairsData pair = pairsList.at(index.row());
    if(index.column() == 2) return getSpinLabel(pair,0);
    else if(index.column() == 3) return pair.lightPi;
    else if(index.column() == 0) return getParticleLabel(pair,0);      
    else if(index.column() == 1) return pair.lightM;
    else if(index.column() == 4) return pair.lightG;
    else if(index.column() == 7) return getSpinLabel(pair,1);
    else if(index.column() == 8) return pair.heavyPi;
    else if(index.column() == 5) return getParticleLabel(pair,1);     
    else if(index.column() == 6) return pair.heavyM;
    else if(index.column() == 9) return pair.heavyG;
    else if(index.column() == 10) return pair.excitationEnergy;
    else if(index.column() == 11) return pair.seperationEnergy;
    else if(index.column() == 12) return pair.channelRadius;
    else if(index.column() == 13) {
      if(pair.pairType == 10) return "Particle, Gamma";
      else if(pair.pairType == 20) return "Beta Decay";
      else return "Particle, Particle";
    } else if(index.column() == 14) return pair.ecMultMask;
    else if(index.column() == 15) return pair.isUPOS;
//    else if(index.column() == 16) return pair.secondaryDecayL;
//    else if(index.column() == 17) return pair.Ic;
  } else if (role == Qt::EditRole) {
    PairsData pair = pairsList.at(index.row());
    if(index.column() == 0) return pair.lightJ;
    else if(index.column() == 1) return pair.lightPi;
    else if(index.column() == 2) return pair.lightZ;
    else if(index.column() == 3) return pair.lightM;
    else if(index.column() == 4) return pair.lightG;
    else if(index.column() == 5) return pair.heavyJ;
    else if(index.column() == 6) return pair.heavyPi;
    else if(index.column() == 7) return pair.heavyZ;
    else if(index.column() == 8) return pair.heavyM;
    else if(index.column() == 9) return pair.heavyG;
    else if(index.column() == 10) return pair.excitationEnergy;
    else if(index.column() == 11) return pair.seperationEnergy;
    else if(index.column() == 12) return pair.channelRadius;
    else if(index.column() == 13) return pair.pairType; 
    else if(index.column() == 14) return pair.ecMultMask;
    else if(index.column() == 15) return pair.isUPOS;
//    else if(index.column() == 16) return pair.secondaryDecayL;
//    else if(index.column() == 17) return pair.Ic;  
  } else if(role == Qt::TextAlignmentRole) return Qt::AlignCenter;

  return QVariant();
}

QVariant PairsModel::headerData(int section, Qt::Orientation orientation, int role) const {
  if(role!= Qt::DisplayRole) return QVariant();
  if(orientation == Qt::Horizontal) {
    switch(section) {
    case 2:
      return tr("Light\nSpin");
    case 3:
      return tr("Light\nParity");
    case 0:
      return tr("Light\n Particle");
    case 1:
      return tr("Light\nM");
    case 4:
      return tr("Light\ng-Factor");
    case 7: 
      return tr("Heavy\nSpin");
    case 8:
      return tr("Heavy\nParity");
    case 5:
      return tr("Heavy\nParticle");
    case 6:
      return tr("Heavy\nM");
    case 9:
      return tr("Heavy\ng-Factor");
    case 10:
      return tr("Excitation\nEnergy");
    case 11:
      return tr("Separation\nEnergy");
    case 12:
      return tr("Channel\nRadius");
    case 13:
      return tr("Pair\nType");
    case 14:
      return tr("EC\nMultipolarities");
    case 15:
      return tr("Unobserved?");
//    case 16:
//      return tr("Final\nL");
//    case 17:
//      return tr("Final\nJ");
    default: 
      return QVariant();
    }
  } else if(orientation == Qt::Vertical) {
    return section+1;
  }
  return QVariant();
}

bool PairsModel::setData(const QModelIndex &index, const QVariant &value, int role) {
  if (index.isValid() && role == Qt::EditRole ) {
    int row = index.row();
    PairsData tempData = pairsList.value(row);
    if (index.column() == 0) tempData.lightJ=value.toDouble();
    else if(index.column() == 1) tempData.lightPi=value.toInt();
    else if(index.column() == 2) tempData.lightZ=value.toInt();
    else if(index.column() == 3) tempData.lightM=value.toDouble();
    else if(index.column() == 4) tempData.lightG=value.toDouble();
    else if(index.column() == 5) tempData.heavyJ=value.toDouble();
    else if(index.column() == 6) tempData.heavyPi=value.toInt();
    else if(index.column() == 7) tempData.heavyZ=value.toInt();
    else if(index.column() == 8) tempData.heavyM=value.toDouble();
    else if(index.column() == 9) tempData.heavyG=value.toDouble();
    else if(index.column() == 10) tempData.excitationEnergy=value.toDouble();
    else if(index.column() == 11) tempData.seperationEnergy=value.toDouble();
    else if(index.column() == 12) tempData.channelRadius=value.toDouble();
    else if(index.column() == 13) tempData.pairType=value.toInt();
    else if(index.column() == 14) tempData.ecMultMask=value.toInt();
    else if(index.column() == 15) tempData.isUPOS=value.toInt();
//    else if(index.column() == 16) tempData.secondaryDecayL=value.toInt();
//    else if(index.column() == 17) tempData.Ic=value.toDouble();
    else return false;

    pairsList.replace(row,tempData);
    emit(dataChanged(index,index));
    return true;
  }
  return false;
}

bool PairsModel::insertRows(int position, int rows, const QModelIndex &index) {
  Q_UNUSED(index);
  if(rows>0) {
    beginInsertRows(QModelIndex(),position,position+rows-1);
    for(int row=0; row<rows; row++) {
      PairsData tempData;
      pairsList.insert(position,tempData);
    }
    endInsertRows();
  }
  return true;
}

bool PairsModel::removeRows(int position, int rows, const QModelIndex &index) {
  Q_UNUSED(index);
  if(rows>0) {
    beginRemoveRows(QModelIndex(),position,position+rows-1);
    for(int row=0; row<rows;++row) {
      pairsList.removeAt(position);
    }
    endRemoveRows();
  }
  return true;
}

int PairsModel::isPair(const PairsData &pair) const {
  int foundPair=-1;
  for(int i=0;i<pairsList.size();i++) {
    PairsData tempPair=pairsList.value(i);
    if(tempPair.lightJ==pair.lightJ&&
       tempPair.lightPi==pair.lightPi&&
       tempPair.lightZ==pair.lightZ&&
       tempPair.lightM==pair.lightM&&
       tempPair.lightG==pair.lightG&&
       tempPair.heavyJ==pair.heavyJ&&
       tempPair.heavyPi==pair.heavyPi&&
       tempPair.heavyZ==pair.heavyZ&&
       tempPair.heavyM==pair.heavyM&&
       tempPair.heavyG==pair.heavyG&&
       tempPair.seperationEnergy==pair.seperationEnergy&&
       tempPair.excitationEnergy==pair.excitationEnergy&&
       tempPair.channelRadius==pair.channelRadius&&
       tempPair.pairType==pair.pairType&&
       tempPair.ecMultMask==pair.ecMultMask&&
       tempPair.isUPOS==pair.isUPOS
//       tempPair.secondaryDecayL==pair.secondaryDecayL&&
//       tempPair.Ic==pair.Ic) 
     ){
      foundPair=i;
      break;
    }
  }
  return foundPair;
}

QString PairsModel::getParticleLabel(const PairsData &pair, int which) const {
  if(which!=-1) {
    if(pair.pairType==10&&which==0) return "<center>&gamma;</center>";
    else if(pair.pairType==20&&which==0) {
      if(pair.lightZ<0) return "<center>&beta;<sup>-</sup></center>";
      else return "<center>&beta;<sup>+</sup></center>";
    } else {
      int tempZ;
      int tempM;
      if(which==0) {
	tempZ=pair.lightZ;
	tempM=round(pair.lightM);
      } else {
	tempZ=pair.heavyZ;
	tempM=round(pair.heavyM);
      }
      std::map<int, QString>::const_iterator it=elementMap.find(tempZ);
      if(it!=elementMap.end()) {
	if(tempM==1) {
	  if(tempZ==1) return "<center><i>p</i></center>";
	  else return QString("<center><i>%1</i></center>").arg(it->second);
	} else if(tempZ==2&&tempM==4) return "<center>&alpha;</center>";
	else return QString("<center><sup>%1</sup>%2</center>").arg(tempM).arg(it->second);
      } else return "?";
    }
  } else {
    QString lightLabel;
    std::map<int, QString>::const_iterator it=elementMap.find(pair.lightZ);
    if(pair.pairType==10) lightLabel="&gamma;";
    else if(pair.pairType==20) {
      if(pair.lightZ<0.) lightLabel="&beta;<sup>-</sup>";
      else lightLabel="&beta;<sup>+</sup>";
    } else if(it!=elementMap.end()) {
      if(round(pair.lightM)==1) {
	if(pair.lightZ==1) lightLabel="<i>p</i>";
	else lightLabel=QString("<i>%1</i>").arg(it->second);
      } else if(pair.lightZ==2&&round(pair.lightM)==4) lightLabel="&alpha;";
      else lightLabel=QString("<sup>%1</sup>%2").arg(round(pair.lightM)).arg(it->second);
    } else lightLabel="?";
    QString heavyLabel;
    it=elementMap.find(pair.heavyZ);
    if(it!=elementMap.end()) {
      if(round(pair.heavyM)==1) {
	if(pair.heavyZ==1) heavyLabel="<i>p</i>";
	else heavyLabel=QString("<i>%1</i>").arg(it->second);
      } else if(pair.heavyZ==2&&round(pair.heavyM)==4) heavyLabel="&alpha;";
      else heavyLabel=QString("<sup>%1</sup>%2").arg(round(pair.heavyM)).arg(it->second);
    } else heavyLabel="?";
    if(pair.pairType==20) return QString("<center>%1(%2) [%3 MeV]</center>").arg(heavyLabel).arg(lightLabel).arg(pair.excitationEnergy,0,'f',3);
    else return QString("<center>%1+%2 [%3 MeV]</center>").arg(heavyLabel).arg(lightLabel).arg(pair.excitationEnergy,0,'f',3);
  }
}

QString PairsModel::getReactionLabel(const PairsData &firstPair, const PairsData &secondPair) {
  QString lightLabel[2];
  std::map<int, QString>::const_iterator it=elementMap.find(firstPair.lightZ);
  if(firstPair.pairType==10) lightLabel[0]="&gamma;";
  else if(firstPair.pairType==20) {
    if(firstPair.lightZ<0) lightLabel[0]="&beta;<sup>-</sup>";
    else lightLabel[0]="&beta;<sup>+</sup>";
  } else if(it!=elementMap.end()) {
    if(round(firstPair.lightM)==1) {
      if(firstPair.lightZ==1) lightLabel[0]="<i>p</i>";
      else lightLabel[0]=QString("<i>%1</i>").arg(it->second);
    } else if(firstPair.lightZ==2&&round(firstPair.lightM)==4) lightLabel[0]="&alpha;";
    else lightLabel[0]=QString("<sup>%1</sup>%2").arg(round(firstPair.lightM)).arg(it->second);
  } else lightLabel[0]="?";
  it=elementMap.find(secondPair.lightZ);
  if(secondPair.pairType==10) lightLabel[1]="&gamma;";
  else if(secondPair.pairType==20) {
    if(secondPair.lightZ<0) lightLabel[1]="&beta;<sup>-</sup>";
    else lightLabel[1]="&beta;<sup>+</sup>";
  } else if(it!=elementMap.end()) {
    if(round(secondPair.lightM)==1) {
      if(secondPair.lightZ==1) lightLabel[1]="<i>p</i>";
      else lightLabel[1]=QString("<i>%1</i>").arg(it->second);
    } else if(secondPair.lightZ==2&&round(secondPair.lightM)==4) lightLabel[1]="&alpha;";
    else lightLabel[1]=QString("<sup>%1</sup>%2").arg(round(secondPair.lightM)).arg(it->second);
  } else lightLabel[1]="?";
  QString heavyLabel[2];
  it=elementMap.find(firstPair.heavyZ);
  if(it!=elementMap.end()) {
    if(round(firstPair.heavyM)==1) {
      if(firstPair.heavyZ==1) heavyLabel[0]="<i>p</i>";
      else heavyLabel[0]=QString("<i>%1</i>").arg(it->second);
    } else if(firstPair.heavyZ==2&&round(firstPair.heavyM)==4) heavyLabel[0]="&alpha;";
    else heavyLabel[0]=QString("<sup>%1</sup>%2").arg(round(firstPair.heavyM)).arg(it->second);
  } else heavyLabel[0]="?";
  it=elementMap.find(secondPair.heavyZ);
  if(it!=elementMap.end()) {
    if(round(secondPair.heavyM)==1) {
      if(secondPair.heavyZ==1) heavyLabel[1]="<i>p</i>";
      else heavyLabel[1]=QString("<i>%1</i>").arg(it->second);
    } else if(secondPair.heavyZ==2&&round(secondPair.heavyM)==4) heavyLabel[1]="&alpha;";
    else heavyLabel[1]=QString("<sup>%1</sup>%2").arg(round(secondPair.heavyM)).arg(it->second);
  } else heavyLabel[1]="?";
  if(firstPair.pairType==20) return QString("%1(%2%3)%4 [%5 MeV]").arg(heavyLabel[0]).arg(lightLabel[0]).arg(lightLabel[1]).arg(heavyLabel[1]).arg(secondPair.excitationEnergy,0,'f',3);
  else return QString("%1(%2,%3)%4 [%5 MeV]").arg(heavyLabel[0]).arg(lightLabel[0]).arg(lightLabel[1]).arg(heavyLabel[1]).arg(secondPair.excitationEnergy,0,'f',3);
}

QString PairsModel::getReactionLabelTotalCapture(const PairsData &firstPair) {
  QString lightLabel[2];
  std::map<int, QString>::const_iterator it=elementMap.find(firstPair.lightZ);
  if(firstPair.pairType==10) lightLabel[0]="&gamma;";
  else if(firstPair.pairType==20) {
    if(firstPair.lightZ<0) lightLabel[0]="&beta;<sup>-</sup>";
    else lightLabel[0]="&beta;<sup>+</sup>";
  } else if(it!=elementMap.end()) {
    if(round(firstPair.lightM)==1) {
      if(firstPair.lightZ==1) lightLabel[0]="<i>p</i>";
      else lightLabel[0]=QString("<i>%1</i>").arg(it->second);
    } else if(firstPair.lightZ==2&&round(firstPair.lightM)==4) lightLabel[0]="&alpha;";
    else lightLabel[0]=QString("<sup>%1</sup>%2").arg(round(firstPair.lightM)).arg(it->second);
  } else lightLabel[0]="?";
  lightLabel[1]="&gamma;";
  QString heavyLabel[2];
  it=elementMap.find(firstPair.heavyZ);
  if(it!=elementMap.end()) {
    if(round(firstPair.heavyM)==1) {
      if(firstPair.heavyZ==1) heavyLabel[0]="<i>p</i>";
      else heavyLabel[0]=QString("<i>%1</i>").arg(it->second);
    } else if(firstPair.heavyZ==2&&round(firstPair.heavyM)==4) heavyLabel[0]="&alpha;";
    else heavyLabel[0]=QString("<sup>%1</sup>%2").arg(round(firstPair.heavyM)).arg(it->second);
  } else heavyLabel[0]="?";
  int i;
  for(i=0;i<pairsList.size();i++) 
    if(pairsList[i].pairType==10) break;
  if(i==pairsList.size()) heavyLabel[1]="?";
  else {
    PairsData secondPair = pairsList[i];
    it=elementMap.find(secondPair.heavyZ);
    if(it!=elementMap.end()) {
      if(round(secondPair.heavyM)==1) {
	if(secondPair.heavyZ==1) heavyLabel[1]="<i>p</i>";
	else heavyLabel[1]=QString("<i>%1</i>").arg(it->second);
      } else if(secondPair.heavyZ==2&&round(secondPair.heavyM)==4) heavyLabel[1]="&alpha;";
      else heavyLabel[1]=QString("<sup>%1</sup>%2").arg(round(secondPair.heavyM)).arg(it->second);
    } else heavyLabel[1]="?";
  }
  if(firstPair.pairType==20) return QString("%1(%2%3)%4 [TOTAL]").arg(heavyLabel[0]).arg(lightLabel[0]).arg(lightLabel[1]).arg(heavyLabel[1]);
  else return QString("%1(%2,%3)%4 [TOTAL]").arg(heavyLabel[0]).arg(lightLabel[0]).arg(lightLabel[1]).arg(heavyLabel[1]);
}

QString PairsModel::getSpinLabel(const PairsData &pair, int which) const {
  double tempJ;
  int tempPi;
  if(which==0) {
    tempJ=pair.lightJ;
    tempPi=pair.lightPi;
  } else {
    tempJ=pair.heavyJ;
    tempPi=pair.heavyPi;
  }
  QString tempSpin;
  if(((int)(tempJ*2))%2!=0&&tempJ!=0.) tempSpin=QString("%1/2").arg((int)(tempJ*2));
  else tempSpin=QString("%1").arg(tempJ);
  if(tempPi==-1) return QString("%1-").arg(tempSpin);
  else return QString("%1+").arg(tempSpin);
}
