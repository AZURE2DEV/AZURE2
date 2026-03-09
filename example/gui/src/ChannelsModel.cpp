#include "ChannelsModel.h"
#include "PairsModel.h"

ChannelsModel::ChannelsModel(QObject *parent) : QAbstractTableModel(parent) {
}

int ChannelsModel::rowCount(const QModelIndex &parent) const {
  Q_UNUSED(parent);
  return channelsList.size();
}

int ChannelsModel::columnCount(const QModelIndex &parent) const {
  Q_UNUSED(parent);
  return ChannelsData::SIZE;
}

QVariant ChannelsModel::data(const QModelIndex &index, int role) const {
  if(!index.isValid()) return QVariant();

  if(index.row() >= channelsList.size() || index.row() < 0)  return QVariant();

  if (role == Qt::DisplayRole) {
    ChannelsData channel = channelsList.at(index.row());
    if(index.column() == 1) return channel.levelIndex;
    else if(index.column() == 2) {
      PairsData pair=pairsModel->getPairs().at(channel.pairIndex);
      if(channel.reducedWidth!=0.) 
	return QString("<center><font style='font-weight:bold;'>%1</font></center>").arg(pairsModel->getParticleLabel(pair));
      else return QString("<center>%1</center>").arg(pairsModel->getParticleLabel(pair));
    } else if(index.column() == 3) {
       if(channel.reducedWidth!=0.) 
	 return QString("<center><font style='font-weight:bold;'>%1</font></center>").arg(getSpinLabel(channel));
       else return QString("<center>%1</center>").arg(getSpinLabel(channel));
    } else if(index.column() == 4) {
      if(channel.radType=='P') {
	if(channel.reducedWidth!=0.) 
	  return QString("<center><font style='font-weight:bold;'>%1</font></center>").arg(channel.lValue);
	else return QString("<center>%1</center>").arg(channel.lValue);
      } else if(channel.radType=='F') {
	if(channel.reducedWidth!=0.) 
	  return QString("<center><font style='font-weight:bold;'>F</font></center>");
	else return QString("<center>F</center>");
      } else if(channel.radType=='G') {
	if(channel.reducedWidth!=0.) 
	  return QString("<center><font style='font-weight:bold;'>GT</font></center>");
	else return QString("<center>GT</center>");
      } else {
	if(channel.reducedWidth!=0.)
	  return QString("<center><font style='font-weight:bold;'>%1%2</font></center>").arg(channel.radType).arg(channel.lValue);
	else return QString("<center>%1%2</center>").arg(channel.radType).arg(channel.lValue);
      }
    } else if(index.column() == 5) return channel.radType;
    else if(index.column() == 6) return channel.reducedWidth;
  } else if(role==Qt::EditRole) {
    ChannelsData channel = channelsList.at(index.row());
    if(index.column() == 1) return channel.levelIndex;
    else if(index.column() == 2) return channel.pairIndex;
    else if(index.column() == 3) return channel.sValue;
    else if(index.column() == 4) return channel.lValue;
    else if(index.column() == 5) return channel.radType;
    else if(index.column() == 6) return channel.reducedWidth;
  } else if(role==Qt::TextAlignmentRole) return Qt::AlignCenter;
  else if (role==Qt::CheckStateRole && index.column()==0) {
    ChannelsData channel = channelsList.at(index.row());
    if(channel.isFixed==1) return Qt::Checked;
    else return Qt::Unchecked;
  }
  return QVariant();
}

QVariant ChannelsModel::headerData(int section, Qt::Orientation orientation, int role) const {
  if(role!= Qt::DisplayRole) return QVariant();
  if(orientation == Qt::Horizontal) {
    switch(section) {
    case 0:
      return tr("Fix?");
    case 1:
      return tr("level");
    case 2:
      return tr("Channel\nPair");
    case 3:
      return tr("s");
    case 4:
      return tr("l");
    case 5:
      return tr("radiation type");
    case 6:
      return tr("reduced width");
    default: 
      return QVariant();
    }
  } else if(orientation == Qt::Vertical) {
    return section+1;
  }
  return QVariant();
}

bool ChannelsModel::setData(const QModelIndex &index, const QVariant &value, int role) {
  if (index.isValid()){
    if(role == Qt::EditRole ) {
      int row = index.row();
      ChannelsData tempData = channelsList.value(row);
      if (index.column() == 0) tempData.isFixed=value.toInt();
      else if(index.column() == 1) tempData.levelIndex=value.toInt();
      else if(index.column() == 2) tempData.pairIndex=value.toInt();
      else if(index.column() == 3) tempData.sValue=value.toDouble();
      else if(index.column() == 4) tempData.lValue=value.toInt();
      else if(index.column() == 5) tempData.radType=value.toChar();
      else if(index.column() == 6) tempData.reducedWidth=value.toDouble();
      else return false;
      
      channelsList.replace(row,tempData);
      if(index.column()!=6) emit(dataChanged(index,index));
      return true;
    } else if(role== Qt::CheckStateRole) {
      int row = index.row();
      ChannelsData tempData = channelsList.value(row);
      if(index.column()==0) {
	if(value==Qt::Checked) tempData.isFixed=1;
	else tempData.isFixed=0;
      } else return false;
      
      channelsList.replace(row,tempData);
      emit(dataChanged(index,index));
      return true;
    }
  }
  return false;
}

bool ChannelsModel::insertRows(int position, int rows, const QModelIndex &index) {
  Q_UNUSED(index);
  if(rows>0) {
    beginInsertRows(QModelIndex(),position,position+rows-1);
    for(int row=0; row<rows; row++) {
      ChannelsData tempData={0,-1,-1,0.0,0,'P',0.0};
      channelsList.insert(position,tempData);
    }
    endInsertRows();
  }
  return true;
}

bool ChannelsModel::removeRows(int position, int rows, const QModelIndex &index) {
  Q_UNUSED(index);
  if(rows>0) {
    beginRemoveRows(QModelIndex(),position,position+rows-1);
    for(int row=0; row<rows;++row) {
      channelsList.removeAt(position);
    }
    endRemoveRows();
  }
  return true;
}

Qt::ItemFlags ChannelsModel::flags(const QModelIndex &index) const {
  if (!index.isValid()) return Qt::ItemIsEnabled;
  if(index.column()==0) return QAbstractTableModel::flags(index) | Qt::ItemIsUserCheckable;
  return QAbstractTableModel::flags(index);
}

bool ChannelsModel::isChannel(const ChannelsData &channel)  const {
  bool foundChannel=false;
  for(int i=0;i<channelsList.size();i++) {
    ChannelsData tempChannel=channelsList.value(i);
    if(tempChannel.levelIndex==channel.levelIndex&&
       tempChannel.pairIndex==channel.pairIndex&&
       tempChannel.sValue==channel.sValue&&
       tempChannel.lValue==channel.lValue&&
       tempChannel.radType==channel.radType) {
      foundChannel=true;
      break;
    }
  }
  return foundChannel;

}

QString ChannelsModel::getSpinLabel(const ChannelsData &channel) const {
  if(((int)(channel.sValue*2))%2!=0&&channel.sValue!=0.) return QString("%1/2").arg((int)(channel.sValue*2));
  else return QString("%1").arg(channel.sValue);
}

void ChannelsModel::setPairsModel(PairsModel* model) {
  pairsModel=model;
}
