#include "LevelsModel.h"

LevelsModel::LevelsModel(QObject *parent) : QAbstractTableModel(parent) {
}

int LevelsModel::rowCount(const QModelIndex &parent) const {
  Q_UNUSED(parent);
  return levelsList.size();
}

int LevelsModel::columnCount(const QModelIndex &parent) const {
  Q_UNUSED(parent);
  return LevelsData::SIZE;
}

QVariant LevelsModel::data(const QModelIndex &index, int role) const {
  if(!index.isValid()) return QVariant();

  if(index.row() >= levelsList.size() || index.row() < 0)  return QVariant();

  if (role == Qt::DisplayRole) {
    LevelsData level = levelsList.at(index.row());
    if(index.column() == 2) return getSpinLabel(level);
    else if(index.column() == 3) return level.piValue;
    else if(index.column() == 4) return level.energy;
  } else if (role == Qt::EditRole) {
    LevelsData level = levelsList.at(index.row());
    if(index.column() == 2) return level.jValue;
    else if(index.column() == 3) return level.piValue;
    else if(index.column() == 4) return level.energy;
  } else if(role==Qt::TextAlignmentRole) return Qt::AlignCenter;
  else if (role==Qt::CheckStateRole && (index.column()==0||index.column()==1)) {
    LevelsData level = levelsList.at(index.row());
    if(index.column()==0) {
      if(level.isActive==1) return Qt::Checked;
      else return Qt::Unchecked;
    } else {
      if(level.isFixed==1) return Qt::Checked;
      else return Qt::Unchecked;
    }
  }
  return QVariant();
}

QVariant LevelsModel::headerData(int section, Qt::Orientation orientation, int role) const {
  if(role!= Qt::DisplayRole) return QVariant();
  if(orientation == Qt::Horizontal) {
    switch(section) {
    case 0:
      return tr("Include?");
    case 1:
      return tr("Fix?");
    case 2:
      return tr("Level\nSpin");
    case 3:
      return tr("Parity");
    case 4:
      return tr("Energy\n[MeV]");
    default: 
      return QVariant();
    }
  } else if(orientation == Qt::Vertical) {
    return section+1;
  }
  return QVariant();
}

bool LevelsModel::setData(const QModelIndex &index, const QVariant &value, int role) {
  if (index.isValid()) {
    if(role == Qt::EditRole ) {
      int row = index.row();
      LevelsData tempData = levelsList.value(row);
      if(index.column() == 0) tempData.isActive=value.toInt();
      else if(index.column() == 1) tempData.isFixed=value.toInt();
      else if(index.column() == 2) tempData.jValue=value.toDouble();
      else if(index.column() == 3) tempData.piValue=value.toInt();
      else if(index.column() == 4) tempData.energy=value.toDouble();
      else return false;
      
      levelsList.replace(row,tempData);
      emit(dataChanged(index,index));
      return true;
    } else if(role== Qt::CheckStateRole) {
      int row = index.row();
      LevelsData tempData = levelsList.value(row);
      if(index.column()==0) {
	if(value==Qt::Checked) tempData.isActive=1;
	else tempData.isActive=0;
      } else if(index.column()==1) {
	if(value==Qt::Checked) tempData.isFixed=1;
	else tempData.isFixed=0;
      } else return false;

      levelsList.replace(row,tempData);
      emit(dataChanged(index,index));
      return true;
    }
  }
  return false;
}

bool LevelsModel::insertRows(int position, int rows, const QModelIndex &index) {
  Q_UNUSED(index);
  if(rows>0) {
    beginInsertRows(QModelIndex(),position,position+rows-1);
    for(int row=0; row<rows; row++) {
      LevelsData tempData={1,0,0.0,-1,0.0};
      levelsList.insert(position,tempData);
    }
    endInsertRows();
  }
  return true;
}

bool LevelsModel::removeRows(int position, int rows, const QModelIndex &index) {
  Q_UNUSED(index);
  if(rows>0) {
    beginRemoveRows(QModelIndex(),position,position+rows-1);
    for(int row=0; row<rows;++row) {
      levelsList.removeAt(position);
    }
    endRemoveRows();
  }
  return true;
}

Qt::ItemFlags LevelsModel::flags(const QModelIndex &index) const {
  if (!index.isValid()) return Qt::ItemIsEnabled;
  if(index.column()==0 || index.column() ==1) return QAbstractTableModel::flags(index) | Qt::ItemIsUserCheckable;
  return QAbstractTableModel::flags(index);
}

int LevelsModel::isLevel(const LevelsData &level) const {
  int foundLevel=-1;
  for(int i=0;i<levelsList.size();i++) {
    LevelsData tempLevel=levelsList.value(i);
    if(tempLevel.jValue==level.jValue&&
       tempLevel.piValue==level.piValue&&
       tempLevel.energy==level.energy) {
      foundLevel=i;
      break;
    }
  }
  return foundLevel;
}

QString LevelsModel::getSpinLabel(const LevelsData &level) const {
  QString tempSpin;
  if(((int)(level.jValue*2))%2!=0&&level.jValue!=0.) tempSpin=QString("%1/2").arg((int)(level.jValue*2));
  else tempSpin=QString("%1").arg(level.jValue);
  if(level.piValue==-1) return QString("%1-").arg(tempSpin);
  else return QString("%1+").arg(tempSpin);
}
