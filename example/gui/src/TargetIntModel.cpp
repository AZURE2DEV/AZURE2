#include "TargetIntModel.h"

TargetIntModel::TargetIntModel(QObject *parent) : QAbstractTableModel(parent) {
}

int TargetIntModel::rowCount(const QModelIndex &parent) const {
  Q_UNUSED(parent);
  return targetIntList.size();
}

int TargetIntModel::columnCount(const QModelIndex &parent) const {
  Q_UNUSED(parent);
  return TargetIntData::SIZE;
}

QVariant TargetIntModel::data(const QModelIndex &index, int role) const {
  if(!index.isValid()) return QVariant();
  if(index.row()>=targetIntList.size() || index.row() < 0) return QVariant();
  if(role == Qt::DisplayRole) {
    TargetIntData targetInt=targetIntList.at(index.row());
    if(index.column() == 1) return targetInt.segmentsList;
    else if(index.column() == 2) {
      if(targetInt.isTargetIntegration||targetInt.isConvolution) return targetInt.numPoints;
      else return QString(tr("N/A"));
    }
    else if(index.column() == 3) {
      if(targetInt.isConvolution) return QString(tr("YES"));
      else return QString(tr("NO"));
    } else if (index.column() == 4) return targetInt.sigma;
    else if(index.column() == 5) {
      if(targetInt.isTargetIntegration) return QString(tr("YES"));
      else return QString(tr("NO"));
    } else if(index.column() == 6) return targetInt.density;
    else if(index.column() == 7) return targetInt.stoppingPowerEq;
    else if(index.column() == 8) return targetInt.numParameters;
    else if(index.column() == 9) return QVariant(); 
    else if(index.column() == 10) {
      if(targetInt.isQCoefficients) return QString(tr("YES"));
      else return QString(tr("NO"));
    } else if(index.column() == 11) return QVariant();
  } else if(role == Qt::EditRole) {
    TargetIntData targetInt=targetIntList.at(index.row());
    if(index.column()==1) return targetInt.segmentsList;
    if(index.column()==2) return targetInt.numPoints;
    if(index.column()==3) return targetInt.isConvolution;
    if(index.column()==4) return targetInt.sigma;
    if(index.column()==5) return targetInt.isTargetIntegration;
    if(index.column()==6) return targetInt.density;
    if(index.column()==7) return targetInt.stoppingPowerEq;
    if(index.column()==8) return targetInt.numParameters;
    if(index.column()==9) return QVariant::fromValue<QList<double> >(targetInt.parameters);    
    if(index.column()==10) return  targetInt.isQCoefficients;
    if(index.column()==11) return QVariant::fromValue<QList<double> >(targetInt.qCoefficients);    
  } else if (role==Qt::CheckStateRole && index.column()==0) {
    TargetIntData targetInt=targetIntList.at(index.row());
    if(targetInt.isActive==1) return Qt::Checked;
    else return Qt::Unchecked;
  } else if(role == Qt::TextAlignmentRole)  return Qt::AlignCenter;
  return QVariant();
}

QVariant TargetIntModel::headerData(int section, Qt::Orientation orientation, int role) const {
  if(role!=Qt::DisplayRole) return QVariant();
  if(orientation == Qt::Horizontal) {
    switch(section) {
    case 0:
      return tr("");
    case 1:
      return tr("Segment List");
    case 2:
      return tr("Number of Integration Points");
    case 3:
      return tr("Convolution Active?");
    case 4:
      return tr("Gaussian Sigma");
    case 5:
      return tr("Target Integration Active?");
    case 6:
      return tr("Target Density");
    case 7:
      return tr("Stopping Power Equation");
    case 8:
      return tr("Number of Parameters");
    case 9:
      return tr("Parameters List");
    case 10:
      return tr("Use Q-Coefficients?");
    case 11:
      return tr("Q-Coefficient List");
    default: 
      return QVariant();
    }
  } else if(orientation == Qt::Vertical) return section+1;
  return QVariant();
}

bool TargetIntModel::setData(const QModelIndex &index, const QVariant &value, int role) {
  if (index.isValid() && role == Qt::EditRole ) {
    int row = index.row();
    TargetIntData tempData = targetIntList.value(row);
    if(index.column() == 0) tempData.isActive = value.toInt();
    else if(index.column() == 1) tempData.segmentsList = value.toString();
    else if(index.column() == 2) tempData.numPoints = value.toInt();
    else if(index.column() == 3) tempData.isConvolution = value.toBool();
    else if(index.column() == 4) tempData.sigma = value.toDouble();
    else if(index.column() == 5) tempData.isTargetIntegration = value.toBool();
    else if(index.column() == 6) tempData.density = value.toDouble();
    else if(index.column() == 7) tempData.stoppingPowerEq = value.toString();
    else if(index.column() == 8) tempData.numParameters = value.toInt();
    else if(index.column() == 9) tempData.parameters = qVariantValue<QList<double> >(value);
    else if(index.column() == 10) tempData.isQCoefficients = value.toBool();
    else if(index.column() == 11) tempData.qCoefficients = qVariantValue<QList<double> >(value);
    else return false;
    targetIntList.replace(row,tempData);
    emit(dataChanged(index,index));
    return true;
  } else if(role == Qt::CheckStateRole) {
    int row = index.row();
    TargetIntData tempData = targetIntList.value(row);
    if(index.column() == 0) {
      if(value==Qt::Checked) tempData.isActive=1;
      else tempData.isActive=0;
    } else return false;
    targetIntList.replace(row,tempData);
    emit(dataChanged(index,index));
    return true;
  }
  return false;
}

bool TargetIntModel::insertRows(int position, int rows, const QModelIndex &index) {
  Q_UNUSED(index);
  if(rows>0) {
    beginInsertRows(QModelIndex(),position,position+rows-1);
    for(int row=0; row<rows; row++) {
     TargetIntData tempData;
      targetIntList.insert(position,tempData);
    }
    endInsertRows();
  }
  return true;
}

bool TargetIntModel::removeRows(int position, int rows, const QModelIndex &index) {
  Q_UNUSED(index);
  if(rows>0) {
    beginRemoveRows(QModelIndex(),position,position+rows-1);
    for(int row=0; row<rows;++row) {
      targetIntList.removeAt(position);
    }
    endRemoveRows();
  }
  return true;
}

Qt::ItemFlags TargetIntModel::flags(const QModelIndex &index) const {
  if (!index.isValid()) return Qt::ItemIsEnabled;
  if(index.column()==0) return QAbstractTableModel::flags(index) | Qt::ItemIsUserCheckable;
  return QAbstractTableModel::flags(index);
}
