#ifndef TARGETINTMODEL_H
#define TARGETINTMODEL_H

#include <QAbstractTableModel>
#include <QList>

Q_DECLARE_METATYPE(QList<double>);

struct TargetIntData {
  static const int SIZE = 12;
  int isActive;
  QString segmentsList;
  int numPoints;
  bool isConvolution;
  double sigma;
  bool isTargetIntegration;
  double density;
  QString stoppingPowerEq;
  int numParameters;
  QList<double> parameters;
  bool isQCoefficients;
  QList<double> qCoefficients;
};

class TargetIntModel : public QAbstractTableModel {
 Q_OBJECT

 public:
  TargetIntModel(QObject *parent = 0);
  int rowCount(const QModelIndex &parent) const;
  int columnCount(const QModelIndex &parent) const;
  QVariant data(const QModelIndex &index, int role) const;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const;
  bool setData(const QModelIndex &index, const QVariant &value, int role=Qt::EditRole);
  bool insertRows(int position, int rows, const QModelIndex &index=QModelIndex());
  bool removeRows(int position, int rows, const QModelIndex &index=QModelIndex());
  Qt::ItemFlags flags(const QModelIndex &index) const;
  QList<TargetIntData> getLines() const {return targetIntList;};
 private:
  QList<TargetIntData> targetIntList;
};

#endif
