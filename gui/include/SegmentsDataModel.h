#ifndef SEGMENTSDATAMODEL_H
#define SEGMENTSDATAMODEL_H

#include <QAbstractTableModel>
#include <QList>

class PairsModel;

struct SegmentsDataData {
  static const int SIZE = 24;
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
  int isAdvanced;
  int operationType;
  QString componentsList;
  int isUPOS;
  int secondaryL;
  double finalJ;
  double delta;
};

class SegmentsDataModel : public QAbstractTableModel {
  Q_OBJECT

 public:
  SegmentsDataModel(QObject *parent = 0);

  int rowCount(const QModelIndex &parent) const;
  int columnCount(const QModelIndex &parent) const;
  QVariant data(const QModelIndex &index, int role) const;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const;
  bool setData(const QModelIndex &index, const QVariant &value, int role=Qt::EditRole);
  bool insertRows(int position, int rows, const QModelIndex &index=QModelIndex());
  bool removeRows(int position, int rows, const QModelIndex &index=QModelIndex());
  Qt::ItemFlags flags(const QModelIndex &index) const;
  int isSegDataLine(const SegmentsDataData &line) const;
  QList<SegmentsDataData> getLines() const {return segDataLineList;};
  void setPairsModel(PairsModel* model);
  QString getReactionLabel(const QModelIndex &index);

 signals:
  void normalizationChanged(int segmentIndex, double value);
  void energyShiftChanged(int segmentIndex, double value);
  void normalizationErrorChanged(int segmentIndex, double error);
  void energyShiftErrorChanged(int segmentIndex, double error);
  void normalizationVaryChanged(int segmentIndex, bool vary);
  void energyShiftVaryChanged(int segmentIndex, bool vary);

 private:
  QList<SegmentsDataData> segDataLineList;  
  PairsModel* pairsModel;
};


#endif
