#ifndef SEGMENTSTESTMODEL_H
#define SEGMENTSTESTMODEL_H

#include <QAbstractTableModel>
#include <QList>

class PairsModel;

struct SegmentsTestData {
  static const int SIZE = 17;  // last column (16) is isTHM
  int isActive;
  int entrancePairIndex;
  int exitPairIndex;
  double lowEnergy;
  double highEnergy;
  double energyStep;
  double lowAngle;
  double highAngle;
  double angleStep;
  int dataType;
  double phaseJ;
  int phaseL;
  int maxAngDistOrder;
  int isAdvanced;
  int operationType;
  QString componentsList;
  int isTHM;  // THM (modified R-matrix) flag; written as isDiff offset +10.
};

class SegmentsTestModel : public QAbstractTableModel {
  Q_OBJECT

 public:
  SegmentsTestModel(QObject *parent = 0);

  int rowCount(const QModelIndex &parent) const;
  int columnCount(const QModelIndex &parent) const;
  QVariant data(const QModelIndex &index, int role) const;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const;
  bool setData(const QModelIndex &index, const QVariant &value, int role=Qt::EditRole);
  bool insertRows(int position, int rows, const QModelIndex &index=QModelIndex());
  bool removeRows(int position, int rows, const QModelIndex &index=QModelIndex());
  Qt::ItemFlags flags(const QModelIndex &index) const;
  int isSegTestLine(const SegmentsTestData &line) const;
  QList<SegmentsTestData> getLines() const {return segTestLineList;};
  void setPairsModel(PairsModel* model);
  QString getReactionLabel(const QModelIndex &index);
 private:
  QList<SegmentsTestData> segTestLineList;  
  PairsModel* pairsModel;
};


#endif
