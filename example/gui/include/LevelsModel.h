#ifndef LEVELSMODEL_H
#define LEVELSMODEL_H

#include <QAbstractTableModel>
#include <QList>

struct LevelsData {
  static const int SIZE = 5;
  int isActive;
  int isFixed;
  double jValue;
  int piValue;
  double energy;
};

class LevelsModel : public QAbstractTableModel {
  Q_OBJECT

 public:
  LevelsModel(QObject *parent = 0);

  int rowCount(const QModelIndex &parent) const;
  int columnCount(const QModelIndex &parent) const;
  QVariant data(const QModelIndex &index, int role) const;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const;
  bool setData(const QModelIndex &index, const QVariant &value, int role=Qt::EditRole);
  bool insertRows(int position, int rows, const QModelIndex &index=QModelIndex());
  bool removeRows(int position, int rows, const QModelIndex &index=QModelIndex());
  Qt::ItemFlags flags(const QModelIndex &index) const;
  int isLevel(const LevelsData &level) const;
  QList<LevelsData> getLevels() const {return levelsList;};
  QString getSpinLabel(const LevelsData &level) const;

 private:
  QList<LevelsData> levelsList;  
};

#endif
