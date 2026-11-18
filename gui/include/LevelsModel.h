#ifndef LEVELSMODEL_H
#define LEVELSMODEL_H

#include <QAbstractTableModel>
#include <QList>
#include "PairsModel.h"

struct LevelsData {
  static const int SIZE = 6;
  int isActive;
  int isFixed;
  double jValue;
  int piValue;
  double energy;
  // The compound nucleus the level belongs to (PairsModel::nucleusOf of its
  // channels' pairs).  (0, 0) = not assigned: the level takes channels from
  // every pair and matches a level of any nucleus, which is the behaviour of
  // a single-nucleus project and of the older tests.
  int nucZ = 0;
  int nucA = 0;
  // The level ID the line carried in the file.  An isospin level (<isospin>
  // block) is one ID whose lines span several nuclei; the GUI holds one level
  // per nucleus and uses this to write them back as one ID.  0 = not from a file.
  int fileId = 0;
  Nucleus nucleus() const { return {nucZ, nucA}; }
  bool hasNucleus() const { return nucZ > 0 || nucA > 0; }
};

/*!
 * Table model behind the Levels tab: one row per level, with its J^pi, energy and fixed flag.
 */
class LevelsModel : public QAbstractTableModel {
  Q_OBJECT

 public:
  LevelsModel(QObject *parent = 0);

  int rowCount(const QModelIndex &parent) const;
  int columnCount(const QModelIndex &parent) const;
  QVariant data(const QModelIndex &index, int role) const;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const;
  bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole);
  bool insertRows(int position, int rows, const QModelIndex &index = QModelIndex());
  bool removeRows(int position, int rows, const QModelIndex &index = QModelIndex());
  Qt::ItemFlags flags(const QModelIndex &index) const;
  /*! Row of a level with the same J, parity, energy and nucleus (an unassigned nucleus matches any), or -1. */
  int isLevel(const LevelsData &level) const;
  QList<LevelsData> getLevels() const { return levelsList; };
  QString getSpinLabel(const LevelsData &level) const;
  /*! The read-only column showing the level's compound nucleus (hidden by the tab with one nucleus). */
  static const int NucleusColumn = 5;
  void setNucleus(int row, const Nucleus &nucleus);
  /*! The level ID the level's lines carried in the file (LevelsData::fileId). */
  void setFileId(int row, int fileId);

 private:
  QList<LevelsData> levelsList;
};

#endif
