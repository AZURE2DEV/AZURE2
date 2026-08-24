#ifndef PAIRSMODEL_H
#define PAIRSMODEL_H

#include <QAbstractTableModel>
#include <QList>

struct PairsData {
  static const int SIZE = 16;
  double lightJ;
  int lightPi;
  int lightZ;
  double lightM;
  double lightG;
  double heavyJ;
  int heavyPi;
  int heavyZ;
  double heavyM;
  double heavyG;
  double excitationEnergy;
  double seperationEnergy;
  double channelRadius;
  int pairType;
  int ecMultMask;
  double bindingEnergy;  // THM: transferred-particle binding energy (MeV)
};

/*!
 * Table model behind the Particle Pairs tab.
 *
 * A pair's key is its 1-based row, which is what the .azr stores and what a channel's pair number refers to -- so deleting a row renumbers the pairs above it.
 */
class PairsModel : public QAbstractTableModel {
  Q_OBJECT

 public:
  PairsModel(QObject *parent = 0);

  int rowCount(const QModelIndex &parent) const;
  int columnCount(const QModelIndex &parent) const;
  QVariant data(const QModelIndex &index, int role) const;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const;
  bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole);
  bool insertRows(int position, int rows, const QModelIndex &index = QModelIndex());
  bool removeRows(int position, int rows, const QModelIndex &index = QModelIndex());

  int isPair(const PairsData &pair) const;
  int numPairs() const { return pairsList.size(); };
  QList<PairsData> getPairs() const { return pairsList; };
  QString getParticleLabel(const PairsData &pair, int which = -1) const;
  QString getReactionLabel(const PairsData &firstPair, const PairsData &secondPair);
  QString getReactionLabelTotalCapture(const PairsData &firstPair);
  QString getSpinLabel(const PairsData &pair, int which) const;

 private:
  QList<PairsData> pairsList;
};

#endif
