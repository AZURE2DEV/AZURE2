#ifndef CHANNELSMODEL_H
#define CHANNELSMODEL_H

#include <QAbstractTableModel>
#include <QList>

struct ChannelsData {
  static const int SIZE = 7;
  int isFixed;
  int levelIndex;
  int pairIndex;
  double sValue;
  int lValue;
  QChar radType;
  double reducedWidth;
};

class PairsModel;

class ChannelsModel : public QAbstractTableModel {
  Q_OBJECT

 public:
  ChannelsModel(QObject *parent = 0);

  int rowCount(const QModelIndex &parent) const;
  int columnCount(const QModelIndex &parent) const;
  QVariant data(const QModelIndex &index, int role) const;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const;
  QList<ChannelsData> getChannels() const {return channelsList;};
  bool setData(const QModelIndex &index, const QVariant &value, int role=Qt::EditRole);
  bool insertRows(int position, int rows, const QModelIndex &index=QModelIndex());
  bool removeRows(int position, int rows, const QModelIndex &index=QModelIndex());
  Qt::ItemFlags flags(const QModelIndex &index) const;
  bool isChannel(const ChannelsData &channel) const;
  QString getSpinLabel(const ChannelsData &channel) const;
  void setPairsModel(PairsModel *model);

 private:
  QList<ChannelsData> channelsList;  
  PairsModel *pairsModel;
};

#endif
