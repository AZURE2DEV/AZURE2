#ifndef CHANNELSMODEL_H
#define CHANNELSMODEL_H

#include <QAbstractTableModel>
#include <QStringList>
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
  /// The <levels> row the channel was read from, its tokens by column (azure2::front::LevelRow), so that a value
  /// the GUI does not change is written back as it was written; empty for a channel made in the GUI.
  QStringList fileTokens = QStringList();
};

class PairsModel;

/*!
 * Table model for the channels of the selected level: pair, l, s and the reduced width or ANC.
 */
class ChannelsModel : public QAbstractTableModel {
  Q_OBJECT

 public:
  /// The file tokens of a channel (ChannelsData::fileTokens).
  void setFileTokens(int row, const QStringList &tokens);
  ChannelsModel(QObject *parent = 0);

  int rowCount(const QModelIndex &parent) const;
  int columnCount(const QModelIndex &parent) const;
  QVariant data(const QModelIndex &index, int role) const;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const;
  QList<ChannelsData> getChannels() const { return channelsList; };
  bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole);
  bool insertRows(int position, int rows, const QModelIndex &index = QModelIndex());
  bool removeRows(int position, int rows, const QModelIndex &index = QModelIndex());
  Qt::ItemFlags flags(const QModelIndex &index) const;
  bool isChannel(const ChannelsData &channel) const;
  QString getSpinLabel(const ChannelsData &channel) const;
  void setPairsModel(PairsModel *model);

 signals:
  /// A channel's width or ANC (column 6) was set.
  void widthEdited();

 private:
  QList<ChannelsData> channelsList;
  PairsModel *pairsModel;
};

#endif
