#ifndef CHANNELRADIIMODEL_H
#define CHANNELRADIIMODEL_H

#include <QAbstractTableModel>
#include <QList>
#include <QString>
#include <functional>

class QTextStream;

/// One line of the <channelRadii> block: the channels (pair key, l) -- only channel spin s and/or total spin J when
/// those are set (>= 0) -- take \c radius (fm) instead of their pair's radius.
struct ChannelRadiusRow {
  int pairKey = 1;
  int l = 0;
  double s = -1.0;  ///< channel spin, or < 0 for any
  double J = -1.0;  ///< total spin, or < 0 for any
  double radius = 5.0;
};

/*!
 * The per-channel radii of a project (the <channelRadii> block after </levels>, 2026-11-10).  A channel takes its
 * pair's radius unless a row gives it its own; the most specific row wins (l+s+J > l+J > l+s > l), as in the engine
 * (PPair::GetChRad).  Columns: pair key, l, s (blank = any), J (blank = any), radius in fm.  Spins read and write as
 * "1/2", "3/2", "1" ...; the engine refuses what this model does not check (an unknown or non-particle pair, a channel
 * matched by an s row and a J row, external capture across two radii).
 */
class ChannelRadiiModel : public QAbstractTableModel {
  Q_OBJECT

 public:
  enum Column { PairCol = 0, LCol, SCol, JCol, RadiusCol, NumCols };
  explicit ChannelRadiiModel(QObject *parent = nullptr);

  int rowCount(const QModelIndex &parent = QModelIndex()) const override;
  int columnCount(const QModelIndex &parent = QModelIndex()) const override;
  QVariant data(const QModelIndex &index, int role) const override;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const override;
  Qt::ItemFlags flags(const QModelIndex &index) const override;
  bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole) override;

  const QList<ChannelRadiusRow> &rows() const { return rows_; }
  void setRows(const QList<ChannelRadiusRow> &rows);
  /// Append a row (pair key, l, any s, any J, \p radius); returns its index.
  int addRow(int pairKey, int l, double radius);
  bool removeRowAt(int row);
  void clear();
  bool isEmpty() const { return rows_.isEmpty(); }

  /// The radius the engine gives channel (pair key, l, s, J): the most specific matching row, else \p pairRadius.
  double radiusFor(int pairKey, int l, double s, double J, double pairRadius) const;
  /// The row that decides channel (pair key, l, s, J) -- the most specific match, as the engine -- or -1.
  int matchFor(int pairKey, int l, double s, double J) const;
  /// Give channel (pair key, l, s, J) the radius \p radius: through a row naming l, s and J (added, changed, or
  /// removed when the channel would take \p radius without it, from its pair or from a broader row).
  void setChannelRadius(int pairKey, int l, double s, double J, double radius, double pairRadius);

  /// Read the block body up to and including </channelRadii>; false on a malformed line or a missing end tag.
  bool readBlock(QTextStream &in);
  /// Write the body lines (without the tags).
  void writeBlock(QTextStream &out) const;

  /// "1/2" or "0.5" -> 0.5; false unless a non-negative multiple of 1/2.  An empty string or "any" -> -1.
  static bool parseSpin(const QString &text, double &value);
  /// 0.5 -> "1/2", 1 -> "1", < 0 -> "".
  static QString spinText(double value);
  /// One block line, "<pair> <l> <radius> [s=<s>] [J=<J>]".
  static QString lineFor(const ChannelRadiusRow &row);
  /// Parse one block line (comment stripped); false if malformed.
  static bool parseLine(const QString &line, ChannelRadiusRow &row);

 signals:
  void rowsEdited();

 private:
  QList<ChannelRadiusRow> rows_;
};

/// One particle channel of the calculation: pair key, l, channel spin s, and the J^pi of its level group.
struct ChannelKey {
  int pairKey = 1;
  int l = 0;
  double s = 0.0;
  double J = 0.0;
  int parity = 1;
};

/*!
 * Every particle channel of the calculation with the radius it takes (the Pairs tab's "Per-channel radii" box,
 * 2026-11-15): the pair's radius, or one of its own from the <channelRadii> rows (\c ChannelRadiiModel, which
 * stays the file's content).  Editing a radius writes a row naming the channel's l, s and J
 * (ChannelRadiiModel::setChannelRadius).  Columns: pair, l, s, J^pi, radius, source ("pair", "this channel", or the
 * broader row that decides it).  The channels come from the Levels tab (setChannels).
 */
class ChannelRadiiListModel : public QAbstractTableModel {
  Q_OBJECT

 public:
  enum Column { PairCol = 0, LCol, SCol, JCol, RadiusCol, SourceCol, NumCols };
  explicit ChannelRadiiListModel(ChannelRadiiModel *rows, QObject *parent = nullptr);

  int rowCount(const QModelIndex &parent = QModelIndex()) const override;
  int columnCount(const QModelIndex &parent = QModelIndex()) const override;
  QVariant data(const QModelIndex &index, int role) const override;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const override;
  Qt::ItemFlags flags(const QModelIndex &index) const override;
  bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole) override;

  /// The pair radius (fm) and a label for a pair key, from the Pairs table.
  void setPairInfo(std::function<double(int)> pairRadius, std::function<QString(int)> pairLabel);
  /// Replace the channel list (sorted by pair, J, parity, l, s; duplicates dropped).
  void setChannels(QList<ChannelKey> channels);
  const QList<ChannelKey> &channels() const { return channels_; }
  /// Give channel \p row back its pair's radius (removes its own row; a broader row still applies).
  void resetToPair(int row);
  /// Rows of the file that decide no channel of the list.
  int unusedRowCount() const;

 private:
  ChannelRadiiModel *rows_;
  QList<ChannelKey> channels_;
  std::function<double(int)> pairRadius_;
  std::function<QString(int)> pairLabel_;
  double pairRadiusOf(int pairKey) const { return pairRadius_ ? pairRadius_(pairKey) : 5.0; }
};

#endif
