#ifndef CHANNELRADIIMODEL_H
#define CHANNELRADIIMODEL_H

#include <QAbstractTableModel>
#include <QList>
#include <QString>

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

#endif
