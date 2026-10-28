#ifndef SCATTERINGLENGTHSMODEL_H
#define SCATTERINGLENGTHSMODEL_H

#include <QAbstractTableModel>
#include <QList>
#include <QStyledItemDelegate>
#include <QTextStream>
#include "PairsModel.h"

/*! One line of the <scatteringLengths> block: a measured s-wave scattering length (fm) of a particle pair in one
 *  channel spin s (J = s, parity pi1*pi2), or the coherent spin-weighted combination (s = -1). A sigma of 0
 *  leaves that part unconstrained; the imaginary part is optional on the line (hasIm). */
struct ScatteringLengthEntry {
  bool active = true;
  int pairKey = 1;
  double s = 0.0;
  double value = 0.0;
  double sigma = 0.0;
  double imValue = 0.0;
  double imSigma = 0.0;
  bool hasIm = false;
};

/*!
 * Table model behind the "Measured Scattering Lengths" table of the Segments tab: one row per line of the
 * optional <scatteringLengths> block (see src/ScatteringLength.cpp, ReadConstraints).
 *
 * Rows refer to a pair by its key, the 1-based row of the Particle Pairs tab. They are never dropped when the
 * pairs change: a row whose pair is missing, is not a particle pair, or does not allow the row's channel spin
 * is shown in red with the reason as its tooltip (problem()), and is written as it stands -- the engine
 * reports the error.
 *
 * Columns: 0 active (check box), 1 pair, 2 channel spin, 3 value, 4 sigma, 5 Im value, 6 Im sigma (fm).
 */
class ScatteringLengthsModel : public QAbstractTableModel {
  Q_OBJECT

 public:
  enum { ActiveCol = 0, PairCol = 1, SpinCol = 2, ValueCol = 3, SigmaCol = 4, ImValueCol = 5, ImSigmaCol = 6, NumCols = 7 };

  ScatteringLengthsModel(QObject *parent = 0);
  void setPairsModel(PairsModel *model);
  PairsModel *getPairsModel() const { return pairsModel; }

  int rowCount(const QModelIndex &parent = QModelIndex()) const;
  int columnCount(const QModelIndex &parent = QModelIndex()) const;
  QVariant data(const QModelIndex &index, int role) const;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const;
  bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole);
  Qt::ItemFlags flags(const QModelIndex &index) const;
  bool removeRows(int position, int rows, const QModelIndex &parent = QModelIndex());

  QList<ScatteringLengthEntry> entries() const { return lines; }
  void addEntry(const ScatteringLengthEntry &entry);
  /*! Replaces row \p row (the Edit dialog). */
  void setEntry(int row, const ScatteringLengthEntry &entry);
  /*! "key: pair label" of a pair key, as the Pair column shows it. */
  QString pairText(int key) const { return pairLabel(key); }
  /*! A new row on the first particle pair with its first allowed channel spin, nothing constrained yet. */
  void addDefaultEntry();
  void clear();

  /*! Reads the block's lines up to </scatteringLengths> (the opening tag already read). Comments and blank
   *  lines are skipped and not kept. False on a malformed line. */
  bool read(QTextStream &in);
  /*! Writes the lines (no tags), a header comment first. */
  bool write(QTextStream &out) const;

  /*! The engine's identical-pair rule (PPair): same Z, mass, spin and parity, light partner and residual in
   *  their ground states, particle pair. */
  static bool isIdentical(const PairsData &pair);
  /*! The channel spins an s wave of the pair has: |j1-j2| .. j1+j2, only even ones for an identical pair. */
  static QList<double> allowedChannelSpins(const PairsData &pair);
  /*! "1", "1/2", "coherent" (s < 0). */
  static QString spinLabel(double s);
  /*! Why the engine would refuse the row ("" when it would not). */
  QString problem(int row) const;
  QString problem(const ScatteringLengthEntry &entry) const;

 public slots:
  /*! The pairs changed: re-evaluate the labels and flags of every row. */
  void pairsChanged();

 private:
  QString pairLabel(int key) const;
  QList<ScatteringLengthEntry> lines;
  PairsModel *pairsModel = nullptr;
};

/*!
 * The scattering lengths as rows under the Segments From Data table (SegmentsTab): the columns of
 * SegmentsDataModel, filled for a scattering length -- Reaction "a0 of <pair>", Energy Range "threshold", Angle
 * Range the channel spin, Data Type "Scattering Length", Data File the measured value(s).  Column 0 is the active
 * check box.  Reads and writes through the ScatteringLengthsModel, which stays the <scatteringLengths> block.
 */
class ScatteringLengthRowsModel : public QAbstractTableModel {
  Q_OBJECT

 public:
  ScatteringLengthRowsModel(ScatteringLengthsModel *source, int columns, QObject *parent = 0);
  int rowCount(const QModelIndex &parent = QModelIndex()) const override;
  int columnCount(const QModelIndex &parent = QModelIndex()) const override;
  QVariant data(const QModelIndex &index, int role) const override;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const override;
  bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole) override;
  Qt::ItemFlags flags(const QModelIndex &index) const override;
  /*! Text of the Data File column: "Re a0 = v +- s fm", and the Im part when given. */
  static QString valueText(const ScatteringLengthEntry &e);

 private slots:
  void sourceAboutToBeReset() { beginResetModel(); }
  void sourceReset() { endResetModel(); }
  void sourceRowsAboutToBeInserted(const QModelIndex &, int first, int last) { beginInsertRows(QModelIndex(), first, last); }
  void sourceRowsInserted() { endInsertRows(); }
  void sourceRowsAboutToBeRemoved(const QModelIndex &, int first, int last) { beginRemoveRows(QModelIndex(), first, last); }
  void sourceRowsRemoved() { endRemoveRows(); }
  void sourceDataChanged(const QModelIndex &topLeft, const QModelIndex &bottomRight);
  void sourceLayoutAboutToBeChanged() { emit layoutAboutToBeChanged(); }
  void sourceLayoutChanged() { emit layoutChanged(); }

 private:
  ScatteringLengthsModel *src;
  int ncols;
};

/*!
 * Editors of the scattering-length table: a drop-down of the particle pairs, one of the channel spins the
 * row's pair allows (plus "coherent" when the pair is not identical), and a line edit for the numbers.
 */
class ScatteringLengthsDelegate : public QStyledItemDelegate {
  Q_OBJECT

 public:
  ScatteringLengthsDelegate(QObject *parent = 0) : QStyledItemDelegate(parent) {}
  QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const override;
  void setEditorData(QWidget *editor, const QModelIndex &index) const override;
  void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const override;
};

#endif
