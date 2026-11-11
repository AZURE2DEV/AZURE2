#ifndef SEGMENTSDATAMODEL_H
#define SEGMENTSDATAMODEL_H

#include <QAbstractTableModel>
#include <QList>

class PairsModel;

struct SegmentsDataData {
  /// Columns 0-23 are the fields below; column 24 (SegmentsDataModel::GroupsColumn) is derived from extraTokens.
  static const int SIZE = 25;
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
  // UPOS (Unobserved Primary, Observed Secondary) fields
  int isUPOS;
  int secondaryDecayL;
  double finalJ;
  double delta;
  // Trailing keyword blocks of the line, verbatim: "sqrtshift <b> <bError>
  // <vary>" (sqrt(E) energy-shift term, SegLine.h) and "mixture <w_self> <n>
  // {"file" <w>}..." (natural-abundance targets, MixtureTable.h).  Carried
  // through so a file that uses them survives a GUI load/save.  The data
  // segment dialog edits the sample, covariance, spin and group blocks
  // (SegmentKeywords.h) and leaves every other token as read.  Set with
  // setLineExtras(); the Groups column (read-only) shows the group blocks.
  QString extraTokens;
  // isDiff 10: the spin-observable token (SpinObservable.h).  Not a table column: set with setLineSpin().
  QString spinObservable;
};

/*!
 * Table model for <segmentsData> -- the segments that carry a data file and contribute to chi-squared.
 */
class SegmentsDataModel : public QAbstractTableModel {
  Q_OBJECT

 public:
  SegmentsDataModel(QObject *parent = 0);

  int rowCount(const QModelIndex &parent) const;
  int columnCount(const QModelIndex &parent) const;
  QVariant data(const QModelIndex &index, int role) const;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const;
  bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole);
  bool insertRows(int position, int rows, const QModelIndex &index = QModelIndex());
  bool removeRows(int position, int rows, const QModelIndex &index = QModelIndex());
  /// Relocate one row wholesale.  The whole struct moves, so fields added
  /// later can never be forgotten the way a per-column copy forgets them.
  bool moveLine(int from, int to);
  /// Store the verbatim keyword blocks of a row (no editable column addresses them).
  void setLineExtras(int row, const QString &tokens);
  void setLineSpin(int row, const QString &token);
  /// The read-only column listing a row's groups, "norm: D, shift: A" (normgroup / shiftgroup / shiftfngroup).
  static const int GroupsColumn = 24;
  /*! Row of the segment whose \p keyword ("normgroup", "shiftgroup", "shiftfngroup") value \p row takes: the first
   *  ACTIVE row carrying the same label (EData counts active segments only), when that is another row and \p row
   *  is active itself; else -1 (no group, the group's first segment, or an inactive row). */
  int groupLeader(int row, const QString &keyword) const;
  /// The labels the rows use with \p keyword, each once, sorted.
  QStringList groupLabels(const QString &keyword) const;
  /// "norm: D, shift: A, shift fn: F" for the row's group blocks; "" when none.
  QString groupsText(int row) const;
  /// One line per group block of the row: whose value it shares, or which rows share its value.
  QString groupsToolTip(int row) const;
  Qt::ItemFlags flags(const QModelIndex &index) const;
  int isSegDataLine(const SegmentsDataData &line) const;
  QList<SegmentsDataData> getLines() const { return segDataLineList; };
  void setPairsModel(PairsModel *model);
  QString getReactionLabel(const QModelIndex &index);

 signals:
  void normalizationChanged(int segmentIndex, double value);
  void energyShiftChanged(int segmentIndex, double value);
  void normalizationErrorChanged(int segmentIndex, double error);
  void energyShiftErrorChanged(int segmentIndex, double error);
  void normalizationVaryChanged(int segmentIndex, bool vary);
  void energyShiftVaryChanged(int segmentIndex, bool vary);

 private:
  /// dataChanged over the whole table when any row is in a group (activating a row can change who leads it).
  void emitGroupCellsChanged();
  QList<SegmentsDataData> segDataLineList;
  PairsModel *pairsModel;
};


#endif
