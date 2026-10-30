#ifndef MASSCHECKDIALOG_H
#define MASSCHECKDIALOG_H

#include <QDialog>
#include <QList>

#include "PairsModel.h"

QT_BEGIN_NAMESPACE
class QCheckBox;
class QLabel;
class QTableWidget;
class QPushButton;
QT_END_NAMESPACE

/*!
 * "Masses from AME2020...": every pair of the project against the mass table (MassTable), with a tick box per
 * pair.  Accepting fills the ticked pairs' masses -- and, with "Also set the separation energies", their S --
 * from the table; filledPairs() returns the new pair data by row.
 *
 * A pair is ticked by default when the table can fill it and its values differ from the table's.  Moving S moves
 * that pair's threshold, so the dialog says so; a pair the table cannot identify is shown but cannot be ticked.
 */
class MassCheckDialog : public QDialog {
  Q_OBJECT

 public:
  MassCheckDialog(const PairsModel *model, QWidget *parent = nullptr);
  /*! The pairs to change, (row, new data), for the ticked rows whose values the fill actually changes. */
  QList<QPair<int, PairsData>> filledPairs() const;
  /*! The pair data of one row as the fill would leave it (unchanged when the table cannot fill it). */
  PairsData filled(int row, bool separation) const;
  QTableWidget *table() const { return table_; }
  QCheckBox *separationCheck() const { return separationCheck_; }

 private:
  QList<PairsData> pairs_;
  QTableWidget *table_;
  QCheckBox *separationCheck_;
  QLabel *summary_;
};

#endif
