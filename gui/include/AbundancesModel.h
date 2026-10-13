#ifndef ABUNDANCESMODEL_H
#define ABUNDANCESMODEL_H

#include <QAbstractTableModel>
#include <QList>
#include <QMap>
#include "PairsModel.h"

/*! One compound nucleus's isotopic abundance as the <abundances> block stores it:
 *  the value, whether it is a fit parameter, and the width of its prior (0 = none). */
struct AbundanceEntry {
  Nucleus nucleus;
  double value = 1.0;
  bool vary = false;
  double sigma = 0.0;
};

/*!
 * Table model behind the "Isotopic abundances" table of the Particle Pairs tab:
 * one row per compound nucleus the pairs form, in order of first appearance.
 * The rows follow the PairsModel (sync() on its row and data signals), so a
 * file load, the pair dialogs and a reset are all covered; an entry whose
 * nucleus disappears is kept for the session but not written.
 *
 * Columns: 0 nucleus (read-only), 1 abundance, 2 vary (check box), 3 prior sigma.
 * A project with a single nucleus has one row at (1, fixed, 0), which is not
 * written as a block (isDefault()), so such a file is unchanged.
 */
class AbundancesModel : public QAbstractTableModel {
  Q_OBJECT

 public:
  enum { NucleusCol = 0, ValueCol = 1, VaryCol = 2, SigmaCol = 3, NumCols = 4 };

  AbundancesModel(QObject *parent = 0);
  void setPairsModel(PairsModel *model);

  int rowCount(const QModelIndex &parent = QModelIndex()) const;
  int columnCount(const QModelIndex &parent = QModelIndex()) const;
  QVariant data(const QModelIndex &index, int role) const;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const;
  bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole);
  Qt::ItemFlags flags(const QModelIndex &index) const;

  /*! The entries in nucleus order. */
  QList<AbundanceEntry> entries() const;
  /*! The entry of one nucleus (defaults when it has none). */
  AbundanceEntry entry(const Nucleus &nucleus) const;
  /*! Set an entry; a nucleus the pairs do not form is kept aside until they do. Returns false for a bad value. */
  bool setEntry(const Nucleus &nucleus, double value, bool vary, double sigma);
  /*! True when every nucleus is at (1, fixed, 0): then no <abundances> block is written. */
  bool isDefault() const;
  /*! Should the file carry an <abundances> block: more than one nucleus, or a non-default entry. The engine
   *  emits one abundance parameter per nucleus exactly when the block is present. */
  bool hasBlock() const { return rows.size() > 1 || !isDefault(); }
  /*! Forget every entry (a new file). */
  void clear();

 public slots:
  /*! Rebuild the rows from the PairsModel's nuclei, keeping the entries of surviving nuclei. */
  void sync();

 private:
  PairsModel *pairsModel = nullptr;
  QList<Nucleus> rows;
  QMap<Nucleus, AbundanceEntry> stored;
};

#endif
