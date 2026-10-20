#ifndef ABUNDANCESMODEL_H
#define ABUNDANCESMODEL_H

#include <QAbstractTableModel>
#include <QList>
#include <QMap>
#include <QStringList>
#include "PairsModel.h"

/*! One compound nucleus's isotopic abundance in one sample, as the <abundances>
 *  block stores it: the value, whether it is a fit parameter, and the width of
 *  its prior (0 = none).  The sample is "" for the unnamed (default) sample. */
struct AbundanceEntry {
  Nucleus nucleus;
  double value = 1.0;
  bool vary = false;
  double sigma = 0.0;
  QString sample;
};

/*!
 * Table model behind the "Isotopic abundances" table of the Particle Pairs tab:
 * one row per compound nucleus the pairs form, in order of first appearance,
 * for the sample currently shown (setCurrentSample).  A project may hold
 * several samples (a natural target and an enriched one, say): the unnamed
 * sample "" is what "natural 1" on a segment line means, a named one is
 * chosen with "natural <name>".  Each sample may carry its own sum-to-one
 * width (the block's "<name>: sum <sigma>" line; -1 = none, the config's
 * abundanceSumSigma applies).
 *
 * The rows follow the PairsModel (sync() on its row and data signals), so a
 * file load, the pair dialogs and a reset are all covered; an entry whose
 * nucleus disappears is kept for the session but not written.
 *
 * Columns: 0 nucleus (read-only), 1 abundance, 2 vary (check box), 3 prior sigma.
 * A project with a single nucleus and the unnamed sample at (1, fixed, 0) has
 * no block (hasBlock()), so such a file is unchanged.
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

  // -- samples --------------------------------------------------------------
  /*! The samples in file order, the unnamed one ("") first; never empty. */
  QStringList samples() const;
  /*! The sample the table shows. */
  QString currentSample() const { return current; }
  void setCurrentSample(const QString &sample);
  /*! Add a named sample (one token, no ':' or '#'); returns false for a bad or existing name. */
  bool addSample(const QString &sample);
  /*! Remove a named sample and its entries (the unnamed one stays). */
  bool removeSample(const QString &sample);
  /*! A sample's own sum-to-one width; -1 = none (the config's abundanceSumSigma applies). */
  double sampleSumSigma(const QString &sample) const;
  void setSampleSumSigma(const QString &sample, double sigma);

  // -- entries --------------------------------------------------------------
  /*! The entries of one sample in nucleus order. */
  QList<AbundanceEntry> entries(const QString &sample) const;
  /*! The entries of the unnamed sample (the single-sample case). */
  QList<AbundanceEntry> entries() const { return entries(QString()); }
  /*! Every sample's entries, sample by sample in file order: the engine's parameter order
   *  (abundance[_<sample>]_<Z>_<A>, samples as blocks of nuclei). */
  QList<AbundanceEntry> allEntries() const;
  /*! The entry of one nucleus in one sample (defaults when it has none). */
  AbundanceEntry entry(const QString &sample, const Nucleus &nucleus) const;
  AbundanceEntry entry(const Nucleus &nucleus) const { return entry(QString(), nucleus); }
  /*! Set an entry; a nucleus the pairs do not form is kept aside until they do; a new sample name
   *  creates the sample.  Returns false for a bad value. */
  bool setEntry(const QString &sample, const Nucleus &nucleus, double value, bool vary, double sigma);
  bool setEntry(const Nucleus &nucleus, double value, bool vary, double sigma) {
    return setEntry(QString(), nucleus, value, vary, sigma);
  }
  /*! The engine's parameter name of an entry. */
  static QString parameterName(const AbundanceEntry &e);
  /*! True when there is one sample, unnamed, with every nucleus at (1, fixed, 0) and no sum width. */
  bool isDefault() const;
  /*! Should the file carry an <abundances> block: more than one nucleus, a named sample, or a
   *  non-default entry.  The engine emits one abundance parameter per (sample, nucleus) exactly
   *  when the block is present. */
  bool hasBlock() const { return rows.size() > 1 || !isDefault(); }
  /*! Forget every entry and sample (a new file). */
  void clear();

 public slots:
  /*! Rebuild the rows from the PairsModel's nuclei, keeping the entries of surviving nuclei. */
  void sync();

 signals:
  /*! The sample list or the shown sample changed. */
  void samplesChanged();

 private:
  PairsModel *pairsModel = nullptr;
  QList<Nucleus> rows;
  QStringList sampleOrder;                              // "" first
  QMap<QString, QMap<Nucleus, AbundanceEntry> > stored;  // sample -> nucleus -> entry
  QMap<QString, double> sumSigma;                       // sample -> own width (absent = -1)
  QString current;
};

#endif
