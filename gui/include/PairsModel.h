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
 * A compound nucleus, as the engine infers it from a pair's columns (PPair):
 * Z1+Z2 / round(M1+M2) for a particle pair, Z2 / round(M2) for a gamma pair,
 * Z2-Z1 / round(M2) for a beta pair.  Pairs with the same (Z, A) form one
 * nucleus; a project may hold several.
 */
struct Nucleus {
  int Z = 0;
  int A = 0;
  bool operator==(const Nucleus &o) const { return Z == o.Z && A == o.A; }
  bool operator!=(const Nucleus &o) const { return !(*this == o); }
  bool operator<(const Nucleus &o) const { return Z != o.Z ? Z < o.Z : A < o.A; }
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
  /*! The pair as plain text, "11B+a [0.000 MeV]" with a superscript mass number, for widgets that render no HTML. */
  QString getPairLabel(const PairsData &pair) const;
  /*! An HTML label of this model (getParticleLabel, getReactionLabel) as plain text: tags dropped, entities and
   *  <sup> digits and signs turned into their Unicode characters. */
  static QString plainText(const QString &html);

  /*! The read-only column showing the compound nucleus a pair forms (after the PairsData columns). */
  static const int NucleusColumn = PairsData::SIZE;
  /*! The compound nucleus a pair forms (the engine's rule, PPair). */
  static Nucleus nucleusOf(const PairsData &pair);
  /*! The distinct compound nuclei of the pairs, in order of first appearance. */
  QList<Nucleus> nuclei() const;
  /*! 0-based position of the nucleus in nuclei(), or -1. */
  int nucleusIndex(const Nucleus &nucleus) const;
  /*! "<sup>23</sup>Ne" (html) or the same as plain text for the compound nucleus. */
  static QString nucleusLabel(const Nucleus &nucleus, bool html = true);
  /*! "23Ne (22Ne+n)" as plain text: the compound nucleus and the first particle pair that forms it -- the
   *  target isotope an abundance refers to is the heavy particle of that pair. */
  QString nucleusLongLabel(const Nucleus &nucleus) const;

 private:
  QList<PairsData> pairsList;
};

#endif
