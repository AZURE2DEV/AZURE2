#ifndef MASSTABLE_H
#define MASSTABLE_H

#include <QHash>
#include <QString>

struct PairsData;

/*!
 * One nuclide of the atomic mass table.
 */
struct MassEntry {
  int Z = 0;
  int A = 0;
  QString symbol;
  double massExcess = 0.0;       //!< keV
  double massExcessUnc = 0.0;    //!< keV
  double atomicMass = 0.0;       //!< u, atomic (electrons included)
  double atomicMassUnc = 0.0;    //!< u
  bool estimated = false;        //!< the evaluators' extrapolation ('#' in the file), not a measurement
  /*! "12C"; the neutron is "n". */
  QString label() const { return Z == 0 ? QString("n") : QString("%1%2").arg(A).arg(symbol); }
};

/*!
 * What the mass table says about one particle pair (MassTable::check).
 *
 * A particle is identified by its charge and its mass rounded to a mass number, the engine's rule, so a pair
 * typed with mass numbers or a bare-proton mass is still recognised.  hasS is set for particle pairs whose three
 * nuclides (light, heavy, compound) are all in the table; a gamma pair's 0 and a beta pair's value are
 * conventions and are not compared.
 */
struct PairMassCheck {
  const MassEntry *light = nullptr;   //!< nullptr for a photon, an electron or a nuclide the table lacks
  const MassEntry *heavy = nullptr;
  bool hasS = false;
  double sTable = 0.0;                //!< MeV
  double dSkeV = 0.0;                 //!< the pair's S minus the table's, keV
  /*! Can the table fill this pair (both particles of a particle pair, the heavy one of a gamma or beta pair)? */
  bool fillable = false;
  /*! Largest |pair mass - table mass| over the identified particles, u. */
  double massOffset = 0.0;
};

/*!
 * The AME2020 atomic mass table, read from the copy of pyazr/data/mass_1.mas20 the GUI embeds as a Qt resource
 * (":/ame/mass_1.mas20"), so the GUI and pyazr fill masses and separation energies from identical numbers.
 *
 * A pair's separation energy is S = D1 + D2 - Dc from the mass excesses (the electrons of the atomic masses
 * cancel), the value AZURE2 compares with m1 + m2 - M when it reads a project
 * (CNuc::CheckMassSeparationConsistency).  It does not depend on the excitation energy of the residual nucleus:
 * an excited-state pair's threshold is S + E_x.
 */
class MassTable {
 public:
  static const MassTable &instance();
  static QString version() { return QStringLiteral("AME2020"); }
  /*! The agreement within which a pair's S counts as consistent, keV (the engine's warning threshold). */
  static constexpr double kToleranceKeV = 1.0;

  bool isLoaded() const { return !entries_.isEmpty(); }
  int size() const { return entries_.size(); }
  /*! The row of (Z, A), or nullptr. */
  const MassEntry *lookup(int Z, int A) const;
  /*! The row a pair column describes: Z and the mass rounded to a mass number; nullptr for A <= 0 or Z < 0. */
  const MassEntry *identify(int Z, double mass) const;
  /*! (Z, A) of "p", "n", "d", "t", "a"/"alpha", "helion", "12C", "C12", "c-12" or "12-C"; false if unreadable.
   *  "h" is refused: it reads as hydrogen as easily as the helion. */
  bool parseNuclide(const QString &text, int *Z, int *A) const;
  /*! S of (Z1, A1) from the compound nucleus (Z1+Z2, A1+A2), MeV; false when a nuclide is missing. */
  bool separationEnergy(int Z1, int A1, int Z2, int A2, double *S) const;
  /*! The pair against the table. */
  PairMassCheck check(const PairsData &pair) const;

 private:
  MassTable();
  QHash<int, MassEntry> entries_;      // key Z * 1000 + A
  QHash<QString, int> symbolToZ_;      // lower-case symbol -> Z
};

#endif
