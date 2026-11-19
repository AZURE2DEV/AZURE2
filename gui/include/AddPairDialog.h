#ifndef ADDPAIRDIALOG_H
#define ADDPAIRDIALOG_H

#include <QDialog>

#include "PairsModel.h"

QT_BEGIN_NAMESPACE

class QLineEdit;
class QLabel;
class QCheckBox;
class QGroupBox;
class QComboBox;
class QDoubleSpinBox;
class QPushButton;
class QSpinBox;
class QToolButton;
class QWidget;

QT_END_NAMESPACE

/*!
 * Dialog for adding, editing or duplicating one particle pair (docs/upgrade/qt_gui.md section 15.6).
 *
 * Each particle is chosen as dev's dialog chose it, by its nucleus: the element by its symbol (a drop-down of the
 * mass table's elements, or Z typed) and the mass number A by a spin box; Z and the atomic mass are filled from the
 * table (AME2020), and J and the parity where they are known without a spin table (n, p, d, t, 3He, alpha, and the
 * 0+ of an even-even nucleus); every value stays editable.  J is a spin box of halves, the parity a +/- toggle.  The
 * pair type is a drop-down that hides what that type does not use (a photon pair: no light particle, separation
 * energy or radius; a beta pair: the light particle's charge only).
 */
class AddPairDialog : public QDialog {
  Q_OBJECT

 public:
  AddPairDialog(QWidget *parent = 0);
  /// The pair as the dialog holds it.
  PairsData pair() const;
  /// Show a pair (its g-factors are kept: the dialog does not edit them).
  void setPair(const PairsData &pair);

  QComboBox *pairTypeCombo;
  /// One particle's widgets: nucleus (element and A), Z, mass, J and parity.
  struct Particle {
    QComboBox *element = nullptr;
    QSpinBox *massNumber = nullptr;
    QLineEdit *z = nullptr;
    QLineEdit *mass = nullptr;
    QDoubleSpinBox *j = nullptr;
    QToolButton *parity = nullptr;  ///< checked: +
    QLabel *note = nullptr;
    QWidget *nucleusRow = nullptr;
  };
  Particle light, heavy;
  QLineEdit *excitationEnergyText;
  /*! The light particle's excitation energy [MeV] (column 32 of a <levels> line); particle pairs only. */
  QLineEdit *lightExcitationText;
  QLineEdit *seperationEnergyText;
  QLineEdit *channelRadiusText;
  QCheckBox *e1Check;
  QCheckBox *e2Check;
  QGroupBox *multBox;
  /*! Fills both masses and the separation energy from the mass table (MassTable, AME2020). */
  QPushButton *fillButton;
  /*! The entered masses and S against the mass table, updated as they are typed. */
  QLabel *massCheckLabel;

  /// Choose particle `p`'s nucleus as the user does (element symbol and A): Z, mass, J and parity filled.
  bool chooseNucleus(Particle &p, const QString &symbol, int A);

 public slots:
  void updateLightParticle(int index);
  /*! Masses (and, for a particle pair, S) from the mass table, the particles identified by Z and M. */
  void fillFromMassTable();
  void updateMassCheck();

 private:
  void buildParticle(Particle &p, const QString &title, QGroupBox *&group);
  /// The element and A shown follow Z and the mass (an edit by hand, a pair shown).
  void showNucleus(Particle &p);
  void nucleusChosen(Particle &p);
  QGroupBox *lightGroup = nullptr, *heavyGroup = nullptr;
  QWidget *excitationRow = nullptr, *lightExcitationRow = nullptr, *separationRow = nullptr, *radiusRow = nullptr;
  double lightG_ = 0., heavyG_ = 0.;
  int lastType_ = 0;     ///< the pair type shown before the last change
  PairsData held_{};     ///< the particle pair's values while another type is looked at
  bool haveHeld_ = false;
  QPushButton *okButton;
  QPushButton *cancelButton;
};

#endif
