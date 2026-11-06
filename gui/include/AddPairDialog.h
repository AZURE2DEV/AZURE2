#ifndef ADDPAIRDIALOG_H
#define ADDPAIRDIALOG_H

#include <QDialog>

QT_BEGIN_NAMESPACE

class QLineEdit;
class QLabel;
class QCheckBox;
class QGroupBox;
class QComboBox;
class QPushButton;

QT_END_NAMESPACE

/*!
 * Dialog for adding or editing one particle pair.
 */
class AddPairDialog : public QDialog {
  Q_OBJECT

 public:
  AddPairDialog(QWidget *parent = 0);
  QLineEdit *lightJText;
  QComboBox *lightPiCombo;
  QLineEdit *lightZText;
  QLineEdit *lightMText;
  // QLineEdit *lightGText;
  QLineEdit *heavyJText;
  QComboBox *heavyPiCombo;
  QLineEdit *heavyZText;
  QLineEdit *heavyMText;
  // QLineEdit *heavyGText;
  QLineEdit *excitationEnergyText;
  QLineEdit *seperationEnergyText;
  QLineEdit *channelRadiusText;
  QComboBox *pairTypeCombo;
  QCheckBox *e1Check;
  // QCheckBox *m1Check;
  QCheckBox *e2Check;
  QGroupBox *multBox;
  /*! Optional isotope names ("p", "a", "12C"): an entry the mass table knows fills that particle's Z and M. */
  QLineEdit *lightIsotopeText;
  QLineEdit *heavyIsotopeText;
  /*! Fills both masses and the separation energy from the mass table (MassTable, AME2020). */
  QPushButton *fillButton;
  /*! The entered masses and S against the mass table, updated as they are typed. */
  QLabel *massCheckLabel;

 public slots:
  void updateLightParticle(int index);
  /*! A particle's isotope field was edited: look it up, fill Z and M. */
  void lightIsotopeEdited();
  void heavyIsotopeEdited();
  /*! Masses (and, for a particle pair, S) from the mass table, the particles identified by Z and M. */
  void fillFromMassTable();
  void updateMassCheck();

 private:
  QLabel *lightJLabel;
  QLabel *lightPiLabel;
  QLabel *lightZLabel;
  QLabel *lightMLabel;
  // QLabel *lightGLabel;
  QLabel *heavyJLabel;
  QLabel *heavyPiLabel;
  QLabel *heavyZLabel;
  QLabel *heavyMLabel;
  // QLabel *heavyGLabel;
  QLabel *excitationEnergyLabel;
  QLabel *seperationEnergyLabel;
  QLabel *channelRadiusLabel;
  QLabel *pairTypeLabel;
  QLabel *lightIsotopeLabel;
  QLabel *heavyIsotopeLabel;
  bool isotopeEdited(QLineEdit *isotope, QLineEdit *zText, QLineEdit *mText);

  QPushButton *okButton;
  QPushButton *cancelButton;
};

#endif
