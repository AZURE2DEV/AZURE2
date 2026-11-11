#ifndef EDITOPTIONSDIALOG_H
#define EDITOPTIONSDIALOG_H


#include <QDialog>

QT_BEGIN_NAMESPACE

class QPushButton;
class QGroupBox;
class QCheckBox;
class QComboBox;
class QLineEdit;

QT_END_NAMESPACE

QT_END_NAMESPACE

/*!
 * Dialog for the Runtime Options -- the formalism switches that become the paramMask bits, and the kinematics.
 */
class EditOptionsDialog : public QDialog {
  Q_OBJECT

 public:
  EditOptionsDialog(QWidget *parent = 0);
  /*! The R-matrix parametrization, one choice (2026-11-09; before, two check boxes "Use Brune formalism" and "Use Park
   *  parametrization", the second implying the first).  Item data and indices are the values below.  "Brune-Park" is
   *  the engine's Park mode (USE_PARK_FORMALISM, --use-park): Park's observed reduced widths on Brune's level matrix,
   *  named so in the GUI to keep Brune's part in it visible (DeBoer). */
  enum Parametrization { Formal = 0, Brune = 1, BrunePark = 2 };
  QComboBox *parametrizationCombo;
  int parametrization() const;
  void setParametrization(int choice);
  /*! The choice a paramMask holds: Park bit -> BrunePark, else Brune bit -> Brune, else Formal. */
  static int parametrizationFromMask(unsigned int mask);
  /*! mask with the Brune and Park bits set for choice (BrunePark sets both, as the engine expects); other bits kept. */
  static unsigned int maskWithParametrization(unsigned int mask, int choice);
  QCheckBox *useGSLCoulCheck;
  QCheckBox *ignoreExternalsCheck;
  QCheckBox *useRMCCheck;
  QCheckBox *noTransformCheck;
  QCheckBox *useHybridMethodCheck;
  QCheckBox *useAdaptiveGridCheck;
  /*! Relativistic two-body kinematics: "relativistic=1" in <config>.  Unlike the switches above it is part of the
   *  project, so the command line and pyazr use it too (before 2026-11-01 it sat below the Particle Pairs table). */
  QCheckBox *relativisticCheck;
  /*! Prior parameter covariance file: "priorCovariance=<path>" in <config> (empty = none), a path relative to the
   *  project's directory.  Part of the project like the kinematics: command-line and pyazr runs use it too. */
  QLineEdit *priorCovarianceEdit;
  /*! The project's directory: a file chosen inside it is entered relative to it.  Empty = absolute paths. */
  void setProjectDirectory(const QString &dir) { projectDirectory_ = dir; }
  // QCheckBox* noLongWavelengthCheck;

 private slots:
  void choosePriorCovariance();
  void parametrizationChanged(int);
  void useRMCCheckChanged(int);

 private:
  QPushButton *okButton;
  QPushButton *cancelButton;
  QPushButton *priorCovarianceButton;
  QString projectDirectory_;
};

#endif
