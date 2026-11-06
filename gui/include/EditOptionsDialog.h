#ifndef EDITOPTIONSDIALOG_H
#define EDITOPTIONSDIALOG_H


#include <QDialog>

QT_BEGIN_NAMESPACE

class QPushButton;
class QGroupBox;
class QCheckBox;
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
  QCheckBox *useBruneCheck;
  QCheckBox *useParkCheck;
  QCheckBox *useGSLCoulCheck;
  QCheckBox *ignoreExternalsCheck;
  QCheckBox *useRMCCheck;
  QCheckBox *noTransformCheck;
  QCheckBox *useHybridMethodCheck;
  QCheckBox *useThmCheck;
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
  void useBruneCheckChanged(int);
  void useParkCheckChanged(int);
  void useRMCCheckChanged(int);

 private:
  QPushButton *okButton;
  QPushButton *cancelButton;
  QPushButton *priorCovarianceButton;
  QString projectDirectory_;
};

#endif
