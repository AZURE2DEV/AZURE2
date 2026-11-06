#ifndef EDITOPTIONSDIALOG_H
#define EDITOPTIONSDIALOG_H


#include <QDialog>

QT_BEGIN_NAMESPACE

class QPushButton;
class QGroupBox;
class QCheckBox;

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
  // QCheckBox* noLongWavelengthCheck;

 private slots:
  void useBruneCheckChanged(int);
  void useParkCheckChanged(int);
  void useRMCCheckChanged(int);

 private:
  QPushButton *okButton;
  QPushButton *cancelButton;
};

#endif
