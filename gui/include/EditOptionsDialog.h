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
 * Dialog for the Runtime Options -- the formalism switches that become the paramMask bits.
 */
class EditOptionsDialog : public QDialog {
  Q_OBJECT

 public:
  EditOptionsDialog(QWidget *parent = 0);
  QCheckBox *useBruneCheck;
  QCheckBox *useGSLCoulCheck;
  QCheckBox *ignoreExternalsCheck;
  QCheckBox *useRMCCheck;
  QCheckBox *noTransformCheck;
  QCheckBox *useHybridMethodCheck;
  QCheckBox *useAdaptiveGridCheck;
  // QCheckBox* noLongWavelengthCheck;

 private slots:
  void useBruneCheckChanged(int);
  void useRMCCheckChanged(int);

 private:
  QPushButton *okButton;
  QPushButton *cancelButton;
};

#endif
