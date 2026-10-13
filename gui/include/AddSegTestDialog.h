#ifndef ADDSEGTESTDIALOG_H
#define ADDSEGTESTDIALOG_H

#include <QDialog>
#include <QLineEdit>
#include <QSpinBox>
#include "PairKeyCombo.h"
#include "MixtureEditor.h"
#include <QComboBox>
#include <QLabel>
#include <QCheckBox>
#include <QListWidget>
#include <QGroupBox>
#include <QPushButton>

QT_BEGIN_NAMESPACE

class QLabel;

QT_END_NAMESPACE

/*!
 * Dialog for adding or editing one test segment -- an energy and angle grid with no data file.
 */
class AddSegTestDialog : public QDialog {
  Q_OBJECT

 public:
  AddSegTestDialog(QWidget *parent = 0);
  /*! Call before exec(): the pair drop-downs list this model's pairs by name and hand back the keys. */
  void setPairsModel(PairsModel *model);
  PairKeyCombo *entrancePairCombo;
  PairKeyCombo *exitPairCombo;
  QLineEdit *lowEnergyText;
  QLineEdit *highEnergyText;
  QLineEdit *energyStepText;
  QLineEdit *lowAngleText;
  QLineEdit *highAngleText;
  QLineEdit *angleStepText;
  QComboBox *dataTypeCombo;

  /*! The observable code written to the .azr. It is not the combo box index:
   *  the analyzing power is code 7 in both segment kinds, but the two lists
   *  have different lengths, so the code is carried as item data instead. */
  int dataTypeCode() const;
  void setDataTypeCode(int code);

  QLineEdit *phaseJValueText;
  QLineEdit *phaseLValueText;
  QLabel *phaseJValueLabel;
  QLabel *phaseLValueLabel;
  QLabel *angDistLabel;
  QLabel *totalCaptureLabel;
  QSpinBox *angDistSpin;

  QCheckBox *advancedModeCheck;
  QGroupBox *advancedModeBox;
  QComboBox *operationCombo;
  QListWidget *componentsList;
  PairKeyCombo *componentEntranceCombo;
  PairKeyCombo *componentExitCombo;
  QCheckBox *useFixedAngleCheck;
  QLineEdit *fixedAngleText;
  QLabel *fixedAngleLabel;
  QLabel *componentScalingLabel;
  QLineEdit *componentScalingText;
  QPushButton *addComponentButton;
  QPushButton *removeComponentButton;

  /*! The row's trailing keyword text; the mixture block is edited here. */
  MixtureEditor *mixtureEditor;
  bool mixtureCompatible(QString *why = 0) const;

 public slots:
  void dataTypeChanged(int);
  void tryAccept();
  void advancedModeChanged(int);
  void operationTypeChanged(int);
  void useFixedAngleChanged(int);
  void addComponent();
  void removeComponent();

 private:
  QPushButton *okButton;
  QPushButton *cancelButton;
};

#endif
