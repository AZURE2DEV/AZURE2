#ifndef ADDSEGDATADIALOG_H
#define ADDSEGDATADIALOG_H

#include <QDialog>
#include <QShowEvent>
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
 * Dialog for adding or editing one data segment, including the sum/ratio components of a composite one.
 */
class AddSegDataDialog : public QDialog {
  Q_OBJECT

 public:
  AddSegDataDialog(QWidget *parent = 0);
  /*! Call before exec(): the pair drop-downs list this model's pairs by name and hand back the keys. */
  void setPairsModel(PairsModel *model);
  PairKeyCombo *entrancePairCombo;
  PairKeyCombo *exitPairCombo;
  QLineEdit *lowEnergyText;
  QLineEdit *highEnergyText;
  QLineEdit *lowAngleText;
  QLineEdit *highAngleText;
  QComboBox *dataTypeCombo;

  /*! The observable code written to the .azr. It is not the combo box index:
   *  the analyzing power is code 7 in both segment kinds, but the two lists
   *  have different lengths, so the code is carried as item data instead. */
  int dataTypeCode() const;
  void setDataTypeCode(int code);

  QLineEdit *dataFileText;
  QLineEdit *dataNormText;
  QLineEdit *dataNormErrorText;
  QLabel *dataNormErrorLabel;
  QCheckBox *varyNormCheck;
  QLineEdit *phaseJValueText;
  QLineEdit *phaseLValueText;
  QLabel *phaseLValueLabel;
  QLabel *phaseJValueLabel;
  QLabel *totalCaptureLabel;
  QLineEdit *energyShiftText;
  QLabel *energyShiftLabel;
  QLineEdit *energyShiftErrorText;
  QLabel *energyShiftErrorLabel;
  QCheckBox *varyEnergyShiftCheck;

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

  // UPOS (Unobserved Primary, Observed Secondary) widgets
  QCheckBox *uposCheck;
  QLabel *secondaryLLabel;
  QSpinBox *secondaryLSpin;
  QLabel *finalJLabel;
  QLineEdit *finalJText;
  QLabel *deltaLabel;
  QLineEdit *deltaText;

  /*! The row's trailing keyword text ("sqrtshift ...", "mixture ..."): the
   *  mixture block is edited here, anything else rides along untouched. */
  MixtureEditor *mixtureEditor;
  /*! A mixture is only defined for angle-integrated and differential cross
   *  sections (plain or summed); the reason it is not, otherwise. */
  bool mixtureCompatible(QString *why = 0) const;

 public slots:
  void setChooseFile();
  void tryAccept();
  void dataTypeChanged(int);
  void varyNormChanged(int);
  void advancedModeChanged(int);
  void operationTypeChanged(int);
  void useFixedAngleChanged(int);
  void addComponent();
  void removeComponent();
  void uposChanged(int);

 protected:
  /*! Sized to the content within the screen on every opening (DialogScrollArea.h). */
  void showEvent(QShowEvent *event) override;

 private:
  QPushButton *okButton;
  QPushButton *cancelButton;
};

#endif
