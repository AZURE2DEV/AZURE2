#ifndef ADDSEGDATADIALOG_H
#define ADDSEGDATADIALOG_H

#include <QDialog>
#include <QShowEvent>
#include <QLineEdit>
#include <QSpinBox>
#include "PairKeyCombo.h"
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
  QCheckBox *thmCheck;

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

  /*! "Natural target": the term and each component are weighted by the isotopic abundance of their
   *  entrance pair's compound nucleus in one sample of the <abundances> block (keyword tail
   *  "natural 1" for the default sample, "natural <name>" for a named one).  Shown only when the pairs
   *  form more than one compound nucleus; reset to "No" and disabled for the observables the engine
   *  refuses (ratio, phase shift, analyzing power). */
  QLabel *naturalLabel;
  QComboBox *naturalCombo;
  /*! Fill the drop-down from the model's samples ("" = the default sample first). */
  void setSamples(const QStringList &samples);
  /*! The chosen token: "" (not natural), "1" (default sample) or a sample name. */
  QString naturalToken() const;
  void setNaturalToken(const QString &token);
  bool naturalAllowed() const;
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

 public slots:
  void setChooseFile();
  void dataTypeChanged(int);
  void varyNormChanged(int);
  void advancedModeChanged(int);
  void operationTypeChanged(int);
  void updateNaturalCheck();
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
