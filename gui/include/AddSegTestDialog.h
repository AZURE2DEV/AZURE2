#ifndef ADDSEGTESTDIALOG_H
#define ADDSEGTESTDIALOG_H

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
  /*! Titled "C.M. Angle [degrees]" for a C.M. differential cross section and an analyzing power, "Lab Angle [degrees]" otherwise. */
  QGroupBox *angleBox;
  QLineEdit *angleStepText;
  /*! The observable: Angle Integrated (0), Differential (1), Phase Shift (2), Angular Distribution
   *  Coefficients (3), Analyzing Power (7). Each item carries that base isDiff code as item data; look items
   *  up with findData, never by position. */
  QComboBox *dataTypeCombo;
  QCheckBox *thmCheck;
  /*! Differential only: Lab (isDiff 1) or C.M. (isDiff 5), as item data. */
  QLabel *frameLabel;
  QComboBox *frameCombo;

  /*! The isDiff code written to the .azr, from the observable, its frame and the exit pair (the exit item
   *  "Total capture", key PairKeyCombo::TotalCapture, with Angle Integrated is code 4). */
  int dataTypeCode() const;
  /*! Sets the observable, frame and, for code 4, the total-capture exit item. */
  void setDataTypeCode(int code);
  /*! The base code of the observable list's current item. */
  int observableCode() const;

  QLineEdit *phaseJValueText;
  QLineEdit *phaseLValueText;
  QLabel *phaseJValueLabel;
  QLabel *phaseLValueLabel;
  QLabel *angDistLabel;
  QSpinBox *angDistSpin;

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

 public slots:
  void dataTypeChanged(int);
  void exitPairChanged(int key);
  void advancedModeChanged(int);
  void operationTypeChanged(int);
  void updateNaturalCheck();
  void useFixedAngleChanged(int);
  void addComponent();
  void removeComponent();

 protected:
  /*! Sized to the content within the screen on every opening (DialogScrollArea.h). */
  void showEvent(QShowEvent *event) override;

 private:
  QPushButton *okButton;
  QPushButton *cancelButton;
};

#endif
