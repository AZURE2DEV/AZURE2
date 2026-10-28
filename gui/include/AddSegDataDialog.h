#ifndef ADDSEGDATADIALOG_H
#define ADDSEGDATADIALOG_H

#include <QDialog>
#include <QShowEvent>
#include <QLineEdit>
#include <QSpinBox>
#include "PairKeyCombo.h"
#include "ScatteringLengthsModel.h"
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
  /*! Titled "C.M. Angle [degrees]" for a C.M. differential cross section and an analyzing power, "Lab Angle [degrees]" otherwise. */
  QGroupBox *angleBox;
  /*! The observable: Angle Integrated (0), Differential (1), Phase Shift (2), Analyzing Power (7). Each item
   *  carries that base isDiff code as item data; look items up with findData, never by position. */
  QComboBox *dataTypeCombo;
  /*! Differential only: Lab (isDiff 1) or C.M. (isDiff 4), as item data. */
  QLabel *frameLabel;
  QComboBox *frameCombo;
  /*! Angle Integrated with a gamma exit pair only: All (0), E1 only (5), E2 only (6), as item data.
   *  Hidden otherwise, and then read as All. */
  QLabel *multipolarityLabel;
  QComboBox *multipolarityCombo;

  /*! The isDiff code written to the .azr, from the observable, its frame or multipolarity, and the exit pair
   *  (the exit item "Total capture", key PairKeyCombo::TotalCapture, with Angle Integrated is code 3). */
  int dataTypeCode() const;
  /*! Sets the observable, frame, multipolarity and, for code 3, the total-capture exit item. */
  void setDataTypeCode(int code);
  /*! The base code of the observable list's current item. */
  int observableCode() const;
  /*! True when the multipolarity drop-down applies: Angle Integrated with a gamma exit pair. */
  bool multipolarityApplies() const;

  /*! Data Type "Scattering Length": not a segment line but a line of the <scatteringLengths> block (pair, channel
   *  spin, Re a0 +- sigma, optional Im a0 +- sigma); SegmentsTab lists those under the data segments. */
  static const int ScatteringLengthCode = 100;
  bool isScatteringLength() const { return observableCode() == ScatteringLengthCode; }
  ScatteringLengthEntry scatteringLength() const;
  /*! Edit mode for a scattering length: the type is fixed and the fields filled from \p entry. */
  void setScatteringLength(const ScatteringLengthEntry &entry);
  /*! Edit mode for a segment line: drop the "Scattering Length" type (a segment cannot become one). */
  void disallowScatteringLength();
  QGroupBox *scatLenBox;
  QComboBox *scatLenSpinCombo;
  QLineEdit *scatLenValueText;
  QLineEdit *scatLenSigmaText;
  QLineEdit *scatLenImValueText;
  QLineEdit *scatLenImSigmaText;

  QLineEdit *dataFileText;
  /*! Optional data covariance file (keyword block "covariance <file>"): every data segment naming the same
   *  file shares one covariance matrix.  Empty = uncorrelated errors (no block). */
  QLineEdit *covarianceFileText;
  QLineEdit *dataNormText;
  QLineEdit *dataNormErrorText;
  QLabel *dataNormErrorLabel;
  QCheckBox *varyNormCheck;
  QLineEdit *phaseJValueText;
  QLineEdit *phaseLValueText;
  QLineEdit *phaseSValueText;   ///< channel spin of a phase-shift segment ("spin <s>"), empty = not given
  QLabel *phaseSValueLabel;
  QLabel *phaseLValueLabel;
  QLabel *phaseJValueLabel;
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
  void setChooseCovarianceFile();
  void dataTypeChanged(int);
  void exitPairChanged(int key);
  void varyNormChanged(int);
  void advancedModeChanged(int);
  void operationTypeChanged(int);
  void updateNaturalCheck();
  void useFixedAngleChanged(int);
  void addComponent();
  void removeComponent();
  void uposChanged(int);
  void fillScatteringLengthSpins();

 protected:
  /*! Sized to the content within the screen on every opening (DialogScrollArea.h). */
  void showEvent(QShowEvent *event) override;

 private:
  /*! Show the scattering-length fields and hide the segment ones, or the reverse (restoring the sections the
   *  check boxes and the observable decide). */
  void applyScatteringLengthMode(bool on);
  QList<QWidget *> segmentOnly;
  QLabel *entrancePairLabel;
  QLabel *exitPairLabel;
  bool scatLenMode = false;
  bool scatLenActive = true;
  PairsModel *pairs = nullptr;
  QPushButton *okButton;
  QPushButton *cancelButton;
};

#endif
