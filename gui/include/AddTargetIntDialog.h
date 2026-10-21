#ifndef ADDTARGETINTDIALOG_H
#define ADDTARGETINTDIALOG_H

#include <QDialog>
#include <QRadioButton>
#include <QTableWidget>
#include <functional>
#include <string>
#include <vector>

QT_BEGIN_NAMESPACE

class QLineEdit;
class QSpinBox;
class QDoubleSpinBox;
class QCheckBox;
class QTableWidget;
class QGroupBox;
class QSize;
class QComboBox;
class QPushButton;

QT_END_NAMESPACE

/*!
 * Dialog for one target effect: convolution, target integration, and the attenuation and convolution coefficients.
 */
class AddTargetIntDialog : public QDialog {
  Q_OBJECT

 public:
  AddTargetIntDialog(QWidget *parent = 0);
  QCheckBox *isConvolutionCheck;
  // One Gaussian beam spread. sigmaText holds its width in lab MeV: a number ("0.03") or an expression of
  // the lab energy x ("a0*x^1.5"); the tab stores a number in the constant slot of the line and an expression
  // with the sigmalab keyword. A line from before 2026-10 carries its expression in the c.m. frame; for such
  // a line the frame row appears, with the exact conversion to the lab form.
  QRadioButton *sigmaLabRadio;
  QRadioButton *sigmaCmRadio;
  QPushButton *convertToLabButton;
  QWidget *frameRow;
  /*! True when sigmaText parses as a plain number, i.e. a constant width. */
  bool sigmaIsConstant() const;
  /*! c.m./lab energy ratio of the entrance pair the effect's segments share (0 = unknown or mixed); it
   *  enables the exact rewrite of a legacy c.m. expression into the lab form. */
  void setConversionFactor(double factor);
  QCheckBox *isTargetIntegrationCheck;
  QCheckBox *isQCoefficientCheck;
  QLineEdit *sigmaText;
  QLineEdit *segmentsListText;
  QSpinBox *numPointsSpin;
  QSpinBox *numParametersSpin;
  QSpinBox *numQCoefficientSpin;
  QSpinBox *numConvCoefficientSpin;
  QLineEdit *densityText;
  QLineEdit *stoppingPowerEqText;
  QComboBox *elementComboBox;
  QLineEdit *compoundText;
  QComboBox *activeElementComboBox;
  QPushButton *fetchStoppingPowerButton;
  QCheckBox *isStraggling;
  QLineEdit *stragglingCoefficientText;
  QCheckBox *isBeamProfileCheck;
  QCheckBox *beamPhotodissociationCheck;
  QSpinBox *numBeamComponentSpin;
  QTableWidget *beamProfileTable;
  QLineEdit *beamTpcSigmaText;
  QLineEdit *beamTruncationText;
  QList<double> tempBeamProfile;   // flattened (xi, omega, alpha, weight) quadruples
  QCheckBox *isUdrCheck;
  QLineEdit *udrFileText;
  QPushButton *udrBrowseButton;
  QLineEdit *udrFlightPathText;
  QLineEdit *udrBurstText;
  QLineEdit *udrChannelText;
  QCheckBox *udrCentredCheck;
  // Electron screening of the listed segments: U_e (eV), fitted or not, with an optional prior width.
  QCheckBox *isScreeningCheck;
  QLineEdit *screeningUeText;
  QPushButton *adiabaticButton;
  QCheckBox *screeningVaryCheck;
  QLineEdit *screeningSigmaText;
  /*! The adiabatic-limit U_e (eV) of a segment list ("6,7"), < 0 when it has none (neutral or mixed
   *  entrance pairs); the tab installs it, the dialog asks it when screening is switched on. */
  void setAdiabaticProvider(std::function<double(const QString &)> provider);
  QLineEdit *applyRangesText;
  QLineEdit *transitionWidthText;
  QLineEdit *autoToleranceText;
  QDoubleSpinBox *resonanceWidthMultiplierSpin;
  QDoubleSpinBox *pointsPerWidthSpin;
  QLineEdit *energyText;
  QLineEdit *deltaEText;
  QPushButton *calculateDeltaEButton;
  QTableWidget *parametersTable;
  QTableWidget *qCoefficientTable;
  QTableWidget *convCoefficientTable;
  QList<double> tempParameters;
  QList<double> tempQCoefficients;
  QList<double> tempConvCoefficients;
  QLineEdit *convolutionEqText;
  /*! Sizes a coefficient table to show every row (up to 12 before it scrolls) and lets the dialog grow with it. */
  void fitTableHeight(QTableWidget *table);
  void createParameterItem(int row, double value = 0.0);
  void createQCoefficientItem(int row, double value = 1.0);
  void createConvCoefficientItem(int row, double value = 1.0);
  void createBeamProfileItem(int row, double xi = 0.0, double omega = 0.0, double alpha = 0.0, double weight = 1.0);

 public slots:
  void convolutionCheckChanged(bool checked);

  void targetIntCheckChanged(bool checked);
  void parameterSpinChanged(int newNumber);
  void parameterChanged(int row, int column);

  void qCoefficientCheckChanged(bool checked);
  void qCoefficientSpinChanged(int newNumber);
  void qCoefficientChanged(int row, int column);

  void sigmaModeChanged();
  void convertExpressionToLab();
  void convCoefficientSpinChanged(int newNumber);
  void convCoefficientChanged(int row, int column);

  void beamProfileCheckChanged(bool checked);
  void beamComponentSpinChanged(int newNumber);
  void beamProfileChanged(int row, int column);

  void udrCheckChanged(bool checked);
  void udrBrowse();

  void screeningCheckChanged(bool checked);
  /// Fill U_e with the adiabatic limit of the listed segments' entrance pair.
  void useAdiabaticScreening();

  void elementSelectionChanged(int index);
  void compoundFormulaChanged(const QString &formula);
  void fetchStoppingPowerParameters();
  void calculateDeltaE();

 private:
  QPushButton *okButton;
  QPushButton *cancelButton;
  QGroupBox *stoppingPowerBox;
  QGroupBox *qCoefficientBox;
  QGroupBox *convCoefficientBox;
  QGroupBox *beamProfileBox;
  QGroupBox *udrBox;
  QGroupBox *screeningBox;
  std::function<double(const QString &)> adiabaticProvider_;

  int selectedElement_;
  double conversionFactor_ = 0.0;
  int activeElementNumber() const;
  void populateElementComboBox();
  void updateStoppingPowerFromElement(int element);
  void updateStoppingPowerFromCompound(const std::string &formula);
  void updateStoppingPowerGUI(const QString &equation, const std::vector<double> &parameters, const QString &materialName);
};

#endif
