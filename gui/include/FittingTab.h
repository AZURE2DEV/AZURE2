#ifndef FITTINGTAB_H
#define FITTINGTAB_H

#include <QWidget>
#include "ThmSettings.h"
#include <QSignalMapper>
#include <QPointer>
#include <QTextStream>
#include "LevelsModel.h"
#include "ChannelsModel.h"
#include "SegmentsDataModel.h"

class AZURESetup;

// Forward declarations
class InfoDialog;
class LevelsTab;
class SegmentsTab;
class Config;

struct FittingParameter {
  QString name;
  double value;
  double lowerLimit;
  double upperLimit;
  double error;
  double fitError;  // NEW: Error from fitting (separate from nuisance calculation error)
  bool useAsNuisance;
  QString category;  // "level", "norm", "shift", "cbkg" (THM coherent background)
  int minuitIndex;   // Index in Minuit parameters

  // For level parameters
  int levelIndex;
  int channelIndex;
};

QT_BEGIN_NAMESPACE

class QTabWidget;
class QTableWidget;
class QTableWidgetItem;
class QPushButton;
class QVBoxLayout;
class QHBoxLayout;
class QGroupBox;
class QLabel;
class QCheckBox;

QT_END_NAMESPACE

/*!
 * The Fitting tab: the parameter list Minuit will minimize, and its results.
 *
 * The tables mirror the level scheme, so they are rebuilt whenever the levels or pairs change; refreshFromMinuitParameters brings the fitted values back afterwards.
 */
class FittingTab : public QWidget {
  Q_OBJECT

 public:
  FittingTab(QWidget *parent = 0);
  friend class AZURESetup;
  void reset();
  void updateParameterTables();
  bool writeParameterSettings(QTextStream &outStream);
  bool readParameterSettings(QTextStream &inStream);
  void refreshFromMinuitParameters();
  void populateFromCurrentGUIState();
  void setTabReferences(LevelsTab *levelsTab, SegmentsTab *segmentsTab);
  /// Writes `values` (cbkg_* name -> value) into the cbackground= values of
  /// the THM block (ApplyThmCoherentValues); false if nothing was changed.
  bool applyCoherentValues(const QMap<QString, double> &values);
  /*! Follows a renumbering of the data segments (SegmentsTab::
      dataSegmentsRenumbered): the segment_N_norm / segment_N_energy_shift
      parameters, their saved settings and prior centres take the new
      numbers, those of a deleted segment are dropped, and so are the THM
      background parameters of an experiment that is gone; a line added at
      the end gets its parameters. */
  void followSegments(const QVector<int> &newNumber);

 protected:
  void showEvent(QShowEvent *event) override;
  void setConfig(Config *config) { config_ = config; }

 public slots:
  void showInfo(int which = 0, QString title = "");

 private slots:
  void parameterItemChanged(QTableWidgetItem *item);
  void loadSettings();
  void refreshParameters();
  void clearLimits();
  void onSegmentNormalizationChanged(int segmentIndex, double value);
  void onSegmentEnergyShiftChanged(int segmentIndex, double value);
  void onSegmentNormalizationErrorChanged(int segmentIndex, double error);
  void onSegmentEnergyShiftErrorChanged(int segmentIndex, double error);
  void onSegmentNormalizationVaryChanged(int segmentIndex, bool vary);
  void onSegmentEnergyShiftVaryChanged(int segmentIndex, bool vary);

 private:
  void setupParameterTable(QTableWidget *table, const QString &title, bool priorCentreColumn = false);
  void addParameterRow(QTableWidget *table, const FittingParameter &param);
  void updateParameterFromTable(const QString &paramName, int column, const QVariant &value);
  void syncSegmentVaryStates();
  void assignMinuitIndices();
  void updateParameterInOtherTabs(const QString &paramName, const FittingParameter &param);
  void writeSegmentColumn(const FittingParameter &param, int column, const QVariant &value);
  void updateParameterTableValue(const QString &paramName, double value);
  void updateParameterTableError(const QString &paramName, double error);
  void updateParameterTableCheckbox(const QString &paramName, bool checked);
  void applyParameterSettings();
  QList<int> engineLevelOrder(const QList<int> &fileOrder);
  QString findMatchingParameterKey(const FittingParameter &param, const QStringList &savKeys);
  double convertReducedToPhysical(double reducedWidth, int levelIndex, int channelIndex);
  double convertPhysicalToReduced(double physicalWidth, int levelIndex, int channelIndex);
  double transformRWAParameterToPhysical(const QString &paramName, double rwaValue);
  /// The THM coherent backgrounds (cbackground= of the <thm> experiments):
  /// their free parameters, last, as EData::FillMnParams adds them.
  void appendCoherentParameters();
  /// The segment_<i+1>_norm (norm) or _energy_shift parameter of segment line i.
  static FittingParameter segmentParameter(int i, const SegmentsDataData &segment, bool norm);
  /*! Every cbkg_* parameter name of the experiment `record`, in AZURE2's
      order (EData::BuildThmGroups, mirrored on the Levels and Segments tabs),
      with each one's start value and fixed flag. */
  QStringList coherentNames(const ThmExperimentRecord &record, QList<double> *values = nullptr,
                            QList<bool> *fixed = nullptr) const;
  AZURESetup *setup() const;

 public:
  // Getter for fitting parameters (for MCMCTab access)
  const QList<FittingParameter> &getFittingParameters() const { return fittingParameters; }

 private:
  QTabWidget *paramTabWidget;
  QTableWidget *levelParamsTable;
  QTableWidget *normParamsTable;
  QTableWidget *shiftParamsTable;
  QTableWidget *cbkgParamsTable;  ///< THM coherent background; its tab only when there are some

  QPushButton *refreshButton;
  QPushButton *loadButton;
  QPushButton *clearLimitsButton;

  QSignalMapper *mapper;
  QPushButton *infoButton[3];
  static const std::vector<QString> infoText;
  QPointer<InfoDialog> infoDialog[3];

  QList<FittingParameter> fittingParameters;
  QList<FittingParameter> savedParameterSettings;  // Settings from <parameterSettings> section
  /*! Explicit prior centres ("segment_N_norm prior_centre c" and
      "segment_N_energy_shift prior_centre c" rows of <parameterSettings>,
      EData::ReadPriorCentres), by parameter name.  A norm or shift without one
      has its prior centred on the Segments-tab value, as in every classic
      file.  Kept apart from fittingParameters, which is rebuilt from the tabs
      (and by a .sav load): a centre survives both. */
  QMap<QString, double> priorCentres_;
  /// The rows above, in segment order (norm before shift), as written.
  QStringList priorCentreRows() const;

  // Tab references for reading current GUI state
  LevelsTab *levelsTab_;
  SegmentsTab *segmentsTab_;
  Config *config_ = nullptr;  // for the current .azr path (backup on load)
};

#endif