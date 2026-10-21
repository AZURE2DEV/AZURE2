#ifndef FITTINGTAB_H
#define FITTINGTAB_H

#include <QWidget>
#include <QMap>
#include <QPair>
#include <QSignalMapper>
#include <QPointer>
#include <QTextStream>
#include "LevelsModel.h"
#include "ChannelsModel.h"
#include "SegmentsDataModel.h"

// Forward declarations
class InfoDialog;
class LevelsTab;
class SegmentsTab;
class AbundancesModel;
class TargetIntModel;
class Config;

struct FittingParameter {
  QString name;
  double value;
  double lowerLimit;
  double upperLimit;
  double error;
  double fitError;  // NEW: Error from fitting (separate from nuisance calculation error)
  bool useAsNuisance;
  QString category;  // "level", "abundance", "norm", "shift"
  int minuitIndex;   // Index in Minuit parameters

  // For level parameters: the model rows of the level and channel (the Levels tab's models)
  int levelIndex;
  int channelIndex;
  // ... and the engine's numbering of them, energy_<engineLevel> / width_<engineLevel>_<engineChannel>
  // (LevelsTab::engineNumbering of the order the project is written in); 0 = not a parameter of the engine.
  int engineLevel = 0;
  int engineChannel = 0;
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
  /*! The isotopic abundances of the Particle Pairs tab: one parameter abundance_<Z>_<A> per compound nucleus when
   *  the project carries an <abundances> block (Value = abundance, Error = prior sigma, Nuisance = fitted). */
  void setAbundancesModel(AbundancesModel *model);
  /// The Experimental Effects tab's lines: their screening potentials are parameters too.
  void setTargetIntModel(TargetIntModel *model);

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
  /*! Mirror the Particle Pairs tab's abundance table into the abundance rows (value, sigma, fitted). */
  void syncAbundances();
  /*! Mirror the Experimental Effects tab's screening potentials into the screening rows. */
  void syncScreening();

 private:
  void setupParameterTable(QTableWidget *table, const QString &title);
  void addParameterRow(QTableWidget *table, const FittingParameter &param);
  void updateParameterFromTable(const QString &paramName, int column, const QVariant &value);
  void syncSegmentVaryStates();
  void writeAbundance(const FittingParameter &param);
  /// Append screening_<n> (eV) for every active <targetInt> line with screening, in the engine's order.
  void appendScreeningParameters();
  /// U_e, its prior sigma and its fit flag go back to the Experimental Effects tab.
  void writeScreening(const FittingParameter &param);
  void assignMinuitIndices();
  void updateParameterInOtherTabs(const QString &paramName, const FittingParameter &param);
  void writeSegmentColumn(const FittingParameter &param, int column, const QVariant &value);
  void updateParameterTableValue(const QString &paramName, double value);
  void updateParameterTableError(const QString &paramName, double error);
  void updateParameterTableCheckbox(const QString &paramName, bool checked);
  void applyParameterSettings();
  /// Rename the level entries just read into the engine numbering of the written order.
  void translateLevelSettings(bool marked);
  QString findMatchingParameterKey(const FittingParameter &param, const QStringList &savKeys);
  double convertReducedToPhysical(double reducedWidth, int levelIndex, int channelIndex);
  double convertPhysicalToReduced(double physicalWidth, int levelIndex, int channelIndex);
  double transformRWAParameterToPhysical(const QString &paramName, double rwaValue);

 public:
  // Getter for fitting parameters (for MCMCTab access)
  const QList<FittingParameter> &getFittingParameters() const { return fittingParameters; }

 private:
  QTabWidget *paramTabWidget;
  QTableWidget *levelParamsTable;
  QTableWidget *abundanceParamsTable;
  QTableWidget *screeningParamsTable;
  QTableWidget *normParamsTable;
  QTableWidget *shiftParamsTable;

  QPushButton *refreshButton;
  QPushButton *loadButton;
  QPushButton *clearLimitsButton;

  QSignalMapper *mapper;
  QPushButton *infoButton[3];
  static const std::vector<QString> infoText;
  QPointer<InfoDialog> infoDialog[3];

  QList<FittingParameter> fittingParameters;
  QList<FittingParameter> savedParameterSettings;  // Settings from <parameterSettings> section
  // loadSettings: the engine's names (energy_N, width_N_M) of the project file on disk, by model rows
  QMap<QPair<int, int>, QString> diskEngineNames_;

  // Tab references for reading current GUI state
  LevelsTab *levelsTab_;
  SegmentsTab *segmentsTab_;
  AbundancesModel *abundancesModel_ = nullptr;
  TargetIntModel *targetIntModel_ = nullptr;
  Config *config_ = nullptr;  // for the current .azr path (backup on load)
};

#endif