#ifndef FITTINGTAB_H
#define FITTINGTAB_H

#include <QWidget>
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

struct FittingParameter {
  QString name;
  double value;
  double lowerLimit;
  double upperLimit;
  double error;
  double fitError; // NEW: Error from fitting (separate from nuisance calculation error)
  bool useAsNuisance;
  QString category; // "level", "norm", "shift"
  int minuitIndex; // Index in Minuit parameters
  
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

QT_END_NAMESPACE

class FittingTab : public QWidget {
  Q_OBJECT

 public:
  FittingTab(QWidget* parent=0);
  friend class AZURESetup;
  void reset();
  void updateParameterTables();
  bool writeParameterSettings(QTextStream& outStream);
  bool readParameterSettings(QTextStream& inStream);
  void refreshFromMinuitParameters();
  void populateFromCurrentGUIState();
  void setTabReferences(LevelsTab* levelsTab, SegmentsTab* segmentsTab);

 public slots:
  void showInfo(int which=0, QString title="");

 private slots:
  void parameterItemChanged(QTableWidgetItem* item);
  void loadSettings();
  void refreshParameters();
  void populateWignerLimits();
  void clearLimits();
  void onSegmentNormalizationChanged(int segmentIndex, double value);
  void onSegmentEnergyShiftChanged(int segmentIndex, double value);
  void onSegmentNormalizationErrorChanged(int segmentIndex, double error);
  void onSegmentEnergyShiftErrorChanged(int segmentIndex, double error);
  void onSegmentNormalizationVaryChanged(int segmentIndex, bool vary);
  void onSegmentEnergyShiftVaryChanged(int segmentIndex, bool vary);

 private:
  void setupParameterTable(QTableWidget* table, const QString& title);
  void addParameterRow(QTableWidget* table, const FittingParameter& param);
  void updateParameterFromTable(const QString& paramName, int column, const QVariant& value);
  void updateParameterInOtherTabs(const QString& paramName, const FittingParameter& param);
  void updateParameterTableValue(const QString& paramName, double value);
  void updateParameterTableError(const QString& paramName, double error);
  void updateParameterTableCheckbox(const QString& paramName, bool checked);
  void applyParameterSettings();
  QString findMatchingParameterKey(const FittingParameter& param, const QStringList& savKeys);
  double convertReducedToPhysical(double reducedWidth, int levelIndex, int channelIndex);
  double convertPhysicalToReduced(double physicalWidth, int levelIndex, int channelIndex);
  double transformRWAParameterToPhysical(const QString& paramName, double rwaValue);

 public:
  // Getter for fitting parameters (for MCMCTab access)
  const QList<FittingParameter>& getFittingParameters() const { return fittingParameters; }

 private:
  QTabWidget* paramTabWidget;
  QTableWidget* levelParamsTable;
  QTableWidget* normParamsTable;
  QTableWidget* shiftParamsTable;
  
  QPushButton* refreshButton;
  QPushButton* loadButton;
  QPushButton* populateWignerLimitsButton;
  QPushButton* clearLimitsButton;
  
  QSignalMapper* mapper;
  QPushButton *infoButton[3];
  static const std::vector<QString> infoText;
  QPointer<InfoDialog> infoDialog[3];
  
  QList<FittingParameter> fittingParameters;
  QList<FittingParameter> savedParameterSettings; // Settings from <parameterSettings> section
  
  // Tab references for reading current GUI state
  LevelsTab* levelsTab_;
  SegmentsTab* segmentsTab_;
};

#endif