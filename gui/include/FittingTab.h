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
  void resetToDefaults();
  void loadSettings();
  void refreshParameters();

 private:
  void setupParameterTable(QTableWidget* table, const QString& title);
  void addParameterRow(QTableWidget* table, const FittingParameter& param);
  void updateParameterFromTable(const QString& paramName, int column, const QVariant& value);
  void updateParameterInOtherTabs(const QString& paramName, const FittingParameter& param);
  void applyParameterSettings();
  double convertReducedToPhysical(double reducedWidth, int levelIndex, int channelIndex);
  double convertPhysicalToReduced(double physicalWidth, int levelIndex, int channelIndex);

  QTabWidget* paramTabWidget;
  QTableWidget* levelParamsTable;
  QTableWidget* normParamsTable;
  QTableWidget* shiftParamsTable;
  
  QPushButton* resetButton;
  QPushButton* refreshButton;
  QPushButton* loadButton;
  
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