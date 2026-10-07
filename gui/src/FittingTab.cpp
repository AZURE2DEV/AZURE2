#include <QTabWidget>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QPushButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QLabel>
#include <QHeaderView>
#include <QMessageBox>
#include <QSignalMapper>
#include <QCheckBox>
#include <QFileDialog>
#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QRegExp>
#include <QShowEvent>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "FittingTab.h"
#include <QSignalBlocker>
#include <QSet>
#include "RoundTripNumber.h"
#include "InfoDialog.h"
#include "LevelsTab.h"
#include "SegmentsTab.h"
#include "AZURESetup.h"
#include "CNuc.h"
#include "Config.h"
#include "AZUREParams.h"

FittingTab::FittingTab(QWidget *parent) :
  QWidget(parent),
  levelsTab_(nullptr),
  segmentsTab_(nullptr) {
  // Create main layout
  QVBoxLayout *mainLayout = new QVBoxLayout;

  // Create parameter tables tab widget
  paramTabWidget = new QTabWidget;

  // Create tables for different parameter types
  levelParamsTable = new QTableWidget;
  setupParameterTable(levelParamsTable, "Level Parameters");
  paramTabWidget->addTab(levelParamsTable, "Level Parameters");

  normParamsTable = new QTableWidget;
  setupParameterTable(normParamsTable, "Normalization Parameters", true);
  paramTabWidget->addTab(normParamsTable, "Normalization");

  shiftParamsTable = new QTableWidget;
  setupParameterTable(shiftParamsTable, "Energy Shift Parameters", true);
  paramTabWidget->addTab(shiftParamsTable, "Energy Shifts");

  // THM coherent backgrounds (cbackground=): a tab only when there are some.
  cbkgParamsTable = new QTableWidget;
  setupParameterTable(cbkgParamsTable, "THM Coherent Background");

  mainLayout->addWidget(paramTabWidget);

  // Create button group
  QHBoxLayout *buttonLayout = new QHBoxLayout;

  refreshButton = new QPushButton("Refresh from Current");
  loadButton = new QPushButton("Load from .sav file");
  clearLimitsButton = new QPushButton("Clear Limits");

  buttonLayout->addWidget(refreshButton);
  buttonLayout->addWidget(loadButton);
  buttonLayout->addWidget(clearLimitsButton);
  buttonLayout->addStretch();

  mainLayout->addLayout(buttonLayout);

  // Connect signals
  connect(refreshButton, SIGNAL(clicked()), this, SLOT(refreshParameters()));
  connect(loadButton, SIGNAL(clicked()), this, SLOT(loadSettings()));
  connect(clearLimitsButton, SIGNAL(clicked()), this, SLOT(clearLimits()));

  setLayout(mainLayout);
}

void FittingTab::setupParameterTable(QTableWidget *table, const QString &title, bool priorCentreColumn) {
  // Set up columns (removed Fixed column - show only unfixed parameters)
  QStringList headers;
  headers << "Parameter" << "Value" << "Lower Limit" << "Upper Limit"
          << "Error" << "Fit Error" << "Use as Nuisance";
  // Normalizations and energy shifts: the centre of the prior, when it is not
  // the Value (an explicit prior_centre row; empty = the Value, as classic).
  if (priorCentreColumn) headers << "Prior Centre";

  table->setColumnCount(headers.size());
  table->setHorizontalHeaderLabels(headers);

  // Set column widths
  table->horizontalHeader()->setStretchLastSection(true);
  table->setColumnWidth(0, 150);  // Parameter name
  table->setColumnWidth(1, 100);  // Value (reduced width for levels)
  table->setColumnWidth(2, 100);  // Lower limit
  table->setColumnWidth(3, 100);  // Upper limit
  table->setColumnWidth(4, 100);  // Error
  table->setColumnWidth(5, 100);  // Fit Error
  table->setColumnWidth(6, 120);  // Nuisance checkbox
  if (priorCentreColumn) {
    table->setColumnWidth(7, 100);  // Prior centre
    table->horizontalHeaderItem(7)->setToolTip(
        "Centre of the prior on this parameter. Empty: the value in the Segments tab "
        "is both the start value and the centre (classic). A number keeps the prior "
        "there whatever value the fit starts from (e.g. 1 for a normalization whose "
        "fitted value is stored in the Segments tab).");
  }

  table->setAlternatingRowColors(true);
  table->setSelectionBehavior(QAbstractItemView::SelectRows);

  // Connect item change signal
  connect(table, SIGNAL(itemChanged(QTableWidgetItem *)),
          this, SLOT(parameterItemChanged(QTableWidgetItem *)));
}

void FittingTab::addParameterRow(QTableWidget *table, const FittingParameter &param) {
  // Filling a row is not an edit.  Unblocked, each setItem emits itemChanged,
  // and parameterItemChanged read the six-digit display text back as the new
  // value (or limit, or error) and pushed it into the levels model: merely
  // opening a project rounded every free level parameter.
  const QSignalBlocker blocker(table);
  int row = table->rowCount();
  table->insertRow(row);

  // Parameter name (read-only)
  QTableWidgetItem *nameItem = new QTableWidgetItem(param.name);
  nameItem->setFlags(nameItem->flags() & ~Qt::ItemIsEditable);
  table->setItem(row, 0, nameItem);

  // Value (reduced width for levels)
  table->setItem(row, 1, new QTableWidgetItem(QString::number(param.value, 'g', 6)));

  // Lower limit
  table->setItem(row, 2, new QTableWidgetItem(QString::number(param.lowerLimit, 'g', 6)));

  // Upper limit
  table->setItem(row, 3, new QTableWidgetItem(QString::number(param.upperLimit, 'g', 6)));

  // Error (for nuisance calculations)
  table->setItem(row, 4, new QTableWidgetItem(QString::number(param.error, 'g', 6)));

  // Fit Error (from fitting results, read-only)
  QTableWidgetItem *fitErrorItem = new QTableWidgetItem(QString::number(param.fitError, 'g', 6));
  fitErrorItem->setFlags(fitErrorItem->flags() & ~Qt::ItemIsEditable);
  table->setItem(row, 5, fitErrorItem);

  // Nuisance checkbox
  QTableWidgetItem *nuisanceItem = new QTableWidgetItem();
  nuisanceItem->setCheckState(param.useAsNuisance ? Qt::Checked : Qt::Unchecked);
  nuisanceItem->setFlags(nuisanceItem->flags() & ~Qt::ItemIsEditable);
  table->setItem(row, 6, nuisanceItem);

  if (table->columnCount() > 7) {
    QMap<QString, double>::const_iterator c = priorCentres_.constFind(param.name);
    table->setItem(row, 7, new QTableWidgetItem(c == priorCentres_.constEnd() ? QString() : roundTripNumber(*c)));
  }
}

void FittingTab::updateParameterTables() {
  // Clear existing tables
  levelParamsTable->setRowCount(0);
  normParamsTable->setRowCount(0);
  shiftParamsTable->setRowCount(0);
  cbkgParamsTable->setRowCount(0);

  // Populate tables with fitting parameters (only non-fixed ones)
  for (const FittingParameter &param : fittingParameters) {
    QTableWidget *targetTable = nullptr;

    if (param.category == "level") {
      targetTable = levelParamsTable;
    } else if (param.category == "norm") {
      targetTable = normParamsTable;
    } else if (param.category == "shift") {
      targetTable = shiftParamsTable;
    } else if (param.category == "cbkg") {
      targetTable = cbkgParamsTable;
    }

    if (targetTable) {
      addParameterRow(targetTable, param);
    }
  }
  showCoherentTab();
}

void FittingTab::reset() {
  fittingParameters.clear();
  savedParameterSettings.clear();
  priorCentres_.clear();
  updateParameterTables();
}

void FittingTab::refreshFromMinuitParameters() {
  fittingParameters.clear();

  // Create some example parameters to test the interface
  // In a real implementation, this would read from the current tab models
  FittingParameter param1 = {"Level_1_Energy", 0.5, 0.0, 2.0, 0.01, 0.01, false, "level", 0, 0, -1};
  FittingParameter param2 = {"Level_1_Width_Ch1", 0.1, 0.0, 1.0, 0.005, 0.005, false, "level", 1, 0, 0};
  FittingParameter param3 = {"Segment_1_Norm", 1.0, 0.5, 2.0, 0.05, 0.05, true, "norm", 2, -1, -1};
  FittingParameter param4 = {"Segment_1_Shift", 0.0, -0.1, 0.1, 0.001, 0.001, true, "shift", 3, -1, -1};

  fittingParameters.append(param1);
  fittingParameters.append(param2);
  fittingParameters.append(param3);
  fittingParameters.append(param4);

  updateParameterTables();
}

void FittingTab::refreshParameters() {
  populateFromCurrentGUIState();
}

void FittingTab::setTabReferences(LevelsTab *levelsTab, SegmentsTab *segmentsTab) {
  levelsTab_ = levelsTab;
  segmentsTab_ = segmentsTab;

  // Connect signals from SegmentsDataModel to update FittingTab when values change
  if (segmentsTab_) {
    SegmentsDataModel *segmentsModel = segmentsTab_->getSegmentsDataModel();
    if (segmentsModel) {
      connect(segmentsModel, SIGNAL(normalizationChanged(int, double)),
              this, SLOT(onSegmentNormalizationChanged(int, double)));
      connect(segmentsModel, SIGNAL(energyShiftChanged(int, double)),
              this, SLOT(onSegmentEnergyShiftChanged(int, double)));
      connect(segmentsModel, SIGNAL(normalizationErrorChanged(int, double)),
              this, SLOT(onSegmentNormalizationErrorChanged(int, double)));
      connect(segmentsModel, SIGNAL(energyShiftErrorChanged(int, double)),
              this, SLOT(onSegmentEnergyShiftErrorChanged(int, double)));
      connect(segmentsModel, SIGNAL(normalizationVaryChanged(int, bool)),
              this, SLOT(onSegmentNormalizationVaryChanged(int, bool)));
      connect(segmentsModel, SIGNAL(energyShiftVaryChanged(int, bool)),
              this, SLOT(onSegmentEnergyShiftVaryChanged(int, bool)));
    }
  }

  // Populate parameters from current tab state
  populateFromCurrentGUIState();
}

// The channels of one level, as indices into the channels model, in the order
// they are written to (and read back from) the file.  AZURE2 numbers a level's
// widths by this position: width_N_k is the k-th of them.
static QList<int> channelsOfLevel(const QList<ChannelsData> &channels, int levelIndex) {
  QList<int> result;
  for (int ch = 0; ch < channels.size(); ch++)
    if (channels.at(ch).levelIndex == levelIndex) result.append(ch);
  return result;
}

// The levels AZURE2 builds from a file whose <levels> block lists the given
// model levels in the given order, in the order it numbers them ("Level N",
// energy_N): CNuc::Fill keeps only active lines, opens a J-group at the first
// level of each (J, parity) and appends later levels to their group, and
// FillMnParams counts group by group.  A level without channels writes no line
// and so does not exist for the engine.
//
// These numbers, not the model's row indices, are what the "Level N ..." names
// in <parameterSettings> mean.  The model keeps the order of the file that was
// opened while the GUI writes its levels sorted by J, parity and energy, so a
// name built from a row index described one level in the file just saved and
// another one after reopening it.
QList<int> FittingTab::engineLevelOrder(const QList<int> &fileOrder) {
  QList<LevelsData> levels;
  QList<ChannelsData> channels;
  if (levelsTab_) {
    if (LevelsModel *m = levelsTab_->getLevelsModel()) levels = m->getLevels();
    if (ChannelsModel *m = levelsTab_->getChannelsModel()) channels = m->getChannels();
  }
  QList<QList<int>> groups;
  for (int la : fileOrder) {
    if (la < 0 || la >= levels.size() || levels[la].isActive == 0) continue;
    if (channelsOfLevel(channels, la).isEmpty()) continue;
    int g = 0;
    while (g < groups.size() && !(levels[groups[g].first()].jValue == levels[la].jValue &&
                                  levels[groups[g].first()].piValue == levels[la].piValue))
      g++;
    if (g == groups.size()) groups.append(QList<int>());
    groups[g].append(la);
  }
  QList<int> order;
  for (const QList<int> &group : groups) order.append(group);
  return order;
}

FittingParameter FittingTab::segmentParameter(int i, const SegmentsDataData &segment, bool norm) {
  FittingParameter param;
  // Use params.sav naming convention: segment_x_norm, segment_x_energy_shift
  param.name = QString(norm ? "segment_%1_norm" : "segment_%1_energy_shift").arg(i + 1);
  param.value = norm ? segment.dataNorm : segment.energyShift;
  param.lowerLimit = 0;
  param.upperLimit = 0;
  param.error = norm ? segment.dataNormError : segment.energyShiftError;
  param.fitError = 0.0;  // No fit error initially
  // Active only when the parameter varies
  param.useAsNuisance = norm ? (segment.varyNorm == 1) : (segment.varyEnergyShift == 1);
  param.category = norm ? "norm" : "shift";
  param.levelIndex = -1;
  param.channelIndex = i;  // Store segment index for reverse lookup
  return param;
}

void FittingTab::populateFromCurrentGUIState() {
  if (!levelsTab_ || !segmentsTab_) return;

  fittingParameters.clear();

  // Get level and channel data
  LevelsModel *levelsModel = levelsTab_->getLevelsModel();
  ChannelsModel *channelsModel = levelsTab_->getChannelsModel();

  if (levelsModel && channelsModel) {
    QList<LevelsData> levels = levelsModel->getLevels();
    QList<ChannelsData> channels = channelsModel->getChannels();

    // AZURE2's order and numbering: the levels of the file this tab's rows
    // will be saved into, energy first, then that level's widths.
    const QList<int> order = engineLevelOrder(levelsTab_->writeOrder());
    for (int n = 0; n < order.size(); n++) {
      const int levelIndex = order.at(n);
      const LevelsData &level = levels[levelIndex];

      // Add energy parameter if not fixed (isFixed == 0 means not fixed)
      if (level.isFixed == 0) {
        FittingParameter energyParam;
        energyParam.name = QString("Level %1 Energy (MeV)").arg(n + 1);
        energyParam.value = level.energy;
        energyParam.lowerLimit = 0;  // Default limits
        energyParam.upperLimit = 0;
        energyParam.error = 0.01;    // Default error
        energyParam.fitError = 0.0;  // No fit error initially
        energyParam.useAsNuisance = false;
        energyParam.category = "level";
        energyParam.levelIndex = levelIndex;
        energyParam.channelIndex = -1;  // Energy parameter

        fittingParameters.append(energyParam);
      }

      // Now add all width parameters for this level
      const QList<int> levelChannels = channelsOfLevel(channels, levelIndex);
      for (int k = 0; k < levelChannels.size(); k++) {
        const int channelIndex = levelChannels.at(k);
        const ChannelsData &channel = channels[channelIndex];

        if (channel.isFixed == 0 && channel.reducedWidth != 0.0) {
          FittingParameter widthParam;
          // Keep "Width" in the name (parameter matching keys on it);
          // the unit reflects the channel's input convention.  The channel
          // is numbered within its level, as AZURE2 does (width_N_k).
          widthParam.name = QString("Level %1 Channel %2 Width (%3)")
                                .arg(n + 1)
                                .arg(k + 1)
                                .arg(channel.gammaIsRWA == 1 ? "MeV^(1/2)" : "eV");

          // LevelsTab stores the value in the channel's input convention:
          // physical width/ANC, or the reduced width amplitude if flagged
          widthParam.value = channel.reducedWidth;
          widthParam.lowerLimit = 0;
          widthParam.upperLimit = 0;                      // Default upper limit
          widthParam.error = channel.reducedWidth * 0.1;  // Default 10% error
          widthParam.fitError = 0.0;                      // No fit error initially
          widthParam.useAsNuisance = false;
          widthParam.category = "level";
          widthParam.levelIndex = levelIndex;
          widthParam.channelIndex = channelIndex;

          fittingParameters.append(widthParam);
        }
      }
    }
  }

  // Get normalization and shift parameters from SegmentsDataModel
  // IMPORTANT: Follow same order as EData::FillMnParams() - ALL norms first, then ALL shifts
  SegmentsDataModel *segmentsModel = segmentsTab_->getSegmentsDataModel();
  if (segmentsModel) {
    QList<SegmentsDataData> segments = segmentsModel->getLines();

    // First pass: Add normalization parameters (only for ACTIVE segments with varyNorm, like EData::FillMnParams)
    for (int i = 0; i < segments.size(); i++) {
      const SegmentsDataData &segment = segments[i];

      // Only process active segments
      // Always show norm parameter for active segments, but only add to Minuit if varyNorm=1 (like EData)
      if (segment.isActive == 1) fittingParameters.append(segmentParameter(i, segment, true));
    }

    // Second pass: Add energy shift parameters (only for ACTIVE segments, like EData::FillMnParams)
    for (int i = 0; i < segments.size(); i++) {
      const SegmentsDataData &segment = segments[i];

      // Only process active segments
      // Always add energy shift parameter to GUI and Minuit for active segments (like EData), but useAsNuisance depends on varyEnergyShift
      if (segment.isActive == 1) fittingParameters.append(segmentParameter(i, segment, false));
    }
  }

  appendCoherentParameters();
  assignMinuitIndices();

  // Apply parameter settings from saved configuration (limits, errors, etc.)
  applyParameterSettings();

  updateParameterTables();
}

void FittingTab::showEvent(QShowEvent *event) {
  QWidget::showEvent(event);
  // Reflect free/fix changes to normalizations or energy shifts made in the
  // Segments tab while this tab was hidden.
  syncSegmentVaryStates();
}

// Number the parameters exactly the way AZURE2 does, so the minuit_index written
// into <parameterSettings> names the parameter it was meant for.
//
// AZURE2 builds its Minuit vector as CNuc::FillMnParams() followed by
// EData::FillMnParams(), then enumerates the entries it will apply settings to as
// "not fixed, or a segment parameter".  That means: a level energy or width counts
// only when it is free, a normalization exists at all only when its segment varies
// it, and an energy shift is always present (fixed when not varied) and always
// counted.  Anything else -- in particular skipping the shifts of segments that do
// not vary theirs -- shifts every later index by one and silently applies one
// parameter's limits and nuisance prior to another.
void FittingTab::assignMinuitIndices() {
  QList<LevelsData> levels;
  QList<ChannelsData> channels;
  if (levelsTab_) {
    if (LevelsModel *m = levelsTab_->getLevelsModel()) levels = m->getLevels();
    if (ChannelsModel *m = levelsTab_->getChannelsModel()) channels = m->getChannels();
  }
  QList<SegmentsDataData> segments;
  if (segmentsTab_) {
    if (SegmentsDataModel *m = segmentsTab_->getSegmentsDataModel()) segments = m->getLines();
  }

  int paramIndex = 0;
  for (int i = 0; i < fittingParameters.size(); i++) {
    FittingParameter &param = fittingParameters[i];
    bool counted = false;
    if (param.category == "level") {
      // Levels excluded from the calculation never reach CNuc, so they
      // occupy no parameter slot at all.
      bool levelActive = (param.levelIndex >= 0 && param.levelIndex < levels.size() &&
                          levels[param.levelIndex].isActive != 0);
      if (!levelActive)
        counted = false;
      else if (param.channelIndex < 0) {
        counted = (levels[param.levelIndex].isFixed == 0);
      } else {
        counted = (param.channelIndex < channels.size() &&
                   channels[param.channelIndex].isFixed == 0 &&
                   channels[param.channelIndex].reducedWidth != 0.0);
      }
    } else if (param.category == "norm" || param.category == "shift") {
      bool segmentActive = (param.channelIndex >= 0 && param.channelIndex < segments.size() &&
                            segments[param.channelIndex].isActive == 1);
      counted = segmentActive &&
          (param.category == "shift" || segments[param.channelIndex].varyNorm == 1);
    } else if (param.category == "cbkg") {
      counted = true;  // only the free ones are listed (appendCoherentParameters)
    }
    param.minuitIndex = counted ? paramIndex++ : -1;
  }
}

// Sync the "Use as Nuisance" state of the norm/shift parameters to the segments'
// current vary flags, updating the tables in place (no rebuild, so any Fitting-tab
// limit edits are preserved).
void FittingTab::syncSegmentVaryStates() {
  if (!segmentsTab_) return;
  SegmentsDataModel *segmentsModel = segmentsTab_->getSegmentsDataModel();
  if (!segmentsModel) return;
  QList<SegmentsDataData> segments = segmentsModel->getLines();
  for (int i = 0; i < fittingParameters.size(); i++) {
    FittingParameter &p = fittingParameters[i];
    if (p.category != "norm" && p.category != "shift") continue;
    int seg = p.channelIndex;  // segment index is stored in channelIndex
    if (seg < 0 || seg >= segments.size()) continue;
    bool vary = (p.category == "norm") ? (segments[seg].varyNorm == 1)
                                       : (segments[seg].varyEnergyShift == 1);
    if (p.useAsNuisance != vary) {
      p.useAsNuisance = vary;
      updateParameterTableCheckbox(p.name, vary);
    }
  }
  // Freeing or fixing a normalization changes which parameters AZURE2 puts in
  // its Minuit vector, so the recorded indices have to be renumbered with it.
  assignMinuitIndices();
}

double FittingTab::transformRWAParameterToPhysical(const QString &paramName, double rwaValue) {
  // Transform RWA parameter to physical using proper R-Matrix transformation via AZURESetup

  // Find the parent AZURESetup widget
  AZURESetup *azureSetup = nullptr;
  QWidget *parent = this->parentWidget();
  while (parent != nullptr) {
    azureSetup = qobject_cast<AZURESetup *>(parent);
    if (azureSetup != nullptr) {
      break;
    }
    parent = parent->parentWidget();
  }

  if (azureSetup != nullptr) {
    // Use the proper RWA to Physical conversion from AZURESetup
    return azureSetup->ConvertRWAToPhysical(paramName, rwaValue);
  } else {
    // Fallback: return the original value if we can't find AZURESetup
    return rwaValue;
  }
}

void FittingTab::applyParameterSettings() {
  // Apply saved parameter settings (limits, errors, etc.) to current parameters
  for (int i = 0; i < fittingParameters.size(); i++) {
    const FittingParameter &current = fittingParameters[i];
    for (const FittingParameter &saved : savedParameterSettings) {
      // Level entries were tied to a level and channel of the model when the
      // file was read (see readParameterSettings); their names are numbers
      // that depend on the order the levels happen to be in.
      bool same = (current.category == "level" || saved.category == "level")
                      ? (current.category == saved.category && current.levelIndex == saved.levelIndex &&
                         current.channelIndex == saved.channelIndex)
                      : (current.name == saved.name);
      if (same) {
        // Apply saved settings but keep current value from models
        fittingParameters[i].lowerLimit = saved.lowerLimit;
        fittingParameters[i].upperLimit = saved.upperLimit;
        fittingParameters[i].error = saved.error;  // Keep saved error for nuisance calculations
        // fitError is maintained from fitting results, not overwritten
        fittingParameters[i].useAsNuisance = saved.useAsNuisance;
        break;
      }
    }
  }
}

void FittingTab::parameterItemChanged(QTableWidgetItem *item) {
  // Handle parameter changes
  int row = item->row();
  int col = item->column();
  QTableWidget *table = qobject_cast<QTableWidget *>(item->tableWidget());
  if (!table) return;

  // Get parameter name from first column
  QTableWidgetItem *nameItem = table->item(row, 0);
  if (!nameItem) return;
  QString paramName = nameItem->text();

  // Find the parameter in our settings
  int paramIndex = -1;
  for (int i = 0; i < fittingParameters.size(); i++) {
    if (fittingParameters[i].name == paramName) {
      paramIndex = i;
      break;
    }
  }

  if (paramIndex == -1) return;  // Parameter not found

  // Column mapping: 0=Name, 1=Value, 2=Lower, 3=Upper, 4=Error, 5=Fit Error (read-only), 6=Nuisance
  if (col == 1) {  // Value (reduced width) changed
    bool ok;
    double value = item->text().toDouble(&ok);
    if (!ok) {
      QMessageBox::warning(this, "Invalid Input", "Please enter a valid number.");
      item->setText(QString::number(fittingParameters[paramIndex].value, 'g', 6));
      return;
    }

    fittingParameters[paramIndex].value = value;

    // Update the corresponding tab with the new value
    if (fittingParameters[paramIndex].category == "cbkg") {
      QMap<QString, double> one;
      one[paramName] = value;
      applyCoherentValues(one);  // the cbackground= of its experiment
    } else {
      updateParameterInOtherTabs(paramName, fittingParameters[paramIndex]);
    }

  } else if (col >= 2 && col <= 4) {  // Lower limit, upper limit, or error changed
    bool ok;
    double value = item->text().toDouble(&ok);
    if (!ok) {
      QMessageBox::warning(this, "Invalid Input", "Please enter a valid number.");
      // Reset to previous value
      if (col == 2)
        item->setText(QString::number(fittingParameters[paramIndex].lowerLimit, 'g', 6));
      else if (col == 3)
        item->setText(QString::number(fittingParameters[paramIndex].upperLimit, 'g', 6));
      else if (col == 4)
        item->setText(QString::number(fittingParameters[paramIndex].error, 'g', 6));
      return;
    }

    // Update parameter setting
    if (col == 2)
      fittingParameters[paramIndex].lowerLimit = value;
    else if (col == 3)
      fittingParameters[paramIndex].upperLimit = value;
    else if (col == 4) {
      fittingParameters[paramIndex].error = value;
      // The user edited the uncertainty itself, so write that column and
      // nothing else: dataNormError (10) or energyShiftError (15).
      const FittingParameter &p = fittingParameters[paramIndex];
      writeSegmentColumn(p, p.category == "norm" ? 10 : 15, value);
    }

  } else if (col == 7) {  // Prior centre (norms and shifts): empty = the Value
    const QString text = item->text().trimmed();
    if (text.isEmpty()) {
      priorCentres_.remove(paramName);
      return;
    }
    bool ok;
    const double centre = text.toDouble(&ok);
    const bool norm = fittingParameters[paramIndex].category == "norm";
    if (!ok || !std::isfinite(centre) || (norm && centre <= 0.0)) {
      QMessageBox::warning(this, "Invalid Input",
                           norm ? "Enter a positive number, or leave the cell empty to centre the prior on the Value."
                                : "Enter a number, or leave the cell empty to centre the prior on the Value.");
      const QSignalBlocker blocker(table);
      QMap<QString, double>::const_iterator c = priorCentres_.constFind(paramName);
      item->setText(c == priorCentres_.constEnd() ? QString() : roundTripNumber(*c));
      return;
    }
    priorCentres_[paramName] = centre;
  } else if (col == 6) {  // Nuisance checkbox
    fittingParameters[paramIndex].useAsNuisance = (item->checkState() == Qt::Checked);
    // Likewise, the user changed the fit freedom, so write only the vary
    // column: varyNorm (11) or varyEnergyShift (16).
    const FittingParameter &p = fittingParameters[paramIndex];
    writeSegmentColumn(p, p.category == "norm" ? 11 : 16, p.useAsNuisance ? 1 : 0);
  }
}

void FittingTab::loadSettings() {
  // The RWA<->physical transformation uses the current GUI state (levels,
  // channels), so refresh from it first -- a stale state gives a wrong result.
  populateFromCurrentGUIState();

  QString filename = QFileDialog::getOpenFileName(this,
                                                  "Load Parameters from .sav file", "", "AZURE2 Parameter Files (*.sav);;All Files (*)");

  if (!filename.isEmpty()) {
    // Back up the current .azr before applying the transformation.
    QString backupNote;
    if (config_ && !config_->configfile.empty()) {
      QString azr = QString::fromStdString(config_->configfile);
      QString bk = azr + ".bk";
      QFile::remove(bk);
      if (QFile::copy(azr, bk))
        backupNote = QString("\n\nThe original file was backed up as:\n%1\n\n"
                             "The parameter transformation can be unstable when the parameters are "
                             "unphysical; restore this backup if the loaded values look wrong.")
                         .arg(bk);
    }

    QFile file(filename);
    if (file.open(QIODevice::ReadOnly | QIODevice::Text)) {
      QTextStream in(&file);

      // Parse .sav file and create lookup map
      QMap<QString, QPair<double, double>> savParams;  // paramName -> (value, error)

      while (!in.atEnd()) {
        QString line = in.readLine().trimmed();
        if (line.isEmpty()) continue;

        QStringList parts = line.split(QRegExp("\\s+"), Qt::SkipEmptyParts);
        if (parts.size() >= 3) {
          QString paramName = parts[0];
          double value = parts[1].toDouble();
          double error = parts[2].toDouble();

          // If "_rwa" in param.name, remove it
          if (paramName.endsWith("_rwa")) {
            paramName.chop(4);  // Remove last 4 characters
          }

          savParams[paramName] = qMakePair(value, error);
        }
      }
      file.close();

      // THM coherent backgrounds: their values go into the cbackground= of
      // their experiments (the <thm> block), not into a tab's model.
      {
        QMap<QString, double> coherent;
        for (const QString &key : savParams.keys())
          if (key.startsWith("cbkg_")) {
            coherent[key] = savParams[key].first;
            savParams.remove(key);
          }
        if (!coherent.isEmpty()) applyCoherentValues(coherent);
      }

      // Now populate ALL parameters from the current GUI state (including fixed/non-varying)
      // This is needed for proper RWA conversion of all level parameters
      if (levelsTab_ && segmentsTab_) {
        fittingParameters.clear();

        // Get level and channel data
        LevelsModel *levelsModel = levelsTab_->getLevelsModel();
        ChannelsModel *channelsModel = levelsTab_->getChannelsModel();

        if (levelsModel && channelsModel) {
          QList<LevelsData> levels = levelsModel->getLevels();
          QList<ChannelsData> channels = channelsModel->getChannels();

          // Add ALL level parameters (including fixed) of the levels AZURE2
          // has, in its order and with its numbers -- the .sav file's
          // energy_N and width_N_k refer to those.
          const QList<int> order = engineLevelOrder(levelsTab_->writeOrder());
          for (int n = 0; n < order.size(); n++) {
            const int levelIndex = order.at(n);
            const LevelsData &level = levels[levelIndex];

            // Add ALL energy parameters
            FittingParameter energyParam;
            energyParam.name = QString("Level %1 Energy (MeV)").arg(n + 1);
            energyParam.value = level.energy;
            energyParam.lowerLimit = 0;
            energyParam.upperLimit = 0;
            energyParam.error = 0.01;
            energyParam.fitError = 0.0;
            energyParam.useAsNuisance = false;
            energyParam.category = "level";
            energyParam.levelIndex = levelIndex;
            energyParam.channelIndex = -1;

            fittingParameters.append(energyParam);

            // Add ALL width parameters for this level (including those with zero width or fixed)
            const QList<int> levelChannels = channelsOfLevel(channels, levelIndex);
            for (int k = 0; k < levelChannels.size(); k++) {
              const int channelIndex = levelChannels.at(k);
              const ChannelsData &channel = channels[channelIndex];
              FittingParameter widthParam;
              widthParam.name = QString("Level %1 Channel %2 Width (%3)")
                                    .arg(n + 1)
                                    .arg(k + 1)
                                    .arg(channel.gammaIsRWA == 1 ? "MeV^(1/2)" : "eV");

              widthParam.value = channel.reducedWidth;
              widthParam.lowerLimit = 0;
              widthParam.upperLimit = 0;
              widthParam.error = (channel.reducedWidth != 0.0) ? channel.reducedWidth * 0.1 : 0.01;
              widthParam.fitError = 0.0;
              widthParam.useAsNuisance = false;
              widthParam.category = "level";
              widthParam.levelIndex = levelIndex;
              widthParam.channelIndex = channelIndex;

              fittingParameters.append(widthParam);
            }
          }
        }

        // Add ALL normalization and shift parameters from all active segments
        SegmentsDataModel *segmentsModel = segmentsTab_->getSegmentsDataModel();
        if (segmentsModel) {
          QList<SegmentsDataData> segments = segmentsModel->getLines();

          // Add ALL norm parameters for active segments
          for (int i = 0; i < segments.size(); i++) {
            const SegmentsDataData &segment = segments[i];

            if (segment.isActive == 1) {
              FittingParameter normParam;
              normParam.name = QString("segment_%1_norm").arg(i + 1);
              normParam.value = segment.dataNorm;
              normParam.lowerLimit = 0;
              normParam.upperLimit = 0;
              normParam.error = segment.dataNormError;
              normParam.fitError = 0.0;
              normParam.useAsNuisance = (segment.varyNorm == 1);
              normParam.category = "norm";
              normParam.levelIndex = -1;
              normParam.channelIndex = i;

              fittingParameters.append(normParam);
            }
          }

          // Add ALL energy shift parameters for active segments
          for (int i = 0; i < segments.size(); i++) {
            const SegmentsDataData &segment = segments[i];

            if (segment.isActive == 1) {
              FittingParameter shiftParam;
              shiftParam.name = QString("segment_%1_energy_shift").arg(i + 1);
              shiftParam.value = segment.energyShift;
              shiftParam.lowerLimit = 0;
              shiftParam.upperLimit = 0;
              shiftParam.error = segment.energyShiftError;
              shiftParam.fitError = 0.0;
              shiftParam.useAsNuisance = (segment.varyEnergyShift == 1);
              shiftParam.category = "shift";
              shiftParam.levelIndex = -1;
              shiftParam.channelIndex = i;

              fittingParameters.append(shiftParam);
            }
          }
        }
        appendCoherentParameters();  // with the values just loaded
        assignMinuitIndices();
      }

      // Separate handling: RWA parameters (levels) vs direct assignment (norms/shifts)
      QMap<QString, QPair<double, double>> rwaParamMap;
      QList<QPair<int, QString>> levelParameterMapping;  // For level parameters that need RWA conversion

      // Parse norms and shifts using segment names from .sav file structure
      QMap<QString, QPair<double, double>> normsFromSav;   // segment_name -> (value, error)
      QMap<QString, QPair<double, double>> shiftsFromSav;  // segment_name -> (value, error)

      // Extract norms and shifts using segment names/numbers from .sav file
      for (const QString &key : savParams.keys()) {
        QPair<double, double> savData = savParams[key];

        // Match norm parameters and extract segment identifier
        if (key.contains("norm", Qt::CaseInsensitive) && !key.contains("_rwa")) {
          // Try multiple patterns to extract segment identifier
          QRegExp segmentNumRegex("segment_(\\d+)_norm", Qt::CaseInsensitive);
          QRegExp normNumRegex("norm_(\\d+)", Qt::CaseInsensitive);
          QRegExp segmentNameRegex("([^_]+)_norm", Qt::CaseInsensitive);

          QString segmentId;
          if (segmentNumRegex.indexIn(key) != -1) {
            segmentId = segmentNumRegex.cap(1);  // Just the number part
          } else if (normNumRegex.indexIn(key) != -1) {
            segmentId = normNumRegex.cap(1);  // Just the number part
          } else if (segmentNameRegex.indexIn(key) != -1) {
            segmentId = segmentNameRegex.cap(1);  // Segment name part
          }

          if (!segmentId.isEmpty()) {
            normsFromSav[segmentId] = savData;
          }
        }
        // Match energy shift parameters and extract segment identifier
        else if (key.contains("shift", Qt::CaseInsensitive) && !key.contains("_rwa")) {
          // Try multiple patterns to extract segment identifier
          QRegExp segmentNumRegex("segment_(\\d+)_energy_shift", Qt::CaseInsensitive);
          QRegExp segmentShiftRegex("segment_(\\d+)_shift", Qt::CaseInsensitive);
          QRegExp shiftNumRegex("shift_(\\d+)", Qt::CaseInsensitive);
          QRegExp segmentNameRegex("([^_]+)_(?:energy_)?shift", Qt::CaseInsensitive);

          QString segmentId;
          if (segmentNumRegex.indexIn(key) != -1) {
            segmentId = segmentNumRegex.cap(1);  // Just the number part
          } else if (segmentShiftRegex.indexIn(key) != -1) {
            segmentId = segmentShiftRegex.cap(1);  // Just the number part
          } else if (shiftNumRegex.indexIn(key) != -1) {
            segmentId = shiftNumRegex.cap(1);  // Just the number part
          } else if (segmentNameRegex.indexIn(key) != -1) {
            segmentId = segmentNameRegex.cap(1);  // Segment name part
          }

          if (!segmentId.isEmpty()) {
            shiftsFromSav[segmentId] = savData;
          }
        }
      }

      // First pass: collect level parameters that need RWA conversion
      // CRITICAL FIX: Extract OLD level structure from .sav RWA parameter names
      // This maps: energyIndex -> number of channels in that level from the .sav file
      QMap<int, int> oldLevelChannelCounts;  // energyIndex -> channelCount
      for (const QString &key : savParams.keys()) {
        // Match width parameters like "width_1_1", "width_1_2", "width_2_1", etc.
        QRegExp widthRegex("^width_(\\d+)_(\\d+)$");
        if (widthRegex.indexIn(key) != -1) {
          int energyIndex = widthRegex.cap(1).toInt();
          int widthIndex = widthRegex.cap(2).toInt();

          // Track the maximum width index for each energy level
          if (!oldLevelChannelCounts.contains(energyIndex) || oldLevelChannelCounts[energyIndex] < widthIndex) {
            oldLevelChannelCounts[energyIndex] = widthIndex;
          }
        }
      }

      for (int i = 0; i < fittingParameters.size(); i++) {
        FittingParameter &param = fittingParameters[i];

        if (param.category == "level") {
          QString matchKey = findMatchingParameterKey(param, savParams.keys());
          if (!matchKey.isEmpty() && savParams.contains(matchKey)) {
            QPair<double, double> savData = savParams[matchKey];
            levelParameterMapping.append(qMakePair(i, matchKey));
            rwaParamMap[matchKey] = savData;
          }
        }
      }

      // Perform batch transformation using AZURESetup
      QMap<QString, QPair<double, double>> physicalParamMap;
      if (!rwaParamMap.isEmpty()) {
        // Find the parent AZURESetup widget for batch conversion
        AZURESetup *azureSetup = nullptr;
        QWidget *parent = this->parentWidget();
        while (parent != nullptr) {
          azureSetup = qobject_cast<AZURESetup *>(parent);
          if (azureSetup != nullptr) break;
          parent = parent->parentWidget();
        }

        if (azureSetup != nullptr) {
          // Use AZUREAPI-style batch transformation with OLD compound structure
          // Build complete RWA parameter vector
          QStringList rwaParamNames = rwaParamMap.keys();
          std::vector<double> rwaValues, rwaValuesPlusError, rwaValuesMinusError;

          for (const QString &paramName : rwaParamNames) {
            QPair<double, double> rwaData = rwaParamMap[paramName];
            rwaValues.push_back(rwaData.first);
            rwaValuesPlusError.push_back(rwaData.first + rwaData.second);
            rwaValuesMinusError.push_back(rwaData.first - rwaData.second);
          }

          // CRITICAL FIX: Use the OLD level structure for transformation
          // Pass the old level channel counts to ensure correct RWA-to-physical conversion
          std::vector<double> physicalValues = azureSetup->BatchConvertRWAToPhysicalWithOldStructure(
              rwaParamNames, rwaValues, oldLevelChannelCounts);
          std::vector<double> physicalValuesPlusError = azureSetup->BatchConvertRWAToPhysicalWithOldStructure(
              rwaParamNames, rwaValuesPlusError, oldLevelChannelCounts);
          std::vector<double> physicalValuesMinusError = azureSetup->BatchConvertRWAToPhysicalWithOldStructure(
              rwaParamNames, rwaValuesMinusError, oldLevelChannelCounts);

          // Store the transformed results
          for (int i = 0; i < rwaParamNames.size(); i++) {
            QString paramName = rwaParamNames[i];
            double physicalValue = physicalValues[i];
            double physicalErrorUp = std::abs(physicalValuesPlusError[i] - physicalValue);
            double physicalErrorDown = std::abs(physicalValue - physicalValuesMinusError[i]);
            double physicalError = std::max(physicalErrorUp, physicalErrorDown);

            // Check for NaN and replace with 0 if found
            if (std::isnan(physicalValue)) {
              physicalValue = 0.0;
              QMessageBox::warning(this, "Invalid Value",
                                   QString("Transformed physical value for parameter '%1' is NaN. Setting to 0.")
                                       .arg(paramName));
            }
            if (std::isnan(physicalError)) {
              physicalError = 0.0;
            }

            physicalParamMap[paramName] = qMakePair(physicalValue, physicalError);
          }
        } else {
          // Fallback: use RWA values as physical values if transformation unavailable
          for (const QString &paramName : rwaParamMap.keys()) {
            physicalParamMap[paramName] = rwaParamMap[paramName];
          }
        }
      }

      // Second pass: apply transformed values and direct assignments
      int updatedCount = 0;

      // Apply transformed level parameters
      for (const QPair<int, QString> &mapping : levelParameterMapping) {
        int paramIndex = mapping.first;
        QString matchKey = mapping.second;
        FittingParameter &param = fittingParameters[paramIndex];

        if (physicalParamMap.contains(matchKey)) {
          QPair<double, double> physicalData = physicalParamMap[matchKey];

          // Check for NaN and replace with 0 if found
          double physicalValue = physicalData.first;
          double physicalError = physicalData.second;
          if (std::isnan(physicalValue)) {
            physicalValue = 0.0;
            // Output a warning if needed
            QMessageBox::warning(this, "Invalid Value",
                                 QString("Transformed physical value for parameter '%1' is NaN. Setting to 0.")
                                     .arg(param.name));
          }
          if (std::isnan(physicalError)) {
            physicalError = 0.0;
          }

          // Update parameter value and fit error (keep original error for nuisance calculations)
          param.value = physicalValue;
          param.fitError = physicalError;  // Store fit error separately from nuisance error

          // Update the underlying models with converted values
          if (param.name.contains("Width") && param.channelIndex >= 0) {
            if (levelsTab_) {
              ChannelsModel *channelsModel = levelsTab_->getChannelsModel();
              if (channelsModel) {
                QList<ChannelsData> channels = channelsModel->getChannels();
                for (int j = 0; j < channels.size(); j++) {
                  if (j == param.channelIndex) {
                    QModelIndex modelIndex = channelsModel->index(j, 6);
                    channelsModel->setData(modelIndex, param.value, Qt::EditRole);
                    break;
                  }
                }
              }
            }
          }

          // Propagate changes to other tabs
          updateParameterInOtherTabs(param.name, param);
          updatedCount++;
        }
      }

      // Apply norms directly from .sav file using segment names/identifiers
      if (segmentsTab_) {
        SegmentsDataModel *segmentsModel = segmentsTab_->getSegmentsDataModel();
        if (segmentsModel) {
          QList<SegmentsDataData> segments = segmentsModel->getLines();

          // Temporarily disconnect signals to prevent double updates
          segmentsModel->blockSignals(true);

          // Match segments by name/identifier from .sav file
          for (auto it = normsFromSav.begin(); it != normsFromSav.end(); ++it) {
            QString segmentId = it.key();  // Segment identifier from .sav file
            QPair<double, double> normData = it.value();

            // Try to find matching segment by different methods
            int segmentIndex = -1;

            for (int i = 0; i < segments.size(); i++) {
              const SegmentsDataData &segment = segments[i];

              // Method 1: Match by segment number (convert segmentId to number and check 1-based index)
              bool segmentIdIsNumber = false;
              int segmentNumber = segmentId.toInt(&segmentIdIsNumber);
              if (segmentIdIsNumber && (segmentNumber == i + 1)) {
                segmentIndex = i;
                break;
              }

              // Method 2: Match by filename (if segment has a dataFile)
              if (!segment.dataFile.isEmpty()) {
                QString fileBaseName = QFileInfo(segment.dataFile).baseName();
                if (fileBaseName == segmentId) {
                  segmentIndex = i;
                  break;
                }

                // Also try matching with file extension
                QString fileName = QFileInfo(segment.dataFile).fileName();
                if (fileName == segmentId || fileName.startsWith(segmentId + ".")) {
                  segmentIndex = i;
                  break;
                }
              }

              // Method 3: Match by constructed segment name pattern
              QString constructedName = segment.dataFile.isEmpty() ? QString("%1").arg(i + 1) : QFileInfo(segment.dataFile).baseName();
              if (constructedName == segmentId) {
                segmentIndex = i;
                break;
              }
            }

            if (segmentIndex >= 0 && segmentIndex < segments.size()) {
              // Update the segments model directly
              // Check for NaN and replace with 0 if found
              double normValue = normData.first;
              if (std::isnan(normValue)) {
                normValue = 0.0;
              }

              QModelIndex normValueIndex = segmentsModel->index(segmentIndex, 9);  // dataNorm column
              segmentsModel->setData(normValueIndex, normValue, Qt::EditRole);
              // dataNormError (column 10) is deliberately left alone.  The .sav's
              // second column is Minuit's uncertainty on the fitted value; the
              // segment's error is the assumed systematic that sets the width of
              // the normalization penalty, and EData::CalcNormChiSquared reads it
              // as a *percentage* of the nominal norm.  The two are neither the
              // same quantity nor the same units, so copying one onto the other
              // silently destroys the researched systematic on every .sav load.
            }
          }

          // Re-enable signals
          segmentsModel->blockSignals(false);

          // Now manually update ALL FittingParameters from the updated model
          // Also update parameter names to match params.sav exactly
          for (int j = 0; j < fittingParameters.size(); j++) {
            FittingParameter &param = fittingParameters[j];
            if (param.category == "norm") {
              int segmentIndex = param.channelIndex;
              if (segmentIndex >= 0 && segmentIndex < segments.size()) {
                // Re-read the updated segments data
                QList<SegmentsDataData> updatedSegments = segmentsModel->getLines();
                param.value = updatedSegments[segmentIndex].dataNorm;
                // The fit uncertainty comes from the .sav itself, below, rather
                // than from the segment's error column: that column holds the
                // assumed systematic and is no longer overwritten on load.
                param.fitError = 0.0;

                // Find the corresponding params.sav name for this segment
                for (auto it = normsFromSav.begin(); it != normsFromSav.end(); ++it) {
                  QString segmentId = it.key();

                  // Check if this segmentId matches this segment
                  bool matches = false;
                  const SegmentsDataData &segment = updatedSegments[segmentIndex];

                  // Same matching logic as above
                  bool segmentIdIsNumber = false;
                  int segmentNumber = segmentId.toInt(&segmentIdIsNumber);
                  if (segmentIdIsNumber && (segmentNumber == segmentIndex + 1)) {
                    matches = true;
                  } else if (!segment.dataFile.isEmpty()) {
                    QString fileBaseName = QFileInfo(segment.dataFile).baseName();
                    QString fileName = QFileInfo(segment.dataFile).fileName();
                    if (fileBaseName == segmentId || fileName == segmentId || fileName.startsWith(segmentId + ".")) {
                      matches = true;
                    }
                  }

                  if (matches) {
                    param.fitError = it.value().second;
                    // Find the exact parameter name from params.sav
                    for (const QString &savKey : savParams.keys()) {
                      if (savKey.contains("norm", Qt::CaseInsensitive) && !savKey.contains("_rwa")) {
                        // Check if this savKey corresponds to this segmentId
                        QRegExp segmentNumRegex("segment_(\\d+)_norm", Qt::CaseInsensitive);
                        QRegExp normNumRegex("norm_(\\d+)", Qt::CaseInsensitive);
                        QRegExp segmentNameRegex("([^_]+)_norm", Qt::CaseInsensitive);

                        QString extractedId;
                        if (segmentNumRegex.indexIn(savKey) != -1) {
                          extractedId = segmentNumRegex.cap(1);
                        } else if (normNumRegex.indexIn(savKey) != -1) {
                          extractedId = normNumRegex.cap(1);
                        } else if (segmentNameRegex.indexIn(savKey) != -1) {
                          extractedId = segmentNameRegex.cap(1);
                        }

                        if (extractedId == segmentId) {
                          // Update the parameter name to match params.sav exactly
                          param.name = savKey;
                          break;
                        }
                      }
                    }
                    break;
                  }
                }

                updatedCount++;
              }
            }
          }
        }
      }

      // Apply energy shifts directly from .sav file using segment names/identifiers
      if (segmentsTab_) {
        SegmentsDataModel *segmentsModel = segmentsTab_->getSegmentsDataModel();
        if (segmentsModel) {
          QList<SegmentsDataData> segments = segmentsModel->getLines();

          // Temporarily disconnect signals to prevent double updates
          segmentsModel->blockSignals(true);

          // Match segments by name/identifier from .sav file
          for (auto it = shiftsFromSav.begin(); it != shiftsFromSav.end(); ++it) {
            QString segmentId = it.key();  // Segment identifier from .sav file
            QPair<double, double> shiftData = it.value();

            // Try to find matching segment by different methods (same logic as norms)
            int segmentIndex = -1;

            for (int i = 0; i < segments.size(); i++) {
              const SegmentsDataData &segment = segments[i];

              // Method 1: Match by segment number (convert segmentId to number and check 1-based index)
              bool segmentIdIsNumber = false;
              int segmentNumber = segmentId.toInt(&segmentIdIsNumber);
              if (segmentIdIsNumber && (segmentNumber == i + 1)) {
                segmentIndex = i;
                break;
              }

              // Method 2: Match by filename (if segment has a dataFile)
              if (!segment.dataFile.isEmpty()) {
                QString fileBaseName = QFileInfo(segment.dataFile).baseName();
                if (fileBaseName == segmentId) {
                  segmentIndex = i;
                  break;
                }

                // Also try matching with file extension
                QString fileName = QFileInfo(segment.dataFile).fileName();
                if (fileName == segmentId || fileName.startsWith(segmentId + ".")) {
                  segmentIndex = i;
                  break;
                }
              }

              // Method 3: Match by constructed segment name pattern
              QString constructedName = segment.dataFile.isEmpty() ? QString("%1").arg(i + 1) : QFileInfo(segment.dataFile).baseName();
              if (constructedName == segmentId) {
                segmentIndex = i;
                break;
              }
            }

            if (segmentIndex >= 0 && segmentIndex < segments.size()) {
              // Update the segments model directly
              // Check for NaN and replace with 0 if found
              double shiftValue = shiftData.first;
              if (std::isnan(shiftValue)) {
                shiftValue = 0.0;
              }

              QModelIndex shiftValueIndex = segmentsModel->index(segmentIndex, 14);  // energyShift column
              segmentsModel->setData(shiftValueIndex, shiftValue, Qt::EditRole);
              // energyShiftError (column 15) left alone for the same reason as
              // dataNormError above: the .sav carries a fit uncertainty, the
              // segment carries the assumed systematic that the shift penalty
              // is built around.
            }
          }

          // Re-enable signals
          segmentsModel->blockSignals(false);

          // Now manually update ALL FittingParameters from the updated model
          // Also update parameter names to match params.sav exactly
          for (int j = 0; j < fittingParameters.size(); j++) {
            FittingParameter &param = fittingParameters[j];
            if (param.category == "shift") {
              int segmentIndex = param.channelIndex;
              if (segmentIndex >= 0 && segmentIndex < segments.size()) {
                // Re-read the updated segments data
                QList<SegmentsDataData> updatedSegments = segmentsModel->getLines();
                param.value = updatedSegments[segmentIndex].energyShift;
                // Fit uncertainty taken from the .sav below, not from the
                // segment's systematic-error column.
                param.fitError = 0.0;

                // Find the corresponding params.sav name for this segment
                for (auto it = shiftsFromSav.begin(); it != shiftsFromSav.end(); ++it) {
                  QString segmentId = it.key();

                  // Check if this segmentId matches this segment
                  bool matches = false;
                  const SegmentsDataData &segment = updatedSegments[segmentIndex];

                  // Same matching logic as above
                  bool segmentIdIsNumber = false;
                  int segmentNumber = segmentId.toInt(&segmentIdIsNumber);
                  if (segmentIdIsNumber && (segmentNumber == segmentIndex + 1)) {
                    matches = true;
                  } else if (!segment.dataFile.isEmpty()) {
                    QString fileBaseName = QFileInfo(segment.dataFile).baseName();
                    QString fileName = QFileInfo(segment.dataFile).fileName();
                    if (fileBaseName == segmentId || fileName == segmentId || fileName.startsWith(segmentId + ".")) {
                      matches = true;
                    }
                  }

                  if (matches) {
                    param.fitError = it.value().second;
                    // Find the exact parameter name from params.sav
                    for (const QString &savKey : savParams.keys()) {
                      if (savKey.contains("shift", Qt::CaseInsensitive) && !savKey.contains("_rwa")) {
                        // Check if this savKey corresponds to this segmentId
                        QRegExp segmentNumRegex("segment_(\\d+)_energy_shift", Qt::CaseInsensitive);
                        QRegExp segmentShiftRegex("segment_(\\d+)_shift", Qt::CaseInsensitive);
                        QRegExp shiftNumRegex("shift_(\\d+)", Qt::CaseInsensitive);
                        QRegExp segmentNameRegex("([^_]+)_(?:energy_)?shift", Qt::CaseInsensitive);

                        QString extractedId;
                        if (segmentNumRegex.indexIn(savKey) != -1) {
                          extractedId = segmentNumRegex.cap(1);
                        } else if (segmentShiftRegex.indexIn(savKey) != -1) {
                          extractedId = segmentShiftRegex.cap(1);
                        } else if (shiftNumRegex.indexIn(savKey) != -1) {
                          extractedId = shiftNumRegex.cap(1);
                        } else if (segmentNameRegex.indexIn(savKey) != -1) {
                          extractedId = segmentNameRegex.cap(1);
                        }

                        if (extractedId == segmentId) {
                          // Update the parameter name to match params.sav exactly
                          param.name = savKey;
                          break;
                        }
                      }
                    }
                    break;
                  }
                }

                updatedCount++;
              }
            }
          }
        }
      }

      // Robust verification: Re-read params.sav directly and check assignments
      QString verificationReport;

      // Re-read the params.sav file to get exact parameter names from first column
      QFile verificationFile(filename);
      QMap<QString, QPair<double, double>> directSavParams;

      if (verificationFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QTextStream verificationIn(&verificationFile);

        while (!verificationIn.atEnd()) {
          QString line = verificationIn.readLine().trimmed();
          if (line.isEmpty()) continue;

          QStringList parts = line.split(QRegExp("\\s+"), Qt::SkipEmptyParts);
          if (parts.size() >= 3) {
            QString exactParamName = parts[0];  // Exact parameter name from first column
            double value = parts[1].toDouble();
            double error = parts[2].toDouble();
            directSavParams[exactParamName] = qMakePair(value, error);
          }
        }
        verificationFile.close();
      }

      // Check normalizations against current assignments
      if (segmentsTab_) {
        SegmentsDataModel *segmentsModel = segmentsTab_->getSegmentsDataModel();
        if (segmentsModel) {
          QList<SegmentsDataData> segments = segmentsModel->getLines();
          verificationReport += "=== ROBUST NORMALIZATION VERIFICATION ===\n";

          // Find all norm parameters in the .sav file
          for (const QString &paramName : directSavParams.keys()) {
            if (paramName.contains("norm", Qt::CaseInsensitive) && !paramName.contains("_rwa")) {
              QPair<double, double> expectedData = directSavParams[paramName];

              // Extract segment number using the same regex as assignment
              QRegExp normRegex("(?:segment_|norm_?)(\\d+)(?:_norm)?", Qt::CaseInsensitive);
              if (normRegex.indexIn(paramName) != -1) {
                int segmentNumber = normRegex.cap(1).toInt();
                int segmentIndex = segmentNumber - 1;

                if (segmentIndex >= 0 && segmentIndex < segments.size()) {
                  const SegmentsDataData &segment = segments[segmentIndex];
                  bool isMatch = (qAbs(expectedData.first - segment.dataNorm) < 1e-10);

                  verificationReport += QString("%1: %2 → Segment[%3] Expected=%4, Found=%5 %6\n")
                                            .arg(paramName)
                                            .arg(segmentNumber)
                                            .arg(segmentIndex)
                                            .arg(expectedData.first, 0, 'e', 6)
                                            .arg(segment.dataNorm, 0, 'e', 6)
                                            .arg(isMatch ? "✓" : "✗ MISMATCH!");

                  // If there's a mismatch, try to find where this value actually went
                  if (!isMatch) {
                    verificationReport += QString("  → Searching for value %1 in other segments...\n")
                                              .arg(expectedData.first, 0, 'e', 6);
                    for (int i = 0; i < segments.size(); i++) {
                      if (qAbs(expectedData.first - segments[i].dataNorm) < 1e-10) {
                        verificationReport += QString("  → FOUND at Segment[%1] (should be [%2])!\n")
                                                  .arg(i)
                                                  .arg(segmentIndex);
                      }
                    }
                  }
                } else {
                  verificationReport += QString("%1: Invalid segment index %2 (out of range)\n")
                                            .arg(paramName)
                                            .arg(segmentIndex);
                }
              } else {
                verificationReport += QString("%1: Could not extract segment number from name\n")
                                          .arg(paramName);
              }
            }
          }
          verificationReport += "\n";

          // Also check energy shifts
          verificationReport += "=== ROBUST ENERGY SHIFT VERIFICATION ===\n";
          for (const QString &paramName : directSavParams.keys()) {
            if (paramName.contains("shift", Qt::CaseInsensitive) && !paramName.contains("_rwa")) {
              QPair<double, double> expectedData = directSavParams[paramName];

              // Extract segment number using the same regex as assignment
              QRegExp shiftRegex("(?:segment_|shift_?)(\\d+)(?:_(?:energy_)?shift)?", Qt::CaseInsensitive);
              if (shiftRegex.indexIn(paramName) != -1) {
                int segmentNumber = shiftRegex.cap(1).toInt();
                int segmentIndex = segmentNumber - 1;

                if (segmentIndex >= 0 && segmentIndex < segments.size()) {
                  const SegmentsDataData &segment = segments[segmentIndex];
                  bool isMatch = (qAbs(expectedData.first - segment.energyShift) < 1e-10);

                  verificationReport += QString("%1: %2 → Segment[%3] Expected=%4, Found=%5 %6\n")
                                            .arg(paramName)
                                            .arg(segmentNumber)
                                            .arg(segmentIndex)
                                            .arg(expectedData.first, 0, 'e', 6)
                                            .arg(segment.energyShift, 0, 'e', 6)
                                            .arg(isMatch ? "✓" : "✗ MISMATCH!");

                  // If there's a mismatch, try to find where this value actually went
                  if (!isMatch) {
                    verificationReport += QString("  → Searching for value %1 in other segments...\n")
                                              .arg(expectedData.first, 0, 'e', 6);
                    for (int i = 0; i < segments.size(); i++) {
                      if (qAbs(expectedData.first - segments[i].energyShift) < 1e-10) {
                        verificationReport += QString("  → FOUND at Segment[%1] (should be [%2])!\n")
                                                  .arg(i)
                                                  .arg(segmentIndex);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }

      // Write the level parameters -- including the fixed and non-varying ones --
      // back to their models, so the RWA-to-physical values just read from the .sav
      // are what the .azr carries when it is saved.
      //
      // Norm and shift parameters are deliberately excluded.  Their values have
      // already gone into the segments model above, value columns only.  Sending
      // them through updateParameterInOtherTabs as well would repeat that write and
      // additionally push param.error and param.useAsNuisance into the uncertainty
      // and "vary" columns -- fit freedom and assumed systematics that the user set
      // and that a .sav does not describe at all.  Loading fitted *values* must not
      // silently redefine which parameters are free or how tightly they are
      // constrained.
      for (const FittingParameter &param : fittingParameters) {
        if (param.category == "level") updateParameterInOtherTabs(param.name, param);
      }

      // Refresh the parameter tables with ALL parameters (including fixed and non-varying)
      updateParameterTables();

      // Show verification results
      QMessageBox msgBox;
      msgBox.setWindowTitle("Load Settings - Verification Report");
      msgBox.setText(QString("Updated %1 fitting parameters from: %2%3")
                         .arg(updatedCount)
                         .arg(filename)
                         .arg(backupNote));
      msgBox.setDetailedText(verificationReport);
      msgBox.exec();

      // Re-read everything from the (now updated) GUI models, i.e. the same
      // "Refresh from Current" behaviour, so the tables end up fully in sync.
      populateFromCurrentGUIState();
    } else {
      QMessageBox::warning(this, "Load Error",
                           "Could not load parameter file.");
    }
  }
}

QString FittingTab::findMatchingParameterKey(const FittingParameter &param, const QStringList &savKeys) {
  // Try to match current parameter with .sav file parameter names

  // First try direct match
  if (savKeys.contains(param.name)) {
    return param.name;
  }

  if (param.category == "norm") {
    // For normalization parameters: look for patterns like "segment_1_norm", "norm_1", etc.
    int segmentIndex = param.channelIndex + 1;  // Convert 0-based to 1-based

    QStringList possibleNames;
    possibleNames << QString("segment_%1_norm").arg(segmentIndex);
    possibleNames << QString("norm_%1").arg(segmentIndex);
    possibleNames << QString("norm%1").arg(segmentIndex);
    possibleNames << QString("segment%1_norm").arg(segmentIndex);

    for (const QString &name : possibleNames) {
      if (savKeys.contains(name)) {
        return name;
      }
    }

  } else if (param.category == "shift") {
    // For energy shift parameters: look for patterns like "segment_1_energy_shift", "shift_1", etc.
    int segmentIndex = param.channelIndex + 1;  // Convert 0-based to 1-based

    QStringList possibleNames;
    possibleNames << QString("segment_%1_energy_shift").arg(segmentIndex);  // Main pattern from .sav files
    possibleNames << QString("segment_%1_shift").arg(segmentIndex);
    possibleNames << QString("shift_%1").arg(segmentIndex);
    possibleNames << QString("shift%1").arg(segmentIndex);
    possibleNames << QString("segment%1_energy_shift").arg(segmentIndex);
    possibleNames << QString("segment%1_shift").arg(segmentIndex);

    for (const QString &name : possibleNames) {
      if (savKeys.contains(name)) {
        return name;
      }
    }

  } else if (param.category == "level") {
    // For level parameters: match with .sav file patterns like "energy_1", "width_1_2"
    // Current GUI names are like "Level 1 Energy (MeV)" and "Level 1 Channel 2 Width (eV)"

    if (param.name.contains("Energy") && param.channelIndex == -1) {
      // Energy parameter: "Level N Energy (MeV)" -> "energy_N", N in AZURE2's numbering
      int levelIndex = engineLevelOrder(levelsTab_ ? levelsTab_->writeOrder() : QList<int>()).indexOf(param.levelIndex) + 1;
      if (levelIndex < 1) return QString();

      QStringList possibleNames;
      possibleNames << QString("energy_%1").arg(levelIndex);
      possibleNames << QString("energy%1").arg(levelIndex);
      possibleNames << QString("level_%1_energy").arg(levelIndex);
      possibleNames << QString("level%1_energy").arg(levelIndex);

      for (const QString &name : possibleNames) {
        if (savKeys.contains(name)) {
          return name;
        }
      }

    } else if (param.name.contains("Width") && param.channelIndex >= 0) {
      // Width parameter: "Level N Channel M Width (eV)" -> "width_N_M"
      // Need to count channels per level, not global channel index
      int levelIndex = engineLevelOrder(levelsTab_ ? levelsTab_->writeOrder() : QList<int>()).indexOf(param.levelIndex) + 1;
      if (levelIndex < 1) return QString();

      // Count which channel this is within this specific level
      int channelWithinLevel = 1;  // Start counting from 1 for .sav file format

      if (levelsTab_ && segmentsTab_) {
        ChannelsModel *channelsModel = levelsTab_->getChannelsModel();
        if (channelsModel) {
          QList<ChannelsData> channels = channelsModel->getChannels();
          for (int i = 0; i < channels.size(); i++) {
            const ChannelsData &channel = channels[i];
            if (channel.levelIndex == param.levelIndex) {
              if (i == param.channelIndex) {
                // Found our channel - channelWithinLevel is correct
                break;
              }
              channelWithinLevel++;  // Count channels for this level
            }
          }
        }
      }

      QStringList possibleNames;
      possibleNames << QString("width_%1_%2").arg(levelIndex).arg(channelWithinLevel);
      possibleNames << QString("width%1_%2").arg(levelIndex).arg(channelWithinLevel);
      possibleNames << QString("level_%1_channel_%2_width").arg(levelIndex).arg(channelWithinLevel);
      possibleNames << QString("level%1_channel%2_width").arg(levelIndex).arg(channelWithinLevel);

      for (const QString &name : possibleNames) {
        if (savKeys.contains(name)) {
          return name;
        }
      }
    }
  }

  // No match found
  return QString();
}

// Write one column of one segment, for the two fields that only a deliberate
// user edit may change: the assumed systematic uncertainty and the vary flag.
// Kept separate from updateParameterInOtherTabs so that propagating a value can
// never carry either of them along with it.
void FittingTab::writeSegmentColumn(const FittingParameter &param, int column, const QVariant &value) {
  if (param.category != "norm" && param.category != "shift") return;
  if (!segmentsTab_) return;
  SegmentsDataModel *segmentsModel = segmentsTab_->getSegmentsDataModel();
  if (!segmentsModel) return;
  const int segmentIndex = param.channelIndex;  // segment index is stored in channelIndex
  if (segmentIndex < 0 || segmentIndex >= segmentsModel->getLines().size()) return;
  segmentsModel->setData(segmentsModel->index(segmentIndex, column), value, Qt::EditRole);
}

void FittingTab::updateParameterInOtherTabs(const QString &paramName, const FittingParameter &param) {
  if (!levelsTab_ || !segmentsTab_) return;

  if (param.category == "level") {
    // Update level parameters in LevelsModel - same pattern as LevelsTab::editLevel()
    LevelsModel *levelsModel = levelsTab_->getLevelsModel();
    if (levelsModel) {
      // The model row and channel this parameter was built from.  Its name
      // carries AZURE2's numbering, which is not the model's.
      int levelIndex = param.levelIndex;
      int channelIndex = param.channelIndex;

      if (levelIndex >= 0) {
        QList<LevelsData> levels = levelsModel->getLevels();
        if (levelIndex < levels.size()) {
          if (channelIndex == -1) {
            // Energy parameter - update column 4 (energy column)
            QModelIndex index = levelsModel->index(levelIndex, 4);
            levelsModel->setData(index, param.value, Qt::EditRole);
          } else {
            // Width parameter - update channels model
            ChannelsModel *channelsModel = levelsTab_->getChannelsModel();
            if (channelsModel) {
              QList<ChannelsData> channels = channelsModel->getChannels();
              // Find the channel with matching level and channel indices
              for (int i = 0; i < channels.size(); i++) {
                if (channels[i].levelIndex == levelIndex && i == channelIndex) {
                  // Update reducedWidth column (column 6)
                  QModelIndex index = channelsModel->index(i, 6);
                  channelsModel->setData(index, param.value, Qt::EditRole);
                  break;
                }
              }
            }
          }
        }
      }
    }
  } else if (param.category == "norm" || param.category == "shift") {
    // Update normalization or shift parameters in SegmentsDataModel
    SegmentsDataModel *segmentsModel = segmentsTab_->getSegmentsDataModel();
    if (segmentsModel) {
      QList<SegmentsDataData> segments = segmentsModel->getLines();

      // Use channelIndex which now stores the segment index
      int segmentIndex = param.channelIndex;

      // Only the value.  The segment's uncertainty is the experimental
      // systematic the user entered to constrain the normalization, and its
      // vary flag is their choice of what the fit may move; neither is
      // implied by a new value, so neither is written here.  Both have a
      // single deliberate editing path -- the Segments tab dialog, and the
      // Error / Nuisance cells of this tab, which write them directly (see
      // parameterItemChanged) -- and nothing else may touch them.
      if (segmentIndex >= 0 && segmentIndex < segments.size()) {
        if (param.category == "norm") {
          // dataNorm column (column 9)
          segmentsModel->setData(segmentsModel->index(segmentIndex, 9), param.value, Qt::EditRole);
        } else if (param.category == "shift") {
          // energyShift column (column 14)
          segmentsModel->setData(segmentsModel->index(segmentIndex, 14), param.value, Qt::EditRole);
        }
      }
    }
  }
}

void FittingTab::updateParameterTableValue(const QString &paramName, double value) {
  // Determine which table the parameter belongs to
  QTableWidget *targetTable = nullptr;
  if (paramName.contains("Normalization")) {
    targetTable = normParamsTable;
  } else if (paramName.contains("Energy Shift")) {
    targetTable = shiftParamsTable;
  } else {
    targetTable = levelParamsTable;
  }

  if (!targetTable) return;

  // Find the row with this parameter name
  for (int row = 0; row < targetTable->rowCount(); row++) {
    QTableWidgetItem *nameItem = targetTable->item(row, 0);
    if (nameItem && nameItem->text() == paramName) {
      // Update the value column (column 1)
      QTableWidgetItem *valueItem = targetTable->item(row, 1);
      if (valueItem) {
        // Temporarily disconnect signals to avoid recursion
        targetTable->blockSignals(true);
        valueItem->setText(QString::number(value, 'g', 6));
        targetTable->blockSignals(false);
      }
      break;
    }
  }
}

void FittingTab::showInfo(int which, QString title) {
  if (!infoDialog[which]) {
    // QString construct info dialog
    infoDialog[which] = new InfoDialog("", this, title);
    if (title.isEmpty()) {
      switch (which) {
        case 0: infoDialog[which]->setWindowTitle("Level Parameters Help"); break;
        case 1: infoDialog[which]->setWindowTitle("Normalization Parameters Help"); break;
        case 2: infoDialog[which]->setWindowTitle("Energy Shift Parameters Help"); break;
      }
    } else {
      infoDialog[which]->setWindowTitle(title);
    }
    // infoDialog[which]->setInfoText(infoText[which]);
  }
  infoDialog[which]->show();
  infoDialog[which]->raise();
}

bool FittingTab::writeParameterSettings(QTextStream &outStream) {
  // Write current parameter settings to AZURE2 file
  outStream << "# Fitting parameter settings (only non-fixed parameters shown)\n";
  outStream << "# Format: name value lower_limit upper_limit error fit_error nuisance category minuit_index\n";

  for (const FittingParameter &param : fittingParameters) {
    outStream << param.name << " "
              << roundTripNumber(param.value) << " "
              << roundTripNumber(param.lowerLimit) << " "
              << roundTripNumber(param.upperLimit) << " "
              << roundTripNumber(param.error) << " "
              << roundTripNumber(param.fitError) << " "
              << (param.useAsNuisance ? 1 : 0) << " "
              << param.category << " "
              << param.minuitIndex << "\n";
  }
  const QStringList priorRows = priorCentreRows();
  if (!priorRows.isEmpty()) {
    outStream << "# Prior centres: segment_N_norm|segment_N_energy_shift prior_centre value\n";
    for (const QString &row : priorRows) outStream << row << "\n";
  }
  return true;
}

QStringList FittingTab::priorCentreRows() const {
  QRegExp entry("^segment_(\\d+)_(norm|energy_shift)$");
  QMap<QPair<int, int>, QString> ordered;  // (segment, 0 norm / 1 shift) -> row
  for (QMap<QString, double>::const_iterator c = priorCentres_.constBegin(); c != priorCentres_.constEnd(); ++c) {
    if (entry.indexIn(c.key()) == -1) continue;
    ordered[qMakePair(entry.cap(1).toInt(), entry.cap(2) == "norm" ? 0 : 1)] =
        c.key() + " prior_centre " + roundTripNumber(c.value());
  }
  return ordered.values();
}

bool FittingTab::readParameterSettings(QTextStream &inStream) {
  savedParameterSettings.clear();

  // Settings for normalizations and energy shifts are keyed by segment number,
  // so an entry that outlives the segment it described will silently attach
  // itself to whatever segment later takes that number -- handing a new data
  // set someone else's nuisance prior and freeing a normalization the user
  // never freed. An entry numbered beyond the segments that actually exist is
  // provably stale, so drop it here rather than carry it forward. Segments are
  // read before this point, so the count is already known.
  int numSegments = 0;
  if (segmentsTab_ && segmentsTab_->getSegmentsDataModel())
    numSegments = segmentsTab_->getSegmentsDataModel()->getLines().size();
  QRegExp segmentEntry("^segment_(\\d+)_(norm|energy_shift)$");

  // Level entries are named the way AZURE2 numbers this file: "Level N" is the
  // N-th level it builds from <levels> (read just before this section, into the
  // model in file order), "Channel k" the k-th channel of that level.  Each one
  // is tied to the model's level and channel here; the names written back on
  // save follow the order the levels are saved in, which need not be this one.
  //
  // Files written by earlier versions of the GUI numbered differently: "Level N"
  // was the N-th row of the levels model and "Channel k" the k-th row of the
  // channels model, whatever level it belonged to.  Such a file is recognised
  // by the values it records, which are those of the parameters it meant.
  QList<LevelsData> levels;
  QList<ChannelsData> channels;
  if (levelsTab_) {
    if (LevelsModel *m = levelsTab_->getLevelsModel()) levels = m->getLevels();
    if (ChannelsModel *m = levelsTab_->getChannelsModel()) channels = m->getChannels();
  }
  QList<int> fileOrder;
  for (int la = 0; la < levels.size(); la++) fileOrder.append(la);
  const QList<int> engineOrder = engineLevelOrder(fileOrder);
  // False when the entry names a level or channel the file does not have.
  auto locate = [&](FittingParameter &param, int levelNumber, int channelNumber, bool legacy) {
    param.levelIndex = -1;
    param.channelIndex = -1;
    if (legacy) {
      if (levelNumber < 1 || levelNumber > levels.size()) return false;
      param.levelIndex = levelNumber - 1;
      if (channelNumber == 0) return true;
      if (channelNumber < 1 || channelNumber > channels.size() ||
          channels.at(channelNumber - 1).levelIndex != param.levelIndex)
        return false;
      param.channelIndex = channelNumber - 1;
      return true;
    }
    if (levelNumber < 1 || levelNumber > engineOrder.size()) return false;
    param.levelIndex = engineOrder.at(levelNumber - 1);
    if (channelNumber == 0) return true;
    const QList<int> levelChannels = channelsOfLevel(channels, param.levelIndex);
    if (channelNumber < 1 || channelNumber > levelChannels.size()) return false;
    param.channelIndex = levelChannels.at(channelNumber - 1);
    return true;
  };
  // Does the recorded value belong to the parameter the entry was located at?
  auto recordedValueFits = [&](const FittingParameter &param) {
    double current = param.channelIndex < 0 ? levels.at(param.levelIndex).energy
                                            : channels.at(param.channelIndex).reducedWidth;
    return std::fabs(param.value - current) <= 1e-5 * std::max(std::fabs(current), 1e-30);
  };
  struct LevelEntry {
    FittingParameter param;
    int levelNumber;
    int channelNumber;
  };
  QList<LevelEntry> levelEntries;

  priorCentres_.clear();

  QString line;
  while (!inStream.atEnd()) {
    line = inStream.readLine().trimmed();
    if (line.startsWith("<")) break;                       // Next section started
    if (line.isEmpty() || line.startsWith("#")) continue;  // Skip empty lines and comments

    QStringList parts = line.split(" ", Qt::SkipEmptyParts);
    if (parts.size() == 3 && parts[1] == "prior_centre") {
      // An explicit prior centre (EData::ReadPriorCentres); a stale or
      // malformed one is dropped, as the entries below.
      bool ok;
      const double centre = parts[2].toDouble(&ok);
      if (!ok || segmentEntry.indexIn(parts[0]) == -1) continue;
      if (numSegments > 0 && segmentEntry.cap(1).toInt() > numSegments) continue;
      priorCentres_[parts[0]] = centre;
    } else if (parts.size() == 9) {
      FittingParameter param;
      param.name = parts[0];
      param.value = parts[1].toDouble();  // This will be overridden by current model values
      param.lowerLimit = parts[2].toDouble();
      param.upperLimit = parts[3].toDouble();
      param.error = parts[4].toDouble();
      param.fitError = parts[5].toDouble();
      param.useAsNuisance = (parts[6].toInt() == 1);
      param.category = parts[7];
      param.minuitIndex = parts[8].toInt();
      param.levelIndex = -1;
      param.channelIndex = -1;
      if (numSegments > 0 && segmentEntry.indexIn(param.name) != -1 &&
          segmentEntry.cap(1).toInt() > numSegments) continue;
      savedParameterSettings.append(param);
    } else if (parts.size() == 12) {
      // Backward compatibility: old format without fitError
      FittingParameter param;
      param.name = parts[0] + " " + parts[1] + " " + parts[2] + " " + parts[3];
      param.value = parts[4].toDouble();
      param.lowerLimit = parts[5].toDouble();
      param.upperLimit = parts[6].toDouble();
      param.error = parts[7].toDouble();
      param.fitError = parts[8].toDouble();  // Default for old format
      param.useAsNuisance = (parts[9].toInt() == 1);
      param.category = parts[10];
      param.minuitIndex = parts[11].toInt();

      levelEntries.append({param, parts[1].toInt(), 0});
    } else if (parts.size() == 14) {
      // Backward compatibility: old format without fitError
      FittingParameter param;
      param.name = parts[0] + " " + parts[1] + " " + parts[2] + " " + parts[3] + " " + parts[4] + " " + parts[5];
      param.value = parts[6].toDouble();
      param.lowerLimit = parts[7].toDouble();
      param.upperLimit = parts[8].toDouble();
      param.error = parts[9].toDouble();
      param.fitError = parts[10].toDouble();  // Default for old format
      param.useAsNuisance = (parts[11].toInt() == 1);
      param.category = parts[12];
      param.minuitIndex = parts[13].toInt();

      levelEntries.append({param, parts[1].toInt(), parts[3].toInt()});
    }
  }

  // Read the entries the way that accounts for more of the values they record,
  // then the way that finds more of them a parameter; on a tie, the way AZURE2
  // itself reads them.
  int fitsEngine = 0, fitsLegacy = 0, foundEngine = 0, foundLegacy = 0;
  for (LevelEntry &entry : levelEntries) {
    if (locate(entry.param, entry.levelNumber, entry.channelNumber, false)) {
      foundEngine++;
      if (recordedValueFits(entry.param)) fitsEngine++;
    }
    if (locate(entry.param, entry.levelNumber, entry.channelNumber, true)) {
      foundLegacy++;
      if (recordedValueFits(entry.param)) fitsLegacy++;
    }
  }
  const bool legacy = fitsLegacy > fitsEngine || (fitsLegacy == fitsEngine && foundLegacy > foundEngine);
  for (LevelEntry &entry : levelEntries)
    if (locate(entry.param, entry.levelNumber, entry.channelNumber, legacy))
      savedParameterSettings.append(entry.param);

  return true;
}

void FittingTab::onSegmentNormalizationChanged(int segmentIndex, double value) {
  // Find normalization parameter by segment index stored in channelIndex
  for (int i = 0; i < fittingParameters.size(); i++) {
    if (fittingParameters[i].category == "norm" && fittingParameters[i].channelIndex == segmentIndex) {
      fittingParameters[i].value = value;

      // Update the corresponding table cell
      updateParameterTableValue(fittingParameters[i].name, value);
      break;
    }
  }
}

void FittingTab::onSegmentEnergyShiftChanged(int segmentIndex, double value) {
  // Find energy shift parameter by segment index stored in channelIndex
  for (int i = 0; i < fittingParameters.size(); i++) {
    if (fittingParameters[i].category == "shift" && fittingParameters[i].channelIndex == segmentIndex) {
      fittingParameters[i].value = value;

      // Update the corresponding table cell
      updateParameterTableValue(fittingParameters[i].name, value);
      break;
    }
  }
}

void FittingTab::updateParameterTableError(const QString &paramName, double error) {
  // Determine which table the parameter belongs to
  QTableWidget *targetTable = nullptr;
  if (paramName.contains("Normalization")) {
    targetTable = normParamsTable;
  } else if (paramName.contains("Energy Shift")) {
    targetTable = shiftParamsTable;
  } else {
    targetTable = levelParamsTable;
  }

  if (!targetTable) return;

  // Find the row with this parameter name
  for (int row = 0; row < targetTable->rowCount(); row++) {
    QTableWidgetItem *nameItem = targetTable->item(row, 0);
    if (nameItem && nameItem->text() == paramName) {
      // Update the error column (column 4)
      QTableWidgetItem *errorItem = targetTable->item(row, 4);
      if (errorItem) {
        // Temporarily disconnect signals to avoid recursion
        targetTable->blockSignals(true);
        errorItem->setText(QString::number(error, 'g', 6));
        targetTable->blockSignals(false);
      }
      break;
    }
  }
}

void FittingTab::updateParameterTableCheckbox(const QString &paramName, bool checked) {
  // Route to the table by parameter name.  Norm/shift parameters are named
  // "segment_N_norm" / "segment_N_energy_shift"; check "shift" before "norm"
  // so the energy-shift name is not mistaken for a normalization.
  QTableWidget *targetTable = nullptr;
  if (paramName.contains("shift", Qt::CaseInsensitive)) {
    targetTable = shiftParamsTable;
  } else if (paramName.contains("norm", Qt::CaseInsensitive)) {
    targetTable = normParamsTable;
  } else {
    targetTable = levelParamsTable;
  }

  if (!targetTable) return;

  // Find the row with this parameter name
  for (int row = 0; row < targetTable->rowCount(); row++) {
    QTableWidgetItem *nameItem = targetTable->item(row, 0);
    if (nameItem && nameItem->text() == paramName) {
      // Update the checkbox column (column 6)
      QTableWidgetItem *checkboxItem = targetTable->item(row, 6);
      if (checkboxItem) {
        // Temporarily disconnect signals to avoid recursion
        targetTable->blockSignals(true);
        checkboxItem->setCheckState(checked ? Qt::Checked : Qt::Unchecked);
        targetTable->blockSignals(false);
      }
      break;
    }
  }
}

void FittingTab::onSegmentNormalizationErrorChanged(int segmentIndex, double error) {
  // Find normalization parameter by segment index stored in channelIndex
  for (int i = 0; i < fittingParameters.size(); i++) {
    if (fittingParameters[i].category == "norm" && fittingParameters[i].channelIndex == segmentIndex) {
      fittingParameters[i].error = error;

      // Update the corresponding table cell (error column is column 4)
      updateParameterTableError(fittingParameters[i].name, error);
      break;
    }
  }
}

void FittingTab::onSegmentEnergyShiftErrorChanged(int segmentIndex, double error) {
  // Find energy shift parameter by segment index stored in channelIndex
  for (int i = 0; i < fittingParameters.size(); i++) {
    if (fittingParameters[i].category == "shift" && fittingParameters[i].channelIndex == segmentIndex) {
      fittingParameters[i].error = error;

      // Update the corresponding table cell (error column is column 4)
      updateParameterTableError(fittingParameters[i].name, error);
      break;
    }
  }
}

void FittingTab::onSegmentNormalizationVaryChanged(int segmentIndex, bool vary) {
  // Find normalization parameter by segment index stored in channelIndex
  for (int i = 0; i < fittingParameters.size(); i++) {
    if (fittingParameters[i].category == "norm" && fittingParameters[i].channelIndex == segmentIndex) {
      fittingParameters[i].useAsNuisance = vary;

      // Update the checkbox in the table
      updateParameterTableCheckbox(fittingParameters[i].name, vary);
      break;
    }
  }
  // A normalization enters or leaves AZURE2's Minuit vector with this flag, so
  // every recorded index after it moves.
  assignMinuitIndices();
}

void FittingTab::onSegmentEnergyShiftVaryChanged(int segmentIndex, bool vary) {
  // Find energy shift parameter by segment index stored in channelIndex
  for (int i = 0; i < fittingParameters.size(); i++) {
    if (fittingParameters[i].category == "shift" && fittingParameters[i].channelIndex == segmentIndex) {
      fittingParameters[i].useAsNuisance = vary;

      // Update the checkbox in the table
      updateParameterTableCheckbox(fittingParameters[i].name, vary);
      break;
    }
  }
}

void FittingTab::clearLimits() {
  // First, populate the fittingParameters from current GUI state
  populateFromCurrentGUIState();

  // Set all limits to 0 for all parameters
  int clearedCount = 0;
  for (int i = 0; i < fittingParameters.size(); i++) {
    if (fittingParameters[i].lowerLimit != 0.0 || fittingParameters[i].upperLimit != 0.0) {
      fittingParameters[i].lowerLimit = 0.0;
      fittingParameters[i].upperLimit = 0.0;
      clearedCount++;
    }
  }

  // Refresh the parameter tables to show the cleared limits
  updateParameterTables();

  if (clearedCount > 0) {
    QMessageBox::information(this, "Success",
                             QString("Cleared limits for %1 parameters. All limits are now set to 0 (unlimited).").arg(clearedCount));
  } else {
    QMessageBox::information(this, "No Changes",
                             "All parameter limits were already set to 0.");
  }
}

// Static info text - this would be defined in InTabDocs.cpp in the real implementation
const std::vector<QString> FittingTab::infoText = {
    QString("Level parameters control the R-matrix level energies and widths. "
            "Set limits to constrain parameter values during fitting. "
            "Enable 'Use as Nuisance' to include parameter uncertainty in chi-squared."),
    QString("Normalization parameters adjust the overall scale of data segments. "
            "These are typically varied during fitting to account for experimental uncertainties."),
    QString("Energy shift parameters correct for energy calibration offsets in data segments. "
            "These parameters shift the energy scale of experimental data points.")};

void FittingTab::followSegments(const QVector<int> &newNumber) {
  QRegExp entry("^segment_(\\d+)_(norm|energy_shift)$");
  // The new name of a segment-keyed name, "" if its segment is gone; other names are kept.
  auto renamed = [&](const QString &name) {
    if (entry.indexIn(name) == -1) return name;
    const int k = entry.cap(1).toInt();
    if (k < 1 || k > newNumber.size()) return name;
    return newNumber[k - 1] ? QString("segment_%1_%2").arg(newNumber[k - 1]).arg(entry.cap(2)) : QString();
  };
  const QSet<QString> coherent = coherentParameterNames();  // the cbkg_* names that still exist
  // In the order populateFromCurrentGUIState builds: levels, norms and
  // shifts by segment, then the THM background.
  QList<FittingParameter> levels, norms, shifts, rest;
  for (FittingParameter p : fittingParameters) {
    if (p.category == "norm" || p.category == "shift") {
      p.name = renamed(p.name);
      if (p.name.isEmpty()) continue;
      if (p.channelIndex >= 0 && p.channelIndex < newNumber.size()) p.channelIndex = newNumber[p.channelIndex] - 1;
      (p.category == "norm" ? norms : shifts) << p;
    } else if (p.category == "cbkg") {
      if (coherent.contains(p.name)) rest << p;
    } else if (p.category == "level") {
      levels << p;
    } else {
      rest << p;
    }
  }
  // A line added at the end (newNumber covers the lines before it) gets its
  // parameters as populateFromCurrentGUIState gives them, with any saved settings.
  QList<SegmentsDataData> segments;
  if (segmentsTab_ && segmentsTab_->getSegmentsDataModel()) segments = segmentsTab_->getSegmentsDataModel()->getLines();
  for (int i = newNumber.size(); i < segments.size(); i++) {
    if (segments[i].isActive != 1) continue;
    for (int norm = 1; norm >= 0; norm--) {
      FittingParameter p = segmentParameter(i, segments[i], norm);
      for (const FittingParameter &saved : savedParameterSettings)
        if (saved.name == p.name) {
          p.lowerLimit = saved.lowerLimit;
          p.upperLimit = saved.upperLimit;
          p.error = saved.error;
          p.useAsNuisance = saved.useAsNuisance;
          break;
        }
      (norm ? norms : shifts) << p;
    }
  }
  auto bySegment = [](const FittingParameter &a, const FittingParameter &b) { return a.channelIndex < b.channelIndex; };
  std::stable_sort(norms.begin(), norms.end(), bySegment);
  std::stable_sort(shifts.begin(), shifts.end(), bySegment);
  fittingParameters = levels + norms + shifts + rest;

  QList<FittingParameter> saved;
  for (FittingParameter p : savedParameterSettings) {
    p.name = renamed(p.name);
    if (!p.name.isEmpty()) saved << p;
  }
  savedParameterSettings = saved;
  QMap<QString, double> centres;
  for (auto c = priorCentres_.constBegin(); c != priorCentres_.constEnd(); ++c)
    if (!renamed(c.key()).isEmpty()) centres[renamed(c.key())] = c.value();
  priorCentres_ = centres;

  assignMinuitIndices();
  updateParameterTables();
}
