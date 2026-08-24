#include "MCMCTab.h"
#include "InfoDialog.h"
#include "LevelsTab.h"
#include "SegmentsTab.h"
#include "FittingTab.h"
#include "LevelsModel.h"
#include "ChannelsModel.h"
#include "SegmentsDataModel.h"
#include "AZURESetup.h"

#include <QGroupBox>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QLabel>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QProgressBar>
#include <QTextEdit>
#include <QTabWidget>
#include <QTableWidget>
#include <QCheckBox>
#include <QComboBox>
#include <QHeaderView>
#include <QFileDialog>
#include <QTableWidgetItem>
#include <QSignalMapper>
#include <QMessageBox>
#include <QBrush>
#include <QColor>
#include <QApplication>
#include <QThread>
#include <QTimer>
#include <QRandomGenerator>
#include <QFileInfo>
#include <QDir>
#include <QStringList>
#include <QTextStream>
#include <QRegExp>
#include <QIODevice>
#include <QtCore>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <thread>
#include <chrono>

#ifdef USE_MCMC
#include "AZURECalcMCMC.h"
#include "ParameterLabel.h"
#include "Config.h"
#include "CNuc.h"
#include "EData.h"
#include "AZUREParams.h"
#endif

// Forward declarations and structure definitions from AZURE2.cpp
struct SegPairs {
  int firstPair;
  int secondPair;
};
extern bool readSegmentFile(const Config &configure, std::vector<SegPairs> &segPairs);
extern bool checkExternalCapture(Config &configure, const std::vector<SegPairs> &segPairs);

const std::vector<QString> MCMCTab::infoText = {
    "Parameter Configuration:\n\n"
    "Configure parameters for MCMC sampling. Each parameter can have:\n"
    "- Prior mean: Central value for Gaussian prior\n"
    "- Prior std: Standard deviation for Gaussian prior\n"
    "- Use Gaussian Prior: Whether to apply Gaussian prior constraint\n\n"
    "Parameters without Gaussian priors use uniform priors over the fitted range.",

    "Sampling Settings:\n\n"
    "- Walkers: Number of parallel MCMC chains (recommended: 2-10 times number of parameters)\n"
    "- Steps: Number of MCMC steps per walker\n"
    "- Burn-in: Number of initial steps to discard\n"
    "- Chain Spread: Percentage spread for initial walker positions around parameter values",

    "Progress Monitoring:\n\n"
    "Track the progress of MCMC sampling:\n"
    "- Progress bar shows completion percentage\n"
    "- Status shows current step and acceptance rate\n"
    "- Log shows detailed sampling information\n\n"
    "Use Stop button to halt sampling early if needed.",

    "Results Analysis:\n\n"
    "After sampling completes:\n"
    "- Statistics: Parameter means, std dev, credible intervals\n"
    "- Chains: Raw MCMC chain data\n"
    "- Convergence diagnostics and autocorrelation analysis\n\n"
    "Results are automatically saved to output directory."};

MCMCTab::MCMCTab(QWidget *parent) :
  QWidget(parent),
  levelsTab_(nullptr),
  segmentsTab_(nullptr),
  isRunning(false),
  stepOffset(0),
  mcmcCompletedSuccessfully(false) {
  // Create main layout
  QVBoxLayout *mainLayout = new QVBoxLayout(this);

  // Create tab widget for organization
  QTabWidget *tabWidget = new QTabWidget();

  // Parameters tab
  QWidget *parametersWidget = new QWidget();
  QVBoxLayout *parametersLayout = new QVBoxLayout(parametersWidget);

  setupParameterTable();
  parametersLayout->addWidget(parametersTable);

  QHBoxLayout *paramButtonLayout = new QHBoxLayout();
  loadButton = new QPushButton("Load Physical Parameters");
  loadSavButton = new QPushButton("Load RWA Parameters");
  paramButtonLayout->addWidget(loadButton);
  paramButtonLayout->addWidget(loadSavButton);
  paramButtonLayout->addStretch();
  parametersLayout->addLayout(paramButtonLayout);

  tabWidget->addTab(parametersWidget, "Parameters");

  // Sampling settings tab
  QWidget *samplingWidget = new QWidget();
  setupSamplingControls(samplingWidget);
  tabWidget->addTab(samplingWidget, "Sampling Settings");

  // Results tab
  QWidget *resultsWidget = new QWidget();
  QVBoxLayout *resultsLayout = new QVBoxLayout(resultsWidget);

  resultsTabWidget = new QTabWidget();

  // Statistics subtab
  statisticsTable = new QTableWidget();
  statisticsTable->setColumnCount(6);
  QStringList statsHeaders;
  statsHeaders << "Parameter" << "Mean" << "Std Dev" << "2.5%" << "50%" << "97.5%";
  statisticsTable->setHorizontalHeaderLabels(statsHeaders);
  statisticsTable->horizontalHeader()->setStretchLastSection(true);
  resultsTabWidget->addTab(statisticsTable, "Statistics");

  // Chains subtab
  chainsTextEdit = new QTextEdit();
  chainsTextEdit->setReadOnly(true);
  chainsTextEdit->setFont(QFont("Courier", 9));
  resultsTabWidget->addTab(chainsTextEdit, "Chain Data");

  // Add refresh button for results
  QHBoxLayout *resultsButtonLayout = new QHBoxLayout();
  QPushButton *refreshResultsButton = new QPushButton("Refresh Results from File");
  refreshResultsButton->setToolTip("Load and calculate statistics from existing samples.mcmc file");
  connect(refreshResultsButton, &QPushButton::clicked, this, &MCMCTab::refreshResultsFromFile);
  resultsButtonLayout->addWidget(refreshResultsButton);
  resultsButtonLayout->addStretch();

  resultsLayout->addWidget(resultsTabWidget);
  resultsLayout->addLayout(resultsButtonLayout);
  tabWidget->addTab(resultsWidget, "Results");

  mainLayout->addWidget(tabWidget);

  // Progress and control section
  setupProgressControls();
  mainLayout->addWidget(progressControlsGroup);

  // Connect signals
  connect(loadButton, SIGNAL(clicked()), this, SLOT(loadFromPhysical()));
  connect(loadSavButton, SIGNAL(clicked()), this, SLOT(loadFromReduced()));
  // Note: runButton and stopButton are connected by AZURESetup, not here

  // Setup info system
  mapper = new QSignalMapper(this);
  connect(mapper, SIGNAL(mapped(int)), this, SLOT(showInfo(int)));

  for (int i = 0; i < 4; i++) {
    infoDialog[i] = 0;
    infoButton[i] = new QPushButton("?");
    infoButton[i]->setMaximumSize(30, 30);
    mapper->setMapping(infoButton[i], i);
    connect(infoButton[i], SIGNAL(clicked()), mapper, SLOT(map()));
  }

  reset();
}

bool MCMCTab::isAutoPriorCategory(const QString &category) {
  // Normalizations and energy shifts carry an experimental error in the
  // segment definition, so AZURECalcMCMC::BuildAutoPriors() derives their
  // priors from the data and overrides whatever this table holds.
  return category == "norm" || category == "shift";
}

QString MCMCTab::categoryFromStoredName(const QString &name) {
  // Names persisted in the .azr file describe the parameter physically
  // ("normalization of segment 3 (data/x.dat)"). Older files carry the
  // previous forms, "Segment N Normalization" or "Parameter
  // (segment_N_norm)"; all of them are recognised here.
  if (name.contains("norm", Qt::CaseInsensitive)) return "norm";
  if (name.contains("shift", Qt::CaseInsensitive)) return "shift";
  return "loaded";
}

int MCMCTab::paramIndexForRow(int row) const {
  if (row < 0 || row >= tableRowToParam.size()) return -1;
  return tableRowToParam[row];
}

void MCMCTab::refreshParameterTable() {
  // Suppress itemChanged while rebuilding, or the handler writes the
  // half-populated row back into mcmcParameters.
  const bool wasBlocked = parametersTable->blockSignals(true);

  // Only the parameters whose priors the user is responsible for: the level
  // energies and widths. Normalizations and energy shifts are still sampled --
  // AZURECalcMCMC derives their priors from the errors quoted in the segment
  // definitions -- so there is nothing here for the user to fill in and a row
  // for them would just be read-only clutter.
  tableRowToParam.clear();
  for (int i = 0; i < mcmcParameters.size(); i++) {
    if (!mcmcParameters[i].autoPrior) tableRowToParam.append(i);
  }

  parametersTable->setRowCount(tableRowToParam.size());
  for (int row = 0; row < tableRowToParam.size(); row++) {
    const int i = tableRowToParam[row];
    const MCMCParameter &param = mcmcParameters[i];

    parametersTable->setItem(row, 0, new QTableWidgetItem(param.name));
    parametersTable->item(row, 0)->setFlags(Qt::ItemIsEnabled);

    parametersTable->setItem(row, 1, new QTableWidgetItem(QString::number(param.value, 'g', 6)));
    parametersTable->item(row, 1)->setFlags(Qt::ItemIsEnabled);

    QTableWidgetItem *meanItem = new QTableWidgetItem(QString::number(param.priorMean, 'g', 6));
    QTableWidgetItem *stdItem = new QTableWidgetItem(QString::number(param.priorStd, 'g', 6));

    QTableWidgetItem *checkItem = new QTableWidgetItem();
    checkItem->setCheckState(param.useGaussianPrior ? Qt::Checked : Qt::Unchecked);

    checkItem->setFlags(checkItem->flags() | Qt::ItemIsUserCheckable);

    parametersTable->setItem(row, 2, meanItem);
    parametersTable->setItem(row, 3, stdItem);
    parametersTable->setItem(row, 4, checkItem);

    parametersTable->setItem(row, 5, new QTableWidgetItem(param.category));
    parametersTable->item(row, 5)->setFlags(Qt::ItemIsEnabled);
  }

  parametersTable->blockSignals(wasBlocked);
}

void MCMCTab::setupParameterTable() {
  parametersTable = new QTableWidget();
  parametersTable->setColumnCount(6);

  QStringList headers;
  headers << "Parameter" << "Current Value" << "Prior Mean" << "Prior Std" << "Use Gaussian Prior" << "Category";
  parametersTable->setHorizontalHeaderLabels(headers);

  parametersTable->horizontalHeader()->setStretchLastSection(true);
  parametersTable->setSelectionBehavior(QAbstractItemView::SelectRows);
  parametersTable->setAlternatingRowColors(true);

  connect(parametersTable, SIGNAL(itemChanged(QTableWidgetItem *)),
          this, SLOT(parameterItemChanged(QTableWidgetItem *)));
}

void MCMCTab::setupSamplingControls(QWidget *samplingWidget) {
  QVBoxLayout *samplingLayout = new QVBoxLayout(samplingWidget);

  // Basic sampling parameters
  QGroupBox *basicGroup = new QGroupBox("Basic Settings");
  QGridLayout *basicLayout = new QGridLayout(basicGroup);

  basicLayout->addWidget(new QLabel("Number of Walkers:"), 0, 0);
  nWalkersSpinBox = new QSpinBox();
  nWalkersSpinBox->setRange(2, 1000);
  nWalkersSpinBox->setValue(50);
  basicLayout->addWidget(nWalkersSpinBox, 0, 1);
  if (infoButton[1]) basicLayout->addWidget(infoButton[1], 0, 2);

  basicLayout->addWidget(new QLabel("Number of Steps:"), 1, 0);
  nStepsSpinBox = new QSpinBox();
  nStepsSpinBox->setRange(100, 1000000);
  nStepsSpinBox->setValue(10000);
  basicLayout->addWidget(nStepsSpinBox, 1, 1);


  basicLayout->addWidget(new QLabel("Initial Chain Spread (%):"), 2, 0);
  chainSpreadSpinBox = new QDoubleSpinBox();
  chainSpreadSpinBox->setRange(0.001, 50.0);
  chainSpreadSpinBox->setValue(5.0);
  chainSpreadSpinBox->setSingleStep(0.5);
  chainSpreadSpinBox->setDecimals(3);
  chainSpreadSpinBox->setSuffix("%");
  chainSpreadSpinBox->setToolTip("Initial walker spread as a percentage of each parameter's value.\n"
                                 "Applies to widths; level energies use the absolute spread below.");
  basicLayout->addWidget(chainSpreadSpinBox, 2, 1);

  basicLayout->addWidget(new QLabel("Level Energy Spread:"), 3, 0);
  energySpreadSpinBox = new QDoubleSpinBox();
  energySpreadSpinBox->setRange(0.001, kMaxEnergySpreadKeV);
  energySpreadSpinBox->setValue(kMaxEnergySpreadKeV);
  energySpreadSpinBox->setSingleStep(0.1);
  energySpreadSpinBox->setDecimals(3);
  energySpreadSpinBox->setSuffix(" keV");
  energySpreadSpinBox->setToolTip(
      QString("Initial walker spread of the level energies, as an absolute\n"
              "energy rather than a percentage. Capped at %1 keV: a percentage\n"
              "of a several-MeV level energy scatters walkers across unrelated\n"
              "resonance structures and the chain does not converge.")
          .arg(kMaxEnergySpreadKeV));
  basicLayout->addWidget(energySpreadSpinBox, 3, 1);

  basicLayout->addWidget(new QLabel("Number of Threads:"), 4, 0);
  nThreadsSpinBox = new QSpinBox();
  nThreadsSpinBox->setRange(1, std::thread::hardware_concurrency());
  nThreadsSpinBox->setValue(std::min(4, (int)std::thread::hardware_concurrency()));  // Default to 4 or max available
  nThreadsSpinBox->setToolTip(QString("Number of parallel threads (1-%1 available)").arg(std::thread::hardware_concurrency()));
  basicLayout->addWidget(nThreadsSpinBox, 4, 1);

  samplingLayout->addWidget(basicGroup);

  // Note: Advanced settings (temperature, parallel sampling, proposal distribution)
  // were removed as per user request since they are not used with the numcmc library
  samplingLayout->addStretch();
}

void MCMCTab::setupProgressControls() {
  progressControlsGroup = new QGroupBox("MCMC Control");
  QVBoxLayout *progressLayout = new QVBoxLayout(progressControlsGroup);

  // Control buttons
  QHBoxLayout *buttonLayout = new QHBoxLayout();
  runButton = new QPushButton("Run MCMC Sampling");
  stopButton = new QPushButton("Stop Sampling");
  stopButton->setEnabled(false);

  buttonLayout->addWidget(runButton);
  buttonLayout->addWidget(stopButton);

  freshStartCheckBox = new QCheckBox("Fresh Start");
  freshStartCheckBox->setToolTip("When checked, ignores existing samples and starts from iteration 1");
  buttonLayout->addWidget(freshStartCheckBox);

  useReducedWidthsCheckBox = new QCheckBox("Use Reduced Widths");
  useReducedWidthsCheckBox->setToolTip("When checked, uses reduced width amplitudes (RWA) for fitting instead of physical parameters");
  buttonLayout->addWidget(useReducedWidthsCheckBox);

  if (infoButton[2]) buttonLayout->addWidget(infoButton[2]);
  buttonLayout->addStretch();

  progressLayout->addLayout(buttonLayout);

  // Progress bar and status
  progressBar = new QProgressBar();
  progressBar->setVisible(false);
  progressLayout->addWidget(progressBar);

  statusLabel = new QLabel("Ready to run MCMC sampling");
  progressLayout->addWidget(statusLabel);

  // Current iteration display
  QHBoxLayout *iterLayout = new QHBoxLayout();
  iterLayout->addWidget(new QLabel("Current Iteration:"));
  currentIterationLabel = new QLabel("0");
  currentIterationLabel->setStyleSheet("font-weight: bold; color: white;");
  iterLayout->addWidget(currentIterationLabel);
  iterLayout->addStretch();
  progressLayout->addLayout(iterLayout);

  // Log probability display
  QHBoxLayout *probLayout = new QHBoxLayout();
  probLayout->addWidget(new QLabel("Current Log Probability:"));
  logProbLabel = new QLabel("N/A");
  logProbLabel->setStyleSheet("font-weight: bold; color: white;");
  probLayout->addWidget(logProbLabel);

  probLayout->addWidget(new QLabel("Log Likelihood:"));
  logLikelihoodLabel = new QLabel("N/A");
  logLikelihoodLabel->setStyleSheet("font-weight: bold; color: white;");
  probLayout->addWidget(logLikelihoodLabel);

  probLayout->addWidget(new QLabel("Log Prior:"));
  logPriorLabel = new QLabel("N/A");
  logPriorLabel->setStyleSheet("font-weight: bold; color: white;");
  probLayout->addWidget(logPriorLabel);

  probLayout->addStretch();
  progressLayout->addLayout(probLayout);

  // Log area
  QGroupBox *logGroup = new QGroupBox("Sampling Log");
  QVBoxLayout *logLayout = new QVBoxLayout(logGroup);
  logTextEdit = new QTextEdit();
  logTextEdit->setMaximumHeight(150);
  logTextEdit->setReadOnly(true);
  logTextEdit->setFont(QFont("Courier", 8));
  logLayout->addWidget(logTextEdit);
  if (infoButton[3]) logLayout->addWidget(infoButton[3]);

  progressLayout->addWidget(logGroup);
}

void MCMCTab::reset() {
  mcmcParameters.clear();
  parametersTable->setRowCount(0);
  statisticsTable->setRowCount(0);
  chainsTextEdit->clear();
  logTextEdit->clear();
  stepOffset = 0;                     // Reset step offset
  mcmcCompletedSuccessfully = false;  // Reset completion status

  // Reset controls to defaults
  nWalkersSpinBox->setValue(50);
  nStepsSpinBox->setValue(10000);
  chainSpreadSpinBox->setValue(5.0);
  energySpreadSpinBox->setValue(kMaxEnergySpreadKeV);
  nThreadsSpinBox->setValue(std::min(4, (int)std::thread::hardware_concurrency()));

  statusLabel->setText("Ready to run MCMC sampling");
  progressBar->setVisible(false);
  progressBar->setValue(0);

  // Reset log probability displays
  currentIterationLabel->setText("0");
  logProbLabel->setText("N/A");
  logLikelihoodLabel->setText("N/A");
  logPriorLabel->setText("N/A");

  isRunning = false;
  runButton->setEnabled(true);
  stopButton->setEnabled(false);
}

void MCMCTab::setTabReferences(LevelsTab *levelsTab, SegmentsTab *segmentsTab) {
  levelsTab_ = levelsTab;
  segmentsTab_ = segmentsTab;
}


void MCMCTab::mcmcFinished() {
  // UI cleanup is now handled by AZURESetup::DeleteMCMCThread
  // Reset local state
  isRunning = false;
  stepOffset = 0;  // Reset step offset for next run
  // Note: mcmcCompletedSuccessfully is set by onMCMCSamplingComplete/Error callbacks

  // Update status
  if (mcmcCompletedSuccessfully) {
    statusLabel->setText("MCMC sampling completed successfully");
    logTextEdit->append("--- MCMC session completed successfully ---");
  } else {
    statusLabel->setText("MCMC sampling stopped or failed");
    logTextEdit->append("--- MCMC session ended (stopped early or failed) ---");
  }
}

void MCMCTab::resetToDefaults() {
  reset();
}

void MCMCTab::loadFromReduced() {
  if (!levelsTab_ || !segmentsTab_) {
    QMessageBox::warning(this, "Tab References Not Set",
                         "Cannot load parameters - tab references not properly set.");
    return;
  }

  mcmcParameters.clear();
  parametersTable->setRowCount(0);

  loadFromAZUREParams(true);
}

void MCMCTab::loadFromPhysical() {
  if (!levelsTab_ || !segmentsTab_) {
    QMessageBox::warning(this, "Tab References Not Set",
                         "Cannot load parameters - tab references not properly set.");
    return;
  }

  mcmcParameters.clear();
  parametersTable->setRowCount(0);

  loadFromAZUREParams(false);
}

void MCMCTab::loadFromAZUREParams(bool isRWA, std::string filename) {
  // Load parameters from current AZUREParams in the exact order expected by AZURECalcMCMC
  // This is used when "Use Reduced Widths" is checked

  AZURESetup *setup = qobject_cast<AZURESetup *>(window());
  if (!setup) {
    logTextEdit->append("Warning: Cannot access AZURESetup - falling back to GUI-based loading");
    return;
  }

  try {
    // Get current parameters from AZURESetup
    // This reads from the current .azr file and respects the parameter settings
    Config &config = setup->GetConfig();

    // CRITICAL: Initialize compound nucleus and data objects first
    CNuc *compound = new CNuc();
    EData *data = new EData();

    // Fill compound nucleus from current GUI configuration
    if (compound->Fill(config) == -1) {
      delete compound;
      delete data;
      logTextEdit->append("Error: Could not initialize compound nucleus from current configuration.");
      return;
    }

    // Fill data object if needed for parameter context
    if (config.paramMask & Config::CALCULATE_WITH_DATA) {
      if (data->Fill(config, compound) == -1) {
        delete compound;
        delete data;
        logTextEdit->append("Error: Could not initialize data structures from current configuration.");
        return;
      }
    } else {
      // Create theoretical points if no experimental data
      if (data->MakePoints(config, compound) == -1) {
        delete compound;
        delete data;
        logTextEdit->append("Error: Could not create theoretical data points.");
        return;
      }
    }

    AZUREParams azureParams;
    if (isRWA) {
      compound->TransformIn(config);
      compound->FillMnParams(azureParams.GetMinuitParams(), &config);
      data->FillMnParams(azureParams.GetMinuitParams());
      // Read the param.par file to get RWA parameters
      if (filename.empty()) {
        filename = config.outputdir + "param.par";
      }
      // Open and read param.par file
      std::ifstream paramFile(filename);
      // The second column contains the RWA values
      if (!paramFile) {
        delete compound;
        delete data;
        logTextEdit->append(QString("Error: Could not open RWA parameter file: %1")
                                .arg(QString::fromStdString(filename)));
        return;
      }
      // The second column contains the RWA values
      std::string line;
      int i = 0;
      while (std::getline(paramFile, line)) {
        std::istringstream iss(line);
        std::string paramName;
        double rwaValue;
        if (!(iss >> paramName >> rwaValue)) {
          logTextEdit->append(QString("Error: Invalid line in RWA parameter file: %1")
                                  .arg(QString::fromStdString(line)));
          continue;
        }
        if (i >= (int)azureParams.GetMinuitParams().Params().size()) {
          logTextEdit->append(QString("Warning: param.par has more entries than expected parameters (%1). Stopping at line %2.")
                                  .arg(azureParams.GetMinuitParams().Params().size())
                                  .arg(i));
          break;
        }
        azureParams.GetMinuitParams().SetValue(i, rwaValue);
        i++;
      }
    } else {
      compound->FillMnParams(azureParams.GetMinuitParams(), &config);
      data->FillMnParams(azureParams.GetMinuitParams());
    }

    // Automatic priors for the normalizations and energy shifts, keyed by
    // segment. These mirror exactly what AZURECalcMCMC::BuildAutoPriors()
    // derives at run time -- and the chi-squared penalties AZURECalc adds
    // during a Minuit fit -- so what the table shows is what gets sampled.
    // Collected before `data` is released.
    QMap<int, QPair<double, double>> autoNormPrior;  // key -> (mean, sigma)
    QMap<int, QPair<double, double>> autoShiftPrior;

    // Physical descriptions of every parameter, gathered here because the
    // compound nucleus and data are released a few lines below.
    QVector<QString> paramLabels;
    for (unsigned int i = 0; i < azureParams.GetMinuitParams().Params().size(); i++) {
      paramLabels.append(QString::fromStdString(
          AZURELabel::Parameter(compound, data, (int)i)));
    }
    {
      std::vector<ESegment> &segments = data->GetSegments();
      for (size_t s = 0; s < segments.size(); s++) {
        const int key = segments[s].GetSegmentKey();
        // GetNormError() is a percentage of the nominal normalization.
        const double nominalNorm = segments[s].GetNominalNorm();
        autoNormPrior[key] = qMakePair(
            nominalNorm, nominalNorm / 100.0 * segments[s].GetNormError());
        autoShiftPrior[key] = qMakePair(
            segments[s].GetNominalEnergyShift(),
            segments[s].GetEnergyShiftError());
      }
    }

    // Clean up temporary objects
    delete compound;
    delete data;

    const ROOT::Minuit2::MnUserParameters &minuitParams = azureParams.GetMinuitParams();

    int mcmcParamIndex = 0;

    // Iterate through ALL parameters in AZUREParams order
    for (unsigned int i = 0; i < minuitParams.Params().size(); i++) {
      // ONLY include non-fixed parameters (same logic as AZURECalcMCMC)
      if (!minuitParams.Parameter(i).IsFixed()) {
        QString paramName = QString::fromStdString(minuitParams.GetName(i));
        double value = minuitParams.Value(i);

        MCMCParameter param;
        // Real quantum numbers rather than a name like "width_1_2":
        // which level (J^pi, energy) and which channel (pair, L, S).
        param.name = (i < (unsigned int)paramLabels.size() && !paramLabels[i].isEmpty())
            ? paramLabels[i]
            : paramName;
        param.value = value;  // Use current RWA value
        param.priorMean = value;
        param.priorStd = std::abs(value) * 0.1;  // Use error if available, else 10%
        param.useGaussianPrior = false;
        param.autoPrior = false;

        // Determine category based on parameter name
        if (paramName.contains("norm", Qt::CaseInsensitive)) {
          param.category = "norm";
        } else if (paramName.contains("shift", Qt::CaseInsensitive)) {
          param.category = "shift";
        } else {
          if (isRWA) {
            param.category = "level_rwa";
          } else {
            param.category = "level";
          }
        }

        // Parameter names are "segment_<key>_norm" and
        // "segment_<key>_energy_shift" (see EData::FillMnParams).
        if (isAutoPriorCategory(param.category)) {
          const QStringList tokens = paramName.split('_');
          bool haveKey = false;
          int segmentKey = tokens.size() > 1 ? tokens[1].toInt(&haveKey) : 0;

          const QMap<int, QPair<double, double>> &source =
              (param.category == "norm") ? autoNormPrior : autoShiftPrior;

          if (haveKey && source.contains(segmentKey)) {
            param.priorMean = source[segmentKey].first;
            param.priorStd = source[segmentKey].second;
            // A segment with no quoted error gets no prior rather
            // than an invented width.
            param.useGaussianPrior = (param.priorStd > 0.0);
            param.autoPrior = true;
          }
        }

        param.minuitIndex = mcmcParamIndex++;  // Sequential index for MCMC array
        param.levelIndex = -1;
        param.channelIndex = -1;

        mcmcParameters.append(param);
      }
    }

    // Tell the user which parameters are still their responsibility.
    int nAuto = 0, nUserNeeded = 0, nUserSet = 0;
    for (const MCMCParameter &param : mcmcParameters) {
      if (param.autoPrior) {
        nAuto++;
      } else {
        nUserNeeded++;
        if (param.useGaussianPrior && param.priorStd > 0.0) nUserSet++;
      }
    }
    logTextEdit->append(QString("Automatic priors set for %1 normalization/energy-shift "
                                "parameter(s) from the quoted experimental errors")
                            .arg(nAuto));
    if (nUserSet < nUserNeeded) {
      logTextEdit->append(QString("Note: %1 of %2 level-energy/width parameter(s) have no "
                                  "prior. Define them in this table if the chain wanders.")
                              .arg(nUserNeeded - nUserSet)
                              .arg(nUserNeeded));
    }

    // Update parameter table display
    refreshParameterTable();

    logTextEdit->append(QString("Loaded %1 non-fixed parameters from current AZUREParams (RWA mode)")
                            .arg(mcmcParameters.size()));

    if (mcmcParameters.isEmpty()) {
      logTextEdit->append("Warning: No non-fixed parameters found in current configuration.");
    }

  } catch (const std::exception &e) {
    logTextEdit->append(QString("Error reading current AZUREParams: %1")
                            .arg(QString::fromStdString(e.what())));
    logTextEdit->append("Falling back to GUI-based parameter loading");
  }
}

void MCMCTab::updateParameterValuesFromFittingTab() {
  // Use the new method to read parameter values directly from the .azr file
  // This matches the approach used by ParameterLimitsManager
  if (!readParameterValuesFromConfigFile()) {
    logTextEdit->append("Warning: Could not read parameter values from config file, using default values");
  }
}

int MCMCTab::calculateStepOffset() const {
  // Check if Fresh Start is enabled - if so, offset is 0
  if (freshStartCheckBox->isChecked()) {
    return 0;
  }

  // Try to get output directory from parent AZURESetup
  QString outputDir;
  AZURESetup *setup = qobject_cast<AZURESetup *>(window());
  if (setup && !setup->GetConfig().configfile.empty()) {
    // Use the same directory as the config file as default
    QString configPath = QString::fromStdString(setup->GetConfig().outputdir);
    QFileInfo configInfo(configPath);
    outputDir = configInfo.absolutePath();
  } else {
    // Default to current working directory
    outputDir = QDir::currentPath() + QStringLiteral("/output/");
  }

  QString filePath = QDir(outputDir).absoluteFilePath("samples.mcmc");

  QFile file(filePath);
  if (!file.exists() || !file.open(QIODevice::ReadOnly | QIODevice::Text)) {
    return 0;  // No existing samples
  }

  QTextStream in(&file);
  int lineCount = 0;
  bool hasHeader = false;

  while (!in.atEnd()) {
    QString line = in.readLine().trimmed();
    if (line.isEmpty()) continue;

    // Skip header line
    if (!hasHeader && line.contains("step,walker,logprob,loglikelihood,logprior")) {
      hasHeader = true;
      continue;
    }

    lineCount++;
  }

  file.close();

  // Each step produces nWalkers lines, so step count = lineCount / nWalkers
  int nWalkers = nWalkersSpinBox->value();
  int existingSteps = (nWalkers > 0) ? lineCount / nWalkers : 0;

  return existingSteps;
}

void MCMCTab::parameterItemChanged(QTableWidgetItem *item) {
  // Handle parameter table changes - use the actual changed item
  if (!item) return;

  updateParameterFromTable(parametersTable->row(item));
}

void MCMCTab::updateParameterFromTable(int row) {
  const int index = paramIndexForRow(row);
  if (index < 0 || index >= mcmcParameters.size()) return;

  MCMCParameter &param = mcmcParameters[index];

  // Update prior mean (column 2) if item exists
  QTableWidgetItem *priorMeanItem = parametersTable->item(row, 2);
  if (priorMeanItem) {
    bool ok;
    double value = priorMeanItem->text().toDouble(&ok);
    if (ok) param.priorMean = value;
  }

  // Update prior std (column 3) if item exists
  QTableWidgetItem *priorStdItem = parametersTable->item(row, 3);
  if (priorStdItem) {
    bool ok;
    double value = priorStdItem->text().toDouble(&ok);
    if (ok) param.priorStd = value;
  }

  // Update Gaussian prior checkbox (column 4) if item exists
  QTableWidgetItem *gaussianItem = parametersTable->item(row, 4);
  if (gaussianItem) {
    param.useGaussianPrior = (gaussianItem->checkState() == Qt::Checked);
  }
}

bool MCMCTab::writeMCMCSettings(QTextStream &outStream) {
  outStream << "<mcmc>\n";
  outStream << "nwalkers " << nWalkersSpinBox->value() << "\n";
  outStream << "nsteps " << nStepsSpinBox->value() << "\n";
  outStream << "chainspread " << chainSpreadSpinBox->value() << "\n";
  outStream << "energyspread " << energySpreadSpinBox->value() << "\n";
  outStream << "nthreads " << nThreadsSpinBox->value() << "\n";
  outStream << "usereducedwidths " << (useReducedWidthsCheckBox->isChecked() ? "1" : "0") << "\n";

  outStream << "<parameters>\n";
  for (const MCMCParameter &param : mcmcParameters) {
    outStream << param.name << " "
              << param.value << " "
              << param.priorMean << " "
              << param.priorStd << " "
              << (param.useGaussianPrior ? "1" : "0") << "\n";
  }
  outStream << "</parameters>\n";
  outStream << "</mcmc>\n";

  return true;
}

bool MCMCTab::readMCMCSettings(QTextStream &inStream) {
  QString line;
  bool inMCMCSection = true;  // We're already in the MCMC section when this is called
  bool inParametersSection = false;

  // Clear existing parameters
  mcmcParameters.clear();

  while (!inStream.atEnd()) {
    line = inStream.readLine().trimmed();

    if (line == "</mcmc>") {
      inMCMCSection = false;
      break;
    }

    if (!inMCMCSection) continue;

    if (line == "<parameters>") {
      inParametersSection = true;
      continue;
    } else if (line == "</parameters>") {
      inParametersSection = false;
      continue;
    }

    if (inParametersSection) {
      // Parse parameter line: "name currentValue priorMean priorStd useGaussianPrior"
      // Parameter names can contain spaces, so we need to parse from the right
      QStringList parts = line.split(' ', Qt::SkipEmptyParts);
      if (parts.size() >= 5) {
        MCMCParameter param;

        // The last 4 parts are currentValue, priorMean, priorStd, useGaussianPrior
        int numParts = parts.size();
        param.value = parts[numParts - 4].toDouble();
        param.priorMean = parts[numParts - 3].toDouble();
        param.priorStd = parts[numParts - 2].toDouble();
        param.useGaussianPrior = (parts[numParts - 1].toInt() == 1);

        // Everything except the last 4 parts is the parameter name
        QStringList nameParts = parts.mid(0, numParts - 4);
        param.name = nameParts.join(" ");

        // The stored name still identifies norm/shift parameters, whose
        // priors come from the data rather than from the user. The
        // values stored here are the automatic ones written at save
        // time; a run re-derives them from the segments regardless.
        param.category = categoryFromStoredName(param.name);
        param.autoPrior = isAutoPriorCategory(param.category);
        param.minuitIndex = -1;
        param.levelIndex = -1;
        param.channelIndex = -1;

        mcmcParameters.append(param);
      } else if (parts.size() >= 4) {
        // Handle legacy format for backward compatibility
        MCMCParameter param;

        // The last 3 parts are priorMean, priorStd, useGaussianPrior
        int numParts = parts.size();
        param.priorMean = parts[numParts - 3].toDouble();
        param.priorStd = parts[numParts - 2].toDouble();
        param.useGaussianPrior = (parts[numParts - 1].toInt() == 1);

        // Everything except the last 3 parts is the parameter name
        QStringList nameParts = parts.mid(0, numParts - 3);
        param.name = nameParts.join(" ");

        param.value = param.priorMean;  // Legacy: set current value to prior mean
        param.category = categoryFromStoredName(param.name);
        param.autoPrior = isAutoPriorCategory(param.category);
        param.minuitIndex = -1;
        param.levelIndex = -1;
        param.channelIndex = -1;

        mcmcParameters.append(param);
      }
    } else {
      // Parse settings lines: "key value"
      QStringList parts = line.split(' ', Qt::SkipEmptyParts);
      if (parts.size() >= 2) {
        QString key = parts[0];
        QString value = parts[1];

        if (key == "nwalkers") {
          nWalkersSpinBox->setValue(value.toInt());
        } else if (key == "nsteps") {
          nStepsSpinBox->setValue(value.toInt());
        } else if (key == "chainspread") {
          chainSpreadSpinBox->setValue(value.toDouble());
        } else if (key == "energyspread") {
          // Older .azr files predate this setting; the spin box keeps
          // its default in that case.
          energySpreadSpinBox->setValue(
              std::min(value.toDouble(), kMaxEnergySpreadKeV));
        } else if (key == "nthreads") {
          nThreadsSpinBox->setValue(value.toInt());
        } else if (key == "usereducedwidths") {
          useReducedWidthsCheckBox->setChecked(value.toInt() == 1);
        }
      }
    }
  }

  // Update the parameters table display
  refreshParameterTable();

  // After reading MCMC settings, update parameter values from fitting tab
  // This ensures we have the most recent parameter values from the parameterSettings section
  if (mcmcParameters.isEmpty()) {
    loadFromPhysical();
  } else {
    // Update parameter values from fitting tab for existing parameters
    updateParameterValuesFromFittingTab();
  }

  logTextEdit->append(QString("Loaded MCMC settings with %1 parameters").arg(mcmcParameters.size()));

  return true;
}

void MCMCTab::showInfo(int which, QString title) {
  if (which < 0 || which >= 4) return;

  if (!infoDialog[which]) {
    infoDialog[which] = new InfoDialog(infoText[which], this);
  }
  infoDialog[which]->show();
  infoDialog[which]->raise();
  infoDialog[which]->activateWindow();
}


// MCMC worker thread slot implementations


void MCMCTab::onMCMCProgressUpdated(int currentStep, int totalSteps, double logProbability,
                                    double logLikelihood, double logPrior) {
  // Initialize timing on first call
  if (stepOffset == 0 && currentStep == 1) {
    stepOffset = calculateStepOffset();
    startTime = std::chrono::steady_clock::now();
    lastUpdateTime = startTime;
  }

  // Adjust step numbers for display
  int displayCurrentStep = currentStep + stepOffset;
  int displayTotalSteps = nStepsSpinBox->value();  // Use the actual target from UI, not remaining steps

  // Update progress bar
  progressBar->setMaximum(displayTotalSteps);
  progressBar->setValue(displayCurrentStep);
  progressBar->setVisible(true);

  // Calculate estimated time remaining
  QString timeEstimate = "";
  auto currentTime = std::chrono::steady_clock::now();
  if (displayCurrentStep > stepOffset + 1) {  // Need at least 2 steps to estimate
    auto elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime);
    int stepsCompleted = displayCurrentStep - stepOffset;
    int stepsRemaining = displayTotalSteps - displayCurrentStep;

    if (stepsCompleted > 0 && stepsRemaining > 0 && elapsedTime.count() > 0) {
      double stepsPerSecond = static_cast<double>(stepsCompleted) / elapsedTime.count();
      int estimatedSecondsRemaining = static_cast<int>(stepsRemaining / stepsPerSecond);

      // Format time estimate
      if (estimatedSecondsRemaining < 60) {
        timeEstimate = QString(" - Estimated time %1s").arg(estimatedSecondsRemaining);
      } else if (estimatedSecondsRemaining < 3600) {
        int minutes = estimatedSecondsRemaining / 60;
        int seconds = estimatedSecondsRemaining % 60;
        timeEstimate = QString(" - Estimated time %1m %2s").arg(minutes).arg(seconds);
      } else {
        int hours = estimatedSecondsRemaining / 3600;
        int minutes = (estimatedSecondsRemaining % 3600) / 60;
        timeEstimate = QString(" - Estimated time %1h %2m").arg(hours).arg(minutes);
      }
    }
  }

  // Update status label with time estimation
  double percentage = (double)displayCurrentStep / displayTotalSteps * 100.0;
  statusLabel->setText(QString("MCMC Step %1 of %2 (%3%)%4")
                           .arg(displayCurrentStep)
                           .arg(displayTotalSteps)
                           .arg(percentage, 0, 'f', 1)
                           .arg(timeEstimate));

  // Update current iteration display
  currentIterationLabel->setText(QString::number(displayCurrentStep));

  // Update log probability displays with exponential notation for better precision
  // Use 'e' format for scientific notation with 6 decimal places

  logProbLabel->setText(QString::number(logProbability, 'e', 6));
  logLikelihoodLabel->setText(QString::number(logLikelihood, 'e', 6));
  logPriorLabel->setText(QString::number(logPrior, 'e', 6));

  // Update log with periodic messages
  if (currentStep % 10 == 0 || currentStep == totalSteps) {
    logTextEdit->append(QString("Step %1: LogProb=%2 (LogL=%3, LogPrior=%4)")
                            .arg(currentStep)
                            .arg(logProbability, 0, 'f', 3)
                            .arg(logLikelihood, 0, 'f', 3)
                            .arg(logPrior, 0, 'f', 3));
  }
}

void MCMCTab::onMCMCLogMessage(const QString &message) {
  logTextEdit->append(message);
}

void MCMCTab::onMCMCSamplingComplete(const std::vector<std::vector<double>> &samples) {
#ifdef USE_MCMC
  // Mark as completed successfully
  mcmcCompletedSuccessfully = true;

  // Store the samples
  mcmcSamples = samples;

  statusLabel->setText("MCMC sampling completed successfully");
  logTextEdit->append(QString("MCMC sampling finished. Generated %1 samples.").arg(samples.size()));

  // Refresh statistics from file to ensure accuracy
  refreshResultsFromFile();
#endif
}


void MCMCTab::onMCMCSamplingError(const QString &errorMessage) {
  // Mark as failed/incomplete
  mcmcCompletedSuccessfully = false;

  statusLabel->setText("MCMC sampling failed");
  logTextEdit->append(QString("ERROR: %1").arg(errorMessage));

  // Refresh statistics from file even when there's an error
  refreshResultsFromFile();

  QMessageBox::critical(this, "MCMC Sampling Error", errorMessage);
}

bool MCMCTab::readParameterValuesFromConfigFile(const QString &configFilePath) {
  QString configFile = configFilePath;

  // If no config file path provided, try to get it from parent AZURESetup
  if (configFile.isEmpty()) {
    // Try to get the current config file from the parent widget
    AZURESetup *setup = qobject_cast<AZURESetup *>(window());
    if (setup) {
      configFile = QString::fromStdString(setup->GetConfig().configfile);
    }

    if (configFile.isEmpty()) {
      logTextEdit->append("No config file path available");
      return false;
    }
  }

  std::ifstream file(configFile.toStdString().c_str());
  if (!file.is_open()) {
    logTextEdit->append(QString("Cannot open config file: %1").arg(configFile));
    return false;
  }

  std::string line;
  bool inParameterSettings = false;
  int updatedParameters = 0;

  while (std::getline(file, line)) {
    // Trim whitespace
    size_t start = line.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) continue;
    line = line.substr(start);
    size_t end = line.find_last_not_of(" \t\r\n");
    if (end != std::string::npos) line = line.substr(0, end + 1);

    if (line == "<parameterSettings>") {
      inParameterSettings = true;
      continue;
    } else if (line == "</parameterSettings>" || line.substr(0, 1) == "<") {
      inParameterSettings = false;
      continue;
    }

    if (inParameterSettings && !line.empty() && line.substr(0, 1) != "#") {
      std::istringstream iss(line);
      std::string paramName, category;
      double value, lowerLimit, upperLimit, error;
      int nuisance;

      if (iss >> paramName >> value >> lowerLimit >> upperLimit >> error >> nuisance >> category) {
        // Find the corresponding MCMC parameter and update its value
        for (MCMCParameter &mcmcParam : mcmcParameters) {
          if (mcmcParam.name == QString::fromStdString(paramName)) {
            mcmcParam.value = value;
            updatedParameters++;
            break;
          }
        }
      }
    }
  }

  file.close();

  if (updatedParameters > 0) {
    logTextEdit->append(QString("Updated %1 parameter values from config file").arg(updatedParameters));
    // Update the parameter table display. Rebuilding is simpler and keeps the
    // row-to-parameter mapping in one place.
    refreshParameterTable();
    return true;
  } else {
    logTextEdit->append("No matching parameters found in config file");
    return false;
  }
}

void MCMCTab::addIncompleteResultsWarning() {
  if (!mcmcCompletedSuccessfully) {
    // Add warning message to statistics table
    int rowCount = statisticsTable->rowCount();
    statisticsTable->insertRow(0);

    QTableWidgetItem *warningItem = new QTableWidgetItem("WARNING: MCMC was stopped early. Statistics may not be reliable.");
    warningItem->setBackground(QColor(255, 255, 0, 100));  // Light yellow background
    warningItem->setFont(QFont("Arial", 10, QFont::Bold));
    statisticsTable->setItem(0, 0, warningItem);
    statisticsTable->setSpan(0, 0, 1, statisticsTable->columnCount());
  }
}

void MCMCTab::onMCMCResultsUpdated(int currentStep, int totalSteps, const std::vector<std::vector<double>> &samples) {
#ifdef USE_MCMC
  // Update results every 100 iterations
  if (samples.empty() || samples.size() < 10) return;

  // Update statistics table with current samples
  statisticsTable->setRowCount(mcmcParameters.size());

  // Calculate statistics for each parameter
  for (int paramIdx = 0; paramIdx < mcmcParameters.size() && paramIdx < (int)samples[0].size(); paramIdx++) {
    const MCMCParameter &param = mcmcParameters[paramIdx];

    // Extract values for this parameter from all samples
    std::vector<double> values;
    for (const auto &sample : samples) {
      if (paramIdx < (int)sample.size()) {
        values.push_back(sample[paramIdx]);
      }
    }

    if (!values.empty()) {
      // Calculate statistics
      std::sort(values.begin(), values.end());
      double mean = 0.0;
      for (double v : values) mean += v;
      mean /= values.size();

      double variance = 0.0;
      for (double v : values) variance += (v - mean) * (v - mean);
      variance /= std::max(1, (int)values.size() - 1);
      double std_dev = sqrt(variance);

      // Percentiles
      double p2_5 = values[values.size() * 0.025];
      double p50 = values[values.size() * 0.5];
      double p97_5 = values[values.size() * 0.975];

      // Populate table
      statisticsTable->setItem(paramIdx, 0, new QTableWidgetItem(param.name));
      statisticsTable->setItem(paramIdx, 1, new QTableWidgetItem(QString::number(mean, 'g', 6)));
      statisticsTable->setItem(paramIdx, 2, new QTableWidgetItem(QString::number(std_dev, 'g', 6)));
      statisticsTable->setItem(paramIdx, 3, new QTableWidgetItem(QString::number(p2_5, 'g', 6)));
      statisticsTable->setItem(paramIdx, 4, new QTableWidgetItem(QString::number(p50, 'g', 6)));
      statisticsTable->setItem(paramIdx, 5, new QTableWidgetItem(QString::number(p97_5, 'g', 6)));
    }
  }

  // Update chains display with summary information
  if (currentStep % 500 == 0) {  // Update chain display less frequently
    QString chainStatus;
    chainStatus += QString("Chains Status at Step %1/%2:\n").arg(currentStep).arg(totalSteps);
    chainStatus += QString("Total samples: %1\n").arg(samples.size());

    // Add warning if MCMC was stopped early or failed
    if (!mcmcCompletedSuccessfully && currentStep < totalSteps) {
      chainStatus += "⚠️ WARNING: MCMC sampling was stopped early or failed.\n";
      chainStatus += "Results may not be reliable for statistical inference.\n\n";
    }

    // Show last few samples for each parameter
    if (!samples.empty()) {
      chainStatus += "Recent samples:\n";
      int showSamples = std::min(5, (int)samples.size());
      for (int i = samples.size() - showSamples; i < (int)samples.size(); i++) {
        const auto &sample = samples[i];
        chainStatus += QString("Sample %1:").arg(i);
        for (int j = 0; j < std::min((int)sample.size(), mcmcParameters.size()); j++) {
          chainStatus += QString(" %1=%2").arg(mcmcParameters[j].name).arg(sample[j], 0, 'g', 4);
        }
        chainStatus += "\n";
      }
    }

    chainsTextEdit->setText(chainStatus);
  }
#endif
}

void MCMCTab::refreshResultsFromFile() {
  std::vector<std::vector<double>> samples;
  QString filePath;

  if (!loadSamplesFromFile(samples, filePath)) {
    return;  // Error message already shown by loadSamplesFromFile
  }

  if (samples.empty()) {
    QMessageBox::information(this, "No Samples", "The samples file exists but contains no data.");
    return;
  }

  // Calculate and display statistics from loaded samples
  calculateStatisticsFromSamples(samples, true);  // Assume complete since loaded from file

  logTextEdit->append(QString("Loaded %1 samples from %2").arg(samples.size()).arg(filePath));

  // Update chains display with summary
  QString chainSummary;
  chainSummary += QString("Loaded from file: %1\n").arg(filePath);
  chainSummary += QString("Total samples: %1\n").arg(samples.size());

  if (!samples.empty()) {
    chainSummary += QString("Parameters per sample: %1\n\n").arg(samples[0].size());
    chainSummary += "Sample preview (last 10 samples):\n";

    int startSample = std::max(0, (int)samples.size() - 10);
    for (int i = startSample; i < (int)samples.size(); i++) {
      const auto &sample = samples[i];
      QString line = QString("Sample %1:").arg(i);
      for (int j = 0; j < std::min(5, (int)sample.size()); j++) {  // Show first 5 parameters
        line += QString("\t%1").arg(sample[j], 0, 'g', 6);
      }
      if (sample.size() > 5) {
        line += QString("\t... (%1 more)").arg(sample.size() - 5);
      }
      chainSummary += line + "\n";
    }
  }

  chainsTextEdit->setText(chainSummary);
}

bool MCMCTab::loadSamplesFromFile(std::vector<std::vector<double>> &samples, QString &filePath) {
  samples.clear();

  // Try to get output directory from parent AZURESetup
  QString outputDir;
  AZURESetup *setup = qobject_cast<AZURESetup *>(window());
  if (setup && !setup->GetConfig().configfile.empty()) {
    // Use the same directory as the config file as default
    QString configPath = QString::fromStdString(setup->GetConfig().outputdir);
    QFileInfo configInfo(configPath);
    outputDir = configInfo.absolutePath();
  } else {
    // Default to current working directory
    outputDir = QDir::currentPath() + QStringLiteral("/output/");
  }

  filePath = QDir(outputDir).absoluteFilePath("samples.mcmc");

  QFile file(filePath);
  if (!file.exists()) {
    QMessageBox::information(this, "No Samples File",
                             QString("No samples.mcmc file found at:\n%1\n\n"
                                     "Run MCMC sampling first to generate samples.")
                                 .arg(filePath));
    return false;
  }

  if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
    QMessageBox::critical(this, "File Error",
                          QString("Cannot open samples file:\n%1\n\nError: %2")
                              .arg(filePath)
                              .arg(file.errorString()));
    return false;
  }

  QTextStream in(&file);
  int lineNumber = 0;
  int expectedParameterCount = -1;

  while (!in.atEnd()) {
    QString line = in.readLine().trimmed();
    lineNumber++;

    if (line.isEmpty() || line.startsWith("#") || line.contains("step,walker,logprob,loglikelihood,logprior")) {
      // Skip empty lines, comments, and header line
      continue;
    }

    QStringList values = line.split(",", Qt::SkipEmptyParts);

    if (values.size() < 6) {  // Need at least step,walker,logprob,loglikelihood,logprior + 1 parameter
      continue;
    }

    // Skip the first 5 columns (step,walker,logprob,loglikelihood,logprior)
    // and convert only the parameter values to doubles
    std::vector<double> sample;
    bool conversionOk = true;
    for (int i = 5; i < values.size(); i++) {
      bool ok;
      double value = values[i].trimmed().toDouble(&ok);
      if (!ok) {
        logTextEdit->append(QString("Warning: Invalid number '%1' at line %2, skipping line")
                                .arg(values[i].trimmed())
                                .arg(lineNumber));
        conversionOk = false;
        break;
      }
      sample.push_back(value);
    }

    if (!conversionOk) {
      continue;
    }

    // Validate parameter count consistency
    if (expectedParameterCount == -1) {
      expectedParameterCount = sample.size();
      if (expectedParameterCount != mcmcParameters.size()) {
        logTextEdit->append(QString("Warning: Sample file has %1 parameters but MCMC tab expects %2")
                                .arg(expectedParameterCount)
                                .arg(mcmcParameters.size()));
      }
    } else if ((int)sample.size() != expectedParameterCount) {
      logTextEdit->append(QString("Warning: Line %1 has %2 parameters, expected %3, skipping line")
                              .arg(lineNumber)
                              .arg(sample.size())
                              .arg(expectedParameterCount));
      continue;
    }

    samples.push_back(sample);
  }

  file.close();
  return true;
}

void MCMCTab::calculateStatisticsFromSamples(const std::vector<std::vector<double>> &samples, bool isComplete) {
  if (samples.empty() || mcmcParameters.isEmpty()) {
    return;
  }

  // Set completion status
  mcmcCompletedSuccessfully = isComplete;

  // Clear existing results
  statisticsTable->setRowCount(mcmcParameters.size());

  // Calculate statistics for each parameter
  for (int paramIdx = 0; paramIdx < mcmcParameters.size() && paramIdx < (int)samples[0].size(); paramIdx++) {
    const MCMCParameter &param = mcmcParameters[paramIdx];

    // Extract values for this parameter from all samples
    std::vector<double> values;
    for (const auto &sample : samples) {
      if (paramIdx < (int)sample.size()) {
        values.push_back(sample[paramIdx]);
      }
    }

    if (!values.empty()) {
      // Calculate statistics
      std::sort(values.begin(), values.end());
      double mean = 0.0;
      for (double v : values) mean += v;
      mean /= values.size();

      double variance = 0.0;
      for (double v : values) variance += (v - mean) * (v - mean);
      variance /= std::max(1, (int)values.size() - 1);
      double std_dev = sqrt(variance);

      // Percentiles
      double p2_5 = values[values.size() * 0.025];
      double p50 = values[values.size() * 0.5];
      double p97_5 = values[values.size() * 0.975];

      // Populate table
      statisticsTable->setItem(paramIdx, 0, new QTableWidgetItem(param.name));
      statisticsTable->setItem(paramIdx, 1, new QTableWidgetItem(QString::number(mean, 'g', 6)));
      statisticsTable->setItem(paramIdx, 2, new QTableWidgetItem(QString::number(std_dev, 'g', 6)));
      statisticsTable->setItem(paramIdx, 3, new QTableWidgetItem(QString::number(p2_5, 'g', 6)));
      statisticsTable->setItem(paramIdx, 4, new QTableWidgetItem(QString::number(p50, 'g', 6)));
      statisticsTable->setItem(paramIdx, 5, new QTableWidgetItem(QString::number(p97_5, 'g', 6)));
    }
  }

  // Add warning if results are incomplete
  addIncompleteResultsWarning();
}
