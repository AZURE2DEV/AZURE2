#ifndef MCMCTAB_H
#define MCMCTAB_H

#include <string>
#include <QWidget>
#include <QSignalMapper>
#include <QPointer>
#include <QTextStream>

#ifdef USE_MCMC
#include <vector>
class AZURECalcMCMC;
#endif

// Forward declarations
class InfoDialog;
class LevelsTab;
class SegmentsTab;
class Config;
class CNuc;
class EData;
class LevelsModel;
class ChannelsModel;
class SegmentsDataModel;

QT_BEGIN_NAMESPACE

class QGroupBox;
class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QLabel;
class QSpinBox;
class QDoubleSpinBox;
class QPushButton;
class QProgressBar;
class QTextEdit;
class QTabWidget;
class QTableWidget;
class QTableWidgetItem;
class QCheckBox;
class QComboBox;

QT_END_NAMESPACE

// Upper bound on the initial spread of a level-energy parameter, in keV. Level
// energies cannot be scattered by a percentage of their value the way widths
// are: tens of keV of scatter puts walkers on unrelated resonance structures
// and the ensemble never contracts.
const double kMaxEnergySpreadKeV = 1.0;

struct MCMCParameter {
    QString name;
    double value = 0.0;
    double priorMean = 0.0;
    double priorStd = 0.0;
    bool useGaussianPrior = false;
    QString category; // "level", "level_rwa", "norm", "shift"
    int minuitIndex = -1;

    // Normalizations and energy shifts carry a prior derived from the error
    // quoted in the segment definition, so the user does not edit theirs.
    bool autoPrior = false;

    // For level parameters
    int levelIndex = -1;
    int channelIndex = -1;
};

class MCMCTab : public QWidget {
    Q_OBJECT

public:
    MCMCTab(QWidget* parent = 0);
    friend class AZURESetup;
    void reset();
    bool writeMCMCSettings(QTextStream& outStream);
    bool readMCMCSettings(QTextStream& inStream);
    void setTabReferences(LevelsTab* levelsTab, SegmentsTab* segmentsTab);
    
    // Public member variables for AZURESetup access (like RunTab pattern)
    QSpinBox* nWalkersSpinBox;
    QSpinBox* nStepsSpinBox;
    QDoubleSpinBox* chainSpreadSpinBox;
    QDoubleSpinBox* energySpreadSpinBox;
    QSpinBox* nThreadsSpinBox;
    QPushButton* runButton;
    QPushButton* stopButton;
    QCheckBox* freshStartCheckBox;
    QCheckBox* useReducedWidthsCheckBox;
    QTableWidget* parametersTable;
    QTextEdit* logTextEdit;
    QLabel* currentIterationLabel;

    /*! Every varying parameter, in AZUREParams order.
     *
     * This -- not the table -- is the parameter list. The table shows only the
     * level energies and widths, the ones whose priors the user sets; the
     * normalizations and energy shifts are still sampled, with priors that
     * AZURECalcMCMC derives from the segment errors. The sampled vector has to
     * cover *all* varying parameters in order, so anything building it must read
     * this list rather than counting table rows. */
    const QList<MCMCParameter>& parameters() const {return mcmcParameters;};

public slots:
    void showInfo(int which = 0, QString title = "");
    
    // MCMC worker thread slots - public so AZUREMCMCThread can call them
    void onMCMCProgressUpdated(int currentStep, int totalSteps, double logProbability, double logLikelihood, double logPrior);
    void onMCMCLogMessage(const QString& message);
    void onMCMCSamplingComplete(const std::vector<std::vector<double>>& samples);
    void onMCMCSamplingError(const QString& errorMessage);
    void onMCMCResultsUpdated(int currentStep, int totalSteps, const std::vector<std::vector<double>>& samples);

private slots:
    void mcmcFinished();
    void parameterItemChanged(QTableWidgetItem* item);
    void resetToDefaults();
    void loadFromPhysical();
    void loadFromReduced(); // Load parameters from .sav file without transformation
    void refreshResultsFromFile(); // Load and calculate statistics from samples.mcmc file

private:
    void setupParameterTable();
    void refreshParameterTable(); // Rebuild the table rows from mcmcParameters
    int paramIndexForRow(int row) const; // Table row -> index into mcmcParameters
    static bool isAutoPriorCategory(const QString& category); // norm / shift
    static QString categoryFromStoredName(const QString& name); // for .azr reload
    void setupSamplingControls(QWidget* samplingWidget);
    void setupProgressControls();
    void updateParameterFromTable(int row);
    void updateParameterValuesFromFittingTab();
    bool readParameterValuesFromConfigFile(const QString& configFilePath = QString()); // Read parameter values from .azr file parameterSettings section
    int calculateStepOffset() const; // Calculate starting step from existing samples
    void addIncompleteResultsWarning(); // Add warning to statistics when MCMC was stopped early
    bool loadSamplesFromFile(std::vector<std::vector<double>>& samples, QString& filePath); // Load samples from samples.mcmc file
    void calculateStatisticsFromSamples(const std::vector<std::vector<double>>& samples, bool isComplete = false); // Calculate and display statistics
    QString createPrettyParameterName(const QString& paramName, class CNuc* compound, class EData* data, int paramIndex) const; // Create user-friendly parameter names
    void loadFromAZUREParams(bool isRWA, std::string filename=""); // Load parameters from current AZUREParams in reduced widths mode

    // Control buttons (private ones not accessed by AZURESetup)
    QPushButton* resetButton;
    QPushButton* loadButton;
    QPushButton* loadSavButton;
    
    // Progress tracking
    QProgressBar* progressBar;
    QLabel* statusLabel;
    QLabel* logProbLabel;
    QLabel* logLikelihoodLabel;
    QLabel* logPriorLabel;
    
    // Results display
    QTabWidget* resultsTabWidget;
    QTableWidget* statisticsTable;
    QTextEdit* chainsTextEdit;
    
    // Progress controls group
    QGroupBox* progressControlsGroup;
    
    // Info system
    QSignalMapper* mapper;
    QPushButton* infoButton[4]; // 4 sections: Parameters, Sampling, Progress, Results
    static const std::vector<QString> infoText;
    QPointer<InfoDialog> infoDialog[4];
    
    // Data storage
    QList<MCMCParameter> mcmcParameters;
    // Table row -> mcmcParameters index. The table omits the auto-prior
    // (normalization / energy shift) parameters, so the two are not 1:1.
    QList<int> tableRowToParam;
    
    // Tab references
    LevelsTab* levelsTab_;
    SegmentsTab* segmentsTab_;
    
    // MCMC state
    bool isRunning;
    int stepOffset; // Starting step number when resuming from existing samples
    bool mcmcCompletedSuccessfully; // Track if MCMC finished or was stopped early
    
    // Timing for estimated time remaining
    std::chrono::steady_clock::time_point startTime;
    std::chrono::steady_clock::time_point lastUpdateTime;
    
#ifdef USE_MCMC
    // MCMC results storage
    std::vector<std::vector<double>> mcmcSamples;
#endif
};

#endif // MCMCTAB_H