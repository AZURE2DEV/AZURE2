#ifndef MCMCTAB_H
#define MCMCTAB_H

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
class MCMCWorkerThread;
struct MCMCSamplingParameters;

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
class QCheckBox;
class QComboBox;

QT_END_NAMESPACE

struct MCMCParameter {
    QString name;
    double value;
    double priorMean;
    double priorStd;
    bool useGaussianPrior;
    QString category; // "level", "norm", "shift"
    int minuitIndex;
    
    // For level parameters
    int levelIndex;
    int channelIndex;
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
    QSpinBox* nThreadsSpinBox;
    QPushButton* runButton;
    QPushButton* stopButton;
    QCheckBox* freshStartCheckBox;
    QTableWidget* parametersTable;
    QTextEdit* logTextEdit;
    QLabel* currentIterationLabel;

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
    void parameterItemChanged();
    void resetToDefaults();
    void loadFromFittingTab();
    void refreshResultsFromFile(); // Load and calculate statistics from samples.mcmc file

private:
    void setupParameterTable();
    void setupSamplingControls(QWidget* samplingWidget);
    void setupProgressControls();
    void updateParameterFromTable(int row);
    void updateParameterValuesFromFittingTab();
    bool readParameterValuesFromConfigFile(const QString& configFilePath = QString()); // Read parameter values from .azr file parameterSettings section
    int calculateStepOffset() const; // Calculate starting step from existing samples
    void addIncompleteResultsWarning(); // Add warning to statistics when MCMC was stopped early
    bool loadSamplesFromFile(std::vector<std::vector<double>>& samples, QString& filePath); // Load samples from samples.mcmc file
    void calculateStatisticsFromSamples(const std::vector<std::vector<double>>& samples, bool isComplete = false); // Calculate and display statistics

    // Control buttons (private ones not accessed by AZURESetup)
    QPushButton* resetButton;
    QPushButton* loadButton;
    
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
    
    // MCMC worker thread
    MCMCWorkerThread* mcmcWorkerThread;
    
#ifdef USE_MCMC
    // MCMC results storage
    std::vector<std::vector<double>> mcmcSamples;
#endif
};

#endif // MCMCTAB_H