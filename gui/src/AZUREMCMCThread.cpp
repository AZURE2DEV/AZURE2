#include "AZUREMCMCThread.h"
#include "MCMCTab.h"
#include "ECAmplitudeCache.h"

#ifdef USE_MCMC
#include "AZURECalcMCMC.h"
#include "CNuc.h"
#include "EData.h"
#include "AZUREParams.h"
#include <algorithm>
#include <numeric>
#include <fstream>
#include <random>
#endif

#include <QMutexLocker>
#include <QTableWidget>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QTableWidgetItem>
#include <QTextEdit>
#include <QLabel>
#include <QCheckBox>
#include <QFile>

// Global worker pointer for static callback bridge
static AZUREMCMCWorker* g_currentWorker = nullptr;

// AZUREMCMCWorker implementation - does the actual MCMC work
AZUREMCMCWorker::AZUREMCMCWorker(MCMCTab* mcmcTab, const Config& configure) 
    : config_(configure), mcmcTab_(mcmcTab), stopRequested_(false) {
}

void AZUREMCMCWorker::run() {
#ifndef USE_MCMC
    emit samplingError("MCMC functionality not enabled. Rebuild with USE_MCMC=ON.");
    return;
#else
    try {
        stopRequested_ = false;
        AZURECalcMCMC::ClearStop(); // Clear any previous stop flags
        emit logMessage("Initializing MCMC calculation...");
        
        // Initialize caching systems
        if (config_.paramMask & Config::USE_EXTERNAL_CAPTURE) {
            // Cleanup any existing caches
            CleanupECAmplitudeCache();
            
            // Initialize ECAmplitude cache
            InitializeECAmplitudeCache();
        }
        
        // Create fresh data and compound objects (like AZUREMain does)
        EData* data = new EData();
        CNuc* compound = new CNuc();
        
        emit logMessage("Filling compound nucleus...");
        if(compound->Fill(config_) == -1) {
            delete compound;
            delete data;
            emit samplingError("Could not fill compound nucleus from file.");
            return;
        } else if(compound->NumPairs() == 0 || compound->NumJGroups() == 0) {
            delete compound;
            delete data;
            emit samplingError("No nuclear data exists. Calculation not possible.");
            return;
        }
        
        if(!(config_.paramMask & Config::CALCULATE_REACTION_RATE)) {
            emit logMessage("Filling data structures...");
            if(config_.paramMask & Config::CALCULATE_WITH_DATA) {
                if(data->Fill(config_, compound) == -1) {
                    delete compound;
                    delete data;
                    emit samplingError("Could not fill data object from file.");
                    return;
                } else if(data->NumSegments() == 0) {
                    delete compound;
                    delete data;
                    emit samplingError("There is no data provided.");
                    return;
                }
            } else {
                if(data->MakePoints(config_, compound) == -1) {
                    delete compound;
                    delete data;
                    emit samplingError("Could not fill data object from file.");
                    return;
                } else if(data->NumSegments() == 0) {
                    delete compound;
                    delete data;
                    emit samplingError("Extrapolation segments produce no data.");
                    return;
                }
            }
        } else {
            if(!compound->IsPairKey(config_.rateParams.entrancePair)||
               !compound->IsPairKey(config_.rateParams.exitPair)) {
                delete compound;
                delete data;
                emit samplingError("Reaction rate pairs do not exist in compound nucleus.");
                return;
            } else {
                compound->GetPair(compound->GetPairNumFromKey(config_.rateParams.entrancePair))->SetEntrance();
            }
        }
        
        // Initialize compound nucleus object
        emit logMessage("Initializing compound nucleus...");
        try {
            compound->Initialize(config_);
        } catch (const std::exception& e) {
            delete compound;
            delete data;
            emit samplingError(QString("Error during compound initialization: %1").arg(e.what()));
            return;
        }
        
        // Create parameters (same as AZUREMain does after compound initialization)
        emit logMessage("Creating parameters...");
        AZUREParams params;
        compound->FillMnParams(params.GetMinuitParams());
        data->FillMnParams(params.GetMinuitParams());
        
        // Initialize data object (after parameter creation, same as AZUREMain)
        emit logMessage("Initializing data object...");
        if(data->Initialize(compound, config_) == -1) {
            delete compound;
            delete data;
            emit samplingError("Failed to initialize data object");
            return;
        }
        
        // Create MCMC calculator with freshly initialized objects
        emit logMessage("Creating MCMC calculator...");
        AZURECalcMCMC* mcmcCalculator = new AZURECalcMCMC(NULL, NULL, config_);
        int out = mcmcCalculator->Initialize( );
        std::cout << "MCMC Calculator initialized with status: " << out << std::endl;

        // Check if we should use reduced widths (RWA)
        bool useReducedWidths = mcmcTab_->useReducedWidthsCheckBox->isChecked();
        emit logMessage(QString("Using %1 for MCMC fitting").arg(useReducedWidths ? "reduced width amplitudes (RWA)" : "physical parameters"));
        
        // Get MCMC parameters from MCMCTab
        std::vector<double> initialParams;
        std::vector<double> priorMeans;
        std::vector<double> priorStds;
        std::vector<bool> usePriors;
        
        for(int i = 0; i < mcmcTab_->parametersTable->rowCount(); i++) {
            // Get parameter value
            double value = mcmcTab_->parametersTable->item(i, 1)->text().toDouble();
            initialParams.push_back(value);
            
            // Get prior information
            double priorMean = mcmcTab_->parametersTable->item(i, 2)->text().toDouble();
            double priorStd = mcmcTab_->parametersTable->item(i, 3)->text().toDouble();
            bool usePrior = (mcmcTab_->parametersTable->item(i, 4)->checkState() == Qt::Checked);
            
            priorMeans.push_back(priorMean);
            priorStds.push_back(priorStd);
            usePriors.push_back(usePrior);
        }
        
        // Set priors in MCMC calculator
        mcmcCalculator->SetPriors(priorMeans, priorStds, usePriors);
        
        // Get sampling parameters from MCMCTab
        int nWalkers = mcmcTab_->nWalkersSpinBox->value();
        int nSteps = mcmcTab_->nStepsSpinBox->value();
        double chainSpread = mcmcTab_->chainSpreadSpinBox->value();
        int nThreads = mcmcTab_->nThreadsSpinBox->value();
        bool freshStart = mcmcTab_->freshStartCheckBox->isChecked();
        
        // Handle Fresh Start option - delete existing samples file if checked
        if (freshStart) {
            std::string samplesFile = config_.outputdir + "samples.mcmc";
            QFile existingFile(QString::fromStdString(samplesFile));
            if (existingFile.exists()) {
                existingFile.remove();
                emit logMessage("Fresh start requested - removed existing samples file");
            }
        }
        
        emit logMessage(QString("Starting MCMC sampling with %1 walkers, %2 steps, %3 threads")
                       .arg(nWalkers).arg(nSteps).arg(nThreads));
        
        // Run MCMC sampling using existing method with progress monitoring
        emit logMessage("Starting MCMC sampling...");
        
        std::vector<std::vector<double>> samples;
        
        // Use the enhanced MCMC sampling with progress callbacks
        emit logMessage(QString("Running %1 MCMC steps with %2 walkers...").arg(nSteps).arg(nWalkers));
        
        // Create a custom progress-tracking MCMC runner
        emit logMessage("Starting custom MCMC runner with progress callbacks...");
        
        // Set up static callbacks to forward progress to this worker
        // Store worker instance for static callback access
        g_currentWorker = this;
        
        // Set up progress callback (every 10 steps) for log probability updates
        AZURECalcMCMC::SetGUIProgressCallback([](int step, int totalSteps, double logProb, double logLikelihood, double logPrior) {
            if(g_currentWorker) {
                emit g_currentWorker->progressUpdated(step, totalSteps, logProb, logLikelihood, logPrior);
            }
        });
        
        // Set up iteration callback (every step) for iteration number updates
        AZURECalcMCMC::SetGUIIterationCallback([](int step, int totalSteps) {
            if(g_currentWorker) {
                // Emit a separate signal for iteration updates without log probability data
                emit g_currentWorker->iterationUpdated(step, totalSteps);
            }
        });
        
        // Set up results callback (every 100 steps) for results tab updates
        AZURECalcMCMC::SetGUIResultsCallback([](int step, int totalSteps, const std::vector<std::vector<double>>& samples) {
            if(g_currentWorker) {
                emit g_currentWorker->resultsUpdated(step, totalSteps, samples);
            }
        });
        
        try {
            if (useReducedWidths) {
                // Use reduced width amplitudes (RWA) for fitting
                emit logMessage("Running MCMC with reduced width amplitudes...");
                mcmcCalculator->RunMCMCSampling(nWalkers, nSteps, initialParams, samples, chainSpread, nThreads, true);
            } else {
                // Use physical parameters for fitting
                emit logMessage("Running MCMC with physical parameters...");
                mcmcCalculator->RunMCMCSampling(nWalkers, nSteps, initialParams, samples, chainSpread, nThreads, false);
            }
            
            emit logMessage(QString("MCMC sampling completed! Generated %1 samples")
                           .arg(samples.size()));
            
            emit samplingComplete(samples);
        } catch(const std::exception& e) {
            emit samplingError(QString("MCMC sampling failed: %1").arg(e.what()));
            return;
        }
        
        // Clear callbacks and global worker pointer
        AZURECalcMCMC::SetGUIProgressCallback(nullptr);
        AZURECalcMCMC::SetGUIIterationCallback(nullptr);
        AZURECalcMCMC::SetGUIResultsCallback(nullptr);
        g_currentWorker = nullptr;
        
        // Cleanup
        delete mcmcCalculator;
        delete compound;
        delete data;
        
        // Cleanup caching systems
        if (config_.paramMask & Config::USE_EXTERNAL_CAPTURE) {
            CleanupECAmplitudeCache();
        }
        
    } catch (const std::exception& e) {
        emit samplingError(QString("MCMC sampling failed with exception: %1").arg(e.what()));
    } catch (...) {
        emit samplingError("Unknown error during MCMC sampling");
    }
    
    emit done();
#endif
}


// AZUREMCMCThread implementation - simple thread manager like AZUREMainThread
AZUREMCMCThread::AZUREMCMCThread(MCMCTab* mcmcTab, const Config& config, QObject* parent)
    : QThread(parent), config_(config), mcmcTab_(mcmcTab), worker_(mcmcTab, config), stopRequested_(false) {
    
    // Connect worker signals to our signals
    connect(this, &AZUREMCMCThread::readyToRun, &worker_, &AZUREMCMCWorker::run);
    connect(&worker_, &AZUREMCMCWorker::done, this, &AZUREMCMCThread::quit);
    connect(&worker_, &AZUREMCMCWorker::logMessage, this, [this](const QString& message) {
        // Forward to MCMCTab
        mcmcTab_->logTextEdit->append(message);
    });
    connect(&worker_, &AZUREMCMCWorker::samplingError, this, [this](const QString& error) {
        mcmcTab_->logTextEdit->append(QString("Error: %1").arg(error));
    });
    connect(&worker_, &AZUREMCMCWorker::samplingComplete, this, [this](const std::vector<std::vector<double>>& samples) {
        mcmcTab_->logTextEdit->append(QString("MCMC sampling completed with %1 samples").arg(samples.size()));
        // Forward to MCMCTab for results processing
        mcmcTab_->onMCMCSamplingComplete(samples);
    });
    connect(&worker_, &AZUREMCMCWorker::progressUpdated, this, [this](int currentStep, int totalSteps, double logProbability, double logLikelihood, double logPrior) {
        // Forward enhanced progress updates to MCMCTab
        mcmcTab_->onMCMCProgressUpdated(currentStep, totalSteps, logProbability, logLikelihood, logPrior);
    });
    connect(&worker_, &AZUREMCMCWorker::iterationUpdated, this, [this](int currentStep, int totalSteps) {
        // currentStep now contains the actual function evaluation count from g_sample_count
        mcmcTab_->currentIterationLabel->setText(QString::number(currentStep));
    });
    connect(&worker_, &AZUREMCMCWorker::resultsUpdated, this, [this](int currentStep, int totalSteps, const std::vector<std::vector<double>>& samples) {
        // Forward results updates to MCMCTab for results display
        mcmcTab_->onMCMCResultsUpdated(currentStep, totalSteps, samples);
    });
    
    worker_.moveToThread(this);
}

void AZUREMCMCThread::stop() {
    stopRequested_ = true;
    worker_.requestStop();
}

void AZUREMCMCThread::stopMCMC() {
    stop();
    // Force thread termination if needed
    if(!wait(5000)) {  // Wait 5 seconds for clean shutdown
        terminate();
        wait();
    }
}

void AZUREMCMCThread::run() {
    emit readyToRun();
    exec();
}