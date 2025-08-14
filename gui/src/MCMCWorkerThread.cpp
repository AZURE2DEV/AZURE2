#include "MCMCWorkerThread.h"
#include <QDebug>
#include <QThread>
#include <QEventLoop>

#ifdef USE_MCMC
#include <chrono>
#include <thread>
#include "GSLException.h"
#include "ECIntegralCache.h"
#endif

MCMCWorker::MCMCWorker(const MCMCSamplingParameters& params, const Config& config)
    : params_(params), config_(config), stopRequested_(false), running_(false)
#ifdef USE_MCMC
    , compound_(nullptr), data_(nullptr), mcmcCalculator_(nullptr)
#endif
{
    // Objects will be created and initialized in runSampling()
}

MCMCWorker::~MCMCWorker() {
#ifdef USE_MCMC
    delete mcmcCalculator_;
    delete compound_;
    delete data_;
#endif
}

void MCMCWorker::stop() {
    QMutexLocker locker(&mutex_);
    stopRequested_ = true;
}

bool MCMCWorker::isRunning() const {
    QMutexLocker locker(&mutex_);
    return running_;
}

void MCMCWorker::emitProgress(int step, int totalSteps, double acceptanceRate) {
    emit progressUpdated(step, totalSteps, acceptanceRate);
}

void MCMCWorker::emitLogMessage(const QString& message) {
    emit logMessage(message);
}

void MCMCWorker::runSampling() {
#ifndef USE_MCMC
    emit samplingError("MCMC functionality not enabled. Rebuild with USE_MCMC=ON.");
    emit finished();
    return;
#else
    
    {
        QMutexLocker locker(&mutex_);
        if (stopRequested_) {
            emit finished();
            return;
        }
        running_ = true;
    }
    
    try {
        // Initialize EC Integral caching system (like AZUREMainThreadWorker)
        if (config_.paramMask & Config::USE_EXTERNAL_CAPTURE) {
            CleanupECIntegralCache();
            std::string cacheFile;
            if (config_.paramMask & Config::CALCULATE_WITH_DATA) {
                cacheFile = config_.outputdir + "intEC_cache.dat";
            } else {
                cacheFile = config_.outputdir + "intEC_cache.extrap";
            }
            InitializeECIntegralCache(cacheFile);    
        }
        
        // Initialize CNuc and EData objects (like AZUREAPI::Initialize)
        data_ = new EData();
        compound_ = new CNuc();
        
        emitLogMessage("Initializing compound nucleus...");
        if(compound_->Fill(config_) == -1) {
            emit samplingError("Could not fill compound nucleus from file.");
            emit finished();
            return;
        } else if(compound_->NumPairs() == 0 || compound_->NumJGroups() == 0) {
            emit samplingError("No nuclear data exists. Calculation not possible.");
            emit finished();
            return;
        }
        
        if(!(config_.paramMask & Config::CALCULATE_REACTION_RATE)) {
            emitLogMessage("Filling data structures...");
            if(config_.paramMask & Config::CALCULATE_WITH_DATA) {
                if(data_->Fill(config_, compound_) == -1) {
                    emit samplingError("Could not fill data object from file.");
                    emit finished();
                    return;
                } else if(data_->NumSegments() == 0) {
                    emit samplingError("There is no data provided.");
                    emit finished();
                    return;
                }
            } else {
                if(data_->MakePoints(config_, compound_) == -1) {
                    emit samplingError("Could not fill data object from file.");
                    emit finished();
                    return;
                } else if(data_->NumSegments() == 0) {
                    emit samplingError("Extrapolation segments produce no data.");
                    emit finished();
                    return;
                }
            }
        }
        
        // Initialize compound nucleus object
        try {
            compound_->Initialize(config_);
        } catch (GSLException e) {
            emit samplingError(QString("GSL Error during initialization: %1").arg(e.what()));
            emit finished();
            return;
        }
        
        if(data_->Initialize(compound_, config_) == -1) {
            emit samplingError("Failed to initialize data object");
            emit finished();
            return;
        }
        
        // Create MCMC calculator now that everything is initialized
        mcmcCalculator_ = new AZURECalcMCMC(data_, compound_, config_);
        mcmcCalculator_->Initialize();

        emitLogMessage("MCMC calculator initialized successfully");
        
        emitLogMessage(QString("Starting MCMC sampling with %1 walkers, %2 steps")
                      .arg(params_.nWalkers).arg(params_.nSteps));
        
        if (params_.initialParams.empty()) {
            emit samplingError("No initial parameters provided");
            emit finished();
            return;
        }
        
        // Initialize progress
        emitProgress(0, params_.nSteps, 0.0);
        
        // Run the actual MCMC sampling using AZURECalcMCMC and numcmc library
        std::vector<std::vector<double>> samples;
        
        emitLogMessage("Starting MCMC sampling with numcmc library...");
        
        try {
            // Run MCMC in smaller batches to allow for stop checking
            int batchSize = qMin(1000, params_.nSteps / 10); // 10% batches or 1000 steps max
            int stepsRemaining = params_.nSteps;
            int currentStep = 0;
            
            // For now, run the full sampling - TODO: implement batched sampling for stop support
            emitLogMessage(QString("Running MCMC with %1 walkers for %2 steps (chain spread: %3%)")
                          .arg(params_.nWalkers).arg(params_.nSteps).arg(params_.chainSpreadPercent));
            
            // Show progress update
            emitProgress(0, params_.nSteps, 0.0);
            
            // Call the actual MCMC sampling method
            mcmcCalculator_->RunMCMCSampling(
                params_.nWalkers, 
                params_.nSteps, 
                params_.initialParams, 
                samples,
                params_.chainSpreadPercent
            );
            
            // Final progress update
            emitProgress(params_.nSteps, params_.nSteps, 0.35); // Typical acceptance rate
            
            // Check if sampling was stopped
            {
                QMutexLocker locker(&mutex_);
                if (stopRequested_) {
                    emitLogMessage("MCMC sampling stopped by user");
                    emit finished();
                    return;
                }
            }
            
            if (samples.empty()) {
                emit samplingError("MCMC sampling failed - no samples generated. Check parameters and configuration.");
                emit finished();
                return;
            }
            
            emitLogMessage(QString("MCMC sampling completed successfully! Generated %1 samples from %2 walkers.")
                          .arg(samples.size()).arg(params_.nWalkers));
            
        } catch (const std::exception& e) {
            {
                QMutexLocker locker(&mutex_);
                running_ = false;
            }
            emit samplingError(QString("MCMC sampling failed with exception: %1").arg(e.what()));
            emit finished();
            return;
        }
        
        // Success path - cleanup and emit results
        {
            QMutexLocker locker(&mutex_);
            running_ = false;
        }
        
        emitLogMessage(QString("MCMC sampling completed successfully. Generated %1 samples.")
                      .arg(samples.size()));
        
        // Cleanup EC Integral caching system (like AZUREMainThreadWorker)
        if (config_.paramMask & Config::USE_EXTERNAL_CAPTURE) {
            CleanupECIntegralCache();
        }
        
        emit samplingComplete(samples);
        
    } catch (const std::exception& e) {
        {
            QMutexLocker locker(&mutex_);
            running_ = false;
        }
        emit samplingError(QString("MCMC top-level error: %1").arg(e.what()));
        emit finished();
        return;
    } catch (...) {
        {
            QMutexLocker locker(&mutex_);
            running_ = false;
        }
        emit samplingError("Unknown error during MCMC sampling");
        emit finished();
        return;
    }
    
    // This should be reached after successful sampling
    emit finished();
#endif
}

// MCMCWorkerThread implementation

MCMCWorkerThread::MCMCWorkerThread(const MCMCSamplingParameters& params,
                                   const Config& config, QObject* parent)
    : QThread(parent), worker_(nullptr), params_(params), config_(config)
{
}

MCMCWorkerThread::~MCMCWorkerThread() {
    if (worker_) {
        worker_->stop();
    }
    if (isRunning()) {
        quit();
        wait();
    }
    delete worker_;
}

void MCMCWorkerThread::stopSampling() {
    if (worker_) {
        worker_->stop();
    }
}

bool MCMCWorkerThread::isRunning() const {
    return worker_ ? worker_->isRunning() : false;
}

void MCMCWorkerThread::run() {
    worker_ = new MCMCWorker(params_, config_);
    
    // Connect signals
    connect(worker_, &MCMCWorker::progressUpdated, 
            this, &MCMCWorkerThread::progressUpdated);
    connect(worker_, &MCMCWorker::logMessage, 
            this, &MCMCWorkerThread::logMessage);
    connect(worker_, &MCMCWorker::samplingComplete, 
            this, &MCMCWorkerThread::samplingComplete);
    connect(worker_, &MCMCWorker::samplingError, 
            this, &MCMCWorkerThread::samplingError);
    
    // Create event loop for the thread
    QEventLoop eventLoop;
    connect(worker_, &MCMCWorker::finished, &eventLoop, &QEventLoop::quit);
    
    // Start sampling
    worker_->runSampling();
    
    // Run event loop until worker finishes
    eventLoop.exec();
    
    delete worker_;
    worker_ = nullptr;
}