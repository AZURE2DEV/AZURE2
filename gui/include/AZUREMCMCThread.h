#ifndef AZUREMCMCTHREAD_H
#define AZUREMCMCTHREAD_H

#include <QThread>
#include "Config.h"

#ifdef USE_MCMC
#include <vector>
#include "AZURECalcMCMC.h"
class MCMCTab;
#endif

// Simple worker class like AZUREMainThreadWorker
class AZUREMCMCWorker : public QObject {
    Q_OBJECT

public:
    AZUREMCMCWorker(MCMCTab* mcmcTab, const Config& configure);
    const Config& configure() const { return config_; }
    void requestStop() { 
        stopRequested_ = true; 
        AZURECalcMCMC::RequestStop();
    }
    bool isStopRequested() const { return stopRequested_; }

signals:
    void done();
    void logMessage(const QString& message);
    void samplingError(const QString& errorMessage);
    void samplingComplete(const std::vector<std::vector<double>>& samples);
    void progressUpdated(int currentStep, int totalSteps, double logProbability, double logLikelihood, double logPrior);
    void iterationUpdated(int currentStep, int totalSteps);
    void resultsUpdated(int currentStep, int totalSteps, const std::vector<std::vector<double>>& samples);

public slots:
    void run();

private:
    
    Config config_;
    MCMCTab* mcmcTab_;
    bool stopRequested_;
};

// Simple thread class like AZUREMainThread
class AZUREMCMCThread : public QThread {
    Q_OBJECT

public:
    AZUREMCMCThread(MCMCTab* mcmcTab, const Config& config, QObject* parent = nullptr);
    
    const Config& configure() const { return worker_.configure(); }
    void stop();
    AZUREMCMCWorker* getWorker() { return &worker_; }

signals:
    void readyToRun();

public slots:
    void stopMCMC();

protected:
    void run() override;

private:
    Config config_;
    MCMCTab* mcmcTab_;
    AZUREMCMCWorker worker_;
    bool stopRequested_;
};

#endif // AZUREMCMCTHREAD_H