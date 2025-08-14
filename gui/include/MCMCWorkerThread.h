#ifndef MCMCWORKERTHREAD_H
#define MCMCWORKERTHREAD_H

#include <QThread>
#include <QObject>
#include <QMutex>
#include <QWaitCondition>
#include <vector>

#ifdef USE_MCMC
#include "AZURECalcMCMC.h"
#include "Config.h"
#include "CNuc.h"
#include "EData.h"
#endif

struct MCMCSamplingParameters {
    int nWalkers;
    int nSteps;
    int burnIn;
    double chainSpreadPercent;
    double temperature;
    bool useParallel;
    std::vector<double> initialParams;
    std::vector<QString> paramNames;
};

class MCMCWorker : public QObject {
    Q_OBJECT

public:
    MCMCWorker(const MCMCSamplingParameters& params, const Config& config);
    ~MCMCWorker();
    
    void stop();
    bool isRunning() const;

signals:
    void progressUpdated(int step, int totalSteps, double acceptanceRate);
    void logMessage(const QString& message);
    void samplingComplete(const std::vector<std::vector<double>>& samples);
    void samplingError(const QString& errorMessage);
    void finished();

public slots:
    void runSampling();

private:
    MCMCSamplingParameters params_;
    Config config_;
    
#ifdef USE_MCMC
    CNuc* compound_;
    EData* data_;
#endif
    
    mutable QMutex mutex_;
    bool stopRequested_;
    bool running_;
    
#ifdef USE_MCMC
    AZURECalcMCMC* mcmcCalculator_;
#endif
    
    void emitProgress(int step, int totalSteps, double acceptanceRate = 0.0);
    void emitLogMessage(const QString& message);
};

class MCMCWorkerThread : public QThread {
    Q_OBJECT

public:
    MCMCWorkerThread(const MCMCSamplingParameters& params,
                     const Config& config, QObject* parent = nullptr);
    ~MCMCWorkerThread();
    
    void stopSampling();
    bool isRunning() const;

signals:
    void progressUpdated(int step, int totalSteps, double acceptanceRate);
    void logMessage(const QString& message);
    void samplingComplete(const std::vector<std::vector<double>>& samples);
    void samplingError(const QString& errorMessage);

protected:
    void run() override;

private:
    MCMCWorker* worker_;
    MCMCSamplingParameters params_;
    Config config_;
};

#endif // MCMCWORKERTHREAD_H