#ifndef AZUREMAINTHREAD_H
#define AZUREMAINTHREAD_H

#include <QThread>
#include <QEventLoop>
#include <QPushButton>

#include "TextEditBuffer.h"
#include "FilteredTextEdit.h"
#include "RunTab.h"
#include "AZUREMain.h"
#include "Config.h"
#include "ECAmplitudeCache.h"
#include <chrono>
#include <cstdio>
#include <iostream>

class AZUREMainThreadWorker : public QObject {
Q_OBJECT

 public:
  AZUREMainThreadWorker(const Config& configure) :
  azureMain_(configure) {};
 signals:
  void done();
 public slots:
  void run() {
    // Initialize and cleanup caches for new calculations
    if (azureMain_.configure().paramMask & Config::USE_EXTERNAL_CAPTURE) {
      // Cleanup any existing caches
      CleanupECAmplitudeCache();
      
      // Initialize ECIntegral cache with proper file
      std::string cacheFile;
      if (azureMain_.configure().paramMask & Config::CALCULATE_WITH_DATA) {
        cacheFile = azureMain_.configure().outputdir + "intEC_cache.dat";
      } else {
        cacheFile = azureMain_.configure().outputdir + "intEC_cache.extrap";
      }
      
      // Initialize ECAmplitude cache
      InitializeECAmplitudeCache();
    }
    
    auto calcStart = std::chrono::steady_clock::now();
    azureMain_();
    int totalSec = (int)std::chrono::duration_cast<std::chrono::seconds>(
                         std::chrono::steady_clock::now() - calcStart).count();
    char timeStr[16];
    std::snprintf(timeStr, sizeof(timeStr), "%d:%02d:%02d",
                  totalSec / 3600, (totalSec % 3600) / 60, totalSec % 60);
    azureMain_.configure().outStream << std::endl
        << "Calculation completed in " << timeStr << "." << std::endl;

    // Cleanup caches after calculations
    if (azureMain_.configure().paramMask & Config::USE_EXTERNAL_CAPTURE) {
      CleanupECAmplitudeCache();
    }
    
    emit done();
  };
 private:
  AZUREMain azureMain_;
};

class AZUREMainThread : public QThread {
  Q_OBJECT
 public:
  AZUREMainThread(RunTab *tab, const Config& configure) :
  stream_(&buffer_), configure_(stream_), worker_(configure_) {
    configure_.configfile = configure.configfile;
    configure_.paramMask = configure.paramMask;
    configure_.screenCheckMask = configure.screenCheckMask;
    configure_.fileCheckMask = configure.fileCheckMask;
    configure_.chiVariance = configure.chiVariance;
    configure_.outputdir = configure.outputdir;
    configure_.checkdir = configure.checkdir;
    configure_.paramfile = configure.paramfile;
    configure_.integralsfile = configure.integralsfile;
    configure_.rateParams = configure.rateParams;
    configure_.nloptAlgorithm = configure.nloptAlgorithm;
    connect(&buffer_,SIGNAL(updateLog(QString)),tab->runtimeText,SLOT(write(QString)));
    connect(tab->stopAZUREButton,SIGNAL(clicked()),this,SLOT(stopAZURE()));
    connect(this,SIGNAL(readyToRun()),&worker_,SLOT(run()));
    connect(&worker_,SIGNAL(done()),this,SLOT(quit()));
    worker_.moveToThread(this);
  };
  const Config& configure() const {return configure_;};
 signals:
  void readyToRun();
 public slots:
  void stopAZURE() {
    configure_.stopFlag=true;
  };
 protected:
  void run() {
    emit readyToRun();
    exec();
  };
 private:
  TextEditBuffer buffer_;
  std::ostream stream_;
  Config configure_;
  AZUREMainThreadWorker worker_;
};

#endif
