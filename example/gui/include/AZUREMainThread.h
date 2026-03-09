#ifndef AZUREMAINTHREAD_H
#define AZUREMAINTHREAD_H

#include "TextEditBuffer.h"
#include "FilteredTextEdit.h"
#include "RunTab.h"
#include "AZUREMain.h"
#include "Config.h"
#include <QThread>
#include <QEventLoop>
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
    azureMain_();
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
