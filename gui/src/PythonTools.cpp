#include "PythonTools.h"

#include <QDesktopServices>
#include <QFileDialog>
#include <QFileInfo>
#include <QInputDialog>
#include <QJsonDocument>
#include <QJsonObject>
#include <QMenu>
#include <QMessageBox>
#include <QSettings>
#include <QUrl>

PythonTools::PythonTools(QWidget *window) : QObject(window), window_(window) {}

PythonTools::~PythonTools() { stopWebPage(); }

QString PythonTools::interpreter() { return QSettings().value("pythonInterpreter", "python3").toString(); }

void PythonTools::setInterpreter(const QString &python) {
  QSettings().setValue("pythonInterpreter", python.trimmed().isEmpty() ? QString("python3") : python.trimmed());
}

bool PythonTools::available(QString &why) {
  const QString python = interpreter();
  if (!checked_ || checkedInterpreter_ != python) {
    QProcess probe;
    probe.start(python, {"-c", "import pyazr, pyazr.cli"});
    if (!probe.waitForStarted(10000)) unavailable_ = tr("%1 does not start (Configure > Directories)").arg(python);
    else if (!probe.waitForFinished(60000) || probe.exitCode() != 0)
      unavailable_ = tr("pyazr is not importable with %1: %2")
                         .arg(python, QString::fromLocal8Bit(probe.readAllStandardError()).trimmed().split('\n').last());
    else unavailable_.clear();
    checked_ = true;
    checkedInterpreter_ = python;
  }
  why = unavailable_;
  return unavailable_.isEmpty();
}

bool PythonTools::start(const QStringList &arguments, const QString &directory, QString &why) {
  if (running()) {
    why = tr("a pyazr command is running");
    return false;
  }
  QProcess *p = new QProcess(this);
  process_ = p;
  p->setWorkingDirectory(directory);
  p->setProcessChannelMode(QProcess::MergedChannels);
  connect(p, &QProcess::readyRead, this, [this, p]() {
    while (p->canReadLine()) emit line(QString::fromLocal8Bit(p->readLine()).chopped(1), false);
  });
  connect(p, &QProcess::finished, this, [this, p](int code) {
    const QString rest = QString::fromLocal8Bit(p->readAll());
    if (!rest.isEmpty()) emit line(rest, false);
    emit line(tr("[pyazr ended with exit %1]").arg(code), false);
    p->deleteLater();
    emit finished(code);
  });
  const QStringList args = QStringList{"-m", "pyazr"} + arguments;
  emit line(QString("%1 %2").arg(interpreter(), args.join(' ')), false);
  p->start(interpreter(), args);
  if (!p->waitForStarted(10000)) {
    why = tr("%1 does not start").arg(interpreter());
    p->deleteLater();
    return false;
  }
  return true;
}

bool PythonTools::openWebPage(const QString &projectFile, QString &why) {
  if (web_ && web_->state() != QProcess::NotRunning) {
    if (!webUrl_.isEmpty()) QDesktopServices::openUrl(QUrl(webUrl_));
    return true;
  }
  QProcess *p = new QProcess(this);
  web_ = p;
  webUrl_.clear();
  webOutput_.clear();
  p->setWorkingDirectory(QFileInfo(projectFile).absolutePath());
  connect(p, &QProcess::readyReadStandardOutput, this, [this, p]() {
    webOutput_ += QString::fromLocal8Bit(p->readAllStandardOutput());
    // --json: one JSON document with the URL as soon as the server listens on 127.0.0.1
    const QJsonDocument doc = QJsonDocument::fromJson(webOutput_.toUtf8());
    if (doc.isObject() && webUrl_.isEmpty()) {
      webUrl_ = doc.object().value("url").toString();
      if (!webUrl_.isEmpty()) {
        emit line(tr("The web page of this project: %1").arg(webUrl_), false);
        emit webPageReady(webUrl_);
      }
    }
  });
  connect(p, &QProcess::readyReadStandardError, this,
          [this, p]() { emit line(QString::fromLocal8Bit(p->readAllStandardError()).trimmed(), false); });
  p->start(interpreter(), {"-m", "pyazr", "ui", QFileInfo(projectFile).absoluteFilePath(), "--json"});
  if (!p->waitForStarted(10000)) {
    why = tr("%1 does not start").arg(interpreter());
    return false;
  }
  return true;
}

void PythonTools::stopWebPage() {
  if (web_ && web_->state() != QProcess::NotRunning) {
    web_->terminate();
    if (!web_->waitForFinished(5000)) web_->kill();
  }
}

void PythonTools::fillMenu(QMenu *menu, const QString &projectFile) {
  menu->clear();
  QString why;
  const bool ok = available(why);
  const QString dir = QFileInfo(projectFile).absolutePath();
  auto add = [&](const QString &title, const QString &help, std::function<void()> act) {
    QAction *a = menu->addAction(title);
    a->setEnabled(ok);
    a->setToolTip(ok ? help : why);
    a->setStatusTip(ok ? help : why);
    QObject::connect(a, &QAction::triggered, this, [this, act]() { act(); });
  };
  auto warn = [this](bool started, const QString &reason) {
    if (!started) QMessageBox::warning(window_, tr("pyazr"), reason);
  };
  menu->setToolTipsVisible(true);
  add(tr("Export ENDF-6..."), tr("python -m pyazr endf write: the model at its start values as an ENDF-6 material"),
      [=]() {
        const QString to = QFileDialog::getSaveFileName(window_, tr("ENDF-6 file"), dir, tr("ENDF-6 (*.endf *.txt *)"));
        QString reason;
        if (!to.isEmpty()) warn(start({"endf", "write", projectFile, to}, dir, reason), reason);
      });
  add(tr("Import an ENDF-6 Evaluation as a New Project..."),
      tr("python -m pyazr endf to-project: the resonance parameters of MF=2 (LRF=7) as a format-3 project"), [=]() {
        const QString from = QFileDialog::getOpenFileName(window_, tr("ENDF-6 file"), dir);
        if (from.isEmpty()) return;
        const QString to = QFileDialog::getSaveFileName(window_, tr("New project"), dir, tr("AZURE2 projects (*.azr)"));
        QString reason;
        if (!to.isEmpty()) warn(start({"endf", "to-project", from, to}, dir, reason), reason);
      });
  add(tr("ENDF-6 File Summary..."), tr("python -m pyazr endf summary"), [=]() {
    const QString from = QFileDialog::getOpenFileName(window_, tr("ENDF-6 file"), dir);
    QString reason;
    if (!from.isEmpty()) warn(start({"endf", "summary", from}, dir, reason), reason);
  });
  menu->addSeparator();
  add(tr("Import an EXFOR Data Set..."), tr("python -m pyazr import exfor ID: the data set's columns and values"), [=]() {
    const QString id = QInputDialog::getText(window_, tr("EXFOR"), tr("Data set (e.g. O2599004):"));
    QString reason;
    if (!id.trimmed().isEmpty()) warn(start({"import", "exfor", id.trimmed()}, dir, reason), reason);
  });
  add(tr("ENSDF Levels of a Nuclide..."), tr("python -m pyazr import ensdf NUCLIDE"), [=]() {
    const QString n = QInputDialog::getText(window_, tr("ENSDF"), tr("Nuclide (e.g. 13N):"));
    QString reason;
    if (!n.trimmed().isEmpty()) warn(start({"import", "ensdf", n.trimmed()}, dir, reason), reason);
  });
  menu->addSeparator();
  add(tr("Build the Coulomb Table"), tr("python -m pyazr coulomb-table build (about 5 min on 4 threads, 88 MB)"),
      [=]() {
        QString reason;
        warn(start({"coulomb-table", "build"}, dir, reason), reason);
      });
  add(tr("Check the Coulomb Table"), tr("python -m pyazr coulomb-table check"), [=]() {
    QString reason;
    warn(start({"coulomb-table", "check"}, dir, reason), reason);
  });
  menu->addSeparator();
  add(tr("Open the Web Page of This Project"),
      tr("python -m pyazr ui FILE: the web page on 127.0.0.1, with the samplers (NUTS, emcee) and the band jobs"),
      [=]() {
        QString reason;
        if (!openWebPage(projectFile, reason)) warn(false, reason);
      });
  QAction *samplers = menu->addAction(tr("Samplers (NUTS, emcee) and Band Jobs"));
  samplers->setEnabled(false);
  samplers->setToolTip(tr("Not commands of pyazr's command line: run them in the web page (the item above) or from "
                          "Python (pyazr.examples.fit_emcee, Session.objective().log_prob, pyazr.bands)."));
  if (!ok) {
    QAction *reason = menu->addAction(why);
    reason->setEnabled(false);
  }
}
