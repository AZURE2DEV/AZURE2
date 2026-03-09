#include "AZURESetup.h"
#include "EditChecksDialog.h"
#include "EditDirsDialog.h"
#include "RunTab.h"
#include "Config.h"
#include "EditOptionsDialog.h"
#include "AZUREMainThread.h"
#include "AboutAZURE2Dialog.h"
#ifdef USE_QWT
#include "PlotTab.h"
#endif
#include <iostream>

struct SegPairs {int firstPair; int secondPair;};
extern bool readSegmentFile(const Config& configure,std::vector<SegPairs>& segPairs);
extern bool checkExternalCapture(Config& configure, const std::vector<SegPairs>& segPairs);
extern void startMessage(const Config& configure);
extern void exitMessage(const Config& configure);


AZURESetup::AZURESetup() : config(std::cout) {
  setMinimumSize(1000,640);

  tabWidget=new QTabWidget();  

  pairsTab=new PairsTab;

  levelsTab = new LevelsTab;
  levelsTab->setPairsModel(pairsTab->getPairsModel());
  connect(pairsTab,SIGNAL(pairAdded(int)),levelsTab,SLOT(updateChannelsPairAddedEdited()));
  connect(pairsTab,SIGNAL(pairEdited(int)),levelsTab,SLOT(updateChannelsPairAddedEdited()));
  connect(pairsTab,SIGNAL(pairRemoved(int)),levelsTab,SLOT(updateChannelsPairRemoved(int)));
  connect(levelsTab,SIGNAL(readNewPair(PairsData,int,bool)),pairsTab,SLOT(addPair(PairsData,int,bool)));
  connect(levelsTab,SIGNAL(readExistingPair(PairsData,int,bool)),pairsTab,SLOT(editPair(PairsData,int,bool)));

  segmentsTab = new SegmentsTab;
  segmentsTab->setPairsModel(pairsTab->getPairsModel());

  targetIntTab=new TargetIntTab;

  runTab = new RunTab();
  connect(runTab->calcButton,SIGNAL(clicked()),this,SLOT(SaveAndRun()));

#ifdef USE_QWT
  plotTab = new PlotTab(config,segmentsTab->getSegmentsDataModel(),segmentsTab->getSegmentsTestModel());
#endif 

  tabWidget->addTab(pairsTab,tr("&Particle Pairs"));
  tabWidget->addTab(levelsTab,tr("&Levels and Channels"));
  tabWidget->addTab(segmentsTab,tr("&Segments"));
  tabWidget->addTab(targetIntTab,tr("&Experimental Effects"));
  tabWidget->addTab(runTab,tr("&Calculate"));
#ifdef USE_QWT
  tabWidget->addTab(plotTab,tr("Pl&ot"));
#endif

  setCentralWidget(tabWidget);

  createActions();
  createMenus();


  setWindowTitle(tr("AZURE2 -- untitled"));
}

Config& AZURESetup::GetConfig() {
  return config;
}

void AZURESetup::createActions() {
  aboutAction = new QAction(tr("&About AZURE2..."),this);
  connect(aboutAction,SIGNAL(triggered()),this,SLOT(showAbout()));

  resetAction = new QAction(tr("&New Project"),this);
  resetAction->setShortcuts(QKeySequence::New);
  connect(resetAction,SIGNAL(triggered()),this,SLOT(reset()));

  quitAction = new QAction(tr("&Quit"),this);
  quitAction->setShortcuts(QKeySequence::Quit);
  connect(quitAction,SIGNAL(triggered()),this,SLOT(close()));

  openAction = new QAction(tr("&Open..."),this);
  openAction->setShortcuts(QKeySequence::Open);
  connect(openAction,SIGNAL(triggered()),this,SLOT(open()));

  saveAction = new QAction(tr("&Save"),this);
  saveAction->setShortcuts(QKeySequence::Save);
  connect(saveAction,SIGNAL(triggered()),this,SLOT(save()));

  saveAsAction = new QAction(tr("Save &As..."),this);
  saveAsAction->setShortcuts(QKeySequence::SaveAs);
  connect(saveAsAction,SIGNAL(triggered()),this,SLOT(saveAs()));

  for(int i=0;i<numRecent;i++) {
    recentFileActions[i] = new QAction(this);
    recentFileActions[i] -> setVisible(false);
    connect(recentFileActions[i],SIGNAL(triggered()),this,SLOT(openRecent()));
  }
  clearRecentAction = new QAction(tr("&Clear"),this);
  clearRecentAction->setVisible(false);
  connect(clearRecentAction,SIGNAL(triggered()),this,SLOT(clearRecent()));

  copyAction = new QAction(tr("&Copy"),this);
  copyAction->setShortcuts(QKeySequence::Copy);

  matrixActionGroup = new QActionGroup(this);
  aMatrixAction = new QAction(tr("&A-Matrix"),this);
  aMatrixAction->setCheckable(true);
  matrixActionGroup->addAction(aMatrixAction);
  rMatrixAction = new QAction(tr("&R-Matrix"),this);
  rMatrixAction->setCheckable(true);
  matrixActionGroup->addAction(rMatrixAction);
  aMatrixAction->setChecked(true);
  connect(matrixActionGroup,SIGNAL(triggered(QAction*)),this,SLOT(matrixChanged(QAction*)));

  editChecksAction = new QAction(tr("&Checks..."),this);
  connect(editChecksAction,SIGNAL(triggered()),this,SLOT(editChecks()));

  editDirsAction = new QAction(tr("&Directories..."),this);
  connect(editDirsAction,SIGNAL(triggered()),this,SLOT(editDirs()));
  
  editOptionsAction = new QAction(tr("&Runtime Options..."),this);
  connect(editOptionsAction,SIGNAL(triggered()),this,SLOT(editOptions()));

  showTabInfoAction = new QAction(tr("Show Documentation For Current Tab"),this);
  showTabInfoAction->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_D));
  connect(showTabInfoAction,SIGNAL(triggered()),this,SLOT(showTabInfo()));

  openAZURESiteAction = new QAction(tr("Open AZURE Website..."),this);
  connect(openAZURESiteAction,SIGNAL(triggered()),this,SLOT(openWebsite()));
}

void AZURESetup::createMenus() {
  fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(aboutAction);
  fileMenu->addSeparator();
  fileMenu->addAction(resetAction);
  fileMenu->addAction(openAction);
  recentFileMenu = fileMenu->addMenu(tr("Open &Recent..."));
  for(int i = 0; i < numRecent; i++ ) recentFileMenu->addAction(recentFileActions[i]);
  recentSeparator = recentFileMenu->addSeparator();
  recentFileMenu->addAction(clearRecentAction);
  updateRecent();
  fileMenu->addAction(saveAction);
  fileMenu->addAction(saveAsAction);
  fileMenu->addSeparator();
  fileMenu->addAction(quitAction);
  
  configMenu = menuBar()->addMenu(tr("Co&nfigure"));
  formalismMenu = configMenu->addMenu(tr("&Formalism"));
  formalismMenu->addAction(aMatrixAction);
  formalismMenu->addAction(rMatrixAction);
  configMenu->addAction(editChecksAction);
  configMenu->addAction(editDirsAction);
  configMenu->addAction(editOptionsAction);

  helpMenu = menuBar()->addMenu(tr("&Documentation"));
  helpMenu->addAction(showTabInfoAction);
  helpMenu->addAction(openAZURESiteAction);
}

void AZURESetup::updateRecent() {
  QSettings settings;
  QStringList files = settings.value("recentFileList").toStringList();

  int numFiles = qMin(files.size(),(int)numRecent);

  for(int i = 0; i<numFiles; i++) {
    recentFileActions[i]->setText(tr("&%1 %2").arg(i+1).arg(QFileInfo(files[i]).fileName()));
    recentFileActions[i]->setData(files[i]);
    recentFileActions[i]->setVisible(true);
  }
  
  for(int i = numFiles; i<numRecent; i++) recentFileActions[i]->setVisible(false);
  recentSeparator->setVisible(numFiles>0);
  clearRecentAction->setVisible(numFiles>0);
}

void AZURESetup::open() {
  QString filename = QFileDialog::getOpenFileName(this);
  if(!filename.isEmpty()) {
    if(!this->readFile(filename)) {
      reset();
      QMessageBox::information(this,
			       tr("Can't Access File"),
			       tr("An error was encountered while reading the file."));
    }
  }
}

void AZURESetup::open(QString filename) {
  if(!filename.isEmpty()) {
    if(!this->readFile(filename)) {
      reset();
      QMessageBox::information(this,
			       tr("Can't Access File"),
			       tr("An error was encountered while reading the file."));
    }
  }
}

void AZURESetup::openRecent() {
  QString filename = qobject_cast<QAction*>(sender())->data().toString();
  open(filename);
}

void AZURESetup::clearRecent() {
  QSettings settings;
  QStringList files = settings.value("recentFileList").toStringList();
  files.clear();
  settings.setValue("recentFileList",files);
  updateRecent();
}

bool AZURESetup::readFile(QString filename) {
  QFile file(filename);
  if(!file.open(QIODevice::ReadOnly)) return false;
  QFileInfo info(file);
  QString directory=info.absolutePath();

  reset();

  QTextStream in(&file);
  QString line("");
  
  while(line.trimmed()!=QString("<config>")&&!in.atEnd()) line = in.readLine();
  if(in.atEnd()) return false;
  if(!this->readConfig(in)) return false;

  line=QString("");
  while(line.trimmed()!=QString("<levels>")&&!in.atEnd()) line = in.readLine();
  if(in.atEnd()) return false;
  if(!levelsTab->readNuclearFile(in)) return false;
  
  line=QString("");
  while(line.trimmed()!=QString("<segmentsData>")&&!in.atEnd()) line = in.readLine();
  if(in.atEnd()) return false;
  if(!segmentsTab->readSegDataFile(in)) return false;

  line=QString("");
  while(line.trimmed()!=QString("<segmentsTest>")&&!in.atEnd()) line = in.readLine();
  if(in.atEnd()) return false;
  if(!segmentsTab->readSegTestFile(in)) return false;

  line=QString("");
  while(line.trimmed()!=QString("<targetInt>")&&!in.atEnd()) line = in.readLine();
  if(in.atEnd()) return false;
  if(!targetIntTab->readFile(in)) return false;

  line=QString("");
  while(line.trimmed()!=QString("<lastRun>")&&!in.atEnd()) line = in.readLine();
  if(!in.atEnd()) 
    if(!this->readLastRun(in)) return false;

  file.close();

  QFile file2(filename);
  if(!file2.open(QIODevice::ReadOnly)) return false;
  QTextStream in2(&file2);
  line=QString("");
  while(line.trimmed()!=QString("<externalCapture>")&&!in2.atEnd()) line = in2.readLine();
  if(!in2.atEnd()) {
    if(!pairsTab->parseOldECSection(in2)) return false;
  }
  file2.close();
  
  GetConfig().configfile=QDir::fromNativeSeparators(info.absoluteFilePath()).toStdString();
  setWindowTitle(QString("AZURE2 -- %1").arg(QString::fromStdString(GetConfig().configfile)));
  QDir::setCurrent(directory);

  QSettings settings;
  QStringList files = settings.value("recentFileList").toStringList();
  QString fullFileName = QDir::fromNativeSeparators(info.absoluteFilePath());
  files.removeAll(fullFileName);
  files.prepend(fullFileName);
  while(files.size()>numRecent) files.removeLast();

  settings.setValue("recentFileList",files);
  updateRecent();

  return true;
}

bool AZURESetup::readLastRun(QTextStream& inStream) {
  unsigned int paramMask;
  unsigned int useTempFile;
  unsigned int rateEntrancePair;
  unsigned int rateExitPair;
  QString paramFile;
  QString integralsFile;
  QString temperatureFile;
  QString dummyString;
  double minTemp;
  double maxTemp;
  double tempStep;

  inStream >> paramMask;dummyString=inStream.readLine();
  dummyString=inStream.readLine();paramFile=dummyString.trimmed();
  dummyString=inStream.readLine();integralsFile=dummyString.trimmed();
  inStream >> rateEntrancePair >> rateExitPair;dummyString=inStream.readLine();
  inStream >> useTempFile;dummyString=inStream.readLine();temperatureFile=dummyString.trimmed();
  inStream >> minTemp >> maxTemp >> tempStep;
  
  QString line("");
  while(line.trimmed()!=QString("</lastRun>")&&!inStream.atEnd()) 
    line=inStream.readLine();
  if(line.trimmed()!=QString("</lastRun>")) return false;

  if(paramMask &  Config::CALCULATE_WITH_DATA) {
    if(paramMask & Config::PERFORM_FIT) runTab->calcType->setCurrentIndex(1);
    else if(paramMask & Config::PERFORM_ERROR_ANALYSIS) runTab->calcType->setCurrentIndex(3);
    else runTab->calcType->setCurrentIndex(0);
  } else {
    if(paramMask &  Config::CALCULATE_REACTION_RATE) runTab->calcType->setCurrentIndex(4);
    else  runTab->calcType->setCurrentIndex(2);
  }

  if(paramMask & Config::USE_GSL_COULOMB_FUNC) GetConfig().paramMask |= Config::USE_GSL_COULOMB_FUNC;
  else GetConfig().paramMask &= ~Config::USE_GSL_COULOMB_FUNC;

  if(paramMask & Config::USE_BRUNE_FORMALISM) GetConfig().paramMask |= Config::USE_BRUNE_FORMALISM;
  else GetConfig().paramMask &= ~Config::USE_BRUNE_FORMALISM;

  if(paramMask & Config::IGNORE_ZERO_WIDTHS) GetConfig().paramMask |= Config::IGNORE_ZERO_WIDTHS;
  else GetConfig().paramMask &= ~Config::IGNORE_ZERO_WIDTHS;
  
  if(paramMask & Config::USE_RMC_FORMALISM) GetConfig().paramMask |= Config::USE_RMC_FORMALISM;
  else GetConfig().paramMask &= ~Config::USE_RMC_FORMALISM;

  if(paramMask & Config::TRANSFORM_PARAMETERS) GetConfig().paramMask |= Config::TRANSFORM_PARAMETERS;
  else GetConfig().paramMask &= ~Config::TRANSFORM_PARAMETERS;

  if(paramMask & Config::USE_LONGWAVELENGTH_APPROX) GetConfig().paramMask |= Config::USE_LONGWAVELENGTH_APPROX;
  else GetConfig().paramMask &= ~Config::USE_LONGWAVELENGTH_APPROX;
  
  if(rateEntrancePair!=0) runTab->rateEntranceKey->setText(QString("%1").arg(rateEntrancePair));
  if(rateExitPair!=0) runTab->rateExitKey->setText(QString("%1").arg(rateExitPair));

  if(minTemp!=-1.) runTab->minTempText->setText(QString("%1").arg(minTemp));
  if(maxTemp!=-1.) runTab->maxTempText->setText(QString("%1").arg(maxTemp));
  if(tempStep!=-1.) runTab->tempStepText->setText(QString("%1").arg(tempStep));
  
  if(useTempFile==1) runTab->fileTempButton->setChecked(true);
  else runTab->gridTempButton->setChecked(true);
  if(temperatureFile[0]==QChar('"')) temperatureFile.remove(0,1);
  if(temperatureFile[temperatureFile.length()-1]==QChar('"')) temperatureFile.remove(temperatureFile.length()-1,1);
  if(!temperatureFile.trimmed().isEmpty()) runTab->fileTempText->setText(temperatureFile.trimmed());
  
  if(paramMask & Config::USE_PREVIOUS_PARAMETERS) runTab->oldParamFileButton->setChecked(true);
  else runTab->newParamFileButton->setChecked(true);
  if(paramFile[0]==QChar('"')) paramFile.remove(0,1);
  if(paramFile[paramFile.length()-1]==QChar('"')) paramFile.remove(paramFile.length()-1,1);
  if(!paramFile.trimmed().isEmpty()) runTab->paramFileText->setText(paramFile.trimmed());

  if(paramMask & Config::USE_PREVIOUS_INTEGRALS) runTab->oldIntegralsFileButton->setChecked(true);
  else runTab->newIntegralsFileButton->setChecked(true);
  if(integralsFile[0]==QChar('"')) integralsFile.remove(0,1);
  if(integralsFile[integralsFile.length()-1]==QChar('"')) integralsFile.remove(integralsFile.length()-1,1);
  if(!integralsFile.trimmed().isEmpty()) runTab->integralsFileText->setText(integralsFile.trimmed());
		      
  return true;
}

bool AZURESetup::readConfig(QTextStream& inStream) {
  
  QString isAMatrix;
  QString outputDirectory;
  QString checksDirectory;
  QString compoundCheck;
  QString boundaryCheck;
  QString dataCheck;
  QString lMatrixCheck;
  QString legendreCheck;
  QString coulAmpCheck;
  QString pathwaysCheck;
  QString angDistsCheck;
  QString dummyString;

  inStream >> isAMatrix;dummyString=inStream.readLine();
  dummyString=inStream.readLine();
  int poundSignPos = dummyString.lastIndexOf('#');
  if(poundSignPos==-1) outputDirectory=dummyString.trimmed();
  else outputDirectory=dummyString.left(poundSignPos).trimmed();
  dummyString=inStream.readLine();
  poundSignPos = dummyString.lastIndexOf('#');
  if(poundSignPos==-1) checksDirectory=dummyString.trimmed();
  else checksDirectory=dummyString.left(poundSignPos).trimmed();
  inStream >> compoundCheck;dummyString=inStream.readLine();
  inStream >> boundaryCheck;dummyString=inStream.readLine();
  inStream >> dataCheck;dummyString=inStream.readLine();
  inStream >> lMatrixCheck;dummyString=inStream.readLine();
  inStream >> legendreCheck;dummyString=inStream.readLine();
  inStream >> coulAmpCheck;dummyString=inStream.readLine();
  inStream >> pathwaysCheck;dummyString=inStream.readLine();
  inStream >> angDistsCheck;dummyString=inStream.readLine();
 
  QString line("");
  while(line.trimmed()!=QString("</config>")&&!inStream.atEnd()) 
    line=inStream.readLine();
  if(line.trimmed()!=QString("</config>")) return false;

  if(isAMatrix=="false") rMatrixAction->activate(QAction::Trigger);
  else aMatrixAction->activate(QAction::Trigger);
  GetConfig().outputdir=outputDirectory.toStdString();
  GetConfig().checkdir=checksDirectory.toStdString();
  if(compoundCheck=="file") GetConfig().fileCheckMask |= Config::CHECK_COMPOUND_NUCLEUS;
  else if(compoundCheck=="screen") GetConfig().screenCheckMask |= Config::CHECK_COMPOUND_NUCLEUS;
  if(boundaryCheck=="file") GetConfig().fileCheckMask |= Config::CHECK_BOUNDARY_CONDITIONS;
  else if(boundaryCheck=="screen")  GetConfig().screenCheckMask |= Config::CHECK_BOUNDARY_CONDITIONS;
  if(dataCheck=="file")  GetConfig().fileCheckMask |= Config::CHECK_DATA;
  else if(dataCheck=="screen") GetConfig().screenCheckMask |= Config::CHECK_DATA;
  if(lMatrixCheck=="file") GetConfig().fileCheckMask |= Config::CHECK_ENERGY_DEP;
  else if(lMatrixCheck=="screen") GetConfig().screenCheckMask |= Config::CHECK_ENERGY_DEP;
  if(legendreCheck=="file") GetConfig().fileCheckMask |= Config::CHECK_LEGENDRE;
  else if(legendreCheck=="screen") GetConfig().screenCheckMask |= Config::CHECK_LEGENDRE;
  if(coulAmpCheck=="file")  GetConfig().fileCheckMask |= Config::CHECK_COUL_AMPLITUDES;
  else if(coulAmpCheck=="screen")  GetConfig().screenCheckMask |= Config::CHECK_COUL_AMPLITUDES;
  if(pathwaysCheck=="file") GetConfig().fileCheckMask |= Config::CHECK_PATHWAYS;
  else if(pathwaysCheck=="screen") GetConfig().screenCheckMask |= Config::CHECK_PATHWAYS;
  if(angDistsCheck=="file") GetConfig().fileCheckMask |= Config::CHECK_ANGULAR_DISTS;
  else if(angDistsCheck=="screen") GetConfig().screenCheckMask |= Config::CHECK_ANGULAR_DISTS;

  return true;
}

void AZURESetup::save() {
  if(!GetConfig().configfile.empty()) {
    if(!this->writeFile(QString::fromStdString(GetConfig().configfile))) 
      QMessageBox::information(this,
			       tr("Can't Access File"),
			       tr("An error occured while writing the file."));
  } else saveAs();
}

void AZURESetup::saveAs() {
  QString filename = QFileDialog::getSaveFileName(this);
  if(!filename.isEmpty()) {
    if(!this->writeFile(filename)) 
      QMessageBox::information(this,
			       tr("Can't Access File"),
			       tr("An error occured while writing the file."));
    else {
      QSettings settings;
      QStringList files = settings.value("recentFileList").toStringList();
      QFile file(filename);
      QFileInfo info(file);
      QString fullFileName = QDir::fromNativeSeparators(info.absoluteFilePath());
      files.removeAll(fullFileName);
      files.prepend(fullFileName);
      while(files.size()>numRecent) files.removeLast();     
      settings.setValue("recentFileList",files);
      updateRecent();
    }
  }
}

bool AZURESetup::writeFile(QString filename) {
  QFile file(filename);
  if(!file.open(QIODevice::WriteOnly)) return false;
  QFileInfo info(file);
  QString directory=info.absolutePath();

  QTextStream out(&file);
  out << "<config>" << endl;
  if(!this->writeConfig(out,directory)) return false;
  out << "</config>" << endl;
  
  out << "<levels>" << endl;
  if(!levelsTab->writeNuclearFile(out)) return false;
  out << "</levels>" << endl;

  out << "<segmentsData>" << endl;
  if(!segmentsTab->writeSegDataFile(out)) return false;
  out << "</segmentsData>" << endl;

  out << "<segmentsTest>" << endl;
  if(!segmentsTab->writeSegTestFile(out)) return false;
  out << "</segmentsTest>" << endl;

  out << "<targetInt>" << endl;
  if(!targetIntTab->writeFile(out)) return false;
  out << "</targetInt>" << endl;  
 
  out << "<lastRun>" << endl;
  if(!writeLastRun(out)) return false;
  out << "</lastRun>" << endl;

  GetConfig().configfile=QDir::fromNativeSeparators(info.absoluteFilePath()).toStdString();
  setWindowTitle(QString("AZURE2 -- %1").arg(QString::fromStdString(GetConfig().configfile)));
  QDir::setCurrent(directory);

  out.flush();
  file.close();
  return true;
}

bool AZURESetup::writeConfig(QTextStream& outStream, QString directory) {
  QString isAMatrix;
  QString outputDirectory;
  QString checksDirectory;
  QString compoundCheck;
  QString boundaryCheck;
  QString dataCheck;
  QString lMatrixCheck;
  QString legendreCheck;
  QString coulAmpCheck;
  QString pathwaysCheck;
  QString angDistsCheck;

  if(GetConfig().paramMask & Config::USE_AMATRIX) isAMatrix="true";
  else isAMatrix="false";
  bool emptyCheckDir=false;
  bool emptyOutputDir=false;
  if(!GetConfig().outputdir.empty())  
    outputDirectory=QString::fromStdString(GetConfig().outputdir);
  else {
    outputDirectory=QDir::fromNativeSeparators(directory)+'/';
    GetConfig().outputdir=outputDirectory.toStdString();
    emptyOutputDir=true;
  }
  if(!GetConfig().checkdir.empty()) 
    checksDirectory=QString::fromStdString(GetConfig().checkdir);
  else {
    checksDirectory=QDir::fromNativeSeparators(directory)+'/';
    GetConfig().checkdir=checksDirectory.toStdString();
    emptyCheckDir=true;
  }
  if(emptyCheckDir&&emptyOutputDir) { 
    QMessageBox::information(this,tr("Unspecified Directories"),
			     QString("The output and checks directories are unspecified. "
				     "They will be set to %1.").arg(outputDirectory.trimmed()));
  } else if(emptyCheckDir) {
    QMessageBox::information(this,tr("Unspecified Directory"),
			     QString("The checks directory is unspecified. "
				     "It will be set to %1.").arg(checksDirectory.trimmed()));
  } else if(emptyOutputDir) {
    QMessageBox::information(this,tr("Unspecified Directory"),
			     QString("The output directory is unspecified. "
				     "It will be set to %1.").arg(outputDirectory.trimmed()));
  }
  if(GetConfig().fileCheckMask & Config::CHECK_COMPOUND_NUCLEUS) compoundCheck="file";
  else if(GetConfig().screenCheckMask & Config::CHECK_COMPOUND_NUCLEUS) compoundCheck="screen";
  else compoundCheck="none";
  if(GetConfig().fileCheckMask & Config::CHECK_BOUNDARY_CONDITIONS) boundaryCheck="file";
  else if(GetConfig().screenCheckMask & Config::CHECK_BOUNDARY_CONDITIONS) boundaryCheck="screen";
  else boundaryCheck="none";
  if(GetConfig().fileCheckMask & Config::CHECK_DATA) dataCheck="file";
  else if(GetConfig().screenCheckMask & Config::CHECK_DATA) dataCheck="screen";
  else dataCheck="none";
  if(GetConfig().fileCheckMask & Config::CHECK_ENERGY_DEP) lMatrixCheck="file";
  else if(GetConfig().screenCheckMask & Config::CHECK_ENERGY_DEP) lMatrixCheck="screen";
  else lMatrixCheck="none";
  if(GetConfig().fileCheckMask & Config::CHECK_LEGENDRE) legendreCheck="file";
  else if(GetConfig().screenCheckMask & Config::CHECK_LEGENDRE) legendreCheck="screen";
  else legendreCheck="none";
  if(GetConfig().fileCheckMask & Config::CHECK_COUL_AMPLITUDES) coulAmpCheck="file";
  else if(GetConfig().screenCheckMask & Config::CHECK_COUL_AMPLITUDES) coulAmpCheck="screen";
  else coulAmpCheck="none";
  if(GetConfig().fileCheckMask & Config::CHECK_PATHWAYS) pathwaysCheck="file";
  else if(GetConfig().screenCheckMask & Config::CHECK_PATHWAYS) pathwaysCheck="screen";
  else pathwaysCheck="none";
  if(GetConfig().fileCheckMask & Config::CHECK_ANGULAR_DISTS) angDistsCheck="file";
  else if(GetConfig().screenCheckMask & Config::CHECK_ANGULAR_DISTS) angDistsCheck="screen";
  else angDistsCheck="none";

  outStream.setFieldAlignment(QTextStream::AlignLeft);
  outStream << qSetFieldWidth(100) << isAMatrix << qSetFieldWidth(0) << "#Perform A-Matrix Calculation" << endl;
  outStream << qSetFieldWidth(100) << outputDirectory << qSetFieldWidth(0) << "#Full Path to Output Directory" << endl;
  outStream << qSetFieldWidth(100) << checksDirectory << qSetFieldWidth(0) << "#Full Path to Checks Directory" << endl;
  outStream << qSetFieldWidth(100) << compoundCheck << qSetFieldWidth(0) << "#Compond Nucleus Check" << endl;
  outStream << qSetFieldWidth(100) << boundaryCheck << qSetFieldWidth(0) << "#Boundary Condition Check" << endl;
  outStream << qSetFieldWidth(100) << dataCheck << qSetFieldWidth(0) << "#Data Check" << endl;
  outStream << qSetFieldWidth(100) << lMatrixCheck << qSetFieldWidth(0) << "#Lo-Matrix and Penetrability Check" << endl;
  outStream << qSetFieldWidth(100) << legendreCheck << qSetFieldWidth(0) << "#Legendre Polynomial Check" << endl;
  outStream << qSetFieldWidth(100) << coulAmpCheck << qSetFieldWidth(0) << "#Coulomb Amplitudes Check" << endl;
  outStream << qSetFieldWidth(100) << pathwaysCheck << qSetFieldWidth(0) << "#Reaction Pathway Check" << endl;
  outStream << qSetFieldWidth(100) << angDistsCheck << qSetFieldWidth(0) << "#Angular Distributions Check" << endl;

  return true;
}

bool AZURESetup::writeLastRun(QTextStream& outStream) {
  unsigned int paramMask = GetConfig().paramMask;

  if(runTab->calcType->currentIndex()==1 ||
     runTab->calcType->currentIndex()==3) paramMask |= Config::PERFORM_FIT;
  else paramMask &= ~Config::PERFORM_FIT;
  if(runTab->calcType->currentIndex()==2||
     runTab->calcType->currentIndex()==4) paramMask &= ~Config::CALCULATE_WITH_DATA;
  else paramMask |= Config::CALCULATE_WITH_DATA;
  if(runTab->calcType->currentIndex()==3) paramMask |= Config::PERFORM_ERROR_ANALYSIS;
  else paramMask &= ~Config::PERFORM_ERROR_ANALYSIS;
  if(runTab->calcType->currentIndex()==4) paramMask |= Config::CALCULATE_REACTION_RATE;
  else paramMask &= ~Config::CALCULATE_REACTION_RATE;

  if(runTab->oldParamFileButton->isChecked())
    paramMask |= Config::USE_PREVIOUS_PARAMETERS;
  else paramMask &= ~Config::USE_PREVIOUS_PARAMETERS;
  if(runTab->oldIntegralsFileButton->isChecked())
    paramMask |= Config::USE_PREVIOUS_INTEGRALS;
  else paramMask &= ~Config::USE_PREVIOUS_INTEGRALS;

  outStream << paramMask << endl;
  outStream << '"' << runTab->paramFileText->text() << '"' << endl;
  outStream << '"' << runTab->integralsFileText->text() << '"' << endl;
  if(!runTab->rateEntranceKey->text().isEmpty()) outStream << runTab->rateEntranceKey->text() << ' ';
  else outStream << "0 "; 
  if(!runTab->rateExitKey->text().isEmpty()) outStream << runTab->rateExitKey->text();
  else outStream << 0; 
  outStream << endl;
  if(runTab->fileTempButton->isChecked()) outStream << "1 "; 
  else outStream << "0 "; 
  outStream << '"' << runTab->fileTempText->text() << '"' << endl;
  if(!runTab->minTempText->text().isEmpty()) outStream << runTab->minTempText->text() << ' ';
  else outStream << "-1. ";
  if(!runTab->maxTempText->text().isEmpty()) outStream << runTab->maxTempText->text() << ' ';
  else outStream << "-1. ";
  if(!runTab->tempStepText->text().isEmpty()) outStream << runTab->tempStepText->text();
  else outStream << "-1.";
  outStream << endl;
	 
  return true;
}

void AZURESetup::matrixChanged(QAction *action) {
  if(action==aMatrixAction) GetConfig().paramMask |= Config::USE_AMATRIX;
  else {
    if(GetConfig().paramMask & Config::IGNORE_ZERO_WIDTHS) {
      QMessageBox::information(this,tr("Incompatible Option"),
			       tr("The option to ignore external widths is not possible for R-Matrix formalism. Remove option to use R-Matrix formalism."));
      aMatrixAction->activate(QAction::Trigger);
    } else GetConfig().paramMask &= ~Config::USE_AMATRIX;
  }
}

void AZURESetup::editChecks() {
  EditChecksDialog aDialog;
  if(GetConfig().fileCheckMask & Config::CHECK_COMPOUND_NUCLEUS) 
    aDialog.compoundCheckCombo->setCurrentIndex(2);
  else if(GetConfig().screenCheckMask & Config::CHECK_COMPOUND_NUCLEUS)
    aDialog.compoundCheckCombo->setCurrentIndex(1);
  else aDialog.compoundCheckCombo->setCurrentIndex(0);

  if(GetConfig().fileCheckMask & Config::CHECK_BOUNDARY_CONDITIONS) 
    aDialog.boundaryCheckCombo->setCurrentIndex(2);
  else if(GetConfig().screenCheckMask & Config::CHECK_BOUNDARY_CONDITIONS)
    aDialog.boundaryCheckCombo->setCurrentIndex(1);
  else aDialog.boundaryCheckCombo->setCurrentIndex(0);

  if(GetConfig().fileCheckMask & Config::CHECK_DATA) 
    aDialog.dataCheckCombo->setCurrentIndex(2);
  else if(GetConfig().screenCheckMask & Config::CHECK_DATA)
    aDialog.dataCheckCombo->setCurrentIndex(1);
  else aDialog.dataCheckCombo->setCurrentIndex(0);

  if(GetConfig().fileCheckMask & Config::CHECK_ENERGY_DEP) 
    aDialog.lMatrixCheckCombo->setCurrentIndex(2);
  else if(GetConfig().screenCheckMask & Config::CHECK_ENERGY_DEP)
    aDialog.lMatrixCheckCombo->setCurrentIndex(1);
  else aDialog.lMatrixCheckCombo->setCurrentIndex(0);

  if(GetConfig().fileCheckMask & Config::CHECK_LEGENDRE) 
    aDialog.legendreCheckCombo->setCurrentIndex(2);
  else if(GetConfig().screenCheckMask & Config::CHECK_LEGENDRE)
    aDialog.legendreCheckCombo->setCurrentIndex(1);
  else aDialog.legendreCheckCombo->setCurrentIndex(0);

  if(GetConfig().fileCheckMask & Config::CHECK_COUL_AMPLITUDES)
    aDialog.coulAmpCheckCombo->setCurrentIndex(2);
  else if(GetConfig().screenCheckMask & Config::CHECK_COUL_AMPLITUDES)
    aDialog.coulAmpCheckCombo->setCurrentIndex(1);
  else aDialog.coulAmpCheckCombo->setCurrentIndex(0);

  if(GetConfig().fileCheckMask & Config::CHECK_PATHWAYS) 
    aDialog.pathwaysCheckCombo->setCurrentIndex(2);
  else if(GetConfig().screenCheckMask & Config::CHECK_PATHWAYS)
    aDialog.pathwaysCheckCombo->setCurrentIndex(1);
  else aDialog.pathwaysCheckCombo->setCurrentIndex(0);

  if(GetConfig().fileCheckMask & Config::CHECK_ANGULAR_DISTS) 
    aDialog.angDistsCheckCombo->setCurrentIndex(2);
  else if(GetConfig().screenCheckMask & Config::CHECK_ANGULAR_DISTS)
    aDialog.angDistsCheckCombo->setCurrentIndex(1);
  else aDialog.angDistsCheckCombo->setCurrentIndex(0);

  if(aDialog.exec()) {
    GetConfig().fileCheckMask=0;
    GetConfig().screenCheckMask=0;
    if(aDialog.compoundCheckCombo->currentIndex()==1) GetConfig().screenCheckMask |= Config::CHECK_COMPOUND_NUCLEUS;
    else if(aDialog.compoundCheckCombo->currentIndex()==2) GetConfig().fileCheckMask |= Config::CHECK_COMPOUND_NUCLEUS;
    if(aDialog.boundaryCheckCombo->currentIndex()==1) GetConfig().screenCheckMask |= Config::CHECK_BOUNDARY_CONDITIONS;
    else if(aDialog.boundaryCheckCombo->currentIndex()==2)  GetConfig().fileCheckMask |= Config::CHECK_BOUNDARY_CONDITIONS;
    if(aDialog.dataCheckCombo->currentIndex()==1) GetConfig().screenCheckMask |= Config::CHECK_DATA;
    else if(aDialog.dataCheckCombo->currentIndex()==2) GetConfig().fileCheckMask |= Config::CHECK_DATA;
    if(aDialog.lMatrixCheckCombo->currentIndex()==1) GetConfig().screenCheckMask |= Config::CHECK_ENERGY_DEP;
    else if(aDialog.lMatrixCheckCombo->currentIndex()==2) GetConfig().fileCheckMask |= Config::CHECK_ENERGY_DEP;
    if(aDialog.legendreCheckCombo->currentIndex()==1) GetConfig().screenCheckMask |= Config::CHECK_LEGENDRE;
    else if(aDialog.legendreCheckCombo->currentIndex()==2) GetConfig().fileCheckMask |= Config::CHECK_LEGENDRE;
    if(aDialog.coulAmpCheckCombo->currentIndex()==1) GetConfig().screenCheckMask |= Config::CHECK_COUL_AMPLITUDES;
    else if(aDialog.coulAmpCheckCombo->currentIndex()==2) GetConfig().fileCheckMask |= Config::CHECK_COUL_AMPLITUDES;
    if(aDialog.pathwaysCheckCombo->currentIndex()==1) GetConfig().screenCheckMask |= Config::CHECK_PATHWAYS;
    else if(aDialog.pathwaysCheckCombo->currentIndex()==2) GetConfig().fileCheckMask |= Config::CHECK_PATHWAYS;
    if(aDialog.angDistsCheckCombo->currentIndex()==1) GetConfig().screenCheckMask |= Config::CHECK_ANGULAR_DISTS;
    else if(aDialog.angDistsCheckCombo->currentIndex()==2) GetConfig().fileCheckMask |= Config::CHECK_ANGULAR_DISTS;
  }
}

void AZURESetup::editDirs() {
  EditDirsDialog aDialog;

  aDialog.outputDirectoryText->setText(QString::fromStdString(GetConfig().outputdir));
  aDialog.checksDirectoryText->setText(QString::fromStdString(GetConfig().checkdir));

  if(aDialog.exec()) {
    GetConfig().outputdir=aDialog.outputDirectoryText->text().toStdString();
    GetConfig().checkdir=aDialog.checksDirectoryText->text().toStdString();
  }
}

void AZURESetup::editOptions() {
  EditOptionsDialog aDialog;

  if(GetConfig().paramMask & Config::USE_GSL_COULOMB_FUNC) aDialog.useGSLCoulCheck->setChecked(true);
  else aDialog.useGSLCoulCheck->setChecked(false);

  if(GetConfig().paramMask & Config::USE_BRUNE_FORMALISM) aDialog.useBruneCheck->setChecked(true);
  else aDialog.useBruneCheck->setChecked(false);
  
  if(GetConfig().paramMask & Config::IGNORE_ZERO_WIDTHS) aDialog.ignoreExternalsCheck->setChecked(true);
  else aDialog.ignoreExternalsCheck->setChecked(false);

  if(GetConfig().paramMask & Config::USE_RMC_FORMALISM) aDialog.useRMCCheck->setChecked(true);
  else aDialog.useRMCCheck->setChecked(false);

  if(!(GetConfig().paramMask & Config::TRANSFORM_PARAMETERS)) aDialog.noTransformCheck->setChecked(true);
  else aDialog.noTransformCheck->setChecked(false);

  if(!(GetConfig().paramMask & Config::USE_LONGWAVELENGTH_APPROX)) aDialog.noLongWavelengthCheck->setChecked(true);
  else aDialog.noLongWavelengthCheck->setChecked(false);

  if(aDialog.exec()) {
    if(aDialog.useGSLCoulCheck->isChecked()) GetConfig().paramMask |= Config::USE_GSL_COULOMB_FUNC;
    else GetConfig().paramMask &= ~Config::USE_GSL_COULOMB_FUNC;

    if(aDialog.useBruneCheck->isChecked()) 
      GetConfig().paramMask |= Config::USE_BRUNE_FORMALISM;
    else GetConfig().paramMask &= ~Config::USE_BRUNE_FORMALISM;

    if(aDialog.ignoreExternalsCheck->isChecked()) {
      GetConfig().paramMask |= Config::IGNORE_ZERO_WIDTHS;
      if(!(GetConfig().paramMask & Config::USE_AMATRIX)) {
	QMessageBox::information(this,tr("Incompatible Option"),
				 tr("The option to ignore external widths is not possible for R-Matrix formalism.  The formalism will be changed to A-Matrix. "));
	aMatrixAction->activate(QAction::Trigger);
      }
    } else  GetConfig().paramMask &= ~Config::IGNORE_ZERO_WIDTHS;
    
    if(aDialog.useRMCCheck->isChecked()) GetConfig().paramMask |= Config::USE_RMC_FORMALISM;
    else GetConfig().paramMask &= ~Config::USE_RMC_FORMALISM;
    
    if(aDialog.noTransformCheck->isChecked()) GetConfig().paramMask &= ~Config::TRANSFORM_PARAMETERS;
    else GetConfig().paramMask |= Config::TRANSFORM_PARAMETERS;

    if(aDialog.noLongWavelengthCheck->isChecked()) GetConfig().paramMask &= ~Config::USE_LONGWAVELENGTH_APPROX;
    else GetConfig().paramMask |= Config::USE_LONGWAVELENGTH_APPROX;
  }
}

void AZURESetup::SaveAndRun() {
  save();
  if(GetConfig().configfile.empty()) return;
  runTab->runtimeText->clear();
  QFile file(QString::fromStdString(GetConfig().configfile));
  QFileInfo info(file);
  QString directory=info.absolutePath();
  if(GetConfig().outputdir.empty()) GetConfig().outputdir=QDir::fromNativeSeparators(directory).toStdString()+'/';
  if(GetConfig().checkdir.empty()) GetConfig().checkdir=QDir::fromNativeSeparators(directory).toStdString()+'/';
  
  GetConfig().chiVariance=runTab->chiVarianceText->text().toDouble();

  if(runTab->calcType->currentIndex()==1 ||
     runTab->calcType->currentIndex()==3) GetConfig().paramMask |= Config::PERFORM_FIT;
  else GetConfig().paramMask &= ~Config::PERFORM_FIT;
  if(runTab->calcType->currentIndex()==2||
     runTab->calcType->currentIndex()==4) GetConfig().paramMask &= ~Config::CALCULATE_WITH_DATA;
  else GetConfig().paramMask |= Config::CALCULATE_WITH_DATA;
  if(runTab->calcType->currentIndex()==3) GetConfig().paramMask |= Config::PERFORM_ERROR_ANALYSIS;
  else GetConfig().paramMask &= ~Config::PERFORM_ERROR_ANALYSIS;
  if(runTab->calcType->currentIndex()==4) GetConfig().paramMask |= Config::CALCULATE_REACTION_RATE;
  else GetConfig().paramMask &= ~Config::CALCULATE_REACTION_RATE;

  if(runTab->oldParamFileButton->isChecked()) {
    GetConfig().paramMask |= Config::USE_PREVIOUS_PARAMETERS;
    GetConfig().paramfile=runTab->paramFileText->text().toStdString();
  } else GetConfig().paramMask &= ~Config::USE_PREVIOUS_PARAMETERS;

  std::vector<SegPairs> segPairs;
  if(!(GetConfig().paramMask & Config::CALCULATE_REACTION_RATE)) {
    if(!readSegmentFile(GetConfig(),segPairs)) return;
  } else {
    GetConfig().rateParams.entrancePair=runTab->rateEntranceKey->text().toInt();
    GetConfig().rateParams.exitPair=runTab->rateExitKey->text().toInt();
    if(GetConfig().rateParams.entrancePair==GetConfig().rateParams.exitPair) {
      QMessageBox::information(this,tr("No Scattering Rates"),
			       tr("Reaction rates cannot be calculated for elastic scattering."));
      return;
    }
    if(!runTab->fileTempButton->isChecked()) {
      GetConfig().rateParams.useFile=false;
      GetConfig().rateParams.minTemp = runTab->minTempText->text().toDouble();
      GetConfig().rateParams.maxTemp = runTab->maxTempText->text().toDouble();
      GetConfig().rateParams.tempStep = runTab->tempStepText->text().toDouble();
    } else {
      GetConfig().rateParams.useFile=true;
      GetConfig().rateParams.temperatureFile = runTab->fileTempText->text().toStdString();
    }
    SegPairs tempPair = {runTab->rateEntranceKey->text().toInt(),
			 runTab->rateExitKey->text().toInt()}; 
    segPairs.push_back(tempPair);
  }
  if(segPairs.size()==0) {
    QMessageBox::information(this,tr("Empty Segments"),tr("No active segments have been found."));
    return;
  }
  int maxPairs=pairsTab->getPairsModel()->getPairs().size();
  for(std::vector<SegPairs>::const_iterator it = segPairs.begin();it<segPairs.end();it++) {
    if(it->secondPair==-1) {
      QList<PairsData> pairsList = pairsTab->getPairsModel()->getPairs();
      int i;
      for(i = 0; i<pairsList.size();i++) 
	if(pairsList[i].pairType==10) break;
      if(i==pairsList.size()) {
	QMessageBox::information(this,tr("No Capture Pairs"),
				 tr("Total capture is specified, but no capture pairs exist."));
	return;
      }
    } else if(it->firstPair>maxPairs||it->secondPair>maxPairs||it->firstPair<1||it->secondPair<1) {
      QMessageBox::information(this,tr("Undefined Key"),tr("An undefined pair key is specified."));
      return;
    }
  }

  GetConfig().paramMask &= ~Config::USE_EXTERNAL_CAPTURE;
  if(!checkExternalCapture(GetConfig(),segPairs)) return;
  if(GetConfig().paramMask &Config::USE_EXTERNAL_CAPTURE) {
    if(runTab->oldIntegralsFileButton->isChecked() && 
       !(GetConfig().paramMask & Config::CALCULATE_REACTION_RATE)) {
      GetConfig().paramMask |= Config::USE_PREVIOUS_INTEGRALS;
      GetConfig().integralsfile=runTab->integralsFileText->text().toStdString();
    } else GetConfig().paramMask &= ~Config::USE_PREVIOUS_INTEGRALS;
  }

  if(!QDir(QString::fromStdString(GetConfig().outputdir)).exists()) {
    QMessageBox::information(this,tr("Directory Doesn't Exist"),
			     tr("The specified output directory doesn't exist."));
    return;
  }
  if(!QDir(QString::fromStdString(GetConfig().checkdir)).exists()) {
    QMessageBox::information(this,tr("Directory Doesn't Exist"),
			     tr("The specified checks directory doesn't exist."));
    return;
  }
  if((GetConfig().paramMask & Config::USE_PREVIOUS_PARAMETERS) &&
     !QFile(QString::fromStdString(GetConfig().paramfile)).exists()) {
    QMessageBox::information(this,tr("File Doesn't Exist"),
			     tr("The specified parameter file doesn't exist."));
    return;
  }
  if(((GetConfig().paramMask & Config::USE_PREVIOUS_INTEGRALS) &&
      (GetConfig().paramMask & Config::USE_EXTERNAL_CAPTURE)) &&
     !QFile(QString::fromStdString(GetConfig().integralsfile)).exists()) {
    QMessageBox::information(this,tr("File Doesn't Exist"),
			     tr("The specified integrals file doesn't exist."));
    return;
  }
  if((GetConfig().paramMask & Config::CALCULATE_REACTION_RATE &&
      GetConfig().rateParams.useFile) &&
     !QFile(QString::fromStdString(GetConfig().rateParams.temperatureFile)).exists()) {
    QMessageBox::information(this,tr("File Doesn't Exist"),
			     tr("The specified rate temperature file doesn't exist."));
    return;
  }

  if(!(GetConfig().paramMask & Config::CALCULATE_WITH_DATA) &&
     !(GetConfig().paramMask & Config::CALCULATE_REACTION_RATE)) {
    QList<TargetIntData> targetIntData = targetIntTab->getTargetIntModel()->getLines();
    QList<SegmentsTestData> segmentsTestData=segmentsTab->getSegmentsTestModel()->getLines();
    for(unsigned int i=0;i<targetIntData.size();i++) {
      if(targetIntData.at(i).isActive==1&&
	 (targetIntData.at(i).isTargetIntegration||targetIntData.at(i).isConvolution)) {
	unsigned int j=0;
	unsigned int lastSegNum=0;
	bool inclusive=false;
	QList<unsigned int> tempList;
	QString segmentsList = targetIntData.at(i).segmentsList;
	while(j<segmentsList.length()) {
	  if(segmentsList[j]>='0'&&segmentsList[j]<='9') {
	    QString tempString;
	    while(segmentsList[j]!=','&&segmentsList[j]!='-'&&
		  j<segmentsList.length()) {
	      tempString+=segmentsList[j];
	      j++;
	    }
	    QTextStream stm(&tempString);
	    unsigned int tempSegNum;stm>>tempSegNum;
	    if(inclusive==true) for(int k=lastSegNum+1;k<=tempSegNum;k++) 
				  tempList.push_back(k);
	    else tempList.push_back(tempSegNum);
	    lastSegNum=tempSegNum;
	  }
	  if(segmentsList[j]=='-') inclusive=true;
	  else inclusive =false;
	  j++;
	}     
	bool isAngularDistribution=false;
	for(j=0;j<tempList.size();j++) {
	  if(tempList.at(j)<=segmentsTestData.size()) {
	    for(int k = 0; k<segmentsTestData.size(); k++) {
	      if(segmentsTestData.at(k).isActive==1&&
		 tempList.at(j)-1==k&&
		 segmentsTestData.at(k).dataType==3) {
		isAngularDistribution=true;
		break;
	      }
	    }
	  }
	  if(isAngularDistribution) break;
	}
	if(isAngularDistribution) {
	  QMessageBox::information(this,tr("Incompatable Options"),
				   tr("Angular distribution coefficients cannot be used with convolution or target integration."));
	  return;	
	}
      }
    }
  }

  azureMain = new AZUREMainThread(runTab,GetConfig());
  connect(azureMain,SIGNAL(finished()),this,SLOT(DeleteThread()));
  setWindowTitle(QString("AZURE2 -- %1 -- Running").arg(QString::fromStdString(GetConfig().configfile)));
  runTab->calcButton->setEnabled(false);
  runTab->stopAZUREButton->setEnabled(true);
  runTab->runtimeText->SetMouseFiltered(true);
  startMessage(azureMain->configure());
  azureMain->start();
}

void AZURESetup::DeleteThread() {
  exitMessage(azureMain->configure());
  QScrollBar *sb = runTab->runtimeText->verticalScrollBar();
  sb->setValue(sb->maximum());

  setWindowTitle(QString("AZURE2 -- %1").arg(QString::fromStdString(GetConfig().configfile)));
  runTab->calcButton->setEnabled(true);
  runTab->stopAZUREButton->setEnabled(false);
  runTab->runtimeText->SetMouseFiltered(false);
  delete azureMain;
}

void AZURESetup::showAbout() {
  AboutAZURE2Dialog aboutDialog;
  aboutDialog.exec();
}

void AZURESetup::reset() {
  GetConfig().Reset();
  aMatrixAction->activate(QAction::Trigger);  
  levelsTab->reset();
  segmentsTab->reset();
  targetIntTab->reset();
  runTab->reset();
#ifdef USE_QWT
  plotTab->reset();
#endif
  setWindowTitle(tr("AZURE2 -- untitled"));
  GetConfig().configfile="";
}

void AZURESetup::showTabInfo() {
  QString tabTitle = tabWidget->tabText(tabWidget->currentIndex()).remove(QChar('&'));
  if(tabWidget->currentIndex()==0) pairsTab->showInfo(0,tabTitle);
  if(tabWidget->currentIndex()==1) levelsTab->showInfo(0,tabTitle);
  if(tabWidget->currentIndex()==2) segmentsTab->showInfo(0,tabTitle);
  if(tabWidget->currentIndex()==3) targetIntTab->showInfo(0,tabTitle);
  if(tabWidget->currentIndex()==4) runTab->showInfo(0,tabTitle);
#ifdef USE_QWT
  if(tabWidget->currentIndex()==5) plotTab->showInfo(0,tabTitle);
#endif
}

void AZURESetup::openWebsite() {
  if(!QDesktopServices::openUrl(QUrl("https://azure.nd.edu")))
    QMessageBox::information(this,
			     tr("Can't Open Browser"),
			     tr("AZURE2 could not access your web browser.  "
				"Please navitgate to https://azure.nd.edu/ "
				"to visit the website."));
}
