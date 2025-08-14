#include <QTabWidget>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QPushButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QLabel>
#include <QHeaderView>
#include <QMessageBox>
#include <QSignalMapper>
#include <QFileDialog>
#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QRegExp>

#include "FittingTab.h"
#include "InfoDialog.h"
#include "LevelsTab.h"
#include "SegmentsTab.h"

FittingTab::FittingTab(QWidget *parent) : QWidget(parent), 
    levelsTab_(nullptr), segmentsTab_(nullptr) {
    
    // Create main layout
    QVBoxLayout *mainLayout = new QVBoxLayout;
    
    // Create parameter tables tab widget
    paramTabWidget = new QTabWidget;
    
    // Create tables for different parameter types  
    levelParamsTable = new QTableWidget;
    setupParameterTable(levelParamsTable, "Level Parameters");
    paramTabWidget->addTab(levelParamsTable, "Level Parameters");
    
    normParamsTable = new QTableWidget;
    setupParameterTable(normParamsTable, "Normalization Parameters");
    paramTabWidget->addTab(normParamsTable, "Normalization");
    
    shiftParamsTable = new QTableWidget;
    setupParameterTable(shiftParamsTable, "Energy Shift Parameters");
    paramTabWidget->addTab(shiftParamsTable, "Energy Shifts");
    
    mainLayout->addWidget(paramTabWidget);
    
    // Create button group
    QHBoxLayout *buttonLayout = new QHBoxLayout;
    
    resetButton = new QPushButton("Reset to Defaults");
    refreshButton = new QPushButton("Refresh from Current");
    loadButton = new QPushButton("Load from .sav file");
    
    buttonLayout->addWidget(resetButton);
    buttonLayout->addWidget(refreshButton);
    buttonLayout->addWidget(loadButton);
    buttonLayout->addStretch();
    
    // Info buttons
    mapper = new QSignalMapper(this);
    for(int i = 0; i < 3; i++) {
        infoButton[i] = new QPushButton("?");
        infoButton[i]->setMaximumSize(20, 20);
        buttonLayout->addWidget(infoButton[i]);
        mapper->setMapping(infoButton[i], i);
        connect(infoButton[i], SIGNAL(clicked()), mapper, SLOT(map()));
    }
    connect(mapper, SIGNAL(mapped(int)), this, SLOT(showInfo(int)));
    
    mainLayout->addLayout(buttonLayout);
    
    // Connect signals
    connect(resetButton, SIGNAL(clicked()), this, SLOT(resetToDefaults()));
    connect(refreshButton, SIGNAL(clicked()), this, SLOT(refreshParameters()));
    connect(loadButton, SIGNAL(clicked()), this, SLOT(loadSettings()));
    
    setLayout(mainLayout);
}

void FittingTab::setupParameterTable(QTableWidget* table, const QString& title) {
    // Set up columns (removed Fixed column - show only unfixed parameters)  
    QStringList headers;
    headers << "Parameter" << "Value" << "Lower Limit" << "Upper Limit" 
            << "Error" << "Use as Nuisance";
    
    table->setColumnCount(headers.size());
    table->setHorizontalHeaderLabels(headers);
    
    // Set column widths
    table->horizontalHeader()->setStretchLastSection(true);
    table->setColumnWidth(0, 150); // Parameter name
    table->setColumnWidth(1, 100); // Value (reduced width for levels)
    table->setColumnWidth(2, 100); // Lower limit
    table->setColumnWidth(3, 100); // Upper limit
    table->setColumnWidth(4, 100); // Error
    table->setColumnWidth(5, 120); // Nuisance checkbox
    
    table->setAlternatingRowColors(true);
    table->setSelectionBehavior(QAbstractItemView::SelectRows);
    
    // Connect item change signal
    connect(table, SIGNAL(itemChanged(QTableWidgetItem*)), 
            this, SLOT(parameterItemChanged(QTableWidgetItem*)));
}

void FittingTab::addParameterRow(QTableWidget* table, const FittingParameter& param) {
    int row = table->rowCount();
    table->insertRow(row);
    
    // Parameter name (read-only)
    QTableWidgetItem* nameItem = new QTableWidgetItem(param.name);
    nameItem->setFlags(nameItem->flags() & ~Qt::ItemIsEditable);
    table->setItem(row, 0, nameItem);
    
    // Value (reduced width for levels)
    table->setItem(row, 1, new QTableWidgetItem(QString::number(param.value, 'g', 6)));
    
    // Lower limit
    table->setItem(row, 2, new QTableWidgetItem(QString::number(param.lowerLimit, 'g', 6)));
    
    // Upper limit
    table->setItem(row, 3, new QTableWidgetItem(QString::number(param.upperLimit, 'g', 6)));
    
    // Error
    table->setItem(row, 4, new QTableWidgetItem(QString::number(param.error, 'g', 6)));
    
    // Nuisance checkbox
    QTableWidgetItem* nuisanceItem = new QTableWidgetItem();
    nuisanceItem->setCheckState(param.useAsNuisance ? Qt::Checked : Qt::Unchecked);
    nuisanceItem->setFlags(nuisanceItem->flags() & ~Qt::ItemIsEditable);
    table->setItem(row, 5, nuisanceItem);
}

void FittingTab::updateParameterTables() {
    // Clear existing tables
    levelParamsTable->setRowCount(0);
    normParamsTable->setRowCount(0);
    shiftParamsTable->setRowCount(0);
    
    // Populate tables with fitting parameters (only non-fixed ones)
    for(const FittingParameter& param : fittingParameters) {
        QTableWidget* targetTable = nullptr;
        
        if(param.category == "level") {
            targetTable = levelParamsTable;
        } else if(param.category == "norm") {
            targetTable = normParamsTable;
        } else if(param.category == "shift") {
            targetTable = shiftParamsTable;
        }
        
        if(targetTable) {
            addParameterRow(targetTable, param);
        }
    }
}

void FittingTab::reset() {
    fittingParameters.clear();
    savedParameterSettings.clear();
    updateParameterTables();
}

void FittingTab::refreshFromMinuitParameters() {
    fittingParameters.clear();
    
    // Create some example parameters to test the interface
    // In a real implementation, this would read from the current tab models
    FittingParameter param1 = {"Level_1_Energy", 0.5, 0.0, 2.0, 0.01, false, "level", 0, 0, -1};
    FittingParameter param2 = {"Level_1_Width_Ch1", 0.1, 0.0, 1.0, 0.005, false, "level", 1, 0, 0};
    FittingParameter param3 = {"Segment_1_Norm", 1.0, 0.5, 2.0, 0.05, true, "norm", 2, -1, -1};
    FittingParameter param4 = {"Segment_1_Shift", 0.0, -0.1, 0.1, 0.001, true, "shift", 3, -1, -1};
    
    fittingParameters.append(param1);
    fittingParameters.append(param2);
    fittingParameters.append(param3);
    fittingParameters.append(param4);
    
    updateParameterTables();
}

void FittingTab::refreshParameters() {
    populateFromCurrentGUIState();
}

void FittingTab::setTabReferences(LevelsTab* levelsTab, SegmentsTab* segmentsTab) {
    levelsTab_ = levelsTab;
    segmentsTab_ = segmentsTab;
    // Populate parameters from current tab state
    populateFromCurrentGUIState();
}

void FittingTab::populateFromCurrentGUIState() {
    if(!levelsTab_ || !segmentsTab_) return;
    
    fittingParameters.clear();
    int paramIndex = 0;
    
    // Get level and channel data
    LevelsModel* levelsModel = levelsTab_->getLevelsModel();
    ChannelsModel* channelsModel = levelsTab_->getChannelsModel();
    
    if(levelsModel && channelsModel) {
        QList<LevelsData> levels = levelsModel->getLevels();
        QList<ChannelsData> channels = channelsModel->getChannels();
        
        // CORRECT ORDER: For each level, add energy first, then all widths for that level
        for(int levelIndex = 0; levelIndex < levels.size(); levelIndex++) {
            const LevelsData& level = levels[levelIndex];
            
            // Add energy parameter if not fixed (isFixed == 0 means not fixed)
            if(level.isFixed == 0) {
                FittingParameter energyParam;
                energyParam.name = QString("Level %1 Energy (MeV)").arg(levelIndex+1);
                energyParam.value = level.energy;
                energyParam.lowerLimit = 0;  // Default limits
                energyParam.upperLimit = 0;
                energyParam.error = 0.01;  // Default error
                energyParam.useAsNuisance = false;
                energyParam.category = "level";
                energyParam.minuitIndex = paramIndex++;
                energyParam.levelIndex = levelIndex;
                energyParam.channelIndex = -1;  // Energy parameter
                
                fittingParameters.append(energyParam);
            }
            
            // Now add all width parameters for this level
            for(int channelIndex = 0; channelIndex < channels.size(); channelIndex++) {
                const ChannelsData& channel = channels[channelIndex];
                
                // Only add widths that belong to this level
                if(channel.levelIndex == levelIndex && channel.isFixed == 0 && channel.reducedWidth != 0.0) {
                    FittingParameter widthParam;
                    widthParam.name = QString("Level %1 Channel %2 Width (eV)")
                                     .arg(levelIndex + 1)
                                     .arg(channelIndex + 1);
                    widthParam.value = channel.reducedWidth;
                    widthParam.lowerLimit = 0;
                    widthParam.upperLimit = 0;  // Default upper limit
                    widthParam.error = channel.reducedWidth * 0.1;  // Default 10% error
                    widthParam.useAsNuisance = false;
                    widthParam.category = "level";
                    widthParam.minuitIndex = paramIndex++;
                    widthParam.levelIndex = levelIndex;
                    widthParam.channelIndex = channelIndex;
                    
                    fittingParameters.append(widthParam);
                }
            }
        }
    }
    
    // Get normalization and shift parameters from SegmentsDataModel
    SegmentsDataModel* segmentsModel = segmentsTab_->getSegmentsDataModel();
    if(segmentsModel) {
        QList<SegmentsDataData> segments = segmentsModel->getLines();
        for(int i = 0; i < segments.size(); i++) {
            const SegmentsDataData& segment = segments[i];
            
            // Add normalization parameter if varied (varyNorm == 1)
            if(segment.varyNorm == 1) {
                FittingParameter normParam;
                // Use data file name if available, otherwise use index
                QString segmentName = segment.dataFile.isEmpty() ? 
                                     QString("Segment %1").arg(i + 1) :
                                     QFileInfo(segment.dataFile).baseName();
                normParam.name = QString("%1 Normalization").arg(segmentName);
                normParam.value = segment.dataNorm;
                normParam.lowerLimit = 0;
                normParam.upperLimit = 0;
                normParam.error = segment.dataNormError;
                normParam.useAsNuisance = true;  // Norms are typically nuisance parameters
                normParam.category = "norm";
                normParam.minuitIndex = paramIndex++;
                normParam.levelIndex = -1;
                normParam.channelIndex = -1;
                
                fittingParameters.append(normParam);
            }
            
            // Add energy shift parameter if varied (varyEnergyShift == 1)
            if(segment.varyEnergyShift == 1) {
                FittingParameter shiftParam;
                // Use data file name if available, otherwise use index
                QString segmentName = segment.dataFile.isEmpty() ? 
                                     QString("Segment %1").arg(i + 1) :
                                     QFileInfo(segment.dataFile).baseName();
                shiftParam.name = QString("%1 Energy Shift (keV)").arg(segmentName);
                shiftParam.value = segment.energyShift;
                shiftParam.lowerLimit = 0;
                shiftParam.upperLimit = 0;
                shiftParam.error = segment.energyShiftError;
                shiftParam.useAsNuisance = true;  // Shifts are typically nuisance parameters
                shiftParam.category = "shift";
                shiftParam.minuitIndex = paramIndex++;
                shiftParam.levelIndex = -1;
                shiftParam.channelIndex = -1;
                
                fittingParameters.append(shiftParam);
            }
        }
    }
    
    // Apply parameter settings from saved configuration (limits, errors, etc.)
    applyParameterSettings();
    
    updateParameterTables();
}

double FittingTab::convertReducedToPhysical(double reducedWidth, int levelIndex, int channelIndex) {
    // Width conversion requires access to the compound nucleus and pair data
    // which is not available at the GUI level. The actual conversion happens
    // during calculations when the CNuc object is created.
    // For display purposes, we return a simplified conversion.
    
    // This is a placeholder conversion - the real calculation is complex
    // and involves penetrability factors, channel radii, etc.
    return reducedWidth * 100.0; // Simple scaling for display
}

double FittingTab::convertPhysicalToReduced(double physicalWidth, int levelIndex, int channelIndex) {
    // Inverse of the above conversion for display purposes
    return physicalWidth / 100.0;
}

void FittingTab::applyParameterSettings() {
    // Apply saved parameter settings (limits, errors, etc.) to current parameters
    for(int i = 0; i < fittingParameters.size(); i++) {
        for(const FittingParameter& saved : savedParameterSettings) {
            if(fittingParameters[i].name == saved.name) {
                // Apply saved settings but keep current value from models
                fittingParameters[i].lowerLimit = saved.lowerLimit;
                fittingParameters[i].upperLimit = saved.upperLimit;
                fittingParameters[i].error = saved.error;
                fittingParameters[i].useAsNuisance = saved.useAsNuisance;
                break;
            }
        }
    }
}

void FittingTab::parameterItemChanged(QTableWidgetItem* item) {
    // Handle parameter changes
    int row = item->row();
    int col = item->column();
    QTableWidget* table = qobject_cast<QTableWidget*>(item->tableWidget());
    if(!table) return;
    
    // Get parameter name from first column
    QTableWidgetItem* nameItem = table->item(row, 0);
    if(!nameItem) return;
    QString paramName = nameItem->text();
    
    // Find the parameter in our settings
    int paramIndex = -1;
    for(int i = 0; i < fittingParameters.size(); i++) {
        if(fittingParameters[i].name == paramName) {
            paramIndex = i;
            break;
        }
    }
    
    if(paramIndex == -1) return; // Parameter not found
    
    // Column mapping: 0=Name, 1=Value, 2=Lower, 3=Upper, 4=Error, 5=Nuisance
    if(col == 1) { // Value (reduced width) changed
        bool ok;
        double value = item->text().toDouble(&ok);
        if(!ok) {
            QMessageBox::warning(this, "Invalid Input", "Please enter a valid number.");
            item->setText(QString::number(fittingParameters[paramIndex].value, 'g', 6));
            return;
        }
        
        fittingParameters[paramIndex].value = value;
        
        // Update the corresponding tab with the new value
        updateParameterInOtherTabs(paramName, fittingParameters[paramIndex]);
        
    } else if(col >= 2 && col <= 4) { // Lower limit, upper limit, or error changed
        bool ok;
        double value = item->text().toDouble(&ok);
        if(!ok) {
            QMessageBox::warning(this, "Invalid Input", "Please enter a valid number.");
            // Reset to previous value
            if(col == 2) item->setText(QString::number(fittingParameters[paramIndex].lowerLimit, 'g', 6));
            else if(col == 3) item->setText(QString::number(fittingParameters[paramIndex].upperLimit, 'g', 6));
            else if(col == 4) item->setText(QString::number(fittingParameters[paramIndex].error, 'g', 6));
            return;
        }
        
        // Update parameter setting
        if(col == 2) fittingParameters[paramIndex].lowerLimit = value;
        else if(col == 3) fittingParameters[paramIndex].upperLimit = value;
        else if(col == 4) fittingParameters[paramIndex].error = value;
        
    } else if(col == 5) { // Nuisance checkbox
        fittingParameters[paramIndex].useAsNuisance = (item->checkState() == Qt::Checked);
    }
}

void FittingTab::resetToDefaults() {
    int ret = QMessageBox::question(this, "Reset Parameters",
                                   "Are you sure you want to reset all parameter settings to defaults?",
                                   QMessageBox::Yes | QMessageBox::No);
    if(ret == QMessageBox::Yes) {
        updateParameterTables();
    }
}



void FittingTab::loadSettings() {
    QString filename = QFileDialog::getOpenFileName(this, 
        "Load Parameters from .sav file", "", "AZURE2 Parameter Files (*.sav);;All Files (*)");
    
    if(!filename.isEmpty()) {
        QFile file(filename);
        if(file.open(QIODevice::ReadOnly | QIODevice::Text)) {
            QTextStream in(&file);
            
            // Parse .sav file format: name value error
            fittingParameters.clear();
            int paramIndex = 0;
            
            while(!in.atEnd()) {
                QString line = in.readLine().trimmed();
                if(line.isEmpty()) continue;
                
                QStringList parts = line.split(QRegExp("\\s+"), Qt::SkipEmptyParts);
                if(parts.size() >= 3) {
                    FittingParameter param;
                    param.name = parts[0];
                    param.value = parts[1].toDouble();
                    param.error = parts[2].toDouble();
                    param.minuitIndex = paramIndex++;

                    // If "_rwa" in param.name, remove it
                    if(param.name.endsWith("_rwa")) {
                        param.name.chop(4); // Remove last 4 characters
                    }
                    
                    // Set default limits (can be adjusted by user)
                    param.lowerLimit = param.value - 5 * param.error;
                    param.upperLimit = param.value + 5 * param.error;
                    param.useAsNuisance = false;
                    
                    // Determine category from parameter name
                    if(param.name.contains("segment") && param.name.contains("norm")) {
                        param.category = "norm";
                        param.useAsNuisance = true; // Norms are typically nuisance parameters
                        param.levelIndex = -1;
                        param.channelIndex = -1;
                    } else if(param.name.contains("segment") && param.name.contains("shift")) {
                        param.category = "shift";
                        param.useAsNuisance = true; // Shifts are typically nuisance parameters
                        param.levelIndex = -1;
                        param.channelIndex = -1;
                    } else {
                        param.category = "level";
                        param.levelIndex = 0; // TODO: Parse from name
                        param.channelIndex = 0; // TODO: Parse from name
                        
                        // Level width parameters are already in reduced width units
                    }
                    
                    fittingParameters.append(param);
                }
            }
            
            file.close();
            updateParameterTables();
            QMessageBox::information(this, "Load Settings", 
                                   QString("Loaded %1 parameters from: %2").arg(fittingParameters.size()).arg(filename));
        } else {
            QMessageBox::warning(this, "Load Error", 
                               "Could not load parameter file.");
        }
    }
}

void FittingTab::updateParameterInOtherTabs(const QString& paramName, const FittingParameter& param) {
    if(!levelsTab_ || !segmentsTab_) return;
    
    if(param.category == "level") {
        // Update level parameters in LevelsModel - same pattern as LevelsTab::editLevel()
        LevelsModel* levelsModel = levelsTab_->getLevelsModel();
        if(levelsModel) {
            // Parse level index from parameter name
            int levelIndex = -1;
            int channelIndex = -1;
            
            if(paramName.contains("_energy")) {
                // Energy parameter: j=%d_la=%d_energy
                QRegExp rx("j=\\d+_la=(\\d+)_energy");
                if(rx.indexIn(paramName) != -1) {
                    levelIndex = rx.cap(1).toInt() - 1; // Convert 1-based to 0-based
                }
            } else if(paramName.contains("_ch=")) {
                // Width parameter: j=%d_la=%d_ch=%d
                QRegExp rx("j=\\d+_la=(\\d+)_ch=(\\d+)");
                if(rx.indexIn(paramName) != -1) {
                    levelIndex = rx.cap(1).toInt() - 1; // Convert 1-based to 0-based
                    channelIndex = rx.cap(2).toInt() - 1; // Convert 1-based to 0-based
                }
            }
            
            if(levelIndex >= 0) {
                QList<LevelsData> levels = levelsModel->getLevels();
                if(levelIndex < levels.size()) {
                    if(channelIndex == -1) {
                        // Energy parameter - update column 4 (energy column)
                        QModelIndex index = levelsModel->index(levelIndex, 4);
                        levelsModel->setData(index, param.value, Qt::EditRole);
                    } else {
                        // Width parameter - update channels model
                        ChannelsModel* channelsModel = levelsTab_->getChannelsModel();
                        if(channelsModel) {
                            QList<ChannelsData> channels = channelsModel->getChannels();
                            // Find the channel with matching level and channel indices
                            for(int i = 0; i < channels.size(); i++) {
                                if(channels[i].levelIndex == levelIndex && i == channelIndex) {
                                    // Update reducedWidth column (column 6)
                                    QModelIndex index = channelsModel->index(i, 6);
                                    channelsModel->setData(index, param.value, Qt::EditRole);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    } else if(param.category == "norm" || param.category == "shift") {
        // Update normalization or shift parameters in SegmentsDataModel
        SegmentsDataModel* segmentsModel = segmentsTab_->getSegmentsDataModel();
        if(segmentsModel) {
            QList<SegmentsDataData> segments = segmentsModel->getLines();
            
            // Parse segment index from parameter name (e.g., "segment_1_norm" -> index = 0)
            QRegExp rx("segment_(\\d+)_");
            if(rx.indexIn(paramName) != -1) {
                int segmentIndex = rx.cap(1).toInt() - 1; // Convert 1-based to 0-based
                
                if(segmentIndex >= 0 && segmentIndex < segments.size()) {
                    QModelIndex index;
                    
                    if(param.category == "norm") {
                        // Update dataNorm column (column 9)
                        index = segmentsModel->index(segmentIndex, 9);
                        segmentsModel->setData(index, param.value, Qt::EditRole);
                    } else if(param.category == "shift") {
                        // Update energyShift column (column 14)
                        index = segmentsModel->index(segmentIndex, 14);
                        segmentsModel->setData(index, param.value, Qt::EditRole);
                    }
                }
            }
        }
    }
}

void FittingTab::showInfo(int which, QString title) {
    if(!infoDialog[which]) {
        // QString construct info dialog
        infoDialog[which] = new InfoDialog( "", this, title);
        if(title.isEmpty()) {
            switch(which) {
                case 0: infoDialog[which]->setWindowTitle("Level Parameters Help"); break;
                case 1: infoDialog[which]->setWindowTitle("Normalization Parameters Help"); break;
                case 2: infoDialog[which]->setWindowTitle("Energy Shift Parameters Help"); break;
            }
        } else {
            infoDialog[which]->setWindowTitle(title);
        }
        //infoDialog[which]->setInfoText(infoText[which]);
    }
    infoDialog[which]->show();
    infoDialog[which]->raise();
}

bool FittingTab::writeParameterSettings(QTextStream& outStream) {
    // Write current parameter settings to AZURE2 file
    outStream << "# Fitting parameter settings (only non-fixed parameters shown)\n";
    outStream << "# Format: name value lower_limit upper_limit error nuisance category minuit_index\n";
    
    for(const FittingParameter& param : fittingParameters) {
        outStream << param.name << " " 
                  << param.value << " "
                  << param.lowerLimit << " "
                  << param.upperLimit << " "
                  << param.error << " "
                  << (param.useAsNuisance ? 1 : 0) << " "
                  << param.category << " "
                  << param.minuitIndex << "\n";
    }
    return true;
}

bool FittingTab::readParameterSettings(QTextStream& inStream) {
    savedParameterSettings.clear();
    
    QString line;
    while(!inStream.atEnd()) {
        line = inStream.readLine().trimmed();
        if(line.startsWith("<")) break; // Next section started
        if(line.isEmpty() || line.startsWith("#")) continue; // Skip empty lines and comments
        
        QStringList parts = line.split(" ", Qt::SkipEmptyParts);
        if(parts.size() >= 7) {
            FittingParameter param;
            param.name = parts[0];
            param.value = parts[1].toDouble(); // This will be overridden by current model values
            param.lowerLimit = parts[2].toDouble();
            param.upperLimit = parts[3].toDouble();
            param.error = parts[4].toDouble();
            param.useAsNuisance = (parts[5].toInt() == 1);
            param.category = parts[6];
            param.minuitIndex = (parts.size() >= 8) ? parts[7].toInt() : -1;
            
            // Set defaults for level-specific parameters
            param.levelIndex = -1;
            param.channelIndex = -1;
            
            savedParameterSettings.append(param);
        }
    }
    
    return true;
}

// Static info text - this would be defined in InTabDocs.cpp in the real implementation
const std::vector<QString> FittingTab::infoText = {
    QString("Level parameters control the R-matrix level energies and widths. "
           "Set limits to constrain parameter values during fitting. "
           "Enable 'Use as Nuisance' to include parameter uncertainty in chi-squared."),
    QString("Normalization parameters adjust the overall scale of data segments. "
           "These are typically varied during fitting to account for experimental uncertainties."),
    QString("Energy shift parameters correct for energy calibration offsets in data segments. "
           "These parameters shift the energy scale of experimental data points.")
};