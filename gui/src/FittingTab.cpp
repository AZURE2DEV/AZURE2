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
#include <cmath>
#include <algorithm>

#include "FittingTab.h"
#include "InfoDialog.h"
#include "LevelsTab.h"
#include "SegmentsTab.h"
#include "AZURESetup.h"

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
    
    refreshButton = new QPushButton("Refresh from Current");
    loadButton = new QPushButton("Load from .sav file");
    
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
    
    // Connect signals from SegmentsDataModel to update FittingTab when values change
    if(segmentsTab_) {
        SegmentsDataModel* segmentsModel = segmentsTab_->getSegmentsDataModel();
        if(segmentsModel) {
            connect(segmentsModel, SIGNAL(normalizationChanged(int, double)),
                    this, SLOT(onSegmentNormalizationChanged(int, double)));
            connect(segmentsModel, SIGNAL(energyShiftChanged(int, double)),
                    this, SLOT(onSegmentEnergyShiftChanged(int, double)));
            connect(segmentsModel, SIGNAL(normalizationErrorChanged(int, double)),
                    this, SLOT(onSegmentNormalizationErrorChanged(int, double)));
            connect(segmentsModel, SIGNAL(energyShiftErrorChanged(int, double)),
                    this, SLOT(onSegmentEnergyShiftErrorChanged(int, double)));
            connect(segmentsModel, SIGNAL(normalizationVaryChanged(int, bool)),
                    this, SLOT(onSegmentNormalizationVaryChanged(int, bool)));
            connect(segmentsModel, SIGNAL(energyShiftVaryChanged(int, bool)),
                    this, SLOT(onSegmentEnergyShiftVaryChanged(int, bool)));
        }
    }
    
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
                    
                    // LevelsTab already stores physical widths, use them directly
                    widthParam.value = channel.reducedWidth; // This is actually physical width in the GUI
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
            
            // Add normalization parameter (always show, but useAsNuisance depends on varyNorm flag)
            FittingParameter normParam;
            // Use same naming pattern as MCMC tab
            QString segmentName = segment.dataFile.isEmpty() ? 
                                 QString("Segment %1").arg(i + 1) :
                                 QFileInfo(segment.dataFile).baseName();
            normParam.name = QString("%1 Normalization").arg(segmentName);
            normParam.value = segment.dataNorm;
            normParam.lowerLimit = 0;
            normParam.upperLimit = 0;
            normParam.error = segment.dataNormError;
            normParam.useAsNuisance = (segment.varyNorm == 1);  // Active only when parameter varies
            normParam.category = "norm";
            normParam.minuitIndex = paramIndex++;
            normParam.levelIndex = -1;
            normParam.channelIndex = i;  // Store segment index for reverse lookup
            
            fittingParameters.append(normParam);
            
            // Add energy shift parameter (always show, but useAsNuisance depends on varyEnergyShift flag)
            FittingParameter shiftParam;
            shiftParam.name = QString("%1 Energy Shift (keV)").arg(segmentName);
            shiftParam.value = segment.energyShift;
            shiftParam.lowerLimit = 0;
            shiftParam.upperLimit = 0;
            shiftParam.error = segment.energyShiftError;
            shiftParam.useAsNuisance = (segment.varyEnergyShift == 1);  // Active only when parameter varies
            shiftParam.category = "shift";
            shiftParam.minuitIndex = paramIndex++;
            shiftParam.levelIndex = -1;
            shiftParam.channelIndex = i;  // Store segment index for reverse lookup
            
            fittingParameters.append(shiftParam);
        }
    }
    
    // Apply parameter settings from saved configuration (limits, errors, etc.)
    applyParameterSettings();
    
    updateParameterTables();
}

double FittingTab::transformRWAParameterToPhysical(const QString& paramName, double rwaValue) {
    // Transform RWA parameter to physical using proper R-Matrix transformation via AZURESetup
    
    // Find the parent AZURESetup widget
    AZURESetup* azureSetup = nullptr;
    QWidget* parent = this->parentWidget();
    while(parent != nullptr) {
        azureSetup = qobject_cast<AZURESetup*>(parent);
        if(azureSetup != nullptr) {
            break;
        }
        parent = parent->parentWidget();
    }
    
    if(azureSetup != nullptr) {
        // Use the proper RWA to Physical conversion from AZURESetup
        return azureSetup->ConvertRWAToPhysical(paramName, rwaValue);
    } else {
        // Fallback: return the original value if we can't find AZURESetup
        return rwaValue;
    }
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
        else if(col == 4) {
            fittingParameters[paramIndex].error = value;
            // Update the corresponding tab with the new error value
            updateParameterInOtherTabs(paramName, fittingParameters[paramIndex]);
        }
        
    } else if(col == 5) { // Nuisance checkbox
        fittingParameters[paramIndex].useAsNuisance = (item->checkState() == Qt::Checked);
        // Update the corresponding tab with the new useAsNuisance value
        updateParameterInOtherTabs(paramName, fittingParameters[paramIndex]);
    }
}

void FittingTab::loadSettings() {
    QString filename = QFileDialog::getOpenFileName(this, 
        "Load Parameters from .sav file", "", "AZURE2 Parameter Files (*.sav);;All Files (*)");
    
    if(!filename.isEmpty()) {
        QFile file(filename);
        if(file.open(QIODevice::ReadOnly | QIODevice::Text)) {
            QTextStream in(&file);
            
            // Parse .sav file and create lookup map
            QMap<QString, QPair<double, double>> savParams; // paramName -> (value, error)
            
            while(!in.atEnd()) {
                QString line = in.readLine().trimmed();
                if(line.isEmpty()) continue;
                
                QStringList parts = line.split(QRegExp("\\s+"), Qt::SkipEmptyParts);
                if(parts.size() >= 3) {
                    QString paramName = parts[0];
                    double value = parts[1].toDouble();
                    double error = parts[2].toDouble();
                    
                    // If "_rwa" in param.name, remove it
                    if(paramName.endsWith("_rwa")) {
                        paramName.chop(4); // Remove last 4 characters
                    }
                    
                    savParams[paramName] = qMakePair(value, error);
                }
            }
            file.close();
            
            // Create lists of parameters and values to convert in batch
            QStringList paramNamesToConvert;
            QList<double> rwaValuesToConvert;
            QList<QPair<int, QString>> parameterMapping; // (fittingParameters index, matchKey)
            
            // First pass: collect all parameters that need conversion
            for(int i = 0; i < fittingParameters.size(); i++) {
                FittingParameter& param = fittingParameters[i];
                QString matchKey = findMatchingParameterKey(param, savParams.keys());
                
                if(!matchKey.isEmpty() && savParams.contains(matchKey)) {
                    QPair<double, double> savData = savParams[matchKey];
                    double rwaValue = savData.first;
                    double rwaError = savData.second;
                    
                    parameterMapping.append(qMakePair(i, matchKey));
                    
                    // Add main value
                    paramNamesToConvert.append(matchKey);
                    rwaValuesToConvert.append(rwaValue);
                    
                    // For error calculation, add +/- error values for level parameters
                    if(param.category == "level") {
                        paramNamesToConvert.append(matchKey + "_plus");
                        rwaValuesToConvert.append(rwaValue + rwaError);
                        paramNamesToConvert.append(matchKey + "_minus");
                        rwaValuesToConvert.append(rwaValue - rwaError);
                    }
                }
            }
            
            // Batch conversion of all RWA parameters to physical
            QList<double> convertedValues;
            if(!paramNamesToConvert.isEmpty()) {
                // Find the parent AZURESetup widget for batch conversion
                AZURESetup* azureSetup = nullptr;
                QWidget* parent = this->parentWidget();
                while(parent != nullptr) {
                    azureSetup = qobject_cast<AZURESetup*>(parent);
                    if(azureSetup != nullptr) break;
                    parent = parent->parentWidget();
                }
                
                if(azureSetup != nullptr) {
                    // Batch convert all parameters at once
                    for(int i = 0; i < paramNamesToConvert.size(); i++) {
                        QString paramName = paramNamesToConvert[i];
                        // Remove suffix markers for error calculation
                        if(paramName.endsWith("_plus") || paramName.endsWith("_minus")) {
                            paramName = paramName.left(paramName.lastIndexOf("_"));
                        }
                        double physicalValue = azureSetup->ConvertRWAToPhysical(paramName, rwaValuesToConvert[i]);
                        convertedValues.append(physicalValue);
                    }
                } else {
                    // Fallback: use original values if conversion unavailable
                    convertedValues = rwaValuesToConvert;
                }
            }
            
            // Second pass: apply converted values to parameters
            int updatedCount = 0;
            int valueIndex = 0;
            
            for(const QPair<int, QString>& mapping : parameterMapping) {
                int paramIndex = mapping.first;
                QString matchKey = mapping.second;
                FittingParameter& param = fittingParameters[paramIndex];
                QPair<double, double> savData = savParams[matchKey];
                
                if(param.category == "level" && (param.name.contains("Width") || param.name.contains("Energy"))) {
                    // Use converted values with error calculation
                    if(valueIndex + 2 < convertedValues.size()) {
                        param.value = convertedValues[valueIndex];
                        double physicalValuePlusError = convertedValues[valueIndex + 1];
                        double physicalValueMinusError = convertedValues[valueIndex + 2];
                        
                        // Calculate error from the converted values
                        double errorUp = std::abs(physicalValuePlusError - param.value);
                        double errorDown = std::abs(param.value - physicalValueMinusError);
                        param.error = std::max(errorUp, errorDown);
                        
                        valueIndex += 3; // Move past the three values (main, plus, minus)
                    }
                } else {
                    // For other parameters, use direct conversion
                    if(valueIndex < convertedValues.size()) {
                        param.value = convertedValues[valueIndex];
                        param.error = savData.second; // Keep original error for non-level parameters
                        valueIndex += 1;
                    }
                }
                
                // Update the underlying models with converted values
                if(param.category == "level" && param.name.contains("Width") && param.channelIndex >= 0) {
                    if(levelsTab_) {
                        ChannelsModel* channelsModel = levelsTab_->getChannelsModel();
                        if(channelsModel) {
                            QList<ChannelsData> channels = channelsModel->getChannels();
                            for(int j = 0; j < channels.size(); j++) {
                                if(j == param.channelIndex) {
                                    QModelIndex modelIndex = channelsModel->index(j, 6);
                                    channelsModel->setData(modelIndex, param.value, Qt::EditRole);
                                    break;
                                }
                            }
                        }
                    }
                }
                
                // Propagate changes to other tabs
                updateParameterInOtherTabs(param.name, param);
                updatedCount++;
            }
            
            // Refresh the parameter tables
            updateParameterTables();
            
            QMessageBox::information(this, "Load Settings", 
                                   QString("Updated %1 of %2 fitting parameters from: %3")
                                   .arg(updatedCount).arg(fittingParameters.size()).arg(filename));
        } else {
            QMessageBox::warning(this, "Load Error", 
                               "Could not load parameter file.");
        }
    }
}

QString FittingTab::findMatchingParameterKey(const FittingParameter& param, const QStringList& savKeys) {
    // Try to match current parameter with .sav file parameter names
    
    // First try direct match
    if(savKeys.contains(param.name)) {
        return param.name;
    }
    
    if(param.category == "norm") {
        // For normalization parameters: look for patterns like "segment_1_norm", "norm_1", etc.
        int segmentIndex = param.channelIndex + 1; // Convert 0-based to 1-based
        
        QStringList possibleNames;
        possibleNames << QString("segment_%1_norm").arg(segmentIndex);
        possibleNames << QString("norm_%1").arg(segmentIndex);
        possibleNames << QString("norm%1").arg(segmentIndex);
        possibleNames << QString("segment%1_norm").arg(segmentIndex);
        
        for(const QString& name : possibleNames) {
            if(savKeys.contains(name)) {
                return name;
            }
        }
        
    } else if(param.category == "shift") {
        // For energy shift parameters: look for patterns like "segment_1_energy_shift", "shift_1", etc.
        int segmentIndex = param.channelIndex + 1; // Convert 0-based to 1-based
        
        QStringList possibleNames;
        possibleNames << QString("segment_%1_energy_shift").arg(segmentIndex);  // Main pattern from .sav files
        possibleNames << QString("segment_%1_shift").arg(segmentIndex);
        possibleNames << QString("shift_%1").arg(segmentIndex);
        possibleNames << QString("shift%1").arg(segmentIndex);
        possibleNames << QString("segment%1_energy_shift").arg(segmentIndex);
        possibleNames << QString("segment%1_shift").arg(segmentIndex);
        
        for(const QString& name : possibleNames) {
            if(savKeys.contains(name)) {
                return name;
            }
        }
        
    } else if(param.category == "level") {
        // For level parameters: match with .sav file patterns like "energy_1", "width_1_2"
        // Current GUI names are like "Level 1 Energy (MeV)" and "Level 1 Channel 2 Width (eV)"
        
        if(param.name.contains("Energy") && param.channelIndex == -1) {
            // Energy parameter: "Level N Energy (MeV)" -> "energy_N"
            int levelIndex = param.levelIndex + 1; // Convert 0-based to 1-based
            
            QStringList possibleNames;
            possibleNames << QString("energy_%1").arg(levelIndex);
            possibleNames << QString("energy%1").arg(levelIndex);
            possibleNames << QString("level_%1_energy").arg(levelIndex);
            possibleNames << QString("level%1_energy").arg(levelIndex);
            
            for(const QString& name : possibleNames) {
                if(savKeys.contains(name)) {
                    return name;
                }
            }
            
        } else if(param.name.contains("Width") && param.channelIndex >= 0) {
            // Width parameter: "Level N Channel M Width (eV)" -> "width_N_M"
            // Need to count channels per level, not global channel index
            int levelIndex = param.levelIndex + 1; // Convert 0-based to 1-based
            
            // Count which channel this is within this specific level
            int channelWithinLevel = 1; // Start counting from 1 for .sav file format
            
            if(levelsTab_ && segmentsTab_) {
                ChannelsModel* channelsModel = levelsTab_->getChannelsModel();
                if(channelsModel) {
                    QList<ChannelsData> channels = channelsModel->getChannels();
                    for(int i = 0; i < channels.size(); i++) {
                        const ChannelsData& channel = channels[i];
                        if(channel.levelIndex == param.levelIndex) {
                            if(i == param.channelIndex) {
                                // Found our channel - channelWithinLevel is correct
                                break;
                            }
                            channelWithinLevel++; // Count channels for this level
                        }
                    }
                }
            }
            
            QStringList possibleNames;
            possibleNames << QString("width_%1_%2").arg(levelIndex).arg(channelWithinLevel);
            possibleNames << QString("width%1_%2").arg(levelIndex).arg(channelWithinLevel);
            possibleNames << QString("level_%1_channel_%2_width").arg(levelIndex).arg(channelWithinLevel);
            possibleNames << QString("level%1_channel%2_width").arg(levelIndex).arg(channelWithinLevel);
            
            for(const QString& name : possibleNames) {
                if(savKeys.contains(name)) {
                    return name;
                }
            }
        }
    }
    
    // No match found
    return QString();
}

void FittingTab::updateParameterInOtherTabs(const QString& paramName, const FittingParameter& param) {
    if(!levelsTab_ || !segmentsTab_) return;
    
    if(param.category == "level") {
        // Update level parameters in LevelsModel - same pattern as LevelsTab::editLevel()
        LevelsModel* levelsModel = levelsTab_->getLevelsModel();
        if(levelsModel) {
            // Parse level index from parameter name using current naming scheme
            // Parameter names are like "Level 1 Energy (MeV)" and "Level 1 Channel 2 Width (eV)"
            int levelIndex = -1;
            int channelIndex = -1;
            
            if(paramName.contains("Energy") && paramName.contains("Level")) {
                // Energy parameter: "Level N Energy (MeV)"
                QRegExp rx("Level (\\d+) Energy");
                if(rx.indexIn(paramName) != -1) {
                    levelIndex = rx.cap(1).toInt() - 1; // Convert 1-based to 0-based
                }
            } else if(paramName.contains("Width") && paramName.contains("Level") && paramName.contains("Channel")) {
                // Width parameter: "Level N Channel M Width (eV)"
                QRegExp rx("Level (\\d+) Channel (\\d+) Width");
                if(rx.indexIn(paramName) != -1) {
                    levelIndex = rx.cap(1).toInt() - 1; // Convert 1-based to 0-based
                    channelIndex = rx.cap(2).toInt() - 1; // Convert 1-based to 0-based
                }
            }
            
            // Alternative: use the stored indices from the FittingParameter structure
            if(levelIndex == -1 && param.levelIndex >= 0) {
                levelIndex = param.levelIndex;
                channelIndex = param.channelIndex;
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
            
            // Use channelIndex which now stores the segment index
            int segmentIndex = param.channelIndex;
            
            if(segmentIndex >= 0 && segmentIndex < segments.size()) {
                QModelIndex index;
                
                if(param.category == "norm") {
                    // Update dataNorm column (column 9)
                    index = segmentsModel->index(segmentIndex, 9);
                    segmentsModel->setData(index, param.value, Qt::EditRole);
                    // Also update dataNormError column (column 10)
                    QModelIndex errorIndex = segmentsModel->index(segmentIndex, 10);
                    segmentsModel->setData(errorIndex, param.error, Qt::EditRole);
                    // Update varyNorm column (column 11) based on useAsNuisance
                    QModelIndex varyIndex = segmentsModel->index(segmentIndex, 11);
                    segmentsModel->setData(varyIndex, param.useAsNuisance ? 1 : 0, Qt::EditRole);
                } else if(param.category == "shift") {
                    // Update energyShift column (column 14)
                    index = segmentsModel->index(segmentIndex, 14);
                    segmentsModel->setData(index, param.value, Qt::EditRole);
                    // Also update energyShiftError column (column 15)
                    QModelIndex errorIndex = segmentsModel->index(segmentIndex, 15);
                    segmentsModel->setData(errorIndex, param.error, Qt::EditRole);
                    // Update varyEnergyShift column (column 16) based on useAsNuisance
                    QModelIndex varyIndex = segmentsModel->index(segmentIndex, 16);
                    segmentsModel->setData(varyIndex, param.useAsNuisance ? 1 : 0, Qt::EditRole);
                }
            }
        }
    }
}

void FittingTab::updateParameterTableValue(const QString& paramName, double value) {
    // Determine which table the parameter belongs to
    QTableWidget* targetTable = nullptr;
    if(paramName.contains("Normalization")) {
        targetTable = normParamsTable;
    } else if(paramName.contains("Energy Shift")) {
        targetTable = shiftParamsTable;
    } else {
        targetTable = levelParamsTable;
    }
    
    if(!targetTable) return;
    
    // Find the row with this parameter name
    for(int row = 0; row < targetTable->rowCount(); row++) {
        QTableWidgetItem* nameItem = targetTable->item(row, 0);
        if(nameItem && nameItem->text() == paramName) {
            // Update the value column (column 1)
            QTableWidgetItem* valueItem = targetTable->item(row, 1);
            if(valueItem) {
                // Temporarily disconnect signals to avoid recursion
                targetTable->blockSignals(true);
                valueItem->setText(QString::number(value, 'g', 6));
                targetTable->blockSignals(false);
            }
            break;
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

void FittingTab::onSegmentNormalizationChanged(int segmentIndex, double value) {
    // Find normalization parameter by segment index stored in channelIndex
    for(int i = 0; i < fittingParameters.size(); i++) {
        if(fittingParameters[i].category == "norm" && fittingParameters[i].channelIndex == segmentIndex) {
            fittingParameters[i].value = value;
            
            // Update the corresponding table cell
            updateParameterTableValue(fittingParameters[i].name, value);
            break;
        }
    }
}

void FittingTab::onSegmentEnergyShiftChanged(int segmentIndex, double value) {
    // Find energy shift parameter by segment index stored in channelIndex
    for(int i = 0; i < fittingParameters.size(); i++) {
        if(fittingParameters[i].category == "shift" && fittingParameters[i].channelIndex == segmentIndex) {
            fittingParameters[i].value = value;
            
            // Update the corresponding table cell
            updateParameterTableValue(fittingParameters[i].name, value);
            break;
        }
    }
}

void FittingTab::updateParameterTableError(const QString& paramName, double error) {
    // Determine which table the parameter belongs to
    QTableWidget* targetTable = nullptr;
    if(paramName.contains("Normalization")) {
        targetTable = normParamsTable;
    } else if(paramName.contains("Energy Shift")) {
        targetTable = shiftParamsTable;
    } else {
        targetTable = levelParamsTable;
    }
    
    if(!targetTable) return;
    
    // Find the row with this parameter name
    for(int row = 0; row < targetTable->rowCount(); row++) {
        QTableWidgetItem* nameItem = targetTable->item(row, 0);
        if(nameItem && nameItem->text() == paramName) {
            // Update the error column (column 4)
            QTableWidgetItem* errorItem = targetTable->item(row, 4);
            if(errorItem) {
                // Temporarily disconnect signals to avoid recursion
                targetTable->blockSignals(true);
                errorItem->setText(QString::number(error, 'g', 6));
                targetTable->blockSignals(false);
            }
            break;
        }
    }
}

void FittingTab::updateParameterTableCheckbox(const QString& paramName, bool checked) {
    // Determine which table the parameter belongs to
    QTableWidget* targetTable = nullptr;
    if(paramName.contains("Normalization")) {
        targetTable = normParamsTable;
    } else if(paramName.contains("Energy Shift")) {
        targetTable = shiftParamsTable;
    } else {
        targetTable = levelParamsTable;
    }
    
    if(!targetTable) return;
    
    // Find the row with this parameter name
    for(int row = 0; row < targetTable->rowCount(); row++) {
        QTableWidgetItem* nameItem = targetTable->item(row, 0);
        if(nameItem && nameItem->text() == paramName) {
            // Update the checkbox column (column 5)
            QTableWidgetItem* checkboxItem = targetTable->item(row, 5);
            if(checkboxItem) {
                // Temporarily disconnect signals to avoid recursion
                targetTable->blockSignals(true);
                checkboxItem->setCheckState(checked ? Qt::Checked : Qt::Unchecked);
                targetTable->blockSignals(false);
            }
            break;
        }
    }
}

void FittingTab::onSegmentNormalizationErrorChanged(int segmentIndex, double error) {
    // Find normalization parameter by segment index stored in channelIndex
    for(int i = 0; i < fittingParameters.size(); i++) {
        if(fittingParameters[i].category == "norm" && fittingParameters[i].channelIndex == segmentIndex) {
            fittingParameters[i].error = error;
            
            // Update the corresponding table cell (error column is column 4)
            updateParameterTableError(fittingParameters[i].name, error);
            break;
        }
    }
}

void FittingTab::onSegmentEnergyShiftErrorChanged(int segmentIndex, double error) {
    // Find energy shift parameter by segment index stored in channelIndex
    for(int i = 0; i < fittingParameters.size(); i++) {
        if(fittingParameters[i].category == "shift" && fittingParameters[i].channelIndex == segmentIndex) {
            fittingParameters[i].error = error;
            
            // Update the corresponding table cell (error column is column 4)
            updateParameterTableError(fittingParameters[i].name, error);
            break;
        }
    }
}

void FittingTab::onSegmentNormalizationVaryChanged(int segmentIndex, bool vary) {
    // Find normalization parameter by segment index stored in channelIndex
    for(int i = 0; i < fittingParameters.size(); i++) {
        if(fittingParameters[i].category == "norm" && fittingParameters[i].channelIndex == segmentIndex) {
            fittingParameters[i].useAsNuisance = vary;
            
            // Update the checkbox in the table
            updateParameterTableCheckbox(fittingParameters[i].name, vary);
            break;
        }
    }
}

void FittingTab::onSegmentEnergyShiftVaryChanged(int segmentIndex, bool vary) {
    // Find energy shift parameter by segment index stored in channelIndex
    for(int i = 0; i < fittingParameters.size(); i++) {
        if(fittingParameters[i].category == "shift" && fittingParameters[i].channelIndex == segmentIndex) {
            fittingParameters[i].useAsNuisance = vary;
            
            // Update the checkbox in the table
            updateParameterTableCheckbox(fittingParameters[i].name, vary);
            break;
        }
    }
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