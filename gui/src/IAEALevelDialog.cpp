#include "IAEALevelDialog.h"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QMessageBox>
#include <QCheckBox>
#include <cmath>

IAEALevelDialog::IAEALevelDialog(QWidget *parent)
    : QDialog(parent)
    , iaeaData_(new IAEALevelData(this))
    , massNumber_(0)
    , atomicNumber_(0)
{
    setWindowTitle(tr("Import Levels from IAEA Database"));
    resize(700, 500);

    setupUI();

    // Connect IAEA data signals
    connect(iaeaData_, &IAEALevelData::levelsRetrieved,
            this, &IAEALevelDialog::onLevelsRetrieved);
    connect(iaeaData_, &IAEALevelData::queryError,
            this, &IAEALevelDialog::onQueryError);
}

void IAEALevelDialog::setupUI() {
    QVBoxLayout *mainLayout = new QVBoxLayout(this);

    // Query section
    QHBoxLayout *queryLayout = new QHBoxLayout();
    nuclideLabel_ = new QLabel(tr("Nuclide:"));
    nuclideEdit_ = new QLineEdit();
    nuclideEdit_->setPlaceholderText(tr("e.g., 15o, 12c"));
    queryButton_ = new QPushButton(tr("Query IAEA Database"));
    connect(queryButton_, &QPushButton::clicked,
            this, &IAEALevelDialog::onQueryButtonClicked);

    queryLayout->addWidget(nuclideLabel_);
    queryLayout->addWidget(nuclideEdit_);
    queryLayout->addWidget(queryButton_);
    mainLayout->addLayout(queryLayout);

    // Progress bar
    progressBar_ = new QProgressBar();
    progressBar_->setRange(0, 0); // Indeterminate mode
    progressBar_->hide();
    mainLayout->addWidget(progressBar_);

    // Status label
    statusLabel_ = new QLabel();
    statusLabel_->setWordWrap(true);
    mainLayout->addWidget(statusLabel_);

    // Levels table
    levelsTable_ = new QTableWidget();
    levelsTable_->setColumnCount(5);
    levelsTable_->setHorizontalHeaderLabels(
        QStringList() << tr("Select") << tr("Energy (MeV)")
                      << tr("J^π (IAEA)") << tr("J Value") << tr("Parity"));
    levelsTable_->horizontalHeader()->setStretchLastSection(true);
    levelsTable_->setSelectionMode(QAbstractItemView::NoSelection);
    levelsTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    levelsTable_->setColumnWidth(0, 60);
    levelsTable_->setColumnWidth(1, 120);
    levelsTable_->setColumnWidth(2, 100);
    levelsTable_->setColumnWidth(3, 150);
    levelsTable_->setColumnWidth(4, 150);
    mainLayout->addWidget(levelsTable_);

    // Buttons
    QHBoxLayout *buttonLayout = new QHBoxLayout();
    addSelectedButton_ = new QPushButton(tr("Add Selected Levels"));
    addSelectedButton_->setEnabled(false);
    connect(addSelectedButton_, &QPushButton::clicked,
            this, &IAEALevelDialog::onAddSelectedClicked);

    cancelButton_ = new QPushButton(tr("Cancel"));
    connect(cancelButton_, &QPushButton::clicked, this, &QDialog::reject);

    buttonLayout->addStretch();
    buttonLayout->addWidget(addSelectedButton_);
    buttonLayout->addWidget(cancelButton_);
    mainLayout->addLayout(buttonLayout);
}

void IAEALevelDialog::setCompoundNucleus(int massNumber, int atomicNumber) {
    massNumber_ = massNumber;
    atomicNumber_ = atomicNumber;

    QString elementSymbol = getElementSymbol(atomicNumber);
    QString nuclide = QString("%1%2").arg(massNumber).arg(elementSymbol.toLower());

    nuclideEdit_->setText(nuclide);
    statusLabel_->setText(tr("Compound nucleus: %1 (A=%2, Z=%3)")
                         .arg(nuclide.toUpper())
                         .arg(massNumber)
                         .arg(atomicNumber));
}

void IAEALevelDialog::onQueryButtonClicked() {
    QString nuclide = nuclideEdit_->text().trimmed();
    if (nuclide.isEmpty()) {
        QMessageBox::warning(this, tr("Input Error"),
                           tr("Please enter a nuclide identifier (e.g., 15o)"));
        return;
    }

    progressBar_->show();
    queryButton_->setEnabled(false);
    statusLabel_->setText(tr("Querying IAEA database for %1...").arg(nuclide));
    levelsTable_->setRowCount(0);
    currentLevels_.clear();

    iaeaData_->queryLevels(nuclide);
}

void IAEALevelDialog::onLevelsRetrieved(const QList<IAEALevel> &levels) {
    progressBar_->hide();
    queryButton_->setEnabled(true);

    if (levels.isEmpty()) {
        statusLabel_->setText(tr("No levels found for the specified nuclide."));
        return;
    }

    currentLevels_ = levels;
    populateTable(levels);
    statusLabel_->setText(tr("Retrieved %1 levels. Select levels to import and specify J and parity.")
                         .arg(levels.size()));
}

void IAEALevelDialog::onQueryError(const QString &errorMessage) {
    progressBar_->hide();
    queryButton_->setEnabled(true);
    statusLabel_->setText(tr("Error: %1").arg(errorMessage));
    QMessageBox::warning(this, tr("Query Error"), errorMessage);
}

void IAEALevelDialog::populateTable(const QList<IAEALevel> &levels) {
    levelsTable_->setRowCount(levels.size());

    for (int i = 0; i < levels.size(); ++i) {
        const IAEALevel &level = levels[i];

        // Checkbox for selection
        QCheckBox *checkBox = new QCheckBox();
        connect(checkBox, &QCheckBox::stateChanged,
                this, &IAEALevelDialog::onLevelSelectionChanged);
        levelsTable_->setCellWidget(i, 0, checkBox);

        // Energy
        QTableWidgetItem *energyItem = new QTableWidgetItem(QString::number(level.energy, 'f', 4));
        energyItem->setTextAlignment(Qt::AlignCenter);
        levelsTable_->setItem(i, 1, energyItem);

        // JP from IAEA
        QTableWidgetItem *jpItem = new QTableWidgetItem(level.jp);
        jpItem->setTextAlignment(Qt::AlignCenter);
        levelsTable_->setItem(i, 2, jpItem);

        // Get possible J and parity combinations
        QList<QPair<double, int>> combinations = level.getPossibleJPCombinations();

        // J value combo box
        QComboBox *jCombo = new QComboBox();
        if (!combinations.isEmpty()) {
            for (const auto &combo : combinations) {
                double j = combo.first;
                QString jStr;
                // Format J as integer or half-integer
                if (std::fabs(j - std::round(j)) < 0.01) {
                    jStr = QString::number(static_cast<int>(std::round(j)));
                } else {
                    jStr = QString::number(j, 'f', 1);
                }
                jCombo->addItem(jStr, j);
            }
        } else {
            // If parsing failed, allow manual input
            jCombo->setEditable(true);
            jCombo->setPlaceholderText(tr("Enter J"));
        }
        levelsTable_->setCellWidget(i, 3, jCombo);

        // Parity combo box
        QComboBox *piCombo = new QComboBox();
        piCombo->addItem(tr("-"), -1);
        piCombo->addItem(tr("+"), 1);

        // Set default parity from first combination if available
        if (!combinations.isEmpty()) {
            int parity = combinations.first().second;
            piCombo->setCurrentIndex(parity == 1 ? 1 : 0);
        }

        levelsTable_->setCellWidget(i, 4, piCombo);
    }
}

void IAEALevelDialog::onLevelSelectionChanged() {
    updateAddButton();
}

void IAEALevelDialog::updateAddButton() {
    int selectedCount = 0;
    for (int i = 0; i < levelsTable_->rowCount(); ++i) {
        QCheckBox *checkBox = qobject_cast<QCheckBox*>(levelsTable_->cellWidget(i, 0));
        if (checkBox && checkBox->isChecked()) {
            selectedCount++;
        }
    }
    addSelectedButton_->setEnabled(selectedCount > 0);
}

void IAEALevelDialog::onJPComboChanged(int row, int index) {
    // This can be used for future enhancements
    Q_UNUSED(row);
    Q_UNUSED(index);
}

void IAEALevelDialog::onAddSelectedClicked() {
    selectedLevels_.clear();

    for (int i = 0; i < levelsTable_->rowCount(); ++i) {
        QCheckBox *checkBox = qobject_cast<QCheckBox*>(levelsTable_->cellWidget(i, 0));
        if (!checkBox || !checkBox->isChecked()) {
            continue;
        }

        // Get energy
        double energy = levelsTable_->item(i, 1)->text().toDouble();

        // Get J value
        QComboBox *jCombo = qobject_cast<QComboBox*>(levelsTable_->cellWidget(i, 3));
        double jValue = 0.0;
        if (jCombo) {
            if (jCombo->isEditable()) {
                jValue = jCombo->currentText().toDouble();
            } else {
                jValue = jCombo->currentData().toDouble();
            }
        }

        // Get parity
        QComboBox *piCombo = qobject_cast<QComboBox*>(levelsTable_->cellWidget(i, 4));
        int piValue = -1;
        if (piCombo) {
            piValue = piCombo->currentData().toInt();
        }

        // Create LevelsData structure
        LevelsData level;
        level.isActive = 1;
        level.isFixed = 0;
        level.jValue = jValue;
        level.piValue = piValue;
        level.energy = energy;

        selectedLevels_.append(level);
    }

    if (selectedLevels_.isEmpty()) {
        QMessageBox::warning(this, tr("No Selection"),
                           tr("Please select at least one level to import."));
        return;
    }

    accept();
}

QList<LevelsData> IAEALevelDialog::getSelectedLevels() const {
    return selectedLevels_;
}

QString IAEALevelDialog::getElementSymbol(int atomicNumber) const {
    // Simple lookup table for common elements
    static const QMap<int, QString> elements = {
        {1, "H"}, {2, "He"}, {3, "Li"}, {4, "Be"}, {5, "B"}, {6, "C"}, {7, "N"}, {8, "O"},
        {9, "F"}, {10, "Ne"}, {11, "Na"}, {12, "Mg"}, {13, "Al"}, {14, "Si"}, {15, "P"},
        {16, "S"}, {17, "Cl"}, {18, "Ar"}, {19, "K"}, {20, "Ca"}, {21, "Sc"}, {22, "Ti"},
        {23, "V"}, {24, "Cr"}, {25, "Mn"}, {26, "Fe"}, {27, "Co"}, {28, "Ni"}, {29, "Cu"},
        {30, "Zn"}, {31, "Ga"}, {32, "Ge"}, {33, "As"}, {34, "Se"}, {35, "Br"}, {36, "Kr"}
    };

    return elements.value(atomicNumber, "X");
}
