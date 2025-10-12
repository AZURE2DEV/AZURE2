#ifndef IAEALEVELDIALOG_H
#define IAEALEVELDIALOG_H

#include <QDialog>
#include <QTableWidget>
#include <QPushButton>
#include <QLabel>
#include <QComboBox>
#include <QLineEdit>
#include <QProgressBar>
#include "IAEALevelData.h"
#include "LevelsModel.h"

/**
 * @brief Dialog for querying and selecting nuclear levels from IAEA database
 *
 * This dialog allows users to:
 * 1. Automatically detect compound nucleus from pairs
 * 2. Query IAEA database for level schemes
 * 3. Select which levels to import
 * 4. Manually specify J and parity for ambiguous cases
 */
class IAEALevelDialog : public QDialog {
    Q_OBJECT

public:
    explicit IAEALevelDialog(QWidget *parent = nullptr);

    /**
     * @brief Set the compound nucleus to query
     * @param massNumber Mass number (A)
     * @param atomicNumber Atomic number (Z)
     */
    void setCompoundNucleus(int massNumber, int atomicNumber);

    /**
     * @brief Get the list of selected levels with user-specified J and parity
     * @return List of LevelsData ready to add to the model
     */
    QList<LevelsData> getSelectedLevels() const;

private slots:
    void onQueryButtonClicked();
    void onLevelsRetrieved(const QList<IAEALevel> &levels);
    void onQueryError(const QString &errorMessage);
    void onLevelSelectionChanged();
    void onJPComboChanged(int row, int index);
    void onAddSelectedClicked();

private:
    void setupUI();
    void populateTable(const QList<IAEALevel> &levels);
    QString getElementSymbol(int atomicNumber) const;
    void updateAddButton();

    IAEALevelData *iaeaData_;

    QLabel *nuclideLabel_;
    QLineEdit *nuclideEdit_;
    QPushButton *queryButton_;
    QProgressBar *progressBar_;
    QLabel *statusLabel_;

    QTableWidget *levelsTable_;
    QPushButton *addSelectedButton_;
    QPushButton *cancelButton_;

    QList<IAEALevel> currentLevels_;
    QList<LevelsData> selectedLevels_;

    int massNumber_;
    int atomicNumber_;
};

#endif // IAEALEVELDIALOG_H
