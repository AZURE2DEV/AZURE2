#ifndef EXFORDIALOG_H
#define EXFORDIALOG_H

#include <QDialog>
#include "ExforData.h"

class PairsModel;
class QComboBox;
class QLineEdit;
class QPushButton;
class QTableWidget;
class QPlainTextEdit;
class QLabel;

/**
 * @brief Dialog to search and download experimental data from EXFOR.
 *
 * The user picks an entrance and an exit particle pair (taken from the project's
 * PairsModel); the EXFOR target / reaction / quantity strings are built
 * automatically but remain editable. A search lists the matching datasets; the
 * selected dataset can be previewed in full (raw EXFOR/CSV and the converted
 * AZURE2 columns) and then optionally saved to the data folder.
 */
class ExforDialog : public QDialog {
  Q_OBJECT

 public:
  explicit ExforDialog(PairsModel* pairsModel, QWidget* parent = nullptr);

 private slots:
  void updateQueryFromPairs();
  void doSearch();
  void onSearchFinished(const QList<ExforDataset>& datasets);
  void onDownloadFinished(const QString& rawCsv,
                          const QList<ExforPoint>& points, bool differential);
  void onError(const QString& message);
  void onDatasetSelectionChanged();
  void downloadSelected();
  void saveToDataFolder();

 private:
  void populatePairCombos();
  QString suggestedFileName() const;

  // EXFOR particle / nucleus code helpers built from a particle pair.
  static QString nucleusCode(int z, double mass);
  static QString lightParticleCode(int z, double mass, int pairType);

  PairsModel* pairsModel_;
  ExforData* exfor_;

  QComboBox* entranceCombo_;
  QComboBox* exitCombo_;
  QLineEdit* targetEdit_;
  QLineEdit* reactionEdit_;
  QComboBox* quantityCombo_;
  QPushButton* searchButton_;
  QLabel* statusLabel_;

  QTableWidget* resultsTable_;
  QLabel* datasetInfoLabel_;
  QPushButton* downloadButton_;

  QPlainTextEdit* rawView_;
  QPlainTextEdit* azureView_;
  QPushButton* saveButton_;

  QString currentDatasetId_;
  bool haveData_;
};

#endif  // EXFORDIALOG_H
