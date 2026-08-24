#ifndef IAEALEVELDATA_H
#define IAEALEVELDATA_H

#include <QObject>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QString>
#include <QList>
#include <QJsonObject>

/**
 * @brief Structure to hold nuclear level data from IAEA database
 */
struct IAEALevel {
  double energy;          // Level energy in MeV
  QString jp;             // Spin-parity value (e.g., "1/2+", "3-", "2")
  QString jpUncertainty;  // Uncertainty or alternative values for jp

  // Helper method to parse jp string into J value and parity
  bool parseJP(double &jValue, int &piValue) const;

  // Get list of possible J and pi combinations if jp is ambiguous
  QList<QPair<double, int>> getPossibleJPCombinations() const;
};

/**
 * @brief Class to retrieve nuclear level data from IAEA API
 *
 * This class handles querying the IAEA Nuclear Data Services API
 * to retrieve level scheme information for specified nuclides.
 */
class IAEALevelData : public QObject {
  Q_OBJECT

 public:
  explicit IAEALevelData(QObject *parent = nullptr);
  ~IAEALevelData();

  /**
   * @brief Query IAEA database for level data
   * @param nuclide Nuclide identifier (e.g., "15o", "12c")
   */
  void queryLevels(const QString &nuclide);

  /**
   * @brief Get the list of retrieved levels
   * @return List of IAEA levels
   */
  QList<IAEALevel> getLevels() const { return levels_; }

  /**
   * @brief Check if a query is currently in progress
   * @return true if query is in progress
   */
  bool isQueryInProgress() const { return queryInProgress_; }

  /**
   * @brief Get the last error message
   * @return Error message string
   */
  QString getLastError() const { return lastError_; }

 signals:
  /**
   * @brief Emitted when level data has been successfully retrieved
   * @param levels List of retrieved levels
   */
  void levelsRetrieved(const QList<IAEALevel> &levels);

  /**
   * @brief Emitted when an error occurs during query
   * @param errorMessage Description of the error
   */
  void queryError(const QString &errorMessage);

 private slots:
  void handleNetworkReply(QNetworkReply *reply);

 private:
  QNetworkAccessManager *networkManager_;
  QList<IAEALevel> levels_;
  bool queryInProgress_;
  QString lastError_;

  bool parseJsonResponse(const QByteArray &data);
};

#endif  // IAEALEVELDATA_H
