#ifndef EXFORDATA_H
#define EXFORDATA_H

#include <QObject>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QString>
#include <QList>

/**
 * @brief One dataset descriptor returned by the EXFOR x4list search.
 */
struct ExforDataset {
  QString id;            // EXFOR DatasetID, e.g. "D4426002"
  QString reactionCode;  // RC, e.g. "6-C-12(P,G)7-N-13,,SIG"
  QString author;        // A1, e.g. "2023,Gy.Gyurky+"
  QString reference;     // ref string
  int npts;              // number of data points
  double enMin;          // minimum incident energy (eV)
  double enMax;          // maximum incident energy (eV)
  ExforDataset() : npts(0), enMin(0.0), enMax(0.0) {}
};

/**
 * @brief One data point converted to AZURE2 conventions.
 *
 * AZURE2 data files have four columns: energy (MeV, lab), angle (deg),
 * cross section (barn for integrated, barn/sr for differential) and its error.
 */
struct ExforPoint {
  double energy;   // MeV, lab frame
  double angle;    // degrees (0 for angle-integrated data)
  double cross;    // barn (SIG) or barn/sr (DA)
  double error;    // same units as cross (0 if not provided)
  ExforPoint() : energy(0.0), angle(0.0), cross(0.0), error(0.0) {}
};

/**
 * @brief Lightweight client for the IAEA EXFOR Web API.
 *
 * Mirrors the Python reference client (exfor_client.py): it searches for
 * datasets with x4list and downloads the data of a single dataset (in the
 * computational CSV format) with x4get, then converts it to the AZURE2 data
 * format. Only the energy/angle/cross/error needed by AZURE2 are extracted,
 * so the focus is on cross sections (SIG) and angular distributions (DA).
 */
class ExforData : public QObject {
  Q_OBJECT

 public:
  explicit ExforData(QObject* parent = nullptr);

  /// Search datasets matching target / reaction / quantity (EXFOR codes).
  void searchDatasets(const QString& target, const QString& reaction,
                       const QString& quantity);
  /// Download (and parse) the data of a single dataset by its DatasetID.
  /// The entrance-channel charges and masses (projectile 1, target 2, in amu)
  /// are used to convert astrophysical S-factor data back to cross sections.
  void downloadDataset(const QString& datasetId,
                       double z1 = 0.0, double z2 = 0.0,
                       double m1 = 0.0, double m2 = 0.0);

  QList<ExforDataset> datasets() const { return datasets_; }
  QList<ExforPoint> lastPoints() const { return lastPoints_; }
  QString lastRawCsv() const { return lastRawCsv_; }
  bool lastWasDifferential() const { return lastWasDifferential_; }
  QString lastError() const { return lastError_; }
  bool isBusy() const { return pending_ != NoRequest; }

  /// Convert parsed points to AZURE2 data-file text (energy angle cross error).
  static QString toAzureText(const QList<ExforPoint>& points);

 signals:
  void searchFinished(const QList<ExforDataset>& datasets);
  void downloadFinished(const QString& rawCsv,
                        const QList<ExforPoint>& points, bool differential);
  void errorOccurred(const QString& message);

 private slots:
  void handleReply(QNetworkReply* reply);

 private:
  enum RequestType { NoRequest, SearchRequest, DownloadRequest };

  static bool parseSearch(const QByteArray& data, QList<ExforDataset>& out,
                          QString& error);
  bool parseCsv(const QString& csv, QList<ExforPoint>& out,
                bool& differential, QString& error) const;

  QNetworkAccessManager* manager_;
  RequestType pending_;
  QList<ExforDataset> datasets_;
  QList<ExforPoint> lastPoints_;
  QString lastRawCsv_;
  bool lastWasDifferential_;
  QString lastError_;

  // Entrance-channel data for S-factor -> cross-section conversion.
  double cz1_, cz2_, cm1_, cm2_;
};

#endif  // EXFORDATA_H
