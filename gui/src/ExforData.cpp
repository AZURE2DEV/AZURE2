#include "ExforData.h"

#include <QNetworkRequest>
#include <QUrl>
#include <QUrlQuery>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QStringList>
#include <cmath>

// Base URL of the IAEA EXFOR Web API (see https://nds.iaea.org/exfor/x4guide/API/)
static const char* kExforBase = "https://nds.iaea.org/exfor";

// Physical constants kept identical to AZURE2's include/Constants.h so the
// astrophysical S-factor <-> cross-section conversion matches what AZURE2 does
// internally (EPoint::CalcEDependentValues / AZUREPlot).
static const double kPi     = 3.141592650;
static const double kUconv  = 931.4940880;
static const double kFstruc = 1.00 / 137.0359996790;

ExforData::ExforData(QObject* parent)
    : QObject(parent),
      manager_(new QNetworkAccessManager(this)),
      pending_(NoRequest),
      lastWasDifferential_(false),
      cz1_(0.0), cz2_(0.0), cm1_(0.0), cm2_(0.0) {
  connect(manager_, &QNetworkAccessManager::finished,
          this, &ExforData::handleReply);
}

void ExforData::searchDatasets(const QString& target, const QString& reaction,
                               const QString& quantity) {
  QUrl url(QString("%1/x4list").arg(kExforBase));
  QUrlQuery query;
  if (!target.isEmpty())   query.addQueryItem("Target", target);
  if (!reaction.isEmpty()) query.addQueryItem("Reaction", reaction);
  if (!quantity.isEmpty()) query.addQueryItem("Quantity", quantity);
  query.addQueryItem("json", QString());
  url.setQuery(query);

  pending_ = SearchRequest;
  lastError_.clear();
  QNetworkRequest request(url);
  request.setHeader(QNetworkRequest::UserAgentHeader,
                    "AZURE2-EXFOR-client/1.0");
  manager_->get(request);
}

void ExforData::downloadDataset(const QString& datasetId,
                                double z1, double z2, double m1, double m2) {
  cz1_ = z1; cz2_ = z2; cm1_ = m1; cm2_ = m2;
  QUrl url(QString("%1/x4get").arg(kExforBase));
  QUrlQuery query;
  query.addQueryItem("DatasetID", datasetId);
  query.addQueryItem("op", "csv");  // native-unit CSV.
  // NB: we deliberately do *not* request the computational ("plus") format.
  // That format silently converts energies to the lab frame and rescales the
  // columns to base units, which collided with the conversions done here. The
  // plain CSV keeps the author's units and the true centre-of-mass energy.
  url.setQuery(query);

  pending_ = DownloadRequest;
  lastError_.clear();
  QNetworkRequest request(url);
  request.setHeader(QNetworkRequest::UserAgentHeader,
                    "AZURE2-EXFOR-client/1.0");
  manager_->get(request);
}

void ExforData::handleReply(QNetworkReply* reply) {
  RequestType type = pending_;
  pending_ = NoRequest;
  reply->deleteLater();

  if (reply->error() != QNetworkReply::NoError) {
    lastError_ = reply->errorString();
    emit errorOccurred(lastError_);
    return;
  }

  QByteArray body = reply->readAll();

  if (type == SearchRequest) {
    datasets_.clear();
    QString error;
    if (!parseSearch(body, datasets_, error)) {
      lastError_ = error;
      emit errorOccurred(error);
      return;
    }
    emit searchFinished(datasets_);
  } else if (type == DownloadRequest) {
    lastRawCsv_ = QString::fromUtf8(body);
    lastPoints_.clear();
    QString error;
    if (!parseCsv(lastRawCsv_, lastPoints_, lastWasDifferential_, error)) {
      lastError_ = error;
      emit errorOccurred(error);
      return;
    }
    emit downloadFinished(lastRawCsv_, lastPoints_, lastWasDifferential_);
  }
}

bool ExforData::parseSearch(const QByteArray& data, QList<ExforDataset>& out,
                            QString& error) {
  QJsonParseError jsonError;
  QJsonDocument doc = QJsonDocument::fromJson(data, &jsonError);
  if (doc.isNull() || !doc.isObject()) {
    QString snippet = QString::fromUtf8(data).trimmed().left(200);
    error = QString("Unexpected response from EXFOR server: %1").arg(snippet);
    return false;
  }
  QJsonArray arr = doc.object().value("x4Datasets").toArray();
  for (const QJsonValue& v : arr) {
    QJsonObject o = v.toObject();
    ExforDataset ds;
    ds.id = o.value("id").toString();
    ds.reactionCode = o.value("RC").toString();
    ds.author = o.value("A1").toString();
    ds.reference = o.value("ref").toString();
    ds.npts = o.value("npts").toInt();
    ds.enMin = o.value("enMin").toDouble();
    // enMax is absent for single-point datasets; fall back to enMin
    ds.enMax = o.contains("enMax") ? o.value("enMax").toDouble() : ds.enMin;
    if (!ds.id.isEmpty()) out.append(ds);
  }
  return true;
}

// Split a CSV line honouring double-quoted fields (the trailing Reacode field
// contains commas).
static QStringList splitCsvLine(const QString& line) {
  QStringList fields;
  QString cur;
  bool inQuotes = false;
  for (int i = 0; i < line.size(); ++i) {
    ushort c = line.at(i).unicode();
    if (c == '"') {
      inQuotes = !inQuotes;
    } else if (c == ',' && !inQuotes) {
      fields.append(cur);
      cur.clear();
    } else {
      cur.append(c);
    }
  }
  fields.append(cur);
  return fields;
}

// Multiplier to convert an energy column to MeV based on its unit token.
static double energyToMeV(const QString& header) {
  QString h = header.toUpper();
  if (h.contains("(MEV)")) return 1.0;
  if (h.contains("(KEV)")) return 1.0e-3;
  if (h.contains("(EV)"))  return 1.0e-6;
  return 1.0e-6;  // EXFOR computational default is eV
}

// Multiplier to convert a cross-section column to barn (or barn/sr).
static double crossToBarn(const QString& header) {
  QString h = header.toUpper();
  if (h.contains("(MICRO-B")) return 1.0e-6;
  if (h.contains("(MB"))      return 1.0e-3;  // MB or MB/SR
  if (h.contains("(B"))       return 1.0;     // B or B/SR
  return 1.0;
}

// An astrophysical S-factor column has barn*energy units, e.g. "(B*MEV)",
// "(B*KEV)", "(MB*MEV)". Returns true and the factor that converts the value
// to barn*MeV; returns false for ordinary (non S-factor) columns. The
// energy-unit part (->MeV) is also returned separately: it is the unit at
// which S(E) was tabulated and therefore the natural unit of the incident
// energy axis as well (see parseCsv).
static bool sFactorBarnMeV(const QString& header, double& factor,
                           double& energyFactor) {
  QString h = header.toUpper();
  // Must look like cross-section * energy.
  if (!h.contains("*") || !h.contains("EV")) return false;
  double barn = 1.0;
  if (h.contains("MICRO-B")) barn = 1.0e-6;
  else if (h.contains("MB")) barn = 1.0e-3;
  else if (h.contains("B"))  barn = 1.0;
  double energy = 1.0e-6;  // EV -> MeV by default
  if (h.contains("MEV"))      energy = 1.0;
  else if (h.contains("KEV")) energy = 1.0e-3;
  else if (h.contains("EV"))  energy = 1.0e-6;
  factor = barn * energy;  // -> barn*MeV
  energyFactor = energy;
  return true;
}

bool ExforData::parseCsv(const QString& csv, QList<ExforPoint>& out,
                         bool& differential, QString& error) const {
  QStringList lines = csv.split('\n', QString::SkipEmptyParts);
  if (lines.isEmpty()) {
    error = "Empty dataset returned by EXFOR server.";
    return false;
  }
  // A wrong request returns a short error marker rather than CSV.
  if (lines.first().contains("wrong request") || !lines.first().contains(',')) {
    error = QString("EXFOR server could not return this dataset: %1")
                .arg(lines.first().trimmed());
    return false;
  }

  // Default relative uncertainty assumed when a point has no error given.
  const double kDefaultRelError = 0.05;  // 5%

  QStringList headers = splitCsvLine(lines.first());
  int eCol = -1, eCmCol = -1, dataCol = -1, angCol = -1;
  int errAbsCol = -1, errPctCol = -1;  // absolute (barn) and percent error cols
  bool angIsCosine = false;
  bool absIsTotal = false, pctIsTotal = false;
  bool isSFactor = false;    // data column is an astrophysical S-factor
  bool energyIsCM = false;   // energy column is centre-of-mass (EN-CM)
  double eFactor = 1.0e-6, eCmFactor = 1.0e-6;
  double crossFactor = 1.0, errAbsFactor = 1.0;
  double sFactor = 1.0;      // S-factor value -> barn*MeV
  double sFactorEnergyFactor = 0.0;  // S-factor's energy unit -> MeV

  for (int i = 0; i < headers.size(); ++i) {
    QString h = headers.at(i).trimmed();
    QString up = h.toUpper();

    // Incident energy column. Prefer the lab energy "EN", but also accept a
    // centre-of-mass energy "EN-CM" (converted to lab below). Skip EN-RSL,
    // EN-ERR, ... The unit may be EV, KEV or MEV - energyToMeV() reads it.
    if (eCol < 0 && (up.startsWith("EN ") || up.startsWith("EN(")) &&
        up.contains("EV")) {
      eCol = i;
      eFactor = energyToMeV(h);
    } else if (eCmCol < 0 && up.startsWith("EN-CM") && up.contains("EV")) {
      eCmCol = i;
      eCmFactor = energyToMeV(h);
    }

    // Cross-section (or S-factor) value column "DATA (B)" / "DATA (B/SR)" /
    // "DATA (B*MEV)" / ... (exclude "DATA-ERR ...").
    if (dataCol < 0 && (up.startsWith("DATA ") || up.startsWith("DATA("))) {
      dataCol = i;
      if (sFactorBarnMeV(h, sFactor, sFactorEnergyFactor)) {
        isSFactor = true;
      } else {
        crossFactor = crossToBarn(h);
      }
    }

    // Angle column.
    if (angCol < 0 && up.startsWith("ANG") && up.contains("ADEG")) {
      angCol = i;
      angIsCosine = false;
    } else if (angCol < 0 && up.startsWith("COS")) {
      angCol = i;
      angIsCosine = true;
    }

    // Error columns: keep the best absolute (in barn) and percent candidate,
    // preferring the total error (ERR-T) over partial components.
    if (up.startsWith("ERR")) {
      bool isTotal = up.startsWith("ERR-T");
      bool isPct = up.contains("PER-CENT") || up.contains("PERCENT") ||
                   up.contains("(PER");
      bool isAbs = up.contains("(B");  // matches (B), (B/SR), (MB), ...
      if (isAbs && (errAbsCol < 0 || (isTotal && !absIsTotal))) {
        errAbsCol = i;
        errAbsFactor = crossToBarn(h);
        absIsTotal = isTotal;
      } else if (isPct && (errPctCol < 0 || (isTotal && !pctIsTotal))) {
        errPctCol = i;
        pctIsTotal = isTotal;
      }
    }
  }

  // If no lab energy was found, fall back to the centre-of-mass energy.
  if (eCol < 0 && eCmCol >= 0) {
    eCol = eCmCol;
    eFactor = eCmFactor;
    energyIsCM = true;
  }

  // EXFOR occasionally mislabels the incident-energy unit of S-factor
  // datasets (e.g. O2597002 stores keV values in a column tagged "MEV", even
  // though its S-factor column is correctly tagged "B*KEV"). The energy unit
  // of the S-factor itself is the unit at which S(E) was tabulated, so it is
  // the self-consistent choice for the energy axis too. Trust it over the
  // (possibly wrong) unit of the energy column.
  if (isSFactor && sFactorEnergyFactor > 0.0) {
    eFactor = sFactorEnergyFactor;
  }

  if (eCol < 0 || dataCol < 0) {
    error = QString("Could not find energy/cross-section columns. "
                    "Headers: %1").arg(headers.join(", "));
    return false;
  }
  differential = (angCol >= 0);

  for (int li = 1; li < lines.size(); ++li) {
    QStringList f = splitCsvLine(lines.at(li));
    int needed = qMax(qMax(eCol, dataCol),
                      qMax(qMax(errAbsCol, errPctCol), angCol));
    if (f.size() <= needed) continue;

    bool okE = false, okD = false;
    double e = f.at(eCol).trimmed().toDouble(&okE);
    double d = f.at(dataCol).trimmed().toDouble(&okD);
    if (!okE || !okD) continue;  // skip header echoes / blank rows

    ExforPoint p;
    // Energy in the column's frame (MeV); derive both lab and cm. AZURE2 wants
    // the lab energy, while the S-factor formula needs the cm energy.
    double eRaw = e * eFactor;
    double totM = cm1_ + cm2_;
    double eLab, eCm;
    if (energyIsCM) {
      eCm = eRaw;
      eLab = (cm2_ > 0.0) ? eRaw * totM / cm2_ : eRaw;
    } else {
      eLab = eRaw;
      eCm = (totM > 0.0) ? eRaw * cm2_ / totM : eRaw;
    }
    p.energy = eLab;

    // Conversion factor applied to the data value (and to an absolute error in
    // the same units). For a plain cross section this is just the unit factor;
    // for an S-factor it also undoes the astrophysical parametrisation:
    //   sigma(E) = S(E)/E_cm * exp(-2*pi*eta).
    double crossConv = crossFactor;
    double errConv = errAbsFactor;
    if (isSFactor) {
      if (eCm <= 0.0) continue;  // S-factor undefined at zero energy
      // Reduced mass in amu, exactly as PPair::GetRedMass() computes it.
      double mu = (totM > 0.0) ? cm1_ * cm2_ / totM : 0.0;
      // sigma(E) = S(E) / (E_cm * exp(2*pi*eta)) with the Sommerfeld parameter
      //   2*pi*eta = 2*pi*sqrt(uconv/2)*fstruc * Z1*Z2 * sqrt(mu/E_cm),
      // identical to AZURE2's S-factor conversion in
      // EPoint::CalcEDependentValues (mu in amu, E_cm in MeV).
      double twoPiEta = 2.0 * kPi * std::sqrt(kUconv / 2.0) * kFstruc *
                        cz1_ * cz2_ * std::sqrt(mu / eCm);
      double sigmaPerS = std::exp(-twoPiEta) / eCm;  // (barn*MeV) -> barn
      crossConv = sFactor * sigmaPerS;
      errConv = sFactor * sigmaPerS;  // abs error assumed in same S-factor units
    }
    p.cross = d * crossConv;

    // Error: prefer absolute (converted to barn), else a percent error
    // applied to the value, else fall back to a default relative uncertainty.
    bool haveErr = false;
    if (errAbsCol >= 0) {
      bool ok = false;
      double er = f.at(errAbsCol).trimmed().toDouble(&ok);
      if (ok && er != 0.0) { p.error = er * errConv; haveErr = true; }
    }
    if (!haveErr && errPctCol >= 0) {
      bool ok = false;
      double pct = f.at(errPctCol).trimmed().toDouble(&ok);
      if (ok && pct != 0.0) {
        p.error = std::fabs(p.cross) * pct / 100.0;
        haveErr = true;
      }
    }
    if (!haveErr) p.error = kDefaultRelError * std::fabs(p.cross);

    if (angCol >= 0) {
      bool okA = false;
      double a = f.at(angCol).trimmed().toDouble(&okA);
      if (okA) {
        p.angle = angIsCosine ? std::acos(qBound(-1.0, a, 1.0)) * 180.0 / M_PI
                              : a;
      }
    }
    out.append(p);
  }

  if (out.isEmpty()) {
    error = "No numeric data points could be parsed from this dataset.";
    return false;
  }
  return true;
}

QString ExforData::toAzureText(const QList<ExforPoint>& points) {
  QString text;
  for (const ExforPoint& p : points) {
    text += QString("%1\t%2\t%3\t%4\n")
                .arg(p.energy, 0, 'g', 8)
                .arg(p.angle, 0, 'g', 8)
                .arg(p.cross, 0, 'g', 8)
                .arg(p.error, 0, 'g', 8);
  }
  return text;
}
