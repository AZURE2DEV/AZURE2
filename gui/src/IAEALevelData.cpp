#include "IAEALevelData.h"
#include <QNetworkRequest>
#include <QUrl>
#include <QJsonDocument>
#include <QJsonArray>
#include <QRegularExpression>
#include <QDebug>
#include <QStringList>

// Helper function to parse jp string into J value and parity
bool IAEALevel::parseJP(double &jValue, int &piValue) const {
    if (jp.isEmpty()) {
        return false;
    }

    // Get possible combinations
    QList<QPair<double, int>> combinations = getPossibleJPCombinations();

    if (combinations.isEmpty()) {
        return false;
    }

    // Return the first (most likely) combination
    jValue = combinations.first().first;
    piValue = combinations.first().second;
    return true;
}

// Parse jp string to get all possible J and parity combinations
QList<QPair<double, int>> IAEALevel::getPossibleJPCombinations() const {
    QList<QPair<double, int>> combinations;

    if (jp.isEmpty()) {
        return combinations;
    }

    // Remove whitespace and convert to lowercase for easier parsing
    QString jpClean = jp.trimmed();

    // Regular expression to match spin-parity patterns
    // Matches patterns like: "1/2+", "3-", "2", "(3/2)+", "1/2(-)", etc.
    QRegularExpression re(R"((\d+(?:/\d+)?)\s*([\+\-\(\)\?]*))");
    QRegularExpressionMatch match = re.match(jpClean);

    if (!match.hasMatch()) {
        return combinations;
    }

    QString spinStr = match.captured(1);
    QString parityStr = match.captured(2);

    // Parse spin value
    double spin = 0.0;
    if (spinStr.contains('/')) {
        QStringList parts = spinStr.split('/');
        if (parts.size() == 2) {
            spin = parts[0].toDouble() / parts[1].toDouble();
        }
    } else {
        spin = spinStr.toDouble();
    }

    // Parse parity
    // If parity is ambiguous or uncertain, add both possibilities
    bool hasPositive = parityStr.contains('+');
    bool hasNegative = parityStr.contains('-');
    bool isUncertain = parityStr.contains('?') || parityStr.contains('(');

    if (hasPositive && !hasNegative) {
        combinations.append(qMakePair(spin, 1));
    } else if (hasNegative && !hasPositive) {
        combinations.append(qMakePair(spin, -1));
    } else if (isUncertain || (!hasPositive && !hasNegative)) {
        // If uncertain or no parity specified, add both possibilities
        combinations.append(qMakePair(spin, 1));
        combinations.append(qMakePair(spin, -1));
    } else if (hasPositive && hasNegative) {
        // Both parities mentioned (rare case)
        combinations.append(qMakePair(spin, 1));
        combinations.append(qMakePair(spin, -1));
    }

    return combinations;
}

IAEALevelData::IAEALevelData(QObject *parent)
    : QObject(parent)
    , networkManager_(new QNetworkAccessManager(this))
    , queryInProgress_(false)
{
    connect(networkManager_, &QNetworkAccessManager::finished,
            this, &IAEALevelData::handleNetworkReply);
}

IAEALevelData::~IAEALevelData() {
    // QObject parent-child relationship handles cleanup
}

void IAEALevelData::queryLevels(const QString &nuclide) {
    if (queryInProgress_) {
        emit queryError("A query is already in progress");
        return;
    }

    if (nuclide.isEmpty()) {
        emit queryError("Nuclide identifier cannot be empty");
        return;
    }

    levels_.clear();
    lastError_.clear();
    queryInProgress_ = true;

    // Construct the IAEA API URL
    QString urlString = QString("https://nds.iaea.org/relnsd/v1/data?fields=levels&nuclides=%1")
                        .arg(nuclide.toLower());

    QUrl url(urlString);
    QNetworkRequest request(url);
    request.setHeader(QNetworkRequest::ContentTypeHeader, "application/json");

    //qDebug() << "Querying IAEA API:" << urlString;
    networkManager_->get(request);
}

void IAEALevelData::handleNetworkReply(QNetworkReply *reply) {
    queryInProgress_ = false;

    if (reply->error() != QNetworkReply::NoError) {
        lastError_ = QString("Network error: %1").arg(reply->errorString());
        qWarning() << lastError_;
        emit queryError(lastError_);
        reply->deleteLater();
        return;
    }

    QByteArray responseData = reply->readAll();
    reply->deleteLater();

    //qDebug() << "Received response:" << responseData.left(500); // Log first 500 chars

    if (parseJsonResponse(responseData)) {
        emit levelsRetrieved(levels_);
    } else {
        emit queryError(lastError_);
    }
}

bool IAEALevelData::parseJsonResponse(const QByteArray &data) {
    // The IAEA API returns CSV format, not JSON
    // Parse CSV with format: z,n,symbol,idx,energy_shift,energy,unc_e,ripl_shift,jp,...

    QString csvData = QString::fromUtf8(data);
    QStringList lines = csvData.split('\n', Qt::SkipEmptyParts);

    if (lines.isEmpty()) {
        lastError_ = "Empty response from IAEA";
        qWarning() << lastError_;
        return false;
    }

    // First line is the header
    if (lines.size() < 2) {
        lastError_ = "No data rows in response";
        qWarning() << lastError_;
        return false;
    }

    // Parse header to find column indices
    QStringList headers = lines[0].split(',');
    int energyCol = headers.indexOf("energy");
    int jpCol = headers.indexOf("jp");

    if (energyCol == -1 || jpCol == -1) {
        lastError_ = "Could not find required columns (energy, jp) in CSV";
        qWarning() << lastError_;
        return false;
    }

    // Parse data rows (skip header)
    for (int i = 1; i < lines.size(); ++i) {
        QString line = lines[i].trimmed();
        if (line.isEmpty()) {
            continue;
        }

        QStringList fields = line.split(',');
        if (fields.size() <= qMax(energyCol, jpCol)) {
            qWarning() << "Skipping malformed line:" << line;
            continue;
        }

        IAEALevel level;

        // Extract energy
        QString energyStr = fields[energyCol].trimmed();
        if (energyStr.isEmpty()) {
            continue; // Skip levels without energy
        }

        // Remove uncertainty information (e.g., "0.0+X", "1.234(5)")
        energyStr = energyStr.split('+').first().split('(').first().trimmed();
        bool ok;
        level.energy = energyStr.toDouble(&ok);
        if (!ok) {
            qWarning() << "Could not parse energy:" << energyStr;
            continue;
        }

        // Energy is in keV, convert to MeV
        level.energy /= 1000.0;

        // Extract jp (spin-parity)
        level.jp = fields[jpCol].trimmed();
        if (level.jp.isEmpty()) {
            // Some levels don't have known spin-parity
            level.jp = "?";
        }

        levels_.append(level);
    }

    if (levels_.isEmpty()) {
        lastError_ = "No valid levels found in response";
        qWarning() << lastError_;
        return false;
    }

    //qDebug() << "Parsed" << levels_.size() << "levels";
    return true;
}
