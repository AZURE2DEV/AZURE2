#include "ChannelRadiiModel.h"

#include <QStringList>
#include <QTextStream>
#include <cmath>

ChannelRadiiModel::ChannelRadiiModel(QObject *parent) : QAbstractTableModel(parent) {}

int ChannelRadiiModel::rowCount(const QModelIndex &parent) const { return parent.isValid() ? 0 : rows_.size(); }

int ChannelRadiiModel::columnCount(const QModelIndex &parent) const { return parent.isValid() ? 0 : NumCols; }

QVariant ChannelRadiiModel::data(const QModelIndex &index, int role) const {
  if (!index.isValid() || index.row() >= rows_.size()) return QVariant();
  const ChannelRadiusRow &r = rows_.at(index.row());
  if (role == Qt::DisplayRole || role == Qt::EditRole) {
    switch (index.column()) {
      case PairCol: return r.pairKey;
      case LCol: return r.l;
      case SCol: return (role == Qt::DisplayRole && r.s < 0.0) ? QString(tr("any")) : spinText(r.s);
      case JCol: return (role == Qt::DisplayRole && r.J < 0.0) ? QString(tr("any")) : spinText(r.J);
      case RadiusCol: return r.radius;
    }
  }
  if (role == Qt::TextAlignmentRole) return int(Qt::AlignCenter);
  if (role == Qt::ToolTipRole) {
    switch (index.column()) {
      case PairCol: return tr("Pair key of a particle pair (the Pairs table's first column).");
      case LCol: return tr("Orbital angular momentum l of the channels.");
      case SCol: return tr("Channel spin s (e.g. 1/2), or blank for every s.");
      case JCol: return tr("Total spin J of the level's J^pi group (e.g. 3/2), or blank for every J.");
      case RadiusCol: return tr("Channel radius in fm for the channels this row matches.");
    }
  }
  return QVariant();
}

QVariant ChannelRadiiModel::headerData(int section, Qt::Orientation orientation, int role) const {
  if (role != Qt::DisplayRole || orientation != Qt::Horizontal) return QVariant();
  switch (section) {
    case PairCol: return tr("Pair");
    case LCol: return tr("l");
    case SCol: return tr("s");
    case JCol: return tr("J");
    case RadiusCol: return tr("Radius [fm]");
  }
  return QVariant();
}

Qt::ItemFlags ChannelRadiiModel::flags(const QModelIndex &index) const {
  if (!index.isValid()) return Qt::ItemIsEnabled;
  return QAbstractTableModel::flags(index) | Qt::ItemIsEditable;
}

bool ChannelRadiiModel::setData(const QModelIndex &index, const QVariant &value, int role) {
  if (!index.isValid() || role != Qt::EditRole || index.row() >= rows_.size()) return false;
  ChannelRadiusRow &r = rows_[index.row()];
  bool ok = false;
  switch (index.column()) {
    case PairCol: {
      int v = value.toString().trimmed().toInt(&ok);
      if (!ok || v < 1) return false;
      r.pairKey = v;
      break;
    }
    case LCol: {
      int v = value.toString().trimmed().toInt(&ok);
      if (!ok || v < 0) return false;
      r.l = v;
      break;
    }
    case SCol:
    case JCol: {
      double v;
      if (!parseSpin(value.toString(), v)) return false;
      (index.column() == SCol ? r.s : r.J) = v;
      break;
    }
    case RadiusCol: {
      double v = value.toString().trimmed().toDouble(&ok);
      if (!ok || !(v > 0.0)) return false;
      r.radius = v;
      break;
    }
    default: return false;
  }
  emit dataChanged(index, index);
  emit rowsEdited();
  return true;
}

void ChannelRadiiModel::setRows(const QList<ChannelRadiusRow> &rows) {
  beginResetModel();
  rows_ = rows;
  endResetModel();
  emit rowsEdited();
}

int ChannelRadiiModel::addRow(int pairKey, int l, double radius) {
  beginInsertRows(QModelIndex(), rows_.size(), rows_.size());
  ChannelRadiusRow r;
  r.pairKey = pairKey;
  r.l = l;
  r.radius = radius;
  rows_.append(r);
  endInsertRows();
  emit rowsEdited();
  return rows_.size() - 1;
}

bool ChannelRadiiModel::removeRowAt(int row) {
  if (row < 0 || row >= rows_.size()) return false;
  beginRemoveRows(QModelIndex(), row, row);
  rows_.removeAt(row);
  endRemoveRows();
  emit rowsEdited();
  return true;
}

void ChannelRadiiModel::clear() { setRows(QList<ChannelRadiusRow>()); }

static bool sameSpin(double a, double b) { return std::fabs(a - b) < 1.0e-6; }

double ChannelRadiiModel::radiusFor(int pairKey, int l, double s, double J, double pairRadius) const {
  int best = 0;
  double radius = pairRadius;
  for (const ChannelRadiusRow &r : rows_) {
    if (r.pairKey != pairKey || r.l != l) continue;
    if (r.s >= 0.0 && !sameSpin(r.s, s)) continue;
    if (r.J >= 0.0 && !sameSpin(r.J, J)) continue;
    const int rank = 1 + (r.s >= 0.0 ? 1 : 0) + (r.J >= 0.0 ? 2 : 0);
    if (rank > best) {
      best = rank;
      radius = r.radius;
    }
  }
  return radius;
}

bool ChannelRadiiModel::parseSpin(const QString &text, double &value) {
  const QString t = text.trimmed();
  if (t.isEmpty() || t.compare("any", Qt::CaseInsensitive) == 0) {
    value = -1.0;
    return true;
  }
  bool ok = false;
  double v;
  if (t.contains('/')) {
    const QStringList parts = t.split('/');
    if (parts.size() != 2) return false;
    bool okNum = false, okDen = false;
    const double num = parts.at(0).toDouble(&okNum), den = parts.at(1).toDouble(&okDen);
    if (!okNum || !okDen || den == 0.0) return false;
    v = num / den;
  } else {
    v = t.toDouble(&ok);
    if (!ok) return false;
  }
  if (v < 0.0 || std::fabs(2.0 * v - std::round(2.0 * v)) > 1.0e-9) return false;
  value = 0.5 * std::round(2.0 * v);
  return true;
}

QString ChannelRadiiModel::spinText(double value) {
  if (value < 0.0) return QString();
  const int twice = (int)std::lround(2.0 * value);
  return (twice % 2 == 0) ? QString::number(twice / 2) : QString("%1/2").arg(twice);
}

QString ChannelRadiiModel::lineFor(const ChannelRadiusRow &row) {
  QString line = QString("%1 %2 %3").arg(row.pairKey).arg(row.l).arg(row.radius, 0, 'g', 10);
  if (row.s >= 0.0) line += " s=" + spinText(row.s);
  if (row.J >= 0.0) line += " J=" + spinText(row.J);
  return line;
}

bool ChannelRadiiModel::parseLine(const QString &line, ChannelRadiusRow &row) {
  const QString body = line.section('#', 0, 0).trimmed();
  const QStringList t = body.split(QRegExp("\\s+"), Qt::SkipEmptyParts);
  if (t.size() < 3) return false;
  bool okKey = false, okL = false, okR = false;
  row.pairKey = t.at(0).toInt(&okKey);
  row.l = t.at(1).toInt(&okL);
  row.radius = t.at(2).toDouble(&okR);
  if (!okKey || !okL || !okR || row.l < 0 || !(row.radius > 0.0)) return false;
  row.s = row.J = -1.0;
  for (int i = 3; i < t.size(); i++) {
    const QString tok = t.at(i);
    double v;
    if (tok.startsWith("s=") && row.s < 0.0 && parseSpin(tok.mid(2), v) && v >= 0.0)
      row.s = v;
    else if (tok.startsWith("J=") && row.J < 0.0 && parseSpin(tok.mid(2), v) && v >= 0.0)
      row.J = v;
    else
      return false;
  }
  return true;
}

bool ChannelRadiiModel::readBlock(QTextStream &in) {
  QList<ChannelRadiusRow> rows;
  while (!in.atEnd()) {
    const QString line = in.readLine();
    if (line.trimmed() == QString("</channelRadii>")) {
      setRows(rows);
      return true;
    }
    if (line.section('#', 0, 0).trimmed().isEmpty()) continue;
    ChannelRadiusRow row;
    if (!parseLine(line, row)) return false;
    rows.append(row);
  }
  return false;
}

void ChannelRadiiModel::writeBlock(QTextStream &out) const {
  for (const ChannelRadiusRow &r : rows_) out << lineFor(r) << Qt::endl;
}
