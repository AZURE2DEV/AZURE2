#include "AbundancesModel.h"

AbundancesModel::AbundancesModel(QObject *parent) : QAbstractTableModel(parent) {
  sampleOrder << QString();
}

void AbundancesModel::setPairsModel(PairsModel *model) {
  pairsModel = model;
  if (pairsModel) {
    connect(pairsModel, SIGNAL(rowsInserted(QModelIndex, int, int)), this, SLOT(sync()));
    connect(pairsModel, SIGNAL(rowsRemoved(QModelIndex, int, int)), this, SLOT(sync()));
    connect(pairsModel, SIGNAL(dataChanged(QModelIndex, QModelIndex)), this, SLOT(sync()));
    connect(pairsModel, SIGNAL(modelReset()), this, SLOT(sync()));
  }
  sync();
}

void AbundancesModel::sync() {
  QList<Nucleus> now = pairsModel ? pairsModel->nuclei() : QList<Nucleus>();
  if (now == rows) return;
  beginResetModel();
  rows = now;
  endResetModel();
}

int AbundancesModel::rowCount(const QModelIndex &parent) const {
  Q_UNUSED(parent);
  return rows.size();
}

int AbundancesModel::columnCount(const QModelIndex &parent) const {
  Q_UNUSED(parent);
  return NumCols;
}

// -- samples ------------------------------------------------------------------

QStringList AbundancesModel::samples() const { return sampleOrder; }

void AbundancesModel::setCurrentSample(const QString &sample) {
  if (!sampleOrder.contains(sample)) return;
  if (sample == current) return;
  beginResetModel();
  current = sample;
  endResetModel();
  emit samplesChanged();
}

bool AbundancesModel::addSample(const QString &sample) {
  QString s = sample.trimmed();
  if (s.isEmpty() || s.contains(':') || s.contains('#') || s.contains(QRegExp("\\s"))) return false;
  if (sampleOrder.contains(s)) return false;
  sampleOrder << s;
  emit samplesChanged();
  return true;
}

bool AbundancesModel::removeSample(const QString &sample) {
  if (sample.isEmpty() || !sampleOrder.contains(sample)) return false;
  sampleOrder.removeAll(sample);
  stored.remove(sample);
  sumSigma.remove(sample);
  if (current == sample) {
    beginResetModel();
    current = QString();
    endResetModel();
  }
  emit samplesChanged();
  return true;
}

double AbundancesModel::sampleSumSigma(const QString &sample) const {
  return sumSigma.value(sample, -1.0);
}

void AbundancesModel::setSampleSumSigma(const QString &sample, double sigma) {
  if (!sampleOrder.contains(sample)) sampleOrder << sample;
  if (sigma < 0.0)
    sumSigma.remove(sample);
  else
    sumSigma[sample] = sigma;
}

// -- entries ------------------------------------------------------------------

AbundanceEntry AbundancesModel::entry(const QString &sample, const Nucleus &nucleus) const {
  if (stored.contains(sample) && stored.value(sample).contains(nucleus)) return stored.value(sample).value(nucleus);
  AbundanceEntry e;
  e.nucleus = nucleus;
  e.sample = sample;
  return e;
}

QList<AbundanceEntry> AbundancesModel::entries(const QString &sample) const {
  QList<AbundanceEntry> out;
  for (int i = 0; i < rows.size(); i++) out.append(entry(sample, rows.at(i)));
  return out;
}

QList<AbundanceEntry> AbundancesModel::allEntries() const {
  QList<AbundanceEntry> out;
  for (int s = 0; s < sampleOrder.size(); s++) out += entries(sampleOrder.at(s));
  return out;
}

bool AbundancesModel::setEntry(const QString &sample, const Nucleus &nucleus, double value, bool vary, double sigma) {
  if (value < 0.0 || sigma < 0.0) return false;
  if (!sampleOrder.contains(sample)) {
    if (!addSample(sample)) return false;
  }
  AbundanceEntry e;
  e.nucleus = nucleus;
  e.value = value;
  e.vary = vary;
  e.sigma = sigma;
  e.sample = sample;
  stored[sample][nucleus] = e;
  int row = rows.indexOf(nucleus);
  if (row >= 0 && sample == current) emit dataChanged(index(row, ValueCol), index(row, SigmaCol));
  return true;
}

QString AbundancesModel::parameterName(const AbundanceEntry &e) {
  if (e.sample.isEmpty()) return QString("abundance_%1_%2").arg(e.nucleus.Z).arg(e.nucleus.A);
  return QString("abundance_%1_%2_%3").arg(e.sample).arg(e.nucleus.Z).arg(e.nucleus.A);
}

bool AbundancesModel::isDefault() const {
  if (sampleOrder.size() > 1 || !sumSigma.isEmpty()) return false;
  for (int i = 0; i < rows.size(); i++) {
    AbundanceEntry e = entry(QString(), rows.at(i));
    if (e.value != 1.0 || e.vary || e.sigma != 0.0) return false;
  }
  return true;
}

void AbundancesModel::clear() {
  beginResetModel();
  stored.clear();
  sumSigma.clear();
  sampleOrder = QStringList() << QString();
  current = QString();
  rows.clear();
  endResetModel();
  emit samplesChanged();
  sync();
}

// -- the table ------------------------------------------------------------------

QVariant AbundancesModel::data(const QModelIndex &index, int role) const {
  if (!index.isValid() || index.row() < 0 || index.row() >= rows.size()) return QVariant();
  const Nucleus &nucleus = rows.at(index.row());
  AbundanceEntry e = entry(current, nucleus);
  if (role == Qt::DisplayRole || role == Qt::EditRole) {
    switch (index.column()) {
      case NucleusCol:
        return pairsModel ? pairsModel->nucleusLongLabel(nucleus) : PairsModel::nucleusLabel(nucleus, false);
      case ValueCol:
        return e.value;
      case SigmaCol:
        return e.sigma;
      default:
        return QVariant();
    }
  } else if (role == Qt::CheckStateRole && index.column() == VaryCol) {
    return e.vary ? Qt::Checked : Qt::Unchecked;
  } else if (role == Qt::TextAlignmentRole) {
    return Qt::AlignCenter;
  } else if (role == Qt::ToolTipRole) {
    switch (index.column()) {
      case NucleusCol:
        return tr("The compound nucleus and the pair that forms it; the abundance is that of the pair's target isotope.");
      case ValueCol:
        return tr("Isotopic (or stoichiometric) abundance of this nucleus in the sample shown.");
      case VaryCol:
        return tr("Fit the abundance (parameter %1).").arg(parameterName(e));
      case SigmaCol:
        return tr("Width of the Gaussian prior on a fitted abundance; 0 = no prior.");
      default:
        return QVariant();
    }
  }
  return QVariant();
}

QVariant AbundancesModel::headerData(int section, Qt::Orientation orientation, int role) const {
  if (role != Qt::DisplayRole) return QVariant();
  if (orientation == Qt::Vertical) return section + 1;
  switch (section) {
    case NucleusCol:
      return tr("Compound Nucleus\n(target pair)");
    case ValueCol:
      return tr("Abundance");
    case VaryCol:
      return tr("Fit?");
    case SigmaCol:
      return tr("Prior\nSigma");
    default:
      return QVariant();
  }
}

bool AbundancesModel::setData(const QModelIndex &index, const QVariant &value, int role) {
  if (!index.isValid() || index.row() < 0 || index.row() >= rows.size()) return false;
  const Nucleus &nucleus = rows.at(index.row());
  AbundanceEntry e = entry(current, nucleus);
  if (role == Qt::CheckStateRole && index.column() == VaryCol) {
    e.vary = (value.toInt() == Qt::Checked);
  } else if (role == Qt::EditRole && index.column() == ValueCol) {
    bool ok;
    double v = value.toDouble(&ok);
    if (!ok || v < 0.0) return false;
    e.value = v;
  } else if (role == Qt::EditRole && index.column() == SigmaCol) {
    bool ok;
    double v = value.toDouble(&ok);
    if (!ok || v < 0.0) return false;
    e.sigma = v;
  } else {
    return false;
  }
  e.sample = current;
  stored[current][nucleus] = e;
  emit dataChanged(index, index);
  return true;
}

Qt::ItemFlags AbundancesModel::flags(const QModelIndex &index) const {
  if (!index.isValid()) return Qt::NoItemFlags;
  Qt::ItemFlags f = Qt::ItemIsEnabled | Qt::ItemIsSelectable;
  if (index.column() == VaryCol) return f | Qt::ItemIsUserCheckable;
  if (index.column() == ValueCol || index.column() == SigmaCol) return f | Qt::ItemIsEditable;
  return f;
}
