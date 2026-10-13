#include "AbundancesModel.h"

AbundancesModel::AbundancesModel(QObject *parent) : QAbstractTableModel(parent) {}

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

AbundanceEntry AbundancesModel::entry(const Nucleus &nucleus) const {
  if (stored.contains(nucleus)) return stored.value(nucleus);
  AbundanceEntry e;
  e.nucleus = nucleus;
  return e;
}

QList<AbundanceEntry> AbundancesModel::entries() const {
  QList<AbundanceEntry> out;
  for (int i = 0; i < rows.size(); i++) out.append(entry(rows.at(i)));
  return out;
}

bool AbundancesModel::setEntry(const Nucleus &nucleus, double value, bool vary, double sigma) {
  if (value < 0.0 || sigma < 0.0) return false;
  AbundanceEntry e;
  e.nucleus = nucleus;
  e.value = value;
  e.vary = vary;
  e.sigma = sigma;
  stored[nucleus] = e;
  int row = rows.indexOf(nucleus);
  if (row >= 0) emit dataChanged(index(row, ValueCol), index(row, SigmaCol));
  return true;
}

bool AbundancesModel::isDefault() const {
  for (int i = 0; i < rows.size(); i++) {
    AbundanceEntry e = entry(rows.at(i));
    if (e.value != 1.0 || e.vary || e.sigma != 0.0) return false;
  }
  return true;
}

void AbundancesModel::clear() {
  beginResetModel();
  stored.clear();
  rows.clear();
  endResetModel();
  sync();
}

QVariant AbundancesModel::data(const QModelIndex &index, int role) const {
  if (!index.isValid() || index.row() < 0 || index.row() >= rows.size()) return QVariant();
  const Nucleus &nucleus = rows.at(index.row());
  AbundanceEntry e = entry(nucleus);
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
        return tr("Isotopic abundance (weight) of this nucleus in a natural-target segment.");
      case VaryCol:
        return tr("Fit the abundance (parameter abundance_<Z>_<A>).");
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
  AbundanceEntry e = entry(nucleus);
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
  stored[nucleus] = e;
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
