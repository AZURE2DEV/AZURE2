#include "ScatteringLengthsModel.h"

#include <QBrush>
#include <QColor>
#include <QComboBox>
#include <QLineEdit>
#include <QLocale>
#include <QRegularExpression>
#include <QStringList>
#include <cmath>

namespace {

QString number(double v) { return QString::number(v, 'g', QLocale::FloatingPointShortest); }

bool sameSpin(double a, double b) { return std::fabs(a - b) < 1.0e-6; }

// "1", "1/2", "3/2" ...
QString plainSpin(double s) {
  int twice = (int)std::floor(2.0 * s + 0.5);
  if (twice % 2 == 0) return QString::number(twice / 2);
  return QString("%1/2").arg(twice);
}

}  // namespace

ScatteringLengthsModel::ScatteringLengthsModel(QObject *parent) : QAbstractTableModel(parent) {}

void ScatteringLengthsModel::setPairsModel(PairsModel *model) {
  if (pairsModel) disconnect(pairsModel, 0, this, 0);
  pairsModel = model;
  if (pairsModel) {
    connect(pairsModel, SIGNAL(rowsInserted(QModelIndex, int, int)), this, SLOT(pairsChanged()));
    connect(pairsModel, SIGNAL(rowsRemoved(QModelIndex, int, int)), this, SLOT(pairsChanged()));
    connect(pairsModel, SIGNAL(dataChanged(QModelIndex, QModelIndex)), this, SLOT(pairsChanged()));
    connect(pairsModel, SIGNAL(modelReset()), this, SLOT(pairsChanged()));
  }
  pairsChanged();
}

void ScatteringLengthsModel::pairsChanged() {
  if (!lines.isEmpty()) emit dataChanged(index(0, 0), index(lines.size() - 1, NumCols - 1));
}

int ScatteringLengthsModel::rowCount(const QModelIndex &parent) const {
  Q_UNUSED(parent);
  return lines.size();
}

int ScatteringLengthsModel::columnCount(const QModelIndex &parent) const {
  Q_UNUSED(parent);
  return NumCols;
}

// -- pairs and channel spins --------------------------------------------------

bool ScatteringLengthsModel::isIdentical(const PairsData &p) {
  return p.pairType == 0 && p.lightZ == p.heavyZ && std::fabs(p.lightM - p.heavyM) < 1.0e-6 &&
         std::fabs(p.lightJ - p.heavyJ) < 1.0e-6 && p.lightPi == p.heavyPi && std::fabs(p.excitationEnergy) < 1.0e-9;
}

QList<double> ScatteringLengthsModel::allowedChannelSpins(const PairsData &p) {
  QList<double> out;
  if (p.pairType != 0) return out;
  const bool identical = isIdentical(p);
  for (double s = std::fabs(p.lightJ - p.heavyJ); s <= p.lightJ + p.heavyJ + 1.0e-6; s += 1.0) {
    // Bose/Fermi symmetry of an identical pair: l + s even, so s even for l = 0.
    if (identical && ((int)(s + 0.5)) % 2 != 0) continue;
    out << s;
  }
  return out;
}

QString ScatteringLengthsModel::spinLabel(double s) {
  if (s < 0.0) return tr("coherent");
  return plainSpin(s);
}

QString ScatteringLengthsModel::pairLabel(int key) const {
  QList<PairsData> pairs = pairsModel ? pairsModel->getPairs() : QList<PairsData>();
  if (key < 1 || key > pairs.size()) return tr("%1: (no such pair)").arg(key);
  const PairsData &p = pairs.at(key - 1);
  QString label = QString("%1: %2").arg(key).arg(pairsModel->getPairLabel(p));
  if (p.pairType != 0) label += tr(" (not a particle pair)");
  return label;
}

QString ScatteringLengthsModel::problem(int row) const {
  if (row < 0 || row >= lines.size()) return QString();
  return problem(lines.at(row));
}

QString ScatteringLengthsModel::problem(const ScatteringLengthEntry &e) const {
  if (!pairsModel) return QString();
  QList<PairsData> pairs = pairsModel->getPairs();
  if (e.pairKey < 1 || e.pairKey > pairs.size())
    return tr("There is no pair %1 in the Particle Pairs tab; the calculation will stop on this line.").arg(e.pairKey);
  const PairsData &p = pairs.at(e.pairKey - 1);
  if (p.pairType != 0)
    return tr("Pair %1 is not a particle pair; a scattering length needs two particles.").arg(e.pairKey);
  if (e.s < 0.0) {
    if (isIdentical(p))
      return tr("The coherent combination is not defined for an identical pair; choose a channel spin.");
    return QString();
  }
  QList<double> allowed = allowedChannelSpins(p);
  for (double s : allowed)
    if (sameSpin(s, e.s)) return QString();
  QStringList names;
  for (double s : allowed) names << plainSpin(s);
  return tr("Channel spin %1 is not allowed for pair %2 (allowed: %3).").arg(plainSpin(e.s)).arg(e.pairKey).arg(names.join(", "));
}

// -- the table ------------------------------------------------------------------

QVariant ScatteringLengthsModel::data(const QModelIndex &index, int role) const {
  if (!index.isValid() || index.row() < 0 || index.row() >= lines.size()) return QVariant();
  const ScatteringLengthEntry &e = lines.at(index.row());
  const int col = index.column();
  if (role == Qt::DisplayRole) {
    switch (col) {
      case PairCol:
        return pairLabel(e.pairKey);
      case SpinCol: {
        if (e.s < 0.0) return spinLabel(e.s);
        // J = s for l = 0, with the parity of the pair: show it as J^pi.
        QList<PairsData> pairs = pairsModel ? pairsModel->getPairs() : QList<PairsData>();
        if (e.pairKey >= 1 && e.pairKey <= pairs.size() && pairs.at(e.pairKey - 1).pairType == 0) {
          const PairsData &p = pairs.at(e.pairKey - 1);
          return plainSpin(e.s) + QString(QChar(p.lightPi * p.heavyPi > 0 ? 0x207A : 0x207B));
        }
        return plainSpin(e.s);
      }
      case ValueCol:
        return e.value;
      case SigmaCol:
        return e.sigma;
      case ImValueCol:
        return e.hasIm ? QVariant(e.imValue) : QVariant(QString());
      case ImSigmaCol:
        return e.hasIm ? QVariant(e.imSigma) : QVariant(QString());
      default:
        return QVariant();
    }
  } else if (role == Qt::EditRole) {
    switch (col) {
      case PairCol:
        return e.pairKey;
      case SpinCol:
        return e.s;
      case ValueCol:
        return e.value;
      case SigmaCol:
        return e.sigma;
      case ImValueCol:
        return e.imValue;
      case ImSigmaCol:
        return e.imSigma;
      default:
        return QVariant();
    }
  } else if (role == Qt::CheckStateRole && col == ActiveCol) {
    return e.active ? Qt::Checked : Qt::Unchecked;
  } else if (role == Qt::TextAlignmentRole) {
    return col == PairCol ? QVariant(Qt::AlignLeft | Qt::AlignVCenter) : QVariant(Qt::AlignCenter);
  } else if (role == Qt::ForegroundRole) {
    if (!problem(e).isEmpty()) return QBrush(QColor(200, 0, 0));
  } else if (role == Qt::ToolTipRole) {
    QString why = problem(e);
    if (!why.isEmpty()) return why;
    switch (col) {
      case ActiveCol:
        return tr("Include this measurement in chi2.");
      case PairCol:
        return tr("The particle pair whose s-wave scattering length was measured.");
      case SpinCol:
        return tr("Channel spin s (J = s, parity pi1*pi2 for l = 0), or \"coherent\": the spin-weighted sum "
                  "of the channel spins, sum (2s+1)/((2j1+1)(2j2+1)) a_s.");
      case ValueCol:
        return tr("Measured Re a0 (fm).");
      case SigmaCol:
        return tr("Uncertainty of Re a0 (fm); 0 = Re a0 is not constrained.");
      case ImValueCol:
        return tr("Measured Im a0 (fm), when the pair has open reaction channels at threshold; blank = none.");
      case ImSigmaCol:
        return tr("Uncertainty of Im a0 (fm); 0 or blank = Im a0 is not constrained.");
      default:
        return QVariant();
    }
  }
  return QVariant();
}

QVariant ScatteringLengthsModel::headerData(int section, Qt::Orientation orientation, int role) const {
  if (role != Qt::DisplayRole) return QVariant();
  if (orientation == Qt::Vertical) return section + 1;
  switch (section) {
    case ActiveCol:
      return QString();
    case PairCol:
      return tr("Pair");
    case SpinCol:
      return tr("Spin s (J%1)").arg(QChar(0x03C0));
    case ValueCol:
      return tr("Re a0 (fm)");
    case SigmaCol:
      return tr("Re %1 (fm)").arg(QChar(0x03C3));
    case ImValueCol:
      return tr("Im a0 (fm)");
    case ImSigmaCol:
      return tr("Im %1 (fm)").arg(QChar(0x03C3));
    default:
      return QVariant();
  }
}

bool ScatteringLengthsModel::setData(const QModelIndex &index, const QVariant &value, int role) {
  if (!index.isValid() || index.row() < 0 || index.row() >= lines.size()) return false;
  ScatteringLengthEntry &e = lines[index.row()];
  const int col = index.column();
  if (role == Qt::CheckStateRole && col == ActiveCol) {
    e.active = (value.toInt() == Qt::Checked);
    emit dataChanged(index, index);
    return true;
  }
  if (role != Qt::EditRole) return false;
  if (col == ActiveCol) {
    e.active = value.toBool();
  } else if (col == PairCol) {
    bool ok;
    int key = value.toInt(&ok);
    if (!ok) return false;
    e.pairKey = key;
    // Keep the channel spin when the new pair allows it, otherwise take its first one.
    QList<PairsData> pairs = pairsModel ? pairsModel->getPairs() : QList<PairsData>();
    if (key >= 1 && key <= pairs.size() && pairs.at(key - 1).pairType == 0) {
      const PairsData &p = pairs.at(key - 1);
      QList<double> allowed = allowedChannelSpins(p);
      bool fine = (e.s < 0.0 && !isIdentical(p));
      for (double s : allowed) fine |= sameSpin(s, e.s);
      if (!fine && !allowed.isEmpty()) e.s = allowed.first();
    }
    emit dataChanged(this->index(index.row(), 0), this->index(index.row(), NumCols - 1));
    return true;
  } else {
    // Numbers; a blank imaginary cell removes that part.
    QString text = value.toString().trimmed();
    bool ok = true;
    double v = text.isEmpty() ? 0.0 : value.toDouble(&ok);
    if (!ok) v = QLocale().toDouble(text, &ok);
    if (!ok) return false;
    if (text.isEmpty() && col != ImValueCol && col != ImSigmaCol) return false;
    switch (col) {
      case SpinCol:
        e.s = (v < 0.0) ? -1.0 : v;
        break;
      case ValueCol:
        e.value = v;
        break;
      case SigmaCol:
        if (v < 0.0) return false;
        e.sigma = v;
        break;
      case ImValueCol:
      case ImSigmaCol:
        if (col == ImSigmaCol && v < 0.0) return false;
        (col == ImValueCol ? e.imValue : e.imSigma) = v;
        e.hasIm = !text.isEmpty() || e.imValue != 0.0 || e.imSigma != 0.0;
        break;
      default:
        return false;
    }
  }
  emit dataChanged(this->index(index.row(), 0), this->index(index.row(), NumCols - 1));
  return true;
}

Qt::ItemFlags ScatteringLengthsModel::flags(const QModelIndex &index) const {
  if (!index.isValid()) return Qt::NoItemFlags;
  Qt::ItemFlags f = Qt::ItemIsEnabled | Qt::ItemIsSelectable;
  if (index.column() == ActiveCol) return f | Qt::ItemIsUserCheckable;
  return f | Qt::ItemIsEditable;
}

bool ScatteringLengthsModel::removeRows(int position, int rows, const QModelIndex &parent) {
  if (rows <= 0 || position < 0 || position + rows > lines.size()) return false;
  beginRemoveRows(parent, position, position + rows - 1);
  for (int i = 0; i < rows; i++) lines.removeAt(position);
  endRemoveRows();
  return true;
}

void ScatteringLengthsModel::addEntry(const ScatteringLengthEntry &entry) {
  beginInsertRows(QModelIndex(), lines.size(), lines.size());
  lines.append(entry);
  endInsertRows();
}

void ScatteringLengthsModel::addDefaultEntry() {
  ScatteringLengthEntry e;
  QList<PairsData> pairs = pairsModel ? pairsModel->getPairs() : QList<PairsData>();
  for (int k = 0; k < pairs.size(); k++) {
    if (pairs.at(k).pairType != 0) continue;
    e.pairKey = k + 1;
    QList<double> allowed = allowedChannelSpins(pairs.at(k));
    if (!allowed.isEmpty()) e.s = allowed.first();
    break;
  }
  addEntry(e);
}

void ScatteringLengthsModel::clear() {
  beginResetModel();
  lines.clear();
  endResetModel();
}

// -- the block ------------------------------------------------------------------

bool ScatteringLengthsModel::read(QTextStream &in) {
  while (!in.atEnd()) {
    QString line = in.readLine();
    if (line.contains("</scatteringLengths>")) return true;
    QString entry = line.section('#', 0, 0).trimmed();
    if (entry.isEmpty()) continue;
    QStringList t = entry.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
    bool ok;
    int active = t.at(0).toInt(&ok);
    if (!ok) continue;  // as the engine: a line that does not start with a number is skipped
    if (t.size() != 5 && t.size() != 7) return false;
    ScatteringLengthEntry e;
    e.active = (active != 0);
    bool good = true;
    e.pairKey = t.at(1).toInt(&ok);
    good &= ok;
    e.s = t.at(2).toDouble(&ok);
    good &= ok;
    e.value = t.at(3).toDouble(&ok);
    good &= ok;
    e.sigma = t.at(4).toDouble(&ok);
    good &= ok;
    if (t.size() == 7) {
      e.hasIm = true;
      e.imValue = t.at(5).toDouble(&ok);
      good &= ok;
      e.imSigma = t.at(6).toDouble(&ok);
      good &= ok;
    }
    if (!good) return false;
    addEntry(e);
  }
  return false;
}

bool ScatteringLengthsModel::write(QTextStream &out) const {
  out << "# active pairKey s value sigma [imValue imSigma]   (fm; s = -1: coherent combination)" << Qt::endl;
  for (const ScatteringLengthEntry &e : lines) {
    out << (e.active ? 1 : 0) << " " << e.pairKey << " " << number(e.s) << " " << number(e.value) << " " << number(e.sigma);
    if (e.hasIm || e.imValue != 0.0 || e.imSigma != 0.0) out << " " << number(e.imValue) << " " << number(e.imSigma);
    out << Qt::endl;
  }
  return true;
}

// -- the editors ------------------------------------------------------------------

QWidget *ScatteringLengthsDelegate::createEditor(QWidget *parent, const QStyleOptionViewItem &option,
                                                 const QModelIndex &index) const {
  const int col = index.column();
  if (col == ScatteringLengthsModel::ActiveCol) return nullptr;
  if (col == ScatteringLengthsModel::PairCol || col == ScatteringLengthsModel::SpinCol) {
    QComboBox *box = new QComboBox(parent);
    box->setSizeAdjustPolicy(QComboBox::AdjustToContents);
    ScatteringLengthsDelegate *self = const_cast<ScatteringLengthsDelegate *>(this);
    connect(box, QOverload<int>::of(&QComboBox::activated), self, [self, box](int) {
      emit self->commitData(box);
      emit self->closeEditor(box);
    });
    return box;
  }
  Q_UNUSED(option);
  return new QLineEdit(parent);
}

void ScatteringLengthsDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const {
  const ScatteringLengthsModel *model = qobject_cast<const ScatteringLengthsModel *>(index.model());
  PairsModel *pairsModel = model ? model->getPairsModel() : nullptr;
  QList<PairsData> pairs = pairsModel ? pairsModel->getPairs() : QList<PairsData>();
  const int col = index.column();
  if (QComboBox *box = qobject_cast<QComboBox *>(editor)) {
    box->clear();
    if (col == ScatteringLengthsModel::PairCol) {
      int key = index.data(Qt::EditRole).toInt();
      for (int k = 1; k <= pairs.size(); k++)
        if (pairs.at(k - 1).pairType == 0) box->addItem(QString("%1: %2").arg(k).arg(pairsModel->getPairLabel(pairs.at(k - 1))), k);
      if (box->findData(key) < 0) box->addItem(index.data(Qt::DisplayRole).toString(), key);
      box->setCurrentIndex(box->findData(key));
    } else {
      double s = index.data(Qt::EditRole).toDouble();
      int key = model ? model->index(index.row(), ScatteringLengthsModel::PairCol).data(Qt::EditRole).toInt() : 0;
      if (key >= 1 && key <= pairs.size() && pairs.at(key - 1).pairType == 0) {
        const PairsData &p = pairs.at(key - 1);
        const QString parity(QChar(p.lightPi * p.heavyPi > 0 ? 0x207A : 0x207B));
        for (double a : ScatteringLengthsModel::allowedChannelSpins(p))
          box->addItem(ScatteringLengthsModel::spinLabel(a) + parity, a);
        if (!ScatteringLengthsModel::isIdentical(p)) box->addItem(ScatteringLengthsModel::spinLabel(-1.0), -1.0);
      }
      int current = -1;
      for (int i = 0; i < box->count(); i++)
        if (std::fabs(box->itemData(i).toDouble() - s) < 1.0e-6) current = i;
      if (current < 0) {
        box->addItem(ScatteringLengthsDelegate::tr("%1 (not allowed)").arg(ScatteringLengthsModel::spinLabel(s)), s);
        current = box->count() - 1;
      }
      box->setCurrentIndex(current);
    }
    return;
  }
  if (QLineEdit *edit = qobject_cast<QLineEdit *>(editor)) {
    QVariant shown = index.data(Qt::DisplayRole);
    edit->setText(shown.toString().isEmpty() ? QString() : number(index.data(Qt::EditRole).toDouble()));
  }
}

void ScatteringLengthsDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const {
  if (QComboBox *box = qobject_cast<QComboBox *>(editor)) {
    if (box->currentIndex() >= 0) model->setData(index, box->currentData(), Qt::EditRole);
    return;
  }
  if (QLineEdit *edit = qobject_cast<QLineEdit *>(editor)) model->setData(index, edit->text(), Qt::EditRole);
}
