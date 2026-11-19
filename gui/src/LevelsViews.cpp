// Views of the Levels tab (docs/upgrade/qt_gui.md sections 11.2, 11.3): the spin box of J in halves, the parity
// toggle, the filter and sort of the levels, and the editor that hands a typed value to the tab for checking.
#include "LevelsViews.h"

#include <QEvent>
#include <QHBoxLayout>
#include <QKeyEvent>
#include <QLineEdit>
#include <QLocale>
#include <QRegularExpression>
#include <cmath>

#include "LevelsModel.h"

HalfSpinBox::HalfSpinBox(QWidget *parent) : QDoubleSpinBox(parent) {
  setDecimals(1);
  setSingleStep(0.5);
  setRange(0.0, 20.0);
  setAccelerated(true);
  setToolTip(tr("J in steps of 1/2: type 3/2 or 1.5, or use the arrows"));
}

QString HalfSpinBox::jText(double j) {
  const int twoJ = (int)std::lround(2.0 * j);
  return twoJ % 2 != 0 ? QString("%1/2").arg(twoJ) : QString::number(twoJ / 2);
}

QString HalfSpinBox::textFromValue(double value) const { return jText(value); }

double HalfSpinBox::valueFromText(const QString &text) const {
  const QString t = text.trimmed();
  const int slash = t.indexOf('/');
  bool ok = false;
  if (slash > 0) {
    const int num = t.left(slash).toInt(&ok);
    if (ok && t.mid(slash + 1).trimmed() == "2") return num / 2.0;
    return value();
  }
  const double v = QLocale::c().toDouble(t, &ok);
  return ok ? std::round(2.0 * v) / 2.0 : value();
}

QValidator::State HalfSpinBox::validate(QString &text, int &pos) const {
  Q_UNUSED(pos);
  const QRegularExpression whole("^\\s*\\d*\\s*$"), half("^\\s*\\d+\\s*/\\s*2?\\s*$"),
      decimal("^\\s*\\d*(\\.(0|5)?)?\\s*$");
  if (text.trimmed().isEmpty()) return QValidator::Intermediate;
  if (whole.match(text).hasMatch() || decimal.match(text).hasMatch()) return QValidator::Acceptable;
  if (half.match(text).hasMatch()) return text.trimmed().endsWith('2') ? QValidator::Acceptable : QValidator::Intermediate;
  return QValidator::Invalid;
}

ParityToggle::ParityToggle(QWidget *parent) : QToolButton(parent) {
  setCheckable(true);
  setText("+");
  setToolTip(tr("Parity: click to switch between + and -"));
  setMinimumWidth(28);
  connect(this, &QToolButton::toggled, this, [this](bool negative) { setText(negative ? QString(QChar(0x2212)) : "+"); });
}

// ---- the filter and the sort ----

LevelsFilterProxy::LevelsFilterProxy(QObject *parent) : QSortFilterProxyModel(parent) { setDynamicSortFilter(true); }

// Qt 6.10 replaced invalidateFilter by begin/endFilterChange
#if QT_VERSION >= QT_VERSION_CHECK(6, 10, 0)
#define LEVELS_FILTER_CHANGE(statement) \
  do { beginFilterChange(); statement; endFilterChange(); } while (0)
#else
#define LEVELS_FILTER_CHANGE(statement) \
  do { statement; invalidateFilter(); } while (0)
#endif

void LevelsFilterProxy::setNucleus(const QString &label) {
  LEVELS_FILTER_CHANGE(nucleus_ = label);
}

void LevelsFilterProxy::setJpi(int twoJ, int parity) {
  LEVELS_FILTER_CHANGE(twoJ_ = twoJ; parity_ = parity);
}

void LevelsFilterProxy::setSearch(const QString &text) {
  LEVELS_FILTER_CHANGE(search_ = text.trimmed());
}

void LevelsFilterProxy::setOnlyIncluded(bool on) {
  LEVELS_FILTER_CHANGE(onlyIncluded_ = on);
}

void LevelsFilterProxy::setHighlight(int twoJ, int parity) {
  if (twoJ == highlightTwoJ_ && parity == highlightParity_) return;
  highlightTwoJ_ = twoJ;
  highlightParity_ = parity;
  if (rowCount() > 0) emit dataChanged(index(0, 2), index(rowCount() - 1, 2), {Qt::FontRole});
}

QVariant LevelsFilterProxy::data(const QModelIndex &index, int role) const {
  if (role == Qt::FontRole && index.column() == 2 && highlightTwoJ_ >= 0) {
    const QModelIndex s = mapToSource(index);
    const QAbstractItemModel *m = sourceModel();
    const int twoJ = (int)std::lround(2.0 * m->data(m->index(s.row(), 2), Qt::EditRole).toDouble());
    const int pi = m->data(m->index(s.row(), 3), Qt::EditRole).toInt();
    if (twoJ == highlightTwoJ_ && pi == highlightParity_) {
      QFont bold;
      bold.setBold(true);
      return bold;
    }
  }
  return QSortFilterProxyModel::data(index, role);
}

bool LevelsFilterProxy::filterAcceptsRow(int row, const QModelIndex &parent) const {
  const QAbstractItemModel *m = sourceModel();
  if (!nucleus_.isEmpty() && m->data(m->index(row, LevelsModel::NucleusColumn, parent), Qt::DisplayRole).toString() != nucleus_)
    return false;
  const LevelsModel *levels = qobject_cast<const LevelsModel *>(m);
  if (onlyIncluded_ && levels && row < levels->getLevels().size() && levels->getLevels().at(row).isActive != 1) return false;
  const int twoJ = (int)std::lround(2.0 * m->data(m->index(row, 2, parent), Qt::EditRole).toDouble());
  const int pi = m->data(m->index(row, 3, parent), Qt::EditRole).toInt();
  if (twoJ_ >= 0 && (twoJ != twoJ_ || pi != parity_)) return false;
  if (!search_.isEmpty()) {
    // J^pi as written ("3/2-", with an ASCII sign) or the energy's digits
    const QString jpi = HalfSpinBox::jText(twoJ / 2.0) + (pi < 0 ? "-" : "+");
    const QString energy = QString::number(m->data(m->index(row, 4, parent), Qt::EditRole).toDouble(), 'g', 12);
    QString s = search_;
    s.replace(QChar(0x2212), '-');
    if (!jpi.startsWith(s, Qt::CaseInsensitive) && !energy.contains(s)) return false;
  }
  return true;
}

bool LevelsFilterProxy::lessThan(const QModelIndex &left, const QModelIndex &right) const {
  const QAbstractItemModel *m = sourceModel();
  auto value = [m](const QModelIndex &i, int column) { return m->data(m->index(i.row(), column), Qt::EditRole); };
  const double el = value(left, 4).toDouble(), er = value(right, 4).toDouble();
  if (left.column() == 2 || left.column() == 3) {
    const double jl = value(left, 2).toDouble(), jr = value(right, 2).toDouble();
    if (jl != jr) return jl < jr;
    const int pl = value(left, 3).toInt(), pr = value(right, 3).toInt();
    if (pl != pr) return pl > pr;  // + before -
    return el < er;
  }
  if (left.column() == 4) return el < er;
  return QSortFilterProxyModel::lessThan(left, right);
}

// ---- the editor ----

namespace {

/// The J^pi editor: J in halves and the parity toggle.
class JpiEditor : public QWidget {
 public:
  explicit JpiEditor(QWidget *parent) : QWidget(parent) {
    auto *h = new QHBoxLayout(this);
    h->setContentsMargins(0, 0, 0, 0);
    h->setSpacing(1);
    j = new HalfSpinBox;
    pi = new ParityToggle;
    h->addWidget(j, 1);
    h->addWidget(pi);
    setFocusProxy(j);
    setAutoFillBackground(true);
  }
  HalfSpinBox *j;
  ParityToggle *pi;
};

}  // namespace

CheckedEditDelegate::CheckedEditDelegate(Kind kind, Commit commit, QObject *parent) :
  QStyledItemDelegate(parent), kind_(kind), commit_(std::move(commit)) {}

QWidget *CheckedEditDelegate::createEditor(QWidget *parent, const QStyleOptionViewItem &option,
                                           const QModelIndex &index) const {
  Q_UNUSED(option);
  QWidget *w;
  if (kind_ == Kind::Jpi) {
    // the view filters the editor's events through this delegate; the spin box leaves Enter to its parent
    w = new JpiEditor(parent);
  } else {
    auto *e = new QLineEdit(parent);
    e->setAlignment(Qt::AlignCenter);
    w = e;
  }
  w->setProperty("cellIndex", QVariant::fromValue(QModelIndex(index)));
  return w;
}

void CheckedEditDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const {
  if (kind_ == Kind::Jpi) {
    auto *e = static_cast<JpiEditor *>(editor);
    const QAbstractItemModel *m = index.model();
    e->j->setValue(m->data(index.sibling(index.row(), 2), Qt::EditRole).toDouble());
    e->pi->setParity(m->data(index.sibling(index.row(), 3), Qt::EditRole).toInt());
  } else {
    // the number as held, in full (the shortest text that reads back to it)
    const double v = index.data(Qt::EditRole).toDouble();
    static_cast<QLineEdit *>(editor)->setText(QString::number(v, 'g', QLocale::FloatingPointShortest));
  }
}

void CheckedEditDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const {
  Q_UNUSED(model);
  QVariant value;
  if (kind_ == Kind::Jpi) {
    auto *e = static_cast<JpiEditor *>(editor);
    e->j->interpretText();
    value = QVariantList{e->j->value(), e->pi->parity()};
  } else {
    value = static_cast<QLineEdit *>(editor)->text();
  }
  QString why;
  if (!commit_(index, value, why)) emit refused(index, why);
}

bool CheckedEditDelegate::eventFilter(QObject *object, QEvent *event) {
  if (event->type() == QEvent::KeyPress) {
    auto *key = static_cast<QKeyEvent *>(event);
    if (key->key() == Qt::Key_Return || key->key() == Qt::Key_Enter) {
      QWidget *editor = qobject_cast<QWidget *>(object);
      if (editor && editor->property("cellIndex").isValid()) {
        const QModelIndex index = editor->property("cellIndex").toModelIndex();
        emit commitData(editor);
        emit closeEditor(editor, QAbstractItemDelegate::NoHint);
        emit moveDown(index);
        return true;
      }
    }
  }
  return QStyledItemDelegate::eventFilter(object, event);
}
