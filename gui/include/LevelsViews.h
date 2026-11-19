#ifndef LEVELSVIEWS_H
#define LEVELSVIEWS_H

#include <functional>
#include <QDoubleSpinBox>
#include <QSortFilterProxyModel>
#include <QStyledItemDelegate>
#include <QToolButton>

/*!
 * A spin box of J in steps of 1/2, written as the level scheme writes it ("3/2", "2"); typing "3/2" or "1.5" works.
 */
class HalfSpinBox : public QDoubleSpinBox {
  Q_OBJECT
 public:
  explicit HalfSpinBox(QWidget *parent = nullptr);
  QString textFromValue(double value) const override;
  double valueFromText(const QString &text) const override;
  QValidator::State validate(QString &text, int &pos) const override;
  /// "3/2", "2": J as text (the labels of the Levels tab).
  static QString jText(double j);
};

/// A +/- toggle of a parity: checked = negative.
class ParityToggle : public QToolButton {
  Q_OBJECT
 public:
  explicit ParityToggle(QWidget *parent = nullptr);
  int parity() const { return isChecked() ? -1 : 1; }
  void setParity(int pi) { setChecked(pi < 0); }
};

/*!
 * The levels as the tab shows them (docs/upgrade/qt_gui.md section 11.3): filtered by compound nucleus, J^pi, a search
 * text (J^pi as written, digits of the energy) and "only included"; sorted by energy or by J^pi (J, parity, then
 * energy, so a J^pi group stands together); the J^pi of the current level's group in bold.
 */
class LevelsFilterProxy : public QSortFilterProxyModel {
  Q_OBJECT
 public:
  explicit LevelsFilterProxy(QObject *parent = nullptr);
  /// The nucleus label (LevelsModel::NucleusColumn's display text), "" for every nucleus.
  void setNucleus(const QString &label);
  /// One J^pi (twoJ >= 0), or every J^pi (twoJ < 0).
  void setJpi(int twoJ, int parity);
  void setSearch(const QString &text);
  void setOnlyIncluded(bool on);
  /// The J^pi shown in bold (twoJ < 0: none).
  void setHighlight(int twoJ, int parity);
  QVariant data(const QModelIndex &index, int role) const override;

 protected:
  bool filterAcceptsRow(int row, const QModelIndex &parent) const override;
  bool lessThan(const QModelIndex &left, const QModelIndex &right) const override;

 private:
  QString nucleus_, search_;
  int twoJ_ = -1, parity_ = 0;
  int highlightTwoJ_ = -1, highlightParity_ = 0;
  bool onlyIncluded_ = false;
};

/*!
 * An editor of a cell whose value the tab checks before it is set (an energy, a width, a J^pi): the delegate hands the
 * typed value to `commit`, which makes the edit (through the tab, with its undo entry) or refuses it with the reason;
 * a refusal is shown at the cell and leaves the model as it was.  Enter commits and moves one row down.
 */
class CheckedEditDelegate : public QStyledItemDelegate {
  Q_OBJECT
 public:
  using Commit = std::function<bool(const QModelIndex &index, const QVariant &value, QString &why)>;
  enum class Kind { Number, Jpi };
  CheckedEditDelegate(Kind kind, Commit commit, QObject *parent = nullptr);
  QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const override;
  void setEditorData(QWidget *editor, const QModelIndex &index) const override;
  void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const override;
  bool eventFilter(QObject *object, QEvent *event) override;

 signals:
  void refused(const QModelIndex &index, const QString &why) const;
  void moveDown(const QModelIndex &index) const;

 private:
  Kind kind_;
  Commit commit_;
};

#endif
