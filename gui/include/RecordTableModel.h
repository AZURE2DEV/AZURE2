#ifndef RECORDTABLEMODEL_H
#define RECORDTABLEMODEL_H

#include <QAbstractTableModel>
#include <QStringList>
#include <QStyledItemDelegate>
#include <vector>

#include "azure2/front/Engine.h"

class ProjectHolder;

/*!
 * The records of one section (of some of its kinds) as a table with real columns (docs/upgrade/qt_gui.md section
 * 11.5): the include switch of the kind (KindInfo::includeField) as a check box, the name, the kind when there are
 * several, then the summary fields of the kinds (KindInfo::summaryFields), each cell the field's short text
 * (ShortValueText), a nuisance followed by a "fit" check box; the whole record is the row's tooltip.  Cells of an
 * enum, unit or reference field, of a number and of a nuisance's value are edited in place; every edit is written
 * through the facade (only the changed line moves) and checked by the engine's parser first.
 */
class RecordTableModel : public QAbstractTableModel {
  Q_OBJECT

 public:
  struct Column {
    enum Type { Include, Name, Kind, Field, Fit };
    Type type = Field;
    std::string key;
    azure2::front::FieldInfo field;
    QString title;
  };

  RecordTableModel(ProjectHolder *holder, const QString &section, const QStringList &kinds, QObject *parent = nullptr);

  int rowCount(const QModelIndex &parent = QModelIndex()) const override;
  int columnCount(const QModelIndex &parent = QModelIndex()) const override;
  QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override;
  QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;
  Qt::ItemFlags flags(const QModelIndex &index) const override;
  bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole) override;

  /// Read the section again.
  void reload();
  const std::vector<azure2::front::Record> &records() const { return records_; }
  const Column &column(int c) const { return columns_[c]; }
  int columnOf(Column::Type type, const std::string &key = std::string()) const;
  int rowOf(const QString &name) const;
  /// The words a cell offers in place (enum, unit, reference), "" first for an optional field; empty for a free cell.
  QStringList choices(const QModelIndex &index) const;
  /// The text a cell's editor starts with (a nuisance: its number; else the value as written).
  QString editText(const QModelIndex &index) const;
  /// Switch records on or off (the include field); false with the engine's reason.
  bool setIncluded(const QStringList &names, bool on, QString &why);
  bool included(int row) const;
  bool writing() const { return writing_; }
  bool hasIncludeColumn() const { return columnOf(Column::Include) >= 0; }

 signals:
  /// A record was written from a cell (its form shows it again).
  void recordWritten(const QString &name);
  /// A cell's value was refused: the parser's words.
  void refused(const QString &message);

 private:
  bool writeField(int row, const std::string &key, const std::string &value, QString &why);
  const azure2::front::KindInfo *kindOf(int row) const;
  std::string quantityOf(int row) const;

  ProjectHolder *holder_;
  QString section_;
  QStringList kinds_;
  std::vector<Column> columns_;
  std::vector<azure2::front::Record> records_;
  std::vector<std::vector<QString>> cells_;  ///< short texts, per row and column
  bool writing_ = false;
};

/// The editors of a table's cells in place: a drop-down of the cell's choices, else a line.
class RecordCellDelegate : public QStyledItemDelegate {
  Q_OBJECT
 public:
  using QStyledItemDelegate::QStyledItemDelegate;
  QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const override;
  void setEditorData(QWidget *editor, const QModelIndex &index) const override;
  void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const override;
};

#endif
