#include "RecordsPanel.h"

#include <QAction>
#include <QCheckBox>
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QPainter>
#include <QLineEdit>
#include <QListWidget>
#include <QMenu>
#include <QMessageBox>
#include <QPushButton>
#include <QScrollArea>
#include <QSortFilterProxyModel>
#include <QSplitter>
#include <QTableView>
#include <QTimer>
#include <QToolButton>
#include <QVBoxLayout>

#include "RecordTableModel.h"
#include "ValueEditor.h"

namespace front = azure2::front;

namespace {

QString Q(const std::string &s) { return QString::fromStdString(s); }

QString ToolTip(const front::FieldInfo &f) {
  QString t = "<b>" + Q(f.title).toHtmlEscaped() + "</b> (<tt>" + Q(f.key).toHtmlEscaped() + "</tt>) &mdash; " +
              Q(f.doc).toHtmlEscaped() + "<br><tt>" + Q(f.grammar).toHtmlEscaped() + "</tt>";
  if (f.required) t += "<br>Required.";
  if (!f.requiredWhen.empty()) t += "<br>Required when " + Q(f.requiredWhen).toHtmlEscaped() + ".";
  if (!f.appliesWhen.empty()) t += "<br>Applies when " + Q(f.appliesWhen).toHtmlEscaped() + ".";
  if (f.byQuantity) t += "<br>The record's quantity decides whether it applies.";
  if (!f.defaultText.empty()) t += "<br>Left out: <tt>" + Q(f.defaultText).toHtmlEscaped() + "</tt>.";
  return t;
}

/// A box of the form that opens and closes: a field and the fields that belong with it, its doc on one line.
class FoldBox : public QWidget {
 public:
  FoldBox(const QString &title, const QString &doc, bool open) {
    QVBoxLayout *v = new QVBoxLayout(this);
    v->setContentsMargins(0, 2, 0, 2);
    v->setSpacing(1);
    header_ = new QToolButton;
    header_->setText(title);
    header_->setCheckable(true);
    header_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    header_->setAutoRaise(true);
    header_->setStyleSheet("QToolButton { font-weight: bold; }");
    content_ = new QWidget;
    QVBoxLayout *c = new QVBoxLayout(content_);
    c->setContentsMargins(18, 0, 0, 4);
    c->setSpacing(2);
    if (!doc.isEmpty()) {
      QLabel *d = new QLabel(doc);
      d->setWordWrap(true);
      d->setStyleSheet("color: palette(mid);");
      c->addWidget(d);
    }
    form_ = new QFormLayout;
    form_->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
    c->addLayout(form_);
    v->addWidget(header_);
    v->addWidget(content_);
    QObject::connect(header_, &QToolButton::toggled, this, [this](bool on) { setOpen(on); });
    header_->setChecked(open);
    setOpen(open);
  }
  void setOpen(bool on) {
    header_->setArrowType(on ? Qt::DownArrow : Qt::RightArrow);
    content_->setVisible(on);
    if (header_->isChecked() != on) header_->setChecked(on);
  }
  bool isOpen() const { return content_->isVisible(); }
  QFormLayout *form() const { return form_; }

 private:
  QToolButton *header_;
  QWidget *content_;
  QFormLayout *form_;
};

/// A small curve of a record's kernel (front::PreviewShape), redrawn as the form is edited.
class ShapeView : public QWidget {
 public:
  ShapeView() {
    setMinimumHeight(130);
    setMaximumHeight(160);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
  }
  void setShape(const front::ShapeCurve &shape) {
    shape_ = shape;
    update();
  }

 protected:
  void paintEvent(QPaintEvent *) override {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);
    const QRect area = rect().adjusted(50, 20, -12, -26);
    p.fillRect(rect(), palette().color(QPalette::Base));
    p.setPen(palette().color(QPalette::Mid));
    p.drawRect(area);
    p.setPen(palette().color(QPalette::Text));
    p.drawText(QRect(area.left(), 2, area.width(), 18), Qt::AlignCenter, QString::fromStdString(shape_.title));
    if (shape_.x.size() < 2) return;
    double lo = shape_.x.front(), hi = shape_.x.back(), top = 0;
    for (double y : shape_.y) top = std::max(top, y);
    if (!(hi > lo) || !(top > 0)) return;
    p.drawText(QRect(area.left(), area.bottom() + 2, area.width(), 20), Qt::AlignCenter, QString::fromStdString(shape_.xLabel));
    p.drawText(QRect(0, area.bottom() + 2, area.left() + 40, 20), Qt::AlignCenter, QString::number(lo, 'g', 4));
    p.drawText(QRect(area.right() - 60, area.bottom() + 2, 70, 20), Qt::AlignRight, QString::number(hi, 'g', 4));
    QPolygonF line;
    for (std::size_t i = 0; i < shape_.x.size(); i++)
      line << QPointF(area.left() + (shape_.x[i] - lo) / (hi - lo) * area.width(), area.bottom() - shape_.y[i] / top * area.height());
    p.setPen(QPen(palette().color(QPalette::Highlight), 2));
    p.drawPolyline(line);
  }

 private:
  front::ShapeCurve shape_;
};

bool IsSimple(const front::FieldInfo &f) {
  const std::string &t = f.value.type;
  return t != "structured" && !f.repeatable;
}

}  // namespace

RecordsPanel::RecordsPanel(ProjectHolder *holder, const QString &section, const QString &help, QWidget *parent,
                           const QStringList &kinds)
    : QWidget(parent), holder_(holder), section_(section), kinds_(kinds.isEmpty() ? holder->kindsOfSection(section) : kinds) {
  const QString id = section + (kinds.isEmpty() ? QString() : "_" + kinds.join("_"));
  setObjectName("records_" + id);

  model_ = new RecordTableModel(holder_, section_, kinds_, this);
  proxy_ = new QSortFilterProxyModel(this);
  proxy_->setSourceModel(model_);
  proxy_->setFilterCaseSensitivity(Qt::CaseInsensitive);
  proxy_->setFilterKeyColumn(-1);
  proxy_->setSortCaseSensitivity(Qt::CaseInsensitive);
  table_ = new QTableView;
  table_->setObjectName("recordList_" + id);
  table_->setModel(proxy_);
  table_->setItemDelegate(new RecordCellDelegate(table_));
  table_->setSortingEnabled(true);
  table_->sortByColumn(-1, Qt::AscendingOrder);
  table_->setSelectionBehavior(QAbstractItemView::SelectRows);
  table_->setSelectionMode(QAbstractItemView::ExtendedSelection);
  table_->setEditTriggers(QAbstractItemView::DoubleClicked | QAbstractItemView::EditKeyPressed | QAbstractItemView::SelectedClicked);
  table_->verticalHeader()->setVisible(false);
  table_->verticalHeader()->setDefaultSectionSize(fontMetrics().height() + 8);
  table_->horizontalHeader()->setSectionResizeMode(QHeaderView::Interactive);
  table_->horizontalHeader()->setStretchLastSection(true);
  table_->horizontalHeader()->setHighlightSections(false);
  table_->setAlternatingRowColors(true);
  table_->setWordWrap(false);
  table_->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(table_, &QTableView::customContextMenuRequested, this, &RecordsPanel::contextMenu);
  connect(table_->selectionModel(), &QItemSelectionModel::selectionChanged, this, &RecordsPanel::selectionChanged);
  connect(model_, &RecordTableModel::refused, this, [this](const QString &m) {
    status_->setText(m);
    status_->setVisible(true);
  });
  connect(model_, &RecordTableModel::recordWritten, this, [this](const QString &name) {
    emit edited();
    if (haveShown_ && (name == Q(shown_.name) || name == shownName()))
      QTimer::singleShot(0, this, [this, name]() {
        for (const front::Record &r : holder_->records(section_))
          if (Q(r.name) == name) buildForm(r);
      });
  });

  filter_ = new QLineEdit;
  filter_->setPlaceholderText(tr("Filter: any text of a row"));
  filter_->setClearButtonEnabled(true);
  connect(filter_, &QLineEdit::textChanged, proxy_, &QSortFilterProxyModel::setFilterFixedString);

  auto button = [this](const QString &text, const QString &tip) {
    QPushButton *b = new QPushButton(text);
    b->setToolTip(tip);
    return b;
  };
  addButton_ = button(tr("Add..."), tr("Add a record (Ins)"));
  duplicateButton_ = button(tr("Duplicate"), tr("A copy of the selected record under a new name (Ctrl+D)"));
  removeButton_ = button(tr("Remove"), tr("Remove the selected records (Del)"));
  upButton_ = button(tr("Up"), tr("Move the record one place up (Ctrl+Up)"));
  downButton_ = button(tr("Down"), tr("Move the record one place down (Ctrl+Down)"));
  includeButton_ = button(tr("Use"), tr("Use the selected records in the run"));
  excludeButton_ = button(tr("Leave out"), tr("Leave the selected records out of the run (they stay in the file)"));
  connect(addButton_, &QPushButton::clicked, this, &RecordsPanel::addRecord);
  connect(duplicateButton_, &QPushButton::clicked, this, &RecordsPanel::duplicateRecord);
  connect(removeButton_, &QPushButton::clicked, this, &RecordsPanel::removeRecord);
  connect(upButton_, &QPushButton::clicked, this, [this]() {
    QString why;
    if (!moveShown(-1, why)) status_->setText(why);
  });
  connect(downButton_, &QPushButton::clicked, this, [this]() {
    QString why;
    if (!moveShown(+1, why)) status_->setText(why);
  });
  connect(includeButton_, &QPushButton::clicked, this, [this]() {
    QString why;
    if (!includeSelected(true, why)) QMessageBox::warning(this, tr("Not changed"), why);
  });
  connect(excludeButton_, &QPushButton::clicked, this, [this]() {
    QString why;
    if (!includeSelected(false, why)) QMessageBox::warning(this, tr("Not changed"), why);
  });
  auto action = [this](const QString &text, const QKeySequence &key, QPushButton *b) {
    QAction *a = new QAction(text, this);
    a->setShortcut(key);
    a->setShortcutContext(Qt::WidgetWithChildrenShortcut);
    connect(a, &QAction::triggered, b, &QPushButton::click);
    table_->addAction(a);
    return a;
  };
  addAction_ = action(tr("Add..."), QKeySequence(Qt::Key_Insert), addButton_);
  duplicateAction_ = action(tr("Duplicate"), QKeySequence(Qt::CTRL | Qt::Key_D), duplicateButton_);
  removeAction_ = action(tr("Remove"), QKeySequence(QKeySequence::Delete), removeButton_);
  upAction_ = action(tr("Move up"), QKeySequence(Qt::CTRL | Qt::Key_Up), upButton_);
  downAction_ = action(tr("Move down"), QKeySequence(Qt::CTRL | Qt::Key_Down), downButton_);

  QHBoxLayout *buttons = new QHBoxLayout;
  for (QPushButton *b : {addButton_, duplicateButton_, removeButton_, upButton_, downButton_}) buttons->addWidget(b);
  buttons->addSpacing(12);
  buttons->addWidget(includeButton_);
  buttons->addWidget(excludeButton_);
  includeButton_->setVisible(model_->hasIncludeColumn());
  excludeButton_->setVisible(model_->hasIncludeColumn());
  buttons->addStretch(1);

  QWidget *top = new QWidget;
  QVBoxLayout *topLayout = new QVBoxLayout(top);
  topLayout->setContentsMargins(0, 0, 0, 0);
  if (!help.isEmpty()) {
    QLabel *helpLabel = new QLabel(help);
    helpLabel->setWordWrap(true);
    topLayout->addWidget(helpLabel);
  }
  topLayout->addWidget(filter_);
  topLayout->addWidget(table_, 1);
  topLayout->addLayout(buttons);

  kindLabel_ = new QLabel;
  kindLabel_->setWordWrap(true);
  usedBy_ = new QLabel;
  usedBy_->setWordWrap(true);
  usedBy_->setObjectName("usedBy");
  assignButton_ = new QPushButton(tr("Assign to..."));
  assignButton_->setToolTip(tr("Choose the records (data sets, extrapolations) that use this one"));
  connect(assignButton_, &QPushButton::clicked, this, &RecordsPanel::assignDialog);
  QHBoxLayout *usedRow = new QHBoxLayout;
  usedRow->addWidget(usedBy_, 1);
  usedRow->addWidget(assignButton_);
  everyField_ = new QCheckBox(tr("Every field"));
  everyField_->setToolTip(tr("Show the fields set rarely too (each is shown anyway when the record sets it)."));
  connect(everyField_, &QCheckBox::toggled, this, &RecordsPanel::showEveryField);
  applyButton_ = new QPushButton(tr("&Apply"));
  applyButton_->setObjectName("applyRecord_" + id);
  applyButton_->setToolTip(tr("Check the record with the engine and write it into the project (only its changed lines move)."));
  connect(applyButton_, &QPushButton::clicked, this, &RecordsPanel::apply);
  revertButton_ = new QPushButton(tr("Revert"));
  revertButton_->setToolTip(tr("Show the record as the project holds it"));
  connect(revertButton_, &QPushButton::clicked, this, &RecordsPanel::revert);
  status_ = new QLabel;
  status_->setWordWrap(true);
  status_->setStyleSheet("color: #c0392b;");
  QHBoxLayout *formButtons = new QHBoxLayout;
  formButtons->addWidget(everyField_);
  formButtons->addWidget(status_, 1);
  formButtons->addWidget(revertButton_);
  formButtons->addWidget(applyButton_);

  formHost_ = new QWidget;
  formLayout_ = new QVBoxLayout(formHost_);
  formLayout_->setContentsMargins(4, 4, 4, 4);
  QScrollArea *scroll = new QScrollArea;
  scroll->setWidgetResizable(true);
  scroll->setWidget(formHost_);
  QWidget *bottom = new QWidget;
  QVBoxLayout *bottomLayout = new QVBoxLayout(bottom);
  bottomLayout->setContentsMargins(0, 0, 0, 0);
  bottomLayout->addWidget(kindLabel_);
  bottomLayout->addLayout(usedRow);
  bottomLayout->addWidget(scroll, 1);
  bottomLayout->addLayout(formButtons);

  QSplitter *splitter = new QSplitter(Qt::Vertical);
  splitter->addWidget(top);
  splitter->addWidget(bottom);
  splitter->setStretchFactor(0, 2);
  splitter->setStretchFactor(1, 5);
  splitter->setSizes({220, 550});
  QVBoxLayout *layout = new QVBoxLayout(this);
  layout->setContentsMargins(0, 0, 0, 0);
  layout->addWidget(splitter);

  connect(holder_, &ProjectHolder::opened, this, &RecordsPanel::reload);
  connect(holder_, &ProjectHolder::sectionEdited, this, [this](const QString &s) {
    if (reloading_ || model_->writing()) return;
    if (s == section_) reload();
    else if (haveShown_ && s != section_) showUsedBy();  // a data set now names this target, ...
  });
  reload();
}

int RecordsPanel::recordCount() const { return model_->rowCount(); }

QString RecordsPanel::shownName() const { return haveShown_ ? Q(shown_.name) : QString(); }

void RecordsPanel::reload() {
  const QString keep = shownName();
  const QStringList selected = selectedNames();
  reloading_ = true;
  model_->reload();
  table_->resizeColumnsToContents();
  for (int c = 0; c < model_->columnCount(); c++)
    table_->setColumnWidth(c, std::min(table_->columnWidth(c), fontMetrics().horizontalAdvance(QString(40, 'x'))));
  reloading_ = false;
  haveShown_ = false;
  QStringList again;
  for (const QString &n : selected)
    if (model_->rowOf(n) >= 0) again << n;
  if (model_->rowOf(keep) >= 0) {
    again.removeAll(keep);
    again.prepend(keep);
  }
  if (again.isEmpty() && model_->rowCount() > 0) again << Q(model_->records()[0].name);
  if (again.isEmpty()) {
    clearForm();
    updateButtons();
  } else {
    selectRecords(again);
  }
  addButton_->setEnabled(!kinds_.isEmpty());
}

void RecordsPanel::selectRecords(const QStringList &names) {
  QItemSelection selection;
  QModelIndex first;
  for (const QString &n : names) {
    const int row = model_->rowOf(n);
    if (row < 0) continue;
    const QModelIndex at = proxy_->mapFromSource(model_->index(row, 0));
    if (!at.isValid()) continue;
    if (!first.isValid()) first = at;
    selection.select(at, at.siblingAtColumn(model_->columnCount() - 1));
  }
  if (!first.isValid()) return;
  table_->selectionModel()->setCurrentIndex(first, QItemSelectionModel::NoUpdate);
  table_->selectionModel()->select(selection, QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
  table_->scrollTo(first);
  selectionChanged();
}

bool RecordsPanel::selectRecord(const QString &name) {
  if (model_->rowOf(name) < 0) {
    // a filtered row: show it
    if (!filter_->text().isEmpty()) filter_->clear();
    if (model_->rowOf(name) < 0) return false;
  }
  if (!proxy_->mapFromSource(model_->index(model_->rowOf(name), 0)).isValid()) filter_->clear();
  selectRecords({name});
  return shownName() == name;
}

QStringList RecordsPanel::selectedNames() const {
  QStringList out;
  for (const QModelIndex &i : table_->selectionModel()->selectedRows()) {
    const QModelIndex s = proxy_->mapToSource(i);
    if (s.isValid() && s.row() < (int)model_->records().size()) out << Q(model_->records()[s.row()].name);
  }
  return out;
}

void RecordsPanel::updateButtons() {
  const int n = (int)table_->selectionModel()->selectedRows().size();
  const bool one = haveShown_;
  removeButton_->setEnabled(n > 0);
  duplicateButton_->setEnabled(one);
  upButton_->setEnabled(one && model_->rowOf(shownName()) > 0);
  downButton_->setEnabled(one && model_->rowOf(shownName()) + 1 < model_->rowCount());
  applyButton_->setEnabled(one);
  revertButton_->setEnabled(one);
  includeButton_->setEnabled(n > 0);
  excludeButton_->setEnabled(n > 0);
}

void RecordsPanel::selectionChanged() {
  if (reloading_) return;
  QModelIndex current = table_->selectionModel()->currentIndex();
  const QModelIndexList rows = table_->selectionModel()->selectedRows();
  if (!current.isValid() || !table_->selectionModel()->isRowSelected(current.row(), QModelIndex()))
    current = rows.isEmpty() ? QModelIndex() : rows.first();
  const QModelIndex source = proxy_->mapToSource(current);
  QString name;
  if (source.isValid() && source.row() < (int)model_->records().size()) name = Q(model_->records()[source.row()].name);
  if (haveShown_ && name == shownName()) {
    updateButtons();
    return;
  }
  // the record shown before: its edits are applied when the engine takes them, else it stays shown
  if (haveShown_ && formEdited()) {
    QString why;
    const QString was = shownName();
    if (!applyForm(why, true)) {
      status_->setText(tr("%1 not left: %2").arg(was, why));
      QTimer::singleShot(0, this, [this, was]() { selectRecords({was}); });
      return;
    }
    // applyForm read the section again and selected the record it wrote; go on to the one chosen
    if (!name.isEmpty() && name != shownName()) selectRecords({name});
    return;
  }
  status_->clear();
  for (const front::Record &r : model_->records())
    if (Q(r.name) == name) {
      buildForm(r);
      updateButtons();
      return;
    }
  clearForm();
  updateButtons();
}

void RecordsPanel::clearForm() {
  editors_.clear();
  boxes_.clear();
  rareRows_.clear();
  nameEdit_ = nullptr;
  preview_ = nullptr;
  haveShown_ = false;
  while (QLayoutItem *item = formLayout_->takeAt(0)) {
    if (QWidget *w = item->widget()) delete w;
    delete item;
  }
  kindLabel_->clear();
  usedBy_->clear();
  usedBy_->setVisible(false);
  assignButton_->setVisible(false);
}

QString RecordsPanel::quantityText() const {
  if (quantityKey_.isEmpty()) return QString();
  const QList<FieldEditor *> q = editors_.value(quantityKey_);
  if (!q.isEmpty()) return q.first()->text();
  return Q(shown_.Value(quantityKey_.toStdString()));
}

void RecordsPanel::quantityChanged() {
  for (const QList<FieldEditor *> &list : editors_)
    for (FieldEditor *e : list)
      if (e->field().value.unitOfQuantity || e->field().unitOfQuantity) e->refresh();
}

FieldEditor *RecordsPanel::makeEditor(const front::FieldInfo &field, const QString &text) {
  EditorContext context;
  context.holder = holder_;
  context.kind = shown_.kind;
  context.quantity = [this]() { return quantityText(); };
  FieldEditor *e = new FieldEditor(field, context);
  e->setToolTip(ToolTip(field));
  e->setText(text);
  editors_[Q(field.key)] << e;
  connect(e, &FieldEditor::changed, this, [this]() {
    status_->clear();
    updatePreview();
  });
  if (Q(field.key) == quantityKey_) connect(e, &FieldEditor::changed, this, &RecordsPanel::quantityChanged);
  return e;
}

void RecordsPanel::addFieldRows(const front::KindInfo &kind, const front::Record &record, QFormLayout *form,
                                const front::FieldInfo &field, bool inBox) {
  Q_UNUSED(kind);
  Q_UNUSED(inBox);
  const QString key = Q(field.key);
  QStringList values;
  for (const front::Field &x : record.fields)
    if (x.key == field.key) values << Q(x.value);
  QLabel *label = new QLabel(Q(field.title) + (field.required ? " *" : ""));
  label->setToolTip(ToolTip(field));
  QWidget *row = nullptr;
  if (field.repeatable) {
    QWidget *host = new QWidget;
    QVBoxLayout *v = new QVBoxLayout(host);
    v->setContentsMargins(0, 0, 0, 0);
    v->setSpacing(2);
    auto addEntry = [this, field, v, key](const QString &text, bool fresh) {
      QWidget *entry = new QWidget;
      QHBoxLayout *h = new QHBoxLayout(entry);
      h->setContentsMargins(0, 0, 0, 0);
      FieldEditor *e = makeEditor(field, text);
      if (fresh) e->setUsed(true);
      QToolButton *remove = new QToolButton;
      remove->setText(QString::fromUtf8("−"));
      remove->setToolTip(tr("Remove this %1").arg(Q(field.title)));
      h->addWidget(e, 1);
      h->addWidget(remove, 0, Qt::AlignTop);
      v->insertWidget(v->count() - 1, entry);
      connect(remove, &QToolButton::clicked, this, [this, entry, e, key]() {
        editors_[key].removeAll(e);
        entry->deleteLater();
        status_->clear();
      });
    };
    QPushButton *more = new QPushButton(tr("Add %1").arg(Q(field.title).toLower()));
    more->setObjectName("add_" + key);
    QHBoxLayout *moreRow = new QHBoxLayout;
    moreRow->addWidget(more);
    moreRow->addStretch(1);
    v->addLayout(moreRow);
    editors_[key];  // the key is the form's even with no entry
    for (const QString &value : values) addEntry(value, false);
    connect(more, &QPushButton::clicked, this, [addEntry]() { addEntry(QString(), true); });
    row = host;
  } else {
    row = makeEditor(field, values.isEmpty() ? QString() : values.first());
  }
  form->addRow(label, row);
  if (field.rare && values.isEmpty()) {
    rareRows_ << label << row;
    label->setVisible(everyField_->isChecked());
    row->setVisible(everyField_->isChecked());
  }
}

void RecordsPanel::buildForm(const front::Record &record) {
  clearForm();
  shown_ = record;
  haveShown_ = true;
  const front::KindInfo *kind = holder_->kind(record.kind);
  kindLabel_->setText(kind ? "<b>" + Q(kind->title).toHtmlEscaped() + "</b> &nbsp;" + Q(kind->doc).toHtmlEscaped()
                           : tr("<b>%1</b>: a kind the schema does not know").arg(Q(record.kind)));
  quantityKey_.clear();
  if (kind)
    for (const front::FieldInfo &f : kind->fields)
      if (IsQuantityField(f, kind->kind)) quantityKey_ = Q(f.key);

  QFormLayout *basic = new QFormLayout;
  basic->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
  QWidget *basicHost = new QWidget;
  basicHost->setLayout(basic);
  formLayout_->addWidget(basicHost);
  nameEdit_ = new QLineEdit(Q(record.name));
  nameEdit_->setObjectName("recordName");
  nameEdit_->setToolTip(tr("The record's name: other records refer to it by this name."));
  basic->addRow(tr("Name"), nameEdit_);
  if (kind) {
    auto has = [&record](const std::string &key) {
      for (const front::Field &x : record.fields)
        if (x.key == key) return true;
      return false;
    };
    for (const front::FieldInfo &f : kind->fields) {
      bool grouped = false;
      if (!f.group.empty())
        for (const front::FieldInfo &g : kind->fields) grouped = grouped || g.key == f.group;
      if (grouped) continue;
      std::vector<const front::FieldInfo *> members;
      for (const front::FieldInfo &m : kind->fields)
        if (m.group == f.key) members.push_back(&m);
      if (IsSimple(f) && members.empty()) {
        addFieldRows(*kind, record, basic, f, false);
        continue;
      }
      bool used = has(f.key);
      for (const front::FieldInfo *m : members) used = used || has(m->key);
      QString title = Q(f.title);
      if (!f.value.units.empty() && f.value.type != "structured") title += "  [" + Q(f.value.units.front()) + "]";
      FoldBox *box = new FoldBox(title, Q(f.doc), used || f.required);
      box->setObjectName("box_" + Q(f.key));
      box->setToolTip(ToolTip(f));
      addFieldRows(*kind, record, box->form(), f, true);
      for (const front::FieldInfo *m : members) addFieldRows(*kind, record, box->form(), *m, true);
      formLayout_->addWidget(box);
      boxes_.insert(Q(f.key), box);
      if (f.rare && !used) {
        rareRows_ << box;
        box->setVisible(everyField_->isChecked());
      }
    }
  }
  // fields the schema does not list for the kind: shown as text, kept on Apply (the engine names them)
  for (const front::Field &x : record.fields) {
    bool known = false;
    if (kind)
      for (const front::FieldInfo &f : kind->fields) known = known || f.key == x.key;
    if (known || editors_.contains(Q(x.key))) continue;
    front::FieldInfo unknown;
    unknown.key = x.key;
    unknown.title = x.key + " (not a field of " + record.kind + ")";
    unknown.type = unknown.value.type = "text";
    addFieldRows(kind ? *kind : front::KindInfo(), record, basic, unknown, false);
  }
  formLayout_->insertWidget(1, preview_ = new ShapeView);
  preview_->setObjectName("shapePreview");
  updatePreview();
  formLayout_->addStretch(1);
  showUsedBy();
}

void RecordsPanel::updatePreview() {
  if (!preview_) return;
  front::ShapeCurve shape;
  std::string why;
  const bool drawn = front::PreviewShape(formRecord(), 200, shape, why);
  static_cast<ShapeView *>(preview_)->setShape(shape);
  preview_->setVisible(drawn);
}

void RecordsPanel::showUsedBy() {
  if (!haveShown_) return;
  bool referable = false;
  for (const front::KindInfo &k : holder_->kinds())
    for (const front::FieldInfo &f : k.fields)
      if (f.value.type == "reference" &&
          std::find(f.value.values.begin(), f.value.values.end(), shown_.kind) != f.value.values.end())
        referable = true;
  usedBy_->setVisible(referable);
  assignButton_->setVisible(referable);
  if (!referable) return;
  const QStringList users = usedBy();
  QStringList names;
  for (const QString &u : users) names << u.section('/', 1);
  usedBy_->setText(users.isEmpty() ? tr("Used by no record.")
                                   : tr("Used by %n record(s): %1", nullptr, users.size()).arg(names.join(", ")));
}

QStringList RecordsPanel::usedBy() const {
  QStringList out;
  if (!haveShown_) return out;
  for (const auto &section : front::RecordSections())
    for (const front::Record &r : holder_->records(Q(section.first))) {
      const front::KindInfo *k = holder_->kind(r.kind);
      if (!k) continue;
      for (const front::FieldInfo &f : k->fields) {
        if (f.value.type != "reference" ||
            std::find(f.value.values.begin(), f.value.values.end(), shown_.kind) == f.value.values.end())
          continue;
        const QString value = Q(r.Value(f.key));
        if (!value.isEmpty() && value.section('@', 0, 0) == Q(shown_.name)) out << Q(section.first) + "/" + Q(r.name);
      }
    }
  return out;
}

bool RecordsPanel::assignShownTo(const QStringList &names, QString &why) {
  if (!haveShown_) {
    why = tr("no record is shown");
    return false;
  }
  for (const auto &section : front::RecordSections())
    for (front::Record r : holder_->records(Q(section.first))) {
      const front::KindInfo *k = holder_->kind(r.kind);
      if (!k) continue;
      for (const front::FieldInfo &f : k->fields) {
        if (f.value.type != "reference" ||
            std::find(f.value.values.begin(), f.value.values.end(), shown_.kind) == f.value.values.end())
          continue;
        const QString value = Q(r.Value(f.key));
        const bool names_ = value.section('@', 0, 0) == Q(shown_.name);
        const bool wanted = names.contains(Q(r.name));
        if (wanted == names_) continue;
        std::vector<front::Field> fields;
        bool done = false;
        for (const front::Field &x : r.fields) {
          if (x.key == f.key) {
            if (wanted && !done) fields.push_back({x.key, shown_.name, x.comment});
            done = true;
            continue;
          }
          fields.push_back(x);
        }
        if (wanted && !done) fields.push_back({f.key, shown_.name, ""});
        r.fields = fields;
        if (!holder_->setRecord(Q(section.first), r, why)) return false;
      }
    }
  showUsedBy();
  emit edited();
  return true;
}

void RecordsPanel::assignDialog() {
  QDialog dialog(this);
  dialog.setWindowTitle(tr("Records that use %1").arg(Q(shown_.name)));
  QVBoxLayout *v = new QVBoxLayout(&dialog);
  v->addWidget(new QLabel(tr("Check the records that use %1 (a record names one %2 at most):").arg(Q(shown_.name), Q(shown_.kind))));
  QListWidget *list = new QListWidget;
  const QStringList users = usedBy();
  for (const auto &section : front::RecordSections())
    for (const front::Record &r : holder_->records(Q(section.first))) {
      const front::KindInfo *k = holder_->kind(r.kind);
      if (!k) continue;
      for (const front::FieldInfo &f : k->fields)
        if (f.value.type == "reference" &&
            std::find(f.value.values.begin(), f.value.values.end(), shown_.kind) != f.value.values.end()) {
          QListWidgetItem *item = new QListWidgetItem(Q(r.name) + "   (" + Q(r.Value("reaction")) + " " + Q(r.Value(f.key)) + ")");
          item->setData(Qt::UserRole, Q(r.name));
          item->setCheckState(users.contains(Q(section.first) + "/" + Q(r.name)) ? Qt::Checked : Qt::Unchecked);
          list->addItem(item);
          break;
        }
    }
  v->addWidget(list);
  QDialogButtonBox *box = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
  connect(box, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
  connect(box, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
  v->addWidget(box);
  if (dialog.exec() != QDialog::Accepted) return;
  QStringList names;
  for (int i = 0; i < list->count(); i++)
    if (list->item(i)->checkState() == Qt::Checked) names << list->item(i)->data(Qt::UserRole).toString();
  QString why;
  if (!assignShownTo(names, why)) QMessageBox::warning(this, tr("Not assigned"), why);
}

void RecordsPanel::showEveryField(bool on) {
  for (QWidget *w : rareRows_) w->setVisible(on);
}

QWidget *RecordsPanel::fieldWidget(const QString &key) const { return fieldEditor(key); }

FieldEditor *RecordsPanel::fieldEditor(const QString &key, int occurrence) const {
  const QList<FieldEditor *> list = editors_.value(key);
  return occurrence < list.size() ? list[occurrence] : nullptr;
}

bool RecordsPanel::setFieldText(const QString &key, const QString &text) {
  FieldEditor *e = fieldEditor(key);
  if (!e) return false;
  e->typeText(text);
  return true;
}

front::Record RecordsPanel::formRecord() const {
  front::Record r;
  r.kind = shown_.kind;
  r.name = shown_.name;
  QMap<QString, QStringList> texts;
  for (auto it = editors_.begin(); it != editors_.end(); ++it)
    for (FieldEditor *e : it.value()) {
      const QString t = e->text().trimmed();
      if (!t.isEmpty()) texts[it.key()] << t;
    }
  QStringList done;
  for (const front::Field &x : shown_.fields) {
    const QString key = Q(x.key);
    if (done.contains(key)) continue;
    done << key;
    if (!editors_.contains(key)) {
      for (const front::Field &y : shown_.fields)
        if (y.key == x.key) r.fields.push_back(y);
      continue;
    }
    std::vector<std::string> comments;
    for (const front::Field &y : shown_.fields)
      if (y.key == x.key) comments.push_back(y.comment);
    const QStringList values = texts.value(key);
    for (int i = 0; i < values.size(); i++)
      r.fields.push_back({x.key, values[i].toStdString(), i < (int)comments.size() ? comments[i] : std::string()});
  }
  if (const front::KindInfo *kind = holder_->kind(shown_.kind))
    for (const front::FieldInfo &f : kind->fields) {
      const QString key = Q(f.key);
      if (done.contains(key)) continue;
      for (const QString &v : texts.value(key)) r.fields.push_back({f.key, v.toStdString(), ""});
    }
  return r;
}

bool RecordsPanel::formEdited() const {
  if (!haveShown_) return false;
  if (nameEdit_ && nameEdit_->text().trimmed() != Q(shown_.name)) return true;
  for (const QList<FieldEditor *> &list : editors_)
    for (FieldEditor *e : list)
      if (e->isEdited()) return true;
  // an entry of a repeatable field added or removed
  const front::Record r = formRecord();
  if (r.fields.size() != shown_.fields.size()) return true;
  for (std::size_t i = 0; i < r.fields.size(); i++)
    if (r.fields[i].key != shown_.fields[i].key || r.fields[i].value != shown_.fields[i].value) return true;
  return false;
}

QMap<QString, QStringList> RecordsPanel::errorsWith(const front::Record &record) const {
  QMap<QString, QStringList> out;
  front::Project copy = front::Project::FromText(holder_->project().Text(), holder_->project().Path());
  std::string why;
  if (!copy.SetRecord(section_.toStdString(), record, why)) {
    out[QString()] << Q(why);
    return out;
  }
  for (const front::Diagnostic &d : copy.Validate().diagnostics)
    if (d.severity == "error" && d.record == record.name) out[Q(d.field)] << Q(d.message);
  return out;
}

bool RecordsPanel::applyForm(QString &why, bool refuseErrors) {
  if (!haveShown_ || !nameEdit_) {
    why = tr("no record is shown");
    return false;
  }
  const QString newName = nameEdit_->text().trimmed();
  if (!formEdited()) return true;  // nothing changed: not a byte moves
  const front::Record r = formRecord();
  const QMap<QString, QStringList> errors = errorsWith(r);
  QStringList said;
  for (auto it = errors.begin(); it != errors.end(); ++it) {
    if (FieldEditor *e = fieldEditor(it.key())) e->showMessage(it.value().join(" "));
    said << (it.key().isEmpty() ? QString() : it.key() + ": ") + it.value().join(" ");
  }
  if (refuseErrors && !errors.isEmpty()) {
    why = said.join("; ");
    status_->setText(why);
    return false;
  }
  reloading_ = true;
  bool ok = holder_->setRecord(section_, r, why);
  if (ok && newName != Q(shown_.name)) ok = holder_->renameRecord(section_, Q(shown_.name), newName, why);
  reloading_ = false;
  haveShown_ = false;
  reload();
  if (ok) {
    selectRecord(newName);
    emit edited();
  }
  if (!said.isEmpty()) {
    status_->setText(tr("Written; the engine refuses: %1").arg(said.join("; ")));
    for (auto it = errors.begin(); it != errors.end(); ++it)
      if (FieldEditor *e = fieldEditor(it.key())) e->showMessage(it.value().join(" "));
  }
  return ok;
}

void RecordsPanel::apply() {
  if (!haveShown_) return;
  const QMap<QString, QStringList> errors = errorsWith(formRecord());
  if (!errors.isEmpty()) {
    QStringList lines;
    for (auto it = errors.begin(); it != errors.end(); ++it) {
      QString title = it.key();
      if (FieldEditor *e = fieldEditor(it.key())) {
        e->showMessage(it.value().join(" "));
        title = Q(e->field().title);
      }
      lines << (title.isEmpty() ? QString() : "<b>" + title.toHtmlEscaped() + "</b>: ") + it.value().join(" ").toHtmlEscaped();
    }
    if (QMessageBox::question(this, tr("The engine refuses this record"),
                              tr("<p>%1</p><p>Write it anyway? The diagnostics will name it until it is fixed.</p>")
                                  .arg(lines.join("<br>")),
                              QMessageBox::Yes | QMessageBox::No, QMessageBox::No) != QMessageBox::Yes)
      return;
  }
  QString why;
  if (!applyForm(why)) QMessageBox::warning(this, tr("Record not written"), why);
}

void RecordsPanel::revert() {
  if (!haveShown_) return;
  const QString name = shownName();
  for (const front::Record &r : holder_->records(section_))
    if (Q(r.name) == name) buildForm(r);
  status_->clear();
}

QString RecordsPanel::uniqueName(const QString &stem) const {
  QStringList names;
  for (const front::Record &r : holder_->records(section_)) names << Q(r.name);
  QString name = stem;
  for (int i = 2; names.contains(name); i++) name = QString("%1_%2").arg(stem).arg(i);
  return name;
}

void RecordsPanel::addRecord() {
  if (kinds_.isEmpty()) return;
  QDialog dialog(this);
  dialog.setWindowTitle(tr("Add a record"));
  dialog.setObjectName("addRecordDialog");
  QVBoxLayout *v = new QVBoxLayout(&dialog);
  QFormLayout *top = new QFormLayout;
  QComboBox *kind = new QComboBox;
  for (const QString &k : kinds_) {
    const front::KindInfo *info = holder_->kind(k.toStdString());
    kind->addItem(info ? Q(info->title) : k, k);
  }
  QLineEdit *name = new QLineEdit;
  top->addRow(tr("Kind"), kind);
  top->addRow(tr("Name"), name);
  kind->setVisible(kinds_.size() > 1);
  top->labelForField(kind)->setVisible(kinds_.size() > 1);
  v->addLayout(top);
  QLabel *doc = new QLabel;
  doc->setWordWrap(true);
  v->addWidget(doc);
  QWidget *fieldsHost = new QWidget;
  QFormLayout *fields = new QFormLayout(fieldsHost);
  fields->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
  v->addWidget(fieldsHost, 1);
  QList<FieldEditor *> editors;
  QString quantityKey;
  auto quantity = [&editors, &quantityKey]() {
    for (FieldEditor *e : editors)
      if (Q(e->field().key) == quantityKey) return e->text();
    return QString();
  };
  auto build = [&]() {
    for (FieldEditor *e : editors) delete e;
    editors.clear();
    while (fields->rowCount()) fields->removeRow(0);
    const std::string k = kind->currentData().toString().toStdString();
    const front::KindInfo *info = holder_->kind(k);
    name->setText(uniqueName(kind->currentData().toString()));
    doc->setText(info ? Q(info->doc) : QString());
    if (!info) return;
    quantityKey.clear();
    for (const front::FieldInfo &f : info->fields)
      if (IsQuantityField(f, k)) quantityKey = Q(f.key);
    EditorContext context;
    context.holder = holder_;
    context.kind = k;
    context.quantity = quantity;
    for (const front::FieldInfo &f : info->fields) {
      if (!f.required || !f.appliesWhen.empty()) continue;  // required wherever it applies, and it always applies
      FieldEditor *e = new FieldEditor(f, context);
      e->setToolTip(ToolTip(f));
      e->setText(QString());
      editors << e;
      fields->addRow(Q(f.title) + " *", e);
      e->prefill();
      if (Q(f.key) == quantityKey)
        connect(e, &FieldEditor::changed, &dialog, [&editors]() {
          for (FieldEditor *x : editors)
            if (x->field().value.unitOfQuantity) {
              x->refresh();
              x->prefill();
            }
        });
    }
  };
  build();
  connect(kind, QOverload<int>::of(&QComboBox::activated), &dialog, build);
  QLabel *note = new QLabel(tr("Fields left empty are written as '?': the diagnostics name them until they are given."));
  note->setWordWrap(true);
  note->setStyleSheet("color: palette(mid);");
  v->addWidget(note);
  QDialogButtonBox *box = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
  connect(box, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
  connect(box, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
  v->addWidget(box);
  dialog.resize(700, 420);
  if (dialog.exec() != QDialog::Accepted || name->text().trimmed().isEmpty()) return;
  front::Record r;
  r.kind = kind->currentData().toString().toStdString();
  r.name = name->text().trimmed().toStdString();
  for (FieldEditor *e : editors) {
    const QString t = e->text().trimmed();
    r.fields.push_back({e->field().key, t.isEmpty() ? "?" : t.toStdString(), ""});
  }
  QString why;
  if (!holder_->setRecord(section_, r, why)) {
    QMessageBox::warning(this, tr("Record not added"), why);
    return;
  }
  selectRecord(Q(r.name));
  emit edited();
}

bool RecordsPanel::addNamedRecord(const QString &kind, const QString &name, QString &why) {
  front::Record r;
  r.kind = kind.toStdString();
  r.name = name.trimmed().toStdString();
  // the required fields as '?' placeholders: the engine names each one until it is given (E_PLACEHOLDER)
  if (const front::KindInfo *k = holder_->kind(r.kind))
    for (const front::FieldInfo &f : k->fields)
      if (f.required && f.appliesWhen.empty()) r.fields.push_back({f.key, "?", ""});
  if (!holder_->setRecord(section_, r, why)) return false;
  selectRecord(name.trimmed());
  emit edited();
  return true;
}

void RecordsPanel::duplicateRecord() {
  QString name, why;
  if (haveShown_ && !duplicateShown(name, why)) QMessageBox::warning(this, tr("Record not added"), why);
}

bool RecordsPanel::duplicateShown(QString &newName, QString &why) {
  if (!haveShown_) {
    why = tr("no record is shown");
    return false;
  }
  front::Record r = shown_;
  r.name = uniqueName(Q(shown_.name) + "_copy").toStdString();
  if (!holder_->setRecord(section_, r, why)) return false;
  newName = Q(r.name);
  selectRecord(newName);
  emit edited();
  return true;
}

void RecordsPanel::removeRecord() {
  const QStringList names = selectedNames();
  if (names.isEmpty()) return;
  if (names.size() > 1 &&
      QMessageBox::question(this, tr("Remove records"), tr("Remove %1 records: %2?").arg(names.size()).arg(names.join(", "))) !=
          QMessageBox::Yes)
    return;
  for (const QString &n : names) {
    QString why;
    if (!holder_->removeRecord(section_, n, why)) {
      QMessageBox::warning(this, tr("Record not removed"), why);
      return;
    }
  }
  emit edited();
}

bool RecordsPanel::removeShown(QString &why) {
  if (!haveShown_) {
    why = tr("no record is shown");
    return false;
  }
  const QString name = shownName();
  haveShown_ = false;
  if (!holder_->removeRecord(section_, name, why)) return false;
  emit edited();
  return true;
}

bool RecordsPanel::renameShown(const QString &newName, QString &why) {
  if (!nameEdit_) {
    why = tr("no record is shown");
    return false;
  }
  nameEdit_->setText(newName);
  return applyForm(why);
}

bool RecordsPanel::moveShown(int delta, QString &why) {
  if (!haveShown_) {
    why = tr("no record is shown");
    return false;
  }
  if (formEdited() && !applyForm(why, true)) return false;
  const QString name = shownName();
  // the neighbour in the file among this panel's kinds
  const std::vector<front::Record> all = holder_->records(section_);
  int at = -1;
  for (int i = 0; i < (int)all.size(); i++)
    if (Q(all[i].name) == name) at = i;
  int steps = 0;
  for (int i = at + delta; i >= 0 && i < (int)all.size(); i += delta) {
    steps++;
    if (kinds_.contains(Q(all[i].kind))) break;
  }
  if (at < 0 || steps == 0 || at + delta * steps < 0 || at + delta * steps >= (int)all.size()) {
    why = tr("%1 is already the %2 record").arg(name, delta < 0 ? tr("first") : tr("last"));
    return false;
  }
  std::string reason;
  reloading_ = true;
  bool ok = true;
  for (int i = 0; i < steps && ok; i++) ok = holder_->project().MoveRecord(section_.toStdString(), name.toStdString(), delta, reason);
  reloading_ = false;
  if (!ok) {
    why = Q(reason);
    return false;
  }
  emit holder_->sectionEdited(section_);
  selectRecord(name);
  emit edited();
  return true;
}

bool RecordsPanel::includeSelected(bool on, QString &why) {
  const QStringList names = selectedNames();
  const bool ok = model_->setIncluded(names, on, why);
  selectRecords(names);
  return ok;
}

void RecordsPanel::setFilter(const QString &text) { filter_->setText(text); }

void RecordsPanel::contextMenu(const QPoint &at) {
  QMenu menu(this);
  menu.addAction(addAction_);
  menu.addAction(duplicateAction_);
  menu.addAction(removeAction_);
  menu.addSeparator();
  menu.addAction(upAction_);
  menu.addAction(downAction_);
  if (model_->hasIncludeColumn()) {
    menu.addSeparator();
    menu.addAction(tr("Use"), includeButton_, &QPushButton::click);
    menu.addAction(tr("Leave out"), excludeButton_, &QPushButton::click);
  }
  duplicateAction_->setEnabled(duplicateButton_->isEnabled());
  removeAction_->setEnabled(removeButton_->isEnabled());
  upAction_->setEnabled(upButton_->isEnabled());
  downAction_->setEnabled(downButton_->isEnabled());
  menu.exec(table_->viewport()->mapToGlobal(at));
}

QStringList RecordsPanel::recordNames() const {
  QStringList out;
  for (const front::Record &r : model_->records()) out << Q(r.name);
  return out;
}
