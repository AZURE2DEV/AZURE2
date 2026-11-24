#include "RecordsPanel.h"
#include "ExpertMode.h"

#include <QAction>
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleValidator>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QEvent>
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
#include <QRegularExpression>
#include <QScrollBar>
#include <QSettings>
#include <QSignalBlocker>
#include <QToolButton>
#include <QToolTip>
#include <algorithm>
#include <QTableView>
#include <QTimer>
#include <QToolButton>
#include <QVBoxLayout>
#include <cmath>
#include <functional>

#include "RecordDialog.h"
#include "RecordTableModel.h"
#include "ValueEditor.h"

namespace front = azure2::front;

namespace {

QString Q(const std::string &s) { return QString::fromStdString(s); }

QString ToolTip(const front::FieldInfo &f) {
  // the file's key and the value's grammar in the expert view; in standard view the title and what it is
  // (qt_gui.md section 19.9)
  const bool expert = expert::Enabled();
  QString t = "<b>" + Q(f.title).toHtmlEscaped() + "</b>" +
              (expert ? " (<tt>" + Q(f.key).toHtmlEscaped() + "</tt>)" : QString()) + " &mdash; " + Q(f.doc).toHtmlEscaped() +
              (expert ? "<br><tt>" + Q(f.grammar).toHtmlEscaped() + "</tt>" : QString());
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

/// A small curve of a record's kernel, its half maximum marked, redrawn as the form is edited.
class ShapeView : public QWidget {
 public:
  ShapeView() {
    setMinimumHeight(120);
    setMaximumHeight(150);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
  }
  void setShape(const front::KernelCurve &shape) {
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
    p.drawText(QRect(0, area.bottom() + 2, area.left() + 40, 20), Qt::AlignCenter, QString::number(lo, 'g', 6));
    p.drawText(QRect(area.right() - 80, area.bottom() + 2, 90, 20), Qt::AlignRight, QString::number(hi, 'g', 6));
    auto X = [&](double x) { return area.left() + (x - lo) / (hi - lo) * area.width(); };
    QPolygonF line;
    for (std::size_t i = 0; i < shape_.x.size(); i++) line << QPointF(X(shape_.x[i]), area.bottom() - shape_.y[i] / top * area.height());
    // the half maximum, between the curve's two crossings (the FWHM written under the curve)
    if (shape_.fwhm_keV == shape_.fwhm_keV && shape_.fwhm_keV > 0) {
      std::size_t peak = 0;
      for (std::size_t i = 1; i < shape_.y.size(); i++)
        if (shape_.y[i] > shape_.y[peak]) peak = i;
      double left = shape_.x[peak], right = shape_.x[peak];
      for (std::size_t i = peak; i > 0 && shape_.y[i - 1] >= 0.5 * top; i--) left = shape_.x[i - 1];
      for (std::size_t i = peak; i + 1 < shape_.y.size() && shape_.y[i + 1] >= 0.5 * top; i++) right = shape_.x[i + 1];
      p.setPen(QPen(palette().color(QPalette::Mid), 1, Qt::DashLine));
      const double y = area.bottom() - 0.5 * area.height();
      p.drawLine(QPointF(X(left), y), QPointF(X(right), y));
    }
    p.setPen(QPen(palette().color(QPalette::Highlight), 2));
    p.drawPolyline(line);
  }

 private:
  front::KernelCurve shape_;
};

/// The kernels a record of <targets> selects (front::KernelPreviewOf): for which data set and at which energy, then
/// one curve per kernel with its FWHM in keV under it.
class KernelView : public QWidget {
 public:
  explicit KernelView(std::function<void(const QString &, double)> pointChosen) : pointChosen_(std::move(pointChosen)) {
    setObjectName("kernelPreview");
    QVBoxLayout *v = new QVBoxLayout(this);
    v->setContentsMargins(0, 4, 0, 4);
    v->setSpacing(2);
    QHBoxLayout *top = new QHBoxLayout;
    QLabel *title = new QLabel(tr("<b>What it does to a point</b>"));
    users_ = new QComboBox;
    users_->setObjectName("kernelUser");
    users_->setSizeAdjustPolicy(QComboBox::AdjustToContents);
    users_->setToolTip(tr("The data set or extrapolation whose entrance pair (masses) the kernels are drawn for"));
    energy_ = new QLineEdit;
    energy_->setObjectName("kernelEnergy");
    energy_->setValidator(new QDoubleValidator(0.0, 1e9, 12, energy_));
    energy_->setMaximumWidth(energy_->fontMetrics().horizontalAdvance("0.000000000"));
    energy_->setToolTip(tr("The beam's lab energy E_0 at which the kernels are drawn (default: the middle of the data)"));
    pair_ = new QLabel;
    pair_->setStyleSheet("color: palette(mid);");
    top->addWidget(title);
    top->addSpacing(8);
    top->addWidget(new QLabel(tr("for")));
    top->addWidget(users_);
    top->addWidget(new QLabel(tr("at E_lab")));
    top->addWidget(energy_);
    top->addWidget(new QLabel(tr("MeV")));
    top->addWidget(pair_, 1);
    v->addLayout(top);
    curves_ = new QVBoxLayout;
    curves_->setSpacing(0);
    v->addLayout(curves_);
    notes_ = new QLabel;
    notes_->setWordWrap(true);
    notes_->setStyleSheet("color: palette(mid);");
    v->addWidget(notes_);
    connect(users_, QOverload<int>::of(&QComboBox::activated), this, [this]() { chosen(); });
    connect(energy_, &QLineEdit::editingFinished, this, [this]() { chosen(); });
  }
  void setPreview(const front::KernelPreview &k) {
    users_->blockSignals(true);
    users_->clear();
    for (const std::string &u : k.users) users_->addItem(QString::fromStdString(u), QString::fromStdString(u));
    if (k.users.empty()) users_->addItem(tr("(no data set names it)"), QString());
    const int at = users_->findData(QString::fromStdString(k.user));
    users_->setCurrentIndex(at >= 0 ? at : 0);
    users_->setEnabled(k.users.size() > 1);
    users_->blockSignals(false);
    if (!energy_->hasFocus()) energy_->setText(QString::number(k.labEnergy, 'g', 6));
    pair_->setText(tr("%1, E_c.m. = %2 MeV%3")
                       .arg(QString::fromStdString(k.pair))
                       .arg(k.cmEnergy, 0, 'g', 6)
                       .arg(k.highEnergy > 0 ? tr("; data %1 - %2 MeV").arg(k.lowEnergy, 0, 'g', 4).arg(k.highEnergy, 0, 'g', 4)
                                             : QString()));
    while (QLayoutItem *item = curves_->takeAt(0)) {
      delete item->widget();
      delete item;
    }
    for (const front::KernelCurve &c : k.curves) {
      ShapeView *view = new ShapeView;
      view->setObjectName("kernel_" + QString::fromStdString(c.kernel));
      view->setShape(c);
      QLabel *width = new QLabel;
      width->setObjectName("fwhm_" + QString::fromStdString(c.kernel));
      width->setAlignment(Qt::AlignCenter);
      width->setText(c.fwhm_keV == c.fwhm_keV ? tr("FWHM %1 keV (lab), %2 keV (c.m.)")
                                                    .arg(c.fwhm_keV, 0, 'g', 4)
                                                    .arg(c.fwhmCm_keV, 0, 'g', 4)
                                              : tr("no half maximum on both sides"));
      curves_->addWidget(view);
      curves_->addWidget(width);
    }
    QStringList notes;
    for (const std::string &n : k.notes) notes << QString::fromStdString(n);
    notes_->setText(notes.join("; "));
    notes_->setVisible(!notes.isEmpty());
  }

 private:
  void chosen() {
    bool ok = false;
    const double e = energy_->text().toDouble(&ok);
    pointChosen_(users_->currentData().toString(), ok && e > 0 ? e : 0.0);
  }
  std::function<void(const QString &, double)> pointChosen_;
  QComboBox *users_;
  QLineEdit *energy_;
  QLabel *pair_, *notes_;
  QVBoxLayout *curves_;
};

/// A title of the schema as the form shows it: J^pi with a superscript pi, as the Levels tab writes it.
QString DisplayTitle(const std::string &title) { return QString::fromStdString(title).replace("J^pi", QString::fromUtf8("Jπ")); }

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

  expert_ = expert::Enabled();
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
  // a column the user drags keeps its width (fitColumns no longer runs on a resize); the last column, which the header
  // stretches, does not count
  table_->viewport()->installEventFilter(this);
  connect(table_->horizontalHeader(), &QHeaderView::sectionResized, this, [this](int section) {
    // only a drag of the header by the user (the layout and the model resize sections too)
    if (fitting_ || reloading_ || !(QApplication::mouseButtons() & Qt::LeftButton)) return;
    QHeaderView *h = table_->horizontalHeader();
    int last = -1;  // the last shown column, the one the header stretches
    for (int v = h->count() - 1; v >= 0 && last < 0; v--)
      if (!h->isSectionHidden(h->logicalIndex(v))) last = h->logicalIndex(v);
    if (section != last) userWidths_ = true;
  });
  table_->setSortingEnabled(true);
  table_->sortByColumn(-1, Qt::AscendingOrder);
  table_->setSelectionBehavior(QAbstractItemView::SelectRows);
  table_->setSelectionMode(QAbstractItemView::ExtendedSelection);
  // a double click opens the record's pop-up (dev's segment dialogs); a cell is edited in place by F2 or a click on
  // the selected cell
  table_->setEditTriggers(QAbstractItemView::EditKeyPressed | QAbstractItemView::SelectedClicked);
  connect(table_, &QTableView::doubleClicked, this, [this](const QModelIndex &at) {
    if (at.isValid()) QTimer::singleShot(0, this, &RecordsPanel::editRecord);
  });
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
  editButton_ = button(tr("Edit..."), tr("The selected record in its pop-up (a double click on its row)"));
  editButton_->setObjectName("editRecord");
  connect(editButton_, &QPushButton::clicked, this, &RecordsPanel::editRecord);
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
  for (QPushButton *b : {addButton_, editButton_, duplicateButton_, removeButton_, upButton_, downButton_}) buttons->addWidget(b);
  buttons->addSpacing(12);
  buttons->addWidget(includeButton_);
  buttons->addWidget(excludeButton_);
  includeButton_->setVisible(model_->hasIncludeColumn());
  excludeButton_->setVisible(model_->hasIncludeColumn());
  buttons->addStretch(1);
  // what the section is: a "?" with its sentence as the tip, not a paragraph above the table (qt_gui.md 19.6)
  if (!help.isEmpty()) {
    QToolButton *about = new QToolButton;
    about->setObjectName("sectionHelp");
    about->setText("?");
    about->setToolTip(help);
    about->setAutoRaise(true);
    connect(about, &QToolButton::clicked, this, [about, help]() { QToolTip::showText(about->mapToGlobal(QPoint(0, about->height())), help, about); });
    buttons->addWidget(about);
  }
  // the form below the table: closed and opened here, or by dragging the splitter
  formToggle_ = new QToolButton;
  formToggle_->setObjectName("formToggle");
  formToggle_->setText(tr("Form"));
  formToggle_->setCheckable(true);
  formToggle_->setChecked(true);
  formToggle_->setToolTip(tr("Show or hide the form of the selected record below the table (a double click on a row "
                             "opens the record in its own window)."));
  connect(formToggle_, &QToolButton::toggled, this, &RecordsPanel::setFormShown);
  buttons->addWidget(formToggle_);
  // "?" and "Form" on the buttons' line at the buttons' height (they were smaller and sat higher)
  for (QToolButton *b : findChildren<QToolButton *>())
    if (b->objectName() == "sectionHelp" || b == formToggle_) b->setFixedHeight(addButton_->sizeHint().height());

  QWidget *top = new QWidget;
  QVBoxLayout *topLayout = new QVBoxLayout(top);
  topLayout->setContentsMargins(0, 0, 0, 0);
  topLayout->addWidget(filter_);
  topLayout->addWidget(table_, 1);
  topLayout->addLayout(buttons);
  emptyHint_ = new QLabel(tr("No record yet: <b>Add...</b> (Ins) makes one, its required fields filled with the project's "
                             "first choices."));
  emptyHint_->setObjectName("emptyHint");
  emptyHint_->setWordWrap(true);
  emptyHint_->setVisible(false);
  topLayout->addWidget(emptyHint_);
  topLayout->addStretch(0);

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
  formScroll_ = scroll;
  QWidget *bottom = new QWidget;
  bottom_ = bottom;
  QVBoxLayout *bottomLayout = new QVBoxLayout(bottom);
  bottomLayout->setContentsMargins(0, 0, 0, 0);
  bottomLayout->addWidget(kindLabel_);
  bottomLayout->addLayout(usedRow);
  bottomLayout->addWidget(scroll, 1);
  bottomLayout->addLayout(formButtons);

  // The table and the form share the height by a splitter (qt_gui.md section 19.3): the table keeps at least five
  // rows, the form may be closed (dragged to the bottom, or the Form button), and the sizes a user sets are kept per
  // panel.  Before, the form took 5/7 of the height and the table of a project's data sets showed two rows at 720 px.
  splitter_ = new QSplitter(Qt::Vertical);
  splitter_->setObjectName("split_" + id);
  splitter_->addWidget(top);
  splitter_->addWidget(bottom);
  splitter_->setStretchFactor(0, 1);
  splitter_->setStretchFactor(1, 1);
  splitter_->setCollapsible(0, false);
  splitter_->setCollapsible(1, true);
  connect(splitter_, &QSplitter::splitterMoved, this, [this]() {
    const QList<int> sizes = splitter_->sizes();
    if (sizes.size() == 2 && sizes[1] > 0) formHeight_ = sizes[1];
    QSignalBlocker block(formToggle_);
    formToggle_->setChecked(sizes.size() == 2 && sizes[1] > 0);
    QSettings().setValue("layout/" + splitter_->objectName(), splitter_->saveState());
  });
  QVBoxLayout *layout = new QVBoxLayout(this);
  layout->setContentsMargins(0, 0, 0, 0);
  layout->addWidget(splitter_);

  connect(holder_, &ProjectHolder::opened, this, &RecordsPanel::reload);
  connect(holder_, &ProjectHolder::sectionEdited, this, [this](const QString &s) {
    if (reloading_ || model_->writing()) return;
    if (s == section_) reload();
    else if (haveShown_ && s != section_) showUsedBy();  // a data set now names this target, ...
  });
  reload();
}

int RecordsPanel::recordCount() const { return model_->rowCount(); }

int RecordsPanel::tableRowsHeight(int rows) const {
  return table_->horizontalHeader()->sizeHint().height() + rows * table_->verticalHeader()->defaultSectionSize() +
         2 * table_->frameWidth() + table_->horizontalScrollBar()->sizeHint().height();
}

int RecordsPanel::visibleRows() const {
  return table_->viewport()->height() / std::max(1, table_->verticalHeader()->defaultSectionSize());
}

void RecordsPanel::showEvent(QShowEvent *event) {
  QWidget::showEvent(event);
  if (!sized_) initialSizes();
  if (!columnsFitted_) fitColumns();
}

bool RecordsPanel::eventFilter(QObject *watched, QEvent *event) {
  // the table's width changed (the window resized, the panel laid out): the columns fitted again, unless the user
  // set a column's width
  if (watched == table_->viewport() && event->type() == QEvent::Resize && !userWidths_ &&
      std::abs(table_->viewport()->width() - fittedWidth_) > 8) {
    columnsFitted_ = false;
    fitColumns();
  }
  return QWidget::eventFilter(watched, event);
}

void RecordsPanel::fitColumns() {
  // The columns as wide as their contents (at most 40 characters), the widest ones narrowed to one width when they
  // would not fit side by side (qt_gui.md section 19.8): no sideways scrolling at 1280 px while a column can give room.
  // A cell's whole text is in the row's tip.  Done once the table has its width; a column the user widens stays so.
  const int available = table_->viewport()->width();
  if (!table_->isVisible() || available < 200) return;
  columnsFitted_ = true;
  fittedWidth_ = available;
  fitting_ = true;
  // the contents' widths (the last column measured unstretched: it is stretched to the table's edge)
  const bool stretch = table_->horizontalHeader()->stretchLastSection();
  table_->horizontalHeader()->setStretchLastSection(false);
  table_->resizeColumnsToContents();
  QVector<int> widths;
  int total = 0;
  for (int c = 0; c < model_->columnCount(); c++) {
    const int w = table_->isColumnHidden(c) ? 0 : std::min(table_->columnWidth(c), fontMetrics().horizontalAdvance(QString(40, 'x')));
    widths << w;
    total += w;
  }
  int cap = fontMetrics().horizontalAdvance(QString(40, 'x'));
  const int floor = fontMetrics().horizontalAdvance(QString(12, 'x'));
  while (total > available && cap > floor) {
    cap -= 8;
    total = 0;
    for (int w : widths) total += std::min(w, cap);
  }
  for (int c = 0; c < model_->columnCount(); c++)
    if (!table_->isColumnHidden(c)) table_->setColumnWidth(c, std::min(widths[c], cap));
  table_->horizontalHeader()->setStretchLastSection(stretch);
  fitting_ = false;
}

void RecordsPanel::initialSizes() {
  const int total = splitter_->height();
  if (total < 100) return;  // not laid out yet
  sized_ = true;
  const QByteArray saved = QSettings().value("layout/" + splitter_->objectName()).toByteArray();
  if (!saved.isEmpty() && splitter_->restoreState(saved)) {
    QSignalBlocker block(formToggle_);
    formToggle_->setChecked(splitter_->sizes().value(1) > 0);
    return;
  }
  // the table: its rows up to ten (at least four), with the filter and the buttons; at most 60 % of the height
  const int rows = std::clamp(model_->rowCount(), 4, 10);
  const int chrome = filter_->sizeHint().height() + addButton_->sizeHint().height() + 24;
  const int top = std::min(tableRowsHeight(rows) + chrome, total * 3 / 5);
  splitter_->setSizes({top, total - top});
}

void RecordsPanel::setFormShown(bool on) {
  const QList<int> sizes = splitter_->sizes();
  if (sizes.size() != 2) return;
  const int total = sizes[0] + sizes[1];
  if (!on) {
    if (sizes[1] > 0) formHeight_ = sizes[1];
    splitter_->setSizes({total, 0});
  } else if (sizes[1] == 0) {
    const int form = formHeight_ > 0 ? formHeight_ : total / 2;
    splitter_->setSizes({std::max(total - form, tableRowsHeight(5)), form});
  }
  {
    QSignalBlocker block(formToggle_);
    formToggle_->setChecked(on);
  }
  QSettings().setValue("layout/" + splitter_->objectName(), splitter_->saveState());
}

bool RecordsPanel::formShown() const { return splitter_->sizes().value(1) > 0; }

void RecordsPanel::setExpertMode(bool on) {
  expert_ = on;
  fitting_ = true;  // columns hidden or shown here are not the user's resizing
  updateColumns();
  fitting_ = false;
  columnsFitted_ = false;
  userWidths_ = false;
  fitColumns();
}

void RecordsPanel::updateColumns() {
  // standard view: the Kind column when every row is of one kind, and the columns empty in every row (a "Fit?" column
  // without a box in any row too), hidden; a column that holds a value in any row stays (nothing in effect is hidden)
  const int rows = model_->rowCount();
  for (int c = 0; c < model_->columnCount(); c++) {
    const RecordTableModel::Column &col = model_->column(c);
    bool hide = false;
    if (!expert_ && rows > 0 && col.type == RecordTableModel::Column::Kind) {
      hide = true;
      const QString first = model_->data(model_->index(0, c)).toString();
      for (int r = 1; r < rows && hide; r++) hide = model_->data(model_->index(r, c)).toString() == first;
    } else if (!expert_ && rows > 0 && col.type == RecordTableModel::Column::Field) {
      hide = true;
      for (int r = 0; r < rows && hide; r++) hide = model_->data(model_->index(r, c)).toString().trimmed().isEmpty();
    } else if (!expert_ && rows > 0 && col.type == RecordTableModel::Column::Fit) {
      // a "Fit?" column whose field no row has (no box in any row)
      hide = true;
      for (int r = 0; r < rows && hide; r++) hide = !model_->data(model_->index(r, c), Qt::CheckStateRole).isValid();
    }
    table_->setColumnHidden(c, hide);
  }
}

QString RecordsPanel::shownName() const { return haveShown_ ? Q(shown_.name) : QString(); }

void RecordsPanel::reload() {
  const QString keep = shownName();
  const QStringList selected = selectedNames();
  reloading_ = true;
  model_->reload();
  table_->resizeColumnsToContents();
  for (int c = 0; c < model_->columnCount(); c++)
    table_->setColumnWidth(c, std::min(table_->columnWidth(c), fontMetrics().horizontalAdvance(QString(40, 'x'))));
  updateColumns();
  columnsFitted_ = false;
  userWidths_ = false;
  fitColumns();
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
  editButton_->setEnabled(one);
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
  kernels_ = front::KernelPreview();
  quantityRows_.clear();
  fieldRows_.clear();
  notApplying_.clear();
  hiddenByQuantity_.clear();
  haveShown_ = false;
  while (QLayoutItem *item = formLayout_->takeAt(0)) {
    if (QWidget *w = item->widget()) delete w;
    delete item;
  }
  // an empty section: the hint under the table, a table of a few rows, no form
  const bool empty = model_->rowCount() == 0;
  kindLabel_->setText(empty ? QString() : tr("Choose a record in the table."));
  emptyHint_->setVisible(empty);
  bottom_->setVisible(!empty);
  table_->setMaximumHeight(empty ? table_->horizontalHeader()->sizeHint().height() + 3 * table_->verticalHeader()->defaultSectionSize() + 4
                                 : QWIDGETSIZE_MAX);
  // a table with records keeps five rows in sight, whatever the splitter (qt_gui.md section 19.3)
  table_->setMinimumHeight(empty ? 0 : tableRowsHeight(std::min(5, std::max(1, model_->rowCount()))));
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
  applyQuantityRules();
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
    if (!previewPending_) {
      previewPending_ = true;
      QTimer::singleShot(0, this, [this]() {
        previewPending_ = false;
        updatePreview();
      });
    }
  });
  if (Q(field.key) == quantityKey_) connect(e, &FieldEditor::changed, this, &RecordsPanel::quantityChanged);
  connect(e, &FieldEditor::changed, this, &RecordsPanel::applyConditions);
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
  QLabel *label = new QLabel(DisplayTitle(field.title) + (field.required ? " *" : ""));
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
      remove->setFixedSize(22, 22);
      // the remove button left of the entry: every field's text button stays in one column at the right edge
      // (qt_gui.md section 19.9)
      h->addWidget(remove, 0, Qt::AlignTop);
      h->addWidget(e, 1);
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
  if (!fieldRows_.contains(key)) fieldRows_.insert(key, {form, label, row, field});
  if (field.byQuantity) quantityRows_.insert(key, {form, label, row, field});
  if (field.rare && values.isEmpty()) {
    rareRows_ << label << row;
    label->setVisible(everyField_->isChecked());
    row->setVisible(everyField_->isChecked());
  }
}

void RecordsPanel::buildForm(const front::Record &record) {
  clearForm();
  if (record.name != shown_.name) {
    previewUser_.clear();
    previewEnergy_ = 0;
  }
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
  basic->setContentsMargins(0, 0, 0, 0);  // the same right edge as the boxes' forms (the text buttons in one column)
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
  formLayout_->insertWidget(1, preview_ = new KernelView([this](const QString &user, double e) { setPreviewPoint(user, e); }));
  updatePreview();
  formLayout_->addStretch(1);
  applyQuantityRules();
  applyConditions();
  showUsedBy();
}

namespace {

/// A field's schema condition in the words experiment::Describe writes for the two simple kinds: "<field> is[ not]
/// <word> or <word> ..." and "<field> is[ not] given".  False for any other condition (it is then taken to hold).
bool ParseCondition(const std::string &words, QString &field, QStringList &values, bool &negate, bool &present) {
  const QRegularExpression re("^([A-Za-z_][A-Za-z0-9_]*) is( not)? ([A-Za-z0-9_.+\\- ]+)$");
  const QRegularExpressionMatch m = re.match(Q(words));
  if (!m.hasMatch()) return false;
  field = m.captured(1);
  negate = !m.captured(2).isEmpty();
  const QString rest = m.captured(3).trimmed();
  present = rest == "given";
  values = present ? QStringList() : rest.split(" or ", Qt::SkipEmptyParts);
  return present || !values.isEmpty();
}

}  // namespace

void RecordsPanel::applyConditions() {
  // A field the schema requires only for some records ("Active nucleus: kind is layer, infinite or sample") is not
  // marked missing on a record it does not apply to, and its row is hidden while it holds no value (qt_gui.md
  // section 19.8).  Before, a valid target of kind none showed red "Required:" lines.
  if (!haveShown_) return;
  const front::KindInfo *kind = holder_->kind(shown_.kind);
  if (!kind) return;
  auto valueOf = [this](const QString &key) {
    for (FieldEditor *e : editors_.value(key)) {
      const QString t = e->text().trimmed();
      if (!t.isEmpty()) return t.section(' ', 0, 0);
    }
    return QString();
  };
  auto requiredField = [kind](const QString &key) {
    for (const front::FieldInfo &f : kind->fields)
      if (Q(f.key) == key) return f.required && f.appliesWhen.empty();
    return false;
  };
  notApplying_.clear();
  for (const front::FieldInfo &f : kind->fields) {
    if (f.appliesWhen.empty()) continue;
    QString field;
    QStringList values;
    bool negate = false, present = false;
    if (!ParseCondition(f.appliesWhen, field, values, negate, present)) continue;
    const QString value = valueOf(field);
    // the validator's rule (experiment Validate.cpp, Holds): a field not given is "not is", unless it is a required
    // field still missing, whose answer is unknown (taken to hold here: the field stays marked)
    bool holds;
    if (present) holds = value.isEmpty() == negate;
    else if (value.isEmpty()) holds = requiredField(field) || negate;
    else holds = values.contains(value.startsWith('"') ? value.mid(1).chopped(1) : value) != negate;
    const QString key = Q(f.key);
    bool given = false;
    for (FieldEditor *e : editors_.value(key)) {
      e->setApplies(holds);
      given = given || !e->text().trimmed().isEmpty();
    }
    if (!holds) notApplying_ << key;
    if (fieldRows_.contains(key)) {
      const QuantityRow &r = fieldRows_[key];
      r.label->setText(DisplayTitle(f.title) + (f.required && holds ? " *" : ""));
      const bool rareHidden = rareRows_.contains(r.row) && !everyField_->isChecked() && !given;
      const bool hide = !holds && !given && !everyField_->isChecked();
      if (!quantityRows_.contains(key)) r.form->setRowVisible(r.row, !rareHidden && !hide);
    }
    if (QWidget *box = boxes_.value(key, nullptr)) {
      bool used = given;
      for (const front::FieldInfo &m : kind->fields)
        if (m.group == f.key)
          for (FieldEditor *e : editors_.value(Q(m.key))) used = used || !e->text().trimmed().isEmpty();
      const bool rareHidden = rareRows_.contains(box) && !everyField_->isChecked() && !used;
      box->setVisible(!rareHidden && (holds || used || everyField_->isChecked()));
    }
  }
}

void RecordsPanel::updatePreview() {
  if (!preview_ || !haveShown_) return;
  // only a record other records name by reference (a target, a beam, a detector) folds their points
  bool referable = false;
  for (const front::KindInfo &k : holder_->kinds())
    for (const front::FieldInfo &f : k.fields)
      referable = referable || (f.value.type == "reference" &&
                                std::find(f.value.values.begin(), f.value.values.end(), shown_.kind) != f.value.values.end());
  std::string why;
  const bool drawn = referable && front::KernelPreviewOf(holder_->project(), section_.toStdString(), formRecord(),
                                                         previewUser_.toStdString(),
                                                         previewEnergy_ > 0 ? previewEnergy_ : std::nan(""), 301, kernels_, why);
  previewWhy_ = drawn ? QString() : Q(why);
  if (!drawn) kernels_ = front::KernelPreview();
  // a record that selects a kernel the engine cannot draw now (the project does not build): the reason in its place
  // (only a reason about this record: an engine message about another record, as a data set's field that does not
  // apply, is that record's matter, shown on its form and in the Diagnostics, not here; qt_gui.md section 19.8)
  const bool reason = referable && !drawn && !why.empty() && kernels_.notes.empty() &&
                      why.find("selects no kernel") == std::string::npos &&
                      why.find(shown_.name) != std::string::npos;
  if (reason) kernels_.notes.push_back(tr("Not drawn: %1").arg(Q(why)).toStdString());
  static_cast<KernelView *>(preview_)->setPreview(kernels_);
  preview_->setVisible(drawn || reason);
}

void RecordsPanel::setPreviewPoint(const QString &user, double labEnergy) {
  previewUser_ = user;
  previewEnergy_ = labEnergy > 0 ? labEnergy : 0.0;
  updatePreview();
}

bool RecordsPanel::quantityEdited() const {
  const FieldEditor *e = quantityKey_.isEmpty() ? nullptr : fieldEditor(quantityKey_);
  return e && e->isEdited();
}

void RecordsPanel::applyQuantityRules() {
  hiddenByQuantity_.clear();
  if (quantityRows_.isEmpty()) return;
  std::vector<front::QuantityFieldRule> rules;
  const bool known = front::QuantityFieldRules(quantityText().trimmed().toStdString(), rules);
  const bool changed = quantityEdited();
  for (auto it = quantityRows_.begin(); it != quantityRows_.end(); ++it) {
    const front::QuantityFieldRule *rule = nullptr;
    for (const front::QuantityFieldRule &r : rules)
      if (QString::fromStdString(r.field) == it.key()) rule = &r;
    bool value = false;
    for (FieldEditor *e : editors_.value(it.key())) value = value || !e->text().trimmed().isEmpty();
    // unknown quantity (being typed): nothing hidden; a value the quantity does not use stays shown (and the
    // engine's refusal beside it) until the quantity is changed in the form, then it is dropped on Apply
    const bool allowed = !known || rule;
    const bool shown = allowed || (value && !changed);
    if (!allowed && changed) hiddenByQuantity_ << it.key();
    const bool required = rule && rule->required && !rule->conditional;
    it->label->setText(DisplayTitle(it->field.title) + (required || it->field.required ? " *" : ""));
    QString tip = ToolTip(it->field);
    if (rule && rule->conditional)
      tip += "<br>" + tr("Quantity '%1' uses it only when %2.").arg(quantityText().toHtmlEscaped(), QString::fromStdString(rule->when).toHtmlEscaped());
    else if (rule && rule->required)
      tip += "<br>" + tr("Quantity '%1' requires it.").arg(quantityText().toHtmlEscaped());
    it->label->setToolTip(tip);
    const bool rare = rareRows_.contains(it->row) && !everyField_->isChecked() && !value;
    it->form->setRowVisible(it->row, shown && !rare);
    // a field of its own box (a structured value): the box goes with it
    if (QWidget *box = boxes_.value(it.key(), nullptr))
      box->setVisible(shown && !(rareRows_.contains(box) && !everyField_->isChecked() && !value));
  }
}

bool RecordsPanel::fieldRowShown(const QString &key) const {
  if (quantityRows_.contains(key)) {
    const QuantityRow &r = quantityRows_[key];
    return r.row && r.row->isVisibleTo(formHost_);
  }
  const FieldEditor *e = fieldEditor(key);
  return e && e->isVisibleTo(formHost_);
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
  applyQuantityRules();
  applyConditions();
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
      if (!t.isEmpty() && !hiddenByQuantity_.contains(it.key())) texts[it.key()] << t;
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
  QStringList dropped;
  for (const QString &key : hiddenByQuantity_)
    for (FieldEditor *e : editors_.value(key))
      if (!e->text().trimmed().isEmpty() && !dropped.contains(key)) dropped << key;
  const QString droppedQuantity = quantityText();
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
  if (ok && !dropped.isEmpty())
    status_->setText(tr("Dropped %1: quantity '%2' does not use %3.")
                         .arg(dropped.join(", "), droppedQuantity, dropped.size() > 1 ? tr("them") : tr("it")));
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
  // the complete pop-up of dev's segment dialogs: every usual field; OK and "Add another" write after the engine's check
  RecordDialog dialog(holder_, section_, kinds_, [this](const front::Record &r, QString &why) {
    if (!holder_->setRecord(section_, r, why)) return false;
    selectRecord(Q(r.name));
    emit edited();
    return true;
  }, this);
  dialog.setWindowTitle(tr("Add a record"));
  dialog.exec();
}

void RecordsPanel::editRecord() {
  if (!haveShown_) return;
  // the form's own edits first (else the pop-up would show the record without them)
  QString why;
  if (formEdited() && !applyForm(why, true)) {
    status_->setText(why);
    return;
  }
  const front::Record shown = shown_;
  RecordDialog dialog(holder_, section_, kinds_, [this, shown](const front::Record &r, QString &why) {
    front::Record same = r;
    same.name = shown.name;
    reloading_ = true;
    bool ok = holder_->setRecord(section_, same, why);
    if (ok && r.name != shown.name) ok = holder_->renameRecord(section_, Q(shown.name), Q(r.name), why);
    reloading_ = false;
    haveShown_ = false;
    reload();
    selectRecord(Q(ok ? r.name : shown.name));
    if (ok) emit edited();
    return ok;
  }, this);
  dialog.setRecord(shown);
  dialog.exec();
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
