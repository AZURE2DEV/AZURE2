#include "DiagnosticsPanel.h"

#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QPushButton>
#include <QStyle>
#include <QTreeWidget>
#include <QVBoxLayout>

DiagnosticsPanel::DiagnosticsPanel(ProjectHolder *holder, QWidget *parent) : QWidget(parent), holder_(holder) {
  setObjectName("diagnostics");
  tree_ = new QTreeWidget;
  tree_->setObjectName("diagnosticsTree");
  tree_->setHeaderLabels({tr("Severity"), tr("Code"), tr("Line"), tr("Record"), tr("Message")});
  tree_->setRootIsDecorated(false);
  tree_->header()->setSectionResizeMode(4, QHeaderView::Stretch);
  tree_->setToolTip(tr("Double-click a diagnostic to go to its tab and record."));
  connect(tree_, &QTreeWidget::itemDoubleClicked, this, [this](QTreeWidgetItem *item) {
    emit showRecord(item->data(0, Qt::UserRole).toString(), item->data(1, Qt::UserRole).toString());
  });
  summary_ = new QLabel(tr("Not validated yet."));
  QPushButton *button = new QPushButton(tr("&Validate"));
  button->setToolTip(tr("Check the project as it stands with the engine (AZURE2 --validate)."));
  connect(button, &QPushButton::clicked, this, &DiagnosticsPanel::validate);
  QHBoxLayout *top = new QHBoxLayout;
  top->addWidget(summary_, 1);
  top->addWidget(button);
  QVBoxLayout *layout = new QVBoxLayout(this);
  layout->setContentsMargins(4, 4, 4, 4);
  layout->addLayout(top);
  layout->addWidget(tree_);
}

QString DiagnosticsPanel::summary() const { return summary_->text(); }

void DiagnosticsPanel::activateRow(int row) {
  if (QTreeWidgetItem *item = tree_->topLevelItem(row)) emit tree_->itemDoubleClicked(item, 0);
}

int DiagnosticsPanel::rowCount() const { return tree_->topLevelItemCount(); }

QStringList DiagnosticsPanel::codes() const {
  QStringList out;
  for (int i = 0; i < tree_->topLevelItemCount(); i++) out << tree_->topLevelItem(i)->text(1);
  return out;
}

bool DiagnosticsPanel::validate() {
  tree_->clear();
  errors_ = 0;
  warnings_ = 0;
  const azure2::front::Validation v = holder_->validate();
  if (!v.format3) {
    summary_->setText(QString::fromStdString(v.refusal));
    errors_ = 1;
    emit validated(false, errors_, 0);
    return false;
  }
  int &warnings = warnings_;
  const QIcon errorIcon = style()->standardIcon(QStyle::SP_MessageBoxCritical);
  const QIcon warningIcon = style()->standardIcon(QStyle::SP_MessageBoxWarning);
  for (const azure2::front::Diagnostic &d : v.diagnostics) {
    QTreeWidgetItem *item = new QTreeWidgetItem(
        {QString::fromStdString(d.severity), QString::fromStdString(d.code), d.line ? QString::number(d.line) : QString(),
         QString::fromStdString(d.record + (d.field.empty() ? "" : "." + d.field)), QString::fromStdString(d.message)});
    item->setToolTip(4, QString::fromStdString(d.text));
    item->setData(0, Qt::UserRole, QString::fromStdString(d.section));
    item->setData(1, Qt::UserRole, QString::fromStdString(d.record));
    if (d.severity == "error") {
      errors_++;
      item->setForeground(0, QColor(0xd0, 0x30, 0x30));
      item->setIcon(0, errorIcon);
    } else {
      warnings++;
      item->setIcon(0, warningIcon);
    }
    tree_->addTopLevelItem(item);
  }
  summary_->setText(v.ok ? tr("Valid: %1 warning(s).").arg(warnings)
                         : tr("Not valid: %1 error(s), %2 other diagnostic(s).").arg(errors_).arg(warnings));
  emit validated(v.ok, errors_, warnings);
  return v.ok;
}
