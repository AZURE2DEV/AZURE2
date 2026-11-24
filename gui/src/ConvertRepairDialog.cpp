// A line of an old project the converter could not read: shown, edited, a reading taken (ConvertRepairDialog.h).
#include "ConvertRepairDialog.h"

#include <QDialogButtonBox>
#include <QFileInfo>
#include <QFontDatabase>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QVBoxLayout>

ConvertRepairDialog::ConvertRepairDialog(const QString &file, const azure2::front::ConvertLineProblem &problem, QWidget *parent)
    : QDialog(parent) {
  setObjectName("convertRepairDialog");
  setWindowTitle(tr("A line the converter cannot read"));
  resize(820, 360);
  QVBoxLayout *v = new QVBoxLayout(this);
  QLabel *where = new QLabel(tr("<b>%1</b>, line %2, in &lt;%3&gt;:")
                                 .arg(QFileInfo(file).fileName().toHtmlEscaped())
                                 .arg(problem.line)
                                 .arg(QString::fromStdString(problem.section)));
  v->addWidget(where);
  edit_ = new QLineEdit(QString::fromStdString(problem.text).simplified());
  edit_->setObjectName("convertRepairLine");
  edit_->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
  edit_->setToolTip(tr("The line as the converted copy will hold it: edit it, or take a reading below."));
  v->addWidget(edit_);
  QLabel *why = new QLabel(QString::fromStdString(problem.problem).section('\n', 0, 0));
  why->setWordWrap(true);
  why->setStyleSheet("color: #c0392b;");
  v->addWidget(why);
  QLabel *expected = new QLabel(tr("Expected: %1").arg(QString::fromStdString(problem.expected).toHtmlEscaped()));
  expected->setWordWrap(true);
  expected->setProperty("role", "hint");
  v->addWidget(expected);
  QGroupBox *box = new QGroupBox(problem.readings.empty() ? tr("No reading can be proposed: edit the line")
                                                          : tr("Possible readings (a click puts it in the line above)"));
  QVBoxLayout *b = new QVBoxLayout(box);
  for (const auto &r : problem.readings) {
    QPushButton *button = new QPushButton(QString::fromStdString(r.first));
    button->setToolTip(QString::fromStdString(r.second));
    button->setStyleSheet("text-align: left; padding: 4px 8px;");
    const QString line = QString::fromStdString(r.second);
    connect(button, &QPushButton::clicked, this, [this, line]() { edit_->setText(line); });
    b->addWidget(button);
    readings_ << button;
  }
  v->addWidget(box);
  v->addStretch(1);
  QLabel *kept = new QLabel(tr("The old file is kept: the line goes into a copy beside it, and the copy is converted."));
  kept->setProperty("role", "hint");
  v->addWidget(kept);
  QDialogButtonBox *buttons = new QDialogButtonBox(QDialogButtonBox::Cancel);
  retry_ = buttons->addButton(tr("Convert again"), QDialogButtonBox::AcceptRole);
  retry_->setObjectName("convertRepairRetry");
  retry_->setProperty("primary", true);
  connect(buttons, &QDialogButtonBox::accepted, this, &QDialog::accept);
  connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
  v->addWidget(buttons);
}

QString ConvertRepairDialog::lineText() const { return edit_->text(); }
