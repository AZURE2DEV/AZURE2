#include "SectionTextPanel.h"

#include <QFontDatabase>
#include <QHBoxLayout>
#include <QLabel>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QVBoxLayout>

SectionTextPanel::SectionTextPanel(ProjectHolder *holder, const QString &section, const QString &help, QWidget *parent)
    : QWidget(parent), holder_(holder), section_(section) {
  setObjectName("section_" + section);
  edit_ = new QPlainTextEdit;
  edit_->setObjectName("sectionText_" + section);
  edit_->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
  edit_->setLineWrapMode(QPlainTextEdit::NoWrap);
  QPushButton *revert = new QPushButton(tr("Revert"));
  QPushButton *applyButton = new QPushButton(tr("&Apply"));
  connect(revert, &QPushButton::clicked, this, &SectionTextPanel::reload);
  connect(applyButton, &QPushButton::clicked, this, &SectionTextPanel::apply);
  QHBoxLayout *buttons = new QHBoxLayout;
  // what the section is: the tip of its name, not a paragraph above the text (qt_gui.md section 19.6)
  QLabel *name = new QLabel(QString("&lt;%1&gt; ... &lt;/%1&gt;").arg(section));
  name->setToolTip(help);
  edit_->setToolTip(help);
  buttons->addWidget(name);
  buttons->addStretch(1);
  buttons->addWidget(revert);
  buttons->addWidget(applyButton);
  QVBoxLayout *layout = new QVBoxLayout(this);
  layout->addWidget(edit_, 1);
  layout->addLayout(buttons);
  connect(holder_, &ProjectHolder::opened, this, &SectionTextPanel::reload);
  connect(holder_, &ProjectHolder::sectionEdited, this, [this](const QString &s) {
    if (s == section_ && !applying_) reload();
  });
  reload();
}

QString SectionTextPanel::text() const { return edit_->toPlainText(); }
void SectionTextPanel::setText(const QString &text) { edit_->setPlainText(text); }

void SectionTextPanel::reload() { edit_->setPlainText(holder_->sectionBody(section_).join("\n")); }

void SectionTextPanel::apply() {
  QStringList lines = edit_->toPlainText().split('\n');
  if (lines.size() == 1 && lines.first().isEmpty()) lines.clear();
  if (lines == holder_->sectionBody(section_)) return;
  applying_ = true;
  holder_->setSectionBody(section_, lines);
  applying_ = false;
}
