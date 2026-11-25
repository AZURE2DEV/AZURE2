// The notes about a project just opened (NotesBar.h).
#include "NotesBar.h"

#include <QAction>
#include <QDialog>
#include <QDialogButtonBox>
#include <QLabel>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QStyle>
#include <QVBoxLayout>

NotesBar::NotesBar(QWidget *parent) : QToolBar(tr("Notes about the project"), parent) {
  setObjectName("notesBar");
  setMovable(false);
  icon_ = new QLabel;
  icon_->setObjectName("notesIcon");
  icon_->setPixmap(style()->standardIcon(QStyle::SP_MessageBoxInformation).pixmap(16, 16));
  icon_->setContentsMargins(8, 0, 4, 0);
  addWidget(icon_);
  text_ = new QLabel;
  text_->setObjectName("notesText");
  text_->setContentsMargins(0, 0, 12, 0);
  addWidget(text_);
  locateAction_ = addAction(tr("Locate Data Files..."));
  locateAction_->setObjectName("notesLocate");
  locateAction_->setToolTip(tr("Choose a folder: each data file the project names but this computer does not have at "
                               "that path is looked for there by its name, and the data set is pointed to it."));
  locateAction_->setVisible(false);
  connect(locateAction_, &QAction::triggered, this, &NotesBar::locateRequested);
  detailsAction_ = addAction(tr("Details"));
  detailsAction_->setObjectName("notesDetails");
  connect(detailsAction_, &QAction::triggered, this, &NotesBar::showDetails);
  dismissAction_ = addAction(tr("Dismiss"));
  dismissAction_->setObjectName("notesDismiss");
  connect(dismissAction_, &QAction::triggered, this, &QWidget::hide);
  hide();
}

void NotesBar::showLocate(bool on) { locateAction_->setVisible(on); }

void NotesBar::clearNotes() {
  locateAction_->setVisible(false);
  summaries_.clear();
  texts_.clear();
  aboutDiagnostics_ = false;
  warning_ = false;
  icon_->setPixmap(style()->standardIcon(QStyle::SP_MessageBoxInformation).pixmap(16, 16));
  if (details_) details_->hide();
  hide();
}

void NotesBar::addNote(const QString &summary, const QString &details, bool diagnostics, bool warning) {
  // a warning first (the one line the bar shows), the others in the order they came
  const int at = warning ? 0 : summaries_.size();
  summaries_.insert(at, summary);
  texts_.insert(at, details.isEmpty() ? summary : details);
  aboutDiagnostics_ = aboutDiagnostics_ || diagnostics;
  if (warning && !warning_) {
    warning_ = true;
    icon_->setPixmap(style()->standardIcon(QStyle::SP_MessageBoxWarning).pixmap(16, 16));
  }
  updateText();
  show();
  if (details_ && details_->isVisible()) showDetails();
}

void NotesBar::updateText() {
  const QString first = summaries_.isEmpty() ? QString() : summaries_.first();
  text_->setText(summaries_.size() <= 1 ? first : tr("%1 notes: %2").arg(summaries_.size()).arg(first));
  text_->setToolTip(summaries_.join("\n"));
}

void NotesBar::showDetails() {
  if (!details_) {
    // a window of its own, not modal: the main window stays usable while it is open
    details_ = new QDialog(window());
    details_->setObjectName("notesDetailsWindow");
    details_->setWindowTitle(tr("Notes about the project"));
    details_->setModal(false);
    details_->resize(640, 360);
    QVBoxLayout *v = new QVBoxLayout(details_);
    detailsText_ = new QPlainTextEdit;
    detailsText_->setReadOnly(true);
    v->addWidget(detailsText_);
    QDialogButtonBox *buttons = new QDialogButtonBox(QDialogButtonBox::Close);
    diagnosticsButton_ = buttons->addButton(tr("Show Diagnostics"), QDialogButtonBox::ActionRole);
    diagnosticsButton_->setObjectName("notesShowDiagnostics");
    connect(diagnosticsButton_, &QPushButton::clicked, this, &NotesBar::diagnosticsRequested);
    connect(buttons, &QDialogButtonBox::rejected, details_, &QDialog::hide);
    v->addWidget(buttons);
  }
  diagnosticsButton_->setVisible(aboutDiagnostics_);
  detailsText_->setPlainText(texts_.join("\n\n"));
  details_->show();
  details_->raise();
}
