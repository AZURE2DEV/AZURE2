#ifndef SECTIONTEXTPANEL_H
#define SECTIONTEXTPANEL_H

#include <QWidget>

#include "ProjectHolder.h"

QT_BEGIN_NAMESPACE
class QPlainTextEdit;
QT_END_NAMESPACE

/*!
 * A section the engine reads with a grammar of its own (<thm>, <parameterSettings>, <mcmc>): its lines as text, written
 * back with Apply (only that section's body changes; a section the file has not is added at its end when the text is
 * not empty).  The engine's diagnostics name a wrong line (Validate, the window's panel); nothing is parsed here.
 */
class SectionTextPanel : public QWidget {
  Q_OBJECT

 public:
  SectionTextPanel(ProjectHolder *holder, const QString &section, const QString &help, QWidget *parent = nullptr);
  QString text() const;
  void setText(const QString &text);
  /// Write the text into the project (Save calls it); no change when the text is the section's.
  void apply();

 public slots:
  void reload();

 private:
  ProjectHolder *holder_;
  QString section_;
  QPlainTextEdit *edit_;
  bool applying_ = false;
};

#endif
