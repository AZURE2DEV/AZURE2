#ifndef DIALOGSCROLLAREA_H
#define DIALOGSCROLLAREA_H

#include <QScrollArea>
#include <QScrollBar>
#include <QGuiApplication>
#include <QScreen>
#include <QDialog>
#include <QFrame>
#include <QEvent>
#include <QLayout>

/*!
 * A scroll area that asks for its content's size, so a dialog built around it
 * opens at the content's natural size and only grows a scroll bar when the
 * content is taller than the screen.  QScrollArea's own sizeHint ignores the
 * content, which would open the dialog as a small box.
 */
class DialogScrollArea : public QScrollArea {
 public:
  explicit DialogScrollArea(QWidget *content, QWidget *parent = 0) : QScrollArea(parent) {
    setWidget(content);
    setWidgetResizable(true);
    setFrameShape(QFrame::NoFrame);
    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  }
  QSize sizeHint() const override {
    if (!widget()) return QScrollArea::sizeHint();
    QSize s = widget()->sizeHint();
    // Room for the vertical bar when it appears, so the content never has to shrink sideways.
    s.rwidth() += verticalScrollBar()->sizeHint().width() + 2 * frameWidth();
    s.rheight() += 2 * frameWidth();
    return s;
  }
  QSize minimumSizeHint() const override {
    QSize s = sizeHint();
    s.setHeight(120);
    return s;
  }

 protected:
  // A scroll area absorbs its content's layout changes (its own hint does not
  // depend on them), so the dialog's layout would keep the hint it cached when
  // the sections were still hidden.  Pass the change on.
  bool eventFilter(QObject *o, QEvent *e) override {
    const bool handled = QScrollArea::eventFilter(o, e);
    if (o == widget() && e->type() == QEvent::LayoutRequest) updateGeometry();
    return handled;
  }
};

/*!
 * Keep a dialog within the screen that shows it: the segment dialogs gained
 * UPOS, composite and mixture sections and outgrew small displays (DeBoer,
 * 2026-10-12).  Call after the layout is set; adjustSize() then stops at this
 * height and the scroll area takes the rest.
 */
inline void capDialogToScreen(QDialog *dialog) {
  const QScreen *screen = dialog->screen() ? dialog->screen() : QGuiApplication::primaryScreen();
  if (!screen) return;
  const QRect avail = screen->availableGeometry();
  dialog->setMaximumHeight(avail.height() - 80);
  dialog->setMaximumWidth(avail.width() - 40);
}

/*!
 * Size a dialog to its content, within the cap above.  Used instead of
 * QWidget::adjustSize(): that caps a window at two thirds of the screen, and
 * it trusts the layout's cached hint, which is stale right after a section
 * was shown or hidden (the scroll area absorbs the content's layout request).
 */
inline void fitDialogToContent(QDialog *dialog) {
  // Every layout keeps a cache, and every widget's entry in its parent's
  // layout keeps another; a section shown or hidden refreshes them only as
  // the event loop delivers layout requests level by level.  Refresh them all
  // now, bottom to top, so the hint read below is the content's current one.
  // Nested layouts (a QGridLayout inside another, as the segment dialogs build them) are not any widget's
  // layout() and keep their own cached hints, so a section hidden inside one left the old height: the dialog
  // stayed tall and the spare height went to the energy and angle boxes (DeBoer, 2026-10-25).
  const QList<QLayout *> layouts = dialog->findChildren<QLayout *>();
  for (int i = layouts.size() - 1; i >= 0; i--) layouts[i]->invalidate();
  const QList<QWidget *> all = dialog->findChildren<QWidget *>();
  for (int i = all.size() - 1; i >= 0; i--) {
    if (QLayout *l = all[i]->layout()) l->invalidate();
    all[i]->updateGeometry();
  }
  if (QLayout *l = dialog->layout()) l->invalidate();
  QSize s = dialog->sizeHint();
  s = s.boundedTo(QSize(dialog->maximumWidth(), dialog->maximumHeight()));
  dialog->resize(s);
}

#endif
