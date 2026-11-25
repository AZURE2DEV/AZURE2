#ifndef EXPERTMODE_H
#define EXPERTMODE_H

#include <QSettings>

/*!
 * Standard or expert view (docs/upgrade/qt_gui.md section 19.4): "Show expert options" in the View menu and on the
 * start page, off by default, kept in QSettings ("view/expertOptions").  It changes what the window shows, never the
 * project nor what the engine computes: in standard view the options an occasional user does not need are hidden, and
 * an expert option that holds a value other than its default (a record, a key set, a section not empty) is shown
 * anyway, marked, so that nothing hidden is in effect.  Each tab answers setExpertMode(bool); the window calls them.
 */
namespace expert {

inline bool Enabled() { return QSettings().value("view/expertOptions", false).toBool(); }
inline void SetEnabled(bool on) { QSettings().setValue("view/expertOptions", on); }
/// The tip added to an expert option shown in standard view because the project sets it.
inline QString ShownBecauseSet() {
  return QObject::tr("An expert option, shown because the project sets it (View > Show Expert Options shows them all).");
}

}  // namespace expert

#endif
