#ifndef RECORDTABS_H
#define RECORDTABS_H

#include <QList>
#include <QWidget>

#include "ProjectHolder.h"

class RecordsPanel;
class SectionTextPanel;

QT_BEGIN_NAMESPACE
class QTabWidget;
QT_END_NAMESPACE

/*!
 * A tab of the window made of format-3 sections (docs/upgrade/qt_gui.md section 2): one RecordsPanel per section of
 * records and one SectionTextPanel per section of its own grammar, as sub-tabs when there are several.  The tabs the
 * old program had keep their names: Segments (Segments From Data, Segments Without Data, and the model data),
 * Experimental Effects (targets, beams, detectors), Nuclear Potential (the THM tab is ThmTab).
 */
class SectionsTab : public QWidget {
  Q_OBJECT

 public:
  struct Part {
    QString section;
    QString title;
    QString help;
    bool text = false;  ///< a section of its own grammar (SectionTextPanel)
  };
  SectionsTab(ProjectHolder *holder, const QList<Part> &parts, QWidget *parent = nullptr);
  /// The panel of a section (nullptr when this tab has none).
  RecordsPanel *records(const QString &section) const;
  SectionTextPanel *text(const QString &section) const;
  /// Show the sub-tab of a section and select a record in it (from a diagnostic).
  bool show(const QString &section, const QString &record);
  /// Write the text panels into the project (before a save).
  void applyText();

 private:
  QTabWidget *tabs_ = nullptr;
  QList<RecordsPanel *> records_;
  QList<SectionTextPanel *> texts_;
};

/// The parts of the window's tabs, in the old program's order.
QList<SectionsTab::Part> SegmentsParts();
QList<SectionsTab::Part> EffectsParts();
QList<SectionsTab::Part> PotentialParts();

#endif
