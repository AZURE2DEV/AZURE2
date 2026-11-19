#include "RecordTabs.h"

#include <QTabWidget>
#include <QVBoxLayout>

#include "RecordsPanel.h"
#include "SectionTextPanel.h"

SectionsTab::SectionsTab(ProjectHolder *holder, const QList<Part> &given, QWidget *parent) : QWidget(parent) {
  QVBoxLayout *layout = new QVBoxLayout(this);
  // a part of one sub-tab per kind: the kinds of its section, in the schema's order, titled by the schema
  QList<Part> parts;
  QList<QStringList> kindsOf;
  for (const Part &p : given) {
    if (!p.perKind) {
      parts << p;
      kindsOf << QStringList();
      continue;
    }
    for (const QString &kind : holder->kindsOfSection(p.section)) {
      Part one = p;
      const azure2::front::KindInfo *k = holder->kind(kind.toStdString());
      one.title = k ? QString::fromStdString(k->title) + "s" : kind;
      one.help = p.help;
      parts << one;
      kindsOf << QStringList{kind};
    }
  }
  if (parts.size() > 1) {
    tabs_ = new QTabWidget;
    layout->addWidget(tabs_);
  }
  for (int i = 0; i < parts.size(); i++) {
    const Part &p = parts[i];
    QWidget *panel = nullptr;
    if (p.text) {
      SectionTextPanel *t = new SectionTextPanel(holder, p.section, p.help);
      texts_ << t;
      panel = t;
    } else {
      RecordsPanel *r = new RecordsPanel(holder, p.section, p.help, nullptr, kindsOf[i]);
      records_ << r;
      panel = r;
    }
    if (tabs_) tabs_->addTab(panel, p.title);
    else layout->addWidget(panel);
  }
}

RecordsPanel *SectionsTab::records(const QString &section) const {
  for (RecordsPanel *p : records_)
    if (p->section() == section) return p;
  return nullptr;
}

RecordsPanel *SectionsTab::recordsOfKind(const QString &section, const QString &kind) const {
  for (RecordsPanel *p : records_)
    if (p->section() == section && p->kinds().contains(kind)) return p;
  return nullptr;
}

RecordsPanel *SectionsTab::panelFor(const QString &section, const QString &name) const {
  for (RecordsPanel *p : records_)
    if (p->section() == section && p->recordNames().contains(name)) return p;
  return records(section);
}

SectionTextPanel *SectionsTab::text(const QString &section) const {
  for (SectionTextPanel *p : texts_)
    if (p->objectName() == "section_" + section) return p;
  return nullptr;
}

bool SectionsTab::show(const QString &section, const QString &record) {
  QWidget *panel = panelFor(section, record);
  if (!panel) panel = text(section);
  if (!panel) return false;
  if (tabs_) tabs_->setCurrentWidget(panel);
  if (RecordsPanel *r = panelFor(section, record)) r->selectRecord(record);
  return true;
}

void SectionsTab::applyText() {
  for (SectionTextPanel *t : texts_) t->apply();
}

QList<SectionsTab::Part> SegmentsParts() {
  return {{"data", QObject::tr("Segments From Data"),
           QObject::tr("The data sets (<data>): one row per data file; shared normalizations and energy scales.  Check "
                       "\"Use\" to include a data set, \"Fit?\" to vary its normalization; double-click a cell to "
                       "change it in place, or edit the record in the form below (Prior... for a Gaussian prior).")},
          {"extrapolations", QObject::tr("Segments Without Data"),
           QObject::tr("The extrapolations (<extrapolations>): a quantity of a reaction on a grid of energies and "
                       "angles, calculated without data (the command extrapolate).")},
          {"modelData", QObject::tr("Model Data"),
           QObject::tr("Model data (<modelData>): measured numbers of the model itself (a resonance strength, a width, "
                       "an ANC, a scattering length) compared in the chi-squared.")}};
}

QList<SectionsTab::Part> EffectsParts() {
  // one list per kind of <targets>, each titled and described by the schema
  return {{"targets", QString(), QObject::tr("A data set or an extrapolation names a target, a beam or a detector by "
                                             "its name: \"Used by\" lists them, \"Assign to...\" chooses them."),
           false, true}};
}

QList<SectionsTab::Part> PotentialParts() {
  return {{"potentials", QObject::tr("Nuclear Potential"),
           QObject::tr("Hybrid-potential channels (<potentials>): the nuclear potential of a pair whose channels are "
                       "computed with it instead of hard-sphere functions.")}};
}

