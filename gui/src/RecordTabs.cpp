#include "RecordTabs.h"

#include <QTabWidget>
#include <QVBoxLayout>

#include "RecordsPanel.h"
#include "SectionTextPanel.h"

SectionsTab::SectionsTab(ProjectHolder *holder, const QList<Part> &parts, QWidget *parent) : QWidget(parent) {
  QVBoxLayout *layout = new QVBoxLayout(this);
  if (parts.size() > 1) {
    tabs_ = new QTabWidget;
    layout->addWidget(tabs_);
  }
  for (const Part &p : parts) {
    QWidget *panel = nullptr;
    if (p.text) {
      SectionTextPanel *t = new SectionTextPanel(holder, p.section, p.help);
      texts_ << t;
      panel = t;
    } else {
      RecordsPanel *r = new RecordsPanel(holder, p.section, p.help);
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

SectionTextPanel *SectionsTab::text(const QString &section) const {
  for (SectionTextPanel *p : texts_)
    if (p->objectName() == "section_" + section) return p;
  return nullptr;
}

bool SectionsTab::show(const QString &section, const QString &record) {
  QWidget *panel = records(section);
  if (!panel) panel = text(section);
  if (!panel) return false;
  if (tabs_) tabs_->setCurrentWidget(panel);
  if (RecordsPanel *r = records(section)) r->selectRecord(record);
  return true;
}

void SectionsTab::applyText() {
  for (SectionTextPanel *t : texts_) t->apply();
}

QList<SectionsTab::Part> SegmentsParts() {
  return {{"data", QObject::tr("Segments From Data"),
           QObject::tr("The data sets (<data>): one record per data file, with its reaction, quantity, unit, columns, "
                       "windows and normalization; shared normalizations and energy scales.  Spin observables, line "
                       "windows, priors (\"norm 1 fit prior 1 0.05\") and nuisance parameters are fields of the record.")},
          {"extrapolations", QObject::tr("Segments Without Data"),
           QObject::tr("The extrapolations (<extrapolations>): a quantity of a reaction on a grid of energies and "
                       "angles, calculated without data (the command extrapolate).")},
          {"modelData", QObject::tr("Model Data"),
           QObject::tr("Model data (<modelData>): measured numbers of the model itself (a resonance strength, a width, "
                       "an ANC, a scattering length) compared in the chi-squared.")}};
}

QList<SectionsTab::Part> EffectsParts() {
  return {{"targets", QObject::tr("Targets, Beams and Detectors"),
           QObject::tr("The apparatus (<targets>): targets (solid, implanted, gas, layers, isotopes and their "
                       "abundances, stopping and straggling, fitted target parameters), beams (spread, profile) and "
                       "detectors (resolution kernels).  A data set names them by name.")}};
}

QList<SectionsTab::Part> PotentialParts() {
  return {{"potentials", QObject::tr("Nuclear Potential"),
           QObject::tr("Hybrid-potential channels (<potentials>): the nuclear potential of a pair whose channels are "
                       "computed with it instead of hard-sphere functions.")}};
}

