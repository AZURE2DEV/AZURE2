// The THM part of the Fitting tab: the coherent-background parameters
// (cbackground= of the <thm> experiments), on a tab of their own that is
// shown only while THM is on (AZURESetup::setThmEnabled) and there are some.

#include <QSet>
#include <QTabWidget>
#include <QTableWidget>
#include <algorithm>
#include <cmath>
#include <map>

#include "AZURESetup.h"
#include "FittingTab.h"
#include "LevelsTab.h"
#include "SegmentsTab.h"
#include "ThmExperiment.h"
#include "ThmSettings.h"

void FittingTab::setThmEnabled(bool on) {
  thmEnabled_ = on;
  showCoherentTab();
}

void FittingTab::showCoherentTab() {
  const int cbkgTab = paramTabWidget->indexOf(cbkgParamsTable);
  const bool shown = thmEnabled_ && cbkgParamsTable->rowCount() > 0;
  // The names (cbkg_<experiment>_<J^pi>_<exit>_<s,l,s',l'>_re0 ...) are long: the column fits them.
  if (shown) {
    const int was = cbkgParamsTable->columnWidth(0);
    cbkgParamsTable->resizeColumnToContents(0);
    cbkgParamsTable->setColumnWidth(0, std::max(was, cbkgParamsTable->columnWidth(0)));
  }
  if (shown && cbkgTab < 0)
    paramTabWidget->addTab(cbkgParamsTable, "THM Background");
  else if (!shown && cbkgTab >= 0)
    paramTabWidget->removeTab(cbkgTab);
}

QSet<QString> FittingTab::coherentParameterNames() const {
  QSet<QString> names;
  AZURESetup *s = setup();
  ThmSettings thm;
  if (s && s->thmSettings(thm))
    for (const ThmExperimentRecord &r : ThmExperimentRecord::read(thm.experimentLines))
      if (!r.cbackground.isEmpty())
        for (const QString &n : coherentNames(r)) names.insert(n);
  return names;
}

AZURESetup *FittingTab::setup() const {
  for (QWidget *p = parentWidget(); p; p = p->parentWidget())
    if (AZURESetup *s = qobject_cast<AZURESetup *>(p)) return s;
  return nullptr;
}

QStringList FittingTab::coherentNames(const ThmExperimentRecord &record, QList<double> *values,
                                      QList<bool> *fixed) const {
  QStringList names;
  std::vector<ThmExperiment::CoherentTerm> terms;
  if (!levelsTab_ || !segmentsTab_ ||
      !ParseThmCoherentBackground(record.cbackground.toStdString(), terms).empty())
    return names;
  const QList<LevelsData> levels = levelsTab_->getLevelsModel()->getLevels();
  const QList<ChannelsData> channels = levelsTab_->getChannelsModel()->getChannels();
  const QList<SegmentsDataData> segments = segmentsTab_->getSegmentsDataModel()->getLines();
  int entrance = 0;
  for (int k : record.segments)
    if (k >= 1 && k <= segments.size() && segments.at(k - 1).isActive) {
      entrance = segments.at(k - 1).entrancePairIndex;
      break;
    }
  // The channels of a J^pi group as CNuc::Fill makes them: the levels in
  // AZURE2's order, each level's channels in file order, (pair, s, l) once.
  struct Channel {
    int pair;
    double s;
    int l;
  };
  const QList<int> order = const_cast<FittingTab *>(this)->engineLevelOrder(levelsTab_->writeOrder());
  static const char *parts[4] = {"re0", "im0", "re1", "im1"};
  for (const ThmExperiment::CoherentTerm &t : terms) {
    QList<Channel> group;
    for (int la : order) {
      const LevelsData &lv = levels.at(la);
      if (std::fabs(lv.jValue - t.J) > 1e-6 || lv.piValue != t.parity) continue;
      for (int ch = 0; ch < channels.size(); ch++) {
        const ChannelsData &c = channels.at(ch);
        if (c.levelIndex != la || c.radType != QChar('P')) continue;
        bool seen = false;
        for (const Channel &g : group) seen = seen || (g.pair == c.pairIndex + 1 && g.s == c.sValue && g.l == c.lValue);
        if (!seen) group.append(Channel{c.pairIndex + 1, c.sValue, c.lValue});
      }
    }
    for (const Channel &in : group) {
      if (in.pair != entrance) continue;
      if (t.hasChannels && (std::fabs(in.s - t.s) > 1e-6 || in.l != t.l)) continue;
      for (const Channel &out : group) {
        if (out.pair != t.exitKey) continue;
        if (t.hasChannels && (std::fabs(out.s - t.sp) > 1e-6 || out.l != t.lp)) continue;
        const QString stem = QString::fromStdString(
            ThmCoherentParamStem(record.name.toStdString(), t.J, t.parity, t.exitKey, in.s, in.l, out.s, out.l));
        for (int k = 0; k < 2 * t.form; k++) {
          names << stem + parts[k];
          if (values) values->append(t.hasValues ? t.value[k] : 0.0);
          if (fixed) fixed->append(t.hasValues && t.fixed[k]);
        }
      }
    }
  }
  return names;
}

void FittingTab::appendCoherentParameters() {
  AZURESetup *s = setup();
  ThmSettings thm;
  if (!s || !s->thmSettings(thm)) return;
  for (const ThmExperimentRecord &r : ThmExperimentRecord::read(thm.experimentLines)) {
    if (r.cbackground.isEmpty()) continue;
    QList<double> values;
    QList<bool> fixed;
    const QStringList names = coherentNames(r, &values, &fixed);
    for (int k = 0; k < names.size(); k++) {
      if (fixed.at(k)) continue;  // shown, like a level's, only when free
      FittingParameter p;
      p.name = names.at(k);
      p.value = values.at(k);
      p.lowerLimit = 0;
      p.upperLimit = 0;
      p.error = values.at(k) != 0.0 ? 0.1 * std::fabs(values.at(k)) : 0.1;
      p.fitError = 0.0;
      p.useAsNuisance = false;
      p.category = "cbkg";
      p.minuitIndex = -1;
      p.levelIndex = -1;
      p.channelIndex = -1;
      fittingParameters.append(p);
    }
  }
}

bool FittingTab::applyCoherentValues(const QMap<QString, double> &values) {
  AZURESetup *s = setup();
  ThmSettings thm;
  if (!s || !s->thmSettings(thm)) return false;
  const QList<ThmExperimentRecord> before = ThmExperimentRecord::read(thm.experimentLines);
  QList<ThmExperimentRecord> after = before;
  std::map<std::string, double> named;
  for (auto it = values.begin(); it != values.end(); ++it) named[it.key().toStdString()] = it.value();
  bool changed = false;
  for (ThmExperimentRecord &r : after) {
    if (r.cbackground.isEmpty()) continue;
    std::vector<std::string> names;
    for (const QString &n : coherentNames(r)) names.push_back(n.toStdString());
    bool mine = false;
    for (const std::string &n : names) mine = mine || named.count(n);
    if (!mine) continue;
    std::string out;
    if (!ApplyThmCoherentValues(r.name.toStdString(), r.cbackground.toStdString(), names, named, out).empty())
      continue;
    if (QString::fromStdString(out) != r.cbackground) {
      r.cbackground = QString::fromStdString(out);
      changed = true;
    }
  }
  if (!changed) return false;
  thm.experimentLines = ThmExperimentRecord::compose(thm.experimentLines, before, after);
  s->setThmSettings(thm);
  return true;
}
