#include "ThmSettings.h"

#include <QDir>
#include <QLocale>
#include <QSet>
#include <sstream>
#include <vector>

#include "Config.h"
#include "ThmExperiment.h"

namespace {

// Shortest text that reads back as the same double.
QString numberText(double x) { return QString::number(x, 'g', QLocale::FloatingPointShortest); }

// The engine reads numbers with operator>> on an istringstream; so does this,
// so that "0.5 " or "5" are taken exactly as AZURE2 takes them.
bool readDouble(const QString &text, double &x) {
  std::istringstream s(text.toStdString());
  return !!(s >> x);
}
bool readInt(const QString &text, int &x, bool strict) {
  std::istringstream s(text.toStdString());
  if (!(s >> x)) return false;
  std::string rest;
  return !strict || !(s >> rest);
}

// The code part of a line: comment stripped, trimmed (as the engine reads it).
QString codeOf(const QString &line) {
  QString code = line;
  int hash = code.indexOf('#');
  if (hash >= 0) code.truncate(hash);
  return code.trimmed();
}

}  // namespace

// ---------------------------------------------------------------------------
// ThmSettings

bool ThmSettings::operator==(const ThmSettings &o) const {
  return entranceL == o.entranceL && vertex == o.vertex && kinematics == o.kinematics &&
         coulombIntegral == o.coulombIntegral && spectatorEnergy == o.spectatorEnergy &&
         spectatorByPair == o.spectatorByPair && weight == o.weight && weightTest == o.weightTest &&
         experimentLines == o.experimentLines;
}

QList<QPair<QString, QString>> ThmSettings::keyValues() const {
  const ThmSettings d;
  QList<QPair<QString, QString>> kv;
  if (entranceL != d.entranceL) kv << qMakePair(QString("entranceL"), entranceL);
  if (vertex != d.vertex) kv << qMakePair(QString("vertex"), vertex);
  if (kinematics != d.kinematics) kv << qMakePair(QString("kinematics"), kinematics);
  if (coulombIntegral != d.coulombIntegral) kv << qMakePair(QString("coulombIntegral"), QString("1"));
  if (spectatorEnergy != d.spectatorEnergy)
    kv << qMakePair(QString("spectatorEnergy"), numberText(spectatorEnergy));
  for (auto it = spectatorByPair.begin(); it != spectatorByPair.end(); ++it)
    kv << qMakePair(QString("spectatorEnergy[%1]").arg(it.key()), numberText(it.value()));
  for (auto it = weight.begin(); it != weight.end(); ++it)
    kv << qMakePair(QString("weight[%1]").arg(it.key()), it.value());
  for (auto it = weightTest.begin(); it != weightTest.end(); ++it)
    kv << qMakePair(QString("weightTest[%1]").arg(it.key()), it.value());
  return kv;
}

QString ThmSettings::experimentName(const QString &rawLine) {
  QString code = codeOf(rawLine);
  if (!code.startsWith("experiment[")) return QString();
  int close = code.indexOf(']');
  return close < 0 ? QString() : code.mid(11, close - 11);
}

bool ThmSettings::parseLine(const QString &rawLine, QString &key, QString &value, ThmSettings &s) {
  key.clear();
  value.clear();
  QString line = codeOf(rawLine);
  if (line.isEmpty()) return true;
  if (line.startsWith("experiment[")) {
    // A THM experiment record: kept raw (an empty key, as a comment); parse()
    // runs the engine's parser on the set.
    s.experimentLines << rawLine;
    return true;
  }
  int eq = line.indexOf('=');
  if (eq < 0) return false;
  QString k = line.left(eq);
  while (!k.isEmpty() && (k.endsWith(' ') || k.endsWith('\t'))) k.chop(1);
  QString v = line.mid(eq + 1);
  while (!v.isEmpty() && (v.startsWith(' ') || v.startsWith('\t'))) v.remove(0, 1);

  if (k == "vertex") {
    if (v == "real") v = "perlevel";
    if (v != "onshell" && v != "constant" && v != "perlevel") return false;
    s.vertex = v;
  } else if (k == "kinematics") {
    if (v != "lacognata" && v != "triple" && v != "kf3body" && v != "lambda32") return false;
    s.kinematics = v;
  } else if (k == "entranceL") {
    if (v != "coherent" && v != "incoherent") return false;
    s.entranceL = v;
  } else if (k == "coulombIntegral") {
    if (v == "1" || v == "true" || v == "on")
      s.coulombIntegral = true;
    else if (v == "0" || v == "false" || v == "off")
      s.coulombIntegral = false;
    else
      return false;
    v = s.coulombIntegral ? "1" : "0";
  } else if (k.startsWith("spectatorEnergy")) {
    double x;
    if (!readDouble(v, x) || !(x >= 0.0)) return false;
    v = numberText(x);
    if (k == "spectatorEnergy") {
      s.spectatorEnergy = x;
    } else if (k.size() > 17 && k[15] == '[' && k.endsWith(']')) {
      int pair;
      if (!readInt(k.mid(16, k.size() - 17), pair, false)) return false;
      s.spectatorByPair[pair] = x;
      k = QString("spectatorEnergy[%1]").arg(pair);
    } else {
      return false;
    }
  } else if (k.startsWith("weight")) {
    bool test = k.startsWith("weightTest[");
    int open = test ? 10 : 6;
    int segment;
    if (!(k.size() > open + 2 && k[open] == '[' && k.endsWith(']') && !v.isEmpty())) return false;
    if (!readInt(k.mid(open + 1, k.size() - open - 2), segment, true) || segment < 1) return false;
    (test ? s.weightTest : s.weight)[segment] = v;
    k = QString(test ? "weightTest[%1]" : "weight[%1]").arg(segment);
  } else {
    return false;
  }
  key = k;
  value = v;
  return true;
}

QString ThmSettings::checkExperimentLines(const QStringList &lines) {
  std::vector<ThmExperiment> experiments;
  for (const QString &line : lines) {
    QString code = codeOf(line);
    if (code.isEmpty()) continue;
    std::string why = ParseThmExperimentLine(code.toStdString(), experiments);
    if (!why.empty()) return "<thm> " + QString::fromStdString(why);
  }
  std::string why = CheckThmExperiments(experiments);
  return why.empty() ? QString() : "<thm> " + QString::fromStdString(why);
}

bool ThmSettings::parse(const QStringList &lines, ThmSettings &out, QString *error) {
  ThmSettings s;
  for (const QString &line : lines) {
    QString key, value;
    if (!parseLine(line, key, value, s)) {
      if (error) *error = QString("<thm> line not understood: '%1'").arg(line.trimmed());
      return false;
    }
  }
  QString why = checkExperimentLines(s.experimentLines);
  if (!why.isEmpty()) {
    if (error) *error = why;
    return false;
  }
  out = s;
  return true;
}

QString ThmSettings::validate(const QString &projectDir) const {
  const QMap<int, QString> *maps[2] = {&weight, &weightTest};
  const char *names[2] = {"weight", "weightTest"};
  for (int m = 0; m < 2; m++) {
    for (auto it = maps[m]->begin(); it != maps[m]->end(); ++it) {
      QString name = it.value().trimmed();
      if (name.isEmpty()) return QString("%1[%2]: no file given.").arg(names[m]).arg(it.key());
      if (name.contains('#'))
        return QString("%1[%2]: the path contains '#', which starts a comment in the .azr.")
            .arg(names[m])
            .arg(it.key());
      QString path = QDir::isAbsolutePath(name) || projectDir.isEmpty()
                         ? name
                         : QDir(projectDir).filePath(name);
      ThmWeightTable table;
      std::string why = table.Read(path.toStdString());
      if (!why.empty())
        return QString("%1[%2]: %3").arg(names[m]).arg(it.key()).arg(QString::fromStdString(why));
    }
  }
  for (auto it = spectatorByPair.begin(); it != spectatorByPair.end(); ++it)
    if (!(it.value() >= 0.0)) return QString("spectatorEnergy[%1] must be >= 0.").arg(it.key());
  if (!(spectatorEnergy >= 0.0)) return QString("spectatorEnergy must be >= 0.");
  return QString();
}

QStringList ThmSettings::compose(const QStringList &oldLines) const {
  const QList<QPair<QString, QString>> kv = keyValues();
  QMap<QString, QString> wanted;
  for (const auto &p : kv) wanted[p.first] = p.second;
  QSet<QString> written;
  QList<bool> used;  // experimentLines already placed
  for (int i = 0; i < experimentLines.size(); i++) used << false;
  QStringList out;
  for (const QString &line : oldLines) {
    const QString name = experimentName(line);
    if (!name.isEmpty()) {
      int same = -1;
      for (int i = 0; i < experimentLines.size() && same < 0; i++)
        if (!used[i] && experimentLines[i] == line) same = i;
      if (same >= 0) {  // unchanged: in place
        used[same] = true;
        out << line;
        continue;
      }
      // Changed or gone: the unplaced new lines of this experiment go here.
      for (int i = 0; i < experimentLines.size(); i++)
        if (!used[i] && experimentName(experimentLines[i]) == name) {
          used[i] = true;
          out << experimentLines[i];
        }
      continue;
    }
    ThmSettings scratch;
    QString key, value;
    if (!parseLine(line, key, value, scratch) || key.isEmpty()) {
      out << line;  // comment, blank (or a line the engine would refuse: left alone)
      continue;
    }
    if (!wanted.contains(key) || written.contains(key)) continue;  // back to default, or repeated
    written.insert(key);
    if (value == wanted[key]) {
      out << line;
      continue;
    }
    // Rewrite the value, keeping the indentation and any inline comment.
    int hash = line.indexOf('#');
    QString code = hash >= 0 ? line.left(hash) : line;
    QString comment = hash >= 0 ? line.mid(hash) : QString();
    int start = 0;
    while (start < code.size() && code[start].isSpace()) start++;
    int end = code.size();
    while (end > start && code[end - 1].isSpace()) end--;
    out << code.left(start) + key + "=" + wanted[key] + code.mid(end) + comment;
  }
  for (const auto &p : kv)
    if (!written.contains(p.first)) out << p.first + "=" + p.second;
  for (int i = 0; i < experimentLines.size(); i++)
    if (!used[i]) out << experimentLines[i];
  return out;
}

// ---------------------------------------------------------------------------
// ThmExperimentRecord

bool ThmExperimentRecord::sameAs(const ThmExperimentRecord &o) const {
  return name == o.name && segments == o.segments && background == o.background &&
         hasBackgroundKey == o.hasBackgroundKey && beam == o.beam && target == o.target &&
         spectator == o.spectator && beamEnergy == o.beamEnergy && lineshape == o.lineshape && ps == o.ps &&
         psNodes == o.psNodes && distortion == o.distortion && opticalAA == o.opticalAA && opticalSF == o.opticalSF &&
         spectatorAngle == o.spectatorAngle && distortionRef == o.distortionRef &&
         distortionRatio == o.distortionRatio && boundState == o.boundState && extraTokens == o.extraTokens;
}

QString ThmExperimentRecord::segmentsListText(const QList<int> &segments) {
  QStringList parts;
  for (int i = 0; i < segments.size();) {
    int j = i;
    while (j + 1 < segments.size() && segments[j + 1] == segments[j] + 1) j++;
    if (j >= i + 2)
      parts << QString("%1-%2").arg(segments[i]).arg(segments[j]);
    else
      for (int k = i; k <= j; k++) parts << QString::number(segments[k]);
    i = j + 1;
  }
  return parts.join(",");
}

bool ThmExperimentRecord::expandSegments(const QString &text, QList<int> &out) {
  out.clear();
  std::vector<ThmExperiment> x;
  if (!ParseThmExperimentLine("experiment[x] segments=" + text.toStdString(), x).empty() || x.empty())
    return false;
  for (int k : x[0].segments) out << k;
  return true;
}

QString ThmExperimentRecord::line() const {
  QStringList tokens;
  if (!segments.isEmpty())
    tokens << "segments=" + (segmentsText.isEmpty() ? segmentsListText(segments) : segmentsText);
  if (hasBackgroundKey || background != "none") tokens << "background=" + background;
  if (!beam.isEmpty()) tokens << "beam=" + beam;
  if (!target.isEmpty()) tokens << "target=" + target;
  if (!spectator.isEmpty()) tokens << "spectator=" + spectator;
  if (!beamEnergy.isEmpty()) tokens << "Ebeam=" + beamEnergy;
  if (lineshape) tokens << "lineshape=on";
  if (!ps.isEmpty()) tokens << "ps=" + ps;
  if (!psNodes.isEmpty()) tokens << "psNodes=" + psNodes;
  const QString distortionKeys[7][2] = {{"distortion", distortion},         {"opticalAA", opticalAA},
                                        {"opticalSF", opticalSF},           {"spectatorAngle", spectatorAngle},
                                        {"distortionRef", distortionRef},   {"distortionRatio", distortionRatio},
                                        {"boundState", boundState}};
  for (const auto &kv : distortionKeys)
    if (!kv[1].isEmpty()) tokens << kv[0] + "=" + kv[1];
  tokens << extraTokens;
  return QString("experiment[%1] %2").arg(name, tokens.join(' '));
}

QList<ThmExperimentRecord> ThmExperimentRecord::read(const QStringList &lines) {
  QList<ThmExperimentRecord> records;
  for (const QString &raw : lines) {
    const QString name = ThmSettings::experimentName(raw);
    if (name.isEmpty()) continue;
    int at = -1;
    for (int i = 0; i < records.size(); i++)
      if (records[i].name == name) at = i;
    if (at < 0) {
      ThmExperimentRecord r;
      r.name = r.originName = name;
      records << r;
      at = records.size() - 1;
    }
    ThmExperimentRecord &r = records[at];
    QString code = codeOf(raw);
    const QStringList tokens = code.mid(code.indexOf(']') + 1).split(QRegExp("[ \t]+"), Qt::SkipEmptyParts);
    for (const QString &token : tokens) {
      int eq = token.indexOf('=');
      QString key = token.left(eq), value = token.mid(eq + 1);
      if (key == "segments") {
        r.segmentsText = value;
        expandSegments(value, r.segments);
      } else if (key == "background") {
        r.background = value;
        r.hasBackgroundKey = true;
      } else if (key == "beam") {
        r.beam = value;
      } else if (key == "target") {
        r.target = value;
      } else if (key == "spectator") {
        r.spectator = value;
      } else if (key == "Ebeam") {
        r.beamEnergy = value;
      } else if (key == "lineshape") {
        r.lineshape = value == "on";  // on | off, checked by the engine's parser
      } else if (key == "ps") {
        r.ps = value;
      } else if (key == "psNodes") {
        r.psNodes = value;
      } else if (key == "distortion") {
        r.distortion = value;
      } else if (key == "opticalAA") {
        r.opticalAA = value;
      } else if (key == "opticalSF") {
        r.opticalSF = value;
      } else if (key == "spectatorAngle") {
        r.spectatorAngle = value;
      } else if (key == "distortionRef") {
        r.distortionRef = value;
      } else if (key == "distortionRatio") {
        r.distortionRatio = value;
      } else if (key == "boundState") {
        r.boundState = value;
      } else {
        r.extraTokens << token;
      }
    }
  }
  return records;
}

QStringList ThmExperimentRecord::compose(const QStringList &oldLines, const QList<ThmExperimentRecord> &oldRecords,
                                         const QList<ThmExperimentRecord> &records) {
  auto current = [&](const QString &origin) -> const ThmExperimentRecord * {
    for (const ThmExperimentRecord &r : records)
      if (!r.originName.isEmpty() && r.originName == origin) return &r;
    return nullptr;
  };
  auto original = [&](const QString &origin) -> const ThmExperimentRecord * {
    for (const ThmExperimentRecord &r : oldRecords)
      if (r.originName == origin) return &r;
    return nullptr;
  };
  QStringList out;
  QSet<QString> written;
  for (const QString &raw : oldLines) {
    const QString name = ThmSettings::experimentName(raw);
    if (name.isEmpty()) continue;
    const ThmExperimentRecord *r = current(name);
    if (!r) continue;  // removed
    const ThmExperimentRecord *o = original(name);
    if (o && r->sameAs(*o)) {
      out << raw;
    } else if (!written.contains(name)) {
      written.insert(name);
      out << r->line();
    }
  }
  for (const ThmExperimentRecord &r : records)
    if (r.originName.isEmpty()) out << r.line();
  return out;
}
