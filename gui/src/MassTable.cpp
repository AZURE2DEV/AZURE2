#include "MassTable.h"

#include <QFile>
#include <QRegularExpression>
#include <QTextStream>
#include <algorithm>
#include <cmath>

#include "PairsModel.h"

// The resource lives in the AZURESetup library; a program that links the library without running the GUI's
// main() (a test) has to register it itself.  Q_INIT_RESOURCE must be used outside any namespace.
static void initMassResource() { Q_INIT_RESOURCE(AZURESetup); }

namespace {

const char *kResource = ":/ame/mass_1.mas20";

// A fixed-width AME field; '#' marks an estimated value in place of the decimal point.  Empty or '*' -> false.
bool amedouble(const QString &line, int start, int width, double *value, bool *estimated) {
  QString t = line.mid(start, width).trimmed();
  if (t.isEmpty() || t == "*") return false;
  if (t.contains('#')) {
    if (estimated) *estimated = true;
    t.replace('#', '.');
  }
  bool ok = false;
  *value = t.toDouble(&ok);
  return ok;
}

}  // namespace

const MassTable &MassTable::instance() {
  static const MassTable table;
  return table;
}

MassTable::MassTable() {
  if (!QFile::exists(kResource)) initMassResource();
  QFile file(kResource);
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
  QTextStream in(&file);
  // format a1,i3,i5,i5,i5,1x,a3,a4,1x,f14.6,f12.6,f13.5,1x,f10.5,1x,a2,f13.5,f11.5,1x,i3,1x,f13.6,f12.6
  while (!in.atEnd()) {
    QString line = in.readLine();
    if (line.size() < 123) continue;
    bool okZ = false, okA = false, okW = false;
    int Z = line.mid(9, 5).trimmed().toInt(&okZ);
    int A = line.mid(14, 5).trimmed().toInt(&okA);
    QString symbol = line.mid(20, 3).trimmed();
    if (!okZ || !okA || symbol.isEmpty()) continue;
    bool alpha = true;
    for (QChar c : symbol) alpha = alpha && c.isLetter();
    if (!alpha) continue;
    MassEntry e;
    e.Z = Z;
    e.A = A;
    e.symbol = Z == 0 ? QString("n") : symbol.left(1).toUpper() + symbol.mid(1).toLower();
    if (!amedouble(line, 28, 14, &e.massExcess, &e.estimated)) continue;
    amedouble(line, 42, 12, &e.massExcessUnc, nullptr);
    int whole = line.mid(106, 3).trimmed().toInt(&okW);
    double micro = 0.0;
    if (!okW || !amedouble(line, 110, 13, &micro, &e.estimated)) continue;
    // the file's 1e-12 u resolution, without float noise
    e.atomicMass = std::round((whole + micro * 1e-6) * 1e12) / 1e12;
    double microUnc = 0.0;
    if (amedouble(line, 123, 12, &microUnc, nullptr)) e.atomicMassUnc = microUnc * 1e-6;
    entries_.insert(Z * 1000 + A, e);
    if (Z > 0) symbolToZ_.insert(e.symbol.toLower(), Z);
  }
}

const MassEntry *MassTable::lookup(int Z, int A) const {
  auto it = entries_.constFind(Z * 1000 + A);
  return it == entries_.constEnd() ? nullptr : &it.value();
}

const MassEntry *MassTable::identify(int Z, double mass) const {
  int A = static_cast<int>(std::lround(mass));
  if (A <= 0 || Z < 0) return nullptr;
  return lookup(Z, A);
}

bool MassTable::parseNuclide(const QString &text, int *Z, int *A) const {
  const QString s = text.trimmed().toLower();
  static const QHash<QString, QPair<int, int>> shorthand = {
      {"n", {0, 1}}, {"neutron", {0, 1}}, {"p", {1, 1}}, {"proton", {1, 1}}, {"d", {1, 2}}, {"deuteron", {1, 2}},
      {"t", {1, 3}}, {"triton", {1, 3}}, {"helion", {2, 3}}, {"a", {2, 4}}, {"alpha", {2, 4}}};
  auto it = shorthand.constFind(s);
  if (it != shorthand.constEnd()) {
    *Z = it.value().first;
    *A = it.value().second;
    return true;
  }
  static const QRegularExpression massFirst("^(\\d+)\\s*-?\\s*([a-z]{1,3})$");
  static const QRegularExpression symbolFirst("^([a-z]{1,3})\\s*-?\\s*(\\d+)$");
  QString sym, mass;
  QRegularExpressionMatch m = massFirst.match(s);
  if (m.hasMatch()) {
    mass = m.captured(1);
    sym = m.captured(2);
  } else {
    m = symbolFirst.match(s);
    if (!m.hasMatch()) return false;
    sym = m.captured(1);
    mass = m.captured(2);
  }
  const int a = mass.toInt();
  if (sym == "n" && a == 1) {   // "1n" / "n1" is the neutron
    *Z = 0;
    *A = 1;
    return true;
  }
  auto z = symbolToZ_.constFind(sym);
  if (z == symbolToZ_.constEnd() || a <= 0) return false;
  *Z = z.value();
  *A = a;
  return true;
}

bool MassTable::separationEnergy(int Z1, int A1, int Z2, int A2, double *S) const {
  const MassEntry *l = lookup(Z1, A1), *h = lookup(Z2, A2), *c = lookup(Z1 + Z2, A1 + A2);
  if (!l || !h || !c) return false;
  // mass excesses carry 1e-6 keV, so S carries 1e-9 MeV
  *S = std::round((l->massExcess + h->massExcess - c->massExcess) * 1e6) / 1e9;
  return true;
}

PairMassCheck MassTable::check(const PairsData &pair) const {
  PairMassCheck r;
  if (pair.pairType == 0) r.light = identify(pair.lightZ, pair.lightM);
  r.heavy = identify(pair.heavyZ, pair.heavyM);
  if (pair.pairType == 0) {
    r.fillable = r.light && r.heavy;
    if (r.fillable && separationEnergy(r.light->Z, r.light->A, r.heavy->Z, r.heavy->A, &r.sTable)) {
      r.hasS = true;
      r.dSkeV = (pair.seperationEnergy - r.sTable) * 1e3;
    }
  } else {
    r.fillable = r.heavy != nullptr;
  }
  if (r.light) r.massOffset = std::max(r.massOffset, std::fabs(pair.lightM - r.light->atomicMass));
  if (r.heavy) r.massOffset = std::max(r.massOffset, std::fabs(pair.heavyM - r.heavy->atomicMass));
  return r;
}
