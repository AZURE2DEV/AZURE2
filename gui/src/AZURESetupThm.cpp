// The THM part of the main window: the opt-in, the <thm> block, which the
// GUI keeps verbatim unless the THM workspace changes it, and what follows a
// renumbering of the segment lines in it.

#include <QAction>
#include <QFile>
#include <QMessageBox>
#include <QTextStream>

#include "AZURESetup.h"
#include "FittingTab.h"
#include "ThmWorkspace.h"

void AZURESetup::setThmEnabled(bool on) {
  thmEnabled_ = on;
  thmWorkspaceAction->setVisible(on);
  segmentsTab->setThmEnabled(on);
  fittingTab->setThmEnabled(on);
}

QString AZURESetup::thmContent() const {
  QStringList parts;
  if (hasThmBlock) parts << tr("a <thm> block");
  int segments = 0;
  for (const SegmentsDataData &s : segmentsTab->getSegmentsDataModel()->getLines()) segments += s.isTHM ? 1 : 0;
  for (const SegmentsTestData &s : segmentsTab->getSegmentsTestModel()->getLines()) segments += s.isTHM ? 1 : 0;
  if (segments) parts << tr("%n THM segment(s)", "", segments);
  bool binding = false, amplitude = false;
  for (const PairsData &p : pairsTab->getPairsModel()->getPairs()) binding = binding || p.bindingEnergy != 0.0;
  for (const ChannelsData &c : levelsTab->getChannelsModel()->getChannels())
    amplitude = amplitude || (c.radType == QChar('P') && c.gammaIsRWA != 0);
  if (binding) parts << tr("THM binding energies");
  if (amplitude) parts << tr("widths entered as reduced width amplitudes");
  return parts.join(", ");
}

void AZURESetup::applyThmOption(bool on) {
  if (!on && thmEnabled_) {
    const QString content = thmContent();
    if (!content.isEmpty()) {
      // Nothing is deleted: the content is saved as read and AZURE2 uses it.
      QMessageBox box(QMessageBox::Warning, tr("Trojan Horse Method"),
                      tr("This project has %1.\nThey stay in the project and AZURE2 still uses them; "
                         "only their controls are hidden.")
                          .arg(content),
                      QMessageBox::Ok | QMessageBox::Cancel, this);
      box.setTextFormat(Qt::PlainText);
      box.setDefaultButton(QMessageBox::Cancel);
      if (box.exec() != QMessageBox::Ok) return;
    }
  }
  setThmEnabled(on);
}

bool AZURESetup::readThmContent(const QString &filename) {
  // The <thm> block may sit anywhere (the engine searches the whole file), so
  // it gets a pass of its own from the top.
  QFile file(filename);
  if (!file.open(QIODevice::ReadOnly)) return false;
  QTextStream in(&file);
  if (!readThmBlock(in, thmBlockLines, hasThmBlock)) return false;
  // The rest of the project is read: THM is on if any of it is THM.
  setThmEnabled(!thmContent().isEmpty());
  return true;
}

bool AZURESetup::readThmBlock(QTextStream &in, QStringList &lines, bool &present) {
  lines.clear();
  present = false;
  QString line;
  // Same tests as Config::ReadThmBlock: the block opens at the first line whose
  // first non-blank characters are <thm>, and closes at a line that is </thm>
  // once a # comment and surrounding blanks are removed.
  while (!in.atEnd()) {
    line = in.readLine();
    if (line.trimmed().startsWith(QString("<thm>"))) {
      present = true;
      break;
    }
  }
  if (!present) return true;
  while (!in.atEnd()) {
    line = in.readLine();
    QString code = line;
    int hash = code.indexOf('#');
    if (hash >= 0) code.truncate(hash);
    if (code.trimmed() == QString("</thm>")) return true;
    lines << line;
  }
  lines.clear();
  present = false;
  return false;
}

void AZURESetup::writeThmBlock(QTextStream &out, const QStringList &lines) {
  out << "<thm>" << Qt::endl;
  for (const QString &line : lines) out << line << Qt::endl;
  out << "</thm>" << Qt::endl;
}

bool AZURESetup::thmSettings(ThmSettings &settings, QString *error) const {
  if (!hasThmBlock) {
    settings = ThmSettings();
    return true;
  }
  return ThmSettings::parse(thmBlockLines, settings, error);
}

void AZURESetup::setThmSettings(const ThmSettings &settings) {
  ThmSettings current;
  if (thmSettings(current) && current == settings) return;  // untouched: keep the block verbatim
  if (settings.isDefault()) {
    hasThmBlock = false;
    thmBlockLines.clear();
    return;
  }
  thmBlockLines = settings.compose(hasThmBlock ? thmBlockLines : QStringList());
  hasThmBlock = true;
}

QStringList AZURESetup::followThmSegments(const QVector<int> &newNumber, bool test) {
  QStringList notes;
  ThmSettings thm;
  if (hasThmBlock && thmSettings(thm)) {
    const QStringList removed = thm.renumberSegments(newNumber, test);
    setThmSettings(thm);
    if (!removed.isEmpty() && !test)
      notes << tr("THM experiment %1 removed: no segments left.").arg(removed.join(", "));
  }
  return notes;
}

void AZURESetup::followTestSegments(const QVector<int> &newNumber) { followThmSegments(newNumber, true); }

void AZURESetup::editThmWorkspace() {
  if (!thmEnabled_) return;
  ThmSettings current;
  QString error;
  if (!thmSettings(current, &error)) {
    QMessageBox::warning(this, tr("THM Workspace"),
                         tr("The <thm> block of this project has a line the engine would refuse:\n%1\n"
                            "Correct it in the .azr file; the block is kept as it is.")
                             .arg(error));
    return;
  }
  ThmWorkspace workspace(this, current, this);
  workspace.exec();  // Accept installs the pages' values (ThmWorkspace::apply)
}
