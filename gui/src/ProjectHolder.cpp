#include "ProjectHolder.h"

ProjectHolder::ProjectHolder(QObject *parent)
    : QObject(parent), project_(azure2::front::Project::Empty("untitled.azr")), kinds_(azure2::front::RecordKinds()),
      sections_(azure2::front::RecordSections()) {}

void ProjectHolder::setProject(azure2::front::Project &&project) {
  project_ = std::move(project);
  emit opened();
}

const azure2::front::KindInfo *ProjectHolder::kind(const std::string &name) const {
  for (const azure2::front::KindInfo &k : kinds_)
    if (k.kind == name) return &k;
  return nullptr;
}

QStringList ProjectHolder::kindsOfSection(const QString &section) const {
  QStringList out;
  for (const auto &s : sections_)
    if (QString::fromStdString(s.first) == section)
      for (const std::string &k : s.second) out << QString::fromStdString(k);
  return out;
}

std::vector<azure2::front::Record> ProjectHolder::records(const QString &section) const {
  return project_.Records(section.toStdString());
}

bool ProjectHolder::setRecord(const QString &section, const azure2::front::Record &record, QString &why) {
  std::string w;
  if (!project_.SetRecord(section.toStdString(), record, w)) {
    why = QString::fromStdString(w);
    return false;
  }
  emit sectionEdited(section);
  return true;
}

bool ProjectHolder::removeRecord(const QString &section, const QString &name, QString &why) {
  std::string w;
  if (!project_.RemoveRecord(section.toStdString(), name.toStdString(), w)) {
    why = QString::fromStdString(w);
    return false;
  }
  emit sectionEdited(section);
  return true;
}

bool ProjectHolder::renameRecord(const QString &section, const QString &from, const QString &to, QString &why) {
  std::string w;
  if (!project_.RenameRecord(section.toStdString(), from.toStdString(), to.toStdString(), w)) {
    why = QString::fromStdString(w);
    return false;
  }
  emit sectionEdited(section);
  return true;
}

QString ProjectHolder::configValue(const QString &name, const QString &absent) const {
  return QString::fromStdString(project_.ConfigValue(name.toStdString(), absent.toStdString()));
}

bool ProjectHolder::setConfig(const QString &name, const QString &value, QString &why, bool remove) {
  std::string w;
  if (!project_.SetConfig(name.toStdString(), value.toStdString(), w, remove)) {
    why = QString::fromStdString(w);
    return false;
  }
  emit sectionEdited("config");
  return true;
}

QStringList ProjectHolder::sectionBody(const QString &name) const {
  QStringList out;
  for (const std::string &l : project_.SectionBody(name.toStdString())) out << QString::fromStdString(l);
  return out;
}

void ProjectHolder::setSectionBody(const QString &name, const QStringList &lines) {
  std::vector<std::string> body;
  for (const QString &l : lines) body.push_back(l.toStdString());
  project_.SetSectionBody(name.toStdString(), body);
  emit sectionEdited(name);
}
