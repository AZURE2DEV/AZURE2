#ifndef THMTAB_H
#define THMTAB_H

#include <QMap>
#include <QWidget>

#include "ProjectHolder.h"

class SectionTextPanel;

QT_BEGIN_NAMESPACE
class QFormLayout;
class QLabel;
class QListWidget;
QT_END_NAMESPACE

/*!
 * The THM tab (docs/upgrade/qt_gui.md section 8): <thm> as two forms built from the engine's tables of its keys
 * (azure2::front::ThmKinds: the global options, and the keys of an experiment[<name>] line: beam, target, spectator,
 * Ebeam, the momentum distribution, the distortion, the vertex model with vertexModel=ff and its well, volume term and
 * nodes, the background, ...), each change written by the engine's edit of the section's lines
 * (azure2::front::SetThmOption, SetThmExperimentKey: a refused value is named by the engine's reader), and a "text"
 * view of the whole section.  A key whose row says it applies when another is given is enabled only then.
 */
class ThmTab : public QWidget {
  Q_OBJECT

 public:
  explicit ThmTab(ProjectHolder *holder, QWidget *parent = nullptr);
  /// Set a global option or a key of the selected experiment through the form (as an edit of its editor); "" removes.
  bool setOption(const QString &key, const QString &value, QString &why);
  bool setExperimentKey(const QString &experiment, const QString &key, const QString &value, QString &why);
  bool addExperiment(const QString &name, QString &why);
  bool removeExperiment(const QString &name, QString &why);
  bool selectExperiment(const QString &name);
  QStringList experiments() const;
  /// The editor of a key in the experiment form (nullptr when none).
  QWidget *experimentEditor(const QString &key) const { return experimentEditors_.value(key, nullptr); }
  /// Write the text view into the project (before a save).
  void applyText();
  /// Standard or expert view (qt_gui.md section 19.4): in standard view the experiment's beam, target, spectator,
  /// Ebeam, ps and background are shown, every other key (and every <thm> option) only when the project gives it.
  void setExpertMode(bool on);
  /// Whether a key's row is shown on the form (the tests).
  bool rowShown(const QString &key, bool experiment) const;

 public slots:
  void reload();

 private:
  QWidget *editorFor(const azure2::front::FieldInfo &field, bool experiment);
  void fill();
  QString currentExperiment() const;

  ProjectHolder *holder_;
  std::vector<azure2::front::KindInfo> kinds_;
  QMap<QString, QWidget *> optionEditors_, experimentEditors_;
  QMap<QString, QString> presentWith_;  ///< experiment key -> the key it needs
  QListWidget *experimentList_;
  QLabel *message_;
  SectionTextPanel *text_;
  bool filling_ = false;
  bool expert_ = false;
  QFormLayout *optionsForm_ = nullptr, *experimentForm_ = nullptr;
};

#endif
