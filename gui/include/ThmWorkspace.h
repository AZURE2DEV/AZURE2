#ifndef THMWORKSPACE_H
#define THMWORKSPACE_H

#include <QDialog>

#include "ThmSettings.h"

QT_BEGIN_NAMESPACE
class QLabel;
class QPushButton;
class QTabWidget;
class QUrl;
QT_END_NAMESPACE

class AZURESetup;
class ThmChannelsPage;
class ThmDiagnosticsPage;
class ThmExperimentsPage;
class ThmModelPage;

/*!
 * The THM workspace (Configure > THM Workspace...): everything that belongs
 * to the Trojan Horse (HOES) observable, in one window, apart from the classic
 * editor.  Pages:
 *   Model        the options of the <thm> block (ThmModelPage);
 *   Experiments  its experiment[<name>] lines (ThmExperimentsPage);
 *   Channels     the THM columns of the <levels> lines: binding energy B
 *                (field 32) and width input flag (field 33) (ThmChannelsPage);
 *   Diagnostics  read-only plots the engine computes on request for one THM
 *                segment: entrance vertex, HOES and on-shell cross sections,
 *                line shape, weight, spectator-momentum window
 *                (ThmDiagnosticsPage; built with USE_API,
 *                as the Plot tab is with USE_QWT).
 * It writes nothing else, and nothing at all when its values are left as
 * they were read: such a project saves byte for byte as before.  Without a
 * THM segment the pages are disabled under a one-line note.  The pages carry
 * no explanatory text of their own: tooltips, and a Help button that opens
 * the user guide.
 */
class ThmWorkspace : public QDialog {
  Q_OBJECT

 public:
  /// `settings`: the project's <thm> block as read (AZURESetup::thmSettings).
  ThmWorkspace(AZURESetup *setup, const ThmSettings &settings, QWidget *parent = 0);

  /// True if the project has a THM data or test segment (the pages are enabled).
  bool hasThmSegments() const { return hasThm_; }
  /// The first thing AZURE2 would refuse, page by page ("" if none); shows that page.
  QString validate();
  /// Installs the pages' values into the project (no check; see validate).
  void apply();
  /*! The project as Accept followed by a save would write it: the <thm> block
      as the Model and Experiments pages set it, the Channels page's values
      installed for the time of writing only.  The project is left as it was.
      False (and the reason) for a project that has no file yet. */
  bool projectSnapshot(QString &text, QString *error = nullptr);
  /*! What the pages hold, as text: changes when an edit does.  Cheap, unlike
      projectSnapshot, which writes the whole project (the workspace is
      modal, so nothing else changes the project while it is open). */
  QString editState() const;

  QLabel *noThmLabel;
  QTabWidget *pages;
  ThmModelPage *modelPage;
  ThmExperimentsPage *experimentsPage;
  ThmChannelsPage *channelsPage;
  /// Null in a build without the engine API (USE_API=OFF): the page needs it.
  ThmDiagnosticsPage *diagnosticsPage = nullptr;
  QPushButton *acceptButton;
  QPushButton *helpButton;  ///< opens helpUrl()
  /// The user guide's THM workspace section: a local build of the docs if found, else online.
  static QUrl helpUrl();

 public slots:
  void accept() override;

 private:
  AZURESetup *setup_;
  bool hasThm_;
};

#endif
