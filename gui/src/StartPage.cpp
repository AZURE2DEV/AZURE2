// The page shown while no project is open (StartPage.h).
#include "StartPage.h"

#include <QCheckBox>
#include <QCoreApplication>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QFrame>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QIcon>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLabel>
#include <QListWidget>
#include <QPushButton>
#include <QScrollArea>
#include <QSet>
#include <QSignalBlocker>
#include <QVBoxLayout>

namespace {

QLabel *Heading(const QString &text) {
  QLabel *l = new QLabel(text);
  QFont f = l->font();
  f.setPointSizeF(f.pointSizeF() * 1.15);
  f.setBold(true);
  l->setFont(f);
  return l;
}

QFrame *Card(QLayout *content) {
  QFrame *card = new QFrame;
  card->setProperty("role", "card");
  card->setLayout(content);
  return card;
}

QListWidget *ProjectList(const QString &name) {
  QListWidget *l = new QListWidget;
  l->setObjectName(name);
  l->setIconSize(QSize(16, 16));
  l->setUniformItemSizes(true);
  l->setAlternatingRowColors(false);
  return l;
}

/// The project roots: a checkout or an installed share/azure2 (a folder with examples/).
bool IsRoot(const QDir &d) { return d.exists("examples"); }

/// The first heading of the folder's README.md ("# Neutrons on natural neon: a worked example"), "" when it has none.
QString ReadmeTitle(const QDir &folder) {
  QFile f(folder.filePath("README.md"));
  if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) return QString();
  for (int n = 0; n < 20 && !f.atEnd(); n++) {
    const QString line = QString::fromUtf8(f.readLine()).trimmed();
    if (line.startsWith("# ")) return line.mid(2).trimmed();
    if (!line.isEmpty() && !line.startsWith('#')) break;  // text before any title: no title
  }
  return QString();
}

}  // namespace

QStringList StartPage::examplePlaces() {
  QStringList out;
  const QString env = qEnvironmentVariable("AZURE2_EXAMPLES");
  if (!env.isEmpty()) out << QDir(env).absolutePath();
  QDir d(QCoreApplication::applicationDirPath());
  for (int up = 0; up < 6; up++) {
    if (IsRoot(d)) out << d.absolutePath();
    const QDir share(d.filePath("share/azure2"));
    if (IsRoot(share)) out << share.absolutePath();
    if (!d.cdUp()) break;
  }
#ifdef AZURE2_SOURCE_DIR
  if (IsRoot(QDir(AZURE2_SOURCE_DIR))) out << QDir(AZURE2_SOURCE_DIR).absolutePath();
#endif
  out.removeDuplicates();
  return out;
}

QList<StartPage::Example> StartPage::readManifest(const QString &root, QString *why) {
  QList<Example> out;
  QFile f(QDir(root).filePath("examples/index.json"));
  if (!f.open(QIODevice::ReadOnly)) return out;
  QJsonParseError error;
  const QJsonDocument doc = QJsonDocument::fromJson(f.readAll(), &error);
  if (doc.isNull()) {
    if (why) *why = tr("%1: %2").arg(f.fileName(), error.errorString());
    return out;
  }
  for (const QJsonValue &g : doc.object().value("groups").toArray()) {
    const QJsonObject group = g.toObject();
    Example heading;
    heading.group = group.value("title").toString();
    out << heading;
    for (const QJsonValue &e : group.value("entries").toArray()) {
      const QJsonObject entry = e.toObject();
      Example x;
      x.group = heading.group;
      x.reaction = entry.value("reaction").toString();
      x.source = entry.value("source").toString();
      x.what = entry.value("what").toString();
      x.project = QDir(root).absoluteFilePath(entry.value("project").toString());
      out << x;
      for (const QJsonValue &v : entry.value("variants").toArray()) {
        Example y = x;
        y.what = v.toObject().value("what").toString();
        y.project = QDir(root).absoluteFilePath(v.toObject().value("project").toString());
        y.variant = true;
        out << y;
      }
    }
  }
  return out;
}

QList<StartPage::Example> StartPage::readManifest(QString *why) {
  for (const QString &place : examplePlaces()) {
    const QList<Example> m = readManifest(place, why);
    if (!m.isEmpty()) return m;
  }
  return {};
}

QStringList StartPage::findExamples() {
  QStringList out;
  const QList<Example> manifest = readManifest();
  for (const Example &e : manifest)
    if (!e.project.isEmpty()) out << e.project;
  if (!manifest.isEmpty()) return out;
  // no manifest: the projects of examples/*/ (and of the folders of $AZURE2_EXAMPLES)
  QSet<QString> seen;
  for (const QString &place : examplePlaces())
    for (const QString &base : {QString("examples"), QString(".")}) {
      if (base == "." && place != QDir(qEnvironmentVariable("AZURE2_EXAMPLES")).absolutePath()) continue;
      const QDir dir(QDir(place).filePath(base));
      for (const QString &folder : dir.entryList(QDir::Dirs | QDir::NoDotAndDotDot, QDir::Name)) {
        const QDir project(dir.filePath(folder));
        for (const QString &f : project.entryList({"*.azr"}, QDir::Files, QDir::Name)) {
          const QString canonical = QFileInfo(project.filePath(f)).canonicalFilePath();
          if (canonical.isEmpty() || seen.contains(canonical)) continue;
          seen.insert(canonical);
          out << QFileInfo(project.filePath(f)).absoluteFilePath();
        }
      }
    }
  return out;
}

StartPage::StartPage(QWidget *parent) : QWidget(parent) {
  setObjectName("startPage");
  QLabel *logo = new QLabel;
  logo->setPixmap(QIcon(":/azure-icon.png").pixmap(QSize(56, 56)));
  QLabel *title = new QLabel(tr("AZURE2"));
  title->setProperty("role", "title");
  QFont tf = title->font();
  tf.setPointSizeF(tf.pointSizeF() * 1.9);
  tf.setBold(true);
  title->setFont(tf);
  QLabel *subtitle = new QLabel(tr("Multi-level, multi-channel R-matrix analysis of nuclear reactions"));
  subtitle->setProperty("role", "hint");
  QVBoxLayout *titles = new QVBoxLayout;
  titles->addWidget(title);
  titles->addWidget(subtitle);
  QHBoxLayout *header = new QHBoxLayout;
  header->addWidget(logo);
  header->addLayout(titles, 1);

  new_ = new QPushButton(tr("New Project"));
  new_->setObjectName("startNew");
  new_->setProperty("primary", true);
  new_->setToolTip(tr("An empty project in project format 3 (Ctrl+N)."));
  open_ = new QPushButton(tr("Open Project..."));
  open_->setObjectName("startOpen");
  open_->setProperty("primary", true);
  open_->setToolTip(tr("Open an .azr project (Ctrl+O); an older project is offered the conversion to format 3."));
  connect(new_, &QPushButton::clicked, this, &StartPage::newRequested);
  connect(open_, &QPushButton::clicked, this, &StartPage::openRequested);
  QLabel *drop = new QLabel(tr("or drop an .azr file on this window"));
  drop->setProperty("role", "hint");
  expert_ = new QCheckBox(tr("Show expert options"));
  expert_->setObjectName("startExpert");
  expert_->setToolTip(tr("Show the options, columns, tabs and sections an occasional user does not need (also View > "
                         "Show Expert Options).  It changes what is shown, never the project."));
  connect(expert_, &QCheckBox::toggled, this, &StartPage::expertToggled);
  QHBoxLayout *buttons = new QHBoxLayout;
  buttons->addWidget(new_);
  buttons->addWidget(open_);
  buttons->addWidget(drop);
  buttons->addStretch(1);
  buttons->addWidget(expert_);

  recent_ = ProjectList("startRecent");
  examples_ = ProjectList("startExamples");
  for (QListWidget *l : {recent_, examples_}) {
    // a heading has no file: nothing to open
    connect(l, &QListWidget::itemActivated, this, [this](QListWidgetItem *item) {
      if (!item->data(Qt::UserRole).toString().isEmpty()) emit openFile(item->data(Qt::UserRole).toString());
    });
    connect(l, &QListWidget::itemClicked, this, [this](QListWidgetItem *item) {
      if (!item->data(Qt::UserRole).toString().isEmpty()) emit openFile(item->data(Qt::UserRole).toString());
    });
  }
  // the manifest's groups, entries and variants (qt_gui.md section 19.9): a heading per group, one line per reaction,
  // its variants indented under it; without a manifest, the projects found by name
  const QList<Example> manifest = readManifest();
  if (!manifest.isEmpty()) {
    for (const Example &e : manifest) {
      QListWidgetItem *item = new QListWidgetItem;
      if (e.project.isEmpty()) {
        item->setText(e.group);
        QFont f = item->font();
        f.setBold(true);
        item->setFont(f);
        item->setFlags(Qt::NoItemFlags);
        item->setForeground(palette().color(QPalette::WindowText));
      } else if (e.variant) {
        item->setText(QString("      ") + QChar(0x21B3) + " " + e.what);
        item->setData(Qt::UserRole, e.project);
        item->setToolTip(QDir::toNativeSeparators(e.project));
      } else {
        item->setIcon(QIcon(":/actions/open.svg"));
        item->setText(QString("%1   %2 %3   %4").arg(e.reaction, QChar(0x2014), e.source, QChar(0x00B7)) + " " + e.what);
        item->setData(Qt::UserRole, e.project);
        item->setToolTip(QString("%1 (%2): %3\n%4").arg(e.reaction, e.source, e.what, QDir::toNativeSeparators(e.project)));
      }
      examples_->addItem(item);
    }
  } else {
    for (const QString &f : findExamples()) {
      const QFileInfo info(f);
      QString text = info.completeBaseName();
      if (info.dir().dirName() != text) text += QString("   (%1)").arg(info.dir().dirName());
      const QString title = ReadmeTitle(info.dir());
      if (!title.isEmpty()) text += QString("  ") + QChar(0x2014) + "  " + title;
      QListWidgetItem *item = new QListWidgetItem(QIcon(":/actions/open.svg"), text);
      item->setData(Qt::UserRole, f);
      item->setToolTip(f);
      examples_->addItem(item);
    }
  }
  QVBoxLayout *recentBox = new QVBoxLayout;
  recentBox->addWidget(Heading(tr("Recent projects")));
  recentBox->addWidget(recent_, 1);
  QVBoxLayout *exampleBox = new QVBoxLayout;
  exampleBox->addWidget(Heading(tr("Example projects")));
  QLabel *exampleNote = new QLabel(examples_->count() ? tr("Projects of this installation; open one, then Save As "
                                                           "to keep your changes elsewhere.")
                                                      : tr("None found (the examples folder of an AZURE2 checkout "
                                                           "or installation, or the folder AZURE2_EXAMPLES names)."));
  exampleNote->setProperty("role", "hint");
  exampleNote->setWordWrap(true);
  exampleBox->addWidget(exampleNote);
  exampleBox->addWidget(examples_, 1);
  QHBoxLayout *lists = new QHBoxLayout;
  lists->addWidget(Card(recentBox), 1);
  lists->addWidget(Card(exampleBox), 1);

  // what each tab is for: the tabs' tool tips (qt_gui.md section 19.6), not eleven lines on this page
  tabNotes_ = new QWidget;
  tabNotes_->hide();

  QWidget *content = new QWidget;
  content->setMaximumWidth(1100);
  QVBoxLayout *layout = new QVBoxLayout(content);
  layout->setSpacing(14);
  layout->addLayout(header);
  layout->addLayout(buttons);
  notices_ = new QVBoxLayout;
  notices_->setContentsMargins(0, 0, 0, 0);
  layout->addLayout(notices_);
  layout->addLayout(lists, 1);
  QHBoxLayout *centred = new QHBoxLayout;
  centred->addStretch(1);
  centred->addWidget(content, 12);
  centred->addStretch(1);
  QWidget *inner = new QWidget;
  inner->setLayout(centred);
  QScrollArea *scroll = new QScrollArea;
  scroll->setWidget(inner);
  scroll->setWidgetResizable(true);
  scroll->setFrameShape(QFrame::NoFrame);
  QVBoxLayout *outer = new QVBoxLayout(this);
  outer->setContentsMargins(0, 0, 0, 0);
  outer->addWidget(scroll);
}

void StartPage::setExpertChecked(bool on) {
  QSignalBlocker block(expert_);
  expert_->setChecked(on);
}

void StartPage::addNotice(QWidget *notice) { notices_->addWidget(notice); }

void StartPage::setRecent(const QStringList &files) {
  recent_->clear();
  for (const QString &f : files) {
    if (!QFileInfo::exists(f)) continue;
    const QFileInfo info(f);
    QListWidgetItem *item =
        new QListWidgetItem(QIcon(":/actions/open.svg"), QString("%1   (%2)").arg(info.fileName(), info.absolutePath()));
    item->setData(Qt::UserRole, f);
    item->setToolTip(f);
    recent_->addItem(item);
  }
  if (recent_->count() == 0) {
    QListWidgetItem *none = new QListWidgetItem(tr("No recent project yet."));
    none->setFlags(Qt::NoItemFlags);
    recent_->addItem(none);
  }
}

void StartPage::setTabNotes(const QList<QPair<QString, QString>> &notes) {
  delete tabNotes_->layout();
  for (QWidget *w : tabNotes_->findChildren<QWidget *>(QString(), Qt::FindDirectChildrenOnly)) delete w;
  QGridLayout *grid = new QGridLayout(tabNotes_);
  grid->setContentsMargins(0, 0, 0, 0);
  grid->setHorizontalSpacing(16);
  grid->setColumnStretch(1, 1);
  for (int i = 0; i < notes.size(); i++) {
    QLabel *name = new QLabel(notes[i].first);
    QFont f = name->font();
    f.setBold(true);
    name->setFont(f);
    QLabel *what = new QLabel(notes[i].second);
    what->setWordWrap(true);
    grid->addWidget(name, i, 0, Qt::AlignTop);
    grid->addWidget(what, i, 1);
  }
}
