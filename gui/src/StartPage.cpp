// The page shown while no project is open (StartPage.h).
#include "StartPage.h"

#include <QCheckBox>
#include <QCoreApplication>
#include <QDir>
#include <QFileInfo>
#include <QFrame>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QIcon>
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

/// The project roots: a checkout (examples/ and tests/upgrade/) or an installed share/azure2.
bool IsRoot(const QDir &d) { return d.exists("examples") || d.exists("tests/upgrade"); }

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

QStringList StartPage::findExamples() {
  QStringList out;
  QSet<QString> seen;
  auto take = [&](const QString &file) {
    const QString canonical = QFileInfo(file).canonicalFilePath();
    if (canonical.isEmpty() || seen.contains(canonical)) return;
    seen.insert(canonical);
    out << QFileInfo(file).absoluteFilePath();
  };
  for (const QString &place : examplePlaces()) {
    const QDir root(place);
    // the env's directory may itself hold project folders
    for (const QString &base : {QString("examples"), QString("tests/upgrade"), QString("tests"), QString(".")}) {
      const QDir dir(root.filePath(base));
      if (!dir.exists()) continue;
      if (base == "tests") {  // the 13N project of the test suite only
        if (dir.exists("13N/13N.azr")) take(dir.filePath("13N/13N.azr"));
        continue;
      }
      if (base == "." && place != QDir(qEnvironmentVariable("AZURE2_EXAMPLES")).absolutePath()) continue;
      for (const QString &folder : dir.entryList(QDir::Dirs | QDir::NoDotAndDotDot, QDir::Name)) {
        if (base == "tests/upgrade" && !(folder.size() == 2 && folder[0] == 'T')) continue;  // T1, T2, T3
        const QDir project(dir.filePath(folder));
        for (const QString &f : project.entryList({"*.azr"}, QDir::Files, QDir::Name)) take(project.filePath(f));
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
    connect(l, &QListWidget::itemActivated, this,
            [this](QListWidgetItem *item) { emit openFile(item->data(Qt::UserRole).toString()); });
    connect(l, &QListWidget::itemClicked, this,
            [this](QListWidgetItem *item) { emit openFile(item->data(Qt::UserRole).toString()); });
  }
  for (const QString &f : findExamples()) {
    const QFileInfo info(f);
    QListWidgetItem *item = new QListWidgetItem(QIcon(":/actions/open.svg"),
                                                QString("%1   (%2)").arg(info.fileName(), info.dir().dirName()));
    item->setData(Qt::UserRole, f);
    item->setToolTip(f);
    examples_->addItem(item);
  }
  QVBoxLayout *recentBox = new QVBoxLayout;
  recentBox->addWidget(Heading(tr("Recent projects")));
  recentBox->addWidget(recent_, 1);
  QVBoxLayout *exampleBox = new QVBoxLayout;
  exampleBox->addWidget(Heading(tr("Example projects")));
  QLabel *exampleNote = new QLabel(examples_->count() ? tr("Projects of this installation; open one, then Save As "
                                                           "to keep your changes elsewhere.")
                                                      : tr("None found (examples/ of an AZURE2 checkout, or "
                                                           "AZURE2_EXAMPLES=<directory>)."));
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
