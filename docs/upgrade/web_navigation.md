# The web page's navigation against the Qt program's

The PI, who has used both, found the web page (`python -m pyazr ui PROJECT.azr`) "in-intuitive and difficult to
navigate into, with respect to the old Qt". This note records what a new user meets on the page as it was at
`98a4cf60` (walked on 2026-11-18 in headless Chrome 124 over its DevTools protocol, on copies of `tests/upgrade/T1`
and `T2`, from the printed URL, in a fresh profile; screenshots in `work/web-shots/before/`, outside the
repository), what was changed, and the number of clicks each task of the Qt program takes before and after.
The reference for "easy" is the Qt program of `gh/dev:gui/src` (`AZURESetup.cpp`): one row of tabs in the order of
the work (Particle Pairs, Levels and Channels, Segments, Experimental Effects, [Nuclear Potential], Fitting
Settings, Calculate, [MCMC], Plot), each tab one table with its detail beside it and + / − at the table, the menus
File, Configure and Documentation, a status bar.

## 1. What was hard (the page of `98a4cf60`)

- **Two navigations.** The page opened in a "first view" (Project, Particle Pairs, Levels and Channels, Segments,
  Experimental Effects, Fit, Plot, ENDF, Reaction Rate) and an "expert" check box replaced the whole tab row by
  another one with other names (Project, Nucleus, Data, Fit, Results, Physics, Rate, ENDF). In the expert's row the
  pairs were at the bottom of the Nucleus tab, under the levels and the γ-decay records.
- **Hidden things.** The Physics tab was not in the first view of a project without THM content (T1): two clicks
  through the expert switch. The covariance and prior builders had no tab (reached only from a data set's
  "Covariance…" and from a link above the Fit tab's parameter table). Channel radii and isospin records appeared as
  groups of Segments under the schema's names `channelRadii` and `isospin`. There was no control to remove a level
  (the route `DELETE /api/project/levels/{label}` existed), none to open another project or to save under another
  name (routes `POST /api/project/open` and `save` with a path existed).
- **An extra step.** In "+ level" the channels a J^π allows were shown only after a press of "Channels"; "Add
  level" without it gave a level of an existing J^π its group's channels unseen, and for a new J^π the engine's
  refusal. After "Add level" the page went back to the list, the new level not selected.
- **Detail below, not beside.** A level's channels, a pair's numbers and a record's form came under the table, the
  form of a data set 2200 px down the page with its plot at the bottom.
- **State not shown.** A running job was visible only on the Fit tab (the header said "working…" for any request);
  the validator's messages only in the Project tab's "Problems" section; the Plot tab stayed an empty page for the
  3 to 40 s its curves take, with "working…" in the header.
- **Names.** "Fit" for the Qt program's Fitting Settings and Calculate; "+ new" at the far right of a group bar,
  away from the table it adds to; the expert's Nucleus, Data, Results.
- What already worked: every view had its address (`#/data/el_lab090`), so Reload and Back kept the place; tables
  were reachable by keyboard (↑ ↓ Enter, `/`); number keys switched tabs; the unsaved state was in the header.

## 2. What was changed

- One tab row in the Qt order and words, the same in both depths: Particle Pairs, Levels and Channels, Segments,
  Experimental Effects (the Qt name; the brief's "Target Effects" was not used), Calculate and Fit, Plot; then three
  groups with a second row: Rate and ENDF (Reaction Rate, ENDF), Physics and THM (Physics functions; THM where the
  project has THM content), Covariances (Data covariance, Parameter prior); and Project. The expert's depth now
  changes fields and columns only. The page opens on Particle Pairs, as the Qt program did.
- Menus File (Open…, New project from this one…, Save, Save as…, Write the output files, ENDF import and export,
  About), Edit (Undo, Redo, History, the text), Configure (Run settings, Session options, Prior, Engine), View
  (every field, theme, the validator's messages, the job), Help (overview, documentation of this tab, keys, About);
  each opens by a click, ↓ or Enter, never by hovering.
- The top bar: name, saved state, Undo, Redo, Save, Validate, Calculate, χ². The status line: the engine's last
  message, the job that runs (progress, Stop, a link to its panel) or ran, the validator's verdict opening its list,
  each message a link to its record or level.
- Particle Pairs, Levels and Channels: the table on the left with Add / Remove at it, the selected row beside it; a
  new level's channels listed as J and parity are typed, the new level opened after Add; Remove level. Channel radii
  and isospin records with the γ decays and potentials under the levels' "More". Segments: "Add" at the list.
- A first visit from the printed URL shows an overview of the tabs in the order of the work (dismissible; Help
  reopens it); a tab that waits for the engine says so; ← → Home End move along a tab row.
- Nothing of the `/api` routes changed. Code: `pyazr/ui/static/layout.js` (the tab table), `shell.js` (menus, status
  marks, overview), `app.js`, `tabs/nucleus.js` and `tabs/levels.js`. Test: `tests/pyazr/ui/navigation_test.py`.

## 3. Clicks from the opened page to each task, before and after

Clicks only (typing into a form's boxes not counted), from the page as the printed URL opens it, with T1.

| Task (the Qt program's) | before | after |
|---|---|---|
| add a particle pair | 3 (Particle Pairs, + pair, Add pair) | 3 (Particle Pairs, Add pair, Add pair) |
| add a level and see its channels | 5 (tab, + level, Channels, Add level, the new row) | 3 (tab, Add level, Add level: the level opens with its channels) |
| remove a level | no control (edit the whole text) | 3 + confirm (tab, the row, Remove level) |
| add a data segment | 3 (Segments, + new, Add) | 3 (Segments, Add to Segments From Data, Add) |
| add a target and name it on a data set | 7 | 7 |
| run a fit | 2 (Fit, Run) | 2 (Calculate and Fit, Run) |
| see that a fit runs, from another tab | 1 (Fit) | 0 (the status line) |
| stop a running fit, from another tab | 2 (Fit, Stop) | 1 (Stop in the status line) |
| see the plot | 1 (Plot or Calculate) | 1 |
| see the parameters after the fit | 1 (Fit; the table below the job) | 1 (Calculate and Fit, or the job in the status line) |
| see whether the project is valid, from any tab | 1 (Project; scroll to Problems) | 0 (the status line) |
| go from a validator message to its row | 2 (Project, → go to) | 2 (the status mark, the message's link) |
| run settings (Qt: Configure › Runtime Options) | 1 (Project; scroll past three sections) | 2 (Configure, Run settings: it scrolls there) |
| open another project / save as | no control | 3 (File, Open… or Save as…, OK) |
| Physics functions on T1 | 2 (expert, Physics) | 1 (Physics and THM) |
| THM on T2 | 1 | 2 (Physics and THM, THM) |
| Reaction Rate / ENDF | 1 / 1 | 1 / 2 |
| covariance builder of a data set | 3 (Segments, the data set, Covariance…) | 3, or 2 (Covariances, the data set) |

Two tasks take one click more than before (THM, ENDF: they sit in a group's second row); the tasks of the Qt
program's own tabs take as many clicks or fewer, and four (remove a level, open, save as, Physics on a project
without THM) became possible or visible. The walk of `navigation_test.py` (the printed URL, the tab row, each main
tab's table, a level chosen into the address, Reload, Back, Calculate, the validator's list, a fit followed and
stopped from the status line on another tab) passes on this branch and fails on `98a4cf60` at its first assertion.
Screenshots of every tab after the change: `work/web-shots/after/`.
