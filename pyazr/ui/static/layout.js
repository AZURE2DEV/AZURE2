// Layout hints of the generated forms (docs/upgrade/web_ui.md section 1.3): which fields of a record kind sit on one
// line.  The only hand-written knowledge of record kinds in the page; a kind without a row here still gets its full
// form, one field per line in the schema's order (forms.js, arrange).  A field named here that the schema lacks is
// skipped, so a renamed field costs a line of this table, never a broken form.
export const LAYOUT = {
  dataset: [["reaction", "quantity", "unit"], ["entrance_pair", "exit_pair"], ["file", "covariance"],
    ["energy_window", "angle_window"], ["target", "beam", "detector"], ["source", "exfor", "doi"]],
  extrapolation: [["reaction", "quantity", "unit"], ["entrance_pair", "exit_pair"], ["target", "beam", "detector"]],
  target: [["kind", "active", "phase"], ["composition", "compound_correction"]],
  datum: [["quantity", "level", "J"], ["value", "error"], ["imag", "imag_error"], ["source", "exfor", "doi"]],
};

// The navigation (docs/upgrade/web_navigation.md): one row of main tabs in the old setup program's order and words
// (Particle Pairs, Levels and Channels, Segments, Experimental Effects, Calculate and Fit, Plot), then the few tabs of
// what the old program did not have, each a group with a second row of its own (Rate and ENDF, Physics and THM,
// Covariances), and the Project tab (the old program's Configure menu: the run's settings, the problems, the history,
// the text).  The same tabs in both depths: the expert's depth shows every field and column, never other tabs.
//
// Each entry is a view of one of the page's tabs (app.js TABS): `tab` the page's tab, `hash` the route it opens,
// `part` the first word of the routes it stands for (the Nucleus tab's "pair"), `sections` the Data tab's groups it
// shows (null: every group no other view of the tab names), `group` the main tab it sits under (its title; none: a
// main tab of its own), `when` a condition of the project under which it is shown ("thm": the project has THM
// content, the THM tab's own rule), `doc` what it is for (the overview and the link's title).
export const FIRST_VIEW = [
  { title: "Particle Pairs", tab: "nucleus", part: "pair", hash: "#/nucleus/pair",
    doc: "the entrance and exit pairs: partners, masses, spins, separation energies, channel radii" },
  { title: "Levels and Channels", tab: "nucleus", hash: "#/nucleus/level",
    doc: "the levels of the compound nucleus, each with its channels and their widths beside it" },
  { title: "Segments", tab: "data", sections: null, hash: "#/data?section=data",
    doc: "the data sets (segments from data), the extrapolations (without data) and model data" },
  { title: "Experimental Effects", tab: "data", sections: ["targets"], hash: "#/data?section=targets",
    doc: "targets, beams and detectors the data sets name (energy loss, resolution)" },
  { title: "Calculate and Fit", tab: "fit", hash: "#/fit",
    doc: "run a fit (and MINOS, the sampler, the bands) and see every parameter with its error" },
  { title: "Plot", tab: "results", hash: "#/results",
    doc: "every data set with the model and the residuals, the χ² per data set" },
  { title: "Reaction Rate", group: "Rate and ENDF", tab: "rate", hash: "#/rate",
    doc: "the thermonuclear reaction rate N_A<σv> with its band" },
  { title: "ENDF", group: "Rate and ENDF", tab: "endf", hash: "#/endf",
    doc: "ENDF-6 evaluations over the data, their check, a project to and from ENDF-6" },
  { title: "Physics functions", group: "Physics and THM", tab: "physics", hash: "#/physics",
    doc: "the engine's functions (penetrabilities, Coulomb functions, stopping) plotted on a range" },
  { title: "THM", group: "Physics and THM", tab: "thm", hash: "#/thm", when: "thm",
    doc: "the Trojan Horse model options, experiments and diagnostics (a project with THM content)" },
  { title: "Data covariance", group: "Covariances", tab: "covariance", hash: "#/covariance",
    doc: "a data set's covariance built from parts (also: a data set's \"Covariance…\")" },
  { title: "Parameter prior", group: "Covariances", tab: "prior", hash: "#/prior",
    doc: "a correlated Gaussian prior on the fit parameters" },
  { title: "Project", tab: "project", hash: "#/project",
    doc: "the file, the next step, the problems, the run's settings, the history and the text" },
];

/** The route of a page opened without one (the printed URL): the first tab of the old program's order. */
export const START = FIRST_VIEW[0].hash;

/** The old program's words for the groups of the Data tab and of the Nucleus tab's further records (a group without
 * a row is shown by its section's name). */
export const SECTION_TITLES = { data: "Segments From Data", extrapolations: "Segments Without Data",
  modelData: "Model data", targets: "Experimental Effects", potentials: "Potentials", gammaDecays: "γ decays",
  channelRadii: "Channel radii", isospin: "Isospin" };

/** The tabs, each {title, hash, tab, part?, sections?, group?, row}: the views of FIRST_VIEW, each with the row of
 * `tabs` (app.js TABS) it is a view of.  The same in both depths (`expert` is kept for the callers). */
export function shownTabs(expert, tabs) {
  return FIRST_VIEW.map((v) => ({ ...v, row: tabs.find((t) => t.id === v.tab) }));
}

/** The main row: one entry per main tab, a group as {title, group: true, members, hash of its first member}. */
export function mainTabs(shown) {
  const out = [];
  for (const v of shown) {
    if (!v.group) { out.push({ ...v, members: [v] }); continue; }
    const g = out.find((m) => m.group && m.title === v.group);
    if (g) g.members.push(v);
    else out.push({ title: v.group, group: true, hash: v.hash, tab: v.tab, members: [v] });
  }
  return out;
}

/** The main tab a view sits under (its group's title, or its own). */
export const mainTitleOf = (view) => (view ? view.group || view.title : null);

/** The shown tab a route belongs to: the one of its tab whose part is the route's first word (a pair) or whose
 * sections hold the group the tab drew (`section`, the Data tab's); else the view of its tab without either. */
export function tabOfRoute(shown, route, section) {
  const mine = shown.filter((t) => t.tab === route.tab);
  return mine.find((t) => t.part && t.part === route.args[0])
    || mine.find((t) => t.sections && section && t.sections.includes(section))
    || mine.find((t) => !t.part && !t.sections) || mine[0] || null;
}

/** The tab a key goes to among the tabs shown (`links`, in the bar's order): 1-9 the first nine, 0 the tenth, [ and ]
 * the one before and after the active tab (around the ends). */
export function tabByKey(links, key) {
  if (/^[1-9]$/.test(key)) return links[+key - 1] || null;
  if (key === "0") return links[9] || null;
  const at = links.findIndex((a) => a.classList?.contains("active"));
  const step = key === "]" ? 1 : -1;
  return links.length ? links[((at < 0 ? 0 : at + step) + links.length) % links.length] : null;
}

/** The link of the tab row an arrow key moves the focus to (the ARIA tab pattern: ← → around the ends, Home, End). */
export function tabByArrow(links, at, key) {
  if (!links.length) return null;
  if (key === "Home") return links[0];
  if (key === "End") return links[links.length - 1];
  const step = key === "ArrowRight" ? 1 : key === "ArrowLeft" ? -1 : 0;
  return step ? links[(at + step + links.length) % links.length] : null;
}
