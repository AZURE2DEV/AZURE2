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

// The first view (docs/upgrade/web_ui.md section 6, the two depths): the old setup program's path, in its order and
// words, each tab a view of one of the page's tabs (app.js TABS): `tab` the page's tab, `hash` the route it opens,
// `part` the first word of the routes it stands for (the Nucleus tab's "pair" or "level"), `sections` the Data tab's
// groups it shows (null: every group no other view of the tab names), `when` a condition of the project under which
// the view is shown ("thm": the project has THM content, the THM tab's own rule).  The expert's view shows TABS as
// they are.
export const FIRST_VIEW = [
  { title: "Project", tab: "project", hash: "#/project" },
  { title: "Particle Pairs", tab: "nucleus", part: "pair", hash: "#/nucleus/pair" },
  { title: "Levels and Channels", tab: "nucleus", hash: "#/nucleus/level" },
  { title: "Segments", tab: "data", sections: null, hash: "#/data?section=data" },
  { title: "Experimental Effects", tab: "data", sections: ["targets"], hash: "#/data?section=targets" },
  { title: "Fit", tab: "fit", hash: "#/fit" },
  { title: "Plot", tab: "results", hash: "#/results" },
  // shown only when the project has THM content (`when`: app.js shows such a view where the condition's call says so)
  { title: "Physics", tab: "physics", hash: "#/physics", when: "thm" },
  { title: "THM", tab: "thm", hash: "#/thm", when: "thm" },
  { title: "Reaction Rate", tab: "rate", hash: "#/rate" },
];

/** The old program's words for the groups of the Data tab (a group without a row is shown by its section's name). */
export const SECTION_TITLES = { data: "Segments From Data", extrapolations: "Segments Without Data",
  modelData: "Model data", targets: "Experimental Effects" };

/** The tabs of a depth, each {title, hash, tab, part?, sections?, row}: the first view's (FIRST_VIEW, each with the
 * row of `tabs` it is a view of) or the expert's (`tabs`, app.js TABS, as they are). */
export function shownTabs(expert, tabs) {
  if (expert) return tabs.map((t) => ({ title: t.title, hash: `#/${t.id}`, tab: t.id, row: t }));
  return FIRST_VIEW.map((v) => ({ ...v, row: tabs.find((t) => t.id === v.tab) }));
}

/** The shown tab a route belongs to: the one of its tab whose part is the route's first word (a pair) or whose
 * sections hold the group the tab drew (`section`, the Data tab's); else the view of its tab without either. */
export function tabOfRoute(shown, route, section) {
  const mine = shown.filter((t) => t.tab === route.tab);
  return mine.find((t) => t.part && t.part === route.args[0])
    || mine.find((t) => t.sections && section && t.sections.includes(section))
    || mine.find((t) => !t.part && !t.sections) || mine[0] || null;
}
