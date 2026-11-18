// The page's navigation (docs/upgrade/web_navigation.md): its tabs are the old program's path, then a few groups,
// each a view of a tab of the page, a route marks the view it belongs to; a form folds the fields the engine's
// table marks rare unless the record sets them; the next step of a project; Calculate's comparison of two
// calculations (the engine's numbers, subtracted).
import { test } from "node:test";
import assert from "node:assert/strict";
import { engineSchema, mod } from "./common.mjs";

const layout = await mod("layout.js");
const calc = await mod("calculate.js");
const forms = await mod("forms.js");
const data = await mod("tabs/data.js");
const project = await mod("tabs/project.js");
const schema = engineSchema();
const skipNoSchema = schema ? false : "the engine's schema was not given (AZURE2_UI_SCHEMA)";

const TABS = ["project", "nucleus", "data", "fit", "results", "physics", "thm", "endf", "rate", "covariance", "prior"]
  .map((id) => ({ id, title: id }));

test("the tabs are the old program's path, then a few groups; the same in both depths", () => {
  const shown = layout.shownTabs(false, TABS);
  assert.deepEqual(layout.mainTabs(shown).map((t) => t.title), ["Particle Pairs", "Levels and Channels", "Segments",
    "Experimental Effects", "Calculate and Fit", "Plot", "Rate and ENDF", "Physics and THM", "Covariances",
    "Project"]);
  assert.deepEqual(layout.mainTabs(shown).filter((t) => t.group).map((t) => [t.title, t.members.map((m) => m.title)]),
    [["Rate and ENDF", ["Reaction Rate", "ENDF"]], ["Physics and THM", ["Physics functions", "THM"]],
      ["Covariances", ["Data covariance", "Parameter prior"]]]);
  // THM only where the project has THM content (app.js shows a view with `when` by its condition)
  assert.deepEqual(shown.filter((t) => t.when).map((t) => [t.title, t.when]), [["THM", "thm"]]);
  // every tab of the page has a view, each view a tab of the page, and a sentence for the overview
  assert.deepEqual([...new Set(shown.map((t) => t.tab))].sort(), TABS.map((t) => t.id).sort());
  for (const t of shown) assert.ok(t.row && t.hash.startsWith(`#/${t.tab}`) && t.doc, t);
  assert.deepEqual(layout.shownTabs(true, TABS), shown);
  // a page opened without a route opens the first tab of the work
  assert.equal(layout.START, "#/nucleus/pair");
  assert.equal(layout.mainTitleOf(shown.find((t) => t.tab === "thm")), "Physics and THM");
  assert.equal(layout.mainTitleOf(shown.find((t) => t.tab === "results")), "Plot");
});

test("a route marks the view it belongs to", () => {
  const first = layout.shownTabs(false, TABS);
  const on = (tab, args, section) => layout.tabOfRoute(first, { tab, args }, section)?.title;
  assert.equal(on("nucleus", ["pair", "2"]), "Particle Pairs");
  assert.equal(on("nucleus", ["level", "1- 12.4"]), "Levels and Channels");
  assert.equal(on("nucleus", []), "Levels and Channels");
  assert.equal(on("nucleus", ["some_gamma_decay"]), "Levels and Channels");
  assert.equal(on("data", ["tgt"], "targets"), "Experimental Effects");
  assert.equal(on("data", ["s01"], "data"), "Segments");
  assert.equal(on("data", [], "a_new_section"), "Segments");
  assert.equal(on("fit", ["job", "3"]), "Calculate and Fit");
  assert.equal(on("results", []), "Plot");
  assert.equal(on("physics", []), "Physics functions");
  assert.equal(on("endf", ["check"]), "ENDF");
  assert.equal(on("rate", []), "Reaction Rate");
  assert.equal(on("covariance", ["s01"]), "Data covariance");
  assert.equal(on("prior", []), "Parameter prior");
});

test("the arrow keys move along a tab row, around the ends", () => {
  const links = ["a", "b", "c"];
  assert.equal(layout.tabByArrow(links, 0, "ArrowRight"), "b");
  assert.equal(layout.tabByArrow(links, 0, "ArrowLeft"), "c");
  assert.equal(layout.tabByArrow(links, 2, "ArrowRight"), "a");
  assert.equal(layout.tabByArrow(links, 1, "Home"), "a");
  assert.equal(layout.tabByArrow(links, 1, "End"), "c");
  assert.equal(layout.tabByArrow(links, 1, "x"), null);
  assert.equal(layout.tabByArrow([], 0, "Home"), null);
});

test("the Data tab's views split its groups; a group the schema gains goes to Segments", () => {
  const sections = [["data", ["a"]], ["targets", ["b"]], ["extrapolations", ["c"]], ["new_section", ["d"]]];
  const views = layout.FIRST_VIEW.filter((v) => v.tab === "data");
  const effects = views.find((v) => v.sections);
  const segments = views.find((v) => !v.sections);
  assert.deepEqual(data.sectionsOfView(sections, effects).map(([s]) => s), ["targets"]);
  assert.deepEqual(data.sectionsOfView(sections, segments).map(([s]) => s), ["data", "extrapolations",
    "new_section"]);
  assert.deepEqual(data.sectionsOfView(sections, null), sections);
});

test("a form folds the rare fields the record does not set, never a set or required one", { skip: skipNoSchema },
  () => {
    const kinds = Object.entries(schema.records).filter(([, r]) => Object.values(r.fields).some((f) => f.rare));
    assert.ok(kinds.length >= 4, "the engine's table marks rare fields in at least four record kinds");
    for (const [kind, rec] of kinds) {
      const form = forms.describeForm(schema, kind, {}, { records: [], foldRare: true });
      const more = form.more.flat().map((e) => e.key);
      const lines = form.lines.flat().map((e) => e.key);
      for (const k of more) assert.ok(rec.fields[k].rare && !rec.fields[k].required, `${kind}.${k}`);
      for (const k of lines) assert.ok(!rec.fields[k].rare || rec.fields[k].required, `${kind}.${k} not folded`);
      // the expert's form (no folding) has every field on its lines
      const all = forms.describeForm(schema, kind, {}, { records: [] });
      assert.deepEqual(all.more.flat(), []);
      assert.equal(all.lines.flat().length, lines.length + more.length);
      // a rare field the record sets is shown, and an edit of it is still found
      if (more.length) {
        const k = more[0];
        const set = forms.describeForm(schema, kind, { [k]: "x" }, { records: [], foldRare: true });
        assert.ok(set.lines.flat().some((e) => e.key === k), `${kind}.${k} set but folded`);
        const changed = forms.changedFields(form, { [k]: ["y"] });
        assert.deepEqual(changed, { [k]: rec.fields[k].repeatable ? ["y"] : "y" }, "a folded field's edit is lost");
      }
    }
  });

test("the next step follows the old program's order", () => {
  const step = (counts, err) => project.nextStep(counts, err).label;
  assert.equal(step({ pairs: 0, levels: 0, data: 0 }), "Particle Pairs");
  assert.equal(step({ pairs: 2, levels: 0, data: 3 }), "Levels and Channels");
  assert.equal(step({ pairs: 2, levels: 4, data: 0 }), "Segments");
  assert.equal(step({ pairs: 2, levels: 4, data: 3 }, { what: "refused" }), "Problems");
  assert.equal(step({ pairs: 2, levels: 4, data: 3 }), "Calculate");
  for (const s of [{ pairs: 0 }, { pairs: 1, levels: 0 }, { pairs: 1, levels: 1, data: 0 }]) {
    const hash = project.nextStep(s).hash;
    assert.ok(layout.FIRST_VIEW.some((v) => v.hash === hash), hash);
  }
});

test("Calculate compares two calculations: the engine's numbers subtracted, the largest change first", () => {
  const ev = (chi2, sets) => ({ chi2, objective: chi2 + 1, datasets: Object.entries(sets).map(([name, c]) => (
    { name, chi2: c, n: 10 })) });
  const a = calc.snapshotOf(ev(30, { s1: 10, s2: 20 }), 4);
  const b = calc.snapshotOf(ev(25.5, { s1: 10.5, s2: 15 }), 5);
  let seen = calc.remember([], a);
  seen = calc.remember(seen, a);
  assert.equal(seen.length, 1, "the same calculation twice is one");
  assert.equal(calc.baselineOf(seen, a), null);
  seen = calc.remember(seen, b);
  assert.equal(calc.baselineOf(seen, b), a);
  const c = calc.changesOf(a, b);
  assert.equal(c.chi2, 25.5 - 30);
  assert.deepEqual(c.sets.map((s) => [s.name, s.change]), [["s2", 15 - 20], ["s1", 10.5 - 10]]);
  assert.deepEqual([c.from, c.to], [4, 5]);
  // a data set only in one of the two comes first, its change unknown
  const d = calc.changesOf(a, calc.snapshotOf(ev(31, { s1: 10, s2: 20, s3: 1 }), 6));
  assert.equal(d.sets[0].name, "s3");
  assert.equal(d.sets[0].change, null);
  assert.deepEqual([calc.changeText(-4.5), calc.changeText(0.5), calc.changeText(0), calc.changeText(null),
    calc.changeText(120000), calc.changeText(1.5e-7)], ["−4.5", "+0.5", "0", "—", "+120000", "+1.5e-7"]);
});

test("keys: 1-9 the first nine tabs shown, 0 the tenth, [ and ] the one before and after the active one", () => {
  const links = Array.from({ length: 11 }, (_, i) => ({ id: i, classList: { contains: (c) => c === "active"
    && i === 10 } }));
  assert.equal(layout.tabByKey(links, "1").id, 0);
  assert.equal(layout.tabByKey(links, "9").id, 8);
  assert.equal(layout.tabByKey(links, "0").id, 9);
  assert.equal(layout.tabByKey(links, "]").id, 0, "after the last, the first");
  assert.equal(layout.tabByKey(links, "[").id, 9);
  assert.equal(layout.tabByKey(links.slice(0, 7), "8"), null);
  assert.equal(layout.tabByKey(links.slice(0, 7), "0"), null);
});
