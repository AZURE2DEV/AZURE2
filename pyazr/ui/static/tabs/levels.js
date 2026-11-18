// The levels of the Levels and Channels tab (docs/upgrade/web_ui.md section 6.2): one level's panel (its energy,
// its fixed and active flags, its channels with width, fixed flag, amplitude flag, prior and limits, each editable)
// and the form of a new level (J, parity, E; the channels its J^π allows with the project's pairs, offered as they
// are typed).  Every edit is a Project call through the project's routes; the status line names the χ² before and
// after it (the engine's numbers).

import { pathOf } from "../api.js";
import { formatNumber, renderTable } from "../tables.js";
import { button, h } from "../ui.js";

/** A number typed as a number; any other text as typed (pyazr converts it or refuses it with its words). */
export const valueOf = (t) => (t !== "" && Number.isFinite(Number(t)) ? Number(t) : t);

/** The body of POST /api/project/levels: J, parity, E (MeV), whether its energy is fixed, and the channels checked
 * among those allowed_channels gave ({pair, l, s}, width 0, fixed) -- or "allowed" when every one is checked
 * (add_level then gives a J^π the project has its group's channels). */
export function newLevelRequest({ J, parity, energy, fixed }, allowed, checked) {
  const chosen = allowed.filter((c, i) => checked[i]);
  return { J: valueOf(String(J).trim()), parity: Number(parity), energy: valueOf(String(energy).trim()),
    fixed: Boolean(fixed), channels: chosen.length === allowed.length && allowed.length ? "allowed"
      : chosen.map((c) => ({ pair: c.pair, l: c.l, s: c.s, width: 0, fixed: true })) };
}

/** The body of PUT /api/project/parameter_settings from the four boxes of a parameter's prior and limits (texts as
 * typed): the prior [centre, sigma] when both are typed, the limits [lower, upper] when both are, null otherwise
 * (both null remove the row); `refused` names a pair of boxes typed half.  The numbers are the engine's to check. */
export function settingRequest(level, channel, { centre, sigma, lower, upper }) {
  const pair = (a, b, what) => {
    const [x, y] = [String(a ?? "").trim(), String(b ?? "").trim()];
    if (!x && !y) return { value: null };
    if (!x || !y) return { refused: `${what}: both boxes or neither` };
    return { value: [Number(x), Number(y)] };
  };
  const prior = pair(centre, sigma, "the prior (centre and σ)");
  const limits = pair(lower, upper, "the limits (lower and upper)");
  const refused = prior.refused || limits.refused;
  if (refused) return { refused };
  return { level, channel: channel ?? null, prior: prior.value, limits: limits.value };
}

/** The row of `settings` (GET /api/project/parameter_settings) that names the energy of `level` (channel null) or its
 * channel [pair, l, s]; undefined when there is none. */
export function settingOf(settings, level, channel) {
  return (settings || []).find((r) => r.level === level && (channel ? Array.isArray(r.channel)
    && r.channel[0] === channel[0] && r.channel[1] === channel[1] && r.channel[2] === channel[2] : !r.channel));
}

/** A row's prior and limits in words ("prior 1.2 ± 0.1; limits 0 … 2"), or "" without either. */
export function settingText(row) {
  if (!row) return "";
  const parts = [];
  if (row.nuisance) parts.push(`prior ${formatNumber(row.value, 6)} ± ${formatNumber(row.error, 4)}`);
  if (!(row.lower === 0 && row.upper === 0)) {
    parts.push(`limits ${formatNumber(row.lower, 6)} … ${formatNumber(row.upper, 6)}`);
  }
  return parts.join("; ");
}


/** One level: its energy, fixed and active flags, and its channels with width and fixed flag, each editable. */
export function levelPanel(ctx, level, pairLabel, msg, settings) {
  const path = pathOf("/api/project/levels", level.label);
  const box = h("section.panel.level", { "aria-label": `level ${level.label}` });
  const settingForm = h("div.setting-form");
  const energy = h("input.num", { type: "text", value: level.energy_text, size: 12, "aria-label": "energy (MeV)" });
  const fixed = h("input", { type: "checkbox", "aria-label": "energy fixed" });
  fixed.checked = level.fixed;
  const active = h("input", { type: "checkbox", "aria-label": "active" });
  active.checked = level.active;
  const applyLevel = () => {
    const changed = {};
    if (energy.value.trim() !== level.energy_text) changed.energy = Number(energy.value);
    if (fixed.checked !== level.fixed) changed.fixed = fixed.checked;
    if (active.checked !== level.active) changed.active = active.checked;
    if (Object.keys(changed).length) chi2Edit(ctx, `set level ${level.label}`, () => ctx.api.patch(path, changed), msg);
  };
  energy.addEventListener("keydown", (e) => { if (e.key === "Enter") applyLevel(); });
  box.append(h("h2", {}, `Level ${level.jpi} at ${level.energy_text} MeV`),
    h("div.inline-form", {}, h("label.ctl", {}, h("span.ctl-label", {}, "E_λ (MeV)"), energy),
      h("label.ctl", {}, fixed, h("span.ctl-label", {}, "energy fixed")),
      h("label.ctl", {}, active, h("span.ctl-label", {}, "active")),
      button("Apply", "set_level with the changed values (Enter in the energy box)", applyLevel, "primary")),
    settingLine(null, level.fixed, "MeV"));

  const rows = level.channels.map((c) => ({ ...c, pair_label: pairLabel[c.pair] ?? c.pair,
    key: `${c.pair}/${c.l}/${c.s}` }));
  const columns = [
    { key: "channel", title: "ch", numeric: true }, { key: "pair_label", title: "pair" },
    { key: "l", title: "l", numeric: true }, { key: "s", title: "s", numeric: true },
    { key: "width", title: "width as written", numeric: true, sortable: false,
      doc: "the row's convention: a partial width Γ (eV), an ANC (fm^-1/2) or an amplitude (MeV^1/2)",
      render: (r) => widthEditor(r) },
    { key: "fixed", title: "fixed", sortable: false, render: (r) => fixedEditor(r) },
    { key: "amplitude", title: "amplitude", sortable: false, render: (r) => amplitudeEditor(r),
      doc: "whether the width column is a reduced-width amplitude (MeV^1/2) instead of a partial width or ANC (a "
        + "particle channel only; the value is not converted: type the width in the new convention)" },
    { key: "setting", title: "prior / limits", sortable: false, render: (r) => settingCell(r),
      doc: "a prior on the width (a term ((Γ - centre)/σ)² of the objective) or limits in a fit and the sampler; "
        + "a free channel only" },
    { key: "line", title: "line", numeric: true },
  ];
  const table = renderTable(document, { columns, rows, rowKey: (r) => r.key, withFilter: false, label: "channels",
    caption: "Channels: Enter in a width box sends the width; the header shows the new χ²" });
  box.appendChild(table.element);
  box.appendChild(settingForm);
  // "?setting=energy" or "?setting=<pair>,<l>,<s>" opens that parameter's prior and limits (a link, a screenshot)
  const asked = ctx.route?.params?.setting;
  if (asked) {
    const channel = asked === "energy" ? null : asked.split(",").map(Number);
    openSetting(channel, settingOf(settings, level.label, channel), channel ? "the width's convention" : "MeV");
  }
  box.appendChild(h("p.hint", {}, "A width is written as the file has it: Γ in eV for a particle channel above "
    + "threshold, an ANC in fm^-1/2 below it, an amplitude in MeV^1/2 where the row says so. The engine's "
    + "physical values and Wigner limits come with the Fit tab."));
  return box;

  function channelPath(r) {
    return `${path}/channels/${r.pair}/${r.l}/${r.s}`;
  }
  // a prior or limits on the energy or a channel's width: <parameterSettings> (set_parameter_setting), one form
  // below the channels opened by the energy's or a channel's "Prior / limits…"
  function settingLine(channel, isFixed, unit) {
    const row = settingOf(settings, level.label, channel);
    return h("div.inline-form.setting-line", {}, h("span.hint", {}, channel ? "" : `Energy: ${settingText(row)
      || "no prior, no limits"} `), button("Prior / limits…", isFixed ? "The engine applies a prior or limits to a "
      + "free parameter only: free it first" : "A prior on this parameter (a term of the objective) or limits in a "
      + "fit and the sampler (set_parameter_setting)", () => openSetting(channel, row, unit)));
  }
  function settingCell(r) {
    const channel = [r.pair, r.l, r.s];
    const row = settingOf(settings, level.label, channel);
    const open = button(row ? "edit" : "set…", r.fixed ? "The engine applies a prior or limits to a free parameter "
      + "only: free the channel first" : "A prior on this width or limits (set_parameter_setting)",
    () => openSetting(channel, row, "the width's convention"));
    return h("span", {}, settingText(row), row ? " " : "", open);
  }
  function openSetting(channel, row, unit) {
    settingForm.textContent = "";
    const what = channel ? `channel ${pairLabel[channel[0]] ?? channel[0]} l=${channel[1]} s=${channel[2]}`
      : "the energy";
    const box = (label, value) => h("input.num", { type: "text", size: 12, "aria-label": label,
      value: value === undefined || value === null ? "" : String(value) });
    const centre = box("prior centre", row?.nuisance ? row.value : "");
    const sigma = box("prior sigma", row?.nuisance ? row.error : "");
    const limited = row && !(row.lower === 0 && row.upper === 0);
    const lower = box("lower limit", limited ? row.lower : "");
    const upper = box("upper limit", limited ? row.upper : "");
    const send = (values) => {
      const body = settingRequest(level.label, channel, values);
      if (body.refused) { ctx.status(body.refused, "error"); return; }
      chi2Edit(ctx, `prior / limits of ${level.label} ${what}`,
        () => ctx.api.put("/api/project/parameter_settings", body), msg);
    };
    settingForm.append(h("h3", {}, `Prior and limits of ${what} (${unit})`),
      h("div.inline-form", {}, h("label.ctl", {}, h("span.ctl-label", {}, "prior centre"), centre),
        h("label.ctl", {}, h("span.ctl-label", {}, "σ"), sigma),
        h("label.ctl", {}, h("span.ctl-label", {}, "lower limit"), lower),
        h("label.ctl", {}, h("span.ctl-label", {}, "upper limit"), upper),
        button("Apply", "Write the prior and the limits (empty boxes: none)", () => send({ centre: centre.value,
          sigma: sigma.value, lower: lower.value, upper: upper.value }), "primary"),
        row ? button("Remove", "Remove the prior and the limits", () => send({})) : null),
      h("p.hint", {}, "A prior adds ((p − centre)/σ)² to the objective; limits bound the parameter in a fit and "
        + "the sampler. Leave a pair of boxes empty for none. The engine applies them to a free parameter only."));
  }
  function widthEditor(r) {
    const input = h("input.num", { type: "text", value: r.width_text, size: 14,
      "aria-label": `width of channel ${r.channel}` });
    let last = r.width_text;  // Enter and the change that follows it send once
    const send = () => {
      if (input.value.trim() === last) return;
      last = input.value.trim();
      chi2Edit(ctx, `width of ${level.label} channel ${r.channel}`,
        () => ctx.api.patch(channelPath(r), { width: Number(input.value) }), msg);
    };
    input.addEventListener("keydown", (e) => {
      if (e.key === "Enter") { e.preventDefault(); send(); } else if (e.key === "Escape") input.value = r.width_text;
    });
    input.addEventListener("change", send);
    return input;
  }
  function amplitudeEditor(r) {
    const flag = h("input", { type: "checkbox", "aria-label": `channel ${r.channel} width is an amplitude` });
    flag.checked = Boolean(r.amplitude);
    flag.addEventListener("change", () => chi2Edit(ctx, `amplitude flag of ${level.label} channel ${r.channel}`,
      () => ctx.api.patch(channelPath(r), { width_is_amplitude: flag.checked }), msg));
    return flag;
  }
  function fixedEditor(r) {
    const box2 = h("input", { type: "checkbox", "aria-label": `channel ${r.channel} fixed` });
    box2.checked = r.fixed;
    box2.addEventListener("change", () => chi2Edit(ctx, `fix ${level.label} channel ${r.channel}`,
      () => ctx.api.patch(channelPath(r), { fixed: box2.checked }), msg));
    return box2;
  }
}

/** A new level: J, parity, E and whether its energy is fixed; the channels its J^π allows with the project's pairs
 * (GET /api/project/allowed_channels: the engine's rule), each to be checked; Add sends add_level. */
export function newLevelPanel(ctx, pairLabel, msg) {
  const box = h("section.panel.new-level", { "aria-label": "new level" }, h("h2", {}, "New level"));
  const J = h("input.num", { type: "text", value: "0.5", size: 6, "aria-label": "J",
    title: "the spin as a number: 0.5, 1, 1.5, ..." });
  const parity = h("select", { "aria-label": "parity" }, h("option", { value: "1" }, "+"),
    h("option", { value: "-1" }, "−"));
  const energy = h("input.num", { type: "text", value: "", size: 12, "aria-label": "E_λ (MeV)" });
  const lmax = h("input.num", { type: "text", value: "2", size: 4, "aria-label": "l max" });
  const fixed = h("input", { type: "checkbox", "aria-label": "energy fixed" });
  fixed.checked = true;
  const list = h("div.allowed-channels");
  let allowed = [];
  let checks = [];
  const jValue = () => valueOf(J.value.trim());
  // the channels this J^π allows are asked as J, parity and l max are typed (the button asks again)
  let asked = 0;
  const offer = async () => {
    const mine = ++asked;
    list.textContent = "";
    try {
      const answer = await ctx.api.get(`/api/project/allowed_channels?J=${encodeURIComponent(jValue())}&parity=`
        + `${parity.value}&l_max=${encodeURIComponent(lmax.value.trim())}`, { quiet: true });
      if (mine !== asked) return;  // a newer question was asked meanwhile
      allowed = answer;
    } catch (e) {
      if (mine !== asked) return;
      allowed = [];
      ctx.showError(e, list);
      return;
    }
    checks = allowed.map((c) => {
      const input = h("input", { type: "checkbox", "aria-label": `channel pair ${c.pair} l ${c.l} s ${c.s}` });
      input.checked = true;
      return input;
    });
    list.append(h("p.hint", {}, `${allowed.length} channel(s), each added with width 0, fixed; with every one `
      + "checked, a J^π the project has gets its group's channels."),
    ...allowed.map((c, i) => h("label.choice", {}, checks[i], ` pair ${c.pair} (${pairLabel[c.pair] ?? ""}), `
      + `l ${c.l}, s ${c.s}`)));
  };
  const channels = button("Channels", "The channels this J^π has with the project's pairs (allowed_channels)", offer);
  channels.dataset.action = "allowed-channels";
  let timer = null;
  for (const input of [J, parity, lmax]) {
    input.addEventListener(input.tagName === "SELECT" ? "change" : "input", () => {
      clearTimeout(timer);
      timer = setTimeout(offer, 300);
    });
  }
  offer();
  const add = button("Add level", "add_level (POST /api/project/levels)", () => {
    const body = newLevelRequest({ J: jValue(), parity: parity.value, energy: energy.value, fixed: fixed.checked },
      allowed, checks.map((c) => c.checked));
    const before = new Set(ctx.levelLabels || []);
    // the new level opens with its channels beside the list (its label: the one the list did not have)
    ctx.edit("add a level", () => ctx.api.post("/api/project/levels", body), msg)
      .then(async () => {
        const after = await ctx.api.get("/api/project/levels");
        const label = after.map((l) => l.label).find((l) => !before.has(l));
        ctx.navigate(label ? `#/nucleus/level/${encodeURIComponent(label)}` : "#/nucleus/level");
      }, () => {});
  }, "primary");
  add.dataset.action = "add-level";
  box.append(h("div.inline-form", {}, h("label.ctl", {}, h("span.ctl-label", {}, "J"), J),
    h("label.ctl", {}, h("span.ctl-label", {}, "parity"), parity),
    h("label.ctl", {}, h("span.ctl-label", {}, "E_λ (MeV)"), energy),
    h("label.ctl", {}, fixed, h("span.ctl-label", {}, "energy fixed")),
    h("label.ctl", {}, h("span.ctl-label", {}, "l max"), lmax), channels), list, h("div.form-actions", {}, add));
  return box;
}

/** An edit after which the status line names the χ² before and after (both the engine's, as sent). */
export function chi2Edit(ctx, what, run, msg) {
  const before = ctx.evaluation?.chi2;
  return ctx.edit(what, run, msg).then(() => {
    const after = ctx.evaluation?.chi2;
    if (before !== undefined && after !== undefined) {
      ctx.status(`${what}: χ² ${formatNumber(before, 10)} → ${formatNumber(after, 10)} (Undo: Ctrl+Z)`, "ok");
    }
  }, () => {});
}
