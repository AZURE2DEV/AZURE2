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
