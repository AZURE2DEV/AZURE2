# Where each feature of `dev` lands in this branch

For the owner and his agents. This branch does **not** have the features `dev` gained after the branch point
(`aab1de0d`), and by the PI's decision of 2026-11-15 they are not ported by the upgrade: the owner writes them
here. This page says, feature by feature, which named places of the new structure each one touches, so that the
work starts from a list of files and recipes and not from a search. The recipes are the sections of the Developer
Guide's page `docs/source/developer/extending.rst`, whose titles are quoted below; `CLAUDE.md` ("Recipes: where a
new feature goes") gives the same places with their consistency checks.

**What this page rests on.** A reading of `dev` at `3e502625` (65 commits beyond the branch point) against this
branch, done on 2026-11-14 and 2026-11-15; a measured plain merge of that commit (93 conflicted files, 21 merged
silently, at least eight of them wrongly: do not merge, re-write); runs where a line says "run". Everything else
was read, not run. `dev` has moved since (to `09e5c1a8` on 2026-11-15, with the branches `feature/rmfit`,
`feature/ec-unbound-final`, `feature/ec-unbound-broad`, `fix/formal-ec-zero-anc`, `gui/channel-radii-list`): those
commits were not read. The hours are estimates of an agent that did not do the work.

## 1. Already here (nothing to port but a test)

- **Coulomb functions of heavy ions just above threshold**: repaired in `physics/coulomb` (a listed local change;
  `FOR_THE_OWNER.md` section 3).
- **Coincident levels of one J^π in the Brune and Park transforms**:
  `physics/include/azure2/physics/LevelPairShift.h`, a gap of 1e-6 MeV with the mean S′ (`FOR_THE_OWNER.md`
  section 3, row 4; the gap is a question to the owner).
- **The scattering length**: `src/ScatteringLength.cpp`, `tests/scattering_length`.
- **Capture in a project with two compound nuclei** (dev's decay-index repair): right here (run: dev's three models
  converted, the merged project's two capture curves equal the standalone ones). dev's test `tests/two_nuclei_ec`
  is not ported: convert it with `AZURE2 --convert` and add it.
- **A final state of external capture in the formal parametrization**: this branch has its own version (the boundary
  condition at the lowest level of a group; the final states' amplitudes carried by `CompoundSystem`). Run: dev's
  four models agree across level orders and between the two bases. **Not run:** a final state that is not the
  lowest level of its group, `--no-transform` input, the ANC-sign and static-moment reports, a zero ANC. Port dev's
  `tests/formal_ec_final_state` (51 checks) and repair what it finds.
- **dev's new blocks in an old project file** (`<channelRadii>`, `<isospin>`, spin-operator data types): refused by
  the converter and by the reader with a message, never read past. Each refusal is lifted where its feature lands.
- **The flight path of a UDR record** is in metres (it was mis-scaled here; repaired).

## 2. Per-channel radii

*On dev:* a channel takes its pair's radius unless a line gives its own; every channel quantity reads the channel's
radius; the cache of external-capture integrals knows the entries.

*Here:* the model of the system already holds the radius per channel (`physics/include/azure2/physics/System.h`,
`radius_fm` of a channel). The engine does not use it: `PPair::GetChRad` is called 39 times in 28 files of `src/`,
`include/` and `api/` (counted on 2026-11-15). **This feature cuts across**; the way that keeps it small:

1. One accessor of the channel's radius, filled once at set-up from the system's channel ("A set-up fact of the
   compound system"), and the call sites changed to it. A refactor of its own: every output file identical.
2. A record kind `channel_radius` (pair, l, radius, optional s and J) of format 3 ("A record kind or a key of the
   project format": `experiment/src/Schema.cpp`, its build, the schema page); the converter reads the old block.
3. `AZUREAPI::SetRadius` (the session's radius override, kept across rebuilds) gains a channel form through the same
   accessor, so that a rebuild does not undo a channel's own radius.
4. `WorkingCopy::ECSignature` and the key of the external-capture table include the entries. The Coulomb table is a
   lookup in (η, ρ) and is not affected.
5. dev's `tests/channel_radii`, converted. Nothing for the Qt Pairs tab: the web page's forms are generated from the
   schema.

Estimate 10–20 h (the refactor 4–6). An unfinished attempt exists (section 8).

## 3. General spin observables

*On dev:* Re(c Tr(F A F† B))/Tr(F F†) for an operator token, on the plumbing of A_y, with its adjoint.

*Here* it is the two recipes "A data quantity" and "An observable of a reaction at a mesh energy":

- a quantity `spin_observable` in `experiment/src/Quantities.cpp` with a required operator field (a value type with
  its parser, so that `--validate` and the web page check the token), its tree in `experiment/src/Build.cpp`;
- a `SiteObservable` kind and its row in `src/SiteObservables.cpp` (the analyzing-power kind already holds the
  channel-spin amplitudes per slot, which is what the product-basis F needs);
- the formula called from `GenMatrixFunc::CalculateBareCrossSection`, its adjoint in `AMatrixFunc::PointTBar`; dev's
  token parser and operator matrices (`SpinObservable.h/.cpp`) as they are, better under `physics/` (pure
  angular-momentum algebra);
- the converter: the old numbered data type becomes the named quantity (dev renumbered the type once; format 3
  names quantities, so that trap is gone);
- the table of literature names (A_yy, ...) beside the quantity, so that pyazr computes nothing.

Test: dev's three projects converted, with their independent data; the gradient against finite differences.
Estimate 12–24 h. An unfinished attempt exists (section 8).

## 4. Time-of-flight resolution at low energy

*On dev:* four absolute floors written for charged particles broke the UDR kernel at thermal energies.

*Here* the same floors exist in other places, and the defect with them (run: zeros and a bias on dev's test project):

- `experiment/src/KernelUdr.cpp`, `UdrSupport`: `max(..., 0.001)`;
- `experiment/src/KernelDoppler.cpp`, `DopplerSupport`: the same floor; dev keeps it for a Doppler window without a
  UDR kernel, so the support must learn whether a UDR kernel is folded with it;
- `include/MeshPlanner.h`, `MeshZeroWidth()` = 1e-6 MeV absolute, used by `src/MeshPlanner.cpp` and
  `src/ShiftedFold.cpp`. **Shared by every kernel**: making it relative to the energy is a named step with its own
  table of all projects before and after (expected to move nothing in `tests/`, not measured);
- `experiment/src/KernelUdr.cpp`, `FoldUdr`: a one-node window returns 0, dev returns the point's value.

Test: dev's `tests/udr_low_energy` converted (`AZURE2 --convert` writes it with no notes, run). Estimate 5 h. An
unfinished attempt exists (section 8).

## 5. Isospin-labelled levels across compound nuclei

*On dev:* about 2,700 lines; a level ID spanning several compound nuclei is one isospin level, the member in a
reference nucleus free and the others derived, with shift and mixing parameters, priors and soft ties.

*Here*, six recipes cover most of it:

- the tie table is "A set-up fact of the compound system": a constant member of `CompoundSystem` built in
  `src/CompoundSystemBuild.cpp` from the level IDs; the boundary condition of tied groups in
  `CompoundSystem::CalcBoundaryConditions`;
- shifts and mixes are two rows of "A parameter kind" (`include/ParameterKind.h`), with their fields in
  `api/src/pybind_records.cpp`, so that pyazr's `Parameter` reads them with no Python code;
- the derived slots follow the one precedent of a derived parameter, the `gamma_decay` record
  (`src/GammaDecayMap.cpp`, `GradAccum::Scatter` in `include/AZUREGrad.h`, the fill in `src/LevelStateFill.cpp`);
- priors and soft ties are two rows of "A term of the objective" (`src/Objective.cpp`): they then reach the fit, the
  least-squares fitters, the sampler and pyazr's `Objective` at once (dev wrote each of the four by hand);
- the block's line kinds are record kinds of format 3; the converter reads the old block;
- `isospin.out` is an output file (`docs/source/reference/output_files.rst`).

**It cuts across in three places**, which is where the risk is: every fill path must honour the derived slots (from
a fit vector, from an observed vector, the transforms in and out, the parameter file, the sampler's walkers); in the
Brune basis the ties act on observed energies and amplitudes, which here are transformed per group
(`src/LevelTransforms.cpp`), so the tie must sit on the right side of that transform in each basis; one
`CompoundSystem` serves one project. Estimate 30–60 h, after the channel radii (the partner-radius check reads
them). Nothing was started.

## 6. Smaller items

- **EXFOR requests with the IAEA's user agent**: `pyazr/sources/nds.py`; dev's test imports the old module name.
- **dev's gradient-based sampler** (`pyazr.sampling`, written on pyazr 2): here `Objective` gives the log
  probability and its gradient in the fitter's coordinates, with the limits; a NUTS script on it is in the
  uncertainty skill. No package for it was installed.
- **One control for Formal / Brune / Brune-Park**: the keys stay `use_brune` and `use_park`; a control on the page's
  Project tab is "A route and a tab of the web UI".
- **dev's skill notes**: each belongs to the owner skill of its feature and is written when the feature lands.

## 7. The features announced for `dev`

rmfit, capture into an unbound final state and fitted secondary γ transitions are traced through this structure in
`FOR_THE_OWNER.md` section 6 (the second table). Two things were kept open for them:

- **Capture into an unbound final state**: *stage 1 (dev 2026-11-13, the narrow resonance) is written* on
  `feature/upgrade-5-dev-ec-unbound` (2026-11-16, the owner's session). The rule `UnboundFinalStates`
  (`experiment/src/FinalStates.cpp`) no longer refuses a final state above a threshold: it refuses one at a
  threshold, or too close above a charged pair's for its outgoing wave (eta > 100), with the same diagnostic
  `E_FINAL_STATE_UNBOUND`. The places that took a final state as bound (the table of the 2026-11-15 entry of
  `HANDOFF.md`) are changed as on dev: `src/ExternalCaptureFinalStates.cpp` (N_f and the ANC conversion),
  `src/ECIntegral.cpp` (`IntegrateUnbound`), `src/ECTable.cpp` (direct integration), `src/UnitExternalCapture.cpp` and
  `src/ChannelTables.cpp` (no photon below the final state: closed). `tests/ec_unbound_final` is dev's test.
  *Stages 2 to 4 are written too* (2026-11-16): the field `line_window` of a capture data set or extrapolation (dev's
  `linewindow W`; the converter reads the old keyword), the line-shape integration of a broad final state
  (`src/LineWindow.cpp`, one routine called by the energy pass; the T-matrix builders read its node through
  `GenMatrixFunc::ActiveLineNode`), its form for internal capture alone, and the line centre that follows the final
  level's current energy. `tests/ec_unbound_broad` is dev's test. *Not written yet:* the capture cascade (UPOS after
  a capture).
- **rmfit** is expected as its own top-level package beside `pyazr/`, on pyazr 3.0's names only (`Session`,
  `Objective`, `Project`, `Parameters`, `pyazr.sources`). It touches `tools/check_skills.py` (a second package's
  skills) and the CI. The ported copy of the `rmfit-campaign` skill here is a question to the owner
  (`FOR_THE_OWNER.md` section 8).

## 8. Three unfinished attempts (not in this branch)

Ports of sections 2, 3 and 4 were started by the upgrade's agents on 2026-11-15 and stopped the same morning by the
PI's decision. They exist as local branches with the PI (`feature/upgrade-5-channel-radii`,
`feature/upgrade-5-spin-observables`, `feature/upgrade-5-dev-resolution-capture`), partly uncommitted, **not
reviewed, not verified, not pushed**. They are not a base to trust; ask the PI if they are wanted as a starting
point.
