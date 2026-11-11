# Defining what a data set measures: totals, shared files, natural targets

A segment states which quantity the model is compared with. A segment that measures a sum over channels, or a
target of several isotopes, must say so in the `.azr`; otherwise the fit compares a partial cross section with a
total and the failure looks like missing physics. Behaviour: `docs/source/user_guide/segments.rst`,
`docs/source/user_guide/particle_pairs.rst`, `docs/source/examples/natural_neon.rst`.

## 1. A total is summed over every channel it measures

A neutron-total or total-yield segment is the sum of all channels open in its energy range. A segment defined with
too few exits looks like a model failure above the next threshold. Example: a 4π (α,n) total-yield set defined as
(α,n₀) only; adding the missing exits lowered χ² by 3.35 × 10⁴. A composite tail is
`1 0 <nComp> <ent> <exit> -999 ... 0` after the file name (operation 0 = sum; the segment's own exit included).

Two cheap tests before any refit:
1. Evaluate the corrected definition at the current parameters (CLI mode 1 with the `.sav`). The right definition
   usually wins already.
2. Plot two total measurements against each other across a threshold. Two sets that track each other with no step
   where a channel opens measure the same total quantity.

## 2. Two lines sharing one data file are not summed

Older archive files define a total as several plain lines naming the same data file under different exit pairs,
with no tail tokens. The builds that wrote them summed such lines (outside this repository's history); current
builds read them as independent segments, so every line but the first compares one partial cross section with the
total. No warning is printed. In one case the (n,α₀) line alone contributed 5.13 × 10⁶ of a 5.95 × 10⁶ χ² where
the stored output said 8.4 × 10⁵.

- **The tell is in the old output**: a `chiSquared.out` with exactly `0` for every line of a same-file group but the
  last, which carries the whole χ².
- **Find them** (prints `file[ent>exit ...]` for every data file on more than one active line):

  ```
  awk '/<segmentsData>/{f=1;next}/<\/segmentsData>/{f=0}f&&$1==1{n=split($0,t," ");for(i=1;i<=n;i++)if(t[i]~/\.dat$/){c[t[i]]++;e[t[i]]=e[t[i]]" "$2">"$3}}END{for(k in c)if(c[k]>1)print k"["e[k]" ]"}' model.azr
  ```

  The same file under one exit (energy windows, one line per angle) is legitimate; under different exits it is this
  trap.
- **Fix**: one composite line; set the others inactive rather than deleting them (deletion renumbers keys). An
  old-format line has no shift tokens and the tail is parsed only after them, so insert `0 0 0` before the file name:
  `... 1 0 0  0 0 0  data/ntotal.dat 1 0 2 2 1 -999 2 3 -999 0` (own exit n₀ plus components α₀ and n₁).
- **Check in two steps**: first the sum the old file defined, which must reproduce the old number (n₀ + α₀ gave
  698,808 against 203.497 × 3434 = 698,809 stored); then add the channels the old definition left out (n₁, open
  above E_n = 1.716 MeV: 723,978).
- A stored `output/` in such a directory cannot be reproduced from the unedited `.azr` with a current build; say so
  before comparing χ² across the edit.

### Neutron totals: prefer `total 1`

`... data/ntot.dat 0 0 total 1` on an n+X → n+X angle-integrated line is the optical-theorem total
σ_tot = (2π/k²) Σ g_J (1 − Re U_cc), i.e. 2 Re T_cc, from one T-matrix. It equals the sum over every particle exit of
the compound nucleus (1e-11 relative on a 9-exit case) and cannot leave a channel out. Capture is not included
(photon channels carry no width in the level matrix): keep a capture pair as a SUM component if it matters. In a
natural-target SUM every elastic component becomes its nucleus's total. Refused for charged or identical pairs, other
exits, differential data, `--use-rmc`. Engines that predate the keyword ignore it and compute the elastic cross
section alone: check the build before trusting a `total 1` file. All terms of a SUM, RATIO or total-capture segment
share one level-matrix solution per point; `AZURE_NO_COMPONENT_MAPPING=1` restores the per-term route (same numbers).

## 3. Several compound nuclei: natural targets and samples

- **Identity is the compound nucleus, abundance is the target isotope's.** The nucleus of a pair is inferred from its
  columns (particle pair Z = z₁ + z₂, A = round(m₁ + m₂); γ pair the heavy particle; β pair the daughter). n + ²²Ne
  lines are ²³Ne lines, and the parameter for natural neon's 9.25 % of ²²Ne is `abundance_10_23`.
- **A natural data set is one segment line**: a SUM composite whose components are the members' terms (one per
  nucleus) plus the `sample 1` keyword (`natural 1` in older files). Each term is multiplied by the abundance of its
  entrance pair's compound nucleus, taken from the `<abundances>` block after `</levels>` (`Z A value vary sigma`;
  a parameter exists only when the block does, emitted after the level parameters and before the norms) and
  `abundanceSumSigma=s` in `<config>`.
- **Samples** (targets of different composition): `<name>: Z A value vary sigma` lines form a named sample,
  `<name>: sum s` its own sum-to-one width (0 = off, e.g. a stoichiometric Al₂O₃ sample with weights 2 and 3);
  `sample <name>` picks it on the segment line; parameters `abundance_<name>_Z_A`. Refused on RATIO, phase shifts
  (2), analysing powers (7), test angular-coefficient lines (3), in a project without the block, and for an
  undefined sample name. Output `AZUREOut_sample_seg=<key>.out/.extrap` (total plus one column per nucleus).
- **Build such a project** (format 3) with `pyazr.merge_projects(projects, out, isotopes={"20Ne": 0.9048, ...},
  shared=[{0: "vaughn", 1: "vaughn"}])` rather than by hand: the engine merges the `<levels>` rows (pair keys and
  level IDs offset), every record is kept (colliding names renamed), and `isotopes` adds a `target` record with one
  `isotope` field per target nucleus; member 0's `<parameterSettings>` and `<mcmc>` are kept, the others' dropped.
  Edit the isotopes as record fields (`Project.set("natural", isotope=[...], isotope_sum=...)`); read
  `parameters.abundances` and `penalties()["abundance"]`.
- **Cross-check a merge against its members**: each nucleus alone (other nuclei's components removed, abundances
  1/0/0) reproduces its standalone project's χ² to 1e-10 (`tests/natural_target/check.sh`). The real-data reference
  model is `tests/natural_neon` (three nuclei, EXFOR data).
- rmfit: keys of kind `abundance`, prior rows from the block; the sum-to-one term is not a row, so fix the abundances
  or accept priors only.
- The `<targetInt>` rules for extrapolation apply per segment as usual (azure2-calculate).
