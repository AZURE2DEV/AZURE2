# The data structure of the engine: what was measured, the data set, the working copy

This note replaces the plan of many small steps (`data_model.md` S7, S12 and Phase 6's share) by the structure the
code should have. `EData`, `ESegment` and `EPoint` mix what was measured, the settings of a data set, kinematics,
energy-dependent tables and apparatus state in one object graph that every working copy clones (`data_model.md`
§1.1). The structure below separates them by who writes them and when. The physics functions (channel tables,
folds, target effects, external capture, THM, covariance, χ²) keep their formulas and the order of their
floating-point operations; what changes is who owns the numbers they read and write. Every output byte of every
project stays as it is.

## 1. Three kinds of state, three types

| state | type | built | written | shared |
|---|---|---|---|---|
| what was measured | `experiment::Measurement` | once at set-up, from the spec and the data file | never after set-up | `shared_ptr<const>` by every working copy |
| the data set: measurements, settings, fixed apparatus, covariance, parameter layout | `DataSet` | once at set-up | never after set-up | const pointer from every working copy |
| what depends on the parameters | `WorkingCopy` | copied from the set-up copy | the single-threaded hand-over | none: a copy is a plain copy |
| what an evaluation produces | `Results` (`include/Results.h`, unchanged) | by the caller | the point pass, each slot by one task | none |

### 1.1 `experiment::Measurement` (`experiment/include/azure2/experiment/Measurement.h`, no engine headers)

```
struct MeasurementRows {                 // one entry per row kept (inside the data set's window), file order
  std::vector<double> energy;            // MeV: E_lab of the beam, or the lab energy of the light decay particle
  std::vector<double> angle;             // deg lab, deg c.m. (c.m. observables), or the Legendre order k
  std::vector<double> value, error;      // the quantity's engine unit (b, b/sr, 1, deg); error after a covariance file
  std::vector<double> binLow, binHigh;   // MeV lab, NaN without columns 5-6
  std::vector<int> record;               // 1-based data row in the file (covariance files index by it)
};
struct Measurement {
  int key;                               // data set number in the file (1-based), the label of its outputs
  int entranceKey, exitKey;              // <pairs> keys of the base term (exit -1: total capture)
  bool extrapolation;                    // rows from a grid (no file; value 0, error 0.1)
  AngleColumn angleColumn;               // Lab | CM | LegendreOrder | Unused
  bool labDifferential;                  // the lab -> c.m. angle and Jacobian apply
  bool convertValue;                     // value and error carry dOmega_lab/dOmega_cm (false for a ratio)
  bool coefficientAbsolute;              // LegendreOrder rows: A_k (b/sr), not a_k
  MeasurementRows rows;
  Provenance provenance;                 // file as opened, rows read, rows kept, kept E_lab range (MeV)
  std::shared_ptr<const Prediction> prediction;   // the expression tree (experiment/Prediction.h)
};
```
`ReadMeasurementRows` (experiment/) reads the file and applies the window (E_lab in [minE, maxE], and θ_lab in
[minA, maxA] for differential data), the rule of `ESegment::IsInSegment`. The engine and `--describe` call it; the
engine then replaces errors from a covariance file (sqrt C_ii) before the measurement is frozen.

### 1.2 `experiment::MeasurementFrame` (`experiment/Frames.h`) and the one conversion

```
struct FrameRow { double eLab, eCM, eX, thetaCM, jacobian, yCM, dyCM, binLowCM, binHighCM; };
struct MeasurementFrame { std::vector<FrameRow> rows; };
FrameRow ConvertRow(const FrameRule &, double eLab, double thetaLab, row values, const FrameRow *atLoad, AngleLimit &);
```
`FrameRule` holds the partitions of the term (entrance, exit, exit kind: none / particle / photon), elastic,
β-delayed, relativistic, `labDifferential`, `convertValue`. `ConvertRow` is the sequence of
`ESegment::ConvertPointToCM` and the `EPoint::Convert*` functions, operation for operation. At load `atLoad` is null
and a lab angle beyond the kinematic maximum stops the run; after a shift the angle and Jacobian fall back to
`atLoad` (H4 of `data_model.md` §1.3: the load-time frame, not an earlier evaluation's). The frame at the nominal
shift is computed once at set-up and shared; the current frame is the working copy's.

### 1.3 `DataSet` (immutable)
The measurements; per measurement the settings each data set had as `ESegment` members (norm, shift and thickness
specification and groups, THM attachments, target-effect number, keys of its leaves); the fixed part of the target
effects (lines, kernels, stopping); covariance groups; the parameter layout (offsets); the frames at the nominal
shift. Shared by every working copy through a const pointer.

### 1.4 `WorkingCopy` (small, flat)
- per measurement: its current `MeasurementFrame`, norm, shift state (constant, sqrt, function values; last applied),
  thickness, channel tables of the last three shift states;
- the calculation units, **one vector** of plain records: a row of a measurement for one leaf, or a sub-point; each
  names its measurement, leaf and row by index, holds its kinematics (E_lab, E_c.m., E_x, θ), its tables through the
  copy-on-write `SharedTable`, its sub-points as a range of the same vector, its apparatus values per ν (density
  thickness, U_e), its result offset;
- the calculation plan by unit index; the apparatus values that depend on ν (densities, Doppler T, U_e).
No pointer from one unit to another or to a parent; a copy is a plain copy.

### 1.5 Components are leaves
A SUM, RATIO, total-capture or sample term is a leaf of the measurement's prediction: entrance, exit, scaling, fixed
angle, total / yield-total flags, and a key from the data set's counter (the old −1000 counter: the keys stay, they
are part of the cache keys and of `intEC.dat`). Its units are units of the base measurement's rows with another
exit; they hold no data. The order of the units (base first, leaves in term order) keeps the `intEC.dat` layout.

## 2. Where every member goes (groups of `data_model.md` §1.1)

| group | members now | goes to |
|---|---|---|
| input datum | `original_energy_`, `lab_angle_` (as read), `lab_crosssection_`, `lab_dcrosssection_`, `bin_low_lab_`, `bin_high_lab_`, `fileRecord_`, `coefficient_order_`, `coefficient_absolute_` | `Measurement::rows` (read by row) |
| segment settings in every point | keys, `is_*` flags, `j_value_`, `l_value_`, UPOS, `thm_*` pointers, `crossSectionComponent_`, `max_ang_dist_order_` | the measurement / leaf, named by index (milestone 2) |
| kinematics from the datum | `cm_energy_`, `excitation_energy_`, `cm_angle_`, bins c.m. | the frame (rows) and the unit (E and θ it is calculated at) |
| | `cm_crosssection_`, `cm_dcrosssection_`, `crossSectionKinFactor_`, `loadFrame_` | the frame only |
| energy-dependent tables | `edep_`, `legendreP_`, `thm_ps_`, `geofactor_`, `sfactorconv_`, Coulomb, EC | the unit (copy-on-write, unchanged) |
| apparatus | `targetEffectNum_`, stopping, thickness, grid spec, UDR, blend, fallback, screening | fixed part: data set; per-ν values: the unit |
| topology | `is_mapped_`, `energy_map_`, `own_e_dependent_`, `plan_task_`, `integrationPoints_`, `parentData_` | the plan (by index); sub-points: a range of the unit vector |
| results | (already in `Results`) | unchanged |
| `ESegment` settings | keys, ranges, flags, norm/shift/thickness specification, groups, THM | the data set |
| `ESegment` state | `dataNorm_`, `energyShift*`, `shiftFnValues_`, `thickness_`, `lastEnergyShift*`, `shiftTables_` | the working copy, per measurement |
| `ESegment` components | `componentSegments_` (pointers), `componentScaling_` | leaves of the prediction (§1.5) |
| `EData` | segments, component deque, target effects, covariance, offsets, THM groups, plan, grids cache | data set (fixed) and working copy (state) |

## 3. The files and the picture (milestone 3, done)

```
  io::Document --specs--> set-up (WorkingCopy::Fill / MakePoints / Initialize, writing the draft)
                                |
  experiment::Measurement  <----+   rows, provenance, prediction           (immutable, shared_ptr<const>)
             |
             v
  DataSet  (include/DataSet.h)      segments (Segment: settings, start values, priors, terms, frame at load),
             |                      target effects, leaves, PointSetUp records (grid layout, UDR kernel),
             |                      covariance groups, THM experiments, parameter layout, prior covariance
             |                      -- frozen at the end of Initialize; shared_ptr<const DataSet> in every copy
             v
  WorkingCopy (include/WorkingCopy.h)  per segment a SegmentState (norm, shifts, shift function, thickness; per term
             |                         the unit range, the present frame, the shift tables), TargetEffectValues
             |                         (U_e, density, Doppler T), THM coherent values, the calculation plan, one
             |                         vector of Units (include/Unit.h) -- copied by a plain copy (Clone)
             v
  Results  (include/Results.h)      leaf records, model values, objective values -- the caller's
```

`include/`: `DataSet.h`, `SegmentState.h`, `WorkingCopy.h`, `Segment.h`, `Unit.h`, `UnitTables.h`, `UnitRange.h`,
`SubPointAccess.h`, `TermSettings.h`.  `src/`: `DataSet.cpp` (lookups), `WorkingCopy.cpp` (set-up from the specs with
its refusals, Initialize, the copy, the result layout, the mapping, the parallel point pass), `WorkingCopyUnits.cpp`,
`HandOver.cpp`, `DataParameters.cpp`, `TargetSetup.cpp`, `ExternalCaptureSetup.cpp`, `DataCovarianceGroups.cpp`
(`DataCovariance.cpp` is the file reader of another branch), `ThmNorm.cpp`, `ThmReports.cpp`, `DataOutputFiles.cpp`,
`DataCheckFiles.cpp`, `Segment.cpp`, `SegmentParameters.cpp`, `ShiftFrames.cpp`, `Unit.cpp`, `PointPass.cpp`,
`ChannelTables.cpp`, `TargetFold.cpp`, `UnitLegendre.cpp`, `UnitExternalCapture.cpp`.  No `EData*`, `EPoint*`,
`ESegment*` file or class (consistency check 31).

Shared and written after the set-up, on purpose, and nothing else: a covariance group's frozen relvector model (under
the group's mutex; DataSet.h) and the node tables of re-placed sub-points (`NodeTableCache.h`, under its shards'
mutexes).  A setting changed after the set-up (a prior covariance from pyazr, a parameter layout laid out again) makes
a new DataSet (`WorkingCopy::ReviseDataSet`).

## 4. Order of the commits

Milestone 1 (what was measured exists once; done, commits 0439701a and 9db6b59a):
1. `Measurement`, `ReadMeasurementRows`, `ConvertRow` in `experiment/` with unit tests by hand.
2. Segments hold a `shared_ptr<const Measurement>`; points read the datum by row; the input members go.
3. The frame: `ESegment::Fill`, the shift step and the component load call one function; y_c.m., δy_c.m. and the
   Jacobian live in the frame; `loadFrame_` goes (the nominal frame takes its place).
4. Components as leaves of the base segment: the component deque and its copies go.  In milestone 1 a leaf lives in
   its segment (`ESegment::Leaf`: pairs, key, fixed angle, scale, role, target effect, units); the segment's own
   points are term 0, a leaf term k, and the shift step, the frame, the grids and the target set-up take
   (segment, term) (`TermRef`); `EData::LeafRefs` lists the leaves in creation order (plan, layout, intEC file).
5. `--describe` shows rows read, kept and the kept range.
Milestone 2 (the working copy is flat), done: the settings of a term are one shared immutable record
(`TermSettings`), which a unit names by a plain `const` pointer (the term owns the shared pointer); a point has no
rows pointer (its row index, Legendre order, nominal energy and angle only) and no pointer to its data set (the
folds, `Calculate`, the plan lookup, the model pass and the gradient take the data set as an argument); every point
of every term and every sub-point is a unit of one vector of the working copy (`EData::Units`,
`src/WorkingCopyUnits.cpp`), a term's points and a point's sub-points are ranges of it, and the accessors
(`GetPoints`, `GetPoint`, `TermPoints`, `GetSubPoint`, `GetSubPoints`) are views that take the data set; sub-points
given at set-up and re-placed grids are placed single-threaded as new ranges and the vector is compacted (points
first, in data order; they keep their indices after the set-up); the plan names units by index (`PointRef` is an
int); `EData::Clone` is the copy constructor, after which the copy's points drop the set-up's external-capture
amplitudes (as before; they recompute them at their energies); the S factor has one formula.
Milestone 3 (the old names are gone), done: `DataSet` (immutable, shared) and `WorkingCopy` apart, the
segment's state in `SegmentState`, the target-effect values and THM coherent values in the copy, the unit 608 -> 480
bytes (its grid layout and UDR kernel are data-set records, the angle-limit message is returned, not kept), the data
iterator gone, the files and classes renamed (`WorkingCopy`, `Segment`, `Unit`).  Still in the unit: the external-capture
amplitudes of the set-up, which a copy drops in `Clone` (the next pathway branch computes them at set-up).

## 5. Bits

Moving a value keeps its double; `ConvertRow` is the old code with the same operations in the same order; the units
are created and calculated in the old order, so the plan, the mapping and `intEC.dat` keep their order. Each commit
is checked against the reference outputs on the projects closest to it; each milestone on all of them.
