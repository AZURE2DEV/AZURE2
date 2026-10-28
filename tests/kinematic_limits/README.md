# kinematic_limits

A lab angle the forward-branch kinematics cannot reach is an error, not a NaN
(audit 2026-10-25, 01_kinematics_frames C1-C4, checklist item 31).

`EPoint::ConvertLabAngle` maps a lab angle to the c.m. on the forward branch
only. Where two c.m. angles reach one lab angle (a reaction near the threshold
of an endothermic exit, or elastic scattering with a projectile heavier than
the target) a lab angle past the kinematic maximum gave NaN: the reaction
branch printed one line and asserted (compiled out in release builds), the
elastic branch said nothing, and a point below the exit threshold went NaN
through the Marion-Young map with no message. `chiSquared.out` and the
extrapolation files then carried `nan` and the run exited 0; in a fit the NaN
chi-squared silently spoiled everything. Now the converters return false with a
reason, the four readers (data segments, test segments, SUM/RATIO components,
and the energy-shift update during a fit, which warns and keeps the last valid
angle) report the point, and the run stops before any output. A UPOS segment
given as lab differential draws a warning that the conversion uses the
kinematics of particle b.

## Checks (`check.sh path/to/AZURE2`)

- (a) `tests/7Li_p_ay` plus a lab-differential test segment 7Li(p,n)7Be at
  60 deg for E_p = 1.89-1.90 MeV (theta_max 22-39 deg), mode 3: refused, no
  extrapolation file;
- (b) the same as a data segment, mode 1: refused, no `chiSquared.out`;
- (c) `tests/13N` with pair 1 written as 12C + p (theta_max 4.8 deg), segment
  3 at 84.3 deg: refused;
- (d) control: the (p,n) segment at 10 deg runs with finite numbers;
- (e) a lab-differential point below the (p,n) threshold (E_p = 1.88 MeV,
  threshold 1.885): refused ("below the threshold of the exit pair").

dev 265c259: (a), (b), (c) and (e) exit 0 with `nan` in the output files.
