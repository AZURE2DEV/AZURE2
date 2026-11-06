# shift_fn_groups fixture

tests/shift_fn's p+p model with every level energy and width fixed, and two segments in one
shift-function group (`shiftfngroup F`): the 40 deg segment declares
`shiftfn "a0*(x-1)+a1*(x-1)^2"` (a0, a1 free, penalty width 1), and both segments have their
own FREE constant shift (penalty width 1 MeV).  `xs_90_lin_mislabelled_p3.dat` is
tests/shift_fn's `xs_90_lin_mislabelled.dat` with every energy label moved by +0.003 MeV, so
the truth is a1(const, 40 deg) = -0.004, a2(const, 90 deg) = -0.007018, a0 = 0.006, a1 = 0
(the "gradcheck" model of tests/shift_fn_groups/check.sh).  For checks of evaluators that
must keep a shift-function group together (rmfit Sharded = Single).  Needs an engine with
`shiftfngroup` (feature/shift-fn-groups, 2026-10-22).
