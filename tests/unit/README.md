# Unit tests

Small standalone checks of individual routines, compiled directly against a
couple of source files rather than the whole engine. The physics regression
suite lives one level up in `tests/`.

Build and run one by hand, e.g. on macOS with the dependencies in a conda
prefix:

```bash
PREFIX=$CONDA_PREFIX          # or wherever gsl lives
c++ -O1 -std=c++17 -Iinclude -isystem $PREFIX/include \
    tests/unit/test_spherical_harmonics.cpp src/AngCoeff.cpp \
    -L$PREFIX/lib -lgsl -lgslcblas -Wl,-rpath,$PREFIX/lib \
    -o /tmp/test_spherical_harmonics
/tmp/test_spherical_harmonics
```

Each returns a non-zero exit code on failure.

- `test_spherical_harmonics.cpp` — `AngCoeff::SphericalHarmonic` against closed
  forms for l ≤ 2, the Condon-Shortley relation Y_l^-m = (-1)^m conj(Y_l^m) up
  to l = 4, and out-of-range behaviour. The Condon-Shortley phase is the one
  convention that silently flips the sign of an analyzing power, so it is
  checked rather than assumed.
