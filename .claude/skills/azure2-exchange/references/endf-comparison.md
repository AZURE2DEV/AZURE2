# Comparing an AZURE2 fit with an ENDF/B evaluation

## Getting the file

The IAEA mirror serves static directory listings that `curl` can walk (NNDC's own pages are rendered by
JavaScript): `https://www-nds.iaea.org/public/download-endf/<version>/<sublibrary>/`, for example
`ENDF-B-VIII.1/p/` for incident protons. Files are named `<projectile>_<Z>-<Elem>-<A>_<MAT>.zip`
(`p_006-C-12_0625.zip` for p + ¹²C); sublibraries are `n`, `p`, `d`, `t`, `he3`, `he4`, `g`, `decay`, ...
(`000-NSUB-index.htm` in the version root lists them). The mirror sometimes answers 402 to a first request; a plain
`curl` retry has succeeded. Fetching needs network access, which an engine run never has: fetch once, keep the file
with its version and MAT number.

## Parsing

`pip install endf` (pure Python, numpy only). Its `Material.interpret()` does not cover the charged-particle
sublibraries (`NotImplementedError: No class implemented for NSUB=10010` for protons in endf 0.1.12); parse at the
section level:

```python
import io, endf
mat = endf.Material("p_006-C-12_0625.dat")
print(mat.sections)                       # [(MF, MT), ...]
sec = endf.mf6.parse_mf6(io.StringIO(mat.section_text[6, 2]))
```

## Charged-particle elastic scattering (MF=6/MT=2, LAW=5)

Coulomb scattering has no finite angle-integrated cross section, so MF=3/MT=2 holds either a placeholder (1.0 at
every energy, LTP = 1 or 2) or, for LTP = 12 or 14, the nuclear-plus-interference integral σ_NI(E) in b, which can
be negative. The angular dependence is in MF=6/MT=2, LAW=5 (ENDF-102, section 6.2.7; pages 144–147 of the 2023
manual). Check `dist['distribution'][i]['LTP']` at each energy before choosing the reconstruction:

- LTP = 1: nuclear amplitude expansion, complex a_l(E) and real b_l(E) Legendre coefficients (ENDF-102 eqs. 6.13,
  6.14).
- LTP = 2: residual cross section as Legendre coefficients (eqs. 6.15–6.18); poor at forward angles and low energy,
  where the interference term is not well represented by a Legendre series.
- LTP = 12 / 14: tabulated (μ, P_NI), linear in μ (12) or in ln P_NI (14), normalised so that ∫P_NI dμ = 1 over the
  tabulated μ range (the representation of ENDF/B-VIII.1 p + ¹²C). Reconstruct
  σ(μ, E) = σ_R(μ, E) + σ_NI(E) P_NI(μ, E) inside the tabulated range and σ_R outside it, with σ_NI from MF=3/MT=2
  and P_NI linearly interpolated in μ (eqs. 6.19, 6.20).

Rutherford term (ENDF-102 eq. 6.9) in the c.m. frame, with T_c.m. = E_lab m_target/(m_beam + m_target) (MeV) and
e² = 1.43996 MeV fm:

```python
sigma_ruth_fm2 = (Z1 * Z2 * 1.43996 / (4 * T_cm))**2 / np.sin(theta_cm / 2)**4
sigma_ruth_b = sigma_ruth_fm2 / 100.0   # 1 b = 100 fm^2
```

## Energy coverage

`dist['distribution'][i]['E']` is the incident lab-energy grid in eV. Light-target proton sublibraries are built for
transport and dosimetry and often start at 1 MeV, missing the resolved-resonance region of a low-energy R-matrix
fit. Check the overlap with the fitted data before building the comparison.
