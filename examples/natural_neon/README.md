# Neutrons on natural neon: a worked example

One AZURE2 project with three compound nuclei, n+²⁰Ne → ²¹Ne, n+²²Ne → ²³Ne and n+²¹Ne → ²²Ne. It is fitted to the
eleven measured neon neutron total cross sections below 2.5 MeV, taken on natural and enriched samples (6,328 points).
The walkthrough is the documentation page *Worked Example: Neutrons on Natural Neon*
(`docs/source/examples/natural_neon.rst`).

| File | Contents |
|---|---|
| `natural_neon.azr` | the project; levels at their fitted values, normalizations at their nominal value 1 |
| `natural_neon_fit.sav` | the fitted parameters, including the normalizations |
| `data/` | the data sets (lab MeV, b, b) |
| `plot_natural_neon.py` | the two figures of the documentation page |

Calculate with the fitted parameters:

```
mkdir -p output checks
printf '1\nnatural_neon_fit.sav\n\n7\n' | AZURE2 --no-gui --no-readline natural_neon.azr
```

Data sources:
- Junghans *et al.* (2020), EXFOR 23755005.
- Vaughn *et al.* (1960), EXFOR 11106003.
- Almeida and Käppeler (1983), EXFOR 21808.
- Sikkema (1958), EXFOR 20203 and 20261.
- Hinterberger *et al.* (KfK), digitized by R. J. deBoer from the original run printouts. It is not in EXFOR.

`tests/natural_neon` is a cut-down copy of this model, used as a regression test.
