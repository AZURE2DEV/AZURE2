"""Names an evaluator recognises, as ENDF-102 (2023) defines them: the files (Table 4, p. 13), the reactions of
the charged-particle and neutron sublibraries (Appendix B, pp. 329-341), the sublibraries (Table 3, p. 7), the
resolved-resonance formats (p. 63) and the LRP flag (pp. 38-39).  Immutable tables."""

FILES = {1: "general information", 2: "resonance parameters", 3: "reaction cross sections",
         4: "angular distributions", 5: "energy distributions", 6: "energy-angle distributions",
         7: "thermal scattering law", 8: "radioactivity and fission-product yields",
         9: "multiplicities for radioactive nuclide production",
         10: "cross sections for radioactive nuclide production", 12: "photon production multiplicities",
         13: "photon production cross sections", 14: "photon angular distributions",
         15: "photon energy distributions", 23: "photo- or electro-atomic cross sections",
         26: "electro-atomic angle and energy distributions", 27: "atomic form factors", 28: "atomic relaxation",
         30: "covariances from parameter covariances", 31: "covariances of nu-bar",
         32: "covariances of resonance parameters", 33: "covariances of cross sections",
         34: "covariances of angular distributions", 35: "covariances of energy distributions",
         39: "covariances of radionuclide yields", 40: "covariances of radionuclide production"}

SUBLIBRARIES = {0: "photo-nuclear", 3: "photo-atomic", 4: "radioactive decay", 10: "incident neutron",
                10010: "incident proton", 10020: "incident deuteron", 10030: "incident triton",
                20030: "incident helion (3He)", 20040: "incident alpha"}

#: Symbols of the light particles by ZA (the projectile of a sublibrary, NSUB = 10 ZA, Table 3 p. 7).
PARTICLES = {1: "n", 1001: "p", 1002: "d", 1003: "t", 2003: "h", 2004: "a", 0: "g"}

#: The families of discrete-level reactions (Appendix B, pp. 332, 339-341): first MT, the emitted particle's ZA,
#: its symbol; MT = first + k for the residual's k-th level, first + 49 the continuum.
FAMILIES = ((50, 1, "n"), (600, 1001, "p"), (650, 1002, "d"), (700, 1003, "t"), (750, 2003, "h"),
            (800, 2004, "a"))

LRP = {-1: "no File 2", 0: "File 2 holds the scattering radius only",
       1: "File 2 parameters, their cross sections added to File 3",
       2: "File 2 parameters given, File 3 already holds their cross sections (derived files only, p. 42)"}

RESOLVED_FORMATS = {1: "SLBW", 2: "MLBW", 3: "Reich-Moore", 4: "Adler-Adler", 7: "R-Matrix Limited"}
RML_FORMULAS = {1: "SLBW", 2: "MLBW", 3: "Reich-Moore", 4: "full R-matrix"}


def reaction_name(mt):
    """The reaction of an MT number, as Appendix B writes it ("(z,p1)"), or "MT=<n>"."""
    simple = {1: "(z,total)", 2: "(z,z0) elastic", 4: "(z,n) total inelastic", 5: "(z,anything)",
              102: "(z,g) capture", 103: "(z,p) total", 104: "(z,d) total", 105: "(z,t) total", 106: "(z,h) total",
              107: "(z,a) total", 151: "resonance parameters", 451: "descriptive data and directory"}
    if mt in simple:
        return simple[mt]
    for first, _, symbol in FAMILIES:
        if first <= mt <= first + 49:
            k = mt - first
            return f"(z,{symbol}{'c' if k == 49 else k})" + (" continuum" if k == 49 else "")
    return f"MT={mt}"


def family_of(particle_za):
    """The first MT of the discrete-level family of an emitted particle (by its ZA), or None."""
    return next((first for first, za, _ in FAMILIES if za == particle_za), None)
