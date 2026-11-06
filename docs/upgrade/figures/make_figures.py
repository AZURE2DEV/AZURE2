"""Two figures for the AZURE2 structure page (draft 7): module dependencies and the evaluation chain."""
import html, pathlib
D = pathlib.Path(__file__).parent
def esc(s): return html.escape(s, quote=False)
MK = ('<defs><marker id="ah{k}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
      'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" class="ahp"/></marker></defs>')

def fig(name, W, H, label, body_fn):
    o = []
    k = name
    def box(x, y, w, h, cls, title, subs=()):
        o.append(f'<g class="bx {cls}"><rect x="{x}" y="{y}" width="{w}" height="{h}" rx="6"/>')
        o.append(f'<text x="{x+12}" y="{y+22}" class="bt">{esc(title)}</text>')
        for i, s in enumerate(subs):
            o.append(f'<text x="{x+12}" y="{y+40+15*i}" class="bs">{esc(s)}</text>')
        o.append('</g>')
    def ar(d, lab=None, lx=0, ly=0, cls=""):
        o.append(f'<path d="{d}" class="ar {cls}" marker-end="url(#ah{k})"/>')
        if lab: o.append(f'<text x="{lx}" y="{ly}" class="al">{esc(lab)}</text>')
    def txt(x, y, s, cls="bl"):
        o.append(f'<text x="{x}" y="{y}" class="{cls}">{esc(s)}</text>')
    body_fn(box, ar, txt)
    svg = (f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="{esc(label)}">'
           + MK.format(k=k) + "\n".join(o) + '</svg>')
    (D / f"{name}.svg").write_text(svg)

def deps(box, ar, txt):
    box(250, 20, 500, 62, "fe", "Front ends", ["azure2 CLI · pyazr · web UI · your scripts"])
    box(150, 130, 330, 70, "in", "inference/", ["parameters θ, ν · χ² and priors", "fitter"])
    box(520, 130, 330, 70, "io", "io/", ["project file · ENSDF, AME, EXFOR readers", "measured stopping data"])
    box(150, 248, 530, 70, "ex", "experiment/", ["measurements · materials and targets", "apparatus kernels K · energy grid E_n"])
    box(150, 366, 700, 70, "ph", "physics/", ["system: partitions, channels, Jπ blocks · channel functions",
                                             "collision matrix U · observables"])
    ar("M380,82 V128"); ar("M620,82 V128")
    ar("M315,200 V246"); ar("M600,200 V246")
    ar("M800,200 V364")
    ar("M415,318 V364")
    ar("M150,165 H118 V401 H148")
    txt(880, 160, "arrow = may use", "al"); txt(880, 178, "nothing points up", "al"); txt(880, 196, "physics uses nothing", "al")

def chain(box, ar, txt):
    X = [20, 250, 480, 710, 940]; W = 210
    def cx(i): return X[i] + W // 2
    txt(20, 24, "fixed for a fit stage · built once, then read-only")
    box(X[1], 34, W, 64, "ex", "Energy grid E_n", ["channel functions P, S, φ per (α, l)", "capture integrals per E_n"])
    box(X[2], 34, W, 64, "ph", "Angular coupling", ["coefficients per reaction", "Coulomb amplitude f_C(θ)"])
    box(X[3], 34, W, 64, "ex", "Kernels K", ["targets, stopping, straggling,", "resolution, Doppler, shifts"])
    box(X[4], 34, W, 64, "in", "Data d, covariance C", ["from the project"])
    for i in (1, 2, 3, 4): ar(f"M{cx(i)},98 V158")
    txt(20, 150, "every evaluation · depends only on θ, ν")
    box(X[0], 160, W, 72, "ph", "Level parameters", ["E_λ, γ_λc from θ", "Brune, Park; conventions"])
    box(X[1], 160, W, 72, "ph", "Collision matrix U(E_n)", ["level matrix per Jπ block,", "one solve per grid energy"])
    box(X[2], 160, W, 72, "ph", "Observables O_n", ["σ, Legendre coefficients,", "nuclear amplitudes"])
    box(X[3], 160, W, 72, "ex", "Apparatus", ["m = s · K(ν) · O", "per measurement"])
    box(X[4], 160, W, 72, "in", "χ² and priors", ["(m − d)ᵀ C⁻¹ (m − d)", "+ priors on θ, ν"])
    for i in range(1, 5): ar(f"M{X[i-1]+W},196 H{X[i]-2}")
    box(X[4], 272, W, 48, "in", "Fitter", ["Minuit2 or LM"])
    ar(f"M{cx(4)},232 V270")
    ar(f"M{X[4]},296 H{cx(0)} V234", "new θ", 300, 290)
    ar(f"M{cx(3)},296 V234", "new ν", cx(3) + 6, 262)

fig("deps", 1010, 450, "AZURE2 modules and which may use which", deps)
fig("chain", 1170, 332, "One AZURE2 evaluation: fixed tables feed a chain from parameters to chi-squared", chain)
print("ok")

if __name__ == "__main__":
    # Rebuild structure.html from the template and the two figures.
    page = (D / "page_template.html").read_text()
    page = page.replace("%%DEPS%%", (D / "deps.svg").read_text()).replace("%%CHAIN%%", (D / "chain.svg").read_text())
    (D / "structure.html").write_text("<!doctype html><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n" + page)
    print("structure.html written")
