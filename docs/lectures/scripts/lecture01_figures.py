#!/usr/bin/env python3
"""Figures of Lecture 1 (the external region: Coulomb and Whittaker functions).

Every function is computed with mpmath (coulombf, coulombg, whitw), independently of AZURE2, so the figures show the
mathematics the engine implements, not the engine's output.  Systems: p + 12C (open channel) and p + 16O -> 17F
(bound states 5/2+ at B = 0.600 MeV, l = 2, and 1/2+ at B = 0.105 MeV, l = 0).

Usage: python3 lecture01_figures.py OUTDIR
"""
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np

sys.path.insert(0, os.path.expanduser("~/.claude/skills/prc-figures/scripts"))
import prcstyle as ps  # noqa: E402

mp.mp.dps = 20
HBARC = 197.3269804      # MeV fm
UC2 = 931.49410242       # MeV / u
ALPHA = 1 / 137.035999084
M_H, M_C12, M_O16 = 1.00782503, 12.0, 15.99491462   # atomic masses (u)


class Pair:
    """A two-body channel: Z1 Z2, reduced mass mu (u), channel radius a (fm)."""

    def __init__(self, z1z2, m1, m2, a):
        self.z1z2, self.mu, self.a = z1z2, m1 * m2 / (m1 + m2), a
        self.muc2 = self.mu * UC2

    def eta(self, e):      # Sommerfeld parameter at c.m. energy e > 0 (MeV) or binding e > 0
        return self.z1z2 * ALPHA * math.sqrt(self.muc2 / (2.0 * e))

    def k(self, e):        # wave number (fm^-1): k at energy e, or kappa at binding e
        return math.sqrt(2.0 * self.muc2 * e) / HBARC


PC = Pair(6, M_H, M_C12, 1.25 * (1 + 12 ** (1 / 3)))     # a = 4.11 fm
PO = Pair(8, M_H, M_O16, 1.25 * (1 + 16 ** (1 / 3)))     # a = 4.40 fm


def coul(l, eta, rho):
    """F, F', G, G' (primes d/drho)."""
    F = mp.coulombf(l, eta, rho)
    G = mp.coulombg(l, eta, rho)
    Fp = mp.diff(lambda r: mp.coulombf(l, eta, r), rho)
    Gp = mp.diff(lambda r: mp.coulombg(l, eta, r), rho)
    return float(F), float(Fp), float(G), float(Gp)


def whit(l, eta_b, rho):
    """W_l(rho) = W_{-eta_b, l+1/2}(2 rho) and rho W'/W (DLMF 13.15.23)."""
    W = mp.whitw(-eta_b, l + 0.5, 2 * rho)
    W1 = mp.whitw(1 - eta_b, l + 0.5, 2 * rho)
    return W, float(rho + eta_b - W1 / W)


def shift(pair, l, e):
    """The shift function S_l at c.m. energy e (MeV, relative to threshold): Coulomb functions above threshold,
    the Whittaker function below, the zero-energy limit at e = 0."""
    if e > 0:
        F, Fp, G, Gp = coul(l, pair.eta(e), pair.k(e) * pair.a)
        rho = pair.k(e) * pair.a
        return rho * (F * Fp + G * Gp) / (F * F + G * G)
    if e < 0:
        return whit(l, pair.eta(-e), pair.k(-e) * pair.a)[1]
    c = pair.z1z2 * ALPHA * pair.muc2 / HBARC          # eta k (fm^-1)
    x = math.sqrt(8 * c * pair.a)
    return float(-l - 0.5 * x * mp.besselk(2 * l, x) / mp.besselk(2 * l + 1, x))


def fig_coulomb(out):
    fig, axs = plt.subplots(1, 2, figsize=(ps.COL_2, 2.6), sharey=True)
    rho = np.linspace(0.02, 20, 600)
    for ax, eta, tag in zip(axs, (0.0, 2.0), ("(a)", "(b)")):
        F = [float(mp.coulombf(0, eta, r)) for r in rho]
        G = [float(mp.coulombg(0, eta, r)) for r in rho]
        ax.plot(rho, F, color=ps.C["gs"], ls="-", label="$F_0$ (regular)")
        ax.plot(rho, G, color=ps.C["sum"], ls="--", label="$G_0$ (irregular)")
        if eta > 0:
            tp = 2 * eta                                  # outer turning point for l = 0: rho_tp = 2 eta
            ax.axvline(tp, color=ps.C["bkg"], ls=":", lw=0.8)
            ax.text(tp + 0.3, -2.1, r"$\rho_{\rm tp} = 2\eta$", fontsize=8,
                    bbox=dict(facecolor="white", edgecolor="none", pad=1))
        ax.axhline(0, color="black", lw=0.5)
        ax.set_xlim(0, 20)
        ax.set_ylim(-2.5, 2.5)
        ax.set_xlabel(r"$\rho = kr$")
        ps.minor(ax)
        ps.label(ax, f"{tag} $l = 0$, $\\eta = {eta:g}$")
    axs[0].set_ylabel("Coulomb functions")
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(os.path.join(out, "fig1_coulomb_functions.png"), dpi=200)
    plt.close(fig)


def fig_penetrability(out):
    fig, ax = plt.subplots(figsize=(ps.COL_1, 2.8))
    E = np.geomspace(0.05, 4.0, 120)
    styles = [("-", ps.C["gs"]), ("--", ps.C["sum"]), ("-.", ps.C["casc"]), (":", ps.C["alt"])]
    for l, (ls, col) in enumerate(styles):
        P = []
        for e in E:
            F, _, G, _ = coul(l, PC.eta(e), PC.k(e) * PC.a)
            P.append(PC.k(e) * PC.a / (F * F + G * G))
        ax.plot(E, P, ls=ls, color=col, label=f"$l = {l}$")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel(r"$E$ (MeV, c.m.)")
    ax.set_ylabel(r"$P_l(E)$ at $a = 4.11$ fm")
    ax.set_ylim(1e-12, 3)
    ps.minor(ax)
    ps.label(ax, r"p + $^{12}$C", loc="lower right")
    ps.legend_above(ax, ncol=4)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig2_penetrability.png"), dpi=200)
    plt.close(fig)


def fig_shift(out):
    fig, ax = plt.subplots(figsize=(ps.COL_1, 2.9))
    styles = [("-", ps.C["gs"]), ("--", ps.C["sum"]), ("-.", ps.C["casc"])]
    E = np.concatenate([np.linspace(-3.0, -0.002, 70), np.linspace(0.002, 3.0, 70)])
    for l, (ls, col) in enumerate(styles):
        S = [shift(PO, l, e) for e in E]
        ax.plot(E, S, ls=ls, color=col, label=f"$l = {l}$")
        ax.plot([0.0], [shift(PO, l, 0.0)], "o", ms=3, color=col)
    ax.axvline(0, color="black", lw=0.5)
    for b, txt in ((0.600, r"$\frac{5}{2}^+$"), (0.105, r"$\frac{1}{2}^+$")):
        ax.axvline(-b, color=ps.C["bkg"], ls=":", lw=0.8)
        ax.text(-b - 0.05, 0.25, txt, fontsize=8, ha="right", bbox=dict(facecolor="white", edgecolor="none", pad=1))
    ax.set_xlabel(r"$E - E_{\rm thr}$ (MeV, c.m.)")
    ax.set_ylabel(r"$S_l(E)$ at $a = 4.40$ fm")
    ax.set_xlim(-3, 3)
    ax.set_ylim(-3.6, 0.6)
    ps.minor(ax)
    ps.label(ax, r"p + $^{16}$O", loc="lower right")
    ps.legend_above(ax, ncol=3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig3_shift_function.png"), dpi=200)
    plt.close(fig)


# 17F bound states: binding (MeV), l, ANC C (fm^-1/2; C^2 = 1.08 and 6490 fm^-1, Gagliardi et al., PRC 59, 1149 (1999))
STATES = [(r"$^{17}$F $\frac{5}{2}^+$ ($B = 0.600$ MeV, $l = 2$)", 0.600, 2, math.sqrt(1.08), ps.C["sum"], "--"),
          (r"$^{17}$F $\frac{1}{2}^+$ ($B = 0.105$ MeV, $l = 0$)", 0.105, 0, math.sqrt(6490.0), ps.C["gs"], "-")]


def fig_whittaker(out):
    fig, axs = plt.subplots(1, 2, figsize=(ps.COL_2, 2.8))
    r = np.linspace(PO.a, 80, 400)
    for name, b, l, C, col, ls in STATES:
        kap, eb = PO.k(b), PO.eta(b)
        u = np.array([C * float(mp.whitw(-eb, l + 0.5, 2 * kap * x)) for x in r])
        axs[0].plot(r, u * u, color=col, ls=ls, label=name)
        # fraction of the external norm beyond r: int_r^inf u^2 / int_a^inf u^2
        rr = np.linspace(PO.a, 400, 4000)
        uu = np.array([float(mp.whitw(-eb, l + 0.5, 2 * kap * x)) ** 2 for x in rr])
        cum = np.cumsum(uu[::-1])[::-1]
        axs[1].plot(rr, cum / cum[0], color=col, ls=ls)
    axs[0].set_yscale("log")
    axs[0].set_ylim(1e-4, 1e2)
    axs[0].set_xlabel(r"$r$ (fm)")
    axs[0].set_ylabel(r"$u^2(r) = C^2 W_l^2(\kappa r)$ (fm$^{-1}$)")
    axs[1].set_xlabel(r"$r$ (fm)")
    axs[1].set_ylabel(r"fraction of $\int_a^\infty u^2\,dr$ beyond $r$")
    axs[1].set_xlim(PO.a, 80)
    axs[1].set_ylim(0, 1.02)
    axs[0].set_xlim(PO.a, 80)
    for ax, tag in zip(axs, ("(a)", "(b)")):
        ax.axvline(PO.a, color=ps.C["bkg"], ls=":", lw=0.8)
        ps.minor(ax)
        ps.label(ax, tag, loc="upper right")
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(os.path.join(out, "fig4_whittaker_tails.png"), dpi=200)
    plt.close(fig)


def fig_ec_integrand(out):
    """E1 radial integrand of external capture, p + 16O (p wave, l_i = 1) -> 17F, at E = 0.3 MeV:
    u_hs(r) r W_lf(kappa r), u_hs the hard-sphere wave F(kr) G(ka) - G(kr) F(ka), zero at r = a."""
    e, li = 0.3, 1
    k, eta = PO.k(e), PO.eta(e)
    Fa, Ga = float(mp.coulombf(li, eta, k * PO.a)), float(mp.coulombg(li, eta, k * PO.a))
    rr = np.linspace(PO.a, 250, 2500)
    uhs = np.array([float(mp.coulombf(li, eta, k * x)) * Ga - float(mp.coulombg(li, eta, k * x)) * Fa for x in rr])
    fig, axs = plt.subplots(1, 2, figsize=(ps.COL_2, 2.8))
    for name, b, l, C, col, ls in STATES:
        kap, eb = PO.k(b), PO.eta(b)
        w = np.array([float(mp.whitw(-eb, l + 0.5, 2 * kap * x)) for x in rr])
        f = uhs * rr * w
        axs[0].plot(rr, f / np.max(np.abs(f)), color=col, ls=ls, label=name)
        cum = np.cumsum(f) * (rr[1] - rr[0])
        axs[1].plot(rr, cum / cum[-1], color=col, ls=ls)
    axs[0].set_xlim(PO.a, 120)
    axs[1].set_xlim(PO.a, 120)
    axs[0].axhline(0, color="black", lw=0.5)
    axs[0].set_xlabel(r"$r$ (fm)")
    axs[0].set_ylabel(r"E1 integrand $u_{\rm hs}(r)\, r\, W_{l_f}(\kappa r)$ (scaled)")
    axs[1].set_xlabel(r"upper limit $r$ (fm)")
    axs[1].set_ylabel(r"$\int_a^r$ / $\int_a^{250\,{\rm fm}}$")
    axs[1].axhline(1, color=ps.C["bkg"], ls=":", lw=0.8)
    for ax, tag in zip(axs, ("(a)", "(b)")):
        ps.minor(ax)
        ps.label(ax, f"{tag} $E = 0.3$ MeV, $l_i = 1$", loc="lower right")
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(os.path.join(out, "fig5_ec_integrand.png"), dpi=200)
    plt.close(fig)


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(out, exist_ok=True)
    ps.apply()
    for f in (fig_coulomb, fig_penetrability, fig_shift, fig_whittaker, fig_ec_integrand):
        f(out)
        print("wrote", f.__name__, flush=True)
    # numbers quoted in the text
    for name, b, l, C, *_ in STATES:
        print(f"{name}: kappa = {PO.k(b):.4f} fm^-1 (1/kappa = {1 / PO.k(b):.1f} fm), eta_b = {PO.eta(b):.3f}")
    print(f"p+12C: a = {PC.a:.2f} fm, mu = {PC.mu:.5f} u; p+16O: a = {PO.a:.2f} fm, mu = {PO.mu:.5f} u")
    for e in (0.1, 0.5, 2.0):
        print(f"p+12C at E = {e} MeV: eta = {PC.eta(e):.3f}, ka = {PC.k(e) * PC.a:.3f}")


if __name__ == "__main__":
    main()
