#!/usr/bin/env python
"""Window non-representability certificate from a hybrid run's per-state masks.

One column per run, two rows sharing the energy axis:
  top    -- the frozen ceiling every k-point needs (k sorted by ceiling), the
            unfrozen states that sit inside the global frozen span
            ("interlopers"), and the two candidate global windows: "reach all"
            (freezes every interloper) and "stay clean" (abandons frozen
            states). The box counts what the best possible global window
            still misassigns.
  bottom -- projectability P_nk against the same energy axis, frozen against
            released states; the grey band marks conduction energies where
            both coexist, so no energy window separates them.

Accepted inputs (any mix):
  lcao2wannier  --dump-masks  <run>_masks.npz   (frozen_mask, eps_nk, P_nk)
  w2w_hybrid    <out>.masks                      (text; the wien2wannier port)

Usage:
  plot_certificate.py RUN1 [RUN2 ...] [--labels L1 L2 ...]
                             [--out certificate.pdf] [--title TEXT]
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

C_FROZEN, C_GOLD, RED = "#648FFF", "#FFB000", "#D62728"


def load(path):
    """-> frozen (nk, npool) bool, eps (nk, npool) eV rel. E_F, P (nk, npool)."""
    if path.endswith(".npz"):
        d = np.load(path)
        if "frozen_mask" in d.files:                      # lcao2wannier dump
            return d["frozen_mask"].astype(bool), d["eps_nk"], d["P_nk"]
        eig = d["eig"] - float(d["e_fermi"]) if "e_fermi" in d.files else d["eig"]
        fr, P = d["frozen"].astype(bool), d["P"]          # bridge dump
        if "pool" in d.files and eig.shape[1] != fr.shape[1]:
            eig, P = eig[:, d["pool"]], P[:, d["pool"]]   # masks are over the pool only
        return fr, eig, P
    rows = np.loadtxt(path, comments="#")
    k = rows[:, 0].astype(int)
    nk, npool = k.max(), int((k == 1).sum())
    if rows.shape[0] != nk * npool:
        raise SystemExit(f"{path}: {rows.shape[0]} rows is not nk*npool = {nk}*{npool}")
    eps = rows[:, 2].reshape(nk, npool)
    P = rows[:, 3].reshape(nk, npool)
    fr = rows[:, 4].reshape(nk, npool).astype(bool)
    return fr, eps, P


def best_window_misassigned(fr, eps):
    """min over a global ceiling E of (unfrozen with E_lo<=eps<=E) + (frozen with eps>E), E_lo the
    lowest frozen state. The count changes only at state energies, so evaluating it at E_lo - 1 eV
    and at every interloper and frozen energy gives the exact minimum; a grid of E can miss it
    (Sc-COF beta: 80 on the former 1200-point grid, 74 exactly)."""
    E_lo, E_hi = float(eps[fr].min()), float(eps[fr].max())
    inter = (~fr) & (eps >= E_lo) & (eps <= E_hi)
    e_int, e_fro = np.sort(eps[inter]), np.sort(eps[fr])
    Es = np.concatenate([[E_lo - 1.0], e_int, e_fro])
    wrong = np.searchsorted(e_int, Es, side="right")
    lost = e_fro.size - np.searchsorted(e_fro, Es, side="right")
    return int((wrong + lost).min())


def style(ax):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="major", direction="out", length=6, width=1.4)
    ax.tick_params(which="minor", direction="in", length=3, width=1.0,
                   top=True, right=True)
    for s in ax.spines.values():
        s.set_linewidth(1.4)
    ax.set_box_aspect(1.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="mask files (see above)")
    ap.add_argument("--labels", nargs="*", default=None, help="panel titles")
    ap.add_argument("--out", default="window_certificate.pdf")
    ap.add_argument("--title", default=None)
    ap.add_argument("--dpi", type=int, default=200)
    a = ap.parse_args()
    labels = a.labels or [os.path.basename(r).split("_masks")[0].split(".masks")[0]
                          for r in a.runs]
    if len(labels) != len(a.runs):
        ap.error("--labels must match the number of runs")

    plt.rcParams.update({"font.size": 12, "axes.labelsize": 13,
                         "legend.frameon": False, "pdf.fonttype": 42})
    m = len(a.runs)
    fig, axes = plt.subplots(2, m, figsize=(5.2 * m, 9.6), sharex="col",
                             constrained_layout=True, squeeze=False)
    for ax in axes.ravel():
        style(ax)

    for j, (path, label) in enumerate(zip(a.runs, labels)):
        fr, eps, P = load(path)
        nk = fr.shape[0]
        if not fr.any():
            print(f"[{label}] no frozen states; nothing to certify")
            continue
        noninterval = sum(1 for k in range(nk)
                          if fr[k].any() and np.any(np.diff(np.where(fr[k])[0]) != 1))
        if noninterval:
            print(f"[{label}] F(k) is not an energy interval at {noninterval}/{nk} k "
                  f"(a stronger obstruction than the certificate assumes)")
        E_lo, E_hi = float(eps[fr].min()), float(eps[fr].max())
        xlo, xhi = E_lo - 1.0, E_hi + 2.0
        top_e = np.array([eps[k][fr[k]].max() if fr[k].any() else np.nan
                          for k in range(nk)])
        order = np.argsort(np.nan_to_num(top_e, nan=-np.inf))
        rank = np.empty(nk, int)
        rank[order] = np.arange(nk)
        inter = (~fr) & (eps >= E_lo) & (eps <= E_hi)
        kk, nn = np.where(inter)
        u_top = float(eps[inter].min()) if kk.size else E_hi
        n_bad = best_window_misassigned(fr, eps)
        print(f"[{label}] frozen span [{E_lo:+.3f}, {E_hi:+.3f}] eV, "
              f"{int(inter.sum())} interlopers, best window misassigns {n_bad}")

        axt = axes[0, j]
        axt.fill_betweenx(np.arange(nk), xlo - 5.0, top_e[order],
                          color=C_FROZEN, alpha=0.16, lw=0)
        axt.plot(top_e[order], np.arange(nk), color=C_FROZEN, lw=2.5,
                 label="required ceiling (per $\\mathbf{k}$)")
        if kk.size:
            axt.plot(eps[kk, nn], rank[kk], "o", ms=3.5, mec="none", color=RED,
                     alpha=0.55, label="unfrozen interlopers")
        axt.axvline(E_hi, color=C_GOLD, lw=2.2, ls="--")
        axt.axvline(u_top, color=C_GOLD, lw=2.2, ls=":")
        if j == 0:
            axt.text(E_hi, 0.985, " reach all", transform=axt.get_xaxis_transform(),
                     rotation=90, va="top", ha="left", fontsize=11, color="0.2")
            axt.text(u_top, 0.985, "stay clean ", transform=axt.get_xaxis_transform(),
                     rotation=90, va="top", ha="right", fontsize=11, color="0.2")
            axt.set_ylabel("$\\mathbf{k}$-points\n(sorted by ceiling)")
            axt.legend(loc="lower left", fontsize=11)
        axt.set_ylim(0, nk - 1)
        axt.set_xlim(xlo, xhi)
        axt.set_title(label, fontsize=13)
        axt.text(0.03, 0.96, f"best window misassigns {n_bad}",
                 transform=axt.transAxes, ha="left", va="top", fontsize=11,
                 bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))

        axb = axes[1, j]
        ef, pf, eu, pu = eps[fr], P[fr], eps[~fr], P[~fr]
        su, sf = (eu >= xlo) & (eu <= xhi), (ef >= xlo) & (ef <= xhi)
        axb.scatter(eu[su], pu[su], s=16, color=RED, alpha=0.35, lw=0,
                    label="released (per-state)")
        axb.scatter(ef[sf], pf[sf], s=16, color=C_FROZEN, alpha=0.45, lw=0,
                    label="frozen (per-state)")
        mc = eps > 0.0
        efc, euc = eps[mc & fr], eps[mc & ~fr]
        if efc.size and euc.size:
            lo, hi = max(efc.min(), euc.min()), min(efc.max(), euc.max())
            if hi > lo:
                axb.axvspan(lo, hi, color="0.85", alpha=0.6, zorder=0)
                if j == 0:
                    axb.text((lo + hi) / 2, 0.03, "frozen & released\ncoexist",
                             ha="center", va="bottom", fontsize=10.5, color="0.25")
        axb.set_ylim(-0.02, 1.03)
        axb.set_xlim(xlo, xhi)
        axb.set_xlabel("$E - E_F$ (eV)")
        if j == 0:
            axb.set_ylabel("projectability  $P_{n\\mathbf{k}}$")
            axb.legend(loc="lower left", fontsize=11, markerscale=2.0)

    if a.title:
        fig.suptitle(a.title, fontsize=14)
    fig.savefig(a.out, dpi=a.dpi)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
