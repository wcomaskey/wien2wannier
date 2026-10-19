#!/usr/bin/env python
"""Auto-trust anatomy: the derived thresholds drawn over the sorted band-averaged
projectabilities of the disentanglement pool.

Each panel: bands sorted by Pbar_n (step curve), coloured by the class the rule
assigned (trusted / gray / junk; open markers = bands whose bottom lies above
the fidelity ceiling, so they never enter the frozen set), with the derived
band gate p_gate (dashed) and per-state threshold p_froz (dotted) in gold and
the gray plateau shaded between the plateau top and the gate. Use it to see
why a run froze what it froze before touching --hybrid-p-froz / --froz-emax.

Accepted inputs (any mix):
  lcao2wannier  --hybrid-report-json <run>_report.json
  w2w_hybrid    <out>.anatomy

Usage:
  plot_anatomy.py RUN1 [RUN2 ...] [--labels L1 L2 ...]
                             [--out anatomy.pdf] [--cols N]
"""
import argparse
import json
import math
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

CLS_COLOR = {"trusted": "#648FFF", "gray": "#FE6100", "junk": "#B0B0B0"}
C_GOLD, C_GREY = "#FFB000", "#7F7F7F"


def load(path):
    """-> (Pbar, cls, relevant) arrays over the pool bands, thresholds dict."""
    if path.endswith(".json"):
        with open(path) as fh:
            d = json.load(fh)
        rows = d["bands"]
        thr = dict(d["thresholds"])
        return (np.array([r["Pbar"] for r in rows]), np.array([r["cls"] for r in rows]),
                np.array([int(r["relevant"]) for r in rows], bool), thr)
    thr = {}
    P, cls, rel = [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                # header tokens look like "p_gate=  0.6078" or "plateau=none (...)"
                for key, val in re.findall(r"(\w+)=\s*(\S+)", line[1:]):
                    try:
                        thr[key] = float(val)
                    except ValueError:
                        thr[key] = None if val == "none" else val
                continue
            f = line.split()
            if len(f) < 6:
                continue
            P.append(float(f[1])); rel.append(f[4] == "1"); cls.append(f[5])
    return np.array(P), np.array(cls), np.array(rel, bool), thr


def panel(ax, P, cls, rel, thr, label):
    order = np.argsort(P)[::-1]
    P, cls, rel = P[order], cls[order], rel[order]
    x = np.arange(len(P))
    ax.step(x, P, where="mid", color=C_GREY, lw=2.0, zorder=1)
    ms = 8 if len(P) <= 60 else (5 if len(P) <= 150 else 3.5)
    mew = 1.0 if len(P) <= 150 else 0.5
    for c in ("trusted", "gray", "junk"):
        m = cls == c
        ax.plot(x[m & rel], P[m & rel], "o", ms=ms, mec="black", mew=mew,
                color=CLS_COLOR[c], ls="none", zorder=3, label=c)
        ax.plot(x[m & ~rel], P[m & ~rel], "o", ms=ms, mec=CLS_COLOR[c],
                mew=max(mew, 0.8), mfc="white", ls="none", zorder=3)
    p_gate, p_froz = thr.get("p_gate"), thr.get("p_froz")
    if p_gate is not None:
        ax.axhline(p_gate, color=C_GOLD, lw=2.2, ls="--", zorder=2)
        ax.text(0.98, p_gate, f"  $p_{{gate}}={p_gate:.3f}$",
                transform=ax.get_yaxis_transform(), ha="right", va="bottom", fontsize=11)
    if p_froz is not None:
        ax.axhline(p_froz, color=C_GOLD, lw=2.2, ls=":", zorder=2)
        ax.text(0.98, p_froz, f"  $p_{{froz}}={p_froz:.3f}$",
                transform=ax.get_yaxis_transform(), ha="right", va="top", fontsize=11)
    plateau = thr.get("plateau")
    if plateau is not None and p_gate is not None and plateau > -1e6:
        ax.axhspan(plateau, p_gate, color=C_GOLD, alpha=0.25, zorder=0)
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlabel("bands, sorted by $\\bar{P}_n$")
    ax.set_ylabel("band-averaged projectability $\\bar{P}_n$")
    ax.set_title(label, fontsize=13)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="major", direction="out", length=6, width=1.4)
    ax.tick_params(which="minor", direction="in", length=3, width=1.0, top=True, right=True)
    for s in ax.spines.values():
        s.set_linewidth(1.4)
    ax.set_box_aspect(1.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="report .json or .anatomy files")
    ap.add_argument("--labels", nargs="*", default=None)
    ap.add_argument("--out", default="auto_trust_anatomy.pdf")
    ap.add_argument("--cols", type=int, default=3, help="panels per row")
    ap.add_argument("--dpi", type=int, default=200)
    a = ap.parse_args()
    labels = a.labels or [os.path.basename(r).split("_report")[0].split(".anatomy")[0]
                          for r in a.runs]
    if len(labels) != len(a.runs):
        ap.error("--labels must match the number of runs")
    plt.rcParams.update({"font.size": 12, "axes.labelsize": 13,
                         "legend.frameon": False, "pdf.fonttype": 42})
    n = len(a.runs)
    ncols = min(a.cols, n)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 5.2 * nrows),
                             constrained_layout=True, squeeze=False)
    flat = axes.ravel()
    for ax, path, label in zip(flat, a.runs, labels):
        P, cls, rel, thr = load(path)
        panel(ax, P, cls, rel, thr, label)
        print(f"[{label}] {len(P)} pool bands: "
              + ", ".join(f"{c} {int((cls == c).sum())}" for c in ("trusted", "gray", "junk"))
              + f"; p_gate={thr.get('p_gate')}, p_froz={thr.get('p_froz')}, "
                f"plateau={thr.get('plateau')}")
    for ax in flat[n:]:
        ax.set_axis_off()
    flat[0].legend(loc="lower left", fontsize=11)
    fig.savefig(a.out, dpi=a.dpi)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
