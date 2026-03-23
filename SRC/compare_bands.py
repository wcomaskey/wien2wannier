#!/usr/bin/env python3
"""
wien2wannier/SRC/compare_bands.py

General-purpose band structure comparison: DFT (Wien2k spaghetti) vs
Wannier-interpolated bands (Wannier90).

Uses INDEX-BASED segment mapping — no BZ convention assumptions needed.
Both DFT and W90 paths share the same high-symmetry labels; within each
segment, DFT k-points are linearly mapped onto the W90 k-axis by their
fractional position within the segment. This is convention-independent
and works for all Bravais lattice types.

Usage:
    compare_bands.py CASE [options]

    -dir DIR     working directory (default: .)
    -out FILE    output PNG file (default: CASE_band_comparison.png)
    -emin E      energy min relative to E_F (default: -15 eV)
    -emax E      energy max relative to E_F (default: +10 eV)
    -ezoom E     zoom window half-width around E_F (default: 4 eV)
"""

import numpy as np
import sys, os, re, argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_w90_bands(bandfile):
    """Read Wannier90 _band.dat file → (k_axis, bands[nk, nb])"""
    data = np.loadtxt(bandfile)
    k, e = data[:, 0], data[:, 1]
    nk = np.argmax(np.diff(k) < -1e-6) + 1
    if nk <= 1:
        nk = len(k)
    nb = len(k) // nk
    return k[:nk], e.reshape(nb, nk).T


def read_w90_labels(labelfile, kaxis):
    """Read _band.labelinfo.dat → (labels[], k_positions[])"""
    labels, kpos = [], []
    with open(labelfile) as f:
        for line in f:
            parts = line.split()
            labels.append(parts[0])
            idx = int(parts[1]) - 1  # 1-indexed → 0-indexed
            kpos.append(kaxis[min(idx, len(kaxis) - 1)])
    return labels, np.array(kpos)


def read_spaghetti(spagfile):
    """Read Wien2k spaghetti_ene → list of arrays, each (nk, 2) = [k, E]"""
    bands = []
    current = []
    with open(spagfile) as f:
        for line in f:
            if 'bandindex' in line:
                if current:
                    bands.append(np.array(current))
                current = []
                continue
            parts = line.split()
            if len(parts) >= 5:
                current.append([float(parts[3]), float(parts[4])])
    if current:
        bands.append(np.array(current))
    return bands


def read_klist_band_labels(klistfile):
    """Read klist_band, return (labels[], indices[]) of labeled k-points."""
    labels, indices = [], []
    idx = 0
    with open(klistfile) as f:
        for line in f:
            s = line.strip()
            if s == 'END' or not s:
                continue
            label = line[0:10].strip()
            if label and not label.replace('-', '').replace('.', '').isdigit():
                labels.append(label)
                indices.append(idx)
            idx += 1
    return labels, indices


def map_dft_to_w90_kaxis(dft_bands, dft_hsym_idx, w90_kpos, nseg):
    """
    Map DFT k-points onto W90 k-axis by index within each segment.
    Convention-independent: uses fractional position within segment.
    """
    nk_dft = len(dft_bands[0])
    k_mapped = np.full(nk_dft, np.nan)

    for iseg in range(nseg):
        i_lo = dft_hsym_idx[iseg]
        i_hi = dft_hsym_idx[iseg + 1]
        w_lo = w90_kpos[iseg]
        w_hi = w90_kpos[iseg + 1]
        n_seg = i_hi - i_lo

        for j in range(i_lo, i_hi + (1 if iseg == nseg - 1 else 0)):
            frac = (j - i_lo) / n_seg if n_seg > 0 else 0
            k_mapped[j] = w_lo + frac * (w_hi - w_lo)

    return k_mapped


def find_band_offset(dft_bands, w90_bands):
    """
    Find which DFT band index corresponds to W90 band 1.
    Matches by eigenvalue at Gamma (k-index 0).
    Returns offset: W90 band i = DFT band (i + offset).
    """
    w90_e0 = w90_bands[0, 0]
    dft_e0 = [b[0, 1] for b in dft_bands]
    return min(range(len(dft_e0)), key=lambda i: abs(dft_e0[i] - w90_e0))


def compute_energy_shift(dft_bands, w90_bands, offset, n_match=5):
    """
    Compute constant energy shift between DFT and W90 at Gamma.
    Uses median over first n_match bands for robustness.
    """
    nb_w90 = w90_bands.shape[1]
    n = min(n_match, nb_w90, len(dft_bands) - offset)
    shifts = [dft_bands[offset + i][0, 1] - w90_bands[0, i] for i in range(n)]
    return np.median(shifts)


def main():
    parser = argparse.ArgumentParser(description='Compare DFT vs Wannier bands')
    parser.add_argument('case', help='Case name (e.g. WANN, GaAs-WANN)')
    parser.add_argument('-dir', default='.', help='Working directory')
    parser.add_argument('-out', default=None, help='Output PNG file')
    parser.add_argument('-emin', type=float, default=-15, help='Energy min (eV rel to E_F)')
    parser.add_argument('-emax', type=float, default=+10, help='Energy max (eV rel to E_F)')
    parser.add_argument('-ezoom', type=float, default=4, help='Zoom half-width (eV)')
    args = parser.parse_args()

    d = args.dir
    case = args.case
    outfile = args.out or os.path.join(d, f'{case}_band_comparison.png')

    # --- Read data ---
    with open(os.path.join(d, f'{case}.fermi')) as f:
        ef_ev = float(f.read().strip()) * 13.605693

    kaxis_w90, bands_w90 = read_w90_bands(os.path.join(d, f'{case}_band.dat'))
    w90_labels, w90_kpos = read_w90_labels(
        os.path.join(d, f'{case}_band.labelinfo.dat'), kaxis_w90)

    dft_bands = read_spaghetti(os.path.join(d, f'{case}.spaghetti_ene'))
    dft_labels, dft_hsym_idx = read_klist_band_labels(
        os.path.join(d, f'{case}.klist_band'))

    # --- Validate path consistency ---
    if len(dft_labels) != len(w90_labels):
        print(f"WARNING: path mismatch: DFT has {len(dft_labels)} labels, "
              f"W90 has {len(w90_labels)}")
        print(f"  DFT: {' -> '.join(dft_labels)}")
        print(f"  W90: {' -> '.join(w90_labels)}")

    nseg = min(len(dft_labels), len(w90_labels)) - 1

    # --- Map DFT k-axis to W90 k-axis (index-based, convention-free) ---
    k_mapped = map_dft_to_w90_kaxis(dft_bands, dft_hsym_idx, w90_kpos, nseg)
    valid = ~np.isnan(k_mapped)

    # --- Align energy reference ---
    offset = find_band_offset(dft_bands, bands_w90)
    eshift = compute_energy_shift(dft_bands, bands_w90, offset)

    # --- Report ---
    nb_w90 = bands_w90.shape[1]
    nk_w90 = bands_w90.shape[0]
    print(f"Fermi energy: {ef_ev:.4f} eV")
    print(f"Wannier: {nb_w90} bands, {nk_w90} k-pts")
    print(f"DFT: {len(dft_bands)} bands, {np.sum(valid)}/{len(dft_bands[0])} mapped k-pts")
    print(f"Path: {' -> '.join(dft_labels[:nseg+1])}")
    print(f"DFT band offset: {offset} (W90 band 1 = DFT band {offset+1})")
    print(f"Energy shift: {eshift:.4f} eV ({eshift*1000:.1f} meV)")

    # --- Plot ---
    pretty = lambda s: s.replace('GAMMA', r'$\Gamma$')
    fig, axes = plt.subplots(1, 2, figsize=(14, 7),
                             gridspec_kw={'width_ratios': [2, 1]})

    for ax, (ylo, yhi), title in [
        (axes[0], (ef_ev + args.emin, ef_ev + args.emax),
         'DFT vs PDWF Wannier Bands'),
        (axes[1], (ef_ev - args.ezoom, ef_ev + args.ezoom),
         'Zoom: Near Fermi Level'),
    ]:
        # DFT bands
        for ib, band in enumerate(dft_bands):
            lbl = 'DFT (Wien2k)' if ib == 0 else None
            ax.plot(k_mapped[valid], band[valid, 1] - eshift,
                    'b-', lw=1.2, alpha=0.7, label=lbl)

        # Wannier bands
        for ib in range(nb_w90):
            lbl = 'PDWF Wannier' if ib == 0 else None
            ax.plot(kaxis_w90, bands_w90[:, ib],
                    'r-.', lw=1.5, alpha=0.8, label=lbl)

        # Fermi level and grid
        ax.axhline(ef_ev, color='green', lw=0.8, ls=':', alpha=0.5)
        for pos in w90_kpos:
            ax.axvline(pos, color='gray', lw=0.5, ls='--')
        ax.set_xticks(w90_kpos)
        ax.set_xticklabels([pretty(l) for l in w90_labels], fontsize=11)
        ax.set_xlim(kaxis_w90[0], kaxis_w90[-1])
        ax.set_ylim(ylo, yhi)
        ax.set_ylabel('Energy (eV)', fontsize=12)
        ax.set_title(title, fontsize=14)

    axes[0].legend(fontsize=11, loc='lower right')
    plt.tight_layout()
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"Plot saved to {outfile}")


if __name__ == '__main__':
    main()
