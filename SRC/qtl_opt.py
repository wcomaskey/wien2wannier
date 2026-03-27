#!/usr/bin/env python3
"""
wien2wannier/SRC/qtl_opt.py

QTL-guided optimization of Wannier90 frozen/disentangle windows.

Loss function balances projectability, window coverage, and constraints:
  L = -(w_p * P + w_c * C + w_d * D + w_f * F) + penalties

where:
  P = sum(target_char[frozen]) / num_wann     (projectability)
  C = (fmax - fmin) / ref_width               (coverage)
  D = min(1, n_disentangle / num_wann)         (gauge freedom)
  F = n_frozen / num_wann                      (frozen count)

Constraints:
  - Frozen window inside outer window
  - At every k: bands_in_frozen <= num_wann
  - At every k: bands_in_outer > num_wann (disentanglement required)
  - Frozen window >= minimum floor around E_F
  - Frozen window wider below E_F than above (valence priority)

Usage:
    qtl_opt.py CASE -dir DIR [options]
"""

import numpy as np
import argparse, os, sys, re, json


def parse_qtl(qtlfile):
    """Parse Wien2k .qtl → (qtl_char[nb,nat,4], qtl_int[nb], meta)."""
    with open(qtlfile) as f:
        lines = f.readlines()

    ef_ry = float(re.search(r'FERMI ENERGY=\s*([\d.+-]+)', lines[2]).group(1))
    nat = int(re.search(r'NAT=\s*(\d+)', lines[3]).group(1))
    spin = int(re.search(r'SPIN=\s*(\d+)', lines[3]).group(1))
    isplit = [int(re.search(r'ISPLIT=\s*(\d+)', lines[4+ja]).group(1))
              for ja in range(nat)]

    band_starts = [i for i, l in enumerate(lines) if l.strip().startswith('BAND')]
    nbands = len(band_starts)
    lines_per_k = nat + 1
    nk = (band_starts[1] - band_starts[0] - 1) // lines_per_k if nbands >= 2 else 1

    qtl_char = np.zeros((nbands, nat, 4))
    qtl_int = np.zeros(nbands)
    for ib, bs in enumerate(band_starts):
        for ik in range(nk):
            off = bs + 1 + ik * lines_per_k
            for ja in range(nat):
                p = lines[off + ja].split()
                qtl_char[ib, ja, 0] += float(p[3])
                qtl_char[ib, ja, 1] += float(p[4])
                qtl_char[ib, ja, 2] += float(p[5])
                if len(p) > 8:
                    qtl_char[ib, ja, 3] += float(p[8])
            ip = lines[off + nat].split()
            qtl_int[ib] += float(ip[2]) if len(ip) > 2 else 0
        qtl_char[ib] /= nk
        qtl_int[ib] /= nk

    return qtl_char, qtl_int, {
        'ef_ry': ef_ry, 'nat': nat, 'spin': spin,
        'isplit': isplit, 'nbands': nbands, 'nk_qtl': nk,
    }


def parse_eig(eigfile):
    """Parse .eig → eigvals[nb, nk]."""
    eig = {}
    with open(eigfile) as f:
        for line in f:
            p = line.split()
            eig[(int(p[0]), int(p[1]))] = float(p[2])
    nb = max(ib for ib, ik in eig)
    nk = max(ik for ib, ik in eig)
    ev = np.zeros((nb, nk))
    for (ib, ik), e in eig.items():
        ev[ib-1, ik-1] = e
    return ev


def valence_config(Z):
    """Target l-channels for element Z."""
    if Z <= 2: return [0]
    if Z <= 10: return [0, 1]
    if Z <= 20: return [0, 1]
    if Z <= 30: return [0, 2]
    if Z <= 36: return [0, 1]
    if Z <= 38: return [0, 1]
    if Z <= 48: return [0, 2]
    if Z <= 54: return [0, 1]
    if Z <= 56: return [0, 1]
    if Z == 57: return [0, 1, 2]
    if Z <= 71: return [0, 2, 3]
    if Z <= 80: return [0, 2]
    return [0, 1]


def compute_targets(qtl_char, struct_file):
    """Compute target character and num_wann from struct + QTL."""
    with open(struct_file) as f:
        slines = f.readlines()
    nneq = int(re.search(r'ATOMS[:\s]+(\d+)', slines[1]).group(1))

    Z_list, mult_list = [], []
    i = 4
    for iat in range(nneq):
        i += 1  # position line
        mult = int(slines[i].split('=')[1].split()[0])
        mult_list.append(mult)
        i += 1  # mult line
        for _ in range(mult - 1):
            i += 1  # equivalent atom position lines
        Z_list.append(int(float(re.search(r'Z:\s*([\d.]+)', slines[i]).group(1))))
        i += 1  # element/Z line
        i += 3  # rot matrix (3 lines)

    targets = []
    num_wann = 0
    for ja, Z in enumerate(Z_list):
        for l in valence_config(Z):
            targets.append((ja, l))
            num_wann += (2 * l + 1) * mult_list[ja]

    nb = qtl_char.shape[0]
    tc = np.zeros(nb)
    for ib in range(nb):
        for ja, l in targets:
            tc[ib] += qtl_char[ib, ja, l]

    return tc, targets, num_wann


def loss_function(fmin, fmax, dmin, dmax, eigvals, target_char,
                  num_wann, exclude, min_fwin):
    """
    Evaluate window config. Returns (loss, details).
    Lower loss = better. Returns 1e6 for invalid configs.
    """
    nb, nk = eigvals.shape
    active = [ib for ib in range(nb) if ib not in exclude]

    # Hard constraints
    if fmin >= fmax or dmin >= dmax:
        return 1e6, {}
    if fmin < dmin or fmax > dmax:
        return 1e6, {}
    # Minimum frozen window: HARD constraint.
    # The floor is computed adaptively from the band structure to ensure
    # Fermi-crossing bands and gap-edge bands are fully inside.
    floor_penalty = 0
    if fmin > -min_fwin[0] or fmax < min_fwin[1]:
        return 1e6, {}

    # Band counts at each k
    min_outer = nb
    max_frozen = 0
    for ik in range(nk):
        n_out = sum(1 for ib in active
                    if dmin <= eigvals[ib, ik] <= dmax)
        n_frz = sum(1 for ib in active
                    if fmin <= eigvals[ib, ik] <= fmax)
        min_outer = min(min_outer, n_out)
        max_frozen = max(max_frozen, n_frz)

    if min_outer <= num_wann:  # need bands_in_outer > num_wann
        return 1e6, {}
    if max_frozen > num_wann:
        return 1e6, {}

    # Frozen bands (fully inside frozen window at all k)
    frozen = [ib for ib in active
              if eigvals[ib].min() >= fmin and eigvals[ib].max() <= fmax]

    # ---- Scoring components ----

    nf = len(frozen)
    n_dis = min_outer - nf
    e_avg_band = eigvals.mean(axis=1)

    # 1. Projectability: Fermi-proximity-weighted average.
    #    Split Gaussian: sigma_below = |E_deepest_active|,
    #    sigma_above = 0.5 * sigma_below (valence > conduction).
    e_deepest = min(e_avg_band[ib] for ib in active if e_avg_band[ib] < 0)
    sigma_below = max(3.0, abs(e_deepest))
    sigma_above = max(3.0, 0.5 * sigma_below)
    if nf > 0:
        wsum, wpsum = 0.0, 0.0
        for ib in frozen:
            sig = sigma_below if e_avg_band[ib] < 0 else sigma_above
            w = np.exp(-0.5 * (e_avg_band[ib] / sig) ** 2)
            wsum += w
            wpsum += target_char[ib] * w
        avg_proj = wpsum / wsum
    else:
        avg_proj = 0

    # 2. Conduction restraint: penalize fmax extending above E_F.
    #    Empirically, tighter frozen windows above E_F give better
    #    band structures because conduction bands benefit from gauge
    #    freedom more than from being frozen. Use a steep penalty
    #    that increases linearly with fmax above a material-dependent
    #    threshold (half the valence bandwidth).
    conduction_width = max(0, fmax)
    valence_width = max(0, -fmin)
    cond_threshold = 0.5 * valence_width  # half the valence side
    overextend = max(0, conduction_width - cond_threshold) * 0.05

    # 3. Coverage: only reward valence-side coverage (below E_F).
    #    Conduction coverage is NOT rewarded — it's handled by the
    #    outer window (disentanglement), not the frozen window.
    coverage = min(1.0, valence_width / 15.0)

    # 4. Gauge freedom: strongly reward having disentangle bands.
    gauge = np.log1p(n_dis) / np.log1p(num_wann)

    # 5. Frozen cap: freeze at most ~65% of num_wann.
    max_frozen_frac = 0.65
    if nf > max_frozen_frac * num_wann:
        return 1e6, {}

    # Minimum frozen: must freeze at least 2 bands (otherwise
    # the frozen window is meaningless)
    if nf < 2:
        return 1e6, {}

    # 6. Parsimony: mild preference for fewer frozen bands.
    #    Uses sqrt to compress the range — going from 6→5 frozen
    #    matters less than going from 2→1.
    parsimony = np.sqrt(1.0 - nf / (max_frozen_frac * num_wann))

    # 7. Weak band penalty
    active_median = np.median([target_char[ib] for ib in active])
    weak_frozen = sum(1 for ib in frozen if target_char[ib] < active_median)
    weak_penalty = 0.1 * weak_frozen / max(1, nf)

    # 8. Asymmetry: valence side should be at least 2x the conduction side
    asym = 0
    if conduction_width > 0 and valence_width > 0:
        ratio = valence_width / conduction_width
        if ratio < 2.0:
            asym = 0.1 * (2.0 - ratio)

    # 9. Target coverage: fraction of the num_wann highest-character
    #    bands that are frozen. These are the bands most likely to
    #    form the Wannier manifold. For semiconductors, these are the
    #    valence+gap-edge bands. For metals, these are the bands near E_F.
    #    We pick the top num_wann bands by target character — if all of
    #    them are frozen, fermi_cov = 1.0.
    tc_ranked = sorted(active, key=lambda ib: target_char[ib], reverse=True)
    target_bands = tc_ranked[:num_wann]  # top num_wann by character
    frozen_targets = [ib for ib in frozen if ib in target_bands]
    fermi_cov = len(frozen_targets) / max(1, len(target_bands))

    # Combine — Fermi coverage and projectability are dominant
    score = (1.0 * avg_proj      # quality of frozen bands
             + 0.3 * coverage    # valence coverage
             + 0.4 * gauge       # gauge freedom
             + 0.6 * fermi_cov   # Fermi-level band coverage (key!)
             + 0.1 * parsimony   # very mild fewer-frozen preference
             - weak_penalty
             - asym
             - floor_penalty
             - overextend)

    return -score, {
        'frozen': [ib + 1 for ib in frozen],
        'nf': nf, 'nd': n_dis,
        'min_outer': min_outer, 'max_frozen': max_frozen,
        'avg_proj': avg_proj, 'coverage': coverage,
        'gauge': gauge, 'parsimony': parsimony, 'overextend': overextend,
        'fermi_cov': fermi_cov,
        'weak_penalty': weak_penalty, 'asym': asym,
        'score': score,
    }


def find_gap_edges(eigvals, active):
    """Find energy gap midpoints between consecutive active bands.
    Returns sorted list of gap midpoint energies suitable for window edges."""
    nb = eigvals.shape[0]
    # Sort active bands by k-averaged energy
    e_avg = eigvals.mean(axis=1)
    sorted_active = sorted(active, key=lambda ib: e_avg[ib])

    edges = []
    for i in range(len(sorted_active) - 1):
        ib_lo = sorted_active[i]
        ib_hi = sorted_active[i + 1]
        gap_min = eigvals[ib_hi].min() - eigvals[ib_lo].max()
        if gap_min > 0.5:  # only use real gaps (> 0.5 eV)
            midpoint = 0.5 * (eigvals[ib_lo].max() + eigvals[ib_hi].min())
            edges.append(midpoint)

    # Also add edges just outside the active range
    e_lo = min(eigvals[ib].min() for ib in active)
    e_hi = max(eigvals[ib].max() for ib in active)
    edges.append(e_lo - 0.5)
    edges.append(e_hi + 0.5)

    return sorted(set(edges))


def grid_sweep(eigvals, target_char, num_wann, exclude,
               min_fwin=(4.0, 2.0), step=0.5):
    """Sweep frozen window using gap-snapped edges + regular grid."""
    active = [ib for ib in range(eigvals.shape[0]) if ib not in exclude]
    e_lo = min(eigvals[ib].min() for ib in active)
    e_hi = max(eigvals[ib].max() for ib in active)

    dmin = e_lo - 2.0
    dmax = e_hi + 2.0

    # Build candidate fmin/fmax values:
    # Regular grid + gap-snapped edges (ensures we test mid-gap positions)
    gap_edges = find_gap_edges(eigvals, active)
    fmin_candidates = sorted(set(
        list(np.arange(e_lo - 1, 1, step)) + [e for e in gap_edges if e < 1]
    ))
    fmax_candidates = sorted(set(
        list(np.arange(0, e_hi + 1, step)) + [e for e in gap_edges if e > 0]
    ))

    results = []
    for fmin in fmin_candidates:
        for fmax in fmax_candidates:
            if fmax <= fmin:
                continue
            loss, det = loss_function(
                fmin, fmax, dmin, dmax, eigvals, target_char,
                num_wann, exclude, min_fwin)
            if loss < 1e5:
                results.append((loss, fmin, fmax, dmin, dmax, det))

    results.sort()
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('case')
    ap.add_argument('-dir', default='.')
    ap.add_argument('-qtl', default=None)
    ap.add_argument('-nwann', type=int, default=None)
    ap.add_argument('-efmin', type=float, default=4.0)
    ap.add_argument('-efmax', type=float, default=2.0)
    ap.add_argument('-step', type=float, default=0.5)
    ap.add_argument('-top', type=int, default=15)
    ap.add_argument('-bmin', type=int, default=None,
                    help='bmin used in w2w run (for QTL-EIG alignment)')
    ap.add_argument('-write', action='store_true')
    args = ap.parse_args()

    d = args.dir
    case = args.case

    qtl_char, qtl_int, meta = parse_qtl(args.qtl or f'{d}/{case}.qtl')
    eigvals = parse_eig(f'{d}/{case}.eig')
    target_char, targets, nw_auto = compute_targets(qtl_char, f'{d}/{case}.struct')
    num_wann = args.nwann or nw_auto
    nb_eig, nk = eigvals.shape
    nb_qtl = qtl_char.shape[0]

    # Determine bmin: which Wien2k band does EIG band 1 correspond to?
    # QTL always starts from Wien2k band 1. EIG starts from Wien2k band bmin.
    # Auto-detect by matching k-averaged eigenvalues between QTL and EIG.
    # QTL energies are in Ry (absolute), EIG are in eV (Fermi-shifted by w2w).
    # We Fermi-shift the QTL energies and convert to eV, then match.

    # Determine bmin by matching band energy PATTERNS between QTL and EIG.
    # Both have k-averaged energies but may use different Fermi references.
    # Strategy: compute energy GAPS between consecutive bands (reference-free),
    # then slide the EIG pattern along the QTL pattern to find best alignment.

    # QTL k-averaged energies (in Ry, then convert to eV — no Fermi shift needed
    # since we only use gaps)
    qtl_lines = open(args.qtl or f'{d}/{case}.qtl').readlines()
    qtl_band_starts = [i for i, l in enumerate(qtl_lines) if l.strip().startswith('BAND')]
    qtl_kavg = []
    for bs in qtl_band_starts:
        lpk = meta['nat'] + 1
        nk_q = meta['nk_qtl']
        esum = sum(float(qtl_lines[bs + 1 + ik * lpk].split()[0])
                   for ik in range(nk_q))
        qtl_kavg.append(esum / nk_q * 13.605693)  # Ry -> eV, absolute

    # EIG k-averaged energies (in eV, Fermi-shifted)
    eig_kavg = eigvals.mean(axis=1).tolist()

    # Compute gap patterns (differences between consecutive bands)
    qtl_gaps = [qtl_kavg[i+1] - qtl_kavg[i] for i in range(len(qtl_kavg)-1)]
    eig_gaps = [eig_kavg[i+1] - eig_kavg[i] for i in range(len(eig_kavg)-1)]

    # Slide EIG gap pattern along QTL gap pattern, find best alignment
    best_score = 1e10
    bmin = 1
    n_match = min(len(eig_gaps), 5)  # match first 5 gaps
    for offset in range(len(qtl_gaps) - n_match + 1):
        score = sum(abs(qtl_gaps[offset + i] - eig_gaps[i]) for i in range(n_match))
        if score < best_score:
            best_score = score
            bmin = offset + 1  # 1-indexed Wien2k band

    # Verify with absolute energy matching (using Fermi-shifted QTL)
    ef_ry = meta['ef_ry']
    qtl_kavg_shifted = [(e - ef_ry * 13.605693) for e in qtl_kavg]
    abs_diff = abs(qtl_kavg_shifted[bmin - 1] - eig_kavg[0])
    if abs_diff > 1.0:
        # Gap pattern might be ambiguous, try fermi file
        fermi_file = f'{d}/{case}.fermi'
        if os.path.exists(fermi_file):
            ef_w2w = float(open(fermi_file).read().strip()) * 13.605693
            qtl_kavg_w2w = [(e - ef_w2w) for e in qtl_kavg]
            diffs = [abs(qtl_kavg_w2w[q] - eig_kavg[0]) for q in range(len(qtl_kavg_w2w))]
            bmin_alt = np.argmin(diffs) + 1
            if diffs[bmin_alt - 1] < abs_diff:
                bmin = bmin_alt

    print(f"System: EIG={nb_eig} bands, QTL={nb_qtl} bands, {nk} k-pts")
    print(f"bmin={bmin} (EIG band 1 = Wien2k band {bmin})")
    print(f"num_wann={num_wann}, Targets: {targets}")
    print(f"Min frozen: [-{args.efmin}, +{args.efmax}] eV\n")

    # Align QTL to EIG: only use QTL bands that are in the EIG window
    # QTL band q (0-indexed) -> EIG band (q - bmin + 1) (0-indexed)
    # We need target_char and exclude indexed by EIG band number
    nb = nb_eig
    tc_aligned = np.zeros(nb)
    qtl_aligned = np.zeros((nb, meta['nat'], 4))
    for eig_ib in range(nb):
        qtl_ib = eig_ib + bmin - 1  # 0-indexed QTL band
        if 0 <= qtl_ib < nb_qtl:
            tc_aligned[eig_ib] = target_char[qtl_ib]
            qtl_aligned[eig_ib] = qtl_char[qtl_ib]
    target_char = tc_aligned
    qtl_char_aligned = qtl_aligned

    # Exclude core/semicore via energy gaps + orbital character
    e_avg = eigvals.mean(axis=1)
    exclude = set()

    # 1. Energy gap detection: exclude ALL bands below the largest
    #    gap > 5 eV below E_F. These are core/semicore states
    #    regardless of their orbital character (e.g. Mg 2p has
    #    high sp character but is semicore).
    largest_gap, largest_gap_idx = 0, -1
    for ib in range(nb - 1):
        gap = e_avg[ib + 1] - e_avg[ib]
        if e_avg[ib] < 0 and gap > 5.0 and gap > largest_gap:
            largest_gap = gap
            largest_gap_idx = ib
    if largest_gap_idx >= 0:
        for jb in range(largest_gap_idx + 1):
            exclude.add(jb)

    # 2. Pure d/f semicore (> 90% single-orbital character)
    for ib in range(nb):
        for ja in range(meta['nat']):
            if qtl_char_aligned[ib, ja, 2] > 0.90 or qtl_char_aligned[ib, ja, 3] > 0.90:
                exclude.add(ib)

    # 3. Bands with negligible target character (< 1% of max target)
    #    These are free-electron-like or belong to non-target orbitals.
    max_tc = max(target_char) if len(target_char) > 0 else 1.0
    for ib in range(nb):
        if ib not in exclude and target_char[ib] < 0.01 * max_tc:
            exclude.add(ib)

    # Align qtl_int too
    qtl_int_aligned = np.zeros(nb)
    for eig_ib in range(nb):
        qtl_ib = eig_ib + bmin - 1
        if 0 <= qtl_ib < nb_qtl:
            qtl_int_aligned[eig_ib] = qtl_int[qtl_ib]

    print(f"{'Band':>4s} {'W2k':>4s} {'E_avg':>7s} {'target':>7s} {'int':>5s} {'Status':>8s}")
    print('-' * 42)
    for ib in range(nb):
        st = 'EXCL' if ib in exclude else ''
        w2k = ib + bmin
        print(f"{ib+1:4d} {w2k:4d} {e_avg[ib]:7.1f} {target_char[ib]:7.4f} "
              f"{qtl_int_aligned[ib]:5.3f} {st:>8s}")

    print(f"\nExcluded: {sorted(ib+1 for ib in exclude)}")

    # Compute adaptive minimum frozen window from band structure.
    # fmin_floor: must cover all active bands below E_F that have
    #   significant target character (> 10% of max).
    # fmax_floor: must cover all Fermi-crossing bands (metals) or
    #   the first conduction band above the gap (semiconductors).
    active_list = [ib for ib in range(nb) if ib not in exclude]
    tc_thresh = 0.1 * max(target_char[ib] for ib in active_list)

    # fmin: deepest target band below E_F
    target_below = [ib for ib in active_list
                    if e_avg[ib] < 0 and target_char[ib] > tc_thresh]
    if target_below:
        fmin_floor = abs(min(eigvals[ib].min() for ib in target_below)) + 0.5
    else:
        fmin_floor = args.efmin

    # fmax: detect Fermi-crossing bands or first conduction band
    crossing = [ib for ib in active_list
                if eigvals[ib].min() < 0 and eigvals[ib].max() > 0]
    if crossing:
        # Metal: fmax must cover the top of all Fermi-crossing bands
        fmax_floor = max(eigvals[ib].max() for ib in crossing) + 0.5
    else:
        # Semiconductor: fmax must cover first conduction band
        cond_bands = [ib for ib in active_list if e_avg[ib] > 0]
        if cond_bands:
            first_cond = min(cond_bands, key=lambda ib: e_avg[ib])
            fmax_floor = eigvals[first_cond].max() + 0.5
        else:
            fmax_floor = args.efmax

    # Use the larger of adaptive and user-specified floors
    fmin_floor = max(fmin_floor, args.efmin)
    fmax_floor = max(fmax_floor, args.efmax)

    print(f"Adaptive frozen floor: [-{fmin_floor:.1f}, +{fmax_floor:.1f}] eV")
    if crossing:
        print(f"  (metal: {len(crossing)} Fermi-crossing bands, "
              f"max extent = +{max(eigvals[ib].max() for ib in crossing):.1f} eV)")
    else:
        print(f"  (semiconductor: first conduction band to +{fmax_floor-0.5:.1f} eV)")

    print(f"\nGrid sweep (step={args.step} eV)...")

    results = grid_sweep(eigvals, target_char, num_wann, exclude,
                         (fmin_floor, fmax_floor), args.step)

    if not results:
        print("ERROR: No valid configs found!")
        return 1

    n = min(args.top, len(results))
    print(f"\n{'='*110}")
    print(f"TOP {n} CONFIGURATIONS (of {len(results)} valid)")
    print(f"{'='*110}")
    print(f"{'#':>3s} {'Score':>7s} {'Frozen':>14s} {'Outer':>14s} "
          f"{'Nf':>3s} {'Nd':>3s} {'AvgP':>6s} {'Gauge':>5s} {'FeCov':>5s} {'OvEx':>5s} "
          f"{'Frozen bands':>30s}")
    print('-' * 120)
    for i, (loss, fmin, fmax, dmin, dmax, det) in enumerate(results[:n]):
        print(f"{i+1:3d} {det['score']:7.4f} [{fmin:5.1f},{fmax:5.1f}] "
              f"[{dmin:5.1f},{dmax:5.1f}] "
              f"{det['nf']:3d} {det['nd']:3d} "
              f"{det['avg_proj']:6.4f} {det['gauge']:5.3f} "
              f"{det.get('fermi_cov',0):5.3f} "
              f"{det.get('overextend',0):5.3f} "
              f"{str(det['frozen']):>30s}")

    if args.write:
        import shutil

        excl_1idx = sorted(ib + 1 for ib in exclude)
        outdir = f'{d}/runs_opt'
        os.makedirs(outdir, exist_ok=True)

        # Read base .win for structural blocks
        base_win = ''
        win_path = f'{d}/{case}.win'
        if os.path.exists(win_path):
            with open(win_path) as f:
                base_win = f.read()
            # Strip all parameter lines (we'll write fresh ones)
            for k in ['num_wann','num_bands','dis_froz_min','dis_froz_max',
                       'dis_win_min','dis_win_max','exclude_bands',
                       'dis_num_iter','dis_conv_tol','dis_conv_window',
                       'num_iter','conv_tol','conv_window','restart']:
                base_win = re.sub(rf'^{k}\s*=.*\n', '', base_win, flags=re.MULTILINE)
            base_win = re.sub(r'! ===.*?! === end.*?===\s*\n?', '', base_win, flags=re.DOTALL)

            # Replace projections block with correct target orbitals
            struct_path = f'{d}/{case}.struct'
            if os.path.exists(struct_path):
                with open(struct_path) as sf:
                    slines = sf.readlines()
                # Parse atom positions (fractional)
                nneq_s = int(re.search(r'ATOMS[:\s]+(\d+)', slines[1]).group(1))
                atom_pos = []
                ii = 4
                for iat in range(nneq_s):
                    pos = slines[ii]
                    x = float(pos[12:22])
                    y = float(pos[25:35])
                    z = float(pos[38:48])
                    ii += 1
                    mult = int(slines[ii].split('=')[1].split()[0])
                    ii += 1
                    for _ in range(mult - 1):
                        ii += 1
                    z_val = int(float(re.search(r'Z:\s*([\d.]+)', slines[ii]).group(1)))
                    ii += 1 + 3
                    atom_pos.append((x, y, z, mult, z_val))

                # Generate projections matching num_wann
                l_names = {0: 's', 1: 'p', 2: 'd', 3: 'f'}
                proj_lines = ['begin projections']
                for x, y, z, mult, Z in atom_pos:
                    for l in valence_config(Z):
                        for mu in range(mult):
                            proj_lines.append(f'  f={x:.8f},{y:.8f},{z:.8f}:{l_names[l]}')
                proj_lines.append('end projections')
                proj_block = '\n'.join(proj_lines)

                base_win = re.sub(r'begin projections.*?end projections',
                                  proj_block, base_win, flags=re.DOTALL)

        config_names = []
        for i, (loss, fmin, fmax, dmin, dmax, det) in enumerate(results[:n]):
            dirname = f'opt_{i+1:02d}'
            config_names.append(dirname)
            od = f'{outdir}/{dirname}'
            os.makedirs(od, exist_ok=True)

            # Copy shared files
            for fn in [f'{case}.amn', f'{case}.mmn', f'{case}.eig',
                       f'{case}.nnkp', f'{case}.struct', f'{case}.fermi',
                       f'{case}.klist_band', f'{case}.spaghetti_ene']:
                src_path = f'{d}/{fn}'
                if os.path.exists(src_path):
                    shutil.copy(src_path, od)

            # Write .win with optimizer-determined parameters
            header = f"""! === {dirname}: {det['nf']} frozen {det['frozen']}, score={det['score']:.4f} ===
num_wann        = {num_wann}
num_bands       = {nb}
dis_froz_min    = {fmin:.6f}
dis_froz_max    = {fmax:.6f}
dis_win_min     = {dmin:.6f}
dis_win_max     = {dmax:.6f}
exclude_bands   = {','.join(str(b) for b in excl_1idx)}
dis_num_iter    = 10000
dis_conv_tol    = 1.0e-12
dis_conv_window = 5
num_iter        = 10000
conv_tol        = 1.0e-12
conv_window     = 5
! === end ===
"""
            with open(f'{od}/{case}.win', 'w') as f:
                f.write(header + '\n' + base_win)

            with open(f'{od}/APPROACH.txt', 'w') as f:
                f.write(f'{dirname}: {det["nf"]} frozen {det["frozen"]}\n')
                f.write(f'frozen=[{fmin:.1f},{fmax:.1f}] outer=[{dmin:.1f},{dmax:.1f}]\n')
                f.write(f'exclude={excl_1idx}\n')
                f.write(f'score={det["score"]:.4f}\n')

            print(f"  {dirname}: frozen={det['frozen']} [{fmin:.1f},{fmax:.1f}] "
                  f"excl={excl_1idx} score={det['score']:.4f}")

        # Write SLURM script
        names_str = ' '.join(config_names)
        with open(f'{outdir}/run_opt.sh', 'w') as f:
            f.write(f"""#!/bin/bash --login
#SBATCH -J qtl_opt
#SBATCH -o qtl_opt-%J.o
#SBATCH --ntasks=25
#SBATCH --cpus-per-task=1
#SBATCH -N 1
#SBATCH --mem-per-cpu=4G
#SBATCH -p general,mendoza_q
#SBATCH -t 4:00:00

module purge
module load Wannier90/3.1.0-intel-2024a

DIR=$SLURM_SUBMIT_DIR

for config in {names_str}; do
    echo "========================================"
    echo "Running config: $config"
    cat "$DIR/$config/APPROACH.txt"
    echo "========================================"
    cd "$DIR/$config"
    mpirun -n $SLURM_NTASKS wannier90.x {case}
    echo "Finished: $config (exit code: $?)"
    echo ""
done

echo "All configs complete."
echo "=== Summary ==="
for config in {names_str}; do
    omega=$(grep "Omega Total" "$DIR/$config/{case}.wout" 2>/dev/null | tail -1 | awk '{{print $NF}}')
    delta=$(grep "CONV" "$DIR/$config/{case}.wout" 2>/dev/null | tail -1 | awk '{{print $2}}')
    echo "$config: Omega_Total=$omega  Delta=$delta"
done
""")
        print(f"\n  SLURM script: {outdir}/run_opt.sh")
        print(f"  Transfer {outdir}/ and: cd runs_opt && sbatch run_opt.sh")

    return 0


if __name__ == '__main__':
    sys.exit(main())
