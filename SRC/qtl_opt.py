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

    # Scoring
    proj = sum(target_char[ib] for ib in frozen) / max(1, num_wann)
    # Coverage: fraction of active energy range covered by frozen window
    active_emin = min(eigvals[ib].min() for ib in active)
    active_emax = max(eigvals[ib].max() for ib in active)
    active_range = active_emax - active_emin
    coverage = min(1.0, (fmax - fmin) / max(1.0, active_range))
    n_dis = min_outer - len(frozen)
    disent = min(1.0, n_dis / num_wann)
    froz_frac = len(frozen) / num_wann

    # Asymmetry penalty: prefer valence side wider
    asym = 0
    if fmax > 0 and fmin < 0:
        ratio = abs(fmin) / max(0.1, fmax)
        if ratio < 1.5:
            asym = 0.1 * (1.5 - ratio)

    score = (1.0 * proj + 0.5 * coverage + 0.3 * disent
             + 0.4 * froz_frac - asym)

    return -score, {
        'frozen': [ib + 1 for ib in frozen],
        'nf': len(frozen), 'nd': n_dis,
        'min_outer': min_outer, 'max_frozen': max_frozen,
        'proj': proj, 'coverage': coverage,
        'disent': disent, 'froz_frac': froz_frac,
        'asym': asym, 'score': score,
    }


def grid_sweep(eigvals, target_char, num_wann, exclude,
               min_fwin=(4.0, 2.0), step=0.5):
    """Sweep frozen window grid, return sorted results."""
    active = [ib for ib in range(eigvals.shape[0]) if ib not in exclude]
    e_lo = min(eigvals[ib].min() for ib in active)
    e_hi = max(eigvals[ib].max() for ib in active)

    dmin = e_lo - 2.0
    dmax = e_hi + 2.0

    results = []
    for fmin in np.arange(e_lo - 1, 1, step):
        for fmax in np.arange(0, e_hi + 1, step):
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
    ap.add_argument('-write', action='store_true')
    args = ap.parse_args()

    d = args.dir
    case = args.case

    qtl_char, qtl_int, meta = parse_qtl(args.qtl or f'{d}/{case}.qtl')
    eigvals = parse_eig(f'{d}/{case}.eig')
    target_char, targets, nw_auto = compute_targets(qtl_char, f'{d}/{case}.struct')
    num_wann = args.nwann or nw_auto
    nb, nk = eigvals.shape

    print(f"System: {nb} bands, {nk} k-pts, num_wann={num_wann}")
    print(f"Targets: {targets}")
    print(f"Min frozen: [-{args.efmin}, +{args.efmax}] eV\n")

    # Exclude core/semicore via energy gaps + orbital character
    e_avg = eigvals.mean(axis=1)
    exclude = set()
    
    # 1. Energy gap detection: exclude bands below gaps > 5 eV
    for ib in range(nb - 1):
        gap = e_avg[ib + 1] - e_avg[ib]
        if e_avg[ib] < 0 and gap > 5.0:
            for jb in range(ib + 1):
                exclude.add(jb)
    
    # 2. Pure d/f semicore (> 90% single-orbital character)
    for ib in range(nb):
        for ja in range(meta['nat']):
            if qtl_char[ib, ja, 2] > 0.90 or qtl_char[ib, ja, 3] > 0.90:
                exclude.add(ib)

    print(f"{'Band':>4s} {'E_avg':>7s} {'target':>7s} {'int':>5s} {'Status':>8s}")
    print('-' * 35)
    for ib in range(nb):
        st = 'EXCL' if ib in exclude else ''
        print(f"{ib+1:4d} {e_avg[ib]:7.1f} {target_char[ib]:7.4f} "
              f"{qtl_int[ib]:5.3f} {st:>8s}")

    print(f"\nExcluded: {sorted(ib+1 for ib in exclude)}")
    print(f"\nGrid sweep (step={args.step} eV)...")

    results = grid_sweep(eigvals, target_char, num_wann, exclude,
                         (args.efmin, args.efmax), args.step)

    if not results:
        print("ERROR: No valid configs found!")
        return 1

    n = min(args.top, len(results))
    print(f"\n{'='*110}")
    print(f"TOP {n} CONFIGURATIONS (of {len(results)} valid)")
    print(f"{'='*110}")
    print(f"{'#':>3s} {'Score':>7s} {'Frozen':>14s} {'Outer':>14s} "
          f"{'Nf':>3s} {'Nd':>3s} {'MxF':>4s} {'Proj':>6s} {'Cov':>5s} "
          f"{'Frozen bands':>30s}")
    print('-' * 110)
    for i, (loss, fmin, fmax, dmin, dmax, det) in enumerate(results[:n]):
        print(f"{i+1:3d} {det['score']:7.4f} [{fmin:5.1f},{fmax:5.1f}] "
              f"[{dmin:5.1f},{dmax:5.1f}] "
              f"{det['nf']:3d} {det['nd']:3d} {det['max_frozen']:4d} "
              f"{det['proj']:6.3f} {det['coverage']:5.3f} "
              f"{str(det['frozen']):>30s}")

    if args.write:
        import shutil
        os.makedirs(f'{d}/runs_opt', exist_ok=True)
        for i, (loss, fmin, fmax, dmin, dmax, det) in enumerate(results[:n]):
            dirname = f'opt_{i+1:02d}'
            od = f'{d}/runs_opt/{dirname}'
            os.makedirs(od, exist_ok=True)
            for fn in [f'{case}.amn', f'{case}.mmn', f'{case}.eig',
                       f'{case}.nnkp', f'{case}.struct', f'{case}.fermi',
                       f'{case}.klist_band', f'{case}.spaghetti_ene']:
                src = f'{d}/{fn}'
                if os.path.exists(src):
                    shutil.copy(src, od)

            excl_1idx = sorted(ib + 1 for ib in exclude)
            with open(f'{od}/config.json', 'w') as f:
                json.dump({
                    'rank': i+1, 'score': det['score'],
                    'fmin': float(fmin), 'fmax': float(fmax),
                    'dmin': float(dmin), 'dmax': float(dmax),
                    'frozen_bands': det['frozen'],
                    'exclude_bands': excl_1idx,
                    'num_wann': num_wann, 'num_bands': nb,
                    'details': {k: (float(v) if isinstance(v, (np.floating, float)) else v)
                                for k, v in det.items()},
                }, f, indent=2)
            print(f"  {dirname}: [{fmin:.1f},{fmax:.1f}] nf={det['nf']} "
                  f"score={det['score']:.4f}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
