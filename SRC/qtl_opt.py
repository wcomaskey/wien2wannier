#!/usr/bin/env python3
"""
wien2wannier/SRC/qtl_opt.py

QTL-guided window optimization for Wannier90.

Finds optimal frozen/outer energy windows by optimizing a loss function
that balances projectability (QTL target character), window coverage,
and Wannier90 constraint satisfaction.

Phase 1: Grid sweep over (fmin, fmax) — predicted score, no Wannier90
Phase 2: Exclude-band exploration for constraint-violating candidates
Phase 3: Generate top-N Wannier90 configs for cluster evaluation
Phase 4: (post-cluster) Analyze results, select optimal

Usage:
    qtl_opt.py CASE [options]
    qtl_opt.py CASE --analyze   (after cluster run)
"""

import numpy as np
import argparse
import json
import os
import re
import shutil
import sys

# ============================================================
#  QTL Parser
# ============================================================

def parse_qtl(qtlfile):
    """Parse Wien2k .qtl file → (qtl_char[band,atom,l], qtl_int[band], eig_ry[band], metadata)."""
    with open(qtlfile) as f:
        lines = f.readlines()

    ef_ry = float(re.search(r'FERMI ENERGY=\s*([\d.+-]+)', lines[2]).group(1))
    nat = int(re.search(r'NAT=\s*(\d+)', lines[3]).group(1))
    spin = int(re.search(r'SPIN=\s*(\d+)', lines[3]).group(1))

    # Parse JATOM metadata
    jatom_info = []
    for ja in range(nat):
        line = lines[4 + ja]
        mult = int(re.search(r'MULT=\s*(\d+)', line).group(1))
        isplit = int(re.search(r'ISPLIT=\s*(\d+)', line).group(1))
        jatom_info.append({'mult': mult, 'isplit': isplit})

    # Find band blocks
    band_starts = [i for i, l in enumerate(lines) if l.strip().startswith('BAND')]
    nbands = len(band_starts)
    lines_per_k = nat + 1
    if nbands >= 2:
        nk = (band_starts[1] - band_starts[0] - 1) // lines_per_k
    else:
        nk = (len(lines) - band_starts[0] - 1) // lines_per_k

    # Parse per-band orbital character
    # qtl_char[band, atom, l] where l = 0(s), 1(p), 2(d), 3(f)
    qtl_char = np.zeros((nbands, nat, 4))
    qtl_int = np.zeros(nbands)
    eig_ry = np.zeros((nbands, nk))

    for ib, bs in enumerate(band_starts):
        for ik in range(nk):
            offset = bs + 1 + ik * lines_per_k
            for ja in range(nat):
                parts = lines[offset + ja].split()
                if ik == 0 and ja == 0:
                    eig_ry[ib, ik] = float(parts[0])
                elif ja == 0:
                    eig_ry[ib, ik] = float(parts[0])
                # ISPLIT=2: cols are E, jatom, tot, s, p, d, D-eg, D-t2g, f
                qtl_char[ib, ja, 0] += float(parts[3])   # s
                qtl_char[ib, ja, 1] += float(parts[4])   # p
                qtl_char[ib, ja, 2] += float(parts[5])   # d
                if len(parts) > 8:
                    qtl_char[ib, ja, 3] += float(parts[8])  # f
            # Interstitial line
            int_parts = lines[offset + nat].split()
            qtl_int[ib] += float(int_parts[2]) if len(int_parts) > 2 else 0

        qtl_char[ib] /= nk
        qtl_int[ib] /= nk

    eig_ev = (eig_ry - ef_ry) * 13.605693

    return {
        'qtl_char': qtl_char,
        'qtl_int': qtl_int,
        'eig_ev': eig_ev,
        'ef_ry': ef_ry,
        'nbands': nbands,
        'nk': nk,
        'nat': nat,
        'jatom_info': jatom_info,
    }


# ============================================================
#  Eigenvalue reader
# ============================================================

def parse_eig(eigfile):
    """Parse Wannier90 .eig file → eigvals[band, k]."""
    data = {}
    with open(eigfile) as f:
        for line in f:
            parts = line.split()
            ib, ik, e = int(parts[0]), int(parts[1]), float(parts[2])
            data[(ib, ik)] = e
    nb = max(ib for ib, ik in data)
    nk = max(ik for ib, ik in data)
    eigvals = np.zeros((nb, nk))
    for (ib, ik), e in data.items():
        eigvals[ib - 1, ik - 1] = e
    return eigvals


# ============================================================
#  Valence configuration
# ============================================================

def valence_targets(Z):
    """Return target l-channels for element Z."""
    if Z <= 2:
        return [0]
    elif Z <= 4:
        return [0]
    elif Z <= 10:
        return [0, 1]
    elif Z <= 12:
        return [0, 1]
    elif Z <= 18:
        return [0, 1]
    elif Z <= 20:
        return [0, 1]
    elif Z <= 30:
        return [0, 2]  # transition metals: s + d
    elif Z <= 36:
        return [0, 1]
    elif Z <= 38:
        return [0, 1]
    elif Z <= 48:
        return [0, 2]
    elif Z <= 54:
        return [0, 1]
    elif Z <= 56:
        return [0, 1]
    elif Z <= 71:
        return [0, 2, 3]  # lanthanides: s + d + f
    elif Z <= 80:
        return [0, 2]
    elif Z <= 86:
        return [0, 1]
    else:
        return [0, 1]


def get_target_atoms_l(structfile):
    """Parse struct, return list of (jatom_idx, l) target channels and Z values."""
    with open(structfile) as f:
        lines = f.readlines()

    nneq = int(re.search(r'ATOMS[:\s]+(\d+)', lines[1]).group(1))
    targets = []
    Z_list = []
    mult_list = []
    i = 4  # first ATOM line
    for ja in range(nneq):
        # i = ATOM position line
        # i+1 = MULT line
        mult = int(re.search(r'MULT=\s*(\d+)', lines[i + 1]).group(1))
        mult_list.append(mult)
        # Element line at i + 1 + mult (after MULT + equiv positions)
        elem_line = lines[i + 1 + mult]
        Z = float(re.search(r'Z:\s*([\d.]+)', elem_line).group(1))
        Z_list.append(int(round(Z)))
        for l in valence_targets(int(round(Z))):
            targets.append((ja, l))
        # Advance: position + MULT + (mult-1) equiv positions + element + 3 rot
        i += 1 + 1 + (mult - 1) + 1 + 3

    return targets, Z_list, nneq, mult_list


# ============================================================
#  Scoring functions
# ============================================================

def compute_target_char(qtl_char, targets):
    """Compute per-band target character from QTL data."""
    nbands = qtl_char.shape[0]
    tc = np.zeros(nbands)
    for ja, l in targets:
        tc += qtl_char[:, ja, l]
    return tc


def projectability_score(fmin, fmax, eigvals, target_char, excl_set):
    """Average QTL target character of active bands within frozen window, k-averaged."""
    nb, nk = eigvals.shape
    nb_qtl = len(target_char)
    total_p = 0.0
    total_count = 0
    for ik in range(nk):
        for ib in range(min(nb, nb_qtl)):
            if (ib + 1) in excl_set:
                continue
            if eigvals[ib, ik] >= fmin and eigvals[ib, ik] <= fmax:
                total_p += target_char[ib]
                total_count += 1
    return total_p / max(total_count, 1)


def coverage_score(fmin, fmax, e_deepest, e_highest):
    """Asymmetric coverage score: values coverage below E_F more."""
    val_below = abs(fmin) / max(abs(e_deepest), 1e-6)
    val_above = fmax / max(abs(e_highest), 1e-6)
    return val_below + 0.5 * val_above


def constraint_penalty(fmin, fmax, dmin, dmax, excl_set, eigvals, num_wann,
                       W_min_below, W_min_above, alpha=10.0, beta=50.0,
                       gamma=50.0, delta=20.0):
    """Compute soft constraint penalty."""
    nb, nk = eigvals.shape
    V = 0.0

    # Minimum floor
    if fmin > -W_min_below:
        V += alpha * (fmin + W_min_below) ** 2
    if fmax < W_min_above:
        V += alpha * (W_min_above - fmax) ** 2

    # num_wann constraint (at every k)
    for ik in range(nk):
        n_in = sum(1 for ib in range(nb) if (ib + 1) not in excl_set
                   and eigvals[ib, ik] >= fmin and eigvals[ib, ik] <= fmax)
        if n_in > num_wann:
            V += beta * (n_in - num_wann) ** 2

    # Disentanglement required: num_bands_in_outer > num_wann at every k
    for ik in range(nk):
        n_outer = sum(1 for ib in range(nb) if (ib + 1) not in excl_set
                      and eigvals[ib, ik] >= dmin and eigvals[ib, ik] <= dmax)
        if n_outer <= num_wann:
            V += gamma * (num_wann + 1 - n_outer) ** 2

    # Frozen strictly inside outer
    if fmin <= dmin:
        V += delta * (dmin - fmin + 1) ** 2
    if fmax >= dmax:
        V += delta * (fmax - dmax + 1) ** 2

    return V


def evaluate_config(fmin, fmax, dmin, dmax, excl_set, eigvals, target_char,
                    num_wann, e_deepest, e_highest, W_min_below, W_min_above,
                    w_p=1.0, w_c=1.0, w_v=100.0):
    """Evaluate a single (fmin, fmax) configuration."""
    P = projectability_score(fmin, fmax, eigvals, target_char, excl_set)
    C = coverage_score(fmin, fmax, e_deepest, e_highest)
    V = constraint_penalty(fmin, fmax, dmin, dmax, excl_set, eigvals,
                           num_wann, W_min_below, W_min_above)
    L = -w_p * P - w_c * C + w_v * V

    # Count frozen bands (fully inside window)
    nb = eigvals.shape[0]
    nb_tc = len(target_char)
    frozen = [ib + 1 for ib in range(min(nb, nb_tc)) if (ib + 1) not in excl_set
              and eigvals[ib].min() >= fmin and eigvals[ib].max() <= fmax]

    return {
        'fmin': fmin, 'fmax': fmax,
        'dmin': dmin, 'dmax': dmax,
        'P': P, 'C': C, 'V': V, 'L': L,
        'frozen': frozen, 'n_frozen': len(frozen),
        'excl': sorted(excl_set),
    }


# ============================================================
#  Phase 1: Grid sweep
# ============================================================

def phase1_grid(eigvals, target_char, excl_set, num_wann,
                e_deepest, e_highest, W_min_below, W_min_above,
                fmin_range=None, fmax_range=None, step=0.5,
                w_p=1.0, w_c=1.0, w_v=100.0):
    """Sweep (fmin, fmax) grid, return sorted results."""
    nb, nk = eigvals.shape

    # Default ranges
    if fmin_range is None:
        fmin_range = (e_deepest - 2.0, -W_min_below)
    if fmax_range is None:
        fmax_range = (W_min_above, e_highest + 2.0)

    # Outer window: fixed at full active band range ± 2 eV
    active = [i for i in range(nb) if (i + 1) not in excl_set]
    dmin = eigvals[active].min() - 2.0
    dmax = eigvals[active].max() + 2.0

    results = []
    fmin_vals = np.arange(fmin_range[0], fmin_range[1] + step, step)
    fmax_vals = np.arange(fmax_range[0], fmax_range[1] + step, step)

    for fmin in fmin_vals:
        for fmax in fmax_vals:
            if fmax <= fmin + 1.0:
                continue  # too narrow
            r = evaluate_config(fmin, fmax, dmin, dmax, excl_set, eigvals,
                                target_char, num_wann, e_deepest, e_highest,
                                W_min_below, W_min_above, w_p, w_c, w_v)
            results.append(r)

    results.sort(key=lambda r: r['L'])
    return results, dmin, dmax


# ============================================================
#  Phase 2: Exclude-band exploration
# ============================================================

def phase2_exclude(top_configs, eigvals, target_char, excl_base, num_wann,
                   e_deepest, e_highest, W_min_below, W_min_above,
                   w_p=1.0, w_c=1.0, w_v=100.0):
    """For configs with num_wann violations, try excluding high-energy bands."""
    nb, nk = eigvals.shape
    e_avg = eigvals.mean(axis=1)
    improved = []

    for cfg in top_configs:
        if cfg['V'] < 1e-6:
            improved.append(cfg)
            continue

        # Try excluding bands from the high-energy end
        active_bands = sorted([ib + 1 for ib in range(nb) if (ib + 1) not in excl_base],
                              key=lambda b: -e_avg[b - 1])

        for n_excl in range(1, min(4, len(active_bands) - num_wann)):
            new_excl = excl_base | set(active_bands[:n_excl])

            # Recompute outer window
            active = [i for i in range(nb) if (i + 1) not in new_excl]
            if len(active) <= num_wann:
                continue
            dmin = eigvals[active].min() - 2.0
            dmax = eigvals[active].max() + 2.0

            r = evaluate_config(cfg['fmin'], cfg['fmax'], dmin, dmax,
                                new_excl, eigvals, target_char, num_wann,
                                e_deepest, e_highest, W_min_below, W_min_above,
                                w_p, w_c, w_v)
            if r['V'] < cfg['V']:
                improved.append(r)

    improved.sort(key=lambda r: r['L'])
    return improved


# ============================================================
#  Phase 3: Config generation
# ============================================================

def generate_configs(top_results, case, base_win, eigvals, num_wann, excl_base,
                     outdir, slurm=False, ntasks=25):
    """Generate Wannier90 .win files for top configs."""
    os.makedirs(outdir, exist_ok=True)
    nb = eigvals.shape[0]

    config_names = []
    summary = []

    for i, r in enumerate(top_results):
        name = f"config_{i + 1:03d}"
        d = os.path.join(outdir, name)
        os.makedirs(d, exist_ok=True)
        config_names.append(name)

        # Copy shared files
        for fn in [f'{case}.amn', f'{case}.mmn', f'{case}.eig', f'{case}.nnkp',
                   f'{case}.struct', f'{case}.fermi', f'{case}.klist_band',
                   f'{case}.spaghetti_ene']:
            if os.path.exists(fn):
                shutil.copy(fn, d)

        # Build .win
        content = base_win
        content = re.sub(r'! === .*?! === end.*?===\s*', '', content, flags=re.DOTALL)
        content = re.sub(r'^restart.*\n', '', content, flags=re.MULTILINE)
        for k in ['num_wann', 'num_bands', 'dis_froz_min', 'dis_froz_max',
                   'dis_win_min', 'dis_win_max', 'exclude_bands',
                   'dis_num_iter', 'dis_conv_tol', 'dis_conv_window',
                   'num_iter', 'conv_tol', 'conv_window']:
            content = re.sub(rf'^{k}\s*=.*\n', '', content, flags=re.MULTILINE)

        excl_str = ','.join(str(b) for b in sorted(r['excl']))
        header = f"""! === {name}: L={r['L']:.4f} P={r['P']:.4f} C={r['C']:.4f} V={r['V']:.4f} ===
num_wann        = {num_wann}
num_bands       = {nb}
dis_froz_min    = {r['fmin']:.6f}
dis_froz_max    = {r['fmax']:.6f}
dis_win_min     = {r['dmin']:.6f}
dis_win_max     = {r['dmax']:.6f}
exclude_bands   = {excl_str}
dis_num_iter    = 10000
dis_conv_tol    = 1.0e-10
dis_conv_window = 5
num_iter        = 10000
conv_tol        = 1.0e-8
conv_window     = 5
! === end ===
"""
        with open(os.path.join(d, f'{case}.win'), 'w') as f:
            f.write(header + '\n' + content)

        with open(os.path.join(d, 'APPROACH.txt'), 'w') as f:
            f.write(f"{name}: frozen=[{r['fmin']:.1f},{r['fmax']:.1f}] "
                    f"outer=[{r['dmin']:.1f},{r['dmax']:.1f}]\n")
            f.write(f"frozen_bands={r['frozen']} excl={r['excl']}\n")
            f.write(f"L={r['L']:.4f} P={r['P']:.4f} C={r['C']:.4f} V={r['V']:.4f}\n")

        summary.append({
            'name': name, **{k: (v if not isinstance(v, np.floating) else float(v))
                             for k, v in r.items()}
        })

    # Summary JSON
    with open(os.path.join(outdir, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    # SLURM script
    with open(os.path.join(outdir, 'run_all.sh'), 'w') as f:
        f.write(f"""#!/bin/bash --login
#SBATCH -J qtl_opt
#SBATCH -o qtl_opt-%J.o
#SBATCH --ntasks={ntasks}
#SBATCH --cpus-per-task=1
#SBATCH -N 1
#SBATCH --mem-per-cpu=4G
#SBATCH -p general,mendoza_q
#SBATCH -t 8:00:00

module purge
module load Wannier90/3.1.0-intel-2024a

DIR=$SLURM_SUBMIT_DIR

for config in {' '.join(config_names)}; do
    echo "========================================"
    echo "Running: $config"
    cat "$DIR/$config/APPROACH.txt"
    echo "========================================"
    cd "$DIR/$config"
    mpirun -n $SLURM_NTASKS wannier90.x {case}
    echo "Exit: $?"
    echo ""
done

echo "=== Summary ==="
for config in {' '.join(config_names)}; do
    omega=$(grep "Omega Total" "$DIR/$config/{case}.wout" 2>/dev/null | tail -1 | awk '{{print $NF}}')
    delta=$(grep "CONV" "$DIR/$config/{case}.wout" 2>/dev/null | tail -1 | awk '{{print $2}}')
    echo "$config: Omega=$omega Delta=$delta"
done
""")

    return config_names, summary


# ============================================================
#  Plotting
# ============================================================

def plot_landscape(results, outfile, num_wann):
    """Plot 2D loss landscape from Phase 1 results."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not available, skipping plot")
        return

    # Extract data
    fmins = np.array([r['fmin'] for r in results])
    fmaxs = np.array([r['fmax'] for r in results])
    losses = np.array([r['L'] for r in results])
    violations = np.array([r['V'] for r in results])

    valid = violations < 1e-6
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Loss landscape
    ax = axes[0]
    sc = ax.scatter(fmins[valid], fmaxs[valid], c=losses[valid],
                    cmap='RdYlGn_r', s=15, alpha=0.8)
    ax.scatter(fmins[~valid], fmaxs[~valid], c='gray', s=5, alpha=0.3, label='violated')
    plt.colorbar(sc, ax=ax, label='Loss L')
    ax.set_xlabel('fmin (eV)')
    ax.set_ylabel('fmax (eV)')
    ax.set_title('Loss Landscape (valid configs)')
    ax.legend()

    # Projectability
    ax = axes[1]
    ps = np.array([r['P'] for r in results])
    sc = ax.scatter(fmins[valid], fmaxs[valid], c=ps[valid],
                    cmap='YlGn', s=15, alpha=0.8)
    plt.colorbar(sc, ax=ax, label='P_score')
    ax.set_xlabel('fmin (eV)')
    ax.set_ylabel('fmax (eV)')
    ax.set_title('Projectability Score')

    # Coverage
    ax = axes[2]
    cs = np.array([r['C'] for r in results])
    sc = ax.scatter(fmins[valid], fmaxs[valid], c=cs[valid],
                    cmap='YlOrRd', s=15, alpha=0.8)
    plt.colorbar(sc, ax=ax, label='C_score')
    ax.set_xlabel('fmin (eV)')
    ax.set_ylabel('fmax (eV)')
    ax.set_title('Coverage Score')

    plt.suptitle(f'QTL Optimization Landscape (num_wann={num_wann})', fontsize=14)
    plt.tight_layout()
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"  Landscape plot: {outfile}")


# ============================================================
#  Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(description='QTL window optimization')
    parser.add_argument('case', help='Case name')
    parser.add_argument('--qtl', default=None, help='QTL file (default: CASE.qtl)')
    parser.add_argument('--eig', default=None, help='EIG file (default: CASE.eig)')
    parser.add_argument('--struct', default=None, help='Struct file (default: CASE.struct)')
    parser.add_argument('--num-wann', type=int, default=0, help='num_wann (default: auto)')
    parser.add_argument('--excl', default='', help='Semicore bands to exclude (comma-separated)')
    parser.add_argument('--min-below', type=float, default=4.0, help='Min frozen below E_F (eV)')
    parser.add_argument('--min-above', type=float, default=2.0, help='Min frozen above E_F (eV)')
    parser.add_argument('--top-n', type=int, default=8, help='Top configs for Phase 3')
    parser.add_argument('--step', type=float, default=0.5, help='Grid step (eV)')
    parser.add_argument('--w-p', type=float, default=1.0, help='Projectability weight')
    parser.add_argument('--w-c', type=float, default=1.0, help='Coverage weight')
    parser.add_argument('--w-v', type=float, default=100.0, help='Violation penalty weight')
    parser.add_argument('--analyze', action='store_true', help='Post-cluster analysis mode')
    args = parser.parse_args()

    case = args.case
    qtlfile = args.qtl or f'{case}.qtl'
    eigfile = args.eig or f'{case}.eig'
    structfile = args.struct or f'{case}.struct'
    outdir = f'{case}_opt'

    if args.analyze:
        print("Post-cluster analysis not yet implemented")
        return

    # Parse inputs
    print(f"QTL Window Optimization: {case}")
    print(f"{'=' * 60}")

    print(f"\nParsing QTL: {qtlfile}")
    qtl = parse_qtl(qtlfile)
    print(f"  {qtl['nbands']} bands, {qtl['nk']} k-points, {qtl['nat']} atoms")

    print(f"Parsing EIG: {eigfile}")
    eigvals = parse_eig(eigfile)
    nb, nk = eigvals.shape
    print(f"  {nb} bands, {nk} k-points")

    # Target channels
    targets, Z_list, nneq, mult_list = get_target_atoms_l(structfile)
    target_char = compute_target_char(qtl['qtl_char'][:nb], targets)
    print(f"Targets: {[(f'atom{ja+1}(Z={Z_list[ja]})', 'spdf'[l]) for ja, l in targets]}")

    # num_wann
    if args.num_wann > 0:
        num_wann = args.num_wann
    else:
        num_wann = sum((2 * l + 1) * mult_list[ja]
                       for ja, l in targets)
    print(f"num_wann = {num_wann}")

    # Semicore exclusion
    if args.excl:
        excl_set = set(int(b) for b in args.excl.split(','))
    else:
        # Auto-detect from QTL: any single (atom, l) character > 0.90
        # This catches d-semicore (Ga 3d), p-semicore (Mg 2p), s-core (B 1s)
        excl_set = set()
        for ib in range(min(nb, qtl['nbands'])):
            for ja in range(qtl['nat']):
                for l in range(4):
                    if qtl['qtl_char'][ib, ja, l] > 0.90:
                        excl_set.add(ib + 1)
            # Also exclude by interstitial fraction (very low = core state)
            if qtl['qtl_int'][ib] < 0.02 and ib + 1 not in excl_set:
                excl_set.add(ib + 1)
    print(f"Excluded (semicore): {sorted(excl_set)}")

    # Energy ranges — use min/max over ALL k-points, not just averages
    active = [i for i in range(nb) if (i + 1) not in excl_set]
    e_avg = eigvals.mean(axis=1)
    e_deepest = eigvals[active].min()   # deepest eigenvalue at any k
    e_highest = eigvals[active].max()   # highest eigenvalue at any k
    print(f"Active energy range: [{e_deepest:.1f}, {e_highest:.1f}] eV")
    print(f"Min frozen floor: [-{args.min_below}, +{args.min_above}] eV")

    # Phase 1
    print(f"\n{'=' * 60}")
    print(f"Phase 1: Grid sweep (step={args.step} eV)")
    results, dmin, dmax = phase1_grid(
        eigvals, target_char, excl_set, num_wann,
        e_deepest, e_highest, args.min_below, args.min_above,
        step=args.step, w_p=args.w_p, w_c=args.w_c, w_v=args.w_v)

    n_valid = sum(1 for r in results if r['V'] < 1e-6)
    print(f"  {len(results)} configs evaluated, {n_valid} valid (no violations)")

    # Top valid results
    valid_results = [r for r in results if r['V'] < 1e-6]
    print(f"\n  Top {min(10, len(valid_results))} valid configs:")
    print(f"  {'#':>3s} {'L':>8s} {'P':>6s} {'C':>6s} {'fmin':>6s} {'fmax':>6s} {'#F':>3s} {'Frozen':>25s}")
    for i, r in enumerate(valid_results[:10]):
        print(f"  {i + 1:3d} {r['L']:8.4f} {r['P']:6.4f} {r['C']:6.4f} "
              f"{r['fmin']:6.1f} {r['fmax']:6.1f} {r['n_frozen']:3d} {str(r['frozen']):>25s}")

    # Phase 2
    violated = [r for r in results[:50] if r['V'] > 1e-6]
    if violated:
        print(f"\n{'=' * 60}")
        print(f"Phase 2: Exclude-band exploration ({len(violated)} violated configs)")
        improved = phase2_exclude(
            violated, eigvals, target_char, excl_set, num_wann,
            e_deepest, e_highest, args.min_below, args.min_above,
            args.w_p, args.w_c, args.w_v)
        new_valid = [r for r in improved if r['V'] < 1e-6]
        print(f"  {len(new_valid)} new valid configs found via band exclusion")
        if new_valid:
            print(f"\n  Top from Phase 2:")
            for i, r in enumerate(new_valid[:5]):
                extra = set(r['excl']) - excl_set
                print(f"  {i + 1:3d} L={r['L']:.4f} P={r['P']:.4f} C={r['C']:.4f} "
                      f"[{r['fmin']:.1f},{r['fmax']:.1f}] +excl={sorted(extra)}")
        valid_results.extend(new_valid)
        valid_results.sort(key=lambda r: r['L'])

    # Phase 3: Generate configs
    print(f"\n{'=' * 60}")
    top_n = min(args.top_n, len(valid_results))
    print(f"Phase 3: Generating top {top_n} Wannier90 configs")

    # Deduplicate by (fmin, fmax, excl)
    seen = set()
    deduped = []
    for r in valid_results:
        key = (round(r['fmin'], 1), round(r['fmax'], 1), tuple(r['excl']))
        if key not in seen:
            seen.add(key)
            deduped.append(r)
    deduped = deduped[:top_n]

    base_win = open(f'{case}.win').read() if os.path.exists(f'{case}.win') else ''
    config_names, summary = generate_configs(
        deduped, case, base_win, eigvals, num_wann, excl_set, outdir)

    print(f"  Generated {len(config_names)} configs in {outdir}/")

    # Plot
    plot_landscape(results, os.path.join(outdir, 'phase1_landscape.png'), num_wann)

    # Final summary
    print(f"\n{'=' * 60}")
    print(f"Summary")
    print(f"{'=' * 60}")
    print(f"\n  {'Config':<12s} {'L':>8s} {'P':>6s} {'C':>6s} {'Frozen':>16s} {'#F':>3s} {'Excl':>15s}")
    print(f"  {'-' * 70}")
    for s in summary:
        print(f"  {s['name']:<12s} {s['L']:8.4f} {s['P']:6.4f} {s['C']:6.4f} "
              f"[{s['fmin']:.1f},{s['fmax']:.1f}]{'':<4s} {s['n_frozen']:3d} {str(s['excl']):>15s}")

    print(f"\nNext: transfer {outdir}/ to cluster and run:")
    print(f"  cd {outdir} && sbatch run_all.sh")
    print(f"Then: python3 qtl_opt.py {case} --analyze")


if __name__ == '__main__':
    main()
