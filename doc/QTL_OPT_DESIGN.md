# QTL Window Optimization via Subspace Search

## Motivation

The current QTL pipeline uses percentile-based thresholds to classify bands, then derives energy windows from the classification. This works but has limitations:

1. The `num_wann` constraint on the frozen window forces aggressive shrinking, often producing tiny frozen windows (e.g. p80 for MgB2 gives [-0, +2] eV)
2. Window selection is decoupled from the actual optimization target (band accuracy)
3. The percentile approach doesn't naturally balance window coverage against projectability

The real problem is an **optimization over the window parameters** that directly maximizes band structure accuracy within the region of interest.

## Problem Formulation

### Decision Variables
- `dis_froz_min`, `dis_froz_max`: frozen window edges (continuous, eV)
- `dis_win_min`, `dis_win_max`: outer window edges (continuous, eV)
- `exclude_bands`: set of excluded band indices (discrete)

### Loss Function

```
L = -w_p * P_score - w_c * C_score + w_v * V_penalty
```

**Minimize L** (equivalently, maximize projectability and coverage while satisfying constraints).

**P_score (Projectability)**: Average QTL target character of bands within the frozen window, weighted by k-point sampling:

```python
P_score = mean over k-points of:
    sum(qtl_target[b,k] for b in frozen_bands_at_k) / n_frozen_at_k
```

Higher P_score means the frozen bands have stronger target orbital character — the Wannier functions will be better localized on the intended atoms/orbitals.

**C_score (Coverage)**: Frozen window width normalized by total bandwidth of interest. Uses an asymmetric weighting that values coverage below E_F more than above (occupied states matter more):

```python
C_score = (|fmin| + 0.5 * fmax) / (|E_deepest_valence| + 0.5 * E_highest_conduction)
```

This prevents the optimizer from choosing tiny windows with perfect projectability but no physical coverage.

**V_penalty (Constraint violations)**: Soft penalties for violating Wannier90 requirements:

```python
V_penalty = 0
# Hard floor: frozen window must be at least [-W_min_below, +W_min_above]
if fmin > -W_min_below: V_penalty += alpha * (fmin + W_min_below)**2
if fmax < W_min_above:  V_penalty += alpha * (W_min_above - fmax)**2

# num_wann constraint: at every k, bands in frozen ≤ num_wann
for each k-point:
    n_in = count(active bands in [fmin, fmax] at k)
    if n_in > num_wann:
        V_penalty += beta * (n_in - num_wann)**2

# Must have disentanglement: num_bands_in_outer > num_wann
for each k-point:
    n_outer = count(active bands in [dmin, dmax] at k)
    if n_outer <= num_wann:
        V_penalty += gamma * (num_wann + 1 - n_outer)**2

# Outer must strictly contain frozen
if fmin <= dmin: V_penalty += delta * (dmin - fmin + 1)**2
if fmax >= dmax: V_penalty += delta * (fmax - dmax + 1)**2
```

### Hard Constraints (non-negotiable)
1. **Minimum frozen floor**: `fmin ≤ -W_min_below` and `fmax ≥ +W_min_above` (user-configurable, default: 4 eV below, 2 eV above E_F)
2. **Disentanglement required**: At every k-point, the number of active bands in the outer window must be **strictly greater than** `num_wann`. No exceptions — this ensures gauge freedom exists.
3. **Frozen ⊂ Outer**: The frozen window must be strictly inside the outer window with margin ≥ 0.5 eV on each side.

### Soft Constraints (penalized)
1. **num_wann in frozen**: At every k-point, active bands in frozen ≤ `num_wann`. Violated configurations are heavily penalized but not immediately rejected — the optimizer can explore whether `exclude_bands` modifications resolve the violation.

## Optimization Strategy

### Phase 1: Predicted Score Grid (no Wannier90, seconds)

Sweep `fmin` and `fmax` on a grid (0.5 eV steps) within the eigenvalue range. For each (fmin, fmax):

```python
def evaluate_config(fmin, fmax, dmin, dmax, excl, eigvals, qtl, num_wann):
    # Compute P_score from QTL data
    P = projectability_score(fmin, fmax, eigvals, qtl, excl)

    # Compute C_score from window geometry
    C = coverage_score(fmin, fmax, E_deepest, E_highest)

    # Compute V_penalty from constraint checks
    V = constraint_penalty(fmin, fmax, dmin, dmax, excl, eigvals, num_wann)

    # Loss
    L = -w_p * P - w_c * C + w_v * V
    return L, P, C, V
```

This produces a 2D landscape of L(fmin, fmax). Find the global minimum and top-N local minima.

### Phase 2: Exclude-Band Exploration (no Wannier90, seconds)

For each top candidate from Phase 1 that has num_wann violations:
- Try excluding high-energy bands one at a time
- Check if the violation is resolved
- If yes, record the modified config with its score

This discovers configs where excluding 1-2 high-energy bands enables a much wider frozen window.

### Phase 3: Wannier90 Evaluation (parallel, minutes-hours)

Generate `.win` files for the top 5-10 configs from Phases 1-2. Run Wannier90 in parallel on the cluster. After completion:
- Read Omega_Total from each `.wout`
- Compute band accuracy from `_band.dat` vs `spaghetti_ene` within the frozen window
- Select the config that minimizes the actual (not predicted) cost function

### Phase 4: Local Refinement (optional)

Take the best Phase 3 result. Perturb fmin and fmax by ±0.25, ±0.5, ±1.0 eV. Run Wannier90 for each perturbation. Select the best.

## Default Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `w_p` | 1.0 | Projectability weight |
| `w_c` | 1.0 | Coverage weight |
| `w_v` | 100.0 | Violation penalty weight |
| `W_min_below` | 4.0 eV | Minimum frozen window below E_F |
| `W_min_above` | 2.0 eV | Minimum frozen window above E_F |
| `alpha` | 10.0 | Floor violation penalty |
| `beta` | 50.0 | num_wann violation penalty |
| `gamma` | 50.0 | Disentanglement violation penalty |
| `delta` | 20.0 | Outer containment penalty |
| `fmin_step` | 0.5 eV | Grid step for Phase 1 |
| `fmax_step` | 0.5 eV | Grid step for Phase 1 |

## Implementation

### File: `SRC/qtl_opt.py`

```bash
python3 qtl_opt.py CASE [options]

Options:
  --qtl FILE        QTL file (default: CASE.qtl)
  --eig FILE        EIG file (default: CASE.eig)
  --num-wann N      number of Wannier functions (default: auto from struct)
  --excl BANDS      semicore bands to exclude (default: auto from QTL)
  --min-below E     minimum frozen below E_F in eV (default: 4.0)
  --min-above E     minimum frozen above E_F in eV (default: 2.0)
  --top-n N         number of configs for Phase 3 (default: 8)
  --slurm           generate SLURM submission script
```

### Outputs
```
CASE_opt/
├── phase1_landscape.png       # 2D loss landscape plot
├── phase1_scores.json         # All grid scores
├── config_001/CASE.win        # Top Wannier90 configs
├── config_002/CASE.win
├── ...
├── run_all.sh                 # SLURM script for Phase 3
└── summary.json               # Config details + predicted scores
```

### Post-Phase 3 Analysis
```bash
python3 qtl_opt.py CASE --analyze
# Reads .wout files from each config, computes actual scores,
# generates comparison plots, recommends optimal config
```

## Design Principles

1. **Projectability and coverage are both essential** — a tiny window with perfect projectability is useless; a huge window with poor projectability gives bad Wannier functions.

2. **Shell completeness is NOT part of the scoring** — it fails for metals/semimetals where bands cross E_F and orbital character is mixed. The QTL projectability score already captures orbital quality without assuming complete shells.

3. **Disentanglement is always required** — `num_bands_in_outer > num_wann` at every k-point. This ensures Wannier90 always has gauge freedom. Configurations without disentanglement are excluded even if they have perfect frozen windows.

4. **The minimum window floor is non-negotiable** — the optimizer cannot choose windows smaller than the floor. This guarantees the Fermi surface is always protected regardless of the projectability landscape.

5. **Phase 1 is cheap** — evaluating the predicted score for thousands of (fmin, fmax) pairs takes seconds. Only Phase 3 (Wannier90 runs) is expensive, and we limit it to the top candidates.
