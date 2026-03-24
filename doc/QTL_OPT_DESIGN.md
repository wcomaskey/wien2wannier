# QTL Window Optimization via Subspace Search

## Motivation

The current QTL pipeline treats `num_wann` as a hard constraint on the frozen window: at every k-point, frozen bands ≤ num_wann. This limits the frozen window size and prevents us from finding the optimal subspace.

The real optimization problem is: **find the energy windows and band selection that maximize band structure accuracy within the frozen region while maintaining well-localized Wannier functions.**

## Problem Formulation

### Variables
- `dis_froz_min`, `dis_froz_max`: frozen window edges (continuous)
- `dis_win_min`, `dis_win_max`: outer window edges (continuous)
- `exclude_bands`: set of excluded band indices (discrete)

### Objective
Minimize a combined cost function:

```
L = w_spread * Omega_Total + w_accuracy * E_rms_frozen + w_coverage * (1 / frozen_width)
```

Where:
- `Omega_Total`: Wannier spread (from Wannier90)
- `E_rms_frozen`: RMS band deviation within frozen window (from band comparison)
- `frozen_width`: width of frozen window in eV
- `w_spread`, `w_accuracy`, `w_coverage`: tunable weights

### Constraints
- `num_wann` = number of target Wannier functions (from valence config)
- At each k-point: bands in frozen window ≤ bands in outer window
- Frozen window ⊂ outer window
- Minimum frozen floor: guaranteed minimum coverage around E_F

## Optimization Strategy

### Approach: Grid Search + Local Refinement

Full SGD is impractical because each evaluation requires running Wannier90 (~minutes to hours). Instead:

**Phase 1: Coarse grid (no Wannier90 runs)**
- Enumerate frozen window candidates using QTL data
- For each candidate, compute:
  - QTL projectability score within frozen window
  - Band manifold completeness (shell completeness)
  - num_wann constraint satisfaction
  - Estimated Omega_I from disentanglement theory
- Rank candidates by predicted quality
- Select top 5-10 candidates

**Phase 2: Wannier90 evaluation (parallel)**
- Run Wannier90 for each candidate (can be parallelized on cluster)
- Compute actual Omega and band accuracy
- Identify the Pareto frontier (spread vs accuracy tradeoff)

**Phase 3: Local refinement**
- Take the best candidate from Phase 2
- Perturb frozen window edges by ±0.5 eV
- Re-run Wannier90 for the perturbed configs
- Select the configuration that minimizes the cost function

### Predicted Quality Score (Phase 1)

Without running Wannier90, estimate quality from:

```python
def predicted_score(fmin, fmax, eigvals, qtl_char, num_wann):
    # 1. Average QTL target character of bands in frozen window
    frozen_bands = bands_fully_in_window(fmin, fmax)
    avg_qtl = mean(qtl_target[frozen_bands])

    # 2. Manifold completeness: do frozen bands form complete shells?
    shell_score = count_complete_shells(frozen_bands) / count_partial_shells(frozen_bands)

    # 3. Window width relative to total bandwidth
    coverage = (fmax - fmin) / total_bandwidth

    # 4. Constraint margin: how close to violating num_wann?
    max_bands = max_bands_in_window(fmin, fmax)
    margin = (num_wann - max_bands) / num_wann  # negative = violation

    # 5. Disentangle band count: more = more gauge freedom
    n_disent = count_bands_in_outer_not_frozen()
    freedom = n_disent / num_wann

    # Combined score (higher = better)
    if margin < 0:
        return -inf  # hard constraint violation
    return avg_qtl * shell_score * coverage * (1 + freedom)
```

### Relaxing the num_wann Constraint

Instead of treating `max_bands_in_frozen ≤ num_wann` as a hard wall, treat it as a soft penalty:

```python
if max_bands > num_wann:
    # Shrink frozen window from conduction side until satisfied
    # BUT: record the "ideal" window before shrinking
    # The optimization can explore whether excluding specific
    # high-energy bands (via exclude_bands) allows a wider
    # frozen window that would otherwise violate the constraint
    penalty = (max_bands - num_wann) * penalty_weight
```

This allows the optimizer to discover that excluding a few high-energy bands can enable a much wider frozen window — a tradeoff the current hard constraint prevents.

## Implementation Plan

### File: `SRC/qtl_opt.py`

Standalone Python script (not part of the pipeline — called separately):

```bash
# After running qtl_lapw -norun to generate base files:
python3 qtl_opt.py CASE --qtl case.qtl --eig case.eig --num-wann 12

# Outputs:
#   case.opt_configs/config_001/WANN.win
#   case.opt_configs/config_002/WANN.win
#   ...
#   case.opt_configs/run_all.sh  (SLURM script)
#   case.opt_configs/summary.json
```

### Inputs
- `case.qtl`: orbital character
- `case.eig`: eigenvalues
- `case.amn`, `case.mmn`, `case.nnkp`: shared across all configs
- `num_wann`: from valence config or user override

### Outputs
- Directory of Wannier90 configs ready for parallel cluster execution
- Summary JSON with predicted scores for each config
- After cluster run: comparison script to analyze results and select optimal

## Future Extensions

- **Bayesian optimization**: Replace grid search with GP-based acquisition function to minimize Wannier90 evaluations
- **Transfer learning**: Use results from one material to warm-start optimization for similar materials
- **Active learning**: Iteratively select the most informative config to evaluate next
