# PDWF Future Plans

## Current State (auto_pdwf branch)

The PDWF module automatically generates Wannier functions from a Wien2k calculation using:
- Split-Gaussian relevance weighting (σ_below, σ_above from eigenvalue distribution)
- Percentile-based frozen/disentangle/excluded classification
- Gap-based semicore detection + .in1 STOP LO diagnostics
- Automatic bmin from deep core gap detection (>20 eV)
- Band manifold detection (75th percentile gap threshold)

**Best result**: MgB2 with 7 frozen bands, Omega_Total = 21.20 Ang²

---

## Planned: Shell Completeness at Frozen Boundary

### Motivation
At the frozen/disentangle boundary, borderline bands should be included or excluded based on whether they complete an atomic shell. In systems with well-separated orbital character (transition metal oxides, rare earths), this would prevent splitting a d or f manifold across the boundary.

### Rules
- **No spin**: s=1, p=3, d=5, f=7 orbitals per shell
- **Spin-polarized**: s=2, p=6, d=10, f=14 orbitals per shell
- At the frozen edge: if including a borderline band **completes** a (atom, l) shell → promote to frozen
- At the frozen edge: if a band **starts a new incomplete** shell → demote to disentangle

### Implementation
- New Phase 3c in `classify_bands`, after gap bridging
- Uses per-(atom, l) projectability decomposition from `A_full` matrix
- Only fires when the shell decomposition shows a clear preference (>70% of shell weight on one side)
- For highly hybridized metals (like MgB2), the rule would be inactive since no clean shell boundary exists

### Data Required
- Per-orbital AMN projections (already computed in `pdwf_compute`)
- Shell size lookup from valence_config (already exists)

### Testing
- MgB2 (metal, hybridized): rule should NOT change the boundary
- SrVO3 (correlated oxide, clear d-manifold): rule should complete the t2g/eg shells
- Fe (transition metal): rule should keep full 3d shell together

---

## Planned: Wien2k QTL Integration

### Motivation
Wien2k's `x lapw2 -qtl` computes orbital-resolved partial charges (fat bands) using the full LAPW basis, including interstitial contributions. These are more accurate than our muffin-tin-only projectabilities and could serve as:
1. **Validation** of PDWF projectabilities
2. **Direct input** for band classification (replace Alm-based projectability)
3. **Shell completeness** data for the boundary rule above

### Approach
- Run `x lapw2 -qtl` on the **uniform k-mesh** (same grid as SCF/Wannier)
- Parse `case.qtl` file: partial charges q(atom, l, band, k)
- Use as projectability input to `classify_bands`
- Compare with Alm-based projectability to validate/calibrate thresholds

### Why Uniform Grid (not -band)
- The `-qtl` on uniform grid gives orbital character at the **exact same k-points** used by w2w for AMN/MMN
- Direct k-by-k correlation with PDWF projectability
- Band path only samples 1D slice, misses k-points with different character
- The uniform grid QTL can also be used for Fermi surface orbital analysis

### Integration Options
1. **Optional enrichment**: if `case.qtl` exists, use it to refine thresholds
2. **Primary mode**: add `-qtl` flag to `pdwf_lapw` that runs `x lapw2 -qtl` before w2w
3. **Validation mode**: compare QTL partial charges with Alm projectabilities, report discrepancies

### File Format
```
case.qtl:
  header (atom info, l-channels)
  per k-point: band index, E(Ry), q(atom1,l=0), q(atom1,l=1), ..., q(atom2,l=0), ...
```

---

## Planned: Automatic Configuration Selection

### Problem
The current defaults (p_high=0.60, p_low=0.15) work well for MgB2 but may not be optimal for all materials. Need a method to automatically determine the best percentile thresholds.

### Approaches Under Consideration

1. **QTL-guided**: Use Wien2k QTL data to identify the "natural" frozen manifold where the target orbital character dominates, then set thresholds to match.

2. **Manifold-aware**: Use the band manifold detector to identify connected groups, then choose thresholds that align the frozen boundary with manifold gaps.

3. **Iterative refinement**: Run with default thresholds, check Omega convergence, tighten/loosen if needed. Expensive but guaranteed to converge.

4. **Material-class templates**: Pre-tuned parameters for common material classes (simple metals, semiconductors, transition metal oxides, rare earths).

---

## Planned: Split Gaussian Outer Window

### Motivation
Currently the outer (disentangle) window uses `exclude_bands` to remove low-projectability bands. An alternative: define the outer window edges directly from the split Gaussian weight threshold, eliminating `exclude_bands` entirely.

### Approach
- Outer window: E where weight(E) > w_threshold (e.g., 0.2)
- This gives: dis_win_min = -σ_below × sqrt(-2 ln(w)), dis_win_max = σ_above × sqrt(-2 ln(w))
- The split Gaussian naturally makes the conduction side narrower
- Must verify that ≥ num_wann bands exist at every k-point within the window

### Constraint
- σ_above for the **outer** window must be larger than σ_above for the **frozen** window
- Recommendation: use σ_above_outer = max(σ_above_frozen, value giving ≥ num_wann bands everywhere)

---

## Known Issues

1. **2× sigma regression**: The sigma formula `sigma = |E_deepest|` (not `2×|E_deepest|`) produces the correct classification. The 2× version was tested and reverted — it makes the weighting too flat.

2. **bmin sensitivity**: Using .in1 STOP LO energies for bmin excludes too many bands (semicore). The current approach uses only the first gap > 20 eV (deep core) and lets `classify_bands` handle semicore via its gap detector.

3. **Convergence**: Wannier90 can get stuck in local minima at ~27-29 Ang² before jumping to the global minimum at ~21 Ang². This jump typically occurs around iteration 3000-4000. The `num_iter = 10000` setting accommodates this.
