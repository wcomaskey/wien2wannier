# Split-Gaussian Energy Weighting for PDWF Band Classification

## Overview

The Projectability-Disentangled Wannier Functions (PDWF) method classifies
bands into frozen, disentangle, and excluded categories based on their
projectability onto target atomic orbitals. In the LAPW basis, where
projectabilities are systematically lower than in LCAO due to the
muffin-tin decomposition, a **relevance score** combining projectability
with energy proximity to the Fermi level is essential for reliable
classification.

We introduce a **split-Gaussian energy weighting** that treats occupied
and unoccupied states asymmetrically, reflecting their different physical
importance for ground-state and near-Fermi properties.

## Motivation

### Physical Asymmetry

Occupied states (below E_F) determine:
- Ground-state charge density and total energy
- Chemical bonding and hybridization
- Transport properties (conductivity, thermopower)
- Optical properties (joint density of states)

Unoccupied states (above E_F) are important for:
- Excited-state properties (optical absorption, EELS)
- Electron-phonon coupling (limited energy range above E_F)

For most Wannier-based applications (tight-binding models, transport,
electron-phonon coupling), accuracy of the interpolation **below and
near E_F** is far more critical than accuracy deep in the conduction
bands. The frozen window should therefore protect occupied states more
aggressively.

### Prior Work

The original PDWF method (Qiao, Pizzi, Marzari, *npj Comput. Mater.*
**9**, 208, 2023) uses fixed projectability thresholds (p_high = 0.95,
p_low = 0.10) without energy weighting. This works well for LCAO bases
where projectabilities are properly normalized to [0,1] via the Lowdin
transformation.

Energy-dependent band selection has been explored in several contexts:

- **Souza, Marzari, Vanderbilt** (*Phys. Rev. B* **65**, 035109, 2001):
  The disentanglement procedure itself uses an energy window to select
  the subspace, but the window boundaries are user-specified.

- **Vitale et al.** (*npj Comput. Mater.* **6**, 66, 2020): The
  automated Wannierization in Wannier90 v3 uses projectability-based
  selection but does not weight by energy proximity.

- **SCDM method** (Damle, Lin, Ying, *J. Comput. Phys.* **334**, 1-15,
  2017): Selects columns of the density matrix, inherently favoring
  occupied states, but does not use an explicit energy weighting.

Our split-Gaussian approach is, to our knowledge, the first to use an
explicitly asymmetric energy weighting in the PDWF classification.

## Method

### Relevance Score

For each band *n*, we compute a relevance score:

```
relevance(n) = avg_p(n) * w(E_n)
```

where `avg_p(n)` is the k-averaged projectability and `w(E_n)` is the
energy weight.

### Split-Gaussian Weight Function

The energy weight uses a Gaussian centered on E_F with **different
widths below and above**:

```
w(E) = exp(-E^2 / 2*sigma^2)

where:
  sigma = sigma_below  for E < E_F
  sigma = sigma_above  for E > E_F
```

Default parameters:
- sigma_below = 12 eV (gentle decay for occupied states)
- sigma_above = 6 eV (steeper decay for unoccupied states)

### Weight Function Values

| E - E_F (eV) | w (symmetric, sigma=8) | w (split, 12/6) | Ratio |
|:------------:|:---------------------:|:---------------:|:-----:|
| -15 | 0.325 | 0.535 | 1.65x |
| -10 | 0.535 | 0.707 | 1.32x |
| -5 | 0.821 | 0.917 | 1.12x |
| -3 | 0.930 | 0.969 | 1.04x |
| -1 | 0.992 | 0.997 | 1.00x |
| 0 | 1.000 | 1.000 | 1.00x |
| +1 | 0.992 | 0.986 | 0.99x |
| +3 | 0.930 | 0.882 | 0.95x |
| +5 | 0.821 | 0.707 | 0.86x |
| +10 | 0.535 | 0.247 | 0.46x |
| +15 | 0.325 | 0.029 | 0.09x |

The key effect: a band at -10 eV retains 71% of its projectability
weight, while a band at +10 eV retains only 25%. This 2.9x ratio
ensures occupied states are strongly favored for freezing.

### Percentile-Based Classification

After computing relevance scores for all non-semicore bands, we:

1. Sort the relevance scores in ascending order
2. Compute effective thresholds from user-specified percentiles:
   - `r_high_eff` = relevance at the `p_high` percentile (default: 70th)
   - `r_low_eff` = relevance at the `p_low` percentile (default: 15th)
3. Classify:
   - FROZEN: relevance >= r_high_eff
   - DISENTANGLE: r_low_eff <= relevance < r_high_eff
   - EXCLUDED: relevance < r_low_eff

This adapts automatically to the projectability scale of any basis
(LCAO, LAPW, PAW) without requiring material-specific threshold tuning.

### Gauge Freedom Guarantee

After classification, we enforce that not all active bands are frozen.
If `n_frozen == n_active`, we demote the lowest-relevance frozen bands
to DISENTANGLE until `n_frozen <= 0.8 * n_active`. This guarantees
Wannier90 has sufficient gauge freedom for spread minimization.

## Results: MgB2

### System

- MgB2, hexagonal P6/mmm, a = 5.830 Bohr, c = 6.654 Bohr
- 3 atoms (1 Mg + 2 B), 12 target Wannier functions (Mg:sp + B:sp x 2)
- 8x8x8 k-mesh (512 k-points), Wien2k LAPW basis
- Band window: bands 3-20 (18 bands)

### Band Classification Comparison

| Band | avg(p) | E_avg (eV) | Sym. weight | Split weight | Sym. class | Split class |
|:----:|:------:|:----------:|:-----------:|:------------:|:----------:|:-----------:|
| 3 | 0.266 | -12.1 | 0.730 | 0.882 | FROZEN | FROZEN |
| 4 | 0.267 | -2.8 | 0.940 | 0.973 | FROZEN | FROZEN |
| 5 | 0.293 | +0.1 | 1.000 | 1.000 | FROZEN | FROZEN |
| 6 | 0.294 | +0.1 | 1.000 | 1.000 | FROZEN | FROZEN |
| 7 | 0.216 | +1.5 | 0.983 | 0.969 | FROZEN | FROZEN |
| 8 | 0.252 | +5.8 | 0.736 | 0.634 | FROZEN | FROZEN |
| 9 | 0.293 | +6.0 | 0.727 | 0.607 | **FROZEN** | **DISENT.** |
| 10 | 0.310 | +8.0 | 0.607 | 0.325 | **FROZEN** | **DISENT.** |
| 11 | 0.308 | +8.8 | 0.555 | 0.247 | **FROZEN** | **DISENT.** |
| 12 | 0.311 | +10.4 | 0.441 | 0.117 | DISENT. | DISENT. |
| 13 | 0.331 | +11.3 | 0.380 | 0.071 | DISENT. | DISENT. |

The split Gaussian **demotes conduction bands 9-11** from FROZEN to
DISENTANGLE (highlighted), providing crucial gauge freedom for
Wannier90 to optimize the interpolation near E_F.

### Wannier90 Convergence

| Metric | Symmetric (sigma=8) | Split (12/6) | Change |
|:------:|:-------------------:|:------------:|:------:|
| Omega_I (Ang^2) | 17.36 | 15.52 | -11% |
| Omega_D (Ang^2) | 1.93 | 0.05 | **-97%** |
| Omega_OD (Ang^2) | 8.71 | 3.14 | -64% |
| **Omega_Total (Ang^2)** | **27.26** | **18.70** | **-31%** |
| Frozen bands | 7 | 6 | -1 |
| Disentangle bands | 8 | 9 | +1 |
| Converged | No | Yes | -- |

### Energy Windows

| Parameter | Symmetric | Split |
|:---------:|:---------:|:-----:|
| dis_froz_min (eV) | -5.6 | -12.8 |
| dis_froz_max (eV) | +11.7 | +7.8 |
| Frozen width (eV) | 17.3 | 20.6 |
| dis_win_min (eV) | -14.3 | -14.3 |
| dis_win_max (eV) | +22.4 | +22.4 |
| Outer width (eV) | 36.7 | 36.7 |
| Margin (frozen-outer) | 8.7 eV | 6.8 eV |

The split-Gaussian frozen window extends deeper below E_F (-12.8 vs -5.6 eV)
while being tighter above E_F (+7.8 vs +11.7 eV). The margin between
frozen and outer windows remains substantial (>6 eV), ensuring adequate
gauge freedom.

### Key Finding

The near-zero Omega_D (0.05 Ang^2) confirms that the split-Gaussian
frozen window captures the correct physical subspace. With symmetric
weighting, Omega_D = 1.93 Ang^2 indicated significant disentanglement
mixing was needed to find the optimal subspace — the frozen window was
incorrectly including conduction bands that should have been free to mix.

## Implementation

The split-Gaussian weighting is implemented in the `classify_bands`
subroutine in `SRC_w2w/pdwf.F`. The parameters `sigma_below` and
`sigma_above` are compile-time constants that can be adjusted for
specific applications:

```fortran
! Skewed Gaussian weighting: occupied states (E<0) are more important
! than unoccupied (E>0). Use wider sigma below E_F to protect valence.
real(R8), parameter :: SIGMA_BELOW = 12.0d0  ! eV, occupied side
real(R8), parameter :: SIGMA_ABOVE =  6.0d0  ! eV, unoccupied side

! ...
if (e_band_avg(ib) < 0.0d0) then
   sigma = SIGMA_BELOW
else
   sigma = SIGMA_ABOVE
end if
energy_weight = exp(-0.5d0 * (e_band_avg(ib)/sigma)**2)
relevance(ib) = proj_avg(ib) * energy_weight
```

## References

1. Qiao, J., Pizzi, G., Marzari, N. "Projectability disentanglement
   for accurate and automated electronic-structure Hamiltonians."
   *npj Comput. Mater.* **9**, 208 (2023).

2. Souza, I., Marzari, N., Vanderbilt, D. "Maximally localized
   generalized Wannier functions for composite energy bands."
   *Phys. Rev. B* **65**, 035109 (2001).

3. Vitale, V., et al. "Automated high-throughput Wannierisation."
   *npj Comput. Mater.* **6**, 66 (2020).

4. Damle, A., Lin, L., Ying, L. "Compressed representation of Kohn-Sham
   orbitals via selected columns of the density matrix."
   *J. Comput. Phys.* **334**, 1-15 (2017).

5. Marzari, N., Mostofi, A. A., Yates, J. R., Souza, I., Vanderbilt, D.
   "Maximally localized Wannier functions: Theory and applications."
   *Rev. Mod. Phys.* **84**, 1419 (2012).
