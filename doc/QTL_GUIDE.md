# QTL-Guided Wannier90 Pipeline

## Overview

The QTL pipeline uses Wien2k's orbital-resolved band character (`x lapw2 -qtl`) to automatically determine optimal energy windows and projections for Wannier90. It produces better results than the PDWF method because QTL provides properly normalized orbital character including interstitial contributions.

## Prerequisites

A completed Wien2k SCF calculation with:
- `case.struct`, `case.klist`, `case.energy`, `case.vector`, `case.vsp`, `case.fermi`
- `case.qtl` from `x lapw2 -qtl` (run on the same k-mesh as the SCF)

Optional:
- `case.klist_band` and `case.spaghetti_ene` for band structure comparison

## Quick Start

```bash
# On your Wien2k cluster, after SCF convergence:
x lapw2 -qtl

# Transfer files to your working directory, then:
qtl_lapw -f CASE -w90exe /path/to/wannier90.x
```

The script handles everything: parsing QTL, classifying bands, generating `.inwf` and `.win`, running `wannier90 -pp`, `w2w`, and `wannier90`.

## How It Works

### Band Classification

1. **Semicore detection**: Bands where any single (atom, l) has QTL character > 0.90 are excluded (e.g., Ga 3d at 0.98)

2. **Target character**: For each band, sum the QTL character of the target orbital channels (determined from the periodic table: sp for main-group, sd for transition metals, etc.)

3. **Relevance weighting**: `relevance = target_char * exp(-E²/2σ²)` where σ is derived from the eigenvalue distribution. This prioritizes bands near E_F.

4. **Percentile classification**: Sort bands by relevance, apply percentile cutoffs:
   - Above the frozen percentile → FROZEN
   - Above the exclusion percentile → DISENTANGLE
   - Below → EXCLUDED

### Energy Windows

- **Frozen window**: Covers all frozen bands with ±0.5 eV padding. A minimum floor (e.g., [-6, +3] eV around E_F) guarantees Fermi surface protection.
- **Outer window**: Covers all non-excluded bands with ±2.0 eV padding. Always strictly wider than the frozen window.
- **num_wann constraint**: At every k-point, bands in the frozen window ≤ num_wann. The frozen window is shrunk from the conduction side if needed.

## Options

```
qtl_lapw [FLAGS] [OPTIONS]

FLAGS:
  -h          help
  -v          verbose
  -norun      generate input files only

OPTIONS:
  -f CASE       case name (default: directory name)
  -bmin N       minimum band index (default: auto from .in1)
  -bmax N       maximum band index (default: auto)
  -phigh X      frozen percentile 0-1 (default: 0.60)
  -plow X       exclusion percentile 0-1 (default: 0.15)
  -nkband N     k-points per band segment (default: 100)
  -w2wdir DIR   w2w/w2wc directory
  -w90exe EXE   wannier90.x path
```

## Examples

### MgB2 (hexagonal metal)

```bash
# On cluster:
x lapw2 -qtl

# Locally:
qtl_lapw -f WANN -w90exe ./wannier90.x
```

Typical result: 7 frozen bands, Omega = 16-22 Ang² depending on percentile.

### GaAs (FCC semiconductor)

```bash
# On cluster (complex calculation):
x lapw2 -qtl -c

# Locally:
qtl_lapw -f GaAs -w90exe ./wannier90.x
```

Typical result: 5-6 frozen bands (As 4s + valence sp), excludes Ga 3d semicore.

## Comparison with PDWF

| Feature | PDWF | QTL |
|---------|------|-----|
| Orbital character | Muffin-tin only (25-40%) | Full (MT + interstitial) |
| Semicore detection | Gap-based heuristic | Direct (QTL > 0.90) |
| Requires cluster step | No | Yes (`x lapw2 -qtl`) |
| Threshold tuning | Needs adaptive sigma | Simpler, more robust |
| Best for | Quick prototyping | Production calculations |

## Band Comparison

After running, generate a comparison plot:

```bash
python3 compare_bands.py CASE -dir .
```

This overlays DFT (spaghetti) and Wannier-interpolated bands using index-based k-path mapping (convention-independent).

## Troubleshooting

- **`FERMI - number of k-points inconsistent`**: The `.energy` and `.kgen` files use different k-meshes. Re-run `x lapw1` to regenerate.
- **`More states in frozen window than target WFs`**: The frozen window is too wide. Reduce `-phigh` or the script will auto-shrink.
- **Poor band fitting near E_F**: Try `-phigh 0.70` for a tighter frozen window with more gauge freedom.
