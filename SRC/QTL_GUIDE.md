# QTL-Guided Wannierization for Wien2Wannier

## Overview

The QTL method uses Wien2k's orbital-resolved band character (`x lapw2 -qtl`) to automatically determine optimal frozen/disentangle energy windows for Wannier90. It replaces the manual trial-and-error process of choosing `dis_froz_min`, `dis_froz_max`, `exclude_bands`, and projection orbitals.

### What it does

1. Reads the QTL file to identify which bands have strong target orbital character
2. Automatically excludes core and semicore bands
3. Computes an adaptive frozen window that covers the Fermi-level manifold
4. Generates a complete Wannier90 `.win` file with optimized parameters
5. Generates correct projection orbitals for all atoms (including equivalent atoms)

### When to use it

- After a converged Wien2k SCF calculation
- When you want automated, reproducible Wannier90 input generation
- For any material type (metals, semiconductors, insulators)

---

## Prerequisites

### Wien2k files needed

| File | Source | Purpose |
|------|--------|---------|
| `case.struct` | Wien2k | Crystal structure |
| `case.energy` | `x lapw1` | Eigenvalues for band count |
| `case.vector` | `x lapw1` | Eigenvectors for w2w |
| `case.vsp` | SCF | Potential for w2w |
| `case.fermi` | SCF | Fermi energy |
| `case.klist` | `x kgen` | K-mesh for Wannier90 |
| **`case.qtl`** | **`x lapw2 -qtl`** | **Orbital character (key input)** |
| `case.eig` | `w2w` | Eigenvalues for optimizer |
| `case.amn` | `w2w` | Projection matrix |
| `case.mmn` | `w2w` | Overlap matrix |
| `case.nnkp` | `wannier90 -pp` | Neighbor list |

### Software

- Wien2k (for SCF, lapw1, lapw2 -qtl)
- wien2wannier w2w/w2wc (for AMN/MMN/EIG)
- Wannier90 (for Wannierization)
- Python 3 with NumPy (for qtl_opt.py)

---

## Step-by-Step Procedure

### Step 1: Run Wien2k SCF

Complete the SCF cycle as normal. You need a converged calculation with the desired k-mesh.

### Step 2: Generate QTL file

```bash
# Edit case.in2 — change first line from TOT to QTL
sed -i 's/^TOT/QTL/' case.in2

# Run lapw2 in QTL mode
x lapw2 -qtl

# Restore in2
sed -i 's/^QTL/TOT/' case.in2
```

This produces `case.qtl` with orbital-resolved band character for each k-point and band.

### Step 3: Run w2w to generate AMN/MMN/EIG

First, generate a minimal `.inwf` and `.win` for the w2w run. The band window should include all bands you want to consider (semicore + valence + some conduction). The optimizer will handle exclusion automatically.

```bash
# Create PDWF .inwf (includes all bands, lets optimizer choose)
cat > case.inwf << EOF
PDWF
 bmin  bmax
 5
 0.60  0.15  0.50
EOF

# Generate .win for wannier90 -pp
# (use pdwf_lapw or create manually)

# Run wannier90 -pp to get .nnkp
wannier90.x -pp case

# Run w2w (use w2wc for complex/SOC calculations)
w2w w2w.def    # or w2wc for complex
```

**Choosing bmin**: Include a few bands below the target valence manifold. The optimizer will exclude semicore bands automatically. A good rule: start from the first band after the deepest core gap (>20 eV).

**Choosing bmax**: Include enough conduction bands for disentanglement — typically 1.5-2x `num_wann`.

### Step 4: Run the QTL optimizer

```bash
python3 qtl_opt.py case -dir . -qtl case.qtl -step 0.5 -top 3 -write
```

This will:
1. Auto-detect `bmin` by matching eigenvalue patterns between QTL and EIG
2. Identify and exclude core/semicore bands
3. Compute target orbital character from the valence configuration
4. Sweep frozen window parameters with the loss function
5. Generate `runs_opt/opt_01/`, `runs_opt/opt_02/`, etc. with complete Wannier90 inputs

### Step 5: Run Wannier90

```bash
cd runs_opt/opt_01
wannier90.x case
# or with MPI:
mpirun -n N wannier90.x case
```

### Step 6: Compare band structures

Generate a DFT band structure for comparison:

```bash
# Create a dense klist_band (100 points per segment)
# Run on cluster:
x lapw1 -band
x spaghetti
```

Then compare:

```bash
python3 compare_bands.py case -dir runs_opt/opt_01
```

---

## Command Reference

### qtl_opt.py

```
Usage: qtl_opt.py CASE [options]

Options:
  -dir DIR      Working directory (default: .)
  -qtl FILE     QTL file path (default: DIR/CASE.qtl)
  -step FLOAT   Grid step size in eV (default: 0.5)
  -top N        Number of top configs to show/write (default: 15)
  -write        Generate Wannier90 input folders for top configs
  -efmin FLOAT  Minimum frozen window below E_F in eV (default: 4.0)
  -efmax FLOAT  Minimum frozen window above E_F in eV (default: 2.0)
  -nwann N      Override num_wann (default: auto from valence config)
  -bmin N       Override bmin (default: auto from eigenvalue matching)
```

### compare_bands.py

```
Usage: compare_bands.py CASE [options]

Options:
  -dir DIR      Working directory (default: .)
  -out FILE     Output PNG (default: CASE_band_comparison.png)
  -emin FLOAT   Energy min relative to E_F (default: -15)
  -emax FLOAT   Energy max relative to E_F (default: +10)
  -ezoom FLOAT  Zoom half-width around E_F (default: 4)
```

---

## How the Optimizer Works

### Band exclusion (automatic)

1. **Energy gaps**: Bands below the largest gap >5 eV below E_F are excluded (core/semicore)
2. **d/f character**: Bands with >90% d or f orbital character on any atom are excluded
3. **Zero target character**: Bands with <1% of the maximum target character are excluded

### Frozen window selection

The loss function scores each candidate (fmin, fmax) configuration:

```
Score = 1.0 × avg_proj       (energy-weighted target character of frozen bands)
      + 0.3 × coverage       (valence-side window width)
      + 0.4 × gauge          (disentanglement freedom: log(1+n_disent))
      + 0.6 × fermi_coverage (fraction of bands near E_F that are frozen)
      + 0.1 × parsimony      (mild preference for fewer frozen bands)
      - penalties             (conduction overextend, asymmetry, weak bands)
```

### Adaptive frozen floor

The minimum frozen window is computed from the band structure:
- **fmin**: Covers the deepest target band below E_F (+ 0.5 eV padding)
- **fmax for metals**: Covers the top of all Fermi-crossing bands (+ 0.5 eV)
- **fmax for semiconductors**: Covers the first conduction band above the gap (+ 0.5 eV)

### Constraints

- At every k-point: bands in frozen window ≤ num_wann
- At every k-point: bands in outer window > num_wann (ensures disentanglement)
- Frozen bands ≤ 65% of num_wann (ensures gauge freedom)
- Minimum 2 frozen bands

---

## Typical Results

### MgB2 (metal, hexagonal)

```
System: 18 bands, 512 k-pts, num_wann=12
Target: Mg(sp) + 2×B(sp) = 4 + 4 + 4 = 12
Excluded: bands 1-2 (Mg 2p semicore, gap=32.4 eV)
Excluded: bands 16-18 (zero target character)
Adaptive floor: [-12.8, +7.2] eV (3 Fermi-crossing bands)

Top config: 6 frozen [3-8], window [-13.3, 7.5]
  Omega_Total = 18.24 Ang²
  Band MAE = 0.004 eV within frozen window
```

### GaAs (semiconductor, FCC)

```
System: 18 bands, 512 k-pts, num_wann=8
Target: Ga(sp) + As(sp) = 4 + 4 = 8
Excluded: bands 1-6 (As 3d + Ga 3d semicore)
Adaptive floor: [-13.3, +5.0] eV (first conduction band to +4.5)

Top config: 5 frozen [7-11], window [-13.8, 5.5]
  Omega_Total = 21.99 Ang²
  Band MAE = 0.008 eV within frozen window
```

---

## Troubleshooting

### "param_get_projections: too few projection functions defined"

The projections block in `.win` doesn't match `num_wann`. Check that:
- The struct file is accessible when running `qtl_opt.py -write`
- Equivalent atoms (mult > 1) have separate projection entries
- The total number of projection functions = num_wann

### "dis_windows: More states in the frozen window than target WFs"

The frozen window contains too many active bands at some k-point. The optimizer checks this constraint, but if you manually edit the `.win`, ensure:
- `exclude_bands` removes all intended bands
- The frozen window doesn't accidentally capture extra dispersive bands

### "num_bands must be greater than or equal to num_wann"

After applying `exclude_bands`, the remaining band count is less than `num_wann`. Include more bands in the w2w run (increase `bmax` in `.inwf`).

### bmin auto-detection warning

If the QTL and EIG files use very different Fermi energies (e.g., QTL from a different SCF), the eigenvalue matching may fail. Solutions:
- Re-run `x lapw2 -qtl` with the same SCF as the EIG
- Specify `-bmin N` manually

### Poor band accuracy despite good Omega

- Check if the frozen window is too wide on the conduction side — tighten fmax
- Check if projections have wrong atom positions (mult > 1 issue)
- Increase convergence: `conv_tol = 1.0e-12`, `num_iter = 10000`
