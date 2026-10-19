# The LAPW models of the paper: one command per cell

This page gives the command behind each row of the LAPW table of the paper
([Comaskey2026] in `REFERENCES.md`) and the commands that follow it. The
WIEN2k cases, the `w2w` products each command reads and the resulting
hand-offs are in the data record of the paper, not in this archive.

## 1. Inputs of one cell

A converged WIEN2k calculation on the full-zone mesh of the table, and from
`x w2w` (see `01-overview-and-quickstart.md`):

| File | Content |
|---|---|
| `<seed>.amn`, `<seed>.mmn`, `<seed>.eig`, `<seed>.nnkp` | projections, overlaps, energies relative to $E_F$, neighbor list |
| `<seed>.amn.tl` | the same projections with the radial functions continued past the muffin-tin radius; the projectability is computed from this file (`W2W_PAMN`) |
| `projmap.dat` | channel map of the projection columns |

## 2. Hand-off

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2
export W2W_SEEDK=opt16 W2W_SEEDLOSS=smooth
W2W_PAMN=$PWD/<seed>.amn.tl \
w2w_hybrid --seed <seed> --nw <N> --out <out> --np-mask <np> \
    --rescale lowdin --regime per-state --froz-emax auto --efermi 0.0 \
    --admit-all --solver <gs|cg> --sweeps 30000 --rtol 1e-10 \
    --check-every 10 --colored --adaptive --relax -1.9 --gauge auto8
```

`--froz-emax auto` walks the frozen ceiling down from $E_F+8$ eV in 1 eV
steps until no k-point has more than `<N>` pool states below it. `--admit-all`
is the default of this release and is included for clarity. The gauge
is transported from the best of up to 16 candidate seed k-points
(`W2W_SEEDK=opt16`).

| Cell | Mesh | `--seed` | `--nw` | `--np-mask` | `--solver` | Ceiling found (eV) |
|---|---|---|---|---|---|---|
| GaAs | $8^3$ | `GaAs-WANN` | 18 | 18 | `gs` | $+8$ |
| MgB2 | $12^2\times8$ | `m12s` | 12 | 12 | `gs` | $+8$ |
| α-Sn (SOC) | $8^3$ | `Wsig` | 20 | 10 | `cg` | $+6$ |
| Bi bilayer (SOC) | $12^2\times1$ | `bi36` | 36 | 18 | `gs` | $+4$ |
| bcc Fe, majority spin | $8^3$ | `Fe_up` | 9 | 9 | `gs` | $+8$ |
| bcc Fe, minority spin | $8^3$ | `Fe_dn` | 9 | 9 | `gs` | $+8$ |
| bcc Fe (SOC) | $8^3$ | `Fe_soc` | 18 | 9 | `gs` | $+8$ |
| SrVO3 | $8^3$ | `SrVO3` | 19 | 19 | `gs` | $+6$ |

On spin-orbit cells `--np-mask` counts the channels of one spin block; the
second block is added automatically. The Bi bilayer's sweep stopped at the
30000-sweep cap with a last change of $1\times10^{-7}$ Å$^2$; every other
minimization met its tolerance. For α-Sn the Gauss–Seidel sweep and
conjugate gradients agree in $\Omega_I$ to $10^{-3}$ Å$^2$, and the table
reports the conjugate-gradient model.

The public driver runs the same rule from a WIEN2k case directory:

```sh
wien2wannier --input /path/to/case --seed <seed> --k-grid <n1 n2 n3> \
    --num-wann <N> --froz-emax auto --solver <gs|cg>
```

It chooses the projections and the band window itself and uses the
integrated HYBR path with that path's defaults for the relaxation and the
gauge, so its hand-off is not byte-identical to the one from the standalone
command above, which is the command that made the table. On MgB2 it gives
a different model from the table's; see Section 7.

## 3. Wannier90

Each hand-off (`<out>.amn`, `.mmn`, `.eig`, with `num_bands = num_wann` and
no `dis_*` keywords) is localized by Wannier90 from its SCDM `.amn`, with
the projection block of the `.win` set to `random`, `conv_tol = 1e-10` and
`conv_window = 5`, and no guiding centres. The Bi bilayer is the exception:
its random-start wannierization did not descend, and it uses guiding centres
on Bi *s*, *p*, *d*. Do not continue a hand-off with `restart = wannierise`;
run it to its criterion in one run (GaAs 59500, Fe with spin-orbit coupling
64500, Bi bilayer 58500 iterations).

## 4. Bands and band distance

`tools/band_compare.py` scores Wannier90's interpolation (with the
Wigner–Seitz correction) against the WIEN2k band path: the band distance
with bands paired by index over each window of the table (`eta_windows`).
`tools/band_compare_gaas.json` is the specification for GaAs; the other
cells use the same form with their own files and windows.

| Cell | Windows (eV from $E_F$) | Band distance (meV) |
|---|---|---|
| GaAs | $[-16, 0]$, $[0, +4]$ | 4.9, 18.0 |
| MgB2 | $[-13, +7.2]$, $[-2, +2]$ | 11.5, 4.4 |
| α-Sn | $[-12, 0]$, $[0, +4]$ | 2.2, 14.2 |
| Bi bilayer | $[-15, 0]$, $[0, +4]$ | 0.2, 17.6 |
| bcc Fe, majority | $[-10, 0]$, $[0, +4]$ | 8.0, 15.1 |
| bcc Fe, minority | $[-10, 0]$, $[0, +4]$ | 11.8, 7.0 |
| bcc Fe (SOC) | $[-10, 0]$, $[0, +4]$ | 4.5, 5.5 |
| SrVO3 | $[-8, 0]$, $[0, +4]$ | 0.7, 1.6 |

## 5. Diagnostic figures

Every hand-off writes `<out>.masks` and `<out>.anatomy`:

```sh
python tools/plot_certificate.py <out>.masks
python tools/plot_anatomy.py <out>.anatomy
```

## 6. Software versions of the paper's runs

wien2wannier 3.0 with WIEN2k 24.1; the Wannier90 development build
`ff955d1`, which v4.0.2 reproduces bitwise.

## 7. What the sweep options of Section 2 change

The seven Gauss–Seidel hand-offs of the table were made with three options
that are not defaults of `w2w_hybrid`: the colored update order
(`--colored`), scheduled over-relaxation that starts at a factor 1.9
(`--relax -1.9`), and adaptive spacing of the convergence checks
(`--adaptive`, starting from `--check-every 10`). Conjugate gradients, the
solver of the α-Sn row and the default of this release, ignores all three.

The table compares each published hand-off with the same command without
`--colored`, that is, with the serial update order. Both were made with the
serial program of this release. The command of Section 2 reproduces the
published `.amn`, `.mmn` and `.eig` byte for byte on GaAs, MgB2, the three
bcc Fe cells and SrVO3 and, with `--solver cg`, on α-Sn; the Bi bilayer
was not rerun with the published options, and its published values are
those of the original run.

| Cell | $\Omega_I$, published (Å$^2$) | Sweeps | $\Omega_I$, serial order (Å$^2$) | Sweeps | Largest change of a subspace eigenvalue | Lowest state that changes |
|---|---|---|---|---|---|---|
| GaAs | 38.97932196 | 120 | 38.97932196 | 120 | 0.004 meV | none |
| bcc Fe, majority | 5.87627656 | 110 | 5.87627656 | 110 | 0.002 meV | none |
| bcc Fe, minority | 6.00223215 | 110 | 6.00223215 | 110 | 0.003 meV | none |
| SrVO3 | 21.96424225 | 592 | 21.96424224 | 1728 | 0.9 meV | none above 1 meV |
| bcc Fe (SOC) | 11.22394220 | 8049 | 11.22393391 | 5239 | 1.1 eV | $+8.2$ eV |
| MgB2 | 12.60297735 | 120 | 12.36906647 | 210 | 3.6 eV | $+8.0$ eV |
| Bi bilayer (SOC) | 80.19172265, at the cap | 30000 | 80.19139918, at the cap | 30000 | 0.19 eV | $+4.3$ eV |
| α-Sn (the sweep; the table reports conjugate gradients) | 55.94369571 | 164 | 55.94286476, at the cap | 30000 | 0.31 eV | $+6.0$ eV |

A subspace eigenvalue is an eigenvalue of the Hamiltonian in the selected
subspace at a mesh k-point, as written to `<out>.eig`. What the comparison
shows:

- Below the frozen ceiling the two hand-offs have the same eigenvalues on
  every cell, as the method guarantees: every state below the ceiling is in
  the subspace whatever the solver does. The update order changes only
  which combination of the states above the ceiling completes the subspace.
- On GaAs, both spin channels of bcc Fe and SrVO3 the two orders give the
  same subspace. On bcc Fe with spin-orbit coupling they give two subspaces
  whose $\Omega_I$ differs by $8\times10^{-6}$ Å$^2$. On MgB2 they stop at
  two different minima; conjugate gradients from the same inputs reaches
  12.60282541 Å$^2$ in 360 gradient evaluations, the minimum of the
  colored order (subspace eigenvalues within 0.12 eV of it, compared to
  3.6 eV for the serial order).
- The two MgB2 minima give different models. Wannierized with the
  protocol of Section 3, the hand-off of the colored order (the model of
  the table) has $\Omega_{\rm tot}=18.153$ Å$^2$ and band distances of 11.5
  and 4.4 meV; the hand-off of the serial order has
  $\Omega_{\rm tot}=14.837$ Å$^2$ and 7.0 and 1.5 meV over the same windows
  (converged after about 42000 iterations; with `num_cg_steps = 200`
  Wannier90 reaches both minima within 1000 iterations). The driver with
  its defaults (conjugate gradients on the projections it chooses) reaches
  the lower minimum on this cell: $\Omega_I=12.369054$ Å$^2$,
  $\Omega_{\rm tot}=14.837$ Å$^2$, 7.1 and 1.5 meV.
- On the Bi bilayer neither order meets the tolerance within 30000
  sweeps; the two hand-offs differ only above the ceiling ($+4$ eV). The
  colored order without the relaxation is further from convergence at
  the cap (80.19400266 Å$^2$).
- On α-Sn the serial order does not meet the tolerance within 30000
  sweeps (55.94286476 Å$^2$ with the options of Section 2, 55.94274898
  Å$^2$ with neither relaxation nor adaptive checks); it is still descending
  towards the value of conjugate gradients (55.94273218 Å$^2$ in 2883
  gradient evaluations), which is the model of the table.
- Over-relaxation and the adaptive checks change the number of sweeps and
  not the result. On MgB2 with the colored order: 120 sweeps with both,
  2210 without the relaxation, 1931 without the adaptive checks, 2211 with
  neither, and the subspace eigenvalues agree to 1.1 meV. With the serial
  order: 210, 2371 (relaxation only) and 2491 (neither), eigenvalues within
  0.6 meV. The saving needs both options together.
