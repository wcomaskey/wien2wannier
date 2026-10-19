# HYBR worked example: GaAs, start to finish

> **v3.0 release note.** This walkthrough uses the named-option `w2w_hybrid`
> CLI (`w2w_hybrid --help` is the authoritative reference).

This example takes fcc GaAs from a WIEN2k SCF to an 18-Wannier-function
valence+conduction model localized by stock Wannier90 with **no `dis_*`
keywords**: the hybrid pipeline does the disentanglement; Wannier90 only
rotates. GaAs is the canonical validated case of this feature line
(band RMS compared to the LAPW bands at the release protocol: 4.9 meV over
[-16, 0] eV and 18.0 meV over [0, 4] eV).

Commands marked *(standard WIEN2k / wien2wannier)* are the stock workflow.
Their options are described in the WIEN2k user guide and
`doc/wien2wannier_userguide.pdf`.

## 0. Prerequisites

- A converged WIEN2k SCF for GaAs (`case = GaAs`): `init_lapw` + `run_lapw`
  *(standard WIEN2k)*.
- This package built: `w2w`/`w2wc` (`make all` in `SRC_w2w/`) and
  `w2w_hybrid` (`make -f Makefile.hybrid`).
- `wannier90.x` (upstream 3.x).

## 1. Wannier k-mesh and projections *(standard wien2wannier)*

From the converged case directory:

```sh
prepare_w2wdir GaAs-WANN && cd GaAs-WANN
x kgen -fbz            # unshifted full-BZ mesh, e.g. 8 8 8
write_inwf             # choose the initial projections: Ga/As s,p (and d)
write_win              # num_wann = 18 for the spd valence+conduction target
x wannier90 -pp        # -> GaAs-WANN.nnkp
```

For the hybrid pipeline the projections in `.inwf` are the *seed* of the
per-state projectability measurement, not the final answer; a plain
atom-centered s/p/d set is enough.

## 2. Bloch matrix elements *(standard wien2wannier)*

```sh
x lapw1
x w2w              # (x w2wc for complex cases; spin-orbit: see guide 02 §6)
```

This writes `GaAs-WANN.amn`, `GaAs-WANN.mmn`, `GaAs-WANN.eig`, the inputs
of the hybrid stage. (The `.eig` written by w2w is already Fermi-shifted and
in eV, so `--efermi 0` (the default) is correct downstream.)

## 3. The hybrid stage: `w2w_hybrid`

Near-default invocation (SCDM gauge, automatic frozen ceiling, per-state or
band-level regime decided automatically):

```sh
path/to/SRC_w2w/w2w_hybrid --seed GaAs-WANN --nw 18
```

or, closer to the regression harness (its full recipe is in guide 03 §6.1):

```sh
path/to/SRC_w2w/w2w_hybrid --seed GaAs-WANN --nw 18 \
    --rescale lowdin --regime per-state --froz-emax 2.0 \
    --sweeps 4000 --relax -1.9 --gauge auto8 --admit-all --solver gs
```

What to read in the output (see guide 01 §6, "Reading the diagnostics", for the full
explanation):

- the **pool** and **capacity** lines: frozen demand vs `--nw`, with
  headroom ≥ 0;
- the per-sweep **Omega_I** trace (set `W2W_VERBOSE=1`) ending in
  convergence, not the sweep cap;
- the **gauge** report; `sigma_min` well away from 0 (a collapsed anchor
  frame is refused);
- the hand-off files `GaAs-WANN-hybF.amn/.mmn/.eig` plus
  `GaAs-WANN-hybF.anchors` (the SCDM anchor columns and gauge-policy stamp).

## 4. Localize with stock Wannier90

Give Wannier90 the hand-off under its own case name:

```sh
mkdir hybF && cd hybF
for e in amn mmn eig; do cp ../GaAs-WANN-hybF.$e GaAs-hybF.$e; done
cp ../GaAs-WANN.nnkp GaAs-hybF.nnkp    # same mesh; regenerate via -pp if renamed
```

`GaAs-hybF.win`: copy the original `.win`, then set

```
num_bands = 18
num_wann  = 18
guiding_centres = F
! no dis_win_*, dis_froz_*, dis_num_iter -- the manifold is already isolated
```

and run:

```sh
wannier90.x GaAs-hybF
```

Because `num_bands = num_wann`, Wannier90 performs rotation-only maximal
localization. The disentanglement is already finished, so no `dis_*` iterations
run, but the localization itself needs tens of thousands of iterations. From the
hand-off of the default CG solver, Wannier90 with its default settings converges
in 24,500 iterations. From the hand-off of `--solver gs` it converges in 59,500
iterations (§5). Once it has converged, the interpolated bands from `bands_plot`
differ from the LAPW reference by 5.6 meV (RMS over [-16, 0] eV) and 17.7 meV
(RMS over [0, 4] eV) for the CG hand-off, and by 4.9 and 18.0 meV for the `gs`
hand-off.

## 5. Verifying the result

- The `.wout` must reach "Wannierisation convergence criteria satisfied".
- Compare the model against the WIEN2k band structure over the model's
  fidelity window (up to `--froz-emax` above E_F); the hybrid stage's masks
  make no promise above that ceiling. `tools/band_compare.py` does this from
  `hybF/GaAs-hybF_hr.dat` and a WIEN2k band path (`x lapw1 -band` in a copy
  of the case directory); `tools/band_compare_gaas.json` is the spec for this
  example and says how to prepare the inputs. Its `.txt` companion lists the
  index-paired RMS over [-16, 0] and [0, 4] eV; the release-protocol GaAs
  model of the paper gives 4.9 and 18.0 meV.
- Read `GaAs-WANN-hybF.dis` (the `wien2wannier` driver copies its two lines to the top of the hand-off `.win`): whether the
  disentanglement converged, and the subspace continuity. A continuity warning
  means the subspace jumps between two neighboring k-points. With the default
  CG solver this GaAs hand-off has continuity 2.8e-2 (no warning), and
  Wannier90 with its default settings converges in 24,500 iterations at 5.6
  and 17.7 meV (Omega_tot 51.40 A^2). From the same hand-off, a
  Stengel-Spaldin stage (`use_ss_functional = true`, `conv_tol = 1e-6`)
  followed by a fresh Marzari-Vanderbilt run from its gauge (written as the
  `.amn`) with `num_cg_steps = 200` reaches a lower minimum, 50.11 A^2 at 5.3
  and 14.1 meV, in 3,500 iterations: which minimum Wannier90 reaches depends
  on its optimizer, not only on the hand-off. `--solver gs` gives continuity
  5.0e-4 (a warning) and still interpolates at 4.9 and 18.0 meV once
  Wannier90 converges (59,500 iterations); that is the model of the paper's
  LAPW table. Compare models only after Wannier90 has converged (one
  uninterrupted run; `conv_tol` met), and treat the warning as a reason to
  compare models, not as evidence that the bands are wrong.
- Re-running step 3 with the same inputs is deterministic; to reproduce a
  gauge exactly across reruns use `W2W_GAUGE0OUT=1` + `--gauge u0` (see
  `w2w_hybrid --help`).

## Variations

- **One-shot** (`MODE = HYBR` in `.inwf`): steps 2–3 fuse; `x w2w` emits
  the hand-off directly; guide 02 lists the per-calculation-type settings.
- **SOC**: use `x w2wc` and the complex binaries; the hybrid stage is
  spinor-aware (sigma-mirrored channel masks; see `--np-mask` in guide 03).
- **Metals / entangled manifolds**: guide 04 covers regime and relaxation
  tuning (`--regime`, `--relax`, `--pool-factor`).
