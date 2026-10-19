# HYBR by calculation type: regular, spin-polarized, spin-orbit

> **Public command.** `wien2wannier` selects and validates the native binary and
> SOC def units for normal workflows; see `00-public-wien2wannier-cli.md` and
> `wien2wannier --help`. Direct `w2w`/`w2wc` and `w2w_hybrid` recipes below are
> retained for expert control. Source line citations (`file.F:NN`) and paths under
> `tests/` in this guide refer to the development tree. The line numbers may
> be offset in the released sources, and the `tests/` cases are in the data
> record of the paper, not in this archive; every input deck quoted from them
> is reproduced in the text. The symmetry-adapted features of the development
> tree are not part of this release. `x_lapw:NN` cites the `x` script of WIEN2k
> 24.1 as distributed (`$WIENROOT/SRC/x_lapw`); this package does not contain
> that script.


This document gives three complete recipes for `MODE = HYBR`: a scalar-relativistic
non-magnetic case, a collinear spin-polarized case, and a spin-orbit case. Each
one is a full input set; `.inwf`, def file, command lines, and the numbers to
check afterwards, taken from a run in the development tree.

Read `01-overview-and-quickstart.md` first if you have not: the single most
important structural fact about HYBR is that the hand-off it writes is an
**isolated `num_bands = num_wann` problem**, so Wannier90 does no
disentanglement and no `dis_*` keyword applies. Everything below assumes that.

Every behavioural claim is anchored with `file:line`. Paths are relative to the
repository root.

---

## 1. Choosing the binary: `w2w` vs `w2wc`

`w2w` is built twice from the same sources, with and without `-D_COMPLEX_`
(`SRC_w2w/Makefile:171`). The only thing the flag changes is the type of the
LAPW eigenvector coefficient array `a`; `real(R8)` in the real build,
`complex(C16)` in the complex build (`SRC_w2w/modules_rc.F:22-33`), and,
following from that, whether `almgen` calls `dgemm` or `zgemm`
(`SRC_w2w/modules_rc.F:165-181`, `818-840`). The native frontend probes the
vector record before physics runs and rejects a known real/complex mismatch
(`SRC_w2w/main.F:232-254`). If the file kind cannot be identified it reports
that the probe was inconclusive before continuing. The public driver chooses
`w2w` or `w2wc` from the calculation contract, and SOC always selects `w2wc`.

| Case | Binary | Why |
|---|---|---|
| Centrosymmetric cell, no SOC | `w2w` | `lapw1` writes real eigenvectors |
| Non-centrosymmetric cell, no SOC | `w2wc` | `lapw1c` writes complex eigenvectors |
| Any SOC calculation | `w2wc` | `lapwso` always writes complex spinor vectors |

`x_lapw` picks the suffix from the presence of the inversion operation in
`case.struct` (`x_lapw:496-500`) and **does not override it for `-so`**.
On a centrosymmetric cell, `x w2w -so -up` can therefore select the *real*
binary, which the native vector probe rejects. Use the public `wien2wannier`
command, or `x w2wc -so -up` (`x_lapw:544-546` forces
`cmplx = c`), or call the binary directly on your own def file, which is what
all the recipes below do.

The standalone `w2w_hybrid_acc` has no real/complex variants: it starts from
`amn`/`mmn` files, which are complex in every case.

---

## 2. The three types at a glance

| | Regular | Spin-polarized collinear | Spin-orbit |
|---|---|---|---|
| Binary | `w2w` or `w2wc` (by inversion symmetry) | `w2w` or `w2wc` | `w2wc` always |
| `w2w` runs | 1 | 2 (`-up`, `-dn`), independent | 1 (sigma loop) |
| Extra def units | — | — | `9` (2nd vector), `19` (2nd vsp, SP+SOC only) |
| `num_wann` parity | free | free, per spin | should be even (Kramers); not enforced by code |
| `.nnkp` files | 1 | 1 *content*, 2 *filenames* | 1 |
| `amn` columns | `Nproj` | `Nproj` per spin | `2 × Nproj` (one block per sigma) |
| Combine step | none | none (two separate models) | none in HYBR; `w2waddsp` only on the legacy route |
| `.win` `spinors` | absent | absent | `spinors = .true.` |
| Hand-off | `<nnkp base>-hyb.*` | one per spin | one |

The `.nnkp` row matters because the hand-off filename is
derived from the *nnkp filename* by stripping `.nnkp` and appending `-hyb`
(`SRC_w2w/main.F:1012-1015`). Two spin runs pointing unit 11 at the same
`case.nnkp` write to the same `case-hyb.{amn,mmn,eig}`, and the second silently
overwrites the first.

---

## 3. What every type needs first

HYBR consumes the output of an ordinary, **converged** WIEN2k run. It does not
run SCF, and it builds a Wannier model from an unconverged
density.

1. **A converged SCF** with the relativistic treatment you want
   (`run_lapw`, `runsp_lapw`, or `runsp_lapw -so`).
2. **A full-BZ k-mesh**, not the IBZ mesh: `x kgen -fbz`. Wannier90 needs every
   k of the Monkhorst–Pack grid.
3. **One `lapw1` (+ `lapwso`) pass on that mesh**, producing `vector` and
   `energy` files for every k.
4. **A working directory with the Fermi energy in it.**
   `SRC/prepare_w2wdir_lapw` copies `struct`, `in1`, `vsp`, `vns` (plus `inso`
   for SOC) into a subdirectory and writes `case.fermi`, or `case.fermiup` and
   `case.fermidn` when spin-polarized, by grepping `:FER` out of `case.scf`,
   falling back to `case.scf2{up,dn}` (`SRC/prepare_w2wdir_lapw:81-127`,
   `291-296`). The file holds one number, in Rydberg.
5. **A `.nnkp`**: write a `.win` with the k-mesh and lattice, then run
   `wannier90.x -pp <seed>`. `init_w2w_lapw` does steps 2–5 interactively
   (`SRC/init_w2w_lapw:217-316`).

Check two things before starting a long run:

* `case.fermiup` and `case.fermidn` must contain the **same** number. A
  spin-polarized calculation has one common Fermi level; the two files exist
  only because the def-file machinery is per-spin.
  `tests/w2w_runs/Fe/Fe.fermi{up,dn}` both hold `0.7151300814`.
* `projmap.dat`, `projectability.dat` and `autoprojection.diag` are written to
  **fixed filenames** in the current directory (`SRC_w2w/main.F:493-509`,
  `SRC_w2w/modules_rc.F:1322-1376`, `1885-1920`). Two HYBR runs in the same
  directory overwrite each other's diagnostics. Give each spin or each variant
  its own directory, or run them sequentially and move the files.

### The `.inwf` skeleton shared by all three

```
HYBR
  0  -1
  3
0 0 T T -1
0 1.0 per-state
froz_emax = 2.0
```

Line by line (parser: `SRC_w2w/util_w2w.F:1544-1723`):

* **`HYBR`**; AUTO channel selection plus in-memory Stage 2
  (`util_w2w.F:1577-1583`).
* **`0 -1`**; `bmin`/`bmax` sentinels. `0` means gap-walk auto-`bmin`
  (`main.F:283-304`); `-1` means "highest band whose bottom lies below
  `E_F + 40 eV`" (`main.F:267-275`). Both fields are mandatory; there is no
  default.
* **`3`**; `LJmax` for the Bessel/Gaunt expansion of the plane-wave overlap
  (`util_w2w.F:1607`). Every validated input in the tree uses 3.
* **`0 0 T T -1`**; `threshold max_proj write_proj p2_weight support_ratio`.
  `threshold = 0` plus `support_ratio = -1` is cliff mode with the algorithmic
  Otsu support cut; both are non-default (0.10 and 0.5) and **both are
  required**, either alone degrades silently. `write_proj = T` (default `F`)
  writes the two diagnostic files. `p2_weight = T` is the default.
* **`0 1.0 per-state`**; `num_wann mix regime`. `num_wann = 0` is the
  capacity-driven auto sentinel and requires a numeric `froz_emax`
  (`main.F:537-538`). `mix = 1.0` switches relaxation off (default `-2.5`). `regime = per-state`
  is non-default (`auto`).
* **`froz_emax = 2.0`**, non-default (`auto`, an 8 eV walk-down). The numeric
  form is what enables the pre-`l2Amn` capacity gate and the `num_wann = 0`
  sentinel (`main.F:531-539`). 2.0 eV is the validated value for the small
  benchmark cells; it is a physics statement about which bands must be
  reproduced exactly, not a convergence parameter.

`pool_factor = 1.5`, `margin = 0`, `gauge = auto8`, `seedk = opt16`,
`seedloss = smooth` and `capacity = strict` are already at their validated
values by default and can be omitted. Do **not** write `rescale = ...`: it is
not a keyword, and unknown
keys are a hard error (`util_w2w.F:1692-1693`).

---

## 4. Regular: scalar-relativistic, non-magnetic

Worked example: **GaAs**, zincblende, 8×8×8 mesh (512 k), `tests/GaAs-WANN`
with the HYBR run staged in `tests/regress_cases/gaas_w2w`.

GaAs has no inversion center, so this is a **`w2wc`** case.

### 4.1 `.inwf`

`tests/regress_cases/gaas_w2w/hybr.inwf`, complete:

```
HYBR
  6  90
  3
0 0 T T -1
0 1.0 per-state
froz_emax = 2.0
```

This file uses an explicit window `[6, 90]` rather than the `0 -1` sentinels. Both work. The explicit form is useful to fix the window once the
gap-walk choice has been inspected. `bmin = 6` trims the **As 3d** semicore complex: bands 1–5, at −2.2753/−2.2745 Ry, i.e. ≈ −35 eV relative to
`E_F = 0.3168583747 Ry` (`tests/GaAs-WANN/GaAs-WANN.fermi`,
`tests/GaAs-WANN/GaAs-WANN.energy:6-15`). The **Ga 3d** complex (bands 6–10,
−0.7760/−0.7702 Ry, ≈ −15 eV) stays *inside* the window and is
the highest-ranked AUTO channel of the run
(`tests/regress_cases/gaas_w2w/w2wH.outputwf:84`, `1 2 0.0075 * [LO-shell]`;
atom 1 = Ga, `tests/GaAs-WANN/GaAs-WANN.struct:7`), selected through its
semicore-LO radial `rf = 3` (`tests/regress_cases/gaas_w2w/projmap.dat:3-7`).

### 4.2 Def file

`tests/regress_cases/gaas_w2w/w2wH.def`, complete (paths shortened):

```
  5,'hybr.inwf',        'old',    'formatted',  0
  6,'w2wH.outputwf',    'unknown','formatted',  0
  7,'w2wH.amn',         'unknown','formatted',  0
  8,'w2wH.mmn',         'unknown','formatted',  0
 10,'GaAs-WANN.vector', 'old',    'unformatted',0
 11,'w2wH.nnkp',        'old',    'formatted',  0
 12,'w2wH.eig',         'unknown','formatted',  0
 18,'GaAs-WANN.vsp',    'old',    'formatted',  0
 20,'GaAs-WANN.struct', 'old',    'formatted',  0
 50,'GaAs-WANN.energy', 'old',    'formatted',  0
 51,'GaAs-WANN.fermi',  'old',    'formatted',  0
```

Unit meanings are fixed by parameter declarations, not by the def reader:
`unit_def=1, unit_in=5, unit_out=6, unit_vector=10` and `unit_vsp=18,
unit_struct=20` in `SRC_w2w/mod_w2w.f:35-36`; `unit_amn=7, unit_mmn=8,
unit_nnkp=11, unit_eig=12` and `unit_ene=50, unit_fermi=51` in
`SRC_w2w/modules.f:46-47`. The def reader at `SRC_w2w/main.F:116-155` dispatches
on those named constants (plus the bare literals `9`, `19` and `50`), never on
the other numbers. `w2wH.nnkp` is a copy of
`GaAs-WANN.nnkp` under a run-specific name, purely so the hand-off comes out as
`w2wH-hyb.*`.

### 4.3 Commands

```sh
cd tests/regress_cases/gaas_w2w
cp ../../GaAs-WANN/GaAs-WANN.nnkp w2wH.nnkp
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
path/to/SRC_w2w/w2wc w2wH.def

python3 ../../make_small_win.py w2wH-hyb 18 "8 8 8" w2wH-hyb.win
wannier90.x w2wH-hyb
```

Set both thread counts. Homebrew OpenBLAS is an OpenMP build, and leaving both
unset nests thread pools.

`make_small_win.py` emits the *structure* the hand-off needs: `num_bands = num_wann`, `guiding_centres = F`, a `random` projections block,
and the lattice and k-points lifted out of the `.nnkp`
(`tests/make_small_win.py:7-22`). It is not byte-identical to the archived
`.win` files: it writes `num_iter = 5000` and **no `write_hr` line at all**
(`make_small_win.py:9`), whereas every checked-in hand-off `.win` begins with
`write_hr=.true.` and carries `num_iter = 10000`
(`tests/regress_cases/gaas_w2w/w2wH-hyb.win:1-4`). Those two edits were applied
by hand afterwards; add them yourself if you want `_hr.dat` out. See §7.1 for
why the stock `write_win` is not a safe substitute.

The equivalent through the WIEN2k wrappers is `x w2wc` with a hand-edited
`case.inwf`; the wrapper's def file (`x_lapw:3088-3100`) is otherwise
identical.

### 4.4 What the run prints, and what it means

From `tests/regress_cases/gaas_w2w/w2wH.outputwf`:

```
 NUM_KPTS=         512
 NNTOT=           8
 MODE= Mmn Amn AUTO
 band window = [6, 90]
 P2 ranking weight: erfc onset E_F+ 5.0 eV, width  2.0 eV
   support spectrum (nsup= 17), otsu ncut= 14
   cliff (support-otsu):  14 mask channels of  16 selected; mask columns =   62
   Channels selected:   16
   Total projections:   84
 HYBR capacity gate: need=11 (froz_emax=  2.000 eV)  margin=4  auto num_wann=18  (6 of 16 ranked channels)
 HYBR: P from 62 physical of 84 spatial columns
 HYBR: Lowdin P (full-cell): max =   1.0000   mean =   0.7294
 HYBR: MV b-weight residual =   2.22E-16
 HYBR: admit-all (validated recipe default)
 HYBR: regime=per-state  pool[1..27] npool=27  frozen/k mean=  9.26  demand(worst k)=11
  disentangle sweeps=  2000  converged=F  |dOmega|_last=  1.05E-04  cap=  2000
  [WARN] disentangle hit the sweep cap WITHOUT converging
 HYBR: Omega_I =      21.83980920
 HYBR: gauge=auto8  seedk=opt16  seedloss=smooth
  SCDM sigma_min(worst k) =   1.4180E-03
  auto gauge (tournament): winner seed k=329  Omega =        90.6728
 HYBR: hand-off max sigma(M) =     0.999044  OK (contractive)
 HYBR: wrote w2wH-hyb.{amn,mmn,eig}
```

### 4.5 How to tell it worked

| Check | GaAs value | Meaning |
|---|---|---|
| `MV b-weight residual` | `2.22E-16` | b-vector weights solve the Marzari–Vanderbilt completeness relation. Anything above ~1e-6 means a broken `.nnkp`/mesh. |
| `demand(worst k)` ≤ `num_wann` | 11 ≤ 18 | The frozen manifold fits. A violation croaks unless `capacity = warn`. |
| `Lowdin P` `max` near 1 | 1.0000 | The projector set spans at least some states completely. A `max` well below 1 means the channel selection missed the physics. |
| `Omega_I` per WF | 21.84/18 = 1.21 Å² | The subspace quality, what Wannier90 cannot improve. Compare only within one material, mesh and `nw`; see the reference table in `04-methods-and-tuning.md` §6.2. |
| disentangler converged | **no** here (`cap= 2000`, `\|dOmega\|_last= 1.05E-04`; a run made at the former fixed cap of 2000) | The cap is the `sweeps` keyword (default 20000). Hitting the cap is a warning, not a failure, but if `\|dOmega\|_last` is large, raise `sweeps`, or rerun Stage 2 from the standalone with a larger **`--sweeps`** (default 20000). The standalone takes named `--key value` options; see `w2w_hybrid --help`. |
| `max sigma(M)` ≤ 1 + 1e-6 | 0.999044 | Hand-off overlaps are between orthonormal subspaces. A `VIOLATION` here means the pipeline is broken, not merely suboptimal. |
| Wannier90 `Omega I` | 21.839810 | **Must equal** the HYBR `Omega_I` to ~7–8 digits. If it does not, Wannier90 is disentangling something, i.e. your `.win` is wrong. |
| Final spread per WF | 23.2997/18 = 1.29 Å² | Total localization. |

The `Omega I` row is the single best end-to-end check. Wannier90 reports
`Omega I = 21.839810407` for this hand-off
(`tests/regress_cases/gaas_w2w/w2wH-hyb.wout:250270`) compared to HYBR's
`21.83980920`, agreement to seven significant figures, because the subspace is
frozen by construction. Fe and Bi2D below agree to seven as well.

Note `SCDM sigma_min = 1.418E-03`, below the hardcoded `2e-3` gate
(`hybrid_core.F:2005`, `real(dp), parameter :: sig_gate = 2.0e-3_dp`, applied at
`:2017` as `if (best_svmin < sig_gate)`), so `auto8` fell through to the seed
tournament instead of keeping the per-k SCDM gauge. Do not confuse `sig_gate`
with `sgate5` (`hybrid_core.F:1886`), a *separate* `2e-3` constant that gates
candidate *stability* in the seed tournament and is settable with
`W2W_SEEDGATE`. Falling through is normal, not an error. It is the only branch
of `auto<N>` in which `seedk`/`seedloss` matter, but the `ptg` and `pts<N>`
gauge modes consume them as well: `ptg` propagates from `k0seed`
(`hybrid_core.F:1992-1994`), and `pts<N>` both runs the tournament with
`lossmode=seed_loss` (`:2079-2086`) and propagates from `k0seed` (`:2093`).

---

## 5. Spin-polarized collinear (magnetic, no SOC)

Worked example: **bcc Fe**, ferromagnetic, `tests/w2w_runs/Fe` with the HYBR
runs staged in `tests/regress_cases/fe_w2w`.

### 5.1 The structure of the calculation

Nothing in Stage 2 is spin-aware. A collinear magnetic calculation is
**two independent HYBR runs** producing **two independent Wannier models**,
which you then use as a pair. There is no combine step, and `w2waddsp` is not
involved, that program exists for SOC only (`doc/wien2wannier_userguide.tex:615-628`).

Per spin you need: `case.vector{up,dn}`, `case.vsp{up,dn}`,
`case.energy{up,dn}`, `case.fermi{up,dn}`. Shared: `case.struct`, and the
b-vector geometry.

**One nnkp, two filenames.** The `.nnkp` content must be identical between the
spins (same k-mesh, same neighbor list, same G-shifts) because it defines
the finite-difference geometry, which is a property of the lattice, not the
spin. But the two def files must name it differently, because the hand-off
filename is derived from it (`SRC_w2w/main.F:1012-1015`). In the reference case
`feupH.nnkp` and `fednH.nnkp` are byte-identical copies of the same file.

### 5.2 `.inwf`, one file, used by both

`tests/regress_cases/fe_w2w/hybr.inwf`, complete:

```
HYBR
  0   0
  3
0 0 T T -1
9 1.0 per-state
froz_emax = 2.0
```

`bmax = 0` resolves to `minNEn`, the number of bands present at *every* k
(`main.F:262-266`); together with the gap-walk `bmin` this gives the window
`[5, 19]` for both spins. `num_wann = 9` is explicit here (Fe s + d = 1 + 5,
plus p, per spin) rather than the auto sentinel.

You may give the two spins **different** `num_wann` if their
capacities differ; use two copies of the file in that case.

### 5.3 Def files, the only things that differ

`tests/regress_cases/fe_w2w/feupH.def`:

```
  5,'hybr.inwf',      'old',    'formatted',  0
  6,'feupH.outputwf', 'unknown','formatted',  0
  7,'feupH.amn',      'unknown','formatted',  0
  8,'feupH.mmn',      'unknown','formatted',  0
 10,'Fe.vectorup',    'old',    'unformatted',0
 11,'feupH.nnkp',     'old',    'formatted',  0
 12,'feupH.eig',      'unknown','formatted',  0
 18,'Fe.vspup',       'old',    'formatted',  0
 20,'Fe.struct',      'old',    'formatted',  0
 50,'Fe.energyup',    'old',    'formatted',  0
 51,'Fe.fermiup',     'old',    'formatted',  0
```

`tests/regress_cases/fe_w2w/fednH.def` is the same file with `up → dn` on units
7, 8, 10, 11, 12, 18, 50 and 51, including the output names and the nnkp name.
**Unit 5 is not renamed**: both def files point at the shared `hybr.inwf`
(`feupH.def:1` and `fednH.def:1` are byte-identical), which is the §5.2 "one
file, used by both".

### 5.4 Commands

```sh
cd tests/regress_cases/fe_w2w
cp ../../w2w_runs/Fe/Fe.nnkp feupH.nnkp
cp ../../w2w_runs/Fe/Fe.nnkp fednH.nnkp
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1

path/to/SRC_w2w/w2w feupH.def     # majority
mv projmap.dat projmap.up.dat
mv autoprojection.diag autoprojection.up.diag
path/to/SRC_w2w/w2w fednH.def     # minority

python3 ../../make_small_win.py feupH-hyb 9 "8 8 8" feupH-hyb.win
python3 ../../make_small_win.py fednH-hyb 9 "8 8 8" fednH-hyb.win
wannier90.x feupH-hyb
wannier90.x fednH-hyb
```

bcc Fe is centrosymmetric, so the **real** build `w2w` is correct here.

The two `mv` lines are the fixed-filename workaround from §3. Skip them and you
lose the majority-spin diagnostics.

Do **not** use `x wannier90 -up` / `-dn` on the hand-off. That wrapper
symlinks `case.amn{up,dn}` to `case.amn` and renames outputs
(`SRC/wannier90_lapw:156-163, 180-186`); it is built for the legacy `case.*`
naming, and the hand-off is `<nnkp base>-hyb.*`. Call `wannier90.x` directly.

### 5.5 The two spin channels differ, and the difference is physical

The reference Fe run gives, for identical inputs apart from the spin files:

| | majority (`up`) | minority (`dn`) |
|---|---|---|
| band window | `[5, 19]` | `[5, 19]` |
| channels selected / columns | 3 / 9 | 3 / 9 |
| `num_wann` | 9 | 9 |
| `need` at `froz_emax = 2.0` | 6 | 6 |
| pool | `[1..14]`, `npool = 14` | `[1..14]`, `npool = 14` |
| **frozen/k mean** | **5.54** | **4.75** |
| **Omega_I** | **6.97607571** | **7.17933349** |
| `SCDM sigma_min` | 4.2601E-01 | 4.1412E-01 |
| `max sigma(M)` | 0.996393 | 0.996527 |
| Wannier90 total spread | 8.0834 Å² (0.898/WF) | 8.3586 Å² (0.929/WF) |

The minority channel freezes ~0.8 fewer states per k and ends up with a 3 %
larger `Omega_I`. **That is the exchange splitting, not a bug.** In FM Fe the
majority d bands sit below `E_F` as a nearly isolated complex, while the
minority d bands straddle `E_F` and hybridize with the sp continuum, so a
larger fraction of the minority window is entangled. Expect:

* different `frozen/k mean` between spins (always);
* different `Omega_I`, typically a few percent, occasionally more for a strong
  magnet;
* possibly a different `demand(worst k)`, and hence a different capacity
  outcome at the same `froz_emax`.

What would be a real problem, as opposed to physics:

* **different band windows** between spins when both used sentinels. The window
  is resolved per run from that spin's `energy` file, so a spin-dependent
  `minNEn` or a spin-dependent semicore gap can move it. If the two windows
  differ, set both explicitly to the same `bmin`/`bmax` before comparing
  anything, or the two models are not describing the same bands.
* **different channel selections** between spins. AUTO ranks channels from that
  spin's projectability; a marginal channel can flip. If it does, fix the
  selection (raise `threshold` above 0, or set `max_proj`) so both spins use
  the same orbital basis, otherwise the two `hr.dat` files are not
  index-compatible and you cannot build a spin-resolved Hamiltonian from them.
* **a capacity violation on one spin only**. Fix it by raising that spin's
  `num_wann` (they may differ), not by lowering `froz_emax` for both.

### 5.6 How to tell it worked

Run the §4.5 checklist on **each** spin separately, then add:

1. `feupH-hyb.eig` and `fednH-hyb.eig` have the same shape (`nw × nk`) and both
   are Fermi-shifted to eV; HYBR passes `e_fermi = 0` to the writer
   (`main.F:1027`) because `eig_h` was already shifted at `main.F:897-901`.
   The two spins are therefore on a **common** energy zero, provided
   `Fe.fermiup` and `Fe.fermidn` held the same number.
2. Wannier90's `Omega I` matches each run's HYBR `Omega_I`: 6.976075920 vs
   6.97607571 (up) and 7.179333707 vs 7.17933349 (dn).
3. The magnetic moment implied by the two models is sane. With
   `num_wann = 9` per spin and `hr_plot = .true.`, integrating the two
   `_hr.dat` densities of states to `E = 0` should reproduce the WIEN2k
   `:MMTOT` to within the model's fidelity.

---

## 6. Spin-orbit coupling

Worked example: **Bi2D** (buckled bismuthene monolayer, non-spin-polarized
SOC), `tests/Bi2D` with HYBR runs in `tests/regress_cases/bi2d_w2w`. For SP+SOC
see `tests/alpha-Sn` and `tests/Fe_soc`.

### 6.1 What SOC changes

`lapwso` produces **spinor** eigenvectors, stored as two files: one per spin
component, at the same k and with the same eigenvalues. With `x lapwso`
(`x_lapw:1016-1054`) these come out as

* non-spin-polarized SOC: `case.vectorso` + `case.vectorsodn`, energies in
  `case.energyso`;
* SP+SOC (`x lapwso -up`): `case.vectorsoup` + `case.vectorsodn`, energies in
  `case.energysoup`.

SOC is therefore **always a complex calculation**; use `w2wc`.

Presenting the second component to `w2w` on **def-file unit 9** turns on the
sigma loop (`SRC_w2w/main.F:122-125`). This is the only switch; nothing in the
`.inwf` mentions SOC. On SP+SOC, add **unit 19** (the second spin's `vsp`) so
the sigma-2 pass uses the right radial potential (`main.F:126-129`,
`713-719`); without it the sigma-2 pass reuses sigma-1's `vsp`, which is only
correct when the two spins share a potential.

The stock `x w2w` def file (`x_lapw:3088-3100`) has **no unit 9 or 19**.
The SOC one-shot needs a hand-written def file.

### 6.2 The sigma loop, and the spinor `amn` layout

With `sigma2 = .true.` the main loop runs twice (`main.F:689-761`):

* **`mmn`**: the overlap is accumulated *in place* across the two passes, so
  the written `mmn` is the spinor overlap
  `M = ⟨u↑|u↑⟩ + ⟨u↓|u↓⟩`, a **sum**.
* **`amn`**: each pass produces its own `Nproj` columns, and they are stored
  into **disjoint column blocks** of a `Nb × 2·Nproj × Nk` array
  (`main.F:744-753`). The header the writer emits is literally
  `Nb  Nk  2*Nproj` (`main.F:654`), and the file is written block by block
  (`main.F:829-837`). This is a **concatenation**, not a sum.
* **`eig`**: written once, from the sigma-1 pass only, the sigma-2 read passes
  `write_eig = .false.` (`main.F:700-707`). Spinor eigenvalues are identical
  between the two components, so this is correct.
* **AUTO projectability**: a spinor *sum* over both components, via the
  `accumulate`/`finalize` arguments to `compute_projectability`
  (`main.F:393-472`).
* **auto `num_wann`**: each candidate channel contributes
  `2·(2l+1)·mult` columns, because `sf = 2` under the sigma loop
  (`main.F:533, 562`), so the auto sentinel always produces an **even** `nw`.
* **mask columns**: the leading `mask_cols` of sigma block 1 are mirrored into
  sigma block 2, so the projectability span is not sigma-asymmetric
  (`main.F:909-923`).

The concatenated layout is needed because the spin-*summed* `amn` that a
naive `w2waddsp` produces is structurally rank-deficient for spinor Wannier
functions: it offers only `Nproj` candidate columns for a manifold whose
natural basis is `orbital × {↑,↓}`. On Rh3Si4 the summed route gave a worst-k
SCDM `sigma_min` of 6e-8 versus 3.7e-4 for the 160-candidate block layout
(`tests/Rh3Si4/STATE.md:22-33`).

### 6.3 `.inwf`, unchanged

`tests/regress_cases/bi2d_w2w/hybr.inwf`, complete:

```
HYBR
  9  120
  3
0 0 T T -1
0 1.0 per-state
froz_emax = 2.0
```

Identical in form to the non-magnetic file. The window indices count **spinor**
states, because the `energy` file from `lapwso` lists spinor states.

### 6.4 Def file, units 9 (and 19)

`tests/regress_cases/bi2d_w2w/w2wBH.def`, complete (non-spin-polarized SOC,
paths shortened):

```
  5,'hybr.inwf',         'old',    'formatted',  0
  6,'w2wBH.outputwf',    'unknown','formatted',  0
  7,'w2wBH.amn',         'unknown','formatted',  0
  8,'w2wBH.mmn',         'unknown','formatted',  0
  9,'Bi2D.vectorsodn',   'old',    'unformatted',0
 10,'Bi2D.vectorso',     'old',    'unformatted',0
 11,'w2wBH.nnkp',        'old',    'formatted',  0
 12,'w2wBH.eig',         'unknown','formatted',  0
 18,'Bi2D.vsp',          'old',    'formatted',  0
 20,'Bi2D.struct',       'old',    'formatted',  0
 50,'Bi2D.energyso',     'old',    'formatted',  0
 51,'Bi2D.fermi',        'old',    'formatted',  0
```

For **SP+SOC**, add unit 19 and use the per-spin potentials and energies.
`tests/alpha-Sn/sigw2w_v2.def`, complete:

```
  5,'W.inwfsig',       'old',    'formatted',  0
  6,'W.outputwf_sig2', 'unknown','formatted',  0
  7,'W.amn_sig2',      'unknown','formatted',  0
  8,'W.mmn_sig2',      'unknown','formatted',  0
 10,'./W.vectorsoup',  'old',    'unformatted',0
 11,'W.nnkp',          'old',    'formatted',  0
 12,'W.eig_sig2',      'unknown','formatted',  0
 18,'W.vspup',         'old',    'formatted',  0
 20,'W.struct',        'old',    'formatted',  0
 50,'W.energysoup',    'old',    'formatted',  0
 51,'W.fermiup',       'old',    'formatted',  0
  9,'./W.vectorsodn',  'old',    'unformatted',0
 19,'W.vspdn',         'old',    'formatted',  0
```

Unit order in the def file does not matter, the reader is a `select case` over
whatever it finds (`main.F:120-152`).

### 6.5 Commands

```sh
cd tests/regress_cases/bi2d_w2w
cp ../../Bi2D/Bi2D.nnkp w2wBH.nnkp
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
path/to/SRC_w2w/w2wc w2wBH.def

python3 ../../make_small_win.py w2wBH-hyb 36 "6 6 1" w2wBH-hyb.win
# add the spinors flag (make_small_win.py does not emit it).
# Insert BEFORE line 1: make_small_win.py's first line is `num_bands = 36`,
# not `write_hr=.true.` as in the archived .win files.
sed -i '' '1i\
spinors=.true.
' w2wBH-hyb.win
wannier90.x w2wBH-hyb
```

Or use `tests/emit_guiding_win.py`, which detects the spinor case from the
anchors file and writes `spinors = .true.` itself (§7.1).

### 6.6 What the run prints

From `tests/regress_cases/bi2d_w2w/w2wBH.outputwf`:

```
 band window = [9, 120]
   cliff (support-otsu):   6 mask channels of   6 selected; mask columns =   44
   Channels selected:    6
   Total projections:   44
 HYBR capacity gate: need=26 (froz_emax=  2.000 eV)  margin=4  auto num_wann=36  (3 of 6 ranked channels)
 SIGMA loop: spinor-summed overlap + per-sigma amn column blocks (Bi2D.vectorso / Bi2D.vectorsodn)
 HYBR: P from all 88 amn columns
 HYBR: Lowdin P (full-cell): max =   1.0000   mean =   0.7857
 HYBR: MV b-weight residual =   3.42E-08
 HYBR: regime=per-state  pool[1..112] npool=112  frozen/k mean= 25.83  demand(worst k)=26
 HYBR: Omega_I =      59.53119273
  SCDM sigma_min(worst k) =   7.2330E-05
  auto gauge (tournament): winner seed k=1  Omega =       144.5568
 HYBR: hand-off max sigma(M) =     0.999693  OK (contractive)
 HYBR: wrote w2wBH-hyb.{amn,mmn,eig}
```

`Total projections: 44` is the **spatial** column count; `P from all 88 amn
columns` is `2 × 44`, the spinor total. Seeing the factor of two in that line is
the quickest confirmation that the sigma loop ran. If the def file
lacks unit 9, the `SIGMA loop:` line is absent and the column count is 44. The
run succeeds and produces a wrong model.

### 6.7 `num_wann` must be even (but nothing enforces it)

There is **no parity check anywhere in the source**. An odd `num_wann` under
the sigma loop will run to completion and produce a hand-off.

It will nonetheless be wrong whenever the bands are Kramers degenerate. With
both time-reversal and inversion symmetry present, every SOC band is doubly
degenerate at every k; an odd `nw` cuts a degenerate pair in half, which no
smooth gauge can do, and the result is a model whose bands cannot follow the
DFT bands. Keep `nw` even. The auto sentinel (`num_wann = 0`) does this for you
because `sf = 2` (`main.F:533, 562`); an explicit value is on you.

Systems without inversion (polar slabs, chiral crystals) do split the pairs, so
the argument is weaker there, but the `amn` still has `2·Nproj` columns and an
even `nw` remains the natural size.

### 6.8 `spinors = .true.` in the `.win`

Add it. In the hand-off setting Wannier90 takes its initial gauge entirely from
the `amn`, so the flag does not change the minimization start. It tells
Wannier90 that a projections block entry stands for a **spin-doubled**
spatial orbital, and downstream spin-resolved post-processing reads it.

With a `random` projections block the flag is accepted directly. Every SOC
hand-off `.win` in `tests/regress_cases/bi2d_w2w/` pairs `spinors=.true.` with a
`random` projections block. The block below is `w2wBH-hyb.win` specifically (the
`nw = 12` wins `w2wB12-hyb.win`, `w2wBI-hyb.win` and `w2wBV-hyb.win` have
`num_bands = num_wann = 12`, and `w2wBV-hyb.win` is a plot restart, prefixed
with `restart = plot` / `bands_plot = true` / `bands_num_points = 40` and
carrying `num_iter = 0`):

```
write_hr=.true.
spinors=.true.
num_bands = 36
num_wann  = 36
num_iter = 10000
conv_window = 5
conv_tol = 1e-8
guiding_centres = F
begin projections
random
end projections
```

Wannier90 3.x runs this file (`w2wBH-hyb.wout`).

If you instead supply **real** projection cards, Wannier90 expects `nw/2` of
them, since it doubles each. `emit_guiding_win.py:107-119` implements exactly
that, emitting each unique `(atom, l, m)` once and falling back to
`random` + `guiding_centres = F` if the anchor set is not closed under the sigma
swap. `tests/check_anchors.py` verifies that closure.

`tests/Rh3Si4/STATE.md:32-34` records dropping the `spinors` flag on the 120-WF
model for the opposite reason: that Wannier90 build *demanded a full explicit
projections block* once `spinors` was set, and only a `random` dummy block was
available. So the flag was dropped and the `random` block kept, because the `amn` carries
the gauge regardless.

### 6.9 Kramers degeneracy: what to look for in the `.wout`

For the 36-WF Bi2D model (`w2wBH-hyb.wout:430268-430305`), the final spreads
pair off exactly. The file lists the WFs in strict index order 1, 2, 3, …, 36;
the excerpt below is **re-ordered** to put each partner next to its mate, so do
not expect a literal `diff` against the file (WF 1 is at `:430269`, its partner
WF 20 at `:430288`):

```
 WF centre and spread    1  (  1.087757,  0.262768,  0.599503 )     2.06678859   <- :430269
 WF centre and spread   20  ( -1.087757, -0.262768, -0.599503 )     2.06678859   <- :430288
 WF centre and spread    2  ( -1.088805, -0.260951, -0.599460 )     2.06362158   <- :430270
 WF centre and spread   19  (  1.088805,  0.260951,  0.599460 )     2.06362158   <- :430287
 ...
 Sum of centres and spreads (  0.000000, -0.000000, -0.000000 )    63.09238661   <- :430305
```

Every spread appears twice to 8 decimals, with partner centers at exactly `−r`,
and the center sum vanishes. That is the Kramers/inversion signature, and it is
the strongest available check that the sigma blocks were built and paired
correctly. Unpaired spreads on a centrosymmetric cell mean something went wrong
upstream, usually an odd `nw`, a missing unit 9, or a mask that was not
mirrored across the sigma blocks.

### 6.10 How to tell it worked (SOC)

Run the §4.5 checklist, plus:

| Check | Bi2D nw=36 | Meaning |
|---|---|---|
| `SIGMA loop:` line present | yes | Unit 9 was read; otherwise this is a single-component run |
| `P from all N amn columns` | 88 = 2 × 44 | Column blocks were assembled |
| `num_wann` even | 36 | Kramers pairs intact |
| Wannier90 `Omega I` | 59.531191949 vs HYBR 59.53119273 | Subspace preserved |
| Final spread per WF | 63.0924/36 = 1.75 Å² | Localization |
| Spreads pair to ~8 decimals | yes | Kramers degeneracy respected |
| `check_anchors.py <out>.anchors projmap.dat` | exit 0 | Anchor set is shell-complete and sigma-paired |

### 6.11 The legacy route: per-spin `w2w` + `w2waddsp`

HYBR does not need this route, because the sigma loop replaces it
(`SRC_w2w/main.F:826-828` says so). It is described here because it is
what the WIEN2k wrappers do and what older data was produced with.

`x w2wc -so -up` then `x w2wc -so -dn` produce `case.amn{up,dn}` and
`case.mmn{up,dn}` (`x_lapw:3061-3100`; note that for non-spin-polarized SOC
the script clears `$updn` for the `inwf`/`vsp`/`fermi` units but keeps it for
`amn`/`mmn`, so the two runs read the **same** `case.inwf`). Then
`x w2waddsp`, which is what `x wannier90 -so` calls
(`SRC/wannier90_lapw:146-152`).

`w2waddsp` does **entrywise summation of both files**:

* `m1 = m1 + m2` for the overlaps (`SRC_trig/w2waddsp.f:71`); correct: this is
  the spinor overlap;
* `a1 = a1 + a2` for the projections (`SRC_trig/w2waddsp.f:104`), and this is
  where the rank deficiency comes from.

The `w2waddsp` built from this release's `SRC_trig` (not copied by `make install`; `x w2waddsp` runs WIEN2k's own 2.0 program, which checks nothing and writes 2E18.5) validates array dimensions on both paths, and on the **mmn path only** also
the neighbor lists and G-shifts: `w2waddsp.f:67-68` (dimensions) plus `:69-70`
(`if (any(neigh1 /= neigh2) .or. any(gsh1 /= gsh2))`). The amn path checks
dimensions alone (`w2waddsp.f:102-103`) and then sums at `:104`; there is no
neighbour/G-shift equivalent there. It writes binary `.bin` sidecars for both
outputs.

**The zero-coefficient block construction.** To make the *sum* behave like a
*concatenation*, write two manual `.inwf` files with `2·N` projection cards
each: in the up file, cards `1..N` carry coefficient `(1.000, 0.000)` and cards
`N+1..2N` carry `(0.000, 0.000)`; in the dn file the pattern is reversed. Since
`l2Amn` multiplies each contribution by the card's coefficient
(`SRC_w2w/modules_rc.F:945, 954, 966, 974`), a zero coefficient yields an
exactly zero column, and the sum reproduces the block layout.

`tests/Rh3Si4/w2w_hybrid/w2w_hybrid.inwf{up,dn}` are a working pair; 160 cards
each, blocks of 80. Head of the up file:

```
AMN
  193  544
  3  160
1
  1  2  -2   1.000 0.000
1
  1  2  -1   1.000 0.000
...            (cards 1-80: coefficient 1)
1
  1  2  -2   0.000 0.000
1
  1  2  -1   0.000 0.000
...            (cards 81-160: coefficient 0)
```

and the dn file has the two blocks swapped. The manual `.inwf` format is
`MODE / bmin bmax / LJmax Nproj / {NY, then NY lines of "iat l m Re Im"} × Nproj`
(`SRC_w2w/util_w2w.F:1699-1721`).

With non-spin-polarized SOC, `x w2w` reads the
plain `case.inwf` for both runs (`x_lapw:3075-3086`), so you must swap the
file between the two invocations, or write your own def files pointing unit 5
at `case.inwfup` / `case.inwfdn`. The summed `case.amn`/`case.mmn` can
then be fed to the **standalone** `w2w_hybrid_acc` for Stage 2. This is how
the `fesoc` and `bi2d` regression cells work
(`tests/regress_hybrid.sh:44, 46`).

---

## 7. Two upstream defects to work around

Both are in the helper scripts that ship with wien2wannier, not in the Stage-2
code. Both are triggered by the documented workflow.

### 7.1 `write_win` emits `guiding_centres = T` with a dummy projections block

`SRC_trig/write_win_backend.f:280-292` writes:

```
!!!! Dummy `projections' block for guiding centres !!!
guiding_centres = <T/F>
begin projections
   1:s
   2:s
   ...
end projections
```

The cards are produced by `print fmt_centers, centers(i)` with
`fmt_centers = '(I4, ":s")'` (`write_win_backend.f:26`, `289-291`), i.e. an
atom index and a literal `s`. The `guiding_centres` value defaults to
`inwf%Nproj > 0` (`write_win_backend.f:102`) (so **true for any real input
file**) and is only turned off if some projection is not uniquely centered
(`write_win_backend.f:104-113`).

Why this is harmful: with an `amn` supplied, Wannier90 takes the initial gauge
from the `amn` and never consults the projections block for it. But
`guiding_centres = T` keeps a guidance term active *during* minimization, and
that term is evaluated against the block's centers. The block says every
Wannier function is an `s` orbital on atom `i` in file order, which bears no
relation to the actual hand-off, whose columns are the QRCP-selected anchors of
a `2·Nproj`-column candidate set. The minimization is then pulled toward
meaningless centers. Compounding this, `num_wann` is set by a `merge`: it takes
`inwf%Nproj` when `Nproj > 0` and falls back to `num_bands` otherwise
(`write_win_backend.f:76-77`,
`merge(inwf%Nproj, num_bands, inwf%Nproj>0)`), i.e. for any real input file it
is the candidate count, not `nw`. A `num_wann` passed explicitly on the command
line does override it, via the `optarg` key/value loop
(`write_win_backend.f:83-95`).

**Correct settings.** Two options exist. The scripts named below are in the development tree, not in this release:

*Safe:* `tests/make_small_win.py`; `num_bands = num_wann = nw`,
`guiding_centres = F`, `begin projections / random / end projections`, lattice
and k-points taken from the `.nnkp`. This is what every validated hand-off
`.win` in `tests/regress_cases/` uses.

*Better:* `tests/emit_guiding_win.py`, reconstructs **real** centers by
composing `<out>.anchors` (one candidate-column index per WF, written by
`SRC_w2w/hybrid_core.F:1810-1817`, the
`open(newunit=ua, file=trim(outp)//'.anchors', ...)` block at `:1812`, inside
`choose_gauge`, which starts at `:1743`), `projmap.dat` (column → global atom, `l`,
`m`, `rf`, written at `SRC_w2w/main.F:493-509`) and the WIEN2k `struct`
positions, then emits genuine `f=x,y,z:l=L,mr=M` cards and sets
`guiding_centres = T` legitimately. It detects the spinor case from whether any
anchor index exceeds the spatial column count, emits `nw/2` cards and
`spinors = .true.` in that case, and falls back to `random` +
`guiding_centres = F` if the card count does not match
(`emit_guiding_win.py:96-125`).

```sh
python3 tests/emit_guiding_win.py w2wH-hyb 18 "8 8 8" w2wH-hyb.win \
        w2wH-hyb.anchors projmap.dat GaAs-WANN.struct --num-iter 10000
```

Never leave `guiding_centres = T` in place in the `write_win` output.
Turning it off is the minimum fix.

### 7.2 `write_inwf`'s hybrid orbital sets are not orthonormal

`SRC/write_inwf_lapw:251` reads:

```python
s2 = np.sqrt(2); s3 = np.sqrt(2); s6 = np.sqrt(2); s12 = np.sqrt(2)
```

Only `s2` is right. `s3`, `s6` and `s12` are meant to be `sqrt(3)`, `sqrt(6)`
and `sqrt(12)`; they are all `sqrt(2)`. Every hybrid built from them is
therefore neither normalized nor mutually orthogonal.

For `sp3d2` (`write_inwf_lapw:298-303`) the resulting Gram matrix over the six
hybrids, computed with the constants as written, is

```
[[1.75 0.75 0.75 0.75 0.   0.  ]
 [0.75 1.75 0.75 0.75 0.   0.  ]
 [0.75 0.75 1.75 0.75 0.   0.  ]
 [0.75 0.75 0.75 1.75 0.   0.  ]
 [0.   0.   0.   0.   1.5  0.5 ]
 [0.   0.   0.   0.   0.5  1.5 ]]
```

instead of the identity. Blast radius, by inspection of which constants each set
uses:

| Set | Constants used | Status |
|---|---|---|
| `sp` (`:284-285`) | `s2` | correct |
| `sp2` (`:286-288`) | `s3`, `s6`, `s2` | **broken** |
| `sp3` (`:289-292`) | literal `/2` only | correct |
| `sp3d` (`:293-297`) | `s3`, `s6`, `s2` | **broken** (its `sp2` part) |
| `sp3d2` (`:298-303`) | `s6`, `s2`, `s12`, `s3` | **broken** |

What it costs you: a non-orthonormal projector set produces an `amn` whose
columns are linearly dependent to within the overlap, which degrades the SCDM
pivot conditioning (`sigma_min`) and biases the Löwdin projectability, since
`S⁻¹` is then far from the identity.

**Workaround.** Do not use `sp2`, `sp3d` or `sp3d2` from `write_inwf`. In HYBR
this is easy, because you should not be writing manual projection cards at all:
AUTO generates pure `Y_lm` channels, which are orthonormal by construction, and
the gauge stage builds whatever hybrids the physics wants out of them via SCDM.
`sp` and `sp3` are safe if you do need a manual file. If you specifically need
one of the broken sets, write the cards by hand with the correct
normalizations. The manual `.inwf` format takes an explicit complex
coefficient per `Y_lm` (`SRC_w2w/util_w2w.F:1712-1719`).

### 7.3 Related: `x wplotc` is a dead branch

Not a HYBR issue, but you will meet it when plotting the Wannier functions.
`x_lapw:547` has a stray `endif` that closes the command-alias `if`-chain
before the `wplotc` branch at `x_lapw:548-551`. csh silently skips the
orphaned `else if`, so `x wplotc` does **not** set `cmplx = c` and runs the real
`wplot`. The `w2wc` branch just above (`:544-546`) is inside the chain and works.
Call `wplotc` directly, or use `x wplot -c`.

---

## 8. One-page summary of the checks

Run these on every HYBR calculation, whatever the type.

**In `case.outputwf`:**

1. `MV b-weight residual` ≲ 1e-6.
2. Band window is what you expected (sentinels are resolved and printed).
3. `Channels selected` / `Total projections`, for SOC, confirm the `P from all
   N amn columns` line shows `2 ×` the spatial count and that `SIGMA loop:`
   appeared.
4. Capacity gate line says `OK` (or you consciously accepted `tight` /
   `capacity = warn`).
5. `demand(worst k)` ≤ `num_wann` in the post-mask line.
6. `Omega_I` is comparable to the recorded value for a similar cell, the
   9-cell reference table is in `04-methods-and-tuning.md` §6.2.
7. The disentangler converged, or hit the cap with a small `|dOmega|_last`.
8. `hand-off max sigma(M)` prints `OK (contractive)`.

**After Wannier90:**

8. Wannier90's `Omega I` equals HYBR's `Omega_I` to 7–8 digits. If not, the
   `.win` is wrong; check for stray `dis_*` keywords and that
   `num_bands = num_wann`.
9. Final spread per WF is not pathological.
10. SOC only: spreads pair to ~8 decimals with partner centers at `−r`, and the
    center sum vanishes.
11. Spin-polarized only: the two spins' windows and channel selections match
    (differences in `frozen/k`, `Omega_I` and spread are expected and physical).

**Fidelity, if the model is going anywhere downstream:** compute the
Vitale-style band distance with `tools/band_compare.py`, which reports `eta` (RMS,
Fermi-weighted to `E_F + 2 eV`) and `eta_max` against the DFT bands on a k-path.
Plot bands against **k-distance**, never band index; `klist_band` steps are
non-uniform and an index axis fabricates jitter.
