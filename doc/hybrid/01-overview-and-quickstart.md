# HYBR: hybrid pre-Wannierization, overview and quickstart

> **Public command.** Use `wien2wannier` for normal integrated AUTO/HYBR
> workflows; its complete workflow contract is in
> `00-public-wien2wannier-cli.md` and `wien2wannier --help`. `w2w_hybrid`
> remains the advanced pre-generated-file backend; its native option reference
> is `w2w_hybrid --help`. Source line citations (`file.F:NN`) and paths under
> `tests/` in this guide refer to the development tree. The line numbers may
> be offset in the released sources, and the `tests/` cases are in the data
> record of the paper, not in this archive; every input deck quoted from them
> is reproduced in the text. The symmetry-adapted features of the development
> tree are not part of this release. `x_lapw:NN` cites the `x` script of WIEN2k
> 24.1 as distributed (`$WIENROOT/SRC/x_lapw`); this package does not contain
> that script.


This document describes `MODE = HYBR`, the hybrid pre-Wannierization mode built into
`w2w`, and its standalone twin `w2w_hybrid`. It assumes you know WIEN2k and
Wannier90 but have never used this mode.

Every statement below is anchored to the source with `file:line`. Paths are
relative to the repository root.

---

## 1. What the mode does

The classical wien2wannier route is: pick a band window and a set of projection
cards by hand, let `w2w` write `case.amn` / `case.mmn` / `case.eig`, and let
Wannier90 do both the disentanglement (`dis_*` keywords) and the gauge
minimization.

HYBR moves the subspace selection out of Wannier90 and into `wien2wannier`,
where the all-electron information (per-state projectability onto muffin-tin
channels, with full-cell tails) is available. (Terminology used
throughout these guides: `wien2wannier` is this modified package, whose
hybrid extension does the disentanglement -- in-pipeline via `MODE = HYBR`,
or standalone via `w2w_hybrid`; `w2w`/`w2wc` names the matrix-element
executable and its invocation. The original `w2w` program never performed
disentanglement.) HYBR runs as one extra
stage inside the ordinary `x w2w` run, immediately after `l2Amn`
(`SRC_w2w/main.F:852-855`, subroutine `hybrid_stage2` at
`SRC_w2w/main.F:874-1034`).

### The Stage-2 pipeline

| # | Step | Where |
|---|------|-------|
| 1 | **AUTO channel selection**, rank every `(atom, l, rf)` channel by projectability-weighted charge, choose the candidate set, write `projmap.dat` | `main.F:362-514`, `modules_rc.F:1322-1920` |
| 2 | **Capacity gate**; `need = max_k count(E − E_F ≤ froz_emax)`; validate or auto-derive `num_wann` *before* any expensive `amn`/`mmn` work | `main.F:516-620` |
| 3 | **Löwdin per-state projectability**; `P_i(k) = diag(A S⁻¹ A†)` on the full-cell (tail-extended) `amn`, exactly bounded in `[0,1]` | `main.F:903-930`, `hybrid_prep.F:25-62` |
| 4 | **Band masks**, pool selection (`pool_factor·nw` contiguous bands), frozen set below `froz_emax`, per-k feasibility repair, pool guard | `hybrid_prep.F:499-636` |
| 5 | **Ω_I disentanglement**, frozen-constrained Gauss–Seidel minimization of the gauge-invariant spread over the pool | `main.F:1001-1003`, `hybrid_core.F:333-662` |
| 6 | **Rotate to the projected-H eigenbasis**, diagonalize `V†HV`, giving `nw` downfolded bands per k | `hybrid_core.F:666-706` |
| 7 | **Gauge choice**; SCDM/QRCP anchors, optional parallel transport and seed tournament, smoothing sweeps | `hybrid_core.F:1743-2108` |
| 8 | **Isolated hand-off**, write `<case>-hyb.{amn,mmn,eig}` plus `<case>-hyb.anchors` | `main.F:1015-1033`, `hybrid_io.F:99-109` |

Steps 1–2 also happen in plain `MODE = AUTO`; steps 3–8 are what HYBR adds.

---

## 2. Why the hand-off is `Nb = Nw`

This is the central conceptual point, and it changes how you write the `.win`.

The disentanglement in step 5 is done in `wien2wannier` (the hybrid stage
of this modified package, which the original `w2w` never did), on the *pool* of
`npool ≈ pool_factor·nw` bands, with a frozen manifold enforced per k
(`hybrid_core.F:333-662`). Its output is a per-k isometry `V_k` of shape
`npool × nw`. Step 6 rotates that subspace into the eigenbasis of the projected
Hamiltonian, producing exactly `nw` downfolded band energies per k
(`hybrid_core.F:666-706`). Step 7 then chooses a gauge `A_k` of shape `nw × nw`.

The writer emits an overlap matrix `M^{(k,b)}` of shape `nw × nw`
(`hybrid_io.F:76-97`), an `amn` whose header is literally `nw  nk  nw`
(`hybrid_io.F:64`), and an `eig` with `nw` entries per k (`hybrid_io.F:49-53`).

So the problem handed to Wannier90 is an **isolated `nw`-band manifold**:

```
num_bands = num_wann = nw
```

Wannier90 therefore has **nothing left to disentangle**. It performs a pure
gauge minimization (the Ω_OD/Ω_D part of the Marzari–Vanderbilt spread).
Ω_I is already fixed by the subspace
`wien2wannier` selected, and no `dis_win_min` / `dis_win_max` / `dis_froz_min` /
`dis_froz_max` / `dis_num_iter` keyword has any meaning for the hand-off files.
Putting them in the `.win` is at best inert and at worst an error.

Two consequences follow:

* The **quality of the model is decided in `wien2wannier`**, not in Wannier90. If the
  bands are wrong, the fix is `froz_emax`, `num_wann`, `pool_factor` or the band
  window, not `dis_*`.
* The hand-off `.eig` holds the **downfolded** eigenvalues (eigenvalues of the
  projected Hamiltonian in the selected subspace), already shifted to eV
  relative to `E_F`; HYBR passes `e_fermi = 0` to the writer
  (`main.F:1027-1028`, and `eig_h` is built Fermi-shifted at
  `main.F:897-901`). They are not raw DFT eigenvalues.

A sanity gate is printed at the end: the largest singular value over a
**sampled** subset of the `M^{(k,b)}` blocks must be ≤ 1 + 1e-6, since the
hand-off overlaps are between orthonormal subspaces (`hybrid_io.F:16-41`,
printed at `main.F:1029-1032`). The sampling is every `max(1, nk/64)`-th
k-point, all `b` at those k (`hybrid_io.F:14-15`, `:27-28`), roughly 64
k-points, so for the 512-point GaAs cell about 12.5% of the blocks. A
contractivity violation at an unsampled k would not be caught.

---

## 3. Two ways to run it

For a normal calculation, the public driver prepares and validates the
integrated route described below:

```sh
wien2wannier --input /path/to/case --seed model --k-grid 8 8 8
```

Its default `all` stage writes the isolated `model-hyb` handoff. It infers
`num_wann` through the integrated `num_wann = 0` capacity route using the
reported `--froz-emax` ceiling, selects `w2w` or `w2wc` including SOC def units,
and refuses a mesh/NNKP mismatch. The direct commands in (a) and (b) remain the
lower-level interfaces for manual cards and backend-specific experiments.

### (a) In-memory one-shot; `MODE = HYBR` inside `w2w`

Put `HYBR` on line 1 of `case.inwf` and run `w2w` as usual. The `amn`, `mmn`
and `eig` never leave memory before Stage 2 consumes them
(`main.F:975-996`), so there is no re-parse of multi-GB text files. Stage 2 runs
on the MPI root only (`main.F:852`).

The `wien2wannier` driver (`00-public-wien2wannier-cli.md`) runs this route by
default: with `--stage all` and `--method hybrid` it writes `HYBR` on line 1 of
the `.inwf`. A normal calculation therefore uses this route through the driver,
and running `w2w` directly with `HYBR` in `case.inwf` is the lower-level form of
the same route. It is also the only route that obtains the tail-extended
projectability directly: `l2Amn` computes the
muffin-tin `amn` and the full-cell `amn` side by side in one plane-wave sweep
and exports both (`modules_rc.F:608-641`), so no second `w2w` run and no
`W2W_PAMN` sidecar file are needed.

### (b) Standalone; `w2w_hybrid` / `w2w_hybrid_acc`

`SRC_w2w/w2w_hybrid.F`, built with `make -f Makefile.hybrid` (`w2w_hybrid`) or `make -f Makefile.hybrid acc` (`w2w_hybrid_acc`)
(`SRC_w2w/Makefile.hybrid:82-85`). It reads `<seed>.amn`, `<seed>.mmn`,
`<seed>.eig`, `<seed>.nnkp` from disk (`w2w_hybrid.F:166-174`; binary
`.bin` sidecars are used transparently when present, `hybrid_files.F:11-32`)
and takes named `--key value` options (`--seed`/`--nw` mandatory; `w2w_hybrid --help` is the authoritative reference; the bare `<case> <N>` shorthand also works).

Use it when you want to:

* sweep Stage-2 parameters (pool factor, gauge mode, mixing, seed policy)
  against a **fixed** `amn`/`mmn`, without recomputing the overlaps;
* reach controls the `.inwf` does not expose: `init_transport`, `nofroz`,
  `--p-floor` and the disentangler's `--check-every` (HYBR fixes it at 10).
  `colored`, `adaptive`, `anderson`, `sweeps`, `rel_tol` and the absolute
  `tol` are HYBR keywords; the standalone passes the literal `1.0e-9` for
  `tol` and has no option for it;
* run Stage 2 on files produced by an older `w2w`.

The two paths call the *same* modules (`build_masks`, `disentangle_frozen`,
`rotate_to_subspace`, `choose_gauge`, `write_handoff`) and were validated to
produce bit-identical hand-offs. The parity harness is
`tests/regress_cases/run_parity_cell.sh`; it tests the equivalence with

```
W2W_PAMN=<tails-amn> W2W_SEEDK=opt16 \
  w2w_hybrid_acc --seed <seed> --nw <nw> --out <out> --rescale lowdin \
                 --regime per-state --froz-emax 2.0 --gauge auto8 \
                 --np-mask <mask_cols> --admit-all
```

The two environment variables are part of that command
(`tests/regress_cases/run_parity_cell.sh:24-26`). `W2W_PAMN` is what makes the
standalone use the tail-extended full-cell `amn` as its projectability source;
without it the tool falls back to the plain MT `amn` (`amnP = amn`,
`w2w_hybrid.F:180-190`), giving a different Löwdin `P`, different masks and a
different hand-off, i.e. not bit-identical. `W2W_SEEDK=opt16` is likewise
needed to match HYBR, which passes `seedk_in=trim(inwf%hyb_seedk)`
(`main.F:1025-1026`) with the type default `opt16` (`util_w2w.F:1533`).

which matches HYBR's defaults `sweeps = 20000`, `solver = cg`,
`rel_tol = 1e-10`, `check_every = 10` (fixed), `colored = false`,
`gauge = auto8`, admit-all (`main.F:949`); HYBR's `adaptive = true` and
`relax = -2.5` act only with `solver = gs`.

The **validated production recipe** for the standalone (this is a tuned choice,
not a default) is the form used by `tests/regress_hybrid.sh`:

```
W2W_PAMN=<seed>.amn.tl W2W_SEEDK=opt16 W2W_SEEDLOSS=smooth \
  w2w_hybrid_acc --seed <seed> --nw <nw> --out <out> --rescale lowdin \
                 --regime per-state --froz-emax 2.0 --sweeps 4000 \
                 --colored --adaptive --relax -1.9 --gauge <gauge> \
                 --np-mask <np> --admit-all --solver gs
```

Note `--pool-factor` is *not* passed there, so it falls back to its default 1.5.

---

## 4. Prerequisites

Before you can run HYBR you need, in order:

1. **A converged WIEN2k SCF run** for the case.
2. **A working directory prepared by `prepare_w2wdir`**, which copies
   `struct`, `in1(c)`, `vsp`, `vns`, … and extracts `case.fermi` from the SCF
   file (`SRC/prepare_w2wdir_lapw:83-125`, `:240-300`). `case.fermi` is
   mandatory, it is unit 51 in the `w2w` def file
   (`x_lapw:3099`) and is read at `main.F:258`.

   **This step comes first.** `prepare_w2wdir` copies only the `must`/`should`/
   `can` suffix lists; `struct`, `in1(c)`, `vsp`, `vns` (plus `klist_band`,
   `output1`, `spaghetti_ene`, `vorb`), and never `case.energy`,
   `case.vector` or `case.klist`
   (`SRC/prepare_w2wdir_lapw:241`, `:291-307`). Everything from here on
   therefore has to happen **inside** the new subdirectory, whose case name is
   the directory name: `x w2w` run there writes a def file pointing at
   `W.vector` / `W.energy` (`x_lapw:3093`, `:3098`), so `lapw1` must have
   been run there too. This is also the order of the user guide
   (`doc/wien2wannier_userguide.tex:165-170`, `:202-210`).
3. **A full-BZ k-mesh**, generated in that subdirectory. Wannier90 needs the
   complete mesh, not the irreducible wedge: `x kgen -fbz` (this is what
   `init_w2w_lapw` does, `SRC/init_w2w_lapw:219`).
4. **`lapw1` (and `lapwso`, if applicable) run on that mesh, in that
   subdirectory**, producing `case.energy` and `case.vector`. Both are read by
   `w2w` (`main.F:175-226`, `main.F:338-346`).
5. **A `case.win`** (from `write_win`), and then the Wannier90 preprocessing
   step that turns it into the b-vector list:

   ```
   x wannier90 -pp        # runs wannier90.x -pp, writes case.nnkp
   ```

   (`SRC/init_w2w_lapw:315`.) `w2w` reads the k-point count and the
   `nnkpts` neighbor list from `case.nnkp` (`main.F:229-248`), and HYBR
   re-reads the same file for the Cartesian b-vectors after the unit is closed
   (`main.F:140-143`, `main.F:933`).

Real or complex: `x w2w` picks `w2w` or `w2wc` from the struct's inversion
symmetry (`x_lapw:496-500`); force it with `x w2w -c`. The variant must
match the one your `lapw1` used.

Spin–orbit and spin-polarized cases need extra def-file units (9 = second-spin
vector, which switches on the in-memory sigma loop; 19 = second-spin `vsp`)
that `x w2w` does **not** write; see `main.F:119-129`, `main.F:689-761`. Those
are outside the scope of this quickstart.

---

## 5. Worked example: GaAs, non-magnetic, no SOC

This is the cell in `tests/regress_cases/gaas_w2w/`, reproduced here in full.
An 8×8×8 mesh, band window 6–90, `num_wann` left to the capacity gate.

### 5.1 Setup

```sh
cd /path/to/GaAs
prepare_w2wdir_lapw W            # copies struct/in1/vsp/vns, writes W/W.fermi
cd W                             # everything below happens HERE, case = W
x kgen -fbz                      # full BZ, 8x8x8 -> 512 k-points
write_win_lapw                   # or hand-write W.win
x wannier90 -pp                  # wannier90.x -pp -> W.nnkp
x lapw1                          # on the full mesh -> W.energy, W.vector
```

`prepare_w2wdir` runs **first**: it does not copy `case.energy` /
`case.vector` / `case.klist` (`SRC/prepare_w2wdir_lapw:241`, `:291-307`), so
running `x lapw1` in the parent directory would leave
`GaAs.energy` / `GaAs.vector` where the later `x w2w`, executed inside `W/`,
with a def file naming `W.vector` / `W.energy` (`x_lapw:3093`, `:3098`) cannot see them.

### 5.2 `W.inwf`, the complete input file

```
HYBR
  6  90
  3
0 0 T T -1
0 1.0 per-state
froz_emax = 2.0
```

Line by line (parser: `util_w2w.F:1544-1723`):

* `HYBR`, mode word (`util_w2w.F:1577-1583`).
* `6  90`, band window `bmin bmax`. Sentinels are available in HYBR:
  `bmin = 0` runs the gap-walk auto-bmin, `bmax = 0` means "all bands present
  at every k", `bmax = -1` means "highest band whose bottom is below
  E_F + 40 eV" (`main.F:262-304`).
* `3`; `LJmax` for the plane-wave overlap expansion (`util_w2w.F:1607`). 3 in
  every validated test input.
* `0 0 T T -1`; AUTO config: `threshold max_proj write_proj p2_weight
  support_ratio`. `threshold = 0` selects cliff mode (breadth for the gauge,
  Otsu cut for the physical/mask boundary); `max_proj = 0` means no hard column
  budget; `write_proj = T` writes the diagnostics; `p2_weight = T` enables the
  erfc energy weight in the channel ranking; `support_ratio = -1` selects the
  algorithmic support cut. This exact line is the validated production choice,
  not a set of defaults. The type defaults are `0.10 / 0 / F / T / 0.5`
  (`util_w2w.F:1497-1502`).
* `0 1.0 per-state`; `num_wann mix regime`. `num_wann = 0` is the
  capacity-driven auto sentinel (`main.F:551-584`); `mix = 1.0` switches the Zd
  relaxation off (default `-2.5`, scheduled; it acts only with `solver = gs`);
  `per-state` is the validated mask regime, versus the default
  `auto`.
* `froz_emax = 2.0`, keyword tail. The default is `auto` (walk down from 8 eV);
  `2.0` eV is the validated choice for this class of cell, and an explicit
  numeric value is **required** whenever `num_wann = 0`
  (`main.F:537-538`, `util_w2w.F:1664-1672`).

Other tail keywords, all optional, with their defaults:
`pool_factor = 1.5`, `margin = 0` (auto: `max(4, ceil(0.04·need))`),
`gauge = auto8`, `seedk = opt16`, `seedloss = smooth`, `capacity = strict`
(`util_w2w.F:1673-1691`). An unknown key is a hard error, and the tail stops at
the first blank line (`util_w2w.F:1652-1653`, `:1692-1693`); do not leave a
blank line in the middle of it.

### 5.3 Run

```sh
x w2w                            # add -c if your case is complex
```

`x w2w` writes the def file itself (`x_lapw:3088-3099`). The test tree uses
an equivalent hand-written def:

```
  5,'hybr.inwf',           'old',    'formatted',  0
  6,'w2wH.outputwf',       'unknown','formatted',  0
  7,'w2wH.amn',            'unknown','formatted',  0
  8,'w2wH.mmn',            'unknown','formatted',  0
 10,'../GaAs-WANN/GaAs-WANN.vector','old','unformatted',0
 11,'w2wH.nnkp',           'old',    'formatted',  0
 12,'w2wH.eig',            'unknown','formatted',  0
 18,'../GaAs-WANN/GaAs-WANN.vsp',   'old',    'formatted',  0
 20,'../GaAs-WANN/GaAs-WANN.struct','old',    'formatted',  0
 50,'../GaAs-WANN/GaAs-WANN.energy','old',    'formatted',  0
 51,'../GaAs-WANN/GaAs-WANN.fermi', 'old',    'formatted',  0
```

**The hand-off name is derived from the `.nnkp` filename** (unit 11), not from
the case name: `main.F:1012-1015` strips `.nnkp` and appends `-hyb`. With
`x w2w` that gives `case-hyb.*`.

### 5.4 Output files

| File | What it is | Written at |
|------|------------|------------|
| `case.outputwf` | the log: window resolution, channel ranking, capacity gate, all HYBR diagnostics | throughout |
| `case.amn`, `case.mmn`, `case.eig` | the **legacy** `Nb × nproj` files (`Nb = 85` here). Under HYBR the text `.amn`/`.mmn` are written only with `text_io = true` (the `.bin` sidecars always); not what you give Wannier90 in HYBR mode | `main.F:651-659`, `:763-824`, `modules_rc.F:2512-2517` |
| `case.mmn.bin` | binary sidecar for fast re-reads by the standalone tool | `main.F:784-822` |
| `case.amn.bin` | ditto, for the `amn`. In the ordinary (non-sigma) case documented here it is written by `l2Amn`; `main.F:826-848` writes it only for the SOC sigma-loop *combined* amn (`main.F:826` guards on `sigma2`) | `modules_rc.F:1008-1023` |
| `projmap.dat` | `amn` column → global atom, `l`, `m`, `rf`; carries `# mask_cols N` in cliff mode | `main.F:493-509` |
| `projectability.dat` | per-(band, k, atom, l) projectability (only if `write_proj = T`) | `modules_rc.F:1322-1376` |
| `autoprojection.diag` | channel ranking, alternative sets, paste-ready `.inwf` templates (only if `write_proj = T`) | `modules_rc.F:1885-1920` |
| **`case-hyb.amn`, `case-hyb.mmn`, `case-hyb.eig`** | **the isolated `Nb = Nw` hand-off; this is what Wannier90 gets** | `hybrid_io.F:99-109` |
| `case-hyb.anchors` | one anchor `amn` column index per Wannier function; join with `projmap.dat` for real guiding centres | `hybrid_core.F:1808-1817` |
| `case-hyb.masks` | the per-state masks behind the hand-off, one row per (k, pool state): `k band eps-E_F P frozen admitted`; `tools/plot_certificate.py` draws the window certificate from it | `hybrid_io.F` (`write_masks`) |
| `case-hyb.anatomy` | the auto-trust anatomy: per pool band `Pbar minP minE_rel_EF relevant cls` with the derived `p_gate`, `p_froz`, `plateau` in the header; `tools/plot_anatomy.py` draws it | `hybrid_io.F` (`write_anatomy`) |

---

## 6. Reading the diagnostics

Below is the actual `w2wH.outputwf` output from the GaAs cell above, annotated.

### 6.1 Window and selection

```
 MODE= Mmn Amn AUTO
 band window = [6, 90]
 P2 ranking weight: erfc onset E_F+ 5.0 eV, width  2.0 eV
```

`MODE=` prints `Mmn Amn AUTO` even in HYBR; HYBR sets all three flags
(`util_w2w.F:1577-1583`); it does not print a separate `HYBR` token
(`main.F:322-330`).

```
 WARNING: Nproj < Nb (entangled bands).
   Nproj =   84, Nb =   85
   Wannier90 will require disentanglement parameters
   (dis_win_min/max, dis_froz_min/max) in the .win file.
   Consider adjusting the threshold or band window.
```

**Ignore this warning in HYBR mode.** It is emitted by the AUTO selection
(`modules_rc.F:1761-1772`) and speaks about the *legacy* `amn`. In HYBR the
disentanglement is done by `wien2wannier`'s hybrid stage, and the hand-off has
no `dis_*` needs at all
(section 2).

### 6.2 The capacity gate

```
 HYBR capacity gate: need=11 (froz_emax=  2.000 eV)  margin=4  auto num_wann=18  (6 of 16 ranked channels)
```

`need` is the worst-k count of window states below `E_F + froz_emax`
(`main.F:541-548`). With `num_wann = 0` the gate walks the ranked channel list
adding whole channels (`(2l+1)·mult` columns each) until the total reaches
`need + margin`, and adopts that as `num_wann` (`main.F:551-584`). Here: 11 + 4
→ 18, using 6 of 16 candidate channels.

With an explicit `num_wann` you get one of three lines instead
(`main.F:592-616`):

* `... >= need+margin= ... OK`, nothing to do.
* `... >= need= ... but inside margin ... (tight -- consider N)`, valid, but the
  frozen manifold nearly fills the model. Raise `num_wann` to the suggested `N`
  if the larger model is acceptable.
* `... < need= ... -> raise num_wann to >= N or lower froz_emax, or set capacity = warn`,
  a hard error under the default `capacity = strict`. Choose one of the
  three actions it names. `capacity = warn` does **not** change any physics: it
  only downgrades the croak, and the frozen budget is then spent on the `nw`
  deepest sub-ceiling states per k (`hybrid_prep.F:596-610`).

On the **auto-nw path only** (`num_wann = 0`) there is one further advisory:

```
 HYBR capacity gate: WARNING window N bands < pool_factor*nw = X -> limited disentanglement slack; consider widening bmax
```

The pool wants `pool_factor·nw` bands of room inside the window
(`main.F:577-583`). Widen `bmax` (or use `bmax = -1`). This check sits inside
the `associate` block of the `if (inwf%num_wann == 0)` branch
(`main.F:550`, `:557`, `:576-584`, `:585`), so with an explicit `num_wann`
(`main.F:586`, the `else` branch) it can never appear; there is no pool-slack
advisory at gate time in that case.

Two gate-time hard errors, both with `num_wann = 0`: *"ranked candidate pool
exhausted below the frozen capacity"* (`main.F:564-566`) and *"auto-nw exceeds
the band window"* (`main.F:567-568`). Both mean the requested `froz_emax` is
incompatible with the candidate set or the window.

When `froz_emax` is `auto` or `none` (`main.F:530-531`, `:533`) the gate
behaves differently depending on `num_wann`:

* `num_wann > 0`, the whole gate block is **skipped**, and only the post-mask
  assertion in section 6.3 applies.
* `num_wann = 0`, it is a **hard error at the gate**:
  `HYBR: num_wann = 0 (auto) requires an explicit froz_emax = <eV> keyword`
  (`main.F:536-537`). The run aborts there and never reaches the post-mask
  assertion.

### 6.3 Projectability and masks

```
 HYBR: P from 62 physical of 84 spatial columns
 HYBR: Lowdin P (full-cell): max =   1.0000   mean =   0.7294
 HYBR: MV b-weight residual =   2.22E-16
 HYBR: admit-all (validated recipe default)
 HYBR: regime=per-state  pool[1..27] npool=27  frozen/k mean=  9.26  demand(worst k)=11  (capacity: nw >= demand)
```

* **`P from N physical of M spatial columns`**, the Otsu mask boundary from
  cliff mode restricted the projectability source; without a boundary you get
  `P from all M amn columns` (`main.F:909-926`).
* **`Lowdin P: max / mean`**; `max` should reach ≈ 1.0. A mean far below ~0.5
  says the selected channels do not span the window; widen the channel set
  (lower `threshold`) or narrow the window.
* **`MV b-weight residual`**, the least-squares residual of the
  Marzari–Vanderbilt b-weight condition (`hybrid_prep.F:98-167`). It should be
  at round-off (here 2e-16). A large residual means the `nnkp` shell set is
  degenerate; regenerate it with `x wannier90 -pp`.
* **`pool[a..b] npool=N`**, the band range handed to the disentangler.
  `max(nw+1, ceil(pool_factor·nw))` is only the intermediate `n_pool_max`
  (`hybrid_prep.F:248`); the pool is then widened to span the `nw`
  most-projectable seed bands, `npool = max(n_pool_max, s1-s0+1)`, and clamped
  by the window, `npool = min(npool, nb-s0+1)` (`hybrid_prep.F:263-265`). So
  `npool` can exceed `ceil(pool_factor·nw)` before the pool guard runs, and can
  also come out smaller near the top of the window. The pool guard, which
  force-includes every sub-ceiling state at both ends, may widen it further
  (`hybrid_prep.F:535-586`).
* **`demand(worst k)`**, the *pre-repair* frozen requirement. This is the
  number the capacity constraint applies to. The returned masks are always
  repaired to fit, so they can never flag the violation themselves
  (`hybrid_prep.F:513-518`).

If the demand exceeds `nw` you get, before the disentangler runs:

```
 HYBR: FROZEN CAPACITY VIOLATION: worst-k demand D > num_wann N -> raise num_wann to >= D or lower froz_emax
```

which croaks under `capacity = strict` (`main.F:961-973`). This is the check
that catches the `froz_emax = auto` case, where no pre-gate ran.

### 6.4 The disentangler

```
  disentangle sweeps=  2000  converged=F  |dOmega|_last=  1.05E-04  cap=  2000
  [WARN] disentangle hit the sweep cap WITHOUT converging
 HYBR: Omega_I =      21.83980920
```

The telemetry line is unconditional (`hybrid_core.F:626-631`): sweep count,
convergence flag, the last `|ΔΩ_I|` between exact checks, and the cap. HYBR's
cap is the `sweeps` keyword (default 20000; the run above used the former
fixed cap of 2000), with `tol = 1e-9`, `rel_tol = 1e-10` (keywords `tol`,
`rel_tol`) and `check_every = 10`.

`Omega_I` is the gauge-invariant spread of the selected subspace, in Å².
Lower is better; it is the one number to compare across HYBR runs of the same
cell, and it is *fixed* for Wannier90, the minimizer cannot improve it.

**What to do when you see the `[WARN]` line** (`hybrid_core.F:634-635`):

1. Look at `|dOmega|_last` first. At 1e-4 on a Ω_I of ≈ 21.8 (a relative change
   of ~5e-6 per check) the subspace is converged for practical purposes. The run
   did not meet the strict `1e-10` relative criterion. This is the GaAs
   case above, and its hand-off was accepted.
2. If `|dOmega|_last` is large (say ≳ 1e-2, or not small compared with `Omega_I`
   times your tolerance), the subspace is unconverged. Raise the
   `sweeps` keyword, or rerun Stage 2 with the standalone tool (`--sweeps`);
   `colored`, `adaptive` and `anderson` are HYBR keywords as well as
   standalone options (`--colored`, `--adaptive`, `--anderson`). Mixing (`--relax`) is *not* standalone-only: HYBR reaches it through
   field 2 of `.inwf` line 5 (`num_wann [mix] [regime]`,
   `util_w2w.F:1636-1645`, default at `:1523`) and passes it straight to the
   disentangler as `mixing=inwf%hyb_mix` (`main.F:1001-1003`).
3. Persistent non-convergence with a large residual usually means the pool is
   too tight or the frozen manifold conflicts with the optimization: raise
   `pool_factor`, widen the window, or lower `froz_emax`.
4. There is no progress trace in HYBR mode. `W2W_VERBOSE=1` is read only by the
   standalone (`w2w_hybrid.F:154-162`); HYBR passes a literal `.false.`
   (`main.F:1003`).

### 6.5 The gauge

```
 HYBR: gauge=auto8  seedk=opt16  seedloss=smooth
  SCDM sigma_min(worst k) =   1.4180E-03
  seed tournament: 17 candidates from 511 of 512 gate-passing k
    seed cand k=    329  Omega =        90.6728   (OD/D =      50.8373     17.9957 )
    ...
  auto gauge (tournament): winner seed k=329  Omega =        90.6728
  gauge: auto8 with 8 smoothing sweeps (in tournament)
  initial spread Omega =        90.6728   (I/OD/D =      21.8398     50.8373     17.9957 )
```

Under the default `gauge = auto8` the dispatcher compares the worst-k SCDM
smallest singular value to a hardcoded gate of 2e-3, the constant
`real(dp), parameter :: sig_gate = 2.0e-3_dp` declared at
`hybrid_core.F:2005`, applied at `:2016` (`if (best_svmin < sig_gate)`), with
the two report lines at `:2035-2039`. (Do not confuse this with the *separate*
gate of the seed-tournament candidate listing, `sgate5 = 2.0e-3_dp` at
`hybrid_core.F:1886`, which the `W2W_SEEDGATE` environment variable can
override, `:1936-1971`.)

* `sigma_min >= 2e-3` → the per-k SCDM gauge is kept, then 8 smoothing sweeps.
  The line reads `auto gauge: sigma_min ... >= ... -> per-k SCDM kept`.
  **In this branch `seedk` and `seedloss` do nothing**, even though the
  candidate list is still printed.
* `sigma_min < 2e-3` (the GaAs case, 1.42e-3) → the anchors are ill-conditioned
  somewhere. What happens next depends on `seedk`
  (`hybrid_core.F:2016-2017`):
  * `seedk` = `opt` or `opt<N>` (the default `opt16`, hence the GaAs output
    above) → the **seed tournament** runs and the winning transported gauge is
    used. This path `exit autotfl`s at `:2031`, so its smoothing happens
    *inside* the tournament and the trailing `smooth_gauge` call at `:2064` is
    skipped.
  * `seedk` = `best`, `median` or an explicit integer → **no tournament**: a
    single `gauge_best_k` seed followed by one `propagate_gauge`
    (`:2033-2037`), printing
    `auto gauge: sigma_min ... < ... -> transport from k=N`.

`initial spread Omega` with its `I/OD/D` decomposition is the spread of the
hand-off *before* Wannier90 touches it. `I` here equals the `Omega_I` from
step 5, as it must. `OD + D` is what Wannier90 will minimize.

```
 HYBR: hand-off max sigma(M) =     0.999044  OK (contractive)
 HYBR: wrote w2wH-hyb.{amn,mmn,eig}
```

If the `max sigma(M)` line reads `VIOLATION`, the hand-off is not a valid isolated model. Do not pass it to Wannier90, and report the case.

---

## 7. Handing off to Wannier90

Copy your original `.win` to the hand-off seedname and edit it. **The unit cell,
`mp_grid` and the `kpoints` block must stay byte-identical to the ones used for
`x wannier90 -pp`**: the hand-off `.mmn` is stored in exactly the neighbor
order of that `.nnkp`, and Wannier90 regenerates the same list from the same
mesh.

Minimal working `case-hyb.win` (this is the file from the GaAs cell, with
the 512 k-points elided):

```
write_hr = .true.
num_bands = 18
num_wann  = 18
num_iter = 10000
conv_window = 5
conv_tol = 1e-8
guiding_centres = F

begin projections
random
end projections

begin unit_cell_cart
   0.0000000   2.8266001   2.8266001
   2.8266001   0.0000000   2.8266001
   2.8266001   2.8266001   0.0000000
end unit_cell_cart

mp_grid = 8 8 8

begin kpoints
    0.00000000    0.00000000    0.00000000
    ... (all 512, unchanged from case.win)
end kpoints
```

Rules:

* `num_bands = num_wann = nw`. Nothing else is consistent with the files.
* **Delete every `dis_*` keyword.** `dis_win_min`, `dis_win_max`,
  `dis_froz_min`, `dis_froz_max`, `dis_num_iter`, all meaningless here.
* The `projections` block is still syntactically required, but Wannier90 takes
  its initial gauge from the supplied `.amn`, so `random` is sufficient.
* **Guiding centres.** With `random` projections you must set
  `guiding_centres = F`; pairing `guiding_centres = T` with dummy projection
  cards guides the minimization toward meaningless centers, the historical
  wien2wannier defect described in `02-calculation-types.md` §7.1. If you
  want the guidance term, generate *real* centers: each hand-off Wannier
  function carries an anchor column in `case-hyb.anchors`, `projmap.dat` maps
  that column to a global `(atom, l, m, rf)`, and the struct gives the position.
  `tests/emit_guiding_win.py` composes the three into a `.win` with genuine
  projection cards; then `guiding_centres = T` is meaningful.
* `fermi_energy = 0.0` if you set it at all, the hand-off `.eig` is already
  referenced to `E_F` in eV.

Then run the minimization. The hand-off files already use the `-hyb` seedname,
so run `wannier90.x` on it directly:

```sh
wannier90.x case-hyb
```

(The `x wannier90` wrapper defaults the seedname to `basename $PWD`
(`SRC/wannier90_lapw:38`), but `-f` overrides it (`:57`, `:104`) and `x` passes
argv straight through (`x_lapw:133`, `:3235`), so
`x wannier90 -f case-hyb` works just as well, no symlink needed. Note the
`-up`/`-dn`/`-so` flags do *not* enter the seedname: the name handed to
wannier90.x is `$file` alone (`SRC/wannier90_lapw:167`); those flags only pick
which files are symlinked to the unsuffixed names beforehand and which outputs
are renamed afterwards.)

Check in `case-hyb.wout` that the initial spread matches the
`initial spread Omega` printed by `wien2wannier` and that `Omega_I` does not move. If
it moves, the `.win` mesh does not match the one the `.mmn` was built on.

---

## 8. Troubleshooting: error messages

`croak` prints `w2w: <message>` to **standard error** and then exits
(`util_w2w.F:463-479`). It writes to neither `case.error` nor `case.outputwf`:
`case.error` is created at startup by `errflg(errfn,'Error in W2W')`
(`main.F:32`) and contains only that generic flag string, which is truncated
again by `ERRCLR` on success (`main.F:859`).

Not every message below is a `croak`. The entries citing `hybrid_core.F`,
`hybrid_prep.F`, `hybrid_files.F` and `w2w_hybrid.F` are Fortran `error stop`
(some preceded by a `write(*,*)`) in the hybrid modules; `hybrid_core.F:154`,
`:155`, `:751`; `hybrid_prep.F:50`, `:145`; `hybrid_files.F:194`;
`w2w_hybrid.F:110-112`. Those also go to stderr and never touch `case.error`.

| Message | Source | Meaning / fix |
|---|---|---|
| `unknown MODE in inwf, XXXX` | `util_w2w.F:1585` | Line 1 is not one of `AMN MMN BOTH AUTO AAMN HYBR`. |
| `invalid band index in inwf` | `util_w2w.F:1594-1595` | **AUTO-family modes only** (AUTO/AAMN/HYBR, guarded by `if (inwf%auto_mode)` at `util_w2w.F:1588`): `bmin < 0` or `bmax < -1`. Only `0` and `-1` are sentinels. The manual modes (AMN/MMN/BOTH) do no sentinel validation at all; `util_w2w.F:1599-1600` checks only `bmax < bmin`, so a negative `bmin` is accepted by the parser and fails later. |
| `bmin > bmax in inwf` | `util_w2w.F:1596-1600` | Both given explicitly and inverted. Raised in both branches. |
| `auto band window is empty` | `main.F:306` | The gap-walk `bmin` landed above the resolved `bmax`. Set `bmin` explicitly. |
| `band window is empty after bmax clamp` | `main.F:317` | `bmax` was clamped to the lowest per-k band count and fell below `bmin`. Your `lapw1` did not produce enough bands at every k; raise `emax` in `case.in1`. |
| `HYBR: num_wann must be >= 0 (0 = capacity-driven auto)` | `util_w2w.F:1646-1647` | Negative `num_wann` on line 5. |
| `HYBR keyword line without "="` | `util_w2w.F:1659-1660` | A tail line has no `=`. Note the tail ends at the first blank line, so stray text after a blank line is silently *ignored* instead. |
| `HYBR: unknown keyword "..."` | `util_w2w.F:1692-1693` | Typo protection. The tail accepts `froz_emax`, `pool_factor`, `margin`, `gauge`, `seedk`, `seedloss`, `capacity`, `edges`, `sweeps`, `relax`, `relax_late`, `taper`, `taper_late`, `kick`, `adaptive`, `colored`, `solver`, `text_io`, `tol`, `rel_tol` and `anderson`. In particular `rescale = ...` is **not** a HYBR keyword; writing it aborts the run. |
| `HYBR froz_emax: bad value ...` / `HYBR pool_factor: bad value` / `HYBR margin: bad value` | `util_w2w.F:1670-1678` | Unparseable number. `froz_emax` accepts `auto`, `none`, or a real. |
| `HYBR capacity: strict\|warn, got ...` | `util_w2w.F:1690` | Only those two values. |
| `HYBR: num_wann = 0 (auto) requires an explicit froz_emax = <eV> keyword` | `main.F:537-538` | Auto `num_wann` needs a numeric ceiling; `auto`/`none` are not enough. |
| `HYBR auto-nw: ranked candidate pool exhausted below the frozen capacity` | `main.F:564-566` | Even all candidate channels cannot supply `need + margin` columns. Widen the window or lower `froz_emax`. |
| `HYBR auto-nw exceeds the band window` | `main.F:567-568` | The derived `nw` is larger than `bmax - bmin + 1`. Widen `bmax` or lower `froz_emax`. |
| `HYBR: num_wann below frozen capacity` | `main.F:599-601` | The pre-gate violation of section 6.2. Raise `num_wann`, lower `froz_emax`, or set `capacity = warn` deliberately. |
| `HYBR: frozen manifold exceeds num_wann (capacity; capacity = warn overrides)` | `main.F:966-967` | The post-mask violation of section 6.3, including the `froz_emax = auto` case. |
| `HYBR: num_wann exceeds the band window` | `main.F:891` | `nw > Nb`. |
| `HYBR: num_wann exceeds the projection count (SCDM needs nproj >= num_wann...)` | `main.F:892-893` | Fewer `amn` columns than Wannier functions. Lower `threshold`, or raise `max_proj`. |
| `HYBR: no amn in memory` / `HYBR: no full-cell amn in memory` | `main.F:889`, `:912` | Internal; Stage 2 ran without `l2Amn` having exported its payload. |
| `HYBR: nnkp dims mismatch` | `main.F:934` | The `.nnkp` re-read disagrees with the first pass. The file changed under the run. |
| `inconsistent numbers of k-points between nnkp and energy files` | `main.F:233-234` | The `.win` mesh and the `lapw1` mesh differ. Regenerate one of them. |
| `inconsistent numbers of k-points between vector and energy files` | `main.F:348-349` | `case.vector` and `case.energy` are from different runs. |
| `build_kinfo: nfroz > nw (infeasible)` | `hybrid_core.F:154` | More frozen states than Wannier functions at some k, after the repair. Should be unreachable under `capacity = strict`; expected if you forced `capacity = warn` past a violation. |
| `build_kinfo: nfroz+ndis < nw (infeasible)` | `hybrid_core.F:155` | The admitted pool at some k is smaller than `nw`. Raise `pool_factor`. |
| `propagate_gauge: k-graph not connected` | `hybrid_core.F:751` | The `nnkp` neighbor graph does not connect the whole mesh. Regenerate the `.nnkp`; check `mp_grid`. |
| `mv_weights: dgelsd failed` | `hybrid_prep.F:145` | The b-vector shell system is singular. Regenerate the `.nnkp`. |
| `lowdin P: zheev failed` | `hybrid_prep.F:50` | The projection overlap matrix `S` is not diagonalisable, usually duplicated `amn` columns. |
| `nnkp: missing <tag>` | `hybrid_files.F:194` | Truncated `.nnkp`. |
| `unknown option '--...'` / `missing value for '--...'` | `w2w_hybrid.F` option parser | Standalone only: mistyped or incomplete named option; `--help` lists the valid set. |

---

## 9. Where to go next

* Full option surface, every `.inwf` keyword, every environment variable and
  every standalone option, with defaults, semantics and the options that
  are silently inert: see the companion reference document in this directory.
* Compile-time constants that are not run-time options: `AUTOW_GAP` (5 eV),
  `AUTOW_BELOW` (1 eV), `AUTOW_ECAP` (40 eV), `P2_ONSET` (5 eV), `P2_WIDTH`
  (2 eV) are `real(r8), parameter` declarations at `main.F:101-103`. They are
  **not** environment variables; changing them means editing `main.F` and
  recompiling.

