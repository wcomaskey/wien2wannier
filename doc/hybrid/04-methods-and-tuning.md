# HYBR: methods and tuning

> **v3.0 release note.** `w2w_hybrid` takes named `--key value` options
> (`w2w_hybrid --help` is the authoritative reference). The old positional
> interface, in which the options were positional arguments, was removed. Only
> the first two bare arguments remain, read as the seed and `num_wann`
> (`w2w_hybrid <case> <N> [options]`). Every further argument must be named, and
> a third bare argument stops the run with `unexpected positional argument`.
> Source line citations (`file.F:NN`) and paths under
> `tests/` in this guide refer to the development tree. The line numbers may
> be offset in the released sources, and the `tests/` cases are in the data
> record of the paper, not in this archive; every input deck quoted from them
> is reproduced in the text. The symmetry-adapted features of the development
> tree are not part of this release.


What each stage of the hybrid pre-Wannierization computes, and which parameter changes it.
This document assumes you already have a HYBR run that completes; it is about making the
result better.

All line references are to `SRC_w2w/` of the development tree. Where a number is a *default* it is
labeled as such; where it is a *validated production choice* (the value used by
`tests/regress_hybrid.sh:29`, the 9-cell regression harness) that is stated explicitly.

Two executables run the same Stage-2 code:

* `MODE=HYBR` inside `w2w`/`w2wc`, the in-memory one-shot (`main.F:852-855`, subroutine
  `hybrid_stage2` at `main.F:875-1034`), configured from the `.inwf` file.
* `w2w_hybrid` / `w2w_hybrid_acc`, the standalone, configured from named
  `--key value` options (`w2w_hybrid --help`).

They enter identical shared routines (`hybrid_prep_m`, `hybrid_core_m`, `hybrid_io_m`) but
they do **not** expose the same options. Where that matters it is flagged inline.

---

## 1. AUTO channel selection

The selection stage decides *which* `(atom, l, rf)` channels become `.amn` columns. It runs
before any Stage-2 work, in `main.F:361-514`.

### 1.1 Two different projectabilities

There are two distinct quantities in the pipeline; conflating them is the most common source
of confusion.

**(a) Channel projectability, for ranking.** `compute_projectability`
(`modules_rc.F:1193-1290`) calls `almgen` per (atom, k), then accumulates
`|A_lm|^2` summed over `m` and over equivalent atoms, with the `u_l`, `u_dot` and LO radial
contributions added (`modules_rc.F:1209-1226`), normalized at the end by
`Nk*(2l+1)*mult` (`modules_rc.F:1310-1317`). This is a muffin-tin quantity: it is not a
partition of unity and does not include interstitial charge
(`modules_rc.F:1319-1325`). Its only job is to rank channels.

Under SOC (def-file unit 9 present) the pass runs twice and the two spinor components are
*summed* before normalization, the `accumulate`/`finalize` flags at `main.F:424-472`.

**(b) Per-state projectability, for the masks.** Stage-2 does not use (a). It uses the
Löwdin form

```
p_i(k) = diag( A S^-1 A† ),    S = A† A
```

over the selected `.amn` columns, computed by `lowdin_pstate`
(`hybrid_prep.F:25-62`; eigen-decomposition of `S` with a `1e-10` relative null cut at
`hybrid_prep.F:53-58`). This is exactly bounded in `[0,1]`, which raw `|amn|^2` is not once
the projectors overlap. In HYBR the source is the **full-cell (tail-extended)** amn held in
memory as `projT_h`, restricted to the physical columns (`main.F:903-930`); in the standalone
it is `rescale=lowdin` (`--rescale`) on the amn named by `W2W_PAMN`
(`w2w_hybrid.F:178-189, 255-261`).

The run prints both the column choice and the resulting distribution:

```
 HYBR: P from 12 physical of 30 spatial columns x 2 sigma blocks
 HYBR: Lowdin P (full-cell): max =   1.0000   mean =   0.6667
```

(`main.F:918-929`). A `max` well below 1 means no state in the window is well represented by
the chosen columns; check the channel selection before anything else.

### 1.2 The gap-walk `bmin`

`bmin = 0` on line 2 of the `.inwf` triggers the gap walk (`main.F:283-304`): scan every
inter-band gap, and keep the **highest** gap that is wider than `AUTOW_GAP` **and** whose
upper band still starts below `E_F - AUTOW_BELOW`; cut there. Defaults are `AUTOW_GAP = 5.0`
eV and `AUTOW_BELOW = 1.0` eV (`main.F:101`). The second condition is what stops a
global-largest-gap rule from slicing into the valence complex (Sr 4s/4p is the case that
motivated it). No qualifying gap gives `bmin = 1`.

`bmax = 0` means "all bands present at every k" (`minNEn`, `main.F:262-266`); `bmax = -1`
means the highest band whose per-band minimum is within `AUTOW_ECAP = 40.0` eV of `E_F`
(`main.F:102, 267-275`). An explicit `bmax` above `minNEn` is clamped with a warning in every
mode (`main.F:312-318`).

`AUTOW_GAP`, `AUTOW_BELOW`, `AUTOW_ECAP` are `real(r8), parameter`, compile-time constants,
**not** environment variables. Changing them means editing `main.F` and rebuilding.

### 1.3 P2: the erfc energy weight (ranking only)

With line-4 field `p2_weight = T` (type default `.true.`, `util_w2w.F:1500`) each state's
contribution to a channel's ranking weight is multiplied by

```
w(E) = 0.5 * erfc( ((E - E_F) - P2_ONSET) / P2_WIDTH )
```

with `P2_ONSET = 5.0` eV and `P2_WIDTH = 2.0` eV (`main.F:103`, applied at `main.F:375-385`).
This suppresses high-energy nearly-free-electron states that otherwise inflate high-`l`
channels. It touches **only** the ranking: `projectability.dat` and the written `.amn` are
unweighted (`main.F:377-379`).

### 1.4 P3 / the support gate, and the Otsu cliff

Line-4 field `support_ratio` (type default `0.5`, `util_w2w.F:1501`) has three regimes
(`modules_rc.F:1507-1581`):

* `support_ratio > 0`, fixed gate. A `u_l` channel (`rf == 1`) whose best per-state
  projectability is below `support_ratio * max` gets weight zero and can never be selected
  (`modules_rc.F:1565-1575`). Across Fe/SrVO3/GaAs/MgB2 the measured dead zone is: physical
  channels `>= 0.56`, tail channels `<= 0.45` (`modules_rc.F:1443-1446`).
* `support_ratio = 0`, no gate.
* `support_ratio < 0`, **algorithmic cliff**, the validated production choice (`-1`). Otsu's
  two-cluster split is run on the log of the sorted per-channel support spectrum
  (`otsu_cut`, `modules_rc.F:1383-1410`), over non-noise channels only (noise floor
  `1e-6 * max weight`, `modules_rc.F:1516-1524`). It **gates nothing out**, because the gauge
  needs breadth. Instead the upper cluster is flagged as the *mask* channel set
  (`modules_rc.F:1546-1560`), and every ranking weight is rewritten as
  `weight * min(support, 1)` (`modules_rc.F:1537-1539`).

The run prints the spectrum and the cut:

```
   support spectrum (nsup= 14), otsu ncut=  4
    sup:   0.9985  0.9421  0.8863  0.7104  0.3312 ...
   cliff (support-otsu):   4 mask channels of  14 selected; mask columns =   12 (recorded in projmap)
```

`mask_cols` is what HYBR later uses to pick the physical columns for the Löwdin `P`
(`main.F:909-925`) and what the standalone reproduces via `--np-mask`.

**Cliff mode needs both halves.** `threshold <= 0` gives eligibility-by-noise-floor
(`modules_rc.F:1634-1636`) and records `mask_cols` (`modules_rc.F:1715-1728`);
`support_ratio < 0` is what sets `mask_flag` at all (`modules_rc.F:1507-1564`). Setting only
one gives `mask_cols = 0` and HYBR silently falls back to "P from all amn columns"
(`main.F:923-925`). The validated line is `0 0 T T -1`.

### 1.5 LO-shell orthogonalized projectors

Semicore shells carried by a local orbital are invisible to a `u_l`-based projector. For each
`(atom, l)` with an LO radial that is *not* parallel to `u_l`, a second candidate channel is
created using the perpendicular component

```
<phi_perp | psi>  =  (A_2 p2 + A_r q2) / sqrt(q2)
```

(`modules_rc.F:1236-1246`, channel creation at `modules_rc.F:1488-1500`). This is
parameter-free: an APW+lo radial parallel to `u` has vanishing perpendicular support and
sinks below the noise floor on its own. The channel carries `rf > 1` in `projmap.dat`. The
run prints the map:

```
 LO-shell atom  1 rf:  0  0  3  0  0  0  maxP:   0.000   0.000   0.913 ...
```

(`main.F:478-483`; the write statement itself is `main.F:480-481`). Do not confuse it with the
`projmap.dat` writer at `main.F:493-509`, whose `# mask_cols` header sits at `main.F:498-501`.

### 1.6 Selection defaults, in one place

| Field | Line | Default | Validated production value |
|---|---|---|---|
| `threshold` | 4.1 | `0.10` (`util_w2w.F:1497`) | `0` (cliff mode) |
| `max_proj` | 4.2 | `0` = auto (`util_w2w.F:1498`) | `0` |
| `write_proj` | 4.3 | `.false.` (`util_w2w.F:1499`) | `T` |
| `p2_weight` | 4.4 | `.true.` (`util_w2w.F:1500`) | `T` |
| `support_ratio` | 4.5 | `0.5` (`util_w2w.F:1501`) | `-1` (Otsu cliff) |
| `LJmax` | 3 | none (mandatory) | `3` in every validated test file |

`max_proj` counts **columns**, `sum (2l+1)*mult`, not channels
(`modules_rc.F:1604-1613`), and when it is `>0` it bypasses the threshold/cliff logic and the
feasibility guard entirely.

---

## 2. The capacity gate

### 2.1 What it computes

`need = max_k count( E_nk - E_F <= froz_emax )` over the band window (`main.F:541-549`). This
is the largest number of states that must be frozen at any single k. The gate runs **before**
`l2Amn`/`l2Mmn` (`main.F:531-620`), because adding channels afterwards would force a full amn
recompute.

`margin` defaults to `0`, which resolves to `max(4, ceiling(0.04*need))`
(`main.F:550-551`).

### 2.2 Why under-capacity is structurally dead

The disentangler fixes the frozen states as unit columns of `V_k` and optimizes only the
remaining `ntake = nw - nfroz` columns (`hybrid_core.F:320-328`). If `nfroz > nw` at any k
there is no subspace: `build_kinfo` hard-stops with `nfroz > nw (infeasible)`
(`hybrid_core.F:154`). To avoid that stop, `select_pool_and_frozen` performs a per-k
*feasibility repair* that silently demotes the excess frozen states
(`hybrid_prep.F:390-401`). The masks handed to the disentangler therefore *always* satisfy
`count <= nw`, and can never report the violation themselves; this is exactly the Rh3Si4
silent-failure mechanism. Hence `build_masks` returns `froz_demand_out`, the **pre-repair**
worst-k demand (`hybrid_prep.F:617-637`), and HYBR asserts on that
(`main.F:961-973`).

Concretely: asking for `nw` Wannier functions when 26 states sit below the ceiling at some k
does not produce a slightly worse model. It produces a model in which 14 of those states were
dropped from the frozen manifold by an energy-ordered rule you did not choose.

### 2.3 Reading the output

Explicit `num_wann`, healthy:

```
 HYBR capacity gate: nw=18 >= need+margin=9+4 at froz_emax=  2.000 eV  OK
```

Explicit `num_wann`, tight (warning only, the margin is a buffer, not a floor):

```
 HYBR capacity gate: nw=9 >= need=6 but inside margin 4 at froz_emax=  2.000 eV  (tight -- consider 10)
```

Explicit `num_wann`, violation (`main.F:592-606`):

```
 HYBR capacity gate: nw=12 < need=26 at froz_emax=  2.000 eV -> raise num_wann to >= 30 or lower froz_emax, or set capacity = warn
```

`num_wann = 0` (auto), which walks the ranked channel list adding whole spatial channels
(`sf*(2l+1)*mult`, `sf = 2` under SOC) until `nw >= need + margin` (`main.F:560-573`):

```
 HYBR capacity gate: need=200 (froz_emax=  2.000 eV)  margin=8  auto num_wann=208  (12 of 26 ranked channels)
```

Auto-`nw` **requires** an explicit numeric `froz_emax` and croaks otherwise
(`main.F:537-538`). It also croaks if the ranked pool is exhausted below capacity, or if the
resulting `nw` exceeds the band window (`main.F:565-570`).

The post-mask assertion is separate and always runs (`main.F:961-973`):

```
 HYBR: FROZEN CAPACITY VIOLATION: worst-k demand 26 > num_wann 12 -> raise num_wann to >= 26 or lower froz_emax
```

### 2.4 strict vs warn

`capacity = strict` (default, `util_w2w.F:1535`) croaks at both points. `capacity = warn`
downgrades both to printed warnings (`main.F:598, 965-972`) and proceeds; the frozen budget is
then spent on the `nw` **deepest** sub-ceiling states per k (`hybrid_prep.F:596-610`) and the
rest are released to the disentangler. This changes no physics elsewhere, but band fidelity near
`froz_emax` is not guaranteed, so check `eta`. Both outcomes are on record: a guarded
under-capacity Bi2D nw=12 model was both best-localized and most faithful, while an
under-capacity Rh3Si4 model failed at 434 meV (`main.F:589-596`).

**Pitfall.** The whole pre-`l2Amn` gate block is skipped when `froz_emax` is `auto` or `none`
(`main.F:531-539`), because the ceiling only exists after the walk-down. With
`froz_emax = auto` and an explicit `num_wann`, `margin` therefore affects nothing, and the
only capacity check you get is the post-mask assertion.

---

## 3. Band masks and the pool

`build_masks` (`hybrid_prep.F:497-637`) turns `P(nb,nk)` and the window eigenvalues into three
things: a contiguous band **pool**, an **admit** mask, and a **frozen** mask.

### 3.1 The pool

`select_pool` (`hybrid_prep.F:236-265`) takes the `nw` most-projectable in-window bands by
k-averaged `P`, then returns the contiguous band range covering them, widened to

```
n_pool_max = max(nw+1, ceiling(pool_factor * nw))
```

`pool_factor` defaults to `1.5` (`util_w2w.F:1530`, `hybrid_prep.F:508-509`); it is the
`.inwf` keyword `pool_factor` and the standalone's `--pool-factor` (which the production recipe does
**not** pass, so the harness also runs at 1.5).

`n_pool_max` is a **floor, not a ceiling**. `select_pool` returns
`npool = max(n_pool_max, s1-s0+1)` where `s0`/`s1` are the extreme band indices of the `nw`
most-projectable in-window bands, then clamps to the window with `npool = min(npool, nb-s0+1)`
(`hybrid_prep.F:248, 263, 265`). A spread-out P-based seed set therefore widens the pool past
`ceiling(pool_factor*nw)` on its own, and the both-end guard of §3.2 widens it further: the
regression logs show `bi2d` nw=12 at `pool[1..34] npool=34` compared to `n_pool_max = 18`
(`tests/regress_out/base/bi2d.log:13`) and `mgb2` nw=12 at `npool=22`
(`tests/regress_out/base/mgb2.log:12`). The practical point stands: pool-ratio
sweeps must set the ratio explicitly (standalone `--pool-factor`, or the `pool_factor` keyword).

The measured Fe **SOC** `nw=18` pool-ratio behavior (windows 9-26 / 9-30 / 9-35; not the
cells of `tests/HYBRID_MASK_RESULTS.md`): an isolated `Nb = Nw` model gives `eta` 67 meV, 1.99 Å²/WF and
`sigma_min` 2.4e-16 with no W90 convergence; 1.22x gives 13.70 meV / 1.19; 1.5x gives 6.79 meV
/ 0.798 with `Omega_I` 13.59. `tests/HYBRID_MASK_RESULTS.md` records the *separate* pf sweep
and the base/lowTLadm pair at pf 1.5 (0.8010 / 5.69 meV and 0.7925 / 6.91 meV,
`tests/HYBRID_MASK_RESULTS.md:109-124`); those are different cells, so do not expect the two
sets of numbers to match. Unprotected
over-width is a failure. A pool_factor sweep on Fe gave `eta` 5.69 / 15.39 / 487
meV at 1.5 / 2.0 / 3.0 with the MT `P`. With the Löwdin full-cell admit mask, 2x and 3x become
fail-soft (10.36 / 10.29 meV) because the mask excludes the widened junk.

### 3.2 The both-end pool guard

Any window state below the frozen ceiling that is **outside** the pool is structurally
unrepresentable. The guard (`hybrid_prep.F:535-586`) walks down from `pool(1)` and up from
`pool(npool)` and force-includes every band that has a sub-ceiling energy at some k, marking
the added rows `admit = .true.` and `frozen = (E - E_F <= fe_used)`. It runs whenever
`froz_emax /= none`.

Both directions were live defects: below `pool(1)` was the MgB2/Rh3Si4 deep-band exclusion
(MgB2 `Omega_I` 11.50 -> 8.59 once fixed), above `pool(npool)` was the pool-factor cap hiding
ceiling states from the disentangler (Bi2D pool top 18 vs true need 26, `eta ~3 eV`).

### 3.3 The coherent deepest-first frozen budget

After the guard, for each k with more sub-ceiling states than `nw`
(`hybrid_prep.F:588-610`):

```fortran
frozen(:,kk)    = .false.
frozen(1:nw,kk) = .true.      ! pool rows are energy-ordered
```

The `nw` **deepest** sub-ceiling states are frozen, giving a deterministic, contiguous
manifold. The alternative (capping the P-based flags row-wise) produces a patchwork that
differs from k to k and gave an incoherent subspace (`eta ~3 eV` vs 3.34 meV on the Bi2D
benchmark).

### 3.4 Regimes, and what HYBR overrides

`regime` (line 5 field 3, default `auto`) picks between per-state masks and band-level
freezing. `auto` decides from the degenerate-pair fraction, going band-level iff the fraction
exceeds 0.7 (`hybrid_prep.F:170-188, 522`); `band-level` additionally symmetrizes `P` inside
0.02 eV degenerate groups and freezes whole bands (`apply_band_level`,
`hybrid_prep.F:454-495`). The validated production value is `per-state`.

`froz_emax = auto` starts at 8 eV and walks down in 1 eV steps until the frozen demand (max over k
of pool states at or below the ceiling) is at most `num_wann`, with a floor of 1 eV (`hybrid_prep.F`, `apply_froz_emax`). `none` returns immediately and
also disables the pool guard and the energy-based demand count
(`hybrid_prep.F:426-429, 542, 625-628`). The validated production value is a literal `2.0`.

**Pitfall.** HYBR overwrites the computed admission mask unconditionally:

```fortran
admit_h = .true.
```

(`main.F:949`). The validated recipe is admit-all, which matches the standalone's `admitall`
(`--admit-all`, the standalone's default since 2026-10-06). The per-state admission that `build_masks` computed is
discarded, and there is no `.inwf` keyword to re-enable it. The standalone sets a floor with `--p-floor <x>`. The standalone's `nofroz`/`allfree` modes have no HYBR
equivalent either.

Typical mask report:

```
 HYBR: regime=per-state  pool[1..27] npool=27  frozen/k mean=  5.86  demand(worst k)=6  (capacity: nw >= demand)
```

### 3.5 The admission floor

The standalone admits every pool state unless `--p-floor <x>` is given (since 2026-10-06; the floor was
0.10 before), and HYBR always does (§3.4). The LCAO frontend (`lcao2wannier`) keeps a floor of 0.10. The
measurements behind the difference, on the eight LAPW cells of the paper at the release protocol
with only the admission changed:

- **A** `--admit-all` (the published models; the rerun reproduces their hand-offs byte for byte);
- **Z** `--p-floor 5e-7`: only the states of exactly zero projectability excluded (alpha-Sn and the Bi bilayer,
  which have none, `3e-4` and `1.6e-3`, at the largest gap of their P distributions);
- **F** `--p-floor 0.10`, the LCAO floor.

Wannier90 runs the published decks with `num_iter` raised to 300,000. Excluded: non-frozen pool states the
admission mask leaves out. Weakest link: the hand-off's subspace-continuity report. Band error: index-paired
window RMS over the two windows of the paper's LAPW table (MgB2: [-13,+7.2] and [-2,+2] eV). Scripts in
the development tree (not in this archive): `tests/regress_out/floor_gap_chain.sh`, `floor_gap_w90.sh`, `bi_floor_chain.sh` (the Bi bilayer with
its published guided deck); scores from `tools/band_compare.py` with the development tree's specs `floor_gap.json`
and `floor_gap_bi2d.json`.

| Cell | Setting | Excluded | Omega_I (A^2) | Weakest link | W90 | Band error (meV) |
|---|---|---|---|---|---|---|
| GaAs | A | 0 | 38.979 | 5.0e-4 | 59,500, met | 4.9 / 18.0 |
|  | Z | 3 | 39.344 | 7.3e-4 | 194,000, met | 8.2 / 22.8 |
|  | F | 26 | 39.623 | 1.6e-3 | 300,000 cap | 7.1 / 21.0 |
| MgB2 | A | 0 | 12.603 | 0.79 | 12,000, met | 11.5 / 4.4 |
|  | Z | 59 | 12.818 | 0.73 | 3,000, met | 5.7 / 2.3 |
|  | F | 2,124 | 14.088 | 0.63 | 6,000, met | 4.8 / 1.5 |
| alpha-Sn | A | 0 | 55.943 | 0.50 | 30,500, met | 2.2 / 14.2 |
|  | Z | 2 | 55.943 | 0.50 | 40,500, met | 2.6 / 13.9 |
|  | F | 557 | 63.702 | 0.24 | 1,500, met | 4.2 / 17.5 |
| Bi bilayer | A | 0 | 80.192 | 0.61 | 58,500, met | 0.2 / 17.6 |
|  | Z | 4 | 80.193 | 0.61 | 56,000, met | 0.2 / 17.6 |
|  | F | 7,292 | 365.03 | 1.2e-3 | 300,000 cap | 2.6 / 457.4 |
| Fe up | A | 0 | 5.876 | 0.65 | met | 8.0 / 15.1 |
|  | Z | 18 | 5.876 | 0.65 | 500, met | 8.2 / 14.9 |
|  | F | 1,666 | 7.086 | 0.64 | 500, met | 14.5 / 19.0 |
| Fe down | A | 0 | 6.002 | 0.61 | 500, met | 11.8 / 7.0 |
|  | Z | 18 | 6.002 | 0.61 | met | 11.8 / 7.0 |
|  | F | 1,420 | 6.941 | 0.50 | 500, met | 11.3 / 4.6 |
| Fe-SOC | A | 0 | 11.224 | 0.74 | 64,500, met | 4.5 / 5.5 |
|  | Z | 118 | 12.015 | 1.9e-5 | 300,000 cap | 54.1 / 14.4 |
|  | F (b) | 3,242 | 17.380 | 3.2e-5 | 300,000 cap | 125.0 / 38.2 |
| SrVO3 | A | 0 | 21.964 | 0.86 | 45,500, met | 0.7 / 1.6 |
|  | Z | 6 | 21.967 | 0.86 | 41,000, met | 0.7 / 1.6 |
|  | F (a) | 3,681 | 22.023 | 0.90 | met | 0.5 / 1.7 |

(a) At the walked ceiling (+6.0 eV) the subspace tore (weakest link 4.8e-19); the pocket skip lowered the
ceiling to +4.62 eV, where the row was measured. (b) The pocket skip tried +6.90 eV, did not repair the
tear, and kept +8.0 eV. "cap": stopped at the iteration limit without meeting the criterion.

What the table shows:

- The floor improves MgB2 (11.5 / 4.4 to 4.8 / 1.5 meV, at an Omega_I 1.5 A^2 higher) and Fe down
  slightly. It tears the Bi bilayer and Fe-SOC, whose wannierizations then do not converge in 300,000 iterations,
  with band errors 7 to 28 times larger, and it worsens alpha-Sn, Fe up and GaAs.
- Excluding only the zero-projectability states improves MgB2, leaves alpha-Sn nearly unchanged (2.6 / 13.9 vs. 2.2 / 14.2 meV), leaves Fe, SrVO3 and the Bi
  bilayer unchanged, worsens GaAs, and tears Fe-SOC: its zero lies at H on a band (+9.86 eV) whose states at the
  neighboring k-points have P ~ 1e-3 and stay admitted, so the subspace loses the band at one k-point.
- Neither Omega_I nor the continuity report identifies the best setting: on MgB2 setting A has the lowest Omega_I and
  a continuous subspace and interpolates worst.
- Setting A is the only one that keeps every tested subspace continuous apart from GaAs, which tears in all
  three settings. That is why it is the default.

**When to set a floor.** Only after a band comparison on the cell itself (§6.3): run the default, then
`--p-floor 0.10`, and keep the floor only if the continuity report stays clean and the band error falls,
as on MgB2. On the LCAO frontend the same scan (h-BN, MgB2, SnTe at floors 0 to 0.3) leaves 0.10 within
0.2 meV of the best floor in every window except two conduction windows; the full tables are in the
paper's Supplemental Material (section "The admission floor") and in `lcao2wannier`'s
`docs/DISENTANGLER_VS_WANNIER90.md`.

---

## 4. Ω_I disentanglement

### 4.1 What is minimized

`disentangle_frozen` (`hybrid_core.F:333-662`) minimizes the gauge-invariant spread

```
Omega_I = (1/Nk) sum_{k,b} w_b ( nw - || V_k† M^(k,b) V_{k+b} ||_F^2 )
```

(`omega_k`, `hybrid_core.F:227-246`) over per-k subspaces `V_k` subject to the frozen
constraint. `w_b` are Marzari-Vanderbilt finite-difference shell weights from the b-vectors in
the `.nnkp` (`mv_weights`, `hybrid_prep.F:96-167`); the printed residual

```
 HYBR: MV b-weight residual =   2.22E-16
```

is `max |sum_b w_b b_a b_c - delta_ac|`. Observed values across the regression cells span
`2.2e-16` to `3.4e-8` (the larger ones come from least-squares shell solves on lower-symmetry
cells). A residual materially above that means the b-shells are not a valid finite-difference
set; fix the mesh, not the disentangler.

### 4.2 Sweep structure

Seeding: per-k SVD of the amn restricted to the disentanglement rows
(`hybrid_core.F:394-407`), optionally replaced by a BFS parallel-transport seed
(`init_transport`, `hybrid_core.F:410-444`; standalone `--init-transport` only).

Each sweep is Gauss-Seidel over k. At each k (`sweep_one_k`, `hybrid_core.F:268-330`):
build `Zd = sum_b w_b (M_dis^(b) V_{k+b})(...)†` on the free sub-block, take its top `ntake`
eigenvectors as the free columns, and fix the frozen rows as unit columns. Neighbors are read
at their current value, so the ordering matters (Jacobi diverges; see the coloring note at
`hybrid_core.F:186-190`).

Sweeps run in blocks. The first block is one sweep and later blocks are `check_every` sweeps,
so exact checks land at sweeps 1, 1+ce, 1+2ce, … (`hybrid_core.F:501-509, 623`).

### 4.3 The convergence test

At each check, `Omega_I` is evaluated exactly and compared to the previous check
(`hybrid_core.F:610-616`):

```fortran
dlast = abs(omp - om)
if (dlast < tol .or. dlast < rel_tol*max(abs(om),1.0)) converged = .true.
```

Telemetry is unconditional (`hybrid_core.F:629-635`):

```
  disentangle sweeps=  2962  converged=T  |dOmega|_last=  5.23E-09  cap=  4000
```

### 4.4 The tunables, and where they are reachable

| Parameter | Standalone option | Value in the recipe | HYBR | `.inwf` route |
|---|---|---|---|---|
| `niter` (sweep cap) | `--sweeps` | **4000** | default `20000` | `sweeps` |
| `tol` (absolute) | — | `1e-9` (literal) | default `1e-9` | `tol` |
| `rel_tol` | `--rtol` | **1e-10** | default `1e-10` | `rel_tol` |
| `check_every` | `--check-every` | **10** | `10` (literal) | none |
| `colored` | `--colored` | on | default `false` | `colored` |
| `adaptive` | `--adaptive` | on | default `true` | `adaptive` |
| `mixing` | `--relax` | **-1.9** | default `-2.5` | `relax` or line 5 field 2 |
| `anderson` | `--anderson` | **0** (off) | default `0` | `anderson` |
| `init_transport` | `--init-transport` | off | not passed -> `.false.` | none |
| `verbose` | env `W2W_VERBOSE` | on **only** for the exact string `1` (`w2w_hybrid.F:157-162`) | literal `.false.` (`main.F:1003`) | none |

In the one-shot, only `check_every` and `init_transport` cannot be set, and `W2W_VERBOSE` has no effect. HYBR passes the
literal `.false.` (`main.F:1003`), so an hour-long disentanglement prints nothing until it
ends. Even on the standalone the trace is enabled only by `W2W_VERBOSE=1`
(`dis_verb = (trim(vbenv) == '1')`, `w2w_hybrid.F:157-162`); `true`, `yes`, `on` and `2` all
leave it off.

`colored` is red-black graph coloring of the k-neighbor graph
(`build_coloring`, `hybrid_core.F:185-223`) so that same-color k can be swept in parallel.
It is an OpenMP path, and **neither** build enables OpenMP by default: `SRC_w2w/Makefile` has
no `-fopenmp` anywhere (`FFLAGS` is `FGEN + FOPT`, `Makefile:71`), and `Makefile.hybrid` adds
it only inside `ifdef OMP` (`Makefile.hybrid:38-40`), which the `acc` target sets. Without that,
`colored` merely re-orders the Gauss-Seidel sweep with no speed benefit, in the standalone,
where the production recipe passes `colored 1`, exactly as much as in HYBR. Color order is a *different* Gauss-Seidel
ordering, and on a non-convex problem it reaches a different fixed point: on Rh3Si4 without
relaxation it converged in 881 sweeps to `Omega_I = 840.98188` where plain serial order took
1071 sweeps and reached `840.05736` (a lower, i.e. better, `Omega_I`, though with 4x worse
anchor conditioning; `SRC_w2w/bench/results.csv:29, 31`). Do not enable colored ordering
without relaxation.

### 4.5 `mixing` in detail

`mixing` blends the Gauss-Seidel Gram matrix, not the subspace
(`hybrid_core.F:309-318`):

```
Zd_eff = mix * Zd_new + (1 - mix) * Zd_eff_prev
```

Three regimes (`hybrid_core.F:462-475`):

* `mixing = 1`, **off**. The type default and the fallback when the optional field is missing
  are `-2.5` (scheduled); every HYBR `.inwf` quoted in these guides writes `1.0` explicitly.
* `mixing > 1`, fixed over-relaxation (SOR-like). `mixing < 1` damps.
* `mixing < 0`, **scheduled relaxation**. `mixv` starts at `|mixing|` for fast basin escape
  and is tapered geometrically toward 1, `mixv <- 1 + 0.5*(mixv - 1)`, after each *failed*
  exact convergence check (`hybrid_core.F:620-622`).

The validated production value is **-1.9** (scheduled, starting at 1.9). On Rh3Si4 H120 it
converged in 142 sweeps / 4.8 s versus 1071 serial sweeps / 187 s at `mixing = 1`, and to a
marginally lower `Omega_I` (`SRC_w2w/bench/results.csv:31, 70`).

**Coupling with `adaptive`.** The taper fires per failed exact check, not per
iterate. With `adaptive = .false.` the check cadence is `check_every` (`hybrid_core.F:623`),
so checks land at sweeps 1, 11, 21, … and `alpha - 1` halves each time: `1.9 -> 1.45 -> 1.225
-> 1.1125 -> 1.056`, i.e. below 1.1 by roughly sweep 31. With `adaptive = .true.` the exact
checks are gated on the free-block overlap plateau with a `10*check_every` fallback
(`hybrid_core.F:601-607`), so alpha stays near 1.9 for the first ~100 sweeps, which is the
regime the Rh3Si4 win was measured in. HYBR's `adaptive` keyword defaults to `true` (it acts with
`solver = gs`); on the standalone pass `--adaptive`. Treat
`mixing = -1.9` and `adaptive = .true.` as one change.

### 4.6 When the cap-hit warning appears

```
  disentangle sweeps=  4000  converged=F  |dOmega|_last=  2.62E-04  cap=  4000
  [WARN] disentangle hit the sweep cap WITHOUT converging
```

The reported `Omega_I` is then not a converged number, and `|dOmega|_last` is only a *lower
bound* on the remaining error (for a linear iteration with rate rho the residual is
`dlast*rho/(1-rho)`). Everything downstream (`rotate_to_subspace`, the gauge, the hand-off) is built on that subspace, and the amn/mmn/eig produced are syntactically valid. The run does not stop at this point.

Response, in order:

1. Judge the magnitude. `|dOmega|_last` of `3.4e-9` compared to `tol = 1e-9` (base cell `fedn`) is
   negligible. `1.05e-4` or `3.9e-2` is not.
2. On the standalone: raise `--sweeps` and make sure `--adaptive` and
   `--relax -1.9` are both on. `--anderson` is a parameter-free alternative
   that cut plain colored sweeps 881 -> 232 on Rh3Si4 with no tuning; it does not compose with
   the relaxation schedule, so use one or the other.
3. In HYBR raise the `sweeps` keyword (default 20000), or re-run Stage-2 with `w2w_hybrid_acc` on
   the files w2w already wrote; the code path is the same.
4. Verify the result rather than trusting it: the `I` component of the initial-spread report
   should equal the disentangler's `Omega_I` (see §6.2), and the bands should be plotted.

Cap hits are on record for the one-shot at its former fixed 2000 cap on several cells, including Rh3Si4 at
`|dOmega| = 3.86e-2` after 2513 s
(`tests/Rh3Si4/w2w_hybrid/outputwf_sig:1475`), and for the production standalone at its 4000
cap on `feup` and `fedn` (`tests/regress_out/base/feup.log`, `fedn.log`). A cap hit is common.
It should be examined before the hand-off is used.

---

## 5. Gauge selection

After `rotate_to_subspace` diagonalizes the projected Hamiltonian
(`hybrid_core.F:666-706`) and `B = A† W` is formed (`main.F:1006-1009`), `choose_gauge`
(`hybrid_core.F:1743-2108`) picks the `nw x nw` gauge that goes into the hand-off `.amn`.

### 5.1 SCDM anchors via QRCP

`scdm_anchors` (`hybrid_core.F:1493-1593`) forms the k-averaged density
`rho = (1/Nk) sum_k B_k B_k†`, then runs column-pivoted QR to select `nw` anchor columns. Two
candidate pivot sets are tried (QRCP on `rho`, and QRCP on `B[Gamma]†`) and each is
evaluated in two gauge frames; the combination with the best worst-k `sigma_min` wins. The
result is reported as

```
  SCDM sigma_min(worst k) =   6.4909E-01
```

This number is the single most informative gauge diagnostic. Order of magnitude on the
validated set: SrVO3 `6.5e-1`, Bi2D `6.9e-1`, MgB2 `3.2e-1`, Fe-up `5.0e-3`, Fe-dn `1.9e-4`,
Fe_soc `1.5e-4`, aSn `1.3e-4`, GaAs `6.0e-16` (`tests/regress_out/base/*.log`). GaAs is
singular: its flat semicore manifold makes the anchor gauge rank-deficient.

Two alternative anchor stages exist. `bscdm`/`bsms<N>` run SCDM per spectrally isolated
manifold (`block_scdm`, `hybrid_core.F:1607-1656`; the 0.5 eV gap literal is passed at the call
site, `hybrid_core.F:1768`) and cure exactly the GaAs singularity (`1.4e-16 -> 1.1e-4` in the
recorded A/B). `lowdin`/`lds<N>` replace QRCP anchors with the leading `nw` physical columns
(`lowdin_gauge`, `hybrid_core.F:1686-1706`, dispatched at `hybrid_core.F:1801-1807`); it is the
standard W90 projection baseline and never beat SCDM in the recorded comparison, regressing
aSn.

`W2W_ANCHORS=shellcomplete` constrains the pivots to whole `(atom,l,rf)` shells so the anchor
set is a representation basis, required by WannSymm/irrep tools, at a cost in spread
(`hybrid_core.F:1774-1796`). **Experimental**: it falls back to unconstrained SCDM with a
warning if `projmap.dat` is unusable (`hybrid_core.F:1788-1791`).

### 5.2 The `auto<N>` switch

Default gauge is `auto8` (`util_w2w.F:1532`). The rule (`hybrid_core.F:2016-2040`):

* worst-k `sigma_min >= sig_gate` (hardcoded `2.0e-3`, `hybrid_core.F:2005`): keep the per-k
  SCDM gauge as-is;
* below the gate: the per-k anchors are corrupt somewhere, so either run the seed tournament
  (if the seed policy is `opt<N>`) or transport from the best-conditioned k;

then `N` global smoothing sweeps either way. `smooth_gauge` (`hybrid_core.F:1196-1244`) is
Gauss-Seidel polar relaxation `U_k <- polar( sum_b w_b M^(k,b) U_{k+b} )`, which distributes
loop obstruction evenly over the BZ rather than piling it on the tree-closing edges of a
transport spanning tree.

Suffix parsing is lenient: unreadable digits silently give 8 sweeps
(`hybrid_core.F:2013-2014, 2076-2077`), and the whole dispatcher has **no else branch**
(`hybrid_core.F:1992-2099`). An unrecognized gauge string is not an error, it ships the bare
per-k SCDM gauge. `gauge = scdm` legitimately means "anchors, no smoothing"; `gauge = auto8x`
is a typo that happens to mean `auto8`.

Other modes: `ptg` (pure transport from the seed), `tfl<N>` (keep SCDM at healthy k,
BFS-transport into sick k only, then smooth), `sms<N>` (smoothing only), `pts<N>` (transport +
smooth).

### 5.3 The seed tournament

`seedk` defaults to `opt16` (`util_w2w.F:1533`), overridable by `W2W_SEEDK`. `opt<N>`
(`hybrid_core.F:1875-1942`) drops k whose per-k `sigma_min` is below `W2W_SEEDGATE`
(default `2.0e-3`), then samples up to `N` candidates evenly across the sigma-sorted passing
pool with both endpoints included, plus Gamma if it passes. `N` is clamped to `[2,64]`,
default 12 when the suffix is absent. Alternatives: `best` (largest sigma), `median`, or an
explicit integer k.

Every candidate then gets its own transport + loss transform + spread, scored by the **true
MV spread** `gauge_spread` (`hybrid_core.F:1062-1186`), and the argmin wins. Fixed
sigma-percentile rules are unsafe: the highest-sigma seeds are atypical k that failed Fe-up,
while median failed Fe_soc catastrophically (`eta` 116.7 meV in the recorded sweep). The
tournament therefore selects by the objective.

`seedloss` (default `smooth`) selects the loss functional:

| Mode | What is scored |
|---|---|
| `smooth` | transport + `nsw` global smoothing sweeps, full MV spread |
| `blind` | raw frame, seam-blind masked spread |
| `local` | raw frame + 8 seam-local sweeps, full spread |
| `abel` | raw frame + one-shot abelian seam removal, full spread |
| `absm` | abelian fix **then** `nsw` smoothing, both before scoring (what is scored is what ships) |

For `blind`/`local`/`abel` the winner additionally receives `nsw` global sweeps if `nsw > 0`,
so the scored frame is not the shipped frame; `absm` is the self-consistent variant. The
regression harness fixes one cell (`fednab`) to `pts8` + `absm` to keep that path
covered.

**Pitfall.** Under the default `auto8`, whenever the worst-k `sigma_min` is at or above `2e-3`,
the candidate list is still computed and printed but never used, the per-k SCDM gauge is
kept and only smoothing runs. The line

```
  seed tournament: 16 candidates from 36 of 36 gate-passing k
  auto gauge: sigma_min  6.92E-01 >=  2.0E-03 -> per-k SCDM kept
```

means the seed policy did nothing this run. The **tournament** itself runs only for `pts<N>`
(`hybrid_core.F:2079-2092`) and for `auto<N>` when the worst-k `sigma_min` is below the gate
(`hybrid_core.F:2018-2031`), in both cases only under an `opt<N>` policy.

`ptg` also ignores the seed policy. `k0seed` is initialized to `1` (`hybrid_core.F:1856`) and
the `opt<N>` branch builds `cand_seeds` without ever reassigning it
(`hybrid_core.F:1865-1942`), while the `ptg` branch calls
`propagate_gauge(..., k0seed, ...)` directly with no tournament (`hybrid_core.F:1992-1994`).
So `gauge = ptg` under the **default** `seedk = opt16` silently transports from k=1, whatever
`W2W_SEEDGATE`/`W2W_SEEDLOSS` say. Only the non-`opt` policies (`best`, `median`, or an
explicit integer k) move the `ptg` seed.

### 5.4 The initial-spread report

Every run ends the gauge stage with (`hybrid_core.F:2101-2107`):

```
  initial spread Omega =        62.5601   (I/OD/D =      19.5589      0.3098     42.6913 )
```

`gauge_spread` (`hybrid_core.F:1251-1305`) computes the three Marzari-Vanderbilt pieces on the
*rotated* mmn with the chosen gauge applied:

* **I**; `(1/Nk) sum_kb w_b sum_n (1 - |M_nn|^2) - OD`, i.e. leakage out of the subspace.
  Gauge-invariant, so it should reproduce the disentangler's `Omega_I` to printing precision.
  Bi2D: `Omega_I = 6.34295913` from the disentangler, `I = 6.3430` here.
* **OD**, off-diagonal weight, the part a good gauge removes.
* **D**, the diagonal/centre term.

This is the *starting point* Wannier90 will see, not the final spread. Its use is comparative:
between gauge modes or seeds on the same subspace, a lower total is a better basin, and W90
usually preserves the ordering. A large OD with a small I means the subspace is right and the
gauge is bad (more smoothing, or a different seed); a large I means the subspace itself is
leaky (pool, masks, `nw`).

The final sanity line is the contractivity check (`main.F:1029-1032`, the `handoff_svmax` call
and its write; the standalone's equivalent is `w2w_hybrid.F:368-370`), backed by
`handoff_svmax` at `hybrid_io.F:16-41`:

```
 HYBR: hand-off max sigma(M) =     0.999358  OK (contractive)
```

Every rotated overlap block must satisfy `sigma_max <= 1 + 1e-6`, but the printed number is a
**spot check, not an exhaustive one**: `handoff_svmax` strides over k with
`stride = max(1, nk/64)` (`hybrid_io.F:27-28`), so it samples at most ~64 k-points (1 in 8 at
`nk = 512`) and takes the max over all `nntot` b-blocks at those k. The source comment says so
(`hybrid_io.F:14`, "over a sampled set of (k,b) blocks"). A `VIOLATION` is therefore conclusive. Stop and investigate, and do not run W90. An `OK` only means that no *sampled* block violated the bound.

---

## 6. Quality metrics

### 6.1 The three numbers

**Spread per Wannier function** (Å²/WF), from the Wannier90 `.wout` after wannierization,
total `Omega` divided by `num_wann`. This is the localization figure and the one that compares
to the literature. It is **not** comparable across different k-meshes: the MV finite-difference
spread underestimates on coarse grids, so spreads *rise* with mesh refinement. It is also not
comparable across different `num_wann`.

**Omega_I**, from the disentangler, the invariant part. Lower is a better subspace. Comparable
only within one material, one mesh, one `nw`.

**eta**, band-interpolation error: the Vitale band distance between the interpolated Wannier
bands and the DFT bands, with an erfc weight at `mu = E_F + 2` eV, `sigma = 0.1`
(`tests/HYBRID_MASK_RESULTS.md:6-7`). This is the fidelity number, and the only one of the
three that can catch a model that is beautifully localized and physically wrong.

The three do not move together. The Fe mask study has cells where `Omega_I` improves while
`eta` degrades, and vice versa (`tests/HYBRID_MASK_RESULTS.md:109-124`).

### 6.2 Reference points

Stage-2 output from the 9-cell production regression
(`tests/regress_out/base/*.log`, standalone at the validated settings: `4000 1e-10 10 1 1
-1.9`, `W2W_SEEDK=opt16`, `W2W_SEEDLOSS=smooth`):

| cell | nw | pool | Omega_I | sweeps | SCDM sigma_min | initial Omega (I/OD/D) |
|---|---|---|---|---|---|---|
| gaas | 18 | 1..20 | 36.4454 | 120 ✓ | 6.0e-16 | 106.51 (36.45 / 57.79 / 12.27) |
| mgb2 | 12 | 1..22 | 8.5937 | 120 ✓ | 3.2e-01 | 37.21 (8.59 / 3.29 / 25.33) |
| feup | 9 | 1..14 | 6.5144 | 4000 ✗ (2.6e-4) | 5.0e-03 | 18.78 (6.51 / 6.66 / 5.61) |
| fedn | 9 | 1..14 | 6.5760 | 4000 ✗ (3.4e-9) | 1.9e-04 | 14.77 (6.58 / 6.44 / 1.75) |
| fesoc | 18 | 1..27 | 13.1924 | 3285 ✓ | 1.5e-04 | 28.15 (13.19 / 13.34 / 1.61) |
| asn | 20 | 1..30 | 52.6598 | 2962 ✓ | 1.3e-04 | 146.95 (52.66 / 58.39 / 35.90) |
| bi2d | 12 | 1..34 | 6.3430 | 20 ✓ | 6.9e-01 | 22.46 (6.34 / 0.01 / 16.11) |
| srvo3 | 19 | 1..27 | 19.5589 | 120 ✓ | 6.5e-01 | 62.56 (19.56 / 0.31 / 42.69) |

Post-W90 figures on record (`tests/HYBRID_MASK_RESULTS.md`), all with the guard-corrected
pools:

| model | Å²/WF | eta |
|---|---|---|
| MgB2 nw=12 (pool 1..22) | 0.8184, converged 985 it | — |
| Bi2D nw=12 (pool 1..34) | 0.5286, converged 12 it | 3.34 meV (max 8.0) |
| Bi2D fully automatic (`0 -1`, `num_wann=0`, fe 2.0) -> nw=16 | 3.1522 | 1.60 meV |
| Fe_soc nw=18 lowTLadm (pf 1.5) | 0.7925 | 6.91 meV |
| Fe_soc nw=18 base MT (pf 1.5) | 0.8010 | 5.69 meV |
| aSn nw=20 (sms gauge) | 3.5026 | — |

GaAs `nw=18` at the `sms8` gauge gave 2.7003 Å²/WF in a separate gauge comparison, not in `tests/HYBRID_MASK_RESULTS.md`.

Caveats attached to those numbers in the record: the Bi2D nw=12 cell
sits below the `froz_emax = 2.0` demand of 26 (capacity-clean Bi is the auto-nw model), and
all coarse-mesh `eta` figures carry a mesh floor; Bi2D `eta` went 13.11 meV on 6x6 to 0.52
meV on 12x12 for the same recipe. Mesh-converge a fidelity claim before quoting it.

### 6.3 Band comparison: the k-distance rule

**Band plots and every per-segment smoothness metric use cumulative
Cartesian k-distance as the x-axis, never the k-point index.**

WIEN2k `klist_band` snaps path points to a fixed denominator. On Bi2D with `den 120`, the
Gamma-M leg alternated steps of 1/120 and 2/120 while K-Gamma was uniform. Plotted against the
index, every band slope becomes a zigzag proportional to that slope, fake "jitter" that
appears on the DFT and the Wannier bands alike. Line
styling, WS-distance interpolation and a 12x12 SCF restart were each investigated before
the k-spacing was checked.

First lines of any band script:

```python
kcart = kfrac @ B                       # B from the .nnkp "begin recip_lattice"
x = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(kcart, axis=0), axis=1))])
# ticks at x[label_index]
```

When a band plot looks noisy, check `np.diff` of the path points **first**. Non-uniform steps
are this bug. Smoothness metrics such as `|second difference|` are likewise invalid on
non-uniform sampling treated as uniform.

A second rule: plot the bands and inspect them before trusting any scalar.
An `eta` weighted fully below `mu` is manifold-blind and can score inert
flat semicore bands as an excellent model. Require the model to straddle `E_F` before quoting
`eta` for a near-`E_F` model.

---

## 7. Tuning playbook

### Spread per WF too large

*Diagnose first:* read the initial-spread `I/OD/D` split.

* **Large OD, small I.** Subspace is fine, gauge is not. Raise the smoothing count
  (`gauge = auto16`), or if the SCDM `sigma_min` is below `2e-3` let the tournament run and
  widen it (`seedk = opt32`). Try `seedloss = absm` so the scored frame is the shipped frame.
  If `sigma_min` is near machine zero (GaAs-like semicore), the mode that helps is
  `gauge = bscdm` (*raw* block-SCDM, no smoothing) which recovered GaAs from 5.2496 (stuck,
  no W90 convergence in 5k iterations) to 2.7018 Å²/WF. `bsms8` is block-SCDM **plus** 8
  smoothing sweeps, and once you are smoothing at all, plain SCDM gives the same or a better result: GaAs `bsms8`
  2.7003 is identical to `sms8`/`auto8` 2.7003, and MgB2 `bsms8` 1.9440 is worse than `sms8`
  1.8917. Block-SCDM's value is removing the *dependence* on
  smoothing, not beating it. The `bs` prefix selects `block_scdm` (`hybrid_core.F:1763-1772`);
  the `bsms` branch applies `smooth_gauge` afterwards (`hybrid_core.F:2069-2098`).
* **Large I.** The subspace is leaky. Raise `pool_factor` (1.5 -> 2.0) so the disentangler has
  more freedom, and check that `nw` is not starved: a model that is too small cannot represent
  the manifold at any gauge. Confirm the disentangler converged.
* **Neither, and W90 does not converge.** Suspect an isolated-model hand-off on an entangled
  manifold. Fe SOC `nw=18` with `Nb = Nw` gives 1.99 Å²/WF and no W90 convergence; the same
  material at 1.5x pool gives 0.798 (see §3.1).

### Disentanglement not converging

See §4.6. In short: standalone, raise `niter`, ensure `adaptive = 1` and `mixing = -1.9`
together; consider Anderson acceleration as a tuning-free alternative (standalone `--anderson`, the
history depth; the recorded runs used `m = 3`, and it does **not** compose with the relaxation
schedule: the recorded AA+schedule probe oscillated to the cap). In HYBR raise `sweeps` (default 20000; `adaptive` is on by default), or finish a cell that does not converge
with `w2w_hybrid_acc` on the same files. Always compare
`|dOmega|_last` to `tol = 1e-9` before deciding it matters.

### The capacity gate refuses

The message names the fix. There are three options:

1. **Raise `num_wann`** to the printed value. The correct option when the demand is physical.
2. **Lower `froz_emax`.** `2.0` is validated for Bi2D-class cells; `froz_emax = auto` walks
   down on its own, at the cost of disabling the pre-`l2Amn` gate.
3. **`capacity = warn`.** Explicit, never silent. The frozen budget goes to the `nw` deepest
   sub-ceiling states per k and the rest are released. Then *check eta*, this configuration
   has produced both the best model in a set and a 434 meV failure.

If the gate refuses with `num_wann = 0` and "ranked candidate pool exhausted", the selection
did not produce enough channels: widen `bmax`, or lower `froz_emax` so `need` shrinks.

### Interpolation error too large on general k, though the symmetry path looks good

This is almost always an **under-resolved DFT k-mesh**, not a subspace defect. The Wannier
model is built from the mmn overlaps on the SCF mesh; a coarse mesh sets a floor that no
amount of gauge or mask work can cross. Bi2D at 6x6x1 gave `eta` 13.11 meV; the identical
recipe on a 12x12x1 DFT restart gave 0.52 meV.

Procedure: rerun the DFT with a denser mesh, regenerate the nnkp/amn/mmn, rerun Stage-2, and
compare `eta` at both meshes before attributing anything to the method. Note that spreads
*rise* with mesh density, so a cross-mesh spread comparison is meaningless in the other
direction too.

The second cause, if the mesh is already converged: the frozen manifold does not cover the
states that carry the error. Wide-window Bi2D at `fe = 2.5` stayed at 14.8 meV because
P-based freezing left the upper 6p free. Raising `froz_emax` to cover them, or moving to
`regime = band-level` on a degenerate manifold, is the lever.

### Symmetry breaking in the final model

The gauge is chosen by a localization objective, which has no reason to respect the point
group; models routinely break their own symmetry at the few-meV level. If you need symmetry
(WannSymm, irrep analysis, symmetric downfolding):

* Set `W2W_ANCHORS=shellcomplete` so the SCDM pivots are closed over whole `(atom,l,rf)`
  shells and both sigma blocks (`hybrid_core.F:1774-1796`). The anchor set is then a
  representation basis. **This costs spread** (it spends WF budget on orbitals SCDM declined)
  and it is flagged experimental here: if `projmap.dat` is missing or unparseable it falls
  back to unconstrained SCDM and prints
  `W2W_ANCHORS=shellcomplete: no usable projmap.dat -- falling back to unconstrained SCDM`
  (`hybrid_core.F:1788-1791`). Check for that line.
* Use `regime = band-level` on manifolds with substantial degeneracy, so `P` is symmetrized
  within 0.02 eV groups and whole bands are frozen together
  (`hybrid_prep.F:454-495`), a per-state mask can break a degenerate multiplet by freezing
  part of it.
* Check the magnitude before acting: a symmetry residual of a few meV near `E_F` is normal
  for a localization-optimal gauge, and "gaps" at that scale in a symmetry analysis are
  typically noise, not physics.

---

## Appendix: minimal and full HYBR inputs

Minimal *validated* input, auto band window, auto `num_wann`, cliff selection, and the two
non-default settings the production recipe needs (`froz_emax = 2.0`, `regime = per-state`):

```
HYBR
  0  -1
  3
0 0 T T -1
0 1.0 per-state
froz_emax = 2.0
```

Full surface, every keyword written at its **validated production value** (this is *not* the
default set; see the table below):

```
HYBR
  0  -1
  3
0 0 T T -1
0 1.0 per-state
froz_emax = 2.0
pool_factor = 1.5
margin = 0
gauge = auto8
seedk = opt16
seedloss = smooth
capacity = strict
```

The **type defaults**, for contrast (`util_w2w.F:1497-1501, 1524-1535`):

| Setting | Type default | Value in the block above |
|---|---|---|
| `threshold` (line 4.1) | `0.10` | `0` (cliff mode) |
| `max_proj` (line 4.2) | `0` | `0` |
| `write_proj` (line 4.3) | `.false.` | `T` |
| `p2_weight` (line 4.4) | `.true.` | `T` |
| `support_ratio` (line 4.5) | `0.5` | `-1` (Otsu cliff) |
| `mix` (line 5.2) | `-2.5` | `1.0` |
| `regime` (line 5.3) | `auto` | `per-state` |
| `froz_emax` | `auto` (numeric-mode value `8.0`) | `2.0` |
| `pool_factor` | `1.5` | `1.5` |
| `margin` | `0` | `0` |
| `gauge` | `auto8` | `auto8` |
| `seedk` | `opt16` | `opt16` |
| `seedloss` | `smooth` | `smooth` |
| `capacity` | `strict` | `strict` |

A file with all defaults would omit the `froz_emax` line entirely and read `0 -2.5 auto` on
line 5. With `froz_emax = auto` the whole pre-`l2Amn` capacity gate is **skipped**
(§2.4), so on such a file the gate is off and `margin` does nothing.

The keyword tail stops at the **first blank line** (`util_w2w.F:1652-1653`). Keywords after a
blank separator are silently ignored. An unrecognized key is a hard error
(`util_w2w.F:1692-1693`); note in particular that `rescale = ...` is **not** a parser case
and will abort the run. HYBR
always uses `rescale = 'none'` on an already-bounded Löwdin `P` (`main.F:943`).

Equivalent standalone invocation (the validated production recipe;
`--pool-factor` is omitted, so it takes its 1.5 default):

```sh
W2W_PAMN=<seed>.amn.tl W2W_SEEDK=opt16 W2W_SEEDLOSS=smooth \
  w2w_hybrid_acc --seed <seed> --nw <nw> --out <out> --rescale lowdin \
    --regime per-state --froz-emax 2.0 --sweeps 4000 \
    --colored --adaptive --relax -1.9 --gauge auto8 \
    --np-mask <np_mask> --admit-all --solver gs
```
