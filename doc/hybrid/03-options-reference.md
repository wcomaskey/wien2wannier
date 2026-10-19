# HYBR options reference

> **Public command.** `wien2wannier --help` is the authoritative workflow
> reference for automatic integrated runs; see also
> `00-public-wien2wannier-cli.md`. This document records lower-level `.inwf`,
> environment, and standalone `w2w_hybrid` controls. Source line citations (`file.F:NN`) and paths under
> `tests/` in this guide refer to the development tree. The line numbers may
> be offset in the released sources, and the `tests/` cases are in the data
> record of the paper, not in this archive; every input deck quoted from them
> is reproduced in the text. The symmetry-adapted features of the development
> tree are not part of this release.


Complete input surface of the hybrid pre-Wannierization mode, extracted from
`SRC_w2w/` of the development tree (line numbers may be offset, see the note above). Every entry cites the
file and line that implements it. Where a value is a *validated production
choice* rather than a code default, that is said explicitly.

There are five input channels:

| Channel | Read by | Notes |
|---|---|---|
| Public command | `SRC/wien2wannier` | Validates and writes the integrated contract; dedicated options plus checked `--hybrid-option KEY=VALUE` and explicit `--bridge-arg` forwarding |
| The `.inwf` file | `inwf_read_unit`, `SRC_w2w/util_w2w.F:1544-1723` | MODE word, band window, `LJmax`, AUTO config line, HYBR line, `key = value` tail |
| The def file (units 9 and 19) | `SRC_w2w/main.F:116-155` | optional second spin component; this is how SOC is switched on |
| Environment variables | see [Environment variables](#4-environment-variables) | the main `W2W_*` variables (`w2w_hybrid --help` lists all of them); not all are honored on both paths |
| Command-line options | `SRC_w2w/w2w_hybrid.F:45-129` | standalone tool `w2w_hybrid` only; named `--key value` options (see `w2w_hybrid --help`) |

`AUTOW_GAP`, `AUTOW_BELOW`, `AUTOW_ECAP`, `P2_ONSET`, `P2_WIDTH` are **not**
environment variables. They are `real(r8), parameter` compile-time constants
(`SRC_w2w/main.F:101-103`); changing them requires editing `main.F` and
recompiling.

---

## 1. The `.inwf` file in HYBR mode

### 1.1 Line-by-line layout

All reads are list-directed unless noted. Lines 4 and 5 are read into a
256-character buffer and parsed with *progressive fallback* (see the pitfalls
in §1.7).

```
line 1   MODE                                                  util_w2w.F:1555-1586
line 2   bmin  bmax                                            util_w2w.F:1588-1601
line 3   LJmax                                                 util_w2w.F:1607
line 4   threshold  max_proj  write_proj  [p2_weight]  [support_ratio]
                                                               util_w2w.F:1608-1628
line 5   num_wann  [mix]  [regime]            (HYBR only)      util_w2w.F:1635-1647
line 6+  key = value                          (HYBR only)      util_w2w.F:1650-1696
```

In the manual modes (`AMN`/`MMN`/`BOTH`) line 3 is instead `LJmax Nproj`
(`util_w2w.F:1699`), and lines 4-6 do not exist. In `AMN` and `BOTH` the
`Nproj` projection blocks follow; in `MMN` they do **not**, the block loop is
guarded by `if (inwf%Amn)` (`util_w2w.F:1701-1721`) and `MMN` sets only
`inwf%Mmn = .true.` (`util_w2w.F:1560-1561`), so an `MMN` file is exactly
three lines and any projection cards written after them are never read.

A complete, working HYBR file is `tests/Bi2D12/fullauto.inwf`, the
minimal fully-automatic case:

```
HYBR
  0  -1
  3
0 0 T T -1
0 1.0 per-state
froz_emax = 2.0
```

The same file with the whole keyword tail written out. Note that this is the
**validated production recipe**, not the all-defaults file: every line is at
its code default *except* `froz_emax`, whose default is `auto` (8 eV
walk-down, `util_w2w.F:1528-1529`). Writing `froz_emax = auto` here would
skip the entire pre-`l2Amn` capacity gate and make `num_wann = 0` a hard
error (`main.F:531-538`).

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

### 1.2 MODE words (line 1)

The word is uppercased before matching (`util_w2w.F:1557`), so lower case
works. An unrecognized word croaks.

| MODE | Writes | Selection | Source |
|---|---|---|---|
| `AMN` | `.amn` | manual projection cards | `util_w2w.F:1558-1559` |
| `MMN` | `.mmn` | — | `util_w2w.F:1560-1561` |
| `BOTH` | `.amn` + `.mmn` | manual projection cards | `util_w2w.F:1562-1564` |
| `AUTO` | `.amn` + `.mmn` | automatic channel selection | `util_w2w.F:1565-1568` |
| `AAMN` | `.amn` only | automatic channel selection | `util_w2w.F:1569-1576` |
| `HYBR` | `.amn.bin` + `.mmn.bin` (text `.amn`/`.mmn` only with `text_io = true`) + `<case>-hyb.{amn,mmn,eig}` | automatic + in-memory Stage-2 | `util_w2w.F:1577-1583`; `main.F:852-855, 874-1034` |

`HYBR` is `AUTO` plus the in-memory Stage-2 pipeline (`hybrid_stage2`,
`main.F:874-1034`), which runs on the root rank only (`main.F:852`).

### 1.3 Band window (line 2)

| Field | Default | Accepted | Meaning | Source |
|---|---|---|---|---|
| `bmin` | none (mandatory) | integer `>= 0`; `0` = auto sentinel (AUTO/AAMN/HYBR); negative croaks | Lowest band index. `0` = gap-walk: cut above the *highest* inter-band gap `> AUTOW_GAP` (5 eV) whose upper band still starts below `E_F - AUTOW_BELOW` (1 eV). No qualifying gap gives `bmin = 1`. | `util_w2w.F:1588,1594-1595`; `main.F:283-304` |
| `bmax` | none (mandatory) | integer `>= -1`; `0` and `-1` are sentinels (AUTO/AAMN/HYBR); `< -1` croaks | Highest band index. `0` = `minNEn`, the number of bands present at *every* k. `-1` = highest band whose per-band minimum lies below `E_F + AUTOW_ECAP` (40 eV). | `util_w2w.F:1588-1601`; `main.F:262-281` |

Additional checks: `bmin > bmax` (both positive) croaks (`util_w2w.F:1596`);
an empty auto window croaks (`main.F:306`); and a literal `bmax` above
`minNEn` is clamped with a warning in **all** modes, not only AUTO
(`main.F:312-318`).

### 1.4 `LJmax` (line 3)

| Field | Default | Accepted | Meaning | Source |
|---|---|---|---|---|
| `LJmax` | none (mandatory); `3` in every validated test input | integer | Maximum `L` in the Bessel/Gaunt expansion of the plane-wave overlap. Used by `l2Mmn`/`radint`/`init_bessel` only; it does **not** affect AUTO channel selection or the `.amn`. | `util_w2w.F:1607`; `main.F:356` |

### 1.5 HYBR line (line 5)

Read with progressive fallback: 3 fields, else `regime = 'auto'` and 2 fields,
else `mix = -2.5` and 1 field (`util_w2w.F:1635-1645`).

| Field | Default | Accepted | Meaning | Source |
|---|---|---|---|---|
| `num_wann` | `0` in the type; mandatory in the file | integer `>= 0`; `0` = capacity-driven auto sentinel; negative croaks | Target Wannier-function count, i.e. the hand-off dimension. See [the capacity gate](#23-num_wann--0-and-the-capacity-gate). Post-`l2Amn` guards: `nw > Nb` croaks, `nw > nproj` croaks. | `util_w2w.F:1522,1646-1647`; `main.F:526-620, 888-893` |
| `mix` | `-2.5` (scheduled; `1` = off) | real. `1` off; `<1` damps; `>1` extrapolates; negative = *scheduled* extrapolation starting at `\|mix\|` | Zd mixing factor handed to `disentangle_frozen` as `mixing`. In scheduled mode the excess over 1 is halved at every convergence check. Standalone equivalent: `--relax`. | `util_w2w.F:1524,1642`; `main.F:1003`; `hybrid_core.F:354-356, 465-475` |
| `regime` | `'auto'` | `auto` \| `per-state` \| `band-level` | Mask regime passed to `build_masks`. `per-state`: no band-level freeze, no degeneracy symmetrization. `band-level`: symmetrize P within 0.02 eV degenerate groups and freeze whole bands. `auto` (or any unrecognized string) lets `auto_regime` decide from the degenerate-pair fraction (band-level iff fraction `> 0.7`). The resolved value is echoed. | `util_w2w.F:1525,1637`; `main.F:944`; `hybrid_prep.F:170-188, 523-527` |

`per-state` is the validated production choice; the code default is `auto`.

### 1.6 The keyword tail (line 6+)

HYBR only, entirely optional, order-free, one `key = value` per line. The key
is everything left of the first `=` (max 64 chars, `adjustl`'d), the value
everything to its right (max 192 chars).

| Keyword | Default | Accepted values | Meaning | Source |
|---|---|---|---|---|
| `froz_emax` | `auto` (with `hyb_froz_val = 8.0` eV) | `auto` \| `none` \| `<real eV>` | Frozen-manifold ceiling relative to `E_F`, in eV. `auto`: start at 8 eV and walk down in 1 eV steps until the frozen demand (max over k of pool states at or below the ceiling) is at most `num_wann`, floor 1 eV. The exit is the demand test alone, and the earlier requirement that at most 15% of k-points be saturated was removed (it never triggered on the eight release cells and misjudges an isolated target manifold, which can saturate every k at `num_wann` states). `none`: no ceiling; `apply_froz_emax` returns immediately, the pool guard is skipped, and the frozen demand falls back to the post-mask frozen count. `<eV>`: literal ceiling. A bad numeric croaks. | `util_w2w.F:1528-1529,1664-1672`; `main.F:944-945`; `hybrid_prep.F:417-451, 542, 625-628` |
| `pool_factor` | `1.5` | any parsable real; unparseable croaks. There is **no** range check: `< 1` is accepted silently and is equivalent to a pool floor of `nw+1` bands | Disentanglement pool **floor**, not a cap: `n_pool_max = max(nw+1, ceiling(pool_factor*nw))` is only one of the two terms in `npool = max(n_pool_max, s1-s0+1)`, where `s0..s1` is the index span of the `nw` most-projectable in-window bands. The pool widens past `pool_factor*nw` whenever that P-selected span is wider, and again if the `froz_emax` pool guard extends it; the only hard cap is the band count (`npool = min(npool, nb-s0+1)`). Standalone equivalent: `--pool-factor`. | `util_w2w.F:1530,1673-1675`; `main.F:947`; `hybrid_prep.F:248, 261-265, 509, 535-586` |
| `margin` | `0`, resolved to `max(4, ceiling(0.04*need))` | integer (bad value croaks); `<= 0` selects the auto rule | Headroom above the frozen demand used by the capacity gate. With `num_wann = 0` it sets the auto target (`need + margin`); with an explicit `num_wann` it only decides between the "OK" line and the "inside margin (tight)" warning. | `util_w2w.F:1531,1676-1678`; `main.F:549-550, 560-573, 607-615` |
| `gauge` | `auto8` | `scdm` \| `ptg` \| `auto[<N>]` \| `tfl<N>` \| `sms<N>` \| `pts<N>` \| `lowdin` \| `lds<N>` \| `bscdm` \| `bsms<N>` (`character(len=16)`) | Gauge dispatcher mode, passed to `choose_gauge`. See [gauge modes](#5-gauge-modes-in-detail). Standalone equivalent: `--gauge`. | `util_w2w.F:1532,1679-1680`; `main.F:1023-1026`; `hybrid_core.F:1743-1807, 1992-2099` |
| `seedk` | `opt16` | `opt[<N>]` \| `best` \| `median` \| `<integer k>` | Transport-seed policy. `opt<N>`: drop k whose per-k `sigma_min` is below `W2W_SEEDGATE`, sample up to N candidates evenly across the sigma-sorted passing pool (endpoints included) plus Gamma if it passes, then score each by true MV spread and keep the argmin. `best`: largest per-k `sigma_min`. `median`: median-sigma k. `<integer>`: explicit seed (out of range warns and falls back to `k = 1`). N is clamped to `[2,64]`, default 12 when the suffix is absent or unparseable. Overridden by `W2W_SEEDK`. | `util_w2w.F:1533,1681-1682`; `hybrid_core.F:1854-1857, 1865-1991` |
| `seedloss` | `smooth` | `smooth` \| `blind` \| `local` \| `abel` \| `absm` | Loss functional used to score tournament candidates. `smooth`: transport + `nsw` global smoothing sweeps, full MV spread. `blind`: raw frame, seam-blind masked spread. `local`: raw frame + 8 seam-local sweeps, full spread. `abel`: raw frame + one-shot abelian seam removal, full spread. `absm`: abelian fix then `nsw` smoothing, both before scoring. For `blind`/`local`/`abel` the winner additionally receives `nsw` global sweeps when `nsw > 0`. Overridden by `W2W_SEEDLOSS`. | `util_w2w.F:1534,1683-1684`; `hybrid_core.F:1072-1077, 1105-1129, 1167-1184` |
| `capacity` | `strict` | `strict` \| `warn` (anything else croaks) | What happens on a frozen-capacity violation. `strict` croaks at both the pre-`l2Amn` gate and the post-mask assertion. `warn` prints and proceeds; the per-k frozen budget is then spent on the `nw` deepest sub-ceiling states and the rest are released to the disentangler. It changes no physics, only whether the run stops. | `util_w2w.F:1535,1685-1691`; `main.F:592-616, 961-973`; `hybrid_prep.F:596-610` |
| `sweeps` | `20000` | positive integer; `< 1` or unparseable croaks | Ω_I disentanglement sweep cap. On reaching it the run does **not** fail: it recomputes Ω_I exactly, continues with the unconverged subspace, and grades the tail by comparing Ω_I with its value one tenth of the sweeps earlier: `flat tail (...); usable, tolerance not met` below a relative change of 1e-5, otherwise `POOR CONVERGENCE: still descending (...)`, which means the cap cut the minimization short (raise `sweeps`, or use `solver = cg`). The classification is printed after `[WARN] disentangle hit the sweep cap WITHOUT converging` and recorded in `<case>-hyb.dis` (§2.6). | `util_w2w.F` (`hyb_sweeps`); `main.F` disentangle call; `hybrid_core.F` (`disentangle_frozen`, `cap_tail_text`) |
| `solver` | `cg` | `cg` \| `gs` \| `auto` (anything else croaks) | Omega_I solver. `cg`: Riemannian conjugate gradients on the Grassmannian, parallel over k with no coloring; it stops when the projected gradient falls to max(1e-10, 1e-6 \|grad0\|) (exit `OK`; `STALL` and `CAP` are not certified minima) and does not read the convergence tolerances. `gs`: the Gauss-Seidel sweep (the default before 2026-10), which the `relax`, `adaptive`, `anderson` and `colored` keywords control. `auto`: `cg` when the largest free-subspace dimension `ndis` reaches 250, else `gs`. The two solvers can stop in different Omega_I minima of an entangled cell; see the continuity diagnostic (§2.6). | `util_w2w.F:1664-1669, 1901-1906`; `main.F:1798-1938`; `hybrid_grass.F` (`grassmann_cg`) |
| `relax` | `-2.5` (with `kick = true`) | any non-zero real; `0` croaks | Zd relaxation factor, and the most important performance parameter. `1` = off. `> 1` = fixed over-relaxation. **Negative = scheduled**: start at `|relax|` and taper geometrically toward 1 (`α ← 1 + (α−1)/2`) after each failed exact convergence check, aggressive early, monotone-safe late. Overrides the positional `mix` field on the `num_wann` line, which is why it is the preferred spelling. Measured on Rh₃Si₄ H120: 1071 sweeps at `relax = 1` vs 50 at `1.7`; fixed `1.9` *overshoots* (Ω_I 843.03 vs 840.06), which is why the shipped default is the scheduled form. The optimum is material-dependent. | `util_w2w.F` (`hyb_mix`); `hybrid_core.F:461-474, 619-622` |
| `adaptive` | `true` | `true`\|`t`\|`yes`\|`y`\|`on`\|`1` / `false`\|`f`\|`no`\|`n`\|`off`\|`0`, case-insensitive; anything else croaks | Defer the exact `omega_invariant` convergence check until the free eigen-sum plateaus, instead of testing every `check_every` sweeps. Pure overhead reduction. It changes when the exact check is made, not the minimum reached. | `util_w2w.F` (`hyb_adaptive`); `hybrid_core.F:601-609` |
| `anderson` | `0` | integer `>= 0`; negative or unparseable croaks | Safeguarded alternating Anderson acceleration: history depth (`0` = off). Extrapolates the free-block subspace *projectors* over the last `m+1` iterates and accepts only on improvement. Implemented and safeguarded but **off in the validated recipe**, so treat non-zero values as experimental. | `util_w2w.F` (`hyb_anderson`); `hybrid_core.F:476-477` |
| `colored` | `false` | as `adaptive` | Red-black k-parallel Gauss-Seidel sweep. **Off by default.** It is a different sweep ordering that converges to a *different, equally valid* fixed point (Rh₃Si₄ H120: serial 840.05736, colored 840.98188), so enabling it invalidates every existing md5 baseline, and `w2w`/`w2wc` are built without `-fopenmp`, so it would change the answer without adding speed. Only turn it on together with an OpenMP build and a BLAS built with `USE_LOCKING`. | `util_w2w.F` (`hyb_colored`); `hybrid_core.F:185-223, 509-528` |

> **History.** Before 2026-07-29 the keywords `sweeps`, `relax`, `adaptive`, `anderson` and `colored` did not exist: the one-shot
> hardcoded `2000` sweeps with `relax` effectively `1.0` (the old default) and
> both `adaptive` and `colored` off, so it ran a configuration that was never
> validated and hit the sweep cap on GaAs and on the Rh₃Si₄ SOC cell. The
> defaults above are the production recipe from `tests/regress_hybrid.sh`.
> Inputs written before that date that pass an explicit `1.0` in the positional
> `mix` field still get `1.0` (an explicit value overrides the default) and the run prints
> `[NOTE] relax=1 disables over-relaxation`. Add `relax = -2.5` to fix them.

> **The standalone `w2w_hybrid` does not default to this recipe.** Its defaults
> (`w2w_hybrid.F:64-80`) are `--rescale occ99`, `--gauge scdm`, `--relax 1` (no
> over-relaxation) and adaptive checking off, while the HYBR keyword that the
> `wien2wannier` driver writes defaults to the full-cell Löwdin projectability,
> `gauge = auto8`, `relax = -2.5` with `kick = true`, and `adaptive = true`
> (`util_w2w.F:1607-1691`). Both default to `solver = cg`, for which `relax` and
> `adaptive` do not apply. A standalone run reproduces a driver run only when
> these options are passed explicitly, as the published LAPW runs did.

The parser also accepts `edges`, `relax_late`, `taper`, `taper_late`, `kick`, `text_io`, `tol` and `rel_tol` (defaults in `wien2wannier --help`). Any other key is a hard error: `HYBR: unknown keyword "..."`
(`util_w2w.F:1692-1693`). A tail line without `=` croaks
(`util_w2w.F:1659-1660`).

### 1.7 Pitfalls in the `.inwf` parser

- **`rescale = ...` is not a keyword.** The parser
  has no such case (`util_w2w.F:1663-1694`), so it hits the unknown-keyword
  croak and the run aborts. HYBR always builds P as the Löwdin per-state
  projectability of the full-cell `amn` and calls `build_masks` with
  `rescale_mode = 'none'` (`main.F:927, 943`). The `rescale` choice exists
  only on the standalone (`--rescale`).
- **The tail stops at the first blank line** (`util_w2w.F:1653`). Keywords
  after a blank separator are silently ignored, no warning, no error.
- **Progressive fallback hides typos on lines 4 and 5.** The fallback is
  *staged*, not all-or-nothing. On line 4: a malformed field 5 fails the
  5-field read, silently resets `support_ratio` to `0.5`, and the 4-field
  retry then succeeds, so `p2_weight` keeps the value written in the file.
  `p2_weight` is forced to `.true.` only when the 4-field read also fails,
  i.e. when field 4 itself is bad (`util_w2w.F:1615-1628`: the outer
  `if (confios /= 0)` sets `support_ratio`, the inner one sets `p2_weight`).
  Same staging on line 5: a bad `regime` resets only `regime` to `'auto'`;
  `mix` reverts to `-2.5` only if field 2 is also bad (`util_w2w.F:1638-1645`).
- **`margin` does nothing unless `froz_emax` is numeric.** The entire
  capacity-gate block is skipped for `froz_emax = auto|none`
  (`main.F:531-539`); only the post-mask assertion at `main.F:961` remains.
- **HYBR discards the per-state admission mask.** `admit_h = .true.` is set
  unconditionally right after `build_masks` (`main.F:949`). There is no
  keyword to re-enable per-state admission, and no HYBR equivalent of the
  standalone's `nofroz`/`allfree` mask modes.
- **HYBR fixes `check_every = 10` and `verbose = .false.`** and never passes
  `init_transport`; `sweeps`, `tol`, `rel_tol`, `relax`, `adaptive`,
  `colored` and `anderson` are keywords (§1.6).

---

## 2. The AUTO projection-selection block

### 2.1 Line 4 fields

Validated production line: `0 0 T T -1`.

| Field | Default | Accepted | Meaning | Source |
|---|---|---|---|---|
| `threshold` | `0.10` (type default; the field itself is mandatory) | real. `>0` = noise-floor mode; `<=0` = **cliff mode** | `>0`: channels with `weight(i)/weight(1) < threshold` are ineligible (legacy relative-weight floor). `<=0`: every channel above a `1e-6` relative numerical-noise floor is eligible; the physical/gauge boundary comes from the support Otsu cut instead, selected channels are stably partitioned mask-first so the mask is a *leading* `amn` block, and `mask_cols` is recorded into `projmap.dat` and the in-process payload. | `util_w2w.F:1497,1610`; `modules_rc.F:1623-1636, 1657-1691, 1715-1728` |
| `max_proj` | `0` (= auto) | integer `>= 0` | When `>0`, channels are taken in ranked order while the running **column** total (`sum of (2l+1)*mult(atom)`) stays `<= max_proj`. This bypasses only what lives in the `else` branch: the eligibility selection (threshold noise floor / cliff drop) and the `Nb` feasibility guard. The support-Otsu mask detection and the cliff-mode ranking-weight rewrite run *before* the branch (keyed on `support_ratio < 0`, `modules_rc.F:1507-1564`) and the `mask_cols` recording runs *after* it (keyed on `threshold <= 0`, `modules_rc.F:1715-1729`) over the `nsel` this branch produced, so with `max_proj > 0`, `threshold <= 0` and `support_ratio < 0` all three still apply. It counts `amn` columns, not channels. | `util_w2w.F:1498,1611`; `modules_rc.F:1604-1613` |
| `write_proj` | `.false.` | Fortran logical `T`/`F` | Write `projectability.dat` (unit 97, per band/k/atom/l) and `autoprojection.diag` (unit 98: channel ranking, alternative sets, ready-to-paste `.inwf` templates). Also allocates `proj_state(Nb,Nk,Nneq,0:Lmax2)`. | `util_w2w.F:1499,1612`; `main.F:391-392,473-476`; `modules_rc.F:1322-1376, 1838-1841` |
| `p2_weight` | `.true.` (also the fallback when the 4-field read fails) | Fortran logical | P2 energy ranking weight: each state's contribution to the channel weight is multiplied by `0.5*erfc((E-E_F-P2_ONSET)/P2_WIDTH)` = `0.5*erfc((E-E_F-5 eV)/2 eV)`, suppressing high-energy nearly-free-electron states. Applied to the `(atom,l)` **ranking only**; the per-state output and the written `.amn` stay unweighted. | `util_w2w.F:1500,1613`; `main.F:375-385` |
| `support_ratio` | `0.5` (also the fallback when the 5-field read fails) | real. `>0` fixed P3 gate ratio; `0` disables; `<0` = algorithmic cliff gate | `>0`: zero the weight of any `u_l` channel (`rf == 1`) whose per-state maximum projectability is below `support_ratio * best-channel max`; those channels sink to the bottom of the ranking and are never eligible. `0`: no-op. `<0`: Otsu two-cluster split of the sorted support spectrum over non-noise channels. It gates nothing out; instead it flags the upper cluster as the *mask* channel set and rewrites every ranking weight to `weight*min(support,1)`. | `util_w2w.F:1501-1502,1614`; `modules_rc.F:1507-1581` |

**Cliff mode needs both halves.** `mask_flag` is only ever set in the
`support_ratio < 0` branch (`modules_rc.F:1507-1564`), while `mask_cols` is
only recorded in the `threshold <= 0` branch (`modules_rc.F:1715-1728`).
`threshold <= 0` alone yields `mask_cols = 0`, and HYBR then falls back to
"P from all `amn` columns" (`main.F:923-925`). `support_ratio < 0` alone still
rewrites the ranking weight by the support factor
(`modules_rc.F:1537-1539`) but records no mask. The validated line
`0 0 T T -1` sets both.

### 2.2 `bmax` sentinels

- `bmax = 0` -> `minNEn`, the number of bands present at every k
  (`main.F:262-266`). Safe upper bound; never reaches `read_vec`'s
  uninitialized-`A` path.
- `bmax = -1` -> the highest band `i` with `(ebmin(i)-E_F)*Ry <= AUTOW_ECAP`
  = `E_F + 40 eV` (`main.F:267-275`). The scan runs over `1..minNEn`, so the
  result is still bounded by the all-k band count.

### 2.3 `num_wann = 0` and the capacity gate

The gate runs *before* `l2Amn`/`l2Mmn` (`main.F:516-620`), because adding
channels later would force an `amn` recompute. It requires an explicit
numeric `froz_emax`; `num_wann = 0` with `froz_emax = auto|none` croaks
(`main.F:537-538`).

```
need = max over k of  count{ i in [bmin,bmax] : (E(i,k) - E_F)*Ry <= froz_emax }   main.F:541-548
marg = hyb_margin  if > 0,  else  max(4, ceiling(0.04*need))                       main.F:549-550
```

`need` is already in WF units; SOC energy files list spinor states.

With `num_wann = 0` the code walks the ranked channel list from the top,
adding `sf*(2l+1)*mult(atom)` columns per channel (`sf = 2` under the SOC
sigma loop), and stops at the first total `>= need + marg`; that total becomes
`num_wann` (`main.F:551-569`). It croaks if the ranked pool is exhausted below
capacity, and if the result exceeds the band window. It warns if the window
holds fewer than `pool_factor*nw` bands (`main.F:577-583`).

With an explicit `num_wann` the gate only classifies: `nw < need` is a
violation (croak, unless `capacity = warn`); `need <= nw < need+marg` prints
a "tight" warning; otherwise "OK" (`main.F:585-617`).

A second, independent assertion runs after the masks are built
(`main.F:961-973`) on the **pre-repair** demand returned by `build_masks`.
This one also covers `froz_emax = auto`, where the ceiling only exists after
the walk-down. It is needed because the per-k feasibility repair inside the
mask builder silently demotes frozen states above `nw`
(`hybrid_prep.F:596-610`), so the returned mask can never flag the violation
itself.

### 2.4 `AAMN`: AUTO selection, `amn` only

`AAMN` sets `Amn = .true.`, `Mmn = .false.`, `auto_mode = .true.`
(`util_w2w.F:1569-1576`). It re-runs channel selection and writes a new
`.amn` against an **existing** `.mmn` computed on the same band window,
skipping the overlap recompute entirely. The comment in the source records
the motivation: on Rh3Si4 the `.mmn` is 12 GB and takes hours, while the
`.amn` takes minutes, which is what makes widening the candidate set
affordable. The band window must be unchanged, since `.amn` and `.mmn` must
agree on `Nb`.

### 2.5 Files the AUTO/HYBR path writes

| File | Written when | Content | Source |
|---|---|---|---|
| `projmap.dat` | always in AUTO/AAMN/HYBR (root) | `amn` column -> `iatom(global) l m rf`; carries `# mask_cols N` in cliff mode | `main.F:493-509` |
| `projectability.dat` | `write_proj = T` | per (band, k, atom, l) raw projectability | `modules_rc.F:1322-1376` |
| `autoprojection.diag` | `write_proj = T` | channel ranking, gap analysis, alternative sets, `.inwf` templates | `modules_rc.F:1885-1920` |
| `<case>-hyb.{amn,mmn,eig}` | HYBR | the isolated `Nb = Nw` hand-off | `main.F:1015, 1027` |
| `<case>-hyb.anchors` | HYBR | one anchor candidate-column index per WF; joins with `projmap.dat` to place real guiding centres in the `.win` | `hybrid_core.F:1810-1817` |
| `<out>.masks` | HYBR (`<case>-hyb.masks`) and standalone `w2w_hybrid` | the per-state masks the hand-off was built from: one row per (k, pool state) with `k band eps-E_F P frozen admitted`, header lines carry `p_froz`, `p_gate`, `froz_emax` and the pool. `tools/plot_certificate.py` draws the window certificate from it | `hybrid_io.F` (`write_masks`) |
| `<out>.anatomy` | HYBR (`<case>-hyb.anatomy`) and standalone `w2w_hybrid` | the auto-trust anatomy: one row per pool band with `Pbar minP minE_rel_EF relevant cls` (`trusted`/`gray`/`junk`) and the derived `p_gate`, `p_froz`, `plateau` in the header. `tools/plot_anatomy.py` draws it | `hybrid_io.F` (`write_anatomy`) |
| `<out>.dis` | HYBR (`<case>-hyb.dis`) and standalone `w2w_hybrid` | two comment lines: the disentanglement outcome (`converged (gs, N sweeps, Omega_I ...)`, `stopped at the sweep cap (...): <verdict>`, or `NOT certified (cg exit ...)` for a CG run that did not reach the gradient tolerance) and the subspace continuity (§2.6). The `wien2wannier` driver puts both lines at the top of the hand-off `.win`. | `hybrid_io.F` (`write_dis_status`) |
| `<file>.mmn.bin` | always (best effort) | binary sidecar read transparently by the standalone; delete to force a text re-parse | `main.F:784-822`; `hybrid_files.F:11-32` |
| `<file>.amn.bin` | always (best effort) | same, for the `.amn`. Written by `l2Amn` itself on the ordinary single-sigma path (`modules_rc.F:995, 1008-1023`); for SOC runs, where `l2Amn` is called with `write_files = .false.` and root assembles the combined block-structured `.amn`, it is written at `main.F:838-847` instead (inside `if (sigma2 .and. inwf%Amn .and. is_root)`, `main.F:829`). | `modules_rc.F:1008-1023`; `main.F:829, 838-847`; `hybrid_files.F:11-32` |

`<case>` is the nnkp filename with `.nnkp` stripped (`main.F:1012-1015`), so a
def file naming unit 11 `w2wH.nnkp` produces `w2wH-hyb.amn` and friends.

### 2.6 Checks that report or refuse

**Subspace continuity (both drivers).** After the disentanglement the hybrid
stage computes the smallest singular value of every hand-off overlap block
W(k)^H M(k,b) W(k+b) and prints
`subspace continuity: smallest singular value of the hand-off overlaps ... at k K -> K'; N of M (k,b) blocks below 0.5`.
For a subspace that varies smoothly in k the values are close to 1. A value
near 0 means that a state of the subspace at k has no partner in the subspace
at k+b: the subspace jumps there, no gauge can smooth the jump, and Wannier90
interpolates worse between those k-points and converges slowly or not at all.
Below 1e-2 the line is followed by `*** WARNING: the disentangled subspace JUMPS
between k K and k K'`. The value is also written to `<case>-hyb.dis`. It is a
diagnostic: nothing written changes. Measured on the regression cells with
`froz_emax = 2.0`: Fe 0.65 (up) and 0.78 (down), Fe-SOC 0.66, MgB2 0.50,
alpha-Sn 0.93, SrVO3 0.97, Bi bilayer 0.99, GaAs 2.3e-3. At the release
protocol (`froz_emax = auto`): alpha-Sn 0.50 (CG solver); GaAs 5.0e-4 with the
Gauss-Seidel solver (4.9 and 18.0 meV over [-16,0] and [0,4] eV, Wannier90
converged in 59,500 iterations) and 2.8e-2 with the CG solver (the default).
The Wannier90 settings, not only the hand-off, select the GaAs CG model. With
Wannier90's defaults the hand-off converges in 24,500 iterations at Omega_tot
51.40 A^2 (5.6 and 17.7 meV); the hand-off of the earlier CG stopping rule,
equal in Omega_I to six digits, converged in 104,500 iterations at 50.12 A^2
(5.2 and 14.0 meV), and with `precond = true` or `num_cg_steps = 200` at
57.57 A^2. A Stengel-Spaldin stage followed by a fresh Marzari-Vanderbilt run
with `num_cg_steps = 200` reaches 50.11 A^2 (5.3 and 14.1 meV) from both
hand-offs in about 3,500 iterations. The threshold
comes from the LCAO frontend (lcao2wannier), where SnTe at 7e-4 interpolated at
90 meV and the CG solver's subspace (1.2e-2) at 11 meV. The value does not rank
models by itself: GaAs warns with a good model. When it warns, compare the
CG solver (`--solver cg`; `solver = cg`) or a different frozen ceiling against
the band distance, e.g. with tools/band_compare.py.

**The sweep cap (both drivers).** See `sweeps` in §1.6; the default cap is
20,000 for HYBR and for `w2w_hybrid --sweeps`.

**A full disk refuses.** gfortran does not report a failed buffered write, so
a `w2w` run that filled the disk used to end normally (empty `case.error`,
`W2W END`) with a short or missing `.amn`, `.mmn`, sidecar or `.eig`. Every
product is closed with a size check; a mismatch stops the run with, e.g.,
`the .mmn text file case.mmn: 73420800 bytes on disk but 127623226 were written ... Refusing`,
rc 1, and `Error in W2W` left in `case.error`. The hand-off files
(`<out>.amn/.mmn/.eig/.anchors/.masks/.anatomy/.dis`, `.gauge0`, and the window
route's `<out>-win.*`) are checked the same way in HYBR mode and in
`w2w_hybrid`: a short one stops the run with
`ERROR STOP hand-off file incomplete (disk full?)`. Free space and rerun.

**`planew:` lines in `case.outputwf`.** When a neighbor k-point has more plane
waves than the head k-point's basis size, `planew` completes its interstitial
mask for that pair and prints one line per pair, e.g.
`planew: head k 3 neighbour k 2 G 0 0 0: 2 neighbour plane waves added to the interstitial mask`.
A `.mmn` written by the first assembly of 3.0 (2026-09, not published) for a case that prints such lines in this release
lost part of the interstitial overlap in exactly those blocks (MgB2 8x8x8:
pairs 3->2 and 7->8, up to 2e-4); no other measured material is affected.

**k-parallel vector files (`w2w case.def N`, `x w2w -p`).** After a k-parallel
`lapw1 -p` the vectors are in `case.vector_1..N`. The real/complex check
opens `case.vector_1` (it was skipped in split runs, so the real `w2w` on
complex split vectors wrote a wrong hand-off without an error), the FFT grid of
the interstitial overlap is the maximum over all parts (it was the last part's,
which left no margin on fcc splits), and
`planew: N plane-wave indices fall outside the interstitial grid ...` refuses
(rc 1) when an index would leave the grid, naming the first offender.

**Threads.** With a BLAS built with OpenMP (Homebrew OpenBLAS is one), the BLAS
thread count follows `OMP_NUM_THREADS` and ignores `OPENBLAS_NUM_THREADS`; unset
means every core. Set `OMP_NUM_THREADS` for every run (the `wien2wannier` driver
does). Any value from 2 up gives identical bytes; 1 can give a `.mmn` that
differs at round-off (up to 4e-16).

**b-vectors of the gauge spreads.** A `.nnkp` lists each k-point's neighbors in
its own order. The gauge policy (the seed tournament and the
`initial spread Omega` line) evaluates centers and Omega_D with the
b-vector of each (k,b) block; the first assembly of 3.0 used k = 1's list at every k, which
printed wrong Omega_D values (MgB2: 29.17 compared to Wannier90's 0.036). Omega_I
and Omega_OD, the hand-off files, and the tournament winners on the release
cells are unchanged.

---

## 3. Def-file input surface (SOC)

Two optional units turn a HYBR run into a SOC run. Nothing in the `.inwf`
controls this.

| Unit | Meaning | Effect | Source |
|---|---|---|---|
| `9` | second spin-component **vector** file | Presence sets `sigma2 = .true.` and enables the sigma loop: one process reads both spin vectors, accumulates the spinor-summed overlap in place, and assembles per-sigma `amn` column blocks. The AUTO projectability pass becomes a spinor sum over both components. | `main.F:122-125, 393-472, 689-761, 829-848` |
| `19` | second spin-component **vsp** (SP+SOC) | Swaps the radial potential for the sigma-2 pass. Without it, the sigma-2 pass reuses sigma-1's `vsp`, which is correct only for non-spin-polarized SOC. | `main.F:126-129, 416-421, 713-719` |

Fixed units: `unit_vector = 10`, `unit_nnkp = 11`, `unit_eig = 12`,
`unit_amn = 7`, `unit_mmn = 8`, `unit_ene = 50`, `unit_fermi = 51`,
`unit_vsp = 18`, `unit_struct = 20`, `unit_in = 5`, `unit_out = 6`.

A complete SP+SOC def file, a cleaned-up rendering of
`tests/alpha-Sn/sigw2w_v2.def` rather than a verbatim copy: the on-disk file
omits the comma after `'W.fermiup'` on line 11. That is harmless, because def
lines are read list-directed (`main.F:117`), where whitespace is a valid
separator.

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

---

## 4. Environment variables

`EXPERIMENTAL` marks options with no validated production use.

### 4.1 Projectability

| Variable | Default | Accepted | Path | Experimental | Meaning | Source |
|---|---|---|---|---|---|---|
| `W2W_TAILS` | unset (off) | exactly `1`; anything else = off (`character(len=8)`) | **HYBR/AUTO only** | no | Full-cell projector tails in the written `.amn`: `u_l` is continued past `R_MT` with a log-derivative-matched decaying exponential, its overlap with the interstitial plane-wave sum is added to the MT matching, and the column is renormalized by `tail_nrm(l)`. | `modules_rc.F:608-641, 947-975` |
| `W2W_PAMN` | unset (P source = the seed `.amn`) | path to an `.amn` (`character(len=256)`); must match the seed `.amn` in `nb` and `nk` | **standalone only** | no | Read a separate `.amn` as the projectability source, so tail-extended full-cell projectors drive P and the masks while the plain MT `.amn` feeds the SCDM gauge and the hand-off. HYBR computes both variants in one plane-wave sweep and uses `projT_h` directly (`main.F:903-930`), so setting this under HYBR does nothing. | `w2w_hybrid.F:178-189` |

`W2W_TAILS` is the exact mirror of `W2W_PAMN`: it is read in exactly one
place, `l2Amn` (`modules_rc.F:631`), which is part of the `w2w`/`w2wc` binary
only. The standalone never sees it; `Makefile.hybrid:56` links only
`hybrid_core.o hybrid_grass.o hybrid_prep.o hybrid_io.o hybrid_files.o
w2w_hybrid.o`, with no `modules_rc.F`; the variables the standalone reads
are listed by `w2w_hybrid --help`.

Under HYBR, `W2W_TAILS` is read inside `l2Amn`
(`modules_rc.F:631-632`). HYBR forces the tails *machinery* on regardless
(`do_tails = env_tails .or. want_full`, `modules_rc.F:637`) so it can export
`projT_h`, but `W2W_TAILS=1` additionally swaps the **written** `.amn` and the
in-memory `proj_h` (which feed the SCDM gauge *and* the hand-off) to
full-cell content (`modules_rc.F:957-976`).

### 4.2 Gauge

| Variable | Default | Accepted | Path | Experimental | Meaning | Source |
|---|---|---|---|---|---|---|
| `W2W_SEEDK` | unset -> falls back to the caller's `seedk_in` (HYBR: `opt16`); if neither, seed `k = 1` | `opt[<N>]` \| `best` \| `median` \| `<integer k>` (`character(len=8)`, so longer values truncate) | both | no | Transport-seed policy; overrides the `.inwf` `seedk` keyword. See the `seedk` row in §1.6. | `hybrid_core.F:1854-1857, 1865-1991` |
| `W2W_SEEDLOSS` | unset -> caller's `seedloss_in` (HYBR: `smooth`); if neither, `smooth` | `smooth` \| `blind` \| `local` \| `abel` \| `absm` (`character(len=8)`) | both | no | Tournament loss functional; overrides the `.inwf` `seedloss` keyword. Unrecognized strings fall into the `case default` (= `smooth` behavior) but still print the `seed loss mode:` notice. | `hybrid_core.F:1858-1863, 1105-1129` |
| `W2W_SEEDGATE` | `2.0e-3` | real (`character(len=16)`); unparseable value silently reverts to the default | both | no | Candidate stability gate for `seedk = opt<N>`: k-points whose per-k anchor `sigma_min` is below this are excluded from the candidate pool before sampling. If nothing passes, the code prints `NO seed passes the stability gate` and samples all k. It does **not** change the separate hardcoded `sig_gate = 2e-3` that drives the `auto<N>`/`tfl<N>` keep-vs-transport switch (`hybrid_core.F:2005`). | `hybrid_core.F:1883-1892, 1919-1928` |
| `W2W_SMTOL` | unset -> `-1` (off): always run the full sweep count | real, e.g. `1e-8` (`character(len=16)`); unparseable value silently reverts to `-1` | both | no | Converged-smoothing mode for `smooth_gauge`: exit as soon as the largest per-sweep `max\|dU\|` drops below the tolerance, treating `nsweep` as a cap. Honored inside every seed-tournament candidate too, so it changes tournament *scores*, not just wall time. | `hybrid_core.F:1206-1215, 1232-1243` |
| `W2W_ANCHORS` | unset (unconstrained SCDM pivots) | exactly `shellcomplete` (`character(len=16)`) | both | **yes** | Constrain the SCDM QRCP pivots to whole `(atom,l,rf)` shells across both sigma blocks, using the channel groups in `projmap.dat`, so the anchor set is a representation basis (required by WannSymm / irrep tools). Costs spread, since it spends WF budget on orbitals SCDM declined. Falls back to unconstrained SCDM with a warning if `projmap.dat` is missing or unusable. | `hybrid_core.F:1774-1796, 1379-1422` |

### 4.3 Diagnostics

| Variable | Default | Accepted | Path | Experimental | Meaning | Source |
|---|---|---|---|---|---|---|
| `W2W_VERBOSE` | unset (off) | exactly `1` (`character(len=8)`) | **standalone only** | no | Per-sweep `Omega_I` trace from the disentangler. Hour-scale runs are otherwise silent, no progress, no sweep count, no ETA. MODE=HYBR passes the literal `.false.` for `verbose` (`main.F:1003`), so this variable is inert in the one-shot. | `w2w_hybrid.F:154-162` |

### 4.4 Compile-time constants often mistaken for variables

| Name | Value | Meaning | Source |
|---|---|---|---|
| `AUTOW_GAP` | 5.0 eV | minimum inter-band gap for the `bmin = 0` gap-walk to treat a gap as a core/valence separator | `main.F:101, 290-292` |
| `AUTOW_BELOW` | 1.0 eV | the band above a candidate gap must still start below `E_F - AUTOW_BELOW` | `main.F:101, 291-292` |
| `AUTOW_ECAP` | 40.0 eV | energy cap for the `bmax = -1` sentinel | `main.F:102, 268-275` |
| `P2_ONSET` | 5.0 eV | onset of the P2 erfc ranking weight above `E_F` | `main.F:103, 380-384` |
| `P2_WIDTH` | 2.0 eV | width of the P2 erfc ranking weight | `main.F:103, 380-384` |
| `sig_gate` | 2.0e-3 | per-k `sigma_min` threshold for the `auto<N>`/`tfl<N>` keep-vs-transport switch (distinct from `W2W_SEEDGATE`) | `hybrid_core.F:2005` |
| block-SCDM gap | 0.5 eV | manifold-isolation gap used by `bscdm`/`bsms<N>` | `hybrid_core.F:1768` |
| `dtol` | 0.02 eV | degeneracy tolerance in `build_masks` | `hybrid_prep.F:519` |
| degenerate-fraction threshold | 0.7 | `auto_regime` picks band-level above this | `hybrid_prep.F:523` |
| `p_floor` | 0.10 | per-state admission threshold in `select_pool_and_frozen` (discarded by HYBR's `admit_h = .true.`) | `hybrid_prep.F:531, 345` |
| `fermi_shell1`, `shell_guard` | 3.0 eV, 0.15 | Fermi-shell and shell-guard parameters of `select_pool_and_frozen` | `hybrid_prep.F:531-532` |
| `thr_floor`, `margin`, `eps`, `fid_emax` | 0.55, 0.015, 0.02, 8.0 eV | `auto_trust` threshold constants | `hybrid_prep.F:276` |
| froz walk-down | start 8 eV, step 1 eV, floor 1 eV | `froz_emax = auto` schedule; exits on the frozen-demand test alone, the k-point-saturation test was removed | `hybrid_prep.F:433, 442-443` |

There is no identifier with a `P3_` prefix anywhere in `SRC_w2w`; "P3" is only
a nickname used in comments for the support gate, whose only parameter is the
line-4 `support_ratio` field.

---

## 5. Gauge modes in detail

The dispatcher is `choose_gauge` (`hybrid_core.F:1743-2108`). It has two
stages.

**Anchor stage** (`hybrid_core.F:1763-1807`). A mode string starting `bs`
selects block-SCDM (per spectrally isolated manifold, 0.5 eV gap,
`hybrid_core.F:1763-1772`). A string starting `low` or `lds` replaces the QRCP
anchors with a Löwdin gauge on the leading `nw` columns
(`hybrid_core.F:1801-1807`). Everything else uses global SCDM/QRCP anchors
(`hybrid_core.F:1794`).

**Propagation stage** (`hybrid_core.F:1992-2099`):

| Mode | Behavior |
|---|---|
| `scdm` | bare per-k SCDM gauge; no propagation, no smoothing |
| `ptg` | parallel transport from the seed k |
| `auto<N>` | if worst-k SCDM `sigma_min >= 2e-3`, keep the per-k SCDM gauge, then N smoothing sweeps; else transport from the best-conditioned k and then N sweeps, **or** (if `seedk` is an `opt` policy) run the seed tournament, in which case the N sweeps happen *inside* the tournament (see below) |
| `tfl<N>` | keep per-k SCDM at healthy k, BFS-transport into sick k, then N sweeps |
| `sms<N>`, `lds<N>`, `bsms<N>` | N smoothing sweeps only |
| `pts<N>` | transport + N sweeps, or the tournament (again with the N sweeps inside it) |

**The tournament sub-case does not get a second smoothing pass.** When
`auto<N>` takes the tournament branch, the code executes `exit autotfl`
immediately after `seed_tournament` returns (`hybrid_core.F:2022-2031`), which
skips the trailing `call smooth_gauge` at `hybrid_core.F:2065`. The N sweeps
are applied *within* the tournament instead, each candidate is smoothed with
`nsw` sweeps before scoring and the winner's frame is rebuilt with them
(`hybrid_core.F:1126, 1178`). The printed line says so: `... smoothing sweeps
(in tournament)`. `pts<N>` has the identical structure; `exit ptsblk` at
`hybrid_core.F:2091`, skipping the `smooth_gauge` at `hybrid_core.F:2095`.

Two pitfalls:

- **Suffix parsing falls back to 8 silently.** Unreadable digits give
  `nsw = 8` (`hybrid_core.F:2013-2014, 2076-2077`), so `gauge = auto` is
  identical to `auto8`, and `gauge = autoX` is too.
- **There is no validation of the gauge string.** Anything that matches none
  of `bs*`/`low*`/`lds*`/`ptg`/`auto*`/`tfl*`/`sms*`/`pts*` falls through the
  dispatcher with no else-branch and silently ships the bare per-k SCDM
  gauge. `gauge = scdm` is a legitimate way to say "no smoothing";
  `gauge = auto8x` is a typo that still means auto with 8 sweeps.

A third consequence: under the default `auto8`, whenever the SCDM worst-k
`sigma_min` is at or above `2e-3`, the per-k SCDM gauge is kept and only
smoothing runs, so `seedk`, `seedloss`, `W2W_SEEDK`, `W2W_SEEDLOSS` and
`W2W_SEEDGATE` have **no effect at all**. The candidate list is still computed
and printed (`hybrid_core.F:1939-1940`), which makes the no-op easy to miss.
These options take effect only for `ptg`, `pts<N>`, or `auto<N>` when `sigma_min`
falls below the gate.

---

## 6. Standalone `w2w_hybrid` options

Built by `make -f Makefile.hybrid` (serial `w2w_hybrid`) or
`make -f Makefile.hybrid acc` (`w2w_hybrid_acc`, OpenMP build; identical
options). In the OpenMP build `OMP_NUM_THREADS` sets the threads of the
conjugate-gradient solver (parallel over k), of the colored sweep, and of
the serial Gauss–Seidel sweep's neighbor products (an OpenMP loop over b
when `nd*nw*npool` reaches 2e5); the conjugate-gradient result and the
serial sweep's result do not depend on the thread count. Options are named `--key value` pairs; `--seed` and `--nw` are
mandatory, everything else has a default. The bare shorthand
`w2w_hybrid <case> <N> [options]` is also accepted. `w2w_hybrid --help`
prints the authoritative reference.

| Option | Default | Meaning | Source |
|---|---|---|---|
| `--seed` | none (mandatory) | Seedname. Reads `<seed>.amn/.mmn/.eig/.nnkp`, each with a transparent `<file>.bin` sidecar (used when present and not older than the text; a missing text file is not staleness, so the `.bin` alone suffices). | `w2w_hybrid.F:49, 166-174`; `hybrid_files.F:11-32` |
| `--nw` | none (mandatory) | Target WF count = hand-off dimension. Error-stops if `nw > nproj`. There is no `0` sentinel here; the auto-capacity sentinel exists only in the `.inwf` HYBR path. | `w2w_hybrid.F:50, 295-299` |
| `--out` | `<seed>-hybF` | Output seedname for `<out>.amn/.mmn/.eig` plus `<out>.anchors`, `<out>.masks`, `<out>.anatomy` and `<out>.dis` (see the output-file table above). | `w2w_hybrid.F:51-52, 362-371` |
| `--window` | off | The automatic window route, in addition to the hand-off: the frozen states' energy hull becomes `dis_froz_min/max` (lowered where a k-point would hold more than `nw` frozen states), the outer window keeps the pool bottom and takes the ceiling that minimizes the frozen-constrained Omega_I among the pool top and fixed offsets above the frozen ceiling; writes `<out>-win.win` (a snippet with the windows and every Omega_I tried), `<out>-win.amn` (the anchor columns of the seed projections over all `nb` bands), `<out>-win.mmn` (the seed overlaps) and `<out>-win.eig`, for stock Wannier90 disentanglement. The public driver's `--method window` completes the snippet with the seed geometry and k-points. | `w2w_hybrid.F` (`window_route_emit`), `hybrid_io.F` (`write_window_route`) |
| `--rescale` | `occ99` | How P is built and normalized. `lowdin`: `p_i(k) = diag(A S^-1 A^H)` over the selected columns, exactly bounded `[0,1]`; the code then rewrites `rescale` to `none` so `build_masks` only clips (so `lowdin` and `share` cannot be combined). `occ99`: divide by the 99th percentile of P over occupied states. `p99`: 99th percentile over all states. `max`: divide by `max(P)`. `none`: no scaling. `share`: mask-column P as a fraction of all-column P (warns if the mask is the full set). | `w2w_hybrid.F:53-54, 255-293`; `hybrid_prep.F:25-62, 191-213` |
| `--regime` | `auto` | `auto` \| `per-state` \| `band-level`; identical to the `.inwf` line-5 `regime` field. | `w2w_hybrid.F:55, 304-306` |
| `--froz-emax` | `auto` (value 8.0) | `auto` \| `none` \| `<real eV>`; identical to the `.inwf` `froz_emax` keyword. | `w2w_hybrid.F:56-59`; `hybrid_prep.F:417-451` |
| `--pocket-skip` | `on` | `on` \| `off`; acts only with `--froz-emax auto`. When the weakest hand-off link (smallest singular value of `W^H M W`) is below 0.5 and touches a k-point where a pool band is frozen that lies below the walked ceiling at fewer than 25% of the k-points and whose frozen states average P >= 0.5 (a mixed-character pocket), the ceiling drops 0.1 eV below that band, the masks are rebuilt and the disentanglement repeats; the result is kept only if the new weakest link is >= 0.5, otherwise the first pass is restored (up to 3 steps; never below 1 eV; the pool must not change). The `.inwf` HYBR path reports the same condition and the `froz_emax` to pass, but does not act on it. Same rule as lcao2wannier `--hybrid-pocket-skip`. | `w2w_hybrid.F` (`pocket_loop`); `hybrid_prep.F` (`pocket_ceiling`); `hybrid_io.F` (`continuity_links`) |
| `--efermi` | `0.0` | Fermi level in the same units as `<seed>.eig`. w2w writes `.eig` already Fermi-shifted and in eV (`modules_rc.F:2515`), so `0` is the correct value for w2w-produced files. | `w2w_hybrid.F:60-62`; `hybrid_io.F:43-55` |
| `--sweeps` | `20000` | Maximum `Omega_I` disentanglement sweeps (or CG gradient evaluations). A Gauss-Seidel run that stops at the cap is graded `flat tail` or `POOR CONVERGENCE: still descending` in the log and in `<out>.dis` (§2.6). HYBR's equivalent is the `sweeps` keyword (default 20000). | `w2w_hybrid.F` (`dis_niter`); `hybrid_core.F` (`cap_tail_text`) |
| `--solver` | `cg` | `cg` \| `gs` \| `auto`; identical to the `.inwf` `solver` keyword (default `cg` since 2026-10, `gs` before). `--rtol`, `--check-every`, `--colored`, `--adaptive` and `--relax` act on `gs` only. | `w2w_hybrid.F:79, 315-459, 723-727` |
| `--rtol` | `1e-10` | Relative `Omega_I` convergence tolerance (the absolute tolerance is hardcoded `1e-9` in the call). | `w2w_hybrid.F:68-70, 347` |
| `--check-every` | `10` | Convergence-check cadence in sweeps. | `w2w_hybrid.F:71-73, 347` |
| `--colored` | off | Red-black / graph-colored sweep ordering; a flag without a value. HYBR keyword `colored`. | `w2w_hybrid.F:74-78`; `hybrid_core.F:449-457` |
| `--adaptive` | off | Adaptive per-k convergence signal driving the check cadence; a flag without a value. HYBR keyword `adaptive` (default `true`). | `w2w_hybrid.F:79-83`; `hybrid_core.F:458-460, 504-508` |
| `--relax` | `1.0` | Zd mixing; negative = scheduled extrapolation from `\|value\|`. Same parameter as the `.inwf` `mix` field. | `w2w_hybrid.F:84-87`; `hybrid_core.F:354-356, 465-475` |
| `--gauge` | `scdm` | Gauge dispatcher mode (see §5). **Note the default differs from HYBR's `auto8`.** | `w2w_hybrid.F:88-89` |
| `--anderson` | `0` (off) | Safeguarded alternating Anderson acceleration history depth. Every 5 sweeps it extrapolates the free-block subspace projectors over the last `m+1` iterates, re-extracts per-k subspaces, and accepts only if the free-block overlap sum improves after one plain sweep, else reverts. Forces `dis_adaptive` on. HYBR keyword `anderson`. | `w2w_hybrid.F:90-93`; `hybrid_core.F:357-363, 476-483, 539-563` |
| `--np-mask` | `0` (full P) | Which `amn` columns feed P (and hence the masks). `>0`: the leading `np_mask` columns; for SOC files (`nproj = 2*rows` in `projmap.dat`) the selection is mirrored into the second sigma block so the P span is not sigma-asymmetric. `<0`: all columns with `l <= \|np_mask\|` (needs `projmap.dat`; also sigma-mirrored). `0`: all columns, it then reads `projmap.dat` only to print a mask hint explicitly annotated `(not adopted; np_mask=0 = full P)`. HYBR's equivalent is automatic: `auto_config%mask_cols` from the cliff-mode Otsu boundary. | `w2w_hybrid.F:94-101, 130-151, 190-254`; `main.F:909-927` |
| `--init-transport` | off | BFS parallel transport of the seed subspace across the k-mesh before disentangling (each child takes the `nt`-dimensional non-frozen subspace of maximal overlap with its parent through the `mmn` link). A flag without a value. Not exposed in HYBR. | `w2w_hybrid.F:102-113`; `hybrid_core.F:349-353, 410-412` |
| `--admit-all` / `--no-froz` / `--all-free` | admit-all (default) | `--admit-all`: admit every pool state (overrides `--p-floor`); `--no-froz`: no frozen constraint; `--all-free`: both | `w2w_hybrid.F:114-120, 312-319` |
| `--p-floor` | none (admit all) | Projectability floor of the per-state admission mask: a pool state outside the frozen set may enter the disentangled subspace only if its P (after `--rescale`) is at least this value. Without it every pool state is admitted, as in HYBR and the paper's LAPW protocol (default since 2026-10-06; the floor was 0.10 before). An explicit `--admit-all` overrides it. On the eight LAPW cells a 0.10 floor excludes 26-7292 states per cell, breaks the subspace continuity of the Bi bilayer and Fe-SOC and raises their band errors 7-28 times, and worsens alpha-Sn (2.2/14.2 -> 4.2/17.5 meV); it helps MgB2 (11.5/4.4 -> 4.8/1.5 meV). | `w2w_hybrid.F` (`p_floor`, `pfloor_set`); `hybrid_prep.F` (`build_masks`, `p_floor_in`) |
| `--pool-factor` | `1.5` | `npool` **floor** `= max(nw+1, ceil(pool_factor*nw))`. Echoed as `  pool factor = ...` when given. It is not a cap: `select_pool` takes `npool = max(n_pool_max, s1-s0+1)`, so P-based selection widens the pool whenever the `nw` most-projectable in-window bands span more indices than that, and the `froz_emax` pool guard can widen it further. The only hard cap is the available band count. | `w2w_hybrid.F:121-129`; `hybrid_prep.F:248, 261-265, 509` |

### 6.1 The validated production recipe

This is the invocation in `tests/regress_hybrid.sh`, with `--solver gs` added:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2

W2W_PAMN=$PWD/<seed>.amn.tl W2W_SEEDK=opt16 W2W_SEEDLOSS=smooth \
  path/to/SRC_w2w/w2w_hybrid_acc --seed <seed> --nw <nw> --out <out> \
    --rescale lowdin --regime per-state --froz-emax 2.0 --sweeps 4000 \
    --colored --adaptive --relax -1.9 --gauge auto8 --np-mask <np> --admit-all \
    --solver gs
```

(The default solver became `cg` in 2026-10. The harness's archived hand-offs
were made with the Gauss-Seidel sweep, and the `--colored`, `--adaptive` and
`--relax` options act on that solver only, so the recipe selects it with
`--solver gs`.)

Option by option, with the reason each non-default value is there:

| Option | Value | Default? | Why |
|---|---|---|---|
| `--rescale` | `lowdin` | non-default (`occ99`) | Löwdin P is exactly bounded `[0,1]`, the correct normalization for overlapping tail-extended projectors |
| `--regime` | `per-state` | non-default (`auto`) | fixes the mask regime instead of letting the degenerate-pair fraction decide |
| `--froz-emax` | `2.0` | non-default (`auto`) | explicit ceiling; also what makes the capacity arithmetic well-defined |
| `--efermi` | `0.0` | default | w2w already writes `.eig` Fermi-shifted in eV |
| `--sweeps` | `4000` | non-default (20000) | the harness's cap; the regression cells converge well below it |
| `--rtol` | `1e-10` | default | — |
| `--check-every` | `10` | default | — |
| `--solver` | `gs` | non-default (`cg`) | the archived hand-offs are Gauss-Seidel results |
| `--colored` | (flag) | non-default (off) | colored sweep ordering |
| `--adaptive` | (flag) | non-default (off) | per-k convergence signal drives the check cadence |
| `--relax` | `-1.9` | non-default (1.0) | scheduled extrapolation from 1.9, halved at each convergence check |
| `--gauge` | `auto8` | non-default (`scdm`) | matches the HYBR default: keep per-k SCDM when it is well conditioned, 8 smoothing sweeps |
| `--anderson` | `0` | default | Anderson acceleration off |
| `--np-mask` | `<np>` | non-default (0) | the physical/mask column count, mirrored across sigma blocks for SOC |
| `--init-transport` | not passed | default | transport initialization off |
| `--admit-all` | (flag) | default since 2026-10-06 (was a floor of 0.10) | every pool state admitted |
| `--pool-factor` | not passed | default 1.5 | omitted, so the default applies |

`<np>` values used by the harness: 18 GaAs, 12 MgB2, 9 Fe (per spin), 9
Fe_soc, 10 aSn, 6 Bi2D, 19 SrVO3. `<nw>`: 18 GaAs, 12 MgB2, 9 Fe per spin, 18
Fe_soc, 20 aSn, 12 Bi2D, 19 SrVO3. One cell (`fednab`) is run with
`gauge = pts8` and `W2W_SEEDLOSS=absm` specifically to keep the abelian/seam
code path covered.

The `<seed>.amn.tl` file passed as `W2W_PAMN` is produced by a separate w2w
run with `W2W_TAILS=1` whose def file points unit 7 at that name (see
`tests/alpha-Sn/lowtl_bench/bench_tl.def`).

### 6.2 Reproducing HYBR's settings from the standalone

To make the standalone follow the same trajectory as the one-shot, mirror
HYBR's disentangler defaults (`sweeps = 20000`, `solver = cg`, no coloring, no
Anderson, no transport init) and pass the cliff-mode
mask boundary explicitly. This is the parity invocation in
`tests/regress_cases/run_parity_cell.sh`:

```sh
mc=$(awk '/# mask_cols/{print $3}' projmap.dat)   # 0 if absent
W2W_PAMN=$PWD/<prefix>T.amn W2W_SEEDK=opt16 \
  path/to/SRC_w2w/w2w_hybrid_acc --seed <prefix>H --nw <nw> --out par-<prefix> \
    --rescale lowdin --regime per-state --froz-emax 2.0 --gauge auto8 \
    --np-mask $mc --admit-all
```

---

## 7. Recommended settings

The three blocks below are the validated production recipe expressed in the
one-shot HYBR path. The only thing that differs between calculation types is
the def file.

### 7.1 Non-magnetic, no SOC

`case.inwfhyb`:

```
HYBR
  0  -1
  3
0 0 T T -1
0 1.0 per-state
froz_emax = 2.0
```

`w2wH.def`:

```
  5,'case.inwfhyb',  'old',    'formatted',  0
  6,'w2wH.outputwf', 'unknown','formatted',  0
  7,'w2wH.amn',      'unknown','formatted',  0
  8,'w2wH.mmn',      'unknown','formatted',  0
 10,'case.vector',   'old',    'unformatted',0
 11,'w2wH.nnkp',     'old',    'formatted',  0
 12,'w2wH.eig',      'unknown','formatted',  0
 18,'case.vsp',      'old',    'formatted',  0
 20,'case.struct',   'old',    'formatted',  0
 50,'case.energy',   'old',    'formatted',  0
 51,'case.fermi',    'old',    'formatted',  0
```

Run it with the build that matches the case; `w2w` for the real build,
`w2wc` for the complex one (`-D_COMPLEX_`), and feed
`w2wH-hyb.{amn,mmn,eig}` to a stock Wannier90 run with
`num_bands = num_wann` and no `dis_*` keywords:

```sh
path/to/SRC_w2w/w2w w2wH.def
```

Non-defaults and why:

- `bmin = 0`, `bmax = -1` (defaults do not exist; both fields are mandatory):
  gap-walk core trimming plus an `E_F + 40 eV` cap, so the window is set by
  the band structure rather than by hand.
- `threshold = 0`, non-default (0.10). Cliff mode: every non-noise channel
  becomes an SCDM candidate, and the physical/gauge boundary is set by the
  support Otsu cut instead of a relative-weight floor.
- `support_ratio = -1`, non-default (0.5). Activates that Otsu cut and
  records `mask_cols`, which is what the P source restricts to. Cliff mode
  needs *both* this and `threshold <= 0`; either one alone silently degrades.
- `write_proj = T`, non-default (F). Produces `projectability.dat` and
  `autoprojection.diag`; the cost is one extra array, and without it a failed
  selection is very hard to diagnose.
- `p2_weight = T`, default. Keep it on: the cliff-mode ranking assumes the
  erfc energy weight has already suppressed high-energy NFE channels.
- `num_wann = 0`, non-default in effect (it is the type default, but every
  file must state it). Capacity-driven auto sizing; requires the numeric
  `froz_emax` below.
- `mix = 1.0`, non-default (`-2.5`): relaxation off; it acts only with `solver = gs`.
- `regime = per-state`, non-default (`auto`).
- `froz_emax = 2.0`, non-default (`auto`, 8 eV walk-down). An explicit
  numeric ceiling is what enables the pre-`l2Amn` capacity gate and the
  `num_wann = 0` sentinel. 2.0 eV is the validated value for the small
  benchmark cells; it is a physics choice about which bands must be
  reproduced exactly, not a numerical parameter.

Everything else (`pool_factor = 1.5`, `margin = 0`, `gauge = auto8`,
`seedk = opt16`, `seedloss = smooth`, `capacity = strict`) is already at the
validated value by default, so those lines can be omitted. The solver
defaults to `cg` since 2026-10; add `solver = gs` to reproduce a hand-off
made before.

### 7.2 Collinear spin-polarized, no SOC

Two independent HYBR runs, one per spin. Nothing in the Stage-2 code is
spin-aware; each spin channel produces its own hand-off and its own Wannier90
run. Use the same `.inwf` for both (or two copies with different `num_wann`)
and two def files that differ only in the spin-dependent files:

```
 10,'case.vectorup','old','unformatted',0     |  10,'case.vectordn','old','unformatted',0
 18,'case.vspup',   'old','formatted',  0     |  18,'case.vspdn',   'old','formatted',  0
 50,'case.energyup','old','formatted',  0     |  50,'case.energydn','old','formatted',  0
 51,'case.fermiup', 'old','formatted',  0     |  51,'case.fermidn', 'old','formatted',  0
```

The two `.nnkp` files should have identical content, so the two models share a
k-mesh and b-vector set, but they need **distinct filenames**: the hand-off
name is derived from the nnkp filename (`main.F:1012-1015`), so reusing one
name would make the second run overwrite the first one's hand-off.
`tests/regress_cases/fe_w2w/{feupH,fednH}.def` are a working pair, same
`hybr.inwf`, per-spin `vector`/`vsp`/`energy`/`fermi`, per-spin nnkp name.
With an explicit `num_wann` (that file uses 9 per spin for bcc Fe) the two
spins can be given different values if their capacities differ.

### 7.3 SOC (with or without spin polarization)

One HYBR run. Add unit 9 (second spin-component vector) to the def file; add
unit 19 as well when the calculation is spin-polarized, so the sigma-2 pass
uses the down-spin radial potential:

```
  5,'case.inwfhyb',    'old',    'formatted',  0
  6,'w2wH.outputwf',   'unknown','formatted',  0
  7,'w2wH.amn',        'unknown','formatted',  0
  8,'w2wH.mmn',        'unknown','formatted',  0
 10,'case.vectorsoup', 'old',    'unformatted',0
 11,'w2wH.nnkp',       'old',    'formatted',  0
 12,'w2wH.eig',        'unknown','formatted',  0
 18,'case.vspup',      'old',    'formatted',  0
 20,'case.struct',     'old',    'formatted',  0
 50,'case.energysoup', 'old',    'formatted',  0
 51,'case.fermiup',    'old',    'formatted',  0
  9,'case.vectorsodn', 'old',    'unformatted',0
 19,'case.vspdn',      'old',    'formatted',  0
```

The `.inwf` is unchanged from §7.1. Note that under the sigma loop:

- the `.amn` holds one column *block* per spin component
  (`main.F:654, 829-837`), and `nproj_h = 2*Nproj`;
- the AUTO projectability is a spinor sum over both components
  (`main.F:393-472`);
- the auto-`num_wann` walk adds `2*(2l+1)*mult` columns per channel
  (`main.F:533, 562`);
- the mask column selection is mirrored into the second sigma block
  (`main.F:916-919`).

For non-spin-polarized SOC, omit unit 19; the sigma-2 pass then reuses
sigma-1's `vsp`, which is correct in that case. SOC runs use the complex
build (`w2wc`).

### 7.4 Post-hoc, on files that already exist

Use the standalone. This is the right choice when the `.amn`/`.mmn` are
expensive and you want to sweep Stage-2 options, or when you need one of the
options HYBR does not expose (`init_transport`, `nofroz`, `--p-floor`,
`W2W_VERBOSE`, `W2W_PAMN`). The full command is
in §6.1.

---

## 8. Options to change only for a specific reason

- **`rescale = ...` in a HYBR `.inwf`.** Not a keyword. It aborts the run.
- **`W2W_ANCHORS=shellcomplete`**, experimental. It constrains SCDM pivots to
  whole shells so the anchor set is a representation basis; it costs spread
  because it spends WF budget on orbitals SCDM declined. Only use it when a
  downstream symmetry tool requires it.
- **`capacity = warn`.** It downgrades exactly two croaks and changes no
  physics. The frozen budget is still spent deepest-first per k
  (`hybrid_prep.F:596-610`), so band fidelity near `froz_emax` is no longer
  guaranteed; check `eta` afterwards.
- **`gauge` strings other than `auto8`.** There is no validation: a typo
  silently ships the bare per-k SCDM gauge. `ptg`, `pts<N>`, `tfl<N>`,
  `lowdin`, `lds<N>`, `bscdm`, `bsms<N>` are all implemented but only `auto8`
  (and `pts8` on one regression cell) is exercised by the harness.
- **`seedk` / `seedloss` and their environment overrides.** Under the default
  `auto8` with a healthy `sigma_min` they do nothing at all, while still
  printing a candidate list. Changing them is only meaningful for `ptg`,
  `pts<N>`, or `auto<N>` below the gate.
- **`W2W_SMTOL`.** It is honored inside every seed-tournament candidate, so
  it changes tournament scores, not only wall time. Results with and without
  it are not comparable.
- **Standalone `--colored`, `--adaptive`, `--anderson`, `--init-transport`.**
  They alter the disentanglement trajectory. The production recipe passes
  `--colored --adaptive` and neither of the others; anything else is unvalidated.
  `anderson > 0` additionally forces `adaptive` on.
- **Standalone `--np-mask 0`** prints a `projmap mask hint: N columns` line that
  is explicitly *not* adopted. It is informational; it does not set anything.
- **Standalone `--no-froz` / `--all-free`.** These remove the frozen
  constraint entirely, giving the pure `Omega_I` optimum with no guarantee
  that any band is reproduced. Diagnostic use only.
- **`AUTOW_*` and `P2_*`.** Compile-time `parameter`s. Editing them changes
  the auto band window and the channel ranking for every material, and
  requires a rebuild.
