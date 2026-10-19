# SRC_w2w INDEX (release v3.0)
Updated: 2026-10-21

The hybrid pre-Wannierization pipeline, without the symmetry-adapted
features of the development branches. `w2w_hybrid --help` is the
authoritative CLI reference; user documentation lives in
`../doc/hybrid/00–06`, and `../REFERENCES.md` lists the literature behind
each file.

## Sources
| File | Purpose / when to use | Updated | Status |
|---|---|---|---|
| `main.F` | w2w main program: stock modes + `MODE = HYBR` one-shot (Stage 1+2 fused); products closed with a size check; k-parallel vector-kind probe; HYBR writes `<case>-hyb.dis` | 2026-10-01 | canonical |
| `w2w_hybrid.F` | Standalone hybrid driver + named-option CLI; reads `<seed>.amn/.mmn/.eig/.nnkp`, writes the `Nb=Nw` hand-off and `<out>.dis` (sweep cap 20000) | 2026-10-01 | canonical |
| `hybrid_prep.F` | Projectability, masks, pool/frozen selection (frozen ceiling: `auto` walks 8 eV down in 1 eV steps until the frozen demand fits num_wann, no k-point-saturation test), MV b-weights, capacity arithmetic | 2026-09-14 | canonical |
| `hybrid_core.F` | Frozen-constrained Ω_I disentanglement (Gauss–Seidel; coloring, mixing/schedules, Anderson; graded classification at the sweep cap; the neighbor products of one k update run as an OpenMP loop over b in the OpenMP build when large, with thread-independent results), subspace rotation, gauge construction/propagation (spreads with per-(k,b) b-vectors) | 2026-10-11 | canonical |
| `hybrid_grass.F` | Grassmann-CG disentangler, the default solver (`--solver cg`; parallel over k, no coloring; stops on the projected gradient; a k-point with no free rows is skipped in the constant terms) | 2026-10-11 | canonical |
| `hybrid_io.F` | `.eig` I/O, hand-off writers + `validate_handoff` final acceptance (full-coverage contractivity, finiteness, SVD success, refuses before writing); `handoff_close` size check; `subspace_continuity`; `.dis` status writer | 2026-10-01 | canonical |
| `hybrid_files.F` | `.amn/.mmn/.nnkp` readers with transparent `.bin` sidecar caching; `nnkp_block_bvectors` (the b-vector of every (k,b) block) | 2026-10-01 | canonical |
| `modules_rc.F` | Real/complex Bloch matrix-element machinery (stock w2w core); `planew_mask_extend`, `planew_check_grid`, grid maximum over k-parallel vector files | 2026-10-01 | canonical (upstream+mods) |
| `modules.f`, `mod_w2w.f`, `util_w2w.F`, `projection_support.F`, `mpi_w2w.F` | Stock w2w support modules (+ `.inwf` HYBR keyword parsing in `util_w2w.F`); two unused modules (`rfovl_m`, `stepfn_m`) removed from `modules.f` | 2026-10-11 | canonical (upstream+mods) |

## Build
| File | Purpose | Updated | Status |
|---|---|---|---|
| `Makefile` | `make real` / `complex` / `all` → `w2w` / `w2wc` (needs `../make.sys`) | 2026-09-07 | canonical |
| `Makefile.hybrid` | `make -f Makefile.hybrid` → `w2w_hybrid` (serial); `acc` target → `w2w_hybrid_acc` (OpenMP), leaving the serial program as it is | 2026-10-12 | canonical |

## Records
| File | Purpose | Updated | Status |
|---|---|---|---|
| `bench/results.csv`, `bench/README.md`, `bench/report.html` | Timings recorded while the performance features were added (colored OpenMP sweep, binary sidecar I/O, relaxation schedules, adaptive monitor, Anderson acceleration, Alm cache, HYBR one-shot, MPI check); the cases they name are development-tree runs | 2026-07-19 | supporting |

## Environment variables
The release keeps a small set (gauge seeding, tolerances, anchor grouping,
CG restart, verbosity); the complete list with defaults is in
`w2w_hybrid --help`, and `../doc/hybrid/03-options-reference.md` §4
documents the main ones.

## Not in this release
The symmetry-adapted features (site-symmetry/.dmn, little-group projection,
sewing dumps, the reduced solver) are on the development branches only.
