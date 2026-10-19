# wien2wannier v3.0 INDEX
Updated: 2026-10-21

## Programs
| Tool | Purpose / when to use | Inputs → Outputs | Updated | Status |
|---|---|---|---|---|
| `SRC_w2w/w2w`, `SRC_w2w/w2wc` | The wien2wannier interface (real / complex): projections, overlaps, eigenvalues from WIEN2k vectors; AUTO projections and the one-shot HYBR pre-Wannierization | `case.vector*`, `case.inwf`, `case.nnkp` → `case.amn/.mmn/.eig`, `case-hyb.*` | 2026-10-11 | canonical |
| `SRC_w2w/w2w_hybrid` | Standalone hybrid pre-Wannierization on existing `.amn/.mmn/.eig`; named options, `--help` is the reference; `--window` adds the automatic window route (`<out>-win.*` for stock Wannier90 disentanglement); reports subspace continuity and grades a capped disentanglement; `--pocket-skip` (default on) lowers a walked ceiling under a mixed-character pocket the weakest link sits on; admits every pool state by default, `--p-floor <x>` sets a projectability floor; `--continuity-chk <seed>` reports the subspace continuity of a finished Wannier90 window-route run from `<seed>.chk` + `<seed>.mmn` | `<seed>.amn/.mmn/.eig/.nnkp` → `<out>.amn/.mmn/.eig`, `.anchors`, `.masks`, `.anatomy`, `.dis`, `-win.*` | 2026-10-06 | canonical |
| `SRC/wien2wannier` | The public driver: `--method hybrid` (default) or `--method window`; prepares, runs the backends, publishes the hand-off (with the `.dis` status at the top of its `.win`) or the window files, optionally localizes; `--mpi-bind`; `--solver` (default `cg`); `--froz-emax <eV>|auto|none` (words need `--num-wann`); after the window route localizes, reports the continuity of Wannier90's disentangled subspace (`w2w_hybrid --continuity-chk`, never fails the run) | WIEN2k case or bridge files → `<seed>-hyb.*` or `<seed>-hyb-win.*` | 2026-10-06 | canonical |
| `SRC_w2w/Makefile.hybrid` | Build `w2w_hybrid` (`make -f Makefile.hybrid`) and the OpenMP program `w2w_hybrid_acc` (`acc` target) | — | 2026-10-12 | canonical |
| `SRC/*_lapw` | Helper scripts of wien2wannier 2.0 (`init_w2w`, `write_inwf`, `write_win`, `prepare_w2wdir`, `wannier90_lapw`, ...), kept for reference; WIEN2k installs its own copies and `make install` does not copy these | — | 2016-12-13 | supporting |
| `Makefile` (`install`, `install-cli`) | Copy `w2w`, `w2wc`, `w2w_hybrid` (and `w2w_hybrid_acc` when built) and the driver, or the driver alone, into `target-dir` (default `$WIENROOT`); no script of WIEN2k is replaced | — | 2026-10-12 | canonical |

## Tools
| Tool | Purpose / when to use | Inputs → Outputs | Updated | Status |
|---|---|---|---|---|
| `tools/plot_certificate.py` | Window non-representability certificate: what any global frozen window must misassign, against the per-state mask (exact count, evaluated at every state energy) | `<out>.masks` (or LCAO `_masks.npz`) → `.pdf/.png` | 2026-10-07 | canonical |
| `tools/band_compare.py` | Score Wannier bands against the WIEN2k band path: Vitale eta / eta_max, index-paired window RMS (`eta_windows`), E_F-convention gate, Wigner-Seitz correction, `.txt` companion; `--dump-npz` exports the curves | JSON spec (`tools/band_compare_gaas.json`), `case.klist_band` + `.energy`, `<seed>_hr.dat` (+ `_wsvec.dat`) → `.png` + `.txt` | 2026-10-01 | canonical |
| `tools/band_common.py` | File readers used by `band_compare.py` (`.klist`, `.klist_band`, `.energy`, `.fermi`, `_hr.dat`, `.eig`, `.win` k-points) | — | 2026-10-01 | supporting |
| `tools/plot_anatomy.py` | Auto-trust anatomy: sorted band-averaged projectability with the derived gate and threshold | `<out>.anatomy` (or LCAO `_report.json`) → `.pdf/.png` | 2026-09-10 | canonical |

## Documentation
| Document | Purpose | Updated | Status |
|---|---|---|---|
| `README` | What this modified wien2wannier is, GPL notice, citation, quick start | 2026-10-21 | canonical |
| `REFERENCES.md` | The literature behind each source file of the hybrid pre-Wannierization, with DOIs; source comments cite its keys | 2026-10-11 | canonical |
| `CITATION.cff` / `.zenodo.json` | How to cite the release; metadata for the Zenodo archive | 2026-10-21 | canonical |
| `NEWS` | Release notes; version 3.0 at the top | 2026-10-21 | canonical |
| `doc/hybrid/00-public-wien2wannier-cli.md` | The public `wien2wannier` command | 2026-10-21 | canonical |
| `doc/hybrid/01-overview-and-quickstart.md` | Method overview, pipeline stages, output files, reading the diagnostics | 2026-10-21 | canonical |
| `doc/hybrid/02-calculation-types.md` | Insulators, metals, SOC, spin-polarized cases | 2026-10-21 | canonical |
| `doc/hybrid/03-options-reference.md` | The `.inwf` keywords, def-file units, standalone options, main environment variables (`w2w_hybrid --help` lists all of them) and output files, including the `solver` keyword and `--solver` (default `cg`); §2.6 the checks that report or refuse | 2026-10-21 | canonical |
| `doc/hybrid/04-methods-and-tuning.md` | What to change when a run needs adjusting; §3.5 has the measured admission-floor settings on the LAPW cells and when to set `--p-floor` | 2026-10-21 | canonical |
| `doc/hybrid/05-worked-example.md` | GaAs end to end, incl. the CG-default model and the two-stage wannierization | 2026-10-21 | canonical |
| `doc/hybrid/06-paper-models.md` | The command behind each LAPW model of the paper, the Wannier90 protocol and the scoring; start here to reproduce a number, or to see what the sweep options of those commands change | 2026-10-21 | canonical |
| `doc/` (upstream) | The original wien2wannier v2.0 user guide | — | supporting |
| `INSTALL`, `WIEN-VERSION` | Build and installation with the current WIEN2k (24.1): which two programs replace WIEN2k's, what `make install` copies | 2026-10-21 | canonical |
| `docs/reference/mathematical-background.md` | Local orbitals in LAPW and the multiple-LO fix in `w2w` | 2025-12-15 | supporting |
