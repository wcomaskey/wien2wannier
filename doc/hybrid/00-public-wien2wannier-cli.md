# Public `wien2wannier` command

`wien2wannier` is the supported workflow entry point for automatic hybrid
pre-Wannierization. It prepares the `.inwf`, `.def`, `.win`, and `.nnkp`
contract, selects the correct real, complex, SOC, serial, threaded-library, or
MPI backend, runs automatic projection and band selection, and publishes an
isolated `<seed>-hyb.{amn,mmn,eig,win}` handoff. The lower-level `w2w`, `w2wc`,
and `w2w_hybrid` commands remain available for compatible expert workflows.

The launcher needs Python 3.10 or newer and only the Python standard library.
NumPy and SciPy are not launcher dependencies. Install only the public command
with:

```sh
make install-cli target-dir="$HOME/bin"
```

The installed script is self-contained. Put the compiled native executables
beside it, in `$WIENROOT`, on `PATH`, or name their directory with
`--backend-dir`. The search order is `--backend-dir`, the script's own
directory, `../SRC_w2w`, `../SRC_trig`, `$WIENROOT`, then `PATH`; the first
directory that holds the program wins.

## Existing calculation contract

The command consumes a completed full-BZ WIEN2k calculation. It never launches
SCF, `kgen`, `lapw1`, or `lapwso`, and it never substitutes a new mesh for the
vectors already on disk. Its native primitive-coordinate conversion follows
the WIEN structure lattice convention before comparing the ordered NNKP points
and neighbors.

Scalar calculations require
`case.{struct,klist,vector,energy,vsp,fermi}`. Collinear calculations use the
matching `vectorup`/`energyup`/`vspup` or `vectordn`/`energydn`/`vspdn` set;
when both exist, select `--spin up` or `--spin dn`. SOC requires `case.inso`,
two vector components (`vectorsoup` and `vectorsodn`, or the supported
`vectorso` layout), the matching energy and potential files, and always uses
the complex backend. Empty `.inso` and `.in1c` placeholders are ignored during
automatic discovery. The SOC def contract contains units 9, 19 when the second
potential is present, and 22, while the first vector remains unit 10.

Supply a compatible `case.win` and `case.nnkp`. If either is missing,
`write_win_backend` creates the WIN geometry/k-point block and `wannier90.x
-pp` creates the neighbor file. The driver validates generated and reused
NNKP coordinates, order, reciprocal lattice, source blocks, neighbor indices,
and displacement shells before computation.

## Normal workflow

```sh
wien2wannier \
  --input /path/to/Si \
  --case Si \
  --seed si-model \
  --output-dir /scratch/si-model-w2w \
  --k-grid 8 8 8 \
  --threads 8
```

The default `--stage all` runs the integrated HYBR frontend and stops after the
complete isolated handoff. Add `--wannier90 /path/to/wannier90.x` to localize
that handoff. Omitted `--num-wann`, or `--num-wann 0`, asks the integrated
frontend to infer the smallest channel-complete target that contains the
frozen capacity plus its documented margin. The public default ceiling is
`--froz-emax 1.0`, in eV above `E_F`, because the automatic target is the
frozen capacity at that ceiling; the resolved plan and native log report the
policy and final target. A positive `--num-wann N` fixes the target. With it,
`--froz-emax auto` walks the ceiling down from `E_F+8` eV until the frozen
demand fits `N`, the protocol of the published LAPW models, and
`--froz-emax none` disables the ceiling; both words are refused without
`--num-wann`.

Stages have these meanings:

| Stage | Work performed |
|---|---|
| `prepare` | Validate inputs and write controls/WIN/NNKP; no projection or overlap calculation. |
| `overlaps` | Run integrated `AUTO` and require fresh raw AMN/MMN/EIG products; no hybrid solver or localization. |
| `hybrid` | Recompute needed integrated prerequisites and publish the isolated hybrid handoff. |
| `all` | Complete preparation and hybrid handoff; this is the default. |
| `localize` | Validate an existing isolated handoff and run the explicitly supplied Wannier90 executable. |

`--check` validates the entire prospective contract without creating the
output directory or running a program. `--dry-run` adds the resolved commands,
def units, scientific controls, and OMP/BLAS allocation, with the same no-write
guarantee. Existing scientific outputs cause an error unless `--force` is
given; a forced run removes the selected products before execution so a
zero-exit/no-output backend cannot republish stale data. Parent inputs and
overlapping input/output directories are never removed by `--force`.
`--quiet` suppresses routine progress while retaining successful-backend
warnings, stderr diagnostics, recognized warning/nonconvergence records from
the native `.outputwf`, and complete failure output.

## Resources

`--parallel threads --threads N` uses a one-process native frontend and sets
`OMP_NUM_THREADS` and every BLAS thread variable (`OPENBLAS_NUM_THREADS`,
`MKL_NUM_THREADS`, `BLIS_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`) to N. The
native w2w programs (built without `-fopenmp`) and Wannier90 3.x have no OpenMP
regions of their own, so `OMP_NUM_THREADS` reaches only an OpenMP-built BLAS.
Homebrew's OpenBLAS (0.3.33) is such a build: it takes its pool size from
`OMP_NUM_THREADS` and ignores `OPENBLAS_NUM_THREADS`. An OpenMP standalone bridge uses N OpenMP threads and
one BLAS thread to avoid nested pools.
`--parallel mpi --mpi-ranks R --threads N` selects `w2w_mpi`, or `w2w_mpic` for
complex and SOC vectors, launches R ranks, and gives each rank the same N-thread
settings. A WIEN2k installation names the complex MPI program `w2wc_mpi`; that
name is accepted in each search directory after `w2w_mpic`. Scheduler per-task,
task-count, and node allocations are checked unless `--allow-oversubscribe` is
explicit. MPI resources apply to the integrated frontend. NNKP preprocessing
and final Wannier90 localization each run as one process.

`--mpi-bind` sets process binding for MPI runs. With more than one thread per
rank, the default `auto` runs the MPI launcher once with `--version`, also
under `--check` and `--dry-run`; that call neither writes nor computes. If the
launcher reports Open MPI, `auto` adds `--bind-to none`. Without it, Open MPI's
documented default binds each rank to one core when it starts two ranks or
fewer, and the rank's threads would share that core. Other launchers keep their own binding. With one
thread per rank nothing is added and the launcher is not called. `--mpi-bind
none` always adds `--bind-to none` (Open MPI syntax); `--mpi-bind launcher`
never adds a binding option. The dry run reports the choice as `MPI binding:`.
Binding changes run time only: on GaAs with Open MPI 5.0.9 the products of all
three settings were byte-identical to a serial run. That machine was macOS,
where Open MPI cannot bind processes at all, so the effect on Linux was not
measured there.

The thread count can change the `.mmn` at round-off. With Homebrew OpenBLAS,
the Rh2B spin-orbit `.mmn` differed by up to 3.7e-16 between one thread and two
or more, while GaAs and MgB2 were identical; every count from 2 to 12 gave
identical bytes (2026-09-25). Compare products byte for byte only between runs
that both use one thread or both use two or more.

Resource selection never changes the solver or update order. `--solver cg`
(Grassmann conjugate gradients, stopping on the projected gradient) is the
default since 2026-10; `--solver gs` (the Gauss-Seidel sweep, the default
before) and `--solver auto` are explicit alternatives.
Colored GS is enabled only by `--colored-gs`.

## The window route (`--method window`)

The second automatic method keeps the selection and hands Wannier90 energy
windows instead of an isolated hand-off. `wien2wannier --method window
--num-wann N ...` runs the raw AUTO export (projections, overlaps,
eigenvalues), then `w2w_hybrid --window` on it: the frozen states' energy hull
becomes `dis_froz_min/max` (lowered where a k-point would hold more than N
frozen states), the outer window keeps the pool bottom and takes the ceiling
that minimizes the frozen-constrained Omega_I, and the anchor columns of the
seed projections become the `.amn`. The driver completes the window snippet
with the seed's geometry and k-points and publishes
`<seed>-hyb-win.{amn,mmn,eig,win}`; `--wannier90` then runs stock Wannier90
disentanglement and localization on `<seed>-hyb-win`. Both methods derive
their thresholds from the same projectability distribution; the hand-off
`<seed>-hyb.*` and its masks/anatomy are written as well. Bridge inputs work
the same way (`--input-kind bridge --method window`). `--num-wann N` is
required.

## Advanced routes

`--hybrid-option KEY=VALUE` forwards the validated integrated HYBR keywords
listed by `wien2wannier --help`. Dedicated options own `froz_emax`,
`pool_factor`, `solver`, and `colored`. Explicit `edges`, `seedk`, and
`seedloss` options take precedence over the matching inherited `W2W_*`
variable. When no matching option is given, native environment controls remain
available. Input paths in `W2W_STRUCT`, `W2W_PAMN`, and `W2W_CGLOAD` must name
existing nonempty files; relative values are resolved from the invocation
directory before the backend changes into `--output-dir`. Backend output paths
such as `W2W_CGDUMP` remain relative to the output directory.

The advanced bridge route starts from `<bridge-seed>.eig/.nnkp/.win` and AMN
and MMN text, native `.amn.bin/.mmn.bin` sidecars, or both:

```sh
wien2wannier \
  --input /path/to/export \
  --input-kind bridge \
  --bridge-seed raw \
  --seed model \
  --output-dir /scratch/model-bridge \
  --num-wann 18
```

The driver preserves fresh binary sidecars, rejects stale local caches, and
requires a positive fixed target. The bridge has no MPI mode. Use
`--parallel threads` only when `w2w_hybrid_acc` is available. A named-option
backend accepts repeated `--bridge-arg=--option --bridge-arg=value`; its
`w2w_hybrid --help` is the complete expert reference. Recognized development
binaries with the documented 20-position interface receive the public
`pool-factor`, `solver`, and `colored-gs` controls through a compatibility
adapter. Those legacy binaries reject `--bridge-arg`, and an unknown interface
fails before any output is written. Existing direct invocations and legacy
binaries are unchanged.
