# Task 2: Profiling Analysis for wien2wannier

## Status: IMPLEMENTED

The following optimizations have been implemented based on this analysis.

## Overview

This document summarizes the profiling analysis of the wien2wannier codebase to identify computational bottlenecks and optimization opportunities.

## Codebase Structure

### Directory Layout

```
wien2wannier/
├── SRC/              - Shell script wrappers (non-computational)
├── SRC_w2w/          - Main w2w computation engine (Fortran) [PRIMARY TARGET]
├── SRC_wplot/        - Wannier function plotting (Fortran)
├── SRC_trig/         - Utility programs (Fortran)
└── SRC_templates/    - Template input files
```

### Main Programs

| Program | Location | Purpose |
|---------|----------|---------|
| `wf` | SRC_w2w/main.f | Main w2w computation (Amn, Mmn matrices) |
| `wf_r` | SRC_wplot/wplot.f | Wannier function plotting |
| `write_win_backend` | SRC_trig/ | Generate .win files |
| `findbands` | SRC_trig/ | Band structure analysis |

## Existing Timing Infrastructure

### ptime Module (util_w2w.F:100-282)

The codebase already has a timing infrastructure:

```fortran
type ptimer
  real(DPk) :: cpu1, cputot        ! CPU time tracking
  integer   :: wall1, walltot      ! Wall time tracking
  integer   :: countrate
end type ptimer

interface ptime
  module procedure ptime_setunit   ! Initialize timer
  module procedure ptime_print     ! Print elapsed time
  module procedure ptick           ! Start timing
end interface
```

**Usage in main.f**:
```fortran
call ptime(unit_out)       ! Initialize
call l2mmn(...)
call ptime('l2Mmn')        ! Print: "Times for l2Mmn (sec): X.XXX wall; X.XXX CPU"
call planew(...)
call ptime('planew')
call l2amn(...)
call ptime('l2Amn')
```

### Internal Profiling in l2Mmn (modules_rc.F:288-396)

The l2Mmn subroutine has detailed internal timing:

```fortran
talm   = 0.  ! Time for almgen() calls
tmeas1 = 0.  ! Time for radint() + YLM()
tmeas2 = 0.  ! Time for overlap accumulation (ZGERU loops)

call cpu_time(t1)
! ... computation ...
call cpu_time(tt0)
call almgen(Blm, ...)
call cpu_time(tt1)
! radint + YLM
call cpu_time(tt2)
! overlap loops
call cpu_time(tt3)

talm   = talm   + tt1 - tt0
tmeas1 = tmeas1 + tt2 - tt1
tmeas2 = tmeas2 + tt3 - tt2
```

**Output format**:
```
CPU time used for atom N = [total] [almgen] [radint+YLM] [overlap]
```

## Computational Hotspots

### Priority 1: l2Mmn (40-60% of runtime)

**File**: modules_rc.F:234-402

**Loop Structure** (7 levels deep):
```
atoms (Nneq)
  └─ k1 (num_kpts)           <- Largest loop, 10s-1000s iterations
     └─ k2 (NNTOT)           <- ~1-12 neighbors
        └─ mu (mult)         <- Atom multiplicities
           └─ L1, L2, LJ     <- Angular momentum (0 to LMAX2=4)
              └─ M1, MJ      <- Magnetic quantum numbers
                 └─ irf1, irf2  <- Radial functions
                    └─ ZGERU    <- BLAS rank-1 update [Nb x Nb]
```

**Computational Complexity**: O(Nk × Nntot × Nat × L^6 × Nrf^2 × Nb^2)

**Key Operations**:
- `almgen()`: DGEMM/ZGEMM matrix multiplications
- `radint()`: Radial integral computation with Bessel functions
- `ZGERU`: Rank-1 outer product updates to overlap matrix

### Priority 2: almgen (20-40% of runtime)

**File**: modules_rc.F:58-230

**Operations**:
- Bessel function computation via `harmon()`
- Spherical harmonics via `YLM()`
- BLAS matrix multiplication:
  ```fortran
  call zgemm('N', 'N', index, Nb, ibb,
             (1.0, 0), h_alyl, lda, a(ii,1,kkk), ldb,
             (1.0, 0), almt(1,1,1,mu), ldc)
  ```

**Blocking Strategy** (line 107):
```fortran
do ii=1,N-(nlo+nlon+nlov),iblock   ! iblock=128 default
   ! Process iblock elements at a time
```

### Priority 3: planew (15-30% of runtime)

**File**: modules_rc.F:668-852

**Operations**:
- G-vector loop: O(Ng^2) where Ng = (2·maxx+1)×(2·maxy+1)×(2·maxz+1)
- ZGERU rank-1 updates

**Memory Intensive**:
```fortran
FAC: 3D array [(2*maxx+1) × (2*maxy+1) × (2*maxz+1)]
A_:  [Nb × Ng × num_kpts]    ! Can be very large
B_:  [Nb × Ng]               ! Temporary
```

### Priority 4: l2Amn (10-20% of runtime)

**File**: modules_rc.F:406-665

Similar structure to l2Mmn but simpler (no k2 neighbor loop).

## BLAS/LAPACK Usage Summary

| Subroutine | Operation | Matrix Size | Call Frequency |
|------------|-----------|-------------|----------------|
| almgen | DGEMM/ZGEMM | [index × Nb] × [Nmat × Nb] | Nk × Nneq × blocking |
| l2Mmn | ZGERU | [Nb] ⊗ [Nb] = Nb^2 | Nk × Nntot × Nat × L^6 × Nrf^2 |
| planew | DGERU/ZGERU | [Nb] ⊗ [Nb] | Nk × Nntot × Ng^2 |
| l2Amn | DGEMM/ZGER | Similar to l2Mmn | Nk × Nneq × Nproj |

## Memory Analysis

### Critical Arrays

| Array | Type | Size | Location |
|-------|------|------|----------|
| `a` | complex/real | [Nmat × Nb × Nk] | xa3 module |
| `ALM`, `BLM` | complex | [Nb × Nrf × (LMAX2+1)^2 × mult] | l2Mmn local |
| `overlap` | complex | [Nb × Nb × (Nk × Nntot)] | Amn_Mmn module |
| `FAC` | complex | [(2·maxx+1)×(2·maxy+1)×(2·maxz+1)] | planew local |
| `A_` | complex | [Nb × Ng × Nk] | planew local |

### Memory Estimates

For typical parameters (Nb=100, Nk=1000, Nntot=12, Ng=1331):
- `overlap`: ~480 MB (complex)
- `A_` in planew: ~16 MB (complex)
- `ALM/BLM`: ~few MB per atom

## Parallelization Opportunities

### OpenMP Targets

1. **k-point loop in l2Mmn** (lines 303-390)
   - Independent iterations
   - Private: ALM, BLM arrays per thread
   - Reduction: overlap accumulation needs care (atomic or per-thread)

2. **k-point loop in l2Amn** (similar structure)

3. **G-vector loops in planew** (lines 774-830)
   - Highly parallel
   - May need atomic updates to overlap

### Challenges

1. **Memory scaling**: ALM/BLM must be thread-private
2. **Overlap reduction**: ZGERU updates shared overlap array
3. **I/O**: Vector file reading is sequential

## Profiling Methodology

### Step 1: Baseline Measurement

Run existing code with ptime output enabled:
```bash
./w2w case.def > profile.log 2>&1
grep "Times for" profile.log
```

### Step 2: Detailed Internal Profiling

The l2Mmn subroutine already outputs:
```
CPU time used for atom N = [total] [almgen] [radint+YLM] [overlap]
```

Parse this to get breakdown per atom.

### Step 3: External Profiling Tools

**gprof** (GNU):
```bash
# Compile with profiling
make FFLAGS="-pg -O2"
./w2w case.def
gprof w2w gmon.out > analysis.txt
```

**perf** (Linux):
```bash
perf record -g ./w2w case.def
perf report
```

**Intel VTune** (if available):
```bash
vtune -collect hotspots ./w2w case.def
vtune -report hotspots
```

### Step 4: BLAS Profiling

Check if optimized BLAS is being used:
```bash
ldd w2w | grep -i blas
# Should show MKL, OpenBLAS, or similar
```

## Recommendations

### Immediate (No Code Changes)

1. **Use optimized BLAS**: Link against MKL, OpenBLAS, or BLIS
2. **Tune iblock parameter**: Current default is 128; profile different values
3. **Enable compiler optimizations**: -O3, -march=native, etc.

### Short-term (Minor Code Changes)

1. **Add more timing points**: Profile individual loops within l2Mmn
2. **Memory pre-allocation**: Avoid repeated ALM/BLM allocation in atom loop
3. **Loop reordering**: Consider swapping loop order for better cache locality

### Medium-term (OpenMP)

1. **Parallelize k1 loop**: See Task 3 for implementation plan
2. **Thread-private arrays**: ALM, BLM per thread
3. **Reduction for overlap**: Use OpenMP reduction or atomic operations

### Long-term (Algorithm)

1. **Batched BLAS**: Replace many small ZGERU with batched operations
2. **Sparse matrix techniques**: For large, sparse overlap matrices
3. **MPI parallelization**: Distribute k-points across nodes

## Test Cases for Profiling

### Small (Baseline)
- System: Si (2 atoms, no LOs)
- Nk: ~50
- Expected: < 1 minute

### Medium (Typical)
- System: Fe with 3d LO
- Nk: ~200
- Expected: 5-15 minutes

### Large (Stress Test)
- System: SrVO3 or similar
- Nk: 1000+
- Expected: 1+ hours

## Output Format

When profiling, collect:

```
=== Profiling Run: [case name] ===
Date: YYYY-MM-DD
Parameters: Nk=[X], Nb=[Y], Nneq=[Z]

High-level timing (ptime):
  l2Mmn:  XX.XXX sec (YY.Y%)
  planew: XX.XXX sec (YY.Y%)
  l2Amn:  XX.XXX sec (YY.Y%)
  Total:  XX.XXX sec

l2Mmn breakdown per atom:
  Atom 1: total=XX.X, almgen=XX.X (YY%), radint=XX.X (YY%), overlap=XX.X (YY%)
  Atom 2: ...

Memory usage: XXX MB peak
```

## Implemented Optimizations

The following optimizations have been implemented based on this analysis:

### 1. Gaunt Coefficient Precomputation (HIGH IMPACT)

**Files Modified**: `modules.f`, `modules_rc.F`, `main.f`

**Problem**: `GAUNT1()` was called inside the innermost loop of l2Mmn, recomputing
the same coefficients millions of times.

**Solution**: Added `gaunt_cache` module that precomputes all Gaunt coefficients
at startup and stores them in a lookup table. The `get_gaunt()` function performs
O(1) lookup instead of O(N) computation.

**Expected Speedup**: 10-50% reduction in l2Mmn time (depends on LMAX).

### 2. Sparse Index Lists in planew (MEDIUM IMPACT)

**Files Modified**: `modules_rc.F`

**Problem**: Double loops over `max` indices (O(max²)) with sparse `a_null` checks
wasted most iterations on empty checks.

**Solution**: Build lists of non-null indices during initialization and iterate
only over those lists. Reduces complexity from O(max²) to O(nnz²) where nnz is
the number of non-null elements.

**Expected Speedup**: 50-90% reduction in planew t_b time for sparse G-vectors.

### 3. Memory Leak Fix in planew (BUG FIX)

**Files Modified**: `modules_rc.F`

**Problem**: Deallocate statements were unreachable (after `return`).

**Solution**: Moved deallocations before the implicit return at end of subroutine.

### 4. Index Hoisting in l2Mmn (LOW IMPACT)

**Files Modified**: `modules_rc.F`

**Problem**: `index1 = L1*(L1+1)+M1+1` was computed inside MJ loop but only
depends on L1 and M1.

**Solution**: Hoisted computation outside the MJ loop.

**Expected Speedup**: Minor (~1-5%), but essentially free.

### 5. LO Coefficient Module Update (CORRECTNESS FIX)

**Files Modified**: `modules.f`

**Problem**: Overlap integrals stored as scalars instead of per-LO arrays.

**Solution**: Added allocatable arrays `pi12lo`, `pe12lo`, `pilolo` with proper
indexing by `(l, jlo)`. Added `init_loabc()` and `cleanup_loabc()` routines.

**Note**: This is a correctness fix (Step 1 of Task 1), not a performance optimization.

## Benchmark Script

A Python benchmark script has been created at `tests/benchmark_w2w.py` to measure
the performance impact of these optimizations.

### Running on Windows (via WSL)

```bash
# Check prerequisites
wsl bash -c "cd 'path/to/wien2wannier/tests' && python3 benchmark_w2w.py --check"

# Test LO coefficient module
wsl bash -c "cd 'path/to/wien2wannier/tests' && python3 benchmark_w2w.py --test-lo"

# Build original and optimized versions
wsl bash -c "cd 'path/to/wien2wannier/tests' && python3 benchmark_w2w.py --build -v"
```

### Running on Linux/macOS

```bash
cd tests
python3 benchmark_w2w.py --check     # Check prerequisites
python3 benchmark_w2w.py --test-lo   # Test LO coefficient fix
python3 benchmark_w2w.py --build     # Build original and optimized versions
python3 benchmark_w2w.py --run       # Run benchmarks (requires WIEN2k test data)
python3 benchmark_w2w.py --compare   # Generate comparison report
python3 benchmark_w2w.py --all       # All of the above
```

### Prerequisites

- gfortran or ifort compiler
- LAPACK/BLAS libraries
- make
- Python 3.6+ with numpy, matplotlib (for plotting)

**Note**: The `--run` option requires actual WIEN2k calculation outputs (vector files,
energy files, etc.) which must be generated by running WIEN2k separately.

## References

- WIEN2k documentation: SRC_lapw1 for similar loop structures
- Wannier90 documentation: AMN/MMN matrix definitions
- BLAS reference: netlib.org/blas
- OpenMP specification: openmp.org

## Change Log

| Date | Author | Change |
|------|--------|--------|
| 2024-12-14 | Analysis | Initial profiling document |
| 2024-12-14 | Implementation | Gaunt cache, sparse planew, index hoisting, memory leak fix |
