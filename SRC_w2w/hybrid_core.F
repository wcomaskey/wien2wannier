!!! hybrid_core.F  --  native-Fortran core of the hybrid pre-Wannierization
!!! Stage-2 recipe (frozen-constrained Omega_I disentanglement -> rotate to
!!! projected-H eigenbasis -> AO/channel SCDM gauge).  Direct LAPACK/BLAS.
!!!
!!! Ported from lcao2wannier/hybrid.py (disentangle_frozen, rotate_to_subspace,
!!! scdm_anchors) + spread.py (omega_invariant), carrying the Python
!!! optimization lessons: eigenproblem restricted to the free dis-subblock,
!!! per-k index/Mdis hoisting, ZHERK Gram accumulation, Gauss-Seidel in place,
!!! throttled Omega_I convergence checks.  Everything is COMPLEX*16 (the .mmn
!!! overlaps are always complex); no _COMPLEX_ branching in the numerical core.
!!!
!!! References (keys of REFERENCES.md in the root directory):
!!!   [Marzari1997]              spread functional, b-vector weights
!!!   [Souza2001]                disentanglement by minimizing Omega_I (the sweep)
!!!   [Frankel1950], [Young1954] over-relaxation of the sweep
!!!   [Anderson1965], [Walker2011]  Anderson acceleration of the sweep
!!!   [Damle2015], [Damle2017], [Damle2018]  SCDM gauge by pivoted QR
!!!   [Green1952], [Schoenemann1966]  orthogonal Procrustes (gauge transport)
!!!   [Lowdin1950]               symmetric orthonormalization of a gauge
module hybrid_core_m
  use, intrinsic :: iso_fortran_env, only: dp => real64, i4 => int32, &
                                           output_unit
!$ use omp_lib, only: omp_in_parallel
  implicit none
  private
  public :: disentangle_frozen, omega_invariant, rotate_to_subspace, &
            scdm_anchors, block_scdm, lowdin_gauge, propagate_gauge, &
            propagate_fill, smooth_gauge, gauge_spread, gauge_best_k, &
            seed_tournament, gauge_spread_masked, seam_smooth, &
            abelian_seam_fix, kinfo_t, build_kinfo, choose_gauge, svd_left

  complex(dp), parameter :: CONE = (1.0_dp, 0.0_dp)
  complex(dp), parameter :: CZERO = (0.0_dp, 0.0_dp)
  ! One Gauss-Seidel k update computes its nntot neighbour products in parallel
  ! only when nd*nw*npool reaches this many complex multiply-adds; below it the
  ! fork costs more than the products (MgB2 with 17 functions: 19*17*19 = 6137).
  integer, parameter :: GS_PAR_MIN = 200000
  ! per-k hoisted invariants for the disentanglement sweep (masks never change)
  type :: kinfo_t
     integer :: nfroz = 0, ndis = 0, ntake = 0
     integer, allocatable :: froz(:)          ! frozen row indices (1..npool)
     integer, allocatable :: dis(:)           ! free (admit & !frozen) rows
     complex(dp), allocatable :: Mdis(:,:,:)  ! (ndis, npool, nntot) dis-rows of M(:,:,b,k)
  end type kinfo_t

contains

  !-------------------------------------------------------------------- helpers
  subroutine hermitize(A, n)   ! A <- 0.5*(A + A^H), in place
    integer, intent(in) :: n
    complex(dp), intent(inout) :: A(n,n)
    integer :: i, j
    complex(dp) :: t
    do j = 1, n
       A(j,j) = cmplx(real(A(j,j),dp), 0.0_dp, dp)
       do i = j+1, n
          t = 0.5_dp*(A(i,j) + conjg(A(j,i)))
          A(i,j) = t
          A(j,i) = conjg(t)
       end do
    end do
  end subroutine hermitize

  ! full Hermitian eigendecomposition (ascending); A overwritten by eigenvectors.
  ! Uses a workspace query (optimal LWORK, matching numpy's zheevd path) so the
  ! result is bit-for-bit reproducible against the Python reference.  The
  ! allocate is not a scaling bottleneck (the kernel is memory/dispatch bound).
  subroutine eigh(A, n, wev)
    integer, intent(in) :: n
    complex(dp), intent(inout) :: A(n,n)
    real(dp), intent(out) :: wev(n)
    complex(dp), allocatable :: work(:)
    real(dp), allocatable :: rwork(:)
    integer, allocatable :: iwork(:)
    complex(dp) :: qw(1)
    real(dp) :: qr(1)
    integer :: qi(1), lwork, lrwork, liwork, info
    if (n == 0) return
    if (n == 1) then
       wev(1) = real(A(1,1),dp); A(1,1) = CONE; return
    end if
    call zheevd('V','U', n, A, n, wev, qw,-1, qr,-1, qi,-1, info)
    lwork = int(real(qw(1),dp)); lrwork = int(qr(1)); liwork = qi(1)
    allocate(work(lwork), rwork(lrwork), iwork(liwork))
    call zheevd('V','U', n, A, n, wev, work,lwork, rwork,lrwork, iwork,liwork, info)
    if (info /= 0) then
       write(*,*) 'zheevd INFO=', info, ' n=', n
       error stop 'eigh failed'
    end if
    deallocate(work, rwork, iwork)
  end subroutine eigh

  ! smallest singular value of an n x n complex matrix (A destroyed)
  function smallest_sv(A, n) result(smin)
    integer, intent(in) :: n
    complex(dp), intent(inout) :: A(n,n)
    real(dp) :: smin
    real(dp), allocatable :: s(:), rwork(:)
    integer, allocatable :: iwork(:)
    complex(dp), allocatable :: work(:)
    complex(dp) :: qw(1), udum(1,1), vtdum(1,1)
    integer :: lwork, info
    allocate(s(n), rwork(7*n), iwork(8*n))
    call zgesdd('N', n, n, A, n, s, udum,1, vtdum,1, qw,-1, rwork, iwork, info)
    lwork = int(real(qw(1),dp)); allocate(work(lwork))
    call zgesdd('N', n, n, A, n, s, udum,1, vtdum,1, work,lwork, rwork, iwork, info)
    if (info /= 0) error stop 'zgesdd (smallest_sv) failed'
    smin = s(n)
    deallocate(s, rwork, iwork, work)
  end function smallest_sv

  ! leading `na` left singular vectors of an m x n matrix A -> U(m,na); A destroyed
  subroutine svd_left(A, m, n, na, U)
    integer, intent(in) :: m, n, na
    complex(dp), intent(inout) :: A(m,n)
    complex(dp), intent(out) :: U(m,na)
    integer :: mn
    real(dp), allocatable :: s(:), rwork(:)
    integer, allocatable :: iwork(:)
    complex(dp), allocatable :: work(:), Uf(:,:), VTf(:,:)
    complex(dp) :: qw(1)
    integer :: lwork, info
    mn = min(m,n)
    allocate(s(mn), Uf(m,mn), VTf(mn,n), iwork(8*mn))
    ! JOBZ='S': thin U(m,mn), VT(mn,n); RWORK per LAPACK 'S' rule
    allocate(rwork(mn*max(5*mn+7, 2*max(m,n)+2*mn+1)))
    call zgesdd('S', m, n, A, m, s, Uf,m, VTf,mn, qw,-1, rwork, iwork, info)
    lwork = int(real(qw(1),dp)); allocate(work(lwork))
    call zgesdd('S', m, n, A, m, s, Uf,m, VTf,mn, work,lwork, rwork, iwork, info)
    if (info /= 0) error stop 'zgesdd (svd_left) failed'
    U(:, 1:na) = Uf(:, 1:na)
    deallocate(s, Uf, VTf, iwork, rwork, work)
  end subroutine svd_left

  ! QR with column pivoting; returns the leading `ntop` pivot columns (1-based)
  subroutine qrcp_pivots(A, m, n, ntop, piv_top)
    integer, intent(in) :: m, n, ntop
    complex(dp), intent(inout) :: A(m,n)     ! destroyed
    integer, intent(out) :: piv_top(ntop)
    integer, allocatable :: jpvt(:)
    complex(dp), allocatable :: tau(:), work(:)
    real(dp), allocatable :: rwork(:)
    complex(dp) :: qw(1)
    integer :: lwork, info
    allocate(jpvt(n), tau(min(m,n)), rwork(2*n))
    jpvt = 0                                  ! all columns free
    call zgeqp3(m, n, A, m, jpvt, tau, qw,-1, rwork, info)
    lwork = int(real(qw(1),dp)); allocate(work(lwork))
    call zgeqp3(m, n, A, m, jpvt, tau, work,lwork, rwork, info)
    if (info /= 0) error stop 'zgeqp3 failed'
    piv_top = jpvt(1:ntop)
    deallocate(jpvt, tau, rwork, work)
  end subroutine qrcp_pivots

  !---------------------------------------------------------------- build_kinfo
  ! Precompute the per-k frozen/dis index sets and hoist the dis-row overlaps.
  subroutine build_kinfo(mmn, neigh, frozen, admit, npool, nntot, nk, nw, ki, &
                         want_mdis)
    integer, intent(in) :: npool, nntot, nk, nw
    complex(dp), intent(in) :: mmn(npool,npool,nntot,nk)
    integer, intent(in) :: neigh(nntot,nk)
    logical, intent(in) :: frozen(npool,nk), admit(npool,nk)
    type(kinfo_t), intent(out) :: ki(nk)
    ! The Mdis hoist is (ndis, npool, nntot) per k -- 2.6 GB at Rh3Si4 scale --
    ! and ONLY the Gauss-Seidel sweep reads it.  The Grassmann solver builds its
    ! own dis-by-dis Mdd, so it asks for the index sets alone.
    logical, intent(in), optional :: want_mdis
    logical :: do_mdis
    integer :: k, m, b, a, nf, nd
    do_mdis = .true.
    if (present(want_mdis)) do_mdis = want_mdis
    do k = 1, nk
       nf = count(frozen(:,k))
       nd = count(admit(:,k) .and. .not. frozen(:,k))
       if (nf > nw) error stop 'build_kinfo: nfroz > nw (infeasible)'
       if (nf + nd < nw) error stop 'build_kinfo: nfroz+ndis < nw (infeasible)'
       ki(k)%nfroz = nf; ki(k)%ndis = nd; ki(k)%ntake = nw - nf
       allocate(ki(k)%froz(nf), ki(k)%dis(nd))
       a = 0
       do m = 1, npool
          if (frozen(m,k)) then
             a = a + 1; ki(k)%froz(a) = m
          end if
       end do
       a = 0
       do m = 1, npool
          if (admit(m,k) .and. .not. frozen(m,k)) then
             a = a + 1; ki(k)%dis(a) = m
          end if
       end do
       ! hoist dis rows of every neighbour block:  Mdis(a,n,b) = M(dis(a),n,b,k)
       if (do_mdis) then
          allocate(ki(k)%Mdis(nd, npool, nntot))
          do b = 1, nntot
             do a = 1, nd
                ki(k)%Mdis(a,:,b) = mmn(ki(k)%dis(a), :, b, k)
             end do
          end do
       end if
    end do
  end subroutine build_kinfo

  !----------------------------------------------------------- graph coloring
  ! Greedy coloring of the (symmetric) k-neighbour graph so that no two
  ! same-colour k are neighbours -> within one colour every V(:,:,k) update is
  ! independent (reads only other-colour neighbours), which is exactly the
  ! Gauss-Seidel ordering restricted to an independent set (Jacobi diverges).
  subroutine build_coloring(neigh, nntot, nk, ncol, korder, cstart)
    integer, intent(in) :: nntot, nk, neigh(nntot,nk)
    integer, intent(out) :: ncol
    integer, allocatable, intent(out) :: korder(:), cstart(:)
    integer, allocatable :: color(:), cnt(:)
    logical, allocatable :: used(:)
    integer :: k, b, kb, bb, c, maxc, found
    ! symmetry check: neigh(:,k) must be the full neighbour set of k
    do k = 1, nk
       do b = 1, nntot
          kb = neigh(b,k); found = 0
          do bb = 1, nntot
             if (neigh(bb,kb) == k) then; found = 1; exit; end if
          end do
          if (found == 0) error stop 'build_coloring: neighbour graph not symmetric'
       end do
    end do
    maxc = 2*nntot + 2
    allocate(color(nk), used(maxc)); color = 0
    do k = 1, nk
       used = .false.
       do b = 1, nntot
          kb = neigh(b,k)
          if (color(kb) > 0) used(color(kb)) = .true.
       end do
       c = 1; do while (c <= maxc); if (.not. used(c)) exit; c = c+1; end do
       color(k) = c
    end do
    ncol = maxval(color)
    allocate(cnt(ncol), cstart(ncol+1)); cnt = 0
    do k = 1, nk; cnt(color(k)) = cnt(color(k)) + 1; end do
    cstart(1) = 1
    do c = 1, ncol; cstart(c+1) = cstart(c) + cnt(c); end do
    allocate(korder(nk)); cnt = 0
    do k = 1, nk
       c = color(k); korder(cstart(c)+cnt(c)) = k; cnt(c) = cnt(c) + 1
    end do
    deallocate(color, used, cnt)
  end subroutine build_coloring

  !------------------------------------------------------------- omega_invariant
  ! Omega_I = (1/nk) sum_{k,b} w_b ( nw - || V_k^H M^{(k,b)} V_{k+b} ||_F^2 )
  function omega_k(k, V, mmn, neigh, w, npool, nntot, nk, nw) result(ok)
    integer, intent(in) :: k, npool, nntot, nk, nw
    complex(dp), intent(in) :: V(npool,nw,nk), mmn(npool,npool,nntot,nk)
    integer, intent(in) :: neigh(nntot,nk)
    real(dp), intent(in) :: w(nntot)
    real(dp) :: ok
    ! heap scratch (see sweep_one_k note: macOS thread-stack safety)
    complex(dp), allocatable :: MV(:,:), Mt(:,:)
    integer :: b, kb
    allocate(MV(npool,nw), Mt(nw,nw))
    ok = 0.0_dp
    do b = 1, nntot
       kb = neigh(b,k)
       call zgemm('N','N', npool, nw, npool, CONE, mmn(1,1,b,k), npool, &
                  V(1,1,kb), npool, CZERO, MV, npool)
       call zgemm('C','N', nw, nw, npool, CONE, V(1,1,k), npool, &
                  MV, npool, CZERO, Mt, nw)
       ok = ok + w(b)*(real(nw,dp) - sum(real(Mt,dp)**2 + aimag(Mt)**2))
    end do
  end function omega_k

  function omega_invariant(V, mmn, neigh, w, npool, nntot, nk, nw) result(om)
    integer, intent(in) :: npool, nntot, nk, nw
    complex(dp), intent(in) :: V(npool,nw,nk), mmn(npool,npool,nntot,nk)
    integer, intent(in) :: neigh(nntot,nk)
    real(dp), intent(in) :: w(nntot)
    real(dp) :: om
    real(dp), allocatable :: omk(:)
    integer :: k
    ! Per-k array then a SERIAL sum in k order, NOT reduction(+:om).  A
    ! reduction sums in a thread-count-dependent order, and this value gates
    ! the convergence test -- so with a reduction, changing OMP_NUM_THREADS (or
    ! merely compiling with -fopenmp) can change the stopping sweep and hence
    ! every hand-off md5.  A compile flag must not change the answer.  The
    ! parallel loop below only fills omk; the accumulation order is fixed.
    allocate(omk(nk))
    !$omp parallel do default(shared) private(k) schedule(guided)
    do k = 1, nk
       omk(k) = omega_k(k, V, mmn, neigh, w, npool, nntot, nk, nw)
    end do
    !$omp end parallel do
    om = 0.0_dp
    do k = 1, nk
       om = om + omk(k)
    end do
    deallocate(omk)
    om = om/real(nk,dp)
  end function omega_invariant

  !------------------------------------------------------ one Gauss-Seidel update
  ! Update V(:,:,kk) from its neighbours: Zd = sum_b w_b (Mdis_b V_{k+b})(..)^H
  ! on the free dis-subblock; top ntake eigenvectors -> free columns; frozen
  ! rows pinned as unit columns.  Scratch is automatic => thread-private.
  subroutine sweep_one_k(kk, kinf, V, neigh, w, npool, nntot, nk, nw, esum, gap, &
                         mix, first, Zh)
    integer, intent(in) :: kk, npool, nntot, nk, nw
    type(kinfo_t), intent(in) :: kinf
    complex(dp), intent(inout) :: V(npool,nw,nk)
    integer, intent(in) :: neigh(nntot,nk)
    real(dp), intent(in) :: w(nntot)
    ! Optional: this k's free-block overlap sum = sum of the top-ntake
    ! eigenvalues of Zd (the quantity the sweep maximizes).  Costs nothing —
    ! evd is already computed.  Used by the adaptive convergence monitor.
    real(dp), intent(out), optional :: esum
    real(dp), intent(out), optional :: gap    ! M1: lambda_nt - lambda_nt+1
    ! Optional Zd mixing (W90 dis_mix_ratio analog, but on the Gauss-Seidel
    ! Gram matrix): Zd_eff = mix*Zd_new + (1-mix)*Zd_eff_prev, history in Zh.
    ! mix < 1 damps, mix > 1 extrapolates (SOR-like).  Hermiticity is
    ! preserved for any real mix; eigh only reads the 'U' triangle zherk
    ! fills.  Zh(:,:) is this k's slab, owner-written only (colored-safe).
    real(dp), intent(in), optional :: mix
    logical, intent(in), optional :: first
    complex(dp), intent(inout), optional :: Zh(:,:)
    ! Heap scratch (allocatable, NOT automatic): the automatic versions
    ! (~1-1.5 MB at npool=180-240) overflow macOS's 512 kB non-main thread
    ! stacks the moment the colored OpenMP path runs.  Heap allocation is
    ! thread-private and safe; the malloc cost is negligible vs the BLAS.
    complex(dp), allocatable :: Gb(:,:,:), Zd(:,:)
    real(dp), allocatable :: evd(:)
    integer :: b, c, a, nd, nt
    logical :: par
    nt = kinf%ntake
    if (present(esum)) esum = 0.0_dp
    if (nt == 0) return
    nd = kinf%ndis
    allocate(Gb(npool,nw,nntot), Zd(npool,npool), evd(npool))
    ! The nntot neighbour products Mdis_b V(k+b) are independent: each zgemm
    ! reads a neighbour's V, which this update does not write, and fills its
    ! own slice of Gb.  They run as an OpenMP loop over b when the products
    ! are large enough to pay for the fork and this update is not already
    ! inside the coloured sweep's parallel region (that one is parallel over
    ! k).  The accumulation below keeps its b order, so Zd, and with it every
    ! hand-off file, is the same at any thread count.
    par = .false.
!$  par = real(nd,dp)*real(nw,dp)*real(npool,dp) >= real(GS_PAR_MIN,dp) &
!$        .and. .not. omp_in_parallel()
    !$omp parallel do default(shared) private(b) schedule(static) if(par)
    do b = 1, nntot
       call zgemm('N','N', nd, nw, npool, CONE, kinf%Mdis(1,1,b), nd, &
                  V(1,1,neigh(b,kk)), npool, CZERO, Gb(1,1,b), npool)
    end do
    !$omp end parallel do
    do b = 1, nntot
       if (b == 1) then
          call zherk('U','N', nd, nw, w(b), Gb(1,1,b), npool, 0.0_dp, Zd, npool)
       else
          call zherk('U','N', nd, nw, w(b), Gb(1,1,b), npool, 1.0_dp, Zd, npool)
       end if
    end do
    if (present(mix)) then
       if (mix /= 1.0_dp) then
          if (first) then
             Zh(1:nd,1:nd) = Zd(1:nd,1:nd)
          else
             Zd(1:nd,1:nd) = mix*Zd(1:nd,1:nd) + (1.0_dp-mix)*Zh(1:nd,1:nd)
             Zh(1:nd,1:nd) = Zd(1:nd,1:nd)
          end if
       end if
    end if
    call eigh(Zd(1:nd,1:nd), nd, evd(1:nd))
    V(:,:,kk) = CZERO
    do c = 1, kinf%nfroz
       V(kinf%froz(c), c, kk) = CONE
    end do
    do c = 1, nt                       ! largest eigenvalue first (descending)
       do a = 1, nd
          V(kinf%dis(a), kinf%nfroz+c, kk) = Zd(a, nd-c+1)
       end do
    end do
    if (present(esum)) esum = sum(evd(nd-nt+1:nd))
    ! M1 -- truncation gap.  evd is sorted ascending and the top nt are kept,
    ! so this is the margin between the last retained and first discarded
    ! state: the subspace stability margin at this k.  When it is small the
    ! two manifolds are near-degenerate and a finite relaxation step can swap
    ! states ACROSS the cut, which is a basin change.  Free -- evd is already
    ! computed and otherwise only summed.
    if (present(gap)) then
       if (nd > nt) then
          gap = evd(nd-nt+1) - evd(nd-nt)
       else
          gap = huge(1.0_dp)             ! nothing discarded: no cut to cross
       end if
    end if
  end subroutine sweep_one_k

  !----------------------------------------------------------- disentangle_frozen
  subroutine disentangle_frozen(mmn, amn, neigh, w, npool, nproj, nntot, nk, nw, &
                                admit, frozen, V, omega, niter, tol, rel_tol, &
                                check_every, verbose, colored, adaptive, &
                                mixing, anderson, init_transport, kick, taper, &
                                relax_late, taper_late, nsweeps, is_converged, cap_verdict)
    integer, intent(in) :: npool, nproj, nntot, nk, nw
    complex(dp), intent(in) :: mmn(npool,npool,nntot,nk), amn(npool,nproj,nk)
    integer, intent(in) :: neigh(nntot,nk)
    real(dp), intent(in) :: w(nntot)
    logical, intent(in) :: admit(npool,nk), frozen(npool,nk)
    complex(dp), intent(out) :: V(npool,nw,nk)
    real(dp), intent(out) :: omega
    integer, intent(in) :: niter, check_every
    real(dp), intent(in) :: tol, rel_tol
    logical, intent(in) :: verbose
    logical, intent(in), optional :: colored
    logical, intent(in), optional :: adaptive
    ! transport init: BFS parallel-transport of the seed subspace across the
    ! k-mesh (child takes the nt-dim non-frozen subspace of maximal overlap
    ! with its parent through the mmn link).  Cures the k-incoherent SVD
    ! seed when projections cannot rank near-degenerate NFE states.
    logical, intent(in), optional :: init_transport
    ! Zd mixing factor (default 1 = off): <1 damps, >1 extrapolates (A5);
    ! negative = scheduled (F2).
    real(dp), intent(in), optional :: mixing
    ! Restart the schedule instead of letting the taper decay alpha to 1 and
    ! finish the run unrelaxed.  Only meaningful when mixing < 0.
    logical, intent(in), optional :: kick
    ! Alpha to restart TO when the taper spends itself.  Defaults to the
    ! starting value; may be set higher because the basin is already fixed
    ! by then (see relate).
    real(dp), intent(in), optional :: relax_late
    ! Geometric down-step of the schedule: alpha <- 1 + taper*(alpha-1) after
    ! each failed check.  Default 0.5 (halving).  The step size is a coarse
    ! GRID in alpha-space: from 1.9 the halving visits 1.9, 1.45, 1.225, ...
    ! and never samples the 1.7 neighbourhood at all, so two starting points
    ! traverse disjoint sequences.  A gentler taper (0.8-0.9) sweeps alpha
    ! CONTINUOUSLY down through the intermediate values instead of jumping
    ! over them.
    real(dp), intent(in), optional :: taper
    ! Taper used AFTER the first restart.  The opening taper must be fast to
    ! bring a super-critical start under alpha=2 (see the stability guard);
    ! the tail has no such constraint, because the restart target is already
    ! below 2.  Keeping them equal is what gives the kick sawtooth its poor
    ! duty cycle: measured on GaAs, the kick phase runs at MEAN alpha 1.236
    ! with alpha < 1.1 half the time, despite restarting to 1.95 each cycle.
    real(dp), intent(in), optional :: taper_late
    ! Safeguarded alternating Anderson acceleration (F4): history depth m
    ! (0/absent = off).  Every AA_PERIOD sweeps, extrapolate the free-block
    ! subspace PROJECTORS (gauge-invariant) over the last m+1 iterates and
    ! re-extract per-k subspaces; accept only if the free-block overlap sum
    ! S improves after one plain sweep, else revert (never slower than the
    ! plain iteration).
    integer, intent(in), optional :: anderson
    ! Outcome, for the caller's log and the hand-off status line: sweeps run,
    ! whether the tolerance was met, and -- when the cap stopped the run -- a
    ! graded verdict on the tail (see cap_tail_text); blank when converged.
    integer, intent(out), optional :: nsweeps
    logical, intent(out), optional :: is_converged
    character(*), intent(out), optional :: cap_verdict

    type(kinfo_t), allocatable :: ki(:)
    ! Omega_I at every exact check (sweep, value), for the cap verdict
    integer, allocatable :: h_it(:)
    real(dp), allocatable :: h_om(:)
    integer :: nh
    character(len=320) :: capv
    complex(dp), allocatable :: Aseed(:,:), Useed(:,:)
    integer, allocatable :: korder(:), cstart(:)
    integer :: k, it, c, nd, nt, a, ic, ncol, ns
    real(dp) :: om, omp, dlast, Sprev, Ssw, mixv, mix0, taperf, taperl
    real(dp), allocatable :: esum_k(:)
    complex(dp), allocatable :: Zhist(:,:,:)
    logical :: do_adapt, need_check, do_mix, do_sched, sw_first
    integer :: last_check, isw
    logical :: converged, do_color
    ! Scheduled-relaxation restart: once the taper has decayed alpha to within
    ! KICK_FLOOR of 1 the acceleration is spent, so re-kick it to its initial
    ! value rather than finish the run unrelaxed.
    real(dp), parameter :: KICK_FLOOR = 1.0e-2_dp
    logical :: do_kick
    integer :: nkick
    real(dp), allocatable :: gap_k(:)      ! M1: per-k truncation gap
    ! Auto-calibration of the late-phase alpha.  The sweeps-vs-alpha curve is
    ! the classical SOR shape -- monotone down to an optimum, then a sharp rise
    ! as the iteration turns oscillatory -- but the optimum is MATERIAL-
    ! DEPENDENT (measured: GaAs 1.99, MgB2/Fe-up/Fe-dn 1.95), so any fixed
    ! value is wrong somewhere: worst-case +168% at 1.95, +41% at 1.99.
    !
    ! EXPERIMENTAL, and currently WORSE than the fixed default -- do not enable
    ! without reading this.  The probe scores candidates by the Omega_I drop
    ! over PROBE_N sweeps, which is the wrong quantity: over-relaxation near
    ! omega_opt is oscillatory-but-fast ASYMPTOTICALLY, so a short window sees
    ! the oscillation as poor progress and rejects the winner.  Measured on
    ! GaAs, where the true optimum is 1.99 (1596 sweeps vs 4274 at 1.95):
    !   PROBE_N   50 -> picks 1.95    100 -> 1.95    200 -> 1.95    400 -> 1.99
    ! It needs a 400-sweep window to see the asymptotic advantage, i.e. 1600
    ! sweeps over four candidates -- more than the ENTIRE run on MgB2 (996) or
    ! Fe (~1030), which would already have picked correctly.  Net: wins on GaAs
    ! against a bad constant, loses badly everywhere else.
    !
    ! The fix is to score the CONTRACTION RATE rather than the total drop, but
    ! over-relaxation is non-monotone so successive decrements change sign and
    ! a naive ratio estimator is unstable; it needs an envelope fit.  Until
    ! then the fixed minimax constant is the better engineering choice.
    logical :: auto_late
    integer :: PROBE_N                     ! env W2W_PROBEN, default 50
    real(dp), parameter :: PROBE_A(4) = [1.99_dp, 1.95_dp, 1.85_dp, 1.70_dp]
    real(dp) :: gapmin, dprev, rho, a_young, relate
    ! Anderson state (F4)
    integer, parameter :: AA_PERIOD = 5
    integer :: aa_m, nhist, aa_i, aa_j, aa_naccept, aa_nreject
    complex(dp), allocatable :: Phist(:,:,:,:)     ! projector iterates ring
    complex(dp), allocatable :: Pacc(:,:,:), Vsave(:,:,:), Qd(:,:)
    real(dp), allocatable :: Gm(:,:), rhs(:), evq(:), esum2(:)
    real(dp) :: Sbefore, Safter
    logical :: do_transp

    allocate(ki(nk))
    call build_kinfo(mmn, neigh, frozen, admit, npool, nntot, nk, nw, ki)

    ! ---- initial subspace: frozen unit columns + project-optimal SVD seed ----
    V = CZERO
    do k = 1, nk
       do c = 1, ki(k)%nfroz
          V(ki(k)%froz(c), c, k) = CONE
       end do
       nt = ki(k)%ntake; nd = ki(k)%ndis
       if (nt > 0) then
          allocate(Aseed(nd,nproj), Useed(nd,nt))
          do a = 1, nd
             Aseed(a,:) = amn(ki(k)%dis(a), :, k)
          end do
          call svd_left(Aseed, nd, nproj, nt, Useed)
          do c = 1, nt
             do a = 1, nd
                V(ki(k)%dis(a), ki(k)%nfroz+c, k) = Useed(a,c)
             end do
          end do
          deallocate(Aseed, Useed)
       end if
    end do

    ! ---- optional: parallel-transport init (BFS spanning tree) ----
    do_transp = .false.
    if (present(init_transport)) do_transp = init_transport
    if (do_transp) then
       block
         integer :: qs(nk), qh, qt, pk, ck, bb, ndc, ntc, a2, c2
         logical :: vis(nk)
         complex(dp), allocatable :: Wfull(:,:), Wd(:,:), Us2(:,:)
         vis = .false.; qh = 1; qt = 1; qs(1) = 1; vis(1) = .true.
         do while (qh <= qt)
            pk = qs(qh); qh = qh + 1
            do bb = 1, nntot
               ck = neigh(bb, pk)
               if (ck < 1 .or. ck > nk) cycle
               if (vis(ck)) cycle
               vis(ck) = .true.; qt = qt + 1; qs(qt) = ck
               ntc = ki(ck)%ntake; ndc = ki(ck)%ndis
               if (ntc <= 0 .or. ndc <= 0) cycle
               allocate(Wfull(npool,nw), Wd(ndc,nw), Us2(ndc,ntc))
               ! image of the parent subspace at the child: W = M^H V_p
               call zgemm('C','N', npool, nw, npool, CONE, mmn(1,1,bb,pk), &
                          npool, V(1,1,pk), npool, CZERO, Wfull, npool)
               do a2 = 1, ndc
                  Wd(a2,:) = Wfull(ki(ck)%dis(a2), :)
               end do
               call svd_left(Wd, ndc, nw, ntc, Us2)
               do c2 = 1, ntc
                  do a2 = 1, ndc
                     V(ki(ck)%dis(a2), ki(ck)%nfroz+c2, ck) = Us2(a2,c2)
                  end do
               end do
               deallocate(Wfull, Wd, Us2)
            end do
         end do
       end block
       if (verbose) write(*,'(a)') &
            '   init: parallel-transport subspace seed (BFS tree)'
    end if

    ! ---- Gauss-Seidel sweeps (in place; neighbours read current-sweep V) ----
    do_color = .false.
    if (present(colored)) do_color = colored
    if (do_color) then
       call build_coloring(neigh, nntot, nk, ncol, korder, cstart)
       if (verbose) write(*,'(a,i0,a,i0,a)') '   red-black: ', ncol, &
            ' colours over ', nk, ' k-points (parallel sweeps)'
    end if
    omp = huge(1.0_dp)
    converged = .false.
    do_adapt = .false.
    if (present(adaptive)) do_adapt = adaptive
    if (do_adapt) allocate(esum_k(nk))
    ! mixing > 1: fixed over-relaxation.  mixing < 0: SCHEDULED relaxation
    ! (SRJ-flavoured): start at |mixing| for fast basin escape, then taper
    ! geometrically toward 1 after each failed exact convergence check --
    ! aggressive early, monotone-safe late.
    mixv = 1.0_dp
    do_sched = .false.
    if (present(mixing)) then
       mixv = mixing
       if (mixv < 0.0_dp) then
          mixv = abs(mixv)
          do_sched = .true.
       end if
    end if
    taperf = 0.5_dp
    if (present(taper)) taperf = taper
    ! alpha0 and the taper are COUPLED, not independent knobs.  A start above
    ! Kahan's bound only survives because the taper drops it back under 2 fast:
    ! from 2.5, taper=0.5 spends ONE step above 2 (2.5 -> 1.75) and converges,
    ! while taper=0.8 spends two (2.5 -> 2.2 -> 1.96) and diverges outright
    ! (measured: Omega_I 181.6 after 200 sweeps).  Refuse the combination
    ! rather than let it blow up thousands of sweeps later.
    if (mixv > 2.0_dp .and. taperf > 0.6_dp) then
       write(*,'(a,f5.2,a,f4.2,a)') '  [FATAL] relax start ', mixv, &
            ' with taper ', taperf, ' stays above the alpha=2 stability'
       write(*,'(a)') '          bound for more than one step and will diverge.'
       write(*,'(a)') '          Use taper <= 0.6 with a super-critical start.'
       error stop 'disentangle_frozen: unstable relax/taper combination'
    end if
    mix0 = mixv                        ! restart target for the schedule
    nkick = 0
    allocate(gap_k(nk)); gap_k = 0.0_dp
    gapmin = 0.0_dp; dprev = 0.0_dp; rho = 0.0_dp; a_young = 0.0_dp
    ! Late-phase relaxation.  Basin selection happens in the OPENING sweeps
    ! while alpha is near its initial value (measured: on GaAs a start of 1.85
    ! or above commits to a different, 0.13%-worse minimum no matter how the
    ! taper is shaped, and 1.80 sits on the separatrix and critically slows).
    ! Once committed, alpha is no longer basin-critical, so the restart can
    ! safely exceed the starting value -- mild early to pick the basin,
    ! aggressive late to converge inside it.
    relate = mix0
    auto_late = .false.
    PROBE_N = 50
    block
      character(len=32) :: pn
      integer :: ios2
      call get_environment_variable('W2W_PROBEN', pn)
      if (len_trim(pn) > 0) then
         read(pn, *, iostat=ios2) PROBE_N
         if (ios2 /= 0 .or. PROBE_N < 5) PROBE_N = 50
      end if
    end block
    if (present(relax_late)) then
       if (relax_late > 1.0_dp) relate = relax_late
       ! negative sentinel = calibrate at the first restart
       if (relax_late < 0.0_dp) auto_late = .true.
    end if
    taperl = taperf
    if (present(taper_late)) then
       if (taper_late > 0.0_dp .and. taper_late < 1.0_dp) taperl = taper_late
    end if
    do_kick = .false.
    if (present(kick)) do_kick = kick .and. do_sched
    do_mix = (mixv /= 1.0_dp) .or. do_sched
    if (do_mix) allocate(Zhist(npool,npool,nk))
    aa_m = 0
    if (present(anderson)) aa_m = max(0, anderson)
    if (aa_m > 0) then
       do_adapt = .true.               ! AA needs the per-sweep S signal
       if (.not. allocated(esum_k)) allocate(esum_k(nk))
       allocate(Phist(npool,npool,nk,aa_m+1), Pacc(npool,npool,nk), &
                Vsave(npool,nw,nk), Qd(npool,npool), evq(npool), &
                Gm(aa_m,aa_m), rhs(aa_m))
       nhist = 0; aa_naccept = 0; aa_nreject = 0
    end if
    ! Sweeps run in blocks; Omega_I is tested between blocks.  The first block is
    ! ONE sweep and later blocks are `check_every` sweeps, so the check points
    ! (sweeps 1, 1+ce, 1+2ce, ...) match the Python reference exactly (=> the
    ! serial path is bit-for-bit reproducible).  The colored path opens ONE
    ! parallel region per block, amortizing fork/join over the block.
    !
    ! ADAPTIVE mode (opt-in): each sweep also collects the free-block overlap
    ! sum S = sum_k (top-ntake eigenvalues of Zd) — free, since the eigen-
    ! decomposition is already done.  The EXPENSIVE exact omega_invariant is
    ! called only when S plateaus (|dS| < 10*rel_tol*|S|), at a fallback
    ! cadence of 10*check_every sweeps, or at the cap; the exact value still
    ! gates the stop, so the converged Omega_I is exact either way.
    it = 0
    ns = 1                                   ! first block: a single sweep
    dlast = huge(1.0_dp)
    Sprev = huge(1.0_dp)
    last_check = 0
    allocate(h_it(niter+1), h_om(niter+1)); nh = 0
    do while (it < niter)
       if (do_adapt) then
          ns = 1                             ! per-sweep S signal
       else
          ns = min(ns, niter-it)
       end if
       if (do_color) then
          ! isw MUST be explicitly private: a BLOCK-declared loop variable
          ! inside the parallel region was bound SHARED by gfortran-16,
          ! desynchronizing the threads' omp-do sequences (NaN + barrier
          ! spin observed on Rh3Si4 2x10x6; Fe happened to survive).
          !$omp parallel default(shared) private(ic, a, k, isw)
          do isw = 1, ns
             do ic = 1, ncol
                ! guided, not static: within one colour no k can observe
                ! another's write (see build_coloring), so thread order cannot
                ! affect the result -- this is bit-identical to static, not
                ! merely close.  Static hands contiguous chunks, and per-k cost
                ! scales with ndis(k), which is spatially correlated, so some
                ! threads got systematically heavier work and the rest waited
                ! at the colour barrier.
                !$omp do schedule(guided)
                do a = cstart(ic), cstart(ic+1)-1
                   k = korder(a)
                   ! `first` MUST be evaluated here, where isw is the
                   ! thread-private copy (host association inside sweep_k
                   ! would reach the shared original).
                   call sweep_k(k, it == 0 .and. isw == 1)
                end do
                !$omp end do
             end do
          end do
          !$omp end parallel
       else
          do isw = 1, ns
             do k = 1, nk
                call sweep_k(k, it == 0 .and. isw == 1)
             end do
          end do
       end if
       it = it + ns

       ! ---- F4: safeguarded alternating Anderson on the subspace projectors
       if (aa_m > 0) then
          if (nhist == aa_m+1) Phist(:,:,:,1:aa_m) = Phist(:,:,:,2:aa_m+1)
          nhist = min(nhist+1, aa_m+1)
          do k = 1, nk
             call zgemm('N','C', npool, npool, nw, CONE, V(1,1,k), npool, &
                        V(1,1,k), npool, CZERO, Phist(1,1,k,nhist), npool)
          end do
          if (nhist == aa_m+1 .and. mod(it, AA_PERIOD) == 0) then
             Sbefore = sum(esum_k)
             ! theta = G^-1 1 / (1^T G^-1 1) with G_ij = Re<f_i, f_j>,
             ! f_i = P_{i+1} - P_i  (minimize ||sum theta_i f_i||, sum=1)
             do aa_i = 1, aa_m
                do aa_j = aa_i, aa_m
                   Gm(aa_i,aa_j) = real(sum(conjg( &
                        Phist(:,:,:,aa_i+1)-Phist(:,:,:,aa_i)) * &
                       (Phist(:,:,:,aa_j+1)-Phist(:,:,:,aa_j))), dp)
                   Gm(aa_j,aa_i) = Gm(aa_i,aa_j)
                end do
                Gm(aa_i,aa_i) = Gm(aa_i,aa_i)*(1.0_dp + 1.0e-10_dp)  ! ridge
             end do
             rhs = 1.0_dp
             call solve_spd(Gm, rhs, aa_m)          ! rhs := G^-1 1
             rhs = rhs / sum(rhs)                   ! theta
             Pacc = CZERO
             do aa_i = 1, aa_m
                Pacc = Pacc + rhs(aa_i)*Phist(:,:,:,aa_i+1)
             end do
             Vsave = V
             do k = 1, nk
                nt = ki(k)%ntake
                if (nt == 0) cycle
                nd = ki(k)%ndis
                do a = 1, nd
                   do c = 1, nd
                      Qd(c,a) = Pacc(ki(k)%dis(c), ki(k)%dis(a), k)
                   end do
                end do
                call eigh(Qd(1:nd,1:nd), nd, evq(1:nd))
                do c = 1, nt
                   V(:, ki(k)%nfroz+c, k) = CZERO
                   do a = 1, nd
                      V(ki(k)%dis(a), ki(k)%nfroz+c, k) = Qd(a, nd-c+1)
                   end do
                end do
             end do
             ! safeguard: one plain sweep on the candidate; keep only if the
             ! ascent objective S improves.
             do k = 1, nk
                call sweep_k(k, .false.)
             end do
             Safter = sum(esum_k)
             if (Safter >= Sbefore) then
                it = it + 1
                aa_naccept = aa_naccept + 1
                nhist = 0                            ! restart history
             else
                V = Vsave
                aa_nreject = aa_nreject + 1
             end if
          end if
       end if

       if (do_adapt) then
          Ssw = sum(esum_k)
          need_check = (abs(Sprev-Ssw) < 10.0_dp*rel_tol*max(abs(Ssw),1.0_dp)) &
                       .or. (it - last_check >= 10*check_every) .or. (it >= niter)
          ! after a failed exact check, wait at least check_every sweeps
          need_check = need_check .and. (it - last_check >= check_every .or. it >= niter)
          Sprev = Ssw
          if (.not. need_check) cycle
       end if

       om = omega_invariant(V, mmn, neigh, w, npool, nntot, nk, nw)
       last_check = it
       if (nh < size(h_it)) then
          nh = nh + 1; h_it(nh) = it; h_om(nh) = om
       end if
       if (verbose) write(*,'(a,i5,a,f18.10)') '   sweep', it, '  Omega_I=', om
       dlast = abs(omp-om)
       ! M1/M2 telemetry at each exact check.
       !   gapmin  = min over k of the truncation gap (subspace safety margin)
       !   rho     = observed contraction ratio dOmega_n / dOmega_{n-1}
       !   a_young = 2/(1+sqrt(1-rho^2)), the classical SOR optimum for a
       !             linear iteration contracting at rho.  For a nonlinear
       !             Grassmann iteration it is an UPPER guide, not a
       !             prescription -- compare it against gapmin before trusting
       !             it, since that is where basin escape comes from.
       gapmin = minval(gap_k)
       if (dprev > 0.0_dp) then
          rho = min(dlast/dprev, 0.999999_dp)
          a_young = 2.0_dp/(1.0_dp + sqrt(max(1.0_dp - rho*rho, 1.0e-12_dp)))
       end if
       dprev = dlast
       if (dlast < tol .or. dlast < rel_tol*max(abs(om),1.0_dp)) then
          converged = .true.; omega = om; exit
       end if
       omp = om
       ! Scheduled relaxation: taper alpha toward 1 after each failed exact
       ! check (mixv is read inside the parallel region but only updated here,
       ! between sweep blocks).
       !
       ! RESTART (do_kick): the plain taper halves alpha-1 every failed check,
       ! so after ~10 checks alpha ~ 1.0007 and the run finishes UNRELAXED,
       ! crawling to the tolerance -- measured on the GaAs one-shot as 10042
       ! sweeps at -1.7 against 1649 at -1.9.  Re-kicking alpha back to its
       ! initial value once the taper has spent itself keeps the acceleration
       ! available for the whole run instead of only its opening.  This is
       ! self-pacing: the kick is driven by the taper having decayed, not by
       ! any per-material tuning.  Note the choice of alpha selects which
       ! minimum you land in (GaAs: -1.7 -> 21.83306, -1.9 -> 21.86151, both
       ! converged to 1e-13), so kicking back to mix0 also keeps the run in
       ! the basin its alpha chose rather than drifting as alpha -> 1.
       if (do_sched) then
          if (mixv > 1.0_dp + KICK_FLOOR) then
             ! fast taper while the opening start is being brought down;
             ! slow taper once we are in the restart cycle, to hold alpha
             ! near the restart target instead of collapsing back toward 1
             if (nkick > 0) then
                mixv = 1.0_dp + taperl*(mixv - 1.0_dp)
             else
                mixv = 1.0_dp + taperf*(mixv - 1.0_dp)
             end if
          else if (do_kick) then
             if (auto_late .and. nkick == 0) call probe_late()
             mixv = relate                 ! re-kick (late alpha may exceed mix0)
             nkick = nkick + 1
          end if
       end if
       if (.not. do_adapt) ns = check_every  ! legacy cadence
    end do
    if (.not. converged) omega = omega_invariant(V, mmn, neigh, w, npool, nntot, nk, nw)
    ! Unconditional telemetry: sweep count, convergence status, last |dOmega|.
    ! (The cap-hit case was previously silent under verbose=.false. — the
    !  Rh3Si4 production runs may have been sweep-limited without any signal.)
    write(*,'(a,i6,a,l1,a,es10.2,a,i6)') '  disentangle sweeps=', it, &
         '  converged=', converged, '  |dOmega|_last=', dlast, &
         '  cap=', niter
    if (aa_m > 0) write(*,'(a,i4,a,i4,a)') '  anderson: ', aa_naccept, &
         ' accepted, ', aa_nreject, ' rejected'
    if (do_kick) write(*,'(a,i0,a,f5.2,a)') '  relax restarts: ', nkick, &
         ' (alpha re-kicked to ', relate, ' when the taper spent itself)'
    write(*,'(a,es11.3,a,f6.3,a,f6.3)') '  measurables: min gap=', gapmin, &
         '  rho=', rho, '  alpha_young=', a_young
    if (allocated(gap_k)) deallocate(gap_k)
    capv = ''
    if (.not. converged) then
       capv = cap_tail_text(h_it(1:nh), h_om(1:nh), it, omega)
       write(*,'(a)') '  [WARN] disentangle hit the sweep cap WITHOUT converging: '//trim(capv)
    end if
    if (present(nsweeps)) nsweeps = it
    if (present(is_converged)) is_converged = converged
    if (present(cap_verdict)) cap_verdict = capv
    deallocate(h_it, h_om)
    deallocate(ki)

  contains

    ! One k update, dispatching the optional esum/mixing plumbing.  Shared
    ! host-associated state only (V, ki, esum_k, Zhist, do_*); everything
    ! thread-varying comes in as arguments.
    subroutine sweep_k(k1, first1)
      integer, intent(in) :: k1
      logical, intent(in) :: first1
      if (do_mix) then
         if (do_adapt) then
            call sweep_one_k(k1, ki(k1), V, neigh, w, npool, nntot, nk, nw, &
                 esum=esum_k(k1), gap=gap_k(k1), mix=mixv, first=first1, &
                 Zh=Zhist(:,:,k1))
         else
            call sweep_one_k(k1, ki(k1), V, neigh, w, npool, nntot, nk, nw, &
                 gap=gap_k(k1), mix=mixv, first=first1, Zh=Zhist(:,:,k1))
         end if
      else if (do_adapt) then
         call sweep_one_k(k1, ki(k1), V, neigh, w, npool, nntot, nk, nw, &
              esum=esum_k(k1), gap=gap_k(k1))
      else
         call sweep_one_k(k1, ki(k1), V, neigh, w, npool, nntot, nk, nw, &
              gap=gap_k(k1))
      end if
    end subroutine sweep_k

    !> Price each candidate late-alpha by running PROBE_N sweeps from the
    !> CURRENT iterate and comparing the Omega_I decrement, then keep the best.
    !> V and the mixing history are saved and restored around every trial, so
    !> the probe leaves the iteration exactly where it found it and only
    !> `relate` changes.  Uses the serial sweep order for all candidates so the
    !> comparison is like-for-like.
    subroutine probe_late()
      real(dp) :: om0, sc(size(PROBE_A))
      complex(dp), allocatable :: Vsav(:,:,:), Zsav(:,:,:)
      integer :: ic, ip, kq
      allocate(Vsav, source=V)
      if (do_mix) allocate(Zsav, source=Zhist)
      om0 = omega_invariant(V, mmn, neigh, w, npool, nntot, nk, nw)
      do ic = 1, size(PROBE_A)
         V = Vsav
         if (do_mix) Zhist = Zsav
         mixv = PROBE_A(ic)
         do ip = 1, PROBE_N
            do kq = 1, nk
               call sweep_k(kq, .false.)
            end do
         end do
         ! decrement achieved in PROBE_N sweeps; larger is better.  A candidate
         ! that diverges gives a negative (or NaN) score and loses on its own.
         sc(ic) = om0 - omega_invariant(V, mmn, neigh, w, npool, nntot, nk, nw)
         if (sc(ic) /= sc(ic)) sc(ic) = -huge(1.0_dp)      ! NaN guard
      end do
      V = Vsav
      if (do_mix) Zhist = Zsav
      relate = PROBE_A(maxloc(sc, 1))
      write(*,'(a,i0,a)') '  late-alpha probe (', PROBE_N, ' sweeps each):'
      do ic = 1, size(PROBE_A)
         write(*,'(a,f5.2,a,es12.4,a)') '    alpha=', PROBE_A(ic), &
              '  dOmega=', sc(ic), merge('<- kept ', '        ', &
              PROBE_A(ic) == relate)   ! both branches 8 chars
      end do
      deallocate(Vsav)
      if (allocated(Zsav)) deallocate(Zsav)
    end subroutine probe_late


  end subroutine disentangle_frozen

  !> Grades a disentanglement that stopped at its sweep cap.  Compares Omega_I at the
  !> end with its value at the last exact check at or before one tenth of the run
  !> earlier: a relative fall above CAP_TAIL_REL (1e-5) means the minimization was
  !> still descending when the cap cut it off; below it the tail is flat and the
  !> result is usable although the tolerance was not met.  The same rule and wording
  !> as lcao2wannier's hand-off (hybrid.cap_tail_verdict).
  function cap_tail_text(h_it, h_om, it_end, om_end) result(txt)
    integer,  intent(in) :: h_it(:), it_end
    real(dp), intent(in) :: h_om(:), om_end
    character(len=320) :: txt
    real(dp), parameter :: CAP_TAIL_REL = 1.0e-5_dp
    integer :: j, jref, tgt
    real(dp) :: rel
    txt = 'too short a history to grade the cap'
    if (size(h_it) < 1) return
    tgt = it_end - max(1, it_end/10)
    jref = 1
    do j = 1, size(h_it)
       if (h_it(j) <= tgt) jref = j
    end do
    if (h_it(jref) >= it_end) return
    rel = (h_om(jref) - om_end)/max(abs(om_end), 1.0e-300_dp)
    if (rel > CAP_TAIL_REL) then
       write(txt,'(a,es8.1,a,i0,a)') 'POOR CONVERGENCE: still descending (Omega_I fell ', rel, &
            ' relative over the last ', it_end - h_it(jref), ' sweeps); raise the sweep cap '// &
            '(w2w_hybrid --sweeps N; HYBR keyword sweeps = N) or use the CG solver '// &
            '(--solver cg; solver = cg)'
    else
       write(txt,'(a,es8.1,a,i0,a)') 'flat tail (Omega_I moved ', rel, ' relative over the last ', &
            it_end - h_it(jref), ' sweeps); usable, tolerance not met'
    end if
  end function cap_tail_text

  !---------------------------------------------------------- rotate_to_subspace
  ! diagonalise the projected Hamiltonian; rotate M into that eigenbasis.
  subroutine rotate_to_subspace(eig, V, mmn, neigh, npool, nntot, nk, nw, &
                                eig_t, W, mmn_rot)
    integer, intent(in) :: npool, nntot, nk, nw
    real(dp), intent(in) :: eig(npool,nk)
    complex(dp), intent(in) :: V(npool,nw,nk), mmn(npool,npool,nntot,nk)
    integer, intent(in) :: neigh(nntot,nk)
    real(dp), intent(out) :: eig_t(nw,nk)
    complex(dp), intent(out) :: W(npool,nw,nk), mmn_rot(nw,nw,nntot,nk)
    complex(dp), allocatable :: tmp(:,:), Hsub(:,:), R(:,:), MW(:,:)
    real(dp), allocatable :: ev(:)
    integer :: k, b, m, c, kb
    allocate(tmp(npool,nw), Hsub(nw,nw), R(nw,nw), MW(npool,nw), ev(nw))
    do k = 1, nk
       ! Hsub = V_k^H diag(eig) V_k
       do c = 1, nw
          do m = 1, npool
             tmp(m,c) = eig(m,k)*V(m,c,k)
          end do
       end do
       call zgemm('C','N', nw, nw, npool, CONE, V(1,1,k), npool, tmp, npool, &
                  CZERO, Hsub, nw)
       call hermitize(Hsub, nw)
       R = Hsub
       call eigh(R, nw, ev)              ! R <- eigenvectors, ev ascending
       eig_t(:,k) = ev
       ! W_k = V_k R
       call zgemm('N','N', npool, nw, nw, CONE, V(1,1,k), npool, R, nw, &
                  CZERO, W(1,1,k), npool)
    end do
    ! mmn_rot(:,:,b,k) = W_k^H M^{(k,b)} W_{k+b}
    do k = 1, nk
       do b = 1, nntot
          kb = neigh(b,k)
          call zgemm('N','N', npool, nw, npool, CONE, mmn(1,1,b,k), npool, &
                     W(1,1,kb), npool, CZERO, MW, npool)
          call zgemm('C','N', nw, nw, npool, CONE, W(1,1,k), npool, MW, npool, &
                     CZERO, mmn_rot(1,1,b,k), nw)
       end do
    end do
    deallocate(tmp, Hsub, R, MW, ev)
  end subroutine rotate_to_subspace

  !------------------------------------------------------------ propagate_gauge
  ! Parallel-transport / Procrustes propagation of the anchor gauge across
  ! the BZ (smooth-frame construction in the spirit of Cances, Levitt,
  ! Panati & Stoltz, PRB 95, 075114): BFS over the k-neighbour graph from
  ! the anchor k0; for each tree edge k'(visited) -> k = neigh(b,k'),
  !   U_k = polar( M(k',b)^H U_k' ),   polar(X) = u v^H from the SVD,
  ! the closest unitary to the transported frame.  The result is exactly
  ! unitary at every k; its virtue over the raw SCDM columns is gauge
  ! SMOOTHNESS far from the anchor -- the single-anchor degradation that
  ! worst-k sigma_min detects is what the transport repairs.
  subroutine propagate_gauge(mmn, neigh, nw, nntot, nk, k0, A_gauge)
    integer, intent(in) :: nw, nntot, nk, k0
    complex(dp), intent(in) :: mmn(nw,nw,nntot,nk)
    integer, intent(in) :: neigh(nntot,nk)
    complex(dp), intent(inout) :: A_gauge(nw,nw,nk)  ! in: valid gauge at k0
    integer, allocatable :: queue(:)
    logical, allocatable :: seen(:)
    complex(dp), allocatable :: X(:,:), Us(:,:), VTs(:,:), work(:)
    real(dp), allocatable :: sv(:), rwork(:)
    complex(dp) :: qw(1)
    integer :: qh, qt, kp, k, b, info, lwork
    allocate(queue(nk), seen(nk), X(nw,nw), Us(nw,nw), VTs(nw,nw), sv(nw), &
             rwork(5*nw))
    call zgesvd('A','A', nw, nw, X, nw, sv, Us, nw, VTs, nw, qw, -1, rwork, info)
    lwork = int(real(qw(1),dp)); allocate(work(lwork))
    seen = .false.; seen(k0) = .true.
    queue(1) = k0; qh = 1; qt = 1
    do while (qh <= qt)
       kp = queue(qh); qh = qh + 1
       do b = 1, nntot
          k = neigh(b,kp)
          if (seen(k)) cycle
          call zgemm('C','N', nw, nw, nw, CONE, mmn(1,1,b,kp), nw, &
                     A_gauge(1,1,kp), nw, CZERO, X, nw)
          call zgesvd('A','A', nw, nw, X, nw, sv, Us, nw, VTs, nw, &
                      work, lwork, rwork, info)
          if (info /= 0) error stop 'propagate_gauge: SVD failed'
          call zgemm('N','N', nw, nw, nw, CONE, Us, nw, VTs, nw, CZERO, &
                     A_gauge(1,1,k), nw)
          seen(k) = .true.
          qt = qt + 1; queue(qt) = k
       end do
    end do
    if (.not. all(seen)) error stop 'propagate_gauge: k-graph not connected'
  end subroutine propagate_gauge

  !------------------------------------------------------------ align_gauge
  ! BFS Procrustes ALIGNMENT.  Unlike propagate_gauge, which REPLACES each
  ! child frame with the transported parent (A_k <- polar(M^H A_p), so every
  ! per-k frame is thrown away), this KEEPS each k's own frame and only
  ! rotates it: A_k <- A_k W with W = polar(A_k^H M^H A_p), the unitary that
  ! best matches the parent's frame carried to k.  Right-multiplication
  ! preserves the column SPAN exactly (Omega_I untouched) and keeps the
  ! frame built from that k's own vectors, while removing the arbitrary
  ! relative orientation that independently-solved per-k frames carry.
  !
  ! MEASURED MOTIVATION (MgB2 12x12x8):
  ! per-k Loewdin frames start at half auto8's spread (34.6 vs 69.1) with
  ! 5x better off-mesh symmetry, but need 2322 W90 iterations vs 437.
  ! Transport-filling the 14 rank-deficient k bought only 21% (-> 1826), so
  ! the incoherence is DISTRIBUTED over all 1152 k, not concentrated in the
  ! sick ones; and smoothing is Jacobi-like (information travels one
  ! neighbour per sweep), so 8 sweeps cannot align a mesh of diameter 12.
  ! A spanning-tree alignment fixes global orientation in ONE pass.
  subroutine align_gauge(mmn, neigh, nw, nntot, nk, k0, A_gauge)
    integer, intent(in) :: nw, nntot, nk, k0
    complex(dp), intent(in) :: mmn(nw,nw,nntot,nk)
    integer, intent(in) :: neigh(nntot,nk)
    complex(dp), intent(inout) :: A_gauge(nw,nw,nk)
    integer, allocatable :: queue(:)
    logical, allocatable :: seen(:)
    complex(dp), allocatable :: T(:,:), O(:,:), W(:,:), Us(:,:), VTs(:,:), work(:)
    real(dp), allocatable :: sv(:), rwork(:)
    complex(dp) :: qw(1)
    integer :: qh, qt, kp, k, b, info, lwork
    allocate(queue(nk), seen(nk), T(nw,nw), O(nw,nw), W(nw,nw), Us(nw,nw), &
             VTs(nw,nw), sv(nw), rwork(5*nw))
    call zgesvd('A','A', nw, nw, O, nw, sv, Us, nw, VTs, nw, qw, -1, rwork, info)
    lwork = int(real(qw(1),dp)); allocate(work(lwork))
    seen = .false.; seen(k0) = .true.
    queue(1) = k0; qh = 1; qt = 1
    do while (qh <= qt)
       kp = queue(qh); qh = qh + 1
       do b = 1, nntot
          k = neigh(b,kp)
          if (seen(k)) cycle
          ! parent frame carried to k
          call zgemm('C','N', nw, nw, nw, CONE, mmn(1,1,b,kp), nw, &
                     A_gauge(1,1,kp), nw, CZERO, T, nw)
          ! Procrustes: W = polar(A_k^H T), then A_k <- A_k W
          call zgemm('C','N', nw, nw, nw, CONE, A_gauge(1,1,k), nw, T, nw, &
                     CZERO, O, nw)
          call zgesvd('A','A', nw, nw, O, nw, sv, Us, nw, VTs, nw, &
                      work, lwork, rwork, info)
          if (info /= 0) error stop 'align_gauge: SVD failed'
          call zgemm('N','N', nw, nw, nw, CONE, Us, nw, VTs, nw, CZERO, W, nw)
          T = A_gauge(:,:,k)
          call zgemm('N','N', nw, nw, nw, CONE, T, nw, W, nw, CZERO, &
                     A_gauge(1,1,k), nw)
          seen(k) = .true.
          qt = qt + 1; queue(qt) = k
       end do
    end do
    if (.not. all(seen)) error stop 'align_gauge: k-graph not connected'
  end subroutine align_gauge

  !------------------------------------------------------------ propagate_fill
  ! Multi-seed variant of propagate_gauge: KEEP the existing gauge at every
  ! k with keep(k) = .true. (healthy per-k anchors) and BFS-transport it
  ! only into the .not.keep region.  This is the middle ground between
  ! per-k SCDM (locally adaptive, but its quality is bounded by the WORST
  ! k) and single-seed transport (globally coherent, but rigid where the
  ! band character rotates): local anchoring wherever it is trustworthy,
  ! transport repair only where it is not.  Sick k adjacent to several
  ! healthy k are entered from whichever healthy front reaches them first;
  ! the subsequent smoothing sweeps reconcile the fronts.
  subroutine propagate_fill(mmn, neigh, nw, nntot, nk, keep, A_gauge)
    integer, intent(in) :: nw, nntot, nk
    complex(dp), intent(in) :: mmn(nw,nw,nntot,nk)
    integer, intent(in) :: neigh(nntot,nk)
    logical, intent(in) :: keep(nk)
    complex(dp), intent(inout) :: A_gauge(nw,nw,nk)  ! in: valid where keep
    integer, allocatable :: queue(:)
    logical, allocatable :: seen(:)
    complex(dp), allocatable :: X(:,:), Us(:,:), VTs(:,:), work(:)
    real(dp), allocatable :: sv(:), rwork(:)
    complex(dp) :: qw(1)
    integer :: qh, qt, kp, k, b, info, lwork
    if (.not. any(keep)) error stop 'propagate_fill: no healthy seeds'
    allocate(queue(nk), seen(nk), X(nw,nw), Us(nw,nw), VTs(nw,nw), sv(nw), &
             rwork(5*nw))
    call zgesvd('A','A', nw, nw, X, nw, sv, Us, nw, VTs, nw, qw, -1, rwork, info)
    lwork = int(real(qw(1),dp)); allocate(work(lwork))
    seen = keep
    qt = 0
    do k = 1, nk
       if (keep(k)) then; qt = qt + 1; queue(qt) = k; end if
    end do
    qh = 1
    do while (qh <= qt)
       kp = queue(qh); qh = qh + 1
       do b = 1, nntot
          k = neigh(b,kp)
          if (seen(k)) cycle
          call zgemm('C','N', nw, nw, nw, CONE, mmn(1,1,b,kp), nw, &
                     A_gauge(1,1,kp), nw, CZERO, X, nw)
          call zgesvd('A','A', nw, nw, X, nw, sv, Us, nw, VTs, nw, &
                      work, lwork, rwork, info)
          if (info /= 0) error stop 'propagate_fill: SVD failed'
          call zgemm('N','N', nw, nw, nw, CONE, Us, nw, VTs, nw, CZERO, &
                     A_gauge(1,1,k), nw)
          seen(k) = .true.
          qt = qt + 1; queue(qt) = k
       end do
    end do
    if (.not. all(seen)) error stop 'propagate_fill: k-graph not connected'
  end subroutine propagate_fill

  !------------------------------------------------------- gauge_spread_masked
  ! Seam-blind MV spread: like gauge_spread, but EXCLUDES edges whose
  ! anti-Hermitian deviation d(k,b) = ||M^U - (M^U)^H||_F exceeds
  ! dcut_rel x median(d).  On a raw transported gauge the tree edges have
  ! M^U exactly Hermitian-positive (unitary polar factor = 1 by
  ! construction), so d is a closed-form seam detector: ~0 on the tree,
  ! O(1) on loop-closing edges.  Scoring with the seams excluded measures
  ! the BASIN content of the frame without the construction artifact --
  ! no smoothing, no diffusion, no information loss.  Kept-weight
  ! renormalised per k so candidates with different trees stay comparable.
  subroutine gauge_spread_masked(mmn, neigh, w, nw, nntot, nk, U, dcut_rel, &
                                 om_score, nexcl)
    integer, intent(in) :: nw, nntot, nk
    complex(dp), intent(in) :: mmn(nw,nw,nntot,nk), U(nw,nw,nk)
    integer, intent(in) :: neigh(nntot,nk)
    real(dp), intent(in) :: w(nntot), dcut_rel
    real(dp), intent(out) :: om_score      ! OD+I-like gauge score, seam-blind
    integer, intent(out) :: nexcl
    complex(dp), allocatable :: MU(:,:), T(:,:)
    real(dp), allocatable :: dev(:,:), devs(:)
    real(dp) :: dmed, acc, wkept, wtot
    integer :: k, b, kb, m, n, nedge, i
    allocate(MU(nw,nw), T(nw,nw), dev(nntot,nk), devs(nntot*nk))
    ! pass 1: anti-Hermitian deviation per edge
    do k = 1, nk
       do b = 1, nntot
          kb = neigh(b,k)
          call zgemm('C','N', nw, nw, nw, CONE, U(1,1,k), nw, &
                     mmn(1,1,b,k), nw, CZERO, T, nw)
          call zgemm('N','N', nw, nw, nw, CONE, T, nw, U(1,1,kb), nw, &
                     CZERO, MU, nw)
          acc = 0.0_dp
          do n = 1, nw
             do m = 1, nw
                acc = acc + (abs(MU(m,n) - conjg(MU(n,m))))**2
             end do
          end do
          dev(b,k) = sqrt(acc)
       end do
    end do
    nedge = nntot*nk
    devs = reshape(dev, (/nedge/))
    ! median via partial insertion on a copy (nedge is modest)
    call rsort(devs, nedge)
    dmed = devs((nedge+1)/2)
    ! pass 2: accumulate the gauge score over kept edges only
    om_score = 0.0_dp; nexcl = 0
    wtot = sum(w)
    do k = 1, nk
       wkept = 0.0_dp; acc = 0.0_dp
       do b = 1, nntot
          if (dev(b,k) > dcut_rel*max(dmed, 1.0e-14_dp)) then
             nexcl = nexcl + 1; cycle
          end if
          kb = neigh(b,k)
          call zgemm('C','N', nw, nw, nw, CONE, U(1,1,k), nw, &
                     mmn(1,1,b,k), nw, CZERO, T, nw)
          call zgemm('N','N', nw, nw, nw, CONE, T, nw, U(1,1,kb), nw, &
                     CZERO, MU, nw)
          do n = 1, nw
             acc = acc + w(b)*(1.0_dp - (real(MU(n,n),dp)**2 + aimag(MU(n,n))**2))
          end do
          wkept = wkept + w(b)
       end do
       if (wkept > 1.0e-14_dp) om_score = om_score + acc*(wtot/wkept)
    end do
    om_score = om_score/real(nk,dp)
    deallocate(MU, T, dev, devs)
  end subroutine gauge_spread_masked

  ! ascending real sort (insertion; helper for the median above)
  subroutine rsort(x, n)
    integer, intent(in) :: n
    real(dp), intent(inout) :: x(n)
    integer :: i, j
    real(dp) :: key
    do i = 2, n
       key = x(i); j = i-1
       do while (j >= 1)
          if (x(j) <= key) exit
          x(j+1) = x(j); j = j-1
       end do
       x(j+1) = key
    end do
  end subroutine rsort

  !--------------------------------------------------------------- seam_smooth
  ! Seam-LOCAL smoothing: flag edges with anti-Hermitian deviation
  ! > dcut_rel x median (the closed-form seam detector), mark both
  ! endpoints, and run Gauss-Seidel polar sweeps ONLY on the marked k.
  ! Diffusion confined to the artifact's support -- the basin-encoding
  ! gauge structure away from the seams is untouched.
  subroutine seam_smooth(mmn, neigh, w, nw, nntot, nk, nsweep, dcut_rel, &
                         A_gauge, nmask)
    integer, intent(in) :: nw, nntot, nk, nsweep
    complex(dp), intent(in) :: mmn(nw,nw,nntot,nk)
    integer, intent(in) :: neigh(nntot,nk)
    real(dp), intent(in) :: w(nntot), dcut_rel
    complex(dp), intent(inout) :: A_gauge(nw,nw,nk)
    integer, intent(out) :: nmask
    logical, allocatable :: kmask(:)
    complex(dp), allocatable :: MU(:,:), T(:,:), Us(:,:), VTs(:,:), work(:)
    real(dp), allocatable :: dev(:,:), devs(:), sv(:), rwork(:)
    complex(dp) :: qw(1)
    real(dp) :: acc, dmed
    integer :: k, b, kb, m, n, is, info, lwork, nedge
    allocate(kmask(nk), MU(nw,nw), T(nw,nw), Us(nw,nw), VTs(nw,nw), &
             dev(nntot,nk), devs(nntot*nk), sv(nw), rwork(5*nw))
    call zgesvd('A','A', nw, nw, T, nw, sv, Us, nw, VTs, nw, qw, -1, rwork, info)
    lwork = int(real(qw(1),dp)); allocate(work(lwork))
    do k = 1, nk
       do b = 1, nntot
          kb = neigh(b,k)
          call zgemm('C','N', nw, nw, nw, CONE, A_gauge(1,1,k), nw, &
                     mmn(1,1,b,k), nw, CZERO, T, nw)
          call zgemm('N','N', nw, nw, nw, CONE, T, nw, A_gauge(1,1,kb), nw, &
                     CZERO, MU, nw)
          acc = 0.0_dp
          do n = 1, nw
             do m = 1, nw
                acc = acc + (abs(MU(m,n) - conjg(MU(n,m))))**2
             end do
          end do
          dev(b,k) = sqrt(acc)
       end do
    end do
    nedge = nntot*nk
    devs = reshape(dev, (/nedge/))
    call rsort(devs, nedge)
    dmed = max(devs((nedge+1)/2), 1.0e-14_dp)
    kmask = .false.
    do k = 1, nk
       do b = 1, nntot
          if (dev(b,k) > dcut_rel*dmed) then
             kmask(k) = .true.; kmask(neigh(b,k)) = .true.
          end if
       end do
    end do
    nmask = count(kmask)
    do is = 1, nsweep
       do k = 1, nk
          if (.not. kmask(k)) cycle
          T = CZERO
          do b = 1, nntot
             kb = neigh(b,k)
             call zgemm('N','N', nw, nw, nw, cmplx(w(b),0.0_dp,dp), &
                        mmn(1,1,b,k), nw, A_gauge(1,1,kb), nw, CONE, T, nw)
          end do
          call zgesvd('A','A', nw, nw, T, nw, sv, Us, nw, VTs, nw, &
                      work, lwork, rwork, info)
          if (info /= 0) error stop 'seam_smooth: SVD failed'
          call zgemm('N','N', nw, nw, nw, CONE, Us, nw, VTs, nw, CZERO, &
                     A_gauge(1,1,k), nw)
       end do
    end do
    deallocate(kmask, MU, T, Us, VTs, dev, devs, sv, rwork, work)
  end subroutine seam_smooth

  !---------------------------------------------------------- abelian_seam_fix
  ! One-shot ABELIAN seam removal: per band n, the diagonal-phase field
  ! phi_n(k,b) = Im ln M^U_nn(k,b) decomposes into gradient + curl parts
  ! on the k-graph.  The gradient part is removable exactly by a diagonal
  ! gauge U_k -> U_k diag(e^{-i theta_n(k)}) with theta the solution of
  ! the weighted graph-Poisson problem  L theta = rho,
  ! rho(k) = sum_b w_b phi_n(k,b)  -- solved here by conjugate gradient
  ! (the k-graph Laplacian is small).  This is the closed-form removal of
  ! the seam's abelian content; the gauge-invariant curl (Berry holonomy)
  ! and non-abelian mixing remain, untouched by any diffusion.
  ! Sign-guarded: applies +theta or -theta, whichever lowers the phase
  ! objective; no-op if neither helps.
  subroutine abelian_seam_fix(mmn, neigh, w, nw, nntot, nk, A_gauge)
    integer, intent(in) :: nw, nntot, nk
    complex(dp), intent(in) :: mmn(nw,nw,nntot,nk)
    integer, intent(in) :: neigh(nntot,nk)
    real(dp), intent(in) :: w(nntot)
    complex(dp), intent(inout) :: A_gauge(nw,nw,nk)
    complex(dp), allocatable :: MU(:,:), T(:,:)
    real(dp), allocatable :: phi(:,:,:), th(:), rr(:), pp(:), Ap(:)
    real(dp) :: rs, rs2, alpha, beta, f0, fp, fm, sgn
    integer :: k, b, kb, n, it, itmax
    allocate(MU(nw,nw), T(nw,nw), phi(nntot,nk,nw), th(nk), rr(nk), &
             pp(nk), Ap(nk))
    do k = 1, nk
       do b = 1, nntot
          kb = neigh(b,k)
          call zgemm('C','N', nw, nw, nw, CONE, A_gauge(1,1,k), nw, &
                     mmn(1,1,b,k), nw, CZERO, T, nw)
          call zgemm('N','N', nw, nw, nw, CONE, T, nw, A_gauge(1,1,kb), nw, &
                     CZERO, MU, nw)
          do n = 1, nw
             phi(b,k,n) = atan2(aimag(MU(n,n)), real(MU(n,n),dp))
          end do
       end do
    end do
    itmax = 400
    do n = 1, nw
       ! rho and CG solve of  L th = rho  (constant mode projected out)
       do k = 1, nk
          rr(k) = 0.0_dp
          do b = 1, nntot
             rr(k) = rr(k) + w(b)*phi(b,k,n)
          end do
       end do
       rr = rr - sum(rr)/real(nk,dp)
       th = 0.0_dp; pp = rr
       rs = dot_product(rr, rr)
       do it = 1, itmax
          if (rs < 1.0e-20_dp) exit
          do k = 1, nk
             Ap(k) = 0.0_dp
             do b = 1, nntot
                Ap(k) = Ap(k) + w(b)*(pp(k) - pp(neigh(b,k)))
             end do
          end do
          alpha = rs/max(dot_product(pp, Ap), 1.0e-300_dp)
          th = th + alpha*pp
          rr = rr - alpha*Ap
          rr = rr - sum(rr)/real(nk,dp)
          rs2 = dot_product(rr, rr)
          beta = rs2/max(rs, 1.0e-300_dp)
          pp = rr + beta*pp
          rs = rs2
       end do
       ! sign guard on the phase objective
       f0 = 0.0_dp; fp = 0.0_dp; fm = 0.0_dp
       do k = 1, nk
          do b = 1, nntot
             kb = neigh(b,k)
             f0 = f0 + w(b)*phi(b,k,n)**2
             fp = fp + w(b)*(phi(b,k,n) - th(k) + th(kb))**2
             fm = fm + w(b)*(phi(b,k,n) + th(k) - th(kb))**2
          end do
       end do
       ! applied phase e^{-i sgn th} on column n realises
       ! phi' = phi + sgn*(th(k) - th(kb)):  sgn=+1 <-> fm, sgn=-1 <-> fp
       sgn = 0.0_dp
       if (fp < f0 .and. fp <= fm) sgn = -1.0_dp
       if (fm < f0 .and. fm < fp)  sgn = 1.0_dp
       if (sgn /= 0.0_dp) then
          do k = 1, nk
             A_gauge(:,n,k) = A_gauge(:,n,k)* &
                  cmplx(cos(sgn*th(k)), -sin(sgn*th(k)), dp)
          end do
       end if
    end do
    deallocate(MU, T, phi, th, rr, pp, Ap)
  end subroutine abelian_seam_fix

  !------------------------------------------------------------ seed_tournament
  ! Loss-based transport-seed selection (2026-07-21 seed sweeps: NO fixed
  ! sigma-percentile is safe -- max-sigma seeds fail Fe-up, median fails
  ! Fe_soc).  For each candidate seed: propagate + nsw smoothing sweeps,
  ! then score with the TRUE MV spread (gauge_spread).  Keep the argmin
  ! gauge.  The high-sigma seed stays in the pool and wins exactly when its
  ! frame really is the best template (GaAs), loses when it is atypical
  ! (Fe-up) -- stability and flexibility balanced by the actual objective.
  subroutine seed_tournament(mmn, neigh, w, bvec, nw, nntot, nk, nsw, &
                             cand, ncand, A_gauge, k_win, om_win, lossmode, bkb)
    integer, intent(in) :: nw, nntot, nk, nsw, ncand, cand(ncand)
    complex(dp), intent(in) :: mmn(nw,nw,nntot,nk)
    integer, intent(in) :: neigh(nntot,nk)
    real(dp), intent(in) :: w(nntot), bvec(3,nntot)
    complex(dp), intent(inout) :: A_gauge(nw,nw,nk)  ! in: anchor gauge; out: winner
    integer, intent(out) :: k_win
    real(dp), intent(out) :: om_win
    character(*), intent(in), optional :: lossmode
    real(dp), intent(in), optional :: bkb(3,nntot,nk)   ! per-(k,b) b-vectors (gauge_spread)
    ! 'smooth' (default): propagate + nsw global sweeps + full spread
    ! 'blind' : raw frame, seam-blind spread (gauge_spread_masked)
    ! 'local' : raw frame + 8 SEAM-LOCAL sweeps (seam_smooth) + full spread
    ! 'abel'  : raw frame + one-shot abelian seam removal + full spread
    ! For blind/local/abel the winner additionally gets nsw GLOBAL sweeps
    ! iff nsw > 0 (so pts0 vs pts8 = without/with final smoothing).
    complex(dp), allocatable :: A0(:,:,:), At(:,:,:)
    real(dp) :: om, omp(3)
    integer :: ic, nexcl, nmsk, icw
    character(len=8) :: lm
    ! per-candidate results, so the parallel phase writes only to its own
    ! slot and the report stays in candidate order
    real(dp), allocatable :: sc(:), sc_od(:), sc_d(:)
    integer,  allocatable :: sc_n(:)
    lm = 'smooth'; if (present(lossmode)) lm = lossmode
    allocate(A0(nw,nw,nk))
    A0 = A_gauge
    allocate(sc(ncand), sc_od(ncand), sc_d(ncand), sc_n(ncand))
    sc = huge(1.0_dp); sc_od = 0; sc_d = 0; sc_n = 0

    ! ---- phase 1: score every candidate INDEPENDENTLY (embarrassingly
    ! parallel -- each seed gets its own transport + loss transform + spread
    ! with zero cross-talk).  `At` is a PRIVATE allocatable: an OpenMP
    ! private copy of an allocatable enters unallocated, so each thread
    ! allocates its own nw*nw*nk scratch (385 MB at nw=224/nk=480 -- the
    ! reason for schedule(dynamic) and modest thread counts).  No shared
    ! writes, so results are independent of thread count and scheduling.
    !$omp parallel do default(shared) schedule(dynamic) &
    !$omp   private(ic, At, om, omp, nexcl, nmsk)
    do ic = 1, ncand
       allocate(At(nw,nw,nk))
       At = A0
       call propagate_gauge(mmn, neigh, nw, nntot, nk, cand(ic), At)
       select case (trim(lm))
       case ('blind')
          call gauge_spread_masked(mmn, neigh, w, nw, nntot, nk, At, &
                                   3.0_dp, om, nexcl)
          sc_n(ic) = nexcl
       case ('local')
          call seam_smooth(mmn, neigh, w, nw, nntot, nk, 8, 3.0_dp, At, nmsk)
          call gauge_spread(mmn, neigh, w, bvec, nw, nntot, nk, At, om, omp, bkb=bkb)
          sc_n(ic) = nmsk; sc_od(ic) = omp(2); sc_d(ic) = omp(3)
       case ('abel')
          call abelian_seam_fix(mmn, neigh, w, nw, nntot, nk, At)
          call gauge_spread(mmn, neigh, w, bvec, nw, nntot, nk, At, om, omp, bkb=bkb)
          sc_od(ic) = omp(2); sc_d(ic) = omp(3)
       case ('absm')
          ! abelian fix THEN smoothing, BOTH before scoring: analytic seam
          ! removal + flow information, self-consistent (ship = scored).
          call abelian_seam_fix(mmn, neigh, w, nw, nntot, nk, At)
          call smooth_gauge(mmn, neigh, w, nw, nntot, nk, nsw, At)
          call gauge_spread(mmn, neigh, w, bvec, nw, nntot, nk, At, om, omp, bkb=bkb)
          sc_od(ic) = omp(2); sc_d(ic) = omp(3)
       case default
          call smooth_gauge(mmn, neigh, w, nw, nntot, nk, nsw, At)
          call gauge_spread(mmn, neigh, w, bvec, nw, nntot, nk, At, om, omp, bkb=bkb)
          sc_od(ic) = omp(2); sc_d(ic) = omp(3)
       end select
       sc(ic) = om
       deallocate(At)
    end do
    !$omp end parallel do

    ! ---- phase 2 (serial): report in candidate order and take the argmin.
    ! Lowest-index tie-break, exactly as the serial loop did.
    om_win = huge(1.0_dp); k_win = cand(1); icw = 1
    do ic = 1, ncand
       select case (trim(lm))
       case ('blind')
          write(*,'(a,i7,a,f14.4,a,i0,a)') '    seed cand k=', cand(ic), &
               '  blind score = ', sc(ic), '   (', sc_n(ic), ' seam edges excluded)'
       case ('local')
          write(*,'(a,i7,a,f14.4,a,i0,a)') '    seed cand k=', cand(ic), &
               '  local score = ', sc(ic), '   (', sc_n(ic), ' seam-adjacent k)'
       case ('abel')
          write(*,'(a,i7,a,f14.4,a,2f12.4,a)') '    seed cand k=', cand(ic), &
               '  abel score = ', sc(ic), '   (OD/D = ', sc_od(ic), sc_d(ic), ' )'
       case ('absm')
          write(*,'(a,i7,a,f14.4,a,2f12.4,a)') '    seed cand k=', cand(ic), &
               '  absm score = ', sc(ic), '   (OD/D = ', sc_od(ic), sc_d(ic), ' )'
       case default
          write(*,'(a,i7,a,f14.4,a,2f12.4,a)') '    seed cand k=', cand(ic), &
               '  Omega = ', sc(ic), '   (OD/D = ', sc_od(ic), sc_d(ic), ' )'
       end select
       if (sc(ic) < om_win) then
          om_win = sc(ic); k_win = cand(ic); icw = ic
       end if
    end do

    ! ---- phase 3: rebuild ONLY the winner's frame (costs 1/ncand of the
    ! tournament; buys a race-free loop and no ncand-sized gauge storage --
    ! holding every candidate would be ncand*385 MB at Rh3Si4 scale).
    allocate(At(nw,nw,nk))
    At = A0
    call propagate_gauge(mmn, neigh, nw, nntot, nk, cand(icw), At)
    select case (trim(lm))
    case ('blind')
       continue                       ! raw frame is what was scored
    case ('local')
       call seam_smooth(mmn, neigh, w, nw, nntot, nk, 8, 3.0_dp, At, nmsk)
    case ('abel')
       call abelian_seam_fix(mmn, neigh, w, nw, nntot, nk, At)
    case ('absm')
       call abelian_seam_fix(mmn, neigh, w, nw, nntot, nk, At)
       call smooth_gauge(mmn, neigh, w, nw, nntot, nk, nsw, At)
    case default
       call smooth_gauge(mmn, neigh, w, nw, nntot, nk, nsw, At)
    end select
    A_gauge = At
    deallocate(sc, sc_od, sc_d, sc_n)
    if ((trim(lm) == 'blind' .or. trim(lm) == 'local' .or. &
         trim(lm) == 'abel') .and. nsw > 0) &
       call smooth_gauge(mmn, neigh, w, nw, nntot, nk, nsw, A_gauge)
    deallocate(A0, At)
  end subroutine seed_tournament

  !------------------------------------------------------------- smooth_gauge
  ! Loop-consistent gauge smoothing: Gauss-Seidel sweeps of the polar update
  !   U_k <- polar( sum_b w_b M(k,b) U_{k+b} ),
  ! the natural relaxation of the gauge-dependent spread term.  Unlike the
  ! spanning-tree transport (propagate_gauge), which concentrates ALL loop
  ! obstruction on tree-closing edges (initial Omega grows with mesh depth),
  ! smoothing distributes it evenly over the BZ.  Seed with any gauge
  ! (SCDM columns recommended); a handful of sweeps suffices.
  subroutine smooth_gauge(mmn, neigh, w, nw, nntot, nk, nsweep, A_gauge)
    integer, intent(in) :: nw, nntot, nk, nsweep
    complex(dp), intent(in) :: mmn(nw,nw,nntot,nk)
    integer, intent(in) :: neigh(nntot,nk)
    real(dp), intent(in) :: w(nntot)
    complex(dp), intent(inout) :: A_gauge(nw,nw,nk)
    complex(dp), allocatable :: X(:,:), Us(:,:), VTs(:,:), work(:), Aold(:,:)
    real(dp), allocatable :: sv(:), rwork(:)
    complex(dp) :: qw(1)
    integer :: is, k, b, kb, info, lwork, ios_t
    real(dp) :: smtol, maxdif
    character(len=16) :: tolenv
    ! env W2W_SMTOL=<tol>: converged-smoothing mode -- exit as soon as the
    ! largest per-sweep gauge change drops below tol (nsweep = cap).
    smtol = -1.0_dp
    call get_environment_variable('W2W_SMTOL', tolenv)
    if (len_trim(tolenv) > 0) then
       read(tolenv, *, iostat=ios_t) smtol
       if (ios_t /= 0) smtol = -1.0_dp
    end if
    allocate(X(nw,nw), Us(nw,nw), VTs(nw,nw), sv(nw), rwork(5*nw), Aold(nw,nw))
    call zgesvd('A','A', nw, nw, X, nw, sv, Us, nw, VTs, nw, qw, -1, rwork, info)
    lwork = int(real(qw(1),dp)); allocate(work(lwork))
    do is = 1, nsweep
       maxdif = 0.0_dp
       do k = 1, nk
          X = CZERO
          do b = 1, nntot
             kb = neigh(b,k)
             ! X += w_b * M(k,b) U_kb   (in-place Gauss-Seidel: freshest U)
             call zgemm('N','N', nw, nw, nw, cmplx(w(b),0.0_dp,dp), &
                        mmn(1,1,b,k), nw, A_gauge(1,1,kb), nw, CONE, X, nw)
          end do
          call zgesvd('A','A', nw, nw, X, nw, sv, Us, nw, VTs, nw, &
                      work, lwork, rwork, info)
          if (info /= 0) error stop 'smooth_gauge: SVD failed'
          if (smtol > 0.0_dp) Aold = A_gauge(:,:,k)
          call zgemm('N','N', nw, nw, nw, CONE, Us, nw, VTs, nw, CZERO, &
                     A_gauge(1,1,k), nw)
          if (smtol > 0.0_dp) &
             maxdif = max(maxdif, maxval(abs(A_gauge(:,:,k) - Aold)))
       end do
       if (smtol > 0.0_dp .and. maxdif < smtol) then
          write(*,'(a,i0,a,es9.2,a)') '    (smoothing converged: ', is, &
               ' sweeps, max change ', maxdif, ')'
          exit
       end if
    end do
  end subroutine smooth_gauge

  !--------------------------------------------------------------- gauge_spread
  ! Full Marzari-Vanderbilt spread of a candidate gauge U (before Wannier90):
  ! M^U(k,b) = U_k^H M(k,b) U_{k+b};  centres r_n = -(1/nk) sum_kb w_b b ImlnM_nn;
  ! Omega = OmegaI + OmegaOD + OmegaD.  Lets gauge variants be scored in
  ! seconds instead of via a Wannier90 startup parse.
  subroutine gauge_spread(mmn, neigh, w, bvec, nw, nntot, nk, U, om_tot, om_parts, bkb)
    integer, intent(in) :: nw, nntot, nk
    complex(dp), intent(in) :: mmn(nw,nw,nntot,nk), U(nw,nw,nk)
    integer, intent(in) :: neigh(nntot,nk)
    real(dp), intent(in) :: w(nntot), bvec(3,nntot)
    real(dp), intent(out) :: om_tot, om_parts(3)   ! I, OD, D
    !> The b-vector of each (k,b) block (hybrid_files_m nnkp_block_bvectors).  A
    !> .nnkp orders each k's neighbours its own way, so bvec(:,b) -- k = 1's -- is
    !> the wrong DIRECTION for block b at most k: without bkb the centres and
    !> Omega_D come out wrong (MgB2 hand-off: Omega_D 29.17 against Wannier90's
    !> 0.036).  Omega_I and Omega_OD use only the per-shell weights w(b) and are
    !> unaffected.  Absent = bvec(:,b) at every k.
    real(dp), intent(in), optional :: bkb(3,nntot,nk)
    complex(dp), allocatable :: MU(:,:), T(:,:), Mdiag(:,:,:)
    real(dp), allocatable :: rbar(:,:), lnim(:,:,:)
    real(dp) :: oi, ood, od, phi
    integer :: k, b, kb, m, n
    allocate(MU(nw,nw), T(nw,nw), Mdiag(nw,nntot,nk), lnim(nw,nntot,nk), &
             rbar(3,nw))
    oi = 0.0_dp; ood = 0.0_dp
    do k = 1, nk
       do b = 1, nntot
          kb = neigh(b,k)
          call zgemm('C','N', nw, nw, nw, CONE, U(1,1,k), nw, &
                     mmn(1,1,b,k), nw, CZERO, T, nw)
          call zgemm('N','N', nw, nw, nw, CONE, T, nw, U(1,1,kb), nw, &
                     CZERO, MU, nw)
          do n = 1, nw
             Mdiag(n,b,k) = MU(n,n)
             lnim(n,b,k) = -aimag(log(MU(n,n)))
             do m = 1, nw
                if (m == n) cycle
                ood = ood + w(b)*(real(MU(m,n),dp)**2 + aimag(MU(m,n))**2)
             end do
             oi = oi + w(b)*(1.0_dp - (real(MU(n,n),dp)**2 + aimag(MU(n,n))**2))
          end do
       end do
    end do
    ! centres from the diagonal phases
    rbar = 0.0_dp
    do k = 1, nk
       do b = 1, nntot
          do n = 1, nw
             if (present(bkb)) then
                rbar(:,n) = rbar(:,n) + w(b)*bkb(:,b,k)*lnim(n,b,k)
             else
                rbar(:,n) = rbar(:,n) + w(b)*bvec(:,b)*lnim(n,b,k)
             end if
          end do
       end do
    end do
    rbar = rbar/real(nk,dp)
    od = 0.0_dp
    do k = 1, nk
       do b = 1, nntot
          do n = 1, nw
             if (present(bkb)) then
                phi = lnim(n,b,k) - dot_product(bkb(:,b,k), rbar(:,n))
             else
                phi = lnim(n,b,k) - dot_product(bvec(:,b), rbar(:,n))
             end if
             od = od + w(b)*phi*phi
          end do
       end do
    end do
    om_parts(1) = (oi - ood)/real(nk,dp)   ! true Omega_I (leakage only)
    om_parts(2) = ood/real(nk,dp)
    om_parts(3) = od/real(nk,dp)
    om_tot = sum(om_parts)
  end subroutine gauge_spread

  !--------------------------------------------------------------- scdm_anchors
  ! AO/channel SCDM gauge from k-independent anchor columns (QRCP), choosing the
  ! (anchor set, gauge) with the best worst-k conditioning.
  !--------------------------------------------------------- group closure
  ! Make a QRCP pivot set closed under whole "groups" of candidate columns.
  ! A group = every m of one (atom, l, rf) shell, in BOTH sigma blocks, so a
  ! single rule delivers two properties of a representation basis:
  !   * SHELL COMPLETENESS -- a shell is taken whole or not at all (tools
  !     that label projections by shell cannot express 4 of 5 d orbitals).
  !   * SIGMA PAIRING -- a SOC basis must be {spatial} x {up,down} for a
  !     spinor rotation to act on it.
  ! Groups are taken in QRCP RANK order, so SCDM still chooses WHICH shells
  ! matter; the constraint only forbids taking a shell in part.  Groups that
  ! would overrun nw are skipped (not truncated).  Any shortfall is filled
  ! with the best remaining individual columns and reported -- a partial fill
  ! means nw is not expressible as a sum of shell sizes.
  subroutine groupify_pivots(piv_ranked, nw, grp, nao, nfull, uout)
    integer, intent(in) :: nw, nao, grp(nao), uout
    integer, intent(inout) :: piv_ranked(nw)
    integer, intent(out) :: nfull
    integer :: i, j, g, s, ngr, cnt, cand, nwx, out(nw)
    logical :: taken(nao)
    integer, allocatable :: gsz(:), gsc(:), gbest(:), tab(:,:), ord(:)
    logical, allocatable :: use_g(:)
    ! Choosing whole shells is a subset-sum, NOT a greedy walk: picking groups
    ! in QRCP rank order and skipping any that would overrun nw can strand
    ! short of nw even when an exact closure exists (alpha-Sn nw=20 stranded at
    ! 16 while both {d,d} = 10+10 and {s,p,d,s} = 2+6+10+2 close exactly).
    ! With one group per (atom,l,rf) shell there are only a handful of groups,
    ! so solve it exactly: 0/1 knapsack over group SIZES constrained to sum to
    ! nw, maximising a rank score so that among exact closures we keep the one
    ! SCDM ranked highest.
    ngr = maxval(grp)
    if (ngr < 1) then
       nfull = 0
       if (uout > 0) write(uout,'(a)') &
            '  anchors: no groups -- unconstrained SCDM'
       return
    end if
    allocate(gsz(ngr), gsc(ngr), gbest(ngr), use_g(ngr), ord(ngr))
    gsz = 0; gsc = 0; gbest = huge(1)
    do j = 1, nao
       g = grp(j)
       if (g > 0) gsz(g) = gsz(g) + 1
    end do
    do i = 1, nw                        ! rank score: earlier pivot = worth more
       g = grp(piv_ranked(i))
       if (g > 0) then
          gsc(g) = gsc(g) + (nw + 1 - i)
          gbest(g) = min(gbest(g), i)
       end if
    end do
    ! tab(g,s) = best score using only groups 1..g summing to exactly s; -1 =
    ! unreachable.  Carried past nw so a failed closure can name the nearest
    ! nw that WOULD close, in both directions.
    nwx = min(nw + 12, nao)
    allocate(tab(0:ngr, 0:nwx))
    tab = -1; tab(0,0) = 0
    do g = 1, ngr
       do s = 0, nwx
          tab(g,s) = tab(g-1,s)                     ! skip group g
          if (gsz(g) > 0 .and. s >= gsz(g)) then
             if (tab(g-1, s-gsz(g)) >= 0) then
                cand = tab(g-1, s-gsz(g)) + gsc(g)
                if (cand > tab(g,s)) tab(g,s) = cand
             end if
          end if
       end do
    end do
    use_g = .false.
    if (tab(ngr,nw) >= 0) then          ! exact closure exists -- reconstruct it
       s = nw
       do g = ngr, 1, -1
          if (tab(g,s) /= tab(g-1,s)) then
             use_g(g) = .true.; s = s - gsz(g)
          end if
       end do
    end if
    ! emit chosen groups ordered by their best pivot rank (deterministic:
    ! ties broken by group id), columns within a group in index order
    do g = 1, ngr
       ord(g) = g
    end do
    do i = 1, ngr-1                     ! insertion sort on (gbest, id)
       do j = i+1, ngr
          if (gbest(ord(j)) < gbest(ord(i)) .or. &
             (gbest(ord(j)) == gbest(ord(i)) .and. ord(j) < ord(i))) then
             g = ord(i); ord(i) = ord(j); ord(j) = g
          end if
       end do
    end do
    taken = .false.; cnt = 0
    do i = 1, ngr
       g = ord(i)
       if (.not. use_g(g)) cycle
       do j = 1, nao
          if (grp(j) == g) then
             cnt = cnt + 1; out(cnt) = j; taken(j) = .true.
          end if
       end do
    end do
    nfull = cnt
    do i = 1, nw                        ! fill any remainder by rank
       if (cnt >= nw) exit
       if (.not. taken(piv_ranked(i))) then
          cnt = cnt + 1; out(cnt) = piv_ranked(i); taken(piv_ranked(i)) = .true.
       end if
    end do
    do j = 1, nao                       ! last resort: any unused column
       if (cnt >= nw) exit
       if (.not. taken(j)) then
          cnt = cnt + 1; out(cnt) = j; taken(j) = .true.
       end if
    end do
    piv_ranked = out
    if (uout > 0) then
       if (nfull == nw) then
          write(uout,'(a,i0,a,i0,a)') '  anchors: shell/sigma-complete (', &
               nw, ' columns in ', count(use_g), ' whole shells)'
       else
          ! genuinely unreachable -- tell the user which nw WOULD close, since
          ! that is the only actionable remedy
          write(uout,'(a,i0,a)') '  anchors: NO exact shell closure for nw=', &
               nw, ' -- nw is not a sum of shell sizes; falling back to rank'
          write(uout,'(a)',advance='no') '    nearest closable nw:'
          do s = max(1, nw-12), nwx
             if (s /= nw .and. tab(ngr,s) >= 0) write(uout,'(1x,i0)',advance='no') s
          end do
          write(uout,'(a)') ''
       end if
    end if
    deallocate(gsz, gsc, gbest, use_g, ord, tab)
  end subroutine groupify_pivots

  !------------------------------------------------- projmap -> group ids
  ! Build a group id per amn column from projmap.dat: columns sharing an
  ! (atom, l, rf) shell get the same id, and the sigma-2 block reuses the
  ! sigma-1 ids, so one grouping enforces shell completeness AND sigma
  ! pairing simultaneously.  grp stays unallocated if projmap is missing or
  ! its row count does not divide nao.
  !--------------------------------------------------------- read_orbit_map
  ! Orbit-closed anchors.  Shell grouping keyed on (atom,l,rf) is
  ! SITE-LOCAL: equivalent atoms get distinct ids, so the knapsack may close
  ! a shell on one atom of an orbit and not on its partners.
  !
  ! The orbit decomposition needs no analysis here: WIEN2k's struct
  ! already groups equivalent positions under each NON-EQUIVALENT atom, and
  ! global atom indices run sequentially through those blocks.  So reading
  ! the MULT values is enough to map global atom -> orbit, with no
  ! structmod dependency (the standalone has no struct access by design).
  !
  ! Enabled by env W2W_STRUCT=<case.struct>; absent => site-local grouping,
  ! bit-identical to before.
  subroutine read_orbit_map(fname, orb, nat, ok)
    character(*), intent(in) :: fname
    integer, allocatable, intent(out) :: orb(:)
    integer, intent(out) :: nat
    logical, intent(out) :: ok
    integer :: u, ios, mult, ineq, g, i
    integer, parameter :: MAXAT = 4096
    integer :: tmp(MAXAT)
    character(len=256) :: line
    ok = .false.; nat = 0
    open(newunit=u, file=trim(fname), status='old', action='read', iostat=ios)
    if (ios /= 0) return
    ineq = 0; g = 0
    do
       read(u,'(a)',iostat=ios) line
       if (ios /= 0) exit
       if (index(line,'MULT=') > 0) then
          read(line(index(line,'MULT=')+5:), *, iostat=ios) mult
          if (ios /= 0) cycle
          ineq = ineq + 1
          do i = 1, mult
             g = g + 1
             if (g > MAXAT) exit
             tmp(g) = ineq
          end do
       end if
    end do
    close(u)
    if (g == 0) return
    nat = g; allocate(orb(g)); orb = tmp(1:g); ok = .true.
  end subroutine read_orbit_map

  subroutine read_projmap_groups(nao, grp, uout)
    integer, intent(in) :: nao, uout
    integer, allocatable, intent(out) :: grp(:)
    integer :: u, ios, n, i, j, iat, l, m, rf, nrows, ng
    integer, allocatable :: at(:), ll(:), rr(:), sid(:)
    character(len=256) :: line
    open(newunit=u, file='projmap.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) return
    allocate(at(nao), ll(nao), rr(nao))
    nrows = 0
    do
       read(u,'(a)',iostat=ios) line
       if (ios /= 0) exit
       if (index(adjustl(line),'#') == 1 .or. len_trim(line) == 0) cycle
       read(line,*,iostat=ios) iat, l, m, rf
       if (ios /= 0) exit
       nrows = nrows + 1
       if (nrows > nao) exit
       at(nrows) = iat; ll(nrows) = l; rr(nrows) = rf
    end do
    close(u)
    ! accept nrows == nao (collinear) or nao == 2*nrows (SOC sigma blocks)
    if (.not. (nrows == nao .or. 2*nrows == nao)) then
       deallocate(at, ll, rr); return
    end if
    ! ARM C: optionally relabel the atom index by its ORBIT, so that a
    ! shell is closed across all symmetry-equivalent atoms at once.
    orbblk: block
      integer, allocatable :: orb(:)
      integer :: natorb, q
      logical :: okorb
      character(len=256) :: sfile
      call get_environment_variable('W2W_STRUCT', sfile)
      if (len_trim(sfile) > 0) then
         call read_orbit_map(trim(sfile), orb, natorb, okorb)
         if (okorb) then
            do q = 1, nrows
               if (at(q) >= 1 .and. at(q) <= natorb) at(q) = -orb(at(q))
            end do
            if (uout > 0) write(uout,'(a,i0,a)') &
                 '  anchor groups: ORBIT-CLOSED (arm C) over ', natorb, &
                 ' atoms from '//trim(sfile)
            deallocate(orb)
         else if (uout > 0) then
            write(uout,'(a)') '  W2W_STRUCT set but unreadable -- '// &
                 'site-local grouping kept'
         end if
      end if
    end block orbblk
    allocate(sid(nrows)); sid = 0; ng = 0
    do i = 1, nrows                       ! unique (atom-or-orbit,l,rf) -> shell
       do j = 1, i-1
          if (at(j) == at(i) .and. ll(j) == ll(i) .and. rr(j) == rr(i)) then
             sid(i) = sid(j); exit
          end if
       end do
       if (sid(i) == 0) then
          ng = ng + 1; sid(i) = ng
       end if
    end do
    allocate(grp(nao))
    do i = 1, nao
       grp(i) = sid(mod(i-1, nrows) + 1)   ! sigma-2 reuses sigma-1 ids
    end do
    if (uout > 0) write(uout,'(a,i0,a,i0,a)') &
         '  anchor groups: ', ng, ' shells over ', nao, ' candidate columns'
    deallocate(at, ll, rr, sid)
  end subroutine read_projmap_groups

  subroutine scdm_anchors(B, nao, nk, nw, gamma_index, anchors, A_gauge, best_svmin, grp, &
                          amode)
    integer, intent(in) :: nao, nk, nw, gamma_index
    complex(dp), intent(in) :: B(nao,nw,nk)
    integer, intent(out) :: anchors(nw)
    complex(dp), intent(out) :: A_gauge(nw,nw,nk)
    real(dp), intent(out) :: best_svmin
    ! optional group id per candidate column (0 = ungrouped); when present the
    ! pivot sets are made group-closed before the gauge bake-off
    integer, intent(in), optional :: grp(nao)
    ! Anchor-selection mode.  '' / absent = plain QRCP (the default).
    !   'pweight' -- PROJECTABILITY-WEIGHTED QRCP: scale each candidate column
    !     by its projectability onto the disentangled subspace before pivoting.
    !     The weight costs nothing: P_j = (1/nk) sum_k ||B(j,:,k)||^2 is
    !     exactly real(rho(j,j)), already computed here.  Pivots then favour
    !     columns that carry the bands, while QRCP still refuses near-parallel
    !     picks -- the failure mode pure P-ranking cannot see (two high-P
    !     columns of the same character on equivalent atoms).  The GAUGE is
    !     built from the UNWEIGHTED columns; weights bias selection only.
    !   'pshell' -- P-RANKED WHOLE SHELLS: rank candidates by P, close under
    !     the (atom,l,rf)-shell groups with the same exact knapsack the
    !     shellcomplete mode uses.  No QRCP in the selection at all; requires
    !     grp (projmap.dat).  Conditioning is still guarded downstream by the
    !     bake-off's sigma_min and choose_gauge's auto gate.
    character(len=*), intent(in), optional :: amode
    complex(dp), allocatable :: rho(:,:), tmp(:,:), Bgh(:,:), Acol(:,:,:), Arow(:,:,:)
    complex(dp), allocatable :: rhoa(:,:), svbuf(:,:)
    real(dp),    allocatable :: pwt(:)
    integer,     allocatable :: prank(:)
    integer :: cand_piv(nw, 2), piv(nw), k, i, j, ic, ig, ncand
    real(dp) :: sv, worst
    integer :: best_set, best_gauge, nfull1, nfull2
    character(len=16) :: am

    ! rho = (1/nk) sum_k B_k B_k^H   (Hermitian)
    allocate(rho(nao,nao))
    rho = CZERO
    do k = 1, nk
       call zherk('U','N', nao, nw, 1.0_dp, B(1,1,k), nao, 1.0_dp, rho, nao)
    end do
    rho = rho/real(nk,dp)
    do j = 1, nao                     ! fill lower triangle
       do i = j+1, nao
          rho(i,j) = conjg(rho(j,i))
       end do
       rho(j,j) = cmplx(real(rho(j,j),dp), 0.0_dp, dp)
    end do

    am = ''
    if (present(amode)) am = amode

    ! per-column projectability onto the disentangled subspace, normalized to
    ! max 1 (used by both P modes; free -- it is rho's diagonal)
    allocate(pwt(nao))
    do j = 1, nao
       pwt(j) = real(rho(j,j), dp)
    end do
    pwt = pwt / max(maxval(pwt), 1.0e-300_dp)

    if (trim(am) == 'pshell') then
       ! ---- option A: P-ranked whole shells, no QRCP in the selection ----
       ! prank = the nw highest-P columns in DESCENDING P order; groupify's
       ! knapsack then closes it under whole shells, maximizing the same
       ! rank score it uses for QRCP order.
       allocate(prank(nao))
       do j = 1, nao
          prank(j) = j
       end do
       do i = 1, nao-1                 ! insertion sort, descending P
          do j = i+1, nao
             if (pwt(prank(j)) > pwt(prank(i))) then
                k = prank(i); prank(i) = prank(j); prank(j) = k
             end if
          end do
       end do
       cand_piv(:,1) = prank(1:nw)
       if (.not. present(grp)) &
            error stop 'scdm_anchors: pshell needs projmap groups (grp)'
       call groupify_pivots(cand_piv(:,1), nw, grp, nao, nfull1, 6)
       call isort(cand_piv(:,1), nw)
       cand_piv(:,2) = cand_piv(:,1)   ! single set; bake-off still tries both gauges
       ncand = 1
       deallocate(prank)
    else
       ! candidate 1: QRCP on rho ; candidate 2: QRCP on B[gamma]^H.
       ! Under pweight, scale COLUMN j of each QRCP input by pwt(j): pivoting
       ! ranks by column norm, so this biases selection toward high-P channels
       ! while the orthogonalization inside QRCP still rejects near-parallel
       ! columns.  Selection only -- rho and B themselves are untouched.
       allocate(tmp(nao,nao))
       tmp = rho
       if (trim(am) == 'pweight') then
          do j = 1, nao
             tmp(:,j) = tmp(:,j) * pwt(j)
          end do
       end if
       call qrcp_pivots(tmp, nao, nao, nw, cand_piv(:,1))
       ! group closure BEFORE the ascending sort: groupify needs QRCP rank order
       if (present(grp)) call groupify_pivots(cand_piv(:,1), nw, grp, nao, nfull1, 6)
       call isort(cand_piv(:,1), nw)
       allocate(Bgh(nw,nao))
       do j = 1, nao
          do i = 1, nw
             Bgh(i,j) = conjg(B(j,i,gamma_index))
          end do
       end do
       if (trim(am) == 'pweight') then
          do j = 1, nao
             Bgh(:,j) = Bgh(:,j) * pwt(j)
          end do
       end if
       call qrcp_pivots(Bgh, nw, nao, nw, cand_piv(:,2))
       if (present(grp)) call groupify_pivots(cand_piv(:,2), nw, grp, nao, nfull2, 0)
       call isort(cand_piv(:,2), nw)
       deallocate(tmp, Bgh)
       ncand = 2
    end if
    deallocate(pwt)

    ! evaluate both gauges for both candidate sets; keep best worst-k sigma_min
    allocate(rhoa(nao,nw), Acol(nw,nw,nk), Arow(nw,nw,nk), svbuf(nw,nw))
    best_svmin = -1.0_dp; best_set = 1; best_gauge = 1
    do ic = 1, ncand
       piv = cand_piv(:,ic)
       do j = 1, nw
          rhoa(:,j) = rho(:, piv(j))         ! rho[:, anchors]
       end do
       do ig = 1, 2
          worst = huge(1.0_dp)
          do k = 1, nk
             if (ig == 1) then
                ! col gauge:  A_k = B_k^H rho[:,anchors]
                call zgemm('C','N', nw, nw, nao, CONE, B(1,1,k), nao, rhoa, nao, &
                           CZERO, svbuf, nw)
             else
                ! row gauge:  A_k(m,n) = conj(B_k[anchors(n), m])
                do j = 1, nw
                   do i = 1, nw
                      svbuf(i,j) = conjg(B(piv(j), i, k))
                   end do
                end do
             end if
             sv = smallest_sv(svbuf, nw)      ! svbuf destroyed
             if (sv < worst) worst = sv
             if (worst <= best_svmin) exit     ! cannot beat current best
          end do
          if (worst > best_svmin) then
             best_svmin = worst; best_set = ic; best_gauge = ig
          end if
       end do
    end do

    ! recompute A_gauge for the winning (set,gauge)
    anchors = cand_piv(:, best_set)
    piv = anchors
    do j = 1, nw
       rhoa(:,j) = rho(:, piv(j))
    end do
    do k = 1, nk
       if (best_gauge == 1) then
          call zgemm('C','N', nw, nw, nao, CONE, B(1,1,k), nao, rhoa, nao, &
                     CZERO, A_gauge(1,1,k), nw)
       else
          do j = 1, nw
             do i = 1, nw
                A_gauge(i,j,k) = conjg(B(piv(j), i, k))
             end do
          end do
       end if
    end do
    deallocate(rho, rhoa, Acol, Arow, svbuf)
  end subroutine scdm_anchors

  ! Block-SCDM: run SCDM independently within each spectrally isolated
  ! manifold of the rotated subspace, assembling a block-diagonal initial
  ! gauge.  A single global QRCP cannot allocate anchor columns across a
  ! flat, near-degenerate sub-manifold (e.g. a semicore shell) and a
  ! dispersive one at once -> the global gauge goes rank-deficient at the
  ! worst k (GaAs Ga-3d: sigma_min ~ 1e-16).  Per-block QRCP gives each
  ! manifold anchors chosen to span ITSELF, curing the singularity and
  ! yielding clean per-shell orbitals.  Valid because the blocks are
  ! gapped at every k (a block boundary is placed only where
  ! min_k(e_{i+1,k} - e_{i,k}) > gap_min); Wannier90 may still mix blocks
  ! afterward, so this only sets a well-conditioned starting frame.
  ! Single-block windows reproduce scdm_anchors bit-for-bit.
  subroutine block_scdm(B, eig_t, nao, nk, nw, gap_min, anchors, A_gauge, &
                        best_svmin, nblock, verbose)
    integer, intent(in) :: nao, nk, nw
    complex(dp), intent(in) :: B(nao,nw,nk)
    real(dp), intent(in) :: eig_t(nw,nk), gap_min
    integer, intent(out) :: anchors(nw), nblock
    complex(dp), intent(out) :: A_gauge(nw,nw,nk)
    real(dp), intent(out) :: best_svmin
    logical, intent(in) :: verbose
    integer :: bstart(nw+1), nb, i, k, i0, i1, nbk, j
    real(dp) :: gmin, svb
    complex(dp), allocatable :: Bb(:,:,:), Ab(:,:,:)
    integer, allocatable :: anb(:)

    ! block boundaries: cut after state i iff the direct gap to i+1 stays
    ! open (> gap_min) at every k -> a genuinely isolated manifold
    nb = 0
    bstart(1) = 1
    do i = 1, nw-1
       gmin = huge(1.0_dp)
       do k = 1, nk
          gmin = min(gmin, eig_t(i+1,k) - eig_t(i,k))
       end do
       if (gmin > gap_min) then
          nb = nb + 1
          bstart(nb+1) = i+1
       end if
    end do
    nb = nb + 1
    bstart(nb+1) = nw + 1
    nblock = nb

    A_gauge = CZERO
    best_svmin = huge(1.0_dp)
    do j = 1, nb
       i0 = bstart(j); i1 = bstart(j+1) - 1; nbk = i1 - i0 + 1
       allocate(Bb(nao,nbk,nk), Ab(nbk,nbk,nk), anb(nbk))
       Bb = B(:, i0:i1, :)
       call scdm_anchors(Bb, nao, nk, nbk, 1, anb, Ab, svb)
       do k = 1, nk
          A_gauge(i0:i1, i0:i1, k) = Ab(:,:,k)
       end do
       anchors(i0:i1) = anb              ! global candidate-column indices
       best_svmin = min(best_svmin, svb)
       if (verbose) write(*,'(a,i0,a,i0,a,i0,a,i0,a,es11.3)') &
            '   block ', j, ': states ', i0, '-', i1, ' (', nbk, &
            ' WF)  sigma_min=', svb
       deallocate(Bb, Ab, anb)
    end do
  end subroutine block_scdm

  ! Index of the best-conditioned k (largest smallest-singular-value of the
  ! per-k gauge).  Seeding parallel transport / smoothing from the healthiest
  ! frame -- rather than the arbitrary k=1 -- gives cleaner propagation: the
  ! transport carries the well-defined frame outward instead of possibly
  ! starting on a near-singular one.  This is the one place SCDM's per-k
  ! conditioning feeds the PTG continuity step new information.
  function gauge_best_k(A, nw, nk) result(k0)
    integer, intent(in) :: nw, nk
    complex(dp), intent(in) :: A(nw,nw,nk)
    integer :: k0, k
    real(dp) :: best, sv
    complex(dp) :: buf(nw,nw)
    k0 = 1; best = -1.0_dp
    do k = 1, nk
       buf = A(:,:,k)
       sv = smallest_sv(buf, nw)        ! buf destroyed
       if (sv > best) then; best = sv; k0 = k; end if
    end do
  end function gauge_best_k

  ! Plain Loewdin projection gauge: trial functions = the leading nw
  ! candidate columns (physical channels, projmap order) rather than a
  ! QRCP-selected anchor set.  A_gauge(m,n,k) = <psi_m^subspace | phi_n> =
  ! conj(B(n,m,k)); Wannier90 orthonormalizes on hand-off (Nb=Nw =>
  ! U = A(A^H A)^-1/2).  This is the classic W90 projection init and the
  ! closest native analogue of the paper's PDWF+gauge-in-AO-space chain.
  ! Reports the worst-k sigma_min of A so the caller can compare
  ! conditioning against the SCDM anchor set.
  subroutine lowdin_gauge(B, nao, nk, nw, sel, A_gauge, worst_svmin)
    integer, intent(in) :: nao, nk, nw
    complex(dp), intent(in) :: B(nao,nw,nk)
    integer, intent(in) :: sel(nw)          ! candidate columns used as trials
    complex(dp), intent(out) :: A_gauge(nw,nw,nk)
    real(dp), intent(out) :: worst_svmin
    complex(dp) :: svb(nw,nw)
    real(dp) :: sv
    integer :: k, m, n
    worst_svmin = huge(1.0_dp)
    do k = 1, nk
       do n = 1, nw
          do m = 1, nw
             A_gauge(m,n,k) = conjg(B(sel(n), m, k))
          end do
       end do
       svb = A_gauge(:,:,k)
       sv = smallest_sv(svb, nw)             ! svb destroyed
       if (sv < worst_svmin) worst_svmin = sv
    end do
  end subroutine lowdin_gauge

  ! Explicit per-k Loewdin orthonormalization: A <- A (A^H A)^{-1/2}.
  ! The symmetric orthonormalization is the unique one commuting with
  ! unitary group action (S transforms covariantly, hence S^{-1/2} does),
  ! so applying it to a group-closed projection gauge preserves exact
  ! covariance -- REQUIRED before smooth_gauge if covariance matters,
  ! because the first polar sweep would otherwise unitarize-and-mix.
  ! Eigenvalue floor mirrors lowdin_pstate's relative cutoff; a floored
  ! eigenvalue means the conditioning gate upstream should have fired.
  subroutine lowdin_orthon(A, nw, nk)
    integer, intent(in) :: nw, nk
    complex(dp), intent(inout) :: A(nw,nw,nk)
    complex(dp) :: S(nw,nw), T(nw,nw), Sinvh(nw,nw)
    real(dp) :: ev(nw), wcut
    integer :: k, i, j
    do k = 1, nk
       call zgemm('C','N', nw, nw, nw, CONE, A(1,1,k), nw, A(1,1,k), nw, &
                  CZERO, S, nw)
       call eigh(S, nw, ev)               ! S overwritten by eigenvectors V
       wcut = max(maxval(ev), 1.0e-300_dp) * 1.0e-12_dp
       do j = 1, nw
          do i = 1, nw
             T(i,j) = S(i,j) / sqrt(max(ev(j), wcut))
          end do
       end do
       ! S^{-1/2} = V diag(ev^{-1/2}) V^H ;  then  A <- A S^{-1/2}
       call zgemm('N','C', nw, nw, nw, CONE, T, nw, S, nw, CZERO, Sinvh, nw)
       T = A(:,:,k)
       call zgemm('N','N', nw, nw, nw, CONE, T, nw, Sinvh, nw, CZERO, &
                  A(1,1,k), nw)
    end do
  end subroutine lowdin_orthon

  ! b := A^-1 b for a small SPD system (Anderson Gram solve); on failure
  ! leaves b as the plain average weights
  subroutine solve_spd(A, b, n)
    integer, intent(in) :: n
    real(dp), intent(inout) :: A(n,n), b(n)
    integer :: info
    call dposv('U', n, 1, A, n, b, n, info)
    if (info /= 0) b = 1.0_dp/real(n,dp)
  end subroutine solve_spd

  ! ascending insertion sort of a short integer vector (anchor pivots)
  subroutine isort(a, n)
    integer, intent(in) :: n
    integer, intent(inout) :: a(n)
    integer :: i, j, key
    do i = 2, n
       key = a(i); j = i-1
       do while (j >= 1)
          if (a(j) <= key) exit
          a(j+1) = a(j); j = j-1
       end do
       a(j+1) = key
    end do
  end subroutine isort

  !--------------------------------------------------------- gauge dispatcher
  ! The full validated gauge-selection pipeline, factored out of the
  ! standalone driver VERBATIM (2026-07-23) so the in-w2w HYBR one-shot
  ! runs the identical code path: anchor gauge (scdm / block-SCDM /
  ! Loewdin trials) -> anchors sidecar + optional per-k sigma dump ->
  ! seed policy (W2W_SEEDK best|median|<int>|opt[<N>] gate-then-sample
  ! tournament, W2W_SEEDGATE, W2W_SEEDLOSS) -> mode dispatch
  ! (ptg / auto<N> / tfl<N> / sms<N> / pts<N> / lds<N> / bsms<N>) ->
  ! initial-spread report.  Env vars are read HERE so both callers behave
  ! identically; uout redirects the narration (default: stdout).
  subroutine choose_gauge(B, eig_t, mmn_rot, neigh, wgt, bcart, nproj, nw, &
                          nntot, nk, gauge_mode, outp, A_gauge, anchors, &
                          best_svmin, uout, seedk_in, seedloss_in, bkb)
    use hybrid_io_m, only: handoff_close
    integer, intent(in) :: nproj, nw, nntot, nk
    complex(dp), intent(in) :: B(nproj,nw,nk), mmn_rot(nw,nw,nntot,nk)
    real(dp), intent(in) :: eig_t(nw,nk), wgt(nntot), bcart(3,nntot)
    integer, intent(in) :: neigh(nntot,nk)
    character(*), intent(in) :: gauge_mode, outp
    complex(dp), intent(out) :: A_gauge(nw,nw,nk)
    integer, intent(out) :: anchors(nw)
    real(dp), intent(out) :: best_svmin
    integer, intent(in), optional :: uout
    ! seed policy defaults for callers that carry them as proper inputs
    ! (the HYBR keyword tail); the W2W_SEEDK/W2W_SEEDLOSS env vars still
    ! OVERRIDE when set, so experiments keep working on both paths
    character(*), intent(in), optional :: seedk_in, seedloss_in
    !> per-(k,b) b-vectors for every spread the gauge policy evaluates (the seed
    !> tournament's ranking and the narration); absent = k = 1's vectors at every k
    real(dp), intent(in), optional :: bkb(3,nntot,nk)
    integer :: uo, i
    !> ---- S3 W1/W2: the A5 stamp.  Collected DURING the run, emitted as '#'
    !> header lines by the .anchors writer (relocated BELOW the seed/transport
    !> policy so the branch actually taken is known at write time -- the list
    !> content is unchanged, so every existing consumer of the integers is
    !> unaffected; consumers are taught to skip '#').
    character(len=64) :: st_branch, st_fallback
    integer  :: st_sweeps, st_nfull
    real(dp) :: st_siggate, st_mineig0
    character(len=16) :: sgenv_l
    integer :: sgst_l
    uo = output_unit
    if (present(uout)) uo = uout
    st_branch = 'none'; st_fallback = 'none'
    st_sweeps = -1; st_nfull = -1; st_mineig0 = -1.0_dp
    ! W2W_SIGGATE: OVERRIDE for the auto/tfl/slwa/slwf conditioning gate.
    ! Default value is the historical hardcoded 2.0e-3, so behaviour is
    ! BIT-IDENTICAL unless the env is set -- the override exists for the A5
    ! forced-branch negative control, and its value is stamped either way.
    st_siggate = 2.0e-3_dp
    call get_environment_variable('W2W_SIGGATE', sgenv_l, status=sgst_l)
    if (sgst_l == 0 .and. len_trim(sgenv_l) > 0) then
       read(sgenv_l, *, iostat=sgst_l) st_siggate
       if (sgst_l /= 0) st_siggate = 2.0e-3_dp
       write(uo,'(a,es10.3,a)') '  W2W_SIGGATE override: sig_gate = ', &
            st_siggate, '   (stamped)'
    end if

    if (trim(gauge_mode) == 'u0') then
       ! ---- S3 W4 (D2): EXTERNAL SEED SOCKET.  gauge_mode='u0' loads the
       ! initial gauge from <outp>.gauge0 (unformatted stream: one header
       ! record nw, nk, neighbor-table hash; then A_gauge in driver k-order).
       ! The hash pins the mesh ORDER, the thing a file cannot otherwise
       ! prove.  Per-k unitarity <= 1e-8 on load or REFUSE -- never Loewdin-
       ! repair here: repair would hide the producer's defect (the aSn
       ! lesson: transport cannot fix a rank-deficient frame, and a loader
       ! that silently fixes files breeds broken producers).
       u0blk: block
         integer :: uu0, ios0, nwf, nkf, hshf, hsh, kk0, bb0
         complex(dp), allocatable :: buf0(:,:)
         real(dp) :: dev0, wdev0
         hsh = 0
         do kk0 = 1, nk
            do bb0 = 1, nntot
               hsh = mod(hsh*31 + neigh(bb0,kk0), 1073741789)
            end do
         end do
         open(newunit=uu0, file=trim(outp)//'.gauge0', status='old', &
              access='stream', form='unformatted', iostat=ios0)
         if (ios0 /= 0) then
            write(uo,'(a)') '  u0: cannot open '//trim(outp)//'.gauge0 -- '// &
                 'the external-seed socket needs its file.  REFUSING.'
            error stop 'u0: missing .gauge0'
         end if
         read(uu0, iostat=ios0) nwf, nkf, hshf
         if (ios0 /= 0 .or. nwf /= nw .or. nkf /= nk .or. hshf /= hsh) then
            write(uo,'(a,i0,a,i0,a,i0,a,i0)') '  u0: header mismatch -- '// &
                 'file nw=', nwf, ' nk=', nkf, '  run nw=', nw, ' nk=', nk
            write(uo,'(a)') '  u0: (or the neighbour-table hash differs: '// &
                 'the file was made for a DIFFERENT mesh/order).  REFUSING.'
            error stop 'u0: .gauge0 does not match this run'
         end if
         read(uu0, iostat=ios0) A_gauge
         close(uu0)
         if (ios0 /= 0) error stop 'u0: truncated .gauge0 stream'
         allocate(buf0(nw,nw))
         wdev0 = 0.0_dp
         do kk0 = 1, nk
            call zgemm('C','N', nw, nw, nw, CONE, A_gauge(1,1,kk0), nw, &
                       A_gauge(1,1,kk0), nw, CZERO, buf0, nw)
            do bb0 = 1, nw
               buf0(bb0,bb0) = buf0(bb0,bb0) - CONE
            end do
            dev0 = maxval(abs(buf0))
            wdev0 = max(wdev0, dev0)
         end do
         deallocate(buf0)
         if (wdev0 > 1.0e-8_dp) then
            write(uo,'(a,es10.3,a)') '  u0: worst per-k |U^H U - I| = ', &
                 wdev0, ' > 1e-8 -- NON-UNITARY external gauge.  REFUSING'
            write(uo,'(a)') '  u0: fix the PRODUCER; this loader never '// &
                 'repairs (repair would hide the defect).'
            error stop 'u0: non-unitary .gauge0'
         end if
         write(uo,'(a,es10.3)') '  u0: external gauge loaded, worst '// &
              'per-k |U^H U - I| = ', wdev0
         do kk0 = 1, nw
            anchors(kk0) = kk0
         end do
         best_svmin = 1.0_dp        ! unitary by the check above
         st_branch = 'u0-external'
       end block u0blk
    else if (gauge_mode(1:2) == 'bs') then
       ! block-SCDM (arg 14 = 'bscdm' | 'bsms<N>'): per-isolated-manifold
       ! anchors, block-diagonal init gauge.  0.5 eV gap => isolated shell.
       block
         integer :: nbl
         call block_scdm(B, eig_t, nproj, nk, nw, 0.5_dp, anchors, A_gauge, &
                         best_svmin, nbl, .true.)
         write(uo,'(a,i0,a,es12.4)') '  block-SCDM: ', nbl, &
              ' block(s), worst-block sigma_min = ', best_svmin
       end block
    else if (gauge_mode(1:5) == 'ident') then
       ! Unanchored control: A_gauge = identity per k, so W90 starts from
       ! rotate_to_subspace's projected-H eigenbasis with no anchoring
       ! information at all.  Every other gauge is SCDM-anchored in some form.
       ! No transport, no smoothing, no seed policy touches it below (prefix
       ! matches none of ptg/auto/tfl/sms/pts/lds).
       block
         integer :: kk, ii
         do ii = 1, nw
            anchors(ii) = ii
         end do
         A_gauge = CZERO
         do kk = 1, nk
            do ii = 1, nw
               A_gauge(ii,ii,kk) = CONE
            end do
         end do
         best_svmin = 1.0_dp
         write(uo,'(a)') '  gauge: identity (UNANCHORED control)'
       end block
    else if (gauge_mode(1:4) == 'slwd' .or. gauge_mode(1:4) == 'slwa' &
         .or. gauge_mode(1:4) == 'slwf' .or. gauge_mode(1:4) == 'slwp') then
       ! Loewdin projectability anchors with site-local closure.  Selection = P-ranked columns closed to whole
       ! (atom,l,rf) shells by the exact knapsack (the pshell rule); gauge =
       ! Loewdin projection columns on the CLOSED set (lowdin_gauge), which
       ! W90 orthonormalizes on hand-off -- Loewdin S^{-1/2} commutes with
       ! unitary transformations of the columns, so on a set closed over
       ! whole shells the gauge transforms with the shells.
       ! 'slwd<N>': explicit per-k S^{-1/2} orthonormalization FIRST, then N
       ! smoothing sweeps (never smooth a non-orthonormal frame: the first
       ! polar sweep would unitarize-and-mix before covariance is checked).
       ! Conditioning gate: min-eig(S) = sigma_min(A)^2 reported per the
       ! plan; on collapse the fallback ladder is slwd -> pts -> auto8
       ! (manual, stamped in .anchors -- no silent auto-fallback here).
       slwdblk: block
         integer, allocatable :: grp(:), prank(:)
         real(dp), allocatable :: pw(:)
         integer :: jj, ii, kk, nfull, nswl, ios_l
         call read_projmap_groups(nproj, grp, uo)
         if (.not. allocated(grp)) then
            write(uo,'(a)') '  slwd: no usable projmap.dat -- the pshell '// &
                 'selection needs shell groups.  FAILING (no silent fallback).'
            error stop 'slwd: projmap.dat with shell groups is required'
         end if
         ! P_j = (1/nk) sum_k ||B(j,:,k)||^2 -- only the diagonal of rho is
         ! needed for ranking, so build just that.
         allocate(pw(nproj), prank(nproj))
         pw = 0.0_dp
         do kk = 1, nk
            do ii = 1, nw
               do jj = 1, nproj
                  pw(jj) = pw(jj) + real(B(jj,ii,kk),dp)**2 + aimag(B(jj,ii,kk))**2
               end do
            end do
         end do
         pw = pw / real(nk,dp)
         do jj = 1, nproj
            prank(jj) = jj
         end do
         do ii = 1, nproj-1              ! insertion sort, descending P
            do jj = ii+1, nproj
               if (pw(prank(jj)) > pw(prank(ii))) then
                  kk = prank(ii); prank(ii) = prank(jj); prank(jj) = kk
               end if
            end do
         end do
         anchors = prank(1:nw)
         call groupify_pivots(anchors, nw, grp, nproj, nfull, uo)
         call isort(anchors, nw)
         write(uo,'(a,i0,a)') '  slwd anchors: P-ranked shell-closed set (', &
              nfull, ' whole shells)'
         call lowdin_gauge(B, nproj, nk, nw, anchors, A_gauge, best_svmin)
         write(uo,'(a,es12.4,a,es12.4)') '  slwd sigma_min(worst k) = ', &
              best_svmin, '   min-eig(S) = ', best_svmin**2
         st_nfull = nfull
         st_mineig0 = best_svmin**2
         ! ---- S3 W1: THE D1 CONDITIONING LADDER (fires on pure slwd only --
         ! slwa/slwf/slwp are explicit manual arms and keep their behaviour).
         ! Threshold 1e-10 = the geometric middle of the measured split
         ! (healthy 2.9e-6 vs collapsed 6.1e-20 / 2.2e-20, three cells; one
         ! future value inside [1e-16,1e-8] re-opens it -- plan risk R1).
         ! Rungs: slwd -> slwf (transport-fill: keep healthy per-k frames,
         ! BFS-fill the sick ones) -> auto8 (fresh SCDM anchors + gate/
         ! transport + 8 sweeps).  Executed AUTONOMOUSLY and STAMPED --
         ! the ladder is the autonomy mechanism (plan :279); the measured
         ! failure it prevents is aSn slwd 256.023138 vs auto8's 129.27
         ! basin.  W2W_ANCHORS_RELAX=1 disables it (forensics only, stamped).
         ladder: block
           real(dp), parameter :: LADDER_MINEIG = 1.0e-10_dp
           character(len=8) :: rlx
           logical, allocatable :: keepl(:)
           complex(dp), allocatable :: bufl(:,:)
           integer :: kl, nheall
           real(dp) :: svl
           if (gauge_mode(1:4) /= 'slwd') exit ladder
           if (best_svmin**2 >= LADDER_MINEIG) exit ladder
           call get_environment_variable('W2W_ANCHORS_RELAX', rlx)
           if (trim(rlx) == '1') then
              st_fallback = 'RELAX(ladder disabled)'
              write(uo,'(a)') '  ladder: min-eig(S) COLLAPSED but '// &
                   'W2W_ANCHORS_RELAX=1 -- keeping the raw slwd frame '// &
                   '(forensics only; stamped)'
              exit ladder
           end if
           write(uo,'(a,es10.3,a)') '  ladder: min-eig(S) ', best_svmin**2, &
                ' < 1.0E-10 -- rung 1: slwf transport-fill'
           allocate(keepl(nk), bufl(nw,nw))
           do kl = 1, nk
              bufl = A_gauge(:,:,kl)
              keepl(kl) = (smallest_sv(bufl, nw) >= st_siggate)
           end do
           nheall = count(keepl)
           if (nheall == 0) then
              kl = gauge_best_k(A_gauge, nw, nk)
              keepl(kl) = .true.; nheall = 1
           end if
           write(uo,'(a,i0,a,i0,a)') '  ladder: transport-fill ', nk-nheall, &
                ' of ', nk, ' rank-deficient k'
           if (nheall < nk) &
              call propagate_fill(mmn_rot, neigh, nw, nntot, nk, keepl, A_gauge)
           ! Re-measure the WORST conditioning BEFORE lowdin_orthon -- the
           ! plan's own D1 warning, which the first version of this ladder
           ! violated: orthon makes every frame unitary, so a post-orthon
           ! measurement returns 1.0 IDENTICALLY and the rung-2 test becomes
           ! unreachable by construction (measured: post-fill "min-eig 1.0"
           ! on aSn regardless of repair quality).
           svl = huge(1.0_dp)
           do kl = 1, nk
              bufl = A_gauge(:,:,kl)
              svl = min(svl, smallest_sv(bufl, nw))
           end do
           best_svmin = svl
           write(uo,'(a,es12.4,a,es12.4)') '  ladder: post-fill (pre-orthon)'// &
                ' sigma_min = ', svl, '   min-eig(S) = ', svl**2
           call lowdin_orthon(A_gauge, nw, nk)
           if (svl**2 >= LADDER_MINEIG) then
              st_fallback = 'ladder:slwd->slwf'
              deallocate(keepl, bufl)
              exit ladder
           end if
           ! rung 2: the frame is beyond repair (broad collapse -- smoothing
           ! and transport cannot fix rank deficiency, MEASURED: min-eig
           ! 2.2e-20 through EVERY sweep count).  Fall back to auto8: fresh
           ! SCDM anchors, the same gate/transport policy, 8 sweeps.
           write(uo,'(a)') '  ladder: still collapsed -- rung 2: auto8 '// &
                '(fresh SCDM anchors + gate/transport + 8 sweeps)'
           call scdm_anchors(B, nproj, nk, nw, 1, anchors, A_gauge, best_svmin)
           if (best_svmin < st_siggate) then
              kl = gauge_best_k(A_gauge, nw, nk)
              call propagate_gauge(mmn_rot, neigh, nw, nntot, nk, kl, A_gauge)
              write(uo,'(a,i0)') '  ladder: transport from k=', kl
              write(st_branch,'(a,i0)') 'transport:k', kl
           else
              st_branch = 'per-k-kept'
           end if
           call smooth_gauge(mmn_rot, neigh, wgt, nw, nntot, nk, 8, A_gauge)
           st_sweeps = 8
           st_fallback = 'ladder:slwd->slwf->auto8'
           deallocate(keepl, bufl)
           deallocate(grp, pw, prank)
           exit slwdblk
         end block ladder
         ! optional 'slwd<N>': explicit Loewdin orthonormalization + N sweeps.
         ! 'slwa<N>' additionally applies auto8's sigma GATE + TRANSPORT.
         !
         ! Measured on MgB2 12x12x8: smoothing
         ! alone takes slwd from 4309 to 2322 W90 iterations and improves
         ! every other axis (iter-0 spread 41.1 -> 34.6, off-mesh eta 136 ->
         ! 49 meV = 5x better than auto8, best basin) -- but cannot reach
         ! auto8's 437, because min-eig(S) stays at 2.2e-20 through EVERY
         ! sweep count: smoothing enforces LOCAL k-to-k coherence and cannot
         ! repair a rank-deficient frame.  auto8's speed comes from the other
         ! half of its recipe -- when conditioning fails it transports the
         ! frame from the best-conditioned k, which establishes GLOBAL
         ! coherence.  slwa = Loewdin anchors + that same repair.
         if (len_trim(gauge_mode) >= 5) then
            nswl = -1
            read(gauge_mode(5:), *, iostat=ios_l) nswl
            if (ios_l == 0 .and. nswl >= 0) then
               if (gauge_mode(1:4) == 'slwa') then
                  ! same gate constant as the auto path
                  if (best_svmin < st_siggate) then
                     kk = gauge_best_k(A_gauge, nw, nk)
                     call propagate_gauge(mmn_rot, neigh, nw, nntot, nk, kk, &
                                          A_gauge)
                     write(uo,'(a,es9.2,a,i0)') '  slwa: sigma_min ', &
                          best_svmin, ' < gate -> transport from k=', kk
                  else
                     write(uo,'(a,es9.2,a)') '  slwa: sigma_min ', best_svmin, &
                          ' >= gate -> per-k Loewdin frames kept'
                  end if
               else if (gauge_mode(1:4) == 'slwp') then
                  kk = gauge_best_k(A_gauge, nw, nk)
                  call align_gauge(mmn_rot, neigh, nw, nntot, nk, kk, A_gauge)
                  write(uo,'(a,i0)') &
                       '  slwp: BFS Procrustes alignment from k=', kk
               else if (gauge_mode(1:4) == 'slwf') then
                  ! SURGICAL repair (measured motivation, MgB2 m12s):
                  ! the Loewdin frames are BETTER conditioned
                  ! than auto8's at the typical k (median sigma_min 0.115 vs
                  ! 0.028, max 0.401 vs 0.097) but 14 of 1152 k (1.2%) are
                  ! rank-deficient (min 1.5e-10).  Global transport (slwa)
                  ! would discard all 1138 good frames to cure 14;
                  ! transport-FILL keeps the healthy ones and BFS-transports
                  ! into the sick ones only.
                  fillblk: block
                    logical, allocatable :: keepk(:)
                    complex(dp), allocatable :: buf(:,:)
                    integer :: kq, nheal
                    allocate(keepk(nk), buf(nw,nw))
                    do kq = 1, nk
                       buf = A_gauge(:,:,kq)
                       keepk(kq) = (smallest_sv(buf, nw) >= st_siggate)
                    end do
                    nheal = count(keepk)
                    if (nheal == 0) then
                       kq = gauge_best_k(A_gauge, nw, nk)
                       keepk(kq) = .true.; nheal = 1
                    end if
                    write(uo,'(a,i0,a,i0,a)') '  slwf: transport-fill ', &
                         nk-nheal, ' of ', nk, ' rank-deficient k'
                    if (nheal < nk) &
                       call propagate_fill(mmn_rot, neigh, nw, nntot, nk, &
                                           keepk, A_gauge)
                    deallocate(keepk, buf)
                  end block fillblk
               end if
               call lowdin_orthon(A_gauge, nw, nk)
               if (nswl > 0) &
                  call smooth_gauge(mmn_rot, neigh, wgt, nw, nntot, nk, nswl, &
                                    A_gauge)
               write(uo,'(a,a,a,i0,a)') '  ', trim(gauge_mode(1:4)), &
                    ': explicit S^{-1/2} + ', nswl, ' smoothing sweeps'
            end if
         end if
         deallocate(grp, pw, prank)
       end block slwdblk
    else
       ! env W2W_ANCHORS=shellcomplete: constrain the SCDM pivots to whole
       ! (atom,l,rf) shells x both sigma blocks, so the anchor set is a
       ! representation basis.  It costs spread, since it spends WF budget
       ! on orbitals SCDM declined.
       block
         character(len=16) :: aenv
         integer, allocatable :: grp(:)
         call get_environment_variable('W2W_ANCHORS', aenv)
         select case (trim(aenv))
         case ('shellcomplete', 'pshell')
            ! both need projmap groups; pshell additionally replaces the QRCP
            ! ranking with projectability ranking (see scdm_anchors)
            call read_projmap_groups(nproj, grp, uo)
            if (allocated(grp)) then
               if (trim(aenv) == 'pshell') then
                  write(uo,'(a)') '  anchors: P-ranked whole shells (pshell)'
                  call scdm_anchors(B, nproj, nk, nw, 1, anchors, A_gauge, &
                                    best_svmin, grp, amode='pshell')
               else
                  call scdm_anchors(B, nproj, nk, nw, 1, anchors, A_gauge, &
                                    best_svmin, grp)
               end if
               deallocate(grp)
            else
               ! S3 W1: the silent fallback was one of the four partial-orbit
               ! slip-throughs -- a run that ASKED for shell closure got
               ! unconstrained SCDM and nothing downstream could tell.
               block
                 character(len=8) :: rlx2
                 call get_environment_variable('W2W_ANCHORS_RELAX', rlx2)
                 if (trim(rlx2) /= '1') then
                    write(uo,'(a)') '  W2W_ANCHORS='//trim(aenv)//': no '// &
                         'usable projmap.dat.  Shell closure was REQUESTED '// &
                         'and cannot be delivered -- REFUSING (no silent '// &
                         'fallback).  Regenerate projmap.dat with the same '// &
                         'selection, or set W2W_ANCHORS_RELAX=1 (forensics '// &
                         'only; stamped).'
                    error stop 'W2W_ANCHORS: projmap.dat required for shell closure'
                 end if
                 write(uo,'(a)') '  W2W_ANCHORS='//trim(aenv)//': no usable '// &
                      'projmap.dat -- RELAX set, falling back to '// &
                      'unconstrained SCDM (stamped)'
                 st_fallback = 'RELAX(no-projmap SCDM)'
                 call scdm_anchors(B, nproj, nk, nw, 1, anchors, A_gauge, best_svmin)
               end block
            end if
         case ('pweight')
            write(uo,'(a)') '  anchors: projectability-weighted QRCP (pweight)'
            call scdm_anchors(B, nproj, nk, nw, 1, anchors, A_gauge, &
                              best_svmin, amode='pweight')
         case default
            call scdm_anchors(B, nproj, nk, nw, 1, anchors, A_gauge, best_svmin)
         end select
       end block
       write(uo,'(a,es12.4)') '  SCDM sigma_min(worst k) = ', best_svmin
    end if
    ! Loewdin projection gauge (arg 14 = 'lowdin' | 'lds<N>'): replace the
    ! QRCP anchor gauge with the leading nw physical channels as trials.
    if (gauge_mode(1:3) == 'low' .or. gauge_mode(1:3) == 'lds') then
       do i = 1, nw
          anchors(i) = i
       end do
       call lowdin_gauge(B, nproj, nk, nw, anchors, A_gauge, best_svmin)
       write(uo,'(a,es12.4)') '  Lowdin sigma_min(worst k) = ', best_svmin
    end if
    ! (the .anchors writer used to sit HERE, before the seed/transport policy;
    ! it is now AFTER that policy so the A5 stamp records the branch taken.
    ! The integer list itself is unchanged.)
    ! PTG seed k: k=1 by default; env W2W_SEEDK=best -> best-conditioned k
    ! (largest per-k sigma_min); W2W_SEEDK=<integer> -> explicit seed k
    ! (seed sweeps).  Only affects the transport (ptg/pts/auto) step.
    block
      integer :: k0seed, ios_k
      logical :: seed_opt
      character(len=8) :: seed_loss
      integer :: cand_seeds(65), ncand_s
      character(len=8) :: seedenv
      call get_environment_variable('W2W_SEEDK', seedenv)
      if (len_trim(seedenv) == 0 .and. present(seedk_in)) seedenv = seedk_in
      k0seed = 1
      seed_opt = (seedenv(1:3) == 'opt')
      call get_environment_variable('W2W_SEEDLOSS', seed_loss)
      if (len_trim(seed_loss) == 0 .and. present(seedloss_in)) &
           seed_loss = seedloss_in
      if (len_trim(seed_loss) == 0) seed_loss = 'smooth'
      if (trim(seed_loss) /= 'smooth') &
         write(uo,'(a)') '  seed loss mode: '//trim(seed_loss)
      ncand_s = 0
      if (seed_opt) then
         ! W2W_SEEDK=opt[<N>]: gate-then-sample candidate selection.
         ! EXCLUDE seeds failing the instability gate (per-k sigma_min <
         ! sig_gate: untrustworthy anchor), then sample up to N (default 12,
         ! max 64) candidates evenly across the sigma-sorted PASSING pool
         ! (endpoints included, so the best-sigma seed is always a
         ! candidate) + Gamma if it passes.  The tournament in hybrid_core
         ! then scores every candidate by the true MV spread and keeps the
         ! argmin -- maximal sampling of the stable pool for the basin
         ! search, no percentile heuristics.
         block
           integer :: kk5, info5, lw5, i5, j5, key5, kc, npass, nsamp, ip5
           real(dp) :: sgate5
           character(len=16) :: gate_env
           integer, allocatable :: idx5(:), pass5(:)
           complex(dp), allocatable :: Ac5(:,:), wk5(:)
           real(dp), allocatable :: sv5(:), rw5(:), svk5(:)
           complex(dp) :: qw5(1)
           ! env W2W_SEEDGATE: candidate stability gate (default 2e-3).
           ! With nsw=0 (no smoothing) the gate carries the full seed-noise
           ! safety burden -- raise it when running raw-transport tournaments.
           sgate5 = 2.0e-3_dp
           call get_environment_variable('W2W_SEEDGATE', gate_env)
           if (len_trim(gate_env) > 0) then
              read(gate_env, *, iostat=info5) sgate5
              if (info5 /= 0) sgate5 = 2.0e-3_dp
              write(uo,'(a,es10.3)') '  seed gate = ', sgate5
           end if
           nsamp = 12
           if (len_trim(seedenv) > 3) then
              read(seedenv(4:), *, iostat=ios_k) nsamp
              if (ios_k /= 0) nsamp = 12
           end if
           nsamp = max(2, min(64, nsamp))
           allocate(Ac5(nw,nw), sv5(nw), rw5(5*nw), svk5(nk), idx5(nk), pass5(nk))
           call zgesvd('N','N', nw, nw, Ac5, nw, sv5, Ac5, nw, Ac5, nw, &
                       qw5, -1, rw5, info5)
           lw5 = int(real(qw5(1),dp)); allocate(wk5(lw5))
           do kk5 = 1, nk
              Ac5 = A_gauge(:,:,kk5)
              call zgesvd('N','N', nw, nw, Ac5, nw, sv5, Ac5, nw, Ac5, nw, &
                          wk5, lw5, rw5, info5)
              if (info5 /= 0) error stop 'seed opt: SVD failed'
              svk5(kk5) = sv5(nw)
           end do
           do i5 = 1, nk; idx5(i5) = i5; end do
           do i5 = 2, nk                       ! sort k by sigma ascending
              key5 = idx5(i5); j5 = i5-1
              do while (j5 >= 1)
                 if (svk5(idx5(j5)) <= svk5(key5)) exit
                 idx5(j5+1) = idx5(j5); j5 = j5-1
              end do
              idx5(j5+1) = key5
           end do
           npass = 0                           ! stability-passing pool, sigma-sorted
           do i5 = 1, nk
              if (svk5(idx5(i5)) >= sgate5) then
                 npass = npass + 1; pass5(npass) = idx5(i5)
              end if
           end do
           if (npass == 0) then                ! pathological: nothing passes;
              npass = nk; pass5 = idx5         ! sample everything instead
              write(uo,'(a)') '  seed opt: NO seed passes the stability gate; sampling all k'
           end if
           do ip5 = 1, min(nsamp, npass)       ! even sample incl. both endpoints
              kc = pass5(1 + nint(real(ip5-1,dp)/real(max(1,min(nsamp,npass)-1),dp) &
                                  *real(npass-1,dp)))
              if (.not. any(cand_seeds(1:ncand_s) == kc)) then
                 ncand_s = ncand_s + 1; cand_seeds(ncand_s) = kc
              end if
           end do
           if (svk5(1) >= sgate5 .and. .not. any(cand_seeds(1:ncand_s) == 1)) then
              ncand_s = ncand_s + 1; cand_seeds(ncand_s) = 1   ! Gamma if stable
           end if
           write(uo,'(a,i0,a,i0,a,i0,a)') '  seed tournament: ', ncand_s, &
                ' candidates from ', npass, ' of ', nk, ' gate-passing k'
           deallocate(Ac5, sv5, rw5, wk5, svk5, idx5, pass5)
         end block
      else if (trim(seedenv) == 'best') then
         k0seed = gauge_best_k(A_gauge, nw, nk)
         write(uo,'(a,i0)') '  PTG seed k (best-conditioned) = ', k0seed
      else if (trim(seedenv) == 'median') then
         ! Fe seed-sweep 2026-07-21: the HIGHEST-sigma seeds are atypical k
         ! (rigid templates -- both Fe-up failures were the top-2 sigma
         ! seeds), while median-sigma seeds converged fast on BOTH spins.
         ! Median = stable enough + representative.
         block
           integer :: kk4, info4, lw4, idx4(nk)
           complex(dp), allocatable :: Ac4(:,:), wk4(:)
           real(dp), allocatable :: sv4(:), rw4(:), svk4(:)
           complex(dp) :: qw4(1)
           integer :: i4, j4, key4
           allocate(Ac4(nw,nw), sv4(nw), rw4(5*nw), svk4(nk))
           call zgesvd('N','N', nw, nw, Ac4, nw, sv4, Ac4, nw, Ac4, nw, &
                       qw4, -1, rw4, info4)
           lw4 = int(real(qw4(1),dp)); allocate(wk4(lw4))
           do kk4 = 1, nk
              Ac4 = A_gauge(:,:,kk4)
              call zgesvd('N','N', nw, nw, Ac4, nw, sv4, Ac4, nw, Ac4, nw, &
                          wk4, lw4, rw4, info4)
              if (info4 /= 0) error stop 'seedk median: SVD failed'
              svk4(kk4) = sv4(nw)
           end do
           do i4 = 1, nk; idx4(i4) = i4; end do
           do i4 = 2, nk        ! insertion sort by svk4
              key4 = idx4(i4); j4 = i4-1
              do while (j4 >= 1)
                 if (svk4(idx4(j4)) <= svk4(key4)) exit
                 idx4(j4+1) = idx4(j4); j4 = j4-1
              end do
              idx4(j4+1) = key4
           end do
           k0seed = idx4((nk+1)/2)
           write(uo,'(a,i0,a,es10.3)') '  PTG seed k (median-conditioned) = ', &
                k0seed, '   sigma = ', svk4(k0seed)
           deallocate(Ac4, sv4, rw4, wk4, svk4)
         end block
      else if (len_trim(seedenv) > 0) then
         read(seedenv, *, iostat=ios_k) k0seed
         if (ios_k /= 0 .or. k0seed < 1 .or. k0seed > nk) then
            write(uo,'(a)') '  WARNING: bad W2W_SEEDK '''//trim(seedenv)// &
                 ''' -> seed k=1'
            k0seed = 1
         else
            write(uo,'(a,i0)') '  PTG seed k (user) = ', k0seed
         end if
      end if
      if (trim(gauge_mode) == 'ptg') then
         call propagate_gauge(mmn_rot, neigh, nw, nntot, nk, k0seed, A_gauge)
         write(uo,'(a,i0)') '  gauge: parallel-transport propagated from k=', k0seed
         write(st_branch,'(a,i0)') 'ptg:k', k0seed
      else if (gauge_mode(1:4) == 'auto' .or. gauge_mode(1:3) == 'tfl') then
         ! 'auto<N>': sigma_min-gated switch (benchmark 2026-07-21, 8/8):
         !   worst-k SCDM sigma_min < sig_gate -> per-k anchors are corrupt
         !   somewhere -> single-seed transport from the best-conditioned k;
         !   else keep the per-k SCDM gauge.  Then N smoothing sweeps.
         ! 'tfl<N>': transport-fill middle ground -- keep per-k SCDM at every
         !   HEALTHY k (per-k sigma_min >= sig_gate), BFS-transport into the
         !   sick k only, then N smoothing sweeps.
         autotfl: block
           integer :: nsw2, ios_g2, kk2, npos2, nhealth, info_s, lws
           real(dp) :: sig_gate
           logical, allocatable :: keepk(:)
           complex(dp), allocatable :: Ac(:,:), wrk(:)
           real(dp), allocatable :: svv(:), rwk(:)
           complex(dp) :: qws(1)
           sig_gate = st_siggate
           nsw2 = 8
           npos2 = merge(5, 4, gauge_mode(1:4) == 'auto')
           if (len_trim(gauge_mode) >= npos2) then
              read(gauge_mode(npos2:), *, iostat=ios_g2) nsw2
              if (ios_g2 /= 0) nsw2 = 8
           end if
           if (gauge_mode(1:4) == 'auto') then
              if (best_svmin < sig_gate) then
                 if (seed_opt) then
                    block
                      integer :: k_win
                      real(dp) :: om_win
                      call seed_tournament(mmn_rot, neigh, wgt, bcart, nw, nntot, &
                           nk, nsw2, cand_seeds(1:ncand_s), ncand_s, A_gauge, &
                           k_win, om_win, lossmode=seed_loss, bkb=bkb)
                      write(uo,'(a,i0,a,f14.4)') &
                           '  auto gauge (tournament): winner seed k=', k_win, &
                           '  Omega = ', om_win
                      write(st_branch,'(a,i0,a,f0.4)') 'tournament:k', &
                           k_win, ':Omega', om_win
                    end block
                    write(uo,'(a,a,a,i0,a)') '  gauge: ', trim(gauge_mode), &
                         ' with ', nsw2, ' smoothing sweeps (in tournament)'
                    exit autotfl
                 end if
                 k0seed = gauge_best_k(A_gauge, nw, nk)
                 call propagate_gauge(mmn_rot, neigh, nw, nntot, nk, k0seed, A_gauge)
                 write(uo,'(a,es9.2,a,es8.1,a,i0)') '  auto gauge: sigma_min ', &
                      best_svmin, ' < ', sig_gate, ' -> transport from k=', k0seed
                 write(st_branch,'(a,i0)') 'transport:k', k0seed
              else
                 write(uo,'(a,es9.2,a,es8.1,a)') '  auto gauge: sigma_min ', &
                      best_svmin, ' >= ', sig_gate, ' -> per-k SCDM kept'
                 st_branch = 'per-k-kept'
              end if
           else
              allocate(keepk(nk), Ac(nw,nw), svv(nw), rwk(5*nw))
              call zgesvd('N','N', nw, nw, Ac, nw, svv, Ac, nw, Ac, nw, &
                          qws, -1, rwk, info_s)
              lws = int(real(qws(1),dp)); allocate(wrk(lws))
              do kk2 = 1, nk
                 Ac = A_gauge(:,:,kk2)
                 call zgesvd('N','N', nw, nw, Ac, nw, svv, Ac, nw, Ac, nw, &
                             wrk, lws, rwk, info_s)
                 if (info_s /= 0) error stop 'tfl: per-k SVD failed'
                 keepk(kk2) = svv(nw) >= sig_gate    ! smallest sv (descending order)
              end do
              nhealth = count(keepk)
              if (nhealth == 0) then
                 k0seed = gauge_best_k(A_gauge, nw, nk)
                 keepk(k0seed) = .true.
                 nhealth = 1
              end if
              write(uo,'(a,i0,a,i0,a)') '  transport-fill: ', nk-nhealth, &
                   ' of ', nk, ' k repaired from healthy anchors'
              if (nhealth < nk) &
                 call propagate_fill(mmn_rot, neigh, nw, nntot, nk, keepk, A_gauge)
              deallocate(keepk, Ac, svv, rwk, wrk)
           end if
           call smooth_gauge(mmn_rot, neigh, wgt, nw, nntot, nk, nsw2, A_gauge)
           write(uo,'(a,a,a,i0,a)') '  gauge: ', trim(gauge_mode), ' with ', nsw2, &
                ' smoothing sweeps'
           st_sweeps = nsw2
           if (gauge_mode(1:3) == 'tfl') st_branch = 'transport-fill'
         end block autotfl
      else if (gauge_mode(1:3) == 'sms' .or. gauge_mode(1:3) == 'pts' &
           .or. gauge_mode(1:3) == 'lds' .or. gauge_mode(1:4) == 'bsms') then
         ptsblk: block
           integer :: nsw, ios_g, npos
           nsw = 8
           npos = merge(5, 4, gauge_mode(1:4) == 'bsms')   ! bsms<N> vs xxx<N>
           if (len_trim(gauge_mode) >= npos) then
              read(gauge_mode(npos:), *, iostat=ios_g) nsw
              if (ios_g /= 0) nsw = 8
           end if
           if (gauge_mode(1:3) == 'pts' .and. seed_opt) then
              block
                integer :: k_win
                real(dp) :: om_win
                call seed_tournament(mmn_rot, neigh, wgt, bcart, nw, nntot, nk, &
                     nsw, cand_seeds(1:ncand_s), ncand_s, A_gauge, k_win, om_win, &
                     lossmode=seed_loss, bkb=bkb)
                write(uo,'(a,i0,a,f14.4)') '  seed tournament winner: k=', k_win, &
                     '  Omega = ', om_win
                write(st_branch,'(a,i0,a,f0.4)') 'tournament:k', k_win, &
                     ':Omega', om_win
              end block
              write(uo,'(a,a,a,i0,a)') '  gauge: ', trim(gauge_mode), ' with ', &
                   nsw, ' smoothing sweeps (in tournament)'
              exit ptsblk
           end if
           if (gauge_mode(1:3) == 'pts') then
                call propagate_gauge(mmn_rot, neigh, nw, nntot, nk, k0seed, A_gauge)
                write(st_branch,'(a,i0)') 'pts:k', k0seed
           end if
           call smooth_gauge(mmn_rot, neigh, wgt, nw, nntot, nk, nsw, A_gauge)
           st_sweeps = nsw
           write(uo,'(a,a,a,i0,a)') '  gauge: ', trim(gauge_mode), ' with ', nsw, &
                ' smoothing sweeps'
         end block ptsblk
      end if
    end block
    ! ---- S3 W4: the ROUND-TRIP PRODUCER.  W2W_GAUGE0OUT=1 writes the final
    ! gauge to <outp>.gauge0 in exactly the format the u0 loader reads, so
    ! A4's round trip is: run any mode with this env, rerun with
    ! gauge_mode=u0, and the hand-off must be BIT-identical (md5 on .amn) --
    ! a stronger instrument than the plan's 1e-9 iter-0 spread.
    g0out: block
      character(len=16) :: g0e
      integer :: ug0, kk9, bb9, hsh9
      call get_environment_variable('W2W_GAUGE0OUT', g0e)
      if (trim(g0e) == '1') then
         hsh9 = 0
         do kk9 = 1, nk
            do bb9 = 1, nntot
               hsh9 = mod(hsh9*31 + neigh(bb9,kk9), 1073741789)
            end do
         end do
         open(newunit=ug0, file=trim(outp)//'.gauge0', status='replace', &
              access='stream', form='unformatted')
         write(ug0) nw, nk, hsh9
         write(ug0) A_gauge
         call handoff_close(ug0, 'the .gauge0 file')
         write(uo,'(a)') '  wrote external-seed gauge -> '//trim(outp)// &
              '.gauge0  (u0 round-trip producer)'
      end if
    end block g0out
    ! ---- S3 W1/W2: the stamped .anchors writer (A5).  '#' lines first, then
    ! the same bare integer list as always.  Anchor -> candidate-column map
    ! joins with projmap.dat to place REAL guiding centres in the .win.
    stampw: block
      integer :: ua
      character(len=64) :: sk_e, sl_e, stf_e
      call get_environment_variable('W2W_SEEDK', sk_e)
      if (len_trim(sk_e) == 0 .and. present(seedk_in)) sk_e = seedk_in
      if (len_trim(sk_e) == 0) sk_e = 'default(k=1)'
      call get_environment_variable('W2W_SEEDLOSS', sl_e)
      if (len_trim(sl_e) == 0 .and. present(seedloss_in)) sl_e = seedloss_in
      if (len_trim(sl_e) == 0) sl_e = 'smooth'
      call get_environment_variable('W2W_STRUCT', stf_e)
      if (len_trim(stf_e) == 0) stf_e = 'none'
      open(newunit=ua, file=trim(outp)//'.anchors', status='replace')
      write(ua,'(a)') '# w2w.anchors v2 (S3 W1 stamp; parsers skip # lines)'
      write(ua,'(a,a,a,i0,a,i0)') '# gauge_mode=', trim(gauge_mode), &
           '  nw=', nw, '  nproj=', nproj
      write(ua,'(a,es13.6,a,es13.6)') '# best_svmin=', best_svmin, &
           '  min_eig_S=', best_svmin**2
      if (st_mineig0 >= 0.0_dp) &
           write(ua,'(a,es13.6)') '# min_eig_S_pre_policy=', st_mineig0
      write(ua,'(a,es10.3,a,a)') '# sig_gate=', st_siggate, &
           '  branch=', trim(st_branch)
      write(ua,'(a,i0)') '# sweeps=', st_sweeps
      write(ua,'(a,a,a,a)') '# seedk=', trim(sk_e), '  seedloss=', trim(sl_e)
      write(ua,'(a,a)') '# fallback=', trim(st_fallback)
      if (st_nfull >= 0) write(ua,'(a,i0,a,i0)') '# nfull=', st_nfull, &
           '  of_nw=', nw
      write(ua,'(a,a)') '# struct=', trim(stf_e)
      do i = 1, nw
         write(ua,'(i6)') anchors(i)
      end do
      call handoff_close(ua, 'the .anchors file')
    end block stampw
    block
      real(dp) :: om_g, om_p(3)
      call gauge_spread(mmn_rot, neigh, wgt, bcart, nw, nntot, nk, A_gauge, &
                        om_g, om_p, bkb=bkb)
      write(uo,'(a,f14.4,a,3f12.4,a)') '  initial spread Omega = ', om_g, &
           '   (I/OD/D = ', om_p, ' )'
    end block
  end subroutine choose_gauge

end module hybrid_core_m
