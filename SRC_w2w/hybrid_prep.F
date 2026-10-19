!!! hybrid_prep.F -- preparation stage for the hybrid Stage-2 port:
!!! Marzari-Vanderbilt finite-difference b-weights (mv_weights) and the
!!! mask-construction pipeline (rescale/degeneracy-symmetrize/auto-regime/
!!! auto-trust/select_pool_and_frozen/froz-emax).  Ported from spread.py
!!! (mv_weights) and hybrid.py (masks) of lcao2wannier.  Consumed by both the
!!! standalone tool and the in-w2w hook.
!!!
!!! References (keys of REFERENCES.md in the root directory): b-vector weights
!!! [Marzari1997]; projectability and its use for freezing [Qiao2023];
!!! Loewdin-orthonormalized projectability [Lowdin1950]; the frozen ceiling,
!!! the thresholds and the pocket rule [Comaskey2026].
module hybrid_prep_m
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: mv_weights, build_masks, rescale_symmetrize, auto_regime, &
            lowdin_pstate, frozen_report, &
            froz_report_line, enforce_frozen_floor, pool_floor, pool_demand, &
            pocket_ceiling

  complex(dp), parameter :: LP_CONE  = (1.0_dp, 0.0_dp)
  complex(dp), parameter :: LP_CZERO = (0.0_dp, 0.0_dp)

  ! froz_emax = auto walk-down ladder.  These ARE HEAD's literals (8.0, 1.0,
  ! 1.0), lifted to named constants only so the refusal message can quote the
  ! ladder it searched instead of hard-coding it a second time.
  real(dp), parameter :: FE_WALK_HI = 8.0_dp, FE_WALK_LO = 1.0_dp, &
                         FE_WALK_STEP = 1.0_dp
  ! froz_emax = auto walk start: .false. = FE_WALK_HI (the fidelity ceiling,
  ! v3.0 behaviour); .true. = max(FE_WALK_HI, top of the trusted frozen set),
  ! the rungs above FE_WALK_HI joining the integer ladder at FE_WALK_HI.
  ! Same rule as lcao2wannier resolve_frozen_ceiling (start 'fid' / 'top').
  logical, parameter :: FE_WALK_FROM_TOP = .false.
  ! Pocket skip (2026-10-04), the rule of lcao2wannier hybrid.pocket_ceiling.
  ! A pool band that crosses the walked ceiling but lies below it at fewer than
  ! POCKET_FRAC of the k-points is a ceiling pocket. When the weakest hand-off
  ! link (smallest singular value of W^H M W) is below POCKET_SV and touches a
  ! k-point where a pocket band of mixed character (mean P of its sub-ceiling
  ! states >= POCKET_P_MIX) is frozen, the ceiling drops POCKET_MARGIN below
  ! that band and the disentanglement is repeated; the caller keeps it only if
  ! the new weakest link is >= POCKET_SV. LCAO SnTe (46 WFs): 0.20 -> 0.92,
  ! 8.4/24.4 -> 3.9/15.4 meV. LAPW alpha-Sn: the weakest link (0.24) sits on a
  ! pocket of P 0.25; skipping it by hand doubled the conduction error, and the
  ! mixed-character test leaves it alone. No other release cell fires.
  real(dp), parameter, public :: POCKET_SV = 0.5_dp, POCKET_FRAC = 0.25_dp, &
                                 POCKET_P_MIX = 0.5_dp, POCKET_MARGIN = 0.1_dp
  integer, parameter, public :: POCKET_MAX_STEPS = 3


contains

  !--------------------------------------------- pocket skip
  !> The lowered ceiling fe_new (fired = .true.) when the weakest hand-off link
  !> sits on a mixed-character ceiling pocket.  eps(npool,nk): pool energies, eV
  !> rel. E_F; Pp(npool,nk): their projectability; fe: the walked ceiling;
  !> sv(nntot,nk): smallest singular value of every hand-off link.  msg names
  !> the pocket bands.  Never lowers below FE_WALK_LO.
  subroutine pocket_ceiling(eps, Pp, npool, nk, fe, sv, neigh, nntot, fe_new, fired, msg)
    integer, intent(in) :: npool, nk, nntot, neigh(nntot,nk)
    real(dp), intent(in) :: eps(npool,nk), Pp(npool,nk), fe, sv(nntot,nk)
    real(dp), intent(out) :: fe_new
    logical, intent(out) :: fired
    character(*), intent(out) :: msg
    integer :: n, k, b, k0, b0, kb0, loc(2), nbl
    logical :: below(nk), weak_k(nk), touch
    real(dp) :: f, pm, emin
    character(len=96) :: part
    fired = .false.; fe_new = fe; msg = ''
    loc = minloc(sv); b0 = loc(1); k0 = loc(2)
    if (sv(b0,k0) >= POCKET_SV) return
    kb0 = neigh(b0,k0)
    weak_k = .false.
    do k = 1, nk
       do b = 1, nntot
          if (sv(b,k) < POCKET_SV) then
             weak_k(k) = .true.; weak_k(neigh(b,k)) = .true.
          end if
       end do
    end do
    touch = .false.
    emin = huge(1.0_dp)
    do n = 1, npool
       below = eps(n,:) <= fe + 1.0e-9_dp
       nbl = count(below)
       if (nbl == 0) cycle
       f = real(nbl,dp) / real(nk,dp)
       if (f >= POCKET_FRAC) cycle
       pm = sum(Pp(n,:), mask=below) / real(nbl,dp)
       if (pm < POCKET_P_MIX) cycle
       if (below(k0) .or. below(kb0)) touch = .true.
       if (any(below .and. weak_k)) then
          emin = min(emin, minval(eps(n,:)))
          write(part,'(a,i0,a,f8.3,a,f5.1,a,f5.2,a)') ' pool band ', n, ' (min E_F', &
               minval(eps(n,:)), ' eV, ', 100.0_dp*f, '% of k, P ', pm, ')'
          msg = trim(msg)//trim(part)
       end if
    end do
    if (.not. touch .or. emin == huge(1.0_dp)) then
       msg = ''; return
    end if
    if (emin - POCKET_MARGIN < FE_WALK_LO) then
       msg = ''; return
    end if
    fe_new = emin - POCKET_MARGIN
    fired = .true.
  end subroutine pocket_ceiling

  !--------------------------------------------- Loewdin per-state projectability
  ! p_i(k) = diag(A S^-1 A^dag),  S = A^dag A over the selected columns.
  ! Exactly bounded [0,1] -- the proper normalisation for overlapping
  ! (tail-extended full-cell) projectors, and the LCAO-native P the mask
  ! thresholds were tuned on.  Shared by the standalone tool
  ! (rescale='lowdin') and the in-w2w HYBR one-shot.
  subroutine lowdin_pstate(amnP, nb, nk, nprojP, pcol, P)
    integer, intent(in) :: nb, nk, nprojP
    complex(dp), intent(in) :: amnP(nb, nprojP, nk)
    logical, intent(in) :: pcol(nprojP)
    real(dp), intent(out) :: P(nb, nk)
    integer :: i, jj, kk, nc, info_l, lworkc
    integer, allocatable :: cidx(:)
    complex(dp), allocatable :: Ak(:,:), Sk(:,:), Yk(:,:), workc(:)
    real(dp), allocatable :: wev(:), rwork(:)
    real(dp) :: wcut
    nc = count(pcol)
    allocate(cidx(nc))
    jj = 0
    do i = 1, nprojP
       if (pcol(i)) then; jj = jj+1; cidx(jj) = i; end if
    end do
    lworkc = max(1, 2*nc-1)*8
    allocate(Ak(nb,nc), Sk(nc,nc), Yk(nb,nc), wev(nc), &
             rwork(max(1,3*nc-2)), workc(lworkc))
    do kk = 1, nk
       do jj = 1, nc
          Ak(:,jj) = amnP(:,cidx(jj),kk)
       end do
       call zgemm('C','N', nc, nc, nb, LP_CONE, Ak, nb, Ak, nb, LP_CZERO, Sk, nc)
       call zheev('V','U', nc, Sk, nc, wev, workc, lworkc, rwork, info_l)
       if (info_l /= 0) error stop 'lowdin P: zheev failed'
       call zgemm('N','N', nb, nc, nc, LP_CONE, Ak, nb, Sk, nc, LP_CZERO, Yk, nb)
       wcut = max(wev(nc), tiny(1.0_dp))*1.0e-10_dp
       do i = 1, nb
          P(i,kk) = 0.0_dp
          do jj = 1, nc
             if (wev(jj) > wcut) P(i,kk) = P(i,kk) + &
                  (real(Yk(i,jj),dp)**2 + aimag(Yk(i,jj))**2)/wev(jj)
          end do
       end do
    end do
    deallocate(cidx, Ak, Sk, Yk, wev, rwork, workc)
  end subroutine lowdin_pstate

  !------------------------------------------------------------------- helpers
  ! ascending-value sort returning the permutation idx (insertion sort; n small)
  subroutine argsort_asc(x, n, idx)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)
    integer, intent(out) :: idx(n)
    integer :: i, j, key
    do i = 1, n; idx(i) = i; end do
    do i = 2, n
       key = idx(i); j = i-1
       do while (j >= 1)
          if (x(idx(j)) <= x(key)) exit
          idx(j+1) = idx(j); j = j-1
       end do
       idx(j+1) = key
    end do
  end subroutine argsort_asc

  ! p-th percentile with linear interpolation (numpy 'linear' default)
  function percentile(x, n, p) result(q)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n), p
    real(dp) :: q, pos, frac, xs(n)
    integer :: idx(n), lo, hi, i
    call argsort_asc(x, n, idx)
    do i = 1, n; xs(i) = x(idx(i)); end do
    pos = (p/100.0_dp)*real(n-1, dp)      ! 0-based position
    lo = int(pos); hi = min(lo+1, n-1); frac = pos - real(lo, dp)
    q = xs(lo+1) + frac*(xs(hi+1) - xs(lo+1))   ! +1 -> 1-based
  end function percentile

  ! MV b-weights: solve the shell least-squares  sum_b w_b b_a b_c = delta_ac.
  ! bcart(3,nntot) are the Cartesian b-vectors (1/Ang).  Returns w(nntot) and
  ! the residual max|sum_b w_b b_a b_c - delta_ac| (diagnostic).
  subroutine mv_weights(bcart, nntot, w, resid, tol)
    integer, intent(in) :: nntot
    real(dp), intent(in) :: bcart(3,nntot)
    real(dp), intent(out) :: w(nntot), resid
    real(dp), intent(in), optional :: tol
    real(dp) :: shtol, bn(nntot), shval(nntot), b(3)
    integer :: shell(nntot), ns, i, s, a, c, mm
    real(dp), allocatable :: Mat(:,:), rhs(:), sv(:), work(:)
    integer, allocatable :: iwork(:)
    real(dp) :: qw(1), mom(3,3), tgt
    integer :: rank, info, lwork, liwork, nlvl

    shtol = 1.0e-6_dp; if (present(tol)) shtol = tol
    ! ---- group b-vectors into shells of equal length ----
    ns = 0
    do i = 1, nntot
       bn(i) = sqrt(sum(bcart(:,i)**2))
       s = 0
       do a = 1, ns
          if (abs(bn(i) - shval(a)) < shtol) then; s = a; exit; end if
       end do
       if (s == 0) then; ns = ns + 1; shval(ns) = bn(i); s = ns; end if
       shell(i) = s
    end do

    ! ---- Mat(6,ns): second-moment components per shell ----
    allocate(Mat(6,ns))
    Mat = 0.0_dp
    do i = 1, nntot
       b = bcart(:,i); s = shell(i)
       Mat(1,s) = Mat(1,s) + b(1)*b(1)
       Mat(2,s) = Mat(2,s) + b(2)*b(2)
       Mat(3,s) = Mat(3,s) + b(3)*b(3)
       Mat(4,s) = Mat(4,s) + b(1)*b(2)
       Mat(5,s) = Mat(5,s) + b(2)*b(3)
       Mat(6,s) = Mat(6,s) + b(3)*b(1)
    end do

    ! ---- least-squares  Mat wshell = [1,1,1,0,0,0]  (dgelsd, SVD; = numpy lstsq)
    mm = max(6, ns)
    allocate(rhs(mm), sv(min(6,ns)))
    rhs = 0.0_dp; rhs(1:3) = 1.0_dp
    nlvl = max(0, int(log(real(min(6,ns),dp)/26.0_dp)/log(2.0_dp)) + 1)
    allocate(iwork(max(1, 3*min(6,ns)*nlvl + 11*min(6,ns))))
    call dgelsd(6, ns, 1, Mat, 6, rhs, mm, sv, -1.0_dp, rank, qw, -1, iwork, info)
    lwork = int(qw(1)); allocate(work(lwork))
    call dgelsd(6, ns, 1, Mat, 6, rhs, mm, sv, -1.0_dp, rank, work, lwork, iwork, info)
    if (info /= 0) error stop 'mv_weights: dgelsd failed'
    do i = 1, nntot
       w(i) = rhs(shell(i))
    end do

    ! ---- residual (diagnostic) ----
    mom = 0.0_dp
    do i = 1, nntot
       do a = 1, 3
          do c = 1, 3
             mom(a,c) = mom(a,c) + w(i)*bcart(a,i)*bcart(c,i)
          end do
       end do
    end do
    resid = 0.0_dp
    do a = 1, 3
       do c = 1, 3
          tgt = merge(1.0_dp, 0.0_dp, a == c)
          resid = max(resid, abs(mom(a,c) - tgt))
       end do
    end do
    deallocate(Mat, rhs, sv, iwork, work)
  end subroutine mv_weights

  !------------------------------------------------------- auto-regime detection
  subroutine auto_regime(eig, nb, nk, dtol, degen_frac_thr, band_level, use_sym, degen_frac)
    integer, intent(in) :: nb, nk
    real(dp), intent(in) :: eig(nb,nk), dtol, degen_frac_thr
    logical, intent(out) :: band_level, use_sym
    real(dp), intent(out) :: degen_frac
    real(dp) :: es(nb)
    integer :: k, i, idx(nb), cnt
    cnt = 0
    do k = 1, nk
       call argsort_asc(eig(:,k), nb, idx)
       do i = 1, nb; es(i) = eig(idx(i),k); end do
       do i = 1, nb-1
          if (es(i+1)-es(i) < dtol) cnt = cnt + 1
       end do
    end do
    degen_frac = real(cnt*2, dp)/real(nk*nb, dp)
    band_level = degen_frac > degen_frac_thr
    use_sym = band_level
  end subroutine auto_regime

  !-------------------------------------------------- rescale + degen-symmetrize
  subroutine rescale_symmetrize(P, eig, nb, nk, e_fermi, rescale_mode, use_sym, dtol)
    integer, intent(in) :: nb, nk
    real(dp), intent(inout) :: P(nb,nk)
    real(dp), intent(in) :: eig(nb,nk), e_fermi, dtol
    character(*), intent(in) :: rescale_mode
    logical, intent(in) :: use_sym
    real(dp) :: scale, occ(nb*nk), allp(nb*nk), es(nb), grpmean
    integer :: k, i, j, m, idx(nb), nocc
    if (rescale_mode == 'occ99') then
       nocc = 0
       do k = 1, nk; do i = 1, nb
          if (eig(i,k)-e_fermi <= 0.0_dp) then; nocc = nocc+1; occ(nocc) = P(i,k); end if
       end do; end do
       scale = percentile(occ(1:nocc), nocc, 99.0_dp)
    else if (rescale_mode == 'p99') then
       m = 0; do k = 1, nk; do i = 1, nb; m = m+1; allp(m) = P(i,k); end do; end do
       scale = percentile(allp, nb*nk, 99.0_dp)
    else if (rescale_mode == 'max') then
       scale = maxval(P)
    else
       scale = 1.0_dp
    end if
    P = max(0.0_dp, min(1.0_dp, P/scale))
    if (use_sym) then
       do k = 1, nk
          call argsort_asc(eig(:,k), nb, idx)
          i = 1
          do while (i <= nb)
             j = i
             do while (j+1 <= nb)
                if (eig(idx(j+1),k) - eig(idx(i),k) < dtol) then; j = j+1; else; exit; end if
             end do
             if (j > i) then
                grpmean = 0.0_dp
                do m = i, j; grpmean = grpmean + P(idx(m),k); end do
                grpmean = grpmean/real(j-i+1, dp)
                do m = i, j; P(idx(m),k) = grpmean; end do
             end if
             i = j+1
          end do
       end do
    end if
  end subroutine rescale_symmetrize

  !------------------------------------------------------- contiguous pool select
  subroutine select_pool(P, eig, nb, nk, nw, e_fermi, pool_factor, seed_lo, seed_hi, &
                         pool, npool)
    integer, intent(in) :: nb, nk, nw
    real(dp), intent(in) :: P(nb,nk), eig(nb,nk), e_fermi, pool_factor, seed_lo, seed_hi
    integer, allocatable, intent(out) :: pool(:)
    integer, intent(out) :: npool
    real(dp) :: avg_p(nb), avg_e(nb), win_p(nb)
    integer :: i, k, nwin, winlist(nb), idx(nb), seed(nw), s0, s1, n_pool_max, j
    do i = 1, nb
       avg_p(i) = sum(P(i,:))/real(nk, dp)
       avg_e(i) = sum(eig(i,:))/real(nk, dp) - e_fermi
    end do
    ! *** DUPLICATED IN main.F Gate B (2/2) `npool_want'.  KEEP IN STEP. ***
    ! This routine is print-free (the standalone shares it), so the ratio
    ! check lives in the caller and restates this expression.  If the sizing
    ! rule here changes, main.F's poolgate block reports a ratio against a
    ! target that no longer exists -- silently.  Cross-referenced 2026-08-26.
    n_pool_max = pool_floor(nw, pool_factor)
    nwin = 0
    do i = 1, nb
       if (avg_e(i) >= seed_lo .and. avg_e(i) <= seed_hi) then
          nwin = nwin+1; winlist(nwin) = i
       end if
    end do
    if (nwin < nw) then
       nwin = nb; do i = 1, nb; winlist(i) = i; end do
    end if
    ! nw most-projectable in-window bands (largest avg_p)
    do i = 1, nwin; win_p(i) = avg_p(winlist(i)); end do
    call argsort_asc(win_p(1:nwin), nwin, idx)   ! ascending
    do j = 1, nw; seed(j) = winlist(idx(nwin-j+1)); end do   ! top nw
    s0 = minval(seed); s1 = maxval(seed)          ! 1-based band indices
    npool = max(n_pool_max, s1-s0+1)
    s0 = max(1, min(s0, nb-npool+1))
    npool = min(npool, nb-s0+1)
    allocate(pool(npool))
    do i = 1, npool; pool(i) = s0 + i - 1; end do
  end subroutine select_pool

  !> W-6: the pool-size FLOOR select_pool builds to, as ONE expression so the
  !> post-exclusion ratio gate (w2w_hybrid.F) cannot drift from it.
  pure integer function pool_floor(nw, pool_factor) result(n)
    integer,  intent(in) :: nw
    real(dp), intent(in) :: pool_factor
    n = max(nw+1, ceiling(pool_factor*real(nw, dp)))
  end function pool_floor

  !> W-6: the pre-repair frozen DEMAND over a pool -- build_masks'
  !> froz_demand_out, shared with the post-exclusion capacity gate so the
  !> two passes count the same thing.
  pure integer function pool_demand(eig, nk, pool, npool, e_fermi, fe_used, &
                                    has_ceiling, frozen) result(d)
    integer,  intent(in) :: nk, npool, pool(npool)
    real(dp), intent(in) :: eig(:,:), e_fermi, fe_used
    logical,  intent(in) :: has_ceiling, frozen(npool,nk)
    integer :: kx, ix, dk
    d = 0
    do kx = 1, nk
       dk = 0
       do ix = 1, npool
          if (has_ceiling) then
             if (eig(pool(ix),kx) - e_fermi <= fe_used) dk = dk + 1
          else
             if (frozen(ix,kx)) dk = dk + 1
          end if
       end do
       d = max(d, dk)
    end do
  end function pool_demand

  !--------------------------------------------------------- auto-trust thresholds
  subroutine auto_trust(P, eig, nb, nk, nw, e_fermi, pool_factor, p_froz, p_froz_band, ok, &
                        avg_out, min_out, bot_out, plateau_out, fid_out)
    integer, intent(in) :: nb, nk, nw
    real(dp), intent(in) :: P(nb,nk), eig(nb,nk), e_fermi, pool_factor
    real(dp), intent(out) :: p_froz, p_froz_band
    logical, intent(out) :: ok
    ! anatomy of the rule (what <out>.anatomy records): per-band average
    ! and minimum projectability, band bottom relative to E_F, the gray
    ! plateau the gate was read from, and the fidelity ceiling used
    real(dp), intent(out), optional :: avg_out(nb), min_out(nb), bot_out(nb)
    real(dp), intent(out), optional :: plateau_out, fid_out
    real(dp), parameter :: fid_emax = 8.0_dp, margin = 0.015_dp, eps = 0.02_dp, thr_floor = 0.55_dp
    integer, allocatable :: pool(:)
    integer :: npool, i, ip
    real(dp) :: avg_p(nb), min_p(nb), bot_e(nb), plateau_top, gate, thr, mn
    logical :: any_irr, any_rel, any_tru
    call select_pool(P, eig, nb, nk, nw, e_fermi, pool_factor, -30.0_dp, 20.0_dp, pool, npool)
    do i = 1, nb
       avg_p(i) = sum(P(i,:))/real(nk, dp)
       min_p(i) = minval(P(i,:))
       bot_e(i) = minval(eig(i,:)) - e_fermi
    end do
    if (present(avg_out)) avg_out = avg_p
    if (present(min_out)) min_out = min_p
    if (present(bot_out)) bot_out = bot_e
    if (present(fid_out)) fid_out = fid_emax
    if (present(plateau_out)) plateau_out = -huge(1.0_dp)
    plateau_top = -huge(1.0_dp); any_irr = .false.; any_rel = .false.
    do ip = 1, npool
       i = pool(ip)
       if (bot_e(i) > fid_emax) then
          any_irr = .true.; plateau_top = max(plateau_top, avg_p(i))
       else
          any_rel = .true.
       end if
    end do
    if (present(plateau_out)) plateau_out = plateau_top
    if (.not. any_irr .or. .not. any_rel) then
       ok = .false.; deallocate(pool); return
    end if
    gate = plateau_top + margin
    any_tru = .false.; mn = huge(1.0_dp)
    do ip = 1, npool
       i = pool(ip)
       if (bot_e(i) <= fid_emax .and. avg_p(i) >= gate) then
          any_tru = .true.; mn = min(mn, min_p(i))
       end if
    end do
    if (.not. any_tru) then
       ok = .false.; deallocate(pool); return
    end if
    thr = max(thr_floor, min(0.95_dp, mn - eps))
    p_froz = thr; p_froz_band = gate; ok = .true.
    deallocate(pool)
  end subroutine auto_trust


  subroutine select_pool_and_frozen(P, eig, nb, nk, nw, p_froz, p_froz_band, p_floor, &
       pool_factor, e_fermi, use_fermi_shell, fermi_shell1, use_shell_guard, shell_guard, &
       pool, npool, admit, frozen)
    integer, intent(in) :: nb, nk, nw
    real(dp), intent(in) :: P(nb,nk), eig(nb,nk), p_froz, p_froz_band, p_floor, &
                            pool_factor, e_fermi, fermi_shell1, shell_guard
    logical, intent(in) :: use_fermi_shell, use_shell_guard
    integer, allocatable, intent(out) :: pool(:)
    integer, intent(out) :: npool
    logical, allocatable, intent(out) :: admit(:,:), frozen(:,:)
    real(dp), allocatable :: pp(:,:), ee(:,:), avgp(:), score(:)
    logical, allocatable :: shell(:,:), fermi_band(:), crossing(:)
    integer :: i, k, j, nf, short, idx0
    real(dp) :: vbm, cbm, ceiling, midgap, mn, mx
    integer, allocatable :: ord(:)
    logical :: insulator

    call select_pool(P, eig, nb, nk, nw, e_fermi, pool_factor, -30.0_dp, 20.0_dp, pool, npool)
    allocate(pp(npool,nk), ee(npool,nk), avgp(npool), admit(npool,nk), frozen(npool,nk), &
             shell(npool,nk), fermi_band(npool), crossing(npool))
    do k = 1, nk
       do i = 1, npool
          pp(i,k) = P(pool(i),k)
          ee(i,k) = eig(pool(i),k) - e_fermi
       end do
    end do
    do i = 1, npool; avgp(i) = sum(pp(i,:))/real(nk, dp); end do
    do k = 1, nk
       do i = 1, npool
          admit(i,k)  = pp(i,k) >= p_floor
          frozen(i,k) = (pp(i,k) >= p_froz) .and. (avgp(i) >= p_froz_band)
       end do
    end do

    shell = .false.
    if (use_fermi_shell) then
       vbm = -huge(1.0_dp); cbm = huge(1.0_dp)
       do k = 1, nk; do i = 1, npool
          if (ee(i,k) <= 0.0_dp) vbm = max(vbm, ee(i,k))
          if (ee(i,k) >  0.0_dp) cbm = min(cbm, ee(i,k))
       end do; end do
       if (vbm == -huge(1.0_dp)) vbm = 0.0_dp
       if (cbm ==  huge(1.0_dp)) cbm = 0.0_dp
       do i = 1, npool
          mn = minval(ee(i,:)); mx = maxval(ee(i,:))
          crossing(i) = (mn < 0.0_dp) .and. (mx > 0.0_dp)
       end do
       insulator = any(ee <= 0.0_dp) .and. any(ee > 0.0_dp) .and. &
                   (.not. any(crossing)) .and. (cbm - vbm > 0.0_dp)
       if (insulator) then
          midgap = 0.5_dp*(vbm+cbm)
          do k = 1, nk; do i = 1, npool; shell(i,k) = ee(i,k) <= midgap; end do; end do
          fermi_band = .false.
       else
          do i = 1, npool; fermi_band(i) = minval(ee(i,:)) <= abs(fermi_shell1); end do
          ceiling = -huge(1.0_dp)
          do k = 1, nk; do i = 1, npool
             if (fermi_band(i)) ceiling = max(ceiling, ee(i,k))
          end do; end do
          do k = 1, nk; do i = 1, npool
             shell(i,k) = (ee(i,k) <= ceiling + 1.0e-6_dp) .or. fermi_band(i)
          end do; end do
       end if
       if (use_shell_guard) then
          do k = 1, nk; do i = 1, npool
             shell(i,k) = shell(i,k) .and. (pp(i,k) >= shell_guard)
          end do; end do
       end if
       do k = 1, nk; do i = 1, npool
          if (shell(i,k)) then; frozen(i,k) = .true.; admit(i,k) = .true.; end if
       end do; end do
    end if

    ! per-k feasibility repair
    allocate(score(npool), ord(npool))
    do k = 1, nk
       nf = count(frozen(:,k))
       if (nf > nw) then
          do i = 1, npool
             score(i) = pp(i,k) + 10.0_dp*merge(1.0_dp, 0.0_dp, shell(i,k))
             if (.not. frozen(i,k)) score(i) = -huge(1.0_dp)   ! only among frozen
          end do
          call argsort_asc(score, npool, ord)
          frozen(:,k) = .false.
          do j = 1, nw; frozen(ord(npool-j+1),k) = .true.; end do
       end if
       short = nw - count(admit(:,k))
       if (short > 0) then
          do i = 1, npool
             score(i) = pp(i,k)
             if (admit(i,k)) score(i) = -huge(1.0_dp)          ! only among non-admitted
          end do
          call argsort_asc(score, npool, ord)
          do j = 1, short; admit(ord(npool-j+1),k) = .true.; end do
       end if
    end do
    frozen = frozen .and. admit
    deallocate(pp, ee, avgp, shell, fermi_band, crossing, score, ord)
  end subroutine select_pool_and_frozen

  !--------------------------------------------------- the frozen DEMAND at fe
  !> demand(fe) = max_k |{ bands b : eig(b,k) - E_F <= fe }|, the quantity the
  !> capacity gate binds (nw >= demand).
  !>
  !> Counted over ALL nb bands, and that is not sloppiness -- it is the only
  !> count that agrees with the one build_masks reports.  build_masks' pool
  !> guard widens the pool to [blo, bhi] with blo the LOWEST and bhi the
  !> HIGHEST band index having any sub-ceiling k, so every band that can
  !> contribute is inside the post-guard pool and the all-band count IS
  !> froz_demand_out.  Counting over the PRE-guard pool (which is what is in
  !> scope here) would under-report, and the walk below would then select a
  !> ceiling the capacity gate goes on to reject -- exactly the failure mode
  !> this routine exists to close.  MEASURED on a 76-band Rh3Si4 pool
  !> (pool[1..76] = all bands): this function reproduces the printed demand
  !> at every rung of the ladder (8.0 -> 76, 3.0 -> 69, 2.0 -> 64, 1.5 -> 60,
  !> 1.0 -> 56, 0.5 -> 53).
  pure function froz_demand_at(eig, nk, e_fermi, fe) result(d)
    real(dp), intent(in) :: eig(:,:), e_fermi, fe
    integer, intent(in) :: nk
    integer :: d, k, b, c
    d = 0
    do k = 1, nk
       c = 0
       do b = 1, size(eig,1)
          if (eig(b,k) - e_fermi <= fe) c = c + 1
       end do
       if (c > d) d = c
    end do
  end function froz_demand_at

  !> The supremum of the ceilings that fit: demand(fe) <= nw for EVERY
  !> fe < fe_sup, and demand(fe) > nw for fe >= fe_sup.  Exact, no ladder:
  !> demand(fe) <= nw at k iff fewer than nw+1 states sit at or below fe, i.e.
  !> iff fe < E_(nw+1)(k); over all k that is the minimum of the (nw+1)-th
  !> smallest eigenvalue.  Returns huge() when nb <= nw (no ceiling can ever
  !> overflow the budget).
  function froz_fe_sup(eig, nk, nw, e_fermi) result(fes)
    real(dp), intent(in) :: eig(:,:), e_fermi
    integer, intent(in) :: nk, nw
    real(dp) :: fes, col(size(eig,1))
    integer :: k, nb, idx(size(eig,1))
    nb = size(eig,1)
    if (nb <= nw) then; fes = huge(1.0_dp); return; end if
    fes = huge(1.0_dp)
    do k = 1, nk
       col = eig(:,k)
       call argsort_asc(col, nb, idx)
       fes = min(fes, col(idx(nw+1)) - e_fermi)
    end do
  end function froz_fe_sup

  !> froz_emax = auto could not find a ceiling that fits: REFUSE, with ranked
  !> alternatives.  Never returns.
  !>
  !> WHY THIS IS A REFUSAL AND NOT A WARNING.  When demand > nw at every rung
  !> of the walk the mask builder's per-k capacity repair pins rows 1..nw
  !> frozen at every k, so frozen/k == nw == the whole budget, there is NOT
  !> ONE free state, the gradient is identically zero and the solver "converges"
  !> after a single evaluation without moving.  MEASURED on the 76-band
  !> Rh3Si4 pool with nw = 60, before this refusal existed: auto stopped at
  !> fe = 2.0 (demand 64), frozen/k mean = 60.00, "cg: LINE SEARCH STALLED
  !> after 1 gradient evaluations", Omega_I 456.85 against 397.68 at
  !> fe = 1.0.  A degenerate problem that reports success is worse than a
  !> failure, so this refuses.
  subroutine refuse_froz_auto(eig, nk, nw, e_fermi, fe_lo, dmd_lo)
    real(dp), intent(in) :: eig(:,:), e_fermi, fe_lo
    integer, intent(in) :: nk, nw, dmd_lo
    real(dp) :: fes, fe_rec
    integer :: dmd_rec
    fes = froz_fe_sup(eig, nk, nw, e_fermi)
    ! highest 0.05 eV grid point strictly below the supremum -- a value a user
    ! can type back in as arg 6 and reproduce
    fe_rec = 0.05_dp * real(ceiling(fes/0.05_dp) - 1, dp)
    dmd_rec = froz_demand_at(eig, nk, e_fermi, fe_rec)
    write(*,'(a,f6.2,a,f6.2,a,f4.2,a)') &
         '  *** froz_emax = auto: NO ceiling in the walk (', FE_WALK_HI, &
         ' -> ', FE_WALK_LO, ' eV, step ', FE_WALK_STEP, ') fits num_wann.'
    write(*,'(a,f7.2,a,i0,a,i0)') &
         '  *** at the walk floor fe = ', fe_lo, ' eV the demand is ', dmd_lo, &
         ' > num_wann ', nw
    write(*,'(a)') &
         '  *** Continuing would freeze num_wann states at EVERY k: no free'
    write(*,'(a)') &
         '  *** states, an identically zero gradient, and a run that reports'
    write(*,'(a)') &
         '  *** success without moving.  REFUSING.  Ranked alternatives:'
    if (fes < huge(1.0_dp)) then
       write(*,'(a,f7.2,a,i0,a,i0,a)') &
            '  ***   (1) froz_emax = ', fe_rec, '   (demand ', dmd_rec, &
            ' <= num_wann ', nw, ') -- the highest 0.05 eV step that fits'
       write(*,'(a,f9.4,a)') &
            '  ***       (exact supremum: any froz_emax < ', fes, &
            ' eV fits; = min_k E_(nw+1)(k) - E_F)'
    end if
    ! (2) is NOT always executable and says so.  MEASURED: GaAs nw=7 -> "raise
    ! to 9" dies at "num_wann = 9 exceeds the 7 projections in the amn", and
    ! SrVO3 nw=9 -> 12 dies at 9 projections.  A ranked alternative that the
    ! next run rejects for a different reason is worse than no alternative, so
    ! the projection precondition is stated with it.
    write(*,'(a,i0,a,f7.2,a)') &
         '  ***   (2) num_wann >= ', dmd_lo, ' -- keep the fe = ', fe_lo, &
         ' eV window instead of shrinking it'
    write(*,'(a,i0,a)') &
         '  ***       (needs at least ', dmd_lo, &
         ' amn projection columns; the run refuses otherwise)'
    write(*,'(a)') &
         '  ***   (3) froz_emax = none -- no window demand at all (the frozen'
    write(*,'(a)') &
         '  ***       constraint then carries no band-fidelity guarantee)'
    error stop 'froz_emax = auto: no ceiling satisfies demand <= num_wann'
  end subroutine refuse_froz_auto

  !------------------------------------------------------------- froz_emax ceiling
  !
  ! DEFECT D2 (fixed here).  The `auto` walk-down tested only
  !         (.not. is_auto) .or. sat_frac <= 0.15 .or. fe <= 1.0
  ! and NEVER tested the capacity constraint demand <= nw, so it could -- and
  ! on rd76 nw=60 did -- halt at a ceiling whose demand overflows the budget
  ! and hand back the degenerate zero-gradient problem described on
  ! refuse_froz_auto.  The ladder is UNCHANGED (8.0 eV down to 1.0 eV in 1.0 eV
  ! steps, HEAD's own); what is added is the demand test in the exit and the
  ! fail-closed refusal when the whole walk is infeasible.
  !
  ! EXPLICIT ceilings are untouched.  mode = '<eV>' is the user's own choice
  ! and exits on the first pass exactly as before, so every archived run
  ! replayed with its recorded ceiling is bit-identical (regression anchor:
  ! rd76 froz_emax = 1.0 -> Omega_I 397.67782615).  A too-high EXPLICIT ceiling
  ! is still reported by the capacity gate and by froz_report_line; it is not
  ! silently overridden, because overriding an explicit argument is its own
  ! kind of dishonesty.
  subroutine apply_froz_emax(frozen, eig, pool, npool, nk, nw, e_fermi, mode, emax_val, fe_out)
    integer, intent(in) :: npool, nk, nw, pool(npool)
    logical, intent(inout) :: frozen(npool,nk)
    real(dp), intent(in) :: eig(:,:), e_fermi, emax_val
    character(*), intent(in) :: mode
    real(dp), intent(out), optional :: fe_out   ! ceiling actually applied
    real(dp) :: fe, ez, etop
    logical :: keep
    integer :: k, i, dmd
    logical :: is_auto
    if (mode == 'none') then
       if (present(fe_out)) fe_out = emax_val
       return
    end if
    is_auto = (mode == 'auto')
    fe = merge(FE_WALK_HI, emax_val, is_auto)
    if (is_auto .and. FE_WALK_FROM_TOP) then
       etop = -huge(1.0_dp)
       do k = 1, nk; do i = 1, npool
          if (frozen(i,k)) etop = max(etop, eig(pool(i),k) - e_fermi)
       end do; end do
       fe = max(FE_WALK_HI, etop)
    end if
    ! Exit on the capacity test alone.  The v3.0 walk also required at most
    ! 15% of k-points to be saturated (trusted frozen states below fe >= nw);
    ! on the eight release cells that test never bound (every sub-8 eV ceiling
    ! was set by demand), and it misjudges an isolated target manifold, which
    ! legitimately freezes nw states at every k (LCAO Bi bilayer: 16 of 16).
    ! demand <= nw still excludes the degenerate truncated problem that
    ! refuse_froz_auto describes.
    do
       if (.not. is_auto) exit
       dmd = froz_demand_at(eig, nk, e_fermi, fe)
       if (dmd <= nw) exit
       if (fe <= FE_WALK_LO) call refuse_froz_auto(eig, nk, nw, e_fermi, fe, dmd)
       if (fe > FE_WALK_HI + 1.0e-9_dp .and. fe - FE_WALK_STEP < FE_WALK_HI) then
          fe = FE_WALK_HI
       else
          fe = fe - FE_WALK_STEP
       end if
    end do
    ! NB: the selected ceiling is NOT printed here.  It is returned in fe_out
    ! and froz_report_line states it ("froz_emax ceiling applied = ... mode =
    ! auto"), which is the one place both the standalone tool and the in-w2w
    ! hook print it -- and this routine is called from an MPI path where only
    ! root may write.
    do k = 1, nk; do i = 1, npool
       ez = eig(pool(i),k) - e_fermi
       keep = (ez <= fe) .or. (ez <= 0.0_dp)
       if (.not. keep) frozen(i,k) = .false.
    end do; end do
    if (present(fe_out)) fe_out = fe
  end subroutine apply_froz_emax

  !------------------------------------------------------------- band-level freeze
  subroutine apply_band_level(frozen, admit, P, eig, pool, npool, nk, nw, e_fermi, &
                              froz_emax_val, p_froz_band, user_set_pfb)
    integer, intent(in) :: npool, nk, nw, pool(npool)
    logical, intent(inout) :: frozen(npool,nk), admit(npool,nk)
    real(dp), intent(in) :: P(:,:), eig(:,:), e_fermi, froz_emax_val, p_froz_band
    logical, intent(in) :: user_set_pfb
    real(dp) :: avgP(npool), ezm(npool), ceil, gate, mx, cavg(npool)
    logical :: below(npool), bandf(npool)
    integer :: i, k, ncand, cand(npool), idx(npool), nkeep, j
    do i = 1, npool
       avgP(i) = sum(P(pool(i),:))/real(nk, dp)
       ezm(i)  = sum(eig(pool(i),:))/real(nk, dp) - e_fermi
    end do
    ceil = froz_emax_val
    below = ezm <= ceil
    if (user_set_pfb) then
       gate = p_froz_band
    else
       mx = -huge(1.0_dp)
       do i = 1, npool; if (below(i)) mx = max(mx, avgP(i)); end do
       gate = merge(0.5_dp*mx, p_froz_band, mx > -huge(1.0_dp))
    end if
    ncand = 0
    do i = 1, npool
       if (below(i) .and. avgP(i) >= gate) then; ncand = ncand+1; cand(ncand) = i; end if
    end do
    if (ncand > nw) then                       ! keep best nw by avgP
       do i = 1, ncand; cavg(i) = avgP(cand(i)); end do
       call argsort_asc(cavg(1:ncand), ncand, idx)
       nkeep = nw
    else
       do i = 1, ncand; idx(i) = i; end do
       nkeep = ncand
    end if
    bandf = .false.
    do j = 1, nkeep; bandf(cand(idx(ncand-j+1))) = .true.; end do
    do k = 1, nk; do i = 1, npool
       frozen(i,k) = bandf(i)
       if (bandf(i)) admit(i,k) = .true.
    end do; end do
  end subroutine apply_band_level

  !------------------------------------------------ S1a: the per-k frozen floor
  ! DEFECT (S1a).  The frozen DEMAND that the capacity gate is checked against
  ! is counted PER K:
  !         F(k) = { pool rows i : eig(pool(i),k) - E_F <= fe },
  !         demand = max_k |F(k)|                       (build_masks, below)
  ! but the mask the optimiser is actually handed CAN be a strict subset of
  ! that, in either regime:
  !   * band-level  frozen(i,k) = bandf(i), one whole-band decision taken on
  !                 the MEAN-over-k energy (apply_band_level).  A band whose
  !                 mean sits above fe but which dips below fe at some k is
  !                 frozen NOWHERE -- so the constraint that made it part of
  !                 the demand is never imposed anywhere.
  !   * per-state   frozen(i,k) survives select_pool_and_frozen's
  !                 projectability test AND apply_froz_emax's ceiling; the
  !                 ceiling only ever REMOVES, so the mask can only shrink
  !                 below F(k).
  ! Measured over 206 archived runs: 199 print a demand
  ! strictly greater than the frozen/k MEAN they enforce (191 band-level, 8
  ! per-state); only 7 agree.  mg8d.hyblog's "demand 5" above "frozen/k mean
  ! 4.00" is the rule, not an outlier, and that one missing state is MgB2's
  ! Mg-derived band -- the state nw=8 silently discards.
  ! That count comparison is an UPPER BOUND on the affected runs, and it is
  ! NOT the criterion.  A first reading of it said "where demand > nw there is
  ! no budget to add anything, so the run is inert"; that repeats the very
  ! max-over-k vs per-k confusion S1a is about.  demand = max_k |F(k)|, so
  ! demand > nw constrains only the WORST k -- at other k |F(k)| < nw and the
  ! floor still adds there, and even where |F(k)| >= nw the enforced set changes
  ! whenever band-level pinned a NON-PREFIX set of nw states.
  !
  ! MEASURED COUNTEREXAMPLE, against the earlier mask path (bit-identical
  ! control, so the change is this fix and nothing else).  GaAs, nw=4, band-level,
  ! admitall, gs, pool 1.5, np_mask 9:
  !     fe = -8.0 eV : demand 6 > nw 4, "+512 / -512 states at 512 of 512 k",
  !                    Omega_I 17.85138106 -> 10.40694414
  !     fe = +4.0 eV : demand 9 > nw 4, same, same Omega_I
  ! At np_mask = 4 the SAME cell is bit-identical at all three ceilings, which
  ! is why the first sweep (over fe only, never over nw or np_mask) missed it.
  !
  ! THE CRITERION IS SET CONTAINMENT: the fix changes a run exactly when the
  ! enforced mask differs from the energy prefix truncated to nw at ANY k.
  ! That is what enforce_frozen_floor and frozen_report test, and nothing is
  ! silent about it -- every such run prints the floor bookkeeping line and,
  ! where containment still fails, the *** S1a block.
  !
  ! WHERE demand > nw NO MASK CONTAINS F(k) and the choice of WHICH nw states
  ! to freeze is a POLICY, not a correctness question.  This code takes the nw
  ! DEEPEST sub-ceiling states, which is HEAD's own stated policy: the pre-
  ! band-level per-k capacity repair already did exactly that and
  ! apply_band_level unconditionally undid it.  Measured on the GaAs arm above
  ! with an independent probe of the band fidelity inside the window, the
  ! substitution is marginally WORSE on band fidelity inside the window
  ! (eta_mean 13014.68 -> 13026.87 meV at fe=4.0; 3713.33 -> 3725.51 at
  ! fe=-8.0) while Omega_I falls 42%.  Both arms are already rejected by the
  ! capacity gate (eta ~ 13 eV: nw=4 cannot carry a 9-state window at all), so
  ! neither policy is meaningful there -- but the substitution is REPORTED
  ! rather than assumed harmless.
  !
  ! A second, independent consequence: the per-k capacity repair in
  ! build_masks (frozen rows 1..nw at an over-demanded k) is UNCONDITIONALLY
  ! UNDONE by apply_band_level, which overwrites every (i,k).  So the log's
  ! "frozen budget spent deepest-first" promise was false in band-level runs.
  !
  ! FIX.  The enforced mask must CONTAIN F(k) at every k, up to the nw budget.
  ! pool() is a contiguous ascending band range (select_pool), so at every k
  ! the pool rows are energy-ordered and F(k) is exactly the prefix 1..nf(k):
  ! honouring it keeps the frozen manifold contiguous in energy, which is the
  ! property the band-level regime was introduced to protect (the incoherence
  ! that cost Bi2D eta ~3 eV came from capping a patchwork of P-BASED flags,
  ! not from an energy prefix).  Band-level's ADDITIONS -- rows it freezes at
  ! k where they sit above fe -- are kept as long as they fit; they are the
  ! whole-band coherence the regime exists for.
  !
  ! WHY NOT THE OTHER REPAIR (redefine the demand as the band-level count):
  ! apply_band_level caps nkeep at nw by construction, so a demand read off
  ! its output satisfies demand <= nw ALWAYS.  The capacity gate would become
  ! a quantity its own subject guarantees, i.e. a test that can never fail --
  ! the self-fulfilling-acceptance-test error.  Rejected on that ground, not
  ! on cost.
  !
  ! nadd/ndrop are counted so the change is attributable, never silent.
  subroutine enforce_frozen_floor(frozen, admit, eig, pool, npool, nk, nw, &
                                  e_fermi, fe, nadd, ndrop, nk_touched)
    integer, intent(in) :: npool, nk, nw, pool(npool)
    logical, intent(inout) :: frozen(npool,nk), admit(npool,nk)
    real(dp), intent(in) :: eig(:,:), e_fermi, fe
    integer, intent(out) :: nadd, ndrop, nk_touched
    integer :: k, i, nf, nextra, nseen, a0, d0
    logical :: want
    nadd = 0; ndrop = 0; nk_touched = 0
    do k = 1, nk
       a0 = nadd; d0 = ndrop
       nf = 0
       do i = 1, npool
          if (eig(pool(i),k) - e_fermi <= fe) nf = nf + 1
       end do
       if (nf >= nw) then
          ! F(k) fills (or overruns) the budget: the nw deepest sub-ceiling
          ! states, exactly the policy the pre-band-level repair intended.
          do i = 1, npool
             want = (i <= nw)
             if (want .and. .not. frozen(i,k)) nadd = nadd + 1
             if (frozen(i,k) .and. .not. want) ndrop = ndrop + 1
             frozen(i,k) = want
             if (want) admit(i,k) = .true.
          end do
       else
          ! honour F(k) first, then keep as many band-level extras as fit,
          ! deepest first (deterministic: pool rows are energy-ordered)
          do i = 1, nf
             if (.not. frozen(i,k)) then
                frozen(i,k) = .true.; nadd = nadd + 1
             end if
             admit(i,k) = .true.
          end do
          nextra = nw - nf
          nseen = 0
          do i = nf+1, npool
             if (frozen(i,k)) then
                nseen = nseen + 1
                if (nseen > nextra) then
                   frozen(i,k) = .false.; ndrop = ndrop + 1
                else
                   admit(i,k) = .true.
                end if
             end if
          end do
       end if
       if (nadd > a0 .or. ndrop > d0) nk_touched = nk_touched + 1
    end do
  end subroutine enforce_frozen_floor

  !------------------------------------------------------- frozen-mask report
  ! One shared implementation of the numbers a caller must print so that the
  ! log states the constraint the optimiser is ACTUALLY held to.  Computed
  ! from the FINAL mask, i.e. after any admit-all / no-froz override, so the
  ! printed value can never describe a mask that was later replaced.
  !   enf_worst   max_k count(frozen(:,k))       -- the binding constraint
  !   enf_mean    mean_k count(frozen(:,k))      -- historical "frozen/k mean"
  !   short_worst max_k |F(k) \ frozen(:,k)|     -- window states NOT imposed
  !   nk_short    #{k : F(k) is not contained}
  ! With has_ceiling = .false. (froz_emax = none) F(k) is undefined and the
  ! shortfall outputs are returned as 0 with nk_short = -1 (= not applicable).
  subroutine frozen_report(frozen, eig, pool, npool, nk, e_fermi, fe, &
                           has_ceiling, enf_worst, enf_mean, short_worst, nk_short)
    integer, intent(in) :: npool, nk, pool(npool)
    logical, intent(in) :: frozen(npool,nk), has_ceiling
    real(dp), intent(in) :: eig(:,:), e_fermi, fe
    integer, intent(out) :: enf_worst, short_worst, nk_short
    real(dp), intent(out) :: enf_mean
    integer :: k, i, c, s, tot
    enf_worst = 0; short_worst = 0; tot = 0
    nk_short = merge(0, -1, has_ceiling)
    do k = 1, nk
       c = 0; s = 0
       do i = 1, npool
          if (frozen(i,k)) c = c + 1
          if (has_ceiling) then
             if (eig(pool(i),k) - e_fermi <= fe .and. .not. frozen(i,k)) s = s + 1
          end if
       end do
       tot = tot + c
       enf_worst = max(enf_worst, c)
       if (has_ceiling) then
          short_worst = max(short_worst, s)
          if (s > 0) nk_short = nk_short + 1
       end if
    end do
    enf_mean = real(tot, dp)/real(max(nk,1), dp)
  end subroutine frozen_report

  !------------------------------------------- the authoritative frozen report
  ! ONE implementation of the lines that state what the optimiser is actually
  ! held to, shared by the standalone tool and the in-w2w hook so the two logs
  ! cannot drift.  The historical pair
  !     "frozen demand(worst k) = D"      (a per-k energy count)
  !     "... frozen/k mean= M"            (a mean over k of the enforced mask)
  ! could not be reconciled by a reader: D and M are different quantities of
  ! different sets, and nothing said which one the capacity gate binds.  This
  ! prints the ENFORCED worst-k count -- the binding constraint -- names the
  ! ceiling that defines the demand, and states the containment relation
  ! explicitly, with ranked alternatives keyed to the actual cause when it
  ! fails.  `iu` is the output unit; callers gate on is_root themselves.
  subroutine froz_report_line(iu, frozen, eig, pool, npool, nk, nw, e_fermi, &
       fe_used, femode, demand, no_froz, froz_legacy, fl_add, fl_drop, fl_nk)
    integer, intent(in) :: iu, npool, nk, nw, pool(npool), demand
    integer, intent(in) :: fl_add, fl_drop, fl_nk
    logical, intent(in) :: frozen(npool,nk), no_froz, froz_legacy
    real(dp), intent(in) :: eig(:,:), e_fermi, fe_used
    character(*), intent(in) :: femode
    integer :: enf_worst, short_worst, nk_short
    real(dp) :: enf_mean
    logical :: has_ceiling
    has_ceiling = (trim(femode) /= 'none') .and. (.not. no_froz)
    call frozen_report(frozen, eig, pool, npool, nk, e_fermi, fe_used, &
                       has_ceiling, enf_worst, enf_mean, short_worst, nk_short)
    write(iu,'(a,i0,a,f7.2,a,i0,a)') &
         '  frozen ENFORCED (this IS the optimiser constraint): worst-k = ', &
         enf_worst, '  mean/k = ', enf_mean, '   (nw = ', nw, ')'
    if (.not. has_ceiling) then
       write(iu,'(a)') '    froz_emax off -- no window demand is defined'
       return
    end if
    write(iu,'(a,f9.4,a,a,a)') '    froz_emax ceiling applied = ', fe_used, &
         ' eV  (mode = ', trim(femode), ');  demand = max_k |{E-E_F <= ceiling}|'
    if (nk_short == 0) then
       write(iu,'(a,i0,a,i0,a)') &
            '    the ceiling set is CONTAINED at all ', nk, &
            ' k -- demand(worst k) ', demand, ' is HELD'
    else
       write(iu,'(a,i0,a,i0,a,i0,a)') &
            '  *** S1a: the ceiling set is NOT contained at ', nk_short, ' of ', &
            nk, ' k (worst k short by ', short_worst, ' state(s))'
       write(iu,'(a)') &
            '  *** the capacity gate above is then checked against a demand the'
       write(iu,'(a)') &
            '  *** optimiser is never made to satisfy.  Ranked alternatives:'
       if (demand > nw) then
          write(iu,'(a,i0,a)') &
               '  ***   (1) raise num_wann to >= ', demand, &
               ' -- the window genuinely needs that many states'
          write(iu,'(a)') &
               '  ***   (2) lower froz_emax until the demand fits num_wann'
          write(iu,'(a)') &
               '  ***   (3) proceed knowing band fidelity near the ceiling is not constrained'
       else
          write(iu,'(a)') &
               '  ***   (1) regime = per-state, so no whole-band mean decides the mask'
          write(iu,'(a)') &
               '  ***   (2) lower froz_emax so no band straddles the ceiling'
       end if
    end if
    if (fl_add > 0 .or. fl_drop > 0) then
       write(iu,'(a,i0,a,i0,a,i0,a,i0,a)') &
            '    per-k frozen floor: +', fl_add, ' / -', fl_drop, &
            ' states at ', fl_nk, ' of ', nk, ' k'
       if (fl_drop > 0 .and. demand > nw) then
          ! No mask can contain F(k) here, so the floor has SUBSTITUTED a
          ! policy (the nw deepest sub-ceiling states) for band-level's
          ! projectability-ranked choice.  That is HEAD's own stated repair
          ! policy, which apply_band_level was undoing -- but it is a policy,
          ! it changes the answer, and it buys no band fidelity.  Say so.
          write(iu,'(a,i0,a,i0,a)') &
               '    NOTE: demand ', demand, ' > num_wann ', nw, &
               ' -- NO mask can contain the ceiling set, so the floor'
          write(iu,'(a)') &
               '    SUBSTITUTED the nw deepest sub-ceiling states for the '// &
               'band-level choice.'
          write(iu,'(a)') &
               '    That is a POLICY (the pre-band-level capacity repair''s '// &
               'own), not a fidelity gain:'
          write(iu,'(a)') &
               '    measured marginally WORSE on window fidelity while '// &
               'Omega_I falls.  Fix the capacity'
          write(iu,'(a,i0,a)') &
               '    first (raise num_wann to >= ', demand, &
               ' or lower froz_emax); this arm is not admissible.'
       end if
    end if
  end subroutine froz_report_line

  !--------------------------------------------------------- full mask orchestrator
  ! Mirrors the bridge main() mask section: regime -> rescale/symmetrize ->
  ! auto-trust -> select_pool_and_frozen -> froz_emax -> band-level ->
  ! per-k frozen floor (S1a).
  subroutine build_masks(P_in, eig, nb, nk, nw, e_fermi, rescale_mode, regime_mode, &
       froz_emax_mode, froz_emax_val, pool, npool, admit, frozen, regime_out, &
       pool_factor_in, froz_demand_out, fe_used_out, froz_legacy_out, &
       floor_add_out, floor_drop_out, floor_nk_out, p_floor_in, &
       P_used_out, anat_avg, anat_min, anat_bot, anat_plateau, anat_fid, &
       anat_pgate, anat_pfroz, anat_ok)
    integer, intent(in) :: nb, nk, nw
    ! the projectability the masks were built from (after rescale/symmetrize)
    ! and the auto-trust anatomy, for <out>.masks / <out>.anatomy
    real(dp), intent(out), optional :: P_used_out(nb,nk)
    real(dp), intent(out), optional :: anat_avg(nb), anat_min(nb), anat_bot(nb)
    real(dp), intent(out), optional :: anat_plateau, anat_fid, anat_pgate, anat_pfroz
    logical, intent(out), optional :: anat_ok
    real(dp), intent(in) :: P_in(nb,nk), eig(nb,nk), e_fermi, froz_emax_val
    character(*), intent(in) :: rescale_mode, regime_mode, froz_emax_mode
    integer, allocatable, intent(out) :: pool(:)
    integer, intent(out) :: npool
    logical, allocatable, intent(out) :: admit(:,:), frozen(:,:)
    character(*), intent(out) :: regime_out
    real(dp), intent(in), optional :: pool_factor_in   ! default 1.5 (historical)
    ! projectability floor of the admission mask: a pool state outside the
    ! frozen set may enter the disentangled subspace only if P >= p_floor
    real(dp), intent(in), optional :: p_floor_in       ! default 0.10
    ! pre-repair frozen DEMAND = max_k count(pool states below the frozen
    ! ceiling).  The per-k feasibility repair in select_pool_and_frozen
    ! silently demotes frozen states when the count exceeds nw (the
    ! Rh3Si4 silent-capacity-failure mechanism), so the returned masks
    ! always satisfy count<=nw; THIS value carries the true requirement.
    integer, intent(out), optional :: froz_demand_out
    ! the ceiling ACTUALLY applied (froz_emax = auto resolves only here, and
    ! it is recorded nowhere else -- without it F(k) is not reconstructible
    ! from the artifacts and condition (iv) can only be bracketed)
    real(dp), intent(out), optional :: fe_used_out
    ! always .false. in this release (the pre-S1a legacy mask path is gone)
    logical, intent(out), optional :: froz_legacy_out
    ! per-k frozen floor bookkeeping (states added / dropped / k touched)
    integer, intent(out), optional :: floor_add_out, floor_drop_out, floor_nk_out
    real(dp), allocatable :: P(:,:)
    real(dp) :: pfl
    real(dp) :: dtol, degen_frac, p_froz, p_froz_band, fe_used, poolf
    logical :: band_level, use_sym, ok, froz_legacy
    integer :: fl_add, fl_drop, fl_nk
    dtol = 0.02_dp
    poolf = 1.5_dp
    if (present(pool_factor_in)) poolf = pool_factor_in
    pfl = 0.10_dp
    if (present(p_floor_in)) pfl = p_floor_in
    froz_legacy = .false.
    fl_add = 0; fl_drop = 0; fl_nk = 0
    allocate(P(nb,nk)); P = P_in
    call auto_regime(eig, nb, nk, dtol, 0.7_dp, band_level, use_sym, degen_frac)
    if (regime_mode == 'per-state') then; band_level = .false.; use_sym = .false.
    else if (regime_mode == 'band-level') then; band_level = .true.; use_sym = .true.
    end if
    regime_out = merge('band-level', 'per-state ', band_level)
    call rescale_symmetrize(P, eig, nb, nk, e_fermi, rescale_mode, use_sym, dtol)
    call auto_trust(P, eig, nb, nk, nw, e_fermi, poolf, p_froz, p_froz_band, ok, &
                    avg_out=anat_avg, min_out=anat_min, bot_out=anat_bot, &
                    plateau_out=anat_plateau, fid_out=anat_fid)
    if (.not. ok) then; p_froz = 0.95_dp; p_froz_band = 0.85_dp; end if
    if (present(anat_ok)) anat_ok = ok
    if (present(anat_pgate)) anat_pgate = p_froz_band
    if (present(anat_pfroz)) anat_pfroz = p_froz
    if (present(P_used_out)) P_used_out = P
    call select_pool_and_frozen(P, eig, nb, nk, nw, p_froz, p_froz_band, pfl, &
         poolf, e_fermi, .true., 3.0_dp, .true., 0.15_dp, pool, npool, admit, frozen)
    call apply_froz_emax(frozen, eig, pool, npool, nk, nw, e_fermi, froz_emax_mode, &
                         froz_emax_val, fe_used)
    ! ---- pool guard: every window state below the frozen ceiling must be
    ! IN the pool -- at BOTH ends.  A sub-froz_emax state outside the pool
    ! is structurally unrepresentable: below pool(1) = the Rh3Si4/MgB2
    ! deep-band exclusion; above pool(npool) = the pool-factor cap hiding
    ! ceiling states from the disentangler (the Bi2D one-shot nw=12
    ! failure: pool top 18 vs true window need 26 -> eta ~3 eV, while the
    ! guard-widened [1..34] benchmark cell gives eta 3.34 meV).
    if (trim(froz_emax_mode) /= 'none') then
       block
         integer :: blo, bhi, bx, kx, add_lo, add_hi, np2
         logical, allocatable :: adm2(:,:), fro2(:,:)
         integer, allocatable :: pool2(:)
         blo = pool(1)
         do bx = 1, pool(1)-1
            if (any(eig(bx,:) - e_fermi <= fe_used)) then
               blo = bx
               exit
            end if
         end do
         bhi = pool(npool)
         do bx = nb, pool(npool)+1, -1
            if (any(eig(bx,:) - e_fermi <= fe_used)) then
               bhi = bx
               exit
            end if
         end do
         add_lo = pool(1) - blo
         add_hi = bhi - pool(npool)
         if (add_lo > 0 .or. add_hi > 0) then
            np2 = npool + add_lo + add_hi
            allocate(pool2(np2), adm2(np2,nk), fro2(np2,nk))
            do bx = 1, add_lo; pool2(bx) = blo+bx-1; end do
            pool2(add_lo+1:add_lo+npool) = pool
            do bx = 1, add_hi; pool2(add_lo+npool+bx) = pool(npool)+bx; end do
            adm2(add_lo+1:add_lo+npool,:) = admit
            fro2(add_lo+1:add_lo+npool,:) = frozen
            do kx = 1, nk
               do bx = 1, add_lo
                  adm2(bx,kx) = .true.
                  fro2(bx,kx) = (eig(blo+bx-1,kx) - e_fermi <= fe_used)
               end do
               do bx = 1, add_hi
                  adm2(add_lo+npool+bx,kx) = .true.
                  fro2(add_lo+npool+bx,kx) = &
                       (eig(pool(npool)+bx,kx) - e_fermi <= fe_used)
               end do
            end do
            call move_alloc(pool2, pool)
            call move_alloc(adm2, admit)
            call move_alloc(fro2, frozen)
            npool = np2
         end if
         ! per-k capacity handling (runs whether or not the pool was
         ! extended).  When the sub-ceiling DEMAND exceeds nw at a k,
         ! the frozen budget is spent on the nw DEEPEST sub-ceiling
         ! states BY ENERGY -- a deterministic, contiguous manifold
         ! (the Bi2D nw=12 benchmark behaviour; a patchwork of P-based
         ! flags capped row-wise gave an incoherent subspace, eta ~3 eV
         ! vs 3.34 meV).  Within budget the P-based frozen set stands
         ! (apply_froz_emax already capped it at the ceiling, so its
         ! count can never exceed the demand).
         block
           integer :: kk, ii, nf
           do kk = 1, nk
              nf = 0
              do ii = 1, npool
                 if (eig(pool(ii),kk) - e_fermi <= fe_used) nf = nf + 1
              end do
              if (nf > nw) then
                 ! pool rows are energy-ordered: the nw deepest
                 ! sub-ceiling states are exactly rows 1..nw
                 frozen(:,kk) = .false.
                 frozen(1:nw,kk) = .true.
              end if
           end do
         end block
       end block
    end if
    if (band_level) &
       call apply_band_level(frozen, admit, P, eig, pool, npool, nk, nw, e_fermi, &
                             fe_used, p_froz_band, .false.)
    ! ---- S1a: make the enforced mask CONTAIN the per-k sub-ceiling set ----
    ! (must run LAST: apply_band_level overwrites every (i,k), including the
    ! per-k capacity repair above)
    if (trim(froz_emax_mode) /= 'none') &
       call enforce_frozen_floor(frozen, admit, eig, pool, npool, nk, nw, &
                                 e_fermi, fe_used, fl_add, fl_drop, fl_nk)
    if (present(fe_used_out))    fe_used_out    = fe_used
    if (present(froz_legacy_out)) froz_legacy_out = froz_legacy
    if (present(floor_add_out))  floor_add_out  = fl_add
    if (present(floor_drop_out)) floor_drop_out = fl_drop
    if (present(floor_nk_out))   floor_nk_out   = fl_nk
    ! pre-repair frozen demand (energy count below the ceiling; for
    ! mode='none' the ceiling concept is off -> post-mask count)
    if (present(froz_demand_out)) froz_demand_out = pool_demand(eig, nk, pool, npool, &
         e_fermi, fe_used, trim(froz_emax_mode) /= 'none', frozen)
    deallocate(P)
  end subroutine build_masks

end module hybrid_prep_m
