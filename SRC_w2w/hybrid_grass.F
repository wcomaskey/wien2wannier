!!! wien2wannier/SRC_w2w/hybrid_grass.F
!!!
!!!    Direct Omega_I minimization on the Grassmann manifold.
!!!
!!! The default solver (`solver = cg'); the Gauss-Seidel sweep of
!!! hybrid_core_m minimizes the same functional (`solver = gs').  Direct
!!! minimization of the entangled-band problem was formulated by Damle,
!!! Levitt and Lin [DamleLevittLin2019]; here it is applied to the
!!! frozen-constrained Omega_I of Souza, Marzari and Vanderbilt [Souza2001]
!!! in the Marzari-Vanderbilt discretization [Marzari1997].  Keys refer to
!!! REFERENCES.md in the root directory.  The same solver is in lcao2wannier
!!! (lcao2wannier/hybrid.py); the per-k timings quoted below were measured
!!! there.
!!!
!!! WHY.  The sweep is a fixed-point iteration: sequential in k, and it forms
!!! an eigendecomposition at every k.  Two measured consequences:
!!!
!!!   * its only k-parallel path is red-black colouring, and colouring reaches
!!!     a DIFFERENT fixed point (Rh3Si4 H120: 840.98 coloured vs 840.06
!!!     serial), so the parallel option changes the answer;
!!!   * eigh dominates at scale.  LCAO-measured per k, single-threaded:
!!!       (ndis, nadd)      eigh     Zd@V + QR
!!!       (261, 114)      12.0 ms      3.1 ms
!!!       (566, 278)      99.9 ms     21.5 ms
!!!       (900, 450)     414.0 ms     38.0 ms
!!!
!!! A gradient method needs only Zd@V and a retraction, and evaluates every k
!!! from the SAME point -- so it parallelizes exactly over k, with no colouring
!!! and no dependency graph, at any nk.  A line search gives monotone descent,
!!! so this is NOT the divergent Jacobi iteration: that divergence is a
!!! property of the fixed-point map, not of a line-searched descent.
!!!
!!! THE OBJECTIVE.  With V_k the FULL (npool,nw) frame whose first nfroz(k)
!!! columns are fixed unit vectors,
!!!
!!!   Omega_I = nw*sum_b w_b - F/nk ,
!!!   F(V) = sum_kb w_b || V_k^H M_kb V_k+b ||_F^2 ,
!!!
!!! so minimizing Omega_I means MAXIMIZING F.  Writing T = V_k^H M V_k+b in
!!! frozen/free COLUMN blocks, the Frobenius norm splits with no cross terms:
!!!
!!!   ||T||^2 = ||T_ff||^2 + ||T_fd||^2 + ||T_df||^2 + ||T_dd||^2
!!!
!!!   T_ff -> aconst                       (constant)
!!!   T_df -> Q_k   += w_b Mf Mf^H,   Mf = M_kb[dis(k), froz(k+b)]
!!!   T_fd -> Q_k+b += w_b X^H X,     X  = M_kb[froz(k), dis(k+b)]
!!!   T_dd -> the free-free term with Mdd = sqrt(w_b) M_kb[dis(k), dis(k+b)]
!!!
!!! giving F = aconst + sum_k Tr(V'^H Q_k V') + sum_kb ||V'^H Mdd V'||^2 with
!!! V' the FREE block only.  This is the same frozen split build_kinfo already
!!! performs -- the index sets come from it, they are not rebuilt here.
!!!
!!! THE GRADIENT, which is where the prototype was first wrong: V'_k appears
!!! in F TWICE, as the CENTRE of ||V_k^H M V_k+b||^2 and as the NEIGHBOUR of
!!! ||V_k'^H M V_k||^2.  The first version had only the centre term.  It ran,
!!! it descended plausibly, and it was caught ONLY by finite differences
!!! (FD/analytic 1.84 and 2.10 -- not a constant, because Q contributes once
!!! while the M term contributes twice).  Hence gcg_check_gradient is a
!!! STARTUP GUARD here, not an optional test.
!!!
!!! PARALLELISM: gather, never scatter.  The neighbour term naturally scatters
!!! into Zd(neigh(k,b)), which races.  rev_* is a reverse neighbour map so the
!!! whole gradient is a loop over TARGET k in which each iteration writes only
!!! its own slice -- race-free with no atomics, and bit-identical at any thread
!!! count.

module hybrid_grass_m
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use hybrid_core_m, only: kinfo_t
  implicit none
  private
  !> grassmann_cg exit reasons.  GCG_OK is the ONLY one that means converged.
  !> GCG_TINYSTEP became UNREACHABLE in R2 A3(d): the vanishing-step failure
  !> rule it belonged to was the lam = 0 half of a lam-gated line search and is
  !> gone, so a failed search now always restarts from a fresh step and gives
  !> up as GCG_STALL.  The constant is kept because callers select on it and a
  !> stale .o linking against a renumbered enum would be silent.
  !> GCG_SMALLDEC is no longer emitted either (2026-10: CG stops on the
  !> gradient, see grassmann_cg); it is kept for the same reason.
  integer, parameter, public :: GCG_OK = 0, GCG_STALL = 1, GCG_TINYSTEP = 2, &
                                GCG_SMALLDEC = 3, GCG_CAP = 4
  !> The stopping rule (2026-10): the projected-gradient norm at or below
  !> max(GCG_GTOL_ABS, GCG_GTOL_REL*|grad0|) ends the solve with GCG_OK.
  !> Until then these two numbers only labelled the exit and the loop stopped
  !> after five consecutive small decrements (GCG_SMALLDEC).
  real(dp), parameter, public :: GCG_GTOL_ABS = 1.0e-10_dp, &
                                 GCG_GTOL_REL = 1.0e-6_dp
  public :: gcg_t, gcg_build, gcg_free, gcg_omega, gcg_check_gradient, &
            grassmann_cg, gcg_seed_from_full, gcg_to_full, gcg_seed_svd

  complex(dp), parameter :: CONE = (1.0_dp, 0.0_dp)
  complex(dp), parameter :: CZERO = (0.0_dp, 0.0_dp)

  !> Run constants for the Grassmann objective.  Built once; Mdd and Q never
  !> change during the solve.
  type :: gcg_t
     integer  :: nk = 0, nntot = 0, dmax = 0, amax = 0
     real(dp) :: aconst = 0.0_dp     ! frozen-frozen block of F
     real(dp) :: nwW = 0.0_dp        ! nw * sum_b w_b
     integer,     allocatable :: ndis(:), ntake(:)
     integer,     allocatable :: neigh(:,:)                 ! (nntot,nk)
     complex(dp), allocatable :: Mdd(:,:,:,:)   ! (dmax,dmax,nntot,nk), sqrt(w) scaled
     complex(dp), allocatable :: Q(:,:,:)       ! (dmax,dmax,nk)
     ! reverse neighbour map: for target m, the (k,b) pairs with neigh(b,k)==m
     integer,     allocatable :: rev_n(:)       ! (nk) how many
     integer,     allocatable :: rev_k(:,:), rev_b(:,:)     ! (maxrev,nk)
  end type gcg_t

contains

  !--------------------------------------------------------------- gcg_build
  !> Build the run constants from the SAME index sets build_kinfo produced.
  subroutine gcg_build(mmn, neigh, w, ki, npool, nntot, nk, nw, P)
    integer, intent(in) :: npool, nntot, nk, nw
    complex(dp), intent(in) :: mmn(npool,npool,nntot,nk)
    integer, intent(in) :: neigh(nntot,nk)
    real(dp), intent(in) :: w(nntot)
    type(kinfo_t), intent(in) :: ki(nk)
    type(gcg_t), intent(out) :: P

    integer :: k, b, kb, a, c, nd, ndb, nf, nfb, maxrev
    real(dp) :: sw
    complex(dp), allocatable :: Mf(:,:), X(:,:)

    P%nk = nk; P%nntot = nntot
    allocate(P%ndis(nk), P%ntake(nk))
    do k = 1, nk
       P%ndis(k) = ki(k)%ndis
       P%ntake(k) = ki(k)%ntake
    end do
    P%dmax = maxval(P%ndis)
    P%amax = maxval(P%ntake)
    P%nwW = real(nw,dp) * sum(w)
    allocate(P%neigh(nntot,nk)); P%neigh = neigh

    allocate(P%Mdd(P%dmax, P%dmax, nntot, nk)); P%Mdd = CZERO
    allocate(P%Q(P%dmax, P%dmax, nk));          P%Q = CZERO
    P%aconst = 0.0_dp

    do k = 1, nk
       nd = ki(k)%ndis; nf = ki(k)%nfroz
       do b = 1, nntot
          kb = neigh(b,k)
          ndb = ki(kb)%ndis; nfb = ki(kb)%nfroz
          sw = sqrt(w(b))

          ! T_dd : free-free, sqrt(w) folded in so the objective is a plain
          ! Frobenius norm of V'^H Mdd V'
          do c = 1, ndb
             do a = 1, nd
                P%Mdd(a,c,b,k) = sw * mmn(ki(k)%dis(a), ki(kb)%dis(c), b, k)
             end do
          end do

          ! T_df : V'_k against the NEIGHBOUR's frozen columns -> constant in k.
          ! A k without free rows (nd = 0: every pool state there is frozen or
          ! excluded by the admission stage) has no such block.  The call would
          ! pass N = 0 with a leading dimension of 0, which BLAS forbids
          ! (LDA >= max(1,N)); Accelerate aborts the process there.
          if (nfb > 0 .and. nd > 0) then
             allocate(Mf(nd, nfb))
             do c = 1, nfb
                do a = 1, nd
                   Mf(a,c) = sw * mmn(ki(k)%dis(a), ki(kb)%froz(c), b, k)
                end do
             end do
             call zherk('U','N', nd, nfb, 1.0_dp, Mf, nd, 1.0_dp, &
                        P%Q(1,1,k), P%dmax)
             deallocate(Mf)
          end if

          ! T_fd : the CENTRE's frozen rows against V'_{k+b} -> constant in k+b
          ! (nothing to add when the neighbour has no free rows)
          if (nf > 0 .and. ndb > 0) then
             allocate(X(nf, ndb))
             do c = 1, ndb
                do a = 1, nf
                   X(a,c) = sw * mmn(ki(k)%froz(a), ki(kb)%dis(c), b, k)
                end do
             end do
             call zherk('U','C', ndb, nf, 1.0_dp, X, nf, 1.0_dp, &
                        P%Q(1,1,kb), P%dmax)
             deallocate(X)
          end if

          ! T_ff : frozen-frozen, wholly constant
          do c = 1, nfb
             do a = 1, nf
                P%aconst = P%aconst + w(b) * &
                     real(conjg(mmn(ki(k)%froz(a), ki(kb)%froz(c), b, k)) * &
                                mmn(ki(k)%froz(a), ki(kb)%froz(c), b, k), dp)
             end do
          end do
       end do
    end do

    ! zherk fills only the upper triangle; mirror so Q is a full Hermitian
    ! matrix (the gradient uses it as a plain matmul operand)
    do k = 1, nk
       nd = P%ndis(k)
       do c = 1, nd
          do a = c+1, nd
             P%Q(a,c,k) = conjg(P%Q(c,a,k))
          end do
       end do
    end do

    ! reverse neighbour map -- GATHER, never scatter (see the header)
    allocate(P%rev_n(nk)); P%rev_n = 0
    do k = 1, nk
       do b = 1, nntot
          P%rev_n(neigh(b,k)) = P%rev_n(neigh(b,k)) + 1
       end do
    end do
    maxrev = maxval(P%rev_n)
    allocate(P%rev_k(maxrev,nk), P%rev_b(maxrev,nk))
    P%rev_n = 0
    do k = 1, nk
       do b = 1, nntot
          kb = neigh(b,k)
          P%rev_n(kb) = P%rev_n(kb) + 1
          P%rev_k(P%rev_n(kb), kb) = k
          P%rev_b(P%rev_n(kb), kb) = b
       end do
    end do
  end subroutine gcg_build

  subroutine gcg_free(P)
    type(gcg_t), intent(inout) :: P
    if (allocated(P%Mdd))   deallocate(P%Mdd)
    if (allocated(P%Q))     deallocate(P%Q)
    if (allocated(P%ndis))  deallocate(P%ndis)
    if (allocated(P%ntake)) deallocate(P%ntake)
    if (allocated(P%neigh)) deallocate(P%neigh)
    if (allocated(P%rev_n)) deallocate(P%rev_n)
    if (allocated(P%rev_k)) deallocate(P%rev_k)
    if (allocated(P%rev_b)) deallocate(P%rev_b)
  end subroutine gcg_free

  !------------------------------------------------------------------ F, Omega
  !> F(V) -- the quantity to MAXIMIZE.
  function gcg_F(P, V) result(F)
    type(gcg_t), intent(in) :: P
    complex(dp), intent(in) :: V(P%dmax, P%amax, P%nk)
    real(dp) :: F
    integer :: k, b, kb, nd, na, ndb, nab
    complex(dp), allocatable :: Y(:,:), T(:,:), QV(:,:)
    real(dp), allocatable :: Fk(:)
    ! Per-k partials + a SERIAL in-k-order sum -- the omega_all pattern
    ! (hybrid_core.F).  F gates line-search acceptance, so a reduction(+:...)
    ! would make the accepted step, and hence every subsequent iterate, depend
    ! on the thread count.  Scratch is allocated INSIDE the loop: an outside
    ! allocation would be shared across threads under the pragma.
    allocate(Fk(P%nk)); Fk = 0.0_dp
    !$omp parallel do default(shared) schedule(dynamic) &
    !$omp   private(k, b, kb, nd, na, ndb, nab, Y, T, QV)
    do k = 1, P%nk
       nd = P%ndis(k); na = P%ntake(k)
       if (na == 0) cycle
       allocate(Y(P%dmax,P%amax), T(P%amax,P%amax), QV(P%dmax,P%amax))
       ! Tr(V^H Q V)
       call zgemm('N','N', nd, na, nd, CONE, P%Q(1,1,k), P%dmax, &
                  V(1,1,k), P%dmax, CZERO, QV, P%dmax)
       Fk(k) = Fk(k) + sum(real(conjg(V(1:nd,1:na,k)) * QV(1:nd,1:na), dp))
       ! sum_b ||V_k^H Mdd V_{k+b}||^2
       do b = 1, P%nntot
          kb = P%neigh(b,k)
          ndb = P%ndis(kb); nab = P%ntake(kb)
          if (nab == 0) cycle
          call zgemm('N','N', nd, nab, ndb, CONE, P%Mdd(1,1,b,k), P%dmax, &
                     V(1,1,kb), P%dmax, CZERO, Y, P%dmax)
          call zgemm('C','N', na, nab, nd, CONE, V(1,1,k), P%dmax, &
                     Y, P%dmax, CZERO, T, P%amax)
          Fk(k) = Fk(k) + sum(real(T(1:na,1:nab))**2 + aimag(T(1:na,1:nab))**2)
       end do
       deallocate(Y, T, QV)
    end do
    !$omp end parallel do
    F = P%aconst + sum(Fk)
    deallocate(Fk)
  end function gcg_F

  !> The TRUE Omega_I, always -- never the penalised value.  Callers report
  !> this, and it stays comparable across lam.
  function gcg_omega(P, V) result(om)
    type(gcg_t), intent(in) :: P
    complex(dp), intent(in) :: V(P%dmax, P%amax, P%nk)
    real(dp) :: om
    om = P%nwW - gcg_F(P, V) / real(P%nk, dp)
  end function gcg_omega

  !------------------------------------------------------------------ gradient
  !> dF/dV* at every k, evaluated from a SINGLE point.
  !>
  !> Loop is over TARGET k and each iteration writes only G(:,:,k), so it is
  !> race-free and bit-identical at any thread count.  Everything -- the
  !> neighbour gather, both Zd terms, and applying Zd to V -- happens inside;
  !> in the prototype, leaving the batched pieces outside capped scaling at
  !> 1.21x on 6 threads (Amdahl-limited by code structure, not the algorithm)
  !> and moving them in took it to 2.30x.
  subroutine gcg_grad(P, V, G)
    type(gcg_t), intent(in) :: P
    complex(dp), intent(in) :: V(P%dmax, P%amax, P%nk)
    complex(dp), intent(out) :: G(P%dmax, P%amax, P%nk)
    ! NB: `gmu`, not `g` -- Fortran is case-insensitive and `G` is the gradient.
    integer :: k, b, kb, r, kp, bp, nd, na, ndb, nab, ndp, nap, a, c, gmu
    complex(dp), allocatable :: Zd(:,:), Y(:,:), Bmat(:,:)

    G = CZERO
    !$omp parallel do default(shared) schedule(dynamic) &
    !$omp   private(k, b, kb, r, kp, bp, nd, na, ndb, nab, ndp, nap, a, c, &
    !$omp           gmu, Zd, Y, Bmat)
    do k = 1, P%nk
       nd = P%ndis(k); na = P%ntake(k)
       if (na == 0) cycle
       allocate(Zd(P%dmax,P%dmax), Y(P%dmax,P%amax), Bmat(P%dmax,P%amax))
       Zd = CZERO

       ! CENTRE: V_k as the LEFT factor of ||V_k^H Mdd V_{k+b}||^2
       do b = 1, P%nntot
          kb = P%neigh(b,k)
          ndb = P%ndis(kb); nab = P%ntake(kb)
          if (nab == 0) cycle
          call zgemm('N','N', nd, nab, ndb, CONE, P%Mdd(1,1,b,k), P%dmax, &
                     V(1,1,kb), P%dmax, CZERO, Y, P%dmax)
          call zgemm('N','C', nd, nd, nab, CONE, Y, P%dmax, Y, P%dmax, &
                     CONE, Zd, P%dmax)
       end do

       ! NEIGHBOUR: V_k as the RIGHT factor, gathered over (k',b) -> k.
       ! This is the term the prototype originally omitted.
       do r = 1, P%rev_n(k)
          kp = P%rev_k(r,k); bp = P%rev_b(r,k)
          ndp = P%ndis(kp); nap = P%ntake(kp)
          if (nap == 0) cycle
          ! B = Mdd(kp,bp)^H V_kp   (rows restricted to dis(k) by construction)
          call zgemm('C','N', nd, nap, ndp, CONE, P%Mdd(1,1,bp,kp), P%dmax, &
                     V(1,1,kp), P%dmax, CZERO, Bmat, P%dmax)
          call zgemm('N','C', nd, nd, nap, CONE, Bmat, P%dmax, Bmat, P%dmax, &
                     CONE, Zd, P%dmax)
       end do

       Zd(1:nd,1:nd) = Zd(1:nd,1:nd) + P%Q(1:nd,1:nd,k)
       call zgemm('N','N', nd, na, nd, CONE, Zd, P%dmax, V(1,1,k), P%dmax, &
                  CZERO, G(1,1,k), P%dmax)
       deallocate(Zd, Y, Bmat)
    end do
    !$omp end parallel do
  end subroutine gcg_grad

  !> Grassmann tangent projection (I - V V^H) G, per k on the true block.
  subroutine gcg_project(P, V, G)
    type(gcg_t), intent(in) :: P
    complex(dp), intent(in) :: V(P%dmax, P%amax, P%nk)
    complex(dp), intent(inout) :: G(P%dmax, P%amax, P%nk)
    integer :: k, nd, na
    complex(dp), allocatable :: VG(:,:)
    ! Own-slice in-place update of G(:,:,k): race-free, no reduction, results
    ! independent of thread count.  Scratch allocated inside the loop
    ! (thread-private); an outside allocation would be shared.
    !$omp parallel do default(shared) schedule(dynamic) private(k, nd, na, VG)
    do k = 1, P%nk
       nd = P%ndis(k); na = P%ntake(k)
       if (na == 0) cycle
       allocate(VG(P%amax,P%amax))
       call zgemm('C','N', na, na, nd, CONE, V(1,1,k), P%dmax, &
                  G(1,1,k), P%dmax, CZERO, VG, P%amax)
       call zgemm('N','N', nd, na, na, -CONE, V(1,1,k), P%dmax, &
                  VG, P%amax, CONE, G(1,1,k), P%dmax)
       deallocate(VG)
    end do
    !$omp end parallel do
  end subroutine gcg_project

  !> QR retraction: W = qf(V + t D), with the sign convention pinned so the
  !> retraction is a function of the point, not of LAPACK's phase choice.
  subroutine gcg_retract(P, V, D, t, W)
    type(gcg_t), intent(in) :: P
    complex(dp), intent(in) :: V(P%dmax,P%amax,P%nk), D(P%dmax,P%amax,P%nk)
    real(dp), intent(in) :: t
    complex(dp), intent(out) :: W(P%dmax,P%amax,P%nk)
    integer :: k, nd, na, c, info, lwork
    complex(dp), allocatable :: A(:,:), tau(:), work(:), sgn(:)
    W = CZERO
    lwork = 64*max(P%amax,1)
    ! Own-slice writes to W(:,:,k); each k's QR is an independent
    ! deterministic LAPACK call, so results are bit-identical at any thread
    ! count.  Workspaces allocated inside the loop (thread-private).
    !$omp parallel do default(shared) schedule(dynamic) &
    !$omp   private(k, nd, na, c, info, A, tau, work, sgn)
    do k = 1, P%nk
       nd = P%ndis(k); na = P%ntake(k)
       if (na == 0) cycle
       allocate(A(P%dmax,P%amax), tau(P%amax), work(lwork), sgn(P%amax))
       A(1:nd,1:na) = V(1:nd,1:na,k) + t * D(1:nd,1:na,k)
       call zgeqrf(nd, na, A, P%dmax, tau, work, lwork, info)
       if (info /= 0) error stop 'gcg_retract: zgeqrf failed'
       ! Capture sign(diag R) BEFORE zungqr overwrites the upper triangle.
       ! Without this the retraction inherits LAPACK's arbitrary column phase,
       ! so the same point could map to differently-signed frames on different
       ! calls -- harmless for Omega_I (which sees only the subspace) but it
       ! makes the line search and the FD check non-reproducible.
       do c = 1, na
          if (real(A(c,c),dp) < 0.0_dp) then
             sgn(c) = -CONE
          else
             sgn(c) = CONE
          end if
       end do
       call zungqr(nd, na, na, A, P%dmax, tau, work, lwork, info)
       if (info /= 0) error stop 'gcg_retract: zungqr failed'
       do c = 1, na
          W(1:nd,c,k) = sgn(c) * A(1:nd,c)
       end do
       deallocate(A, tau, work, sgn)
    end do
    !$omp end parallel do
  end subroutine gcg_retract

  function gcg_ip(P, A, B) result(s)
    type(gcg_t), intent(in) :: P
    complex(dp), intent(in) :: A(P%dmax,P%amax,P%nk), B(P%dmax,P%amax,P%nk)
    real(dp) :: s
    s = sum(real(conjg(A) * B, dp))
  end function gcg_ip

  !------------------------------------------------- finite-difference guard
  !> Checked at EVERY run, not on request.  A missing neighbour term still
  !> descends plausibly and is invisible without this.
  function gcg_check_gradient(P, V) result(rel)
    type(gcg_t), intent(in) :: P
    complex(dp), intent(in) :: V(P%dmax,P%amax,P%nk)
    real(dp) :: rel
    complex(dp), allocatable :: G(:,:,:), E(:,:,:), W(:,:,:), Vc(:,:,:)
    real(dp) :: fd, an, h, f1, f0, rr, ri
    integer :: i, j, k, seed
    allocate(G(P%dmax,P%amax,P%nk), E(P%dmax,P%amax,P%nk), &
             W(P%dmax,P%amax,P%nk), Vc(P%dmax,P%amax,P%nk))
    call gcg_grad(P, V, G)
    call gcg_project(P, V, G)
    ! deterministic pseudo-random tangent direction (no RANDOM_NUMBER, so the
    ! guard is reproducible run to run)
    seed = 1
    do k = 1, P%nk
       do j = 1, P%amax
          do i = 1, P%dmax
             seed = mod(1103515245*seed + 12345, 2147483647)
             rr = real(seed,dp)/2147483647.0_dp - 0.5_dp
             seed = mod(1103515245*seed + 12345, 2147483647)
             ri = real(seed,dp)/2147483647.0_dp - 0.5_dp
             E(i,j,k) = cmplx(rr, ri, dp)
          end do
       end do
    end do
    call gcg_project(P, V, E)
    Vc = V
    h = 1.0e-6_dp
    call gcg_retract(P, Vc, E,  h, W); f1 = gcg_F(P, W)
    call gcg_retract(P, Vc, E, -h, W); f0 = gcg_F(P, W)
    fd = (f1 - f0) / (2.0_dp*h)
    ! d/deps f(V+eps E) = 2 Re<E, df/dV*> for ANY real f, not just a quadratic
    ! one -- so the factor 2 is the Wirtinger convention.
    an = 2.0_dp * gcg_ip(P, G, E)
    rel = abs(fd - an) / max(abs(fd), 1.0e-30_dp)
    write(*,'(a,es16.8,a,es16.8,a,f10.6)') '   [gradchk] FD ', fd, &
         '  analytic ', an, '  an/FD ', an/sign(max(abs(fd),1.0e-30_dp), fd)
    deallocate(G, E, W, Vc)
  end function gcg_check_gradient

  !--------------------------------------------------------------- the solver
  !> EXIT REASON, and why it is not optional.  This routine stops when the
  !> gradient falls below tolerance (GCG_OK), on a line-search stall
  !> (GCG_STALL), or at the iteration cap (GCG_CAP); step underflow and the
  !> small-decrement rule of earlier versions are gone.
  !> Until now it reported none of them, and its caller printed "converged
  !> after N gradient evaluations" unconditionally for all five.  That made
  !> every eval count in the project's cost tables unfalsifiable: a run that
  !> stalled, a run that hit the cap, and a run that genuinely converged were
  !> typographically identical.  The Rh3Si4 lam=2048 row (3907 evals) sits at
  !> 97.7% of its 4000 cap and cannot now be told apart from a capped run.
  !> `why` fixes that; GCG_OK is the only value that means converged.
  subroutine grassmann_cg(P, V, niter, tol, rel_tol, om, nit, verbose, why, &
                          gnorm, gnorm0)
    type(gcg_t), intent(in) :: P
    complex(dp), intent(inout) :: V(P%dmax,P%amax,P%nk)
    integer, intent(in) :: niter
    !> tol, rel_tol: no longer read (CG stops on the gradient since 2026-10);
    !> kept so the interface matches the Gauss-Seidel solver's.
    real(dp), intent(in) :: tol, rel_tol
    real(dp), intent(out) :: om
    integer, intent(out) :: nit
    logical, intent(in) :: verbose
    !> GCG_OK | GCG_STALL | GCG_TINYSTEP | GCG_SMALLDEC | GCG_CAP
    integer, intent(out), optional :: why
    !> R10-4a: the projected-gradient norm at the returned V (gnorm) and at
    !> entry (gnorm0) -- the measured quantity the GCG_OK label tests, for a
    !> caller that wants to print it.
    real(dp), intent(out), optional :: gnorm, gnorm0

    complex(dp), allocatable :: G(:,:,:), Gn(:,:,:), D(:,:,:), W(:,:,:)
    ! NB: `gnv`/`gnv0`, not `gn`/`gn0` -- Fortran is case-insensitive and
    ! `Gn` above is already the next-gradient array.
    real(dp) :: t, gg, ggn, beta, omw, t_init, gnv, gnv0, gtol
    integer :: it, ls, nfail, iwhy
    logical :: ok

    allocate(G(P%dmax,P%amax,P%nk), Gn(P%dmax,P%amax,P%nk), &
             D(P%dmax,P%amax,P%nk), W(P%dmax,P%amax,P%nk))
    om = gcg_omega(P, V)
    call gcg_grad(P, V, G); call gcg_project(P, V, G)
    D = G
    gg = gcg_ip(P, G, G)
    gnv0 = sqrt(gg)             ! R10-4a: the entry gradient norm, for `why`
    t = 1.0_dp / max(maxval(abs(P%Q)), 1.0e-12_dp)
    t_init = t; nfail = 0
    gtol = max(GCG_GTOL_ABS, GCG_GTOL_REL*gnv0)
    iwhy = GCG_CAP              ! falling out of the loop = hit the cap
    nit = 1

    do it = 1, niter
       ok = .false.
       do ls = 1, 20              ! backtracking on the DESCENT of the objective
          call gcg_retract(P, V, D, t, W)
          omw = gcg_omega(P, W)
          if (omw < om - 1.0e-14_dp*abs(om)) then
             ! R10-4a: `iwhy` is NOT set here any more.  It used to be, which
             ! made GCG_CAP dead code (any accepted step, i.e. almost every
             ! real run, latched OK and nothing after this loop could unlatch
             ! it) -- the actual exit reason below is what `iwhy` now tracks.
             ok = .true.; exit
          end if
          t = 0.5_dp*t
       end do
       if (.not. ok) then
          ! ONE line-search failure rule for every lambda (R2 A3(d)).  A
          ! failed search restarts on steepest ascent from a FRESH step
          ! t_init*0.1^nfail -- NOT from the already-halved t, which after 20
          ! backtracks is t/2^20 and drives the next search under 1e-18 after
          ! three failures, quitting the solve early (measured on MgB2 12x12x8
          ! as a spurious Omega_I jump to ~13.0 with the band asymmetry WORSE
          ! than the unpenalised run, at a third of the wall time).
          !
          ! This branch USED TO BE GATED ON lam > 0, with lam = 0 taking the
          ! vanishing-step rule instead.  That made lam = 0 a DIFFERENT
          ! ALGORITHM from lam > 0 and silently confounded every free-vs-
          ! penalised comparison in the record: the Rh3Si4 "80 evaluations,
          ! Omega_I 376.155" free run was compared against a lam = 2048 run
          ! under the other rule and read as a 49x penalty cost, of which
          ! 42.5x was the rule and 1.15x the penalty.
          if (nfail >= 3) then; iwhy = GCG_STALL; exit; end if
          nfail = nfail + 1
          D = G
          t = t_init * (0.1_dp**nfail)
          cycle
       end if
       nfail = 0
       V = W; om = omw
       t = 1.8_dp*t               ! grow the step while it keeps working
       call gcg_grad(P, V, Gn); call gcg_project(P, V, Gn)
       nit = nit + 1
       ggn = gcg_ip(P, Gn, Gn)
       beta = max(0.0_dp, (ggn - gcg_ip(P, Gn, G)) / max(gg, 1.0e-300_dp))
       D = Gn + beta*D
       if (gcg_ip(P, Gn, D) <= 0.0_dp) D = Gn      ! keep it an ascent direction
       G = Gn; gg = ggn
       if (verbose .and. mod(it,25) == 0) then
          write(*,'(a,i6,a,f18.10)') '   cg iter', it, '  Omega_I=', om
       end if
       ! Stopping rule (2026-10): the projected gradient, |grad| <=
       ! max(1e-10, 1e-6 |grad0|).  The earlier rule, five consecutive
       ! decrements below tol / rel_tol (GCG_SMALLDEC), stopped while the
       ! subspace was still moving: on LCAO MgB2 12^3 it exited after 106
       ! evaluations with Omega_I 1.1e-4 A^2 above the minimum and the subspace
       ! 0.14 (max |dP|) away from it; the gradient test stops after 655 within
       ! 1e-3 of it.  tol and rel_tol no longer stop CG (they still govern the
       ! Gauss-Seidel sweep).
       if (sqrt(gg) <= gtol) then; iwhy = GCG_OK; exit; end if
    end do
    ! R10-4a: G/gg are ALREADY the tangent-projected gradient at the V being
    ! returned, on every exit above -- STALL leaves V, G and gg untouched since
    ! the last accepted step, and a CAP fall-through's last accepted step set
    ! them too -- so no extra gcg_grad/gcg_project call is needed here.  A
    ! stall or the cap with the gradient already below gtol is also a
    ! converged minimum.
    gnv = sqrt(gg)
    if (gnv <= gtol) iwhy = GCG_OK
    if (present(gnorm))  gnorm  = gnv
    if (present(gnorm0)) gnorm0 = gnv0
    if (present(why)) why = iwhy
    om = gcg_omega(P, V)        ! report the TRUE Omega_I, never the penalised
    deallocate(G, Gn, D, W)
  end subroutine grassmann_cg

  !------------------------------------------------------------------- seeding
  !> Project-optimal SVD seed for the FREE block -- the SAME seed
  !> disentangle_frozen uses, so a CG-vs-GS comparison starts from one point.
  !>
  !> This must not be skipped.  V is required to be ORTHONORMAL per k: the
  !> tangent projection (I - V V^H) is only a projection when V^H V = I, and
  !> the QR retraction is only first-order accurate there.  Seeding V with
  !> anything else (an uninitialized array, say) silently breaks BOTH, and the
  !> symptom is a finite-difference check that fails by a single common factor
  !> across every term of the objective -- which is exactly how this was found.
  subroutine gcg_seed_svd(P, ki, amn, npool, nproj, V)
    use hybrid_core_m, only: svd_left
    type(gcg_t), intent(in) :: P
    type(kinfo_t), intent(in) :: ki(P%nk)
    integer, intent(in) :: npool, nproj
    complex(dp), intent(in) :: amn(npool, nproj, P%nk)
    complex(dp), intent(out) :: V(P%dmax, P%amax, P%nk)
    integer :: k, a, c, nd, nt
    complex(dp), allocatable :: Aseed(:,:), Useed(:,:)
    V = CZERO
    do k = 1, P%nk
       nd = P%ndis(k); nt = P%ntake(k)
       if (nt == 0) cycle
       allocate(Aseed(nd,nproj), Useed(nd,nt))
       do a = 1, nd
          Aseed(a,:) = amn(ki(k)%dis(a), :, k)
       end do
       call svd_left(Aseed, nd, nproj, nt, Useed)
       do c = 1, nt
          do a = 1, nd
             V(a,c,k) = Useed(a,c)
          end do
       end do
       deallocate(Aseed, Useed)
    end do
  end subroutine gcg_seed_svd

  !------------------------------------------------- full <-> free conversions
  !> Extract the FREE block of the full (npool,nw,nk) frame.
  subroutine gcg_seed_from_full(P, ki, Vfull, npool, nw, V)
    type(gcg_t), intent(in) :: P
    type(kinfo_t), intent(in) :: ki(P%nk)
    integer, intent(in) :: npool, nw
    complex(dp), intent(in) :: Vfull(npool,nw,P%nk)
    complex(dp), intent(out) :: V(P%dmax,P%amax,P%nk)
    integer :: k, a, c
    V = CZERO
    do k = 1, P%nk
       do c = 1, P%ntake(k)
          do a = 1, P%ndis(k)
             V(a,c,k) = Vfull(ki(k)%dis(a), ki(k)%nfroz+c, k)
          end do
       end do
    end do
  end subroutine gcg_seed_from_full

  !> Scatter the free block back, re-pinning the frozen unit columns.
  subroutine gcg_to_full(P, ki, V, npool, nw, Vfull)
    type(gcg_t), intent(in) :: P
    type(kinfo_t), intent(in) :: ki(P%nk)
    integer, intent(in) :: npool, nw
    complex(dp), intent(in) :: V(P%dmax,P%amax,P%nk)
    complex(dp), intent(out) :: Vfull(npool,nw,P%nk)
    integer :: k, a, c
    Vfull = CZERO
    do k = 1, P%nk
       do c = 1, ki(k)%nfroz
          Vfull(ki(k)%froz(c), c, k) = CONE
       end do
       do c = 1, P%ntake(k)
          do a = 1, P%ndis(k)
             Vfull(ki(k)%dis(a), ki(k)%nfroz+c, k) = V(a,c,k)
          end do
       end do
    end do
  end subroutine gcg_to_full

end module hybrid_grass_m
