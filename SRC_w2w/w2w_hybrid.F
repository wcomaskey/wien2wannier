!!! w2w_hybrid.F -- self-contained wien2wannier hybrid pre-Wannierization tool.
!!! Reads the standard W90 files w2w emits (<seed>.amn/.mmn/.eig/.nnkp), builds
!!! per-state projectability from |amn|^2, runs the full validated Stage-2
!!! recipe (masks -> frozen-constrained Omega_I disentanglement -> rotate to
!!! projected-H eigenbasis -> AO/channel SCDM gauge) and writes an isolated
!!! Nb=Nw stock-Wannier90 hand-off <out>.amn/.mmn/.eig.
!!!
!!! Usage:  w2w_hybrid --seed <case> --nw <N> [options]
!!!         (the first two bare arguments are also accepted as seed and nw:
!!!          `w2w_hybrid <case> <N> [options]`)
!!! Run `w2w_hybrid --help` for the complete option and environment-variable
!!! reference.  The method is that of [Comaskey2026]; the interface it extends
!!! is wien2wannier [Kunes2010] for WIEN2k [Blaha2020] (keys of REFERENCES.md
!!! in the root directory).
program w2w_hybrid
  use, intrinsic :: iso_fortran_env, only: dp => real64, output_unit
  use hybrid_files_m
  use hybrid_prep_m
  use hybrid_core_m
  use hybrid_io_m
  use hybrid_grass_m
  implicit none
  ! gauge mode (arg 14): 'scdm' (default) | 'ptg' = SCDM anchor at k=1
  ! then Procrustes/parallel-transport propagation across the BZ;
  ! 'sms<N>'/'pts<N>' = scdm/ptg + N smoothing sweeps; 'lowdin' = plain
  ! Loewdin projection gauge (leading nw physical channels as trials,
  ! W90 orthonormalizes), 'lds<N>' = lowdin + N smoothing sweeps;
  ! 'bscdm' = block-SCDM (per-isolated-manifold anchors), 'bsms<N>' =
  ! block-SCDM + N smoothing sweeps
  complex(dp), parameter :: CONE = (1.0_dp,0.0_dp), CZERO = (0.0_dp,0.0_dp)

  character(len=256) :: seed, outp, arg, pamn_file
  character(len=256) :: cont_chk_seed
  character(len=16) :: rescale, regime, femode, regime_out, gauge_mode
  character(len=8) :: solver_sel
  integer :: whyg
  real(dp) :: gng, gng0
  ! the disentanglement outcome and the subspace continuity, for the log and <out>.dis
  character(len=320) :: dis_status, cont_line
  real(dp) :: cont_smin
  integer :: cont_k, cont_kb
  integer :: nw, nb, nk, nntot, nproj, nb2, nk2, nntot2, nkq, npool, i, j, k, ib
  integer :: dis_niter, dis_check, dis_aa, np_mask, nprojP, froz_dmd
  integer :: fl_add, fl_drop, fl_nk
  real(dp) :: fe_used
  logical :: froz_legacy
  logical :: init_transp, admit_all, no_froz, dis_verb, pfloor_set
  logical :: window_route
  logical :: dis_colored, dis_adaptive
  real(dp) :: e_fermi, femax_val, omega, best_svmin, t0, t1, dis_rtol, dis_mix
  real(dp) :: pool_factor, p_floor
  integer(8) :: wc0, wc1, wrate
  complex(dp), allocatable :: amn(:,:,:), mmn(:,:,:,:), amn_p(:,:,:), mmn_p(:,:,:,:)
  complex(dp), allocatable :: amnP(:,:,:)
  complex(dp), allocatable :: V(:,:,:), W(:,:,:), mmn_rot(:,:,:,:), B(:,:,:), A_gauge(:,:,:)
  real(dp), allocatable :: eig(:,:), eig_p(:,:), eig_t(:,:), P(:,:), bcart(:,:), wgt(:)
  integer, allocatable :: neigh(:,:), gshift(:,:,:), pool(:), anchors(:)
  ! per-(k,b) b-vectors (nnkp_block_bvectors) for the gauge policy's spreads
  real(dp), allocatable :: bkb(:,:,:)
  logical, allocatable :: admit(:,:), frozen(:,:), pcol(:)
  real(dp) :: wresid
  ! anatomy of the mask rule, for <out>.masks / <out>.anatomy
  real(dp), allocatable :: P_used(:,:), an_avg(:), an_min(:), an_bot(:)
  real(dp) :: an_plateau, an_fid, an_pgate, an_pfroz
  logical :: an_ok
  ! pocket skip (hybrid_prep_m::pocket_ceiling): the first pass is saved while a
  ! lowered ceiling is tried, and restored if the try does not clear POCKET_SV
  logical :: pocket_on, pk_try, pk_fired
  integer :: npocket, cont_k_s, cont_kb_s
  real(dp) :: fe_try, fe_s, omega_s, cont_smin_s
  character(len=400) :: pk_msg
  character(len=320) :: dis_status_s
  character(len=16) :: femode_try
  real(dp), allocatable :: sv_links(:,:), eps_pool(:,:), P_pool(:,:), eig_t_s(:,:)
  complex(dp), allocatable :: V_s(:,:,:), W_s(:,:,:), mmn_rot_s(:,:,:,:)
  logical, allocatable :: admit_s(:,:), frozen_s(:,:)
  integer, allocatable :: pool_s(:)

  ! ---- defaults (every option that has one) ----
  seed = ''
  nw = -1
  outp = ''
  rescale = 'occ99'; regime = 'auto'; femode = 'auto'; femax_val = 8.0_dp
  e_fermi = 0.0_dp
  dis_niter = 20000; dis_rtol = 1.0e-10_dp; dis_check = 10
  dis_colored = .false.; dis_adaptive = .false.
  dis_mix = 1.0_dp
  gauge_mode = 'scdm'
  dis_aa = 0
  np_mask = 0
  init_transp = .false.
  admit_all = .false.; no_froz = .false.; pfloor_set = .false.
  pool_factor = 1.5_dp
  p_floor = 0.10_dp
  solver_sel = 'cg'
  window_route = .false.
  pocket_on = .true.
  cont_chk_seed = ''
  call parse_cli()
  if (len_trim(cont_chk_seed) > 0) then
     call window_chk_continuity(trim(cont_chk_seed))
     call exit(0)     ! as util_w2w.F: no STOP, so no IEEE-flag note from the SVDs
  end if
  if (len_trim(outp) == 0) outp = trim(seed)//'-hybF'

  if (np_mask == 0) then
     ! default: full P (all amn columns) -- the validated champion paths
     ! (Fe 1.0002/WF, aSn 3.775/WF) used the complete candidate set.  The
     ! otsu mask boundary recorded by cliff-mode AUTO ('# mask_cols N' in
     ! projmap.dat) is reported as information only; pass --np-mask > 0 to
     ! adopt it explicitly.
     block
       integer :: up, iosj, nmc
       character(len=16) :: t1, t2
       open(newunit=up, file='projmap.dat', status='old', action='read', &
            iostat=iosj)
       if (iosj == 0) then
          read(up,'(a)',iostat=iosj) t1
          if (iosj == 0) then
             read(up,*,iostat=iosj) t1, t2, nmc
             if (iosj == 0 .and. trim(t1) == '#' .and. &
                 trim(t2) == 'mask_cols') write(*,'(a,i0,a)') &
                 '  projmap mask hint: ', nmc, ' columns (not adopted; np_mask=0 = full P)'
          end if
          close(up)
       end if
     end block
  end if

  ! env W2W_VERBOSE=1: per-sweep Omega_I trace from the disentangler.
  ! Hour-scale runs (Rh3Si4 at nw>200) are otherwise completely silent
  ! until they finish -- no progress, no sweep count, no ETA.
  block
    character(len=8) :: vbenv
    call get_environment_variable('W2W_VERBOSE', vbenv)
    dis_verb = (trim(vbenv) == '1')
    if (dis_verb) write(*,'(a)') '  verbose: per-sweep Omega_I trace ON'
  end block

  ! ---- read inputs ----
  call system_clock(wc0, wrate)
  call read_amn(trim(seed)//'.amn', amn, nb, nk, nproj)
  call read_mmn(trim(seed)//'.mmn', mmn, neigh, gshift, nb2, nk2, nntot)
  call read_eig(trim(seed)//'.eig', nb, nk, eig)
  call read_nnkp(trim(seed)//'.nnkp', nkq, nntot2, bcart)
  call system_clock(wc1)
  write(*,'(a,f8.2,a)') '  inputs read in ', real(wc1-wc0,dp)/real(wrate,dp), ' s'
  if (nb2 /= nb .or. nk2 /= nk) error stop 'amn/mmn dims inconsistent'
  if (nntot2 /= nntot) error stop 'nnkp/mmn nntot mismatch'
  write(*,'(a,4(a,i0))') trim(seed), ':  nb=',nb,' nk=',nk,' nntot=',nntot,' nproj=',nproj
  write(*,'(a,i0)') '  num_wann = ', nw

  ! ---- projectability from |amn|^2 ; b-weights ----
  ! env W2W_PAMN=<file>: read a SEPARATE amn as the projectability source
  ! (e.g. tail-extended full-cell projectors for P/masks while the plain
  ! MT amn keeps feeding the SCDM gauge and the hand-off).
  call get_environment_variable('W2W_PAMN', pamn_file)
  if (len_trim(pamn_file) > 0) then
     call read_amn(trim(pamn_file), amnP, i, j, nprojP)
     if (i /= nb .or. j /= nk) error stop 'W2W_PAMN: nb/nk mismatch with seed amn'
     write(*,'(a,a,a,i0,a)') '  P source: ', trim(pamn_file), ' (', nprojP, ' columns)'
  else
     amnP = amn
     nprojP = nproj
  end if
  allocate(P(nb,nk), pcol(nprojP))
  pcol = .true.
  if (np_mask > 0 .and. np_mask < nprojP) then
     ! SOC amn files hold two sigma blocks (nproj = 2*nchan, projmap has
     ! nchan rows): mirror the leading-block selection into the second
     ! block, else the P span is silently sigma-asymmetric.
     block
       integer :: up, nrows, iosj, i1, i2, i3, i4
       character(len=256) :: linep
       nrows = 0
       open(newunit=up, file='projmap.dat', status='old', action='read', iostat=iosj)
       if (iosj == 0) then
          do
             read(up,'(a)',iostat=iosj) linep
             if (iosj /= 0) exit
             if (index(adjustl(linep), '#') == 1 .or. len_trim(linep) == 0) cycle
             read(linep,*,iostat=iosj) i1, i2, i3, i4
             if (iosj /= 0) exit
             nrows = nrows + 1
          end do
          close(up)
       end if
       pcol = .false.
       pcol(1:np_mask) = .true.
       if (nrows > 0 .and. nprojP == 2*nrows .and. np_mask <= nrows) then
          pcol(nrows+1:nrows+np_mask) = .true.
          write(*,'(a,i0,a,i0,a)') '  P from leading ', np_mask, &
               ' channels x 2 sigma blocks of ', nprojP, ' amn columns'
       else
          write(*,'(a,i0,a,i0,a)') '  P from leading ', np_mask, ' of ', &
               nprojP, ' amn columns (mask/gauge decoupled)'
       end if
     end block
  else if (np_mask < 0) then
     ! P from columns with l <= |np_mask| (channel map from projmap.dat);
     ! restores the sp-only mask profile even when d/f candidates are
     ! present for the gauge (the historical paper-model behaviour)
     block
       integer :: up, colj, iatj, lj, mj, iosj, nrows
       character(len=256) :: linep
       open(newunit=up, file='projmap.dat', status='old', action='read', &
            iostat=iosj)
       if (iosj /= 0) error stop 'np_mask<0 needs projmap.dat in cwd'
       nrows = 0
       do
          read(up,'(a)',iostat=iosj) linep
          if (iosj /= 0) exit
          if (index(adjustl(linep), '#') == 1 .or. len_trim(linep) == 0) cycle
          read(linep,*,iostat=iosj) iatj, lj, mj
          if (iosj /= 0) exit
          nrows = nrows + 1
          if (nrows > nprojP) exit
          pcol(nrows) = (lj <= abs(np_mask))
       end do
       close(up)
       ! SOC sigma-block mirroring (projmap lists spatial channels once)
       if (nrows > 0 .and. nprojP == 2*nrows) then
          do colj = 1, nrows
             pcol(nrows+colj) = pcol(colj)
          end do
       end if
     end block
     write(*,'(a,i0,a,i0,a,i0)') '  P from l<=', abs(np_mask), &
          ' channels: ', count(pcol), ' of ', nprojP
  end if
  if (trim(rescale) == 'lowdin') then
     ! Loewdin-orthonormalised projectability (shared with the HYBR
     ! one-shot; see hybrid_prep_m::lowdin_pstate)
     call lowdin_pstate(amnP, nb, nk, nprojP, pcol, P)
     write(*,'(a,f8.4,a,f8.4)') '  Lowdin P: max = ', maxval(P), &
          '   mean = ', sum(P)/real(nb*nk,dp)
     rescale = 'none'   ! already normalised; build_masks just clips
  else
  do k = 1, nk; do i = 1, nb
     P(i,k) = 0.0_dp
     do j = 1, nprojP
        if (pcol(j)) P(i,k) = P(i,k) + real(amnP(i,j,k),dp)**2 + aimag(amnP(i,j,k))**2
     end do
  end do; end do
  end if
  if (trim(rescale) == 'share') then
     ! per-state SHARE statistic: mask-column P as a fraction of the
     ! all-column P.  Scale-invariant per state, so it undoes both the
     ! interstitial suppression of MT-projected P (sp valence in open
     ! structures) and the LO overcounting inflation -- this is the
     ! native analogue of the AO-bridge projectability the paper GaAs
     ! chain used.  Meaningful only when the mask is a proper subset.
     block
       real(dp) :: ptot
       if (count(pcol) == nprojP) write(*,'(a)') &
            '  WARNING: share rescale with full mask gives uniform trust'
       do k = 1, nk; do i = 1, nb
          ptot = 0.0_dp
          do j = 1, nprojP
             ptot = ptot + real(amnP(i,j,k),dp)**2 + aimag(amnP(i,j,k))**2
          end do
          if (ptot > 1.0e-14_dp) then
             P(i,k) = P(i,k) / ptot
          else
             P(i,k) = 0.0_dp
          end if
       end do; end do
     end block
  end if
  deallocate(amnP)
  if (nproj < nw) then
     write(*,'(a,i0,a,i0,a)') 'ERROR: num_wann = ', nw, ' exceeds the ', &
          nproj, ' projections in the amn (SCDM needs nproj >= nw)'
     error stop 'num_wann > nproj'
  end if
  allocate(wgt(nntot)); call mv_weights(bcart, nntot, wgt, wresid)
  write(*,'(a,es10.2)') '  MV b-weight residual = ', wresid

  ! ---- masks ----
  allocate(P_used(nb,nk), an_avg(nb), an_min(nb), an_bot(nb))
  call build_masks(P, eig, nb, nk, nw, e_fermi, rescale, regime, femode, femax_val, &
                   pool, npool, admit, frozen, regime_out, pool_factor_in=pool_factor, p_floor_in=p_floor, &
                   froz_demand_out=froz_dmd, fe_used_out=fe_used, &
                   froz_legacy_out=froz_legacy, floor_add_out=fl_add, &
                   floor_drop_out=fl_drop, floor_nk_out=fl_nk, &
                   P_used_out=P_used, anat_avg=an_avg, anat_min=an_min, anat_bot=an_bot, &
                   anat_plateau=an_plateau, anat_fid=an_fid, anat_pgate=an_pgate, &
                   anat_pfroz=an_pfroz, anat_ok=an_ok)
  write(*,'(a,i0,a)') '  frozen demand(worst k) = ', froz_dmd, &
       '   (states inside froz_emax at the worst k; capacity: nw >= demand)'
  if (froz_dmd > nw) write(*,'(a,i0,a,i0,a)') &
       '  WARNING: frozen capacity violation: demand ', froz_dmd, ' > nw ', &
       nw, ' (the frozen budget is spent on the nw deepest sub-ceiling states)'
  if (admit_all) then
     admit = .true.
     write(*,'(a)') '  admit-all: every pool state admitted (default; --p-floor <x> sets a projectability floor)'
  end if
  if (no_froz) then
     frozen = .false.
     write(*,'(a)') '  no-froz: frozen constraint disabled (pure Omega_I optimum)'
  end if
  write(*,'(a,a,a,i0,a,i0,a,i0,a,f6.2)') '  regime=', trim(regime_out), &
       '  pool[',pool(1),'..',pool(npool),'] npool=',npool, &
       '  frozen/k mean=', real(sum(count(frozen,1)),dp)/real(nk,dp)
  call froz_report_line(6, frozen, eig, pool, npool, nk, nw, e_fermi, fe_used, &
       femode, froz_dmd, no_froz, froz_legacy, fl_add, fl_drop, fl_nk)
  write(*,'(a,f8.4,a,i0,a,i0,a)') '  admit frac = ', &
       real(count(admit),dp)/real(npool*nk,dp), '   (', &
       npool*nk - count(admit), ' of ', npool*nk, ' pool states excluded)'
  call write_masks(trim(outp)//'.masks', pool, npool, frozen, admit, P_used, eig, e_fermi, &
                   nb, nk, nw, an_pfroz, an_pgate, fe_used)
  call write_anatomy(trim(outp)//'.anatomy', pool, npool, an_avg, an_min, an_bot, nb, nw, &
                     an_pgate, an_pfroz, an_plateau, an_fid, an_ok)
  write(*,'(a)') '  wrote '//trim(outp)//'.masks and '//trim(outp)//'.anatomy'

  ! ---- restrict to pool ----
  allocate(mmn_p(npool,npool,nntot,nk), amn_p(npool,nproj,nk), eig_p(npool,nk))
  do k = 1, nk
     do ib = 1, nntot
        do j = 1, npool; do i = 1, npool
           mmn_p(i,j,ib,k) = mmn(pool(i),pool(j),ib,k)
        end do; end do
     end do
     do i = 1, npool
        amn_p(i,:,k) = amn(pool(i),:,k)
        eig_p(i,k) = eig(pool(i),k)
     end do
  end do

  ! ---- disentangle -> rotate -> SCDM ----
  allocate(V(npool,nw,nk), W(npool,nw,nk), eig_t(nw,nk), mmn_rot(nw,nw,nntot,nk), &
           B(nproj,nw,nk), A_gauge(nw,nw,nk), anchors(nw))
  npocket = 0; pk_try = .false.
  pocket_loop: do
  call cpu_time(t0)
  call system_clock(wc0, wrate)
  ! Solver dispatch mirrors main.F: the branch lives here, not inside
  ! disentangle_frozen, because hybrid_grass_m uses hybrid_core_m for kinfo_t.
  solver_pick: block
    character(len=8) :: solv
    integer :: ndis_max, kq, nitg
    type(kinfo_t), allocatable :: kig(:)
    type(gcg_t) :: Pg
    complex(dp), allocatable :: Vg(:,:,:)
    real(dp) :: relg

    solv = trim(solver_sel)
    if (solv == 'auto') then
       ndis_max = 0
       do kq = 1, nk
          ndis_max = max(ndis_max, count(admit(:,kq) .and. .not. frozen(:,kq)))
       end do
       if (ndis_max >= 250) then
          solv = 'cg'
       else
          solv = 'gs'
       end if
       write(*,'(a,i0,a,a,a)') '  solver: auto -> max ndis = ', ndis_max, &
            ' so choosing ', trim(solv), ' (LCAO crossover ~250, unverified here)'
    end if
    if (solv == 'cg') then
       allocate(kig(nk))
       call build_kinfo(mmn_p, neigh, frozen, admit, npool, nntot, nk, nw, &
            kig, want_mdis=.false.)
       call gcg_build(mmn_p, neigh, wgt, kig, npool, nntot, nk, nw, Pg)
       allocate(Vg(Pg%dmax, Pg%amax, nk))
       ! same project-optimal SVD seed the sweep uses, so CG-vs-GS starts
       ! from one point
       call gcg_seed_svd(Pg, kig, amn_p, npool, nproj, Vg)
       ! ---- OPTIONAL RESTART (env W2W_CGLOAD=<file>) ----
       ! Replace the SVD seed with a Vg dumped by an earlier run
       ! (W2W_CGDUMP).  Vg is the REDUCED-domain array, which is why this
       ! lives here and not in disentangle_frozen.
       block
         character(len=256) :: cgl
         integer :: cu, d2, a2, k2, ios3
         call get_environment_variable('W2W_CGLOAD', cgl)
         if (len_trim(cgl) > 0) then
            open(newunit=cu, file=trim(cgl), form='unformatted', &
                 access='stream', status='old', iostat=ios3)
            if (ios3 /= 0) error stop 'W2W_CGLOAD: cannot open file'
            read(cu) d2, a2, k2
            ! checked, not trusted: seeding from a Vg of the wrong shape (or
            ! a different window) would run plausibly and mean nothing
            if (d2 /= Pg%dmax .or. a2 /= Pg%amax .or. k2 /= nk) then
               write(*,'(a,3i7)') '  W2W_CGLOAD file  dmax,amax,nk = ', d2,a2,k2
               write(*,'(a,3i7)') '  W2W_CGLOAD run   dmax,amax,nk = ', &
                    Pg%dmax, Pg%amax, nk
               error stop 'W2W_CGLOAD: dimension mismatch'
            end if
            read(cu) Vg
            close(cu)
            write(*,'(a,a)') '  cg: RESTART from ', trim(cgl)
         end if
       end block
       ! STARTUP GUARD (see main.F): a wrong gradient still descends.
       relg = gcg_check_gradient(Pg, Vg)
       write(*,'(a,es10.2)') '  cg: FD gradient check, rel err = ', relg
       if (relg > 1.0e-4_dp) error stop 'cg: analytic gradient disagrees &
            &with finite differences -- refusing to run'
       call grassmann_cg(Pg, Vg, dis_niter, 1.0e-9_dp, dis_rtol, omega, &
            nitg, dis_verb, whyg, gnorm=gng, gnorm0=gng0)
       block
         character(len=256) :: cgd
         integer :: cu2
         call get_environment_variable('W2W_CGDUMP', cgd)
         if (len_trim(cgd) > 0) then
            open(newunit=cu2, file=trim(cgd), form='unformatted', &
                 access='stream', status='replace')
            write(cu2) Pg%dmax, Pg%amax, nk
            write(cu2) Vg
            close(cu2)
            write(*,'(a,a)') '  cg: converged Vg dumped -> ', trim(cgd)
         end if
       end block
       call gcg_to_full(Pg, kig, Vg, npool, nw, V)
       ! cross-check the Grassmann objective against the independent
       ! full-frame evaluation; disagreement = not the same functional
       relg = omega_invariant(V, mmn_p, neigh, wgt, npool, nntot, nk, nw)
       write(*,'(a,f18.10,a,es10.2)') '  cg: independent Omega_I = ', relg, &
            '   diff ', abs(relg - omega)
       if (abs(relg - omega) > 1.0e-6_dp*max(abs(relg), 1.0_dp)) &
            error stop 'cg: Grassmann objective disagrees with omega_invariant'
       omega = relg
       ! Say WHY it stopped.  This line used to read "converged" for all five
       ! exit paths, which made every eval count in the cost tables
       ! unfalsifiable -- a stalled run, a capped run and a converged run were
       ! typographically identical.  Only GCG_OK is convergence.
       write(*,'(a,es10.3,a,es10.3)') '  cg: |grad| = ', gng, &
            '   |grad0| = ', gng0
       select case (whyg)
       case (GCG_OK)
          write(*,'(a,es10.3,a,es10.3,a,i0,a)') &
               '  cg: CONVERGED: |grad| = ', gng, &
               ' <= max(abs,rel*|grad0|) = ', &
               max(GCG_GTOL_ABS, GCG_GTOL_REL*gng0), ' after ', nitg, &
               ' gradient evaluations'
       case (GCG_SMALLDEC)
          write(*,'(a,i0,a)') '  cg: stopped on small decrement after ', nitg, &
               ' gradient evaluations (five consecutive small steps -- this '// &
               'is NOT a certified converged minimum)'
       case (GCG_STALL)
          write(*,'(a,i0,a)') '  cg: LINE SEARCH STALLED after ', nitg, &
               ' gradient evaluations -- NOT converged'
       case (GCG_TINYSTEP)
          write(*,'(a,i0,a)') '  cg: step underflow after ', nitg, &
               ' gradient evaluations -- NOT converged'
       case default
          write(*,'(a,i0,a,i0,a)') '  cg: HIT THE ITERATION CAP after ', nitg, &
               ' gradient evaluations (niter = ', dis_niter, ') -- NOT converged'
       end select
       if (whyg == GCG_OK) then
          write(dis_status,'(a,i0,a,f0.6,a)') 'converged (cg, ', nitg, &
               ' gradient evaluations, Omega_I ', omega, ' Ang^2)'
       else
          write(dis_status,'(a,a,a,i0,a,es8.2,a,f0.6,a)') 'NOT certified (cg exit ', &
               trim(merge('SMALLDEC', merge('STALL   ', merge('TINYSTEP', 'CAP     ', &
               whyg == GCG_TINYSTEP), whyg == GCG_STALL), whyg == GCG_SMALLDEC)), ', ', nitg, &
               ' gradient evaluations, |grad| ', gng, ', Omega_I ', omega, ' Ang^2)'
       end if
       call gcg_free(Pg)
       deallocate(Vg, kig)
    else
       block
         integer :: nsw_gs
         logical :: conv_gs
         character(len=320) :: capv_gs
         call disentangle_frozen(mmn_p, amn_p, neigh, wgt, npool, nproj, nntot, nk, nw, &
              admit, frozen, V, omega, dis_niter, 1.0e-9_dp, dis_rtol, dis_check, dis_verb, &
              colored=dis_colored, adaptive=dis_adaptive, mixing=dis_mix, &
              anderson=dis_aa, init_transport=init_transp, &
              nsweeps=nsw_gs, is_converged=conv_gs, cap_verdict=capv_gs)
         if (conv_gs) then
            write(dis_status,'(a,i0,a,f0.6,a)') 'converged (gs, ', nsw_gs, ' sweeps, Omega_I ', &
                 omega, ' Ang^2)'
         else
            write(dis_status,'(a,i0,a,f0.6,a)') 'stopped at the sweep cap (', nsw_gs, &
                 ' sweeps, Omega_I ', omega, ' Ang^2): '
            dis_status = trim(dis_status)//' '//trim(capv_gs)
         end if
       end block
    end if
  end block solver_pick
  call system_clock(wc1)
  call cpu_time(t1)
  write(*,'(a,f16.8,a,f8.2,a,f8.2,a)') '  Omega_I = ', omega, '   (disentangle ', &
       real(wc1-wc0,dp)/real(wrate,dp), ' s wall, ', t1-t0, ' s cpu)'
  call rotate_to_subspace(eig_p, V, mmn_p, neigh, npool, nntot, nk, nw, eig_t, W, mmn_rot)
  ! ---- subspace continuity ----
  call subspace_continuity(mmn_rot, neigh, nw, nntot, nk, output_unit, cont_smin, cont_k, cont_kb)
  ! ---- pocket skip (walked ceilings only; hybrid_prep_m::pocket_ceiling) ----
  if (pk_try) then
     if (cont_smin < POCKET_SV) then
        write(*,'(a,es10.3,a,f8.3,a)') '  pocket skip: the weakest link stays below 0.5 (', &
             cont_smin, '); the ceiling ', fe_s, ' eV and its disentanglement are kept'
        V = V_s; W = W_s; eig_t = eig_t_s; mmn_rot = mmn_rot_s; omega = omega_s
        frozen = frozen_s; admit = admit_s; dis_status = dis_status_s; fe_used = fe_s
        cont_smin = cont_smin_s; cont_k = cont_k_s; cont_kb = cont_kb_s
        call write_masks(trim(outp)//'.masks', pool, npool, frozen, admit, P_used, eig, e_fermi, &
                         nb, nk, nw, an_pfroz, an_pgate, fe_used)
        exit pocket_loop
     end if
     npocket = npocket + 1
     write(*,'(a,f8.3,a,i0,a)') '  pocket skip: accepted, ceiling ', fe_used, ' eV (step ', npocket, ')'
     call write_masks(trim(outp)//'.masks', pool, npool, frozen, admit, P_used, eig, e_fermi, &
                      nb, nk, nw, an_pfroz, an_pgate, fe_used)
  end if
  if (.not. pocket_on .or. trim(femode) /= 'auto' .or. no_froz .or. &
      npocket >= POCKET_MAX_STEPS) exit pocket_loop
  if (.not. allocated(sv_links)) allocate(sv_links(nntot,nk), eps_pool(npool,nk), P_pool(npool,nk))
  call continuity_links(mmn_rot, nw, nntot, nk, sv_links)
  do k = 1, nk
     do i = 1, npool
        eps_pool(i,k) = eig_p(i,k) - e_fermi
        P_pool(i,k) = P_used(pool(i),k)
     end do
  end do
  call pocket_ceiling(eps_pool, P_pool, npool, nk, fe_used, sv_links, neigh, nntot, &
                      fe_try, pk_fired, pk_msg)
  if (.not. pk_fired) exit pocket_loop
  write(*,'(a,es10.3,a)') '  pocket skip: the weakest link ', cont_smin, ' touches'//trim(pk_msg)
  write(*,'(a,f8.3,a,f8.3,a)') '  pocket skip: ceiling ', fe_used, ' -> ', fe_try, &
       ' eV, disentangling again'
  if (.not. allocated(V_s)) allocate(V_s(npool,nw,nk), W_s(npool,nw,nk), eig_t_s(nw,nk), &
       mmn_rot_s(nw,nw,nntot,nk), admit_s(npool,nk), frozen_s(npool,nk), pool_s(npool))
  V_s = V; W_s = W; eig_t_s = eig_t; mmn_rot_s = mmn_rot; omega_s = omega
  admit_s = admit; frozen_s = frozen; pool_s = pool; dis_status_s = dis_status; fe_s = fe_used
  cont_smin_s = cont_smin; cont_k_s = cont_k; cont_kb_s = cont_kb
  write(femode_try,'(f10.4)') fe_try
  call build_masks(P, eig, nb, nk, nw, e_fermi, rescale, regime, trim(adjustl(femode_try)), fe_try, &
                   pool, npool, admit, frozen, regime_out, pool_factor_in=pool_factor, p_floor_in=p_floor, &
                   froz_demand_out=froz_dmd, fe_used_out=fe_used, P_used_out=P_used)
  if (npool /= size(pool_s) .or. any(pool /= pool_s)) then
     ! the pool must not move: every array below is in the pool basis
     write(*,'(a)') '  pocket skip: the pool changed at the lowered ceiling; the walked ceiling is kept'
     deallocate(pool, admit, frozen)
     allocate(pool(size(pool_s)), admit(size(pool_s),nk), frozen(size(pool_s),nk))
     npool = size(pool_s); pool = pool_s; admit = admit_s; frozen = frozen_s; fe_used = fe_s
     exit pocket_loop
  end if
  if (admit_all) admit = .true.
  pk_try = .true.
  end do pocket_loop
  write(cont_line,'(es9.3,a,i0,a,i0)') cont_smin, &
       ' (smallest singular value of the hand-off overlaps), at k ', cont_k, ' -> ', cont_kb
  if (cont_smin < SUBSPACE_JUMP_SV) cont_line = trim(cont_line)// &
       '; WARNING: below 1e-2, the subspace jumps there'
  do k = 1, nk
     call zgemm('C','N', nproj, nw, npool, CONE, amn_p(1,1,k), npool, W(1,1,k), npool, &
                CZERO, B(1,1,k), nproj)
  end do
  ! ---- gauge selection (shared dispatcher; env W2W_SEEDK/W2W_SEEDLOSS/
  !      W2W_SEEDGATE are honoured inside choose_gauge) ----
  ! ---- the b-vector of every (k,b) block, from the same .nnkp: a .nnkp lists
  ! each k's neighbours in its own order (511 of 512 k on MgB2), so bcart(:,b)
  ! -- k = 1's -- is the vector of block (k,b) only at k = 1.  The gauge
  ! policy's spreads (seed tournament, narration) need the per-block vectors.
  allocate(bkb(3, nntot, nk))
  block
    character(len=256) :: bvmsg
    logical :: bvok
    call nnkp_block_bvectors(trim(seed)//'.nnkp', nk, nntot, neigh, gshift, bcart, bkb, &
                             bvok, bvmsg)
    if (.not. bvok) then
       write(*,'(a)') 'ERROR: per-block b-vectors: '//trim(bvmsg)
       error stop 'nnkp_block_bvectors refused'
    end if
  end block
  call choose_gauge(B, eig_t, mmn_rot, neigh, wgt, bcart, nproj, nw, &
                    nntot, nk, gauge_mode, trim(outp), A_gauge, anchors, &
                    best_svmin, bkb=bkb)
  ! ---- hand-off; acceptance runs INSIDE write_handoff, before anything
  ! is written: full coverage over every (k,b) block, error stop on
  ! violation ----
  call write_handoff(trim(outp), trim(seed), eig_t, A_gauge, mmn_rot, neigh, gshift, &
                     e_fermi, nw, nk, nntot)
  write(*,'(a)') '  wrote '//trim(outp)//'.{amn,mmn,eig}'
  call write_dis_status(trim(outp)//'.dis', dis_status, cont_line)
  write(*,'(a)') '  disentanglement: '//trim(dis_status)
  if (window_route) call window_route_emit()

contains

  !> The automatic window route: the frozen states' energy hull becomes the
  !> frozen window (lowered where a k-point would hold more than nw frozen
  !> states, which Wannier90 refuses), the outer window keeps the pool bottom
  !> and takes the ceiling that minimises the frozen-constrained Omega_I
  !> among a few candidates above the frozen ceiling. Stock Wannier90 then
  !> disentangles from the seed overlaps and the anchor projections.
  subroutine window_route_emit()
    real(dp) :: fmin, fmax, dlo, dhi, cap_k, cand(5), om_c(5), om_best, d_best, ek
    real(dp), allocatable :: esrt(:)
    logical, allocatable :: frz_w(:,:), adm_w(:,:)
    complex(dp), allocatable :: Vw(:,:,:)
    integer :: i, k, nc, ic, nfz, nad, j
    logical :: ok
    real(dp), parameter :: offs(4) = [8.0_dp, 5.0_dp, 3.0_dp, 1.5_dp]
    fmin = huge(1.0_dp); fmax = -huge(1.0_dp); dlo = huge(1.0_dp); dhi = -huge(1.0_dp)
    do k = 1, nk
       do i = 1, npool
          ek = eig_p(i,k) - e_fermi
          dlo = min(dlo, ek); dhi = max(dhi, ek)
          if (frozen(i,k)) then
             fmin = min(fmin, ek); fmax = max(fmax, ek)
          end if
       end do
    end do
    if (fmax < fmin) then
       write(*,'(a)') '  window route: no frozen states, nothing to express as a window'
       return
    end if
    ! Wannier90 needs at most nw frozen states at every k: cap the frozen
    ! ceiling just below the (nw+1)-th pool energy above fmin at the worst k
    allocate(esrt(npool))
    cap_k = huge(1.0_dp)
    do k = 1, nk
       esrt = eig_p(:,k) - e_fermi
       call sort_asc(esrt, npool)
       j = 0
       do i = 1, npool
          if (esrt(i) >= fmin) then
             j = j + 1
             if (j == nw + 1) then
                cap_k = min(cap_k, esrt(i))
                exit
             end if
          end if
       end do
    end do
    if (fmax >= cap_k) then
       write(*,'(a,f9.3,a,f9.3,a)') '  window route: frozen ceiling ', fmax, &
            ' eV lowered to ', cap_k - 1.0e-4_dp, ' eV so no k-point holds more than num_wann frozen states'
       fmax = cap_k - 1.0e-4_dp
    end if
    ! candidate outer ceilings: the pool top, and fixed offsets above the frozen ceiling
    nc = 1; cand(1) = dhi
    do i = 1, size(offs)
       if (fmax + offs(i) < dhi - 1.0e-6_dp) then
          nc = nc + 1; cand(nc) = fmax + offs(i)
       end if
    end do
    allocate(frz_w(npool,nk), adm_w(npool,nk), Vw(npool,nw,nk))
    om_best = huge(1.0_dp); d_best = dhi; om_c = huge(1.0_dp)
    do ic = 1, nc
       ok = .true.
       do k = 1, nk
          do i = 1, npool
             ek = eig_p(i,k) - e_fermi
             frz_w(i,k) = (ek >= fmin) .and. (ek <= fmax)
             adm_w(i,k) = (ek >= dlo) .and. (ek <= cand(ic))
          end do
          nfz = count(frz_w(:,k)); nad = count(adm_w(:,k))
          if (nfz > nw .or. nad < nw) ok = .false.
       end do
       if (.not. ok) then
          write(*,'(a,f9.3,a)') '  window route: dis_win_max = ', cand(ic), &
               ' eV rejected (a k-point holds more frozen than num_wann or fewer than num_wann admitted states)'
          cycle
       end if
       call disentangle_frozen(mmn_p, amn_p, neigh, wgt, npool, nproj, nntot, nk, nw, &
            adm_w, frz_w, Vw, om_c(ic), min(dis_niter, 2000), 1.0e-9_dp, &
            max(dis_rtol, 1.0e-6_dp), dis_check, .false., &
            colored=dis_colored, adaptive=dis_adaptive, mixing=dis_mix, &
            anderson=dis_aa, init_transport=.false.)
       write(*,'(a,f9.3,a,f16.8)') '  window route: dis_win_max = ', cand(ic), &
            ' eV  Omega_I = ', om_c(ic)
       if (om_c(ic) < om_best) then
          om_best = om_c(ic); d_best = cand(ic)
       end if
    end do
    if (om_best >= huge(1.0_dp)) then
       write(*,'(a)') '  window route: no admissible outer window; nothing written'
       return
    end if
    write(*,'(a,f9.3,a,f9.3,a,f9.3,a,f9.3,a,f14.6)') '  window route: dis_froz = [', &
         fmin, ',', fmax, ']  dis_win = [', dlo, ',', d_best, ']  Omega_I = ', om_best
    call write_window_route(trim(outp), eig, nb, nk, nw, e_fermi, amn, nproj, anchors, &
                            mmn, neigh, gshift, nntot, &
                            fmin, fmax, dlo, d_best, nc, cand(1:nc), om_c(1:nc))
    write(*,'(a)') '  window route: wrote '//trim(outp)//'-win.{win,amn,mmn,eig}; '// &
         'run wannier90 on '//trim(outp)//'-win with those windows'
  end subroutine window_route_emit

  subroutine sort_asc(a, n)
    integer, intent(in) :: n
    real(dp), intent(inout) :: a(n)
    integer :: i, j
    real(dp) :: t
    do i = 2, n
       t = a(i); j = i - 1
       do while (j >= 1)
          if (a(j) <= t) exit
          a(j+1) = a(j); j = j - 1
       end do
       a(j+1) = t
    end do
  end subroutine sort_asc

  !> Named-option command line.  A plain loop over `--flag [value]` pairs;
  !> the first one or two BARE tokens are accepted as seed and nw so the
  !> historical `w2w_hybrid <case> <N>` muscle memory keeps working.
  subroutine parse_cli()
    integer :: ia, na, ios
    character(len=256) :: a, val
    logical :: have_seed, have_nw
    have_seed = .false.; have_nw = .false.
    na = command_argument_count()
    if (na == 0) then
       call print_help()
       stop
    end if
    ia = 1
    do while (ia <= na)
       call get_command_argument(ia, a)
       select case (trim(a))
       case ('--help', '-h')
          call print_help()
          stop
       case ('--seed')
          call optval(ia, na, a, val); seed = val; have_seed = .true.
       case ('--nw', '--num-wann')
          call optval(ia, na, a, val)
          read(val, *, iostat=ios) nw
          if (ios /= 0 .or. nw < 1) error stop '--nw: want a positive integer'
          have_nw = .true.
       case ('--out')
          call optval(ia, na, a, val); outp = val
       case ('--rescale')
          call optval(ia, na, a, val)
          select case (trim(val))
          case ('occ99', 'p99', 'max', 'none', 'lowdin', 'share')
             rescale = trim(val)
          case default
             error stop '--rescale: occ99|p99|max|none|lowdin|share'
          end select
       case ('--regime')
          call optval(ia, na, a, val)
          select case (trim(val))
          case ('auto', 'per-state', 'band-level')
             regime = trim(val)
          case default
             error stop '--regime: auto|per-state|band-level'
          end select
       case ('--froz-emax')
          call optval(ia, na, a, val)
          femode = val
          if (trim(val) /= 'auto' .and. trim(val) /= 'none') then
             read(val, *, iostat=ios) femax_val
             if (ios /= 0) error stop '--froz-emax: auto|none|<eV>'
          end if
       case ('--efermi', '--e-fermi')
          call optval(ia, na, a, val)
          read(val, *, iostat=ios) e_fermi
          if (ios /= 0) error stop '--efermi: want a real (eV)'
       case ('--sweeps')
          call optval(ia, na, a, val)
          read(val, *, iostat=ios) dis_niter
          if (ios /= 0 .or. dis_niter < 1) error stop '--sweeps: want a positive integer'
       case ('--rtol')
          call optval(ia, na, a, val)
          read(val, *, iostat=ios) dis_rtol
          if (ios /= 0 .or. dis_rtol <= 0.0_dp) error stop '--rtol: want a positive real'
       case ('--check-every')
          call optval(ia, na, a, val)
          read(val, *, iostat=ios) dis_check
          if (ios /= 0 .or. dis_check < 1) error stop '--check-every: want a positive integer'
       case ('--colored')
          dis_colored = .true.
       case ('--adaptive')
          dis_adaptive = .true.
       case ('--relax')
          call optval(ia, na, a, val)
          read(val, *, iostat=ios) dis_mix
          if (ios /= 0) error stop '--relax: want a real (negative = scheduled)'
       case ('--gauge')
          call optval(ia, na, a, val); gauge_mode = val
       case ('--anderson')
          call optval(ia, na, a, val)
          read(val, *, iostat=ios) dis_aa
          if (ios /= 0 .or. dis_aa < 0) error stop '--anderson: want an integer >= 0'
       case ('--np-mask')
          call optval(ia, na, a, val)
          read(val, *, iostat=ios) np_mask
          if (ios /= 0) error stop '--np-mask: want an integer'
       case ('--init-transport')
          init_transp = .true.
       case ('--admit-all')
          admit_all = .true.
       case ('--p-floor')
          call optval(ia, na, a, val)
          read(val, *, iostat=ios) p_floor
          if (ios /= 0 .or. p_floor < 0.0_dp .or. p_floor > 1.0_dp) error stop '--p-floor: want a real in [0, 1]'
          pfloor_set = .true.
          write(*,'(a,es10.3)') '  projectability floor of the admission mask = ', p_floor
       case ('--pocket-skip')
          call optval(ia, na, a, val)
          select case (trim(val))
          case ('on')
             pocket_on = .true.
          case ('off')
             pocket_on = .false.
          case default
             error stop '--pocket-skip: on|off'
          end select
       case ('--window')
          window_route = .true.
       case ('--continuity-chk')
          call optval(ia, na, a, val); cont_chk_seed = val
       case ('--no-froz')
          no_froz = .true.
       case ('--all-free')
          admit_all = .true.; no_froz = .true.
       case ('--pool-factor')
          call optval(ia, na, a, val)
          read(val, *, iostat=ios) pool_factor
          if (ios /= 0 .or. pool_factor <= 0.0_dp) error stop '--pool-factor: want a positive real'
          write(*,'(a,f6.3)') '  pool factor = ', pool_factor
       case ('--solver')
          call optval(ia, na, a, val)
          select case (trim(val))
          case ('gs', 'cg', 'auto'); solver_sel = trim(val)
          case default; error stop '--solver: gs|cg|auto'
          end select
       case default
          if (a(1:2) == '--') then
             write(*,'(a)') 'ERROR: unknown option '''//trim(a)// &
                  ''' (see w2w_hybrid --help)'
             error stop 'unknown option'
          else if (.not. have_seed) then
             seed = a; have_seed = .true.
          else if (.not. have_nw) then
             read(a, *, iostat=ios) nw
             if (ios /= 0 .or. nw < 1) &
                  error stop 'second bare argument (num_wann): want a positive integer'
             have_nw = .true.
          else
             write(*,'(a)') 'ERROR: unexpected argument '''//trim(a)// &
                  ''' -- all further options must be named (see w2w_hybrid --help)'
             error stop 'unexpected positional argument'
          end if
       end select
       ia = ia + 1
    end do
    ! admission default: every pool state (as the HYBR path, main.F, and the paper's LAPW protocol);
    ! --p-floor <x> sets a projectability floor, an explicit --admit-all overrides it
    if (.not. pfloor_set) admit_all = .true.
    if (len_trim(cont_chk_seed) > 0) return     ! --continuity-chk needs neither --seed nor --nw
    if (.not. have_seed .or. .not. have_nw) then
       write(*,'(a)') 'usage: w2w_hybrid --seed <case> --nw <N> [options]   (or: w2w_hybrid <case> <N> [options])'
       write(*,'(a)') 'run `w2w_hybrid --help` for the full reference'
       error stop 'missing required arguments (--seed and --nw)'
    end if
  end subroutine parse_cli

  !> Window route: subspace continuity of the subspace Wannier90 disentangled itself,
  !> read from <cseed>.chk and <cseed>.mmn after localization.  The same measure and
  !> thresholds as the hand-off's (subspace_continuity): the smallest singular value of
  !> V(k)^H M(k,b) V(k+b) over all links, links below 0.5 counted, a warning below
  !> SUBSPACE_JUMP_SV.  On the LCAO SnTe window route this is 0.013-0.031 at X against
  !> 0.20 for the hand-off.  A diagnostic only; nothing is written.
  subroutine window_chk_continuity(cseed)
    character(*), intent(in) :: cseed
    integer :: cnb, cnk, cnw, mnb, mnk, mnn, k, b, kmn, kbmn, nhalf
    logical :: cdis
    real(dp) :: om_i, smn
    real(dp), allocatable :: ckp(:,:), sv(:,:)
    complex(dp), allocatable :: cV(:,:,:), cM(:,:,:,:), crot(:,:,:,:), tmp(:,:)
    integer, allocatable :: cng(:,:), cgs(:,:,:)
    call read_w90_chk_subspace(cseed//'.chk', cnb, cnk, cnw, cdis, cV, om_i, ckp)
    if (.not. cdis) then
       write(*,'(a)') '  subspace continuity: '//cseed//'.chk holds no disentangled subspace '// &
            '(num_bands = num_wann): nothing to check'
       return
    end if
    call read_mmn(cseed//'.mmn', cM, cng, cgs, mnb, mnk, mnn)
    if (mnb /= cnb .or. mnk /= cnk) error stop '--continuity-chk: .mmn and .chk dimensions differ'
    allocate(crot(cnw,cnw,mnn,cnk), tmp(cnb,cnw), sv(mnn,cnk))
    do k = 1, cnk
       do b = 1, mnn
          tmp = matmul(cM(:,:,b,k), cV(:,:,cng(b,k)))
          crot(:,:,b,k) = matmul(conjg(transpose(cV(:,:,k))), tmp)
       end do
    end do
    call continuity_links(crot, cnw, mnn, cnk, sv)
    smn = minval(sv); nhalf = count(sv < 0.5_dp)
    kmn = 0; kbmn = 0
    do k = 1, cnk
       do b = 1, mnn
          if (sv(b,k) == smn .and. kmn == 0) then
             kmn = k; kbmn = cng(b,k)
          end if
       end do
    end do
    write(*,'(a,f7.4,a,3f7.3,a,3f7.3,a,i0,a,i0,a,es8.1,a,f12.4,a)') '  subspace continuity: smallest '// &
         'singular value of Wannier90''s disentangled overlaps ', smn, ' at k (', ckp(:,kmn), ') -> (', &
         ckp(:,kbmn), '); ', nhalf, ' of ', cnk*mnn, ' (k,b) blocks below 0.5 (warning below ', &
         SUBSPACE_JUMP_SV, '; Omega_I ', om_i, ' Ang^2)'
    if (smn < SUBSPACE_JUMP_SV) then
       write(*,'(a,es8.1,a)') '  *** WARNING: Wannier90''s disentangled subspace JUMPS between neighbouring '// &
            'k-points (smallest singular value ', smn, ').'
       write(*,'(a)') '  *** No gauge repairs a jump: expect poor interpolation. The hand-off '// &
            '(--method hybrid) keeps its own subspace continuous.'
    else if (nhalf > 0) then
       write(*,'(a,i0,a)') '  note: ', nhalf, ' links below 0.5: the subspace changes sharply there; '// &
            'compare the band distance with the hand-off''s.'
    end if
  end subroutine window_chk_continuity

  !> fetch the value token of `--flag value`; errors if it is missing
  subroutine optval(ia, na, flag, val)
    integer, intent(inout) :: ia
    integer, intent(in)    :: na
    character(len=*), intent(in)  :: flag
    character(len=*), intent(out) :: val
    if (ia + 1 > na) then
       write(*,'(a)') 'ERROR: option '//trim(flag)//' needs a value'
       error stop 'missing option value'
    end if
    ia = ia + 1
    call get_command_argument(ia, val)
  end subroutine optval

  subroutine print_help()
    character(len=*), parameter :: nl = char(10)
    write(*,'(a)') &
'w2w_hybrid -- standalone hybrid pre-Wannierization (wien2wannier release v3.0)'//nl//&
''//nl//&
'Builds per-state projectability from the seed .amn, selects a band pool and'//nl//&
'frozen manifold, minimizes the gauge-invariant spread Omega_I over the pool'//nl//&
'(frozen-constrained disentanglement), rotates to the projected-Hamiltonian'//nl//&
'eigenbasis and fixes an SCDM/Loewdin gauge, then writes an isolated Nb=Nw'//nl//&
'hand-off that stock Wannier90 can localize without dis_* keywords.'//nl//&
''//nl//&
'USAGE'//nl//&
'  w2w_hybrid --seed <case> --nw <N> [options]'//nl//&
'  w2w_hybrid <case> <N> [options]        (bare seed/nw shorthand)'//nl//&
''//nl//&
'CASE FILES (read from the current directory)'//nl//&
'  <seed>.amn   projections A_mn(k)  (text; a fresh binary sidecar <seed>.amn.bin'//nl//&
'               is used/written automatically -- delete the .bin to force a re-parse)'//nl//&
'  <seed>.mmn   overlaps M_mn(k,b)   (text; same .bin sidecar caching)'//nl//&
'  <seed>.eig   band energies (eV); --efermi is subtracted from them'//nl//&
'  <seed>.nnkp  k-mesh + b-vector neighbour list (wannier90 -pp output)'//nl//&
'  projmap.dat  optional: AUTO-projection channel map; enables the --np-mask'//nl//&
'               channel selection and shell-complete anchor grouping'//nl//&
'WRITES'//nl//&
'  <out>.amn/.mmn/.eig   the Nb=Nw hand-off (default <out> = <seed>-hybF)'//nl//&
'  <out>.anchors         the SCDM anchor columns + a gauge-policy stamp'//nl//&
'  <out>.masks           per-state masks: k, band, eps-E_F, P, frozen, admitted'//nl//&
'                        (tools/plot_certificate.py draws the window certificate)'//nl//&
'  <out>.anatomy         auto-trust anatomy: per pool band Pbar, minP, bottom,'//nl//&
'                        relevance and class, plus the derived thresholds'//nl//&
'                        (tools/plot_anatomy.py draws it)'//nl//&
'  <out>.dis             the disentanglement outcome (converged, or stopped at the'//nl//&
'                        cap with a graded tail, or the CG exit) and the subspace'//nl//&
'                        continuity, as comment lines for the hand-off .win'//nl//&
''//nl//&
'REQUIRED'//nl//&
'  --seed <case>        input basename (see CASE FILES)'//nl//&
'  --nw <N>             number of Wannier functions (target subspace dimension)'//nl//&
''//nl//&
'OUTPUT AND MASKS'//nl//&
'  --out <name>         hand-off basename                       [<seed>-hybF]'//nl//&
'  --rescale <mode>     projectability rescaling: occ99|p99|max|none|lowdin|share'//nl//&
'                       (lowdin = Loewdin-orthonormalised P)    [occ99]'//nl//&
'  --regime <mode>      mask regime: auto|per-state|band-level  [auto]'//nl//&
'  --froz-emax <x>      frozen-window ceiling above E_F: auto|none|<eV>  [auto]'//nl//&
'  --pocket-skip <on|off> with --froz-emax auto: when the weakest hand-off link'//nl//&
'                       (< 0.5) touches a mixed-character band (mean P >= 0.5)'//nl//&
'                       frozen at fewer than 25% of the k-points, lower the'//nl//&
'                       ceiling 0.1 eV under it and disentangle again; kept'//nl//&
'                       only if the new weakest link is >= 0.5      [on]'//nl//&
'  --efermi <eV>        Fermi energy subtracted from the .eig energies  [0.0]'//nl//&
'  --np-mask <N>        amn columns used for P/masks: 0 = all; N > 0 = leading N'//nl//&
'                       channels (mirrored over SOC sigma blocks); N < 0 = all'//nl//&
'                       channels with l <= |N| (needs projmap.dat)   [0]'//nl//&
'  --admit-all          admit every pool state (the default); overrides --p-floor'//nl//&
'  --p-floor <x>        admission mask: a pool state outside the frozen set'//nl//&
'                       enters only if P >= x               [none: admit all]'//nl//&
'  --continuity-chk <seed>  report the subspace continuity of a finished Wannier90'//nl//&
'                       window-route run from <seed>.chk and <seed>.mmn, then exit'//nl//&
'  --window             also express the selection as Wannier90 disentanglement'//nl//&
'                       windows (the automatic window route): frozen window ='//nl//&
'                       the frozen states'' energy hull, outer ceiling chosen by'//nl//&
'                       Omega_I, written to <out>-win.{win,amn,mmn,eig}'//nl//&
'  --no-froz            disable the frozen constraint (pure Omega_I optimum)'//nl//&
'  --all-free           both --admit-all and --no-froz'//nl//&
'  --pool-factor <x>    disentanglement pool floor, npool >= max(nw+1, ceil(x*nw)) [1.5]'//nl//&
''//nl//&
'DISENTANGLEMENT'//nl//&
'  --solver <s>         cg = Grassmann CG, stops on the projected gradient'//nl//&
'                       (parallel over k) | gs = Gauss-Seidel sweep | auto'//nl//&
'                       (picks by max ndis; crossover ~250)       [cg]'//nl//&
'  --sweeps <N>         sweep / iteration cap; a run that stops at the cap is'//nl//&
'                       graded (flat or still descending) in the log and in'//nl//&
'                       <out>.dis                                 [20000]'//nl//&
'  --rtol <x>           relative convergence tolerance            [1e-10]'//nl//&
'  --check-every <N>    sweeps between exact convergence checks   [10]'//nl//&
'  --relax <x>          GS mixing (SOR) factor: 1 = off, >1 fixed over-'//nl//&
'                       relaxation, negative = scheduled (start at |x|,'//nl//&
'                       taper toward 1; -1.9 is the production recipe)  [1.0]'//nl//&
'  --colored            red-black colored GS sweeps (OpenMP parallel; NOTE:'//nl//&
'                       reaches a different fixed point than the serial order)'//nl//&
'  --adaptive           adaptive convergence-check throttling'//nl//&
'  --anderson <N>       Anderson-acceleration history depth (0 = off)  [0]'//nl//&
'  --init-transport     seed the disentangler by BFS parallel transport'//nl//&
'                       instead of the per-k SVD seed'//nl//&
''//nl//&
'GAUGE'//nl//&
'  --gauge <mode>       hand-off gauge construction               [scdm]'//nl//&
'                         scdm    per-k SCDM on the projected amn columns'//nl//&
'                         sms<N>  scdm + N smoothing sweeps'//nl//&
'                         ptg     SCDM anchors at one seed k, then Procrustes/'//nl//&
'                                 parallel-transport propagation over the BZ'//nl//&
'                         pts<N>  ptg + N smoothing sweeps'//nl//&
'                         auto<N> conditioning-gated ptg/scdm dispatch + N'//nl//&
'                                 smoothing sweeps (auto8 = production recipe)'//nl//&
'                         lowdin  plain Loewdin projection gauge; lds<N> + smooth'//nl//&
'                         bscdm   block-SCDM (per-isolated-manifold); bsms<N>'//nl//&
'                         slwa/slwd/slwf/slwp/tfl  Loewdin-anchor variants'//nl//&
'                         u0      load <out>.gauge0 written by W2W_GAUGE0OUT=1'//nl//&
''//nl//&
'MISC'//nl//&
'  --help, -h           this reference'//nl//&
''//nl//&
'ENVIRONMENT VARIABLES (all optional)'//nl//&
'  W2W_VERBOSE=1        per-sweep Omega_I trace (hour-scale runs are otherwise'//nl//&
'                       silent until they finish)'//nl//&
'  W2W_PAMN=<file>      separate .amn used as the projectability/mask source'//nl//&
'                       while the seed .amn still feeds the gauge and hand-off'//nl//&
'  W2W_SEEDK=<sel>      transport-gauge seed k: best|median|<int>|opt[<N>]'//nl//&
'                       (opt<N> = gate-then-sample tournament)    [1]'//nl//&
'  W2W_SEEDLOSS=<loss>  loss the seed tournament ranks by (smooth|absm|...)'//nl//&
'                       [smooth]'//nl//&
'  W2W_SEEDGATE=<x>     per-k sigma_min stability gate for seed candidates'//nl//&
'                       [2e-3]'//nl//&
'  W2W_SIGGATE=<x>      sigma_min conditioning gate for the auto/tfl/slw*'//nl//&
'                       gauge dispatch                            [2e-3]'//nl//&
'  W2W_SMTOL=<x>        converged-smoothing exit tolerance (stop smoothing'//nl//&
'                       once the largest per-sweep change drops below it;'//nl//&
'                       unset = run the full sweep count)'//nl//&
'  W2W_PROBEN=<N>       probe sweeps per candidate in the late-alpha ladder'//nl//&
'                       (scheduled relaxation only)               [50]'//nl//&
'  W2W_STRUCT=<file>    case.struct: shell-complete SCDM anchor grouping keyed'//nl//&
'                       off the real structure''s site orbits (unset = site-'//nl//&
'                       local grouping)'//nl//&
'  W2W_ANCHORS=<mode>   shellcomplete|pshell: constrain SCDM pivots to whole'//nl//&
'                       (atom,l,rf) shells x both spin blocks (needs projmap.dat)'//nl//&
'  W2W_ANCHORS_RELAX=1  disable the shell-closure/conditioning refusals on a'//nl//&
'                       collapsed anchor frame (forensics; the run is stamped)'//nl//&
'  W2W_CGDUMP=<file>    dump the converged CG iterate for a later restart'//nl//&
'  W2W_CGLOAD=<file>    restart the CG solver from a W2W_CGDUMP file'//nl//&
'  W2W_GAUGE0OUT=1      write the final gauge to <out>.gauge0 (the input of'//nl//&
'                       --gauge u0; exact-gauge reruns)'//nl//&
''//nl//&
'  (The w2w/w2wc binaries additionally honour W2W_TAILS, W2W_AUTO_DEGTOL'//nl//&
'   and W2W_EDGES; those do not affect this tool.)'
  end subroutine print_help

end program w2w_hybrid
