!!! wien2wannier/SRC_w2w/main.F
!!!
!!!    Main program 'w2w'
!!!
!!! Copyright 2010-2012 Jan Kuneš, Philipp Wissgott
!!!           2013-2016 Elias Assmann

!!/--- Files expected in 'def' ---
!!  5 inwf	'old'	  'formatted'
!!  6 outputwf	'unknown' 'formatted'
!!  7 amn	'unknown' 'formatted'
!!  8 mmn	'unknown' 'formatted'
!!  9 vector	'unknown' 'unformatted'
!! 10 nnkp	'old'	  'formatted'
!! 12 eig	'unknown' 'formatted'
!! 18 vsp	'old'	  'formatted'
!! 20 struct	'old'	  'formatted'
!! 50 energy	'old'     'formatted'
!! 51 fermi	'old'     'formatted'
!!\---

program wf
  use param,     only: unit_def, unit_vector, unit_in, unit_out, unit_struct, &
       &               unit_vsp, wien2wannier_version, Lmax2, Nrad, Nrf
  use w2w,       only: Nmat, unit_nnkp, unit_amn, unit_mmn, unit_ene, &
       unit_inso, soc_hkl, soc_hkl_set, &
       &               unit_fermi, unit_eig
  use const,     only: R8, C16, BUFSZ, Ryd_eV
  use xa,        only: init_xa
  use xa3,       only: init_xa3, A, XK, YK, ZK, GX, GY, GZ, vecsz
  use bessel,    only: init_bessel
  use Amn_Mmn,   only: init_Amn_Mmn, overlap
  use pairs,     only: kp, kpb, bqx,bqy,bqz, bqx1,bqy1,bqz1, init_pairs
  use util,      only: paropen, ptime
  use wien2k,    only: errflg, errclr, gtfnam, vector_kind_probe
  use clio,      only: croak, close_checked
  use structmod, only: struct_t, struct_read
  use inwfmod,   only: inwf_t, inwf_read
  use gaunt_cache, only: init_gaunt_cache, cleanup_gaunt_cache
  use loabc,       only: cleanup_loabc
  use mpi_w2w

  !! in-memory MODE=HYBR Stage-2 (validated hybrid pre-Wannierization modules)
  use hybrid_prep_m,  only: mv_weights, build_masks, lowdin_pstate, &
       &                    froz_report_line, pocket_ceiling
  use hybrid_core_m,  only: disentangle_frozen, rotate_to_subspace, choose_gauge, &
       &                    kinfo_t, build_kinfo, omega_invariant
  use hybrid_grass_m, only: gcg_t, gcg_build, gcg_free, gcg_check_gradient, &
       &                    grassmann_cg, gcg_seed_from_full, gcg_to_full, &
       &                    gcg_seed_svd, &
       &                    GCG_OK, GCG_STALL, GCG_TINYSTEP, GCG_SMALLDEC
  use hybrid_io_m,    only: write_handoff, write_masks, write_anatomy, &
       &                    subspace_continuity, write_dis_status, SUBSPACE_JUMP_SV, &
       &                    continuity_links
  use hybrid_files_m, only: read_nnkp, nnkp_block_bvectors

  !! procedure includes
  use read_vec_m
  use gaunt2_m
  use planew_m
  use l2mmn_m
  use l2amn_m
  use autoprojection_m

  implicit none

  type(inwf_t) :: inwf

  character(len=   11)  :: status,form
  character(len=BUFSZ)  :: deffn, errfn, aline
  character(len=BUFSZ)  :: fname, vecfn, enefn, nnkpfn, iomsg
  complex(C16), allocatable :: proj_h(:,:,:)   ! amn export for MODE=HYBR
  complex(C16), allocatable :: projT_h(:,:,:)  ! full-cell (tails) amn: P source
  ! SOC sigma-loop (def-file unit 9 = second spin-component vector file):
  ! one process reads both spin vectors, accumulates the spinor-summed
  ! overlap in place, and assembles per-sigma amn COLUMN BLOCKS -- the
  ! in-memory replacement for per-spin w2w runs + w2waddsp.  Assumes the
  ! non-spin-polarized SOC layout (shared vsp/inwf/eig); SP+SOC needs
  ! per-sigma vsp and is not handled here.
  ! optional def unit 19 = second spin-component vsp (SP+SOC, e.g. Fe:
  ! the radial potentials genuinely differ per spin; without it the
  ! sigma loop reuses sigma-1's vsp, valid for non-SP SOC only)
  character(len=BUFSZ)  :: vecfn2, vspfn2, vspfn
  character(len=180) :: insofn = ' '
  logical :: sigma2, vsp2
  integer :: isig, nsig, isig_p
  logical :: acc_p, fin_p
  complex(C16), allocatable :: proj_sig(:,:,:), proj_sig_t(:,:,:)

  integer :: iloop, i, j, ios, iproc, irecl, iunit
  integer :: maxx,maxy,maxz, maxg, n, n_pair, nen
  integer :: Nb, Nk, nntot, kkk, iostat
  integer :: pair_index, k1, k2, num1, num2
  integer :: u_sc
  integer, allocatable :: ngh_sc(:,:), gsh_sc(:,:,:)
  complex(C16), allocatable :: Mt_sc(:,:)
  integer :: NEmax, minNEn, bcut, idum

  real(r8) :: efermi, ee, gap_eV, etop_eV
  real(r8), allocatable :: proj_char(:,:,:)
  real(r8), allocatable :: proj_state(:,:,:,:)
  real(r8), allocatable :: proj_state_lo(:,:,:,:)   ! SOC component_sum_max scratch
  integer,  allocatable :: lo_rf_map1(:,:)          ! component-1 LO identity (SOC guard)
  character(len=2*BUFSZ+256) :: socmsg
  real(r8), allocatable :: proj_maxs(:,:)
  real(r8), allocatable :: proj_char_lo(:,:,:), proj_max_lo(:,:)
  integer,  allocatable :: lo_rf_map(:,:)
  real(r8) :: socdegtol, degenv       ! AUTO SOC multiplet-trace degtol (Ry)
  character(len=32) :: degtxt
  integer :: degst, degios
  real(r8), allocatable :: ebmin(:), ebmax(:), eall(:,:), eweight(:,:)
  integer  :: ncut_lo, ncut_hi, bsug
  real(r8) :: gg, gmin_lo, gmin_hi
  ! ---- Gate A (W8) state.  `sent_lo'/`sent_hi' record whether THIS edge was
  ! resolved from an AUTO sentinel (and may therefore be snapped) or typed by
  ! the user (and may not).  nkx_* is the number of k at which the band BEYOND
  ! the edge exists at all -- not Nk, because NEn varies with k and a window
  ! ending at minNEn still truncates at every k that carries more bands.
  integer  :: nkx_lo, nkx_hi, ewcap, emult, bsnap, nsnap
  logical  :: sent_lo, sent_hi, eg_symact, eg_strict, eg_snapany, eg_stuck
  character(len=16)  :: eg_mode
  character(len=256) :: eg_env
  !> degeneracy tolerance for window-edge closure (Ry); 1e-6 Ry ~ 14 ueV,
  !> far below any real splitting and far above the 1e-14 numerical floor
  real(r8), parameter :: DEGTOL = 1.0e-6_r8

  ! auto band-window parameters (energies in eV):
  !   a gap counts as a core/valence separator if it exceeds AUTOW_GAP
  !   and the band above it still starts below E_F - AUTOW_BELOW.
  !   P2 ranking weight: 0.5*erfc(((e - E_F) - P2_ONSET)/P2_WIDTH).
  real(r8), parameter :: AUTOW_GAP = 5.0_R8, AUTOW_BELOW = 1.0_R8
  real(r8), parameter :: AUTOW_ECAP = 40.0_R8
  real(r8), parameter :: P2_ONSET = 5.0_R8, P2_WIDTH = 2.0_R8

  !> AUTO SOC support: multiplet-trace degeneracy tolerance (D5/D4-2,
  !> compute_projectability's component_sum_max mode).  Default 2.0e-6 eV:
  !> measured Kramers splittings on Rh2B (rh2br82) and alpha-Sn are both
  !> <= 9.8e-7 eV, and the smallest physical gap on both is >= 1.03e-6 eV --
  !> 2.0e-6 sits comfortably between.  eV -> Ry with the CODATA constant
  !> (13.605693123, distinct from this file's WIEN-compat Ryd_eV = 13.605698
  !> used for display conversions -- the 3e-7 relative difference is far
  !> below this tolerance's own margin, but the literal is what the default
  !> was calibrated against).  Override: W2W_AUTO_DEGTOL (eV); <= 0 opts out
  !> (compute_projectability falls back to plain state_support_max).
  real(r8), parameter :: SOCDEGTOL_EV2RY   = 1.0_r8 / 13.605693123_r8
  real(r8), parameter :: SOCDEGTOL_DEFAULT = 2.0e-6_r8 * SOCDEGTOL_EV2RY

  type(struct_t) :: stru

!-----------------------------------------------------------------------

  call mpi_w2w_init()

  call gtfnam(deffn,errfn,iproc)
  call errflg(errfn,'Error in W2W')

  sigma2 = .false.
  vsp2   = .false.
  open(unit_def, FILE=deffn, STATUS='old')
  def: do
     read(unit_def, *, end=20) iunit, fname, status, form, irecl

     select case (iunit)
     case (unit_vector)
        vecfn=fname
     case (9)
        ! second spin-component vector file -> enable the sigma loop
        vecfn2=fname
        sigma2=.true.
     case (19)
        ! second spin-component vsp (SP+SOC)
        vspfn2=fname
        vsp2=.true.
     case (unit_inso)
        ! case.inso (optional): the SOC magnetization direction, same file
        ! lapwso reads.  Captured here, parsed after the def loop.
        insofn=fname
     case (50)
        enefn=fname
     case (unit_out, unit_amn, unit_mmn, unit_eig)
        ! Write-only files: only rank 0 opens them
        if (is_root) then
           open(iunit, FILE=fname, STATUS=status, FORM=form, &
                IOSTAT=iostat, IOMSG=iomsg)
           if(iostat /= 0) call croak('error while processing def file `' &
                // trim(deffn) // "': " // trim(iomsg))
        end if
     case (unit_nnkp)
        ! keep the name: MODE=HYBR re-reads the b-vectors after the unit
        ! is closed
        nnkpfn = fname
        open(iunit, FILE=fname, STATUS=status, FORM=form, &
             IOSTAT=iostat, IOMSG=iomsg)
        if(iostat /= 0) call croak('error while processing def file `' &
             // trim(deffn) // "': " // trim(iomsg))
     case default
        open(iunit, FILE=fname, STATUS=status, FORM=form, &
             IOSTAT=iostat, IOMSG=iomsg)

        if(iostat /= 0) call croak('error while processing def file `' &
             // trim(deffn) // "': " // trim(iomsg))
     end select
  end do def
20 close(unit_def)

  ! case.inso: parse the SOC direction if the def named one.  Strict -- a
  ! named-but-unreadable inso croaks rather than silently proceeding without
  ! the classification.  A silent fallback against a stale def is exactly the
  ! failure mode that once produced SOC-free vectors.
  if (len_trim(insofn) > 0) then
     open(unit_inso, FILE=insofn, STATUS='old', IOSTAT=iostat, IOMSG=iomsg)
     if (iostat /= 0) call croak('cannot open inso file `' // trim(insofn) &
          // "': " // trim(iomsg))
     read(unit_inso, *, IOSTAT=iostat)                     ! WFFIL
     if (iostat == 0) read(unit_inso, *, IOSTAT=iostat)    ! llmax, ipr, kpot
     if (iostat == 0) read(unit_inso, *, IOSTAT=iostat)    ! Emin, Emax
     if (iostat == 0) read(unit_inso, *, IOSTAT=iostat) soc_hkl
     if (iostat /= 0) call croak('cannot parse h,k,l from line 4 of `' &
          // trim(insofn) // "'")
     close(unit_inso)
     soc_hkl_set = .true.
  end if

  if (is_root) then
     write(unit_out, '("W2W ", A /)') wien2wannier_version
#ifdef MPI
     write(unit_out, '("MPI enabled: ", I0, " ranks")') nprocs
#endif
  end if

!!!.....READ STRUCT
  call struct_read(unit_struct, stru)

  ! Real/complex is a COMPILE-TIME choice (w2w vs w2wc) with nothing in any
  ! input file to record it, and the wrong choice is silent in one direction
  ! (see vector_kind_probe).  Settle it from the vector file itself, before
  ! any physics runs, and refuse rather than produce a plausible-looking
  ! hand-off built from misread coefficients.
  block
    integer :: vkind
    character(len=BUFSZ+12) :: probefn
    ! With k-parallel vector files (w2w case.def N) the vectors are in vecfn_1 .. vecfn_N
    ! (paropen) and vecfn itself need not exist: probe the first part.
    if (iproc > 0) then
       write(probefn, '(A, "_", I0)') trim(vecfn), 1
    else
       probefn = vecfn
    end if
    vkind = vector_kind_probe(trim(probefn), stru%Nneq)
#ifdef _COMPLEX_
    if (vkind == 1) call croak('this is a REAL vector file but you are &
         &running the COMPLEX binary (w2wc). Rerun with w2w, i.e. drop the &
         &-c from `x w2w`. Real/complex is a compile-time choice and is not &
         &recorded in any input file: '//trim(probefn))
#else
    if (vkind == 2) call croak('this is a COMPLEX vector file but you are &
         &running the REAL binary (w2w). Rerun with w2wc, i.e. `x w2w -c`. &
         &This would otherwise fail SILENTLY -- reading N reals from a record &
         &of N complex numbers is a legal short read -- and yield a hand-off &
         &built from misread coefficients: '//trim(probefn))
#endif
    if (vkind == 0 .and. is_root) write(unit_out,'(a)') &
         ' NOTE: could not probe the vector file kind; proceeding'
  end block

!!!....Find nmat and Nk
!!! Nk could be gotten easier from 'klist' -- can we get nmat
!!! somewhere else?  ('vector'?)
  Nk=0
  NEmax=0
  minNEn=huge(1)

  enefile: do iloop=1,max(iproc,1)
     call paropen(unit_ene, ENEFN, iproc, iloop, STATUS='old')

     header: do I=1,stru%Nneq
        read(unit_ene,*)
        read(unit_ene,*)
     end do header
     kpts: do
        read(unit_ene,'(67X, 2i6)', IOSTAT=ios) N, NEn
        if (ios /= 0) exit kpts
        Nk=Nk+1
        nmat=max(n,nmat)
        NEmax=max(NEn,NEmax)
        minNEn=min(NEn,minNEn)
        do i=1,NEn
           read(unit_ene,*)
        enddo
     end do kpts

     close(unit_ene)
  end do enefile

  ! Second pass: per-band extrema over k (for the auto band window) and
  ! per-(band,k) eigenvalues (for the P2 ranking weight).  Every rank
  ! reads the full set, so no broadcast is needed.
  allocate(ebmin(NEmax), ebmax(NEmax), eall(NEmax, Nk))
  ebmin =  huge(1.0_R8)
  ebmax = -huge(1.0_R8)
  eall  =  huge(1.0_R8)

  kkk=0
  enescan: do iloop=1,max(iproc,1)
     call paropen(unit_ene, ENEFN, iproc, iloop, STATUS='old')

     do I=1,stru%Nneq
        read(unit_ene,*)
        read(unit_ene,*)
     end do
     kpts2: do
        read(unit_ene,'(67X, 2i6)', IOSTAT=ios) N, NEn
        if (ios /= 0) exit kpts2
        kkk=kkk+1
        do i=1,NEn
           read(unit_ene,*) idum, ee
           ebmin(i)=min(ebmin(i),ee)
           ebmax(i)=max(ebmax(i),ee)
           eall(i,kkk)=ee
        enddo
     end do kpts2

     close(unit_ene)
  end do enescan

!!!.....READ INPUT AND POTE
  nnkpt: do
     read(unit_nnkp,'(a80)',end=114) aline
     if(index(aline,'begin kpoints').ne.0) then
        read(unit_nnkp,*) n
        if (n/=Nk) call croak('inconsistent numbers of k-points between&
             & nnkp and energy files')
        cycle nnkpt
     endif

     if(index(aline,'begin nnkpts').ne.0) then
        read(unit_nnkp,*) nntot
        n_pair=Nk*nntot
        call init_pairs(n_pair)
        do i=1,N_pair
           read(unit_nnkp,*) kp(i),kpb(i),bqx1(i),bqy1(i),bqz1(i)
        enddo
        exit nnkpt
     endif
  end do nnkpt
114 close(unit_nnkp)

  if (is_root) then
     write(unit_out,*)'NUM_KPTS=',Nk
     write(unit_out,*)'NNTOT=',NNTOT
     write(unit_out,*)'N_pair=',N_pair
  end if

  call inwf_read(unit_in, inwf)

  read(unit_fermi, *) efermi

  ! Resolve AUTO band-window sentinels from the energy-file scan.
  ! Gate A may only MOVE an edge the code itself chose; a band index the user
  ! typed is never silently overridden (`edges = snap' is the explicit opt-in).
  sent_lo = .false.; sent_hi = .false.
  autowin: if (inwf%auto_mode) then
     sent_lo = (inwf%bmin == 0)
     sent_hi = (inwf%bmax == 0 .or. inwf%bmax == -1)
     if (inwf%bmax == 0) then
        inwf%bmax = minNEn
        if (is_root) write(unit_out, &
             '(" AUTO bmax = ", I0, " (all bands present at every k)")') &
             inwf%bmax
     else if (inwf%bmax == -1) then
        ! energy cap: highest band whose bottom is below E_F + AUTOW_ECAP
        inwf%bmax = 1
        do i = 1, minNEn
           if ((ebmin(i)-efermi)*Ryd_eV <= AUTOW_ECAP) inwf%bmax = i
        end do
        if (is_root) write(unit_out, &
             '(" AUTO bmax = ", I0, " (band bottom <= E_F + ", F5.1, " eV)")') &
             inwf%bmax, AUTOW_ECAP
     else if (inwf%bmax > minNEn) then
        if (is_root) write(unit_out, &
             '(" WARNING: bmax = ", I0, " > ", I0, " bands at every k; clamped")') &
             inwf%bmax, minNEn
        inwf%bmax = minNEn
     end if

     if (inwf%bmin == 0) then
        ! gap-walk: cut above the HIGHEST inter-band gap > AUTOW_GAP whose
        ! upper band still starts below E_F - AUTOW_BELOW.  This trims
        ! core/semicore complexes without ever cutting into the valence
        ! complex (a global-largest-gap rule would, e.g. Sr 4s/4p).
        bcut = 1
        do i = 2, minNEn
           gap_eV  = (ebmin(i) - ebmax(i-1))*Ryd_eV
           etop_eV = (ebmin(i) - efermi)*Ryd_eV
           if (gap_eV > AUTOW_GAP .and. etop_eV < -AUTOW_BELOW) bcut = i
        end do
        inwf%bmin = bcut
        if (is_root) then
           if (bcut > 1) then
              write(unit_out, '(" AUTO bmin = ", I0, " (gap ", F7.2, &
                   & " eV below band ", I0, ")")') bcut, &
                   (ebmin(bcut)-ebmax(bcut-1))*Ryd_eV, bcut
           else
              write(unit_out, '(" AUTO bmin = 1 (no core gap found)")')
           end if
        end if
     end if

     if (inwf%bmin > inwf%bmax) call croak('auto band window is empty')
  end if autowin

  ! Unconditional bmax clamp (previously AUTO-only): a manual-mode bmax
  ! beyond the lowest per-k band count reaches read_vec's uninitialized-A
  ! path.  minNEn comes from the energy-file scan, which runs in all modes.
  if (inwf%bmax > minNEn) then
     if (is_root) write(unit_out, &
          '(" WARNING: bmax = ", I0, " > ", I0, " bands at every k; clamped")') &
          inwf%bmax, minNEn
     inwf%bmax = minNEn
     if (inwf%bmin > inwf%bmax) call croak('band window is empty after bmax clamp')
  end if

  if (is_root) then
     write(unit_out, "(' MODE=')", ADVANCE='no')
     if (inwf%Mmn) write(unit_out, '(" Mmn")', ADVANCE='no')
     if (inwf%Amn) write(unit_out, '(" Amn")', ADVANCE='no')
     if (inwf%auto_mode) write(unit_out, '(" AUTO")', ADVANCE='no')
     write(unit_out,*)
     write(unit_out, '(" band window = [", I0, ", ", I0, "]")') &
          inwf%bmin, inwf%bmax
  end if

  !------------------------------------------------------------------ GATE A
  ! DEGENERACY CLOSURE OF THE WINDOW EDGES.
  !
  ! If a window edge falls INSIDE a degenerate multiplet, the window holds
  ! only part of the multiplet, and no unitary rotation of the states in the
  ! window can restore the degeneracy.  The gate reports such an edge with
  ! the nearest closed edges; `edges = strict' makes it fatal and
  ! `edges = snap' moves the edge to the nearest closed cut.
  !
  ! *** THE TEST IS THE PER-k GAP, NEVER ARITHMETIC ON THE BAND INDEX. ***
  ! Cuts are legal only at multiples of the enforced multiplicity, but NOT
  ! EVERY MULTIPLE IS LEGAL: on MgB2 12x12x8 the cuts after 8, 12, 20 and 24
  ! bands (all multiples of 4) fall inside a multiplet at some k, while the
  ! cuts after 4, 16 and the ODD cut after 21 are clean.
  !
  ! DEGTOL = 1e-6 Ry has to clear the degeneracy floor of the eigensolver
  ! that wrote the energies and stay below every resolved gap:
  !   * lapw1 (scalar), measured on every gap of every band at every k:
  !                       samples   [1e-12,1e-6)   max enforced   min resolved
  !       MgB2 12x12x8    262,517        0          4.44e-14 Ry    1.06e-05 Ry
  !       GaAs  8x8x8     108,057        0          2.90e-13 Ry    1.47e-06 Ry
  !       MgB2  8x8x8      18,546        0          2.18e-14 Ry    1.40e-05 Ry
  !     The bin [1e-12,1e-6) Ry is empty on every cell, so DEGTOL sits in a
  !     void.  The margin above it is thin: on GaAs the nearest resolved gap
  !     is 1.47x DEGTOL, and on MgB2 8x8x8 the cuts 16/17 and 21/22 sit
  !     within two decades.  Do not raise DEGTOL without measuring the gaps
  !     of the cell in hand.
  !   * lapwso (spin-orbit), one dense Pnma cell, not re-measured: the
  !     splitting within a Kramers pair never exceeded 1.462e-07 Ry, while
  !     between pairs the median was 2.8e-4 Ry.  DEGTOL clears that floor by
  !     6.84x; do not lower it toward the lapw1 floor, which would stop
  !     seeing Kramers pairs.  On such a cell 24% of the cuts BETWEEN pairs
  !     are also below DEGTOL and read as inside a multiplet, which is why
  !     `auto' does not make the gate fatal and why the snap walk refuses
  !     rather than walking on.
  !
  ! Hoisted OUT of `is_root': all ranks hold identical eall(NEmax,Nk), so the
  ! decision is the same on every rank.  A croak reached on the root alone
  ! hangs MPI.
  gateA: block
    integer  :: ib, nc1, nkx1, bwant
    real(R8) :: gm1

    ! ---- policy.  Keyword (HYBR tail) first, environment W2W_EDGES on top --
    ! the environment is the only handle a non-HYBR deck has, since the
    ! keyword tail belongs to MODE=HYBR.
    eg_mode = trim(inwf%hyb_edges)
    eg_env = ''
    call get_environment_variable('W2W_EDGES', eg_env)
    if (len_trim(eg_env) > 0) then
       select case (trim(adjustl(eg_env)))
       case ('auto', 'warn', 'strict', 'snap')
          eg_mode = trim(adjustl(eg_env))
       case default
          call croak('W2W_EDGES: auto|warn|strict|snap, got '//trim(eg_env))
       end select
    end if

    ! ---- `auto' never escalates to strict in this release; use
    ! edges=strict/snap (or W2W_EDGES) to make a stuck edge fatal.
    eg_symact = .false.

    select case (trim(eg_mode))
    case ('warn');   eg_strict = .false.; eg_snapany = .false.
    case ('strict'); eg_strict = .true.;  eg_snapany = .false.
    case ('snap');   eg_strict = .true.;  eg_snapany = .true.
    case default                                  ! 'auto'
       eg_strict = eg_symact; eg_snapany = .false.
    end select

    ! ---- walk cap, MEASURED IN THIS RUN, not typed.  Closing a cut through
    ! ONE multiplet of size m costs at most m-1 extra bands, so the largest
    ! multiplet actually present bounds an honest snap.  A walk longer than
    ! that is not closing a multiplet, it is widening the model -- the Bi
    ! band-dense lesson -- and the gate refuses rather than silently growing
    ! the window.  (Measured max multiplet: MgB2 2, GaAs/Fe/SrVO3 3, Bi2D
    ! SOC 4.  The snap-walk length itself has NO measured void -- it runs
    ! 1..61 without a gap -- which is exactly why the cap is derived from
    ! the multiplet size and not from the walk distribution.)
    emult = max_multiplet()
    ewcap = max(1, emult)

    if (is_root) write(unit_out, '(" Gate A: edges = ", A, "  strict = ", &
         & L1, "  symmetry-active = ", L1, "  walk cap = ", I0, &
         & " band(s) (max multiplet ", I0, ")  DEGTOL = ", ES8.1, " Ry")') &
         trim(eg_mode), eg_strict, eg_symact, ewcap, emult, DEGTOL

    nsnap = 0
    call edge_gap_stat(inwf%bmin-1, ncut_lo, gmin_lo, nkx_lo)
    call edge_gap_stat(inwf%bmax,   ncut_hi, gmin_hi, nkx_hi)

    ! ---- snap the TOP edge upward to the nearest closed cut.
    if (ncut_hi > 0 .and. (sent_hi .or. eg_snapany)) then
       bsnap = 0
       do ib = inwf%bmax+1, min(minNEn, inwf%bmax+ewcap)
          call edge_gap_stat(ib, nc1, gm1, nkx1)
          if (nc1 == 0) then; bsnap = ib; exit; end if
       end do
       if (bsnap > 0) then
          bwant = inwf%bmax
          inwf%bmax = bsnap
          nsnap = nsnap + 1
          if (is_root) write(unit_out, '(" Gate A: SNAP top edge bmax ", &
               & I0, " -> ", I0, " (walk ", I0, " <= cap ", I0, &
               & "); the cut at ", I0, "/", I0, " was degenerate at ", I0, &
               & " k")') bwant, bsnap, bsnap-bwant, ewcap, bwant, bwant+1, &
               ncut_hi
          call edge_gap_stat(inwf%bmax, ncut_hi, gmin_hi, nkx_hi)
       else if (is_root) then
          write(unit_out, '(" Gate A: snap REFUSED at the top edge -- no &
               &closed cut within the walk cap of ", I0, " band(s)")') ewcap
       end if
    end if

    ! ---- snap the BOTTOM edge downward.  Band 1 has nothing below it and is
    ! therefore closed by construction, so the walk may always terminate there.
    if (ncut_lo > 0 .and. (sent_lo .or. eg_snapany)) then
       bsnap = 0
       do ib = inwf%bmin-1, max(1, inwf%bmin-ewcap), -1
          if (ib == 1) then; bsnap = 1; exit; end if
          call edge_gap_stat(ib-1, nc1, gm1, nkx1)
          if (nc1 == 0) then; bsnap = ib; exit; end if
       end do
       if (bsnap > 0) then
          bwant = inwf%bmin
          inwf%bmin = bsnap
          nsnap = nsnap + 1
          if (is_root) write(unit_out, '(" Gate A: SNAP bottom edge bmin ", &
               & I0, " -> ", I0, " (walk ", I0, " <= cap ", I0, &
               & "); the cut at ", I0, "/", I0, " was degenerate at ", I0, &
               & " k")') bwant, bsnap, bwant-bsnap, ewcap, bwant-1, bwant, &
               ncut_lo
          call edge_gap_stat(inwf%bmin-1, ncut_lo, gmin_lo, nkx_lo)
       else if (is_root) then
          write(unit_out, '(" Gate A: snap REFUSED at the bottom edge -- no &
               &closed cut within the walk cap of ", I0, " band(s)")') ewcap
       end if
    end if

    if (nsnap > 0) then
       if (inwf%bmin > inwf%bmax) call croak('Gate A: band window is empty &
            &after the edge snap')
       if (is_root) write(unit_out, &
            '(" band window (after Gate A snap) = [", I0, ", ", I0, "]")') &
            inwf%bmin, inwf%bmax
    end if

    eg_stuck = (ncut_lo > 0 .or. ncut_hi > 0)

    if (is_root) then
       if (eg_stuck) then
          write(unit_out, '(" WARNING: the band window CUTS a degenerate", &
               & " multiplet")')
          if (ncut_lo > 0) write(unit_out, &
               '("   bottom edge: bands ", I0, "/", I0, " degenerate at ", &
               & I0, " of ", I0, " k-points  (min gap ", ES9.2, " Ry)")') &
               inwf%bmin-1, inwf%bmin, ncut_lo, nkx_lo, gmin_lo
          if (ncut_hi > 0) write(unit_out, &
               '("   top edge:    bands ", I0, "/", I0, " degenerate at ", &
               & I0, " of ", I0, " k-points  (min gap ", ES9.2, " Ry)")') &
               inwf%bmax, inwf%bmax+1, ncut_hi, nkx_hi, gmin_hi
          ! nearest CLOSED edges, reported whether or not they were taken --
          ! the ranked-alternative pattern of autoprojection.diag, so an
          ! automated retry needs no human.
          !
          ! *** BOTH DIRECTIONS.  CORRECTED 2026-08-26 (R4 finalize). ***
          ! The first implementation searched UPWARD only, so the ranked
          ! alternative was missing in exactly the configuration that needs
          ! it most: the AUTO sentinel `bmax = 0' resolves to minNEn, so
          ! when that last cut is stuck there is no upward closed edge BY
          ! CONSTRUCTION and the gate croaked saying only "no closed top
          ! edge exists".  MEASURED on the shipped mgH deck (bmax = minNEn
          ! = 32, stuck at 2 of 511 k): a downward alternative exists and is
          ! good -- bmax = 29, min gap 5.7075e-04 Ry, 571x DEGTOL.  The SNAP
          ! still refuses to move down unasked (that SHRINKS the model and
          ! is the user's call), but the report must not be silent.
          if (ncut_hi > 0) then
             bsug = 0
             do ib = inwf%bmax+1, minNEn
                call edge_gap_stat(ib, nc1, gm1, nkx1)
                if (nc1 == 0) then; bsug = ib; exit; end if
             end do
             if (bsug > 0) then
                write(unit_out, '("   nearest closed top edge: bmax = ", I0, &
                     & "  (walk ", I0, ", cap ", I0, ")")') &
                     bsug, bsug-inwf%bmax, ewcap
             else
                write(unit_out, '("   no closed top edge exists ABOVE bmax = "&
                     &, I0, " (minNEn = ", I0, ")")') inwf%bmax, minNEn
                bsug = 0
                do ib = inwf%bmax-1, inwf%bmin, -1
                   call edge_gap_stat(ib, nc1, gm1, nkx1)
                   if (nc1 == 0) then; bsug = ib; exit; end if
                end do
                if (bsug > 0) then
                   call edge_gap_stat(bsug, nc1, gm1, nkx1)
                   write(unit_out, '("   nearest closed top edge BELOW bmax: "&
                        &, I0, "  (walk ", I0, " down, min gap ", ES9.2, &
                        & " Ry) -- narrows the model, so the snap will not &
                        &take it unasked; set bmax = ", I0)') &
                        bsug, inwf%bmax-bsug, gm1, bsug
                else
                   write(unit_out, '("   no closed top edge exists anywhere &
                        &in [", I0, ", ", I0, "] -- widen the window or use &
                        &a denser/less degenerate mesh")') inwf%bmin, minNEn
                end if
             end if
          end if
          if (ncut_lo > 0) then
             bsug = 0
             do ib = inwf%bmin-1, 1, -1
                if (ib == 1) then; bsug = 1; exit; end if
                call edge_gap_stat(ib-1, nc1, gm1, nkx1)
                if (nc1 == 0) then; bsug = ib; exit; end if
             end do
             if (bsug > 0) write(unit_out, &
                  '("   nearest closed bottom edge: bmin = ", I0, &
                  & "  (walk ", I0, ", cap ", I0, ")")') &
                  bsug, inwf%bmin-bsug, ewcap
          end if
          write(unit_out, '("   symmetry-adapted Wannier functions are NOT", &
               & " possible with this window")')
          if (.not. eg_strict) write(unit_out, '("   Gate A: PROCEEDING &
               &(edges = ", A, "); set edges = strict, or edges = snap to &
               &let the gate move the edge")') trim(eg_mode)
       else
          ! an edge with no band beyond it cannot be cut; report it as such
          ! rather than printing huge(), which looks like a defect
          write(unit_out, '(" band window is degeneracy-closed at both edges")')
          if (inwf%bmin > 1) then
             write(unit_out, '("   min gap below bmin: ", ES9.2, " Ry  (", &
                  & I0, " k)")') gmin_lo, nkx_lo
          else
             write(unit_out, '("   bottom edge: bmin = 1, nothing below")')
          end if
          if (nkx_hi > 0) then
             write(unit_out, '("   min gap above bmax: ", ES9.2, " Ry  (", &
                  & I0, " k)")') gmin_hi, nkx_hi
          else
             write(unit_out, '("   top edge: bmax is the last band available")')
          end if
       end if
    end if

    ! The gate's verdict must survive whatever the run does next: everything
    ! after this point is expensive and can be killed, run out of disk or die
    ! in the vector read, and a verdict still sitting in the output buffer
    ! would be lost exactly when it is most wanted.
    if (is_root) flush(unit_out)

    ! fail-closed, on EVERY rank
    if (eg_stuck .and. eg_strict) call croak('Gate A: a band-window edge &
         &cuts a degenerate multiplet -- symmetry-adapted Wannier functions &
         &are not possible with this window (see the Gate A block for the &
         &nearest closed edges; `edges = warn'' overrides, `edges = snap'' &
         &lets the gate move the edge)')
  end block gateA

  Nb = inwf%bmax - inwf%bmin + 1

  call init_xa3(Nb,nmat,Nk+1)
  call init_xa(LMAX2,NMAT,NRAD,Nb)
  call init_Amn_Mmn(Nb, N_pair)

  kkk=0
  maxx=0; maxy=0; maxz=0
  vectorfiles: do iloop=1, max(iproc, 1)
     call paropen(unit_vector, vecfn, iproc, iloop, &
          STATUS='old', FORM='unformatted')

     call read_vec(inwf%Bmin, inwf%Bmax, stru%Nneq, &
          &        kkk, maxx, maxy, maxz, Efermi)

     close(unit_vector)
  end do vectorfiles

  if(kkk /= Nk) call croak('inconsistent numbers of k-points between&
       & vector and energy files')

  ! (The historical rank-0 broadcast of maxx/vecsz/XK/GX/A was REDUNDANT:
  !  the vectorfiles loop above has no root guard, so every rank already
  !  read every file.  Deleted: on a 5200-k cell it re-sent 2.37 GB of
  !  data that all ranks already held.)

  call init_bessel(Lmax2, inwf%LJmax, Nrad, Nrf)
  call gaunt2
  call init_gaunt_cache()

  ! AUTO mode: compute projectability and generate projections
  ! Open vector file for atpar/almgen calls in compute_projectability
  auto_proj: if (inwf%auto_mode) then
     call paropen(unit_vector, vecfn, iproc, 1, &
          &       STATUS='old', FORM='unformatted')
     if (is_root) then
        write(unit_out, *)
        write(unit_out, '(A)') ' AUTO mode: computing projectability...'
     end if

     allocate(proj_char(Nb, stru%Nneq, 0:Lmax2))

     ! P2: occupation-like erfc weight suppressing high-energy NFE states
     ! in the (atom, l) channel RANKING only; the per-state output and the
     ! .amn stay unweighted.
     if (inwf%auto_config%p2_weight) then
        allocate(eweight(Nb, Nk))
        do kkk = 1, Nk
           do i = 1, Nb
              ee = (eall(inwf%bmin-1+i, kkk) - efermi)*Ryd_eV
              eweight(i, kkk) = 0.5_R8*erfc((ee - P2_ONSET)/P2_WIDTH)
           end do
        end do
        if (is_root) write(unit_out, '(" P2 ranking weight: erfc onset E_F+", &
             & F4.1, " eV, width ", F4.1, " eV")') P2_ONSET, P2_WIDTH
     end if

     allocate(proj_maxs(stru%Nneq, 0:Lmax2))
     allocate(proj_char_lo(Nb, stru%Nneq, 0:Lmax2))
     allocate(proj_max_lo(stru%Nneq, 0:Lmax2))
     allocate(lo_rf_map(stru%Nneq, 0:Lmax2))
     ! proj_state is the .autoproj payload (write_proj) AND, for SOC, the
     ! component_sum_max scratch the channel support is derived from
     if (inwf%auto_config%write_proj .or. sigma2) &
          allocate(proj_state(Nb, Nk, stru%Nneq, 0:Lmax2))
     if (sigma2) allocate(proj_state_lo(Nb, Nk, stru%Nneq, 0:Lmax2))
     if (sigma2 .and. inwf%auto_config%write_proj .and. is_root) &
          write(unit_out, '(A,F10.1,A)') ' AUTO SOC: per-state scratch ', &
          2.0_r8*Nb*Nk*stru%Nneq*(Lmax2+1)*8.0_r8/2.0_r8**20, ' MiB'
     ! AUTO SOC support: multiplet-trace degeneracy tolerance (D5/D4-2),
     ! passed to compute_projectability's SOC (sigma2) calls only -- the
     ! non-SOC calls below are untouched.  See SOCDEGTOL_DEFAULT above.
     if (sigma2) then
        socdegtol = SOCDEGTOL_DEFAULT
        call get_environment_variable('W2W_AUTO_DEGTOL', degtxt, status=degst)
        if (degst == 0) then
           read(degtxt,*,iostat=degios) degenv
           if (degios == 0) socdegtol = degenv * SOCDEGTOL_EV2RY
        end if
     end if
     ! SOC: the channel-selection projectability must be the SPINOR sum
     ! over both vector-file components (vecfn, vecfn2).  Run one pass per
     ! component, re-filling the module A array between passes and letting
     ! compute_projectability accumulate (|Alm_s1|^2 + |Alm_s2|^2); the MPI
     ! reduce + normalization are deferred to the final pass.  The support
     ! maxima are taken from the component-SUMMED per-state arrays
     ! (component_sum_max), i.e. max_{n,k}[P_1 + P_2], NOT max-merged per
     ! component (D4).  Non-SOC = a single pass, identical to before,
     ! streaming maxima and all.  For SP+SOC the second component
     ! has its own radial basis, so the sigma-2 pass reads vecfn2 (its
     ! linearization energies) and, if unit 19 was given, swaps unit_vsp to
     ! vspfn2 (the down-spin potential) -- mirroring the l2Amn sigma loop.
     if (vsp2) inquire(unit_vsp, name=vspfn)
     sigproj: do isig_p = 1, merge(2, 1, sigma2)
        if (isig_p == 2) then
           close(unit_vector)
           kkk = 0
           do iloop = 1, max(iproc, 1)
              call paropen(unit_vector, vecfn2, iproc, iloop, &
                   STATUS='old', FORM='unformatted')
              call read_vec(inwf%Bmin, inwf%Bmax, stru%Nneq, &
                   kkk, maxx, maxy, maxz, Efermi, write_eig=.false.)
              close(unit_vector)
           end do
           call paropen(unit_vector, vecfn2, iproc, 1, &
                STATUS='old', FORM='unformatted')
           if (vsp2) then
              close(unit_vsp)
              open(unit_vsp, FILE=vspfn2, STATUS='old', FORM='formatted', &
                   IOSTAT=iostat, IOMSG=iomsg)
              if (iostat /= 0) call croak('AUTO sigma-2 vsp: '//trim(iomsg))
           end if
        end if
        acc_p = (isig_p == 2)
        fin_p = (isig_p == merge(2, 1, sigma2))
        if (sigma2) then
           ! SOC: spinor-summed support maxima; proj_state is only MPI-reduced
           ! when it is also the .autoproj payload
           if (inwf%auto_config%p2_weight) then
              call compute_projectability(stru, inwf, Nk, proj_char, &
                   &          proj_state, eweight, proj_max=proj_maxs, &
                   &          proj_char_lo=proj_char_lo, proj_max_lo=proj_max_lo, &
                   &          lo_rf_out=lo_rf_map, accumulate=acc_p, finalize=fin_p, &
                   &          proj_state_lo=proj_state_lo, component_sum_max=.true., &
                   &          reduce_proj_state=inwf%auto_config%write_proj, &
                   &          degtol=socdegtol)
           else
              call compute_projectability(stru, inwf, Nk, proj_char, &
                   &          proj_state, proj_max=proj_maxs, &
                   &          proj_char_lo=proj_char_lo, proj_max_lo=proj_max_lo, &
                   &          lo_rf_out=lo_rf_map, accumulate=acc_p, finalize=fin_p, &
                   &          proj_state_lo=proj_state_lo, component_sum_max=.true., &
                   &          reduce_proj_state=inwf%auto_config%write_proj, &
                   &          degtol=socdegtol)
           end if
           if (isig_p == 1) then
              allocate(lo_rf_map1(stru%Nneq, 0:Lmax2))
              lo_rf_map1 = lo_rf_map
           end if
        else if (inwf%auto_config%write_proj) then
           if (inwf%auto_config%p2_weight) then
              call compute_projectability(stru, inwf, Nk, proj_char, &
                   &          proj_state, eweight, proj_max=proj_maxs, &
                   &          proj_char_lo=proj_char_lo, proj_max_lo=proj_max_lo, &
                   &          lo_rf_out=lo_rf_map, accumulate=acc_p, finalize=fin_p)
           else
              call compute_projectability(stru, inwf, Nk, proj_char, &
                   &          proj_state, proj_max=proj_maxs, &
                   &          proj_char_lo=proj_char_lo, proj_max_lo=proj_max_lo, &
                   &          lo_rf_out=lo_rf_map, accumulate=acc_p, finalize=fin_p)
           end if
        else
           if (inwf%auto_config%p2_weight) then
              call compute_projectability(stru, inwf, Nk, proj_char, &
                   &          eweight=eweight, proj_max=proj_maxs, &
                   &          proj_char_lo=proj_char_lo, proj_max_lo=proj_max_lo, &
                   &          lo_rf_out=lo_rf_map, accumulate=acc_p, finalize=fin_p)
           else
              call compute_projectability(stru, inwf, Nk, proj_char, &
                   &          proj_max=proj_maxs, &
                   &          proj_char_lo=proj_char_lo, proj_max_lo=proj_max_lo, &
                   &          lo_rf_out=lo_rf_map, accumulate=acc_p, finalize=fin_p)
           end if
        end if
     end do sigproj
     if (sigma2) then
        ! LO radial-identity guard: the sigma-2 pass overwrites lo_rf_map
        ! with ITS chosen LO radial per (atom, l).  The summed proj_state_lo
        ! only means something if both components chose the same radial;
        ! two different LO radials cannot be summed under one candidate
        ! identity.  Fail closed, no override.
        do i = 1, stru%Nneq
           do j = 0, Lmax2
              if (lo_rf_map1(i, j) /= lo_rf_map(i, j)) then
                 write(socmsg, '(a,i0,a,i0,a,a,a,i0,a,a,a,i0,a)') &
                      'AUTO SOC: LO radial identity differs between spinor ' // &
                      'components at inequivalent atom ', i, ', l=', j, &
                      ': component 1 (', trim(vecfn), ') rf=', lo_rf_map1(i, j), &
                      ', component 2 (', trim(vecfn2), ') rf=', lo_rf_map(i, j), &
                      ' -- two different LO radials cannot be summed under ' // &
                      'one candidate identity'
                 call croak(trim(socmsg))
              end if
           end do
        end do
        deallocate(lo_rf_map1)
        ! restore A = sigma-1 for the downstream l2Mmn/l2Amn sigma loop
        ! (which relies on the module A array holding component 1)
        close(unit_vector)
        kkk = 0
        do iloop = 1, max(iproc, 1)
           call paropen(unit_vector, vecfn, iproc, iloop, &
                STATUS='old', FORM='unformatted')
           call read_vec(inwf%Bmin, inwf%Bmax, stru%Nneq, &
                kkk, maxx, maxy, maxz, Efermi, write_eig=.false.)
           close(unit_vector)
        end do
        call paropen(unit_vector, vecfn, iproc, 1, &
             STATUS='old', FORM='unformatted')
        if (vsp2) then
           ! restore the up-spin potential for the downstream main path
           close(unit_vsp)
           open(unit_vsp, FILE=vspfn, STATUS='old', FORM='formatted', &
                IOSTAT=iostat, IOMSG=iomsg)
           if (iostat /= 0) call croak('AUTO vsp restore: '//trim(iomsg))
        end if
     end if
     if (inwf%auto_config%write_proj .and. is_root) &
          call write_projstate(stru, inwf, Nb, Nk, proj_state)
     if (allocated(proj_state))    deallocate(proj_state)
     if (allocated(proj_state_lo)) deallocate(proj_state_lo)
     if (allocated(eweight)) deallocate(eweight)
     if (is_root) then
        do i = 1, stru%Nneq
           write(unit_out, '(A,I3,A,6I3,A,6F8.3)') ' LO-shell atom', i, &
                ' rf:', (lo_rf_map(i,j), j=0,5), '  maxP:', (proj_max_lo(i,j), j=0,5)
        end do
     end if
     call generate_projections(stru, proj_char, Nb, Nk, &
          inwf%auto_config, inwf, proj_max=proj_maxs, &
          proj_char_lo=proj_char_lo, proj_max_lo=proj_max_lo, &
          lo_rf_in=lo_rf_map)
     deallocate(proj_char, proj_maxs, proj_char_lo, proj_max_lo, lo_rf_map)

     ! Sidecar map of .amn columns -> (global atom, l, m) so downstream
     ! tools can emit position-specific projection cards (guiding centres)
     ! in hand-off column order.  AUTO projections always have NY = 1.
     if (is_root) then
        open(96, FILE='projmap.dat', STATUS='replace', FORM='formatted')
        write(96, '(A)') '# amn column -> iatom(global) l m rf'
        ! cliff mode records the otsu mask boundary (columns of the
        ! physical cluster); in-process consumers read the payload directly
        if (inwf%auto_config%threshold <= 0 .and. &
             inwf%auto_config%mask_cols > 0) then
           write(96, '(A,I5)') '# mask_cols ', inwf%auto_config%mask_cols
        end if
        do i = 1, inwf%Nproj
           j = 1
           if (allocated(inwf%projections(i)%rf)) j = inwf%projections(i)%rf(1)
           write(96, '(4I6)') inwf%projections(i)%iat(1), &
                inwf%projections(i)%l(1), inwf%projections(i)%m(1), j
        end do
        close(96)
     end if

     close(unit_vector)
     rewind(unit_vsp)
     if (is_root) call ptime('auto-projection')
  end if auto_proj

  ! ---- HYBR capacity gate (BEFORE l2Amn/l2Mmn: decides/validates nw) ----
  ! The frozen-manifold capacity constraint depends only on the window
  ! eigenvalues + froz_emax, both known here -- so the decision is made
  ! before any expensive amn/mmn work (adding channels later would force
  ! an amn recompute).  need = max_k count(E - E_F <= froz_emax), already
  ! in WF units (SOC energy files list spinor states).
  !   num_wann = 0 -> auto: nw = sigma*mask_cols, extended down the ranked
  !     channel list in whole spatial channels until nw >= need + margin
  !   num_wann > 0 -> explicit: hard error below need + margin, with the
  !     required value in the message (raise num_wann or lower froz_emax)
  hyb_gate: if (inwf%hybrid) then
     block
       integer :: need, marg, sf, nwcap, jch, kg, ig, ncols
       logical :: explicit_fe
       real(R8) :: fe
       explicit_fe = trim(inwf%hyb_froz_mode) /= 'auto' .and. &
                     trim(inwf%hyb_froz_mode) /= 'none'
       sf = 1; if (sigma2) sf = 2
       if (.not. explicit_fe) then
          ! auto/none froz_emax resolves only once P exists (walk-down);
          ! the post-mask capacity assertion in hybrid_stage2 still guards
          if (inwf%num_wann == 0) call croak('HYBR: num_wann = 0 (auto) &
               &requires an explicit froz_emax = <eV> keyword')
       else
          fe = inwf%hyb_froz_val
          need = 0
          do kg = 1, Nk
             ncols = 0
             do ig = inwf%bmin, inwf%bmax
                if ((eall(ig,kg) - efermi)*Ryd_eV <= fe) ncols = ncols + 1
             end do
             need = max(need, ncols)
          end do
          marg = inwf%hyb_margin
          if (marg <= 0) marg = max(4, ceiling(0.04_R8*need))
          if (inwf%num_wann == 0) then
             ! auto-nw = SMALLEST channel-aligned nw >= need+margin,
             ! walking the ranked channel list from the top.  (Bi2D
             ! nw=88 lesson: the Otsu mask block is a candidate set,
             ! not a model size -- flooring nw at mask_cols saturated
             ! the window and hurt BOTH localization and eta.)
             associate (ac => inwf%auto_config)
               nwcap = 0
               jch = 0
               do while (nwcap < need + marg .and. jch < ac%nsel)
                  jch = jch + 1
                  nwcap = nwcap + sf*(2*ac%sel_l(jch)+1)*stru%mult(ac%sel_atom(jch))
               end do
               if (nwcap < need + marg) call croak('HYBR auto-nw: ranked &
                    &candidate pool exhausted below the frozen capacity; &
                    &widen the band window or lower froz_emax')
               if (nwcap > inwf%bmax-inwf%bmin+1) call croak('HYBR auto-nw &
                    &exceeds the band window; widen bmax or lower froz_emax')
               inwf%num_wann = nwcap
               if (is_root) write(unit_out, '(" HYBR capacity gate: need=", &
                    & I0, " (froz_emax=", F7.3, " eV)  margin=", I0, &
                    & "  auto num_wann=", I0, "  (", I0, " of ", I0, &
                    & " ranked channels)")') need, fe, marg, nwcap, jch, &
                    ac%nsel
               ! disentanglement-headroom advisory: the pool wants
               ! pool_factor*nw bands of slack inside the window
               if (real(inwf%bmax-inwf%bmin+1,R8) < &
                    inwf%hyb_pool_factor*real(nwcap,R8) .and. is_root) &
                    write(unit_out, '(" HYBR capacity gate: WARNING window ", &
                    & I0, " bands < pool_factor*nw = ", F8.1, &
                    & " -> limited disentanglement slack; consider widening &
                    &bmax")') inwf%bmax-inwf%bmin+1, &
                    inwf%hyb_pool_factor*real(nwcap,R8)
             end associate
          else
             ! explicit nw: below need = violation (croak by default;
             ! `capacity = warn` proceeds -- Bi2D nw=12 shows a guarded
             ! under-capacity model CAN be both best-localized and most
             ! faithful; Rh3Si4 shows it can also fail at 434 meV, so
             ! the override is explicit, never silent); inside the
             ! margin = valid but tight (warn only)
             if (inwf%num_wann < need) then
                if (is_root) write(unit_out, '(" HYBR capacity gate: nw=", &
                     & I0, " < need=", I0, " at froz_emax=", &
                     & F7.3, " eV -> raise num_wann to >= ", I0, &
                     & " or lower froz_emax, or set capacity = warn")') &
                     inwf%num_wann, need, fe, need+marg
                if (inwf%hyb_cap_strict) then
                   call croak('HYBR: num_wann below frozen capacity &
                        &(see the capacity gate line; capacity = warn &
                        &overrides)')
                else if (is_root) then
                   write(unit_out, '(" HYBR capacity gate: PROCEEDING &
                        &under capacity (capacity = warn); band fidelity &
                        &near froz_emax is NOT guaranteed -- check eta")')
                end if
             else if (inwf%num_wann < need + marg) then
                if (is_root) write(unit_out, '(" HYBR capacity gate: nw=", &
                     & I0, " >= need=", I0, " but inside margin ", I0, &
                     & " at froz_emax=", F7.3, " eV  (tight -- consider ", &
                     & I0, ")")') inwf%num_wann, need, marg, fe, need+marg
             else
                if (is_root) write(unit_out, '(" HYBR capacity gate: nw=", &
                     & I0, " >= need+margin=", I0, "+", I0, " at froz_emax=", &
                     & F7.3, " eV  OK")') inwf%num_wann, need, marg, fe
             end if
          end if
       end if
       !--------------------------------------------------------------- GATE B
       ! (1/2) THE POOL-RATIO LINE.  Printed in BOTH the auto-nw and the
       ! explicit-nw branches -- and also when froz_emax is auto|none, where
       ! the capacity gate itself cannot speak because the ceiling only
       ! exists after the walk-down.  Until W8 the explicit-nw path printed
       ! no pool message at all and the ratio npool/nw was never printed
       ! anywhere, while select_pool CLAMPS the pool silently against the
       ! window top (hybrid_prep.F select_pool, the `npool = min(npool,
       ! nb-s0+1)' line).  The standing user requirement is
       ! num_dft >= 1.5 x num_wann.
       if (is_root .and. inwf%num_wann > 0) then
          write(unit_out, '(" HYBR pool: window Nb=", I0, "  num_wann=", I0, &
               & "  ratio=", F6.2, "  pool_factor=", F6.2)') &
               Nb, inwf%num_wann, real(Nb,R8)/real(inwf%num_wann,R8), &
               inwf%hyb_pool_factor
          if (real(Nb,R8) < inwf%hyb_pool_factor*real(inwf%num_wann,R8)) &
               write(unit_out, '("   the WINDOW itself is shorter than &
               &pool_factor*nw = ", F8.1, " bands; the pool will be clamped &
               &-- widen bmax to >= ", I0)') &
               inwf%hyb_pool_factor*real(inwf%num_wann,R8), &
               inwf%bmin - 1 + ceiling(inwf%hyb_pool_factor* &
               real(inwf%num_wann,R8))
       end if
     end block
  end if hyb_gate

  read_proj: if (inwf%Amn .and. is_root) then
     write(unit_out, *)
     if (inwf%auto_mode) then
        write(unit_out,*)'Auto-selected orbital projections:'
     else
        write(unit_out,*)'Initial orbital projections:'
     end if

     do i=1,inwf%Nproj
        write(unit_out, *)
        if (i==1) then
           write(unit_out, "(' atom  l  m    Re      Im')")
           write(unit_out, "(' ', 28('-'))")
        end if

        do j=1, inwf%projections(i)%NY
           write(unit_out, '(i5,2x,i1,1x,i2,"  (",f6.3,", ",f6.3,")"  )') &
                inwf%projections(i)%iat(j), inwf%projections(i)%l    (j), &
                inwf%projections(i)%m  (j), inwf%projections(i)%coeff(j)
        enddo
     enddo
  endif read_proj

  if (is_root) then
     write(unit_out, "(////,30X,50(1H-),/,33X,'S T R U C T U R A L   ', &
            & 'I N F O R M A T I O N',/,30X,50(1H-),//)")
     write(unit_out, "(3X,'SUBSTANCE',20X,'= ',A80,/)") stru%title
  end if

  if (inwf%Amn .and. is_root) then
     write(unit_amn,'(A20)')  stru%title
     ! sigma loop: the combined amn holds one column BLOCK per spin
     write(unit_amn,'(3I12)') Nb, Nk, merge(2,1,sigma2)*inwf%Nproj
  endif
  if (inwf%Mmn .and. is_root) then
     write(unit_mmn,'(A20)')  stru%title
     write(unit_mmn,'(3I12)') Nb, Nk, Nntot
  endif

  if (is_root) then
     write(unit_out, "(3X,'LATTICE',22X,'= ',A4)")                  stru%lattic
     write(unit_out, "(3X,'LATTICE CONSTANTS ARE',8X,'= ',3F12.7)") stru%a
     write(unit_out, "(3X,'NUMBER OF ATOMS IN UNITCELL  = ',I3)")   stru%Nneq
     write(unit_out, "(3X,'MODE OF CALCULATION IS',7X,'= ',A4)")    stru%mode

     write(unit_out,*)' gbas,  stru2frac'
     do i=1,3
        write(unit_out, '(3f10.5,3x,3i3)') &
             stru%conv_rec(i,:), stru%stru2frac(i,:)
     end do
  end if

  !     rotate boundary vectors
  do i=1,N_pair
     bqx(i) = dot_product(stru%stru2frac(1,:), (/bqx1(i), bqy1(i), bqz1(i)/))
     bqy(i) = dot_product(stru%stru2frac(2,:), (/bqx1(i), bqy1(i), bqz1(i)/))
     bqz(i) = dot_product(stru%stru2frac(3,:), (/bqx1(i), bqy1(i), bqz1(i)/))
  enddo

  !.....CALCULATE CHARGE DENSITY CLM(R) IN SPHERES,  PARTIAL CHARGES

  ! l2mmn, l2amn need vector file to be open
  call paropen(unit_vector, vecfn, iproc, 1, &
       &       STATUS='old', FORM='unformatted')

  if (is_root) call ptime(unit_out)

  nsig = 1
  if (sigma2) then
     nsig = 2
     if (is_root) write(unit_out, '(" SIGMA loop: spinor-summed overlap +&
          & per-sigma amn column blocks (", A, " / ", A, ")")') &
          trim(vecfn), trim(vecfn2)
  end if

  sigmaloop: do isig = 1, nsig
     if (isig == 2) then
        ! second spin component: re-read its vectors into xa3 in place
        ! (same basis and eigenvalues in SOC; only spinor parts differ)
        close(unit_vector)
        kkk = 0
        do iloop = 1, max(iproc, 1)
           call paropen(unit_vector, vecfn2, iproc, iloop, &
                STATUS='old', FORM='unformatted')
           call read_vec(inwf%Bmin, inwf%Bmax, stru%Nneq, &
                &        kkk, maxx, maxy, maxz, Efermi, write_eig=.false.)
           close(unit_vector)
        end do
        if (kkk /= Nk) call croak('sigma-2 vector/energy k-count mismatch')
        call paropen(unit_vector, vecfn2, iproc, 1, &
             &       STATUS='old', FORM='unformatted')
        if (vsp2) then
           ! SP+SOC: swap in the sigma-2 radial potentials
           close(unit_vsp)
           open(unit_vsp, FILE=vspfn2, STATUS='old', FORM='formatted', &
                IOSTAT=iostat, IOMSG=iomsg)
           if (iostat /= 0) call croak('sigma-2 vsp: '//trim(iomsg))
        end if
     end if

     if (inwf%Mmn) then
        ! l2mmn reads the vsp ISCF header with NO rewind of its own; after
        ! the AUTO pass (or a previous sigma's l2Amn) vsp is consumed
        rewind(unit_vsp)
        call l2mmn(stru, inwf, Nk, Nntot)
        if (is_root) call ptime('l2Mmn')
        !JXZ: MAXG is not pre-defined
        MAXG = 0
        if (is_root) write(unit_out,*)'MXG=',MAXG
        call planew(stru, Nb, Nk, Nntot, maxx+1, maxy+1, maxz+1)
        if (is_root) call ptime('planew')
     end if

     if (inwf%Amn) then
        if (is_root) call ptime(unit_out)
        if (sigma2) then
           if (inwf%hybrid) then
              call l2amn(stru, inwf, Nk, proj_out=proj_sig, &
                   write_files=.false., proj_full_out=proj_sig_t)
           else
              call l2amn(stru, inwf, Nk, proj_out=proj_sig, write_files=.false.)
           end if
           if (is_root) then
              if (.not. allocated(proj_h)) &
                   allocate(proj_h(Nb, 2*inwf%Nproj, Nk))
              proj_h(:, (isig-1)*inwf%Nproj+1:isig*inwf%Nproj, :) = proj_sig
              if (inwf%hybrid) then
                 if (.not. allocated(projT_h)) &
                      allocate(projT_h(Nb, 2*inwf%Nproj, Nk))
                 projT_h(:, (isig-1)*inwf%Nproj+1:isig*inwf%Nproj, :) = proj_sig_t
              end if
           end if
        else if (inwf%hybrid) then
           call l2amn(stru, inwf, Nk, proj_out=proj_h, proj_full_out=projT_h)
        else
           call l2amn(stru, inwf, Nk)
        end if
        if (is_root) call ptime('l2Amn')
     end if
  end do sigmaloop

  if (inwf%Mmn) then
     ! Reduce overlap from all ranks to rank 0, then write MMN
     ! (after the sigma loop: the accumulated overlap IS the spin sum)
     call mpi_w2w_reduce_sum_complex_inplace(overlap, &
          int(Nb,8)*int(Nb,8)*int(N_pair,8))   ! promote BEFORE multiplying
     if (is_root) then
        ! The TEXT .mmn is an INTERMEDIATE under MODE=HYBR, not a product:
        ! Stage-2 takes its overlaps from the in-memory `overlap` array and
        ! Wannier90 consumes <case>-hyb.{amn,mmn,eig}, so nothing in-process or
        ! in the hand-off chain reads this file.  The binary sidecar written
        ! just below carries the same numbers losslessly and is what every
        ! in-tree reader prefers (hybrid_files.F treats a MISSING text file as
        ! normal, not as staleness).  Writing it costs 1.1 GB / ~11 s on GaAs
        ! and an estimated 12 GB / ~130 s at Rh3Si4 scale.
        !
        ! It is still the only way to feed a stock-Wannier90 FULL-WINDOW run
        ! (the hybrid-vs-window comparison), so it stays available behind
        ! `text_io = true`, and is always written when NOT in HYBR mode.
        if (inwf%hybrid .and. .not. inwf%hyb_text_io) then
           write(unit_out,'(a)') ' HYBR: text .mmn skipped (binary sidecar &
                &only); set text_io = true to keep it'
        else
           pair_index = 0
           do k1=1,Nk
              do k2=1,Nntot
                 pair_index = pair_index + 1
                 write(unit_mmn,'(5i8)') KP(pair_index), KPB(pair_index), &
                      BQX1(pair_index), BQY1(pair_index), BQZ1(pair_index)
                 ! One WRITE per (k,b) block instead of Nb*Nb of them.  FORMAT
                 ! REVERSION restarts the format per complex, so the record
                 ! stream is byte-identical -- but the repeat count MUST become
                 ! 2: with an implied-do, '(32f18.12)' would pack 16 complexes
                 ! onto a single record.
                 write(unit_mmn,'(2f18.12)') &
                      ((conjg(overlap(num2,num1,pair_index)), &
                        num1=1,NB), num2=1,NB)
              end do
           end do
        end if
        call ptime('MMN write')

        ! Binary sidecar (B4): same 'MMNB' stream layout hybrid_files_m
        ! reads, so the combine/hybrid stages skip the multi-minute text
        ! parse entirely.  The reader convention is mmn(m,n,b,k) =
        ! conjg(overlap(n,m,pair)) -- a conjugate TRANSPOSE of each block,
        ! matching the text element order (NUM2 outer, NUM1 fastest).
        ! Close the text unit FIRST or its later flush would out-date the
        ! sidecar's mtime and permanently disable it.
        inquire(unit=unit_mmn, name=fname)
        call close_checked(unit_mmn, 'the .mmn text file')
        open(newunit=u_sc, file=trim(fname)//'.bin', form='unformatted', &
             access='stream', status='replace', action='write', iostat=ios)
        if (ios == 0) then
           allocate(ngh_sc(Nntot,Nk), gsh_sc(3,Nntot,Nk), Mt_sc(Nb,Nb))
           do k1 = 1, Nk
              do k2 = 1, Nntot
                 pair_index = (k1-1)*Nntot + k2
                 ngh_sc(k2,k1)   = KPB(pair_index)
                 gsh_sc(1,k2,k1) = BQX1(pair_index)
                 gsh_sc(2,k2,k1) = BQY1(pair_index)
                 gsh_sc(3,k2,k1) = BQZ1(pair_index)
              end do
           end do
           write(u_sc) 'MMNB', 1, Nb, Nk, Nntot
           write(u_sc) ngh_sc, gsh_sc
           do k1 = 1, Nk
              do k2 = 1, Nntot
                 pair_index = (k1-1)*Nntot + k2
                 do NUM1 = 1, Nb
                    do NUM2 = 1, Nb
                       Mt_sc(NUM2,NUM1) = conjg(overlap(num1,num2,pair_index))
                    end do
                 end do
                 write(u_sc) Mt_sc
              end do
           end do
           call close_checked(u_sc, 'the .mmn binary sidecar')
           deallocate(ngh_sc, gsh_sc, Mt_sc)
           call ptime('MMN sidecar')
        end if
     end if
  endif

  ! sigma loop: root writes the combined block-structured amn (this is
  ! exactly the doubled-inwf + w2waddsp file the SOC workflow produced,
  ! without the zero-coefficient trick or the extra tool)
  if (sigma2 .and. inwf%Amn .and. is_root) then
     do kkk = 1, Nk
        do j = 1, 2*inwf%Nproj
           do i = 1, Nb
              write(unit_amn, '(2i4,1x,i5,1x,2es24.15)') &
                   i, j, kkk, proj_h(i, j, kkk)
           end do
        end do
     end do
     inquire(unit=unit_amn, name=fname)
     call close_checked(unit_amn, 'the combined .amn text file')
     open(newunit=u_sc, file=trim(fname)//'.bin', form='unformatted', &
          access='stream', status='replace', action='write', iostat=ios)
     if (ios == 0) then
        write(u_sc) 'AMNB', 1, Nb, Nk, 2*inwf%Nproj
        write(u_sc) proj_h
        call close_checked(u_sc, 'the combined .amn binary sidecar')
     end if
     call ptime('AMN combined write')
  end if

  ! MODE=HYBR: run the validated Stage-2 in memory on root (the overlap
  ! reduce at the Mmn stage is to root, so only root holds the full sum).
  if (inwf%hybrid .and. is_root) then
     call hybrid_stage2()
     call ptime('hybrid Stage-2')
  end if

  call cleanup_gaunt_cache()  ! Free Gaunt coefficient cache
  call cleanup_loabc()        ! Free LO overlap arrays
  ! The .eig is written by read_vec's first pass and never closed before; close it with the size check.
  if (is_root) call close_checked(unit_eig, 'the .eig file')
  if (is_root) call ERRCLR(ERRFN)
  if (is_root) print "('W2W END')"

  call mpi_w2w_finalize()

contains

  !------------------------------------------------------------- Gate A tools
  !> The cut between band `b' and band `b+1', tested PER k on the DFT
  !> eigenvalues alone -- no model, no projection, nothing a construction of
  !> ours could guarantee.  This is the only right form of the test: cuts are
  !> legal only at multiples of the enforced multiplicity AND not every
  !> multiple is legal, so arithmetic on the band index is wrong in both
  !> directions (Rh2B scalar 8/12/20 and SOC 4/28-44/52 are illegal multiples;
  !> here MgB2's ODD cut after 21 bands is legal while 8, 12, 20, 24 are not).
  !>
  !> `nkx' counts the k at which band b+1 EXISTS -- which is NOT Nk.  NEn
  !> varies with k; eall is initialised to huge() and filled only for
  !> i <= NEn, so a window ending at minNEn still truncates the Hilbert space
  !> at every k that carries more bands.  The pre-W8 diagnostic skipped those
  !> k entirely (`if (inwf%bmax < minNEn)'), which is why Fe [5,19] and MgB2
  !> [1,32] both reported "degeneracy-closed" while cutting a stuck pair at
  !> 1 of 506 and 2 of 511 k respectively (measured 2026-08-26).
  subroutine edge_gap_stat(b, ncut, gmin, nkx)
    integer,  intent(in)  :: b
    integer,  intent(out) :: ncut, nkx
    real(R8), intent(out) :: gmin
    real(R8), parameter   :: EXISTS = 0.5_R8*huge(1.0_R8)
    integer  :: kx
    real(R8) :: g
    ncut = 0; nkx = 0; gmin = huge(1.0_R8)
    if (b < 1 .or. b >= NEmax) return        ! no band beyond this edge at any k
    do kx = 1, Nk
       if (eall(b,   kx) >= EXISTS) cycle    ! band b absent at this k
       if (eall(b+1, kx) >= EXISTS) cycle    ! nothing above the edge here
       nkx = nkx + 1
       g = abs(eall(b+1,kx) - eall(b,kx))
       gmin = min(gmin, g)
       if (g < DEGTOL) ncut = ncut + 1
    end do
  end subroutine edge_gap_stat

  !> One pool edge, reported and counted.  `b' is the ABSOLUTE band index
  !> below the cut, so the cut is b/b+1.  An edge with no band beyond it at
  !> any k is not a truncation and is not counted.
  subroutine pool_edge_report(b, side, ncut_tot, nedge_tot)
    integer,          intent(in)    :: b
    character(len=*), intent(in)    :: side
    integer,          intent(inout) :: ncut_tot, nedge_tot
    integer  :: nc, nkx
    real(R8) :: gm
    call edge_gap_stat(b, nc, gm, nkx)
    if (nkx == 0) return                     ! nothing beyond: not a cut
    nedge_tot = nedge_tot + 1
    if (nc > 0) then
       ncut_tot = ncut_tot + 1
       write(unit_out, '("   pool ", A, " edge: bands ", I0, "/", I0, &
            & " degenerate at ", I0, " of ", I0, " k  (min gap ", ES9.2, &
            & " Ry)")') trim(side), b, b+1, nc, nkx, gm
    end if
  end subroutine pool_edge_report

  !> Largest run of consecutive bands MUTUALLY within DEGTOL at a single k,
  !> over the whole energy scan.  Used ONLY to cap the Gate A snap walk:
  !> closing a cut through one multiplet of size m costs at most m-1 extra
  !> bands, so the largest multiplet actually present is the honest bound.
  !> Measured, per run, from this run's own eigenvalues -- never typed.
  !>
  !> *** MUTUAL, NOT A CHAIN.  CORRECTED 2026-08-26 (R4 finalize). ***
  !> The first implementation tested the ADJACENT gap |E(jb) - E(jb-1)| and
  !> so measured a transitive CHAIN, which is >= the mutual multiplet and is
  !> unbounded: a ladder of adjacent steps just under DEGTOL reads as one
  !> giant multiplet however far apart its endpoints are.  MEASURED on a
  !> synthetic MgB2 file with bands 25..31 laddered at 9e-7 Ry (endpoints
  !> 5.4e-6 Ry apart = 5.4x DEGTOL, so NOT a multiplet by any definition):
  !> the chain form returned cap 8, and `edges = snap' then grew the window
  !> 24 -> 31 and declared it closed.  That is verbatim the failure the cap
  !> exists to prevent -- and the inflation is driven by the very band
  !> density (dense SOC metals, where 24% of BETWEEN-pair cuts are already
  !> below DEGTOL) that makes the cap load-bearing.  Fail-open, in the one
  !> place that had to fail closed.
  !>
  !> The fix is a two-pointer sliding window over the sorted eigenvalues at
  !> each k: E is ascending, so `mutually within DEGTOL' for a contiguous
  !> run is exactly `last - first < DEGTOL', and the widest such window is
  !> the true maximum.  A greedy left-to-right partition (restart at the
  !> first band that leaves the current run) is NOT equivalent -- it
  !> under-counts, e.g. 0, 0.9, 1.0, 1.1 at DEGTOL 1 gives 2 where the true
  !> answer is 3 -- so the window is slid, not restarted.  Still O(NEmax)
  !> per k.  Measured after the fix: MgB2 2, Fe 3, GaAs 3 (unchanged on the
  !> real cells, where chain and mutual coincide), synthetic ladder 8 -> 2.
  integer function max_multiplet() result(m)
    real(R8), parameter :: EXISTS = 0.5_R8*huge(1.0_R8)
    integer :: kx, lo, hi, nb_k
    m = 1
    do kx = 1, Nk
       ! bands actually present at this k (eall is padded with huge() above)
       nb_k = 0
       do hi = 1, NEmax
          if (eall(hi,kx) >= EXISTS) exit
          nb_k = hi
       end do
       lo = 1
       do hi = 1, nb_k
          do while (eall(hi,kx) - eall(lo,kx) >= DEGTOL)
             lo = lo + 1
          end do
          m = max(m, hi - lo + 1)
       end do
    end do
  end function max_multiplet

  !! MODE=HYBR Stage-2, in memory (mirrors the standalone w2w_hybrid tool):
  !! P from |amn|^2 -> build_masks -> pool restriction -> frozen-constrained
  !! disentangle -> rotate to projected-H eigenbasis -> SCDM gauge ->
  !! isolated Nb=Nw hand-off <case>-hyb.{amn,mmn,eig} for stock wannier90.
  !! Data mapped straight from w2w memory: amn = l2Amn's reduced projection;
  !! mmn(m,n,b,k) = conjg(overlap(n,m,(k-1)*Nntot+b)) (conjugate TRANSPOSE
  !! per block -- the same identity the .mmn writer encodes); neigh/gshift
  !! from the UNrotated nnkp pair arrays; eig from eall, Fermi-shifted to eV.
  subroutine hybrid_stage2()
    complex(C16), parameter :: H_CONE = (1.0_R8,0.0_R8), H_CZERO = (0.0_R8,0.0_R8)
    integer :: nw, npool, nproj_h, kh, ibv, ih, jh, pairh, nkq, nnt2, froz_dmd
    integer :: fl_add_h, fl_drop_h, fl_nk_h
    real(R8) :: fe_used_h
    logical  :: froz_legacy_h
    ! the disentanglement outcome and the subspace continuity, for the log and <outp>.dis
    character(len=320) :: dis_status_h, cont_line_h
    real(R8) :: cont_smin_h
    integer  :: cont_k_h, cont_kb_h
    real(R8) :: omega_h, wres_h, svmin_h
    real(R8),     allocatable :: P_h(:,:), eig_h(:,:), bcart_h(:,:), wgt_h(:)
    real(R8),     allocatable :: bkb_h(:,:,:)   ! per-(k,b) b-vectors (nnkp_block_bvectors)
    real(R8),     allocatable :: eig_p(:,:), eig_t(:,:)
    complex(C16), allocatable :: mmn_p(:,:,:,:), amn_p(:,:,:)
    complex(C16), allocatable :: Vh(:,:,:), Wh(:,:,:), mmn_rot(:,:,:,:)
    complex(C16), allocatable :: Bm_h(:,:,:), Ag(:,:,:)
    integer,      allocatable :: neigh_h(:,:), gshift_h(:,:,:), pool_h(:), anch_h(:)
    logical,      allocatable :: admit_h(:,:), frozen_h(:,:)
    character(len=16)    :: regime_out
    character(len=BUFSZ) :: outp, casenm
    ! anatomy of the mask rule, for <case>-hyb.masks / <case>-hyb.anatomy
    real(R8), allocatable :: P_used(:,:), an_avg(:), an_min(:), an_bot(:)
    real(R8) :: an_plateau, an_fid, an_pgate, an_pfroz
    logical  :: an_ok

    nw = inwf%num_wann
    if (.not. allocated(proj_h)) call croak('HYBR: no amn in memory')
    nproj_h = size(proj_h, 2)      ! = 2*Nproj under the sigma loop
    if (nw > Nb) call croak('HYBR: num_wann exceeds the band window')
    if (nw > nproj_h) call croak('HYBR: num_wann exceeds the projection count &
         &(SCDM needs nproj >= num_wann; raise max_proj or the threshold)')

    ! window-relative Fermi-shifted eigenvalues (eV)
    allocate(P_h(Nb,Nk), eig_h(Nb,Nk))
    do kh = 1, Nk
       do ih = 1, Nb
          eig_h(ih,kh) = (eall(inwf%bmin-1+ih,kh) - efermi)*Ryd_eV
       end do
    end do

    ! Loewdin per-state projectability from the FULL-CELL (tails) amn on
    ! the PHYSICAL columns -- the validated lowTLadm P source, computed
    ! in memory (no W2W_PAMN sidecar file).  Physical columns = leading
    ! mask_cols of each sigma block (the Otsu mask boundary exported by
    ! the AUTO selection); no boundary -> all columns (the validated
    ! full-P default).  Replaces the raw |amn|^2 reconstruction.
    block
      logical, allocatable :: pcol(:)
      integer :: nspat, mc
      if (.not. allocated(projT_h)) call croak('HYBR: no full-cell amn in memory')
      allocate(pcol(nproj_h))
      nspat = inwf%Nproj              ! spatial columns per sigma block
      mc = inwf%auto_config%mask_cols
      if (mc > 0 .and. mc < nspat) then
         pcol = .false.
         pcol(1:mc) = .true.
         if (nproj_h == 2*nspat) pcol(nspat+1:nspat+mc) = .true.
         write(unit_out,'(" HYBR: P from ", I0, " physical of ", I0, &
              &  " spatial columns", A)') mc, nspat, &
              merge(' x 2 sigma blocks', '                 ', nproj_h == 2*nspat)
      else
         pcol = .true.
         write(unit_out,'(" HYBR: P from all ", I0, " amn columns")') nproj_h
      end if
      call lowdin_pstate(projT_h, Nb, Nk, nproj_h, pcol, P_h)
      write(unit_out,'(" HYBR: Lowdin P (full-cell): max = ", F8.4, &
           &  "   mean = ", F8.4)') maxval(P_h), sum(P_h)/real(Nb*Nk,R8)
    end block

    ! b-vectors from the nnkp file (same source/units as the standalone tool)
    call read_nnkp(trim(nnkpfn), nkq, nnt2, bcart_h)
    if (nkq /= Nk .or. nnt2 /= Nntot) call croak('HYBR: nnkp dims mismatch')
    allocate(wgt_h(Nntot))
    call mv_weights(bcart_h, Nntot, wgt_h, wres_h)
    write(unit_out,'(" HYBR: MV b-weight residual = ", ES10.2)') wres_h

    ! masks from the Loewdin P (already in [0,1]: rescale='none', the
    ! same convention the standalone uses after its lowdin path).
    ! Validated-recipe defaults: admit-all (per-state admission disabled
    ! -- lowTLadm); froz_emax + pool factor from the HYBR keyword tail.
    allocate(P_used(Nb,Nk), an_avg(Nb), an_min(Nb), an_bot(Nb))
    call build_masks(P_h, eig_h, Nb, Nk, nw, 0.0_R8, 'none', &
         &           trim(inwf%hyb_regime), trim(inwf%hyb_froz_mode), &
         &           inwf%hyb_froz_val, pool_h, npool, &
         &           admit_h, frozen_h, regime_out, &
         &           pool_factor_in=inwf%hyb_pool_factor, &
         &           froz_demand_out=froz_dmd, fe_used_out=fe_used_h, &
         &           froz_legacy_out=froz_legacy_h, floor_add_out=fl_add_h, &
         &           floor_drop_out=fl_drop_h, floor_nk_out=fl_nk_h, &
         &           P_used_out=P_used, anat_avg=an_avg, anat_min=an_min, &
         &           anat_bot=an_bot, anat_plateau=an_plateau, anat_fid=an_fid, &
         &           anat_pgate=an_pgate, anat_pfroz=an_pfroz, anat_ok=an_ok)
    admit_h = .true.
    write(unit_out,'(" HYBR: admit-all (validated recipe default)")')
    write(unit_out,'(" HYBR: regime=", A, "  pool[", I0, "..", I0, &
         &  "] npool=", I0, "  frozen/k mean=", F6.2, "  demand(worst k)=", &
         &  I0, "  (capacity: nw >= demand)")') trim(regime_out), &
         pool_h(1), pool_h(npool), npool, &
         real(sum(count(frozen_h,1)),R8)/real(Nk,R8), froz_dmd
    ! the authoritative statement of what the optimiser is held to (S1a);
    ! same wording as the standalone tool, one shared implementation
    call froz_report_line(unit_out, frozen_h, eig_h, pool_h, npool, Nk, nw, &
         0.0_R8, fe_used_h, trim(inwf%hyb_froz_mode), froz_dmd, .false., &
         froz_legacy_h, fl_add_h, fl_drop_h, fl_nk_h)
    ! capacity assertion on the PRE-REPAIR demand (the per-k feasibility
    ! repair inside the mask builder silently demotes frozen states above
    ! nw -- the Rh3Si4 silent failure -- so the returned mask can never
    ! flag the violation itself).  Also covers froz_emax = auto, where
    ! the ceiling exists only after the walk-down.
    if (froz_dmd > nw) then
       write(unit_out,'(" HYBR: FROZEN CAPACITY VIOLATION: worst-k demand ", &
            & I0, " > num_wann ", I0, " -> raise num_wann to >= ", I0, &
            & " or lower froz_emax")') froz_dmd, nw, froz_dmd
       if (inwf%hyb_cap_strict) then
          call croak('HYBR: frozen manifold exceeds num_wann (capacity; &
               &`capacity = warn'' overrides)')
       else
          write(unit_out,'(" HYBR: proceeding under capacity (capacity = &
               &warn); frozen budget spent deepest-first, ceiling states &
               &released to the disentangler -- check eta")')
       end if
    end if

    !--------------------------------------------------------------- GATE B
    ! (2/2) POST-build_masks ENFORCEMENT.  select_pool sizes the pool as
    ! max(nw+1, ceil(pool_factor*nw)) and then CLAMPS it against the top of
    ! the window without saying so, so the only way npool can fall below the
    ! request is that silent clamp.  Enforced in the CALLER: hybrid_prep.F
    ! stays print-free because the standalone shares it.
    !
    ! *** DUPLICATED FORMULA -- KEEP IN STEP WITH hybrid_prep.F select_pool
    ! (the `n_pool_max = max(nw+1, ceiling(pool_factor*real(nw, dp)))' line,
    ! cross-referenced there).  Verified EXACT on 2026-08-26: select_pool
    ! computes npool = max(n_pool_max, seed span), then clamps s0 into
    ! 1..nb-npool+1, then npool = min(npool, nb-s0+1); since the seed span
    ! is always inside 1..nb, the ONLY path that lowers npool below the
    ! request is that final clamp, which is exactly `npool < npool_want'.
    ! This is a maintenance hazard, not a live defect: the two expressions
    ! live in files owned by different lanes and nothing ties them, so if
    ! select_pool's sizing ever changes, Gate B will silently report a ratio
    ! against a target the pool builder no longer uses.
    poolgate: block
      integer :: npool_want, bmax_need
      npool_want = max(nw+1, ceiling(inwf%hyb_pool_factor*real(nw,R8)))
      if (npool < npool_want) then
         bmax_need = inwf%bmin + pool_h(1) - 1 + npool_want - 1
         write(unit_out,'(" HYBR pool: RATIO VIOLATION: npool ", I0, " < ", &
              & I0, " = pool_factor*nw  (measured ratio ", F6.2, " < ", &
              & F6.2, ")  -- select_pool clamped against the window top; &
              &raise bmax to >= ", I0, " (or lower pool_factor)")') &
              npool, npool_want, real(npool,R8)/real(nw,R8), &
              inwf%hyb_pool_factor, bmax_need
         if (inwf%hyb_cap_strict) then
            call croak('HYBR: pool ratio below pool_factor (see the HYBR &
                 &pool line; `capacity = warn'' overrides)')
         else
            write(unit_out,'(" HYBR pool: PROCEEDING under the ratio &
                 &(capacity = warn); disentanglement slack is limited -- &
                 &check eta")')
         end if
      else
         write(unit_out,'(" HYBR pool: npool ", I0, " >= ", I0, &
              & " = pool_factor*nw  (measured ratio ", F6.2, ")  OK")') &
              npool, npool_want, real(npool,R8)/real(nw,R8)
      end if
    end block poolgate

    !--------------------------------------------------------------- GATE A
    ! (2/2) THE nw TRUNCATION CUT.  The outer band window is not the only
    ! truncation: the POOL is a second one, and the disentangler never sees a
    ! band outside it, so a pool edge inside a stuck multiplet is the same
    ! defect one level down.  It matters more, not less, than the window
    ! edge, because select_pool sizes the pool from pool_factor and knows
    ! nothing about multiplets -- nobody typed this edge and until now nobody
    ! looked at it.  Tested with the SAME per-k gap test on the SAME DFT
    ! eigenvalues, on the pool's ABSOLUTE band indices.
    pooledge: block
      logical  :: inpool(Nb)
      integer  :: ib, nabs, nc1, nkx1, ncut_pool, nedge_pool
      integer  :: ptop, pbot, psug, pnew
      real(R8) :: gm1
      inpool = .false.
      do ih = 1, npool
         inpool(pool_h(ih)) = .true.
      end do
      ncut_pool = 0; nedge_pool = 0
      do ib = 1, Nb
         if (.not. inpool(ib)) cycle
         nabs = inwf%bmin + ib - 1
         ! a cut BELOW this band: it is a pool edge iff the band under it is
         ! outside the pool (or outside the window)
         if (ib == 1) then
            call pool_edge_report(nabs-1, 'bottom', ncut_pool, nedge_pool)
         else if (.not. inpool(ib-1)) then
            call pool_edge_report(nabs-1, 'bottom', ncut_pool, nedge_pool)
         end if
         ! a cut ABOVE this band
         if (ib == Nb) then
            call pool_edge_report(nabs, 'top   ', ncut_pool, nedge_pool)
         else if (.not. inpool(ib+1)) then
            call pool_edge_report(nabs, 'top   ', ncut_pool, nedge_pool)
         end if
      end do
      if (ncut_pool == 0) then
         write(unit_out,'(" HYBR pool: all ", I0, " pool edge(s) are &
              &degeneracy-closed")') nedge_pool
      else
         write(unit_out,'(" HYBR pool: WARNING ", I0, " of ", I0, " pool &
              &edge(s) CUT a degenerate multiplet (see above)")') &
              ncut_pool, nedge_pool

         ! ---- RANKED ALTERNATIVE.  Added 2026-08-26 (R4 finalize) so this
         ! half can fail closed at all: a gate that stops with no way
         ! forward is not autonomous.  The knob is pool_factor (it sets
         ! npool = max(nw+1, ceil(pool_factor*nw))), so walk the top pool
         ! edge up to the nearest closed cut and report the npool and the
         ! pool_factor that reach it -- and whether bmax has to rise too.
         ptop = inwf%bmin + maxval(pool_h(1:npool)) - 1
         pbot = inwf%bmin + minval(pool_h(1:npool)) - 1
         psug = 0
         do ib = ptop+1, minNEn
            call edge_gap_stat(ib, nc1, gm1, nkx1)
            if (nc1 == 0) then; psug = ib; exit; end if
         end do
         if (psug > 0) then
            pnew = npool + (psug - ptop)
            write(unit_out,'("   nearest closed pool TOP edge: band ", I0, &
                 & " -> npool ", I0, " (pool_factor >= ", F6.3, ")")') &
                 psug, pnew, real(pnew,R8)/real(nw,R8)
            if (psug > inwf%bmax) write(unit_out,'("   ...and that needs &
                 &bmax >= ", I0, " (currently ", I0, ")")') psug, inwf%bmax
         else
            write(unit_out,'("   no closed pool TOP edge exists at or below &
                 &minNEn = ", I0)') minNEn
         end if
         psug = 0
         do ib = pbot-1, 1, -1
            if (ib == 1) then; psug = 1; exit; end if
            call edge_gap_stat(ib-1, nc1, gm1, nkx1)
            if (nc1 == 0) then; psug = ib; exit; end if
         end do
         if (psug > 0 .and. psug < pbot) write(unit_out,'("   nearest closed &
              &pool BOTTOM edge: the pool would have to start at band ", I0, &
              & " (currently ", I0, ")")') psug, pbot

         ! ---- A cut through a multiplet at the pool edge is the same defect
         ! as at a window edge, one level down (the disentangler never sees a
         ! band outside the pool), so both follow the one switch eg_strict:
         ! fatal under `edges = strict' or `snap', reported otherwise.
         ! `edges = warn' never stops, and `capacity = warn' does the same
         ! for Gate B.
         if (eg_strict) call croak('Gate A: a POOL edge cuts a degenerate &
              &multiplet -- the disentangler never sees a band outside the &
              &pool, so this truncation breaks symmetry adaptation exactly &
              &as a window edge does (see the nearest closed pool edges &
              &above; `edges = warn'' overrides)')
      end if
    end block pooledge

    ! pool-restricted mmn/amn/eig straight from memory (no full-Nb copy):
    ! the conjugate transpose of each overlap block IS the w90 convention
    allocate(mmn_p(npool,npool,Nntot,Nk), amn_p(npool,nproj_h,Nk), &
         &   eig_p(npool,Nk), neigh_h(Nntot,Nk), gshift_h(3,Nntot,Nk))
    do kh = 1, Nk
       do ibv = 1, Nntot
          pairh = (kh-1)*Nntot + ibv
          neigh_h(ibv,kh)   = KPB(pairh)
          gshift_h(1,ibv,kh) = BQX1(pairh)
          gshift_h(2,ibv,kh) = BQY1(pairh)
          gshift_h(3,ibv,kh) = BQZ1(pairh)
          do jh = 1, npool
             do ih = 1, npool
                mmn_p(ih,jh,ibv,kh) = conjg(overlap(pool_h(jh), pool_h(ih), pairh))
             end do
          end do
       end do
       do ih = 1, npool
          amn_p(ih,:,kh) = proj_h(pool_h(ih),:,kh)
          eig_p(ih,kh)   = eig_h(pool_h(ih),kh)
       end do
    end do

    allocate(Vh(npool,nw,Nk), Wh(npool,nw,Nk), eig_t(nw,Nk), &
         &   mmn_rot(nw,nw,Nntot,Nk), Bm_h(nproj_h,nw,Nk), Ag(nw,nw,Nk), &
         &   anch_h(nw))
    ! Every one of these was a hardcoded literal until 2026-07-29, which left
    ! the one-shot running a configuration that was never validated: half the
    ! production sweep budget with the relaxation and the adaptive convergence
    ! cadence both off.  It hit the cap without converging on GaAs and on the
    ! Rh3Si4 SOC cell (|dOmega| 3.9e-2 after 2513 s) where the standalone
    ! converges in 142 sweeps.  Defaults now ARE the production recipe.
    if (is_root) write(unit_out,'(" HYBR disentangle: sweeps<=", I0, &
         &"  relax=", F6.2, "  adaptive=", L1, "  colored=", L1, &
         &"  anderson=", I0)') inwf%hyb_sweeps, inwf%hyb_mix, &
         inwf%hyb_adaptive, inwf%hyb_colored, inwf%hyb_anderson
    if (inwf%hyb_mix == 1.0_R8 .and. is_root) write(unit_out,'(a)') &
         '  [NOTE] relax=1 disables over-relaxation; expect ~20x more sweeps'
    ! ---- solver dispatch: Gauss-Seidel sweep, or direct Grassmann CG ----
    ! The dispatch lives HERE, not inside disentangle_frozen, because
    ! hybrid_grass_m uses hybrid_core_m for kinfo_t -- putting the branch in
    ! hybrid_core_m would make the dependency circular.

    solver_pick: block
      character(len=8) :: solv
      integer :: ndis_max, kq
      type(kinfo_t), allocatable :: kig(:)
      type(gcg_t) :: Pg
      complex(R8), allocatable :: Vg(:,:,:)
      real(R8) :: relg
      integer :: nitg
      integer :: whyg
      real(R8) :: gng, gng0

      solv = trim(inwf%hyb_solver)
      if (solv == 'auto') then
         ! ndis is what sets the crossover: the sweep is eigh-dominated (n^3)
         ! while CG is matmul + QR.  LCAO measured the crossing at ndis ~ 250;
         ! that number has NOT been reproduced on LAPW inputs, so 'auto' is a
         ! starting hypothesis, not a calibrated threshold.
         ndis_max = 0
         do kq = 1, Nk
            ndis_max = max(ndis_max, count(admit_h(:,kq) .and. .not. frozen_h(:,kq)))
         end do
         if (ndis_max >= 250) then
            solv = 'cg'
         else
            solv = 'gs'
         end if
         if (is_root) write(unit_out,'(a,i0,a,a,a)') &
              ' HYBR solver: auto -> max ndis = ', ndis_max, &
              ' so choosing ', trim(solv), ' (LCAO crossover ~250, unverified here)'
      end if

      if (solv == 'cg') then
         ! Share the frozen split; ask build_kinfo to SKIP the Mdis hoist,
         ! which only the sweep reads.
         allocate(kig(Nk))
         call build_kinfo(mmn_p, neigh_h, frozen_h, admit_h, npool, Nntot, &
              Nk, nw, kig, want_mdis=.false.)
         call gcg_build(mmn_p, neigh_h, wgt_h, kig, npool, Nntot, Nk, nw, Pg)
         allocate(Vg(Pg%dmax, Pg%amax, Nk))
         ! The SAME project-optimal SVD seed the sweep uses, so a CG-vs-GS
         ! comparison starts from one point.  V MUST be orthonormal here or
         ! the tangent projection and the retraction are both invalid.
         call gcg_seed_svd(Pg, kig, amn_p, npool, nproj_h, Vg)

         ! STARTUP GUARD, not an optional test.  The reference implementation
         ! shipped a gradient missing its neighbour term; it descended
         ! plausibly and was caught only here (FD/analytic 1.84 and 2.10).
         relg = gcg_check_gradient(Pg, Vg)
         if (is_root) write(unit_out,'(a,es10.2)') &
              ' HYBR cg: finite-difference gradient check, rel err = ', relg
         if (relg > 1.0e-4_R8) call croak('HYBR cg: the analytic gradient &
              &disagrees with finite differences. Refusing to run -- a wrong &
              &gradient still descends and would silently give a worse minimum.')

         call grassmann_cg(Pg, Vg, inwf%hyb_sweeps, inwf%hyb_tol, &
              inwf%hyb_rel_tol, omega_h, nitg, .false., why=whyg, &
              gnorm=gng, gnorm0=gng0)
         call gcg_to_full(Pg, kig, Vg, npool, nw, Vh)
         ! CROSS-CHECK against an INDEPENDENT implementation of the same
         ! quantity.  gcg_omega evaluates the frozen-split form (aconst + Q +
         ! Mdd); omega_invariant evaluates the full-frame form directly from
         ! mmn.  If they disagree, the CG objective is not the one the sweep
         ! minimizes and any apparent improvement is illusory -- so this is a
         ! guard, not a diagnostic.
         relg = omega_invariant(Vh, mmn_p, neigh_h, wgt_h, npool, Nntot, Nk, nw)
         if (is_root) write(unit_out,'(a,f18.10,a,es10.2)') &
              ' HYBR cg: independent Omega_I = ', relg, '   diff ', &
              abs(relg - omega_h)
         if (abs(relg - omega_h) > 1.0e-6_R8*max(abs(relg),1.0_R8)) &
              call croak('HYBR cg: the Grassmann objective disagrees with &
              &omega_invariant on the converged subspace -- the two are not &
              &minimizing the same functional.')
         omega_h = relg
         ! Say WHY it stopped -- same five-way reason as the standalone
         ! driver (w2w_hybrid.F); only GCG_OK is convergence.
         if (is_root) then
            select case (whyg)
            case (GCG_OK)
               write(unit_out,'(a,es10.3,a,es10.3,a,i0,a)') &
                    ' HYBR cg: CONVERGED, |grad| = ', gng, ' / |grad0| = ', &
                    gng0, ' after ', nitg, ' gradient evaluations'
            case (GCG_SMALLDEC)
               write(unit_out,'(a,i0,a)') &
                    ' HYBR cg: stopped on small decrement after ', nitg, &
                    ' gradient evaluations -- NOT converged'
            case (GCG_STALL)
               write(unit_out,'(a,i0,a)') &
                    ' HYBR cg: LINE SEARCH STALLED after ', nitg, &
                    ' gradient evaluations -- NOT converged'
            case (GCG_TINYSTEP)
               write(unit_out,'(a,i0,a)') &
                    ' HYBR cg: step underflow after ', nitg, &
                    ' gradient evaluations -- NOT converged'
            case default
               write(unit_out,'(a,i0,a)') &
                    ' HYBR cg: HIT THE ITERATION CAP after ', nitg, &
                    ' gradient evaluations -- NOT converged'
            end select
         end if
         if (whyg == GCG_OK) then
            write(dis_status_h,'(a,i0,a,f0.6,a)') 'converged (cg, ', nitg, &
                 ' gradient evaluations, Omega_I ', omega_h, ' Ang^2)'
         else
            write(dis_status_h,'(a,a,a,i0,a,es8.2,a,f0.6,a)') 'NOT certified (cg exit ', &
                 trim(merge('SMALLDEC', merge('STALL   ', merge('TINYSTEP', 'CAP     ', &
                 whyg == GCG_TINYSTEP), whyg == GCG_STALL), whyg == GCG_SMALLDEC)), ', ', nitg, &
                 ' gradient evaluations, |grad| ', gng, ', Omega_I ', omega_h, ' Ang^2)'
         end if
         call gcg_free(Pg)
         deallocate(Vg, kig)
      else
         block
           integer :: nsw_gs
           logical :: conv_gs
           character(len=320) :: capv_gs
           call disentangle_frozen(mmn_p, amn_p, neigh_h, wgt_h, npool, nproj_h, &
                Nntot, Nk, nw, admit_h, frozen_h, Vh, omega_h, inwf%hyb_sweeps, &
                inwf%hyb_tol, inwf%hyb_rel_tol, 10, .false., &
                colored=inwf%hyb_colored, &
                adaptive=inwf%hyb_adaptive, mixing=inwf%hyb_mix, &
                anderson=inwf%hyb_anderson, kick=inwf%hyb_kick, &
                taper=inwf%hyb_taper, relax_late=inwf%hyb_relax_late, &
                taper_late=inwf%hyb_taper_late, &
                nsweeps=nsw_gs, is_converged=conv_gs, cap_verdict=capv_gs)
           if (conv_gs) then
              write(dis_status_h,'(a,i0,a,f0.6,a)') 'converged (gs, ', nsw_gs, &
                   ' sweeps, Omega_I ', omega_h, ' Ang^2)'
           else
              write(dis_status_h,'(a,i0,a,f0.6,a)') 'stopped at the sweep cap (', nsw_gs, &
                   ' sweeps, Omega_I ', omega_h, ' Ang^2): '
              dis_status_h = trim(dis_status_h)//' '//trim(capv_gs)
              write(unit_out,'(a)') ' HYBR: '//trim(dis_status_h)
           end if
         end block
      end if
    end block solver_pick
    write(unit_out,'(" HYBR: Omega_I = ", F16.8)') omega_h
    ! Omega_I is a spread: sum_kb w_b (nw - ||V^H M V||_F^2) with each block
    ! contractive, so it is mathematically NON-NEGATIVE.  A negative value
    ! means the inputs are not what the code thinks they are.
    !
    ! The way this actually happens: running the REAL binary (w2w) on a
    ! COMPLEX case.  Nothing detects that -- real vs complex is a COMPILE-TIME
    ! choice (-D_COMPLEX_) selected by the caller via `x w2w -c`, not a
    ! property the code reads from any input file -- and the mistake is
    ! SILENT, because reading N reals out of a record holding N complex
    ! numbers is a legal short read.  It runs to completion and writes a
    ! hand-off.  Measured on GaAs: Omega_I = -127.61 instead of 21.83.
    ! (The opposite mistake, w2wc on a real case, at least dies loudly with
    ! "I/O past end of record".)
    if (omega_h < 0.0_R8) call croak('HYBR: Omega_I is NEGATIVE, which is &
         &impossible for a spread. The usual cause is running the real binary &
         &(w2w) on a complex case -- rerun with w2wc (x w2w -c). Real/complex &
         &is a compile-time choice and is NOT auto-detected.')
    call rotate_to_subspace(eig_p, Vh, mmn_p, neigh_h, npool, Nntot, Nk, &
         nw, eig_t, Wh, mmn_rot)
    ! subspace continuity (diagnostic: nothing written changes)
    call subspace_continuity(mmn_rot, neigh_h, nw, Nntot, Nk, unit_out, cont_smin_h, &
         cont_k_h, cont_kb_h)
    write(cont_line_h,'(es9.3,a,i0,a,i0)') cont_smin_h, &
         ' (smallest singular value of the hand-off overlaps), at k ', cont_k_h, ' -> ', cont_kb_h
    if (cont_smin_h < SUBSPACE_JUMP_SV) cont_line_h = trim(cont_line_h)// &
         '; WARNING: below 1e-2, the subspace jumps there'
    ! pocket skip (hybrid_prep_m::pocket_ceiling): the standalone w2w_hybrid
    ! lowers a walked ceiling under a mixed-character pocket and disentangles
    ! again; this one-shot reports the ceiling to pass as froz_emax instead.
    if (trim(inwf%hyb_froz_mode) == 'auto') then
       block
         real(R8), allocatable :: svl(:,:), epsp(:,:), Ppk(:,:)
         real(R8) :: fe_pk
         logical :: pk_f
         character(len=400) :: pk_m
         integer :: ip, kp
         allocate(svl(Nntot,Nk), epsp(npool,Nk), Ppk(npool,Nk))
         call continuity_links(mmn_rot, nw, Nntot, Nk, svl)
         do kp = 1, Nk
            do ip = 1, npool
               epsp(ip,kp) = eig_p(ip,kp)
               Ppk(ip,kp) = P_used(pool_h(ip),kp)
            end do
         end do
         call pocket_ceiling(epsp, Ppk, npool, Nk, fe_used_h, svl, neigh_h, Nntot, &
              fe_pk, pk_f, pk_m)
         if (pk_f) then
            write(unit_out,'(a,es10.3,a)') ' HYBR: pocket: the weakest link ', cont_smin_h, &
                 ' touches'//trim(pk_m)
            write(unit_out,'(a,f8.3,a)') ' HYBR: pocket: rerun with froz_emax = ', fe_pk, &
                 ' eV (w2w_hybrid --froz-emax auto applies this automatically)'
         end if
         deallocate(svl, epsp, Ppk)
       end block
    end if
    do kh = 1, Nk
       call zgemm('C','N', nproj_h, nw, npool, H_CONE, amn_p(1,1,kh), &
            npool, Wh(1,1,kh), npool, H_CZERO, Bm_h(1,1,kh), nproj_h)
    end do
    ! hand-off name first (choose_gauge writes <outp>.anchors)
    casenm = nnkpfn
    ih = index(casenm, '.nnkp', back=.true.)
    if (ih > 0) casenm = casenm(1:ih-1)
    outp = trim(casenm)//'-hyb'
    ! the masks and the anatomy behind them (tools/plot_certificate.py,
    ! tools/plot_anatomy.py); admit_h is the admit-all mask applied above
    call write_masks(trim(outp)//'.masks', pool_h, npool, frozen_h, admit_h, &
         P_used, eig_h, 0.0_R8, Nb, Nk, nw, an_pfroz, an_pgate, fe_used_h)
    call write_anatomy(trim(outp)//'.anatomy', pool_h, npool, an_avg, an_min, &
         an_bot, Nb, nw, an_pgate, an_pfroz, an_plateau, an_fid, an_ok)
    write(unit_out,'(" HYBR: wrote ", A, ".masks and ", A, ".anatomy")') &
         trim(outp), trim(outp)

    ! gauge via the SHARED dispatcher (identical code to the standalone):
    ! mode + seed policy from the HYBR keyword tail (defaults auto8 +
    ! opt16 + smooth = the validated recipe); W2W_SEEDK/W2W_SEEDLOSS
    ! env vars still override for experiments.  Narration -> unit_out.
    write(unit_out,'(" HYBR: gauge=", A, "  seedk=", A, "  seedloss=", A)') &
         trim(inwf%hyb_gauge), trim(inwf%hyb_seedk), trim(inwf%hyb_seedloss)
    ! the b-vector of every (k,b) block, from the same .nnkp: each k lists its
    ! neighbours in its own order, so bcart_h(:,b) is right only at k = 1
    allocate(bkb_h(3, Nntot, Nk))
    block
      character(len=256) :: bvmsg
      logical :: bvok
      call nnkp_block_bvectors(trim(nnkpfn), Nk, Nntot, neigh_h, gshift_h, bcart_h, bkb_h, &
                               bvok, bvmsg)
      if (.not. bvok) call croak('HYBR: per-block b-vectors: '//trim(bvmsg))
    end block
    call choose_gauge(Bm_h, eig_t, mmn_rot, neigh_h, wgt_h, bcart_h, &
         nproj_h, nw, Nntot, Nk, trim(inwf%hyb_gauge), trim(outp), Ag, &
         anch_h, svmin_h, uout=unit_out, seedk_in=trim(inwf%hyb_seedk), &
         seedloss_in=trim(inwf%hyb_seedloss), bkb=bkb_h)
    call write_handoff(trim(outp), trim(casenm), eig_t, Ag, mmn_rot, &
         neigh_h, gshift_h, 0.0_R8, nw, Nk, Nntot)
    ! acceptance runs INSIDE write_handoff, before anything is written:
    ! full coverage over every (k,b) block, error stop on violation.
    write(unit_out,'(" HYBR: wrote ", A, ".{amn,mmn,eig}")') trim(outp)
    call write_dis_status(trim(outp)//'.dis', dis_status_h, cont_line_h)
    write(unit_out,'(a)') ' HYBR: disentanglement: '//trim(dis_status_h)
  end subroutine hybrid_stage2

end program wf


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
