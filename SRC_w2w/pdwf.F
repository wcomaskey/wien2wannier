!!! wien2wannier/SRC_w2w/pdwf.F
!!!
!!!    Projectability-Disentangled Wannier Functions (PDWF) module
!!!
!!!    Implements automated band selection and projection generation
!!!    based on the PDWF method (Qiao, Pizzi, Marzari, npj Comput.
!!!    Mater. 9, 208, 2023), adapted for the LAPW basis of Wien2k.
!!!
!!!    The method:
!!!    1. Selects target orbitals from a valence configuration table
!!!    2. Computes per-band projectability using Alm coefficients
!!!    3. Classifies bands into frozen/disentangle/excluded
!!!    4. Determines energy windows automatically
!!!    5. Writes AMN projection matrix for Wannier90

module pdwf_m
  use const,     only: R8, C16, DPk, TAU
  use param,     only: Lmax2, Nrf, LOmax, Nloat, Nrad, &
       &               unit_vsp, unit_out, unit_vector
  use w2w,       only: Nmat, iblock, unit_amn
  use structmod, only: struct_t
  use inwfmod,   only: inwf_t

  implicit none
  private

  public :: pdwf_target_t, build_targets, pdwf_compute

  integer, parameter :: BAND_FROZEN      = 1
  integer, parameter :: BAND_DISENTANGLE = 2
  integer, parameter :: BAND_EXCLUDED    = 3

  type pdwf_target_t
     integer :: Ntarget   ! total number of target orbitals (= num_wann)
     integer, allocatable :: atom(:)   ! (Ntarget) atom index (1..Nat)
     integer, allocatable :: jatom(:)  ! (Ntarget) inequivalent atom index
     integer, allocatable :: mu(:)     ! (Ntarget) equivalent atom index within jatom
     integer, allocatable :: l(:)      ! (Ntarget) angular momentum
     integer, allocatable :: m(:)      ! (Ntarget) magnetic quantum number (-l..l)
     integer, allocatable :: lm(:)     ! (Ntarget) compound index l*(l+1)+m+1
  end type pdwf_target_t

contains

!-----------------------------------------------------------------------
! Valence configuration lookup: returns target l-channels for element Z
!-----------------------------------------------------------------------
subroutine valence_config(Z, l_channels, n_channels)
  integer, intent(in)  :: Z
  integer, intent(out) :: l_channels(0:3)  ! 1 if channel is target, 0 otherwise
  integer, intent(out) :: n_channels       ! number of active channels

  l_channels = 0
  n_channels = 0

  select case(Z)
  ! Period 1
  case(1:2)     ! H, He
     l_channels(0) = 1; n_channels = 1
  ! Period 2
  case(3:4)     ! Li, Be
     l_channels(0) = 1; n_channels = 1
  case(5:10)    ! B-Ne
     l_channels(0:1) = 1; n_channels = 2
  ! Period 3
  case(11)      ! Na
     l_channels(0) = 1; n_channels = 1
  case(12:18)   ! Mg-Ar
     l_channels(0:1) = 1; n_channels = 2
  ! Period 4
  case(19)      ! K
     l_channels(0) = 1; n_channels = 1
  case(20)      ! Ca
     l_channels(0:1) = 1; n_channels = 2
  case(21:30)   ! Sc-Zn (3d transition metals)
     l_channels(0) = 1; l_channels(2) = 1; n_channels = 2
  case(31:36)   ! Ga-Kr
     l_channels(0:1) = 1; n_channels = 2
  ! Period 5
  case(37)      ! Rb
     l_channels(0) = 1; n_channels = 1
  case(38)      ! Sr
     l_channels(0:1) = 1; n_channels = 2
  case(39:48)   ! Y-Cd (4d transition metals)
     l_channels(0) = 1; l_channels(2) = 1; n_channels = 2
  case(49:54)   ! In-Xe
     l_channels(0:1) = 1; n_channels = 2
  ! Period 6
  case(55)      ! Cs
     l_channels(0) = 1; n_channels = 1
  case(56)      ! Ba
     l_channels(0:1) = 1; n_channels = 2
  case(57)      ! La
     l_channels(0:2) = 1; n_channels = 3
  case(58:71)   ! Ce-Lu (lanthanides)
     l_channels(0) = 1; l_channels(2:3) = 1; n_channels = 3
  case(72:80)   ! Hf-Hg (5d transition metals)
     l_channels(0) = 1; l_channels(2) = 1; n_channels = 2
  case(81:86)   ! Tl-Rn
     l_channels(0:1) = 1; n_channels = 2
  ! Period 7
  case(87)      ! Fr
     l_channels(0) = 1; n_channels = 1
  case(88)      ! Ra
     l_channels(0:1) = 1; n_channels = 2
  case(89:103)  ! Ac-Lr (actinides)
     l_channels(0) = 1; l_channels(2:3) = 1; n_channels = 3
  case default
     l_channels(0:1) = 1; n_channels = 2
  end select
end subroutine valence_config


!-----------------------------------------------------------------------
! Build target orbital list from structure and valence configuration
!-----------------------------------------------------------------------
subroutine build_targets(stru, targets)
  type(struct_t),     intent(in)  :: stru
  type(pdwf_target_t), intent(out) :: targets

  integer :: jatom, iZ, l, m, mu_eq, latom, lfirst
  integer :: l_channels(0:3), n_channels
  integer :: Ntarget, idx

  ! First pass: count total target orbitals
  Ntarget = 0
  do jatom = 1, stru%Nneq
     iZ = nint(stru%Z(jatom))
     call valence_config(iZ, l_channels, n_channels)
     do l = 0, min(Lmax2, 3)
        if (l_channels(l) == 1) then
           Ntarget = Ntarget + (2*l+1) * stru%mult(jatom)
        end if
     end do
  end do

  targets%Ntarget = Ntarget
  allocate(targets%atom(Ntarget), targets%jatom(Ntarget), &
       &   targets%mu(Ntarget), targets%l(Ntarget), &
       &   targets%m(Ntarget), targets%lm(Ntarget))

  ! Second pass: fill target arrays
  idx = 0
  lfirst = 1
  do jatom = 1, stru%Nneq
     iZ = nint(stru%Z(jatom))
     call valence_config(iZ, l_channels, n_channels)

     do mu_eq = 1, stru%mult(jatom)
        latom = lfirst - 1 + mu_eq
        do l = 0, min(Lmax2, 3)
           if (l_channels(l) /= 1) cycle
           do m = -l, l
              idx = idx + 1
              targets%atom(idx)  = latom
              targets%jatom(idx) = jatom
              targets%mu(idx)    = mu_eq
              targets%l(idx)     = l
              targets%m(idx)     = m
              targets%lm(idx)    = l*(l+1) + m + 1
           end do
        end do
     end do

     lfirst = lfirst + stru%mult(jatom)
  end do

  ! Print target info
  write(unit_out, *)
  write(unit_out, '(A)')         ' PDWF Target Orbitals:'
  write(unit_out, '(A,I6)')      '   num_wann = ', Ntarget
  write(unit_out, '(A)')         '   atom   Z   l-channels'
  write(unit_out, '(A)')         '   ----  ---  ----------'
  lfirst = 1
  do jatom = 1, stru%Nneq
     iZ = nint(stru%Z(jatom))
     call valence_config(iZ, l_channels, n_channels)
     write(unit_out, '(I6, I5, 2X)', ADVANCE='no') jatom, iZ
     do l = 0, 3
        if (l_channels(l) == 1) then
           select case(l)
           case(0); write(unit_out, '(A)', ADVANCE='no') 's'
           case(1); write(unit_out, '(A)', ADVANCE='no') 'p'
           case(2); write(unit_out, '(A)', ADVANCE='no') 'd'
           case(3); write(unit_out, '(A)', ADVANCE='no') 'f'
           end select
        end if
     end do
     write(unit_out, '(A,I3,A)') '  (x', stru%mult(jatom), ')'
     lfirst = lfirst + stru%mult(jatom)
  end do
  write(unit_out, *)
end subroutine build_targets


!-----------------------------------------------------------------------
! Classify bands into frozen/disentangle/excluded based on projectability
!-----------------------------------------------------------------------
subroutine classify_bands(proj_avg, Nb, p_high, p_low, band_class, &
     &                    n_frozen, n_disent, n_excluded)
  real(R8), intent(in)  :: proj_avg(Nb)
  integer,  intent(in)  :: Nb
  real(R8), intent(in)  :: p_high, p_low
  integer,  intent(out) :: band_class(Nb)
  integer,  intent(out) :: n_frozen, n_disent, n_excluded

  integer :: ib

  n_frozen = 0; n_disent = 0; n_excluded = 0
  do ib = 1, Nb
     if (proj_avg(ib) >= p_high) then
        band_class(ib) = BAND_FROZEN
        n_frozen = n_frozen + 1
     else if (proj_avg(ib) >= p_low) then
        band_class(ib) = BAND_DISENTANGLE
        n_disent = n_disent + 1
     else
        band_class(ib) = BAND_EXCLUDED
        n_excluded = n_excluded + 1
     end if
  end do
end subroutine classify_bands


!-----------------------------------------------------------------------
! Determine energy windows from eigenvalues and band classification
!-----------------------------------------------------------------------
subroutine determine_windows(eigvals, band_class, Nb, Nk, margin, &
     &                       froz_min, froz_max, dis_min, dis_max, &
     &                       has_frozen, has_disent)
  real(R8), intent(in)  :: eigvals(Nb, Nk)
  integer,  intent(in)  :: band_class(Nb)
  integer,  intent(in)  :: Nb, Nk
  real(R8), intent(in)  :: margin
  real(R8), intent(out) :: froz_min, froz_max, dis_min, dis_max
  logical,  intent(out) :: has_frozen, has_disent

  integer  :: ib, ik
  real(R8) :: emin, emax

  has_frozen = .false.
  has_disent = .false.
  froz_min =  1.0d30; froz_max = -1.0d30
  dis_min  =  1.0d30; dis_max  = -1.0d30

  do ib = 1, Nb
     if (band_class(ib) == BAND_FROZEN) then
        has_frozen = .true.
        do ik = 1, Nk
           if (eigvals(ib,ik) < froz_min) froz_min = eigvals(ib,ik)
           if (eigvals(ib,ik) > froz_max) froz_max = eigvals(ib,ik)
        end do
     end if
     if (band_class(ib) /= BAND_EXCLUDED) then
        has_disent = .true.
        do ik = 1, Nk
           if (eigvals(ib,ik) < dis_min) dis_min = eigvals(ib,ik)
           if (eigvals(ib,ik) > dis_max) dis_max = eigvals(ib,ik)
        end do
     end if
  end do

  ! Apply margins
  if (has_frozen) then
     froz_min = froz_min - margin
     froz_max = froz_max + margin
  end if
  if (has_disent) then
     dis_min = dis_min - margin
     dis_max = dis_max + margin
  end if

  ! Validate: check for excluded bands inside frozen window
  if (has_frozen) then
     do ib = 1, Nb
        if (band_class(ib) == BAND_EXCLUDED) then
           do ik = 1, Nk
              if (eigvals(ib,ik) >= froz_min .and. &
                   eigvals(ib,ik) <= froz_max) then
                 write(unit_out, '(A,I4,A,I4,A)') &
                      ' WARNING: Excluded band ', ib, &
                      ' at k-point ', ik, &
                      ' falls inside frozen window!'
              end if
           end do
        end if
     end do
  end if
end subroutine determine_windows


!-----------------------------------------------------------------------
! Write PDWF diagnostic output
!-----------------------------------------------------------------------
subroutine write_pdwf_diagnostics(stru, targets, proj_avg, proj_k, &
     &     band_class, Nb, Nk, froz_min, froz_max, dis_min, dis_max, &
     &     has_frozen, has_disent)
  type(struct_t),      intent(in) :: stru
  type(pdwf_target_t), intent(in) :: targets
  real(R8),            intent(in) :: proj_avg(Nb), proj_k(Nb, Nk)
  integer,             intent(in) :: band_class(Nb)
  integer,             intent(in) :: Nb, Nk
  real(R8),            intent(in) :: froz_min, froz_max, dis_min, dis_max
  logical,             intent(in) :: has_frozen, has_disent

  integer  :: ib
  real(R8) :: pmin, pmax
  character(len=12) :: class_str

  write(unit_out, *)
  write(unit_out, '(A)') ' ===== PDWF Band Classification ====='
  write(unit_out, '(A)') &
       '  Band   avg(p)    min(p)    max(p)  Class'
  write(unit_out, '(A)') &
       '  ----  --------  --------  --------  -----------'

  do ib = 1, Nb
     pmin = minval(proj_k(ib,:))
     pmax = maxval(proj_k(ib,:))
     select case(band_class(ib))
     case(BAND_FROZEN);      class_str = 'FROZEN'
     case(BAND_DISENTANGLE); class_str = 'DISENTANGLE'
     case(BAND_EXCLUDED);    class_str = 'EXCLUDED'
     end select
     write(unit_out, '(I5, 3F10.4, 2X, A)') &
          ib, proj_avg(ib), pmin, pmax, trim(class_str)
  end do

  write(unit_out, *)
  write(unit_out, '(A)') ' ===== Wannier90 Window Parameters ====='
  write(unit_out, '(A,I6)')  '  num_wann         = ', targets%Ntarget
  if (has_frozen) then
     write(unit_out, '(A,F12.6)') '  dis_froz_min     = ', froz_min
     write(unit_out, '(A,F12.6)') '  dis_froz_max     = ', froz_max
  end if
  if (has_disent) then
     write(unit_out, '(A,F12.6)') '  dis_win_min      = ', dis_min
     write(unit_out, '(A,F12.6)') '  dis_win_max      = ', dis_max
  end if

  ! Print exclude_bands list
  write(unit_out, '(A)', ADVANCE='no') '  exclude_bands    = '
  do ib = 1, Nb
     if (band_class(ib) == BAND_EXCLUDED) then
        write(unit_out, '(I0, A)', ADVANCE='no') ib, ','
     end if
  end do
  write(unit_out, *)
  write(unit_out, '(A)') ' ======================================='
  write(unit_out, *)
end subroutine write_pdwf_diagnostics


!-----------------------------------------------------------------------
! Main PDWF computation: projectability + AMN generation
!
! Follows the same atom-loop + k-loop pattern as l2Amn.
! For each atom at each k-point, calls almgen to get Alm, then:
!   - Extracts target orbital columns into A_full
!   - Accumulates |Alm|^2 into proj_k for projectability
! After all atoms: classifies bands, determines windows, writes AMN.
!-----------------------------------------------------------------------
subroutine pdwf_compute(stru, inwf, targets, Nb, Nk)
  use xa,        only: BK, FJ, DFJ, PHS, BKrot, BKRloc
  use xa3,       only: A, XK,YK,ZK, GX,GY,GZ, BKX,BKY,BKZ, vecsz
  use lolog,     only: nlo,nlov,nlon,loor,ilo,lapw
  use loabc,     only: alo
  use atspdt,    only: P, DP
  use radfu,     only: RF1, RF2
  use eigenval_store, only: eigvals

  !! procedure includes
  use rint13_m
  use harmon_m
  use atpar_m
  use Ylm_m

  implicit none

  type(struct_t),      intent(in) :: stru
  type(inwf_t),        intent(in) :: inwf
  type(pdwf_target_t), intent(in) :: targets
  integer,             intent(in) :: Nb, Nk

  ! Local arrays
  complex(C16), allocatable :: A_full(:,:,:) ! (Nb, Ntarget, Nk) AMN matrix
  real(R8),     allocatable :: proj_k(:,:)   ! (Nb, Nk) per-band projectability
  real(R8),     allocatable :: proj_avg(:)   ! (Nb) k-averaged projectability
  integer,      allocatable :: band_class(:) ! (Nb) classification

  ! Alm workspace (same layout as l2Amn)
  complex(C16) :: Yl((Lmax2+1)**2)
  complex(C16) :: Alm((Lmax2+1)**2, Nb, stru%Nat, Nrf)
  complex(C16) :: phshel, cfac

  real(R8) :: h_al(iblock), h_bl(iblock)
  complex(C16) :: h_yl  ((Lmax2+1)**2, iblock)
  complex(C16) :: h_alyl((Lmax2+1)**2, iblock)
  complex(C16) :: h_blyl((Lmax2+1)**2, iblock)

  real(R8) :: pi12lo(0:Lmax2)
  real(R8) :: fac, arg
  real(R8) :: froz_min, froz_max, dis_min, dis_max
  logical  :: has_frozen, has_disent
  integer  :: n_frozen, n_disent, n_excluded

  integer :: i, ii, i3, l, m, m1, mu_eq, n, idx, ibb, iscf
  integer :: jatom, jlo, jneq, irf, lda, ldb, ldc
  integer :: kkk, max_m, itape, jtape, latom, lfirst
  integer :: itgt, ib, ipr

  ! Allocate working arrays
  allocate(A_full(Nb, targets%Ntarget, Nk))
  allocate(proj_k(Nb, Nk))
  allocate(proj_avg(Nb))
  allocate(band_class(Nb))

  A_full  = (0.0d0, 0.0d0)
  proj_k  = 0.0d0

  ! ------ Phase 1-2: Compute Alm and extract target projections ------

  rewind(unit_vsp)
  read(unit_vsp, '(50X,I2,//)') ISCF
  LFIRST = 1

  atoms: do JATOM = 1, stru%Nneq
     if (jatom > 1) lfirst = lfirst + stru%mult(jatom-1)
     itape = unit_vector
     jtape = unit_vsp
     rewind(itape)
     call atpar(stru, jatom, itape, jtape)
     rewind(itape)

     ! Compute LO overlap corrections (same as l2Amn)
     pi12lo = 0
     do l = 0, LOmax
        if (loor(l)) call rint13(stru, jatom, rf1(:,l,1), rf2(:,l,1), &
             &                   rf1(:,l,3), rf2(:,l,3), pi12lo(l))
     end do

     kpoints: do KKK = 1, Nk
        N = vecsz(kkk)

        ! plane wave Q+k
        BKX(1:N) = XK(kkk) + GX(1:N, kkk)
        BKY(1:N) = YK(kkk) + GY(1:N, kkk)
        BKZ(1:N) = ZK(kkk) + GZ(1:N, kkk)

        ! calculate Bessel functions j_l(|Q+k|*RMT)
        call harmon(stru, jatom, N, Bkx, Bky, Bkz, Lmax2, FJ, DFJ)

        fac = 2*TAU * stru%RMT(jatom)**2 / sqrt(stru%vol)
        Alm = 0

        ! ---------- Compute Alm for each equivalent atom ----------
        equiv: do mu_eq = 1, stru%mult(jatom)
           latom = lfirst - 1 + mu_eq

           ! matching of plane wave and atomic functions (blocked)
           do ii = 1, N-(nlo+nlon+nlov), iblock
              i3 = 0
              do i = ii, min(ii+iblock-1, N-(nlo+nlon+nlov))
                 i3 = i3 + 1

                 bkrot = matmul(stru%rotij(:,:,latom), &
                      &         (/bkx(i), bky(i), bkz(i)/))
                 bk = matmul(stru%conv_rec, bkrot)
                 bkrloc = matmul(stru%locrot(:,:,jatom), bk)
                 call YLM(bkrloc, LMAX2, YL)
                 arg = dot_product(bkrot, stru%pos(:, lfirst)) &
                      + dot_product((/bkx(i), bky(i), bkz(i)/), &
                      &             stru%transij(:,latom))
                 phshel = exp((0,1)*arg*TAU)

                 idx = 0
                 do l = 0, Lmax2
                    do m = 1, 2*l+1
                       idx = idx + 1
                       h_yl(idx, i3) = conjg(yl(idx)) * phshel
                    enddo
                 enddo
              enddo

              IDX = 0
              do L = 0, LMAX2
                 i3 = 0
                 do i = ii, min(ii+iblock-1, N-(nlo+nlon+nlov))
                    i3 = i3 + 1
                    if (lapw(l)) then
                       h_AL(i3) = DFJ(L,I)*P(L,2) - FJ(L,I)*DP(L,2)
                       h_BL(i3) = FJ(L,I)*DP(L,1) - DFJ(L,I)*P(L,1)
                    else
                       h_AL(i3) = FJ(L,I)/P(L,1)/stru%RMT(jatom)**2
                       h_BL(i3) = 0
                    endif
                 enddo
                 MAX_m = 2*L+1
                 do M = 1, MAX_m
                    IDX = IDX + 1
                    i3 = 0
                    do i = ii, min(ii+iblock-1, N-(nlo+nlon+nlov))
                       i3 = i3 + 1
                       h_alyl(idx,i3) = h_AL(i3) * h_YL(idx,i3)
                       h_blyl(idx,i3) = h_BL(i3) * h_YL(idx,i3)
                    enddo
                 enddo
              enddo
              ibb = min(iblock, N-(nlo+nlon+nlov)-ii+1)
#ifndef _COMPLEX_
              lda = 2*(LMAX2+1)*(LMAX2+1)
#else
              lda = (LMAX2+1)*(LMAX2+1)
#endif
              ldc = lda
              ldb = nmat

#ifndef _COMPLEX_
              call dgemm('N','N',2*idx,Nb,ibb,1.d0, &
                   h_alyl,lda,a(ii,1,kkk),ldb,1.d0, &
                   alm(1,1,mu_eq,1),ldc)
              call dgemm('N','N',2*idx,Nb,ibb,1.d0, &
                   h_blyl,lda,a(ii,1,kkk),ldb,1.d0, &
                   alm(1,1,mu_eq,2),ldc)
#else
              call zgemm('N','N',idx,Nb,ibb,(1.d0,0.d0), &
                   h_alyl,lda,a(ii,1,kkk),ldb,(1.d0,0.d0), &
                   alm(1,1,mu_eq,1),ldc)
              call zgemm('N','N',idx,Nb,ibb,(1.d0,0.d0), &
                   h_blyl,lda,a(ii,1,kkk),ldb,(1.d0,0.d0), &
                   alm(1,1,mu_eq,2),ldc)
#endif
           enddo

           ! ******** LOCAL ORBITALS *****************************
           lo: if (nlo /= 0) then
              i = n - (nlo+nlon)
              do L = 0, LoMAX
                 do jlo = 1, ilo(l)
                    do jneq = 1, stru%mult(jatom)
                       do m1 = -l, l
                          i = i + 1
                          BK(1) = BKX(I)
                          BK(2) = BKY(I)
                          BK(3) = BKZ(I)
                          bkrot = matmul(stru%rotij(:,:,latom), BK)
                          BK    = matmul(stru%conv_rec, BKROT)
                          bkrloc = matmul(stru%locrot(:,:,jatom), BK)
                          call Ylm(Bkrloc, LOmax, YL)
                          arg = dot_product(bkrot(:), stru%pos(:, lfirst)) &
                               + dot_product((/bkx(i), bky(i), bkz(i)/), &
                               &             stru%transij(:,latom))
                          phshel = exp((0,1)*TAU*arg)
                          PHS(:) = phshel * A(i,:,kkk)

                          do m = -l, l
                             idx = l*(l+1) + m + 1
                             do ib = 1, Nb
                                Alm(idx, ib, mu_eq, :) = &
                                     Alm(idx, ib, mu_eq, :) + &
                                     Alo(l, jlo, :) * conjg(Yl(idx)) * phs(ib)
                             end do
                          end do
                       end do
                    end do
                 end do
              end do
           end if lo

           ! Apply normalization factor
           idx = 0
           do l = 0, Lmax2
              cfac = (0,1)**l
              do m = 1, 2*l+1
                 idx = idx + 1
                 do ib = 1, Nb
                    do irf = 1, Nrf
                       Alm(idx,ib,mu_eq,irf) = Alm(idx,ib,mu_eq,irf)*fac*cfac
                    enddo
                 enddo
              enddo
           enddo
        enddo equiv

        ! ---------- Extract target projections from Alm ----------
        do itgt = 1, targets%Ntarget
           if (targets%jatom(itgt) /= jatom) cycle
           mu_eq = targets%mu(itgt)
           idx = targets%lm(itgt)
           l = targets%l(itgt)

           do ib = 1, Nb
              ! AMN = conjg(Alm) for rf=1 + LO correction (rf=3)
              A_full(ib, itgt, kkk) = conjg(Alm(idx, ib, mu_eq, 1)) &
                   + conjg(Alm(idx, ib, mu_eq, 3)) * pi12lo(l)

              ! Projectability: |projection|^2
              proj_k(ib, kkk) = proj_k(ib, kkk) + &
                   abs(A_full(ib, itgt, kkk))**2
           end do
        end do

     enddo kpoints
  enddo atoms

  ! ------ Phase 3: K-average projectability and classify bands ------

  do ib = 1, Nb
     proj_avg(ib) = sum(proj_k(ib,:)) / dble(Nk)
  end do

  call classify_bands(proj_avg, Nb, inwf%pdwf_p_high, inwf%pdwf_p_low, &
       &              band_class, n_frozen, n_disent, n_excluded)

  write(unit_out, '(A,I4,A)') ' PDWF: ', n_frozen,   ' frozen bands'
  write(unit_out, '(A,I4,A)') ' PDWF: ', n_disent,   ' disentangle bands'
  write(unit_out, '(A,I4,A)') ' PDWF: ', n_excluded,  ' excluded bands'

  ! ------ Phase 4: Determine energy windows ------

  call determine_windows(eigvals, band_class, Nb, Nk, inwf%pdwf_margin, &
       &                 froz_min, froz_max, dis_min, dis_max, &
       &                 has_frozen, has_disent)

  ! ------ Write diagnostics ------

  call write_pdwf_diagnostics(stru, targets, proj_avg, proj_k, &
       &     band_class, Nb, Nk, froz_min, froz_max, dis_min, dis_max, &
       &     has_frozen, has_disent)

  ! ------ Phase 5: Write AMN file ------

  do kkk = 1, Nk
     do ipr = 1, targets%Ntarget
        do ib = 1, Nb
           write(unit_amn, '(2I4,1X,I5,1X,2E18.5)') &
                ib, ipr, kkk, A_full(ib, ipr, kkk)
        enddo
     enddo
  enddo

  ! Cleanup
  deallocate(A_full, proj_k, proj_avg, band_class)
end subroutine pdwf_compute

end module pdwf_m
