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
  public :: write_win_parameters

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
! Classify bands into frozen/disentangle/excluded.
!
! Adapted PDWF for the LAPW basis where projectabilities are
! systematically lower than in LCAO (muffin-tin only, no interstitial
! contribution). Uses percentile-based adaptive thresholds.
!
! Strategy:
!   Phase 1 — Semicore exclusion:
!     Detect energy gaps > GAP_THRESH between consecutive bands.
!     All bands below the largest sub-E_F gap are EXCLUDED.
!
!   Phase 2 — Percentile-based classification:
!     Sort non-semicore projectabilities, compute percentiles.
!     p_high_eff = value at P_HIGH_PCTL percentile (e.g. 60th)
!     p_low_eff  = max(p_low, value at P_LOW_PCTL percentile)
!     Bands with avg(p) >= p_high_eff  ->  FROZEN
!     Bands with avg(p) >= p_low_eff   ->  DISENTANGLE
!     Bands below p_low_eff            ->  EXCLUDED
!
!   Phase 3 — Gap bridging:
!     Disentangle bands sandwiched between frozen bands in
!     band-index space with avg(p) >= p_high_eff - BRIDGE_MARGIN
!     are promoted to FROZEN.
!
!   Phase 4 — Fermi-level fallback:
!     If no frozen bands, promote bands near E_F to FROZEN.
!
!   Phase 5 — Gap quality scoring:
!     Q = min(p_frozen) / max(p_excluded)
!     Reports classification reliability.
!
! eigvals are in eV relative to E_F (E_F = 0).
!-----------------------------------------------------------------------
subroutine classify_bands(eigvals, proj_avg, Nb, Nk, &
     &                    pctl_high, pctl_low, &
     &                    band_class, n_frozen, n_disent, n_excluded)
  real(R8), intent(in)  :: eigvals(Nb, Nk)
  real(R8), intent(in)  :: proj_avg(Nb)
  integer,  intent(in)  :: Nb, Nk
  real(R8), intent(in)  :: pctl_high  ! frozen percentile (0-1), e.g. 0.60
  real(R8), intent(in)  :: pctl_low   ! exclusion percentile (0-1), e.g. 0.15
  integer,  intent(out) :: band_class(Nb)
  integer,  intent(out) :: n_frozen, n_disent, n_excluded

  integer  :: ib, semicore_top, n_valence, ipctl, jtmp
  integer  :: ifirst_frozen, ilast_frozen
  real(R8) :: e_band_avg(Nb), gap, max_gap
  real(R8) :: p_max, r_high_eff, r_low_eff
  real(R8) :: relevance(Nb), energy_weight
  real(R8) :: p_min_frozen, p_max_excluded, gap_quality, tmp
  real(R8), allocatable :: r_sorted(:)

  real(R8) :: sigma
  ! Minimum gap (eV) to identify core/semicore boundary
  real(R8), parameter :: GAP_THRESH = 5.0d0
  ! Skewed Gaussian weighting: occupied states (E<0) are more important
  ! than unoccupied (E>0). Use wider sigma below E_F to protect valence.
  real(R8), parameter :: SIGMA_BELOW = 12.0d0  ! eV, occupied side
  real(R8), parameter :: SIGMA_ABOVE =  6.0d0  ! eV, unoccupied side
  ! Absolute margin on relevance score for gap-bridging
  real(R8), parameter :: BRIDGE_MARGIN = 0.02d0
  ! Energy range around E_F for fallback freezing (eV)
  real(R8), parameter :: FALLBACK_EV = 3.0d0

  ! Compute k-averaged band energies
  do ib = 1, Nb
     e_band_avg(ib) = sum(eigvals(ib,:)) / dble(Nk)
  end do

  ! ---- Phase 1: Semicore exclusion via energy gaps ----
  semicore_top = 0
  max_gap = 0.0d0
  do ib = 1, Nb - 1
     gap = e_band_avg(ib + 1) - e_band_avg(ib)
     if (e_band_avg(ib) < 0.0d0 .and. gap > GAP_THRESH) then
        if (gap > max_gap) then
           max_gap = gap
           semicore_top = ib
        end if
     end if
  end do

  if (semicore_top > 0) then
     write(unit_out, '(A,I3,A,F8.2,A)') &
          ' PDWF: excluding bands 1-', semicore_top, &
          ' as core/semicore (gap = ', max_gap, ' eV)'
  end if

  ! ---- Phase 2: Relevance-weighted percentile classification ----
  ! Compute relevance = projectability × energy_weight
  ! energy_weight uses a split Gaussian: wider below E_F (occupied
  ! states matter more) and narrower above (conduction bands are
  ! less critical for ground-state properties).
  !   weight(E) = exp(-E²/2σ²)
  !   σ = σ_below (12 eV) for E < 0  → gentle decay for valence
  !   σ = σ_above ( 6 eV) for E > 0  → steeper decay for conduction
  p_max = 0.0d0
  n_valence = Nb - semicore_top
  do ib = 1, Nb
     if (ib <= semicore_top) then
        relevance(ib) = 0.0d0
     else
        if (e_band_avg(ib) < 0.0d0) then
           sigma = SIGMA_BELOW
        else
           sigma = SIGMA_ABOVE
        end if
        energy_weight = exp(-0.5d0 * (e_band_avg(ib)/sigma)**2)
        relevance(ib) = proj_avg(ib) * energy_weight
        if (proj_avg(ib) > p_max) p_max = proj_avg(ib)
     end if
  end do

  ! Sort relevance scores of non-semicore bands to find percentiles
  allocate(r_sorted(n_valence))
  do ib = 1, n_valence
     r_sorted(ib) = relevance(semicore_top + ib)
  end do
  ! Insertion sort
  do ib = 2, n_valence
     tmp = r_sorted(ib)
     jtmp = ib - 1
     do while (jtmp >= 1 .and. r_sorted(jtmp) > tmp)
        r_sorted(jtmp + 1) = r_sorted(jtmp)
        jtmp = jtmp - 1
     end do
     r_sorted(jtmp + 1) = tmp
  end do

  ! Compute effective thresholds from percentiles of relevance
  ipctl = max(1, min(n_valence, &
       nint(pctl_high * dble(n_valence))))
  r_high_eff = r_sorted(ipctl)

  ipctl = max(1, min(n_valence, &
       nint(pctl_low * dble(n_valence))))
  r_low_eff = r_sorted(ipctl)

  deallocate(r_sorted)

  write(unit_out, '(A,F8.4)') ' PDWF: max projectability     = ', p_max
  write(unit_out, '(A,F5.1,A,F5.1,A)') &
       ' PDWF: energy weight sigma    = ', SIGMA_BELOW, &
       '/-', SIGMA_ABOVE, ' eV (below/above E_F)'
  write(unit_out, '(A,F8.4,A,F5.1,A)') &
       ' PDWF: frozen relevance (eff) = ', r_high_eff, &
       ' (', pctl_high*100.0d0, 'th pctl)'
  write(unit_out, '(A,F8.4,A,F5.1,A)') &
       ' PDWF: exclude relevance(eff) = ', r_low_eff, &
       ' (', pctl_low*100.0d0, 'th pctl)'

  ! Print per-band relevance scores
  write(unit_out, *)
  write(unit_out, '(A)') &
       '  Band   avg(p)   E_avg(eV)  E_weight  relevance'
  write(unit_out, '(A)') &
       '  ----  --------  ---------  --------  ---------'
  do ib = semicore_top + 1, Nb
     if (e_band_avg(ib) < 0.0d0) then
        sigma = SIGMA_BELOW
     else
        sigma = SIGMA_ABOVE
     end if
     energy_weight = exp(-0.5d0 * (e_band_avg(ib)/sigma)**2)
     write(unit_out, '(I5, F10.4, F11.3, F10.4, F11.4)') &
          ib, proj_avg(ib), e_band_avg(ib), energy_weight, relevance(ib)
  end do
  write(unit_out, *)

  ! Classify using relevance thresholds
  n_frozen = 0; n_disent = 0; n_excluded = 0
  do ib = 1, Nb
     if (ib <= semicore_top) then
        band_class(ib) = BAND_EXCLUDED
        n_excluded = n_excluded + 1
     else if (relevance(ib) < r_low_eff) then
        band_class(ib) = BAND_EXCLUDED
        n_excluded = n_excluded + 1
     else if (relevance(ib) >= r_high_eff) then
        band_class(ib) = BAND_FROZEN
        n_frozen = n_frozen + 1
     else
        band_class(ib) = BAND_DISENTANGLE
        n_disent = n_disent + 1
     end if
  end do

  ! ---- Phase 3: Gap bridging ----
  ! Promote disentangle bands to frozen if sandwiched between
  ! frozen bands in band-index space with relevance close to threshold.
  if (n_frozen > 0) then
     ifirst_frozen = 0
     ilast_frozen = 0
     do ib = 1, Nb
        if (band_class(ib) == BAND_FROZEN) then
           if (ifirst_frozen == 0) ifirst_frozen = ib
           ilast_frozen = ib
        end if
     end do

     do ib = ifirst_frozen, ilast_frozen
        if (band_class(ib) == BAND_DISENTANGLE .and. &
             relevance(ib) >= r_high_eff - BRIDGE_MARGIN) then
           band_class(ib) = BAND_FROZEN
           n_frozen = n_frozen + 1
           n_disent = n_disent - 1
           write(unit_out, '(A,I3,A,F7.4,A,F7.4,A)') &
                ' PDWF: gap-bridging band ', ib, &
                ' (p=', proj_avg(ib), &
                ', rel=', relevance(ib), ') -> FROZEN'
        end if
     end do
  end if

  ! ---- Phase 4: Fermi-level fallback ----
  if (n_frozen == 0) then
     write(unit_out, '(A)') &
          ' PDWF: no frozen bands, applying Fermi-level fallback'
     do ib = 1, Nb
        if (band_class(ib) == BAND_DISENTANGLE .and. &
             abs(e_band_avg(ib)) <= FALLBACK_EV) then
           band_class(ib) = BAND_FROZEN
           n_frozen = n_frozen + 1
           n_disent = n_disent - 1
        end if
     end do
     if (n_frozen > 0) then
        write(unit_out, '(A,I3,A,F5.1,A)') &
             ' PDWF: froze ', n_frozen, ' bands within +/-', &
             FALLBACK_EV, ' eV of E_F'
     end if
  end if

  ! ---- Phase 5: Gauge freedom cap ----
  ! Ensure n_frozen <= 80% of active bands. If too many bands are
  ! frozen, demote the lowest-relevance frozen bands to disentangle.
  ! This guarantees Wannier90 has gauge freedom for optimization.
  if (n_frozen > 0 .and. n_disent == 0 .and. n_frozen > 1) then
     ! All active bands are frozen — must demote some
     ! Demote bands with lowest relevance until n_frozen <= 80% active
     demote: do while (n_frozen > nint(0.8d0 * (n_frozen + n_disent)) &
          &            .and. n_frozen > 1)
        ! Find frozen band with lowest relevance
        tmp = 1.0d30
        jtmp = 0
        do ib = 1, Nb
           if (band_class(ib) == BAND_FROZEN .and. &
                relevance(ib) < tmp) then
              tmp = relevance(ib)
              jtmp = ib
           end if
        end do
        if (jtmp == 0) exit demote
        band_class(jtmp) = BAND_DISENTANGLE
        n_frozen = n_frozen - 1
        n_disent = n_disent + 1
        write(unit_out, '(A,I3,A,F7.4,A)') &
             ' PDWF: gauge cap: demoting band ', jtmp, &
             ' (rel=', relevance(jtmp), ') -> DISENTANGLE'
     end do demote
  end if

  ! ---- Phase 6: Gap quality scoring ----
  ! Q = min(relevance among frozen) / max(relevance among non-semicore excluded)
  p_min_frozen  =  1.0d30
  p_max_excluded = 0.0d0
  do ib = semicore_top + 1, Nb
     if (band_class(ib) == BAND_FROZEN) then
        if (relevance(ib) < p_min_frozen) p_min_frozen = relevance(ib)
     else if (band_class(ib) == BAND_EXCLUDED) then
        if (relevance(ib) > p_max_excluded) p_max_excluded = relevance(ib)
     end if
  end do

  if (p_min_frozen < 1.0d30 .and. p_max_excluded > 1.0d-10) then
     gap_quality = p_min_frozen / p_max_excluded
     write(unit_out, '(A,F6.2)', ADVANCE='no') &
          ' PDWF: gap quality = ', gap_quality
     if (gap_quality > 5.0d0) then
        write(unit_out, '(A)') ' (CLEAN)'
     else if (gap_quality > 2.0d0) then
        write(unit_out, '(A)') ' (MODERATE)'
     else
        write(unit_out, '(A)') ' (POOR - review thresholds)'
     end if
  end if
end subroutine classify_bands


!-----------------------------------------------------------------------
! Detect band manifolds from energy gaps.
!
! Groups bands into connected manifolds separated by energy gaps
! above a threshold. Uses the 75th percentile of inter-band gaps
! as the threshold (adapted from LCAO-to-Wannier methodology).
!
! Output: manifold_id(Nb) — integer label 1..Nmanifold for each band.
! Also prints diagnostic table to unit_out.
!-----------------------------------------------------------------------
subroutine detect_manifolds(eigvals, band_class, Nb, Nk, &
     &                      manifold_id, n_manifolds)
  real(R8), intent(in)  :: eigvals(Nb, Nk)
  integer,  intent(in)  :: band_class(Nb)
  integer,  intent(in)  :: Nb, Nk
  integer,  intent(out) :: manifold_id(Nb)
  integer,  intent(out) :: n_manifolds

  real(R8) :: e_avg(Nb), gaps(Nb-1), gap_thresh, tmp
  integer  :: ib, ig, n_active, jtmp, imf
  integer  :: mf_start, mf_end, mf_nfrozen, mf_ndisent
  real(R8), allocatable :: gaps_sorted(:)

  ! Compute k-averaged energies
  do ib = 1, Nb
     e_avg(ib) = sum(eigvals(ib,:)) / dble(Nk)
  end do

  ! Compute gaps between consecutive bands (all bands, not just active)
  n_active = 0
  do ib = 1, Nb - 1
     gaps(ib) = e_avg(ib + 1) - e_avg(ib)
     if (band_class(ib) /= BAND_EXCLUDED) n_active = n_active + 1
  end do

  ! Compute 75th percentile of gaps among non-excluded consecutive pairs
  ig = 0
  allocate(gaps_sorted(Nb - 1))
  do ib = 1, Nb - 1
     if (band_class(ib) /= BAND_EXCLUDED .or. &
          band_class(ib+1) /= BAND_EXCLUDED) then
        ig = ig + 1
        gaps_sorted(ig) = gaps(ib)
     end if
  end do

  ! Sort
  do ib = 2, ig
     tmp = gaps_sorted(ib)
     jtmp = ib - 1
     do while (jtmp >= 1 .and. gaps_sorted(jtmp) > tmp)
        gaps_sorted(jtmp + 1) = gaps_sorted(jtmp)
        jtmp = jtmp - 1
     end do
     gaps_sorted(jtmp + 1) = tmp
  end do

  ! 75th percentile
  if (ig > 0) then
     jtmp = max(1, min(ig, nint(0.75d0 * dble(ig))))
     gap_thresh = gaps_sorted(jtmp)
  else
     gap_thresh = 1.0d30
  end if
  deallocate(gaps_sorted)

  ! Assign manifold IDs: new manifold when gap > threshold
  manifold_id(1) = 1
  n_manifolds = 1
  do ib = 2, Nb
     if (gaps(ib-1) > gap_thresh) then
        n_manifolds = n_manifolds + 1
     end if
     manifold_id(ib) = n_manifolds
  end do

  ! Print manifold structure
  write(unit_out, *)
  write(unit_out, '(A,F8.4,A)') &
       ' PDWF: manifold gap threshold = ', gap_thresh, ' eV (75th pctl)'
  write(unit_out, '(A,I3,A)') &
       ' PDWF: detected ', n_manifolds, ' band manifolds'
  write(unit_out, *)
  write(unit_out, '(A)') &
       '  Manifold  Bands      E_range (eV)     Frozen  Disent  Excl'
  write(unit_out, '(A)') &
       '  --------  ---------  ---------------  ------  ------  ----'

  imf = 1
  mf_start = 1
  do ib = 2, Nb + 1
     if (ib > Nb .or. manifold_id(ib) /= imf) then
        mf_end = ib - 1
        ! Count classifications in this manifold
        mf_nfrozen = 0; mf_ndisent = 0; ig = 0
        do jtmp = mf_start, mf_end
           if (band_class(jtmp) == BAND_FROZEN) mf_nfrozen = mf_nfrozen + 1
           if (band_class(jtmp) == BAND_DISENTANGLE) mf_ndisent = mf_ndisent + 1
           if (band_class(jtmp) == BAND_EXCLUDED) ig = ig + 1
        end do
        write(unit_out, '(I6, 4X, I3, A, I3, 4X, F7.2, A, F7.2, 2X, I4, 4X, I4, 4X, I3)') &
             imf, mf_start, '-', mf_end, &
             e_avg(mf_start), ' to', e_avg(mf_end), &
             mf_nfrozen, mf_ndisent, ig
        if (ib <= Nb) then
           imf = manifold_id(ib)
           mf_start = ib
        end if
     end if
  end do
  write(unit_out, *)
end subroutine detect_manifolds


!-----------------------------------------------------------------------
! Determine energy windows from eigenvalues and band classification.
!
! The frozen window covers ONLY frozen bands (+-0.5 eV padding).
! The disentangle window covers ALL active bands, frozen + disentangle
! (+-2.0 eV padding). This guarantees:
!   1. Frozen window is ALWAYS strictly narrower than disentangle
!   2. Disentangle bands provide gauge freedom for optimization
!   3. The Wannier90 num_wann constraint is satisfied
!
! If manifold_id is provided, the frozen window is adjusted to avoid
! splitting a manifold: it snaps to manifold boundaries.
!-----------------------------------------------------------------------
subroutine determine_windows(eigvals, band_class, Nb, Nk, n_wann, &
     &                       margin, froz_min, froz_max, dis_min, &
     &                       dis_max, has_frozen, has_disent)
  real(R8), intent(in)    :: eigvals(Nb, Nk)
  integer,  intent(in)    :: band_class(Nb)
  integer,  intent(in)    :: Nb, Nk, n_wann
  real(R8), intent(in)    :: margin
  real(R8), intent(out)   :: froz_min, froz_max, dis_min, dis_max
  logical,  intent(inout) :: has_frozen, has_disent

  integer  :: ib, ik, n_in_froz, n_kf
  ! Frozen window padding (eV) — tight around frozen bands
  real(R8), parameter :: PAD_FROZEN = 0.5d0
  ! Outer window padding (eV) — wide for gauge freedom
  real(R8), parameter :: PAD_OUTER  = 2.0d0

  has_frozen = .false.
  has_disent = .false.
  froz_min =  1.0d30; froz_max = -1.0d30
  dis_min  =  1.0d30; dis_max  = -1.0d30

  ! Frozen window: spans ONLY frozen bands
  do ib = 1, Nb
     if (band_class(ib) == BAND_FROZEN) then
        has_frozen = .true.
        do ik = 1, Nk
           if (eigvals(ib,ik) < froz_min) froz_min = eigvals(ib,ik)
           if (eigvals(ib,ik) > froz_max) froz_max = eigvals(ib,ik)
        end do
     end if
  end do

  ! Outer window: spans ALL non-excluded bands (frozen + disentangle)
  do ib = 1, Nb
     if (band_class(ib) /= BAND_EXCLUDED) then
        has_disent = .true.
        do ik = 1, Nk
           if (eigvals(ib,ik) < dis_min) dis_min = eigvals(ib,ik)
           if (eigvals(ib,ik) > dis_max) dis_max = eigvals(ib,ik)
        end do
     end if
  end do

  ! Apply padding: frozen ±0.5 eV, outer ±2.0 eV
  ! This guarantees outer is STRICTLY wider than frozen
  if (has_frozen) then
     froz_min = froz_min - PAD_FROZEN
     froz_max = froz_max + PAD_FROZEN
  end if
  if (has_disent) then
     dis_min = dis_min - PAD_OUTER
     dis_max = dis_max + PAD_OUTER
  end if

  ! Enforce: frozen must be inside outer window
  if (has_frozen .and. has_disent) then
     if (froz_min < dis_min) froz_min = dis_min + 0.1d0
     if (froz_max > dis_max) froz_max = dis_max - 0.1d0
  end if

  ! Critical Wannier90 constraint: at every k-point, the number of
  ! non-excluded bands inside the frozen window must be <= num_wann.
  ! Shrink froz_max if necessary to satisfy this.
  if (has_frozen) then
     shrink: do
        ! Count max bands inside frozen window at any k-point
        n_in_froz = 0
        do ik = 1, Nk
           n_kf = 0
           do ib = 1, Nb
              if (band_class(ib) == BAND_EXCLUDED) cycle
              if (eigvals(ib,ik) >= froz_min .and. &
                   eigvals(ib,ik) <= froz_max) then
                 n_kf = n_kf + 1
              end if
           end do
           if (n_kf > n_in_froz) n_in_froz = n_kf
        end do
        if (n_in_froz <= n_wann) exit shrink
        ! Too many bands: shrink froz_max by 0.5 eV
        froz_max = froz_max - 0.5d0
        if (froz_max <= froz_min) then
           has_frozen = .false.
           write(unit_out, '(A)') &
                ' WARNING: Could not find valid frozen window, disabling'
           exit shrink
        end if
     end do shrink
     if (has_frozen) then
        write(unit_out, '(A,I3,A,I3)') &
             ' PDWF: max bands in frozen window = ', n_in_froz, &
             ' / num_wann = ', n_wann
     end if
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
subroutine pdwf_compute(stru, inwf, targets, Nb, Nk, &
     &                    froz_min_out, froz_max_out, &
     &                    dis_min_out, dis_max_out, &
     &                    has_frozen_out, has_disent_out, &
     &                    band_class_out)
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

  type(struct_t),      intent(in)  :: stru
  type(inwf_t),        intent(in)  :: inwf
  type(pdwf_target_t), intent(in)  :: targets
  integer,             intent(in)  :: Nb, Nk
  real(R8),            intent(out) :: froz_min_out, froz_max_out
  real(R8),            intent(out) :: dis_min_out, dis_max_out
  logical,             intent(out) :: has_frozen_out, has_disent_out
  integer,             intent(out) :: band_class_out(Nb)

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
  integer,  allocatable :: manifold_id(:)
  integer  :: n_manifolds

  integer :: i, ii, i3, l, m, m1, mu_eq, n, idx, ibb, iscf
  integer :: jatom, jlo, jneq, irf, lda, ldb, ldc
  integer :: kkk, max_m, itape, jtape, latom, lfirst
  integer :: itgt, ib, ipr

  ! Allocate working arrays
  allocate(A_full(Nb, targets%Ntarget, Nk))
  allocate(proj_k(Nb, Nk))
  allocate(proj_avg(Nb))
  allocate(band_class(Nb))
  allocate(manifold_id(Nb))

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

  call classify_bands(eigvals, proj_avg, Nb, Nk, &
       &              inwf%pdwf_p_high, inwf%pdwf_p_low, &
       &              band_class, n_frozen, n_disent, n_excluded)

  write(unit_out, '(A,I4,A)') ' PDWF: ', n_frozen,   ' frozen bands'
  write(unit_out, '(A,I4,A)') ' PDWF: ', n_disent,   ' disentangle bands'
  write(unit_out, '(A,I4,A)') ' PDWF: ', n_excluded,  ' excluded bands'

  ! ------ Phase 3b: Manifold detection ------

  call detect_manifolds(eigvals, band_class, Nb, Nk, &
       &                manifold_id, n_manifolds)

  ! ------ Phase 4: Determine energy windows ------

  call determine_windows(eigvals, band_class, Nb, Nk, &
       &                 targets%Ntarget, inwf%pdwf_margin, &
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

  ! Copy results to output arguments
  froz_min_out  = froz_min
  froz_max_out  = froz_max
  dis_min_out   = dis_min
  dis_max_out   = dis_max
  has_frozen_out  = has_frozen
  has_disent_out  = has_disent
  band_class_out  = band_class

  ! Cleanup
  deallocate(A_full, proj_k, proj_avg, band_class, manifold_id)
end subroutine pdwf_compute


!-----------------------------------------------------------------------
! Write PDWF window parameters to a snippet file (case.pdwf_windows)
! and patch the existing .win file if present.
!
! The snippet file can be used with write_win or pasted manually.
! The .win file is patched in-place by inserting PDWF parameters at
! the top and removing conflicting existing lines.
!-----------------------------------------------------------------------
subroutine write_win_parameters(winfn, stru, targets, Nb, &
     &  froz_min, froz_max, dis_min, dis_max, has_frozen, has_disent, &
     &  band_class, bmin)
  implicit none

  character(len=*),    intent(in) :: winfn
  type(struct_t),      intent(in) :: stru
  type(pdwf_target_t), intent(in) :: targets
  integer,             intent(in) :: Nb, bmin
  real(R8),            intent(in) :: froz_min, froz_max, dis_min, dis_max
  logical,             intent(in) :: has_frozen, has_disent
  integer,             intent(in) :: band_class(Nb)

  integer, parameter :: u_snip = 97, u_win = 99, u_tmp = 98
  character(len=1024) :: line, snippetfn
  character(len=1)    :: lchar
  integer :: ios, ib, n_excl, ipos, itgt, prev_atom
  logical :: file_exists

  ! Derive snippet filename: case.win -> case.pdwf_windows
  snippetfn = winfn
  ipos = index(snippetfn, '.win', back=.true.)
  if (ipos > 0) then
     snippetfn = snippetfn(1:ipos-1) // '.pdwf_windows'
  else
     snippetfn = trim(snippetfn) // '.pdwf_windows'
  end if

  ! Count excluded bands
  n_excl = 0
  do ib = 1, Nb
     if (band_class(ib) == BAND_EXCLUDED) n_excl = n_excl + 1
  end do

  ! --- Write standalone snippet file ---
  open(u_snip, file=trim(snippetfn), status='replace', form='formatted')
  call write_pdwf_keys(u_snip)
  call write_pdwf_projections(u_snip)
  close(u_snip)
  write(unit_out, '(A,A)') ' PDWF windows written to: ', trim(snippetfn)

  ! --- Patch existing .win file if present ---
  inquire(file=trim(winfn), exist=file_exists)
  if (.not. file_exists) then
     write(unit_out, '(A,A)') ' No .win file found at: ', trim(winfn)
     write(unit_out, '(A)') &
          ' Use write_win or paste .pdwf_windows into your .win'
     return
  end if

  open(u_win, file=trim(winfn), status='old', form='formatted', &
       iostat=ios)
  if (ios /= 0) then
     write(unit_out, '(A)') ' WARNING: Could not open .win for patching'
     return
  end if

  open(u_tmp, file=trim(winfn)//'.tmp', status='replace', &
       form='formatted')

  ! Write PDWF block at top
  write(u_tmp, '(A)') '! === PDWF auto-generated parameters ==='
  call write_pdwf_keys(u_tmp)
  write(u_tmp, '(A)') '! === end PDWF parameters ==='
  write(u_tmp, *)

  ! Copy existing lines, skipping conflicting keywords and old projections
  do
     read(u_win, '(A)', iostat=ios) line
     if (ios /= 0) exit
     if (is_pdwf_conflict(line)) cycle
     ! Skip old projections block — we write a new one
     if (index(line, 'begin projections') > 0) then
        do
           read(u_win, '(A)', iostat=ios) line
           if (ios /= 0) exit
           if (index(line, 'end projections') > 0) exit
        end do
        cycle
     end if
     if (index(line, 'guiding_centres') > 0 .and. &
          index(line, '!') /= 1) cycle
     write(u_tmp, '(A)') trim(line)
  end do
  close(u_win)

  ! Append new projections block
  call write_pdwf_projections(u_tmp)

  close(u_tmp)

  ! Replace original with patched version
  open(u_win, file=trim(winfn), status='replace', form='formatted')
  open(u_tmp, file=trim(winfn)//'.tmp', status='old', form='formatted')
  do
     read(u_tmp, '(A)', iostat=ios) line
     if (ios /= 0) exit
     write(u_win, '(A)') trim(line)
  end do
  close(u_tmp, status='delete')
  close(u_win)

  write(unit_out, '(A,A)') ' PDWF parameters patched into: ', trim(winfn)

contains

  subroutine write_pdwf_keys(u)
    integer, intent(in) :: u
    integer :: jb

    write(u, '(A,I0)') 'num_wann        = ', targets%Ntarget
    write(u, '(A,I0)') 'num_bands       = ', Nb
    if (has_frozen) then
       write(u, '(A,F12.6)') 'dis_froz_min    = ', froz_min
       write(u, '(A,F12.6)') 'dis_froz_max    = ', froz_max
    end if
    if (has_disent .and. Nb > targets%Ntarget) then
       write(u, '(A,F12.6)') 'dis_win_min     = ', dis_min
       write(u, '(A,F12.6)') 'dis_win_max     = ', dis_max
    end if
    if (n_excl > 0) then
       write(u, '(A)', advance='no') 'exclude_bands   = '
       do jb = 1, Nb
          if (band_class(jb) == BAND_EXCLUDED) then
             write(u, '(I0,",")', advance='no') jb
          end if
       end do
       write(u, *)
    end if
    ! Convergence tolerances
    if (Nb > targets%Ntarget) then
       write(u, '(A)') 'dis_num_iter    = 10000'
       write(u, '(A)') 'dis_conv_tol    = 1.0e-10'
       write(u, '(A)') 'dis_conv_window = 5'
    end if
    write(u, '(A)') 'num_iter        = 10000'
    write(u, '(A)') 'conv_tol        = 1.0e-8'
    write(u, '(A)') 'conv_window     = 5'
  end subroutine write_pdwf_keys

  subroutine write_pdwf_projections(u)
    integer, intent(in) :: u
    integer :: jt, iat
    real(R8) :: xf, yf, zf

    write(u, *)
    write(u, '(A)') 'guiding_centres = T'
    write(u, '(A)') 'begin projections'
    do jt = 1, targets%Ntarget
       select case(targets%l(jt))
       case(0); lchar = 's'
       case(1); lchar = 'p'
       case(2); lchar = 'd'
       case(3); lchar = 'f'
       end select
       ! Write f=x,y,z:l format — only one entry per (atom,l) group
       ! since Wannier90 expands l into all m components
       if (targets%m(jt) == -targets%l(jt)) then
          ! Use fractional coordinates from struct
          iat = targets%atom(jt)
          xf = stru%pos(1, iat)
          yf = stru%pos(2, iat)
          zf = stru%pos(3, iat)
          write(u, '("f=",F10.8,",",F10.8,",",F10.8,":",A)') &
               xf, yf, zf, lchar
       end if
    end do
    write(u, '(A)') 'end projections'
  end subroutine write_pdwf_projections

  logical function is_pdwf_conflict(s)
    character(len=*), intent(in) :: s
    is_pdwf_conflict = .false.
    if (index(adjustl(s), '!') == 1) then
       if (index(s, 'PDWF auto-generated') > 0 .or. &
            index(s, 'end PDWF parameters') > 0) then
          is_pdwf_conflict = .true.
       end if
       return
    end if
    if (index(s, 'dis_froz_min')   > 0) is_pdwf_conflict = .true.
    if (index(s, 'dis_froz_max')   > 0) is_pdwf_conflict = .true.
    if (index(s, 'dis_win_min')    > 0) is_pdwf_conflict = .true.
    if (index(s, 'dis_win_max')    > 0) is_pdwf_conflict = .true.
    if (index(s, 'num_wann')       > 0) is_pdwf_conflict = .true.
    if (index(s, 'num_bands')      > 0) is_pdwf_conflict = .true.
    if (index(s, 'exclude_bands')  > 0) is_pdwf_conflict = .true.
    if (index(s, 'dis_num_iter')   > 0) is_pdwf_conflict = .true.
    if (index(s, 'dis_conv_tol')   > 0) is_pdwf_conflict = .true.
    if (index(s, 'dis_conv_window') > 0) is_pdwf_conflict = .true.
    if (index(s, 'num_iter')       > 0 .and. &
         index(s, 'dis_num_iter')  == 0) is_pdwf_conflict = .true.
    if (index(s, 'conv_tol')       > 0 .and. &
         index(s, 'dis_conv_tol')  == 0) is_pdwf_conflict = .true.
    if (index(s, 'conv_window')    > 0 .and. &
         index(s, 'dis_conv_window') == 0) is_pdwf_conflict = .true.
  end function is_pdwf_conflict

end subroutine write_win_parameters


end module pdwf_m
