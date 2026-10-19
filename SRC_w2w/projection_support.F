!!! wien2wannier/SRC_w2w/projection_support.F
!!!
!!!    Per-state support reduction behind compute_projectability's
!!!    component_sum_max mode (D4).
!!!
!!!    WHY A SEPARATE FILE: modules_rc.F cannot be linked into a unit test
!!!    (compute_projectability pulls almgen/atpar/xa3/mpi_w2w and the whole
!!!    vector-file machinery), so the one reduction whose semantics were
!!!    wrong -- "max over (band, k, spin) of a per-spin weight" instead of
!!!    "max over (band, k) of the spinor SUM" -- is isolated here where a
!!!    unit test can exercise the PRODUCTION source against analytic inputs.  This module owns no SOC policy: the
!!!    caller decides what state_sum holds; this only reduces it.

module projection_support_m
  use const, only: R8

  implicit none
  private
  public :: state_support_max, state_support_max_trace

contains

  !! max_support(j, l) = max_{n,k} state_sum(n, k, j, l)
  !!
  !! Reduces over the band and k dimensions ONLY.  When state_sum is the
  !! component-summed per-state projectability P_1 + P_2, this is the
  !! spinor-sum channel support max_{n,k} [P_1(n,k) + P_2(n,k)] -- NOT
  !! max(P_1) + max(P_2) and NOT max over (n, k, s) of P_s(n,k).
  !! On a single-component array it reproduces the per-k streaming
  !! max(proj_max, stmp) of compute_projectability bit-for-bit: both are
  !! fed the same val in the same (mu, m) order from the same zero, and
  !! max is exact.
  pure subroutine state_support_max(state_sum, max_support)
    real(R8), intent(in)  :: state_sum(:,:,:,0:)   ! (Nb, Nk, Nneq, 0:Lmax2)
    real(R8), intent(out) :: max_support(:,0:)     ! (Nneq, 0:Lmax2)

    integer :: j, l

    do l = 0, ubound(state_sum, 4)
       do j = 1, size(state_sum, 3)
          max_support(j, l) = maxval(state_sum(:, :, j, l))
       end do
    end do
  end subroutine state_support_max

  !! max_support(j, l) = max_{k, c} (1/|c|) sum_{n in c} state_sum(n, k, j, l),
  !! where c ranges over the DEGENERATE CLUSTERS at k: contiguous runs of
  !! bands (in energy order) whose adjacent gap eb(n+1,k) - eb(n,k) <= degtol.
  !! This is the multiplet TRACE, normalized by cluster size -- the D5/D4-2
  !! fix for state_support_max's basis dependence inside an exactly
  !! degenerate multiplet: a rank-one support operator diag(.8,0) on a
  !! Kramers pair becomes (.4,.4) after a 45-degree rotation of the pair's
  !! gauge, so max() over the two states reads .8 or .4 for the SAME
  !! physics.  Two reasons to normalize by |c| rather than use the raw trace:
  !!   (i) trace/|c| is exactly the per-state support a TIME-REVERSAL-
  !!       SYMMETRIC gauge would assign to each Kramers partner -- TR maps
  !!       partner to partner and the spin-summed |Alm|^2 projector is
  !!       TR-invariant, so the pair's total is gauge-fixed and splits
  !!       equally between the partners in a TR-symmetric gauge;
  !!  (ii) it keeps the support spectrum on the same [0,1] per-state scale
  !!       the SUPPORT_RATIO / Otsu thresholds were calibrated on -- the raw
  !!       trace would double every SOC channel's support relative to a
  !!       scalar cell.
  !! On a nondegenerate state |c| = 1 and the value is unchanged (division by
  !! the exact integer 1.0 and max() are both exact in IEEE arithmetic), so
  !! on nondegenerate input this reduces to state_support_max EXACTLY.
  !!
  !! Clustering is the CONTIGUOUS-run partition (cut wherever an adjacent gap
  !! exceeds degtol) -- NOT main.F's
  !! max_multiplet "mutual window" (last - first < DEGTOL), which bounds a
  !! single worst-case multiplet size and is not well-defined as a partition
  !! (the mutual relation is not transitive).  Bands are assumed ascending in
  !! energy per k, as read_vec fills Eb and as every other consumer of
  !! Eb/eall in this tree assumes; not re-sorted here.
  !!
  !! ncl_min/ncl_max: smallest/largest cluster size seen, over all k.
  !! nclus: total number of clusters, summed over all k.  Both are properties
  !! of the energy clustering alone (independent of j, l), for the census
  !! line compute_projectability prints.  Optional gap_merged_max/
  !! gap_open_min: the largest adjacent gap merged into a cluster and the
  !! smallest gap left open between two clusters, over all k -- what a
  !! reader needs to judge whether degtol was a safe choice for this cell.
  pure subroutine state_support_max_trace(state_sum, eb, degtol, max_support, &
       ncl_min, ncl_max, nclus, gap_merged_max, gap_open_min)
    real(R8), intent(in)  :: state_sum(:,:,:,0:)  ! (Nb, Nk, Nneq, 0:Lmax2)
    real(R8), intent(in)  :: eb(:,:)              ! (Nb, Nk), Ry
    real(R8), intent(in)  :: degtol               ! Ry
    real(R8), intent(out) :: max_support(:,0:)    ! (Nneq, 0:Lmax2)
    integer,  intent(out) :: ncl_min, ncl_max, nclus
    real(R8), intent(out), optional :: gap_merged_max, gap_open_min

    integer  :: Nb, Nk, j, l, k, ib, i0, csize
    real(R8) :: gap, gmax_m, gmin_o

    Nb = size(state_sum, 1);  Nk = size(state_sum, 2)
    max_support = -huge(1.0_R8)
    ncl_min = huge(1);  ncl_max = 0;  nclus = 0
    gmax_m = -1.0_R8;   gmin_o = huge(1.0_R8)

    do k = 1, Nk
       i0 = 1
       do ib = 2, Nb+1
          if (ib <= Nb) gap = eb(ib, k) - eb(ib-1, k)
          if (ib > Nb .or. gap > degtol) then
             ! close the cluster [i0, ib-1]
             csize = ib - i0
             nclus   = nclus + 1
             ncl_min = min(ncl_min, csize)
             ncl_max = max(ncl_max, csize)
             do l = 0, ubound(state_sum, 4)
                do j = 1, size(state_sum, 3)
                   max_support(j, l) = max(max_support(j, l), &
                        sum(state_sum(i0:ib-1, k, j, l)) / real(csize, R8))
                end do
             end do
             if (ib <= Nb) gmin_o = min(gmin_o, gap)
             i0 = ib
          else
             gmax_m = max(gmax_m, gap)
          end if
       end do
    end do

    if (present(gap_merged_max)) gap_merged_max = gmax_m
    if (present(gap_open_min))   gap_open_min   = gmin_o
  end subroutine state_support_max_trace

end module projection_support_m
