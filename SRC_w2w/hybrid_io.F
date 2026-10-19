!!! hybrid_io.F -- Nb=Nw hand-off writers for the hybrid Stage-2 port.
!!! Emits stock-Wannier90 .amn/.mmn/.eig (isolated problem, no dis_* keywords)
!!! in the exact formats the Python bridge (hybrid_pipeline._write_*) produces,
!!! plus the hand-off contractivity (unitarity) check.  Reused by both the
!!! standalone driver and the in-w2w hook.  File formats: Wannier90
!!! [Mostofi2008, Pizzi2020] (keys of REFERENCES.md in the root directory).
module hybrid_io_m
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: write_handoff, handoff_svmax, validate_handoff, write_masks, write_anatomy, &
            write_window_route, handoff_close, subspace_continuity, write_dis_status, &
            SUBSPACE_JUMP_SV, continuity_links

  !> Below this smallest singular value of a hand-off overlap block W(k)^H M(k,b) W(k+b)
  !> the disentangled subspace is treated as jumping between neighbouring k-points.
  !> Calibrated on the LCAO frontend (lcao2wannier, 2026-10-01): usable models 0.012-0.93,
  !> poor ones 0.0062 and 0.0007 (SnTe, 90 meV).
  real(dp), parameter :: SUBSPACE_JUMP_SV = 1.0e-2_dp

contains

  !> Closes a hand-off file and stops the run (error stop) when the size on disk differs from the size written.
  !> gfortran reports no error for buffered writes that hit a full disk -- not on the write, the flush or the
  !> close, with or without iostat -- so a short file would otherwise pass as a complete hand-off.  The same
  !> test as clio's close_checked in the w2w driver; this copy lives here because the standalone w2w_hybrid
  !> does not link clio.
  subroutine handoff_close(u, what)
    use, intrinsic :: iso_fortran_env, only: int64, error_unit
    integer,      intent(in) :: u
    character(*), intent(in) :: what

    character(len=4096) :: fname
    integer(int64)      :: szw, szd
    integer             :: ios
    logical             :: named

    fname = ''; szw = -1; szd = -1
    inquire(unit=u, named=named, name=fname, size=szw)
    close(u, iostat=ios)
    if (named) inquire(file=trim(fname), size=szd)
    if (ios /= 0 .or. (named .and. szw >= 0 .and. szd /= szw)) then
       write(error_unit, '(a,a,a,a,i0,a,i0,a,i0,a)') trim(what), ' ', trim(fname), ': ', szd, &
            ' bytes on disk but ', szw, ' were written (close iostat ', ios, &
            '). The write was lost, most likely because the disk is full.'
       error stop 'hand-off file incomplete (disk full?). Refusing the hand-off. Free disk space and rerun.'
    end if
  end subroutine handoff_close

  ! max singular value over EVERY (k,b) block of mmn_rot -- full coverage,
  ! no sampling: final acceptance must cover every link (a violation at an
  ! unsampled k is exactly the silent failure this exists to stop). The
  ! hand-off is contractive iff this is <= 1 + 1e-6 (paper criterion).
  ! Error-stops on SVD failure or a non-finite singular value, so a NaN
  ! block can never pass by comparing false against the running max.
  function handoff_svmax(mmn_rot, nw, nntot, nk) result(svmax)
    integer, intent(in) :: nw, nntot, nk
    complex(dp), intent(in) :: mmn_rot(nw,nw,nntot,nk)
    real(dp) :: svmax
    complex(dp), allocatable :: Ac(:,:), work(:)
    complex(dp) :: qw(1), udum(1,1), vtdum(1,1)
    real(dp), allocatable :: s(:), rwork(:)
    integer, allocatable :: iwork(:)
    integer :: k, b, lwork, info
    allocate(Ac(nw,nw), s(nw), rwork(7*nw), iwork(8*nw))
    svmax = 0.0_dp
    do k = 1, nk
       do b = 1, nntot
          Ac = mmn_rot(:,:,b,k)
          if (k == 1 .and. b == 1) then
             call zgesdd('N', nw, nw, Ac, nw, s, udum,1, vtdum,1, qw,-1, rwork, iwork, info)
             lwork = int(real(qw(1),dp)); allocate(work(lwork))
          end if
          Ac = mmn_rot(:,:,b,k)
          call zgesdd('N', nw, nw, Ac, nw, s, udum,1, vtdum,1, work,lwork, rwork, iwork, info)
          if (info /= 0) then
             write(*,'(a,i0,a,i0,a,i0)') '  *** handoff_svmax: zgesdd failed (info = ', &
                  info, ') at k = ', k, ', b = ', b
             error stop 'handoff_svmax: SVD failed on a hand-off overlap block'
          end if
          if (s(1) /= s(1)) then
             write(*,'(a,i0,a,i0)') '  *** handoff_svmax: non-finite singular value at k = ', k, ', b = ', b
             error stop 'handoff_svmax: NaN in a hand-off overlap block'
          end if
          if (s(1) > svmax) svmax = s(1)
       end do
    end do
    deallocate(Ac, s, rwork, iwork, work)
  end function handoff_svmax

  !> Subspace continuity of the hand-off: the smallest singular value of every (k,b)
  !> overlap block of the rotated subspace, W(k)^H M(k,b) W(k+b).  For a subspace that
  !> varies smoothly in k all of them are close to 1 (a finite b lowers them a little);
  !> a value near 0 means a state of the subspace at k has no partner in the subspace at
  !> k+b: the subspace jumps there.  No gauge can smooth a jump: Wannier90 then
  !> interpolates worse between those k-points and converges slowly or not at all
  !> (LCAO SnTe: 7e-4, 90 meV, unconverged at 300,000 iterations; LAPW GaAs at the
  !> release protocol: 5e-4, 4.9/18.0 meV, unconverged at 50,000).  A diagnostic
  !> only: it changes nothing written.  Prints one line, and a warning below
  !> SUBSPACE_JUMP_SV, on unit uo.
  subroutine subspace_continuity(mmn_rot, neigh, nw, nntot, nk, uo, smin, kmin, kbmin)
    integer, intent(in) :: nw, nntot, nk, uo
    complex(dp), intent(in) :: mmn_rot(nw,nw,nntot,nk)
    integer, intent(in) :: neigh(nntot,nk)
    real(dp), intent(out) :: smin
    integer, intent(out) :: kmin, kbmin
    complex(dp), allocatable :: Ac(:,:), work(:)
    complex(dp) :: qw(1), udum(1,1), vtdum(1,1)
    real(dp), allocatable :: s(:), rwork(:)
    integer, allocatable :: iwork(:)
    integer :: k, b, lwork, info, nhalf
    allocate(Ac(nw,nw), s(nw), rwork(7*nw), iwork(8*nw))
    Ac = (0.0_dp, 0.0_dp)
    call zgesdd('N', nw, nw, Ac, nw, s, udum,1, vtdum,1, qw,-1, rwork, iwork, info)
    lwork = max(1, int(real(qw(1),dp))); allocate(work(lwork))
    smin = huge(1.0_dp); kmin = 0; kbmin = 0; nhalf = 0
    do k = 1, nk
       do b = 1, nntot
          Ac = mmn_rot(:,:,b,k)
          call zgesdd('N', nw, nw, Ac, nw, s, udum,1, vtdum,1, work,lwork, rwork, iwork, info)
          if (info /= 0) then
             write(uo,'(a,i0,a,i0,a,i0)') '  *** subspace_continuity: zgesdd failed (info = ', &
                  info, ') at k = ', k, ', b = ', b
             error stop 'subspace_continuity: SVD failed on a hand-off overlap block'
          end if
          if (s(nw) < 0.5_dp) nhalf = nhalf + 1
          if (s(nw) < smin) then
             smin = s(nw); kmin = k; kbmin = neigh(b,k)
          end if
       end do
    end do
    deallocate(Ac, s, rwork, iwork, work)
    write(uo,'(a,es10.3,a,i0,a,i0,a,i0,a,i0,a,es8.1,a)') '  subspace continuity: smallest '// &
         'singular value of the hand-off overlaps ', smin, ' at k ', kmin, ' -> ', kbmin, '; ', &
         nhalf, ' of ', nk*nntot, ' (k,b) blocks below 0.5 (warning below ', SUBSPACE_JUMP_SV, ')'
    if (smin < SUBSPACE_JUMP_SV) then
       write(uo,'(a,i0,a,i0,a,es8.1,a)') '  *** WARNING: the disentangled subspace JUMPS between k ', &
            kmin, ' and k ', kbmin, ' (singular value ', smin, ').'
       write(uo,'(a)') '  *** Wannier90 cannot smooth a jump: expect larger interpolation errors between '// &
            'those k-points and a wannierisation that converges slowly or not at all.'
       write(uo,'(a)') '  *** Try the CG solver (w2w_hybrid --solver cg; HYBR keyword solver = cg), which '// &
            'can reach a different minimum, or a different frozen ceiling (--froz-emax; froz_emax =).'
    end if
  end subroutine subspace_continuity

  !> The smallest singular value of every hand-off link W^H(k) M W(k+b), sv(nntot,nk),
  !> without printing (subspace_continuity reports the minimum; the pocket skip
  !> needs every link).
  subroutine continuity_links(mmn_rot, nw, nntot, nk, sv)
    integer, intent(in) :: nw, nntot, nk
    complex(dp), intent(in) :: mmn_rot(nw,nw,nntot,nk)
    real(dp), intent(out) :: sv(nntot,nk)
    complex(dp), allocatable :: Ac(:,:), work(:)
    complex(dp) :: qw(1), udum(1,1), vtdum(1,1)
    real(dp), allocatable :: s(:), rwork(:)
    integer, allocatable :: iwork(:)
    integer :: k, b, lwork, info
    allocate(Ac(nw,nw), s(nw), rwork(7*nw), iwork(8*nw))
    Ac = (0.0_dp, 0.0_dp)
    call zgesdd('N', nw, nw, Ac, nw, s, udum,1, vtdum,1, qw,-1, rwork, iwork, info)
    lwork = max(1, int(real(qw(1),dp))); allocate(work(lwork))
    do k = 1, nk
       do b = 1, nntot
          Ac = mmn_rot(:,:,b,k)
          call zgesdd('N', nw, nw, Ac, nw, s, udum,1, vtdum,1, work,lwork, rwork, iwork, info)
          if (info /= 0) error stop 'continuity_links: SVD failed on a hand-off overlap block'
          sv(b,k) = s(nw)
       end do
    end do
    deallocate(Ac, s, rwork, iwork, work)
  end subroutine continuity_links

  !> <out>.dis: the disentanglement outcome as comment lines for the hand-off .win
  !> (the wien2wannier launcher copies them there), so the record travels with the
  !> files Wannier90 reads.
  subroutine write_dis_status(fname, status_line, continuity_line)
    character(*), intent(in) :: fname, status_line, continuity_line
    integer :: u
    open(newunit=u, file=fname, status='replace')
    write(u,'(a)') '! wien2wannier disentanglement: '//trim(status_line)
    write(u,'(a)') '! wien2wannier subspace continuity: '//trim(continuity_line)
    call handoff_close(u, 'the .dis status file')
  end subroutine write_dis_status

  ! Complete-hand-off acceptance: finite values everywhere, SVD success,
  ! and sigma_max <= 1 + 1e-6 over EVERY (k,b) block. Refuses by error
  ! stop -- overlaps are never clipped or unitarized to pass. The 1e-6
  ! tolerance covers SVD roundoff plus the .amn/.mmn text serialization
  ! (18.12f, relative ~1e-12 per element -- orders below the tolerance).
  subroutine validate_handoff(eig_t, A_gauge, mmn_rot, nw, nk, nntot)
    integer, intent(in) :: nw, nk, nntot
    real(dp), intent(in) :: eig_t(nw,nk)
    complex(dp), intent(in) :: A_gauge(nw,nw,nk), mmn_rot(nw,nw,nntot,nk)
    real(dp) :: svmax
    if (any(eig_t /= eig_t) .or. any(abs(eig_t) > 1.0e30_dp)) &
       error stop 'validate_handoff: non-finite eigenvalue in the hand-off .eig'
    if (any(A_gauge /= A_gauge) .or. any(abs(A_gauge) > 1.0e30_dp)) &
       error stop 'validate_handoff: non-finite entry in the hand-off gauge (.amn)'
    if (any(mmn_rot /= mmn_rot)) &
       error stop 'validate_handoff: non-finite entry in the hand-off overlaps (.mmn)'
    svmax = handoff_svmax(mmn_rot, nw, nntot, nk)
    if (svmax > 1.0_dp + 1.0e-6_dp) then
       write(*,'(a,f12.8)') '  *** hand-off max sigma(M) = ', svmax
       error stop 'validate_handoff: hand-off is NOT contractive '// &
            '(sigma_max > 1 + 1e-6); nothing was written -- fix the '// &
            'gauge/rotation that produced it, do not clip the overlaps'
    end if
    write(*,'(a,f12.6,a)') '  hand-off acceptance: max sigma(M) = ', svmax, &
         '  (full coverage, every (k,b) block; tolerance 1 + 1e-6)'
  end subroutine validate_handoff

  subroutine write_masks(fname, pool, npool, frozen, admit, P, eig, e_fermi, nb, nk, nw, &
                         p_froz, p_gate, fe_used)
    ! <out>.masks: the per-state masks the hand-off was built from, one row
    ! per (k, pool state). Read by tools/plot_certificate.py (window
    ! non-representability certificate) and by any script that wants to
    ! reconstruct F(k). Columns: k  band  eps-E_F(eV)  P  frozen  admitted
    character(*), intent(in) :: fname
    integer, intent(in) :: npool, nb, nk, nw, pool(npool)
    logical, intent(in) :: frozen(npool,nk), admit(npool,nk)
    real(dp), intent(in) :: P(nb,nk), eig(nb,nk), e_fermi, p_froz, p_gate, fe_used
    integer :: u, k, i
    open(newunit=u, file=fname, status='replace')
    write(u,'(a)') '# w2w hybrid per-state masks v1 (parsers skip # lines)'
    write(u,'(a,i0,a,i0,a,i0,a,i0)') '# nk=', nk, ' npool=', npool, ' nw=', nw, ' nb=', nb
    write(u,'(a,f14.8,a,f8.4,a,f8.4,a,f8.4)') '# e_fermi_eV=', e_fermi, &
         ' p_froz=', p_froz, ' p_gate=', p_gate, ' froz_emax_eV=', fe_used
    write(u,'(a)') '# pool bands (1-based, into the seed .eig band set):'
    write(u,'(a,*(1x,i0))') '#', (pool(i), i = 1, npool)
    write(u,'(a)') '# columns: k  band  eps_rel_EF_eV  P  frozen  admitted'
    do k = 1, nk
       do i = 1, npool
          write(u,'(i6,1x,i5,1x,f14.8,1x,f10.6,1x,i1,1x,i1)') k, pool(i), &
               eig(pool(i),k) - e_fermi, P(pool(i),k), &
               merge(1, 0, frozen(i,k)), merge(1, 0, admit(i,k))
       end do
    end do
    call handoff_close(u, 'the .masks file')
  end subroutine write_masks

  subroutine write_anatomy(fname, pool, npool, avg, minp, bot, nb, nw, p_gate, p_froz, &
                           plateau, fid_emax, ok)
    ! <out>.anatomy: the auto-trust rule made visible -- one row per pool
    ! band with its band-averaged and minimum projectability, its bottom
    ! relative to E_F, whether it is fidelity-relevant, and the class the
    ! rule assigned (trusted / gray / junk). Same schema as the LCAO
    ! --hybrid-report-json; read by tools/plot_anatomy.py.
    character(*), intent(in) :: fname
    integer, intent(in) :: npool, nb, nw, pool(npool)
    real(dp), intent(in) :: avg(nb), minp(nb), bot(nb), p_gate, p_froz, plateau, fid_emax
    logical, intent(in) :: ok
    integer :: u, ip, i
    logical :: rel
    character(len=8) :: cls
    open(newunit=u, file=fname, status='replace')
    write(u,'(a)') '# w2w hybrid auto-trust anatomy v1 (parsers skip # lines)'
    if (ok) then
       write(u,'(a,f8.4,a,f8.4,a,f8.4,a,f6.2,a,i0)') '# p_gate=', p_gate, &
            ' p_froz=', p_froz, ' plateau=', plateau, ' fid_emax_eV=', fid_emax, ' nw=', nw
    else
       write(u,'(a,f8.4,a,f8.4,a,f6.2,a,i0)') '# p_gate=', p_gate, &
            ' p_froz=', p_froz, ' plateau=none (auto rule undecidable; fixed defaults used)' // &
            ' fid_emax_eV=', fid_emax, ' nw=', nw
    end if
    write(u,'(a)') '# columns: band  Pbar  minP  minE_rel_EF_eV  relevant  cls'
    do ip = 1, npool
       i = pool(ip)
       rel = bot(i) <= fid_emax
       if (rel .and. avg(i) >= p_gate) then
          cls = 'trusted'
       else if (avg(i) >= p_froz) then
          cls = 'gray'
       else
          cls = 'junk'
       end if
       write(u,'(i5,1x,f10.6,1x,f10.6,1x,f12.6,1x,i1,1x,a)') i, avg(i), minp(i), bot(i), &
            merge(1, 0, rel), trim(cls)
    end do
    call handoff_close(u, 'the .anatomy file')
  end subroutine write_anatomy

  subroutine write_window_route(outp, eig, nb, nk, nw, e_fermi, amn, nproj, anchors, &
                                mmn, neigh, gshift, nntot, &
                                fmin, fmax, dmin, dmax, ncand, cand_d, cand_om)
    ! <out>-win.{win,amn,mmn,eig}: the same automatic selection as the
    ! hand-off, expressed as Wannier90 disentanglement windows for the stock
    ! window route. The .amn holds the anchor columns of the seed projections
    ! over ALL nb bands (num_proj = num_wann, as Wannier90 requires), the
    ! .mmn the seed's own overlaps, the .eig every band relative to E_F.
    ! The .win is a snippet: windows, num_bands, num_wann and the Omega_I of
    ! every ceiling tried.
    character(*), intent(in) :: outp
    integer, intent(in) :: nb, nk, nw, nproj, ncand, nntot, anchors(nw)
    integer, intent(in) :: neigh(nntot,nk), gshift(3,nntot,nk)
    real(dp), intent(in) :: eig(nb,nk), e_fermi, fmin, fmax, dmin, dmax
    real(dp), intent(in) :: cand_d(ncand), cand_om(ncand)
    complex(dp), intent(in) :: amn(nb,nproj,nk), mmn(nb,nb,nntot,nk)
    integer :: u, k, n, m, i, b
    open(newunit=u, file=trim(outp)//'-win.win', status='replace')
    write(u,'(a)') '! automatic window route (w2w_hybrid --window): the hand-off''s own'
    write(u,'(a)') '! selection as Wannier90 windows; dis_win_max chosen by Omega_I'
    do i = 1, ncand
       if (cand_om(i) < huge(1.0_dp)) then
          write(u,'(a,f10.4,a,f14.6)') '!   dis_win_max = ', cand_d(i), &
               ' eV  ->  Omega_I = ', cand_om(i)
       else
          write(u,'(a,f10.4,a)') '!   dis_win_max = ', cand_d(i), &
               ' eV  ->  rejected (frozen > num_wann or admitted < num_wann at a k-point)'
       end if
    end do
    write(u,'(a,i0)') 'num_bands = ', nb
    write(u,'(a,i0)') 'num_wann = ', nw
    write(u,'(a,f12.6)') 'dis_froz_min = ', fmin
    write(u,'(a,f12.6)') 'dis_froz_max = ', fmax
    write(u,'(a,f12.6)') 'dis_win_min = ', dmin
    write(u,'(a,f12.6)') 'dis_win_max = ', dmax
    call handoff_close(u, 'the window-route .win')
    open(newunit=u, file=trim(outp)//'-win.eig', status='replace')
    do k = 1, nk
       write(u,'(2i5,f18.12)') (n, k, eig(n,k) - e_fermi, n = 1, nb)
    end do
    call handoff_close(u, 'the window-route .eig')
    open(newunit=u, file=trim(outp)//'-win.amn', status='replace')
    write(u,'(a)') trim(outp)//' : window route, anchor columns of the seed projections'
    write(u,'(3i12)') nb, nk, nw
    do k = 1, nk
       write(u,'(3i5,2f18.12)') ((m, n, k, real(amn(m,anchors(n),k),dp), &
            aimag(amn(m,anchors(n),k)), m = 1, nb), n = 1, nw)
    end do
    call handoff_close(u, 'the window-route .amn')
    open(newunit=u, file=trim(outp)//'-win.mmn', status='replace')
    write(u,'(a)') trim(outp)//' : window route, seed overlaps'
    write(u,'(3i12)') nb, nk, nntot
    do k = 1, nk
       do b = 1, nntot
          write(u,'(5i5)') k, neigh(b,k), gshift(1,b,k), gshift(2,b,k), gshift(3,b,k)
          write(u,'(f18.12,1x,f18.12)') ((real(mmn(m,n,b,k),dp), &
               aimag(mmn(m,n,b,k)), m = 1, nb), n = 1, nb)
       end do
    end do
    call handoff_close(u, 'the window-route .mmn')
  end subroutine write_window_route

  subroutine write_eig(fname, eig_t, e_fermi, nw, nk)
    character(*), intent(in) :: fname
    integer, intent(in) :: nw, nk
    real(dp), intent(in) :: eig_t(nw,nk), e_fermi
    integer :: u, k, n
    open(newunit=u, file=fname, status='replace')
    ! one WRITE per k; format reversion emits nw records of '(2i5,f18.12)',
    ! byte-identical to the per-element loop it replaces
    do k = 1, nk
       write(u,'(2i5,f18.12)') (n, k, eig_t(n,k) - e_fermi, n = 1, nw)
    end do
    call handoff_close(u, 'the hand-off .eig')
  end subroutine write_eig

  subroutine write_amn(fname, title, A_gauge, nw, nk)
    character(*), intent(in) :: fname, title
    integer, intent(in) :: nw, nk
    complex(dp), intent(in) :: A_gauge(nw,nw,nk)   ! (band m, wann n, k)
    integer :: u, k, n, m
    open(newunit=u, file=fname, status='replace')
    write(u,'(a)') trim(title)//' : hybrid Nb=Nw hand-off'
    write(u,'(3i12)') nw, nk, nw
    ! one WRITE per k (nw*nw records by reversion); m fastest, as before
    do k = 1, nk
       write(u,'(3i5,2f18.12)') ((m, n, k, real(A_gauge(m,n,k),dp), &
            aimag(A_gauge(m,n,k)), m = 1, nw), n = 1, nw)
    end do
    call handoff_close(u, 'the hand-off .amn')
  end subroutine write_amn

  subroutine write_mmn(fname, title, mmn_rot, neigh, gshift, nw, nk, nntot)
    character(*), intent(in) :: fname, title
    integer, intent(in) :: nw, nk, nntot
    complex(dp), intent(in) :: mmn_rot(nw,nw,nntot,nk)
    integer, intent(in) :: neigh(nntot,nk), gshift(3,nntot,nk)
    integer :: u, k, b, n, m
    open(newunit=u, file=fname, status='replace')
    write(u,'(a)') trim(title)//' : hybrid Nb=Nw hand-off'
    write(u,'(3i12)') nw, nk, nntot
    do k = 1, nk
       do b = 1, nntot
          write(u,'(5i5)') k, neigh(b,k), gshift(1,b,k), gshift(2,b,k), gshift(3,b,k)
          ! one WRITE per (k,b) block; column-major, m fastest, as before
          write(u,'(f18.12,1x,f18.12)') ((real(mmn_rot(m,n,b,k),dp), &
               aimag(mmn_rot(m,n,b,k)), m = 1, nw), n = 1, nw)
       end do
    end do
    call handoff_close(u, 'the hand-off .mmn')
  end subroutine write_mmn

  subroutine write_handoff(seed, title, eig_t, A_gauge, mmn_rot, neigh, gshift, &
                           e_fermi, nw, nk, nntot)
    character(*), intent(in) :: seed, title
    integer, intent(in) :: nw, nk, nntot
    real(dp), intent(in) :: eig_t(nw,nk), e_fermi
    complex(dp), intent(in) :: A_gauge(nw,nw,nk), mmn_rot(nw,nw,nntot,nk)
    integer, intent(in) :: neigh(nntot,nk), gshift(3,nntot,nk)
    ! FINAL ACCEPTANCE before anything is written: a hand-off that fails
    ! here must not exist on disk at all (print-only checks after the write
    ! left a published, inconsistent hand-off behind).
    call validate_handoff(eig_t, A_gauge, mmn_rot, nw, nk, nntot)
    call write_eig(trim(seed)//'.eig', eig_t, e_fermi, nw, nk)
    call write_amn(trim(seed)//'.amn', title, A_gauge, nw, nk)
    call write_mmn(trim(seed)//'.mmn', title, mmn_rot, neigh, gshift, nw, nk, nntot)
  end subroutine write_handoff

end module hybrid_io_m
