!!! hybrid_files.F -- readers for the standard Wannier90 files wien2wannier
!!! emits (.amn/.mmn/.eig/.nnkp), for the self-contained standalone driver.
!!! Conventions match the Python bridge readers (spread.read_amn/mmn/eig,
!!! bridge.read_nnkp_full).  All arrays complex(dp)/real(dp), column-major.
module hybrid_files_m
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: read_amn, read_mmn, read_eig, read_nnkp, read_nnkp_kgeom, nnkp_block_bvectors
  public :: read_w90_chk_subspace

  ! Binary sidecars (<file>.bin, unformatted stream): the 3.3 GB Rh3Si4 .mmn
  ! parses at ~2.5 MB/s (22 min cold); the sidecar reads in seconds.  A
  ! sidecar is used only if it is at least as new as its text file and its
  ! magic+dims validate; otherwise the text is parsed and the sidecar
  ! (re)written best-effort.  Delete the .bin to force a re-parse.
  integer, parameter :: SIDECAR_VER = 1

contains

  !> .true. when `path` is a SYMLINK.  Fortran's status='replace' FOLLOWS a
  !> symlink and TRUNCATES its target, so writing a cache sidecar through one
  !> silently rewrites whatever the link points at.  Every test harness in this
  !> tree symlinks archived per-cell .bin files into a scratch run directory, so
  !> a sidecar that fails its magic/version/dims check -- a version bump, a
  !> truncated file -- would rewrite the ARCHIVE.  Proven by controlled
  !> experiment: the link stayed a link and the target went from its original
  !> contents to the new record.  The cache is an optimisation; declining to
  !> write it costs a re-read, so refuse rather than risk the archive.
  logical function is_symlink(path) result(yes)
    character(len=*), intent(in) :: path
    integer :: st, sb(13)
    yes = .false.
    call lstat(trim(path), sb, st)
    if (st /= 0) return
    ! st_mode is element 3 of gfortran's stat buffer; S_IFLNK = 0o120000
    yes = (iand(sb(3), 61440) == 40960)
  end function is_symlink

  ! .true. if `bin` exists and is not older than `txt`
  logical function sidecar_fresh(txt, bin)
    character(*), intent(in) :: txt, bin
    integer :: st(13), sb(13), ios1, ios2
    call stat(txt, st, ios1)
    call stat(bin, sb, ios2)
    ! the sidecar is usable if it exists AND is not older than the text.
    ! A MISSING text file is not staleness: the .bin then holds the only
    ! copy, so use it (lets huge .amn/.mmn text files be reclaimed --
    ! Rh3Si4's mmn text is 12 GB against a 5.3 GB sidecar -- without
    ! disabling the fast path).
    sidecar_fresh = (ios2 == 0 .and. (ios1 /= 0 .or. sb(10) >= st(10)))
  end function sidecar_fresh

  subroutine read_amn(path, amn, nb, nk, nproj)
    character(*), intent(in) :: path
    complex(dp), allocatable, intent(out) :: amn(:,:,:)   ! (nb,nproj,nk)
    integer, intent(out) :: nb, nk, nproj
    integer :: u, m, n, k, i, ntot, ios
    real(dp) :: re, im
    character(len=4) :: magic
    if (sidecar_fresh(path, trim(path)//'.bin')) then
       open(newunit=u, file=trim(path)//'.bin', form='unformatted', &
            access='stream', status='old', action='read', iostat=ios)
       if (ios == 0) then
          read(u, iostat=ios) magic, i, nb, nk, nproj
          if (ios == 0 .and. magic == 'AMNB' .and. i == SIDECAR_VER) then
             allocate(amn(nb,nproj,nk))
             read(u, iostat=ios) amn
             close(u)
             if (ios == 0) then
                write(*,'(a)') '  read '//trim(path)//' via binary sidecar'
                return
             end if
             deallocate(amn)
          else
             close(u)
          end if
       end if
    end if
    open(newunit=u, file=path, status='old', action='read')
    read(u,*)                            ! title
    read(u,*) nb, nk, nproj
    allocate(amn(nb,nproj,nk)); amn = (0.0_dp,0.0_dp)
    ntot = nb*nk*nproj
    do i = 1, ntot
       read(u,*) m, n, k, re, im
       amn(m,n,k) = cmplx(re, im, dp)
    end do
    close(u)
    ios = -1
    if (.not. is_symlink(trim(path)//'.bin')) &
         open(newunit=u, file=trim(path)//'.bin', form='unformatted', &
              access='stream', status='replace', action='write', iostat=ios)
    if (ios == 0) then
       write(u) 'AMNB', SIDECAR_VER, nb, nk, nproj
       write(u) amn
       close(u)
    end if
  end subroutine read_amn

  subroutine read_mmn(path, mmn, neigh, gshift, nb, nk, nntot)
    character(*), intent(in) :: path
    complex(dp), allocatable, intent(out) :: mmn(:,:,:,:)  ! (nb,nb,nntot,nk)
    integer, allocatable, intent(out) :: neigh(:,:), gshift(:,:,:)
    integer, intent(out) :: nb, nk, nntot
    integer :: u, k, b, m, n, kk, kpb, g1, g2, g3, ip, npair, ios
    real(dp) :: re, im
    character(len=4) :: magic
    if (sidecar_fresh(path, trim(path)//'.bin')) then
       open(newunit=u, file=trim(path)//'.bin', form='unformatted', &
            access='stream', status='old', action='read', iostat=ios)
       if (ios == 0) then
          read(u, iostat=ios) magic, ip, nb, nk, nntot
          if (ios == 0 .and. magic == 'MMNB' .and. ip == SIDECAR_VER) then
             allocate(mmn(nb,nb,nntot,nk), neigh(nntot,nk), gshift(3,nntot,nk))
             read(u, iostat=ios) neigh, gshift, mmn
             close(u)
             if (ios == 0) then
                write(*,'(a)') '  read '//trim(path)//' via binary sidecar'
                return
             end if
             deallocate(mmn, neigh, gshift)
          else
             close(u)
          end if
       end if
    end if
    open(newunit=u, file=path, status='old', action='read')
    read(u,*)                            ! title
    read(u,*) nb, nk, nntot
    allocate(mmn(nb,nb,nntot,nk), neigh(nntot,nk), gshift(3,nntot,nk))
    npair = nk*nntot
    do ip = 1, npair
       read(u,*) kk, kpb, g1, g2, g3
       k = (ip-1)/nntot + 1
       b = mod(ip-1, nntot) + 1
       neigh(b,k) = kpb
       gshift(1,b,k) = g1; gshift(2,b,k) = g2; gshift(3,b,k) = g3
       do n = 1, nb          ! column-major block: m fastest
          do m = 1, nb
             read(u,*) re, im
             mmn(m,n,b,k) = cmplx(re, im, dp)
          end do
       end do
    end do
    close(u)
    ios = -1
    if (.not. is_symlink(trim(path)//'.bin')) &
         open(newunit=u, file=trim(path)//'.bin', form='unformatted', &
              access='stream', status='replace', action='write', iostat=ios)
    if (ios == 0) then
       write(u) 'MMNB', SIDECAR_VER, nb, nk, nntot
       write(u) neigh, gshift, mmn
       close(u)
    end if
  end subroutine read_mmn

  !> The disentangled subspace of a Wannier90 checkpoint (unformatted, the 3.x and 4.x
  !> record layout of readwrite.F90): u_matrix_opt placed on the outer-window bands
  !> (lwindow) at each k, V(nb,nw,nk), zero elsewhere.  have_dis = .false. (V not
  !> allocated) when the run did not disentangle, i.e. num_bands = num_wann.
  subroutine read_w90_chk_subspace(path, nb, nk, nw, have_dis, V, omega_i, kpts)
    character(*), intent(in) :: path
    integer, intent(out) :: nb, nk, nw
    logical, intent(out) :: have_dis
    complex(dp), allocatable, intent(out) :: V(:,:,:)
    real(dp), intent(out) :: omega_i
    real(dp), allocatable, intent(out) :: kpts(:,:)
    character(len=33) :: header
    character(len=20) :: label
    integer :: u, nex, nntot, mp(3), i, k, n
    integer, allocatable :: excl(:), ndimwin(:)
    logical, allocatable :: lwin(:,:)
    real(dp) :: latt(3,3)
    complex(dp), allocatable :: uopt(:,:,:)
    open(newunit=u, file=path, form='unformatted', status='old', action='read')
    read(u) header
    read(u) nb
    read(u) nex
    allocate(excl(max(nex,1)))
    read(u) (excl(i), i=1, nex)
    read(u) latt
    read(u) latt
    read(u) nk
    read(u) mp
    allocate(kpts(3,nk))
    read(u) kpts
    read(u) nntot
    read(u) nw
    read(u) label
    read(u) have_dis
    omega_i = 0.0_dp
    if (have_dis) then
       read(u) omega_i
       allocate(lwin(nb,nk), ndimwin(nk), uopt(nb,nw,nk))
       read(u) lwin
       read(u) ndimwin
       read(u) uopt
       allocate(V(nb,nw,nk)); V = (0.0_dp, 0.0_dp)
       do k = 1, nk
          n = 0
          do i = 1, nb
             if (lwin(i,k)) then
                n = n + 1
                V(i,:,k) = uopt(n,:,k)
             end if
          end do
       end do
    end if
    close(u)
  end subroutine read_w90_chk_subspace

  subroutine read_eig(path, nb, nk, eig)
    character(*), intent(in) :: path
    integer, intent(in) :: nb, nk
    real(dp), allocatable, intent(out) :: eig(:,:)   ! (nb,nk)
    integer :: u, n, k, i
    real(dp) :: e
    allocate(eig(nb,nk))
    open(newunit=u, file=path, status='old', action='read')
    do i = 1, nb*nk
       read(u,*) n, k, e
       eig(n,k) = e
    end do
    close(u)
  end subroutine read_eig

  ! read recip lattice, kpoints, neighbour list; return b-vectors for k=1
  subroutine read_nnkp(path, nk, nntot, bcart)
    character(*), intent(in) :: path
    integer, intent(out) :: nk, nntot
    real(dp), allocatable, intent(out) :: bcart(:,:)   ! (3,nntot)
    real(dp) :: recip(3,3)
    real(dp), allocatable :: kpts(:,:)
    integer, allocatable :: nn_kpb(:), nn_g(:,:)
    integer :: u, i, j, k, kpb, g1, g2, g3, ios
    character(len=256) :: line
    real(dp) :: dk(3)
    open(newunit=u, file=path, status='old', action='read')
    ! recip_lattice
    call seek(u, 'begin recip_lattice')
    do i = 1, 3; read(u,*) recip(i,1), recip(i,2), recip(i,3); end do
    ! kpoints
    call seek(u, 'begin kpoints')
    read(u,*) nk
    allocate(kpts(3,nk))
    do k = 1, nk; read(u,*) kpts(1,k), kpts(2,k), kpts(3,k); end do
    ! nnkpts
    call seek(u, 'begin nnkpts')
    read(u,*) nntot
    allocate(nn_kpb(nntot), nn_g(3,nntot))
    ! read the first k's nntot lines (k=1); the b-set is shared on a regular mesh
    do i = 1, nntot
       read(u,*) k, kpb, g1, g2, g3
       nn_kpb(i) = kpb; nn_g(1,i) = g1; nn_g(2,i) = g2; nn_g(3,i) = g3
    end do
    close(u)
    allocate(bcart(3,nntot))
    do i = 1, nntot
       dk = kpts(:,nn_kpb(i)) + real(nn_g(:,i),dp) - kpts(:,1)
       ! bcart = dk . recip   (recip rows are the reciprocal vectors)
       do j = 1, 3
          bcart(j,i) = dk(1)*recip(1,j) + dk(2)*recip(2,j) + dk(3)*recip(3,j)
       end do
    end do
    deallocate(kpts, nn_kpb, nn_g)
  contains
    subroutine seek(unit, tag)
      integer, intent(in) :: unit
      character(*), intent(in) :: tag
      do
         read(unit,'(a)', iostat=ios) line
         if (ios /= 0) then; write(*,*) 'nnkp: missing ', tag; error stop; end if
         if (index(line, tag) > 0) exit
      end do
    end subroutine seek
  end subroutine read_nnkp

  !> The mesh k-points (fractional) and the reciprocal lattice
  !> (rows, as in the .nnkp) so a caller can compute the b-vector of EVERY (k,b)
  !> block: the per-k neighbour order of a .nnkp differs from k = 1's (measured:
  !> 511 of 512 k on MgB2, 95 of 96 on Rh2B -- within shells), so read_nnkp's
  !> bcart(:,b) is k = 1's only.
  subroutine read_nnkp_kgeom(path, nk, kpts, recip)
    character(*), intent(in) :: path
    integer, intent(out) :: nk
    real(dp), allocatable, intent(out) :: kpts(:,:)   ! (3,nk)
    real(dp), intent(out) :: recip(3,3)
    integer :: u, i, k, ios
    character(len=256) :: line
    open(newunit=u, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(*,*) 'nnkp: cannot open ', trim(path); error stop
    end if
    call seek2('begin recip_lattice')
    do i = 1, 3; read(u,*) recip(i,1), recip(i,2), recip(i,3); end do
    call seek2('begin kpoints')
    read(u,*) nk
    allocate(kpts(3,nk))
    do k = 1, nk; read(u,*) kpts(1,k), kpts(2,k), kpts(3,k); end do
    close(u)
  contains
    subroutine seek2(tag)
      character(*), intent(in) :: tag
      do
         read(u,'(a)', iostat=ios) line
         if (ios /= 0) then; write(*,*) 'nnkp: missing ', tag; error stop; end if
         if (index(line, tag) > 0) exit
      end do
    end subroutine seek2
  end subroutine read_nnkp_kgeom

  !> The b-vector of EVERY (k,b) overlap block.  read_nnkp
  !> returns k = 1's neighbour vectors, but a .nnkp lists each k's neighbours in
  !> its own order (measured: 511 of 512 k on MgB2-WANN, 95 of 96 on the Rh2B
  !> band-17 deck), so bcart(:,b) is the vector of block (k,b) only at k = 1.
  !>   bkb(:,b,k) = recip^T ( k_neigh(b,k) + G(b,k) - k )      [1/Angstrom]
  !> with the k-points and reciprocal rows of the same .nnkp.  Gates (ok =
  !> .false., msg set): every block's vector must equal SOME k = 1 vector to
  !> BV_TOL (a mesh/nnkp mismatch otherwise), and must lie in the SAME shell as
  !> k = 1's slot b (|bkb(:,b,k)| = |bcart(:,b)|): the per-slot MV weights
  !> w(b) that every consumer applies at every k are only right if the per-k
  !> reordering stays within shells.
  subroutine nnkp_block_bvectors(path, nk, nntot, neigh, gshift, bcart, bkb, ok, msg)
    character(*), intent(in)  :: path
    integer,      intent(in)  :: nk, nntot
    integer,      intent(in)  :: neigh(nntot, nk), gshift(3, nntot, nk)
    real(dp),     intent(in)  :: bcart(3, nntot)
    real(dp),     intent(out) :: bkb(3, nntot, nk)
    logical,      intent(out) :: ok
    character(*), intent(out) :: msg
    real(dp), parameter :: BV_TOL = 1.0e-6_dp
    real(dp), allocatable :: kpts(:,:)
    real(dp) :: recip(3,3), dk(3), dmin, dj
    integer :: nkf, k, b, j, i
    ok = .false.; msg = ''
    call read_nnkp_kgeom(path, nkf, kpts, recip)
    if (nkf /= nk) then
       write(msg,'(a,i0,a,i0)') 'the .nnkp has ', nkf, ' k-points, the overlaps ', nk
       return
    end if
    do k = 1, nk
       do b = 1, nntot
          dk = kpts(:, neigh(b,k)) + real(gshift(:, b, k), dp) - kpts(:, k)
          do i = 1, 3
             bkb(i, b, k) = dk(1)*recip(1,i) + dk(2)*recip(2,i) + dk(3)*recip(3,i)
          end do
          dmin = huge(1.0_dp)
          do j = 1, nntot
             dj = maxval(abs(bcart(:, j) - bkb(:, b, k)))
             dmin = min(dmin, dj)
          end do
          if (dmin > BV_TOL) then
             write(msg,'(a,i0,a,i0,a,es9.2,a)') 'block (k ', k, ', b ', b, &
                  ') has a b-vector outside the k = 1 set (distance ', dmin, ')'
             return
          end if
          if (abs(sqrt(sum(bkb(:,b,k)**2)) - sqrt(sum(bcart(:,b)**2))) > BV_TOL) then
             write(msg,'(a,i0,a,i0,a)') 'block (k ', k, ', b ', b, ') lies in a different shell '// &
                  'from k = 1''s slot -- the per-slot weights would be wrong'
             return
          end if
       end do
    end do
    ok = .true.
  end subroutine nnkp_block_bvectors

end module hybrid_files_m
