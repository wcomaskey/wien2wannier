!!! wien2wannier/SRC_trig/w2waddsp.f
!!!
!!!    Combines amn,mmn files for up/dn spin from wien2k to one
!!!    amn,mmn by entrywise summation.
!!!
!!!    Rewritten on the hybrid_files_m readers (perf redesign, B4):
!!!    - inputs are read through read_amn/read_mmn, so the transparent
!!!      binary sidecars (<file>.bin) are used when fresh and rewritten
!!!      after a text parse -- the 3.3 GB Rh3Si4 text parse happens at
!!!      most once per spin file;
!!!    - the two spin files are VALIDATED against each other (dims,
!!!      neighbour lists, G-shifts); the old version discarded file-2's
!!!      headers unread and summed blind;
!!!    - the combined .amn keeps full double precision (2es24.15); the
!!!      old 2E18.5 floored the combined file to 5 significant digits,
!!!      undoing the upstream l2Amn precision fix exactly where the
!!!      hybrid pipeline reads it;
!!!    - combined outputs get .bin sidecars (text closed first, so the
!!!      sidecar mtime test stays valid).
!!!
!!! Copyright 2009-2012 Philipp Wissgott
!!!           2014-2015 Elias Assmann

PROGRAM combine_spinfiles
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use const,  only: bufsz
  use wien2k, only: gtfnam, errflg, errclr
  use hybrid_files_m, only: read_amn, read_mmn

  implicit none

  integer, parameter :: unit_mmn=13, unit_mmn1=11, unit_mmn2=12
  integer, parameter :: unit_amn=23, unit_amn1=21, unit_amn2=22
  integer, parameter :: SIDECAR_VER = 1

  character(BUFSZ) :: deffn, errfn, status, form, fname, title
  character(len=8) :: envtxt
  character(BUFSZ) :: fn(30)
  logical          :: seen(30)

  complex(dp), allocatable :: a1(:,:,:), a2(:,:,:)
  complex(dp), allocatable :: m1(:,:,:,:), m2(:,:,:,:)
  integer, allocatable :: neigh1(:,:), gsh1(:,:,:), neigh2(:,:), gsh2(:,:,:)
  integer :: nb1, nk1, nn1, nb2, nk2, nn2, ios
  integer :: iunit, ik, k, b, m, n, ipr, u

  call gtfnam(deffn, errfn, ik)
  call errflg(errfn, 'Error in W2WADDSP')
  seen = .false.; fn = ''
  open (1, FILE=deffn, STATUS='old')
  do
     read(1, *, END=110) iunit, fname, status, form
     if (iunit >= 1 .and. iunit <= 30) then
        seen(iunit) = .true.; fn(iunit) = fname
     end if
     ! only the OUTPUT units are opened as Fortran units; inputs are
     ! read by path through hybrid_files_m
     if (iunit == unit_mmn .or. iunit == unit_amn) &
          open(iunit, FILE=fname, STATUS=status, FORM=form)
  end do
110 continue
  close(1)

!!! case.mmn
  if (seen(unit_mmn1)) then
     call read_mmn(trim(fn(unit_mmn1)), m1, neigh1, gsh1, nb1, nk1, nn1)
     call read_mmn(trim(fn(unit_mmn2)), m2, neigh2, gsh2, nb2, nk2, nn2)
     if (nb1 /= nb2 .or. nk1 /= nk2 .or. nn1 /= nn2) &
          error stop 'w2waddsp: spin mmn files have different dimensions'
     if (any(neigh1 /= neigh2) .or. any(gsh1 /= gsh2)) &
          error stop 'w2waddsp: spin mmn files have different k-neighbour lists'
     m1 = m1 + m2
     deallocate(m2)

     call first_line(trim(fn(unit_mmn1)), title)
     write(unit_mmn, '(A)') trim(title)
     write(unit_mmn, '(3I12)') nb1, nk1, nn1
     ! W2W_MMNTEXT=0 -> header-only text, data lives in the .bin sidecar
     ! alone (readers treat missing/older text + present bin as fresh; the
     ! combined text is ~2.4x the bin and can exceed the free disk at scale)
     call get_environment_variable('W2W_MMNTEXT', envtxt)
     if (trim(envtxt) /= '0') then
        do k = 1, nk1
           do b = 1, nn1
              write(unit_mmn, '(5i8)') k, neigh1(b,k), gsh1(:,b,k)
              do n = 1, nb1
                 do m = 1, nb1
                    write(unit_mmn, '(2F18.12)') m1(m,n,b,k)
                 end do
              end do
           end do
        end do
     end if
     close(unit_mmn)                    ! BEFORE the sidecar: mtime order
     open(newunit=u, file=trim(fn(unit_mmn))//'.bin', form='unformatted', &
          access='stream', status='replace', action='write', iostat=ios)
     if (ios == 0) then
        write(u) 'MMNB', SIDECAR_VER, nb1, nk1, nn1
        write(u) neigh1, gsh1, m1
        close(u)
     end if
     deallocate(m1, neigh1, gsh1, neigh2, gsh2)
  end if

!!! case.amn
  if (seen(unit_amn1)) then
     call read_amn(trim(fn(unit_amn1)), a1, nb1, nk1, nn1)
     call read_amn(trim(fn(unit_amn2)), a2, nb2, nk2, nn2)
     if (nb1 /= nb2 .or. nk1 /= nk2 .or. nn1 /= nn2) &
          error stop 'w2waddsp: spin amn files have different dimensions'
     a1 = a1 + a2
     deallocate(a2)

     call first_line(trim(fn(unit_amn1)), title)
     write(unit_amn, '(A)') trim(title)
     write(unit_amn, '(3I12)') nb1, nk1, nn1
     do k = 1, nk1
        do ipr = 1, nn1
           do m = 1, nb1
              write(unit_amn, '(2i4,1x,i5,1x,2es24.15)') m, ipr, k, a1(m,ipr,k)
           end do
        end do
     end do
     close(unit_amn)
     open(newunit=u, file=trim(fn(unit_amn))//'.bin', form='unformatted', &
          access='stream', status='replace', action='write', iostat=ios)
     if (ios == 0) then
        write(u) 'AMNB', SIDECAR_VER, nb1, nk1, nn1
        write(u) a1
        close(u)
     end if
     deallocate(a1)
  end if

  call errclr(errfn)

contains

  ! title for the combined file = file-1's first line (read '(A)', so blanks
  ! survive -- the old list-directed read truncated at the first space).
  ! Falls back to a fixed tag when only the .bin sidecar exists.
  subroutine first_line(path, line)
    character(*), intent(in)  :: path
    character(*), intent(out) :: line
    integer :: lu, jos
    line = ''
    open(newunit=lu, file=path, status='old', action='read', iostat=jos)
    if (jos == 0) then
       read(lu, '(A)', iostat=jos) line
       close(lu)
    end if
    if (len_trim(line) == 0) line = 'w2waddsp combined'
  end subroutine first_line

END PROGRAM combine_spinfiles


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
