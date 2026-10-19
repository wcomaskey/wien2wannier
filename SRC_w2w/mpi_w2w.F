!!! wien2wannier/SRC_w2w/mpi_w2w.F
!!!
!!!    MPI wrapper module for w2w parallelization.
!!!    All MPI calls are guarded with #ifdef MPI.
!!!    When compiled without -DMPI, provides no-op stubs
!!!    so the code behaves identically to the serial version.
!!!
!!! Usage:
!!!    use mpi_w2w
!!!    call mpi_w2w_init()
!!!    ! ... computation ...
!!!    call mpi_w2w_finalize()

module mpi_w2w
  use const, only: R8, C16
#ifdef MPI
  use mpi
#endif

  implicit none
  public

  integer :: myrank  = 0
  integer :: nprocs  = 1
  logical :: is_root = .true.

contains

  ! Contiguous balanced k-block ownership (D2): rank r owns k-plane block
  ! [r*nk/P+1 .. (r+1)*nk/P].  Replaces the historical round-robin in ALL
  ! FOUR kernels at once (l2Mmn, l2Amn, compute_projectability, planew) --
  ! ownership must stay identical across them or the single overlap reduce
  ! double-counts.  Contiguous blocks shrink the Alm-cache working set and
  ! the planew neighbour halo from ~all-k to block+boundary planes.
  pure integer function k_owner(k, nk)
    integer, intent(in) :: k, nk
    k_owner = ((k-1)*nprocs)/nk
  end function k_owner

  subroutine mpi_w2w_init()
#ifdef MPI
    integer :: ierr
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    is_root = (myrank == 0)
#endif
  end subroutine mpi_w2w_init


  subroutine mpi_w2w_finalize()
#ifdef MPI
    integer :: ierr
    call MPI_Finalize(ierr)
#endif
  end subroutine mpi_w2w_finalize


  subroutine mpi_w2w_barrier()
#ifdef MPI
    integer :: ierr
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
  end subroutine mpi_w2w_barrier


  subroutine mpi_w2w_bcast_int(buf, count)
    integer, intent(inout) :: buf(*)
    integer, intent(in)    :: count
#ifdef MPI
    integer :: ierr
    call MPI_Bcast(buf, count, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
  end subroutine mpi_w2w_bcast_int


  subroutine mpi_w2w_bcast_int_scalar(val)
    integer, intent(inout) :: val
#ifdef MPI
    integer :: ierr
    call MPI_Bcast(val, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
  end subroutine mpi_w2w_bcast_int_scalar


  subroutine mpi_w2w_bcast_real8(buf, count)
    real(R8), intent(inout) :: buf(*)
    integer,  intent(in)    :: count
#ifdef MPI
    integer :: ierr
    call MPI_Bcast(buf, count, MPI_DOUBLE_PRECISION, 0, &
         MPI_COMM_WORLD, ierr)
#endif
  end subroutine mpi_w2w_bcast_real8


  subroutine mpi_w2w_bcast_complex16(buf, count)
    complex(C16), intent(inout) :: buf(*)
    integer,      intent(in)    :: count
#ifdef MPI
    integer :: ierr
    call MPI_Bcast(buf, count, MPI_DOUBLE_COMPLEX, 0, &
         MPI_COMM_WORLD, ierr)
#endif
  end subroutine mpi_w2w_bcast_complex16


  subroutine mpi_w2w_reduce_sum_complex_inplace(buf, count)
    !! In-place reduction: root receives sum into buf, non-root ranks send.
    !!
    !! Two latent defects fixed here:
    !!
    !! (1) ALIASING.  The non-root branch passed `buf` as BOTH sendbuf and
    !!     recvbuf.  MPI forbids aliasing those arguments; it worked by
    !!     implementation accident.  MPI_IN_PLACE is not the fix -- for
    !!     MPI_Reduce it is defined only at the root -- so non-root now passes
    !!     a local dummy.
    !!
    !! (2) SIZE.  `count` was a default INTEGER and the whole buffer went in
    !!     ONE message.  Callers compute Nb*Nb*N_pair, which overflows int32
    !!     above roughly Nb=600 with N_pair=6000; and at Rh3Si4 the message is
    !!     352^2*2880*16 B = 5.32 GiB, 2.66x over the 2 GiB line that several
    !!     MPI builds mishandle.  Count is now 64-bit and the reduction is
    !!     chunked.  Chunking does not re-associate any per-element sum, so the
    !!     result is bit-identical.
    complex(C16),    intent(inout) :: buf(*)
    integer(kind=8), intent(in)    :: count
#ifdef MPI
    integer :: ierr, n
    integer(kind=8) :: off, rem
    integer(kind=8), parameter :: CHUNK = 67108864_8   ! 64 Mi elts = 1 GiB
    complex(C16) :: dummy(1)
    off = 0_8
    rem = count
    do while (rem > 0_8)
       n = int(min(rem, CHUNK))
       if (is_root) then
          call MPI_Reduce(MPI_IN_PLACE, buf(off+1), n, &
               MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       else
          call MPI_Reduce(buf(off+1), dummy, n, &
               MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       end if
       off = off + int(n, 8)
       rem = rem - int(n, 8)
    end do
#endif
  end subroutine mpi_w2w_reduce_sum_complex_inplace


  subroutine mpi_w2w_allreduce_sum_real_inplace(buf, count)
    !! In-place allreduce: all ranks receive sum into buf.
    real(R8), intent(inout) :: buf(*)
    integer,  intent(in)    :: count
#ifdef MPI
    integer :: ierr
    call MPI_Allreduce(MPI_IN_PLACE, buf, count, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
  end subroutine mpi_w2w_allreduce_sum_real_inplace


  subroutine mpi_w2w_allreduce_max_real_inplace(buf, count)
    !! In-place allreduce: all ranks receive elementwise MAX into buf.
    real(R8), intent(inout) :: buf(*)
    integer,  intent(in)    :: count
#ifdef MPI
    integer :: ierr
    call MPI_Allreduce(MPI_IN_PLACE, buf, count, &
         MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif
  end subroutine mpi_w2w_allreduce_max_real_inplace

end module mpi_w2w
