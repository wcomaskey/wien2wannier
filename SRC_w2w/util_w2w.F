!!! wien2wannier/SRC_w2w/util_w2w.F
!!!
!!!    Collection of routines for wien2wannier
!!!
!!! Copyright 2010-2012 Philipp Wissgott
!!!           2013-2016 Elias Assmann


!!/=== Things that are exported here =============
!!
!! const: Constants mathematical and configurational
!!
!!    BUFSZ, DPk, IPk, R4, R8, C16
!!    TAU, ORTHO_TEST, SQ3, Ryd_eV, Ang_Bohr
!!
!!
!! util: Assorted miscellanea
!!
!!    string(), uppercase(), lowercase(), invert3x3()
!!    ptime(), type(ptimer), ptick(), ptock(), ptot(),
!!    line_count(), paropen(), newunit()
!!
!!
!! clio: Helper procedures for command line interface
!!
!!    croak(), carp(), cluck(), type(argstr), fetcharg(), fetchenv()
!!
!!
!! structmod: ‘struct_t’ type and associated procedures
!!
!!    type(struct_t), struct_read()
!!
!!
!! inwfmod: ‘inwf_t’ type and associated procedures
!!
!!    type(inwf_t), inwf_read(), type(proj)
!!
!!
!! kpoints: Reading ‘klist’ files
!!
!!    get_kmesh_klist(), get_kmesh_band(),
!!    count_kmesh_klist(), count_kmesh_band()
!!
!!
!! Wien2k: Procedures nicked from Wien2k
!!
!!    errflg(), errclr(), gtfnam()
!!
!!
!! reallocate: Reallocation of arrays of various types
!!
!!    realloc()
!!
!!\===============================================


!---------------  Constants mathematical and configurational  ---------------
module const
  use iso_fortran_env, only: error_unit

  implicit none
  public

  integer, parameter   :: BUFSZ = 256
  !! Default kinds
  integer, parameter   :: DPk = kind(1.0d0) ! inherited from W90 2.0.0
  integer, parameter   :: IPk = selected_int_kind(8)
  !! Kinds for WIEN-compatibility
  real*4,     private  :: four_real
  real*8,     private  :: eight_real
  complex*16, private  :: sixteen_complex
  integer, parameter   :: R4 = kind(four_real)
  integer, parameter   :: R8 = kind(eight_real)
  integer, parameter   :: C16= kind(sixteen_complex)

  ! The circle constant -- taken from http://tauday.com/tau-digits
  real(DPk), parameter :: TAU        = &
       6.28318530717958647692528676655900576839433879875021_DPk
  real(DPk), parameter :: ORTHO_TEST = 0.0001_DPk ! seems to be WIEN standard
  real(DPk), parameter :: SQ3        = sqrt(3._DPk)
  real(DPk), parameter :: Ryd_eV     = 13.605698_DPk ! value for WIEN compat.
  real(DPk), parameter :: Ang_Bohr   = 1.8897261246_DPk
  real(DPk), parameter :: clight     = 137.0359991_DPk
end module const


!---------------  Assorted miscellanea                        ---------------
module util
  use const, only: DPk, BUFSZ, ERROR_unit
  implicit none

  private
  public :: string, invert3x3, ptime, ptimer, ptick, ptock, ptot, &
       &    uppercase, lowercase, basename, line_count, paropen, newunit

  interface string
     module procedure int2str, ints2str, real2str, reals2str
  end interface string

  interface ptime
!!! Ptime() takes and outputs CPU and wall time measurements.  Typical
!!! usage is
!!!
!!!	ptime([UNIT=…])
!!!
!!! to set output unit and initialize the timers; then
!!!
!!!	… do stuff …
!!!
!!!	call ptime("doing stuff")
!!!
!!!	… do other stuff …
!!!
!!!	call ptime("other stuff")
     module procedure ptime_setunit, ptime_print, ptick
  end interface

  type ptimer
     real(DPk)            :: cpu1 =0, cputot =0
     integer              :: wall1=0, walltot=0, countrate=0
     character(len=BUFSZ) :: &
          fmt = "('Times for ', A, T33, '(sec):', &
          &       F8.3, ' wall;', F9.3, ' CPU')"
     integer              :: unit=ERROR_UNIT
     ! flush after every write? -- initialized to false here because
     ! we do not know if the target unit will support it
     logical              :: flush=.false.
  end type ptimer

  type(ptimer), target, save :: default_timer

contains
  pure character(len=15) function int2str(n)
    integer, intent(in) :: n
    write(int2str,"(I0)") n
  end function int2str

  pure character(len=BUFSZ) function ints2str(na)
    integer, intent(in) :: na(:)

    character(BUFSZ) :: fmt

    ! The generated format MUST carry commas between the literals and the
    ! repeat group: without them gfortran dies AT RUNTIME with "Missing comma
    ! between descriptors" the first time the format executes.  The shipped
    ! version emitted ( '[ 'N(I0,1X)' ]' ) and was runtime-fatal on any call.
    write(fmt, '( "( ''[ '', ", I0, "(I0,1X), '' ]'' )" )') size(na)
    write(ints2str, fmt) na
  end function ints2str

  pure character(len=15) function real2str(r)
    real(DPk), intent(in) :: r
    write(real2str,"(E15.8)") r
  end function real2str

  pure character(len=BUFSZ) function reals2str(ra)
    real(DPk), intent(in) :: ra(:)

    character(BUFSZ) :: fmt

    ! Same missing-comma defect as ints2str above; this one was the MEASURED
    ! crash -- struct_read on a Pnma Rh3Si4 struct file reaches a message
    ! that stringifies reals and died with "Missing comma between descriptors"
    ! (the w2w_hybrid.struct copy does not reach that path, which is why the
    ! unit suite never saw it).
    write(fmt, '( "( ''['', ", I0, "(F0.4,1X), '' ]'' )" )') size(ra)
    write(reals2str, fmt) ra
  end function reals2str

  pure function invert3x3(A) result(Ainv)
    implicit none

    real(DPk), dimension(3,3), intent(in) :: A
    real(DPk), dimension(3,3)             :: Ainv

    real(DPk) :: det

    det = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
         +A(1,3)*A(2,1)*A(3,2) - A(3,1)*A(2,2)*A(1,3) &
         -A(1,1)*A(3,2)*A(2,3) - A(2,1)*A(1,2)*A(3,3)

    Ainv(1,1) = (  A(2,2)*A(3,3) - A(2,3)*A(3,2) ) / det
    Ainv(2,1) = (- A(2,1)*A(3,3) + A(2,3)*A(3,1) ) / det
    Ainv(3,1) = (  A(2,1)*A(3,2) - A(2,2)*A(3,1) ) / det

    Ainv(1,2) = (- A(1,2)*A(3,3) + A(1,3)*A(3,2) ) / det
    Ainv(2,2) = (  A(1,1)*A(3,3) - A(1,3)*A(3,1) ) / det
    Ainv(3,2) = (- A(1,1)*A(3,2) + A(1,2)*A(3,1) ) / det

    Ainv(1,3) = (  A(1,2)*A(2,3) - A(1,3)*A(2,2) ) / det
    Ainv(2,3) = (- A(1,1)*A(2,3) + A(1,3)*A(2,1) ) / det
    Ainv(3,3) = (  A(1,1)*A(2,2) - A(1,2)*A(2,1) ) / det
  end function invert3x3

!!! ptime stuff
  subroutine ptime_setunit(unit, timer, flush)
    integer,              intent(in)              :: unit
    type(ptimer), target, intent(inout), optional :: timer
    logical,              intent(in),    optional :: flush

    integer               :: ios
    type(ptimer), pointer :: t
    t                     => default_timer
    if (present(timer)) t => timer

    t%unit = unit
    if (present(flush)) then
       t%flush = flush
    else                ! by default, we flush if possible
       flush(UNIT=unit, IOSTAT=ios)

       t%flush = ios==0 ! ios>0 signals error, ios<0 means flush not supported
    end if

    call ptick(t)
  end subroutine ptime_setunit

  subroutine ptick(timer)
    type(ptimer), target, intent(inout), optional :: timer

    type(ptimer), pointer :: t
    t                      => default_timer
    if (present(timer)) t  => timer

    call cpu_time    (t% cpu1)
    call system_clock(t%wall1, t%countrate)
  end subroutine ptick

  subroutine ptock(timer)
    type(ptimer), target, intent(inout), optional :: timer

    real(DPk)             :: cpu2
    integer               :: wall2
    type(ptimer), pointer :: t
    t                      => default_timer
    if (present(timer)) t  => timer

    call cpu_time    ( cpu2)
    call system_clock(wall2)

    t% cputot = t% cputot +  cpu2 - t% cpu1
    t%walltot = t%walltot + wall2 - t%wall1

    t% cpu1 =  cpu2
    t%wall1 = wall2
  end subroutine ptock

  subroutine ptot(descr, timer, unit)
    character(len=*),     intent(in)              :: descr
    type(ptimer), target, intent(inout), optional :: timer
    integer,              intent(in),    optional :: unit

    integer               :: lun
    type(ptimer), pointer :: t
    t                      => default_timer
    if (present(timer)) t  => timer
    lun                    =  t%unit
    if (present(unit)) lun =  unit

    write(lun, t%fmt) descr, real(t%walltot, DPk)/t%countrate, &
         &                       (t% cputot)
    if (t%flush) flush(lun)
  end subroutine ptot

  subroutine ptime_print(descr, timer, unit)
    character(len=*),     intent(in)              :: descr
    type(ptimer), target, intent(inout), optional :: timer
    integer,              intent(in),    optional :: unit

    real(DPk)             :: cpu2
    integer               :: wall2, lun
    type(ptimer), pointer :: t

    t                      => default_timer
    if (present(timer)) t  => timer
    lun                    =  t%unit
    if (present(unit)) lun =  unit

    call cpu_time    (cpu2)
    call system_clock(wall2)
    t% cputot = t% cputot +  cpu2 - t% cpu1
    t%walltot = t%walltot + wall2 - t%wall1

    write(lun, t%fmt) descr, real(wall2-t%wall1, DPk)/t%countrate, &
         &                       (cpu2 -t% cpu1)
    if (t%flush) flush(lun)

    t% cpu1 =  cpu2
    t%wall1 = wall2
  end subroutine ptime_print
!!! end ptime stuff

  pure elemental function uppercase(str)
    character(*), intent(in) :: str
    character(len(str))      :: uppercase

    integer :: ic, i

    character(26), parameter :: CAP = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

    uppercase = str
    do i = 1, len_trim(str)
       ic = index(low, str(i:i))
       if (ic > 0) uppercase(i:i) = CAP(ic:ic)
    end do
  end function uppercase

  pure elemental function lowercase(str)
    character(*), intent(in) :: str
    character(len(str))      :: lowercase

    integer :: ic, i

    character(26), parameter :: CAP = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

    lowercase = str
    do i = 1, len_trim(str)
       ic = index(CAP, str(i:i))
       if (ic > 0) lowercase(i:i) = low(ic:ic)
    end do
  end function lowercase

  pure elemental function basename(str)
    character(*), intent(in) :: str
    character(len(str))      :: basename

    basename = str(scan(str, '/', BACK=.true.)+1 : )
  end function basename


  integer function line_count(fid)
    !>Returns the number of lines in a file.
    !>\param fname Name of the file
    !>\return Number of lines in the file.

    implicit none

    !input parameters
    !character(len=*) fname   !filename of the file to count
    integer, intent(in) :: fid

    !local parameters
    character(20) dummy
    integer ioStat
    logical ioEndFlag

    ioEndFlag = .false.
    line_count = 0

    do while (.not. ioEndFlag )
       read(fid,"(A20)", iostat=ioStat ) dummy
       if( iostat .eq. 0) line_count = line_count + 1
       if( iostat < 0 ) ioEndFlag = .TRUE.
    end do
    rewind(fid)
  end function line_count

  subroutine paropen(unit, file, iproc, iloop, iostat, status, form)
    integer,      intent(in)            :: unit
    character(*), intent(in)            :: file
    integer,      intent(in)            :: iloop
    integer,      intent(in)            :: iproc
    integer,      intent(out), optional :: iostat
    character(*), intent(in),  optional :: status
    character(*), intent(in),  optional :: form

    character(len=BUFSZ) :: tmpnam
    character(len=   11) :: f
    character(len=    7) :: s

    f = 'FORMATTED'
    s = 'UNKNOWN'

    if (present(status)) s = status
    if (present(form))   f = form

    if (iproc==0) then
       tmpnam = file
    else
       write(tmpnam, '(A, "_", I0)') trim(file), iloop
    end if

    if (present(iostat)) then
       open(unit, FILE=tmpnam, STATUS=s, FORM=f, IOSTAT=iostat)
    else
       open(unit, FILE=tmpnam, STATUS=s, FORM=f)
    end if
  end subroutine paropen

!!! We would like to use the NEWUNIT keyword to open(), but many
!!! compilers do not support that yet.  Workaround from
!!! <http://fortranwiki.org/fortran/show/newunit>:
  ! This is a simple function to search for an available unit.
  ! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
  ! The UNIT value is returned by the function, and also by the optional
  ! argument. This allows the function to be used directly in an OPEN
  ! statement, and optionally save the result in a local variable.
  ! If no units are available, -1 is returned.
  integer function newunit(unit)
    integer, intent(out), optional :: unit
    ! local
    integer, parameter :: LUN_MIN=1000, LUN_MAX=10000
    logical :: opened
    integer :: lun
    ! begin
    newunit=-1
    do lun=LUN_MIN,LUN_MAX
       inquire(unit=lun,opened=opened)
       if (.not. opened) then
          newunit=lun
          exit
       end if
    end do
    if (present(unit)) unit=newunit
  end function newunit
end module util

!---------------  Helper procedures for command line interface --------------
module clio
  use iso_fortran_env, only: ERROR_UNIT, OUTPUT_UNIT
  use util,            only: string
  use const,           only: BUFSZ, DPk

  implicit none
  private

  public :: croak, carp, cluck, argstr, argstr_set, fetcharg, fetchenv, close_checked

  type argstr
#ifdef HAVE_VARLEN_STR
     character(len=:), allocatable :: s
#else
     character(len=BUFSZ) :: s
#endif
  end type argstr

  interface fetcharg
#ifdef HAVE_VARLEN_STR
     module procedure fetcharg_alloc
#else
     module procedure fetcharg_buf
#endif

     module procedure fetcharg_str, fetcharg_int, fetcharg_real
  end interface fetcharg

  interface fetchenv
#ifdef HAVE_VARLEN_STR
     module procedure fetchenv_alloc
#else
     module procedure fetchenv_buf
#endif

     module procedure fetchenv_str, fetchenv_int, fetchenv_real
  end interface fetchenv

  type(argstr), save :: progname
  logical,      save :: got_progname = .false.
contains
  subroutine get_progname()
    integer s

    call fetcharg(0, progname, status=s)

    if (s /= 0) progname%s='WIEN2WANNIER'

    got_progname = .true.
  end subroutine get_progname

  subroutine croak(message, status)
    character(len=*), intent(in), optional   :: message
    integer,          intent(in), optional   :: status

    integer            :: s
    s=1

    if (present(status)) s=status

    if (present(message)) then
       if (.not. got_progname) call get_progname()

       write(ERROR_UNIT, '(A, ": ", A)') trim(progname%s), message
    end if

    call exit(s)
  end subroutine croak

  subroutine carp(message)
    character(len=*), intent(in) :: message

    if (.not. got_progname) call get_progname()

    write(ERROR_UNIT, '(A, ": ", A)') trim(progname%s), message
  end subroutine carp

  subroutine cluck(message)
    character(len=*), intent(in) :: message

    if (.not. got_progname) call get_progname()

    write(OUTPUT_UNIT, '(A, ": ", A)') trim(progname%s), message
  end subroutine cluck

  subroutine argstr_set(a, c)
    type(argstr),     intent(out) :: a
    character(len=*), intent(in)  :: c

#ifdef HAVE_VARLEN_STR
    if (allocated(a%s)) deallocate(a%s)
    allocate(character(len=len_trim(c)) :: a%s)
#endif

    a%s = c
  end subroutine argstr_set

!!!                               fetcharg ()                              !!!
#ifdef HAVE_VARLEN_STR
  !! This elegant variable-length version of `fetcharg´ does not work
  !! on pre-4.8 gfortran
  subroutine fetcharg_alloc(i, str, message, status)
    integer,          intent(in)            :: i
    type(argstr),     intent(out)           :: str
    character(len=*), intent(in),  optional :: message
    integer,          intent(out), optional :: status

    integer :: s, l

    if (allocated(str%s)) deallocate(str%s)

    call get_command_argument(i, length=l, status=s)
    if (present(status)) then
       status = s
    elseif (s /= 0) then
       if (present(message)) then
          call croak(message)
       else
          call croak("FETCHARG: failed to get command argument #" // &
               &     trim(string(i)) // "'s length: " // trim(string(s)))
       end if
    end if

    allocate(character(len=l) :: str%s)
    !! zero-length arguments seem to need special treatment
    !! (hooray Fortran!)
    if (l==0) return

    call get_command_argument(i, value=str%s, status=s)
    if (present(status)) then
       status = s
    elseif (s /= 0) then
       if (present(message)) then
          call croak(message)
       else
          call croak("FETCHARG: failed to get command argument #" // &
               &     trim(string(i)) // ": " // trim(string(s)))
       end if
    end if
  end subroutine fetcharg_alloc
#else
  subroutine fetcharg_buf(i, str, message, status)
    integer,          intent(in)               :: i
    type(argstr),     intent(out)              :: str
    character(len=*), intent(in),  optional    :: message
    integer,          intent(out), optional    :: status

    integer :: s

    !! here, a zero-length argument is okay
    !! (probably because len(str%s) /= 0)
    call get_command_argument(i, value=str%s, status=s)
    if (present(status)) then
       status = s
    elseif (s /= 0) then
       if (present(message)) then
          call croak(message)
       elseif (s < 0) then
          call croak("FETCHARG_BUF: buffer too small for command argument #"&
               &     // trim(string(i)))
       else
          call croak("FETCHARG_BUF: failed to get command argument #" // &
               &     trim(string(i)) // ": " // trim(string(s)))
       end if
    end if
  end subroutine fetcharg_buf
#endif /* HAVE_VARLEN_STR */

  subroutine fetcharg_str(i, str, message, status)
    integer,          intent(in)               :: i
    character(len=*), intent(out)              :: str
    character(len=*), intent(in),  optional    :: message
    integer,          intent(out), optional    :: status

    integer :: s

    !! here, a zero-length argument is okay
    !! (probably because len(str) /= 0)
    call get_command_argument(i, value=str, status=s)
    if (present(status)) then
       status = s
    elseif (s /= 0) then
       if (present(message)) then
          call croak(message)
       elseif (s < 0) then
          call croak("FETCHARG_STR: buffer too small for command argument #"&
               &     // trim(string(i)))
       else
          call croak("FETCHARG_STR: failed to get command argument #" // &
               &     trim(string(i)) // ": " // trim(string(s)))
       end if
    end if
  end subroutine fetcharg_str

  subroutine fetcharg_int(i, n, message, status)
    integer,          intent(in)               :: i
    integer,          intent(out)              :: n
    character(len=*), intent(in), optional     :: message
    integer,          intent(out), optional    :: status

    integer      :: s
    type(argstr) :: str

    call fetcharg(i, str, STATUS=s)

    if (s /= 0) then
       if (present(status)) then
          status = s
          return
       else
          if (present(message)) then
             call croak(message)
          elseif (s < 0) then
             call croak("FETCHARG_INT: buffer too small &
                  & for command argument #" // trim(string(i)))
          else
             call croak("FETCHARG_INT: failed to get command argument #" // &
                  &     trim(string(i)) // ": " // trim(string(s)))
          end if
       end if
    end if

    read(str%s, *, IOSTAT=s) n

    if (present(status)) then
       status = s
    else if (s /= 0) then
       if (present(message)) then
          call croak(message)
       else
          call croak("FETCHARG_INT: failed to read command argument #" // &
               &  trim(string(i)) // " (`"//trim(str%s)//"') as an integer")
       end if
    end if
  end subroutine fetcharg_int

  subroutine fetcharg_real(i, x, message, status)
    integer,          intent(in)               :: i
    real(DPk),        intent(out)              :: x
    character(len=*), intent(in),  optional    :: message
    integer,          intent(out), optional    :: status

    integer :: s
    type(argstr) :: str

    call fetcharg(i, str, STATUS=s)

    if (s /= 0) then
       if (present(status)) then
          status = s
          return
       else
          if (present(message)) then
             call croak(message)
          elseif (s < 0) then
             call croak("FETCHARG_REAL: buffer too small &
                  & for command argument #" // trim(string(i)))
          else
             call croak("FETCHARG_REAL: failed to get command argument #"// &
                  &     trim(string(i)) // ": " // trim(string(s)))
          end if
       end if
    end if

    read(str%s, *, IOSTAT=s) x

    if (present(status)) then
       status = s
    else if (s /= 0) then
       if (present(message)) then
          call croak(message)
       else
          call croak("FETCHARG_REAL: failed to read command argument #" // &
               &  trim(string(i)) // " (`"//trim(str%s)//"') as a real")
       end if
    end if
  end subroutine fetcharg_real


!!!                               fetchenv ()                              !!!
!!!
!!! If STATUS is supplied, it is returned as set by
!!! get_environment_variable():
!!!
!!!    STATUS= 1: variable not defined
!!!          =-1: VALUE buffer is supplied but too short
!!!          > 1: “other error”
!!!
!!! (This table derived from a comparison of gfortran and ifort docs.)
!!!
!!! Else, for <0 and >1 we raise an error; a non-existent variable
!!! will be returned as blank.
#ifdef HAVE_VARLEN_STR
  !! This elegant variable-length version of `fetchenv´ does not work
  !! on pre-4.8 gfortran
  subroutine fetchenv_alloc(name, str, message, status, trim_name)
    character(len=*), intent(in)            :: name
    type(argstr),     intent(out)           :: str
    character(len=*), intent(in),  optional :: message
    integer,          intent(out), optional :: status
    logical,          intent(in),  optional :: trim_name

    integer :: s, l
    logical :: t

    t = .true.
    if (present(trim_name)) t=trim_name

    if (allocated(str%s)) deallocate(str%s)

    call get_environment_variable(name, length=l, status=s, trim_name=t)
    if (present(status)) then
       status = s
    elseif (s > 1) then
       if (present(message)) then
          call croak(message)
       else
          if (t) then
             call croak("FETCHENV: failed to get environment variable `" // &
                  &     trim(name) // "'s length: " // trim(string(s)))
          else
             call croak("FETCHENV: failed to get environment variable `" // &
                  &     name // "'s length: " // trim(string(s)))
          end if
       end if
    end if

    allocate(character(len=l) :: str%s)
    !! zero-length arguments seem to need special treatment
    !! (hooray Fortran!)
    if (l==0) return

    call get_environment_variable(name, value=str%s, status=s, trim_name=t)
    if (present(status)) then
       status = s
    elseif (s < 0 .or. s > 1) then
       if (present(message)) then
          call croak(message)
       else
          if (t) then
             call croak("FETCHENV: failed to get environment variable `" // &
                  &     trim(name) // "': " // trim(string(s)))
          else
             call croak("FETCHENV: failed to get environment variable `" // &
                  &     name // "': " // trim(string(s)))
          end if
       end if
    end if
  end subroutine fetchenv_alloc
#else
  subroutine fetchenv_buf(name, str, message, status, trim_name)
    character(len=*), intent(in)            :: name
    type(argstr),     intent(out)           :: str
    character(len=*), intent(in),  optional :: message
    integer,          intent(out), optional :: status
    logical,          intent(in),  optional :: trim_name

    integer :: s
    logical :: t

    t = .true.
    if (present(trim_name)) t=trim_name

    !! here, a zero-length argument is okay
    !! (probably because len(str%s) /= 0)
    call get_environment_variable(name, value=str%s, status=s, trim_name=t)
    if (present(status)) then
       status = s
    elseif (s < 0 .or. s > 1) then
       if (present(message)) then
          call croak(message)
       elseif (s < 0) then
          if (t) then
             call croak("FETCHENV_BUF: buffer too small for environment &
                  &variable `" // trim(name) // "'")
          else
             call croak("FETCHENV_BUF: buffer too small for environment &
                  &variable `" // name // "'")
          end if
       else
          if (t) then
             call croak("FETCHENV_BUF: failed to get environment &
                  &variable `" // trim(name) //"': " // trim(string(s)))
          else
             call croak("FETCHENV_BUF: failed to get environment &
                  &variable `" // name // "': " // trim(string(s)))
          end if
       end if
    end if
  end subroutine fetchenv_buf
#endif /* HAVE_VARLEN_STR */

  subroutine fetchenv_str(name, str, message, status, trim_name)
    character(len=*), intent(in)            :: name
    character(len=*), intent(out)           :: str
    character(len=*), intent(in),  optional :: message
    integer,          intent(out), optional :: status
    logical,          intent(in),  optional :: trim_name

    integer :: s
    logical :: t

    t = .true.
    if (present(trim_name)) t=trim_name

    !! here, a zero-length argument is okay
    !! (probably because len(str) /= 0)
    call get_environment_variable(name, value=str, status=s, trim_name=t)
    if (present(status)) then
       status = s
    elseif (s < 0 .or. s > 1) then
       if (present(message)) then
          call croak(message)
       elseif (s < 0) then
          if (t) then
             call croak("FETCHENV_STR: buffer too small for environment &
                  &variable `" // trim(name) // "'")
          else
             call croak("FETCHENV_STR: buffer too small for environment &
                  &variable `" // name // "'")
          end if
       else
          if (t) then
             call croak("FETCHENV_STR: failed to get environment variable `"&
                  // trim(name) // "': " // trim(string(s)))
          else
             call croak("FETCHENV_STR: failed to get environment variable `"&
                  // name // "': " // trim(string(s)))
          end if
       end if
    end if
  end subroutine fetchenv_str

  subroutine fetchenv_int(name, n, message, status, trim_name)
    character(len=*), intent(in)            :: name
    integer,          intent(out)           :: n
    character(len=*), intent(in),  optional :: message
    integer,          intent(out), optional :: status
    logical,          intent(in),  optional :: trim_name

    integer      :: s
    type(argstr) :: str
    logical      :: t

    t = .true.
    if (present(trim_name)) t=trim_name

    call fetchenv(name, str, STATUS=s, TRIM_NAME=t)

    if (s < 0 .or. s > 1) then
       if (present(status)) then
          status = s
          return
       else
          if (present(message)) then
             call croak(message)
          elseif (s < 0) then
             if (t) then
                call croak("FETCHENV_INT: buffer too small &
                     & for environment variable `" // trim(name) // "'")
             else
                call croak("FETCHENV_INT: buffer too small &
                     & for environment variable `" // name // "'")
             end if
          else
             if (t) then
                call croak("FETCHENV_INT: failed to get environment &
                     &variable `" // trim(name) // "': " // trim(string(s)))
             else
                call croak("FETCHENV_INT: failed to get environment &
                     &variable `" // name // "': " // trim(string(s)))
             end if
          end if
       end if
    end if

    read(str%s, *, IOSTAT=s) n

    if (present(status)) then
       status = s
    else if (s /= 0) then
       if (present(message)) then
          call croak(message)
       else if (t) then
          call croak("FETCHENV_INT: failed to read environment variable `" &
               //trim(name)// "' (`"//trim(str%s)//"') as an integer")
       else
          call croak("FETCHENV_INT: failed to read environment variable `" &
               //name// "' (`"//trim(str%s)//"') as an integer")
       end if
    end if
  end subroutine fetchenv_int

  subroutine fetchenv_real(name, x, message, status, trim_name)
    character(len=*), intent(in)            :: name
    real(DPk),        intent(out)           :: x
    character(len=*), intent(in),  optional :: message
    integer,          intent(out), optional :: status
    logical,          intent(in),  optional :: trim_name

    integer      :: s
    type(argstr) :: str
    logical      :: t

    t = .true.
    if (present(trim_name)) t=trim_name

    call fetchenv(name, str, STATUS=s, TRIM_NAME=t)

    if (s < 0 .or. s > 1) then
       if (present(status)) then
          status = s
          return
       else
          if (present(message)) then
             call croak(message)
          elseif (s < 0) then
             if (t) then
                call croak("FETCHENV_INT: buffer too small &
                     & for environment variable `" // trim(name) // "'")
             else
                call croak("FETCHENV_INT: buffer too small &
                     & for environment variable `" // name // "'")
             end if
          else
             if (t) then
                call croak("FETCHENV_INT: failed to get environment &
                     &variable `" // trim(name) // "': " // trim(string(s)))
             else
                call croak("FETCHENV_INT: failed to get environment &
                     &variable `" // trim(name) // "': " // trim(string(s)))
             end if
          end if
       end if
    end if

    read(str%s, *, IOSTAT=s) x

    if (present(status)) then
       status = s
    else if (s /= 0) then
       if (present(message)) then
          call croak(message)
       else if (t) then
          call croak("FETCHENV_REAL: failed to read environment variable `" &
               & // trim(name) // "' (`"//trim(str%s)//"') as a real")
       else
          call croak("FETCHENV_REAL: failed to read environment variable `" &
               & // name // "' (`"//trim(str%s)//"') as a real")
       end if
    end if
  end subroutine fetchenv_real

  !> Closes `unit` and refuses the run (croak) when the file on disk is shorter than what was written to it.
  !> gfortran reports no error for buffered writes that hit a full disk -- not on the write, the flush or the
  !> close, with or without iostat -- and the file simply stops short (measured on this machine, gfortran 16,
  !> macOS).  The size libgfortran believes it wrote, inquired before the close,
  !> is compared with the size on disk after it; on a healthy disk the two are equal.  Unnamed units are just closed.
  subroutine close_checked(unit, what)
    use iso_fortran_env, only: int64
    integer,          intent(in) :: unit
    character(len=*), intent(in) :: what

    character(len=4096) :: fname
    character(len=512)  :: msg
    integer(int64)      :: szw, szd
    integer             :: ios
    logical             :: named

    fname = ''; szw = -1; szd = -1
    inquire(unit=unit, named=named, name=fname, size=szw)
    close(unit, iostat=ios)
    if (named) inquire(file=trim(fname), size=szd)
    if (ios /= 0 .or. (named .and. szw >= 0 .and. szd /= szw)) then
       write(msg, '(a,i0,a,i0,a,i0,a)') ': ', szd, ' bytes on disk but ', szw, &
            ' were written (close iostat ', ios, '). The write was lost, most likely because the disk is full.'
       call croak(trim(what)//' '//trim(fname)//trim(msg)// &
            ' Refusing: the product is incomplete. Free disk space and rerun.')
    end if
  end subroutine close_checked
end module clio


!---------------- ‘struct_t’ type and associated procedures   ---------------
!
! NB: individual wien2wannier programs may still use their own
! struct-related modules
module structmod
  use const, only: DPk

  implicit none

  private
  public :: struct_t, struct_read, radgrid

  interface struct_read
     module procedure struct_read_unit, struct_read_fname
  end interface struct_read

  type intptr
     integer, allocatable :: p(:)
  end type intptr

  type struct_t
     !! Structure to represent a ‘struct’ file

     character(len=80) :: title
     character(len= 4) :: lattic       ! lattice type
     integer           :: Nneq, Nat    ! number of noneq. and total atoms
     character(len= 4) :: mode         ! ‘RELA’ or ‘NREL’
     logical           :: rel          ! rel.?
     real(DPk)         :: a(3)         ! lattice constants
     real(DPk)         :: alpha(3)     ! angles
     logical           :: ortho, prim  ! orthogonal axes?  primitive latt.?
     real(DPk)         :: vol          ! u.c. volume

     ! Primitive and conventional unit cells in direct and reciprocal space
     ! in “column is a vector” convention:
     !
     !    conv_rec(:, i) = b_i   conv_dir(:, i) = a_i
     !    prim_rec(:, i) = b_i'  prim_dir(:, i) = a_i'
     !
     ! We include the factor of τ such that a_i ⋅ b_j = τ δ_ij or
     !       T
     !    rec ⋅ prim = τ I

     real(DPk), dimension(3,3) :: conv_dir, conv_rec, prim_dir, prim_rec
     ! ‘struct’ to “fractional” coordinates
     integer,   dimension(3,3) :: stru2frac

     ! Rotational part of the symmetries in R-space and k-space
     integer,   allocatable :: rsym(:,:,:), ksym(:,:,:)
     ! Translational part of the symmetries in R-space
     real(DPk), allocatable :: rtrans(:,:)
     ! Symmetry transformations that transform the position of one
     ! atom to the corresponding position of an equivalent atom
     integer,   allocatable :: rotij(:,:,:)
     real(DPk), allocatable :: transij(:,:)

     real(DPk),         allocatable :: Z(:)     ! atomic number
     character(len=10), allocatable :: aname(:) ! atom name
     ! positions (3 × Nat); local rotation matrices (3 × 3 × Nneq)
     real(DPk),         allocatable :: pos(:,:), locrot(:,:,:)
     ! multiplicity, ‘isplit’ parameter
     integer,           allocatable :: mult(:), isplit(:)
     ! neq2at(ineq) is the atom numer corresponding to ineq
     integer,           allocatable :: neq2at(:)
     ! at2neq(iat)%p is an array of the ‘ineq’s corresponding to iat
     type(intptr),      allocatable :: at2neq(:)
     ! radial mesh: first point, log. step, muffin-tin radius; #points
     real(DPk),         allocatable :: R0(:), dx(:), RMT(:)
     integer,           allocatable :: Npt(:)
  end type struct_t
contains
  subroutine struct_read_unit(lun, stru)
    use const, only: TAU, SQ3, ORTHO_TEST
    use clio,  only: croak
    use util,  only: invert3x3, lowercase

    integer,        intent(in)  :: lun
    type(struct_t), intent(out) :: stru

    integer   :: iat, ineq, N
    ! abbreviations
    real(DPk) :: pia(3), alpha(3), cosab, cosac, cosbc, sinab, sinbc
    real(DPk) :: rvfac, wurzel

    read(lun, '(A)')          stru%title
    read(lun, '(A, 23X, I3)') stru%lattic, stru%Nneq
    read(lun, '(13X, A)')     stru%mode
    read(lun, '(6F10.6)')     stru%a, stru%alpha

    relativ: select case(lowercase(stru%mode))
    case ('rela')
       stru%rel = .true.
    case ('nrel')
       stru%rel = .false.
    case default
       call croak('MODE should be ‘RELA’ or ‘NREL’, not ‘' &
            //trim(stru%mode)//'’')
    end select relativ

    where (abs(stru%alpha) < 1e-5) stru%alpha = 90

    allocate(stru%Z  (stru%Nneq), stru%mult  (    stru%Nneq), &
         &   stru%RMT(stru%Nneq), stru%isplit(    stru%Nneq), &
         &   stru%Npt(stru%Nneq), stru%at2neq(    stru%Nneq), &
         &   stru%R0 (stru%Nneq), stru%aname (    stru%Nneq), &
         &   stru%dx (stru%Nneq), stru%locrot(3,3,stru%Nneq))

    countat: do iat = 1, stru%Nneq
       read(lun,*) ! first position
       read(lun, '(15X, I2, 17X, I2)') stru%mult(iat), stru%isplit(iat)
       ! rest of positions
       do ineq = 2, stru%mult(iat)
          read(lun,*)
       end do
       read(lun, '(A, 5X, I5, 5X, F10.8, 5X, F10.5, 5X, F5.2)') &
            stru%aname(iat), stru%Npt(iat), stru%R0(iat), stru%RMT(iat), &
            stru%Z(iat)
       read(lun,'(20X, 3F10.7)') stru%locrot(:, 1, iat)
       read(lun,'(20X, 3F10.7)') stru%locrot(:, 2, iat)
       read(lun,'(20X, 3F10.7)') stru%locrot(:, 3, iat)
    end do countat

    ! The logarithmic step is given by the minimal (R0) and maximal
    ! (RMT) radii, and the number of intervals
    stru%dx  = log(stru%RMT / stru%R0) / (stru%Npt-1)
    ! Now we re-set RMT from the inverse of the above — in tests with
    ! typical numbers, the difference is minute (~ 10^-15), but it
    ! seems to be Wien2k convention to do this
    stru%rmt = stru%R0 * exp(stru%dx * (stru%Npt-1))

    stru%Nat = sum(stru%mult)
    allocate( stru%pos  (3,    stru%Nat), stru%transij(3, stru%Nat), &
         &    stru%rotij(3, 3, stru%Nat), stru%neq2at (   stru%Nat))

    N=0
    do iat = 1, stru%Nneq
       stru%neq2at(N+1 : N+stru%mult(iat)) = iat

       allocate(stru%at2neq(iat)%p(stru%mult(iat)))
       stru%at2neq(iat)%p = (/( ineq, ineq = N+1, N+stru%mult(iat) )/)

       N = N+stru%mult(iat)
    end do

    rewind(lun)
    read(lun,*) ! title
    read(lun,*) ! lattic
    read(lun,*) ! mode
    read(lun,*) ! uc
    N=0
    readpos: do iat = 1, stru%Nneq
       N=N+1
       read(lun, '(9X, 3(3X, F10.8))') stru%pos(:, N)
       read(lun,*) ! mult
       do ineq = 2, stru%mult(iat)
          N=N+1
          read(lun, '(9X, 3(3X, F10.8))') stru%pos(:, N)
       end do
       read(lun,*) ! aname
       read(lun,*) ! locrot
       read(lun,*) !
       read(lun,*) !
    end do readpos

    !! Now compute lattice vectors.  This is copied from
    !! SRC_structeditor/module.f (subroutine latgen_struct)
    pia = TAU/stru%a
    alpha = stru%alpha*TAU/360

    cosab = cos(alpha(3)); sinab = sin(alpha(3))
    cosac = cos(alpha(2))
    cosbc = cos(alpha(1)); sinbc = sin(alpha(1))

    stru%conv_rec = 0
    stru%prim_rec = 0

    lattic: select case (stru%lattic(1:1))
    case ('H')
       stru%conv_rec(1,1) = pia(1)*2/SQ3
       stru%conv_rec(1,2) = pia(1)  /SQ3
       stru%conv_rec(2,2) = pia(2)
       stru%conv_rec(3,3) = pia(3)

       stru%prim_rec = stru%conv_rec

       rvfac = 2/SQ3

       stru%ortho = .false.
       stru%prim  = .true.

    case ('S', 'P')             ! what's ‘S’??
       wurzel = sqrt(sinbc**2 - cosac**2 - cosab**2 + 2*cosbc*cosac*cosab)

       stru%conv_rec(1,1) = pia(1) * sinbc/wurzel
       stru%conv_rec(1,2) = pia(2) * (-cosab + cosbc*cosac)/(sinbc*wurzel)
       stru%conv_rec(1,3) = pia(3) * (-cosac + cosbc*cosab)/(sinbc*wurzel)
       stru%conv_rec(2,2) = pia(2) / sinbc
       stru%conv_rec(2,3) =-pia(3) * cosbc/sinbc
       stru%conv_rec(3,3) = pia(3)

       stru%prim_rec = stru%conv_rec

       rvfac = 1/wurzel

       stru%ortho = all(abs(alpha - TAU/4) <= ortho_test)
       stru%prim  = .true.

    case ('F')
       stru%conv_rec(1,1) = pia(1)
       stru%conv_rec(2,2) = pia(2)
       stru%conv_rec(3,3) = pia(3)

       stru%prim_rec(1,:) = pia(1) * (/-1, 1, 1 /)
       stru%prim_rec(2,:) = pia(2) * (/ 1,-1, 1 /)
       stru%prim_rec(3,:) = pia(3) * (/ 1, 1,-1 /)

       rvfac = 4
       stru%ortho = .true.
       stru%prim  = .false.

    case ('B')
       stru%conv_rec(1,1) = pia(1)
       stru%conv_rec(2,2) = pia(2)
       stru%conv_rec(3,3) = pia(3)

       stru%prim_rec(1,:) = pia(1) * (/ 0, 1, 1 /)
       stru%prim_rec(2,:) = pia(2) * (/ 1, 0, 1 /)
       stru%prim_rec(3,:) = pia(3) * (/ 1, 1, 0 /)

       rvfac = 2

       stru%ortho = .true.
       stru%prim  = .false.

    case ('R')
       stru%conv_rec(1, :) = pia(1) * (/  1, 1, -2 /)/SQ3
       stru%conv_rec(2, :) = pia(2) * (/ -1, 1,  0 /)
       stru%conv_rec(3, :) = pia(3)

       stru%prim_rec = stru%conv_rec

       rvfac = 6/SQ3

       stru%ortho = .false.
       stru%prim  = .true.

    case ('C')
       ! “defaults”, to be changed in nonorthogonal XZ-case
       stru%conv_rec(1,1) = pia(1)
       stru%conv_rec(2,2) = pia(2)
       stru%conv_rec(3,3) = pia(3)

       rvfac = 2
       stru%ortho = .true.      ! may be overridden for CXZ
       stru%prim  = .false.

       C: select case (stru%lattic(2:3))
       case ('XY')
          stru%prim_rec(1,:) = pia(1) * (/ 1, 1, 0 /)
          stru%prim_rec(2,:) = pia(2) * (/-1, 1, 0 /)
          stru%prim_rec(3,:) = pia(3) * (/ 0, 0, 1 /)

       case ('YZ')
          stru%prim_rec(1,:) = pia(1) * (/ 1, 0, 0 /)
          stru%prim_rec(2,:) = pia(2) * (/ 0, 1, 1 /)
          stru%prim_rec(3,:) = pia(3) * (/ 0,-1, 1 /)

       case ('XZ')
          ortho: if (abs(alpha(3) - TAU/4) <= ORTHO_TEST) then
             stru%prim_rec(1,:) = pia(1) * (/ 1, 0, 1 /)
             stru%prim_rec(2,:) = pia(2) * (/ 0, 1, 0 /)
             stru%prim_rec(3,:) = pia(3) * (/-1, 0, 1 /)
          else                  ! CXZ monoclinic case
             stru%conv_rec(1,1) = pia(1) / sinab
             stru%conv_rec(1,2) =-pia(2) / sinab * cosab

             stru%prim_rec(1,1) = pia(1) / sinab
             stru%prim_rec(1,2) =-pia(2) / sinab * cosab
             stru%prim_rec(1,3) = pia(1) / sinab
             stru%prim_rec(2,:) = pia(2) * (/ 0, 1, 0 /)
             stru%prim_rec(3,:) = pia(3) * (/-1, 0, 1 /)

             rvfac = 2/sinab
             stru%ortho = .false.
          end if ortho
       end select C

    case default
       rvfac = 0                ! silence warning
       call croak('unknown lattice type `' // trim(stru%lattic) // "'")
    end select lattic

    ! Between the direct and reciprocal lattice: an inversion, a
    ! transposition, and a τ.
    !     _        _       _        _
    !    |          |     |          | -1 T
    !    | a₁ a₂ a₃ | = τ | b₁ b₂ b₃ |
    !    |_        _|     |_        _|

    stru%conv_dir = transpose(invert3x3(stru%conv_rec)) * TAU
    stru%prim_dir = transpose(invert3x3(stru%prim_rec)) * TAU

    ! We don't want no τs in “stru2frac”.
    stru%stru2frac = nint(matmul(invert3x3(stru%conv_rec), &
         &                       stru%prim_rec))

    stru%vol = product(stru%a) / rvfac

    call read_rsym()
    call make_ksym()
    call rotdef()

  contains

    subroutine read_rsym()
!!! This is from SRC_kgen (main.f, addinv.f)
      integer :: Nrsym, i, id

      read(lun, '(I4)') Nrsym
      allocate(stru%rsym  (3,3,Nrsym), &
           &   stru%rtrans(3,  Nrsym))
      do i=1,Nrsym
         read(lun, '(3(3I2,F11.8/), I8)') &
              (stru%rsym(:,id,i), stru%rtrans(id, i), id=1,3)
      end do
    end subroutine read_rsym

    subroutine make_ksym()
!!! This is from SRC_kgen/sdefl.f

      logical                   :: addinv
      integer                   :: i, Nrsym, Nksym

      addinv = .true.
      Nrsym  = size(stru%rsym, 3)
!!! Redefine symmetry operations from struct with
!!! unitary transformation  u(-1) * S * u

      do i=1,Nrsym
         ! if inversion is already present, we do not need to add it
         ! later
         if (all(stru%rsym(:,:,i) == reshape( (/ -1, 0, 0,    &
              &                                   0,-1, 0,    &
              &                                   0, 0,-1 /), &
              &                               (/ 3, 3 /) ))) then
            addinv = .false.
            exit
         end if
      end do

      if (Nrsym>24 .and. addinv) &
           call croak('Nsym > 24 without inversion')

      if (addinv) then
         ! add inversion
         ! FIXME: check if we should actually do this (non-sp non-so)
         Nksym = 2*Nrsym
      else
         Nksym =   Nrsym
      end if
      allocate(stru%ksym(3,3,Nksym))

                  stru%ksym(:,:,       1:Nrsym) =  stru%rsym(:,:, 1:Nrsym)
      if (addinv) stru%ksym(:,:, Nrsym+1:Nksym) = -stru%rsym(:,:, 1:Nrsym)

      if (stru%ortho .or. stru%lattic=='CXZ') then
         ! FIXME: why the check for CXZ??  sdefl.f has it …
         do i=1,Nksym
            stru%ksym(:,:,i) = nint(                             &
                 matmul(matmul(transpose(stru%prim_dir), stru%ksym(:,:,i)), &
                 &      stru%prim_rec) / TAU)
         end do
      end if
    end subroutine make_ksym

    subroutine rotdef
      !     The matrices  rotij(3,3,index)  transform the position of an
      !     atom to its corresponding position of an equivalent atom.
      use const,     only: R8
      use util,      only: string

      implicit none

      real(R8) ::   x(3), x1(3)
      integer  ::   i, m, idx, jdx, jatom, Ncount

      real(R8), parameter :: &
           TOL=1e-7_R8, TOL2=1.5_R8*TOL, ONE=1_R8, HALF=0.5_R8
      integer, parameter :: &
           xy(2) = (/1, 2/), xz(2) = (/1, 3/), yz(2) = (/2, 3/)

      idx=0
      ncount=0

      Nneq: do jatom = 1,stru%Nneq
         jdx = idx+1
         mult: do m = 1,stru%mult(jatom)
            idx = idx+1

            ! Canonical frame for the representative: the identity.  Every
            ! site-group operation maps the representative to itself, so a
            ! first match here is decided by the ORDER of the operations in
            ! the struct file -- identity for SrVO3, but a z-flip plus an
            ! in-plane rotation for MgB2, and arbitrary in general.  That
            ! silently redefines what the inwf (l,m) labels mean, while the
            ! downstream conventions (projmap m-ordering, W90/WannSymm
            ! positional assignment) all assume the rotloc frame.
            ! Measured A/B (MgB2, the case where op 1 is not the identity):
            ! the amn gauge changes by O(1), but mmn is bit-identical and
            ! Omega_I agrees to 8 digits -- a gauge choice, not physics.
            if (idx == jdx) then
               Ncount = Ncount+1
               stru%transij(:, idx) = 0
               stru%rotij(:,:, idx) = 0
               stru%rotij(1,1,idx) = 1
               stru%rotij(2,2,idx) = 1
               stru%rotij(3,3,idx) = 1
               cycle mult
            end if

            sym: do i = 1,size(stru%rsym, 3)
               x = matmul(transpose(stru%rsym(:,:, i)), stru%pos(:, jdx))
               x = x + stru%rtrans(:, i)

               x1 = mod(abs(x - stru%pos(:, idx))+TOL, ONE) - TOL

               if(maxval(abs(x1)) < TOL2) then
                  Ncount = Ncount+1
                  stru%transij(:, idx) = stru%rtrans(:, i)
                  stru%rotij(:,:, idx) = stru%rsym(:,:, i)
                  cycle mult
               end if

               ! check positions for centered lattices
               B: if (stru%lattic == 'B') then
                  x1 = mod(x1 + HALF+TOL, ONE) - TOL

                  if(maxval(abs(x1)) < TOL2) then
                     Ncount = Ncount+1

                     stru%transij(:, idx) = stru%rtrans(:, i)
                     stru%rotij(:,:, idx) = stru%rsym(:,:, i)
                     cycle mult
                  end if
               end if B

               FXY: if (stru%lattic == 'F' .or. stru%lattic == 'CXY') then
                  x1(xy) = mod(x1(xy) + HALF+TOL, ONE) - TOL

                  if(maxval(abs(X1)) < TOL2) then
                     Ncount = Ncount+1
                     stru%transij(:, idx) = stru%rtrans(:, i)
                     stru%rotij(:,:, idx) = stru%rsym(:,:, i)
                     cycle mult
                  end if

                  x1(xy) = mod(x1(xy) + HALF, ONE)
               end if FXY

               FXZ: if (stru%lattic == 'F' .or. stru%lattic == 'CXZ') then
                  x1(xz) = mod(x1(xz) + HALF+TOL, ONE) - TOL

                  if(maxval(abs(x1)) < TOL2) then
                     Ncount = Ncount+1
                     stru%transij(:, idx)=stru%rtrans(:, i)
                     stru%rotij(:,:, idx)=stru%rsym(:,:, i)
                     cycle mult
                  end if
                  x1(xz)=mod(x1(xz)+HALF, ONE)
               end if FXZ

               FYZ: if(stru%lattic == 'F' .or. stru%lattic == 'CYZ') then
                  x1(yz)=mod(x1(yz) + HALF+TOL, ONE) - TOL

                  if(maxval(abs(x1)) < TOL2) then
                     Ncount = Ncount+1
                     stru%transij(:,idx) = stru%rtrans(:, i)
                     stru%rotij(:,:,idx) = stru%rsym(:,:, i)
                     cycle mult
                  end if
               end if FYZ
            end do sym

            call croak("error in ROTDEF: no symmetry operation found for &
                 &atoms "//trim(string(jatom))//" and "//trim(string(idx))  &
                 &//", positions "//trim(string(stru%pos(:, jdx)))//" and " &
                 &//trim(string(stru%pos(:, idx))))
         end do mult
      end do Nneq

      if (Ncount /= idx) call croak('error in ROTDEF: &
           &ROTIJ not defined for all atoms of basis.  NCOUNT = ' &
           &// trim(string(Ncount)))
    end subroutine rotdef
end subroutine struct_read_unit

subroutine struct_read_fname(fname, stru)
  use util, only: newunit

  character(*),   intent(in)  :: fname
  type(struct_t), intent(out) :: stru

  integer :: lun

  open(unit=newunit(lun), file=fname, status='old')
  call struct_read_unit(lun, stru)
  close(lun)
end subroutine struct_read_fname

subroutine struct_read_argstr(arg, stru)
  use clio, only: argstr

  type(argstr),   intent(in)  :: arg
  type(struct_t), intent(out) :: stru

  call struct_read_fname(arg%s, stru)
end subroutine struct_read_argstr

pure elemental real(DPk) function radgrid(stru, iAt, iR)
  type(struct_t), intent(in) :: stru
  integer,        intent(in) :: iAt, iR

  radgrid = stru%R0(iAt) * exp( (iR-1) * stru%dx(iAt) )
end function radgrid
end module structmod


!---------------- ‘inwf_t’ type and associated procedures   -----------------
module inwfmod
  use const, only: DPk

  implicit none

  private
  public :: inwf_t, inwf_read, auto_config_t

  interface inwf_read
     module procedure inwf_read_unit, inwf_read_fname, inwf_read_argstr
  end interface inwf_read

  type proj
     integer                   :: NY
     integer,      allocatable :: iat(:), l(:), m(:)
     complex(DPk), allocatable :: coeff(:)
     ! radial tag per Y-term: 1 = u_l (legacy projector, default);
     ! 3..5 = orthogonalized-shell LO projector using that rf's radial
     ! (semicore shells invisible to the u-based projector, e.g. Ga-3d)
     integer,      allocatable :: rf(:)
  end type proj

  type auto_config_t
     real(DPk) :: threshold  = 0.10_DPk  ! min relative weight (w_i/w_max)
     integer   :: max_proj   = 0         ! max projections (0=auto)
     logical   :: write_proj = .false.   ! write .autoproj diagnostic
     logical   :: p2_weight  = .true.    ! erfc energy weight in ranking (P2)
     real(DPk) :: support_ratio = 0.5_DPk ! P3 gate (0 disables); full windows
                                          ! with a dominant channel need < 0.5
     ! ---- selection payload, set by generate_projections (in-process
     ! consumers read these instead of parsing projmap.dat; replaces the
     ! old negative-max_proj mask_cols smuggle) ----
     integer :: mask_cols = 0            ! cliff mode: Otsu mask boundary in
                                         ! amn columns (leading block)
     integer :: nsel      = 0            ! selected channels (post-feasibility)
     integer,   allocatable :: sel_atom(:), sel_l(:), sel_rf(:)  ! per channel
     logical,   allocatable :: sel_mask(:)     ! physical (Otsu upper) channel
     real(DPk), allocatable :: sel_weight(:), sel_support(:)
  end type auto_config_t

  type inwf_t
     !! Structure to represent a 'inwf' file
     character(len=4)        :: mode
     logical                 :: Mmn  = .false.
     logical                 :: Amn  = .false.
     logical                 :: auto_mode = .false.
     logical                 :: hybrid    = .false.  ! MODE=HYBR in-memory Stage-2
     integer                 :: bmin, bmax, Nproj, LJmax
     integer                 :: num_wann = 0         ! HYBR: target WF count
                                                     ! (0 = capacity-driven auto)
     ! HYBR: Zd relaxation.  1 = off; >1 fixed over-relaxation; <0 SCHEDULED
     ! (start at |mix|, taper geometrically toward 1).
     !
     ! Default -2.5 + kick, chosen by WORST-CASE REGRET across four cells,
     ! not by best-case on one.  Measured Omega_I, best-of-three per cell:
     !   cell     fixed 1.9      -1.7            -2.5+kick
     !   MgB2     9.21800613     9.21800615      9.21800614   (one minimum)
     !   Fe-up    6.51159711     6.97607570 !!   6.51382534
     !   Fe-dn    6.57988366     7.17933349 !!   6.57597181
     !   GaAs    21.86150562    21.83306107     21.82791941
     ! Worst-case excess: -2.5+kick 0.034%, fixed 1.9 0.154%, -1.7 9.2%.
     ! -1.7 converges on Fe in SIXTY sweeps to a 7-9% worse minimum with
     ! converged=T and no warning -- the fastest run and the worst answer.
     ! The ranking REVERSES between materials, so no fixed alpha is optimal
     ! everywhere; this is the minimax choice, not the best choice.
     !
     ! Historical note -- the earlier reasoning, kept because it is still  Measured on the GaAs
     ! one-shot: -1.5 and -1.7 both converge to Omega_I = 21.83306107 (the
     ! true minimum, reached independently), while -1.9 "converges" 6.6x
     ! sooner at 21.86150570 -- 0.13% HIGHER.  That is a false convergence
     ! into a flatter, worse basin, and `converged=T` does not distinguish
     ! the two.  -1.7 is the fastest factor that still reaches the minimum,
     ! and it matches the Rh3Si4 H120 relaxation probes.  At mix=1 the
     ! disentangler needs far more sweeps and still does not get there; the
     ! one-shot silently ran that way until 2026-07-29.
     real(DPk)               :: hyb_mix  = -2.5_DPk
     character(len=16)       :: hyb_regime = 'auto'  ! HYBR: mask regime
     ! HYBR keyword tail (optional `key = value' lines after line 5;
     ! defaults = the validated recipe)
     character(len=16)       :: hyb_froz_mode = 'auto'   ! auto | none | <eV>
     real(DPk)               :: hyb_froz_val  = 8.0_DPk  ! ceiling for numeric mode
     real(DPk)               :: hyb_pool_factor = 1.5_DPk
     integer                 :: hyb_margin = 0           ! 0 = max(4, 4% of need)
     character(len=16)       :: hyb_gauge  = 'auto8'     ! gauge dispatch mode
     character(len=16)       :: hyb_seedk  = 'opt16'     ! seed policy
     character(len=16)       :: hyb_seedloss = 'smooth'  ! tournament loss
     logical                 :: hyb_cap_strict = .true.  ! capacity violations
                                                         ! croak (strict) or
                                                         ! warn (capacity=warn)
     ! ---- Gate A (W8): what to do when a band-window edge or a pool edge
     ! falls INSIDE a degenerate multiplet.  A cut through a multiplet puts
     ! only part of it in the model, and no unitary rotation of the model's
     ! states can restore the degeneracy (MgB2 12x12x8: the window 5-22 cuts
     ! the doublet at bands 22/23 at 8 of 1152 k; the window 5-23 is closed).
     !   auto   -- snap AUTO-sentinel edges to the nearest closed edge, then
     !             warn on a surviving stuck edge.
     !   warn   -- diagnose only; never snap, never croak (pre-W8 behaviour,
     !             and the keyword-level rollback)
     !   strict -- croak on any surviving stuck edge
     !   snap   -- snap ANY stuck edge within the walk cap, sentinel or
     !             explicit, then behave as `strict' on what survives
     ! Overridden by the environment variable W2W_EDGES (which is the only
     ! way to set it for a non-HYBR deck, since the keyword tail is HYBR's).
     character(len=16)       :: hyb_edges = 'auto'
     ! Disentangler controls.  These were hardcoded at the call site until
     ! 2026-07-29, which left the one-shot running a configuration we never
     ! validated: half the sweep budget of the production recipe, with both
     ! the relaxation and the adaptive convergence cadence off.  Defaults
     ! below ARE the production recipe.
     ! Sweep cap 20000, not the recipe's 4000.  A cap costs nothing when
     ! convergence is quick (the loop exits on the tolerance), but 4000 is
     ! demonstrably too few: GaAs needs 10042 sweeps at relax=-1.7 to reach
     ! its minimum, and feup/fedn hit 4000 unconverged in every recorded
     ! baseline; on the LCAO frontend SnTe needs 15889 sweeps.  A run that
     ! still stops at the cap is graded (flat or still descending) in the
     ! log and in <case>-hyb.dis.  Lower it deliberately on very large cells
     ! if you would rather fail fast than wait.
     integer                 :: hyb_sweeps   = 20000    ! disentangle sweep cap
     ! Convergence floor.  The test is `dlast < tol .OR. dlast < rel_tol*|Om|',
     ! so the LOOSER branch decides: at Omega ~ 22 the rel_tol branch fires at
     ! 2.2e-9 and `tol' never gets a say.  Tighten BOTH to move the floor.
     ! Re-kick the relaxation schedule when the taper has spent itself,
     ! instead of finishing the run effectively unrelaxed.  Experimental
     ! until measured across materials -- default off.
     ! The upstream text .amn/.mmn are INTERMEDIATES under HYBR: nothing in
     ! process or in the hand-off chain reads them (Stage-2 works from memory,
     ! Wannier90 gets <case>-hyb.*), and the binary sidecars carry the same
     ! numbers for every in-tree reader.  Off by default; turn on for debugging
     ! or to feed a stock-Wannier90 full-window run.
     ! Omega_I solver.  'cg' (default since 2026-10) = direct Riemannian CG on
     ! the Grassmann manifold, stopping on the projected gradient; it
     ! parallelizes exactly over k with no colouring.  'gs' = the Gauss-Seidel
     ! sweep (the default before).  'auto' picks on ndis and says why
     ! (LCAO-measured per-step crossover ndis ~ 250).
     character(len=8)        :: hyb_solver   = 'cg'
     logical                 :: hyb_text_io  = .false.
     logical                 :: hyb_kick     = .true.
     real(DPk)               :: hyb_taper    = 0.5_DPk   ! schedule down-step
     real(DPk)               :: hyb_taper_late = 0.0_DPk ! 0 = same as taper
     ! Late-phase restart target.  MEASURED optimum per cell (sweeps at the
     ! minimum, opening held at -2.5 so the basin is fixed and this isolates
     ! the local rate):
     !   GaAs  1.99 (1596)   MgB2 1.95 (996)
     !   Fe-up 1.95 (1034)   Fe-dn 1.95 (1020)
     ! The curve is the classical SOR shape -- monotone down to an optimum,
     ! then a sharp rise as the iteration turns oscillatory -- but the optimum
     ! is MATERIAL-DEPENDENT, so any fixed value is wrong somewhere.  Worst-case
     ! extra sweeps for a single hard-coded choice: 1.85 +392%, 1.95 +168%,
     ! 1.99 +41%.  1.99 is therefore the minimax pick, NOT a good pick; a
     ! runtime probe (golden section on alpha at the first kick, ~10% overhead)
     ! would drive this to zero and is the right eventual fix.
     ! At exactly 2.0 GaAs loses quality as well as rate (Omega_I 21.8280846 vs
     ! 21.8279193) -- Kahan's bound showing up in the answer, not just the cost.
     real(DPk)               :: hyb_relax_late = 1.99_DPk ! 0 = same as relax
     real(DPk)               :: hyb_tol      = 1.0e-9_DPk
     real(DPk)               :: hyb_rel_tol  = 1.0e-10_DPk
     logical                 :: hyb_adaptive = .true.   ! defer the exact
                                                        ! Omega_I check until
                                                        ! the free eigen-sum
                                                        ! plateaus
     integer                 :: hyb_anderson = 0        ! AA history depth
                                                        ! (0 = off)
     ! Red-black k-parallel sweep.  OFF by default and deliberately so: it is
     ! a different Gauss-Seidel ordering that lands on a different (equally
     ! valid) fixed point, and w2w/w2wc are built WITHOUT -fopenmp, so
     ! enabling it here would change the answer while buying no speed.
     logical                 :: hyb_colored  = .false.
     type(proj),          allocatable :: projections(:)
     type(auto_config_t) :: auto_config
  end type inwf_t

contains

  subroutine inwf_read_unit(lun, inwf)
    use clio,  only: croak
    use util,  only: uppercase

    integer,      intent(in)  :: lun
    type(inwf_t), intent(out) :: inwf

    integer             :: w, y, NY, confios
    real(DPk)           :: re, im
    character(len=256)  :: confline

    read(lun, *) inwf%mode

    select case (uppercase(inwf%mode))
    case('AMN')
       inwf%Amn = .true.
    case('MMN')
       inwf%Mmn = .true.
    case('BOTH')
       inwf%Amn = .true.
       inwf%Mmn = .true.
    case('AUTO')
       inwf%Amn = .true.
       inwf%Mmn = .true.
       inwf%auto_mode = .true.
    case('AAMN')
       ! AUTO channel selection, amn ONLY: re-select projections against an
       ! EXISTING mmn (same band window) without the overlap recompute.
       ! Rh3Si4: the mmn is 12 GB / hours, the amn minutes -- this is what
       ! makes widening the candidate set affordable.
       inwf%Amn = .true.
       inwf%Mmn = .false.
       inwf%auto_mode = .true.
    case('HYBR')
       ! AUTO + in-memory hybrid Stage-2 (masks -> disentangle -> SCDM ->
       ! isolated Nb=Nw hand-off) run inside w2w after l2Amn.
       inwf%Amn = .true.
       inwf%Mmn = .true.
       inwf%auto_mode = .true.
       inwf%hybrid = .true.
    case default
       call croak('unknown MODE in inwf, '//inwf%mode)
    end select

    read(lun, *) inwf%bmin, inwf%bmax
    if (inwf%auto_mode) then
       ! AUTO sentinels resolved in main after the energy-file scan:
       !   bmin =  0 -> gap-walk auto-bmin (trim core/semicore complexes)
       !   bmax =  0 -> all bands present at every k
       !   bmax = -1 -> highest band with bottom <= E_F + 40 eV
       if (inwf%bmin < 0 .or. inwf%bmax < -1) &
            call croak('invalid band index in inwf')
       if (inwf%bmin > 0 .and. inwf%bmax > 0 .and. inwf%bmax < inwf%bmin) &
            call croak('bmin > bmax in inwf')
    else
       if (inwf%bmax < inwf%bmin) &
            call croak('bmin > bmax in inwf')
    end if

    if (inwf%auto_mode) then
       ! AUTO mode: read LJmax and auto-config parameters.  The config
       ! line is `threshold max_proj write_proj [p2_weight]'; the 4th
       ! field is optional (default T) for back-compatibility.
       read(lun, *) inwf%LJmax
       read(lun, '(A)') confline
       read(confline, *, IOSTAT=confios)          &
            &       inwf%auto_config%threshold,   &
            &       inwf%auto_config%max_proj,    &
            &       inwf%auto_config%write_proj,  &
            &       inwf%auto_config%p2_weight,   &
            &       inwf%auto_config%support_ratio
       if (confios /= 0) then
          inwf%auto_config%support_ratio = 0.5_DPk
          read(confline, *, IOSTAT=confios)          &
               &       inwf%auto_config%threshold,   &
               &       inwf%auto_config%max_proj,    &
               &       inwf%auto_config%write_proj,  &
               &       inwf%auto_config%p2_weight
          if (confios /= 0) then
             inwf%auto_config%p2_weight = .true.
             read(confline, *) inwf%auto_config%threshold, &
                  &            inwf%auto_config%max_proj,  &
                  &            inwf%auto_config%write_proj
          end if
       end if
       inwf%Nproj = 0  ! will be set by generate_projections
       if (inwf%hybrid) then
          ! HYBR extra line: `num_wann [mix] [regime]'.  mix = Zd relaxation:
          ! 1 = off, >1 fixed over-relaxation, NEGATIVE = scheduled (start at
          ! |mix|, taper toward 1); default -2.5 (see hyb_mix: -1.9 reaches
          ! the convergence test sooner but at a worse Omega_I).
          ! regime = auto | per-state | band-level, default auto.
          ! num_wann = 0 means capacity-driven auto (requires an explicit
          ! froz_emax keyword).  Prefer the `relax =' keyword to the
          ! positional field -- it is clearer and does not depend on the
          ! optional-field parse succeeding.
          read(lun, '(A)') confline
          read(confline, *, IOSTAT=confios) inwf%num_wann, inwf%hyb_mix, &
               inwf%hyb_regime
          if (confios /= 0) then
             inwf%hyb_regime = 'auto'
             read(confline, *, IOSTAT=confios) inwf%num_wann, inwf%hyb_mix
             if (confios /= 0) then
                inwf%hyb_mix = -2.5_DPk
                read(confline, *) inwf%num_wann
             end if
          end if
          if (inwf%num_wann < 0) call croak('HYBR: num_wann must be >= 0 &
               &(0 = capacity-driven auto)')
          ! optional keyword tail: `key = value' lines until EOF/blank.
          ! Unknown keys are a hard error (typo protection).
          kwloop: do
             read(lun, '(A)', IOSTAT=confios) confline
             if (confios /= 0) exit kwloop
             if (len_trim(confline) == 0) exit kwloop
             block
               integer :: ieq, iosk
               character(len=64)  :: key
               character(len=192) :: val
               ieq = index(confline, '=')
               if (ieq == 0) call croak('HYBR keyword line without "=": '// &
                    trim(confline))
               key = adjustl(confline(1:ieq-1))
               val = adjustl(confline(ieq+1:))
               select case (trim(key))
               case ('froz_emax')
                  if (trim(val) == 'auto' .or. trim(val) == 'none') then
                     inwf%hyb_froz_mode = trim(val)
                  else
                     inwf%hyb_froz_mode = trim(val)   ! numeric mode: the value
                     read(val, *, IOSTAT=iosk) inwf%hyb_froz_val
                     if (iosk /= 0) call croak('HYBR froz_emax: bad value '// &
                          trim(val))
                  end if
               case ('pool_factor')
                  read(val, *, IOSTAT=iosk) inwf%hyb_pool_factor
                  if (iosk /= 0) call croak('HYBR pool_factor: bad value')
               case ('margin')
                  read(val, *, IOSTAT=iosk) inwf%hyb_margin
                  if (iosk /= 0) call croak('HYBR margin: bad value')
               case ('gauge')
                  inwf%hyb_gauge = trim(val)
               case ('seedk')
                  inwf%hyb_seedk = trim(val)
               case ('seedloss')
                  inwf%hyb_seedloss = trim(val)
               case ('capacity')
                  select case (trim(val))
                  case ('strict'); inwf%hyb_cap_strict = .true.
                  case ('warn');   inwf%hyb_cap_strict = .false.
                  case default
                     call croak('HYBR capacity: strict|warn, got '//trim(val))
                  end select
               case ('edges')
                  select case (trim(val))
                  case ('auto', 'warn', 'strict', 'snap')
                     inwf%hyb_edges = trim(val)
                  case default
                     call croak('HYBR edges: auto|warn|strict|snap, got '// &
                          trim(val))
                  end select
               case ('sweeps')
                  read(val, *, IOSTAT=iosk) inwf%hyb_sweeps
                  if (iosk /= 0 .or. inwf%hyb_sweeps < 1) &
                       call croak('HYBR sweeps: want a positive integer, got '&
                       &//trim(val))
               case ('relax')
                  read(val, *, IOSTAT=iosk) inwf%hyb_mix
                  if (iosk /= 0) call croak('HYBR relax: bad value '//trim(val))
                  if (inwf%hyb_mix == 0.0_DPk) call croak('HYBR relax: 0 is &
                       &not a relaxation factor (1 = off, <0 = scheduled)')
               case ('adaptive')
                  inwf%hyb_adaptive = yesno(val, 'adaptive')
               case ('colored')
                  inwf%hyb_colored = yesno(val, 'colored')
               case ('relax_late')
                  if (trim(val) == 'auto') then
                     inwf%hyb_relax_late = -1.0_DPk   ! probe at the first kick
                  else
                     read(val, *, IOSTAT=iosk) inwf%hyb_relax_late
                     if (iosk /= 0) call croak('HYBR relax_late: want a real &
                          &or "auto", got '//trim(val))
                  end if
               case ('taper_late')
                  read(val, *, IOSTAT=iosk) inwf%hyb_taper_late
                  if (iosk /= 0 .or. inwf%hyb_taper_late <= 0.0_DPk &
                       .or. inwf%hyb_taper_late >= 1.0_DPk) call croak('HYBR &
                       &taper_late: want 0 < taper_late < 1, got '//trim(val))
               case ('taper')
                  read(val, *, IOSTAT=iosk) inwf%hyb_taper
                  if (iosk /= 0 .or. inwf%hyb_taper <= 0.0_DPk &
                       .or. inwf%hyb_taper >= 1.0_DPk) call croak('HYBR &
                       &taper: want 0 < taper < 1, got '//trim(val))
               case ('solver')
                  select case (trim(val))
                  case ('gs', 'cg', 'auto'); inwf%hyb_solver = trim(val)
                  case default
                     call croak('HYBR solver: gs|cg|auto, got '//trim(val))
                  end select
               case ('text_io')
                  inwf%hyb_text_io = yesno(val, 'text_io')
               case ('kick')
                  inwf%hyb_kick = yesno(val, 'kick')
               case ('tol')
                  read(val, *, IOSTAT=iosk) inwf%hyb_tol
                  if (iosk /= 0 .or. inwf%hyb_tol <= 0.0_DPk) &
                       call croak('HYBR tol: want a positive real, got '&
                       &//trim(val))
               case ('rel_tol')
                  read(val, *, IOSTAT=iosk) inwf%hyb_rel_tol
                  if (iosk /= 0 .or. inwf%hyb_rel_tol <= 0.0_DPk) &
                       call croak('HYBR rel_tol: want a positive real, got '&
                       &//trim(val))
               case ('anderson')
                  read(val, *, IOSTAT=iosk) inwf%hyb_anderson
                  if (iosk /= 0 .or. inwf%hyb_anderson < 0) &
                       call croak('HYBR anderson: want a history depth >= 0 &
                       &(0 = off), got '//trim(val))
               case default
                  call croak('HYBR: unknown keyword "'//trim(key)//'"')
               end select
             end block
          end do kwloop
       end if
    else
       read(lun,*) inwf%LJmax, inwf%Nproj

       if (inwf%Amn) then
          allocate(inwf%projections(inwf%Nproj))
          do w = 1, inwf%Nproj
             read(lun, *) NY
             inwf%projections(w)%NY = NY
             allocate( inwf%projections(w)%iat  (NY), &
                  &    inwf%projections(w)%l    (NY), &
                  &    inwf%projections(w)%m    (NY), &
                  &    inwf%projections(w)%coeff(NY), &
                  &    inwf%projections(w)%rf   (NY))
             inwf%projections(w)%rf = 1

             ylm: do y = 1,NY
                ! optional 6th field: radial-function index (1 = u_l,
                ! 3..5 = orthogonalized LO-shell projector, as in projmap)
                read(lun, '(A)') confline
                read(confline, *, IOSTAT=confios)      &
                     &       inwf%projections(w)%iat(y), &
                     &       inwf%projections(w)%l  (y), &
                     &       inwf%projections(w)%m  (y), &
                     re, im, inwf%projections(w)%rf(y)
                if (confios /= 0) then
                   inwf%projections(w)%rf(y) = 1
                   read(confline, *) inwf%projections(w)%iat(y), &
                        &            inwf%projections(w)%l  (y), &
                        &            inwf%projections(w)%m  (y), &
                        re, im
                end if
                inwf%projections(w)%coeff(y) = cmplx(re, im, DPk)
             end do ylm
          end do
       end if
    end if
  end subroutine inwf_read_unit

  subroutine inwf_read_fname(fname, inwf)
    use util, only: newunit

    character(*), intent(in)  :: fname
    type(inwf_t), intent(out) :: inwf

    integer :: lun

    open(unit=newunit(lun), file=fname, status='old')
    call inwf_read_unit(lun, inwf)
    close(lun)
  end subroutine inwf_read_fname

  subroutine inwf_read_argstr(arg, inwf)
    use clio, only: argstr

    type(argstr), intent(in)  :: arg
    type(inwf_t), intent(out) :: inwf

    call inwf_read_fname(arg%s, inwf)
  end subroutine inwf_read_argstr

  !> Parse a HYBR keyword-tail logical.  Accepts the spellings a WIEN2k user
  !> is likely to reach for; anything else is a hard error rather than a
  !> silent false, on the same typo-protection principle as the unknown-key
  !> croak (a mistyped `adaptive = ture' must not quietly halve throughput).
  logical function yesno(val, key)
    use clio, only: croak

    character(len=*), intent(in) :: val, key
    character(len=16) :: v
    integer :: i
    v = adjustl(val)
    do i = 1, len_trim(v)                       ! fold to lower case
       if (v(i:i) >= 'A' .and. v(i:i) <= 'Z') &
            v(i:i) = achar(iachar(v(i:i)) + 32)
    end do
    select case (trim(v))
    case ('t', 'true', '.true.', 'yes', 'y', 'on', '1')
       yesno = .true.
    case ('f', 'false', '.false.', 'no', 'n', 'off', '0')
       yesno = .false.
    case default
       call croak('HYBR '//trim(key)//': want true|false, got '//trim(val))
       yesno = .false.                          ! unreachable; quiets -Wall
    end select
  end function yesno
end module inwfmod


!---------------- Reading ‘klist’ files            --------------------------
module kpoints
  use const,     only: BUFSZ, DPk, TAU
  use structmod, only: struct_t
  use util,      only: string, newunit
  use clio,      only: croak

  implicit none
  private
  public :: get_kmesh_klist, get_kmesh_band, &
       &    count_kmesh_klist, count_kmesh_band

  interface get_kmesh_klist
     module procedure get_kmesh_klist_unit, get_kmesh_klist_fname
  end interface get_kmesh_klist

  interface get_kmesh_band
     module procedure get_kmesh_band_unit, get_kmesh_band_fname
  end interface get_kmesh_band

  interface count_kmesh_klist
     module procedure count_kmesh_klist_unit, count_kmesh_klist_fname
  end interface count_kmesh_klist

  interface count_kmesh_band
     module procedure count_kmesh_band_unit, count_kmesh_band_fname
  end interface count_kmesh_band

  character(len=*), parameter :: &
       fmt_klist_old   = "(A10,4I 5)",         &
       fmt_klist_new   = "(A10,4I10)",         &
       fmt_klist_old_1 = "(10X,4I 5,35X,3I3)", &
       fmt_klist_new_1 = "(10X,4I10,35X,3I3)"
  integer,      parameter :: fmt_klist_len = len(fmt_klist_new)
contains
!!! Read standard ‘klist’
  subroutine get_kmesh_klist_unit(unit, kpoints, stru, mp_grid)
    integer,        intent(in)               :: unit
    real(DPk),      intent(out), allocatable :: kpoints(:,:)
    type(struct_t), intent(in)               :: stru
    integer,        intent(out), optional    :: mp_grid(3)

    integer                  :: nk, ik
    character(fmt_klist_len) :: fmt
    character(10)            :: kname
    integer                  :: kfrac(3), div
    real(DPk)                :: RR(3,3)

    if (stru%ortho) then
       RR(1,:) = stru%prim_dir(:,1) / stru%a
       RR(2,:) = stru%prim_dir(:,2) / stru%a
       RR(3,:) = stru%prim_dir(:,3) / stru%a
    end if

    if (present(mp_grid)) then
       call klist_detect_fmt(unit, fmt, mp_grid)
    else
       call klist_detect_fmt(unit, fmt)
    end if

    nk = count_kmesh_klist(unit)
    allocate(kpoints(nk, 3))

    ik = 0
    readk: do
       read(unit, fmt) kname, kfrac, div
       if (kname == "END") exit

       ik = ik+1

       kpoints(ik,:) = real(kfrac, DPk)/real(div, DPk)

       if (stru%ortho) kpoints(ik,:) = matmul(RR, kpoints(ik,:))
    end do readk

    if (ik /= nk) &
         call croak('get_kmesh_klist: numbers of k-points do not match' // &
         &          trim(string(ik)) // ' /= ' // trim(string(nk)))
  end subroutine get_kmesh_klist_unit

!!! Read ‘klist_band’, taking into account only labeled points
  subroutine get_kmesh_band_unit(unit, kpoints, stru, knames)
!!! Reading the band-structure klist gets its own subroutine because
!!! the required work (k-point labels, #points per segment) as well as
!!! the structure of the file (no running numbers, no mp_grid, ...)
!!! are sufficiently different.

    integer,                intent(in)  :: unit
    real(DPk), allocatable, intent(out) :: kpoints(:,:)
    type(struct_t),         intent(in)  :: stru
    character, allocatable, intent(out) :: knames(:)

    integer                  :: nk, ik
    character(fmt_klist_len) :: fmt
    character(10)            :: kname
    integer                  :: kfrac(3), div
    real(DPk)                :: RR(3,3)

    if (stru%ortho) then
       RR(1,:) = stru%prim_dir(:,1) / stru%a
       RR(2,:) = stru%prim_dir(:,2) / stru%a
       RR(3,:) = stru%prim_dir(:,3) / stru%a
    end if

    call klist_detect_fmt(unit, fmt)
    nk = count_kmesh_band(unit)

    allocate(kpoints(nk, 3), knames(nk))

    ik  = 0
    readk: do
       read(unit, fmt) kname, kfrac, div
       if (kname == "END") exit
       if (kname == "")    cycle

       ! now we have a k-vector we should save
       ik = ik+1
       knames(ik) = kname(1:1)

       kpoints(ik,:) = real(kfrac, DPk)/real(div, DPk)

       if (stru%ortho) &
            kpoints(ik,:) = matmul(RR, kpoints(ik,:))
    end do readk

    if (ik /= nk) &
         call croak('get_kmesh_band: numbers of k-points do not match' // &
         &          trim(string(ik)) // ' /= ' // trim(string(nk)))
  end subroutine get_kmesh_band_unit

!!! Count k-points in standard ‘klist’
  function count_kmesh_klist_unit(unit) result(nk)
    integer, intent(in)  :: unit
    integer              :: nk

    character(len=BUFSZ) :: buf

    nk = 0
    countk: do
       read(unit, '(A)') buf
       if (buf == "END") exit
       nk = nk+1
    end do countk
    rewind(unit)
  end function count_kmesh_klist_unit

!!! Count labelled k-points in ‘klist_band’
  function count_kmesh_band_unit(unit) result(nk)
    integer, intent(in)  :: unit
    integer              :: nk

    character(len=6)     :: buf

    nk = 0
    countk: do
       read(unit, '(A)') buf
       if (buf == "END") exit
       ! we only count k-points with a label, assuming that these
       ! delimit "paths" that will be passed to wannier90
       if (buf /= "") nk = nk+1
    end do countk
    rewind(unit)
  end function count_kmesh_band_unit

!!! Wrappers for calling the above with a filename instead of an open unit
  subroutine get_kmesh_klist_fname(fname, kpoints, stru, mp_grid)
    character(*),   intent(in)               :: fname
    real(DPk),      intent(out), allocatable :: kpoints(:,:)
    type(struct_t), intent(in)               :: stru
    integer,        intent(out), optional    :: mp_grid(3)

    integer :: lun

    open(unit=newunit(lun), file=fname, status='old')
    call get_kmesh_klist_unit(lun, kpoints, stru, mp_grid)
    close(lun)
  end subroutine get_kmesh_klist_fname

  subroutine get_kmesh_band_fname(fname, kpoints, stru, knames)
    character(*),           intent(in)  :: fname
    real(DPk), allocatable, intent(out) :: kpoints(:,:)
    type(struct_t),         intent(in)  :: stru
    character, allocatable, intent(out) :: knames(:)

    integer :: lun

    open(unit=newunit(lun), file=fname, status='old')
    call get_kmesh_band_unit(lun, kpoints, stru, knames)
    close(lun)
  end subroutine get_kmesh_band_fname

  function count_kmesh_klist_fname(fname) result(nk)
    character(*), intent(in)  :: fname
    integer                   :: nk

    integer :: lun

    open(unit=newunit(lun), file=fname, status='old')
    nk = count_kmesh_klist_unit(lun)
    close(lun)
  end function count_kmesh_klist_fname

  function count_kmesh_band_fname(fname) result(nk)
    character(*), intent(in)  :: fname
    integer                   :: nk

    integer :: lun

    open(unit=newunit(lun), file=fname, status='old')
    nk = count_kmesh_band_unit(lun)
    close(lun)
  end function count_kmesh_band_fname

!!! Try to figure out which format our klist file follows, optionally
!!! reading ‘mp_grid’ in the process.
  subroutine klist_detect_fmt(unit, fmt, mp_grid)
    integer,                  intent(in)  :: unit
    character(fmt_klist_len), intent(out) :: fmt
    integer,        optional, intent(out) :: mp_grid(3)

    character(BUFSZ) :: iom
    integer          :: idum(4), ios

    ! First, we try the new format:
    if (present(mp_grid)) then
       read(unit, fmt_klist_new_1, iostat=ios, iomsg=iom) idum, mp_grid
    else
       read(unit, fmt_klist_new_1, iostat=ios, iomsg=iom) idum
    end if
    if (ios == 0) then
       ! Success — new format it is
       fmt = fmt_klist_new
    else if (ios < 0) then
       call croak('get_kmesh_klist(): EOF encountered while reading first&
            & line of ‘klist’ file: '//trim(iom))
    else
       ! Error in read, try old format
       if (present(mp_grid)) then
          read(unit, fmt_klist_old_1, iostat=ios, iomsg=iom) idum, mp_grid
       else
          read(unit, fmt_klist_old_1, iostat=ios, iomsg=iom) idum
       end if
       if (ios == 0) then
          fmt = fmt_klist_old
       else
          call croak('get_kmesh_klist(): something went wrong while reading&
               & first line of ‘klist’ file: '//trim(iom))
       end if
    end if

    rewind(unit)
  end subroutine klist_detect_fmt
end module kpoints


!--------------- Procedures nicked from Wien2k       ------------------------
module Wien2k
  private
  public :: errflg, errclr, gtfnam, vector_kind_probe
contains
  subroutine ERRFLG(FNAME,MSG)
    character(*)      FNAME, MSG
!
!     ..................................................................
!
! 1.     PROGRAM UNIT 'ERRFLG'
!           Notify that an error has occured.
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           Because there is no standard (or even semi-standard) way to
!           generate exit codes in FORTRAN 77, this routine writes a
!           non-empty file to the current subdirectory as an indication
!           that some serious error has occured. Other programs can then
!           check the contents of this file to determine whether an
!           error has occured. The errorflag-file is left opened when
!           returning from this routine to enable writing other
!           errormessages to it.
!
! 3.     USAGE
!           CALL ERRFLG('lapw2.error','Error in OUTWIN')
!
!        ARGUMENT-DESCRIPTION
!           FNAME  - CHARACTER*(*) string                        (input)
!                    The name of the file acting as error-flag.
!           MSG    - CHARACTER*(*) string                        (input)
!                    The (error) message which should be written to the
!                    errorflag-file.
!
!        USED SUBROUTINES (DIRECTLY CALLED)
!           none
!
!        INDIRECTLY CALLED SUBROUTINES
!           none
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           ERRCLR - clears the contents of a file
!
!        INPUT/OUTPUT (READ/WRITE)
!           A message given as argument MSG is written to a the file
!           given as argument FNAME. File FNAME is created if not
!           existing and otherwise overwritten.
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           none
!
! 4.     REMARKS
!           The best way to use this routine is to call ERRFLG at the
!           start of a program writing some message to the errorflag-
!           file and ERRCLR before a successful exit of the program. By
!           checking the contents of the errorflag-file it is possible
!           to determine if the program was successfully completed.
!
!           This method has the advantage of working even if some
!           runtime-error occurs which is not taken care of in the
!           program.
!
! 5.     METHOD
!           - open errorflag-file
!           - write some message to this errorflag-file
!           - exit leaving the file opened
!
! 6.     DATE
!           24. August 1993                                 Version 1.02
!
!        INSTITUT FUER TECHNISCHE ELEKTROCHEMIE            --  TU VIENNA
!     ..................................................................
!
    open (99,FILE=FNAME,ERR=900)
    write (99,9000) MSG
    close (99)
    open (99,FILE=FNAME,ERR=900)
!
    return
!
!        Errors
!
900 stop 'ERRFLG - couldn''t open errorflag-file.'
!
9000 format (A)
!
!        End of 'ERRFLG'
!
  end subroutine ERRFLG


  subroutine ERRCLR(FNAME)
    character(*)      FNAME

    close (99)
    open (99,FILE=FNAME,STATUS='REPLACE')
  end subroutine ERRCLR

  subroutine GTFNAM(DEFFN,ERRFN,IPROC)
    implicit none

    character(*)      DEFFN, ERRFN
    integer            IPROC
!
!     ..................................................................
!
! 1.     PROGRAM UNIT 'GTFNAM'
!           Provide names for definition- and error-files for LAPW.
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           Read the commandline-argument
!           specifying the name of the definition file ('lapw1.def')
!           and generate the name of the error file by replacing the
!           extension of the definition filename with '.error'. If no
!           extension can be found '.error' is appended.
!
!           For the parallel version in the second commandline parameter
!           the number of parallel processes is specified
!
! 3.     USAGE
!           CHARACTER*80 DEFFN, ERRFN
!           CALL GTFNAM(DEFFN,ERRFN)
!
!        ARGUMENT-DESCRIPTION
!           DEFFN  - CHARACTER*(*) string                       (output)
!                    on exit contains the filename of the lapw2-
!                    definition file 'lapw2.def'.
!           ERRFN  - CHARACTER*(*) string                       (output)
!                    on exit contains the filename of the file where
!                    error messages are stored (derived from the file-
!                    name of the definition-file).
!           IPROC  - number of parallel processes (if specified)
!
!
!        USED SUBROUTINES (DIRECTLY CALLED)
!           none
!
!        INDIRECTLY CALLED SUBROUTINES
!           none
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           none
!
!        INPUT/OUTPUT (READ/WRITE)
!           none
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           - Subroutine GETARG for the extraction of command-line
!             arguments is used (the index for referencing command-
!             line argument starts with 0 on some machines and with
!             1 on other machines).
!           - A compiler directive to enable the use of the extension
!             subroutine GETARG is used (HP-Version)
!
! 4.     REMARKS
!           It is assumed that filename-extensions are separated by
!           character '.'.
!
! 5.     METHOD
!           - get commandline-argument (taking into account that
!             comandline-arguments are referenced starting with index
!             0 on some machines and starting with index 1 on others).
!           - 'lapw1.def' := commandline-argument
!           - look for the last occurence of character '.' in the
!             commandline-argument
!           - if found replace all characters after that '.' with
!             'error' giving the error filename
!           - otherwise append '.error' giving the error filename
!
! 6.     DATE
!           25. August 1993                                 Version 1.01
!
!        INSTITUT FUER TECHNISCHE ELEKTROCHEMIE            --  TU VIENNA
!     ..................................................................
!
!        Local Parameters
      !
    character(5)       ERREXT
    parameter          (ERREXT = 'error')
!
!        Local Scalars
!
    integer            I, iarg
!
!        extract the command-line argument
!
    IPROC=0
    iarg=command_argument_count()
    if(iarg.eq.1) then
       CALL get_command_argument(iarg,DEFFN)
    else if(iarg.eq.2) then
       CALL get_command_argument(2,DEFFN)
       READ(DEFFN,*)IPROC
       CALL get_command_argument(1,DEFFN)
    else
       STOP 'GTFNAM - One or two commandline arguments have to be given.'
    endif

!        generate a name for the error-message file
    do I = len(DEFFN), 1, -1
       if (DEFFN(I:I) .eq. '.') then
          if (len(ERRFN) .lt. (I+len(ERREXT))) goto 910
          ERRFN(1:I) = DEFFN(1:I)
          ERRFN(I+1:len(ERRFN)) = ERREXT
          goto 30
       end if
    end do

!        the name of the definition file contains no '.', it is assumed
!        that this name contains no extension - append the extension
!        '.error' to get a name for the error file.
    do I = len(DEFFN), 1, -1
       if (DEFFN(I:I) .ne. ' ') then
          if (len(ERRFN) .lt. (I+1+len(ERREXT))) goto 910
          ERRFN(1:I) = DEFFN(1:I)
          ERRFN(I+1:LEN(ERRFN)) = '.' // ERREXT
          goto 30
       end if
    end do

!        filename contains only spaces
    stop 'GTFNAM - string ERRFN contains just spaces.'
30  continue

    return

!        Errors
910 stop 'GTFNAM - string ERRFN too short to hold filename.'
  end subroutine GTFNAM

  !> Is this WIEN2k vector file COMPLEX or REAL?  2 = complex, 1 = real,
  !> 0 = could not tell (the caller must NOT act on 0).
  !>
  !> Real vs complex is a COMPILE-TIME choice (-D_COMPLEX_ selecting w2w /
  !> w2wc) made by the caller via `x w2w -c`.  No input file records it, and
  !> WIEN2k's "case.in1c exists => complex" convention does not settle it
  !> either (a real MgB2 case can carry BOTH in1 and in1c).  The two
  !> mistakes are badly asymmetric:
  !>
  !>   w2wc on a REAL file   -> dies loudly ("I/O past end of record"): a
  !>                            complex read wants 2N reals from a record of N.
  !>   w2w on a COMPLEX file -> SILENT.  Reading N reals from a record holding
  !>                            N complex numbers is a legal SHORT read, so the
  !>                            run completes and writes a hand-off.  Measured
  !>                            on GaAs: Omega_I = -127.61 vs the correct 21.83.
  !>
  !> The test is exact, not heuristic: a band-coefficient record holds N values
  !> of the file's own kind, hence N reals in a real file and 2N in a complex
  !> one.  Sequential unformatted allows reading FEWER items than a record
  !> holds but errors on MORE, so a 2N-real read succeeds iff the file is
  !> complex.  Probing uses a SEPARATE unit on the same path, so the main read
  !> stream is untouched and no BACKSPACE (undefined after a failed read) is
  !> needed.
  function vector_kind_probe(fname, nneq) result(kindcode)
    use const, only: DPk

    character(len=*), intent(in) :: fname
    integer,          intent(in) :: nneq
    integer :: kindcode

    integer :: u, i, n, ne, ios, num
    real(DPk) :: xk, yk, zk, e1
    character(len=10) :: bname
    real(DPk), allocatable :: probe(:)
    integer,   allocatable :: g(:,:)

    kindcode = 0                       ! unknown until proven otherwise
    open(newunit=u, file=fname, status='old', form='unformatted', &
         action='read', iostat=ios)
    if (ios /= 0) return

    do i = 1, nneq                     ! two per-atom header records each
       read(u, iostat=ios) ; if (ios /= 0) goto 90
       read(u, iostat=ios) ; if (ios /= 0) goto 90
    end do

    read(u, iostat=ios) xk, yk, zk, bname, n, ne
    if (ios /= 0 .or. n < 1) goto 90

    allocate(g(3,n))
    read(u, iostat=ios) (g(1,i), g(2,i), g(3,i), i = 1, n)
    if (ios /= 0) goto 90

    read(u, iostat=ios) num, e1        ! first (band index, eigenvalue)
    if (ios /= 0) goto 90

    allocate(probe(2*n))               ! the discriminating read
    read(u, iostat=ios) (probe(i), i = 1, 2*n)
    if (ios == 0) then
       kindcode = 2                    ! record held 2N reals -> complex
    else
       kindcode = 1                    ! over-read -> only N reals -> real
    end if

90  continue
    close(u, iostat=ios)
  end function vector_kind_probe
end module Wien2k


!----------------  Reallocation of various arrays           -----------------
module reallocate
  use const, only: DPk, BUFSZ
  use clio,  only: croak

  implicit none

  private
  public :: realloc

  interface realloc
     ! Fortran makes this pretty cumbersome.  We need to provide one
     ! specific subroutine for each type and rank.
     module procedure realloc_r_d1, realloc_r_d2, realloc_r_d3
     module procedure realloc_c_d1, realloc_c_d2, realloc_c_d3
     module procedure realloc_i_d1, realloc_i_d2, realloc_i_d3
     module procedure realloc_a_d1, realloc_a_d2, realloc_a_d3
#ifdef HAVE_PTR_ALLOC_GENERICS
     ! In Fortran 2008, generic procedures can distinguish between
     ! ‘pointer’ and ‘allocatable’ variables, but older compilers do
     ! not support this yet, so we only support ‘allocatable’ there.
     module procedure realloc_r_p1, realloc_r_p2, realloc_r_p3
     module procedure realloc_c_p1, realloc_c_p2, realloc_c_p3
     module procedure realloc_i_p1, realloc_i_p2, realloc_i_p3
     module procedure realloc_a_p1, realloc_a_p2, realloc_a_p3
#endif
  end interface realloc

contains

!!! Procedures for ‘allocatable’ variables, which use move_alloc()
#define EXPDIM1(NAME) \
    (1:NAME(1))
#define EXPDIM2(NAME) \
    (1:NAME(1), 1:NAME(2))
#define EXPDIM3(NAME) \
    (1:NAME(1), 1:NAME(2), 1:NAME(3))

#define DIMS1 (:)
#define DIMS2 (:,:)
#define DIMS3 (:,:,:)

#ifdef __GFORTRAN__
#define EXPDIM(NDIM, NAME) EXPDIM/**/NDIM (NAME)
#define DIMS(NDIM) DIMS/**/NDIM
#else
#define EXPDIM(NDIM, NAME) EXPDIM##NDIM (NAME)
#define DIMS(NDIM) DIMS##NDIM
#endif

    !     nur 1 mal kopieren reicht
    !     auch fuer mehrdimensionale Felder schaut die Zuweisung gleich aus
#define MK_REALLOC(NDIM, T, TFATTR, TFINT) (tf, nd)               ;\
    implicit none                                                 ;\
                                                                   \
    integer,          intent(in)    :: nd(NDIM)                   ;\
    T,        TFATTR, intent(TFINT) :: tf DIMS(NDIM)              ;\
                                                                   \
    T,        TFATTR :: hilfsfeld DIMS(NDIM)                      ;\
    integer          :: mins(NDIM)                                ;\
                                                                   \
    allocate( hilfsfeld EXPDIM(NDIM, nd) )                        ;\
    mins = min(nd, shape(tf))                                     ;\
    hilfsfeld EXPDIM(NDIM, mins) = tf EXPDIM(NDIM, mins)

#define MK_REALLOC_MV(NDIM, T)               \
    MK_REALLOC(NDIM, T, allocatable, inout) ;\
    call move_alloc(hilfsfeld, tf)

  subroutine realloc_r_d1 MK_REALLOC_MV(1, real(DPk))
  end subroutine
  subroutine realloc_r_d2 MK_REALLOC_MV(2, real(DPk))
  end subroutine
  subroutine realloc_r_d3 MK_REALLOC_MV(3, real(DPk))
  end subroutine

  subroutine realloc_c_d1 MK_REALLOC_MV(1, complex(DPk))
  end subroutine
  subroutine realloc_c_d2 MK_REALLOC_MV(2, complex(DPk))
  end subroutine
  subroutine realloc_c_d3 MK_REALLOC_MV(3, complex(DPk))
  end subroutine

  subroutine realloc_i_d1 MK_REALLOC_MV(1, integer)
  end subroutine
  subroutine realloc_i_d2 MK_REALLOC_MV(2, integer)
  end subroutine
  subroutine realloc_i_d3 MK_REALLOC_MV(3, integer)
  end subroutine

  subroutine realloc_a_d1 MK_REALLOC_MV(1, character(BUFSZ))
  end subroutine
  subroutine realloc_a_d2 MK_REALLOC_MV(2, character(BUFSZ))
  end subroutine
  subroutine realloc_a_d3 MK_REALLOC_MV(3, character(BUFSZ))
  end subroutine

#ifdef HAVE_PTR_ALLOC_GENERICS
!!! Procedures for ‘pointer’ variables, which use pointer assignment

  ! pointer INTENT refers to association status
#define MK_REALLOC_PTR(NDIM, T)          \
    MK_REALLOC(NDIM, T, pointer, out)   ;\
    deallocate(tf)                      ;\
    tf=>hilfsfeld ! der Zeiger wird nur auf das neue Feld umgebogen,
                  ! nicht neu alloziert

  subroutine realloc_r_p1 MK_REALLOC_PTR(1, real(DPk))
  end subroutine realloc_r_p1
  subroutine realloc_r_p2 MK_REALLOC_PTR(2, real(DPk))
  end subroutine
  subroutine realloc_r_p3 MK_REALLOC_PTR(3, real(DPk))
  end subroutine

  subroutine realloc_c_p1 MK_REALLOC_PTR(1, complex(DPk))
  end subroutine
  subroutine realloc_c_p2 MK_REALLOC_PTR(2, complex(DPk))
  end subroutine
  subroutine realloc_c_p3 MK_REALLOC_PTR(3, complex(DPk))
  end subroutine

  subroutine realloc_i_p1 MK_REALLOC_PTR(1, integer)
  end subroutine
  subroutine realloc_i_p2 MK_REALLOC_PTR(2, integer)
  end subroutine
  subroutine realloc_i_p3 MK_REALLOC_PTR(3, integer)
  end subroutine

  subroutine realloc_a_p1 MK_REALLOC_PTR(1, character(BUFSZ))
  end subroutine
  subroutine realloc_a_p2 MK_REALLOC_PTR(2, character(BUFSZ))
  end subroutine
  subroutine realloc_a_p3 MK_REALLOC_PTR(3, character(BUFSZ))
  end subroutine

#endif /* HAVE_PTR_ALLOC_GENERICS */
end module reallocate


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
