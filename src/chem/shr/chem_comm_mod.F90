module chem_comm_mod

  use mpi

  implicit none

  integer :: mpi_comm_chem
  integer :: chem_comm_rootpe = 0

  integer, parameter :: RC_COMM_SUCCESS = 0
  integer, parameter :: RC_COMM_FAILURE = -1

  interface chem_comm_bcast
    module procedure chem_comm_bcast_i0
    module procedure chem_comm_bcast_i1
    module procedure chem_comm_bcast_i1a
    module procedure chem_comm_bcast_r0
    module procedure chem_comm_bcast_r1
    module procedure chem_comm_bcast_r1a
    module procedure chem_comm_bcast_s0
    module procedure chem_comm_bcast_s1
    module procedure chem_comm_bcast_s1a
  end interface chem_comm_bcast

  private

  public :: chem_comm_rootpe
  public :: RC_COMM_SUCCESS
  public :: RC_COMM_FAILURE
  public :: chem_comm_init
  public :: chem_comm_isroot
  public :: chem_comm_get
  public :: chem_comm_bcast
  public :: chem_comm_abort
  public :: chem_comm_log

contains

  subroutine chem_comm_init(rc, comm, isolate)
    integer,           intent(out) :: rc
    integer, optional, intent(in)  :: comm
    logical, optional, intent(in)  :: isolate

    ! -- local variables
    integer :: ierr
    logical :: flag

    ! -- begin
    rc = RC_COMM_FAILURE

    call mpi_initialized(flag, ierr)
    if (.not.flag) then
      call mpi_init(ierr)
      if (ierr /= MPI_SUCCESS) return
    end if

    flag = .false.
    if (present(isolate)) flag = isolate
      
    if (present(comm)) then
      if (flag) then
        call mpi_comm_dup(comm, mpi_comm_chem, ierr)
        if (ierr /= MPI_SUCCESS) return
      else
        mpi_comm_chem = comm
      end if
    end if

    ! -- change MPI default error handler to return
    call mpi_errhandler_set(mpi_comm_chem, MPI_ERRORS_RETURN, ierr)
    if (ierr == MPI_SUCCESS) rc = RC_COMM_SUCCESS
   
  end subroutine chem_comm_init

  logical function chem_comm_isroot()
    ! -- local variables
    integer :: ierr, rank
    ! -- begin
    chem_comm_isroot = .false.
    call mpi_comm_rank(mpi_comm_chem, rank, ierr)
    if (ierr /= MPI_SUCCESS) return
    
    chem_comm_isroot = (rank == chem_comm_rootpe)

  end function chem_comm_isroot

  subroutine chem_comm_get(localpe, pecount, comm, rc)
    integer, optional, intent(out) :: localpe
    integer, optional, intent(out) :: pecount
    integer, optional, intent(out) :: comm
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: ierr

    ! -- begin
    if (present(localpe)) then
      localpe = -1
      call mpi_comm_rank(mpi_comm_chem, localpe, ierr) 
    end if
    if (present(pecount)) then
      pecount = -1
      call mpi_comm_size(mpi_comm_chem, pecount, ierr) 
    end if
    if (present(comm)) comm = mpi_comm_chem

    if (present(rc)) rc = RC_COMM_SUCCESS

  end subroutine chem_comm_get

  subroutine chem_comm_log(msg, sync)
    character(len=*),  intent(in) :: msg
    logical, optional, intent(in) :: sync

    write(0,'("chem_comm_log:",a)') trim(msg)
    if (present(sync)) then
      if (sync) call flush(0)
    end if

  end subroutine chem_comm_log

  subroutine chem_comm_abort(errcode, msg)
    integer,          optional, intent(in) :: errcode
    character(len=*), optional, intent(in) :: msg

    ! -- local variables
    integer :: ierr, localerrcode

    ! -- begin
    localerrcode = RC_COMM_FAILURE
    if (present(msg)) write(0,'("chem_comm_abort:",a)') trim(msg)
    if (present(errcode)) localerrcode = errcode
    call mpi_abort(mpi_comm_chem, localerrcode, ierr)

  end subroutine chem_comm_abort

  ! -- broadcast

  subroutine chem_comm_bcast_i0(data, rootpe, rc)
    integer, intent(inout) :: data
    integer, intent(in)    :: rootpe
    integer, intent(out)   :: rc

    ! -- local variables
    integer :: buffer(1)

    ! -- begin
    buffer(1) = data
    call chem_comm_bcast(buffer, 1, rootpe, rc)
    if (rc /= RC_COMM_FAILURE) data = buffer(1)
    
  end subroutine chem_comm_bcast_i0

  subroutine chem_comm_bcast_i1(buffer, count, rootpe, rc)
    integer, intent(inout) :: buffer(:)
    integer, intent(in)    :: count
    integer, intent(in)    :: rootpe
    integer, intent(out)   :: rc

    ! -- local variables
    integer :: ierr

    ! -- begin
    rc = RC_COMM_FAILURE
    if (count > size(buffer)) return
    call mpi_bcast(buffer, count, MPI_INTEGER, rootpe, mpi_comm_chem, ierr)
    if (ierr == MPI_SUCCESS) rc = RC_COMM_SUCCESS
    
  end subroutine chem_comm_bcast_i1

  subroutine chem_comm_bcast_i1a(buffer, rootpe, rc)
    integer, intent(inout) :: buffer(:)
    integer, intent(in)    :: rootpe
    integer, intent(out)   :: rc

    ! -- local variables
    integer :: ierr

    ! -- begin
    rc = RC_COMM_FAILURE
    call mpi_bcast(buffer, size(buffer), MPI_INTEGER, rootpe, mpi_comm_chem, ierr)
    if (ierr == MPI_SUCCESS) rc = RC_COMM_SUCCESS
    
  end subroutine chem_comm_bcast_i1a

  subroutine chem_comm_bcast_r0(data, rootpe, rc)
    real,    intent(inout) :: data
    integer, intent(in)    :: rootpe
    integer, intent(out)   :: rc

    ! -- local variables
    real    :: buffer(1)

    ! -- begin
    buffer(1) = data
    call chem_comm_bcast(buffer, rootpe, rc)
    if (rc /= RC_COMM_FAILURE) data = buffer(1)
    
  end subroutine chem_comm_bcast_r0

  subroutine chem_comm_bcast_r1(buffer, count, rootpe, rc)
    real,    intent(inout) :: buffer(:)
    integer, intent(in)    :: count
    integer, intent(in)    :: rootpe
    integer, intent(out)   :: rc

    ! -- local variables
    integer :: ierr

    ! -- begin
    rc = RC_COMM_FAILURE
    if (count > size(buffer)) return
    call mpi_bcast(buffer, count, MPI_REAL, rootpe, mpi_comm_chem, ierr)
    if (ierr == MPI_SUCCESS) rc = RC_COMM_SUCCESS
    
  end subroutine chem_comm_bcast_r1

  subroutine chem_comm_bcast_r1a(buffer, rootpe, rc)
    real,    intent(inout) :: buffer(:)
    integer, intent(in)    :: rootpe
    integer, intent(out)   :: rc

    ! -- local variables
    integer :: ierr

    ! -- begin
    rc = RC_COMM_FAILURE
    call mpi_bcast(buffer, size(buffer), MPI_REAL, rootpe, mpi_comm_chem, ierr)
    if (ierr == MPI_SUCCESS) rc = RC_COMM_SUCCESS
    
  end subroutine chem_comm_bcast_r1a

  subroutine chem_comm_bcast_s0(data, rootpe, rc)
    character(len=*), intent(inout) :: data
    integer,          intent(in)    :: rootpe
    integer,          intent(out)   :: rc

    ! -- local variables
    integer :: ierr

    ! -- begin
    rc = RC_COMM_FAILURE
    call mpi_bcast(data, len(data), MPI_CHARACTER, rootpe, mpi_comm_chem, ierr)
    if (ierr == MPI_SUCCESS) rc = RC_COMM_SUCCESS
    
  end subroutine chem_comm_bcast_s0

  subroutine chem_comm_bcast_s1(buffer, count, rootpe, rc)
    character(len=*), intent(inout) :: buffer(:)
    integer,          intent(in)    :: count
    integer,          intent(in)    :: rootpe
    integer,          intent(out)   :: rc

    ! -- local variables
    integer :: ierr

    ! -- begin
    rc = RC_COMM_FAILURE
    if (count > size(buffer)) return
    call mpi_bcast(buffer, count*len(buffer(1)), MPI_CHARACTER, rootpe, mpi_comm_chem, ierr)
    if (ierr == MPI_SUCCESS) rc = RC_COMM_SUCCESS
    
  end subroutine chem_comm_bcast_s1

  subroutine chem_comm_bcast_s1a(buffer, rootpe, rc)
    character(len=*), intent(inout) :: buffer(:)
    integer,          intent(in)    :: rootpe
    integer,          intent(out)   :: rc

    ! -- local variables
    integer :: ierr

    ! -- begin
    rc = RC_COMM_FAILURE
    call mpi_bcast(buffer, size(buffer)*len(buffer(1)), MPI_CHARACTER, rootpe, mpi_comm_chem, ierr)
    if (ierr == MPI_SUCCESS) rc = RC_COMM_SUCCESS
    
  end subroutine chem_comm_bcast_s1a

end module chem_comm_mod
