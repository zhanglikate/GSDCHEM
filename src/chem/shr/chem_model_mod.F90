module chem_model_mod

  use chem_types_mod
  use chem_rc_mod
  use chem_domain_mod, only : chem_domain_type
  use chem_state_mod,  only : chem_state_type

  implicit none

  type chem_model_type
    type(chem_domain_type) :: domain
    type(chem_state_type)  :: stateIn, stateOut
  end type chem_model_type

  type(chem_model_type), dimension(:), allocatable, target :: chem_model

  private

  public :: chem_model_init
  public :: chem_model_domain_set
  public :: chem_model_domain_coord_set
  public :: chem_model_get

contains

  subroutine chem_model_init(deCount, rc)

    integer, intent(in),  optional :: deCount
    integer, intent(out), optional :: rc

    !-- local variables
    integer :: localDeCount, localrc

    !-- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    localDeCount = 0
    if (present(deCount)) localDeCount = deCount

    allocate(chem_model(0:localDeCount), stat=localrc)
    if (localrc /= 0) then
      if (present(rc)) rc = CHEM_RC_FAILURE
      return
    end if
    
  end subroutine chem_model_init

  subroutine chem_model_domain_set(minIndex, maxIndex, tile, de, rc)

    integer, dimension(2), intent(in)  :: minIndex, maxIndex
    integer, optional,     intent(in)  :: tile, de
    integer, optional,     intent(out) :: rc

    !-- local variables
    integer :: localDe, localTile

    !-- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    localDe = 0
    if (present(de)) localDe = de

    localTile = 0
    if (present(tile)) localTile = tile

    if (.not.allocated(chem_model)) then
      write(0,'("ERROR: model not allocated. file=",a,", line=",i0)') &
        __FILE__, __LINE__
      if (present(rc)) rc = CHEM_RC_FAILURE
      return
    end if

    if (ubound(chem_model, dim=1) < localDe) then
      write(0,'("ERROR: model does not include localDe=",i0,". file=",a,", line=",i0)') &
        de, __FILE__, __LINE__
      if (present(rc)) rc = CHEM_RC_FAILURE
      return
    end if

    chem_model(localDe) % domain % tile = localTile
    chem_model(localDe) % domain % is   = minIndex(1)
    chem_model(localDe) % domain % js   = minIndex(2)
    chem_model(localDe) % domain % ie   = maxIndex(1)
    chem_model(localDe) % domain % je   = maxIndex(2)

  end subroutine chem_model_domain_set

  subroutine chem_model_domain_coord_set(coordDim, coord, de, rc)

    integer,            intent(in)  :: coordDim
!   real(ESMF_KIND_R8), pointer     :: coord(:,:)
    real(CHEM_KIND_R8), pointer     :: coord(:,:)
    integer, optional,  intent(in)  :: de
    integer, optional,  intent(out) :: rc

    !-- local variables
    integer :: localDe

    !-- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    localDe = 0
    if (present(de)) localDe = de

    if (.not.allocated(chem_model)) then
      write(0,'("ERROR: model not allocated. file=",a,", line=",i0)') &
        __FILE__, __LINE__
      if (present(rc)) rc = CHEM_RC_FAILURE
      return
    end if

    if (ubound(chem_model, dim=1) < localDe) then
      write(0,'("ERROR: model does not include localDe=",i0,'//&
             &'"file=",a,", line=",i0)') de, __FILE__, __LINE__
      if (present(rc)) rc = CHEM_RC_FAILURE
      return
    end if

    select case (coordDim)
      case(1)
        chem_model(localDe) % domain % lon => coord
      case(2)
        chem_model(localDe) % domain % lat => coord
      case default
        write(0,'("ERROR: coordDim can only be 1 or 2."'//&
             &'"file=",a,", line=",i0)')  __FILE__, __LINE__
      if (present(rc)) rc = CHEM_RC_FAILURE
    end select

  end subroutine chem_model_domain_coord_set

  subroutine chem_model_fields_set(numLayers, numTracers, rc)

    integer,           intent(in)  :: numLayers, numTracers
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localDe
    
    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    if (.not.allocated(chem_model)) then
      write(0,'("ERROR: model not allocated. file=",a,", line=",i0)') &
        __FILE__, __LINE__
      if (present(rc)) rc = CHEM_RC_FAILURE
      return
    end if

!   do localDe = 0, size(chem_model)-1
!     call chem_state_init(chem_model(localDe) % StateIn
!   end do

  end subroutine chem_model_fields_set

  subroutine chem_model_get(de, stateIn, stateOut, tile, rc)

    integer,               optional,  intent(in)  :: de
    type(chem_state_type), optional,  pointer     :: stateIn, stateOut
    integer,               optional,  intent(out) :: tile
    integer,               optional,  intent(out) :: rc

    !-- local variables
    integer :: localDe

    !-- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    localDe = 0
    if (present(de)) localDe = de

    if (.not.allocated(chem_model)) then
      write(0,'("ERROR: model not allocated. file=",a,", line=",i0)') &
        __FILE__, __LINE__
      if (present(rc)) rc = CHEM_RC_FAILURE
      return
    end if

    if (ubound(chem_model, dim=1) < localDe) then
      write(0,'("ERROR: model does not include localDe=",i0,'//&
             &'"file=",a,", line=",i0)') de, __FILE__, __LINE__
      if (present(rc)) rc = CHEM_RC_FAILURE
      return
    end if

    if (present(tile))     tile     =  chem_model(localDe) % domain % tile
    if (present(stateIn))  stateIn  => chem_model(localDe) % stateIn
    if (present(stateOut)) stateOut => chem_model(localDe) % stateOut

  end subroutine chem_model_get

end module chem_model_mod
