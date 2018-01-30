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

  public :: chem_model_create
  public :: chem_model_destroy
  public :: chem_model_domain_set
  public :: chem_model_domain_coord_set
  public :: chem_model_get
  public :: chem_model_set

contains

  subroutine chem_model_create(deCount, rc)

    integer, intent(in),  optional :: deCount
    integer, intent(out), optional :: rc

    !-- local variables
    integer :: localDeCount, localrc

    !-- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    localDeCount = 1
    if (present(deCount)) localDeCount = deCount

    allocate(chem_model(0:localDeCount-1), stat=localrc)
    if (localrc /= 0) then
      if (present(rc)) rc = CHEM_RC_FAILURE
      return
    end if
    
  end subroutine chem_model_create

  subroutine chem_model_destroy(rc)

    integer, intent(out), optional :: rc

    !-- local variables
    integer :: localrc

    !-- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    if (allocated(chem_model)) then
      ! -- TODO: deallocate underlying types and arrays
      deallocate(chem_model, stat=localrc)
      if (localrc /= 0) rc = CHEM_RC_FAILURE
    end if

  end subroutine chem_model_destroy

  subroutine chem_model_domain_set(minIndex, maxIndex, minIndexPTile, maxIndexPTile, tile, de, rc)

    integer, dimension(2),           intent(in)  :: minIndex, maxIndex
    integer, dimension(2), optional, intent(in)  :: minIndexPTile, maxIndexPTile
    integer,               optional, intent(in)  :: tile, de
    integer,               optional, intent(out) :: rc

    !-- local variables
    integer                        :: localrc, localTile
    type(chem_model_type), pointer :: model

    !-- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_local_model_get(model, de=de, rc=localrc)
    if (localrc /= CHEM_RC_SUCCESS) then
      if (present(rc)) rc = localrc
      return
    end if

    localTile = 0
    if (present(tile)) localTile = tile

    model % domain % tile = localTile
    model % domain % is   = minIndex(1)
    model % domain % js   = minIndex(2)
    model % domain % ie   = maxIndex(1)
    model % domain % je   = maxIndex(2)
    if (present(minIndexPTile)) then
      model % domain % its = minIndexPTile(1)
      model % domain % jts = minIndexPTile(2)
    end if
    if (present(maxIndexPTile)) then
      model % domain % ite = maxIndexPTile(1)
      model % domain % jte = maxIndexPTile(2)
      if (.not.present(minIndexPTile)) then
        model % domain % its = 1
        model % domain % jts = 1
      end if
    end if

  end subroutine chem_model_domain_set

  subroutine chem_model_domain_coord_set(coordDim, coord, de, rc)

    integer,            intent(in)  :: coordDim
    real(CHEM_KIND_R8), pointer     :: coord(:,:)
    integer, optional,  intent(in)  :: de
    integer, optional,  intent(out) :: rc

    !-- local variables
    integer                        :: localrc
    type(chem_model_type), pointer :: model

    !-- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_local_model_get(model, de=de, rc=localrc)
    if (localrc /= CHEM_RC_SUCCESS) then
      if (present(rc)) rc = localrc
      return
    end if

    select case (coordDim)
      case(1)
        model % domain % lon => coord
      case(2)
        model % domain % lat => coord
      case default
        write(0,'("ERROR: coordDim can only be 1 or 2."'//&
             &'"file=",a,", line=",i0)')  __FILE__, __LINE__
      if (present(rc)) rc = CHEM_RC_FAILURE
    end select

  end subroutine chem_model_domain_coord_set

  subroutine chem_model_set(de, numIntLayers, numModLayers, numTracers, rc)

    integer, optional, intent(in)  :: de
    integer, optional, intent(in)  :: numIntLayers
    integer, optional, intent(in)  :: numModLayers
    integer, optional, intent(in)  :: numTracers
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer                        :: localrc
    type(chem_model_type), pointer :: model
    
    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_local_model_get(model, de=de, rc=localrc)
    if (localrc /= CHEM_RC_SUCCESS) then
      if (present(rc)) rc = localrc
      return
    end if

    if (present(numIntLayers)) model % domain % ni = numIntLayers
    if (present(numModLayers)) model % domain % nl = numModLayers
    if (present(numTracers))   model % domain % nt = numTracers

  end subroutine chem_model_set

  subroutine chem_model_get(de, stateIn, stateOut, tile, rc)

    integer,               optional,  intent(in)  :: de
    type(chem_state_type), optional,  pointer     :: stateIn, stateOut
    integer,               optional,  intent(out) :: tile
    integer,               optional,  intent(out) :: rc

    !-- local variables
    integer                        :: localrc
    type(chem_model_type), pointer :: model
    
    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_local_model_get(model, de=de, rc=localrc)
    if (localrc /= CHEM_RC_SUCCESS) then
      if (present(rc)) rc = localrc
      return
    end if

    if (present(tile))     tile     =  model % domain % tile
    if (present(stateIn))  stateIn  => model % stateIn
    if (present(stateOut)) stateOut => model % stateOut

  end subroutine chem_model_get

  subroutine chem_local_model_get(model, de, rc)

    type(chem_model_type), optional,  pointer     :: model
    integer,               optional,  intent(in)  :: de
    integer,               optional,  intent(out) :: rc

    !-- local variables
    integer :: localDe

    !-- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    nullify(model)

    localDe = 0
    if (present(de)) localDe = de

    if (localDe < 0) then
      write(0,'("ERROR: DE must be >= 0. file=",a,", line=",i0)') &
        __FILE__, __LINE__
      if (present(rc)) rc = CHEM_RC_FAILURE
      return
    end if

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

    model => chem_model(localDe)

  end subroutine chem_local_model_get

end module chem_model_mod
