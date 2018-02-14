module CHM

  !-----------------------------------------------------------------------------
  ! CHM Component.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, inheritModel => SetServices

  use chem_methods
  
  implicit none

#if 0
  integer, parameter :: importFieldCount = 9
  character(len=*), dimension(importFieldCount), parameter :: &
    importFieldNames = (/ &
      "air_pressure                   ",  &
      "air_pressure_in_model_layers   ",  &
      "geopotential                   ",  &
      "geopotential_in_model_layers   ",  &
      "air_temperature                ",  &
      "x_wind                         ",  &
      "y_wind                         ",  &
      "omega                          ",  &
      "mass_fraction_of_tracers_in_air"   &
    /)
#else
  integer, parameter :: importFieldCount = 25
  character(len=*), dimension(importFieldCount), parameter :: &
    importFieldNames = (/ &
      "air_pressure                             ", &
      "air_pressure_in_model_layers             ", &
      "geopotential                             ", &
      "geopotential_in_model_layers             ", &
      "air_temperature                          ", &
      "x_wind                                   ", &
      "y_wind                                   ", &
      "omega                                    ", &
      "mass_fraction_of_tracers_in_air          ", &
      "area_type                                ", &
      "atmosphere_boundary_layer_thickness      ", &
      "cell_area                                ", &
      "convective_rainfall_amount               ", &
      "exchange                                 ", &
      "friction_velocity                        ", &
      "rainfall_amount                          ", &
      "soil_moisture_content                    ", &
      "surface_downwelling_shortwave_flux_in_air", &
      "surface_mask                             ", &
      "surface_skin_temperature                 ", &
      "surface_upward_sensible_heat_flux        ", &
      "thickness_of_snowfall_amount             ", &
      "vegetation_type                          ", &
      "vegetation_area_fraction                 ", &
      "z_over_l                                 "  &
    /)
#endif

#if 0
#ifdef ORIG_COORD
  integer, parameter :: importFieldCount = 32
#else
  integer, parameter :: importFieldCount = 30
#endif
  character(len=*), dimension(importFieldCount), parameter :: &
    importFieldNames = (/ &
#ifdef ORIG_COORD
      "deg_lon", &
      "deg_lat", &
#endif
      "grid_index", &
      "grid_perm", &
      "surface_mask", &
      "cell_area", &
      "area_type", &
      "vegetation_type", &
      "vegetation_area_fraction", &
      "thickness_of_snowfall_amount", &
      "atmosphere_boundary_layer_thickness", &
      "atmosphere_optical_thickness_due_to_ambient_aerosol", &
      "convective_rainfall_amount", &
      "accumulated_convective_rainfall_amount", &
      "rainfall_amount", &
      "accumulated_rainfall_amount", &
      "surface_skin_temperature", &
      "surface_downwelling_shortwave_flux_in_air", &
      "surface_upward_sensible_heat_flux", &
      "x_wind_at_10m", &
      "y_wind_at_10m", &
      "friction_velocity", &
      "z_over_l", &
      "air_temperature", &
      "exchange", &
      "omega", &
      "x_wind", &
      "y_wind", &
      "air_pressure", &
      "geopotential", &
      "soil_moisture_content", &
      "tracers" &
      /)
  integer, parameter :: importFieldCount = 2
  character(len=*), dimension(importFieldCount), parameter :: &
    importFieldNames = (/ &
      "inst_zonal_wind_height10m", &
      "inst_merid_wind_height10m"  &
      /)
#endif

  integer :: localDeCount

  private
  
  public SetServices
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine SetServices(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(model, inheritModel, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Provide InitializeP0 to switch from default IPDv00 to IPDv03
    call ESMF_GridCompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!   call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
!     phaseLabelList=(/"IPDv03p3"/), userRoutine=InitializeP3, rc=rc)
!   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!     line=__LINE__, &
!     file=__FILE__)) &
!     return  ! bail out

!   call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
!     phaseLabelList=(/"IPDv03p4"/), userRoutine=InitializeP4, rc=rc)
!   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!     line=__LINE__, &
!     file=__FILE__)) &
!     return  ! bail out

!   call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
!     phaseLabelList=(/"IPDv03p5"/), userRoutine=InitializeP5, rc=rc)
!   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!     line=__LINE__, &
!     file=__FILE__)) &
!     return  ! bail out
    
    ! attach specializing method(s)
    call NUOPC_CompSpecialize(model, specLabel=label_DataInitialize, &
      specRoutine=DataInitialize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(model, specLabel=label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(model, specLabel=label_Finalize, &
      specRoutine=ModelFinalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! attach specializing method(s)
    ! -> NUOPC specializes by default --->>> first need to remove the default
!   if (importFieldCount > 1) then
!     call ESMF_MethodRemove(model, label_CheckImport, rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!       line=__LINE__, &
!       file=__FILE__)) &
!       return  ! bail out
!     call NUOPC_CompSpecialize(model, specLabel=label_CheckImport, &
!       specRoutine=CheckImport, rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!       line=__LINE__, &
!       file=__FILE__)) &
!       return  ! bail out
!   end if

  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP0(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS

    ! Switch to IPDv03 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(model, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv03p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine
  
  !-----------------------------------------------------------------------------
  subroutine InitializeP1(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS

    ! -- advertise imported fields
    if (importFieldCount > 0) then
      call NUOPC_Advertise(importState, importFieldNames, &
        TransferOfferGeomObject="cannot provide", &
        TransferOfferField="cannot provide", &
        SharePolicyField="share", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end if

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine InitializeP3(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    type(ESMF_Field) :: field
    character(len=ESMF_MAXSTR) :: connectedValue, transferAction
    integer :: item

    rc = ESMF_SUCCESS

    do item = 1, importFieldCount
      call ESMF_StateGet(importState, field=field, &
        itemName=trim(importFieldNames(item)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call NUOPC_GetAttribute(field, name="Connected", &
        value=connectedValue, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if (trim(connectedValue) == "false") then
        call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
          msg="Field "//trim(importFieldNames(item))&
          //" in importState is not connected.", &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)
        return ! bail out
      else
        call NUOPC_GetAttribute(field, name="TransferActionGeomObject", &
          value=transferAction, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if (trim(transferAction)=="provide") then
          ! the Connector instructed the Mediator to provide geom object
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg="Cannot fulfill request to provide geom object for "// &
            trim(importFieldNames(item))//" in importState", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
          return ! bail out
        end if
      end if
    end do

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine InitializeP5(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    type(ESMF_Field)               :: field
    type(ESMF_FieldStatus_flag)    :: fieldStatus
    type(ESMF_GeomType_flag)       :: geomtype
    integer                        :: item, itemCount, localDeCount
    integer, dimension(:), pointer :: ugLBound, ugUBound, gridToFieldMap
    character(len=ESMF_MAXSTR)     :: value
    type(ESMF_Array) :: array

    rc = ESMF_SUCCESS

    do item = 1, importFieldCount
      ! -- retrieve field object
      call ESMF_StateGet(importState, field=field, &
        itemName=trim(importFieldNames(item)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      ! -- check field status
      call ESMF_FieldGet(field, status=fieldStatus, localDeCount=localDeCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

!     if (localDeCount < 1) then
!       call ESMF_StateRemove(importState, (/trim(importFieldNames(item))/), rc=rc)
!       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, &
!         file=__FILE__)) &
!         return  ! bail out
!     end if

      if (fieldStatus == ESMF_FIELDSTATUS_GRIDSET) then
        write(6,'("InitializeP5: ",a," status is GRIDSET")') trim(importFieldNames(item))
      else if (fieldStatus == ESMF_FIELDSTATUS_COMPLETE) then
        write(6,'("InitializeP5: ",a," status is COMPLETE - localDeCount=",i0)') trim(importFieldNames(item)), localDeCount
        call ESMF_FieldGet(field, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        write(6,'("InitializeP5: ",a," status is COMPLETE - got Array")') trim(importFieldNames(item))
        write(6,'("InitializeP5: ",a," status is COMPLETE - is Array created?",l8)') trim(importFieldNames(item)), ESMF_ArrayIsCreated(array)
      else if (fieldStatus == ESMF_FIELDSTATUS_EMPTY) then
        write(6,'("InitializeP5: ",a," status is EMPTY")') trim(importFieldNames(item))
        cycle
      else
        write(6,'("InitializeP5: ",a," status is UNKNOWN")') trim(importFieldNames(item))
      end if

      value=""
      call NUOPC_GetAttribute(field, name="SharePolicyGeomObject", &
        value=value, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      write(6,'("InitializeP5: ",a," SharePolicy is ",a)') trim(importFieldNames(item)), trim(value)

      value=""
      call NUOPC_GetAttribute(field, name="ShareStatusGeomObject", &
        value=value, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      write(6,'("InitializeP5: ",a," ShareStatus is ",a)') trim(importFieldNames(item)), trim(value)

      value=""
      call NUOPC_GetAttribute(field, name="TransferActionField", &
        value=value, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      write(6,'("InitializeP5: ",a," TransferActionField is ",a)') trim(importFieldNames(item)), trim(value)

      value=""
      call NUOPC_GetAttribute(field, name="TransferActionGeomObject", &
        value=value, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      write(6,'("InitializeP5: ",a," TransferActionGeomObject is ",a)') trim(importFieldNames(item)), trim(value)

      call ESMF_FieldGet(field, geomtype=geomtype, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      if (fieldStatus == ESMF_FIELDSTATUS_GRIDSET) then
        ! the Connector instructed the Mediator to accept geom object
        ! the transferred geom object is already set, allocate memory
        ! for data by complete
        nullify(ugLBound, ugUBound, gridToFieldMap)
        ! deal with gridToFieldMap
        call ESMF_AttributeGet(field, name="GridToFieldMap", &
          convention="NUOPC", purpose="Instance", &
          itemCount=itemCount, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if (itemCount > 0) then
          allocate(gridToFieldMap(itemCount))
          call ESMF_AttributeGet(field, name="GridToFieldMap", &
            convention="NUOPC", purpose="Instance", &
            valueList=gridToFieldMap, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        endif
        ! deal with ungriddedLBound
        call ESMF_AttributeGet(field, name="UngriddedLBound", &
          convention="NUOPC", purpose="Instance", &
          itemCount=itemCount, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if (itemCount > 0) then
          allocate(ugLBound(itemCount))
          call ESMF_AttributeGet(field, name="UngriddedLBound", &
            convention="NUOPC", purpose="Instance", &
            valueList=ugLBound, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        end if
        ! deal with ungriddedUBound
        call ESMF_AttributeGet(field, name="UngriddedUBound", &
          convention="NUOPC", purpose="Instance", &
          itemCount=itemCount, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if (itemCount > 0) then
          allocate(ugUBound(itemCount))
          call ESMF_AttributeGet(field, name="UngriddedUBound", &
            convention="NUOPC", purpose="Instance", &
            valueList=ugUBound, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        end if

        if (associated(ugLBound) .and. associated(ugUBound)) then
          if (associated(gridToFieldMap)) then
            call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, &
              ungriddedLBound=ugLBound, ungriddedUBound=ugUBound, &
              gridToFieldMap=gridToFieldMap, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out
          else
            call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, &
              ungriddedLBound=ugLBound, ungriddedUBound=ugUBound, &
              rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out
          end if
        else
          call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        end if
      end if
    end do

  end subroutine InitializeP5

  !-----------------------------------------------------------------------------

  subroutine DataInitialize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)              :: importState
    type(ESMF_Field)              :: field
    type(ESMF_Clock)              :: clock
    type(ESMF_Grid)               :: grid
    type(ESMF_VM)                 :: vm
    type(ESMF_GeomType_flag)      :: geomtype
    type(ESMF_DistGrid)           :: distgrid
    type(ESMF_Array)              :: array
    integer                       :: de, item, localrc, localDe, tile
    integer                       :: comm, localPet
    real(ESMF_KIND_R8), dimension(:,:), pointer :: coord

    integer :: dimCount, tileCount, deCount
    integer, dimension(:),   allocatable :: deToTileMap, localDeToDeMap
    integer, dimension(:,:), allocatable :: minIndexPDe, maxIndexPDe, minIndexPTile, maxIndexPTile
    integer, dimension(:,:), allocatable :: computationalLBound, computationalUBound

    integer :: yy, mm, dd, h, m
    real(ESMF_KIND_R8) :: dts
    type(ESMF_Time) :: startTime
    type(ESMF_TimeInterval) :: TimeStep

    ! -- initialize chemistry model
    call ESMF_GridCompGet(model, vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! -- check if import fields are defined
    if (importFieldCount < 1) then 
      call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
        msg="This component requires imported fields to be defined.", &
        line=__LINE__, file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    ! -- query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, importState=importState, modelClock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! get coordinates from Grid object
    ! assume all fields on same grid
    ! use first field 
    call ESMF_StateGet(importState, field=field, &
      itemName=trim(importFieldNames(1)), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_FieldGet(field, geomtype=geomtype, localDeCount=localDeCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write(6,'("-- localDeCount =",i0)') localDeCount
    if (geomtype == ESMF_GEOMTYPE_GRID) then
      call ESMF_FieldGet(field, grid=grid, array=array, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_ArrayGet(array, deCount=deCount, dimCount=dimCount, &
        tileCount=tileCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      allocate(minIndexPDe(dimCount, deCount), maxIndexPDe(dimCount, deCount),  &
        minIndexPTile(dimCount, tileCount), maxIndexPTile(dimCount, tileCount), &
        computationalLBound(dimCount, localDeCount), computationalUBound(dimCount, localDeCount), &
        deToTileMap(deCount), localDeToDeMap(localDeCount), stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_ArrayGet(array, distgrid=distgrid, &
        deToTileMap=deToTileMap, localDeToDeMap=localDeToDeMap, &
        computationalLBound=computationalLBound, &
        computationalUBound=computationalUBound, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_DistGridGet(distgrid, &
        minIndexPDe=minIndexPDe, maxIndexPDe=maxIndexPDe, &
        minIndexPTile=minIndexPTile, maxIndexPTile=maxIndexPTile, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! -- init chemistry model on local DEs
      call chem_model_create(deCount=localDeCount, rc=rc)
      if (chem_rc_check(rc)) then
        call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
          msg="Failed to initialize chemistry model for localDeCount", &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)
        return  ! bail out
      end if

      call ESMF_VMGet(vm, localPet=localPet, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      do localDe = 0, localDeCount-1
        de   = localDeToDeMap(localDe+1) + 1
        tile = deToTileMap(de)

!       write(6,'("CHEM on localDe: ",i0," DE: ",i0)') localDe, de-1
        write(6,'("CHEM: PET:",i04," localDe: ",i4," DE: ",i4, " tile=",i2," minIndexPDe=",2i4,2x," maxIndexPDe=",2i4," minIndexPTile=",2i4," maxIndexPTile=",2i4,4i4)') &
          localPet, localDe, de-1, tile, minIndexPDe(:,de), maxIndexPDe(:,de), minIndexPTile(:,tile), maxIndexPTile(:,tile), &
          computationalLBound(:,localDe+1), computationalUBound(:,localDe+1)
        flush(6)
!       write(6,'(4x,"minIndexPDe=",2i8,2x,"maxIndexPDe=",2i8)') minIndexPDe(:,de), maxIndexPDe(:,de)
!       write(6,'(4x,"tile=",i0)') deToTileMap(de)

        ! -- set model domains for local DEs
        call chem_model_domain_set(minIndexPDe=minIndexPDe(:,de), maxIndexPDe=maxIndexPDe(:,de), &
          minIndexPTile=minIndexPTile(:,tile), maxIndexPTile=maxIndexPTile(:,tile), &
          minIndexLocal=computationalLBound(:,localDe+1), maxIndexLocal=computationalUBound(:,localDe+1), &
          tile=deToTileMap(de), tileCount=tileCount, de=localDe, rc=rc)
        if (chem_rc_check(rc)) then
          call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
            msg="Failed to initialize chemistry model for localDeCount", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
        end if

        ! -- get local coordinate arrays
        do item = 1, 2
          call ESMF_GridGetCoord(grid, coordDim=item, staggerloc=ESMF_STAGGERLOC_CENTER, &
            localDE=localDe, farrayPtr=coord, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          call chem_model_domain_coord_set(item, coord, de=localDe, rc=rc)
          if (chem_rc_check(rc)) then
            call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
              msg="Failed to set coordinates for chemistry model", &
              line=__LINE__, &
              file=__FILE__, &
              rcToReturn=rc)
            return  ! bail out
          end if
        end do

      end do
      deallocate(minIndexPDe, maxIndexPDe, minIndexPTile, maxIndexPTile, &
        computationalLBound, computationalUBound, &
        deToTileMap, localDeToDeMap, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

    else
      call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
        msg="Imported fields can only be defined on Grid objects.", &
        line=__LINE__, file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    ! -- init model communication subsystem over ESMF communicator
    call ESMF_GridCompGet(model, vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_VMGet(vm, mpiCommunicator=comm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! -- initialize model after creation to setup correct communication
    print *, ' -- chemCap: before chem_model_init() '
    call chem_model_init(comm=comm, isolate=.true., rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to initialize chemistry model", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if
    print *, ' -- chemCap: after chem_model_init() '

    ! -- initialize model I/O
    call chem_io_init(rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to initialize I/O model subsystem", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if
    print *, ' -- chemCap: after chem_io_init() '

    ! -- read-in emission and background fields, setup internal parameters
    call chem_model_config_init(rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to initialize model configuration", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if
    print *, ' -- chemCap: after chem_model_config_init() '

    ! -- initialize internal clock
    ! -- get clock information
    call ESMF_ClockGet(clock, startTime=startTime, timeStep=timeStep, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    ! -- get forecast initial time
    call ESMF_TimeGet(startTime, yy=yy, mm=mm, dd=dd, h=h, m=m, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    ! -- get time step
    call ESMF_TimeIntervalGet(timeStep, s_r8=dts, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    ! -- set internal clock
    call chem_model_clock_create(rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to create model clock", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if
    call chem_model_clock_set(yy=yy, mm=mm, dd=dd, h=h, m=m, dts=dts, tz=0, rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to initialize model clock", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if
    print *, ' -- chemCap: after chem_model_clock_init() '

    ! -- allocate memory for internal workspace
    call chem_backgd_init(rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to initialize I/O model subsystem", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if
    print *, ' -- chemCap: after chem_backgd_init() '

    ! -- read-in emission and background fields
    call chem_backgd_read(rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to initialize I/O model subsystem", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if
    print *, ' -- chemCap: after chem_backgd_read() '

#if 0
    ! -- diagnostics: write out emission and background fields
    call chem_backgd_write(rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to initialize I/O model subsystem", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if
    print *, ' -- chemCap: after chem_backgd_read() '
#endif

    ! -- connect import fields to model
    ! -- this can be done only once since remote fields are accessed by reference
    call chem_comp_connect('import', importState, importFieldNames, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    print *, ' -- chemCap: after chem_comp_connect() '

    ! -- initialize internal component (GOCART)
    call chem_comp_init(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    print *, ' -- chemCap: after chem_comp_init() '

    ! -- allocate arrays for tracer output
    call chem_output_init(rc=rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to allocate output arrays", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if
    print *, ' -- chemCap: after chem_output_init() '

    ! indicate that data initialization is complete (breaking out of init-loop)
    call NUOPC_CompAttributeSet(model, &
      name="InitializeDataComplete", value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine DataInitialize

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_State)              :: importState, exportState
    type(ESMF_Time)               :: currTime
    type(ESMF_TimeInterval)       :: timeStep
    type(ESMF_Field)              :: field
    type(ESMF_VM)                 :: vm
    integer                       :: item

    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    
    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in 
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.
    
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing CHM from: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_TimePrint(currTime + timeStep, &
      preString="--------------------------------> to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! -- debug
    call ESMF_GridCompGet(model, vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! -- print field diagnostics
    do item = 1, importFieldCount
      call ESMF_StateGet(importState, field=field, &
        itemName=trim(importFieldNames(item)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call fieldPrintMinMax(field, vm, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end do

    ! -- advance model
    call chem_comp_advance(clock, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! -- write output if time
    call chem_output_write(rc)
    if (chem_rc_check(rc)) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
        msg="Failed to allocate output arrays", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

  end subroutine ModelAdvance


  subroutine ModelFinalize(model, rc)
    type(ESMF_GridComp)   :: model
    integer, intent(out)  :: rc

    rc = ESMF_SUCCESS

    call chem_comp_finalize(rc)

  end subroutine ModelFinalize

  !-----------------------------------------------------------------------------

  subroutine CheckImport(model, rc)
    type(ESMF_GridComp)   :: model
    integer, intent(out)  :: rc
    
    ! local variables
    type(ESMF_Clock)                :: clock
    type(ESMF_Time)                 :: currTime, invalidTime
    type(ESMF_State)                :: importState
    logical                         :: timeCheck
    type(ESMF_Field),       pointer :: fieldList(:)
    integer                         :: i
    
    rc = ESMF_SUCCESS

    ! query the Component for its clock and importState
    call ESMF_GridCompGet(model, clock=clock, importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! get the current time out of the clock
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set up invalid time (by convention)
    call ESMF_TimeSet(invalidTime, yy=99999999, mm=01, dd=01, &
      h=00, m=00, s=00, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! Loop through all the field in the importState, and test whether they
    ! are at invalidTime (ignore them for now), or at currTime. Any other
    ! time coming in would flag an incompatibility.
    
    nullify(fieldList)
    call NUOPC_GetStateMemberLists(importState, fieldList=fieldList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    do i=1, size(fieldList)
      timeCheck = NUOPC_IsAtTime(fieldList(i), invalidTime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if (timeCheck) then
        ! The field is at invalidTime

        ! -> In a real application mark the field with a flag as invalid 
        !    so the actual model code can act accordingly.

        ! Here for purpose of demonstration just log a message and continue on.
      
        call ESMF_LogWrite("CHM: detected import field at invalidTime", &
          ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      else
        ! The field is NOT at invalidTime -> it must then be at currTime or it
        ! is incompatible!
      
        ! check that Fields in the importState show correct timestamp
        timeCheck = NUOPC_IsAtTime(fieldList(i), currTime, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if (.not.timeCheck) then
          !TODO: introduce and use INCOMPATIBILITY return codes!!!!
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg="NUOPC INCOMPATIBILITY DETECTED: "//&
            "Import Field not at current time", &
            line=__LINE__, file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
        endif
      
      endif
    enddo
    
  end subroutine

end module CHM
