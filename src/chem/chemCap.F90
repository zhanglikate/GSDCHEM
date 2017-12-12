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
#else
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

    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p3"/), userRoutine=InitializeP3, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!   call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
!     phaseLabelList=(/"IPDv03p4"/), userRoutine=InitializeP4, rc=rc)
!   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!     line=__LINE__, &
!     file=__FILE__)) &
!     return  ! bail out

    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p5"/), userRoutine=InitializeP5, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
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
    if (importFieldCount > 1) then
      call ESMF_MethodRemove(model, label_CheckImport, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call NUOPC_CompSpecialize(model, specLabel=label_CheckImport, &
        specRoutine=CheckImport, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end if

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
    integer                        :: item, itemCount
    integer, dimension(:), pointer :: ugLBound, ugUBound, gridToFieldMap

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
      call ESMF_FieldGet(field, status=fieldStatus, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

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
            call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R4, &
              ungriddedLBound=ugLBound, ungriddedUBound=ugUBound, &
              gridToFieldMap=gridToFieldMap, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out
          else
            call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R4, &
              ungriddedLBound=ugLBound, ungriddedUBound=ugUBound, &
              rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out
          end if
        else
          call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R4, rc=rc)
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
    type(ESMF_Mesh)               :: mesh
    type(ESMF_Grid)               :: grid
    type(ESMF_VM)                 :: vm
    type(ESMF_GeomType_flag)      :: geomtype
    type(ESMF_CoordSys_Flag)      :: coordSys
    type(ESMF_DistGrid)           :: distgrid
    integer                       :: item, localrc, localDe, numOwnedNodes, spatialDim
    integer                       :: comm, localPet
    integer                       :: recvData(1)
    real, allocatable :: fperm(:)
    real(ESMF_KIND_R8), dimension(:), allocatable :: ownedNodeCoords, lon, lat

    integer :: dimCount, tileCount, deCount, elementCount
    integer, dimension(:),   allocatable :: seqIndexList
    integer, dimension(:,:), allocatable :: minIndexPDe, maxIndexPDe

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, importState=importState, modelClock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! initialize chemistry model
    is_nuopc = .true.

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

    call chemInitialize(clock, phase=0, comm=comm, rc=rc)

    ! -- mark import fields as updated
    do item = 1, importFieldCount
      call ESMF_StateGet(importState, field=field, &
        itemName=trim(importFieldNames(item)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call NUOPC_SetAttribute(field, name="Updated", value="true", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end do

    ! get coordinates from Grid/Mesh object
    ! assume all fields on same grid/mesh
    ! use last field from previous loop, if available
    if (importFieldCount > 0) then
      call ESMF_FieldGet(field, geomtype=geomtype, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if      (geomtype == ESMF_GEOMTYPE_MESH) then
        call ESMF_FieldGet(field, mesh=mesh, localDeCount=localDeCount, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        print *,'-- FIELD distgrid: localDeCount = ',localDeCount
        if (localDeCount > 1) then
          call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
            msg="Imported fields can only have localDeCount = 1",&
            line=__LINE__, file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
        end if
        call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedNodes=numOwnedNodes, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if (spatialDim /= 2) then
          call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
            msg="Imported Mesh can only have spatialDim = 2.", &
            line=__LINE__, file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
        end if
#ifndef ORIG_COORD
        allocate(ownedNodeCoords(spatialDim*numOwnedNodes), stat=localrc)
        if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)) return
        call ESMF_MeshGet(mesh, ownedNodeCoords=ownedNodeCoords, coordSys=coordSys, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        allocate(lon(numOwnedNodes), lat(numOwnedNodes), stat=localrc)
        if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)) return
        if      (coordSys == ESMF_COORDSYS_SPH_RAD) then
          lon = ESMF_COORDSYS_RAD2DEG * ownedNodeCoords(1::2)
          lat = ESMF_COORDSYS_RAD2DEG * ownedNodeCoords(2::2)
        else if (coordSys == ESMF_COORDSYS_SPH_DEG) then
          lon = ownedNodeCoords(1::2)
          lat = ownedNodeCoords(2::2)
        end if
        deallocate(ownedNodeCoords, stat=localrc)
        if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)) return
#endif
        ! -- get total number of nodes
        call ESMF_GridCompGet(model, vm=vm, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        recvData = 0
        call ESMF_VMAllReduce(vm, (/ numOwnedNodes /), recvData, 1, ESMF_REDUCE_SUM, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        nip = recvData(1)
        glvl = int(log((nip-2.)/10.)/log(4.))
        ! -- get index decomposition
        call ESMF_MeshGet(mesh, nodalDistgrid=distgrid, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_DistGridGet(distgrid, dimCount=dimCount, tileCount=tileCount, deCount=deCount, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if (tileCount > 1) then
          call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
            msg="Imported Mesh should have 1 tile.", &
            line=__LINE__, file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
        end if
        if (dimCount > 1) then
          call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
            msg="Imported Mesh should have 1 dimension.", &
            line=__LINE__, file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
        end if
        call ESMF_DistGridGet(distgrid, localDe=0, elementCount=elementCount, rc=rc)
        allocate(seqIndexList(elementCount), stat=localrc)
        seqIndexList = 0
        if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)) return
        call ESMF_DistGridGet(distgrid, localDe=0, seqIndexList=seqIndexList, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        ! -- set horizontal domain decomposition
        call chemDomainSet(js=minval(seqIndexList), je=maxval(seqIndexList), jmax=nip)
        deallocate(seqIndexList, stat=localrc)
        if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)) return
      else if (geomtype == ESMF_GEOMTYPE_GRID) then
        call ESMF_FieldGet(field, grid=grid, localDeCount=localDeCount, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        write(6,'("-- localDeCount =",i0)') localDeCount
        call ESMF_GridGet(grid, distgrid=distgrid, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_DistGridGet(distgrid, dimCount=dimCount, tileCount=tileCount, deCount=deCount, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        write(6,'("-- GRID distgrid: dimCount = ",i0," tileCount = ",i0," deCount = ",i0)') dimCount, tileCount, deCount
        if (dimCount /= 2) then
          call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
            msg="Imported Grid should have 2 dimensions.", &
            line=__LINE__, file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
        end if
        allocate(minIndexPDe(dimCount, deCount), maxIndexPDe(dimCount, deCount), stat=localrc)
        if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
        minIndexPDe = 0
        maxIndexPDe = 0
        call ESMF_DistGridGet(distgrid, &
          minIndexPDe=minIndexPDe, &
          maxIndexPDe=maxIndexPDe, &
          rc=rc)
        write(6,'("-- GRID distgrid: minIndexPDe   = ",20i8)') minIndexPDe
        write(6,'("-- GRID distgrid: maxIndexPDe   = ",20i8)') maxIndexPDe
        deallocate(minIndexPDe, maxIndexPDe, stat=localrc)
        if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      else
        call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
          msg="Imported fields can only be defined on Grid or Mesh objects.", &
          line=__LINE__, file=__FILE__, &
          rcToReturn=rc)
        return  ! bail out
      end if
    end if

    if (importFieldCount > 0) then
      if      (geomtype == ESMF_GEOMTYPE_MESH) then
        ! -- if single DE/PET, set pointers here once and for all
        call importFieldSet(importState, importFieldNames, 0, rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        ! -- set vertical domain
        call chemDomainSet(kmax=size(tr3d, dim=1))

        ! -- get index order and permutation information
        call ESMF_VMGet(vm, localPet=localPet, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        allocate(indx(nip), perm(nip), fperm(nip))
        indx = 0
        perm = 0
        fperm = 0.
        call ESMF_StateGet(importState, field=field, itemName="grid_index", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_FieldGather(field, fperm, 0, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if (localPet == 0) then
          indx = nint(fperm)
        end if
        fperm = 0.
        call ESMF_StateGet(importState, field=field, itemName="grid_perm", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_FieldGather(field, fperm, 0, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if (localPet == 0) then
          perm = nint(fperm)
        end if
        deallocate(fperm)

        ! -- create mapping array for FIM node indices (FIM-to-local)
        allocate(map(nip))
        map = 0
        if (localPet == 0) then
          do item = 1, nip
            map(indx(item)) = item
          end do
          allocate(seqIndexList(nip))
          seqIndexList = 0
          do item = 1, nip
            seqIndexList(item) = perm(map(item))
          end do
          perm = seqIndexList
          deallocate(seqIndexList)
        end if

        call ESMF_VMBroadcast(vm, indx, nip, 0, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        call ESMF_VMBroadcast(vm, map, nip, 0, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        call ESMF_VMBroadcast(vm, perm, nip, 0, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        ! initialize chemistry model (phase 2)
        call chemInitialize(clock, phase=1)

#ifndef ORIG_COORD
        ! set model's grid coordinates
        deg_lon(jms:jme) = lon
        deg_lat(jms:jme) = lat

        deallocate(lon, lat, stat=localrc)
        if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)) return
#endif

        ! prepare additional arrays
        call chem_data_prep()
      end if
    end if

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
    integer                       :: item, localDe

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

#if 0
    ! -- if imported Fields are defined on multiple DEs, reset pointers for each localDE
    if (localDeCount > 1) then
      do localDe = 0, localDeCount - 1
        call importFieldSet(importState, importFieldNames, localDe, rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      end do
    end if
#endif

    return

    ! -- prepare additional arrays
    call chem_data_prep()

    ! -- advance model
    call chemAdvance(clock, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  subroutine ModelFinalize(model, rc)
    type(ESMF_GridComp)   :: model
    integer, intent(out)  :: rc

    rc = ESMF_SUCCESS

    call chemFinalize(rc)

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

  !-----------------------------------------------------------------------------

  subroutine fieldPrintMinMax(field, vm, rc)
    type(ESMF_Field), intent(in) :: field
    type(ESMF_VM),    intent(in) :: vm
    integer, intent(out) :: rc

    ! local variables
    real(ESMF_KIND_R4), pointer :: fp1d(:), fp2d(:,:), fp3d(:,:,:), fp4d(:,:,:,:)
    real(ESMF_KIND_R4)          :: fieldMaxValue, fieldMinValue, maxValue, minValue
    real(ESMF_KIND_R4)          :: globalMaxValue(1), globalMinValue(1)
    integer                     :: localDe, localDeCount, localPet, rank
    integer                     :: i, j
    integer, dimension(2)       :: clb, cub
    character(len=ESMF_MAXSTR)  :: fieldName

    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field, rank=rank, localDeCount=localDeCount, &
      name=fieldName, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail

    fieldMinValue = huge(1.0)
    fieldMaxValue = -fieldMinValue
    do localDe = 0, localDeCount - 1
      select case(rank)
        case(1)
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp1d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          minValue = minval(fp1d)
          maxValue = maxval(fp1d)
#if 0
        case(2)
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp2d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          minValue = minval(fp2d)
          maxValue = maxval(fp2d)
#else
        case(2)
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp2d, &
            computationalLBound=clb, computationalUBound=cub, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          do j = clb(2), cub(2)
            do i = clb(1), cub(1)
              write(6,'("-- data",2i8,4x,g16.6)') i,j,fp2d(i,j)
            end do
          end do
          minValue = minval(fp2d(clb(1):cub(1),clb(2):cub(2)))
          maxValue = maxval(fp2d(clb(1):cub(1),clb(2):cub(2)))
          write(6,'("-- dataend Min/Max = ",2g16.6)') minValue, maxValue
#endif
        case(3)
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp3d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          minValue = minval(fp3d)
          maxValue = maxval(fp3d)
        case(4)
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp4d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          minValue = minval(fp4d)
          maxValue = maxval(fp4d)
        case default
          call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
            msg="Field rank not implemented.", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
          return ! bail out
      end select
      fieldMinValue = min(fieldMinValue, minValue)
      fieldMaxValue = max(fieldMaxValue, maxValue)
      write(6,'(a,":",i0,2x,"DE: ",i0,2x,a," - checking  - min/max = ",2g16.6)') 'PET', &
         localPet, localDe, trim(fieldName), fieldMinValue, fieldMaxValue
    end do

    globalMinValue(1) = 0.
    globalMaxValue(1) = 0.

    call ESMF_VMReduce(vm, (/ fieldMinValue /), globalMinValue, 1, &
      reduceflag=ESMF_REDUCE_MIN, rootPet=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_VMReduce(vm, (/ fieldMaxValue /), globalMaxValue, 1, &
      reduceflag=ESMF_REDUCE_MAX, rootPet=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail

    if (localPet == 0) then
       write(6,'(a,":",a," - checking  - min/max = ",2g16.6)') 'Field', &
         trim(fieldName), globalMinValue, globalMaxValue
    end if

  end subroutine fieldPrintMinMax

  !-----------------------------------------------------------------------------

  subroutine importFieldSet(importState, fieldNames, localDe, rc)
    type(ESMF_State) :: importState
    character(len=*), dimension(:), intent(in) :: fieldNames
    integer, intent(in) :: localDe
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Field) :: field
    integer          :: item

    rc = ESMF_SUCCESS

    do item = 1, size(fieldNames)
      call ESMF_StateGet(importState, field=field, &
        itemName=trim(fieldNames(item)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail
      select case (trim(fieldNames(item)))
#ifdef ORIG_COORD
        case ("deg_lon")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=deg_lon, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("deg_lat")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=deg_lat, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
#endif
        case ("grid_perm")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=loc_perm, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("surface_mask")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=slmsk2d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("cell_area")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=area, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("area_type")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=stype2d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("vegetation_type")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=vtype2d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("vegetation_area_fraction")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=vfrac2d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("thickness_of_snowfall_amount")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=snwdph2d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("atmosphere_boundary_layer_thickness")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=pb2d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
!       case ("atmosphere_optical_thickness_due_to_ambient_aerosol")
!         call ESMF_FieldGet(field, localDe=localDe, farrayPtr=?????, rc=rc)
!         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!           line=__LINE__, &
!           file=__FILE__)) &
!           return  ! bail
        case ("convective_rainfall_amount")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=rc2d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("rainfall_amount")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=rn2d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("surface_skin_temperature")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=ts2d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("surface_downwelling_shortwave_flux_in_air")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=rsds, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("surface_upward_sensible_heat_flux")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=hf2d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("friction_velocity")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=us2d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("z_over_l")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=zorl2d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("air_temperature")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=tk3d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("exchange")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=exch, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("omega")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=ws3d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("x_wind")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=us3d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("y_wind")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=vs3d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("air_pressure")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=pr3d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("geopotential")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=ph3d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("soil_moisture_content")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=sm3d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
        case ("tracers")
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=tr3d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail
      end select
      call NUOPC_SetAttribute(field, name="Updated", value="true", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail
    end do

  end subroutine importFieldSet

  !-----------------------------------------------------------------------------

end module CHM
