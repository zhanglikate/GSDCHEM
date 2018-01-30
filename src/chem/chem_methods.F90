module chem_methods

  use ESMF
  use NUOPC
  use chem_rc_mod
  use chem_comm_mod
  use chem_shr_mod
  use chem_config_mod
  use chem_species_mod
  use chem_domain_mod
  use chem_model_mod
  use chem_vars_mod
  use gocart_mod

  implicit none

  public 

contains

  subroutine chemInitialize(clock, phase, comm, rc)

    type(ESMF_Clock),  intent(in)  :: clock
    integer, optional, intent(in)  :: phase
    integer, optional, intent(in)  :: comm
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: lphase
    integer :: pe, pecount
    integer :: yy, mm, dd, h, m
    integer, dimension(:), allocatable :: pelist
    real(ESMF_KIND_R8) :: dts
    type(ESMF_Time) :: startTime
    type(ESMF_TimeInterval) :: TimeStep
    type(ESMF_CalKind_Flag) :: calkindflag

    ! -- begin
    ! -- start timer
    call chem_comm_init(localrc, comm=comm, isolate=.true.)

    print *,"chem: enter init()"

    ! -- set phase
    lphase = -1
    if (present(phase)) lphase = phase

    if (lphase < 1) then

      call chem_comm_get(pecount=pecount)
      allocate(chem_pelist(pecount), stat=localrc)
      if (localrc /= 0) call chem_comm_abort(msg='CHEM_INIT: failed to allocate memory for pelist')
      do pe = 1, pecount
        chem_pelist(pe) = pe-1
      end do

      ! -- print PE list
      print *,"chem: PEs: ", chem_pelist

      ! -- get configuration
      print *,'chem: reading config namelist..'
      call chem_config_read(chem_config)
      print *,'chem: done reading config namelist'

      ! -- get clock information
      call ESMF_ClockGet(clock, startTime=startTime, timeStep=timeStep, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      ! -- set forecast initial time string
      call ESMF_TimeGet(startTime, yy=yy, mm=mm, dd=dd, h=h, m=m, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      write(yyyymmddhhmm, '(i4.4,4i2.2)') yy, mm, dd, h, m
      ! -- set time step
      call ESMF_TimeIntervalGet(timeStep, s_r8=dts, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      dt = real(dts)
      numphr = nint(3600./dt) ! # of time steps/hr

      call chem_control_setup(chem_config)
      call chem_species_setup(chem_config)

    end if

    if (lphase /= 0) then

      ! -- initialize after field are connected

      call chemAllocate()
      
      call chem_background_init()
      print *,'chem: done chem_background_init'
      print *,'chem: chem_opt = ',chem_config % chem_opt

      select case (chem_config % chem_opt)
        case(CHEM_OPT_GOCART)
          print *,'chem: calling gocart init'
          call gocart_init()
          print *,'chem: done gocart init'
        case default
          call chem_comm_abort(msg='CHEM_INIT: chem_opt value not implemented')
      end select

    end if

    print *,"chem: exit init()"

  end subroutine chemInitialize

  subroutine chemAdvance(clock, rc)

    type(ESMF_Clock), intent(in) :: clock
    integer,         intent(out) :: rc

    ! -- local variables
    integer :: its, julday, month
    integer(ESMF_KIND_I8) :: advanceCount
    real(ESMF_KIND_R8) :: dts
    type(ESMF_Time) :: currTime
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_CalKind_Flag) :: calkindflag

    ! -- begin
    rc = ESMF_SUCCESS

    print *,'chem: run(): pelist = ',chem_pelist

    call ESMF_ClockPrint(clock, &
      preString="chem: run(): time step : ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="chem: run(): time stamp: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, &
      advanceCount=advanceCount, calkindflag=calkindflag, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    its = real(advanceCount) + 1

    call ESMF_TimeIntervalGet(timeStep, s_r8=dts, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    dt = real(dts)
    numphr = nint(3600./dt) ! # of time steps/hr

    call ESMF_TimeGet(currTime, mm=month, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimeGet(currTime, dayOfYear=julday, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    select case (chem_config % chem_opt)
      case (CHEM_OPT_GOCART)
        call gocart_run(its, julday, 0, month, dt)
      case default
    end select

    ! -- write history
    call chem_history_write(chem_config, its)

    ! -- stop timer
!   call mpp_clock_end(clockId)

  end subroutine chemAdvance

  subroutine chemFinalize(rc)
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    call chem_model_destroy()
!   clockId = mpp_clock_id( 'Chemistry finalize')
!   call mpp_clock_begin(clockId)
!   call chem_free_memory()
    if (allocated(chem_pelist)) deallocate(chem_pelist)
!   call mpp_clock_end(clockId)
  end subroutine chemFinalize

  subroutine chemAllocate

    allocate(aod2d   (jms:jme))
#ifndef ORIG_COORD
    allocate(deg_lat (jms:jme))
    allocate(deg_lon (jms:jme))
#endif
    allocate(dp3d(nvl  ,jms:jme)) 
    allocate(trdp(nvl  ,jms:jme,ntra+ntrb))

    allocate(gd_cloud ( 1, nvlp1, jms:jme ))
    allocate(gd_cldfr ( 1, nvlp1, jms:jme ))

    aod2d  = 0.
#ifndef ORIG_COORD
    deg_lat = 0.
    deg_lon = 0.
#endif
    dp3d    = 0.
    trdp    = 0.
    gd_cloud = 0.
    gd_cldfr = 0.

    call chem_alloc_workspace
    call chem_alloc_prep
    call chem_alloc_seas

  end subroutine chemAllocate

  !-----------------------------------------------------------------------------

  subroutine chemFieldsConnect(importState, fieldNames, rc)
    type(ESMF_State),               intent(in) :: importState
    character(len=*), dimension(:), intent(in) :: fieldNames
    integer, intent(out) :: rc

    ! local variables
    type(chem_state_type), pointer :: state
    type(ESMF_Field)               :: field
    integer                        :: item, localDe, localDeCount

    rc = ESMF_SUCCESS

    do item = 1, size(fieldNames)

      write(6,'("CHEM: connecting field ",i0," - ",a)') item, trim(fieldNames(item))
      flush(6)
      call ESMF_StateGet(importState, field=field, &
        itemName=trim(fieldNames(item)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail
      call ESMF_FieldGet(field, localDeCount=localDeCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail

      do localDe = 0, localDeCount-1

        call chem_model_get(stateIn=state, de=localDe, rc=rc)
        if (rc /= CHEM_RC_SUCCESS) then
          call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
            msg="Failed to initialize chemistry component", &
            line=__LINE__, &
            file=__FILE__) 
          return  ! bail out
        end if

        select case (trim(fieldNames(item)))
          case ("air_pressure")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % pr3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
            call chem_model_set(numModLayers=size(state % pr3d,dim=3), de=localDe)
          case ("air_pressure_in_model_layers")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % prl3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
            call chem_model_set(numIntLayers=size(state % prl3d,dim=3), de=localDe)
          case ("air_temperature")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % tk3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("area_type")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % stype2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("atmosphere_boundary_layer_thickness")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % pb2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("cell_area")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % area, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("convective_rainfall_amount")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % rc2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("exchange")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % exch, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("friction_velocity")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % us2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("geopotential")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % ph3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("geopotential_in_model_layers")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % phl3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("mass_fraction_of_tracers_in_air")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % tr3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
            call chem_model_set(numTracers=size(state % tr3d, dim=4), de=localDe)
          case ("omega")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % ws3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("rainfall_amount")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % rn2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("soil_moisture_content")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % sm3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("surface_downwelling_shortwave_flux_in_air")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % rsds, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("surface_mask")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % slmsk2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("surface_skin_temperature")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % ts2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("surface_upward_sensible_heat_flux")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % hf2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("thickness_of_snowfall_amount")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % snwdph2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("vegetation_type")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % vtype2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("vegetation_area_fraction")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % vfrac2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("x_wind")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % us3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("y_wind")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % vs3d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
          case ("z_over_l")
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % zorl2d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail
!         case ("inst_zonal_wind_height10m")
!           call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % u10mi, rc=rc)
!           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!             line=__LINE__, &
!             file=__FILE__)) &
!             return  ! bail
!         case ("inst_merid_wind_height10m")
!           call ESMF_FieldGet(field, localDe=localDe, farrayPtr=state % v10mi, rc=rc)
!           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!             line=__LINE__, &
!             file=__FILE__)) &
!             return  ! bail
        end select
      end do
      call NUOPC_SetAttribute(field, name="Updated", value="true", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail
    end do

  end subroutine chemFieldsConnect

  !-----------------------------------------------------------------------------

  subroutine fieldPrintMinMax(field, vm, rc)
    type(ESMF_Field), intent(in) :: field
    type(ESMF_VM),    intent(in) :: vm
    integer, intent(out) :: rc

    ! local variables
    real(ESMF_KIND_R8), pointer :: fp1d(:), fp2d(:,:), fp3d(:,:,:), fp4d(:,:,:,:)
    real(ESMF_KIND_R8)          :: fieldMaxValue, fieldMinValue, maxValue, minValue, avgValue
    real(ESMF_KIND_R8)          :: globalMaxValue(1), globalMinValue(1)
    integer                     :: localDe, localDeCount, localPet, rank
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

    fieldMinValue = huge(1.0_ESMF_KIND_R8)
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
          avgValue = sum(fp1d)/size(fp1d)
        case(2)
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp2d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          minValue = minval(fp2d)
          maxValue = maxval(fp2d)
          avgValue = sum(fp2d)/size(fp2d)
        case(3)
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp3d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          minValue = minval(fp3d)
          maxValue = maxval(fp3d)
          avgValue = sum(fp3d)/size(fp3d)
        case(4)
          call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp4d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          minValue = minval(fp4d)
          maxValue = maxval(fp4d)
          avgValue = sum(fp4d)/size(fp4d)
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
      write(6,'(a,":",i0,2x,"DE: ",i0,2x,a," - min/max/avg = ",3g16.6)') 'PET', &
         localPet, localDe, trim(fieldName), minValue, maxValue, avgValue
    end do

    globalMinValue(1) = 0._ESMF_KIND_R8
    globalMaxValue(1) = 0._ESMF_KIND_R8

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
       write(6,'(a,":",a," - global min/max = ",2g16.6)') 'Field', &
         trim(fieldName), globalMinValue, globalMaxValue
    end if

  end subroutine fieldPrintMinMax

  !-----------------------------------------------------------------------------

end module chem_methods
