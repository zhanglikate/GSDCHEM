module chem_methods

  use ESMF
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
!   integer :: clockId
    integer :: its, julday, month
    integer(ESMF_KIND_I8) :: advanceCount
    real(ESMF_KIND_R8) :: dts
!   character(len=clock_ts_len) :: datetime
    type(ESMF_Time) :: currTime
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_CalKind_Flag) :: calkindflag

    ! -- begin
    rc = ESMF_SUCCESS

    ! -- start timer
!   clockId = mpp_clock_id( 'Chemistry advance')
!   call mpp_clock_begin(clockId)

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

    integer :: clockId

    rc = ESMF_SUCCESS
!   clockId = mpp_clock_id( 'Chemistry finalize')
!   call mpp_clock_begin(clockId)
!   call chem_free_memory()
    if (allocated(chem_pelist)) deallocate(chem_pelist)
!   call mpp_clock_end(clockId)
  end subroutine chemFinalize

  subroutine chemDomainSet(js, je, jmax, kmax)

    integer, optional, intent(in) :: js, je, jmax, kmax

    ids =   1 ; ide =   1 
    ims =   1 ; ime =   1 
    its =   1 ; ite =   1

    if (present(js)) then
      jts = js
      jms = js
    end if
    if (present(je)) then
      jte = je
      jme = je
    end if
    if (present(jmax)) then
      jds = 1
      jde = jmax
    end if
    if (present(kmax)) then
      nvl = kmax
      nvlp1 = nvl + 1
      kts = 1
      kds = 1
      kms = 1
      kte = nvl
      kde = nvlp1
      kme = nvlp1
    end if

  end subroutine chemDomainSet

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

end module chem_methods
