module gocart_model_mod

  use chem_rc_mod
  use chem_comm_mod
  use chem_model_mod
  use chem_tracers_mod
  use gocart_mod

  implicit none

  private

  public :: gocart_model_init
  public :: gocart_model_advance

contains

  ! -- public methods

  subroutine gocart_model_init(rc)
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: deCount
    type(chem_config_type), pointer :: config

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    ! -- get model config
    call chem_model_get(deCount=deCount, config=config, rc=localrc)
    if (chem_rc_check(localrc, msg="Failed to retrieve model", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    if (deCount < 1) return

    ! -- initialize species pointers for GOCART internals
    call chem_tracers_set(config, rc=localrc)
    if (chem_rc_check(localrc, msg="Failed to set tracer pointers", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- initialize GOCART modules
    call gocart_init(config, rc=localrc)
    if (chem_rc_check(localrc, msg="Failed to initialize GOCART modules", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- print out model configuration
    if (chem_comm_isroot()) then
      write(6,'(28("-"))')
      write(6,'("GOCART configuration:")')
      write(6,'(28("-"))')
      write(6,'("    chem_opt         = ",i0)') config % chem_opt
      write(6,'("    chem_in_opt      = ",i0)') config % chem_in_opt
      write(6,'("    dust_opt         = ",i0)') config % dust_opt
      write(6,'("    dmsemis_opt      = ",i0)') config % dmsemis_opt
      write(6,'("    seas_opt         = ",i0)') config % seas_opt
      write(6,'("    biomass_burn_opt = ",i0)') config % biomass_burn_opt
      write(6,'(28("-"))')
    end if

  end subroutine gocart_model_init


  subroutine gocart_model_advance(rc)
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc,tile
    integer :: de, deCount
    integer :: advanceCount, julday, mm, tz
    integer :: is, ie, js, je, ni, nl
    real(CHEM_KIND_R8) :: dts
    real(CHEM_KIND_R8), dimension(:,:), pointer :: lat, lon
    type(chem_config_type), pointer :: config
    type(chem_data_type),   pointer :: data
    type(chem_state_type),  pointer :: stateIn, stateOut

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(deCount=deCount, rc=localrc) 
    if (chem_rc_check(localrc, msg="Failed to retrieve model", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    if (deCount < 1) return

    call chem_model_clock_get(advanceCount=advanceCount, dts=dts, mm=mm, tz=tz, julday=julday, rc=localrc)
    if (chem_rc_check(localrc, msg="Failed to retrieve model clock on local DE", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- GOCART time steps start from 1, while model time steps start from 0
    advanceCount = advanceCount + 1

    do de = 0, deCount-1
      call chem_model_get(de=de, config=config, data=data, &
        stateIn=stateIn, stateOut=stateOut, tile=tile, rc=localrc)
      if (chem_rc_check(localrc, msg="Failed to retrieve model on local DE", &
        file=__FILE__, line=__LINE__, rc=rc)) return

      call chem_model_domain_get(de=de, ids=is, ide=ie, jds=js, jde=je, ni=ni, nl=nl, &
        lon=lon, lat=lat, rc=localrc)
      if (chem_rc_check(localrc, msg="Failed to retrieve model domain on local DE", &
        file=__FILE__, line=__LINE__, rc=rc)) return

      call gocart_advance(config % readrestart, config % chem_opt,&
        config % chem_in_opt, config % chem_conv_tr,config % biomass_burn_opt, &
        config % seas_opt, config % dust_opt, config % dmsemis_opt, &
        config % aer_ra_feedback, &
        config % call_chemistry, config % call_radiation, &
        config % plumerise_flag, config % plumerisefire_frq, &
        config % kemit, &
        advanceCount, dts, mm, tz, julday, &
        ! -- background data 
        data % p_gocart, &
        data % clayfrac, &
        data % dm0,      &               ! dms reference emissions
        data % emiss_ab, &
        data % emiss_abu, &
        data % emiss_ash_dt, &
        data % emiss_ash_height, &
        data % emiss_ash_mass,   &
        data % emiss_tr_dt, &
        data % emiss_tr_height, &
        data % emiss_tr_mass,   &
        data % ero1, &
        data % ero2, &
        data % ero3, &
        data % ssm,  &
        data % h2o2_backgd, &
        data % no3_backgd, &
        data % oh_backgd, &
        data % plumefrp, &
        data % plumestuff, &
        data % sandfrac, &
        data % th_pvsrf, &
        ! -- imported atmospheric fields
        stateIn % area, &
        stateIn % hf2d, &
        stateIn % pb2d, &
        stateIn % rc2d, &
        stateIn % rn2d, &
        stateIn % rsds, &
        stateIn % slmsk2d, &
        stateIn % snwdph2d, &
        stateIn % stype2d, &
        stateIn % ts2d, &
        stateIn % us2d, &
        stateIn % vtype2d, &
        stateIn % vfrac2d, &
        stateIn % zorl2d, &
        stateIn % exch, &
        stateIn % ph3d, &
        stateIn % phl3d, &
        stateIn % pr3d, &
        stateIn % prl3d, &
        stateIn % sm3d, &
        stateIn % tk3d, &
        stateIn % us3d, &
        stateIn % vs3d, &
        stateIn % ws3d, &
        stateIn % tr3d, &
        ! -- output tracers and tracer diagnostics
        stateOut % tr3d, &
        stateOut % trcm, &
        stateOut % trab, &
        stateOut % truf, &
        stateOut % trdf, &
        data % trdp, &
        data % ext_cof, &
        data % sscal, &
        data % asymp, &
        data % aod2d, &
        data % pm10, &
        data % pm25, &
        data % ebu_oc, &
        data % oh_bg, &
        data % h2o2_bg, &
        data % no3_bg, &
        data % wet_dep, &
        ! -- buffers
        data % rainl, &
        data % rainc, &
        data % eburn, &
        ! -- array size
        nl, ni, &
        config % ntra, config % ntrb, config % nvl_gocart, config % nbands, &
        config % numgas, config % num_ebu, &
        config % num_ebu_in, config % num_soil_layers, config % num_chem, config % num_moist, &
        config % num_emis_vol, config % num_emis_ant, &
        config % num_emis_dust, config % num_emis_seas, &
        config % num_asym_par, config % num_bscat_coef, &
        config % num_ext_coef, &
        ! -- domain
        lon, lat, &
        is, ie, js, je, 1, nl, &
        is, ie, js, je, 1, ni, &
        tile, &
        verbose=chem_comm_isroot(), rc=localrc)

        if (chem_rc_check(localrc, msg="Failure advancing GOCART", &
          file=__FILE__, line=__LINE__, rc=rc)) return
    end do


  end subroutine gocart_model_advance

end module gocart_model_mod
