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
    integer :: localrc
    integer :: de, deCount
    integer :: advanceCount, julday, mm, tz
    integer :: is, ie, js, je, ni, nl
    real(CHEM_KIND_R8) :: dts
    real(CHEM_KIND_R8), dimension(:,:), pointer :: lat, lon
    type(chem_config_type), pointer :: config
    type(chem_core_type),   pointer :: core
    type(chem_data_type),   pointer :: data
    type(chem_state_type),  pointer :: stateIn

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    call chem_model_get(deCount=deCount, rc=localrc)
    if (chem_rc_check(localrc, msg="Failed to retrieve model", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    if (deCount < 1) return

    call chem_model_clock_get(advanceCount=advanceCount, dts=dts, mm=mm, tz=tz, julday=julday, rc=localrc)
    if (chem_rc_check(localrc, msg="Failed to retrieve model clock on local DE", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    do de = 0, deCount-1
      call chem_model_get(de=de, config=config, core=core, data=data, stateIn=stateIn, rc=localrc)
      if (chem_rc_check(localrc, msg="Failed to retrieve model on local DE", &
        file=__FILE__, line=__LINE__, rc=rc)) return

      call chem_model_domain_get(de=de, ids=is, ide=ie, jds=js, jde=je, ni=ni, nl=nl, &
        lon=lon, lat=lat, rc=localrc)
      if (chem_rc_check(localrc, msg="Failed to retrieve model domain on local DE", &
        file=__FILE__, line=__LINE__, rc=rc)) return

      print *,'gocart_model_advance(rc) -- before gocart_advance(), de=',de, deCount
      call gocart_advance(config % readrestart, &
        config % chem_opt, config % chem_in_opt, config % biomass_burn_opt, &
        config % seas_opt, config % dust_opt, config % dmsemis_opt, &
        config % aer_ra_feedback, &
        config % call_biomass, config % call_chemistry, config % call_radiation, &
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
        data % h2o2_backgd, &
        data % no3_backgd, &
        data % oh_backgd, &
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
        ! -- output tracers
        data % tr3d, &
        data % trdp, &
        data % emi_d1, &
        data % emi_d2, &
        data % emi_d3, &
        data % emi_d4, &
        data % emi_d5, &
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
!       real(lon, CHEM_KIND_R4), real(lat, CHEM_KIND_R4), &
        lon, lat, &
        is, ie, js, je, 1, nl, &
        is, ie, js, je, 1, ni, &
        rc=localrc)

        if (chem_rc_check(localrc, msg="Failure advancing GOCART", &
          file=__FILE__, line=__LINE__, rc=rc)) return
    end do


  end subroutine gocart_model_advance

#if 0

  subroutine gocart_model_run(de, rc)
    integer, optional, intent(in)  :: de
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    ! -- pre-step
    call gocart_model_prep(de=de, rc=localrc)
    if (chem_rc_check(localrc, msg="Failure in prep phase to GOCART run", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- step

    ! -- post step


  end subroutine gocart_model_run

! -- private methods

  subroutine gocart_model_prep(de, rc)
    integer, optional, intent(in)  :: de
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: is, ie, js, je, ks, ke
    integer :: i, j, k, kk, kk4, kl, ko, k_initial, k_final, l, ll, nv
    real(CHEM_KIND_R4) :: aln, ashz_above_vent, eh, pl, pu, pwant, threshold
    real(CHEM_KIND_R4) :: factor, factor2
    real(CHEM_KIND_R4) :: so2_mass, thv, x1
    real(CHEM_KIND_R4), dimension(:), allocatable :: vert_mass_dist

    type(chem_core_type),    pointer :: core
    type(chem_config_type),  pointer :: config
    type(chem_clock_type),   pointer :: clock
    type(chem_data_type),    pointer :: data
    type(chem_state_type),   pointer :: state
    type(chem_domain_type),  pointer :: domain
    type(chem_species_type), pointer :: s

    ! -- volcano ashes parameters
    !  + original
    ! real, dimension(6) :: h = (/ 9., 16., 58., 79., 109., 129., huge(1.0) /)
    !  + real-time default (if volcano starts at h = 0)
    real(CHEM_KIND_R4), parameter :: h(7)               = (/ (240., i = 1, 6), huge(1.0) /)
    real(CHEM_KIND_R4), parameter :: emiss_ash_table(6) = (/  5834.,  3834.,  5834.,  3334.,  3334.,  2334. /)
    real(CHEM_KIND_R4), parameter :: eh_table(6)        = (/ 3.11e5, 3.87e4, 3.11e5, 2.17e4, 2.17e4, 4.93e3 /)
    real(CHEM_KIND_R4), parameter :: percen_mass_umbrel = 0.75_CHEM_KIND_R4
    real(CHEM_KIND_R4), parameter :: base_umbrel        = 0.25_CHEM_KIND_R4    ! fraction
    real(CHEM_KIND_R4), parameter :: base_umbrel2       = 1.00_CHEM_KIND_R4    ! evenly distribution
    ! .. Intrinsic Functions ..
    INTRINSIC max, min, float

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    ! -- get model on input DE
    call chem_model_get(de=de, clock=clock, config=config, core=core, data=data, &
      domain=domain, stateIn=state, rc=localrc)
    if (chem_rc_check(localrc, msg="Unable to retrieve model info.", &
      file=__FILE__, line=__LINE__, rc=rc)
    return

    ! -- set pointer to container of species pointers
    s => config % species

    ! -- get domain bounds
    call chem_domain_get(domain, ids=is, ide=ie, jds=js, jde=je)

    ! -- compute heights and thickness of geopotential layers
    do j = js, je
      do i = is, ie
        core % z_at_w(i,ks,j) = max(0._CHEM_KIND_R4, state % ph3d(i,j,ks)/g)
      end do
    end do

    do j = js, je
      do k = ks, ke
        do i = is, ie
          core % dz8w(i,k,j) = abs((state % ph3d(i,j,k+1)-state % ph3d(i,j,k))/g)
          core % z_at_w(i,k+1,j) = core % z_at_w(i,k,j) + core % dz8w(i,k,j)
        end do
      end do
    end do

    ! -- pressure at interface
    do j = js, je
      do k = ks, ke+1
        do i = is, ie
          core % p8w(i,k,j) = state % pr3d(i,j,k)
        end do
      end do
    end do

    ! ... TODO: 2D variables - BEGIN
    do j = js, je
      do i = is, ie
!       raincv_b(i,j)=rcav(j)
!       pbl(i,j)=pb2d(j)
!       dms_0(i,j)=dm0(j)
!       hfx(i,j)=hf2d(j)
!       snowh(i,j)=snwdph2d(j)*.001
        core % erod(i,j,1) = data % ero1(i,j)
        core % erod(i,j,2) = data % ero2(i,j)
        core % erod(i,j,3) = data % ero3(i,j)
!       xlat(i,j)=deg_lat(j)
!       xlong(i,j)=deg_lon(j)
!       ust(i,j)=us2d(j)
!       tsk(i,j)=ts2d(j)
!       ivgtyp(i,j)=VTYPE2d(j)
!       isltyp(i,j)=STYPE2d(j)
!       gsw(i,j)=rsds(j)
!       vegfra(i,j)=VFRAC2d(j)
!       rmol(i,j)=0.
!       znt(i,j)=zorl2d(j)*.01
!SLMSK   - SEA(0),LAND(1),ICE(2) MASK
!       xland(i,j)=1.
!       if (slmsk2d(j) == 0.) then
!         xland(i,j) = 0.
!       else if (slmsk2d(j) == 1.) then
!         xland(i,j) = 1.
!       else if (slmsk2d(j) == 2.) then
!         xland(i,j) = 2.
!       end if
!       dxy(i,j)=area(j)
        core % u10(i,j) = state % us3d(i,j,1)
        core % v10(i,j) = state % vs3d(i,j,1)
!       clayf(i,j) = clayfrac(j)
!       sandf(i,j) = sandfrac(j)
      end do
    end do
    ! ... TODO: 2D variables - END

    ! -- emissions
    if ((config % chem_opt == CHEM_OPT_RACM_SOA_VBS) .or. (s % p_bc2 > 1)) then

      k = 1
      do j = js, je
        do i = is, ie
          core % emis_ant(i,k,j,s % p_e_bc)    = data % emiss_ab(i,j,s % p_e_bc)
          core % emis_ant(i,k,j,s % p_e_oc)    = data % emiss_ab(i,j,s % p_e_oc)
          core % emis_ant(i,k,j,s % p_e_sulf)  = data % emiss_ab(i,j,s % p_e_sulf)
          core % emis_ant(i,k,j,s % p_e_so2)   = data % emiss_ab(i,j,s % p_e_so2)
          core % emis_ant(i,k,j,s % p_e_dms)   = 0._CHEM_KIND_R4 !emiss_ab(j,p_e_dms)
          core % emis_ant(i,k,j,s % p_e_pm_25) = data % emiss_ab(i,j,s % p_e_pm_25)
          core % emis_ant(i,k,j,s % p_e_pm_10) = data % emiss_ab(i,j,s % p_e_pm_10)

          ! -- gas anth emission for chem_opt 301
          if ((config % chem_opt == CHEM_OPT_GOCART_RACM) .or. &
              (config % chem_opt == CHEM_OPT_RACM_SOA_VBS)) then
            core % emis_ant(i,k,j,s % p_e_iso)  = data % emiss_ab(i,j,s % p_e_iso)
            core % emis_ant(i,k,j,s % p_e_no)   = data % emiss_ab(i,j,s % p_e_no)
            core % emis_ant(i,k,j,s % p_e_no2)  = data % emiss_ab(i,j,s % p_e_no2)
            core % emis_ant(i,k,j,s % p_e_co)   = data % emiss_ab(i,j,s % p_e_co)
            core % emis_ant(i,k,j,s % p_e_eth)  = data % emiss_ab(i,j,s % p_e_eth)
            core % emis_ant(i,k,j,s % p_e_hc3)  = data % emiss_ab(i,j,s % p_e_hc3)
            core % emis_ant(i,k,j,s % p_e_hc5)  = data % emiss_ab(i,j,s % p_e_hc5)
            core % emis_ant(i,k,j,s % p_e_hc8)  = data % emiss_ab(i,j,s % p_e_hc8)
            core % emis_ant(i,k,j,s % p_e_xyl)  = data % emiss_ab(i,j,s % p_e_xyl)
            core % emis_ant(i,k,j,s % p_e_olt)  = data % emiss_ab(i,j,s % p_e_olt)
            core % emis_ant(i,k,j,s % p_e_oli)  = data % emiss_ab(i,j,s % p_e_oli)
            core % emis_ant(i,k,j,s % p_e_tol)  = data % emiss_ab(i,j,s % p_e_tol)
            core % emis_ant(i,k,j,s % p_e_csl)  = data % emiss_ab(i,j,s % p_e_csl)
            core % emis_ant(i,k,j,s % p_e_hcho) = data % emiss_ab(i,j,s % p_e_hcho)
            core % emis_ant(i,k,j,s % p_e_ald)  = data % emiss_ab(i,j,s % p_e_ald)
            core % emis_ant(i,k,j,s % p_e_ket)  = data % emiss_ab(i,j,s % p_e_ket)
            core % emis_ant(i,k,j,s % p_e_ora2) = data % emiss_ab(i,j,s % p_e_ora2)
            core % emis_ant(i,k,j,s % p_e_nh3)  = data % emiss_ab(i,j,s % p_e_nh3)
          endif

          ! -- gas biomass burning emission for chem_opt 301
          if (config % biomass_burn_opt > 0) then

            core % ebu_in(i,j,s % p_ebu_in_oc)   = data % emiss_abu(i,j,s % p_e_oc)
            core % ebu_in(i,j,s % p_ebu_in_bc)   = data % emiss_abu(i,j,s % p_e_bc)
            core % ebu_in(i,j,s % p_ebu_in_pm25) = data % emiss_abu(i,j,s % p_e_pm_25)
            core % ebu_in(i,j,s % p_ebu_in_pm10) = data % emiss_abu(i,j,s % p_e_pm_10)
            core % ebu_in(i,j,s % p_ebu_in_so2)  = data % emiss_abu(i,j,s % p_e_so2)
            core % ebu_in(i,j,s % p_ebu_in_dms)  = 0._CHEM_KIND_R4 !emiss_abu(j,p_e_dms)

            core % mean_fct_agtf(i,j) = data % plumestuff(i,j,1)
            core % mean_fct_agef(i,j) = data % plumestuff(i,j,2)
            core % mean_fct_agsv(i,j) = data % plumestuff(i,j,3)
            core % mean_fct_aggr(i,j) = data % plumestuff(i,j,4)
            core % firesize_agtf(i,j) = data % plumestuff(i,j,5)
            core % firesize_agef(i,j) = data % plumestuff(i,j,6)
            core % firesize_agsv(i,j) = data % plumestuff(i,j,7)
            core % firesize_aggr(i,j) = data % plumestuff(i,j,8)

            if ((config % chem_opt == CHEM_OPT_GOCART_RACM) .or. &
                (config % chem_opt == CHEM_OPT_RACM_SOA_VBS)) then
              core % ebu_in(i,j,s % p_ebu_in_iso) = data % emiss_abu(i,j,s % p_e_iso)
              core % ebu_in(i,j,s % p_ebu_in_no)  = data % emiss_abu(i,j,s % p_e_no)
              core % ebu_in(i,j,s % p_ebu_in_no2) = data % emiss_abu(i,j,s % p_e_no2)
              core % ebu_in(i,j,s % p_ebu_in_co)  = data % emiss_abu(i,j,s % p_e_co)
              core % ebu_in(i,j,s % p_ebu_in_eth) = data % emiss_abu(i,j,s % p_e_eth)
              core % ebu_in(i,j,s % p_ebu_in_hc3) = data % emiss_abu(i,j,s % p_e_hc3)
              core % ebu_in(i,j,s % p_ebu_in_hc5) = data % emiss_abu(i,j,s % p_e_hc5)
              core % ebu_in(i,j,s % p_ebu_in_hc8) = data % emiss_abu(i,j,s % p_e_hc8)
              core % ebu_in(i,j,s % p_ebu_in_xyl) = data % emiss_abu(i,j,s % p_e_xyl)
              core % ebu_in(i,j,s % p_ebu_in_olt) = data % emiss_abu(i,j,s % p_e_olt)
              core % ebu_in(i,j,s % p_ebu_in_oli) = data % emiss_abu(i,j,s % p_e_oli)
              core % ebu_in(i,j,s % p_ebu_in_tol) = data % emiss_abu(i,j,s % p_e_tol)
              core % ebu_in(i,j,s % p_ebu_in_csl) = data % emiss_abu(i,j,s % p_e_csl)
              core % ebu_in(i,j,s % p_ebu_in_hcho)= data % emiss_abu(i,j,s % p_e_hcho)
              core % ebu_in(i,j,s % p_ebu_in_ald) = data % emiss_abu(i,j,s % p_e_ald)
              core % ebu_in(i,j,s % p_ebu_in_ket) = data % emiss_abu(i,j,s % p_e_ket)
              core % ebu_in(i,j,s % p_ebu_in_ora2)= data % emiss_abu(i,j,s % p_e_ora2)
              core % ebu_in(i,j,s % p_ebu_in_nh3) = data % emiss_abu(i,j,s % p_e_nh3)
            endif
          endif

        end do
      end do

    else if ( s % p_tr2 > 1 ) then  ! tracer options

      ! -- tracer run
      k = ks
      do j = js, je
        do i = is, ie
          core % emis_ant(i,k,j,s % p_e_tr1) = data % emiss_ab(i,j,s % p_e_tr1)
          core % emis_ant(i,k,j,s % p_e_tr2) = data % emiss_ab(i,j,s % p_e_tr2)
        enddo
      enddo

    else if (( s % p_tr2 > 1 ) .and. ( s % p_bc2 > 1 )) then

      call chem_rc_set(CHEM_RC_FAILURE, msg="Inconsistent options detected.", &
        file=__FILE__, line=__LINE__, rc=rc)
      return

    end if

    ! -- convert atmospheric input
    do k = ks, ke+1
      kk = min(k, ke)
      do j = js, je
        do i = is, ie
          core % zmid(i,k,j)    = .5_CHEM_KIND_R4 * (state % ph3d(i,j,kk+1) + state % ph3d(i,j,kk))/g
          core % dz8w(i,k,j)    = core % z_at_w(i,kk+1,j) - core % z_at_w(i,kk,j)
          core % p_phy(i,k,j)   = .5_CHEM_KIND_R4 * (core % p8w(i,kk,j) + core % p8w(i,kk+1,j))
          thv                   = state % tr3d(i,j,kk,1)/(1.+0.6078*state % tr3d(i,j,kk,2))
          core % t_phy(i,k,j)   = thv * (core % p_phy(i,k,j)/p1000)**(rd/cp)
          core % u_phy(i,k,j)   = state % us3d(i,j,kk)
          core % v_phy(i,k,j)   = state % vs3d(i,j,kk)
          core % exch_h(i,k,j)  = state % exch(i,j,kk)
          core % rho_phy(i,k,j) = core % p_phy(i,k,j)/(RD*core % t_phy(i,k,j)*(1.+.608*state % tr3d(i,j,kk,2)))
          core % rri(i,k,j)     = 1._CHEM_KIND_R4 / core % rho_phy(i,k,j)
          core % vvel(i,k,j)    = -state % ws3d(i,j,kk) * core % rri(i,k,j) / g
          core % convfac(i,k,j) = core % p_phy(i,k,j) / rgasuniv / core % t_phy(i,k,j)
          core % moist(i,k,j,:) = 0._CHEM_KIND_R4
          core % moist(i,k,j,1) = state % tr3d(i,j,kk,2)
          if (core % t_phy(i,k,j) > 265._CHEM_KIND_R4) then
            core % moist(i,k,j,2) = state % tr3d(i,j,kk,3)
            core % moist(i,k,j,3) = 0._CHEM_KIND_R4
            if (core % moist(i,k,j,2) < 1.e-8_CHEM_KIND_R4) core % moist(i,k,j,2) = 0._CHEM_KIND_R4
          else
            core % moist(i,k,j,2) = 0._CHEM_KIND_R4
            core % moist(i,k,j,3) = state % tr3d(i,j,kk,3)
            if (core % moist(i,k,j,3) < 1.e-8_CHEM_KIND_R4) core % moist(i,k,j,3) = 0._CHEM_KIND_R4
          endif
          core % relhum(i,k,j) = MIN( .95_CHEM_KIND_R4, core % moist(i,k,j,1) / &
            (3.80_CHEM_KIND_R4*exp(17.27_CHEM_KIND_R4*(core % t_phy(i,k,j)-273._CHEM_KIND_R4)/ &
            (core % t_phy(i,k,j)-36._CHEM_KIND_R4))/(.01_CHEM_KIND_R4*core % p_phy(i,k,j))))
          core % relhum(i,k,j) = max(0.1_CHEM_KIND_R4,core % relhum(i,k,j))
        end do
      end do
    end do
    
    do j = js, je
      do k = 2, ke
        do i = is, ie
          core % t8w(i,k,j) = .5_CHEM_KIND_R4 * (core % t_phy(i,k,j) + core % t_phy(i,k-1,j)) ! .5*(tk3d(k-1,j)+tk3d(k,j))
        end do
      end do
    end do

    ! -- only used in phtolysis....
    do j = js, je
      do i = is, ie
        core % t8w(i,1,j)    = core % t_phy(i,1,j)
        core % t8w(i,ke+1,j) = core % t_phy(i,ke,j)
      end do
    end do

    do nv = 1, config % num_chem
      do j = js, je
        do k = ks, ke+1
          kk = min(k, ke)
          do i = is, ie
            core % chem(i,k,j,nv) = state % tr3d(i,j,kk,config % ntra + nv)
          end do
        end do
      end do
    end do

    ! -- run mode

    if ( .not.config % readrestart) then
      if (clock % advanceCount <= 1) then
        if ((config % chem_opt == CHEM_OPT_GOCART_RACM) .or. &
            (config % chem_opt == CHEM_OPT_RACM_SOA_VBS)) then
          do k = ks, ke + 1
            kk = min(k, ke)
            do j = js, je
              do i = is, ie
                threshold = min(400._CHEM_KIND_R4, core % th_pvsrf(i,j))
                if (core % chem(i,k,j,1) >= threshold) then
                  !convert kg/kg to ppm
                  core % chem(i,j,k,s % p_o3) = 1.e+06_CHEM_KIND_R4 * (airmw/48._CHEM_KIND_R4) * state % tr3d(i,j,kk,4)
                end if
              end do
            end do
          end do
        end if
      end if
    else
      ! -- if restart, do nothing for now
    end if

    ! -- if GOCART is called, prepare background fields
    if ((config % chem_opt == CHEM_OPT_GOCART) .and. config % call_gocart) then

      do j = js, je
        do k = ks, ke
          do i = is, ie
            do ll = 2, config % nvl_gocart
              l = ll
              if (data % p_gocart(l) < .01_CHEM_KIND_R4 * core % p_phy(i,k,j)) exit
            end do
            pu=alog(data % p_gocart(l))
            pl=alog(data % p_gocart(l-1))
            pwant=alog(.01_CHEM_KIND_R4*core % p_phy(i,k,j))
            if (pwant > pl)then
              core % backg_oh(i,k,j)   = data % oh_backgd(i,j,l)
              core % backg_h2o2(i,k,j) = data % h2o2_backgd(i,j,l)
              core % backg_no3(i,k,j)  = data % no3_backgd(i,j,l)
            else
              aln=(data % oh_backgd(i,j,l)*(pwant-pl)+            &
                data % oh_backgd(i,j,l-1)*(pu-pwant))/(pu-pl)
              core % backg_oh(i,k,j)=aln
              aln=(data % h2o2_backgd(i,j,l)*(pwant-pl)+            &
                data % h2o2_backgd(i,j,l-1)*(pu-pwant))/(pu-pl)
              core % backg_h2o2(i,k,j)=aln
              aln=(data % no3_backgd(i,j,l)*(pwant-pl)+            &
                data % no3_backgd(i,j,l-1)*(pu-pwant))/(pu-pl)
              core % backg_no3(i,k,j)=aln
            endif
          end do
        end do
      end do

    end if

    ! -- include emissions into chem array

    if (s % p_bc2 > 1)then

      if (config % chem_opt == CHEM_OPT_GOCART) then
        do j = js, je
          do i = is, ie
            factor  = clock % dts * core % rri(i,k,j) / core % dz8w(i,k,j)
            factor2 = 4.828e-4_CHEM_KIND_R4 * clock % dts * core % rri(i,k,j) / &
              (60._CHEM_KIND_R4 * core % dz8w(i,k,j))
            core % chem(i,k,j,s % p_bc1)  = core % chem(i,k,j,s % p_bc1)  + core % emis_ant(i,k,j,s % p_e_bc)    * factor
            core % chem(i,k,j,s % p_oc1)  = core % chem(i,k,j,s % p_oc1)  + core % emis_ant(i,k,j,s % p_e_oc)    * factor
            core % chem(i,k,j,s % p_p25)  = core % chem(i,k,j,s % p_p25)  + core % emis_ant(i,k,j,s % p_e_pm_25) * factor
            core % chem(i,k,j,s % p_p10)  = core % chem(i,k,j,s % p_p10)  + core % emis_ant(i,k,j,s % p_e_pm_10) * factor
            core % chem(i,k,j,s % p_sulf) = core % chem(i,k,j,s % p_sulf) + core % emis_ant(i,k,j,s % p_e_sulf)  * factor
            core % chem(i,k,j,s % p_so2)  = core % chem(i,k,j,s % p_so2)  + core % emis_ant(i,k,j,s % p_e_so2)   * factor2
          end do
        end do
      end if

    else if (s % p_tr2 > 1) then    !co2 here

      do j = js, je
        do i = is, ie
          factor2 = 4.828e-4_CHEM_KIND_R4 * clock % dts * core % rri(i,k,j) / &
            (60._CHEM_KIND_R4 * core % dz8w(i,k,j))
          core % chem(i,k,j,s % p_tr1) = core % chem(i,k,j,s % p_tr1) + core % emis_ant(i,k,j,s % p_e_tr1) * factor2
          core % chem(i,k,j,s % p_tr2) = core % chem(i,k,j,s % p_tr2) + core % emis_ant(i,k,j,s % p_e_tr2) * factor2
        end do
      end do

    else if ((s % p_tr2 > 1) .and. (s % p_bc2 > 1))then

      call chem_rc_set(CHEM_RC_FAILURE, msg="Inconsistent options detected.", &
        file=__FILE__, line=__LINE__, rc=rc)
      return

    endif

    ! -- soil moisture
    do nv = 1, config % num_soil_layers
      do j = js, je
        do i = is, ie
          core % smois(i,nv,j) = state % sm3d(i,j,nv)
        end do
      end do
    end do

    ! -- volcanic emissions
!   allocate(so2_mass(is:je, js:je), stat=localrc)
!   if (chem_rc_test((localrc /= 0), msg="Unable to allocate work space.", &
!     file=__FILE__, line=__LINE__, rc=rc)) return

   ! -- overwrite values?
    do j = js, je
      do i = is, ie
        if (data % emiss_ash_dt(i,j) > 0._CHEM_KIND_R4) then
!         so2_mass(i,j)=1.5e4*3600.*1.e9/64./state % area(i,j)
          eh=2600.*(data % emiss_ash_height(i,j)*.0005)**4.1494
          data % emiss_ash_mass(i,j)=eh*1.e9/domain % area(i,j)
        end if
      end do
    end do

    ! -- hardcode for special retro case (set h1 - h6 properly
    do nv = 1, 6
      if ((clock % h  >= h(nv)) .and. (clock % h < h(nv+1))) then
        do j = js, je
          if (data % emiss_ash_dt(i,j) > 0) then
            data % emiss_ash_height(i,j) = emiss_ash_table(nv)
            data % emiss_ash_mass(i,j)   = 1.e+09 * eh_table(nv) / domain % area(i,j)
          end if
        end do
        exit
      end if
    end do

    print *,'chem_prep: stage 2'
    ! -- real-time application, keeping eruption constant
!
    if (clock % advanceCount <= 2) then

      ! -- volcanic emissions
      core % emis_vol(:,:,:,:) = 0._CHEM_KIND_R4

      allocate(vert_mass_dist(ks:ke), stat=localrc)
      if (chem_rc_test((localrc /= 0), msg="Unable to allocate work space.", &
        file=__FILE__, line=__LINE__, rc=rc)) return
      vert_mass_dist = 0._CHEM_KIND_R4

      do j = js, je
        do i = is, ie
          if (data % emiss_ash_dt(i,j)     <= 0) cycle
          if (data % emiss_ash_height(i,j) <= 0) cycle
          ashz_above_vent = data % emiss_ash_height(i,j) + core % z_at_w(i,ks,j)
          do k = ke-1, ks, -1
            if (core % z_at_w(i,k,j) < ashz_above_vent) then
              k_final = k+1
              exit
            end if !inner
          end do
          do k = ke-1, ks, -1
            if (core % z_at_w(i,k,j) < (1.-base_umbrel)*ashz_above_vent) then
              k_initial = k
              exit
            endif !inner
          enddo
          vert_mass_dist=0._CHEM_KIND_R4
             
          ! -- parabolic vertical distribution between k_initial and k_final
          kk4 = k_final-k_initial+2
          do ko=1,kk4-1
            kl=ko+k_initial-1
            vert_mass_dist(kl) = 6.*percen_mass_umbrel* float(ko)/float(kk4)**2 * (1. - float(ko)/float(kk4))
          enddo
          if (sum(vert_mass_dist(kts:kte)) /= percen_mass_umbrel) then
            x1= ( percen_mass_umbrel- sum(vert_mass_dist(kts:kte)) )/float(k_final-k_initial+1)
            do ko=k_initial,k_final
              vert_mass_dist(ko) = vert_mass_dist(ko)+ x1 !- values between 0 and 1.
            enddo
                   !pause
          endif !inner
          !k_final > 0 .and. k_initial >
    
          ! -- linear detrainment from vent to base of umbrella
          do ko=1,k_initial-1
            vert_mass_dist(ko)=float(ko)/float(k_initial-1)
          enddo
          x1=sum(vert_mass_dist(1:k_initial-1))
    
          do ko=1,k_initial-1
            vert_mass_dist(ko)=(1.-percen_mass_umbrel)*vert_mass_dist(ko)/x1
          enddo

          if (config % chem_opt == CHEM_OPT_GOCART) then

            ! -- if applied to gocart we only need finest ash bins, we use the coarse one for so2
            so2_mass = 1.5e4*3600.*1.e9/64./state % area(i,j)
            do ko = 1, k_final
              core % emis_vol(i,ko,j,s % p_e_vash1) =        vert_mass_dist(ko) * so2_mass
              core % emis_vol(i,ko,j,s % p_e_vash2) = .08  * vert_mass_dist(ko) * data % emiss_ash_mass(i,j)
              core % emis_vol(i,ko,j,s % p_e_vash3) = .05  * vert_mass_dist(ko) * data % emiss_ash_mass(i,j)
              core % emis_vol(i,ko,j,s % p_e_vash4) = .035 * vert_mass_dist(ko) * data % emiss_ash_mass(i,j)
            end do
          end if !chem_opt==316 or 317,300,502
 
          do ko = k_final+1, kte
            core % emis_vol(i,ko,j,s % p_e_vash1) = 0._CHEM_KIND_R4
            core % emis_vol(i,ko,j,s % p_e_vash2) = 0._CHEM_KIND_R4
            core % emis_vol(i,ko,j,s % p_e_vash3) = 0._CHEM_KIND_R4
            core % emis_vol(i,ko,j,s % p_e_vash4) = 0._CHEM_KIND_R4
          end do
        end do
      end do

      deallocate(vert_mass_dist, stat=localrc)
      if (chem_rc_test((localrc /= 0), msg="Unable to free up work space.", &
        file=__FILE__, line=__LINE__, rc=rc)) return

    end if ! curr_mins 
    print *,'chem_prep: stage 3'

    ! next is done to scale background oh and no3 in dependence on average zenith angle and day/night for no3
    ! this is done since background values are only available as average/month. It will not be necessary if other
    ! chemistry packages are used that provide oh,no3,h2o2

    if ((config % chem_opt == CHEM_OPT_RACM_SOA_VBS) .or. &
        (config % chem_opt >= CHEM_OPT_GOCART .and. config % chem_opt < CHEM_OPT_MAX)) then

#if 0
      ndystep=86400/ifix(dtstep)

      ! ------ TODO: to be moved to clock_prep()? - BEGIN
      do j = js, je
        do i = is, ie
          tcosz(i,j)=0.
          ttday(i,j)=0.
          rlat=xlat(i,j)*3.1415926535590/180.
          xlonn=xlong(i,j)
          do n=1,ndystep
            xtime=n*dtstep/60.
            ixhour=ifix(gmt+.01)+ifix(xtime/60.)
            xhour=float(ixhour)
            xmin=60.*gmt+(xtime-xhour*60.)
            gmtp=mod(xhour,24.)
            gmtp=gmtp+xmin/60.
            CALL szangle(1, 1, julday, gmtp, sza, cosszax,xlonn,rlat)
            TCOSZ(i,j)=TCOSZ(I,J)+cosszax(1,1)
            if (cosszax(1,1) > 0.) ttday(i,j)=ttday(i,j)+dtstep
          enddo
        enddo
      enddo
      ! ------ TODO: to be moved to clock_prep()? - END
#endif

    endif !chem_opt >= 300 .and. chem_opt <  500
         !endif ! ktau

    print *,'chem_prep: stage 4'

    if(config % chem_opt == CHEM_OPT_GOCART)  then
      ! -- for gocart only lump ash into p25 and p10
      do j = js, je
        if (data % emiss_ash_dt(i,j) <= 0._CHEM_KIND_R4) cycle
        do k = ks, ke-2
          do i = is, ie
            factor = 4.828e-4_CHEM_KIND_R4 * clock % dts * core % rri(i,k,j) / &
              (60._CHEM_KIND_R4 * core % dz8w(i,k,j))
            factor2 = clock % dts * core % rri(i,k,j) / core % dz8w(i,k,j)
            core % chem(i,k,j,s % p_p25) = core % chem(i,k,j,s % p_p25) + core % emis_vol(i,k,j,s % p_e_vash4) * factor2   
            core % chem(i,k,j,s % p_so2) = core % chem(i,k,j,s % p_so2) + core % emis_vol(i,k,j,s % p_e_vash1) * factor
            core % chem(i,k,j,s % p_p10) = core % chem(i,k,j,s % p_p10) + &
              (core % emis_vol(i,k,j,s % p_e_vash3) + .5_CHEM_KIND_R4 * core % emis_vol(i,k,j,s % p_e_vash2)) * factor2
          end do
        end do
      end do
    end if

  end subroutine gocart_model_prep

#endif

end module gocart_model_mod
