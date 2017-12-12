module chem_shr_mod

  use mpp_mod,        only : mpp_broadcast, mpp_pe, mpp_root_pe, mpp_error, stdlog, FATAL
  use mpp_mod, only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use mpp_domains_mod,only : domain2d, mpp_get_compute_domain
  use chem_comm_mod
  use chem_const_mod, only : airmw, RGASUNIV, epsilc, rd, g
  use chem_types_mod, only : CHEM_MAXSTR
  use chem_state_mod
  use chem_config_mod
  use chem_vars_mod
  use chem_domain_mod
  use chem_species_mod, only : chem_species_setup
  use chem_io_mod, maybe_write => chem_data_write
! use time_mod

  implicit none

  integer, parameter :: year_track_data = 0

  public
  

contains

  subroutine chem_config_read(config)

    type(chem_config_type), intent(inout) :: config

    ! -- local variables
    integer :: rc, localrc
    integer :: unit
    integer :: buffer(24)
    real    :: rbuffer(4)
    character(len=CHEM_MAXSTR) :: cbuffer(1)

    ! -- begin

    ! -- set defaults
    chem_opt          = 301 ! chem option, 0=off, 300=on
    chemdt            = 3
    kemit             = 1
    phot_opt          = 1
    photdt            = 60
    depo_fact         = 0.
    DUST_OPT          = 3
    DMSEMIS_OPT       = 1
    SEAS_OPT          = 1
    EMISS_OPT         = 5
    BIO_EMISS_OPT     = 0
    BIOMASS_BURN_OPT  = 1
    PLUMERISEFIRE_FRQ = 60
    EMISS_INPT_OPT    = 1
    GAS_BC_OPT        = 1
    GAS_IC_OPT        = 1
    AER_BC_OPT        = 1
    AER_IC_OPT        = 1
    gaschem_onoff     = 1
    aerchem_onoff     = 1
    wetscav_onoff     = 0
    cldchem_onoff     = 0
    vertmix_onoff     = 1
    chem_conv_tr      = 1
    aer_ra_feedback   = 0
    chem_in_opt       = 0
    archive_step      = 1
    emi_inname         = ""
    fireemi_inname     = ""
    emi_outname        = ""
    fireemi_outname    = ""
    input_chem_inname  = ""
    input_chem_outname = ""
    chem_hist_outname  = "fim_out_"

    ! -- read chem configuration namelist

    if (chem_comm_isroot()) then
      unit = 200
      open(unit, file=chem_file_nml, form='formatted', status='old', iostat=rc)
    end if
    call chem_comm_bcast(rc, chem_comm_rootpe, localrc)
    if (rc /= RC_COMM_SUCCESS) call chem_comm_abort(msg='CHEM_CONFIG_READ: error opening input file: '//chem_file_nml)
    if (chem_comm_isroot()) then
      rewind(unit)
      read( unit, nml=chem_nml, iostat=rc )
      close(unit)
!     write(stdlog(), nml=chem_nml)
      write(6, nml=chem_nml)
    end if
    call chem_comm_bcast(rc, chem_comm_rootpe, localrc)
    if (rc /= RC_COMM_SUCCESS) call chem_comm_abort(msg='CHEM_CONFIG_READ: error reading &chem namelist from '//chem_file_nml)

    buffer = (/ &
      chem_opt,          &
      kemit,             &
      phot_opt,          &
      DUST_OPT,          &
      DMSEMIS_OPT,       &
      SEAS_OPT,          &
      EMISS_OPT,         & 
      BIO_EMISS_OPT,     &
      BIOMASS_BURN_OPT,  &
      PLUMERISEFIRE_FRQ, &
      EMISS_INPT_OPT,    &
      GAS_BC_OPT,        &
      GAS_IC_OPT,        &
      AER_BC_OPT,        &
      AER_IC_OPT,        &
      gaschem_onoff,     &
      aerchem_onoff,     &
      wetscav_onoff,     &
      cldchem_onoff,     &
      vertmix_onoff,     &
      chem_conv_tr,      &
      aer_ra_feedback,   &
      chem_in_opt,       &
      archive_step       &
      /)

    call chem_comm_bcast(buffer, chem_comm_rootpe, rc)
    if (rc /= RC_COMM_SUCCESS) call chem_comm_abort(msg='CHEM_CONFIG_READ: error broadcasting &chem namelist (1)')

    config % chem_opt          = buffer( 1 )
    config % kemit             = buffer( 2 )
    config % phot_opt          = buffer( 3 )
    config % dust_opt          = buffer( 4 )
    config % dmsemis_opt       = buffer( 5 )
    config % seas_opt          = buffer( 6 )
    config % emiss_opt         = buffer( 7 )
    config % bio_emiss_opt     = buffer( 8 )
    config % biomass_burn_opt  = buffer( 9 )
    config % plumerisefire_frq = buffer( 10 )
    config % emiss_inpt_opt    = buffer( 11 )
    config % gas_bc_opt        = buffer( 12 )
    config % gas_ic_opt        = buffer( 13 )
    config % aer_bc_opt        = buffer( 14 )
    config % aer_ic_opt        = buffer( 15 )
    config % gaschem_onoff     = buffer( 16 )
    config % aerchem_onoff     = buffer( 17 )
    config % wetscav_onoff     = buffer( 18 )
    config % cldchem_onoff     = buffer( 19 )
    config % vertmix_onoff     = buffer( 20 )
    config % chem_conv_tr      = buffer( 21 )
    config % aer_ra_feedback   = buffer( 22 )
    config % chem_in_opt       = buffer( 23 )
    config % archive_step      = buffer( 24 )

    rbuffer = (/ bioemdt, photdt, chemdt, depo_fact /)

    call chem_comm_bcast(rbuffer, chem_comm_rootpe, rc)
    if (rc /= RC_COMM_SUCCESS) call chem_comm_abort(msg='CHEM_CONFIG_READ: error broadcasting &chem namelist (2)')

    config % bioemdt = rbuffer(1)
    config % photdt  = rbuffer(2)
    config % chemdt  = rbuffer(3)
    config % depo_fact = rbuffer(4)

    cbuffer(1) = chem_hist_outname
    call chem_comm_bcast(cbuffer, chem_comm_rootpe, rc)
    if (rc /= RC_COMM_SUCCESS) call chem_comm_abort(msg='CHEM_CONFIG_READ: error broadcasting &chem namelist (3)')
    config % chem_hist_outname = cbuffer(1)

    print *,'chem: config read exit'
    
  end subroutine chem_config_read

  subroutine chem_config_set(nvl_in, nip_in, ntra_in, ntrb_in, glvl_in, yyyymmddhhmm_in)

    integer, intent(in) :: nvl_in, nip_in, ntra_in, ntrb_in, glvl_in
    character(len=*), intent(in) :: yyyymmddhhmm_in

    ! -- set chemistry spatial domain
    nip   = nip_in
    nvl   = nvl_in 
    nvlp1 = nvl + 1
    glvl  = glvl_in

    ! -- set tracers config
    ntra = ntra_in
    ntrb = ntrb_in

    ! -- set time stamp for output
    yyyymmddhhmm = yyyymmddhhmm_in
    
  end subroutine chem_config_set

  subroutine chem_control_setup(config)
    type(chem_config_type), intent(in) :: config

    integer, parameter :: mp_physics = 0
    integer :: chem_opt

    chem_opt = config % chem_opt

    ntrb = 0

    numgas=1
    num_ebu=0
    num_ebu_in=0
    if(mp_physics.eq.0 .and. chem_opt .eq. 317)then
       num_chem=17
       num_moist=3
       num_emis_ant=6
       num_emis_vol=4
       ntrb=ntrb+num_moist+num_chem-3           ! # of tracers + num_moist-3, num_chem - no ice variable transported
       num_ebu=30
       num_ebu_in=30
!    nvarp=nvarp+num_chem+num_moist - 3
    elseif (chem_opt .eq. 300 .and. mp_physics.eq.0 )then
       num_ebu=7
       num_ebu_in=7
       num_chem=19
       num_moist=3
       num_emis_ant=7
       num_emis_vol=0
       num_emis_season_bb=7
       num_emis_season_ant=7
       ntrb=ntrb+num_moist+num_chem-3           ! # of tracers + num_moist-3, num_chem - no ice variable transported
    elseif (chem_opt .eq. 301 .and. mp_physics.eq.0 )then
       num_chem=66
       numgas=49
       num_moist=3
       num_emis_ant=25
       num_ebu=25
       num_ebu_in=25
       num_emis_vol=0
       num_emis_season_bb=0
       num_emis_season_ant=0
       ntrb=ntrb+num_moist+num_chem-3           ! # of tracers + num_moist-3, num_chem - no ice variable transported
!    nvarp=nvarp+num_chem+num_moist - 2       ! nvarp includes all variables for
!    pressure level output
     elseif (chem_opt .eq. 108 .and. mp_physics.eq.0 )then
       num_chem=103
       numgas=65
       num_moist=3
       num_emis_ant=25
       num_ebu=25
       num_ebu_in=25
       num_emis_vol=0
       num_emis_season_bb=0
       num_emis_season_ant=0
       ntrb=ntrb+num_moist+num_chem-3           ! # of tracers + num_moist-3,num_chem - no ice variable transported
!    nvarp=nvarp+num_chem+num_moist - 2       ! nvarp includes all variables for
!    pressure level output
    elseif (chem_opt .eq. 0 .and. mp_physics.eq.28 )then
       num_chem=0
       num_moist=12
       num_emis_ant=0
       num_emis_vol=0
       ntrb=ntrb+num_moist+num_chem-2           ! # of tracers + num_moist-2, num_chem - transport of qi
    elseif (chem_opt .eq. 317 .and. mp_physics.eq.3 )then
       num_ebu=30
       num_ebu_in=30
       num_chem=17
       num_moist=5
       num_emis_ant=6
       num_emis_vol=4
       ntrb=ntrb+num_moist+num_chem-2           ! # of tracers + num_moist-2, num_chem - transport of qi
    elseif (chem_opt .eq. 500 .and. mp_physics.eq.0 )then
       num_chem=5
       num_moist=3
       num_emis_ant=5
       num_emis_vol=0
       ntrb=ntrb+num_moist+num_chem-2           ! # of tracers + num_moist-2, num_chem - transport of qi
    end if

    ! -- nbegin is the start address (-1) of the first chem variable in tr3d
    if (num_moist > 3) then
      nbegin = ntra + num_moist - 2
    else
      nbegin = ntra + num_moist - 3
    end if


  end subroutine chem_control_setup

  subroutine chem_setup(config)
    type(chem_config_type), intent(in) :: config
    call chem_control_setup(config)
    call chem_species_setup(config)
    call chem_domain_setup()
  end subroutine chem_setup

  subroutine chem_alloc_input

    print *,'chem_alloc_input: entering ...'
    allocate(rn2d    (jms:jme))
    allocate(rc2d    (jms:jme))
    allocate(ts2d    (jms:jme))
    allocate(us2d    (jms:jme))
    allocate(hf2d    (jms:jme))
    allocate(pb2d    (jms:jme))
    allocate(zorl2d  (jms:jme))
    allocate(vfrac2d (jms:jme))
    allocate(vtype2d (jms:jme))
    allocate(stype2d (jms:jme))
!   allocate(u10m    (jms:jme))
!   allocate(v10m    (jms:jme))
    allocate(snwdph2d(jms:jme))
    allocate(aod2d   (jms:jme))
    allocate(rsds    (jms:jme))
    allocate(area    (jms:jme))
    allocate(deg_lat (jms:jme))
    allocate(deg_lon (jms:jme))
    allocate(slmsk2d (jms:jme))

    allocate(sm3d  (4,jms:jme))

    allocate(tr3d(nvl  ,jms:jme,ntra+ntrb)) ! 1=pot.temp, 2=water vapor, 3=cloud water,4=ozone
    allocate(trdp(nvl  ,jms:jme,ntra+ntrb)) ! 1=pot.temp, 2=water vapor, 3=cloud water,4=ozone
    allocate(us3d(nvl  ,jms:jme))      ! zonal wind (m/s)
    allocate(vs3d(nvl  ,jms:jme))      ! meridional wind (m/s)
    allocate(ws3d(nvl  ,jms:jme))      ! vertical wind (Pa/s)
    allocate(tk3d(nvl  ,jms:jme))      ! temperature, kelvin
    allocate(exch(nvl  ,jms:jme))      ! exchange coeffs
    allocate(dp3d(nvl  ,jms:jme))      ! press.diff. between coord levels (Pa)
    allocate(pr3d(nvlp1,jms:jme))      ! pressure (pascal)
    allocate(ph3d(nvlp1,jms:jme))      ! geopotential (=gz), m^2/s^2

    allocate(gd_cloud ( 1, nvlp1, jms:jme ))
    allocate(gd_cldfr ( 1, nvlp1, jms:jme ))

    allocate(perm(nip))
    perm = 0

    rn2d     = 0.
    rc2d     = 0.
    ts2d     = 0.
    us2d     = 0.
    hf2d     = 0.
    pb2d     = 0.
    zorl2d   = 0.
    vfrac2d  = 0.
    vtype2d  = 0.
    stype2d  = 0.
    snwdph2d = 0.
    aod2d    = 0.
    rsds     = 0.
    area     = 0.
    deg_lat  = 0.
    deg_lon  = 0.
    slmsk2d  = 0.

    sm3d     = 0.

    tr3d = 0.
    trdp = 0.
    us3d = 0.
    vs3d = 0.
    ws3d = 0.
    pr3d = 0.
    ph3d = 0.
    tk3d = 0.
    exch = 0.
    dp3d = 0.

!   aerwrf    = 0.
    gd_cloud  = 0.
!   gd_cloud2 = 0.
    gd_cldfr  = 0.
    print *,'chem_alloc_input: done'

  end subroutine chem_alloc_input

  subroutine chem_alloc_workspace

    print *,'chem_alloc_workspace: entering ...'
    print *,'chem_alloc_workspace: step 0 ...'

    allocate(th_pvsrf(jms:jme))
    th_pvsrf = 0.
    ! -- chem_alloc
    allocate(pm25(nvl,jms:jme))
    pm25=0.
    allocate(p10(nvl,jms:jme))
    p10=0.
    allocate(tr1_tavg(nvl,jms:jme))
    tr1_tavg=0.
    allocate(d1st_ave(nvl,jms:jme))
    d1st_ave=0.
    allocate(d2st_ave(nvl,jms:jme))
    d2st_ave=0.
    allocate(d3st_ave(nvl,jms:jme))
    d3st_ave=0.
    allocate(d4st_ave(nvl,jms:jme))
    d4st_ave=0.
    allocate(d5st_ave(nvl,jms:jme))
    d5st_ave=0.
    allocate(ebu_oc(nvl,jms:jme))
    ebu_oc=0.
    allocate(oh_bg(nvl,jms:jme))
    oh_bg=0.
    allocate(h2o2_bg(nvl,jms:jme))
    h2o2_bg=0.
    allocate(no3_bg(nvl,jms:jme))
    no3_bg=0.
    allocate(oh_backgd(nvl_gocart,jms:jme))
    oh_backgd=0.
    allocate( h2o2_backgd(nvl_gocart,jms:jme) )
    h2o2_backgd = 0.
    allocate( no3_backgd(nvl_gocart,jms:jme) )
    no3_backgd = 0.
    allocate( rcav(jms:jme) )
    rcav = 0.
    allocate( rnav(jms:jme) )
    rnav = 0.
    allocate( ero1(jms:jme) )
    ero1 = 0.
    allocate( ero2(jms:jme) )
    ero2 = 0.
    allocate( ero3(jms:jme) )
    ero3 = 0.
    allocate( clayfrac(jms:jme) )
    clayfrac = 0.
    allocate( sandfrac(jms:jme) )
    sandfrac = 0.
    allocate( ashfall(jms:jme) )
    ashfall  = 0.
!  allocate( aod2d(jms:jme) )
!  aod2d = 0.
    allocate( wet_dep(jms:jme,num_chem) )
    wet_dep = 0.
    allocate( dry_dep(jms:jme,num_chem) )
    dry_dep = 0.
    allocate( plumestuff(jms:jme,8) )
    plumestuff = 0.
    allocate( emiss_ab(jms:jme,num_emis_ant) )
    emiss_ab = 0.
    allocate( emiss_ab1(jms:jme,num_emis_ant) )
    emiss_ab1 = 0.
    allocate( emiss_abu(jms:jme,num_ebu_in) )
    emiss_abu = 0.
    allocate( emiss_ash_mass(jms:jme) )
    emiss_ash_mass = 0.
    allocate( emiss_ash_height(jms:jme) )
    emiss_ash_height = 0.
    allocate( emiss_ash_dt(jms:jme) )
    emiss_ash_dt = 0.
    allocate( emiss_co2(jms:jme) )
    emiss_co2 = 0.
    allocate( emiss_ch4(jms:jme) )
    emiss_ch4 = 0.
    allocate( emiss_sf6(jms:jme) )
    emiss_sf6 = 0.

    print *,'chem_alloc_workspace: step 1 ...'
!if(chem_opt == 500)then
      allocate( emiss_tr_mass(jms:jme) )
      emiss_tr_mass = 0.
      allocate( emiss_tr_height(jms:jme) )
      emiss_tr_height = 0.
      allocate( emiss_tr_dt(jms:jme) )
      emiss_tr_dt = 0.
      ALLOCATE( trfall( jms:jme, num_chem ) )
      trfall = 0.
! endif


    print *,'chem_alloc_workspace: step 2 ...'
    allocate( emiss_oc(jms:jme) )
    emiss_oc = 0.
    allocate( emiss_bc(jms:jme) )
    emiss_bc = 0.
    allocate( emiss_sulf(jms:jme) )
    emiss_sulf = 0.
    allocate( emiss_pm25(jms:jme) )
    emiss_pm25 = 0.
    allocate( emiss_pm10(jms:jme) )
    emiss_pm10 = 0.
    allocate(  dm0(jms:jme) )
    dm0 = 0.
    allocate( emi_d1(jms:jme) )
    emi_d1 = 0.
    allocate( emi_d2(jms:jme) )
    emi_d2 = 0.
    allocate( emi_d3(jms:jme) )
    emi_d3 = 0.
    allocate( emi_d4(jms:jme) )
    emi_d4 = 0.
    allocate( emi_d5(jms:jme) )
    emi_d5 = 0.
    allocate( emid1_ave(jms:jme) )
    emid1_ave = 0.
    allocate( emid2_ave(jms:jme) )
    emid2_ave = 0.
    allocate( emid3_ave(jms:jme) )
    emid3_ave = 0.
    allocate( emid4_ave(jms:jme) )
    emid4_ave = 0.
    allocate( emid5_ave(jms:jme) )
    emid5_ave = 0.
    allocate( aod2d_ave(jms:jme) )
    aod2d_ave = 0.

    ! -- chem_alloc2

    print *,'chem_alloc_workspace: alloc2 ...'
  ALLOCATE( chem( ims:ime, kms:kme, jms:jme, num_chem ) )
  chem = 0.
  ALLOCATE( e_bio( ims:ime, jms:jme, ne_area ) )
  e_bio = 0.
  ALLOCATE( emis_ant( ims:ime, 1:kemit, jms:jme,num_emis_ant) )
emis_ant = 0.
  ALLOCATE( emis_vol( ims:ime, kms:kme, jms:jme,num_emis_vol) )
emis_vol=0.
  ALLOCATE( relhum( ims:ime, kms:kme, jms:jme ) )
  relhum = 0.
  ALLOCATE( dms_0( ims:ime, jms:jme) )
  dms_0 = 0.
  ALLOCATE( erod( ims:ime, jms:jme,3) )
  erod = 0.
  ALLOCATE( emis_dust( ims:ime, 1, jms:jme,num_emis_dust) )
emis_dust = 0.
  ALLOCATE( srce_dust( ims:ime, 1, jms:jme,5) ) 
srce_dust = 0. 
  ALLOCATE( emis_seas( ims:ime, 1, jms:jme,num_emis_seas) )
emis_seas = 0.
  ALLOCATE( backg_oh( ims:ime, kms:kme, jms:jme ) )
backg_oh = 0.
  ALLOCATE( backg_h2o2( ims:ime, kms:kme, jms:jme ) )
backg_h2o2 = 0.
  ALLOCATE( backg_no3( ims:ime, kms:kme, jms:jme ) )
backg_no3 = 0.
  ALLOCATE( oh_t( ims:ime, kms:kme, jms:jme ) )
oh_t = 0.
  ALLOCATE( h2o2_t( ims:ime, kms:kme, jms:jme ) )
h2o2_t = 0.
  ALLOCATE( no3_t( ims:ime, kms:kme, jms:jme ) )
no3_t = 0.
  ALLOCATE( h2oai(ims:ime, kms:kme, jms:jme ) )
  h2oai = 0.
  ALLOCATE( h2oaj(ims:ime, kms:kme, jms:jme ) )
  h2oaj = 0.
  ALLOCATE( nu3(ims:ime, kms:kme, jms:jme ) )
  nu3 = 0.
  ALLOCATE( ac3(ims:ime, kms:kme, jms:jme ) )
  ac3 = 0.
  ALLOCATE( cor3(ims:ime, kms:kme, jms:jme ) )
  cor3 = 0.
  ALLOCATE( asulf(ims:ime, kms:kme, jms:jme ) )
  asulf = 0.
  ALLOCATE( ahno3(ims:ime, kms:kme, jms:jme ) )
  ahno3 = 0.
  ALLOCATE( anh3(ims:ime, kms:kme, jms:jme ) )
  anh3 = 0.
!
  ALLOCATE( ebu_in( ims:ime, jms:jme, num_ebu_in ) )
  ebu_in=0.
  ALLOCATE( ebu( ims:ime, kms:kme, jms:jme,num_ebu ) )
  ebu=0.
  ALLOCATE( mean_fct_agtf( ims:ime,  jms:jme ) )
  mean_fct_agtf = 0.
  ALLOCATE( mean_fct_agef( ims:ime,  jms:jme ) )
  mean_fct_agef = 0.
  ALLOCATE( mean_fct_agsv( ims:ime,  jms:jme ) )
  mean_fct_agsv = 0.
  ALLOCATE( mean_fct_aggr( ims:ime,  jms:jme ) )
  mean_fct_aggr = 0.
  ALLOCATE( firesize_agtf( ims:ime,  jms:jme ) )
  firesize_agtf = 0.
  ALLOCATE( firesize_agef( ims:ime,  jms:jme ) )
  firesize_agef = 0.
  ALLOCATE( firesize_agsv( ims:ime,  jms:jme ) )
  firesize_agsv = 0.
  ALLOCATE( firesize_aggr( ims:ime,  jms:jme ) )
  firesize_aggr = 0.
  ALLOCATE( ash_fall( ims:ime,  jms:jme ) )
  ash_fall=0.
  ALLOCATE( dust_fall( ims:ime,  jms:jme ) )
  dust_fall = 0.
  ALLOCATE( pm2_5_dry( ims:ime , kms:kme , jms:jme ) )
  pm2_5_dry=0.
  ALLOCATE( pm2_5_water( ims:ime , kms:kme , jms:jme ) )
  pm2_5_water=0.
  ALLOCATE( aerwrf( ims:ime , kms:kme , jms:jme ) )
  aerwrf=0.
  ALLOCATE( pm2_5_dry_ec( ims:ime , kms:kme , jms:jme ) )
  pm2_5_dry_ec = 0.
  ALLOCATE( pm10( ims:ime , kms:kme , jms:jme ) )
  pm10 = 0.
  ALLOCATE( tcosz( ims:ime , jms:jme ) )
  tcosz = 0.
  ALLOCATE( ttday( ims:ime , jms:jme ) )
  ttday = 0.

  ALLOCATE( sebio_iso( ims:ime , jms:jme ) )
  ALLOCATE( sebio_oli( ims:ime , jms:jme ) )
  ALLOCATE( sebio_api( ims:ime , jms:jme ) )
  ALLOCATE( sebio_lim( ims:ime , jms:jme ) )
  ALLOCATE( sebio_xyl( ims:ime , jms:jme ) )
  ALLOCATE( sebio_hc3( ims:ime , jms:jme ) )
  ALLOCATE( sebio_ete( ims:ime , jms:jme ) )
  ALLOCATE( sebio_olt( ims:ime , jms:jme ) )
  ALLOCATE( sebio_ket( ims:ime , jms:jme ) )
  ALLOCATE( sebio_ald( ims:ime , jms:jme ) )
  ALLOCATE( sebio_hcho( ims:ime , jms:jme ) )
  ALLOCATE( sebio_eth( ims:ime , jms:jme ) )
  ALLOCATE( sebio_ora2( ims:ime , jms:jme ) )
  ALLOCATE( sebio_co( ims:ime , jms:jme ) )
  ALLOCATE( sebio_nr( ims:ime , jms:jme ) )
  ALLOCATE( noag_grow( ims:ime , jms:jme ) )
  ALLOCATE( noag_nongrow( ims:ime , jms:jme ) )
  ALLOCATE( nononag( ims:ime , jms:jme ) )
  ALLOCATE( slai( ims:ime , jms:jme ) )
  ALLOCATE( ebio_iso( ims:ime , jms:jme ) )
  ALLOCATE( ebio_oli( ims:ime , jms:jme ) )
  ALLOCATE( ebio_api( ims:ime , jms:jme ) )
  ALLOCATE( ebio_lim( ims:ime , jms:jme ) )
  ALLOCATE( ebio_xyl( ims:ime , jms:jme ) )
  ALLOCATE( ebio_hc3( ims:ime , jms:jme ) )
  ALLOCATE( ebio_ete( ims:ime , jms:jme ) )
  ALLOCATE( ebio_olt( ims:ime , jms:jme ) )
  ALLOCATE( ebio_ket( ims:ime , jms:jme ) )
  ALLOCATE( ebio_ald( ims:ime , jms:jme ) )
  ALLOCATE( ebio_hcho( ims:ime , jms:jme ) )
  ALLOCATE( ebio_eth( ims:ime , jms:jme ) )
  ALLOCATE( ebio_ora2( ims:ime , jms:jme ) )
  ALLOCATE( ebio_co( ims:ime , jms:jme ) )
  ALLOCATE( ebio_nr( ims:ime , jms:jme ) )
  ALLOCATE( ebio_no( ims:ime , jms:jme ) )
! ------

  sebio_iso = 0.
  sebio_oli = 0.
  sebio_api = 0.
  sebio_lim = 0.
  sebio_xyl = 0.
  sebio_hc3 = 0.
  sebio_ete = 0.
  sebio_olt = 0.
  sebio_ket = 0.
  sebio_ald = 0.
  sebio_hcho = 0.
  sebio_eth = 0.
  sebio_ora2 = 0.
  sebio_co = 0.
  sebio_nr = 0.
  noag_grow = 0.
  noag_nongrow = 0.
  nononag = 0.
  slai    = 0.
  ebio_iso = 0.
  ebio_oli = 0.
  ebio_api = 0.
  ebio_lim = 0.
  ebio_xyl = 0.
  ebio_hc3 = 0.
  ebio_ete = 0.
  ebio_olt = 0.
  ebio_ket = 0.
  ebio_ald = 0.
  ebio_hcho = 0.
  ebio_eth  = 0.
  ebio_ora2 = 0.
  ebio_co = 0.
  ebio_nr = 0.
  ebio_no = 0.
! ------

  if(chem_config % bio_emiss_opt == 3)then
  ALLOCATE( EFmegan(ims:ime, jms:jme , nmegan) )

  ALLOCATE( msebio_isop(ims:ime, jms:jme ) )
  ALLOCATE( pftp_bt(ims:ime, jms:jme ) )
  ALLOCATE( pftp_nt(ims:ime, jms:jme ) )
  ALLOCATE( pftp_sb(ims:ime, jms:jme ) )
  ALLOCATE( pftp_hb(ims:ime, jms:jme ) )

  ALLOCATE( mlai(ims:ime, jms:jme, 12 ) )
  ALLOCATE( mtsa(ims:ime, jms:jme, 12 ) )
  ALLOCATE( mswdown(ims:ime, jms:jme, 12 ) )

  ALLOCATE( mebio_isop(ims:ime, jms:jme ) )
  ALLOCATE( mebio_apin(ims:ime, jms:jme ) )
  ALLOCATE( mebio_bpin(ims:ime, jms:jme ) )
  ALLOCATE( mebio_bcar(ims:ime, jms:jme ) )
  ALLOCATE( mebio_acet(ims:ime, jms:jme ) )
  ALLOCATE( mebio_mbo(ims:ime, jms:jme ) )
  ALLOCATE( mebio_no(ims:ime, jms:jme ) )

  EFmegan = 0.

  msebio_isop = 0.
  pftp_bt = 0.
  pftp_nt = 0.
  pftp_sb = 0.
  pftp_hb = 0.

  mlai = 0.
  mtsa = 0.
  mswdown = 0.

  mebio_isop = 0.
  mebio_apin = 0.
  mebio_bpin = 0.
  mebio_bcar = 0.
  mebio_acet = 0.
  mebio_mbo  = 0.
  mebio_no   = 0.
  endif


    print *,'chem_alloc_workspace: alloc2 step 2 ...'
  if((chem_config % chem_opt == CHEM_OPT_GOCART_RACM) .or. &
     (chem_config % chem_opt == CHEM_OPT_RACM_SOA_VBS)) then
   allocate( ph_o31d(ims:ime, kms:kme, jms:jme ) )
   ph_o31d=0.
   allocate( ph_o33p(ims:ime, kms:kme, jms:jme ) )
   ph_o33p=0.
   allocate( ph_no2(ims:ime, kms:kme, jms:jme ) )
   ph_no2=0.
   allocate( ph_no3o2(ims:ime, kms:kme, jms:jme ) )
   ph_no3o2=0.
   allocate( ph_no3o(ims:ime, kms:kme, jms:jme ) )
   ph_no3o=0.
   allocate( ph_hno2(ims:ime, kms:kme, jms:jme ) )
   ph_hno2=0.
   allocate( ph_hno3(ims:ime, kms:kme, jms:jme ) )
   ph_hno3=0.
   allocate( ph_hno4(ims:ime, kms:kme, jms:jme ) )
   ph_hno4=0.
   allocate( ph_h2o2(ims:ime, kms:kme, jms:jme ) )
   ph_h2o2=0.
   allocate( ph_ch2or(ims:ime, kms:kme, jms:jme ) )
   ph_ch2or=0.
   allocate( ph_ch2om(ims:ime, kms:kme, jms:jme ) )
   ph_ch2om=0.
   allocate( ph_ch3cho(ims:ime, kms:kme, jms:jme ) )
   ph_ch3cho=0.
   allocate( ph_ch3coch3(ims:ime, kms:kme, jms:jme ) )
   ph_ch3coch3=0.
   allocate( ph_ch3coc2h5(ims:ime, kms:kme, jms:jme ) )
   ph_ch3coc2h5=0.
   allocate( ph_hcocho(ims:ime, kms:kme, jms:jme ) )
   ph_hcocho=0.
   allocate( ph_ch3cocho(ims:ime, kms:kme, jms:jme ) )
   ph_ch3cocho=0.
   allocate( ph_hcochest(ims:ime, kms:kme, jms:jme ) )
   ph_hcochest=0.
   allocate( ph_ch3o2h(ims:ime, kms:kme, jms:jme ) )
   ph_ch3o2h=0.
   allocate( ph_ch3coo2h(ims:ime, kms:kme, jms:jme ) )
   ph_ch3coo2h=0.
   allocate( ph_ch3ono2(ims:ime, kms:kme, jms:jme ) )
   ph_ch3ono2=0.
   allocate( ph_hcochob(ims:ime, kms:kme, jms:jme ) )
   ph_hcochob=0.
   allocate( ph_macr(ims:ime, kms:kme, jms:jme ) )
   ph_macr=0.
   allocate( ph_n2o5(ims:ime, kms:kme, jms:jme ) )
   ph_n2o5=0.
   allocate( ph_o2(ims:ime, kms:kme, jms:jme ) )
   ph_o2=0.
   allocate( ph_pan(ims:ime, kms:kme, jms:jme ) )
   ph_pan=0.
   allocate( ph_acet(ims:ime, kms:kme, jms:jme ) )
   ph_acet=0.
   allocate( ph_mglo(ims:ime, kms:kme, jms:jme ) )
   ph_mglo=0.
   allocate( ph_hno4_2(ims:ime, kms:kme, jms:jme ) )
   ph_hno4_2=0.
   allocate( ph_n2o (ims:ime, kms:kme, jms:jme ) )
   ph_n2o=0.
   allocate( ph_pooh (ims:ime, kms:kme, jms:jme ) )
   ph_pooh=0.
   allocate( ph_mpan (ims:ime, kms:kme, jms:jme ) )
   ph_mpan=0.
   allocate( ph_mvk (ims:ime, kms:kme, jms:jme ) )
   ph_mvk=0.
   allocate( ph_etooh (ims:ime, kms:kme, jms:jme ) )
   ph_etooh=0.
   allocate( ph_prooh (ims:ime, kms:kme, jms:jme ) )
   ph_prooh=0.
   allocate( ph_onitr (ims:ime, kms:kme, jms:jme ) )
   ph_onitr=0.
   allocate( ph_acetol (ims:ime, kms:kme, jms:jme ) )
   ph_acetol=0.
   allocate( ph_glyald (ims:ime, kms:kme, jms:jme ) )
   ph_glyald=0.
   allocate( ph_hyac (ims:ime, kms:kme, jms:jme ) )
   ph_hyac=0.
   allocate( ph_mek (ims:ime, kms:kme, jms:jme ) )
   ph_mek=0.
   allocate( ph_open (ims:ime, kms:kme, jms:jme ) )
   ph_open=0.
   allocate( ph_gly (ims:ime, kms:kme, jms:jme ) )
   ph_gly=0.
   allocate( ph_acetp (ims:ime, kms:kme, jms:jme ) )
   ph_acetp=0.
   allocate( ph_xooh (ims:ime, kms:kme, jms:jme ) )
   ph_xooh=0.
   allocate( ph_isooh (ims:ime, kms:kme, jms:jme ) )
   ph_isooh=0.
   allocate( ph_alkooh (ims:ime, kms:kme, jms:jme ) )
   ph_alkooh=0.
   allocate( ph_mekooh (ims:ime, kms:kme, jms:jme ) )
   ph_mekooh=0.
   allocate( ph_tolooh (ims:ime, kms:kme, jms:jme ) )
   ph_tolooh=0.
   allocate( ph_terpooh (ims:ime, kms:kme, jms:jme ) )
   ph_terpooh=0.
   allocate( ph_cl2 (ims:ime, kms:kme, jms:jme ) )
   ph_cl2=0.
   allocate( ph_hocl (ims:ime, kms:kme, jms:jme ) )
   ph_hocl=0.
   allocate( ph_fmcl (ims:ime, kms:kme, jms:jme ) )
   ph_fmcl=0.

endif


    print *,'chem_alloc_workspace: alloc2 step 3 ...'
  if(chem_config % chem_opt == 2)then
     ALLOCATE( h2oai(ims:ime, kms:kme, jms:jme ) )
     ALLOCATE( h2oaj(ims:ime, kms:kme, jms:jme ) )
     h2oai = 0.
     h2oaj = 0.
  endif
  if(chem_config % aer_ra_feedback == 1)then
     ALLOCATE( extt(ims:ime, kms:kme, jms:jme,nbands) )
     ALLOCATE( ssca(ims:ime, kms:kme, jms:jme,nbands) )
     ALLOCATE( asympar(ims:ime, kms:kme, jms:jme,nbands) )
     ALLOCATE( aod(ims:ime, jms:jme ) )
     ALLOCATE( ext_coeff(ims:ime, kms:kme, jms:jme,1:num_ext_coef ) )
     ALLOCATE( bscat_coeff(ims:ime, kms:kme, jms:jme,1:num_bscat_coef ) )
     ALLOCATE( asym_par(ims:ime, kms:kme, jms:jme,1:num_asym_par ) )
     ALLOCATE( tauaerlw(ims:ime, kms:kme, jms:jme,1:16 ) )
     ALLOCATE( tauaersw(ims:ime, kms:kme, jms:jme,1:4 ) )
     ALLOCATE( gaersw(ims:ime, kms:kme, jms:jme, 1:4 ) )
     ALLOCATE( waersw(ims:ime, kms:kme, jms:jme, 1:4 ) )
     ALLOCATE( bscoefsw(ims:ime, kms:kme, jms:jme, 1:4 ) )
     ALLOCATE( l2aer(ims:ime, kms:kme, jms:jme, 1:4 ) )
     ALLOCATE( l3aer(ims:ime, kms:kme, jms:jme, 1:4 ) )
     ALLOCATE( l4aer(ims:ime, kms:kme, jms:jme, 1:4 ) )
     ALLOCATE( l5aer(ims:ime, kms:kme, jms:jme, 1:4 ) )
     ALLOCATE( l6aer(ims:ime, kms:kme, jms:jme, 1:4 ) )
     ALLOCATE( l7aer(ims:ime, kms:kme, jms:jme, 1:4 ) )

     extt = 0.
     ssca = 0.
     asympar = 0.
     aod = 0.
     ext_coeff = 0.
     bscat_coeff = 0.
     asym_par = 0.
     tauaerlw = 0.
     tauaersw = 0.
     gaersw = 0.
     waersw = 0.
     bscoefsw = 0.
     l2aer = 0.
     l3aer = 0.
     l4aer = 0.
     l5aer = 0.
     l6aer = 0.
     l7aer = 0.
  endif
! if (chem_config % aer_ra_feedback > 0) then
    allocate(ext_cof(kms:kme,jms:jme,nbands))
    allocate(sscal  (kms:kme,jms:jme,nbands))
    allocate(asymp  (kms:kme,jms:jme,nbands))
    ext_cof = 0.
    asymp   = 0.
    sscal   = 1.
! end if
  print *,'chem_alloc_workspace: alloc2 step 4 ...'
  if((chem_config % chem_opt == CHEM_OPT_GOCART_RACM) .or. &
     (chem_config % chem_opt == CHEM_OPT_RACM_SOA_VBS)) then

     ALLOCATE( addt(ims:ime, kms:kme, jms:jme ) )
     addt=0.
     ALLOCATE( addx(ims:ime, kms:kme, jms:jme ) )
     addx=0.
     ALLOCATE( addc(ims:ime, kms:kme, jms:jme ) )
     addc=0.
     ALLOCATE( etep(ims:ime, kms:kme, jms:jme ) )
     etep=0.
     ALLOCATE( oltp(ims:ime, kms:kme, jms:jme ) )
     oltp=0.
     ALLOCATE( olip(ims:ime, kms:kme, jms:jme ) )
     olip=0.
     ALLOCATE( cslp(ims:ime, kms:kme, jms:jme ) )
     cslp=0.
     ALLOCATE( limp(ims:ime, kms:kme, jms:jme ) )
     limp=0.
     ALLOCATE( hc5p(ims:ime, kms:kme, jms:jme ) )
     hc5p=0.
     ALLOCATE( hc8p(ims:ime, kms:kme, jms:jme ) )
     hc8p=0.
     ALLOCATE( tolp(ims:ime, kms:kme, jms:jme ) )
     tolp=0.
     ALLOCATE( xylp(ims:ime, kms:kme, jms:jme ) )
     xylp=0.
     ALLOCATE( apip(ims:ime, kms:kme, jms:jme ) )
     apip=0.
     ALLOCATE( isop(ims:ime, kms:kme, jms:jme ) )
     isop=0.
     ALLOCATE( hc3p(ims:ime, kms:kme, jms:jme ) )
     hc3p=0.
     ALLOCATE( ethp(ims:ime, kms:kme, jms:jme ) )
     ethp=0.
     ALLOCATE( o3p(ims:ime, kms:kme, jms:jme ) )
     o3p=0.
     ALLOCATE( tco3(ims:ime, kms:kme, jms:jme ) )
     tco3=0.
     ALLOCATE( mo2(ims:ime, kms:kme, jms:jme ) )
     mo2=0.
     ALLOCATE( o1d(ims:ime, kms:kme, jms:jme ) )
     o1d=0.
     ALLOCATE( olnn(ims:ime, kms:kme, jms:jme ) )
     olnn=0.
     ALLOCATE( olnd(ims:ime, kms:kme, jms:jme ) )
     olnd=0.
     ALLOCATE( rpho(ims:ime, kms:kme, jms:jme ) )
     rpho=0.
     ALLOCATE( xo2(ims:ime, kms:kme, jms:jme ) )
     xo2=0.
     ALLOCATE( ketp(ims:ime, kms:kme, jms:jme ) )
     ketp=0.
     ALLOCATE( xno2(ims:ime, kms:kme, jms:jme ) )
     xno2=0.
     ALLOCATE( ol2p(ims:ime, kms:kme, jms:jme ) )
     ol2p=0.
     ALLOCATE( oln(ims:ime, kms:kme, jms:jme ) )
     oln=0.
     ALLOCATE( macp(ims:ime, kms:kme, jms:jme ) )
     macp=0.
     ALLOCATE( hocoo(ims:ime, kms:kme, jms:jme ) )
     hocoo=0.
     ALLOCATE( bzno2_o(ims:ime, kms:kme, jms:jme ) )
     bzno2_o=0.
     ALLOCATE( bz_o(ims:ime, kms:kme, jms:jme ) )
     bz_o=0.
     ALLOCATE( tbu_o(ims:ime, kms:kme, jms:jme ) )
     tbu_o=0.
  endif

    print *,'chem_alloc_workspace: alloc2 step 4 ...'
  if (chem_config % chem_opt == CHEM_OPT_RACM_SOA_VBS) then
     ALLOCATE( cvaro1(ims:ime, kms:kme, jms:jme ) )
     cvaro1=0.
     ALLOCATE( cvaro2(ims:ime, kms:kme, jms:jme ) )
     cvaro2=0.
     ALLOCATE( cvalk1(ims:ime, kms:kme, jms:jme ) )
     cvalk1=0.
     ALLOCATE( cvole1(ims:ime, kms:kme, jms:jme ) )
     cvole1=0.
     ALLOCATE( cvapi1(ims:ime, kms:kme, jms:jme ) )
     cvapi1=0.
     ALLOCATE( cvapi2(ims:ime, kms:kme, jms:jme ) )
     cvapi2=0.
     ALLOCATE( cvlim1(ims:ime, kms:kme, jms:jme ) )
     cvlim1=0.
     ALLOCATE( cvlim2(ims:ime, kms:kme, jms:jme ) )
     cvlim2=0.
     ALLOCATE( mob(ims:ime, kms:kme, jms:jme ) )
     mob=0.
     ALLOCATE( cvasoaX(ims:ime, kms:kme, jms:jme ) )
     cvasoaX=0.
     ALLOCATE( cvasoa1(ims:ime, kms:kme, jms:jme ) )
     cvasoa1=0.
     ALLOCATE( cvasoa2(ims:ime, kms:kme, jms:jme ) )
     cvasoa2=0.
     ALLOCATE( cvasoa3(ims:ime, kms:kme, jms:jme ) )
     cvasoa3=0.
     ALLOCATE( cvasoa4(ims:ime, kms:kme, jms:jme ) )
     cvasoa4=0.
     ALLOCATE( cvbsoaX(ims:ime, kms:kme, jms:jme ) )
     cvbsoaX=0.
     ALLOCATE( cvbsoa1(ims:ime, kms:kme, jms:jme ) )
     cvbsoa1=0.
     ALLOCATE( cvbsoa2(ims:ime, kms:kme, jms:jme ) )
     cvbsoa2=0.
     ALLOCATE( cvbsoa3(ims:ime, kms:kme, jms:jme ) )
     cvbsoa3=0.
     ALLOCATE( cvbsoa4(ims:ime, kms:kme, jms:jme ) )
     cvbsoa4=0.
     ALLOCATE( asoa1j(ims:ime, kms:kme, jms:jme ) )
     asoa1j=0.
     ALLOCATE( asoa1i(ims:ime, kms:kme, jms:jme ) )
     asoa1i=0.
     ALLOCATE( asoa2j(ims:ime, kms:kme, jms:jme ) )
     asoa2j=0.
     ALLOCATE( asoa2i(ims:ime, kms:kme, jms:jme ) )
     asoa2i=0.
     ALLOCATE( asoa3j(ims:ime, kms:kme, jms:jme ) )
     asoa3j=0.
     ALLOCATE( asoa3i(ims:ime, kms:kme, jms:jme ) )
     asoa3i=0.
     ALLOCATE( asoa4j(ims:ime, kms:kme, jms:jme ) )
     asoa4j=0.
     ALLOCATE( asoa4i(ims:ime, kms:kme, jms:jme ) )
     asoa4i=0.
     ALLOCATE( bsoa1j(ims:ime, kms:kme, jms:jme ) )
     bsoa1j=0.
     ALLOCATE( bsoa1i(ims:ime, kms:kme, jms:jme ) )
     bsoa1i=0.
     ALLOCATE( bsoa2j(ims:ime, kms:kme, jms:jme ) )
     bsoa2j=0.
     ALLOCATE( bsoa2i(ims:ime, kms:kme, jms:jme ) )
     bsoa2i=0.
     ALLOCATE( bsoa3j(ims:ime, kms:kme, jms:jme ) )
     bsoa3j=0.
     ALLOCATE( bsoa3i(ims:ime, kms:kme, jms:jme ) )
     bsoa3i=0.
     ALLOCATE( bsoa4j(ims:ime, kms:kme, jms:jme ) )
     bsoa4j=0.
     ALLOCATE( bsoa4i(ims:ime, kms:kme, jms:jme ) )
     bsoa4i=0.
  endif
    print *,'chem_alloc_workspace: done'
    
  end subroutine chem_alloc_workspace

  subroutine chem_alloc_prep

    allocate(rri(ims:ime, kms:kme, jms:jme))
    allocate(t_phy(ims:ime, kms:kme, jms:jme))
    allocate(p_phy(ims:ime, kms:kme, jms:jme))
!   allocate(relhum(ims:ime, kms:kme, jms:jme))
    allocate(dz8w(ims:ime, kms:kme, jms:jme))
    allocate(p8w(ims:ime, kms:kme, jms:jme))
    allocate(t8w(ims:ime, kms:kme, jms:jme))
    allocate(z_at_w(ims:ime, kms:kme, jms:jme))
    allocate(zmid(ims:ime, kms:kme, jms:jme))
    allocate(exch_h(ims:ime, kms:kme, jms:jme))
    allocate(u_phy(ims:ime, kms:kme, jms:jme))
    allocate(v_phy(ims:ime, kms:kme, jms:jme))
    allocate(vvel(ims:ime, kms:kme, jms:jme))
    allocate(rho_phy(ims:ime, kms:kme, jms:jme))
    allocate(convfac(ims:ime, kms:kme, jms:jme))
    allocate(raincv_b(ims:ime, jms:jme))
    allocate(pbl(ims:ime, jms:jme))
    allocate(hfx(ims:ime, jms:jme))
    allocate(snowh(ims:ime, jms:jme))
    allocate(ust(ims:ime, jms:jme))
    allocate(tsk(ims:ime, jms:jme))
    allocate(gsw(ims:ime, jms:jme))
    allocate(znt(ims:ime, jms:jme))
    allocate(rmol(ims:ime, jms:jme))
    allocate(dxy(ims:ime, jms:jme))
    allocate(u10(ims:ime, jms:jme))
    allocate(v10(ims:ime, jms:jme))
    allocate(xlong(ims:ime, jms:jme))
    allocate(xlat(ims:ime, jms:jme))
    allocate(xland(ims:ime, jms:jme))
    allocate(ivgtyp(ims:ime, jms:jme))
    allocate(isltyp(ims:ime, jms:jme))
    allocate(vegfra(ims:ime, jms:jme))
    allocate(clayf(ims:ime, jms:jme))
    allocate(sandf(ims:ime, jms:jme))
    allocate(moist(ims:ime, kms:kme, jms:jme, num_moist ))
    allocate(smois( ims:ime, num_soil_layers, jms:jme ))

    rri     = 0.
    t_phy   = 0.
    p_phy   = 0.
    dz8w    = 0.
    p8w     = 0.
    t8w     = 0.
    z_at_w  = 0.
    zmid    = 0.
    exch_h  = 0.
    u_phy   = 0.
    v_phy   = 0.
    vvel    = 0.
    rho_phy = 0.
    convfac = 0.
    raincv_b = 0.
    pbl      = 0.
    hfx      = 0.
    snowh    = 0.
    ust      = 0.
    tsk      = 0.
    gsw      = 0.
    znt      = 0.
    rmol     = 0.
    dxy      = 0.
    u10      = 0.
    v10      = 0.
    xlong    = 0.
    xlat     = 0.
    xland    = 0.
    ivgtyp   = 0.
    isltyp   = 0.
    vegfra   = 0.
    clayf    = 0.
    sandf    = 0.
    moist    = 0.
    smois    = 0.

  end subroutine chem_alloc_prep

  subroutine chem_alloc_seas
    allocate(seashelp(ims:ime, jms:jme))
    seashelp = 0.
  end subroutine chem_alloc_seas

  subroutine chem_allocate
    call chem_alloc_input
    call chem_alloc_workspace
    call chem_alloc_prep
    call chem_alloc_seas
  end subroutine chem_allocate

  subroutine chem_free_memory

    deallocate(rn2d,rc2d,ts2d,us2d,hf2d,pb2d,zorl2d,vfrac2d,vtype2d,stype2d,snwdph2d)
    deallocate(aod2d,rsds,area,deg_lat,deg_lon,sm3d,tr3d)
    deallocate(trdp,us3d,vs3d,ws3d,tk3d,exch,dp3d,pr3d,ph3d,gd_cloud,gd_cldfr)
    deallocate(pm25,p10,tr1_tavg,d1st_ave,d2st_ave,d3st_ave,d4st_ave,d5st_ave)
    deallocate(ebu_oc,oh_bg,h2o2_bg,no3_bg,oh_backgd,h2o2_backgd,no3_backgd,rcav,rnav,ero1)
    deallocate(ero2,ero3,clayfrac,sandfrac,ashfall,aod2d,wet_dep,dry_dep,plumestuff,emiss_ab)
    deallocate(emiss_ab1,emiss_abu,emiss_ash_mass,emiss_ash_height,emiss_ash_dt,emiss_co2,emiss_ch4,emiss_sf6,emiss_tr_mass,emiss_tr_height)
    deallocate(emiss_tr_dt,trfall,emiss_oc,emiss_bc,emiss_sulf,emiss_pm25,emiss_pm10,dm0,emi_d1,emi_d2)
    deallocate(emi_d3,emi_d4,emi_d5,emid1_ave,emid2_ave,emid3_ave,emid4_ave,emid5_ave,aod2d_ave,chem)
    deallocate(e_bio,emis_ant,emis_vol,relhum,dms_0,erod,emis_dust,srce_dust,emis_seas,backg_oh)
    deallocate(backg_h2o2,backg_no3,oh_t,h2o2_t,no3_t,h2oai,h2oaj,nu3,ac3,cor3)
    deallocate(asulf,ahno3,anh3,ebu_in,ebu,mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,firesize_agtf)
    deallocate(firesize_agef,firesize_agsv,firesize_aggr,ash_fall,dust_fall,pm2_5_dry,pm2_5_water,aerwrf,pm2_5_dry_ec,pm10)
    deallocate(tcosz,ttday,sebio_iso,sebio_oli,sebio_api,sebio_lim,sebio_xyl,sebio_hc3,sebio_ete,sebio_olt)
    deallocate(sebio_ket,sebio_ald,sebio_hcho,sebio_eth,sebio_ora2,sebio_co,sebio_nr,noag_grow,noag_nongrow,nononag)
    deallocate(slai,ebio_iso,ebio_oli,ebio_api,ebio_lim,ebio_xyl,ebio_hc3,ebio_ete,ebio_olt,ebio_ket)
    deallocate(ebio_ald,ebio_hcho,ebio_eth,ebio_ora2,ebio_co,ebio_nr,ebio_no,EFmegan,msebio_isop,pftp_bt)
    deallocate(pftp_nt,pftp_sb,pftp_hb,mlai,mtsa,mswdown,mebio_isop,mebio_apin,mebio_bpin,mebio_bcar)
    deallocate(mebio_acet,mebio_mbo,mebio_no,ph_o31d,ph_o33p,ph_no2,ph_no3o2,ph_no3o,ph_hno2,ph_hno3)
    deallocate(ph_hno4,ph_h2o2,ph_ch2or,ph_ch2om,ph_ch3cho,ph_ch3coch3,ph_ch3coc2h5,ph_hcocho,ph_ch3cocho,ph_hcochest)
    deallocate(ph_ch3o2h,ph_ch3coo2h,ph_ch3ono2,ph_hcochob,ph_macr,ph_n2o5,ph_o2,ph_pan,ph_acet,ph_mglo)
    deallocate(ph_hno4_2,ph_n2o,ph_pooh,ph_mpan,ph_mvk,ph_etooh,ph_prooh,ph_onitr,ph_acetol,ph_glyald)
    deallocate(ph_hyac,ph_mek,ph_open,ph_gly,ph_acetp,ph_xooh,ph_isooh,ph_alkooh,ph_mekooh,ph_tolooh)
    deallocate(ph_terpooh,ph_cl2,ph_hocl,ph_fmcl,h2oai,h2oaj,extt,ssca,asympar,aod)
    deallocate(ext_coeff,bscat_coeff,asym_par,tauaerlw,tauaersw,gaersw,waersw,bscoefsw,l2aer,l3aer)
    deallocate(l4aer,l5aer,l6aer,l7aer,ext_cof,sscal,asymp,addt,addx,addc)
    deallocate(etep,oltp,olip,cslp,limp,hc5p,hc8p,tolp,xylp,apip)
    deallocate(isop,hc3p,ethp,o3p,tco3,mo2,o1d,olnn,olnd,rpho)
    deallocate(xo2,ketp,xno2,ol2p,oln,macp,hocoo,bzno2_o,bz_o,tbu_o)
    deallocate(cvaro1,cvaro2,cvalk1,cvole1,cvapi1,cvapi2,cvlim1,cvlim2,mob,cvasoaX)
    deallocate(cvasoa1,cvasoa2,cvasoa3,cvasoa4,cvbsoaX,cvbsoa1,cvbsoa2,cvbsoa3,cvbsoa4,asoa1j)
    deallocate(asoa1i,asoa2j,asoa2i,asoa3j,asoa3i,asoa4j,asoa4i,bsoa1j,bsoa1i,bsoa2j)
    deallocate(bsoa2i,bsoa3j,bsoa3i,bsoa4j,bsoa4i,rri,t_phy,p_phy,relhum,dz8w)
    deallocate(p8w,t8w,z_at_w,zmid,exch_h,u_phy,v_phy,vvel,rho_phy,convfac)
    deallocate(raincv_b,pbl,hfx,snowh,ust,tsk,gsw,znt,rmol,dxy,u10)
    deallocate(v10,xlong,xlat,xland,ivgtyp,isltyp,vegfra,clayf,sandf,slmsk2d)
    deallocate(moist,smois,seashelp)
    deallocate(perm)

  end subroutine chem_free_memory

  subroutine chem_background_init

    ! -- local variables
    integer :: rc
    integer :: p, pos

    integer, parameter :: max_v = 29

    character(len=*), dimension(max_v), parameter :: &
      in_file = (/ 'vash1.in', 'vash2.in', 'vash3.in', 'vash4.in', 'vash5.in',  &
                   'vash6.in', 'vash7.in', 'vash8.in', 'vash9.in', 'vash10.in', &
                   'so2.in', 'sulf.in', 'dms.in', 'msa.in', 'p10.in', &
                   'p25.in', 'bc1.in', 'bc2.in', 'oc1.in', 'oc2.in', &
                   'dust1.in', 'dust2.in', 'dust3.in', 'dust4.in', 'dust5.in', &
                   'seas1.in', 'seas2.in', 'seas3.in', 'seas4.in' /)

    character (len=4), parameter :: chem_301(49) = (/ &
                           'so2','sulf','no2','no','o3',        &
                          'hno3','h2o2','ald','hcho','op1','op2',  &
                          'paa','ora1','ora2','nh3','n2o5','no3',  &
                          'pan','hc3','hc5','hc8','eth','co',  &
                          'ete','olt','oli','tol','xyl','aco3',  &
                          'tpan','hono','hno4','ket','gly','mgly',  &
                          'dcb','onit','csl','iso','co2','ch4',  &
                          'udd','hket','api','lim','dien','macr',  &
                          'ho','ho2'/)

     character (len=4), parameter :: chem_108(103) = (/ &
                           'so2','sulf','no2','no','o3',        &
                          'hno3','h2o2','ald','hcho','op1','op2',  &
                          'paa','ora1','ora2','nh3','n2o5','no3',  &
                          'pan','hc3','hc5','hc8','eth','co',  &
                          'ete','olt','oli','tol','xyl','aco3',  &
                          'tpan','hono','hno4','ket','gly','mgly',  &
                          'dcb','onit','csl','iso','co2','ch4',  &
                          'udd','hket','api','lim','dien','macr',  &
                          'hace','ishp','ison','mahp','mpan','nald',&
                          'sesq','mbo','cva1','cva2','cva3',&
                          'cva4','cvb1','cvb2','cvb3',&
                          'cvb4','ho','ho2','soaj','soai','nhaj',&
                          'nhai','n3aj','n3ai','naaj','naai','claj',&
                          'clai','as1j','as1i','as2j','as2i',&
                          'as3j','as3i','as4j','as4i','bs1j',&
                          'bs1i','bs2j','bs2i','bs3j','bs3i',&
                          'bs4j','bs4i','opaj','opai','ecj','eci',&
                          'p25j','p25i','atha','seas','sila','nu0',&
                          'ac0','corn'/)

    integer, dimension(max_v) :: p_v

    ! -- begin
    p_v = (/ p_vash_1, p_vash_2, p_vash_3, p_vash_4, p_vash_5, &
             p_vash_6, p_vash_7, p_vash_8, p_vash_9, p_vash_10, &
             p_so2, p_sulf, p_dms, p_msa, p_p10, &
             p_p25, p_bc1, p_bc2, p_oc1, p_oc2, &
             p_dust_1, p_dust_2, p_dust_3, p_dust_4, p_dust_5, &
             p_seas_1, p_seas_2, p_seas_3, p_seas_4 /)

    if ((chem_config % chem_opt == CHEM_OPT_RACM_SOA_VBS) .or. &
        (chem_config % chem_opt >= CHEM_OPT_GOCART) .and. &
        (chem_config % chem_opt < 500)) then
      write(6,*)'reading gocart background fields'
      call chem_data_read('oh.dat', oh_backgd, rc, perm=perm, label='oh')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading oh')
      call chem_data_read('h2o2.dat', h2o2_backgd, rc, perm=perm, label='h2o2')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading h2o2')
      call chem_data_read('no3.dat', no3_backgd, rc, perm=perm, label='no3')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading no3')
      call chem_data_read_global('p_gocart.dat', p_gocart, rc, label='p_gocart')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading p_gocart')
      call chem_data_read('e_bc.dat', emiss_ab(:,p_e_bc), rc, perm=perm, label='e_bc')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_bc')
      call chem_data_read('e_oc.dat', emiss_ab(:,p_e_oc), rc, perm=perm, label='e_oc')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_oc')
      call chem_data_read('e_sulf.dat', emiss_ab(:,p_e_sulf), rc, perm=perm, label='e_sulf')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_sulf')
      call chem_data_read('e_pm_25.dat', emiss_ab(:,p_e_pm_25), rc, perm=perm, label='pm25')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_pm_25')
      call chem_data_read('e_pm_10.dat', emiss_ab(:,p_e_pm_10), rc, perm=perm, label='pm10')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_pm_10')
      call chem_data_read('dm0.dat', dm0, rc, perm=perm, label='dm0')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading dm0')
      call chem_data_read('erod1.dat', ero1, rc, perm=perm, label='erod1')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ero1')
      call chem_data_read('erod2.dat', ero2, rc, perm=perm, label='erod2')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ero2')
      call chem_data_read('erod3.dat', ero3, rc, perm=perm, label='erod3')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ero3')
      if (chem_config % dust_opt == DUST_OPT_AFWA) then
        call chem_data_read('sand.dat', sandfrac, rc, perm=perm, label='sandfrac')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading sandfrac')
        call chem_data_read('clay.dat', clayfrac, rc, perm=perm, label='clayfrac')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading clayfrac')
      end if
      call chem_data_read('e_so2.dat', emiss_ab(:,p_e_so2), rc, perm=perm, label='so2')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_so2')
      if ((chem_config % chem_opt == CHEM_OPT_GOCART_RACM) .or. &
          (chem_config % chem_opt == CHEM_OPT_RACM_SOA_VBS)) then
        call chem_data_read('e_ald.dat', emiss_ab(:,p_e_ald), rc, perm=perm, label='e_ald')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_ald')
        call chem_data_read('e_co.dat', emiss_ab(:,p_e_co), rc, perm=perm, label='e_co')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_co')
        call chem_data_read('e_csl.dat', emiss_ab(:,p_e_csl), rc, perm=perm, label='e_csl')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_csl')
        call chem_data_read('e_dms.dat', emiss_ab(:,p_e_dms), rc, perm=perm, label='e_dms')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_dms')
        call chem_data_read('e_eth.dat', emiss_ab(:,p_e_eth), rc, perm=perm, label='e_eth')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_eth')
        call chem_data_read('e_hc3.dat', emiss_ab(:,p_e_hc3), rc, perm=perm, label='e_hc3')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_hc3')
        call chem_data_read('e_hc5.dat', emiss_ab(:,p_e_hc5), rc, perm=perm, label='e_hc5')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_hc5')
        call chem_data_read('e_hc8.dat', emiss_ab(:,p_e_hc8), rc, perm=perm, label='e_hc8')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_hc8')
        call chem_data_read('e_hcho.dat', emiss_ab(:,p_e_hcho), rc, perm=perm, label='e_hcho')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_hcho')
        call chem_data_read('e_iso.dat', emiss_ab(:,p_e_iso), rc, perm=perm, label='e_iso')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_iso')
        call chem_data_read('e_ket.dat', emiss_ab(:,p_e_ket), rc, perm=perm, label='e_ket')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_ket')
        call chem_data_read('e_nh3.dat', emiss_ab(:,p_e_nh3), rc, perm=perm, label='e_nh3')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_nh3')
        call chem_data_read('e_no2.dat', emiss_ab(:,p_e_no2), rc, perm=perm, label='e_no2')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_no2')
        call chem_data_read('e_no.dat', emiss_ab(:,p_e_no), rc, perm=perm, label='e_no')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_no')
        call chem_data_read('e_oli.dat', emiss_ab(:,p_e_oli), rc, perm=perm, label='e_oli')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_oli')
        call chem_data_read('e_olt.dat', emiss_ab(:,p_e_oli), rc, perm=perm, label='e_olt')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_olt')
        call chem_data_read('e_ora2.dat', emiss_ab(:,p_e_ora2), rc, perm=perm, label='e_ora2')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_ora2')
        call chem_data_read('e_tol.dat', emiss_ab(:,p_e_tol), rc, perm=perm, label='e_tol')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_tol')
        call chem_data_read('e_xyl.dat', emiss_ab(:,p_e_xyl), rc, perm=perm, label='e_xyl')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading e_xyl')
      end if
    end if

    if (chem_config % biomass_burn_opt > 0) then
      call chem_data_read('ebu_oc.dat', emiss_abu(:,p_e_oc), rc, perm=perm, label='ebu_oc')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_oc')
      call chem_data_read('ebu_bc.dat', emiss_abu(:,p_e_bc), rc, perm=perm, label='ebu_bc')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_bc')
      call chem_data_read('ebu_so2.dat', emiss_abu(:,p_e_so2), rc, perm=perm, label='ebu_so2')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_so2')
      call chem_data_read('ebu_sulf.dat', emiss_abu(:,p_e_sulf), rc, perm=perm, label='ebu_sulf')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_sulf')
      call chem_data_read('ebu_pm25.dat', emiss_abu(:,p_e_pm_25), rc, perm=perm, label='ebu_pm25')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_pm25')
      call chem_data_read('ebu_pm10.dat', emiss_abu(:,p_e_pm_10), rc, perm=perm, label='ebu_pm10')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_pm10')
      call chem_data_read('plumestuff.dat', plumestuff, 1, 8, rc, perm=perm, label='plumestuff')
      if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading plumestuff')

      if ((chem_config % chem_opt == CHEM_OPT_GOCART_RACM) .or. &
          (chem_config % chem_opt == CHEM_OPT_RACM_SOA_VBS)) then
        call chem_data_read('ebu_ald.dat', emiss_abu(:,p_e_ald), rc, perm=perm, label='ebu_ald')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_ald')
        call chem_data_read('ebu_co.dat', emiss_abu(:,p_e_co), rc, perm=perm, label='ebu_co')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_co')
        call chem_data_read('ebu_csl.dat', emiss_abu(:,p_e_csl), rc, perm=perm, label='ebu_csl')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_csl')
        call chem_data_read('ebu_dms.dat', emiss_abu(:,p_e_dms), rc, perm=perm, label='ebu_dms')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_dms')
        call chem_data_read('ebu_eth.dat', emiss_abu(:,p_e_eth), rc, perm=perm, label='ebu_eth')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_eth')
        call chem_data_read('ebu_hc3.dat', emiss_abu(:,p_e_hc3), rc, perm=perm, label='ebu_hc3')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_hc3')
        call chem_data_read('ebu_hc5.dat', emiss_abu(:,p_e_hc5), rc, perm=perm, label='ebu_hc5')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_hc5')
        call chem_data_read('ebu_hc8.dat', emiss_abu(:,p_e_hc8), rc, perm=perm, label='ebu_hc8')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_hc8')
        call chem_data_read('ebu_hcho.dat', emiss_abu(:,p_e_hcho), rc, perm=perm, label='ebu_hcho')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_hcho')
        call chem_data_read('ebu_iso.dat', emiss_abu(:,p_e_iso), rc, perm=perm, label='ebu_iso')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_iso')
        call chem_data_read('ebu_ket.dat', emiss_abu(:,p_e_ket), rc, perm=perm, label='ebu_ket')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_ket')
        call chem_data_read('ebu_nh3.dat', emiss_abu(:,p_e_nh3), rc, perm=perm, label='ebu_nh3')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_nh3')
        call chem_data_read('ebu_no2.dat', emiss_abu(:,p_e_no2), rc, perm=perm, label='ebu+no2')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_no2')
        call chem_data_read('ebu_no.dat', emiss_abu(:,p_e_no), rc, perm=perm, label='ebu_no')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_no')
        call chem_data_read('ebu_oli.dat', emiss_abu(:,p_e_oli), rc, perm=perm, label='ebu_oli')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_oli')
        call chem_data_read('ebu_olt.dat', emiss_abu(:,p_e_olt), rc, perm=perm, label='ebu_olt')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_olt')
        call chem_data_read('ebu_ora2.dat', emiss_abu(:,p_e_ora2), rc, perm=perm, label='ebu_ora2')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_ora2')
        call chem_data_read('ebu_tol.dat', emiss_abu(:,p_e_tol), rc, perm=perm, label='ebu_tol')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_tol')
        call chem_data_read('ebu_xyl.dat', emiss_abu(:,p_e_xyl), rc, perm=perm, label='ebu_xyl')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading ebu_xyl')
      end if
    end if

    ! -- volcanic stuff

    if ((chem_config % chem_opt == CHEM_OPT_GOCART) .or. &
        (chem_config % chem_opt == 316) .or. &
        (chem_config % chem_opt == 317) .or. &
        (chem_config % chem_opt == 502)) then

      if (ash_mass >- -900.) then
        call chem_data_read('volcanic.dat', nv_g, rc, label='nv_g')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading nv_g')
        call chem_data_read('volcanic.dat', emiss_ash_mass, 4, rc, perm=perm, label='emiss_ash_mass')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading emiss_ash_mass')
        call chem_data_read('volcanic.dat', emiss_ash_height, 5, rc, perm=perm, label='emiss_ash_height')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading emiss_ash_height')
        call chem_data_read('volcanic.dat', emiss_ash_dt, 6, rc, perm=perm, label='emiss_ash_dt')
        if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading emiss_ash_dt')

      end if
      ! -- overwrite ash_mass if nameist value exists
      if(ash_mass > -100.)then
        write(0,*)'using namelist value for ash_mass'
        where(emiss_ash_mass > 0.) emiss_ash_mass = ash_mass
      endif
      ! -- overwrite ash_height if nameist value exists
      if(ash_height > 0.)then
        write(0,*)'using namelist value for ash_height'
        where(emiss_ash_height > 0.) emiss_ash_height = ash_height
        where(emiss_ash_height < 1.) emiss_ash_dt     = 0.
      else if (ash_height < -990.) then
        write(0,*)'resetting all ash variables to zero'
        emiss_ash_mass   = 0.
        emiss_ash_height = 0.
        emiss_ash_dt     = 0.
      endif
    end if

    ! --  Initialize chem arrays
    select case (chem_config % chem_opt)
      case(501, 502)
        tr3d = 0.
      case(500)
        tr3d = 390.
      case default
        tr3d = 1.e-16
    end select

    pm25 = 0.
    p10  = 0.
    rcav = 0.

    chem_config % call_chemistry = max(1,numphr*(int(chem_config % Chemdt+.01)*60)/3600)
    chem_config % call_biomass   = max(1,numphr*(int(chem_config % PLUMERISEFIRE_FRQ+.01)*60)/3600)

    if (chem_config % chem_in_opt == 1) then
      if ((chem_config % chem_opt == 316) .or. (chem_config % chem_opt == 317)) then
        if (ash_mass /= 0.) then
          do p = 1, 4
            pos = nbegin + p_v(p)
            call chem_data_read(trim(in_file(p)), tr3d(:,:,pos), 2, rc, label=trim(in_file(p)))
            if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading '//trim(in_file(p)))
            trdp(:,:,pos) = tr3d(:,:,pos) * dp3d
          end do
        end if
      end if
      if ((chem_config % chem_opt == 316) .or. (chem_config % chem_opt == 317)) then
        do p = 5, 10
          pos = nbegin + p_v(p)
          call chem_data_read(trim(in_file(p)), tr3d(:,:,pos), 2, rc, label=trim(in_file(p)))
          if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading '//trim(in_file(p)))
          trdp(:,:,pos) = tr3d(:,:,pos) * dp3d
        end do
      end if
      ! -- GOCART options
      if ((chem_config % chem_opt == CHEM_OPT_RACM_SOA_VBS) .or. &
          ((chem_config % chem_opt >= 300 ) .and. (chem_config % chem_opt < 500))) then
        do p = 11, 12
          pos = nbegin + p_v(p)
          call chem_data_read(trim(in_file(p)), tr3d(:,:,pos), 2, rc, label=trim(in_file(p)))
          if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading '//trim(in_file(p)))
          trdp(:,:,pos) = tr3d(:,:,pos) * dp3d
        end do
        if ((chem_config % chem_opt >= 300 ) .and. (chem_config % chem_opt < 500)) then
          do p = 13, 25
            pos = nbegin + p_v(p)
            call chem_data_read(trim(in_file(p)), tr3d(:,:,pos), 2, rc, label=trim(in_file(p)))
            if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading '//trim(in_file(p)))
            trdp(:,:,pos) = tr3d(:,:,pos) * dp3d
          end do
          if (chem_config % seas_opt == 1) then
            do p = 26, 29
              pos = nbegin + p_v(p)
              call chem_data_read(trim(in_file(p)), tr3d(:,:,pos), 2, rc, label=trim(in_file(p)))
              if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading '//trim(in_file(p)))
              trdp(:,:,pos) = tr3d(:,:,pos) * dp3d
            end do
          end if
          ! -- add gas phase chemistry input
          if (chem_config % chem_opt == 301) then
            do p = 1, 49
              pos = nbegin + p
              call chem_data_read(trim(chem_301(p))//'.in', tr3d(:,:,pos), 2, rc, label=trim(in_file(p)))
              if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading '//trim(in_file(p)))
              trdp(:,:,pos) = tr3d(:,:,pos) * dp3d
            end do
          end if
        end if
        if (chem_config % chem_opt == 108) then
          do p = 1, 103
            pos = nbegin + p
            call chem_data_read(trim(chem_108(p))//'.in', tr3d(:,:,pos), 2, rc, label=trim(in_file(p)))
            if (rc /= 0) call mpp_error(FATAL, 'chem_fields_read: error reading '//trim(in_file(p)))
            trdp(:,:,pos) = tr3d(:,:,pos) * dp3d
          end do
        end if
      end if
    end if

  end subroutine chem_background_init

  subroutine chem_data_set(area_in, deg_lon_in, deg_lat_in, &
                           u10m_in, v10m_in, rn2d_in, rc2d_in, &
                           tr3d_in, us3d_in, vs3d_in, ws3d_in, tk3d_in, pr3d_in, ph3d_in, &
                           slmsk2d_in, zorl2d_in, sm3d_in, snwdph2d_in, &
                           vtype2d_in, stype2d_in, vfrac2d_in, &
                           ts2d_in, us2d_in, hf2d_in, rsds_in, &
                           rcav_in, rnav_in, pb2d_in, exch_in, &
                           aerwrf_in, gd_cloud_in, gd_cldfr_in, &
                           perm_in, &
                           isc, iec, jsc, jec)

    integer, intent(in) :: isc, iec, jsc, jec
    real, dimension(jsc:jec),             intent(in) :: area_in, deg_lon_in, deg_lat_in, &
                                                        u10m_in, v10m_in, rn2d_in, rc2d_in, &
                                                        slmsk2d_in, zorl2d_in, snwdph2d_in, &
                                                        vtype2d_in, stype2d_in, vfrac2d_in, &
                                                        ts2d_in, us2d_in, hf2d_in, rsds_in, &
                                                        rcav_in, rnav_in, pb2d_in
    real, dimension(isc:iec, jsc:jec, ntra+ntrb), intent(in) :: tr3d_in
    real, dimension(isc:iec,jsc:jec),     intent(in) :: us3d_in, vs3d_in, ws3d_in, tk3d_in, exch_in
    real, dimension(isc:iec+1,jsc:jec),   intent(in) :: pr3d_in, ph3d_in
    real, dimension(4,jsc:jec),           intent(in) :: sm3d_in
    real, dimension(1,isc:iec+1,jsc:jec), intent(in) :: aerwrf_in, gd_cloud_in, gd_cldfr_in
    integer, dimension(nip), intent(in) :: perm_in

    ! -- begin
    area     = area_in
    deg_lon  = deg_lon_in
    deg_lat  = deg_lat_in
!   u10m     = u10m_in
!   v10m     = v10m_in
    rn2d     = rn2d_in
    rc2d     = rc2d_in
    tr3d     = tr3d_in
    us3d     = us3d_in
    vs3d     = vs3d_in
    ws3d     = ws3d_in
    tk3d     = tk3d_in
    pr3d     = pr3d_in
    ph3d     = ph3d_in
    slmsk2d  = slmsk2d_in
    zorl2d   = zorl2d_in
    sm3d     = sm3d_in
    snwdph2d = snwdph2d_in
    vtype2d  = vtype2d_in
    stype2d  = stype2d_in
    vfrac2d  = vfrac2d_in
    ts2d     = ts2d_in
    us2d     = us2d_in
    hf2d     = hf2d_in
    rsds     = rsds_in
    rcav     = rcav_in
    rnav     = rnav_in
    pb2d     = pb2d_in
    exch     = exch_in
    aerwrf   = aerwrf_in
    gd_cloud = gd_cloud_in
    gd_cldfr = gd_cldfr_in

    perm     = perm_in

  end subroutine chem_data_set

! -------------------------------------------------------------------------------

  subroutine chem_data_prep

    ! -- compute additional arrays
    dp3d = pr3d(1:nvl, :) - pr3d(2:, :)

  end subroutine chem_data_prep

! -------------------------------------------------------------------------------

  subroutine chem_history_write(config, its)

    type(chem_config_type), intent(in) :: config
    integer,                intent(in) :: its

    ! -- local variables
    character (len=4), parameter :: chem_names301(49) = (/ &
                          'pso2','sulf','pno2','ppno','ppo3',        &
                          'hno3','h2o2','pald','hcho','pop1','pop2',  &
                          'ppaa','ora1','ora2','pnh3','n2o5','pno3',  &
                          'ppan','phc3','phc5','phc8','peth','ppco',  &
                          'pete','polt','poli','ptol','pxyl','aco3',  &
                          'tpan','hono','hno4','pket','pgly','mgly',  &
                          'pdcb','onit','pcsl','piso','pco2','pch4',  &
                          'pudd','hket','papi','plim','dien','macr',  &
                          'ppho','pho2'/)
      character (len=4), parameter :: chem_names108(103) = (/ &
                          'pso2','sulf','pno2','ppno','ppo3',         &
                          'hno3','h2o2','pald','hcho','pop1','pop2',  &
                          'ppaa','ora1','ora2','pnh3','n2o5','pno3',  &
                          'ppan','phc3','phc5','phc8','peth','ppco',  &
                          'pete','polt','poli','ptol','pxyl','aco3',  &
                          'tpan','hono','hno4','pket','pgly','mgly',  &
                          'pdcb','onit','pcsl','piso','pco2','pch4',  &
                          'pudd','hket','papi','plim','dien','macr',  &
                          'hace','ishp','ison','mahp','mpan','nald',  &
                          'sesq','pmbo','cva1','cva2','cva3','cva4',  &
                          'cvb1','cvb2','cvb3','cvb4','ppho','pho2',  &
                          'soaj','soai','nhaj','nhai','n3aj','n3ai',  &
                          'naaj','naai','claj','clai','as1j','as1i',  &
                          'as2j','as2i','as3j','as3i','as4j','as4i',  &
                          'bs1j','bs1i','bs2j','bs2i','bs3j','bs3i',  &
                          'bs4j','bs4i','opaj','opai','pecj','peci',  &
                          'p25j','p25i','atha','seas','sila','pnu0',  &
                          'pac0','corn'/)

    integer :: ArchvStep, ichem_start, imoist_start, j, k, nv
    real    :: dpsum
    logical, parameter:: flux_avg=.false.        ! T: write fluxes averaged over "ArchvIntvl" !lzhang
    integer, parameter :: mp_physics = 0
    integer, parameter :: cu_physics = 0
    real :: exttsum
    real, dimension(jts:jte)      :: intaer, intash, intbc, intdust, intoc, intsulf, o3dg, o3du
    real, dimension(nvl, jts:jte) :: d1st, d2st, d3st, d4st, d5st, dms1, rho_phys, &
                                     sea1, sea2, sea3, sea4, trco

    integer :: item

    ! -- begin
    print *,'--- chem_history_write: entering ...'

    ichem_start = ntra

    if (config % chem_opt >= 300) then
       d1st(:,:) = tr3d(:,:,ichem_start+p_dust_1)
       d2st(:,:) = tr3d(:,:,ichem_start+p_dust_2)
       d3st(:,:) = tr3d(:,:,ichem_start+p_dust_3)
       d4st(:,:) = tr3d(:,:,ichem_start+p_dust_4)
       d5st(:,:) = tr3d(:,:,ichem_start+p_dust_5)
       if (its .eq. 0) then
           d1st_ave=0.
           d2st_ave=0.
           d3st_ave=0.
           d4st_ave=0.
           d5st_ave=0.
        endif
    endif ! chem_opt>=300

    emid1_ave=0.
    emid2_ave=0.
    emid3_ave=0.
    emid4_ave=0.
    emid5_ave=0.
    aod2d_ave=0.

    if (its > 0) then
#if 1
      call average(pr3d, tk3d)
#else
      do j = jts, jte
        exttsum=0.
        do k = 1, nvl
          rho_phys(k,j)=.5*(pr3d(k,j)+pr3d(k+1,j))/(RD*tk3d(k,j))
!*(1.+.608*qv3d(k,j))
          if (config % chem_opt >= 300) then
            d1st_ave(k,j) = d1st_ave(k,j)+d1st(k,j)
            d2st_ave(k,j) = d2st_ave(k,j)+d2st(k,j)
            d3st_ave(k,j) = d3st_ave(k,j)+d3st(k,j)
            d4st_ave(k,j) = d4st_ave(k,j)+d4st(k,j)
            d5st_ave(k,j) = d5st_ave(k,j)+d5st(k,j)
          endif ! chem_opt>=300
          if (config % aer_ra_feedback == 0) then
              exttsum=exttsum+ext_cof(k,j,p_extcof55)
          endif
        enddo
        if (config % aer_ra_feedback == 0 )then
          aod2d(j)=exttsum
        endif
        emid1_ave(j) = emid1_ave(j)+emi_d1(j)
        emid2_ave(j) = emid2_ave(j)+emi_d2(j)
        emid3_ave(j) = emid3_ave(j)+emi_d3(j)
        emid4_ave(j) = emid4_ave(j)+emi_d4(j)
        emid5_ave(j) = emid5_ave(j)+emi_d5(j)
        aod2d_ave(j) = aod2d_ave(j)+aod2d(j)
      enddo
#endif
    end if

    ! -- check if it is time to write history
    
    ArchvStep = config % archive_step

    if (mod(its, ArchvStep) == 0) then
      ichem_start = ntra
      if (config % chem_opt >= 300 .and. config % chem_opt < 500) then
        dms1(:,:) = tr3d(:,:,ichem_start+p_dms)
        sea1(:,:) = tr3d(:,:,ichem_start+p_seas_1)!*rho_phys(:,:)
        sea2(:,:) = tr3d(:,:,ichem_start+p_seas_2)!*rho_phys(:,:)
        sea3(:,:) = tr3d(:,:,ichem_start+p_seas_3)!*rho_phys(:,:)
        sea4(:,:) = tr3d(:,:,ichem_start+p_seas_4)!*rho_phys(:,:)
      end if !chem_opt >= 300 .and. chem_opt < 500
!
! the floowing requires some type of GOCART or volcanic ash option
! will only work for gocart type mass variables (such as p_bc, .....)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (((config % chem_opt >= 300) .and. (config % chem_opt < 500)) &
          .or. (config % chem_opt == 108)) then
#if 1
        call compute_tracers(pr3d, tk3d, ph3d, tr3d)
#else
        do j=jts, jte
          dpsum=0.
          intash(j)=0.
          intaer(j)=0.
          intbc(j)=0.
          intoc(j)=0.
          intsulf(j)=0.
          intdust(j)=0.
          o3du(j)=0.
          o3dg(j)=0.
          do k = 1, nvl
            dpsum=dpsum+(pr3d(k,j)-pr3d(k+1,j))
            rho_phys(k,j)=.5*(pr3d(k,j)+pr3d(k+1,j))&
              /(RD*tk3d(k,j)) !*(1.+.608*qv3d(k,j))
            if ((config % chem_opt >= 300) .and. (config % chem_opt < 500)) then
              intaer(j)=intaer(j)+1e-6*tr3d(k,j,ichem_start+p_p25)*(pr3d(k,j)-pr3d(k+1,j))/g
              intbc(j)=intbc(j)+(tr3d(k,j,ichem_start+p_bc1)&
                +tr3d(k,j,ichem_start+p_bc2))*1e-6*(pr3d(k,j)-pr3d(k+1,j))/g
              intoc(j)=intoc(j)+(tr3d(k,j,ichem_start+p_oc1)&
                +tr3d(k,j,ichem_start+p_oc2))*1e-6*(pr3d(k,j)-pr3d(k+1,j))/g
              intdust(j)=intdust(j)+(d1st(k,j)&
                +.286*d2st(k,j))*1e-6*(pr3d(k,j)-pr3d(k+1,j))/g
              intsulf(j)=intsulf(j)+1e-6*tr3d(k,j,ichem_start+p_sulf)*(pr3d(k,j)&
                -pr3d(k+1,j))/g
            end if !chem_opt >= 300 .and. chem_opt < 500

            if (config % chem_opt == 108) then
              intaer(j)=intaer(j)+(tr3d(k,j,ichem_start+p_p25j)+&
                tr3d(k,j,ichem_start+p_p25i))*1e-6*(pr3d(k,j)-pr3d(k+1,j))/g
              intbc(j)=intbc(j)+(tr3d(k,j,ichem_start+p_ecj)&
                +tr3d(k,j,ichem_start+p_eci))*1e-6*(pr3d(k,j)-pr3d(k+1,j))/g
              intoc(j)=intoc(j)+(tr3d(k,j,ichem_start+p_orgpaj)&
                +tr3d(k,j,ichem_start+p_orgpai))*1e-6*(pr3d(k,j)-pr3d(k+1,j))/g
              intdust(j)=intdust(j)+(tr3d(k,j,ichem_start+p_soila)&
                )*1e-6*(pr3d(k,j)-pr3d(k+1,j))/g
              intsulf(j)=intsulf(j)+(tr3d(k,j,ichem_start+p_so4aj)+&
                tr3d(k,j,ichem_start+p_so4ai))*1e-6*(pr3d(k,j)&
                -pr3d(k+1,j))/g
            end if !chem_opt == 108

            if ((config % chem_opt == 301) .or. (config % chem_opt == 108))then
              o3du(j)=o3du(j)+(0.001*tr3d(k,j,ichem_start+p_o3)*(ph3d(k+1,j)-ph3d(k,j))&
                *rho_phys(k,j)*6.022*1e23/(48.*9.8*2.69*1e20))
            end if
            o3dg(j)=o3dg(j)+(1e6*0.001*tr3d(k,j,4)*airmw/48.*(ph3d(k+1,j)-ph3d(k,j))&
              *rho_phys(k,j)*6.022*1e23/(48.*9.8*2.69*1e20))

          if (config % chem_opt==304.or.config % chem_opt==316.or.config % chem_opt==317) then
              intdust(j)=intdust(j)+d1st(k,j)*1e-6*(pr3d(k,j)-pr3d(k+1,j))&
                /g
            endif
!
            if (config % chem_opt == 316) then
              intash(j)=intash(j)+(tr3d(k,j,ichem_start+p_vash_1)       &
                + tr3d(k,j,ichem_start+p_vash_2)       &
                + tr3d(k,j,ichem_start+p_vash_3)       &
                + tr3d(k,j,ichem_start+p_vash_4)       &
                + tr3d(k,j,ichem_start+p_vash_5)       &
                + tr3d(k,j,ichem_start+p_vash_6)       &
                + tr3d(k,j,ichem_start+p_vash_7)       &
                + tr3d(k,j,ichem_start+p_vash_8)       &
                + tr3d(k,j,ichem_start+p_vash_9)       &
                + tr3d(k,j,ichem_start+p_vash_10))     &
                *1e-6*(pr3d(k,j)-pr3d(k+1,j))/g
            endif
            if (config % chem_opt == 317 ) then
              intash(j)=intash(j)+(tr3d(k,j,ichem_start+p_vash_1)       &
                + tr3d(k,j,ichem_start+p_vash_2)       &
                + tr3d(k,j,ichem_start+p_vash_3)       &
                + tr3d(k,j,ichem_start+p_vash_4))      &
                *1e-6*(pr3d(k,j)-pr3d(k+1,j))/g
            endif

          end do
          if (config % chem_opt == 316 .or. config % chem_opt == 317 ) intash(j)=intash(j)
          intaer(j)=intaer(j)
          intbc(j)=intbc(j)
          intoc(j)=intoc(j)
          intsulf(j)=intsulf(j)
          intdust(j)=intdust(j)
        end do
#endif
      end if ! chem_opt >= 300 .and. chem_opt < 500 .or. chem_opt=108

      ! --- accumulate fields averaged over "ArchvIntvl"
      do j = jts, jte
        do k = 1, nvl
          if (config % chem_opt >= 300) then
            d1st_ave(k,j) = d1st_ave(k,j)/float(ArchvStep)
            d2st_ave(k,j) = d2st_ave(k,j)/float(ArchvStep)
            d3st_ave(k,j) = d3st_ave(k,j)/float(ArchvStep)
            d4st_ave(k,j) = d4st_ave(k,j)/float(ArchvStep)
            d5st_ave(k,j) = d5st_ave(k,j)/float(ArchvStep)
          end if
        end do
        aod2d_ave(j) = aod2d_ave(j)/float(ArchvStep)
      end do

      if ((mp_physics.ne.0).or.(cu_physics.ne.0)) then
        call maybe_write(its,'trco',trco,nvl)
      end if

      if ((config % chem_opt >= 300) .and. (config % chem_opt < 500)) then
        call maybe_write(its,'ex3D',exch,nvl)
        call maybe_write(its,'pm25',pm25,nvl)
        call maybe_write(its,'pm10',p10,nvl)
        call maybe_write(its,'dms1',dms1,nvl)
        call maybe_write(its,'s1ea',sea1,nvl)
        call maybe_write(its,'s2ea',sea2,nvl)
        call maybe_write(its,'s3ea',sea3,nvl)
        call maybe_write(its,'s4ea',sea4,nvl)
        if (flux_avg) then
          call maybe_write(its,'d1st',d1st_ave,nvl)
          call maybe_write(its,'d2st',d2st_ave,nvl)
          call maybe_write(its,'d3st',d3st_ave,nvl)
          call maybe_write(its,'d4st',d4st_ave,nvl)
          call maybe_write(its,'d5st',d5st_ave,nvl)
          call maybe_write(its,'emd1',emid1_ave, 1, twodfile=.true.)
          call maybe_write(its,'emd2',emid2_ave, 1, twodfile=.true.)
          call maybe_write(its,'emd3',emid3_ave, 1, twodfile=.true.)
          call maybe_write(its,'emd4',emid4_ave, 1, twodfile=.true.)
          call maybe_write(its,'emd5',emid5_ave, 1, twodfile=.true.)
          call maybe_write(its,'ao2D',aod2d_ave ,1, twodfile=.true.)
          d1st_ave=0.
          d2st_ave=0.
          d3st_ave=0.
          d4st_ave=0.
          d5st_ave=0.
          emid1_ave=0.
          emid2_ave=0.
          emid3_ave=0.
          emid4_ave=0.
          emid5_ave=0.
          aod2d_ave=0.
        else
          call maybe_write(its,'d1st',d1st,nvl)
          call maybe_write(its,'d2st',d2st,nvl)
          call maybe_write(its,'d3st',d3st,nvl)
          call maybe_write(its,'d4st',d4st,nvl)
          call maybe_write(its,'d5st',d5st,nvl)

          ! -- output dust emis and source
          call maybe_write(its,'emd1',emi_d1, 1, twodfile=.true.)
          call maybe_write(its,'emd2',emi_d2, 1, twodfile=.true.)
          call maybe_write(its,'emd3',emi_d3, 1, twodfile=.true.)
          call maybe_write(its,'emd4',emi_d4, 1, twodfile=.true.)
          call maybe_write(its,'emd5',emi_d5, 1, twodfile=.true.)
          call maybe_write(its,'ao2D',aod2d,1, twodfile=.true.)
        end if
        call maybe_write(its,'ero1',ero1, 1, twodfile=.true.)
        call maybe_write(its,'ero2',ero2, 1, twodfile=.true.)
        call maybe_write(its,'ero3',ero3, 1, twodfile=.true.)
        call maybe_write(its,'dms0',dm0, 1, twodfile=.true.)

        ! -- output aeroosol wetdepostion
        !call maybe_write(its,'wbc1',wet_dep(:,p_bc1), 1, twodfile=.true.)
        call maybe_write(its,'wbc2',wet_dep(:,p_bc2), 1, twodfile=.true.)
        !call maybe_write(its,'woc1',wet_dep(:,p_oc1), 1, twodfile=.true.)
        call maybe_write(its,'woc2',wet_dep(:,p_oc2), 1, twodfile=.true.)
        call maybe_write(its,'wp25',wet_dep(:,p_p25), 1, twodfile=.true.)
        call maybe_write(its,'wp10',wet_dep(:,p_p10), 1, twodfile=.true.)
        call maybe_write(its,'wso4',wet_dep(:,p_sulf), 1, twodfile=.true.)
        call maybe_write(its,'wdt1',wet_dep(:,p_dust_1), 1, twodfile=.true.)
        call maybe_write(its,'wdt2',wet_dep(:,p_dust_2), 1, twodfile=.true.)
        call maybe_write(its,'wdt3',wet_dep(:,p_dust_3), 1, twodfile=.true.)
        call maybe_write(its,'wdt4',wet_dep(:,p_dust_4), 1, twodfile=.true.)
        call maybe_write(its,'wdt5',wet_dep(:,p_dust_5), 1, twodfile=.true.)
        call maybe_write(its,'wse1',wet_dep(:,p_seas_1), 1, twodfile=.true.)
        call maybe_write(its,'wse2',wet_dep(:,p_seas_2), 1, twodfile=.true.)
        call maybe_write(its,'wse3',wet_dep(:,p_seas_3), 1, twodfile=.true.)
        call maybe_write(its,'wse4',wet_dep(:,p_seas_4), 1, twodfile=.true.)
        call maybe_write(its,'aiso',emiss_ab1(:,p_e_iso), 1, twodfile=.true.)
        call maybe_write(its,'aso2',emiss_ab1(:,p_e_so2), 1, twodfile=.true.)
        call maybe_write(its,'ano2',emiss_ab1(:,p_e_no2), 1, twodfile=.true.)
        call maybe_write(its,'fiso',emiss_abu(:,p_e_iso), 1, twodfile=.true.)
!     call maybe_write(its,'aeth',emiss_ab(:,p_e_eth), 1, twodfile=.true.)
!     call maybe_write(its,'feth',emiss_abu(:,p_e_eth), 1, twodfile=.true.)
        if (config % chem_opt == 300) then
          call maybe_write(its,'ohbg',oh_bg, nvl)
          call maybe_write(its,'hobg',h2o2_bg, nvl)
          call maybe_write(its,'no3b',no3_bg, nvl)
          call maybe_write(its,'ocbb',ebu_oc, nvl)
        end if
        ! -- output air density
        call maybe_write(its,'ts2D',ts2d, 1, twodfile=.true.)
        call maybe_write(its,'aird',rho_phys, nvl)
        call maybe_write(its,'extt',ext_cof,nvl)! out_put ext_cof  

        ! -- output GOCART aerosols variables
        !
        !  change unit from ug/kg to ug/m3, need to multiply rho_phys 
        !
        dms1(:,:) = tr3d(:,:,ichem_start+p_bc1)!*rho_phys(:,:)
        call maybe_write(its,'pbc1',dms1,nvl)
        dms1(:,:) = tr3d(:,:,ichem_start+p_bc2)!*rho_phys(:,:)
        call maybe_write(its,'pbc2',dms1,nvl)
        dms1(:,:) = tr3d(:,:,ichem_start+p_oc1)!*rho_phys(:,:)
        call maybe_write(its,'obc1',dms1,nvl)
        dms1(:,:) = tr3d(:,:,ichem_start+p_oc2)!*rho_phys(:,:)
        call maybe_write(its,'obc2',dms1,nvl)
        dms1(:,:) = tr3d(:,:,ichem_start+p_sulf)

        if (config % chem_opt == 300) call maybe_write(its,'sulf',dms1,nvl)
        dms1(:,:) = tr3d(:,:,ichem_start+p_so2)

        if (config % chem_opt == 300) call maybe_write(its,'pso2',dms1,nvl)
        dms1(:,:) = tr3d(:,:,ichem_start+p_msa)
        call maybe_write(its,'pmsa',dms1,nvl)
        dms1(:,:) = tr3d(:,:,ichem_start+p_p25)!*rho_phys(:,:)
        call maybe_write(its,'pp25',dms1,nvl)
        dms1(:,:) = tr3d(:,:,ichem_start+p_p10)!*rho_phys(:,:)
        call maybe_write(its,'pp10',dms1,nvl)

        ! -- output gas phase chemistry
        if (config % chem_opt == 301) then
          call maybe_write(its,'o3du',o3du,1, twodfile=.true.)
          call maybe_write(its,'o3dg',o3dg,1, twodfile=.true.)
          do nv = 1, 49
            dms1(:,:) = tr3d(:,:,ichem_start+nv)
            call maybe_write(its,chem_names301(nv),dms1,nvl)
          end do
        end if ! chem_opt=301
     
          if (config % chem_opt.eq.316.or.config % chem_opt.eq.317) then
            print *,'p_vash_1,p_vash_4 = ',p_vash_1,p_vash_4
            dms1(:,:) = tr3d(:,:,ichem_start+p_vash_1)
            call maybe_write(its,'ash1',dms1,nvl)
            dms1(:,:) = tr3d(:,:,ichem_start+p_vash_2)
            call maybe_write(its,'ash2',dms1,nvl)
            dms1(:,:) = tr3d(:,:,ichem_start+p_vash_3)
            call maybe_write(its,'ash3',dms1,nvl)
            dms1(:,:) = tr3d(:,:,ichem_start+p_vash_4)
            call maybe_write(its,'ash4',dms1,nvl)
            if (config % chem_opt.eq.316) then
              dms1(:,:) = tr3d(:,:,ichem_start+p_vash_5)
              call maybe_write(its,'ash5',dms1,nvl)
              dms1(:,:) = tr3d(:,:,ichem_start+p_vash_6)
              call maybe_write(its,'ash6',dms1,nvl)
              dms1(:,:) = tr3d(:,:,ichem_start+p_vash_7)
              call maybe_write(its,'ash7',dms1,nvl)
              dms1(:,:) = tr3d(:,:,ichem_start+p_vash_8)
              call maybe_write(its,'ash8',dms1,nvl)
              dms1(:,:) = tr3d(:,:,ichem_start+p_vash_9)
              call maybe_write(its,'ash9',dms1,nvl)
              dms1(:,:) = tr3d(:,:,ichem_start+p_vash_10)
              call maybe_write(its,'ash0',dms1,nvl)
            endif !chem_opt=316
          endif !chem_opt=316 or chem_opt=317

        call maybe_write(its,'ia2D',intaer,1, twodfile=.true.)
        call maybe_write(its,'ib2D',intbc,1, twodfile=.true.)
        call maybe_write(its,'io2D',intoc,1, twodfile=.true.)
        call maybe_write(its,'is2D',intsulf,1, twodfile=.true.)
        call maybe_write(its,'id2D',intdust,1, twodfile=.true.)
      end if !chem_opt.ge.300 .and. chem_opt.lt.500

      ! -- output racm soa vbs chemistry
      if (config % chem_opt == 108) then
        if (flux_avg) then
          call maybe_write(its,'ao2D',aod2d_ave ,1, twodfile=.true.)
          aod2d_ave=0.0
        else
          call maybe_write(its,'ao2D',aod2d,1, twodfile=.true.)
        end if
        call maybe_write(its,'ia2D',intaer,1, twodfile=.true.)
        call maybe_write(its,'ib2D',intbc,1, twodfile=.true.)
        call maybe_write(its,'io2D',intoc,1, twodfile=.true.)
        call maybe_write(its,'is2D',intsulf,1, twodfile=.true.)
        call maybe_write(its,'id2D',intdust,1, twodfile=.true.)
        call maybe_write(its,'aiso',emiss_ab1(:,p_e_iso), 1,twodfile=.true.)
        call maybe_write(its,'aso2',emiss_ab1(:,p_e_so2), 1,twodfile=.true.)
        call maybe_write(its,'ano2',emiss_ab1(:,p_e_no2), 1,twodfile=.true.)
        call maybe_write(its,'fiso',emiss_abu(:,p_e_iso), 1,twodfile=.true.)
        call maybe_write(its,'o3du',o3du,1, twodfile=.true.)
        call maybe_write(its,'o3dg',o3dg,1, twodfile=.true.)
        call maybe_write(its,'pm25',pm25,nvl)
        call maybe_write(its,'pm10',p10,nvl)
        call maybe_write(its,'aird',rho_phys, nvl)
        call maybe_write(its,'extt',ext_cof,nvl)! out_put ext_cof  
        do nv=1,103
          dms1(:,:) = tr3d(:,:,ichem_start+nv)
          call maybe_write(its,chem_names108(nv),dms1,nvl)
        enddo
      end if !chem_opt==108
    end if

    ! -- write meteorological fields
    call maybe_write(its,'pr3D',pr3d, nvlp1)
    call maybe_write(its,'ph3D',ph3d, nvlp1, scalefactor=1./9.8)
    call maybe_write(its,'tk3D',tk3d, nvl)
    call maybe_write(its,'ws3D',ws3d, nvl)
    call maybe_write(its,'us3D',us3d, nvl)
    call maybe_write(its,'vs3D',vs3d, nvl)

    call maybe_write(its,'dp3D',dp3d, nvl)
    call maybe_write(its,'th3D',tr3d(:,:,1),nvl)
    call maybe_write(its,'qv3D',tr3d(:,:,2),nvl, scalefactor=1000.)
    call maybe_write(its,'qw3D',tr3d(:,:,3),nvl, scalefactor=1000.)
    call maybe_write(its,'oz3D',tr3d(:,:,4),nvl, scalefactor=1000.)

    print *,'--- chem_history_write: exiting ...'
    
  contains

    subroutine average(pr3d, tk3d)

      real, intent(in) :: pr3d(nvlp1,jms:jme)    ! pressure (pascal)
      real, intent(in) :: tk3d(nvl,jms:jme)      ! temperature, kelvin
      ! -- local variables
      integer :: j, k
      real    :: exttsum
!     real, dimension(nvl, jts:jte) :: rho_phys

      do j = jts, jte
        exttsum=0.
        do k = 1, nvl
          rho_phys(k,j)=.5*(pr3d(k,j)+pr3d(k+1,j))/(RD*tk3d(k,j))
!*(1.+.608*qv3d(k,j))
          if (config % chem_opt >= 300) then
            d1st_ave(k,j) = d1st_ave(k,j)+d1st(k,j)
            d2st_ave(k,j) = d2st_ave(k,j)+d2st(k,j)
            d3st_ave(k,j) = d3st_ave(k,j)+d3st(k,j)
            d4st_ave(k,j) = d4st_ave(k,j)+d4st(k,j)
            d5st_ave(k,j) = d5st_ave(k,j)+d5st(k,j)
          endif ! chem_opt>=300
          if (config % aer_ra_feedback == 0) then
              exttsum=exttsum+ext_cof(k,j,p_extcof55)
          endif
        enddo
        if (config % aer_ra_feedback == 0 )then
          aod2d(j)=exttsum
        endif
        emid1_ave(j) = emid1_ave(j)+emi_d1(j)
        emid2_ave(j) = emid2_ave(j)+emi_d2(j)
        emid3_ave(j) = emid3_ave(j)+emi_d3(j)
        emid4_ave(j) = emid4_ave(j)+emi_d4(j)
        emid5_ave(j) = emid5_ave(j)+emi_d5(j)
        aod2d_ave(j) = aod2d_ave(j)+aod2d(j)
      enddo

    end subroutine average

    subroutine compute_tracers(pr3d, tk3d, ph3d, tr3d)

      real, intent(in) :: pr3d(nvlp1,jms:jme)    ! pressure (pascal)
      real, intent(in) :: tk3d(nvl,jms:jme)      ! temperature, kelvin
      real, intent(in) :: ph3d(nvlp1,jms:jme)    ! geopotential (=gz), m^2/s^2
      real, intent(in) :: tr3d(nvl,jms:jme,ntra+ntrb)  ! 1=pot.temp, 2=water vapor, 3=cloud water, 4=ozone

      ! -- local variables
      integer :: j, k
      real    :: dpsum

        do j=jts, jte
          dpsum=0.
          intash(j)=0.
          intaer(j)=0.
          intbc(j)=0.
          intoc(j)=0.
          intsulf(j)=0.
          intdust(j)=0.
          o3du(j)=0.
          o3dg(j)=0.
          do k = 1, nvl
            dpsum=dpsum+(pr3d(k,j)-pr3d(k+1,j))
            rho_phys(k,j)=.5*(pr3d(k,j)+pr3d(k+1,j))&
              /(RD*tk3d(k,j)) !*(1.+.608*qv3d(k,j))
            if ((config % chem_opt >= 300) .and. (config % chem_opt < 500)) then
              intaer(j)=intaer(j)+1e-6*tr3d(k,j,ichem_start+p_p25)*(pr3d(k,j)-pr3d(k+1,j))/g
              intbc(j)=intbc(j)+(tr3d(k,j,ichem_start+p_bc1)&
                +tr3d(k,j,ichem_start+p_bc2))*1e-6*(pr3d(k,j)-pr3d(k+1,j))/g
              intoc(j)=intoc(j)+(tr3d(k,j,ichem_start+p_oc1)&
                +tr3d(k,j,ichem_start+p_oc2))*1e-6*(pr3d(k,j)-pr3d(k+1,j))/g
              intdust(j)=intdust(j)+(d1st(k,j)&
                +.286*d2st(k,j))*1e-6*(pr3d(k,j)-pr3d(k+1,j))/g
              intsulf(j)=intsulf(j)+1e-6*tr3d(k,j,ichem_start+p_sulf)*(pr3d(k,j)&
                -pr3d(k+1,j))/g
            end if !chem_opt >= 300 .and. chem_opt < 500

            if (config % chem_opt == 108) then
              intaer(j)=intaer(j)+(tr3d(k,j,ichem_start+p_p25j)+&
                tr3d(k,j,ichem_start+p_p25i))*1e-6*(pr3d(k,j)-pr3d(k+1,j))/g
              intbc(j)=intbc(j)+(tr3d(k,j,ichem_start+p_ecj)&
                +tr3d(k,j,ichem_start+p_eci))*1e-6*(pr3d(k,j)-pr3d(k+1,j))/g
              intoc(j)=intoc(j)+(tr3d(k,j,ichem_start+p_orgpaj)&
                +tr3d(k,j,ichem_start+p_orgpai))*1e-6*(pr3d(k,j)-pr3d(k+1,j))/g
              intdust(j)=intdust(j)+(tr3d(k,j,ichem_start+p_soila)&
                )*1e-6*(pr3d(k,j)-pr3d(k+1,j))/g
              intsulf(j)=intsulf(j)+(tr3d(k,j,ichem_start+p_so4aj)+&
                tr3d(k,j,ichem_start+p_so4ai))*1e-6*(pr3d(k,j)&
                -pr3d(k+1,j))/g
            end if !chem_opt == 108

            if ((config % chem_opt == 301) .or. (config % chem_opt == 108))then
              o3du(j)=o3du(j)+(0.001*tr3d(k,j,ichem_start+p_o3)*(ph3d(k+1,j)-ph3d(k,j))&
                *rho_phys(k,j)*6.022*1e23/(48.*9.8*2.69*1e20))
            end if
            o3dg(j)=o3dg(j)+(1e6*0.001*tr3d(k,j,4)*airmw/48.*(ph3d(k+1,j)-ph3d(k,j))&
              *rho_phys(k,j)*6.022*1e23/(48.*9.8*2.69*1e20))

          if (config % chem_opt==304.or.config % chem_opt==316.or.config % chem_opt==317) then
              intdust(j)=intdust(j)+d1st(k,j)*1e-6*(pr3d(k,j)-pr3d(k+1,j))&
                /g
            endif
!
            if (config % chem_opt == 316) then
              intash(j)=intash(j)+(tr3d(k,j,ichem_start+p_vash_1)       &
                + tr3d(k,j,ichem_start+p_vash_2)       &
                + tr3d(k,j,ichem_start+p_vash_3)       &
                + tr3d(k,j,ichem_start+p_vash_4)       &
                + tr3d(k,j,ichem_start+p_vash_5)       &
                + tr3d(k,j,ichem_start+p_vash_6)       &
                + tr3d(k,j,ichem_start+p_vash_7)       &
                + tr3d(k,j,ichem_start+p_vash_8)       &
                + tr3d(k,j,ichem_start+p_vash_9)       &
                + tr3d(k,j,ichem_start+p_vash_10))     &
                *1e-6*(pr3d(k,j)-pr3d(k+1,j))/g
            endif
            if (config % chem_opt == 317 ) then
              intash(j)=intash(j)+(tr3d(k,j,ichem_start+p_vash_1)       &
                + tr3d(k,j,ichem_start+p_vash_2)       &
                + tr3d(k,j,ichem_start+p_vash_3)       &
                + tr3d(k,j,ichem_start+p_vash_4))      &
                *1e-6*(pr3d(k,j)-pr3d(k+1,j))/g
            endif

          end do
          if (config % chem_opt == 316 .or. config % chem_opt == 317 ) intash(j)=intash(j)
          intaer(j)=intaer(j)
          intbc(j)=intbc(j)
          intoc(j)=intoc(j)
          intsulf(j)=intsulf(j)
          intdust(j)=intdust(j)
        end do

    end subroutine compute_tracers

  end subroutine chem_history_write

! -------------------------------------------------------------------------------

  subroutine chem_prep(config,ktau,dtstep,tr3d,tk3d,sm3d,                              &
                       ts2d,us2d,rsds,pr3d,emiss_ash_mass,emiss_ash_height,          &
                       emiss_ash_dt,dm0,emiss_tr_mass,emiss_tr_height,               &
                       emiss_tr_dt,snwdph2d,VFRAC2d,VTYPE2d,STYPE2d,us3d,vs3d,ws3d,           &
                       slmsk2d,zorl2d,exch,pb2d,hf2d,th_pvsrf,oh_backgd,h2o2_backgd, &
                       no3_backgd,backg_oh,backg_h2o2,backg_no3,p_gocart,            &
                       nvl_gocart, ttday,tcosz,gmt,julday,ph3d,area,ero1,            &
                       ero2,ero3,rcav,raincv_b,deg_lat,deg_lon,nvl,nvlp1,ntra,       &
                       relhum,rri,t_phy,moist,u_phy,v_phy,p_phy,chem,tsk,ntrb,       &
                       g,rd,p1000,cp,erod,emis_ant,emis_vol,e_co,dms_0,              &
                       u10,v10,ivgtyp,isltyp,gsw,vegfra,rmol,ust,znt,xland,dxy,      &
                       t8w,p8w,exch_h,pbl,hfx,snowh,xlat,xlong,convfac,z_at_w,zmid,dz8w,vvel,&
                       rho_phy,smois,num_soil_layers,num_chem,num_moist,             &
                       emiss_abu,ebu_in,emiss_ab,num_ebu_in,num_emis_ant,            &
                       num_emis_vol,kemit,call_gocart,plumestuff,                    &
                       mean_fct_agtf,mean_fct_agef,mean_fct_agsv,                    &
                       mean_fct_aggr,firesize_agtf,firesize_agef,                    &
!                      firesize_agsv,firesize_aggr,readrestart,chem_in_opt,          &
                       firesize_agsv,firesize_aggr,                                  &
                       ids,ide, jds,jde, kds,kde,                                    &
                       ims,ime, jms,jme, kms,kme,                                    &
                       its,ite, jts,jte, kts,kte)

    IMPLICIT NONE

    ! -- input variables
    type(chem_config_type), intent(in) :: config

    INTEGER,      INTENT(IN) :: ktau,nvl,nvlp1,ntra,ntrb,nvl_gocart
    INTEGER,      INTENT(IN) :: num_ebu_in,num_soil_layers,num_chem,num_moist,julday,     &
                                   num_emis_vol,num_emis_ant,ids,ide, jds,jde, kds,kde,      &
                                   kemit,ims,ime, jms,jme, kms,kme,                          &
                                   its,ite, jts,jte, kts,kte
    LOGICAL,      INTENT(IN) :: call_gocart
    REAL,         INTENT(IN) :: g,rd,p1000,cp,dtstep,gmt
    real,      intent(inout) :: tk3d(nvl,jms:jme)      ! temperature, kelvin
    real,         intent(in) :: exch(nvl,jms:jme)      !
    real,         intent(in) :: oh_backgd(nvl_gocart,jms:jme)      !
    real,         intent(in) :: h2o2_backgd(nvl_gocart,jms:jme)      !
    real,         intent(in) :: no3_backgd(nvl_gocart,jms:jme)      !
    real,         intent(in) :: tr3d(nvl,jms:jme,ntra+ntrb)  ! 1=pot.temp, 2=water vapor, 3=cloud water, 4=ozone
    real,         intent(in) :: sm3d(4,jms:jme)        ! soil moisture
    real,         intent(in) :: ts2d(jms:jme)          ! skin temperature
    real,         intent(in) :: us2d(jms:jme)          ! friction velocity/equivalent momentum flux
    real,         intent(in) :: pb2d(jms:jme)          ! 
    real,         intent(in) :: th_pvsrf(jms:jme)          ! 
    real,         intent(in) :: rcav(jms:jme)          ! 
    real,         intent(in) :: hf2d(jms:jme)          ! 
    real,         intent(in) :: snwdph2d(jms:jme)      !
    real,         intent(in) :: rsds(jms:jme)          ! downward short-wave radiation flux
    real,         intent(in) :: pr3d(nvlp1,jms:jme)    ! pressure (pascal)
    real,         intent(in) :: ph3d(nvlp1,jms:jme)    ! geopotential (=gz), m^2/s^2
    real,         intent(in) :: emiss_ab(jms:jme,num_emis_ant)           ! 
    real,         intent(in) :: emiss_abu(jms:jme,num_ebu_in)           ! 
    real,         intent(in) :: plumestuff(jms:jme,8)           ! 
    real,         intent(in) :: ero1(jms:jme)           ! 
    real,         intent(in) :: ero2(jms:jme)           ! 
    real,         intent(in) :: ero3(jms:jme)           ! 
    real,      intent(inout) :: emiss_ash_mass(jms:jme)           ! 
    real,      intent(inout) :: emiss_ash_height(jms:jme)           ! 
    real,         intent(in) :: emiss_ash_dt(jms:jme)           ! 
    real,         intent(in) :: emiss_tr_mass(jms:jme)           ! 
    real,         intent(in) :: emiss_tr_height(jms:jme)           ! 
    real,         intent(in) :: emiss_tr_dt(jms:jme)           ! 
    real,         intent(in) :: dm0(jms:jme)           ! 
    real,         intent(in) :: p_gocart(56)           ! 
    real,         intent(in) :: area(jms:jme)             ! the area of cell polygon (m**2)
    real, dimension (jms:jme), intent(in) :: vfrac2d,VTYPE2d,STYPE2d,zorl2d,slmsk2d
    real, dimension (nvl,jms:jme), intent(in) :: us3d,vs3d,ws3d
    real,         intent(in) :: deg_lat(jms:jme),deg_lon(jms:jme)  ! lat and lon in degrees

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
          INTENT(OUT ) ::                                   moist
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
          INTENT(OUT ) ::                                   chem
    REAL, DIMENSION( ims:ime, kms:kemit, jms:jme, num_emis_ant ),                 &
          INTENT(inout ) ::                                   emis_ant
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_emis_vol ),                 &
          INTENT(inout ) ::                                   emis_vol
    REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
           INTENT(OUT   ) ::                                                 &
                                                         rri,               &
                                                       t_phy,               &
                                                       p_phy,               &
                                               relhum, dz8w,p8w,t8w,        &
                                               z_at_w , zmid ,exch_h,       &
                                               u_phy,v_phy,vvel,rho_phy,    &
                                               convfac
    REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
           INTENT(OUT   ) ::                                                 &
                                    backg_oh,backg_h2o2,backg_no3
    REAL,  DIMENSION( ims:ime , jms:jme )         ,               &
           INTENT(OUT   ) ::                                                 &
                              ttday,tcosz
    REAL,DIMENSION( ims:ime , jms:jme,num_ebu_in )           ,               &
           INTENT(OUT   ) ::                                                 &
          ebu_in
    REAL,DIMENSION( ims:ime , jms:jme )                  ,               &
           INTENT(OUT   ) ::                                                 &
          mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,            &
          firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr
    INTEGER,DIMENSION( ims:ime , jms:jme )                  ,               &
           INTENT(OUT   ) ::                                                 &
                                                      ivgtyp,               &
                                                      isltyp
    REAL, DIMENSION( ims:ime, jms:jme,3)::&
          erod

    REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &
           INTENT(OUT   ) ::                                                 &
                                                      u10,                  &
                                                      v10,                  &
                                                      gsw,                  &
                                                   vegfra,                  &
                                                      rmol,                 &
                                                      ust,                  &
                                                      snowh,                &
                                                      xland,                &
                                                      xlat,e_co,dms_0,      &
                                                      xlong,tsk,raincv_b,   &
                                                      dxy,znt,pbl,hfx
    REAL, DIMENSION( ims:ime, num_soil_layers, jms:jme ) ,      &
       INTENT(OUT) ::                               smois
    integer i,j,k,kk,nv,jmax,jmaxi,l,ll,n,ndystep,ixhour
    real maxv,factor,factor2,pu,pl,aln,pwant,rlat
    real thv,xhour,xmin,gmtp,xlonn,xtime,real_time
    real, DIMENSION (1,1) :: sza,cosszax
    real, DIMENSION (jms:jme) :: so2_mass

    ! -- volcanic stuff
    integer :: ko,k_final,k_initial,kl,kk4,curr_hours,curr_secs
    real :: x1,ashz_above_vent
    real, DIMENSION (kms:kme) :: vert_mass_dist
    real :: eh,h1,h2,h3,h4,h5,h6,maxth
    logical, save :: first_init = .true.
    integer :: clockId

    ! -- volcano ashes parameters
    !  + original
    ! real, dimension(6) :: h = (/ 9., 16., 58., 79., 109., 129., huge(1.0) /)
    !  + real-time default (if volcano starts at h = 0)
    real, dimension(7) :: h = (/ (240., i = 1, 6), huge(1.0) /)
    real, dimension(6) :: emiss_ash_table = (/  5834.,  3834.,  5834.,  3334.,  3334.,  2334. /)
    real, dimension(6) :: eh_table        = (/ 3.11e5, 3.87e4, 3.11e5, 2.17e4, 2.17e4, 4.93e3 /)
    real, parameter :: percen_mass_umbrel = 0.75
    real, parameter :: base_umbrel        = 0.25    ! fraction
    real, parameter :: base_umbrel2       = 1.0     ! evenly distribution
    ! .. Intrinsic Functions ..
    INTRINSIC max, min, float

    ! -- begin
    clockId = mpp_clock_id('chem_prep')
    call mpp_clock_begin(clockId)

    print *,'chem_prep: begin...'

    real_time=float(ktau)*dtstep/60.

    so2_mass = 0.

    if (ktau <= 1) then
      emis_ant = 0.
      emis_vol = 0.
    end if

    e_co = 0.

    do j=jts,jte
      do i=its,ite
         z_at_w(i,kts,j)=max(0.,ph3d(kts,j)/g)
      enddo
    enddo

    do j=jts,jte
      do k=kts,kte
        do i=its,ite
          dz8w(i,k,j)=(ph3d(k+1,j)-ph3d(k,j))/g
          if (dz8w(i,k,j) < 0.) dz8w(i,k,j)=-dz8w(i,k,j)
          z_at_w(i,k+1,j)=z_at_w(i,k,j)+dz8w(i,k,j)
        enddo
      enddo
    enddo

    do j=jts,jte
      do k=kts,kte+1
        do i=its,ite
          p8w(i,k,j)=pr3d(k,j)
        enddo
      enddo
    enddo

    do j=jts,jte
      do i=its,ite
        raincv_b(i,j)=rcav(j)
        pbl(i,j)=pb2d(j)
        dms_0(i,j)=dm0(j)
        hfx(i,j)=hf2d(j)
        snowh(i,j)=snwdph2d(j)*.001
        erod(i,j,1)=ero1(j)
        erod(i,j,2)=ero2(j)
        erod(i,j,3)=ero3(j)
        xlat(i,j)=deg_lat(j)
        xlong(i,j)=deg_lon(j)
        ust(i,j)=us2d(j)
        tsk(i,j)=ts2d(j)
        ivgtyp(i,j)=VTYPE2d(j)
        isltyp(i,j)=STYPE2d(j)
        gsw(i,j)=rsds(j)
        vegfra(i,j)=VFRAC2d(j)
        rmol(i,j)=0.
        znt(i,j)=zorl2d(j)*.01
!SLMSK   - SEA(0),LAND(1),ICE(2) MASK
        xland(i,j)=1.
        if (slmsk2d(j) == 0.) then
          xland(i,j) = 0.
        else if (slmsk2d(j) == 1.) then
          xland(i,j) = 1.
        else if (slmsk2d(j) == 2.) then
          xland(i,j) = 2.
        end if
        dxy(i,j)=area(j)
        u10(i,j)=us3d(1,j)
        v10(i,j)=vs3d(1,j)
        clayf(i,j) = clayfrac(j)
        sandf(i,j) = sandfrac(j)
      enddo
    enddo


    factor=0.
    jmax=0
    jmaxi=0
    k=1
    if ((p_bc2 > 1) .or. (config % chem_opt == 108)) then  ! "regular" chem options
      do j=jts,jte
        do i=its,ite
          k=1
          emis_ant(i,k,j,p_e_bc)=emiss_ab(j,p_e_bc)
          emis_ant(i,k,j,p_e_oc)=emiss_ab(j,p_e_oc)
          emis_ant(i,k,j,p_e_sulf)=emiss_ab(j,p_e_sulf)
          emis_ant(i,k,j,p_e_so2)=emiss_ab(j,p_e_so2)
          emis_ant(i,k,j,p_e_dms)= 0. !emiss_ab(j,p_e_dms)
          emis_ant(i,k,j,p_e_pm_25)=emiss_ab(j,p_e_pm_25)
          emis_ant(i,k,j,p_e_pm_10)=emiss_ab(j,p_e_pm_10)

          ! -- gas anth emission for chem_opt 301
          if ((config % chem_opt == 301) .or. (config % chem_opt == 108)) then
            emis_ant(i,k,j,p_e_iso)=emiss_ab(j,p_e_iso)
            emis_ant(i,k,j,p_e_no)=emiss_ab(j,p_e_no)
            emis_ant(i,k,j,p_e_no2)=emiss_ab(j,p_e_no2)
            emis_ant(i,k,j,p_e_co)=emiss_ab(j,p_e_co)
            emis_ant(i,k,j,p_e_eth)=emiss_ab(j,p_e_eth)
            emis_ant(i,k,j,p_e_hc3)=emiss_ab(j,p_e_hc3)
            emis_ant(i,k,j,p_e_hc5)=emiss_ab(j,p_e_hc5)
            emis_ant(i,k,j,p_e_hc8)=emiss_ab(j,p_e_hc8)
            emis_ant(i,k,j,p_e_xyl)=emiss_ab(j,p_e_xyl)
            emis_ant(i,k,j,p_e_olt)=emiss_ab(j,p_e_olt)
            emis_ant(i,k,j,p_e_oli)=emiss_ab(j,p_e_oli)
            emis_ant(i,k,j,p_e_tol)=emiss_ab(j,p_e_tol)
            emis_ant(i,k,j,p_e_csl)=emiss_ab(j,p_e_csl)
            emis_ant(i,k,j,p_e_hcho)=emiss_ab(j,p_e_hcho)
            emis_ant(i,k,j,p_e_ald)=emiss_ab(j,p_e_ald)
            emis_ant(i,k,j,p_e_ket)=emiss_ab(j,p_e_ket)
            emis_ant(i,k,j,p_e_ora2)=emiss_ab(j,p_e_ora2)
            emis_ant(i,k,j,p_e_nh3)=emiss_ab(j,p_e_nh3)
          endif
         
          ! -- gas biomass burning emission for chem_opt 301
          if ((config % chem_opt == 301) .or. (config % chem_opt == 108)) then
            ebu_in(i,j,p_ebu_in_iso)=emiss_abu(j,p_e_iso)
            ebu_in(i,j,p_ebu_in_no)=emiss_abu(j,p_e_no)
            ebu_in(i,j,p_ebu_in_no2)=emiss_abu(j,p_e_no2)
            ebu_in(i,j,p_ebu_in_co)=emiss_abu(j,p_e_co)
            ebu_in(i,j,p_ebu_in_eth)=emiss_abu(j,p_e_eth)
            ebu_in(i,j,p_ebu_in_hc3)=emiss_abu(j,p_e_hc3)
            ebu_in(i,j,p_ebu_in_hc5)=emiss_abu(j,p_e_hc5)
            ebu_in(i,j,p_ebu_in_hc8)=emiss_abu(j,p_e_hc8)
            ebu_in(i,j,p_ebu_in_xyl)=emiss_abu(j,p_e_xyl)
            ebu_in(i,j,p_ebu_in_olt)=emiss_abu(j,p_e_olt)
            ebu_in(i,j,p_ebu_in_oli)=emiss_abu(j,p_e_oli)
            ebu_in(i,j,p_ebu_in_tol)=emiss_abu(j,p_e_tol)
            ebu_in(i,j,p_ebu_in_csl)=emiss_abu(j,p_e_csl)
            ebu_in(i,j,p_ebu_in_hcho)=emiss_abu(j,p_e_hcho)
            ebu_in(i,j,p_ebu_in_ald)=emiss_abu(j,p_e_ald)
            ebu_in(i,j,p_ebu_in_ket)=emiss_abu(j,p_e_ket)
            ebu_in(i,j,p_ebu_in_ora2)=emiss_abu(j,p_e_ora2)
            ebu_in(i,j,p_ebu_in_nh3)=emiss_abu(j,p_e_nh3)
          endif
          ebu_in(i,j,p_ebu_in_oc)=emiss_abu(j,p_e_oc)
          ebu_in(i,j,p_ebu_in_bc)=emiss_abu(j,p_e_bc)
          ebu_in(i,j,p_ebu_in_pm25)=emiss_abu(j,p_e_pm_25)
          ebu_in(i,j,p_ebu_in_pm10)=emiss_abu(j,p_e_pm_10)
          ebu_in(i,j,p_ebu_in_so2)=emiss_abu(j,p_e_so2)
          ebu_in(i,j,p_ebu_in_dms)= 0. !emiss_abu(j,p_e_dms)
         !ebu_in(i,j,p_ebu_in_sulf)=0. ! for now

          mean_fct_agtf(i,j)=plumestuff(j,1)
          mean_fct_agef(i,j)=plumestuff(j,2)
          mean_fct_agsv(i,j)=plumestuff(j,3)
          mean_fct_aggr(i,j)=plumestuff(j,4)
          firesize_agtf(i,j)=plumestuff(j,5)
          firesize_agef(i,j)=plumestuff(j,6)
          firesize_agsv(i,j)=plumestuff(j,7)
          firesize_aggr(i,j)=plumestuff(j,8)
        enddo
      enddo

    else if (p_tr2 > 1) then  ! tracer options

      ! -- tracer run
      do j=jts,jte
        do i=its,ite
          k=kts
          emis_ant(i,k,j,p_e_tr1)=emiss_ab(j,p_e_tr1)
          emis_ant(i,k,j,p_e_tr2)=emiss_ab(j,p_e_tr2)
        enddo
      enddo

    else if ((p_tr2 > 1) .and. (p_bc2 > 1)) then

      call mpp_error(FATAL, 'in chem_prep_fim, 111')

    endif

    do i=its,ite
      do j=jts,jte
        do k=kts,kte
          thv=tr3d(k,j,1)/(1.+0.6078*tr3d(k,j,2))
          tk3d(k,j)=thv*(.5*(p8w(i,k,j)+p8w(i,k+1,j))/p1000)**(rd/cp)
        enddo
      enddo
    enddo

    do j=jts,jte
      do k=kts,kte+1
        kk=min(k,kte)
        do i=its,ite
          zmid(i,k,j)=.5*(ph3d(kk+1,j)+ph3d(kk,j))/g
          dz8w(i,k,j)=z_at_w(i,kk+1,j)-z_at_w(i,kk,j)
          t_phy(i,k,j)=tk3d(kk,j)
!         relhum(i,k,j)=rh3d(kk,j)
          p_phy(i,k,j)=.5*(p8w(i,kk,j)+p8w(i,kk+1,j))
          u_phy(i,k,j)=us3d(kk,j)
          exch_h(i,k,j)=exch(kk,j)
          v_phy(i,k,j)=vs3d(kk,j)
!       print *,'--> i,k,j, RD, p, T, t3d = ',i,k,j,RD,p_phy(i,k,j),T_phy(i,k,j),tr3d(kk,j,2)
          rho_phy(i,k,j)= p_phy(i,k,j)/(RD*T_phy(i,k,j)*(1.+.608*tr3d(kk,j,2)))
          rri(i,k,j)=1./rho_phy(i,k,j)
!       print *,'--> i,k,j, rho, 1/rho = ',i,k,j,rho_phy(i,k,j),rri(i,k,j)
          vvel(i,k,j)=-ws3d(kk,j)*rri(i,k,j)/g
          convfac(i,k,j)=p_phy(i,k,j)/rgasuniv/t_phy(i,k,j)
          moist(i,k,j,:)=0.
          moist(i,k,j,1)=tr3d(kk,j,2)
          if (t_phy(i,k,j) > 265.) then
            moist(i,k,j,2)=tr3d(kk,j,3)
            moist(i,k,j,3)=0.
            if (moist(i,k,j,2) < 1.e-8) moist(i,k,j,2)=0.
          else
            moist(i,k,j,2)=0.
            moist(i,k,j,3)=tr3d(kk,j,3)
            if(moist(i,k,j,3) < 1.e-8)moist(i,k,j,3)=0.
          endif
          relhum(i,k,j) = .95
          relhum(i,k,j) = MIN( .95, moist(i,k,j,1) / &
            (3.80*exp(17.27*(t_phy(i,k,j)-273.)/ &
            (t_phy(i,k,j)-36.))/(.01*p_phy(i,k,j))))
          relhum(i,k,j)=max(0.1,relhum(i,k,j))
        enddo
      enddo
    enddo

    do j=jts,jte
      do k=2,kte
        do i=its,ite
          t8w(i,k,j)=.5*(t_phy(i,k,j)+t_phy(i,k-1,j)) ! .5*(tk3d(k-1,j)+tk3d(k,j))
        enddo
      enddo
    enddo

    ! -- only used in phtolysis....
    do j=jts,jte
      do i=its,ite
        t8w(i,1,j)=t_phy(i,1,j)
        t8w(i,kte+1,j)=t_phy(i,kte,j)
      enddo
    enddo

    do nv=1,num_chem
      do j=jts,jte
        do k=kts,kte+1
          kk=min(k,kte)
          do i=its,ite
            chem(i,k,j,nv)=tr3d(kk,j,ntra+nv)
          enddo
        enddo
      enddo
    enddo

    print *,'chem_prep: chem array done'
    print *,'chem_prep: readrestart  =',config % readrestart

    if (.NOT. config % readrestart) then
      if (config % chem_in_opt == 0 ) then
        if(ktau.le.1)then
!           if(chem_opt > 0 ) then
          do j=jts,jte
            do k=kts,kte
              do i=its,ite
                if (config % chem_opt == 300) then
                  do n=1,num_chem
                    chem(i,k,j,n)=1.e-12
                  enddo
                endif  ! chem_opt==300
                chem(i,k,j,p_so2)=5.e-6
                chem(i,k,j,p_sulf)=3.e-6
                if ((config % chem_opt >= 300) .and. (config % chem_opt <  500)) then
                  chem(i,k,j,p_msa)=0.1e-6
                  chem(i,k,j,p_dms)=0.1e-6
                  chem(i,k,j,p_bc1)=0.1e-3
                  chem(i,k,j,p_bc2)=0.1e-3
                  chem(i,k,j,p_oc1)=0.1e-3
                  chem(i,k,j,p_oc2)=0.1e-3
                  chem(i,k,j,p_p25)=1.
                  chem(i,k,j,p_p10)=1.
                endif !chem_opt >= 300 .and. chem_opt <  500

                if ((config % chem_opt == 301) .or. (config % chem_opt == 108)) then  !added o3 background !lzhang
                  kk=min(k,kte)
                  ! -- add initial constant into O3,CH4 and CO ect.
                  chem(i,k,j,p_o3)=epsilc
                  maxth=min(400.,th_pvsrf(j))
                  if (tr3d(kk,j,1) > maxth) then 
                    chem(i,k,j,p_o3)=(airmw/48.)*tr3d(kk,j,4)*1e6 !convert kg/kg to ppm
                  else
                    chem(i,k,j,p_o3)=0.03 !ppm
                  endif
                  chem(i,k,j,p_ch4)=1.85 !ppm
                  chem(i,k,j,p_co)=0.06 !ppm
                  chem(i,k,j,p_co2)=380.
                  chem(i,k,j,p_ete)=epsilc
                  chem(i,k,j,p_udd)=chem(i,k,j,p_ete)
                  chem(i,k,j,p_hket)=chem(i,k,j,p_ete)
                  chem(i,k,j,p_api)=chem(i,k,j,p_ete)
                  chem(i,k,j,p_lim)=chem(i,k,j,p_ete)
                  chem(i,k,j,p_dien)=chem(i,k,j,p_ete)
                  chem(i,k,j,p_macr)=chem(i,k,j,p_ete)
                endif !( (chem_opt == 301.or.chem_opt==108))
              enddo
            enddo
          enddo
        endif !(ktau<=1)

      else !(chem_in_opt == 0 )

        if ((ktau<=1).and.((config % chem_opt == 301).or.(config % chem_opt == 108))) then  !added GFS o3 background above 380K!lzhang
          do j=jts,jte
            maxth=min(400.,th_pvsrf(j))
            do k=kts,kte+1
              kk=min(k,kte)
              if (tr3d(kk,j,1) >= maxth) then
                chem(its:ite,k,j,p_o3)=(airmw/48.)*tr3d(kk,j,4)*1e6 !convert kg/kg to ppm
              endif !380K
            enddo
          enddo
        endif ! chem_opt == 301.or.chem_opt==108
        
      endif !(chem_in_opt == 1 )

      print *,'chem_prep: chem in done'
      
#if 0
      if (ktau<=1.and.config % chem_opt == 108) then
           call aerosols_soa_vbs_init(chem,convfac,z_at_w,           &
               pm2_5_dry,pm2_5_water,pm2_5_dry_ec,                   &
               config % chem_in_opt, aer_ic_opt,                              &
               ids,ide, jds,jde, kds,kde,                            &
               ims,ime, jms,jme, kms,kme,                            &
               its,ite, jts,jte, kts,kte                             )
        !!!TUCCELLA (BUG, before it was called in module_aerosols_soa_vbs.F)
       ! initialize pointers used by aerosol-cloud-interaction routines

       if( .not.allocated(is_aerosol) ) then
        allocate (is_aerosol(num_chem))
       else
        if( size(is_aerosol) /= num_chem ) &
          print *, 'The number of chemistry species has changed between nests. Are you trying to mix chem_opt settings between nests? Shame on you!'
       end if
           call aerosols_soa_vbs_init_aercld_ptrs(                    &
                num_chem, is_aerosol )
        !...Convert aerosols to mixing ratio
     !   if (.NOT. readrestart) then
     !   if( chem_in_opt == 0 .and. num_chem.gt.numgas)then
     !   do l=numgas+1,num_chem
     !      do j=jts,jte
     !         do k=kts,kte
     !            kk = min(k,kde-1)
     !            do i=its,ite
     !               chem(i,k,j,l)=chem(i,kk,j,l)*rri(i,kk,j)
     !            enddo
     !         enddo
     !      enddo
     !   enddo
     !   endif 
     !   endif


        chem(its:ite,kts:min(kte,kde-1),jts:jte,:)=max(chem(its:ite,kts:min(kte,kde-1),jts:jte,:),epsilc)

       endif !ktau<=1.and.chem_opt=108
#endif
    else !restart
#if 0
if (first_init .and. config % chem_opt == 108) then
           first_init = .false.
           call aerosols_soa_vbs_init(chem,convfac,z_at_w,           &
               pm2_5_dry,pm2_5_water,pm2_5_dry_ec,                   &
               config % chem_in_opt, aer_ic_opt,                              &
               ids,ide, jds,jde, kds,kde,                            &
               ims,ime, jms,jme, kms,kme,                            &
               its,ite, jts,jte, kts,kte                             )

     if( .not.allocated(is_aerosol) ) then
        allocate (is_aerosol(num_chem))
     else
        if( size(is_aerosol) /= num_chem ) &
          print *, 'The number of chemistry species has changed between nests. Are you trying to mix chem_opt settings between nests? Shame on you!'
     end if
          call aerosols_soa_vbs_init_aercld_ptrs(                    &
                num_chem, is_aerosol )
chem(its:ite,kts:min(kte,kde-1),jts:jte,:)=max(chem(its:ite,kts:min(kte,kde-1),jts:jte,:),epsilc)
endif
#endif
!
    endif ! restart

    !
    ! -- gocart background fields only if gocart is called
    !
    !if (.NOT. readrestart) then
    if (call_gocart .and. (config % chem_opt == 300))then
      print *,'chem_prep: call_gocart enter...'
      do j=jts,jte
        do k=kts,kte
          do i=its,ite
            do ll=2,nvl_gocart
              l=ll
              if (p_gocart(l) < .01*p_phy(i,k,j)) exit
            enddo
            pu=alog(p_gocart(l))
            pl=alog(p_gocart(l-1))
            pwant=alog(.01*p_phy(i,k,j))
            if (pwant > pl)then
              backg_oh(i,k,j)=oh_backgd(l,j)
              backg_h2o2(i,k,j)=h2o2_backgd(l,j)
              backg_no3(i,k,j)=no3_backgd(l,j)
            else
              aln=(oh_backgd(l,j)*(pwant-pl)+            &
                oh_backgd(l-1,j)*(pu-pwant))/(pu-pl)
              backg_oh(i,k,j)=aln
              aln=(h2o2_backgd(l,j)*(pwant-pl)+            &
                h2o2_backgd(l-1,j)*(pu-pwant))/(pu-pl)
              backg_h2o2(i,k,j)=aln
              aln=(no3_backgd(l,j)*(pwant-pl)+            &
                no3_backgd(l-1,j)*(pu-pwant))/(pu-pl)
              backg_no3(i,k,j)=aln
            endif
          enddo
        enddo
      enddo
      print *,'chem_prep: call_gocart exit'
    endif   ! end gocart stuff
    !endif !restart

    print *,'chem_prep: chem + emiss: starting ...'
!   emis_ant=0.
    nv=1
    k=1
    factor2=0.
    factor=0.
    if (p_bc2 > 1)then
      if (config % chem_opt == 300) then
        do j=jts,jte
          do i=its,ite
            factor=dtstep*rri(i,k,j)/dz8w(i,k,j)
            factor2=4.828e-4*dtstep*rri(i,k,j)/(60.*dz8w(i,k,j))
            chem(i,k,j,p_bc1)=chem(i,k,j,p_bc1)+emis_ant(i,k,j,p_e_bc)*factor
            chem(i,k,j,p_oc1)=chem(i,k,j,p_oc1)+emis_ant(i,k,j,p_e_oc)*factor
            chem(i,k,j,p_p25)=chem(i,k,j,p_p25)+emis_ant(i,k,j,p_e_pm_25)*factor
            chem(i,k,j,p_p10)=chem(i,k,j,p_p10)+emis_ant(i,k,j,p_e_pm_10)*factor
            chem(i,k,j,p_sulf)=chem(i,k,j,p_sulf)+emis_ant(i,k,j,p_e_sulf)*factor
            chem(i,k,j,p_so2)=chem(i,k,j,p_so2)+emis_ant(i,k,j,p_e_so2)*factor2
          enddo
        enddo
      endif
    else if (p_tr2 > 1)then    !co2 here
      do j=jts,jte
        do i=its,ite
!           factor2=dtstep*rri(i,k,j)/dz8w(i,k,j)
          factor2=4.828e-4*dtstep*rri(i,k,j)/(60.*dz8w(i,k,j))
          chem(i,k,j,p_tr1)=chem(i,k,j,p_tr1)+emis_ant(i,k,j,p_e_tr1)*factor2
          chem(i,k,j,p_tr2)=chem(i,k,j,p_tr2)+emis_ant(i,k,j,p_e_tr2)*factor2
        enddo
      enddo   
    else if ((p_tr2 > 1) .and. (p_bc2 > 1))then
      call mpp_error(FATAL, 'in chem_prep_fim, 112')
    endif
    print *,'chem_prep: chem + emiss: done ...'

    print *,'chem_prep: smois: start'
    do j=jts,jte
      do nv=1,num_soil_layers
        smois(its:ite,nv,j)=sm3d(nv,j)
      enddo
    enddo
    print *,'chem_prep: smois: end'

    curr_secs=ktau*ifix(dtstep)
    curr_hours=curr_secs/3600
!
!     do volcanoes if avaiable
!
!     if(chem_opt == 502 ) then
!       do j=jts,jte
!       if(emiss_ash_dt(j).le.0)CYCLE
!       emiss_ash_mass(j)=0.
!       emiss_ash_height(j)=0.
!       enddo
!
! default
!
    do j=jts,jte
      if (emiss_ash_dt(j) > 0.) then
        so2_mass(j)=1.5e4*3600.*1.e9/64./area(j)
        eh=2600.*(emiss_ash_height(j)*.0005)**4.1494
        emiss_ash_mass(j)=eh*1.e9/area(j)
      end if
    enddo

    ! -- hardcode for special retro case (set h1 - h6 properly
    do nv = 1, 6
      if ((curr_hours >= h(nv)) .and. (curr_hours < h(nv+1))) then
        do j = jts, jte
          if (emiss_ash_dt(j) > 0) then
            emiss_ash_height(j) = emiss_ash_table(nv)
            emiss_ash_mass(j)   = 1.e+09 * eh_table(nv) / area(j)
          end if
        end do
        exit
      end if
    end do

    print *,'chem_prep: stage 2'
!     endif ! chem_opt = 502
!
    ! -- real-time application, keeping eruption constant
!
    if (ktau <= 2) then
      ! -- volcanic emissions
      emis_vol(:,:,:,:)=0.
!      if(curr_hours.eq.h1 .or. curr_hours.eq.h2 .or. curr_hours.eq.h3 &
!         .or. curr_hours.eq.h4 .or. curr_hours.eq.h5 .or. curr_hours.eq.h6 .or. h1.gt.239)then
!         .or. curr_hours.eq.0)then
!         if(chem_opt == 316 .or. chem_opt == 317 .or. chem_opt == 502) then

      do j=jts,jte
        if (emiss_ash_dt(j)     <= 0) cycle
        if (emiss_ash_height(j) <= 0) cycle
        do i=its,ite
          ashz_above_vent=emiss_ash_height(j) +z_at_w(i,kts,j)
          do k=kte-1,kts,-1
            if (z_at_w(i,k,j) < ashz_above_vent)then
              k_final=k+1
              exit
            endif !inner
          enddo
          do k=kte-1,kts,-1
            if (z_at_w(i,k,j) < (1.-base_umbrel)*ashz_above_vent)then
              k_initial=k
              exit
            endif !inner
          enddo
          vert_mass_dist=0.
!              k_initial=int((k_final+k_initial)*0.5)
             
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
          if (config % chem_opt == 316 ) then 
            do ko=1,k_final
              emis_vol(i,ko,j,p_e_vash1)=.02*vert_mass_dist(ko)*emiss_ash_mass(j)
              emis_vol(i,ko,j,p_e_vash2)=.04*vert_mass_dist(ko)*emiss_ash_mass(j)
              emis_vol(i,ko,j,p_e_vash3)=.11*vert_mass_dist(ko)*emiss_ash_mass(j)
              emis_vol(i,ko,j,p_e_vash4)=.09*vert_mass_dist(ko)*emiss_ash_mass(j)
              emis_vol(i,ko,j,p_e_vash5)=.09*vert_mass_dist(ko)*emiss_ash_mass(j)
              emis_vol(i,ko,j,p_e_vash6)=.13*vert_mass_dist(ko)*emiss_ash_mass(j)
              emis_vol(i,ko,j,p_e_vash7)=.16*vert_mass_dist(ko)*emiss_ash_mass(j)
              emis_vol(i,ko,j,p_e_vash8)=.16*vert_mass_dist(ko)*emiss_ash_mass(j)
              emis_vol(i,ko,j,p_e_vash9)=.1*vert_mass_dist(ko)*emiss_ash_mass(j)
              emis_vol(i,ko,j,p_e_vash10)=.1*vert_mass_dist(ko)*emiss_ash_mass(j)
            enddo
            do ko=k_final+1,kte
              emis_vol(i,ko,j,p_e_vash1)=0.
              emis_vol(i,ko,j,p_e_vash2)=0.
              emis_vol(i,ko,j,p_e_vash3)=0.
              emis_vol(i,ko,j,p_e_vash4)=0.
              emis_vol(i,ko,j,p_e_vash5)=0.
              emis_vol(i,ko,j,p_e_vash6)=0.
              emis_vol(i,ko,j,p_e_vash7)=0.
              emis_vol(i,ko,j,p_e_vash8)=0.
              emis_vol(i,ko,j,p_e_vash9)=0.
              emis_vol(i,ko,j,p_e_vash10)=0.
            enddo
          elseif (config % chem_opt == 317 .or. config % chem_opt == 502) then

            ! -- reduced vocanic ash transport
            do ko=1,k_final
              emis_vol(i,ko,j,p_e_vash1)=.11*vert_mass_dist(ko)*emiss_ash_mass(j)
              emis_vol(i,ko,j,p_e_vash2)=.08*vert_mass_dist(ko)*emiss_ash_mass(j)
              emis_vol(i,ko,j,p_e_vash3)=.05*vert_mass_dist(ko)*emiss_ash_mass(j)
              emis_vol(i,ko,j,p_e_vash4)=.035*vert_mass_dist(ko)*emiss_ash_mass(j)
            enddo
          elseif (config % chem_opt == 300) then

            ! -- if applied to gocart we only need finest ash bins, we use the coarse one for so2
            do ko=1,k_final
              emis_vol(i,ko,j,p_e_vash1)=vert_mass_dist(ko)*so2_mass(j)
              emis_vol(i,ko,j,p_e_vash2)=.08*vert_mass_dist(ko)*emiss_ash_mass(j)
              emis_vol(i,ko,j,p_e_vash3)=.05*vert_mass_dist(ko)*emiss_ash_mass(j)
              emis_vol(i,ko,j,p_e_vash4)=.035*vert_mass_dist(ko)*emiss_ash_mass(j)
            enddo
          endif !chem_opt==316 or 317,300,502
 
          do ko=k_final+1,kte
            emis_vol(i,ko,j,p_e_vash1)=0.
            emis_vol(i,ko,j,p_e_vash2)=0.
            emis_vol(i,ko,j,p_e_vash3)=0.
            emis_vol(i,ko,j,p_e_vash4)=0.
          enddo
        enddo
      enddo
!      endif ! chem_opt == 316 .or. chem_opt == 317 .or. chem_opt == 502
    endif ! curr_mins 
    print *,'chem_prep: stage 3'

    ! next is done to scale background oh and no3 in dependence on average zenith angle and day/night for no3
    ! this is done since background values are only available as average/month. It will not be necessary if other
    ! chemistry packages are used that provide oh,no3,h2o2

    !if(ktau.le.1 .or.readrestart)then
    print *,'-- chem_opt = ',config % chem_opt
    if ((config % chem_opt == 108) .or. (config % chem_opt >= 300 .and. config % chem_opt <  500)) then
      !ndystep=86400/ifix(dtstepc)
      ndystep=86400/ifix(dtstep)
      do j=jts,jte
        do i=its,ite
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
    endif !chem_opt >= 300 .and. chem_opt <  500
         !endif ! ktau

    print *,'chem_prep: stage 4'
    if (config % chem_opt == 316) then
      do j = jts, jte
        if (emiss_ash_dt(j) <= 0.) cycle
        do k = kts, kte-2
          do i = its, ite
            factor2=dtstep*rri(i,k,j)/dz8w(i,k,j)
            chem(i,k,j,p_vash_1)=chem(i,k,j,p_vash_1)+emis_vol(i,k,j,p_e_vash1)*factor2
            chem(i,k,j,p_vash_2)=chem(i,k,j,p_vash_2)+emis_vol(i,k,j,p_e_vash2)*factor2
            chem(i,k,j,p_vash_3)=chem(i,k,j,p_vash_3)+emis_vol(i,k,j,p_e_vash3)*factor2
            chem(i,k,j,p_vash_4)=chem(i,k,j,p_vash_4)+emis_vol(i,k,j,p_e_vash4)*factor2
            chem(i,k,j,p_vash_5)=chem(i,k,j,p_vash_5)+emis_vol(i,k,j,p_e_vash5)*factor2
            chem(i,k,j,p_vash_6)=chem(i,k,j,p_vash_6)+emis_vol(i,k,j,p_e_vash6)*factor2
            chem(i,k,j,p_vash_7)=chem(i,k,j,p_vash_7)+emis_vol(i,k,j,p_e_vash7)*factor2
            chem(i,k,j,p_vash_8)=chem(i,k,j,p_vash_8)+emis_vol(i,k,j,p_e_vash8)*factor2
            chem(i,k,j,p_vash_9)=chem(i,k,j,p_vash_9)+emis_vol(i,k,j,p_e_vash9)*factor2
            chem(i,k,j,p_vash_10)=chem(i,k,j,p_vash_10)+emis_vol(i,k,j,p_e_vash10)*factor2
          enddo
        enddo
      enddo
    endif

    if ((config % chem_opt == 317) .or. (config % chem_opt == 502))  then
      do j=jts,jte
        if (emiss_ash_dt(j) <= 0.) cycle
        do k=kts,kte-2
          do i=its,ite
            factor2=dtstep*rri(i,k,j)/dz8w(i,k,j)
            chem(i,k,j,p_vash_1)=chem(i,k,j,p_vash_1)+emis_vol(i,k,j,p_e_vash1)*factor2
            chem(i,k,j,p_vash_2)=chem(i,k,j,p_vash_2)+emis_vol(i,k,j,p_e_vash2)*factor2
            chem(i,k,j,p_vash_3)=chem(i,k,j,p_vash_3)+emis_vol(i,k,j,p_e_vash3)*factor2
            chem(i,k,j,p_vash_4)=chem(i,k,j,p_vash_4)+emis_vol(i,k,j,p_e_vash4)*factor2
          enddo
        enddo
      enddo
    endif

    if(config % chem_opt == 300)  then
      ! -- for gocart only lump ash into p25 and p10
      do j=jts,jte
        if (emiss_ash_dt(j) <= 0.) cycle
        do k=kts,kte-2
          do i=its,ite
            factor=4.828e-4*dtstep*rri(i,k,j)/(60.*dz8w(i,k,j))
            factor2=dtstep*rri(i,k,j)/dz8w(i,k,j)
            chem(i,k,j,p_p25)=chem(i,k,j,p_p25)                          &
                             +emis_vol(i,k,j,p_e_vash4)*factor2   
            chem(i,k,j,p_so2)=chem(i,k,j,p_so2)                          &
                             +emis_vol(i,k,j,p_e_vash1)*factor   
            chem(i,k,j,p_p10)=chem(i,k,j,p_p10)                          &
!                         +.5* emis_vol(i,k,j,p_e_vash4)*factor2      &  
                             +1.* emis_vol(i,k,j,p_e_vash3)*factor2      &  
                             +.5* emis_vol(i,k,j,p_e_vash2)*factor2
          enddo
        enddo
      enddo
    endif
!
    print *,'chem_prep: stage 5'

    ! -- option 501 was only used for cesium ensemble - Japan 2010
    if (config % chem_opt == 501) then
      ! -- explosive tr emissions
      do j=jts,jte
        do i=its,ite
          if ((emiss_tr_dt(j) <= 0.) .or. (emiss_tr_height(j) <= 0.)) cycle
            ashz_above_vent=emiss_tr_height(j)+z_at_w(i,kts,j)
            do k=kte-1,kts,-1
              if (z_at_w(i,k,j) < ashz_above_vent) then
                k_final=k+1
                exit
              endif
            enddo
            do k=kte-1,kts,-1
              if (z_at_w(i,k,j) < (1.-base_umbrel)*ashz_above_vent) then
                k_initial=k
                exit
              endif
            enddo
            vert_mass_dist=0.
            k_initial=int((k_final+k_initial)*0.5)

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
            endif

            ! -- linear detrainment from vent to base of umbrella
            do ko=1,k_initial-1
              vert_mass_dist(ko)=float(ko)/float(k_initial-1)
            enddo
            x1=sum(vert_mass_dist(1:k_initial-1))

            do ko=1,k_initial-1
              vert_mass_dist(ko)=(1.-percen_mass_umbrel)*vert_mass_dist(ko)/x1
            enddo
            ! -- tr emissions for umbrella (explosive) type emissons
            do ko=1,k_final
              emis_ant(i,ko,j,p_e_tr1)=vert_mass_dist(ko)*emiss_tr_mass(j)
              emis_ant(i,ko,j,p_e_tr2)=1./float(k_final)*emiss_tr_mass(j)
            enddo
            if (emiss_tr_dt(j) <= 0.) cycle
            do k=kts,kte-2
              factor2=dtstep*rri(i,k,j)/dz8w(i,k,j)
              chem(i,k,j,p_tr2)=chem(i,k,j,p_tr2)+emis_ant(i,k,j,p_e_tr2)*factor2
              if(real_time > 360.) chem(i,k,j,p_tr1)=chem(i,k,j,p_tr1)+emis_ant(i,k,j,p_e_tr2)*factor2
            enddo
          enddo 
        enddo       
      endif       
      print *,'chem_prep: stage end'

      call mpp_clock_end(clockId)

  end subroutine chem_prep

! -------------------------------------------------------------------------------

  real function calday(year, month, day, hour, minute, sec)
! For time interpolation; Julian day (0 to 365 for non-leap year)
! input:
    integer, intent(in):: year, month, day, hour
    integer, intent(in):: minute, sec
! Local:
      integer n, m, ds, nday
      real tsec
      integer days(12)
      data days /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

      ds = day - 1

      if( month /= 1 ) then
          do m=1, month-1
            if( m==2  .and. leap_year(year) ) then
                ds = ds + 29
            else
                ds = ds + days(m)
            endif
          enddo
      endif

      if ( leap_year(year_track_data) ) then
           nday = 366
      else
           nday = 365
      endif

      calday = real((year-year_track_data)*nday + ds)  + real(hour*3600 + minute*60 + sec)/86400.

  end function calday


  logical function leap_year(ny)
  integer, intent(in):: ny
!
! Determine if year ny is a leap year
! Author: S.-J. Lin
   integer ny00
!
! No leap years prior to 0000
!
      parameter ( ny00 = 0000 )   ! The threshold for starting leap-year 

      if( ny >= ny00 ) then
         if( mod(ny,100) == 0. .and. mod(ny,400) == 0. ) then
             leap_year = .true.
         elseif( mod(ny,4) == 0. .and. mod(ny,100) /= 0.  ) then
             leap_year = .true.
         else
             leap_year = .false.
         endif
      else
          leap_year = .false.
      endif

  end function leap_year

  SUBROUTINE szangle(imx, jmx, doy, xhour, sza, cossza,xlon,rlat)

!
! ****************************************************************************
! **                                                                        **
! **  This subroutine computes solar zenith angle (SZA):                    **
! **                                                                        **
! **      cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR)         **
! **                                                                        **
! **  where LAT is the latitude angle, DEC is the solar declination angle,  **
! **  and AHR is the hour angle, all in radius.                             **
! **                                                                        **
! **  DOY = day-of-year, XHOUR = UT time (hrs).                             **
! **  XLON = longitude in degree, RLAT = latitude in radian.                **
! ****************************************************************************
!

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: imx, jmx
  INTEGER, INTENT(IN)    :: doy
  REAL,    INTENT(IN)    :: xhour
  REAL,    INTENT(OUT)   :: sza(imx,jmx), cossza(imx,jmx)

  REAL    :: a0, a1, a2, a3, b1, b2, b3, r, dec, timloc, ahr,xlon,rlat
  real, parameter :: pi=3.14
  INTEGER :: i, j

  ! executable statements

  ! ***************************************************************************
  ! *  Solar declination angle:                                               *
  ! ***************************************************************************
  a0 = 0.006918
  a1 = 0.399912
  a2 = 0.006758
  a3 = 0.002697
  b1 = 0.070257
  b2 = 0.000907
  b3 = 0.000148
  r  = 2.0* pi * REAL(doy-1)/365.0
  !
  dec = a0 - a1*COS(  r)   + b1*SIN(  r)   &
           - a2*COS(2.0*r) + b2*SIN(2.0*r) &
           - a3*COS(3.0*r) + b3*SIN(3.0*r)
  !
  DO i = 1,imx
     ! ************************************************************************
     ! *  Hour angle (AHR) is a function of longitude.  AHR is zero at        *
     ! *  solar noon, and increases by 15 deg for every hour before or        *
     ! *  after solar noon.                                                   *
     ! ************************************************************************
     ! -- Local time in hours
     timloc  = xhour + xlon/15.0
     !      IF (timloc < 0.0) timloc = 24.0 + timloc
     IF (timloc > 24.0) timloc = timloc - 24.0
     !
     ! -- Hour angle
     ahr = ABS(timloc - 12.0) * 15.0 * pi/180.0

     DO j = 1,jmx
        ! -- Solar zenith angle      
        cossza(i,j) = SIN(rlat) * SIN(dec) + &
                      COS(rlat) * COS(dec) * COS(ahr)
        sza(i,j)    = ACOS(cossza(i,j)) * 180.0/pi
        IF (cossza(i,j) < 0.0)   cossza(i,j) = 0.0
        !
     END do
  END DO

  END SUBROUTINE szangle

end module chem_shr_mod
