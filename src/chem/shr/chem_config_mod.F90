module chem_config_mod

  use chem_types_mod, only : CHEM_MAXSTR

  implicit none

  ! -- Currently available modules
  integer, parameter :: CHEM_OPT_GOCART       = 300
  integer, parameter :: CHEM_OPT_GOCART_RACM  = 301
  integer, parameter :: CHEM_OPT_RACM_SOA_VBS = 108

  integer, parameter :: DUST_OPT_GOCART = 1
  integer, parameter :: DUST_OPT_AFWA   = 3


  ! -- data structure for configuration options

  type chem_config_type
    sequence
    character(len=CHEM_MAXSTR) :: emi_inname
    character(len=CHEM_MAXSTR) :: fireemi_inname
    character(len=CHEM_MAXSTR) :: emi_outname
    character(len=CHEM_MAXSTR) :: fireemi_outname
    character(len=CHEM_MAXSTR) :: input_chem_inname
    character(len=CHEM_MAXSTR) :: input_chem_outname
    character(len=CHEM_MAXSTR) :: chem_hist_outname = 'fim_out_'
    integer :: io_style_emissions
    real    :: bioemdt
    real    :: photdt
    real    :: chemdt
!   integer :: ne_area
    integer :: kemit
!   integer :: nmegan
    integer :: kfuture
    integer :: chem_conv_tr
    integer :: chem_opt
    integer :: gaschem_onoff
    integer :: aerchem_onoff
    integer :: wetscav_onoff
    integer :: cldchem_onoff
    integer :: vertmix_onoff
    integer :: chem_in_opt
    integer :: phot_opt
    integer :: drydep_opt
    real    :: depo_fact
    integer :: emiss_opt
    integer :: dust_opt
    integer :: dmsemis_opt
    integer :: seas_opt
    integer :: bio_emiss_opt
    integer :: biomass_burn_opt
    integer :: plumerisefire_frq
    integer :: emiss_inpt_opt
    integer :: gas_bc_opt
    integer :: gas_ic_opt
    integer :: aer_bc_opt
    integer :: aer_ic_opt
    logical :: have_bcs_chem
    integer :: aer_ra_feedback    = 0
    integer :: aer_op_opt
    integer :: conv_tr_aqchem
    integer :: call_biomass       = 1
    integer :: call_chemistry     = 1
    integer :: call_radiation     = 1
    logical :: readrestart        = .false.
    integer :: archive_step
  end type chem_config_type

  character(len=CHEM_MAXSTR) :: emi_inname
  character(len=CHEM_MAXSTR) :: fireemi_inname
  character(len=CHEM_MAXSTR) :: emi_outname
  character(len=CHEM_MAXSTR) :: fireemi_outname
  character(len=CHEM_MAXSTR) :: input_chem_inname
  character(len=CHEM_MAXSTR) :: input_chem_outname
  character(len=CHEM_MAXSTR) :: chem_hist_outname
  integer :: io_style_emissions
  real    :: bioemdt
  real    :: photdt
  real    :: chemdt
! integer :: ne_area
  integer :: kemit
! integer :: nmegan
  integer :: kfuture
  integer :: chem_conv_tr
  integer :: chem_opt
  integer :: gaschem_onoff
  integer :: aerchem_onoff
  integer :: wetscav_onoff
  integer :: cldchem_onoff
  integer :: vertmix_onoff
  integer :: chem_in_opt
  integer :: phot_opt
  integer :: drydep_opt
  real    :: depo_fact
  integer :: emiss_opt
  integer :: dust_opt
  integer :: dmsemis_opt
  integer :: seas_opt
  integer :: bio_emiss_opt
  integer :: biomass_burn_opt
  integer :: plumerisefire_frq
  integer :: emiss_inpt_opt
  integer :: gas_bc_opt
  integer :: gas_ic_opt
  integer :: aer_bc_opt
  integer :: aer_ic_opt
  logical :: have_bcs_chem
  integer :: aer_ra_feedback
  integer :: aer_op_opt
  integer :: conv_tr_aqchem
  integer :: archive_step

  integer :: numphr                        ! # of time steps/hr
  integer :: nv_g 

  real :: ash_mass
  real :: ash_height


  namelist /chem_nml/          &
    emi_inname,                &
    fireemi_inname,            &
    emi_outname,               &
    fireemi_outname,           &
    input_chem_inname,         &
    input_chem_outname,        &
    chem_hist_outname,         &
    io_style_emissions,        &
!   bioemdt,                   &
    photdt,                    &
    chemdt,                    &
!   ne_area,                   &
    kemit,                     &
!   nmegan,                    &
    kfuture,                   &
    chem_conv_tr,              &
    chem_opt,                  &
    gaschem_onoff,             &
    aerchem_onoff,             &
    wetscav_onoff,             &
    cldchem_onoff,             &
    vertmix_onoff,             &
    chem_in_opt,               &
    phot_opt,                  &
    drydep_opt,                &
    depo_fact,                 &
    emiss_opt,                 &
    dust_opt,                  &
    dmsemis_opt,               &
    seas_opt,                  &
    bio_emiss_opt,             &
    biomass_burn_opt,          &
    plumerisefire_frq,         &
    emiss_inpt_opt,            &
    gas_bc_opt,                &
    gas_ic_opt,                &
    aer_bc_opt,                &
    aer_ic_opt,                &
    have_bcs_chem,             &
    aer_ra_feedback,           &
    aer_op_opt,                &
    conv_tr_aqchem,            &
    archive_step

  character(len=*), parameter :: chem_file_nml = 'input.nml'


  type(chem_config_type) :: chem_config

  ! -- control variables
  integer :: ntra                =  4      ! # of tracers advected on small dt: 1=theta 2=qv 3=qw 4=O3
  integer :: ntrb                =  0      ! # of tracers advected on large dt: will include chemistry
! integer :: nip                 =  0      ! # of icosaedral cells
! integer :: nvl                   ! Number of vertical native levels
! integer :: nvlp                =  0      ! # of isobaric vertical levels - ex.  1000-25 hPa
! integer :: nvlp1                         ! # of vertical levels ( = layers+1)
  real    :: dt                  = 10      ! model time step (seconds)

  ! -- GFS physics _ gocart very light for fim
  !integer, parameter :: num_moist=2+1
  !integer, parameter :: num_chem=13
  !integer, parameter :: num_emis_ant = 6
  !
  ! following for Lin et al. + regular GOCART
  !
  !integer, parameter :: num_moist=6+1
  !integer, parameter :: num_chem=18
  !integer, parameter :: num_emis_ant = 6
  !integer, parameter :: num_emis_vol = 0
  !
  ! volcanic ash only
  !integer, parameter :: num_moist=2+1
  !integer, parameter :: num_chem=23
  !integer, parameter :: num_emis_ant = 6
  !integer, parameter :: num_emis_vol = 10
  !
  ! light gocart + reduced volcanic ash only (4 size bins)
  !!!!!!!!!! REG TEST SETUP !!!!
  !integer, parameter :: num_moist=2+1
  !integer, parameter :: num_chem=17
  !integer, parameter :: num_emis_ant = 6
  !integer, parameter :: num_emis_vol = 4
  integer            :: num_chem,num_moist,num_emis_ant,   &
                        num_emis_vol,num_ebu,num_ebu_in,numgas
  !integer            :: num_moist =0
  !integer            :: num_emis_ant =0
  !integer            :: num_emis_vol =0
  ! volcanic ash only (4 size bins)
  !integer, parameter :: num_moist=2+1
  !integer, parameter :: num_chem=4
  !integer, parameter :: num_emis_ant = 0
  !integer, parameter :: num_emis_vol = 4
  ! Pure GOCART (volcanic ash included in p25 and p10)
  !!!!! CURRENT REAL-TIME
  !integer, parameter :: num_moist=2+1
  !integer, parameter :: num_chem=19
  !integer, parameter :: num_emis_ant = 6
  !integer, parameter :: num_emis_vol = 4
  ! Pure GOCART (volcanic ash included in p25 and p10) + mp_phys=4
  ! (qv,qc,qr,qi,qs)
  !integer, parameter :: num_moist=4+1
  !integer, parameter :: num_chem=19
  !integer, parameter :: num_emis_ant = 6
  !integer, parameter :: num_emis_vol = 4
  !
  !
  integer            :: num_emis_season_bb=0
  integer            :: num_emis_season_ant=0
  integer, parameter :: nbands=14
  integer, parameter :: nbandlw=16
  integer, parameter :: num_soil_layers=4
  integer, parameter :: num_scalar=1
  integer, parameter :: nvl_gocart=55  ! number of input levels from gocart file
  integer, parameter :: num_ext_coef = 5
  integer, parameter :: num_bscat_coef = 3
  integer, parameter :: num_asym_par = 3
  INTEGER, PARAMETER          :: num_emis_dust = 5
  INTEGER, PARAMETER          :: num_emis_seas = 4
  INTEGER, PARAMETER          :: ne_area = 41
  INTEGER, PARAMETER          :: nmegan = 1

  ! -- configuration variables for output:
  !  . grid level defined in FIM Makefile
  integer :: glvl
  !  . time stamp
  character(len=12) :: yyyymmddhhmm = ""

  logical :: is_nuopc = .false.

  public
  
end module chem_config_mod
