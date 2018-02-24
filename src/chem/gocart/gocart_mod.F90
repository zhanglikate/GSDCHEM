module gocart_mod

  use chem_rc_mod
  use chem_types_mod
  use chem_const_mod,  only : cp, grvity, mwdry, p1000, rd, epsilc
  use chem_tracers_mod
  use chem_config_mod, only : CHEM_OPT_NONE,   &
                              CHEM_OPT_GOCART, &
                              DUST_OPT_GOCART, &
                              DUST_OPT_AFWA

  use gocart_prep_mod
  use gocart_settling_mod
  use gocart_aerosols_mod
  use gocart_dmsemis_mod
  use gocart_chem_mod
  use plume_rise_mod
  use vash_settling_mod
  use dep_mod
  use dust_mod
  use seas_mod
  use opt_mod

  implicit none

  public

contains

  subroutine gocart_init
  end subroutine gocart_init

  subroutine gocart_advance(readrestart, chem_opt, chem_in_opt, &
    biomass_burn_opt, seas_opt, dust_opt, dmsemis_opt, aer_ra_feedback, &
    call_biomass, call_chemistry, call_rad, &
    kemit, ktau, dts, current_month, tz, julday,      &
    p_gocart, clayfrac, dm0, emiss_ab, emiss_abu,                         &
    emiss_ash_dt, emiss_ash_height, emiss_ash_mass, &
    emiss_tr_dt, emiss_tr_height, emiss_tr_mass, ero1, ero2, ero3,     &
    h2o2_backgd, no3_backgd, oh_backgd,  plumestuff, sandfrac, th_pvsrf,  &
    area_in, hf2d_in, pb2d_in, rc2d_in, rn2d_in, rsds_in, slmsk2d_in, snwdph2d_in, stype2d_in,       &
    ts2d_in, us2d_in, vtype2d_in, vfrac2d_in, zorl2d_in, exch_in, ph3d_in, phl3d_in, pr3d_in, prl3d_in, &
    sm3d_in, tk3d_in, us3d_in, vs3d_in, ws3d_in, tr3d_in, tr3d_out, tr3d, trdp, &
    emi_d1, emi_d2, emi_d3, emi_d4, emi_d5, ext_cof, sscal, asymp, aod2d, &
    p10, pm25, ebu_oc, oh_bg, h2o2_bg, no3_bg, wet_dep, &
    nvl, nvi, ntra, ntrb, nvl_gocart, nbands, numgas, num_ebu, num_ebu_in, num_soil_layers, &
    num_chem, num_moist, num_emis_vol, num_emis_ant, num_emis_dust, num_emis_seas, &
    num_asym_par, num_bscat_coef, num_ext_coef, lon, lat, &
    its, ite, jts, jte, kts, kte, &
    ims, ime, jms, jme, kms, kme, rc)

    logical,            intent(in) :: readrestart
    integer,            intent(in) :: chem_opt
    integer,            intent(in) :: chem_in_opt
    integer,            intent(in) :: biomass_burn_opt
    integer,            intent(in) :: seas_opt
    integer,            intent(in) :: dust_opt
    integer,            intent(in) :: dmsemis_opt
    integer,            intent(in) :: aer_ra_feedback
    integer,            intent(in) :: call_biomass
    integer,            intent(in) :: call_chemistry
    integer,            intent(in) :: call_rad
    integer,            intent(in) :: kemit
    integer,            intent(in) :: ktau
    integer,            intent(in) :: current_month
    integer,            intent(in) :: julday
    integer,            intent(in) :: tz
    integer,            intent(in) :: its, ite, jts, jte, kts, kte
    integer,            intent(in) :: ims, ime, jms, jme, kms, kme
    integer,            intent(in) :: nvl, nvi, ntra, ntrb, nvl_gocart, &
                                      nbands, num_ebu, num_ebu_in, num_soil_layers, &
                                      num_chem, num_moist, num_emis_vol, num_emis_ant, &
                                      num_emis_dust, num_emis_seas
    integer,            intent(in) :: num_asym_par, num_bscat_coef, num_ext_coef
    integer,            intent(in) :: numgas
    integer, optional, intent(out) :: rc

    real(CHEM_KIND_R8), intent(in) :: dts

    real(CHEM_KIND_R4), dimension(nvl_gocart+1),     intent(in) :: p_gocart
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: dm0
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: ero1
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: ero2
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: ero3
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: emiss_tr_dt
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: emiss_tr_height
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: emiss_tr_mass
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: emiss_ash_dt
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(inout) :: emiss_ash_height
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(inout) :: emiss_ash_mass
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: clayfrac
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: sandfrac
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: th_pvsrf
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl_gocart),     intent(in) :: h2o2_backgd
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl_gocart),     intent(in) :: no3_backgd
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl_gocart),     intent(in) :: oh_backgd
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, num_emis_ant),   intent(in) :: emiss_ab
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, num_ebu_in),     intent(in) :: emiss_abu
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 8), intent(in) :: plumestuff

    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme), intent(in) :: area_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme), intent(in) :: hf2d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme), intent(in) :: pb2d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme), intent(in) :: rc2d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme), intent(in) :: rn2d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme), intent(in) :: rsds_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme), intent(in) :: slmsk2d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme), intent(in) :: snwdph2d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme), intent(in) :: stype2d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme), intent(in) :: ts2d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme), intent(in) :: us2d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme), intent(in) :: vtype2d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme), intent(in) :: vfrac2d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme), intent(in) :: zorl2d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme), intent(in) :: lat
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme), intent(in) :: lon

    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme, nvl), intent(in) :: exch_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme, nvi), intent(in) :: ph3d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme, nvl), intent(in) :: phl3d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme, nvi), intent(in) :: pr3d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme, nvl), intent(in) :: prl3d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme, num_soil_layers), intent(in) :: sm3d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme, nvl), intent(in) :: tk3d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme, nvl), intent(in) :: us3d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme, nvl), intent(in) :: vs3d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme, nvl), intent(in) :: ws3d_in

    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme, nvl, ntra+ntrb), intent(in)  :: tr3d_in
    real(CHEM_KIND_R8), dimension(ims:ime, jms:jme, nvl, ntra+ntrb), intent(out) :: tr3d_out
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl, ntra+ntrb), intent(out) :: tr3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl, ntra+ntrb), intent(inout) :: trdp

    ! -- output tracers
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: emi_d1
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: emi_d2
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: emi_d3
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: emi_d4
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: emi_d5
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: aod2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:num_chem), intent(out) :: wet_dep
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl), intent(out) :: p10
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl), intent(out) :: pm25
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl), intent(out) :: ebu_oc
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl), intent(out) :: oh_bg
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl), intent(out) :: h2o2_bg
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl), intent(out) :: no3_bg
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl, 1:nbands), intent(out) :: ext_cof
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl, 1:nbands), intent(out) :: sscal
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:nvl, 1:nbands), intent(out) :: asymp

    ! -- local variables

    ! -- work arrays
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: aod
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: ash_fall
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: clayf
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: cu_co_ten
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: dep_vel_o3
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: e_co
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: dms_0
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: dusthelp
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: dxy
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: firesize_agef
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: firesize_aggr
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: firesize_agsv
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: firesize_agtf
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: gsw
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: hfx
    integer,            dimension(ims:ime, jms:jme) :: isltyp
    integer,            dimension(ims:ime, jms:jme) :: ivgtyp
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: mean_fct_agef
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: mean_fct_aggr
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: mean_fct_agsv
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: mean_fct_agtf
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: pbl
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: raincv_b
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: rcav
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: rnav
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: rmol
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: sandf
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: seashelp
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: snowh
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: tcosz
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: tsk
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: ttday
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: u10
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: ust
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: v10
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: vegfra
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: xland
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: xlat
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: xlong
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: znt
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: deg_lat
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: deg_lon

    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, num_ebu_in) :: ebu_in
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:3)        :: erod

    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: ac3
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: ahno3
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: anh3
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: asulf
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: backg_h2o2
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: backg_no3
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: backg_oh
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: cor3
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: dz8w
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: exch_h
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: h2o2_t
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: h2oai
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: h2oaj
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: no3_t
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: nu3
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: oh_t
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: p8w
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: p_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: pm10
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: pm2_5_dry
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: pm2_5_dry_ec
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: relhum
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: rho_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: rri
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: t8w
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: t_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: u_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: v_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: vvel
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: z_at_w
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: zmid
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme) :: convfac
    real(CHEM_KIND_R4), dimension(ims:ime, 1:num_soil_layers, jms:jme) :: smois


    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:num_asym_par)   :: asym_par
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:nbands)         :: asympar
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:num_bscat_coef) :: bscat_coeff
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: bscoefsw
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:num_chem)       :: chem
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:num_ebu)        :: ebu
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kemit, jms:jme, 1:num_emis_ant) :: emis_ant
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:num_emis_vol)   :: emis_vol
    real(CHEM_KIND_R4), dimension(ims:ime,   1:1  , jms:jme, 1:num_emis_dust)  :: emis_dust
    real(CHEM_KIND_R4), dimension(ims:ime,   1:1  , jms:jme, 1:num_emis_seas)  :: emis_seas
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:num_ext_coef)   :: ext_coeff
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:nbands)         :: extt
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: gaersw
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: l2aer
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: l3aer
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: l4aer
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: l5aer
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: l6aer
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: l7aer
    real(CHEM_KIND_R4), dimension(ims:ime,   1:1  , jms:jme, 1:5)              :: srce_dust
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:nbands)         :: ssca

    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:num_chem)       :: var_rmv
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 1:num_chem)       :: tr_fall
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:num_moist)      :: moist
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:16)             :: tauaerlw
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: tauaersw
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, 1:4)              :: waersw

    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl) :: tk3d

    ! -----------------------------------------------------------------------
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: area
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: hf2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: pb2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: rc2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: rn2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: rsds
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: slmsk2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: snwdph2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: stype2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: ts2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: us2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: vtype2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: vfrac2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme) :: zorl2d

    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl) :: exch
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvi) :: ph3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl) :: phl3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvi) :: pr3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl) :: prl3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, num_soil_layers) :: sm3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl) :: us3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl) :: vs3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl) :: ws3d
!   real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl) :: dp3d

    ! -----------------------------------------------------------------------
    ! -- local variables
    logical, save :: firstfire = .true.
    logical :: call_gocart, call_plume, call_radiation, scale_fire_emiss
    logical :: store_arrays

!   integer :: ktau
!   integer :: julday
    integer :: localrc
    integer :: nbegin, nv, nvv
    integer :: i, j, jp, jps, k
    integer :: ids, ide, jds, jde, kds, kde
!   integer :: current_month, current_gmt, current_secs, current_msecs
    real(CHEM_KIND_R4) :: dt
    real(CHEM_KIND_R8) :: curr_secs

    real(CHEM_KIND_R4) :: factor, factor2
    real(CHEM_KIND_R4) :: dtstep, gmt
    real(CHEM_KIND_R4) :: dust_alpha,dust_gamma

    ! -- begin
    print *,'gocart_run: entering ...', chem_opt, ims, ime, jms, jme, num_chem

    if (present(rc)) rc = CHEM_RC_SUCCESS

    if (chem_opt == CHEM_OPT_NONE) return

    ! -----------------------------------------------------------------------
    area = area_in
    hf2d = hf2d_in
    pb2d = pb2d_in
    rc2d = rc2d_in
    rn2d = rn2d_in
    rsds = rsds_in
    slmsk2d = slmsk2d_in
    snwdph2d = snwdph2d_in
    stype2d  = stype2d_in
    ts2d = ts2d_in
    us2d = us2d_in
    vtype2d = vtype2d_in
    vfrac2d = vfrac2d_in
    zorl2d  = zorl2d_in

    exch = exch_in
    ph3d = ph3d_in
    phl3d = phl3d_in
    pr3d  = pr3d_in
    prl3d = prl3d_in
    sm3d  = sm3d_in
    us3d  = us3d_in
    vs3d  = vs3d_in
    ws3d  = ws3d_in
    ! -----------------------------------------------------------------------

    ! -- set domain
    ids = ims
    ide = ime
    jds = jms
    jde = jme
    kds = kms
    kde = kme

    ! -- initialize local arrays
    print *,'gocart_run: initializing local arrays ...'
    dep_vel_o3 = 0.
    e_co       = 0.
    raincv_b   = 0.
    cu_co_ten  = 0.
    dusthelp   = 0.
    seashelp   = 0.
    var_rmv    = 0.
    tr_fall    = 0.
    deg_lon    = real(lon,     CHEM_KIND_R4)
    deg_lat    = real(lat,     CHEM_KIND_R4)
    tk3d       = real(tk3d_in, CHEM_KIND_R4)
    tr3d       = real(tr3d_in, CHEM_KIND_R4)

!   dp3d(:,:,1:nvl) = pr3d(:,:,1:nvl) - pr3d(:,:,2:nvi)

    print *,'gocart_run: longitude min/max = ',minval(deg_lon),maxval(deg_lon)
    print *,'gocart_run: local array initialized ...'
    ! -- set run parameters

    ! -- nbegin is the start address (-1) of the first chem variable in tr3d
!   nbegin = ntra + num_moist - 3 + min(max(num_moist-3,0),1)
    if (num_moist > 3) then
      nbegin = ntra + num_moist - 2
    else
      nbegin = ntra + num_moist - 3
    end if

    ! -- set numerical parameters
    dust_alpha = 1.0
    dust_gamma = 1.6
!   xlv = 2.5e+06

    ! -- get time & time step
    dt = real(dts, kind=CHEM_KIND_R4)
    curr_secs = ktau * dts
    gmt = real(tz)

    ! -- set control flags
    call_plume       = (biomass_burn_opt > 0) .and. &
                      ((mod(ktau, call_biomass  ) == 0) .or. (ktau == 1) .or. firstfire)
    call_gocart      = (mod(ktau, call_chemistry) == 0) .or. (ktau == 1)
    call_radiation   = (mod(ktau, call_rad) == 0) .or. (ktau == 1)
    scale_fire_emiss = .false.
    print *,'gocart_run: control flags set'

    print *,'gocart_run: set control flags ...'
    print *,'gocart_run: set control flags ...', biomass_burn_opt
    print *,'gocart_run: set control flags ...', call_biomass
    print *,'gocart_run: set control flags ...', call_chemistry
    print *,'gocart_run: set control flags ...', call_radiation
    print *,'gocart_run: set control flags ...', ktau, firstfire

    print *,'gocart_run: get domain bounds ...', its, ite, jts, jte
    ! -- start working
    if (ktau <= 1) then
      dtstep = dt
      rcav = rc2d
      rnav = rn2d-rc2d
    else
      dtstep = call_chemistry * dt
      rcav = max(0.,rc2d-rcav)
      rnav = max(0.,rn2d-rc2d-rnav)
    end if

    ! -- add time interval to initial date (to be replaced by ESMF_Clock
    ! -- this may not be needed

    ! -- get ready for chemistry run
    ! -- ttday and tcosz are computed -- should they be replaced with FV3 quantities?
    ! -- gmt is needed to compute quantities above
    ! -- get julday ----

    print *,'gocart_run: entering gocart_prep ...'
    call gocart_prep(readrestart,chem_opt,chem_in_opt,ktau,dt,tr3d,tk3d,sm3d,   &
                   ts2d,us2d,rsds,pr3d,prl3d,ph3d,phl3d,emiss_ash_mass,emiss_ash_height, &
                   emiss_ash_dt,dm0,emiss_tr_mass,emiss_tr_height,      &
                   emiss_tr_dt,snwdph2d,VFRAC2d,VTYPE2d,STYPE2d,us3d,vs3d,ws3d,  &
                   slmsk2d,zorl2d,exch,pb2d,hf2d,clayfrac,clayf,sandfrac,sandf,th_pvsrf, &
                   oh_backgd,h2o2_backgd,no3_backgd,backg_oh,backg_h2o2,backg_no3,p_gocart,   &
                   nvl_gocart,ttday,tcosz,gmt,julday,area,ero1,   &
                   ero2,ero3,rcav,raincv_b,deg_lat,deg_lon,nvl,nvi,ntra, &
                   relhum,rri,t_phy,moist,u_phy,v_phy,p_phy,chem,tsk,ntrb, &
                   grvity,rd,p1000,cp,erod,emis_ant,emis_vol,e_co,dms_0,        &
                   u10,v10,ivgtyp,isltyp,gsw,vegfra,rmol,ust,znt,xland,dxy, &
                   t8w,p8w,exch_h,pbl,hfx,snowh,xlat,xlong,convfac,z_at_w,zmid,dz8w,vvel,&
                   rho_phy,smois,num_soil_layers,num_chem,num_moist,        &
                   emiss_abu,ebu_in,emiss_ab,num_ebu_in,num_emis_ant,       &
                   num_emis_vol,kemit,call_gocart,plumestuff, &
                   mean_fct_agtf,mean_fct_agef,mean_fct_agsv, &
                   mean_fct_aggr,firesize_agtf,firesize_agef, &
                   firesize_agsv,firesize_aggr, &
                   ids,ide, jds,jde, kds,kde, &
                   ims,ime, jms,jme, kms,kme, &
                   its,ite, jts,jte, kts,kte, rc=localrc)
    if (chem_rc_check(localrc, msg="Failure in gocart_prep", &
      file=__FILE__, line=__LINE__, rc=rc)) return

    ! -- compute sea salt
    print *,'gocart_run: entering sea salt ...'
    if (seas_opt == SEAS_OPT_DEFAULT) then
      call gocart_seasalt_driver(ktau,dt,rri,t_phy,moist, &
        u_phy,v_phy,chem,rho_phy,dz8w,u10,v10,p8w,        &
        xland,xlat,xlong,dxy,grvity,emis_seas,           &
        seashelp,num_emis_seas,num_moist,num_chem,        &
        ids,ide, jds,jde, kds,kde,                        &
        ims,ime, jms,jme, kms,kme,                        &
        its,ite, jts,jte, kts,kte)
    endif
    print *,'gocart_run: exit sea salt ...'

    print *,'gocart_run: check dust ...', dust_opt
    store_arrays = .false.
    select case (dust_opt)
      case (DUST_OPT_GOCART)
    print *,'gocart_run: dust is GOCART ...'
        call gocart_dust_driver(chem_opt,ktau,dt,rri,t_phy,moist,u_phy,&
          v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,ivgtyp,isltyp,&
          vegfra,xland,xlat,xlong,gsw,dxy,grvity,emis_dust,srce_dust, &
          dusthelp,num_emis_dust,num_moist,num_chem,num_soil_layers,   &
          current_month,                                               &
          ids,ide, jds,jde, kds,kde,                                   &
          ims,ime, jms,jme, kms,kme,                                   &
          its,ite, jts,jte, kts,kte)
        store_arrays = .true.
      case (DUST_OPT_AFWA)
    print *,'gocart_run: dust is AFWA ...'
        call gocart_dust_afwa_driver(ktau,dt,rri,t_phy,moist,u_phy,    &
          v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,ivgtyp,isltyp,&
          vegfra,xland,xlat,xlong,gsw,dxy,grvity,emis_dust,srce_dust, &
          dusthelp,ust,znt,clayf,sandf,dust_alpha,dust_gamma,    &
          num_emis_dust,num_moist,num_chem,num_soil_layers,            &
          ids,ide, jds,jde, kds,kde,                                   &
          ims,ime, jms,jme, kms,kme,                                   &
          its,ite, jts,jte, kts,kte)
        store_arrays = .true.
    end select
    store_arrays = store_arrays .and. (chem_opt >= CHEM_OPT_GOCART)
    print *,'gocart_run: done dust ...'

    ! -- set output arrays
    if (store_arrays) then
      print *,'gocart_run: storing arrays ...'
      emi_d1(its:ite,jts:jte) = emis_dust(its:ite,1,jts:jte,1)
      emi_d2(its:ite,jts:jte) = emis_dust(its:ite,1,jts:jte,2)
      emi_d3(its:ite,jts:jte) = emis_dust(its:ite,1,jts:jte,3)
      emi_d4(its:ite,jts:jte) = emis_dust(its:ite,1,jts:jte,4)
      emi_d5(its:ite,jts:jte) = emis_dust(its:ite,1,jts:jte,5)
      print *,'gocart_run: storing arrays done'
    end if

    if (call_plume) then
      print *,'gocart_run: calling plume'
      firstfire = .false.
      call plumerise_driver (ktau,dtstep,num_chem,num_moist,num_ebu, &
        num_ebu_in,ebu,ebu_in,mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr, &
        firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr, &
        'GOCART','BIOMASSB', t_phy,moist, &
        rho_phy,vvel,u_phy,v_phy,p_phy,                          &
        z_at_w,scale_fire_emiss,                                 &
        ids,ide, jds,jde, kds,kde,                               &
        ims,ime, jms,jme, kms,kme,                               &
        its,ite, jts,jte, kts,kte                                )
      print *,'gocart_run: calling plume done'
    end if

    if (dmsemis_opt == 1) then
      print *,'gocart_run: calling dmsemis ...'
      call gocart_dmsemis(dt,rri,t_phy,u_phy,v_phy,     &
         chem,rho_phy,dz8w,u10,v10,p8w,dms_0,tsk,       &
         ivgtyp,isltyp,xland,dxy,grvity,mwdry,         &
         num_chem,p_dms,                                &
         ids,ide, jds,jde, kds,kde,                     &
         ims,ime, jms,jme, kms,kme,                     &
         its,ite, jts,jte, kts,kte)
      print *,'gocart_run: calling dmsemis done'
    endif

    if ((dust_opt == DUST_OPT_GOCART) .or. &
        (dust_opt == DUST_OPT_AFWA  ) .or. &
        (seas_opt == SEAS_OPT_DEFAULT)) then
      print *,'gocart_run: calling settling ...'
      call gocart_settling_driver(dt,t_phy,moist,  &
        chem,rho_phy,dz8w,p8w,p_phy,   &
        dusthelp,seashelp,dxy,grvity, &
        num_moist,num_chem,            &
        ids,ide, jds,jde, kds,kde,     &
        ims,ime, jms,jme, kms,kme,     &
        its,ite, jts,jte, kts,kte)
      print *,'gocart_run: calling settling done'
    end if

#if 0
    if (chem_opt == 316) then
      ! -- 10 volcanic size bins
      call vash_settling_driver(dt,t_phy,moist, &
         chem,rho_phy,dz8w,p8w,p_phy,dxy,      &
         ash_fall,grvity,num_moist,num_chem,    &
         ids,ide, jds,jde, kds,kde,             &
         ims,ime, jms,jme, kms,kme,             &
         its,ite, jts,jte, kts,kte)
      ashfall = ash_fall
    else
#endif
      ! -- 4 volcanic size bins
      print *,'gocart_run: calling vashshort_settling_driver ...'
      call vashshort_settling_driver(dt,t_phy,moist, &
           chem,rho_phy,dz8w,p8w,p_phy,dxy,         &
           ash_fall,grvity,num_moist,num_chem,       &
           ids,ide, jds,jde, kds,kde,                &
           ims,ime, jms,jme, kms,kme,                &
           its,ite, jts,jte, kts,kte) 
      print *,'gocart_run: calling vashshort_settling_driver done'
!     ashfall = ash_fall !!!!!!!!!!!!!!!!!! WARNING: rewrites ashfall if chem_opt = 316
      print *,'gocart_run: ashfall done'
#if 0
    end if
#endif

    ! -- add biomass burning emissions at every timestep
    if (biomass_burn_opt == 1) then
      print *,'gocart_run: set chem: biomass_burn_opt ...'
      do j = jts, jte
        do k = kts, kte
          do i = its, ite
            ! -- factor for pm emissions, factor2 for burn emissions
            factor  = dt*rri(i,k,j)/dz8w(i,k,j)
            factor2 = 4.828e-4*dt*rri(i,k,j)/(60.*dz8w(i,k,j))
            chem(i,k,j,p_oc1) = chem(i,k,j,p_oc1) + factor  * ebu(i,k,j,p_ebu_oc  )
            chem(i,k,j,p_bc1) = chem(i,k,j,p_bc1) + factor  * ebu(i,k,j,p_ebu_bc  )
            chem(i,k,j,p_p25) = chem(i,k,j,p_p25) + factor  * ebu(i,k,j,p_ebu_pm25)
            chem(i,k,j,p_p10) = chem(i,k,j,p_p10) + factor  * ebu(i,k,j,p_ebu_pm10)
            chem(i,k,j,p_so2) = chem(i,k,j,p_so2) + factor2 * ebu(i,k,j,p_ebu_so2 )
          end do
        end do
      end do
      print *,'gocart_run: set chem: done'
    end if

#if 0
    ! -- subgrid convective transport
    if (chem_conv_tr == 2 )then
      call grelldrvct(DT,ktau,             &
        rho_phy,RAINCV_b,chem,tr_fall,     &
        U_phy,V_phy,t_phy,moist,dz8w,p_phy,&
        XLV,CP,grvity,rv,z_at_w,cu_co_ten, &
        numgas,chem_opt,                   &
        num_chem,num_moist,                &
        ids,ide, jds,jde, kds,kde,         &
        ims,ime, jms,jme, kms,kme,         &
        its,ite, jts,jte, kts,kte)
     endif
#endif

     print *,'gocart_run: calling dry_dep_driver ...'
     call dry_dep_driver(ktau,dt,julday,current_month,t_phy,p_phy,&
       moist,p8w,rmol,rri,gmt,t8w,rcav,                           &
       chem,rho_phy,dz8w,exch_h,hfx,                              &
       ivgtyp,tsk,gsw,vegfra,pbl,ust,znt,zmid,z_at_w,             &
       xland,xlat,xlong,h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,     &
       anh3,dep_vel_o3,grvity,                                    &
       e_co,kemit,snowh,numgas,                          &
       num_chem,num_moist,                                        &
       ids,ide, jds,jde, kds,kde,                                 &
       ims,ime, jms,jme, kms,kme,                                 &
       its,ite, jts,jte, kts,kte)
     print *,'gocart_run: done dry_dep_driver ...'

     ! -- ls wet deposition
     print *,'gocart_run: calling wetdep_ls ...'
     call wetdep_ls(dt,chem,rnav,moist,rho_phy,var_rmv,num_moist, &
         num_chem,numgas,p_qc,dz8w,vvel,chem_opt,        &
         ids,ide, jds,jde, kds,kde,                               &
         ims,ime, jms,jme, kms,kme,                               &
         its,ite, jts,jte, kts,kte)
     print *,'gocart_run: calling wetdep_ls ...'

    if (call_gocart) then
     print *,'gocart_run: calling GOCART CHEM driver ...'
      call gocart_chem_driver(ktau,dt,dtstep, gmt,julday,t_phy,moist, &
        chem,rho_phy,dz8w,p8w,backg_oh,oh_t,backg_h2o2,h2o2_t,backg_no3,no3_t, &
         dxy,grvity,xlat,xlong,ttday,tcosz, &
         chem_opt,num_chem,num_moist,                                      &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
     print *,'gocart_run: done GOCART CHEM driver'
     print *,'gocart_run: calling GOCART aerosols driver ...'
       call gocart_aerosols_driver(ktau,dtstep,t_phy,moist,  &
         chem,rho_phy,dz8w,p8w,dxy,grvity,         &
         chem_opt,num_chem,num_moist,                                      &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
     print *,'gocart_run: done GOCART aerosols driver'
    endif

    if (call_radiation) then
     print *,'gocart_run: calling radiation ...'
      store_arrays = .false.
      select case (aer_ra_feedback)
        case (1) 
     print *,'gocart_run: calling optical_driver ...'
          call optical_driver(curr_secs,dtstep,          &
               chem,dz8w,rri,relhum,                               &
               h2oai,h2oaj,                                        &
               tauaersw,gaersw,waersw,bscoefsw,tauaerlw,           &
               l2aer,l3aer,l4aer,l5aer,l6aer,l7aer,                &
               num_chem,chem_opt,ids,ide, jds,jde, kds,kde,        &
               ims,ime, jms,jme, kms,kme,                          &
               its,ite, jts,jte, kts,kte)
          call aer_opt_out(aod,dz8w,                           &
               ext_coeff,bscat_coeff,asym_par,                &
               tauaersw,gaersw,waersw,tauaerlw,               &
               num_ext_coef,num_bscat_coef,num_asym_par,      &
               ids,ide, jds,jde, kds,kde,                     &
               ims,ime, jms,jme, kms,kme,                     &
               its,ite, jts,jte, kts,kte)
          call aer_ra(dz8w                                     &
               ,extt,ssca,asympar,nbands                       &
               ,tauaersw,gaersw,waersw,tauaerlw                &
               ,ids,ide, jds,jde, kds,kde                      &
               ,ims,ime, jms,jme, kms,kme                      &
               ,its,ite, jts,jte, kts,kte)
          store_arrays = .true.
        case (2) 
          call aero_opt('sw',dz8w,chem         &
                   ,rri,relhum,aod             &
                   ,extt,ssca,asympar,num_chem &
                   ,ids,ide, jds,jde, kds,kde  &
                   ,ims,ime, jms,jme, kms,kme  &
                   ,its,ite, jts,jte, kts,kte)
          store_arrays = .true.
        case default
          ! -- no feedback
      end select
      
      ! -- store aerosol optical variables for feedback in radiation
      if (store_arrays) then
        do nv = 1, nbands
          do j = jts, jte
            do k = kts, kte
              do i = its, ite
                ext_cof(i,j,k,nv) = extt   (i,k,j,nv)
                sscal  (i,j,k,nv) = ssca   (i,k,j,nv)
                asymp  (i,j,k,nv) = asympar(i,k,j,nv)
              end do
            end do
          end do
        end do
        aod2d(its:ite,jts:jte) = aod(its:ite,jts:jte)
      end if
    endif

    print *,'gocart_run: calling sum_pm_gocart ...'
    call sum_pm_gocart (                              &
         rri, chem,pm2_5_dry, pm2_5_dry_ec, pm10,     &
         num_chem,chem_opt,                  &
         ids,ide, jds,jde, kds,kde,                   &
         ims,ime, jms,jme, kms,kme,                   &
         its,ite, jts,jte, kts,kte)
    print *,'gocart_run: done sum_pm_gocart'

    ! -- pm25 and pm10 for output , not for tracer options
    print *,'gocart_run: setting output arrays ...'
    do j = jts, jte
      do k = kts, kte
        do i = its, ite
          pm25  (i,j,k) = pm2_5_dry(i,k,j)
          p10   (i,j,k) = pm10     (i,k,j)
          ebu_oc(i,j,k) = ebu      (i,k,j,p_ebu_oc)
        end do
      end do
    end do
    print *,'gocart_run: done output arrays'

    print *,'gocart_run: call GOCART?',call_gocart
    if (call_gocart) then
      print *,'gocart_run: set bg fields ...'
      do j = jts, jte
        do k = kts, kte
          do i = its, ite
            oh_bg  (i,j,k) = max(0., oh_t  (i,k,j))
            h2o2_bg(i,j,k) = max(0., h2o2_t(i,k,j))
            no3_bg (i,j,k) = max(0., no3_t (i,k,j))
          end do
        end do
      end do
      print *,'gocart_run: set bg fields done'
    end if

    ! -- put chem stuff back into tracer array
    print *,'gocart_run: put chem stuff back into tracer array ...'
    print *,'gocart_run: ntra, ntrb, ntra+ntrb, nbegin, num_chem, nbegin+num_chem', &
                         ntra, ntrb, ntra+ntrb, nbegin, num_chem, nbegin+num_chem

    do nv = 1, num_chem
      nvv = nbegin + nv
      do j = jts, jte
        do k = kts, kte
          do i = its, ite
            tr3d(i,j,k,nvv) = max(epsilc,chem(i,k,j,nv))
            trdp(i,j,k,nvv) = tr3d(i,j,k,nvv)*(pr3d(i,j,k)-pr3d(i,j,k+1))
            ! -- export tr3d
            tr3d_out(i,j,k,nvv) = real(tr3d(i,j,k,nvv), kind=CHEM_KIND_R8)
          end do
        end do
        wet_dep(its:ite,j,nv) = var_rmv(its:ite,j,nv)
      end do
    end do

  end subroutine gocart_advance

end module gocart_mod
