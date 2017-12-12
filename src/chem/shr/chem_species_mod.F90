module chem_species_mod

  use mpp_mod, only : mpp_error, FATAL
  use chem_config_mod
  use chem_state_mod

  implicit none

  private

  public :: chem_species_setup

contains

  subroutine chem_species_setup(config)

    type(chem_config_type), intent(in) :: config

    ! -- begin

    if (config % aer_ra_feedback == 1) then
      P_extcof3    = 1
      P_extcof55   = 2
      P_extcof106  = 3
      P_extcof3_5  = 4
      P_extcof8_12 = 5
      P_bscof3     = 1 
      P_bscof55    = 2 
      P_bscof106   = 3
      P_asympar3   = 1
      P_asympar55  = 2
      P_asympar106 = 3
    endif

    select case (config % chem_opt)
      case (CHEM_OPT_GOCART)
        ! -- gocart simple
        if (num_chem    /= 19) call mpp_error(FATAL, 'num_chem is not equal to 19')
        if ((num_emis_ant < 6) .and. (config % biomass_burn_opt == 1)) call mpp_error(FATAL, 'num_emis_ant smaller than 6')
        if (num_emis_ant  < 4) call mpp_error(FATAL, 'num_emis_ant smaller than 4')
        p_qv=1
        p_qc=2
        p_qi=3
        p_so2=1
        numgas=4
        p_sulf=2
        p_dms=3
        p_msa=4
        p_p25=5
        p_bc1=6
        p_bc2=7
        p_oc1=8
        p_oc2=9
        p_dust_1=10
        p_dust_2=11
        p_dust_3=12
        p_dust_4=13
        p_dust_5=14
        p_seas_1=15
        p_seas_2=16
        p_seas_3=17
        p_seas_4=18
        p_p10   =19
        p_e_bc  =1
        p_e_oc  =2
        p_e_sulf=3
        p_e_pm_25=4
        p_e_so2=5
        p_e_pm_10=6
        p_e_dms=7

        p_ebu_bc  =1
        p_ebu_oc  =2
        p_ebu_sulf=3
        p_ebu_pm25=4
        p_ebu_so2=5
        p_ebu_pm10=6
        p_ebu_dms=7

        p_ebu_in_bc  =1
        p_ebu_in_oc  =2
        p_ebu_in_sulf=3
        p_ebu_in_pm25=4
        p_ebu_in_so2=5
        p_ebu_in_pm10=6
        p_ebu_in_dms=7

        ! -- diagnostic dust and seasale stuff
        p_edust1=1
        p_edust2=2
        p_edust3=3
        p_edust4=4
        p_edust5=5
        p_eseas1=1
        p_eseas2=2
        p_eseas3=3
        p_eseas4=4

      case (CHEM_OPT_GOCART_RACM)
        ! -- gocart + racm
        if (num_chem    /= 66) call mpp_error(FATAL, 'num_chem is not equal to 66')
        if (num_emis_ant < 25) call mpp_error(FATAL, 'num_emis_ant smaller than 25')
        conv_tr_aqchem = 1
        ! -- initialize pointers for gas phase and aerosol stuff
        call gocartracm_pointers
        ! -- initializing photolysis pointers so far done in chem_alloc (no pointers yet)
        ! -- hydrometeors
        p_qv=1
        p_qc=2
        p_qi=3

      case (CHEM_OPT_RACM_SOA_VBS)
        ! -- racm + soa
        if (num_chem    /= 103) call mpp_error(FATAL, 'num_chem is not equal to 103')
        if (num_emis_ant <  25) call mpp_error(FATAL, 'num_emis_ant smaller than 25')
        ! -- initialize pointers for gas phase and aerosol stuff
        conv_tr_aqchem = 1
        call racmsoavbs_pointers
        ! -- initializing photolysis pointers so far done in chem_alloc (no pointers yet)
        ! -- hydrometeors
        p_qv=1
        p_qc=2
        p_qi=3
#if 0
      case (304)
        ! -- gocart fim light
        if (num_chem /= 13) call mpp_error(FATAL, ' num_chem is not equal 13 for gocart fimlight ')
        if ((num_emis_ant < 6) .and. (config % biomass_burn_opt == 1)) call mpp_error(FATAL, ' num_emis_ant smaller than 6 ')
        if (num_emis_ant < 4) call mpp_error(FATAL, ' num_emis_ant smaller than 4 ')

        ! -- set species
        p_qv=1
        p_qc=2
        p_qi=3
        p_so2=1
        numgas=4
        p_sulf=2
        p_dms=3
        p_msa=4
        p_p25=5
        p_bc1=6
        p_bc2=7
        p_oc1=8
        p_oc2=9
        p_dust_1=10
        p_dust_2=11
        p_seas_1=12
        p_seas_2=13
        p_e_bc  =1
        p_e_oc  =2
        p_e_sulf=3
        p_e_pm_25=4
        p_e_so2=5
        p_e_pm_10=6

        ! -- diagnostic dust and seasale stuff
        p_edust1=1
        p_edust2=2
        p_edust3=3
        p_edust4=4
        p_edust5=5
        p_eseas1=1
        p_eseas2=2
        p_eseas3=3
        p_eseas4=4
      case (500)
        ! -- 2 tracers

        if (num_chem    /= 5) call mpp_error(FATAL, 'num_chem is not equal to 5')
        if (num_emis_ant < 5) call mpp_error(FATAL, 'num_emis_ant smaller than 5')
        p_qv=1
        p_qc=2
        p_qi=3

        ! -- used for tracer only transport (originally implemented for "earth analyzer", June 2015)
        p_tr1=1
        p_tr2=2
        p_tr3=3
        p_tr4=4
        p_tr5=5
        p_e_tr1=1
        p_e_tr2=2
        p_e_tr3=3
        p_e_tr4=4
        p_e_tr5=5
      case (16)
        ! -- volcanic ash
        if (num_chem    /= 10) call mpp_error(FATAL, 'num_chem is not equal to 10 for Volcano run')
        if (num_emis_vol < 10) call mpp_error(FATAL, 'num_emis_vol smaller than 10')
        p_qv=1
        p_qc=2
        p_qi=3
        emiss_opt=7
        p_vash_1 = 1
        p_vash_2 = 2
        p_vash_3 = 3
        p_vash_4 = 4
        p_vash_5 = 5
        p_vash_6 = 6
        p_vash_7 = 7
        p_vash_8 = 8
        p_vash_9 = 9
        p_vash_10 = 10
        p_e_vash1 = 1
        p_e_vash2 = 2
        p_e_vash3 = 3
        p_e_vash4 = 4
        p_e_vash5 = 5
        p_e_vash6 = 6
        p_e_vash7 = 7
        p_e_vash8 = 8
        p_e_vash9 = 9
        p_e_vash10 = 10
        numgas=0
      case (502)
        ! -- volcanoc ash (4 bins) only
        if (num_chem    /= 4) call mpp_error(FATAL, 'num_chem is not equal to 4')
        if (num_emis_vol < 4) call mpp_error(FATAL, 'num_emis_vol smaller than 4')
        p_qv=1
        p_qc=2
        p_qi=3
        numgas=0
        p_vash_1 = 1
        p_vash_2 = 2
        p_vash_3 = 3
        p_vash_4 = 4
        p_e_vash1 = 1
        p_e_vash2 = 2
        p_e_vash3 = 3
        p_e_vash4 = 4
      case (317)
        ! -- gocart simple +volcanic ash simple
        if (num_chem    /= 17) call mpp_error(FATAL, 'num_chem is not equal to 17')
        if (num_emis_vol <  4) call mpp_error(FATAL, 'num_emis_vol smaller than 4')
        if ((num_emis_ant < 6) .and. (config % biomass_burn_opt == 1)) call mpp_error(FATAL, ' num_emis_ant smaller than 6 ')
        if (num_emis_ant < 4) call mpp_error(FATAL, ' num_emis_ant smaller than 4 ')
        p_qv=1
        p_qc=2
        p_qi=3
        p_so2=1
        numgas=4
        p_sulf=2
        p_dms=3
        p_msa=4
        p_p25=5
        p_bc1=6
        p_bc2=7
        p_oc1=8
        p_oc2=9
        p_dust_1=10
        p_dust_2=11
        p_seas_1=12
        p_seas_2=13
        p_e_bc  =1
        p_e_oc  =2
        p_e_sulf=3
        p_e_pm_25=4
        p_e_so2=5
        p_e_pm_10=6
        ! -- diagnostic dust and seasale stuff
        p_edust1=1
        p_edust2=2
        p_edust3=3
        p_edust4=4
        p_edust5=5
        p_eseas1=1
        p_eseas2=2
        p_eseas3=3
        p_eseas4=4
        p_vash_1 = 14
        p_vash_2 = 15
        p_vash_3 = 16
        p_vash_4 = 17
        p_e_vash1 = 1
        p_e_vash2 = 2
        p_e_vash3 = 3
        p_e_vash4 = 4
      case (316)
        ! -- gocart simple +volcanic ash
        if (num_chem    /= 23) call mpp_error(FATAL, 'num_chem is not equal to 23')
        if (num_emis_vol < 10) call mpp_error(FATAL, 'num_emis_vol smaller than 10')
        if ((num_emis_ant < 6) .and. (config % biomass_burn_opt == 1)) call mpp_error(FATAL, 'num_emis_ant smaller than 6')
        if (num_emis_ant  < 4) call mpp_error(FATAL, 'num_emis_ant smaller than 4')
        p_qv=1
        p_qc=2
        p_qi=3
        p_so2=1
        numgas=4
        p_sulf=2
        p_dms=3
        p_msa=4
        p_p25=5
        p_bc1=6
        p_bc2=7
        p_oc1=8
        p_oc2=9
        p_dust_1=10
        p_dust_2=11
        p_seas_1=12
        p_seas_2=13
        p_e_bc  =1
        p_e_oc  =2
        p_e_sulf=3
        p_e_pm_25=4
        p_e_so2=5
        p_e_pm_10=6
        ! -- diagnostic dust and seasale stuff
        p_edust1=1
        p_edust2=2
        p_edust3=3
        p_edust4=4
        p_edust5=5
        p_eseas1=1
        p_eseas2=2
        p_eseas3=3
        p_eseas4=4
        p_vash_1 = 14
        p_vash_2 = 15
        p_vash_3 = 16
        p_vash_4 = 17
        p_vash_5 = 18
        p_vash_6 = 19
        p_vash_7 = 20
        p_vash_8 = 21
        p_vash_9 = 22
        p_vash_10 =23
        p_e_vash1 = 1
        p_e_vash2 = 2
        p_e_vash3 = 3
        p_e_vash4 = 4
        p_e_vash5 = 5
        p_e_vash6 = 6
        p_e_vash7 = 7
        p_e_vash8 = 8
        p_e_vash9 = 9
        p_e_vash10 = 10
#endif
      case default
        call mpp_error(FATAL, 'chem_opt not implemented')
    end select

  end subroutine chem_species_setup

  subroutine tracer_pointers_500

    p_tr1=1
    p_tr2=2
    p_tr3=3
    p_tr4=4
    p_tr5=5
    p_e_tr1=1
    p_e_tr2=2
    p_e_tr3=3
    p_e_tr4=4
    p_e_tr5=5

  end subroutine tracer_pointers_500
!
!package   gocartracm_kpp     chem_opt==301                   -
!(chem:so2,sulf,no2,no,o3,hno3,h2o2,ald,hcho,op1,op2,paa,ora1,ora2,nh3,n2o5,no3,pan,hc3,hc5,hc8,eth,co,ete,olt,oli,tol,xyl,aco3,tpan,hono,hno4,ket,gly,mgly,dcb,onit,csl,iso,co2,ch4,udd,hket,api,lim,dien,macr,ho,ho2,dms,msa,p25,bc1,bc2,oc1,oc2,dust_1,dust_2,dust_3,dust_4,dust_5,seas_1,seas_2,seas_3,seas_4,p10
  subroutine gocartracm_pointers
  ! -- this module will initialize variable pointers for chem_opt "301"
  ! -- June 2015

  ! -- 51 variabkes for RACM gas phase
    p_so2  = 1
    p_sulf  = 2
    p_no2  = 3
    p_no  = 4
    p_o3  = 5
    p_hno3  = 6
    p_h2o2  = 7
    p_ald  = 8
    p_hcho  = 9
    p_op1  = 10
    p_op2  = 11
    p_paa  = 12
    p_ora1  = 13
    p_ora2  = 14
    p_nh3  = 15
    p_n2o5  = 16
    p_no3  = 17
    p_pan  = 18
    p_hc3  = 19
    p_hc5  = 20
    p_hc8  = 21
    p_eth  = 22
    p_co  = 23
    p_ete  = 24
    p_olt  = 25
    p_oli  = 26
    p_tol  = 27
    p_xyl  = 28
    p_aco3  = 29
    p_tpan  = 30
    p_hono  = 31
    p_hno4  = 32
    p_ket  = 33
    p_gly  = 34
    p_mgly  = 35
    p_dcb  = 36
    p_onit  = 37
    p_csl  = 38
    p_iso  = 39
    p_co2  = 40
    p_ch4  = 41
    p_udd  = 42
    p_hket  = 43
    p_api  = 44
    p_lim  = 45
    p_dien  = 46
    p_macr  = 47
    p_ho  = 48
    p_ho2  = 49
    p_dms  = 50
    p_msa  = 51
!
! 15 more for GOCART
!
    p_p25  = 52
    p_bc1  = 53
    p_bc2  = 54
    p_oc1  = 55
    p_oc2  = 56
    p_dust_1  = 57
    p_dust_2  = 58
    p_dust_3  = 59
    p_dust_4  = 60
    p_dust_5  = 61
    p_seas_1  = 62
    p_seas_2  = 63
    p_seas_3  = 64
    p_seas_4  = 65
    p_p10  = 66
! emissions
!package   ecptec          emiss_opt==5                   -
!emis_ant:e_iso,e_so2,e_no,e_no2,e_co,e_eth,e_hc3,e_hc5,e_hc8,e_xyl,e_ol2,e_olt,e_oli,e_tol,e_csl,e_hcho,e_ald,e_ket,e_ora2,e_nh3,e_pm_25,e_pm_10,e_oc,e_sulf,e_bc

    p_e_iso = 1
    p_e_so2 = 2
    p_e_no = 3
    p_e_no2 = 4
    p_e_co = 5
    p_e_eth = 6
    p_e_hc3 = 7
    p_e_hc5 = 8
    p_e_hc8 = 9
    p_e_xyl = 10
    p_e_olt = 11
    p_e_oli = 12
    p_e_tol = 13
    p_e_csl = 14
    p_e_hcho = 15
    p_e_ald = 16
    p_e_ket = 17
    p_e_ora2 = 18
    p_e_nh3 = 19
    p_e_pm_25 = 20
    p_e_pm_10 = 21
    p_e_oc = 22
    p_e_sulf = 23
    p_e_bc = 24
    p_e_dms = 25
! biomass burning
    p_ebu_iso = 1
    p_ebu_so2 = 2
    p_ebu_no = 3
    p_ebu_no2=4
    p_ebu_co =5
    p_ebu_eth=6
    p_ebu_hc3=7
    p_ebu_hc5=8
    p_ebu_hc8=9
    p_ebu_xyl=10
    p_ebu_olt=11
    p_ebu_oli=12
    p_ebu_tol=13
    p_ebu_csl=14
    p_ebu_hcho=15
    p_ebu_ald=16
    p_ebu_ket=17
    p_ebu_ora2=18
    p_ebu_nh3=19
    p_ebu_pm25=20
    p_ebu_pm10=21
    p_ebu_oc=22
    p_ebu_sulf=23
    p_ebu_bc=24
    p_ebu_dms=25


    p_ebu_in_iso = 1
    p_ebu_in_so2 = 2
    p_ebu_in_no = 3
    p_ebu_in_no2=4
    p_ebu_in_co =5
    p_ebu_in_eth=6
    p_ebu_in_hc3=7
    p_ebu_in_hc5=8
    p_ebu_in_hc8=9
    p_ebu_in_xyl=10
    p_ebu_in_olt=11
    p_ebu_in_oli=12
    p_ebu_in_tol=13
    p_ebu_in_csl=14
    p_ebu_in_hcho=15
    p_ebu_in_ald=16
    p_ebu_in_ket=17
    p_ebu_in_ora2=18
    p_ebu_in_nh3=19
    p_ebu_in_pm25=20
    p_ebu_in_pm10=21
    p_ebu_in_oc=22
    p_ebu_in_sulf=23
    p_ebu_in_bc=24
    p_ebu_in_dms=25

!dust and seasalt
! diagnostic dust and seasale stuff
    p_edust1=1
    p_edust2=2
    p_edust3=3
    p_edust4=4
    p_edust5=5
    p_eseas1=1
    p_eseas2=2
    p_eseas3=3
    p_eseas4=4

  end subroutine gocartracm_pointers
!
!package   gocartracm_kpp     chem_opt==108                   -
!(chem:so2,sulf,no2,no,o3,hno3,h2o2,ald,hcho,op1,op2,paa,ora1,ora2,nh3,n2o5,no3,pan,hc3,hc5,hc8,eth,co,ete,olt,oli,tol,xyl,aco3,tpan,hono,hno4,ket,gly,mgly,dcb,onit,csl,iso,co2,ch4,udd,hket,api,lim,dien,macr,ho,ho2,dms,msa,p25,bc1,bc2,oc1,oc2,dust_1,dust_2,dust_3,dust_4,dust_5,seas_1,seas_2,seas_3,seas_4,p10
  subroutine racmsoavbs_pointers

  ! -- this subroutine will initialize variable pointers for chem_opt "108"
  ! -- June 2015

  ! -- 51 variabkes for RACM gas phase
    p_so2  = 1
    p_sulf  = 2
    p_no2  = 3
    p_no  = 4
    p_o3  = 5
    p_hno3  = 6
    p_h2o2  = 7
    p_ald  = 8
    p_hcho  = 9
    p_op1  = 10
    p_op2  = 11
    p_paa  = 12
    p_ora1  = 13
    p_ora2  = 14
    p_nh3  = 15
    p_n2o5  = 16
    p_no3  = 17
    p_pan  = 18
    p_hc3  = 19
    p_hc5  = 20
    p_hc8  = 21
    p_eth  = 22
    p_co  = 23
    p_ete  = 24
    p_olt  = 25
    p_oli  = 26
    p_tol  = 27
    p_xyl  = 28
    p_aco3  = 29
    p_tpan  = 30
    p_hono  = 31
    p_hno4  = 32
    p_ket  = 33
    p_gly  = 34
    p_mgly  = 35
    p_dcb  = 36
    p_onit  = 37
    p_csl  = 38
    p_iso  = 39
    p_co2  = 40
    p_ch4  = 41
    p_udd  = 42
    p_hket  = 43
    p_api  = 44
    p_lim  = 45
    p_dien = 46
    p_macr = 47
    p_hace = 48
    p_ishp = 49
    p_ison = 50
    p_mahp = 51
    p_mpan = 52
    p_nald = 53
    p_sesq = 54
    p_mbo  = 55
    p_cvasoa1 = 56
    p_cvasoa2 = 57
    p_cvasoa3 = 58
    p_cvasoa4 = 59
    p_cvbsoa1 = 60
    p_cvbsoa2 = 61
    p_cvbsoa3 = 62
    p_cvbsoa4 = 63
    p_ho  = 64
    p_ho2 = 65
    p_so4aj = 66
    p_so4ai = 67
    p_nh4aj = 68
    p_nh4ai = 69
    p_no3aj = 70
    p_no3ai = 71
    p_naaj  = 72
    p_naai  = 73
    p_claj  = 74
    p_clai  = 75
    p_asoa1j = 76
    p_asoa1i = 77
    p_asoa2j = 78
    p_asoa2i = 79
    p_asoa3j = 80
    p_asoa3i = 81
    p_asoa4j = 82
    p_asoa4i = 83
    p_bsoa1j = 84
    p_bsoa1i = 85
    p_bsoa2j = 86
    p_bsoa2i = 87
    p_bsoa3j = 88
    p_bsoa3i = 89
    p_bsoa4j = 90
    p_bsoa4i = 91
    p_orgpaj = 92
    p_orgpai = 93
    p_ecj    = 94
    p_eci    = 95
    p_p25j   = 96
    p_p25i   = 97
    p_antha  = 98
    p_seas   = 99
    p_soila  = 100
    p_nu0    = 101
    p_ac0    = 102
    p_corn   = 103
!
! emissions
!package   ecptec          emiss_opt==5                   -
!emis_ant:e_iso,e_so2,e_no,e_no2,e_co,e_eth,e_hc3,e_hc5,e_hc8,e_xyl,e_ol2,e_olt,e_oli,e_tol,e_csl,e_hcho,e_ald,e_ket,e_ora2,e_nh3,e_pm_25,e_pm_10,e_oc,e_sulf,e_bc
    p_e_iso = 1
    p_e_so2 = 2
    p_e_no = 3
    p_e_no2 = 4
    p_e_co = 5
    p_e_eth = 6
    p_e_hc3 = 7
    p_e_hc5 = 8
    p_e_hc8 = 9
    p_e_xyl = 10
    p_e_olt = 11
    p_e_oli = 12
    p_e_tol = 13
    p_e_csl = 14
    p_e_hcho = 15
    p_e_ald = 16
    p_e_ket = 17
    p_e_ora2 = 18
    p_e_nh3 = 19
    p_e_pm_25 = 20
    p_e_pm_10 = 21
    p_e_oc = 22
    p_e_sulf = 23
    p_e_bc = 24
    p_e_dms = 25

! biomass burning
    p_ebu_iso = 1
    p_ebu_so2 = 2
    p_ebu_no = 3
    p_ebu_no2=4
    p_ebu_co =5
    p_ebu_eth=6
    p_ebu_hc3=7
    p_ebu_hc5=8
    p_ebu_hc8=9
    p_ebu_xyl=10
    p_ebu_olt=11
    p_ebu_oli=12
    p_ebu_tol=13
    p_ebu_csl=14
    p_ebu_hcho=15
    p_ebu_ald=16
    p_ebu_ket=17
    p_ebu_ora2=18
    p_ebu_nh3=19
    p_ebu_pm25=20
    p_ebu_pm10=21
    p_ebu_oc=22
    p_ebu_sulf=23
    p_ebu_bc=24
    p_ebu_dms=25

    p_ebu_in_iso = 1
    p_ebu_in_so2 = 2
    p_ebu_in_no = 3
    p_ebu_in_no2=4
    p_ebu_in_co =5
    p_ebu_in_eth=6
    p_ebu_in_hc3=7
    p_ebu_in_hc5=8
    p_ebu_in_hc8=9
    p_ebu_in_xyl=10
    p_ebu_in_olt=11
    p_ebu_in_oli=12
    p_ebu_in_tol=13
    p_ebu_in_csl=14
    p_ebu_in_hcho=15
    p_ebu_in_ald=16
    p_ebu_in_ket=17
    p_ebu_in_ora2=18
    p_ebu_in_nh3=19
    p_ebu_in_pm25=20
    p_ebu_in_pm10=21
    p_ebu_in_oc=22
    p_ebu_in_sulf=23
    p_ebu_in_bc=24
    p_ebu_in_dms=25
!dust and seasalt
! diagnostic dust and seasale stuff
    p_edust1=1
    p_edust2=2
    p_edust3=3
    p_edust4=4
    p_edust5=5
    p_eseas1=1
    p_eseas2=2
    p_eseas3=3
    p_eseas4=4

  end subroutine racmsoavbs_pointers

end module chem_species_mod
