module gocart_prep_mod

  use chem_rc_mod
  use chem_types_mod
  use chem_tracers_mod
  use chem_config_mod, only : CHEM_OPT_GOCART,       &
                              CHEM_OPT_GOCART_RACM,  &
                              CHEM_OPT_RACM_SOA_VBS, &
                              CHEM_OPT_MAX
  use chem_const_mod,  only : airmw, epsilc, rgasuniv

  implicit none

  private

  public :: gocart_prep

contains

  subroutine gocart_prep(readrestart,chem_opt,chem_in_opt,ktau,dtstep,tr3d,tk3d,sm3d,&
                       ts2d,us2d,rsds,pr3d,prl3d,ph3d,phl3d,emiss_ash_mass,emiss_ash_height,          &
                       emiss_ash_dt,dm0,emiss_tr_mass,emiss_tr_height,               &
                       emiss_tr_dt,snwdph2d,VFRAC2d,VTYPE2d,STYPE2d,us3d,vs3d,ws3d,           &
                       slmsk2d,zorl2d,exch,pb2d,hf2d,clayfrac,clayf,sandfrac,sandf,th_pvsrf,&
                       oh_backgd,h2o2_backgd,no3_backgd,backg_oh,backg_h2o2,backg_no3,p_gocart,            &
                       nvl_gocart, ttday,tcosz,gmt,julday,area,ero1,                 &
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
                       firesize_agsv,firesize_aggr,                                  &
                       ids,ide, jds,jde, kds,kde,                                    &
                       ims,ime, jms,jme, kms,kme,                                    &
                       its,ite, jts,jte, kts,kte, rc)

    IMPLICIT NONE

    ! -- input variables

    LOGICAL,      INTENT(IN) :: readrestart
    INTEGER,      INTENT(IN) :: chem_opt,chem_in_opt
    INTEGER,      INTENT(IN) :: ktau,nvl,nvlp1,ntra,ntrb,nvl_gocart
    INTEGER,      INTENT(IN) :: num_ebu_in,num_soil_layers,num_chem,num_moist,julday,     &
                                   num_emis_vol,num_emis_ant,ids,ide, jds,jde, kds,kde,      &
                                   kemit,ims,ime, jms,jme, kms,kme,                          &
                                   its,ite, jts,jte, kts,kte
    LOGICAL,      INTENT(IN) :: call_gocart
    REAL(CHEM_KIND_R4), INTENT(IN) :: g,rd,p1000,cp,dtstep,gmt
!---------------------------
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl, ntra+ntrb), intent(inout) :: tr3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl),            intent(inout) :: tk3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, num_soil_layers), intent(in) :: sm3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: ts2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: us2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: rsds
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvlp1), intent(in) :: pr3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl),   intent(in) :: prl3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvlp1), intent(in) :: ph3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl),   intent(in) :: phl3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(inout) :: emiss_ash_mass
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(inout) :: emiss_ash_height
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: emiss_ash_dt
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: dm0
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: emiss_tr_mass
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: emiss_tr_height
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: emiss_tr_dt
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: snwdph2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: vfrac2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: vtype2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: stype2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl), intent(in) :: us3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl), intent(in) :: vs3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl), intent(in) :: ws3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: slmsk2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: zorl2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl), intent(in) :: exch
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: pb2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: hf2d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: clayfrac
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: sandfrac
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: th_pvsrf
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl_gocart), intent(in) :: oh_backgd
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl_gocart), intent(in) :: h2o2_backgd
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl_gocart), intent(in) :: no3_backgd
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: backg_oh
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: backg_h2o2
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: backg_no3
    real(CHEM_KIND_R4), dimension(nvl_gocart+1),     intent(in) :: p_gocart
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: ttday
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: tcosz
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: area
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: ero1
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: ero2
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: ero3
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: rcav
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: raincv_b
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: deg_lat
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(in) :: deg_lon
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: relhum
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: rri
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: t_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, num_moist), intent(out) :: moist
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: u_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: v_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: p_phy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme, num_chem), intent(out) :: chem
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: tsk
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 3), intent(inout) :: erod
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kemit, jms:jme, num_emis_ant), intent(inout) :: emis_ant
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme,   jms:jme, num_emis_vol), intent(inout) :: emis_vol
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: e_co
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: dms_0
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: u10
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: v10
    integer,            dimension(ims:ime, jms:jme), intent(out) :: ivgtyp
    integer,            dimension(ims:ime, jms:jme), intent(out) :: isltyp
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: gsw
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: vegfra
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: rmol
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: ust
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: znt 
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: xland
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: dxy
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: t8w
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: p8w
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: exch_h
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: pbl
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: hfx
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: snowh
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: clayf
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: sandf
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: xlat
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme), intent(out) :: xlong
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: convfac
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: z_at_w
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: zmid
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: dz8w
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: vvel
    real(CHEM_KIND_R4), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: rho_phy
    real(CHEM_KIND_R4), dimension(ims:ime, num_soil_layers, jms:jme), intent(out) :: smois
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, num_emis_ant),   intent(in) :: emiss_ab
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, num_ebu_in),     intent(in) :: emiss_abu
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, num_ebu_in),     intent(out) :: ebu_in
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 8), intent(in) :: plumestuff
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme),    intent(out) :: mean_fct_agtf
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme),    intent(out) :: mean_fct_agef
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme),    intent(out) :: mean_fct_agsv
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme),    intent(out) :: mean_fct_aggr
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme),    intent(out) :: firesize_agtf
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme),    intent(out) :: firesize_agef
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme),    intent(out) :: firesize_agsv
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme),    intent(out) :: firesize_aggr


    INTEGER, OPTIONAL, INTENT(OUT) :: rc

#if 0
    REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
           INTENT(OUT   ) ::                                                 &
                                                         rri,               &
                                                       t_phy,               &
                                                       p_phy,               &
                                               relhum, dz8w,p8w,t8w,        &
                                               z_at_w , zmid ,exch_h,       &
                                               u_phy,v_phy,vvel,rho_phy,    &
                                               convfac


    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl), intent(in) :: exch
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvlp1), intent(in) :: ph3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl), intent(in) :: phl3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvlp1), intent(in) :: pr3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl), intent(in) :: prl3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, 4),   intent(in) :: sm3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl), intent(inout) :: tk3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl, ntra+ntrb), intent(inout) :: tr3d

!---------------------------
    real,         intent(in) :: oh_backgd(nvl_gocart,jms:jme)      !
    real,         intent(in) :: h2o2_backgd(nvl_gocart,jms:jme)      !
    real,         intent(in) :: no3_backgd(nvl_gocart,jms:jme)      !

    real,         intent(in) :: ts2d(jms:jme)          ! skin temperature
    real,         intent(in) :: us2d(jms:jme)          ! friction velocity/equivalent momentum flux
    real,         intent(in) :: pb2d(jms:jme)          ! 
    real,         intent(in) :: th_pvsrf(jms:jme)          ! 
    real,         intent(in) :: rcav(jms:jme)          ! 
    real,         intent(in) :: hf2d(jms:jme)          ! 
    real,         intent(in) :: snwdph2d(jms:jme)      !
    real,         intent(in) :: rsds(jms:jme)          ! downward short-wave radiation flux

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
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl), intent(in) :: us3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl), intent(in) :: vs3d
    real(CHEM_KIND_R4), dimension(ims:ime, jms:jme, nvl), intent(in) :: ws3d
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

#endif
    ! -- local variables
    integer i,j,k,kk,nv,jmax,jmaxi,l,ll,n,ndystep,ixhour
    real(CHEM_KIND_R4) ::  maxv,factor,factor2,pu,pl,aln,pwant,rlat
    real(CHEM_KIND_R4) ::  thv,xhour,xmin,gmtp,xlonn,xtime,real_time
    real(CHEM_KIND_R4), DIMENSION (1,1) :: sza,cosszax
    real(CHEM_KIND_R4), DIMENSION (ims:ime,jms:jme) :: so2_mass

    ! -- volcanic stuff
    integer :: ko,k_final,k_initial,kl,kk4,curr_hours,curr_secs
    real(CHEM_KIND_R4) :: x1,ashz_above_vent
    real(CHEM_KIND_R4), DIMENSION (kms:kme) :: vert_mass_dist
    real(CHEM_KIND_R4) :: eh,h1,h2,h3,h4,h5,h6,maxth
    logical, save :: first_init = .true.

    ! -- volcano ashes parameters
    !  + original
    ! real, dimension(6) :: h = (/ 9., 16., 58., 79., 109., 129., huge(1.0) /)
    !  + real-time default (if volcano starts at h = 0)
    real(CHEM_KIND_R4), dimension(7) :: h = (/ (240., i = 1, 6), huge(1.0) /)
    real(CHEM_KIND_R4), dimension(6) :: emiss_ash_table = (/  5834.,  3834.,  5834.,  3334.,  3334.,  2334. /)
    real(CHEM_KIND_R4), dimension(6) :: eh_table        = (/ 3.11e5, 3.87e4, 3.11e5, 2.17e4, 2.17e4, 4.93e3 /)
    real(CHEM_KIND_R4), parameter :: percen_mass_umbrel = 0.75
    real(CHEM_KIND_R4), parameter :: base_umbrel        = 0.25    ! fraction
    real(CHEM_KIND_R4), parameter :: base_umbrel2       = 1.0     ! evenly distribution
    ! .. Intrinsic Functions ..
    INTRINSIC max, min, float

    ! -- begin
    print *,'chem_prep: begin...'
    if (present(rc)) rc = CHEM_RC_SUCCESS

    real_time=float(ktau)*dtstep/60.

    so2_mass = 0.

    if (ktau <= 1) then
      emis_ant = 0.
      emis_vol = 0.
    end if

    e_co = 0.

    do j=jts,jte
      do i=its,ite
         z_at_w(i,kts,j)=max(0.,ph3d(i,j,kts)/g)
      enddo
    enddo

    do j=jts,jte
      do k=kts,kte
        do i=its,ite
          dz8w(i,k,j)=(ph3d(i,j,k+1)-ph3d(i,j,k))/g
          if (dz8w(i,k,j) < 0.) dz8w(i,k,j)=-dz8w(i,k,j)
          z_at_w(i,k+1,j)=z_at_w(i,k,j)+dz8w(i,k,j)
        enddo
      enddo
    enddo

    do j=jts,jte
      do k=kts,kte+1
        do i=its,ite
          p8w(i,k,j)=pr3d(i,j,k)
        enddo
      enddo
    enddo

    do j=jts,jte
      do i=its,ite
        raincv_b(i,j)=rcav(i,j)
        pbl(i,j)=pb2d(i,j)
        dms_0(i,j)=dm0(i,j)
        hfx(i,j)=hf2d(i,j)
        snowh(i,j)=snwdph2d(i,j)*.001
        erod(i,j,1)=ero1(i,j)
        erod(i,j,2)=ero2(i,j)
        erod(i,j,3)=ero3(i,j)
        xlat(i,j)=deg_lat(i,j)
        xlong(i,j)=deg_lon(i,j)
        ust(i,j)=us2d(i,j)
        tsk(i,j)=ts2d(i,j)
!       ivgtyp(i,j)=VTYPE2d(i,j)
!       isltyp(i,j)=STYPE2d(i,j)
        gsw(i,j)=rsds(i,j)
        vegfra(i,j)=VFRAC2d(i,j)
        rmol(i,j)=0.
        znt(i,j)=zorl2d(i,j)*.01
!SLMSK   - SEA(0),LAND(1),ICE(2) MASK
!       xland(i,j)=1.
!       if (slmsk2d(i,j) == 0.) then
!         xland(i,j) = 0.
!       else if (slmsk2d(i,j) == 1.) then
!         xland(i,j) = 1.
!       else if (slmsk2d(i,j) == 2.) then
!         xland(i,j) = 2.
!       end if
        xland(i,j) = real(nint(slmsk2d(i,j)), CHEM_KIND_R4)
        if (nint(slmsk2d(i,j)) == 2) then
          isltyp(i,j) = 16
        else
          isltyp(i,j) = int(stype2d(i,j)+0.5_CHEM_KIND_R4)
        end if
        ivgtyp(i,j) = int(vtype2d(i,j)+0.5_CHEM_KIND_R4)
        dxy(i,j)=area(i,j)
        u10(i,j)=us3d(i,j,1)
        v10(i,j)=vs3d(i,j,1)
        clayf(i,j) = clayfrac(i,j)
        sandf(i,j) = sandfrac(i,j)
      enddo
    enddo


    factor=0.
    jmax=0
    jmaxi=0
    k=1
    if ((p_bc2 > 1) .or. (chem_opt == CHEM_OPT_RACM_SOA_VBS)) then  ! "regular" chem options
      do j=jts,jte
        do i=its,ite
          k=1
          emis_ant(i,k,j,p_e_bc)=emiss_ab(i,j,p_e_bc)
          emis_ant(i,k,j,p_e_oc)=emiss_ab(i,j,p_e_oc)
          emis_ant(i,k,j,p_e_sulf)=emiss_ab(i,j,p_e_sulf)
          emis_ant(i,k,j,p_e_so2)=emiss_ab(i,j,p_e_so2)
          emis_ant(i,k,j,p_e_dms)= 0. !emiss_ab(j,p_e_dms)
          emis_ant(i,k,j,p_e_pm_25)=emiss_ab(i,j,p_e_pm_25)
          emis_ant(i,k,j,p_e_pm_10)=emiss_ab(i,j,p_e_pm_10)

          ! -- gas anth emission for chem_opt 301
          if ((chem_opt == CHEM_OPT_GOCART_RACM) .or. (chem_opt == CHEM_OPT_RACM_SOA_VBS)) then
            emis_ant(i,k,j,p_e_iso)=emiss_ab(i,j,p_e_iso)
            emis_ant(i,k,j,p_e_no)=emiss_ab(i,j,p_e_no)
            emis_ant(i,k,j,p_e_no2)=emiss_ab(i,j,p_e_no2)
            emis_ant(i,k,j,p_e_co)=emiss_ab(i,j,p_e_co)
            emis_ant(i,k,j,p_e_eth)=emiss_ab(i,j,p_e_eth)
            emis_ant(i,k,j,p_e_hc3)=emiss_ab(i,j,p_e_hc3)
            emis_ant(i,k,j,p_e_hc5)=emiss_ab(i,j,p_e_hc5)
            emis_ant(i,k,j,p_e_hc8)=emiss_ab(i,j,p_e_hc8)
            emis_ant(i,k,j,p_e_xyl)=emiss_ab(i,j,p_e_xyl)
            emis_ant(i,k,j,p_e_olt)=emiss_ab(i,j,p_e_olt)
            emis_ant(i,k,j,p_e_oli)=emiss_ab(i,j,p_e_oli)
            emis_ant(i,k,j,p_e_tol)=emiss_ab(i,j,p_e_tol)
            emis_ant(i,k,j,p_e_csl)=emiss_ab(i,j,p_e_csl)
            emis_ant(i,k,j,p_e_hcho)=emiss_ab(i,j,p_e_hcho)
            emis_ant(i,k,j,p_e_ald)=emiss_ab(i,j,p_e_ald)
            emis_ant(i,k,j,p_e_ket)=emiss_ab(i,j,p_e_ket)
            emis_ant(i,k,j,p_e_ora2)=emiss_ab(i,j,p_e_ora2)
            emis_ant(i,k,j,p_e_nh3)=emiss_ab(i,j,p_e_nh3)
          endif
         
          ! -- gas biomass burning emission for chem_opt 301
          if ((chem_opt == CHEM_OPT_GOCART_RACM) .or. (chem_opt == CHEM_OPT_RACM_SOA_VBS)) then
            ebu_in(i,j,p_ebu_in_iso)=emiss_abu(i,j,p_e_iso)
            ebu_in(i,j,p_ebu_in_no)=emiss_abu(i,j,p_e_no)
            ebu_in(i,j,p_ebu_in_no2)=emiss_abu(i,j,p_e_no2)
            ebu_in(i,j,p_ebu_in_co)=emiss_abu(i,j,p_e_co)
            ebu_in(i,j,p_ebu_in_eth)=emiss_abu(i,j,p_e_eth)
            ebu_in(i,j,p_ebu_in_hc3)=emiss_abu(i,j,p_e_hc3)
            ebu_in(i,j,p_ebu_in_hc5)=emiss_abu(i,j,p_e_hc5)
            ebu_in(i,j,p_ebu_in_hc8)=emiss_abu(i,j,p_e_hc8)
            ebu_in(i,j,p_ebu_in_xyl)=emiss_abu(i,j,p_e_xyl)
            ebu_in(i,j,p_ebu_in_olt)=emiss_abu(i,j,p_e_olt)
            ebu_in(i,j,p_ebu_in_oli)=emiss_abu(i,j,p_e_oli)
            ebu_in(i,j,p_ebu_in_tol)=emiss_abu(i,j,p_e_tol)
            ebu_in(i,j,p_ebu_in_csl)=emiss_abu(i,j,p_e_csl)
            ebu_in(i,j,p_ebu_in_hcho)=emiss_abu(i,j,p_e_hcho)
            ebu_in(i,j,p_ebu_in_ald)=emiss_abu(i,j,p_e_ald)
            ebu_in(i,j,p_ebu_in_ket)=emiss_abu(i,j,p_e_ket)
            ebu_in(i,j,p_ebu_in_ora2)=emiss_abu(i,j,p_e_ora2)
            ebu_in(i,j,p_ebu_in_nh3)=emiss_abu(i,j,p_e_nh3)
          endif
          ebu_in(i,j,p_ebu_in_oc)=emiss_abu(i,j,p_e_oc)
          ebu_in(i,j,p_ebu_in_bc)=emiss_abu(i,j,p_e_bc)
          ebu_in(i,j,p_ebu_in_pm25)=emiss_abu(i,j,p_e_pm_25)
          ebu_in(i,j,p_ebu_in_pm10)=emiss_abu(i,j,p_e_pm_10)
          ebu_in(i,j,p_ebu_in_so2)=emiss_abu(i,j,p_e_so2)
          ebu_in(i,j,p_ebu_in_dms)= 0. !emiss_abu(j,p_e_dms)
         !ebu_in(i,j,p_ebu_in_sulf)=0. ! for now

          mean_fct_agtf(i,j)=plumestuff(i,j,1)
          mean_fct_agef(i,j)=plumestuff(i,j,2)
          mean_fct_agsv(i,j)=plumestuff(i,j,3)
          mean_fct_aggr(i,j)=plumestuff(i,j,4)
          firesize_agtf(i,j)=plumestuff(i,j,5)
          firesize_agef(i,j)=plumestuff(i,j,6)
          firesize_agsv(i,j)=plumestuff(i,j,7)
          firesize_aggr(i,j)=plumestuff(i,j,8)
        enddo
      enddo

    else if (p_tr2 > 1) then  ! tracer options

      ! -- tracer run
      do j=jts,jte
        do i=its,ite
          k=kts
          emis_ant(i,k,j,p_e_tr1)=emiss_ab(i,j,p_e_tr1)
          emis_ant(i,k,j,p_e_tr2)=emiss_ab(i,j,p_e_tr2)
        enddo
      enddo

    else if ((p_tr2 > 1) .and. (p_bc2 > 1)) then

      call chem_rc_set(CHEM_RC_FAILURE, msg="Inconsistent options detected.", &
        file=__FILE__, line=__LINE__, rc=rc)
      return

    endif

#if 0
    do i=its,ite
      do j=jts,jte
        do k=kts,kte
          thv=tr3d(i,j,k,p_atm_ptem)/(1.+0.6078*tr3d(i,j,k,p_atm_cldq))
          tk3d(i,j,k)=thv*(.5*(p8w(i,k,j)+p8w(i,k+1,j))/p1000)**(rd/cp)
        enddo
      enddo
    enddo
#endif

    do j=jts,jte
      do k=kts,kte+1
        kk=min(k,kte)
        do i=its,ite
!         zmid(i,k,j)=.5*(ph3d(i,j,kk+1)+ph3d(i,j,kk))/g
          zmid(i,k,j)=phl3d(i,j,kk)/g
          dz8w(i,k,j)=z_at_w(i,kk+1,j)-z_at_w(i,kk,j)
          t_phy(i,k,j)=tk3d(i,j,kk)
!         relhum(i,k,j)=rh3d(kk,j)
!         p_phy(i,k,j)=.5*(p8w(i,kk,j)+p8w(i,kk+1,j))
          p_phy(i,k,j)=prl3d(i,j,kk)
          u_phy(i,k,j)=us3d(i,j,kk)
          exch_h(i,k,j)=exch(i,j,kk)
          v_phy(i,k,j)=vs3d(i,j,kk)
!       print *,'--> i,k,j, RD, p, T, t3d = ',i,k,j,RD,p_phy(i,k,j),T_phy(i,k,j),tr3d(kk,j,2)
          rho_phy(i,k,j)=p_phy(i,k,j)/(RD*T_phy(i,k,j)*(1.+.608*tr3d(i,j,kk,p_atm_shum)))
          rri(i,k,j)=1./rho_phy(i,k,j)
!       print *,'--> i,k,j, rho, 1/rho = ',i,k,j,rho_phy(i,k,j),rri(i,k,j)
          vvel(i,k,j)=-ws3d(i,j,kk)*rri(i,k,j)/g
          convfac(i,k,j)=p_phy(i,k,j)/rgasuniv/t_phy(i,k,j)
          moist(i,k,j,:)=0.
          moist(i,k,j,1)=tr3d(i,j,kk,p_atm_shum)
          if (t_phy(i,k,j) > 265.) then
            moist(i,k,j,2)=tr3d(i,j,kk,p_atm_cldq)
            moist(i,k,j,3)=0.
            if (moist(i,k,j,2) < 1.e-8) moist(i,k,j,2)=0.
          else
            moist(i,k,j,2)=0.
            moist(i,k,j,3)=tr3d(i,j,kk,p_atm_cldq)
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
            chem(i,k,j,nv)=tr3d(i,j,kk,ntra+nv)
          enddo
        enddo
      enddo
    enddo

    print *,'chem_prep: chem array done'
    print *,'chem_prep: readrestart  =',readrestart

    if (.NOT. readrestart) then
      if (chem_in_opt == 0 ) then
        if(ktau.le.1)then
!           if(chem_opt > 0 ) then
          do j=jts,jte
            do k=kts,kte
              do i=its,ite
                if (chem_opt == CHEM_OPT_GOCART) then
                  do n=1,num_chem
                    chem(i,k,j,n)=1.e-12
                  enddo
                endif  ! chem_opt==300
                chem(i,k,j,p_so2)=5.e-6
                chem(i,k,j,p_sulf)=3.e-6
                if ((chem_opt >= CHEM_OPT_GOCART) .and. (chem_opt < CHEM_OPT_MAX)) then
                  chem(i,k,j,p_msa)=0.1e-6
                  chem(i,k,j,p_dms)=0.1e-6
                  chem(i,k,j,p_bc1)=0.1e-3
                  chem(i,k,j,p_bc2)=0.1e-3
                  chem(i,k,j,p_oc1)=0.1e-3
                  chem(i,k,j,p_oc2)=0.1e-3
                  chem(i,k,j,p_p25)=1.
                  chem(i,k,j,p_p10)=1.
                endif !chem_opt >= 300 .and. chem_opt <  500

                if ((chem_opt == CHEM_OPT_GOCART_RACM) .or. (chem_opt == CHEM_OPT_RACM_SOA_VBS)) then  !added o3 background !lzhang
                  kk=min(k,kte)
                  ! -- add initial constant into O3,CH4 and CO ect.
                  chem(i,k,j,p_o3)=epsilc
                  maxth=min(400.,th_pvsrf(i,j))
                  if (tr3d(i,j,kk,p_atm_ptem) > maxth) then
                    chem(i,k,j,p_o3)=(airmw/48.)*tr3d(i,j,kk,p_atm_o3mr)*1e6 !convert kg/kg to ppm
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

        if ((ktau<=1).and.((chem_opt == CHEM_OPT_GOCART_RACM).or.(chem_opt == CHEM_OPT_RACM_SOA_VBS))) then  !added GFS o3 background above 380K!lzhang
          do j=jts,jte
            do k=kts,kte+1
              kk=min(k,kte)
              do i=its,ite
                maxth=min(400.,th_pvsrf(i,j))
                if (tr3d(i,j,kk,p_atm_ptem) >= maxth) then
                  chem(i,k,j,p_o3)=(airmw/48.)*tr3d(i,j,kk,p_atm_o3mr)*1e6 !convert kg/kg to ppm
                endif !380K
              enddo
            enddo
          enddo
        endif ! chem_opt == 301.or.chem_opt==108
        
      endif !(chem_in_opt == 1 )

      print *,'chem_prep: chem in done'
      
#if 0
      if (ktau<=1.and.chem_opt == CHEM_OPT_RACM_SOA_VBS) then
           call aerosols_soa_vbs_init(chem,convfac,z_at_w,           &
               pm2_5_dry,pm2_5_water,pm2_5_dry_ec,                   &
               chem_in_opt, aer_ic_opt,                              &
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
if (first_init .and. chem_opt == CHEM_OPT_RACM_SOA_VBS) then
           first_init = .false.
           call aerosols_soa_vbs_init(chem,convfac,z_at_w,           &
               pm2_5_dry,pm2_5_water,pm2_5_dry_ec,                   &
               chem_in_opt, aer_ic_opt,                              &
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
    if (call_gocart .and. (chem_opt == CHEM_OPT_GOCART))then
      print *,'chem_prep: call_gocart enter...'
      do j=jts,jte
        do i=its,ite
          do k=kts,kte
            do ll=2,nvl_gocart
              l=ll
              if (p_gocart(l) < .01*p_phy(i,k,j)) exit
            enddo
            pu=alog(p_gocart(l))
            pl=alog(p_gocart(l-1))
            pwant=alog(.01*p_phy(i,k,j))
            if (pwant > pl)then
              backg_oh(i,k,j)=oh_backgd(i,j,l)
              backg_h2o2(i,k,j)=h2o2_backgd(i,j,l)
              backg_no3(i,k,j)=no3_backgd(i,j,l)
            else
              aln=(oh_backgd(i,j,l)*(pwant-pl)+            &
                oh_backgd(i,j,l-1)*(pu-pwant))/(pu-pl)
              backg_oh(i,k,j)=aln
              aln=(h2o2_backgd(i,j,l)*(pwant-pl)+            &
                h2o2_backgd(i,j,l-1)*(pu-pwant))/(pu-pl)
              backg_h2o2(i,k,j)=aln
              aln=(no3_backgd(i,j,l)*(pwant-pl)+            &
                no3_backgd(i,j,l-1)*(pu-pwant))/(pu-pl)
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
      if (chem_opt == CHEM_OPT_GOCART) then
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
      call chem_rc_set(CHEM_RC_FAILURE, msg="Inconsistent options detected.", &
        file=__FILE__, line=__LINE__, rc=rc)
      return
    endif
    print *,'chem_prep: chem + emiss: done ...'

    print *,'chem_prep: smois: start'
    do j=jts,jte
      do nv=1,num_soil_layers
        do i=its,ite
          smois(i,nv,j)=sm3d(i,j,nv)
        enddo
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
      do i=its,ite
        if (emiss_ash_dt(i,j) > 0.) then
          so2_mass(i,j)=1.5e4*3600.*1.e9/64./area(i,j)
          eh=2600.*(emiss_ash_height(i,j)*.0005)**4.1494
          emiss_ash_mass(i,j)=eh*1.e9/area(i,j)
        end if
      enddo
    enddo

    ! -- hardcode for special retro case (set h1 - h6 properly
    do nv = 1, 6
      if ((curr_hours >= h(nv)) .and. (curr_hours < h(nv+1))) then
        do j = jts, jte
          do i = its, ite
            if (emiss_ash_dt(i,j) > 0) then
              emiss_ash_height(i,j) = emiss_ash_table(nv)
              emiss_ash_mass(i,j)   = 1.e+09 * eh_table(nv) / area(i,j)
            end if
          end do
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
        do i=its,ite
          if (emiss_ash_dt(i,j)     <= 0) cycle
          if (emiss_ash_height(i,j) <= 0) cycle
          ashz_above_vent=emiss_ash_height(i,j) +z_at_w(i,kts,j)
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
          if (chem_opt == 316 ) then 
            do ko=1,k_final
              emis_vol(i,ko,j,p_e_vash1)=.02*vert_mass_dist(ko)*emiss_ash_mass(i,j)
              emis_vol(i,ko,j,p_e_vash2)=.04*vert_mass_dist(ko)*emiss_ash_mass(i,j)
              emis_vol(i,ko,j,p_e_vash3)=.11*vert_mass_dist(ko)*emiss_ash_mass(i,j)
              emis_vol(i,ko,j,p_e_vash4)=.09*vert_mass_dist(ko)*emiss_ash_mass(i,j)
              emis_vol(i,ko,j,p_e_vash5)=.09*vert_mass_dist(ko)*emiss_ash_mass(i,j)
              emis_vol(i,ko,j,p_e_vash6)=.13*vert_mass_dist(ko)*emiss_ash_mass(i,j)
              emis_vol(i,ko,j,p_e_vash7)=.16*vert_mass_dist(ko)*emiss_ash_mass(i,j)
              emis_vol(i,ko,j,p_e_vash8)=.16*vert_mass_dist(ko)*emiss_ash_mass(i,j)
              emis_vol(i,ko,j,p_e_vash9)=.1*vert_mass_dist(ko)*emiss_ash_mass(i,j)
              emis_vol(i,ko,j,p_e_vash10)=.1*vert_mass_dist(ko)*emiss_ash_mass(i,j)
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
          elseif (chem_opt == 317 .or. chem_opt == 502) then

            ! -- reduced vocanic ash transport
            do ko=1,k_final
              emis_vol(i,ko,j,p_e_vash1)=.11*vert_mass_dist(ko)*emiss_ash_mass(i,j)
              emis_vol(i,ko,j,p_e_vash2)=.08*vert_mass_dist(ko)*emiss_ash_mass(i,j)
              emis_vol(i,ko,j,p_e_vash3)=.05*vert_mass_dist(ko)*emiss_ash_mass(i,j)
              emis_vol(i,ko,j,p_e_vash4)=.035*vert_mass_dist(ko)*emiss_ash_mass(i,j)
            enddo
          elseif (chem_opt == CHEM_OPT_GOCART) then

            ! -- if applied to gocart we only need finest ash bins, we use the coarse one for so2
            do ko=1,k_final
              emis_vol(i,ko,j,p_e_vash1)=vert_mass_dist(ko)*so2_mass(i,j)
              emis_vol(i,ko,j,p_e_vash2)=.08*vert_mass_dist(ko)*emiss_ash_mass(i,j)
              emis_vol(i,ko,j,p_e_vash3)=.05*vert_mass_dist(ko)*emiss_ash_mass(i,j)
              emis_vol(i,ko,j,p_e_vash4)=.035*vert_mass_dist(ko)*emiss_ash_mass(i,j)
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
    print *,'-- chem_opt = ',chem_opt
    if ((chem_opt == CHEM_OPT_RACM_SOA_VBS) .or. (chem_opt >= CHEM_OPT_GOCART .and. chem_opt < CHEM_OPT_MAX)) then
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
    if (chem_opt == 316) then
      do j = jts, jte
        do k = kts, kte-2
          do i = its, ite
            if (emiss_ash_dt(i,j) <= 0.) cycle
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

    if ((chem_opt == 317) .or. (chem_opt == 502))  then
      do j=jts,jte
        do k=kts,kte-2
          do i=its,ite
            if (emiss_ash_dt(i,j) <= 0.) cycle
            factor2=dtstep*rri(i,k,j)/dz8w(i,k,j)
            chem(i,k,j,p_vash_1)=chem(i,k,j,p_vash_1)+emis_vol(i,k,j,p_e_vash1)*factor2
            chem(i,k,j,p_vash_2)=chem(i,k,j,p_vash_2)+emis_vol(i,k,j,p_e_vash2)*factor2
            chem(i,k,j,p_vash_3)=chem(i,k,j,p_vash_3)+emis_vol(i,k,j,p_e_vash3)*factor2
            chem(i,k,j,p_vash_4)=chem(i,k,j,p_vash_4)+emis_vol(i,k,j,p_e_vash4)*factor2
          enddo
        enddo
      enddo
    endif

    if(chem_opt == CHEM_OPT_GOCART)  then
      ! -- for gocart only lump ash into p25 and p10
      do j=jts,jte
        do k=kts,kte-2
          do i=its,ite
            if (emiss_ash_dt(i,j) <= 0.) cycle
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
#if 0
    print *,'chem_prep: stage 5'

    ! -- option 501 was only used for cesium ensemble - Japan 2010
    if (chem_opt == 501) then
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
#endif
      print *,'chem_prep: stage end'

  end subroutine gocart_prep

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

end module gocart_prep_mod
