module gocart_mod

  use chem_const_mod
  use chem_shr_mod,    only : chem_config_type, chem_prep
  use chem_config_mod, config => chem_config
  use chem_state_mod
  use chem_vars_mod
  use chem_domain_mod

! use gocart_dust_mod
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
! use time_mod

  implicit none

  public

contains

  subroutine gocart_init

    ! -- local variables

    ! -- begin

#if 0
    write(6,'("GOCART input config:")')
    write(6,'("    config % chem_opt         = ",i0)') config % chem_opt
    write(6,'("    config % chem_in_opt      = ",i0)') config % chem_in_opt
    write(6,'("    config % dust_opt         = ",i0)') config % dust_opt
    write(6,'("    config % dmsemis_opt      = ",i0)') config % dmsemis_opt
    write(6,'("    config % seas_opt         = ",i0)') config % seas_opt
    write(6,'("    config % biomass_burn_opt = ",i0)') config % biomass_burn_opt

    if (config % chem_opt /= 300) return

    if (config % biomass_burn_opt > 0) then
    end if
#endif

  end subroutine gocart_init

  subroutine gocart_run(ktau, julday, tz, current_month, dt)
! subroutine gocart_run(clock)

    ! -- input variables
!   type(clock_type), intent(in) :: clock
    integer, intent(in) :: ktau
    integer, intent(in) :: julday
    integer, intent(in) :: tz
    integer, intent(in) :: current_month
    real,    intent(in) :: dt

    ! -- local variables
    logical, save :: firstfire = .true.
    logical :: call_gocart, call_plume, call_radiation, scale_fire_emiss
    logical :: store_arrays

!   integer :: ktau
!   integer :: julday
    integer :: nbegin, nv, nvv
    integer :: i, j, jp, jps, k
!   integer :: current_month, current_gmt, current_secs, current_msecs
    real(CHEM_KIND_R8) :: curr_secs

    real    :: factor, factor2
    real    :: dtstep, gmt
    real    :: dust_alpha,dust_gamma

    real, dimension(ims:ime, jms:jme)           :: dep_vel_o3, &
                                                   e_co,       &
                                                   raincv_b,   &
                                                   cu_co_ten,  &
                                                   dusthelp,   &
                                                   seashelp
    real, dimension(ims:ime, jms:jme, num_chem) :: var_rmv,    &
                                                   tr_fall

    ! -- begin
    print *,'gocart_run: entering ...', config % chem_opt, ims, ime, jms, jme, num_chem
    if (config % chem_opt == 0) return

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

    print *,'gocart_run: set control flags ...'
    print *,'gocart_run: set control flags ...', config % biomass_burn_opt
    print *,'gocart_run: set control flags ...', config % call_biomass
    print *,'gocart_run: set control flags ...', config % call_chemistry
    print *,'gocart_run: set control flags ...', config % call_radiation
    print *,'gocart_run: set control flags ...', ktau, firstfire

    ! -- get time & time step
    curr_secs = ktau * dt
    gmt = real(tz)

    ! -- set control flags
    call_plume       = (config % biomass_burn_opt > 0) .and. &
                      ((mod(ktau, config % call_biomass  ) == 0) .or. (ktau == 1) .or. firstfire)
    call_gocart      = (mod(ktau, config % call_chemistry) == 0) .or. (ktau == 1)
    call_radiation   = (mod(ktau, config % call_radiation) == 0) .or. (ktau == 1)
    scale_fire_emiss = .false.
    print *,'gocart_run: control flags set'

    print *,'gocart_run: get domain bounds ...', jts, jte
    ! -- start working
    if (ktau <= 1) then
      dtstep = dt
      rcav = rc2d
      rnav = rn2d-rc2d
    else
      dtstep = config % call_chemistry * dt
      rcav = max(0.,rc2d-rcav)
      rnav = max(0.,rn2d-rc2d-rnav)
    end if

    ! -- add time interval to initial date (to be replaced by ESMF_Clock
    ! -- this may not be needed

    ! -- get ready for chemistry run
    ! -- ttday and tcosz are computed -- should they be replaced with FV3 quantities?
    ! -- gmt is needed to compute quantities above
    ! -- get julday ----

    print *,'gocart_run: entering chem_prep ...'
    call chem_prep(config,ktau,dt,tr3d,tk3d,sm3d,                       &
                   ts2d,us2d,rsds,pr3d,emiss_ash_mass,emiss_ash_height, &
                   emiss_ash_dt,dm0,emiss_tr_mass,emiss_tr_height,      &
                   emiss_tr_dt,snwdph2d,VFRAC2d,VTYPE2d,STYPE2d,us3d,vs3d,ws3d,  &
                   slmsk2d,zorl2d,exch,pb2d,hf2d,th_pvsrf,oh_backgd,h2o2_backgd, &
                   no3_backgd,backg_oh,backg_h2o2,backg_no3,p_gocart,   &
                   nvl_chem, ttday,tcosz,gmt,julday,ph3d,area,ero1,   &
                   ero2,ero3,rcav,raincv_b,deg_lat,deg_lon,nvl,nvlp1,ntra, &
                   relhum,rri,t_phy,moist,u_phy,v_phy,p_phy,chem,tsk,ntrb, &
                   grvity,rd,p1000,cp,erod,emis_ant,emis_vol,e_co,dms_0,        &
                   u10,v10,ivgtyp,isltyp,gsw,vegfra,rmol,ust,znt,xland,dxy, &
                   t8w,p8w,exch_h,pbl,hfx,snowh,xlat,xlong,convfac,z_at_w,zmid,dz8w,vvel,&
                   rho_phy,smois,num_soil_layers,num_chem,num_moist,        &
                   emiss_abu,ebu_in,emiss_ab,num_ebu_in,num_emis_ant,       &
                   num_emis_vol,config % kemit,call_gocart,plumestuff, &
                   mean_fct_agtf,mean_fct_agef,mean_fct_agsv, &
                   mean_fct_aggr,firesize_agtf,firesize_agef, &
                   firesize_agsv,firesize_aggr, &
                   ids,ide, jds,jde, kds,kde, &
                   ims,ime, jms,jme, kms,kme, &
                   its,ite, jts,jte, kts,kte)

    ! -- compute sea salt
    print *,'gocart_run: entering sea salt ...'
    if(config % seas_opt == 1 )then
      call gocart_seasalt_driver(ktau,dt,rri,t_phy,moist, &
        u_phy,v_phy,chem,rho_phy,dz8w,u10,v10,p8w,        &
        xland,xlat,xlong,dxy,grvity,emis_seas,           &
        seashelp,num_emis_seas,num_moist,num_chem,        &
        ids,ide, jds,jde, kds,kde,                        &
        ims,ime, jms,jme, kms,kme,                        &
        its,ite, jts,jte, kts,kte)
    endif
    print *,'gocart_run: exit sea salt ...'

    print *,'gocart_run: check dust ...', config % dust_opt
    store_arrays = .false.
    select case (config % dust_opt)
      case (DUST_OPT_GOCART)
    print *,'gocart_run: dust is GOCART ...'
        call gocart_dust_driver(ktau,dt,rri,t_phy,moist,u_phy,         &
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
    store_arrays = store_arrays .and. (config % chem_opt >= 300)
    print *,'gocart_run: done dust ...'

    ! -- set output arrays
    if (store_arrays) then
      print *,'gocart_run: storing arrays ...'
      emi_d1(jts:jte) = emis_dust(its,1,jts:jte,1)
      emi_d2(jts:jte) = emis_dust(its,1,jts:jte,2)
      emi_d3(jts:jte) = emis_dust(its,1,jts:jte,3)
      emi_d4(jts:jte) = emis_dust(its,1,jts:jte,4)
      emi_d5(jts:jte) = emis_dust(its,1,jts:jte,5)
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

    if (config % dmsemis_opt == 1) then
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

    if ((config % dust_opt == DUST_OPT_GOCART) .or. &
        (config % dust_opt == DUST_OPT_AFWA  ) .or. &
        (config % seas_opt == 1)) then
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

    if (config % chem_opt == 316) then
      ! -- 10 volcanic size bins
      call vash_settling_driver(dt,t_phy,moist, &
         chem,rho_phy,dz8w,p8w,p_phy,dxy,      &
         ash_fall,grvity,num_moist,num_chem,    &
         ids,ide, jds,jde, kds,kde,             &
         ims,ime, jms,jme, kms,kme,             &
         its,ite, jts,jte, kts,kte)
      ashfall = ash_fall(its,:)
    else
      ! -- 4 volcanic size bins
      print *,'gocart_run: calling vashshort_settling_driver ...'
      call vashshort_settling_driver(dt,t_phy,moist, &
           chem,rho_phy,dz8w,p8w,p_phy,dxy,         &
           ash_fall,grvity,num_moist,num_chem,       &
           ids,ide, jds,jde, kds,kde,                &
           ims,ime, jms,jme, kms,kme,                &
           its,ite, jts,jte, kts,kte) 
      print *,'gocart_run: calling vashshort_settling_driver done'
      ashfall = ash_fall(its,:) !!!!!!!!!!!!!!!!!! WARNING: rewrites ashfall if chem_opt = 316
      print *,'gocart_run: ashfall done'
    end if

    ! -- add biomass burning emissions at every timestep
    if (config % biomass_burn_opt == 1) then
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
    if (config % chem_conv_tr == 2 )then
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
       e_co,config % kemit,snowh,numgas,                          &
       num_chem,num_moist,                                        &
       ids,ide, jds,jde, kds,kde,                                 &
       ims,ime, jms,jme, kms,kme,                                 &
       its,ite, jts,jte, kts,kte)
     print *,'gocart_run: done dry_dep_driver ...'

     ! -- ls wet deposition
     print *,'gocart_run: calling wetdep_ls ...'
     call wetdep_ls(dt,chem,rnav,moist,rho_phy,var_rmv,num_moist, &
         num_chem,numgas,p_qc,dz8w,vvel,config % chem_opt,        &
         ids,ide, jds,jde, kds,kde,                               &
         ims,ime, jms,jme, kms,kme,                               &
         its,ite, jts,jte, kts,kte)
     print *,'gocart_run: calling wetdep_ls ...'

    if (call_gocart) then
     print *,'gocart_run: calling GOCART CHEM driver ...'
      call gocart_chem_driver(ktau,dt,dtstep, gmt,julday,t_phy,moist, &
        chem,rho_phy,dz8w,p8w,backg_oh,oh_t,backg_h2o2,h2o2_t,backg_no3,no3_t, &
         dxy,grvity,xlat,xlong,ttday,tcosz, &
         config % chem_opt,num_chem,num_moist,                                      &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
     print *,'gocart_run: done GOCART CHEM driver'
     print *,'gocart_run: calling GOCART aerosols driver ...'
       call gocart_aerosols_driver(ktau,dtstep,t_phy,moist,  &
         chem,rho_phy,dz8w,p8w,dxy,grvity,         &
         config % chem_opt,num_chem,num_moist,                                      &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
     print *,'gocart_run: done GOCART aerosols driver'
    endif

    if (call_radiation) then
     print *,'gocart_run: calling radiation ...'
      store_arrays = .false.
      select case (config % aer_ra_feedback)
        case (1) 
     print *,'gocart_run: calling optical_driver ...'
          call optical_driver(curr_secs,dtstep,          &
               chem,dz8w,rri,relhum,                               &
               h2oai,h2oaj,                                        &
               tauaersw,gaersw,waersw,bscoefsw,tauaerlw,           &
               l2aer,l3aer,l4aer,l5aer,l6aer,l7aer,                &
               num_chem,config % chem_opt,ids,ide, jds,jde, kds,kde,        &
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
        ext_cof(kts:kte,jts:jte,1:nbands) = extt   (its,kts:kte,jts:jte,1:nbands)
        sscal  (kts:kte,jts:jte,1:nbands) = ssca   (its,kts:kte,jts:jte,1:nbands)
        asymp  (kts:kte,jts:jte,1:nbands) = asympar(its,kts:kte,jts:jte,1:nbands)
        aod2d  (jts:jte)                  = aod    (its,jts:jte)
      end if
    endif

    print *,'gocart_run: calling sum_pm_gocart ...'
    call sum_pm_gocart (                              &
         rri, chem,pm2_5_dry, pm2_5_dry_ec, pm10,     &
         num_chem,config % chem_opt,                  &
         ids,ide, jds,jde, kds,kde,                   &
         ims,ime, jms,jme, kms,kme,                   &
         its,ite, jts,jte, kts,kte)
    print *,'gocart_run: done sum_pm_gocart'

    ! -- pm25 and pm10 for output , not for tracer options
    print *,'gocart_run: setting output arrays ...'
    pm25  (kts:kte,jts:jte) = pm2_5_dry(its,kts:kte,jts:jte)
    p10   (kts:kte,jts:jte) = pm10     (its,kts:kte,jts:jte)
    ebu_oc(kts:kte,jts:jte) = ebu      (its,kts:kte,jts:jte,p_ebu_oc)
    print *,'gocart_run: done output arrays'

    print *,'gocart_run: call GOCART?',call_gocart
    if (call_gocart) then
        print *,'gocart_run: set bg fields ...'
        oh_bg  (kts:kte,jts:jte) = max(0., oh_t  (its,kts:kte,jts:jte))
        h2o2_bg(kts:kte,jts:jte) = max(0., h2o2_t(its,kts:kte,jts:jte))
        no3_bg (kts:kte,jts:jte) = max(0., no3_t (its,kts:kte,jts:jte))
        print *,'gocart_run: set bg fields done'
    end if

    ! -- put chem stuff back into tracer array
    print *,'gocart_run: put chem stuff back into tracer array ...'
    call update_tracers(tr3d, nvl, jms, jme, ntra, ntrb)
    print *,'gocart_run: done chem stuff back into tracer array'

  contains

    subroutine update_tracers(tr3d, nvl, jms, jme, ntra, ntrb)
      integer, intent(in) :: nvl, jms, jme, ntra, ntrb
      real, dimension(nvl,jms:jme,ntra+ntrb), intent(inout) :: tr3d

      ! -- local variables
      integer :: j, k, nv, nvv

      do nv = 1, num_chem
        nvv = nbegin + nv
        do j = jts, jte
          do k = kts, kte
            tr3d(k,j,nvv) = max(epsilc,chem(its,k,j,nv))
            trdp(k,j,nvv) = tr3d(k,j,nvv)*dp3d(k,j)
          end do
          wet_dep(j,nv) = var_rmv(its,j,nv)
        end do
      end do

    end subroutine update_tracers

  end subroutine gocart_run

end module gocart_mod
