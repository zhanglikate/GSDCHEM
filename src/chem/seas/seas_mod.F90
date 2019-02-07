module seas_mod

  use chem_const_mod,   only : pi => pi_r8
  use chem_types_mod,   only : CHEM_KIND_R8
  use chem_tracers_mod, only : p_seas_1, p_seas_2, p_seas_3, p_seas_4, p_seas_5, &
                               config => chem_config
  use seas_data_mod

  implicit none

  integer, parameter :: SEAS_OPT_DEFAULT = 1

  private

  public :: SEAS_OPT_DEFAULT

  public :: gocart_seasalt_driver

   real, parameter    :: r80fac = 1.65     ! ratio of radius(RH=0.8)/radius(RH=0.) [Gerber]
   real, parameter    :: rhop = 2200.      ! dry seasalt density [kg m-3]
   !real, parameter    :: pi = 3.1415       ! ratio of circumference to diameterof circle
CONTAINS
  subroutine gocart_seasalt_driver(ktau,dt,alt,t_phy,moist,u_phy,  &
         v_phy,chem,rho_phy,dz8w,u10,v10,ustar,p8w,tsk,            &
         xland,xlat,xlong,area,g,emis_seas, &
         seashelp,num_emis_seas,num_moist,num_chem,seas_opt,  &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
!   USE module_initial_chem_namelists
! USE module_configure
! USE module_state_description
! USE module_model_constants, ONLY: mwdry
! IMPLICIT NONE
!  TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

     INTEGER,      INTENT(IN   ) :: ktau,num_emis_seas,num_moist,num_chem,   &
                                    ids,ide, jds,jde, kds,kde,               &
                                    ims,ime, jms,jme, kms,kme,               &
                                    its,ite, jts,jte, kts,kte,seas_opt
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
           INTENT(IN ) ::                                   moist
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
           INTENT(INOUT ) ::                                   chem
     REAL, DIMENSION( ims:ime, 1, jms:jme,num_emis_seas),OPTIONAL,&
           INTENT(INOUT ) ::                                                 &
           emis_seas
     REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &
            INTENT(IN   ) ::                                                 &
                                                       u10,                  &
                                                       v10,                  &
                                                       ustar,tsk,            &
                                                       xland,                &
                                                       xlat,                 &
                                                       xlong,area
     REAL,  DIMENSION( ims:ime , jms:jme ),                        &
            INTENT(OUT   ) :: seashelp
     REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
            INTENT(IN   ) ::                                                 &
                                                          alt,               &
                                                        t_phy,               &
                                                       dz8w,p8w,             &
                                                u_phy,v_phy,rho_phy

    REAL, INTENT(IN   ) :: dt,g
!
! local variables
!
    integer :: ipr,nmx,i,j,k,ndt,imx,jmx,lmx,n
    integer,dimension (1,1) :: ilwi
    real(CHEM_KIND_R8), DIMENSION (4) :: tc,bems
    real(CHEM_KIND_R8), dimension (1,1) ::w10m,uustar,gwet,airden,airmas,tskin
    real(CHEM_KIND_R8), dimension (1) :: dxy
    real(CHEM_KIND_R8) conver,converi,fsstemis,tskin_c

    real(CHEM_KIND_R8), dimension(1,1,1) :: airmas1
    real(CHEM_KIND_R8), dimension(1,1,1,4) :: tc1
    real(CHEM_KIND_R8), dimension(1,1,4) :: bems1
    
    conver=1.d-9
    converi=1.d9
!
! number of dust bins
!
    imx=1
    jmx=1
    lmx=1
    nmx=4
    k=kts
! p_seas_1=1
! p_seas_2=2
! p_seas_3=3
! p_seas_4=4
!    write(6,*)'call seasalt'
    select case (config % chem_opt)
      case (304, 316, 317)
        seashelp(:,:)=0.
        do j=jts,jte
          do i=its,ite
!
! donṫ do dust over water!!!
!
            if(xland(i,j).lt.0.5)then
              ilwi(1,1)=0
              tc(1)=chem(i,kts,j,p_seas_1)*conver
              tc(2)=1.d-30
              tc(3)=chem(i,kts,j,p_seas_2)*conver
              tc(4)=1.d-30
              w10m(1,1)=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
              uustar(1,1)=ustar(i,j)
              tskin(1,1)=tsk(i,j)
              airmas(1,1)=-(p8w(i,kts+1,j)-p8w(i,kts,j))*area(i,j)/g
!
! we donṫ trust the u10,v10 values, is model layers are very thin near surface
!
              if(dz8w(i,kts,j).lt.12.)w10m=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))
!
              dxy(1)=area(i,j)
              ipr=0

!    if(j.eq.2478)write(6,*)'call seasalt ',w10m(1,1),airmas(1,1),tc(1),tc(2)
              airmas1(1,1,1) = airmas(1,1)
              tc1(1,1,1,:) = tc
              bems1(1,1,:) = bems
!             call source_ss( imx,jmx,lmx,nmx, dt, tc,ilwi, dxy, w10m, airmas, bems,ipr)
              call source_ss( imx,jmx,lmx,nmx, dt, tc1,ilwi, dxy, w10m, airmas1, bems1,ipr)
              tc = tc1(1,1,1,:)
!    if(j.eq.2558)write(6,*)'call seasalt after',tc(1),tc(2),bems(1)
!    write(6,*)'call seasalt after',tc(1),tc(2),bems(1)
              chem(i,kts,j,p_seas_1)=(tc(1)+.75*tc(2))*converi
              chem(i,kts,j,p_seas_2)=(tc(3)+.25*tc(2))*converi
              seashelp(i,j)=tc(2)*converi
! for output diagnostics
!    emis_seas(i,1,j,p_seas_1)=bems(1)
!    emis_seas(i,1,j,p_seas_2)=bems(2)
!    emis_seas(i,1,j,p_seas_3)=bems(3)
            endif
          enddo
        enddo

      case default
        do j=jts,jte
          do i=its,ite
!
! donṫ do dust over water!!!
!
            if(xland(i,j).lt.0.5)then
              ilwi(1,1)=0
              tc(1)=chem(i,kts,j,p_seas_1)*conver
              tc(2)=chem(i,kts,j,p_seas_2)*conver
              tc(3)=chem(i,kts,j,p_seas_3)*conver
              tc(4)=chem(i,kts,j,p_seas_4)*conver
              w10m(1,1)=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
              uustar(1,1)=ustar(i,j)
              tskin(1,1)=tsk(i,j)
              airmas(1,1)=-(p8w(i,kts+1,j)-p8w(i,kts,j))*area(i,j)/g
!NGAC SST correction:

              fsstemis = 0.0

              tskin_c  = tskin(1,1) - 273.15

              if(tskin_c < -0.1) tskin_c = -0.1    ! temperature range (0, 36) C 
              if(tskin_c > 36.0) tskin_c = 36.0    !

           fsstemis = (-1.107211 -0.010681*tskin_c -0.002276*tskin_c**2 + &
                      60.288927*1.0/(40.0 - tskin_c))
           if(fsstemis < 0.0) fsstemis = 0.0
           if(fsstemis > 7.0) fsstemis = 7.0

!------------------------------------------------------------------------------
!
! we donṫ trust the u10,v10 values, is model layers are very thin near surface
!
              if(dz8w(i,kts,j).lt.12.)w10m=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))
!
              dxy(1)=area(i,j)
              ipr=0

!    if(j.eq.2478)write(6,*)'call seasalt ',w10m(1,1),airmas(1,1),tc(1),tc(2)
              airmas1(1,1,1) = airmas(1,1)
              tc1(1,1,1,:) = tc
              bems1(1,1,:) = bems
            if (seas_opt == 1) then
              call source_ss( imx,jmx,lmx,nmx, dt, tc1,ilwi, dxy, w10m, airmas1, bems1,ipr)
            endif

            if (seas_opt == 2) then            
            do n=1,nmx
              call SeasaltEmission ( imx,jmx,lmx,dt,tc1(1,1,1,n),ra(n), rb(n), 3, dxy,w10m,uustar, &
                                airmas1, bems1(1,1,n), fsstemis,ipr )
            enddo !nmx
            endif !over sea

              tc = tc1(1,1,1,:)
!    if(j.eq.2558)write(6,*)'call seasalt after',tc(1),tc(2),bems(1)
!    write(6,*)'call seasalt after',tc(1),tc(2),bems(1)
              chem(i,kts,j,p_seas_1)=tc(1)*converi
              chem(i,kts,j,p_seas_2)=tc(2)*converi
              chem(i,kts,j,p_seas_3)=tc(3)*converi
              chem(i,kts,j,p_seas_4)=tc(4)*converi
! for output diagnostics
    emis_seas(i,1,j,p_seas_1)=bems(1)
    emis_seas(i,1,j,p_seas_2)=bems(2)
    emis_seas(i,1,j,p_seas_3)=bems(3)
    emis_seas(i,1,j,p_seas_4)=bems(4)
            endif



          enddo
        enddo

    end select

  end subroutine gocart_seasalt_driver

  SUBROUTINE source_ss(imx,jmx,lmx,nmx, dt1, tc, &
                       ilwi, dxy, w10m, airmas, &
                       bems,ipr)

! ****************************************************************************
! *  Evaluate the source of each seasalt particles size classes  (kg/m3) 
! *  by soil emission.
! *  Input:
! *         SSALTDEN  Sea salt density                               (kg/m3)
! *         DXY       Surface of each grid cell                     (m2)
! *         NDT1      Time step                                     (s)
! *         W10m      Velocity at the anemometer level (10meters)   (m/s)
! *      
! *  Output:
! *         DSRC      Source of each sea salt bins       (kg/timestep/cell) 
! *
! *
! * Number flux density: Original formula by Monahan et al. (1986) adapted
! * by Sunling Gong (JGR 1997 (old) and GBC 2003 (new)).  The new version is
! * to better represent emission of sub-micron sea salt particles.
!
! * dFn/dr = c1*u10**c2/(r**A) * (1+c3*r**c4)*10**(c5*exp(-B**2))
! * where B = (b1 -log(r))/b2
! * see c_old, c_new, b_old, b_new below for the constants.
! * number fluxes are at 80% RH.
! *
! * To calculate the flux:
! * 1) Calculate dFn based on Monahan et al. (1986) and Gong (2003)
! * 2) Assume that wet radius r at 80% RH = dry radius r_d *frh
! * 3) Convert particles flux to mass flux :
! *    dFM/dr_d = 4/3*pi*rho_d*r_d^3 *(dr/dr_d) * dFn/dr
! *             = 4/3*pi*rho_d*r_d^3 * frh * dFn/dr
! *               where rho_p is particle density [kg/m3]
! *    The factor 1.e-18 is to convert in micro-meter r_d^3
! ****************************************************************************
   

    IMPLICIT NONE

    INTEGER, INTENT(IN)    :: nmx,imx,jmx,lmx,ipr
    INTEGER, INTENT(IN)    :: ilwi(imx,jmx)
    REAL(CHEM_KIND_R8),    INTENT(IN)    :: dxy(jmx), w10m(imx,jmx)
    REAL(CHEM_KIND_R8),    INTENT(IN)    :: airmas(imx,jmx,lmx)
    REAL(CHEM_KIND_R8),    INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
    REAL(CHEM_KIND_R8),    INTENT(OUT)   :: bems(imx,jmx,nmx)

    REAL(CHEM_KIND_R8) :: c0(5), b0(2)
!  REAL(CHEM_KIND_R8), PARAMETER :: c_old(5)=(/1.373, 3.41, 0.057, 1.05, 1.190/) 
!  REAL(CHEM_KIND_R8), PARAMETER :: c_new(5)=(/1.373, 3.41, 0.057, 3.45, 1.607/)
    ! Change suggested by MC
    REAL(CHEM_KIND_R8), PARAMETER :: c_old(5)=(/1.373, 3.2, 0.057, 1.05, 1.190/) 
    REAL(CHEM_KIND_R8), PARAMETER :: c_new(5)=(/1.373, 3.2, 0.057, 3.45, 1.607/)
    REAL(CHEM_KIND_R8), PARAMETER :: b_old(2)=(/0.380, 0.650/)
    REAL(CHEM_KIND_R8), PARAMETER :: b_new(2)=(/0.433, 0.433/)
    REAL(CHEM_KIND_R8), PARAMETER :: dr=5.0D-2 ! um   
    REAL(CHEM_KIND_R8), PARAMETER :: theta=30.0
    ! Swelling coefficient frh (d rwet / d rd)
!!!  REAL(CHEM_KIND_R8),    PARAMETER :: frh = 1.65
    REAL(CHEM_KIND_R8),    PARAMETER :: frh = 2.d0
    LOGICAL, PARAMETER :: old=.TRUE., new=.FALSE.
    REAL(CHEM_KIND_R8) :: rho_d, r0, r1, r, r_w, a, b, dfn, r_d, dfm, src
    INTEGER :: i, j, n, nr, ir
    REAL :: dt1,fudge_fac


    REAL(CHEM_KIND_R8)    :: tcmw(nmx), ar(nmx), tcvv(nmx)
    REAL(CHEM_KIND_R8)    :: ar_wetdep(nmx), kc(nmx)
    CHARACTER(LEN=20)     :: tcname(nmx), tcunits(nmx)
    LOGICAL               :: aerosol(nmx)


    REAL(CHEM_KIND_R8) :: tc1(imx,jmx,lmx,nmx)
    REAL(CHEM_KIND_R8), TARGET :: tcms(imx,jmx,lmx,nmx) ! tracer mass (kg; kgS for sulfur case)
    REAL(CHEM_KIND_R8), TARGET :: tcgm(imx,jmx,lmx,nmx) ! g/m3

    !-----------------------------------------------------------------------  
    ! sea salt specific
    !-----------------------------------------------------------------------  
! REAL(CHEM_KIND_R8), DIMENSION(nmx) :: ra, rb
! REAL(CHEM_KIND_R8) :: ch_ss(nmx,12)

    !-----------------------------------------------------------------------  
    ! emissions (input)
    !-----------------------------------------------------------------------  
    REAL(CHEM_KIND_R8) :: e_an(imx,jmx,2,nmx), e_bb(imx,jmx,nmx), &
            e_ac(imx,jmx,lmx,nmx)

    !-----------------------------------------------------------------------  
    ! diagnostics (budget)
    !-----------------------------------------------------------------------
!  ! tendencies per time step and process
!  REAL(CHEM_KIND_R8), TARGET :: bems(imx,jmx,nmx), bdry(imx,jmx,nmx), bstl(imx,jmx,nmx)
!  REAL(CHEM_KIND_R8), TARGET :: bwet(imx,jmx,nmx), bcnv(imx,jmx,nmx)!

!  ! integrated tendencies per process
!  REAL(CHEM_KIND_R8), TARGET :: tems(imx,jmx,nmx), tstl(imx,jmx,nmx)
!  REAL(CHEM_KIND_R8), TARGET :: tdry(imx,jmx,nmx), twet(imx,jmx,nmx), tcnv(imx,jmx,nmx)

    ! global mass balance per time step 
    REAL(CHEM_KIND_R8) :: tmas0(nmx), tmas1(nmx)
    REAL(CHEM_KIND_R8) :: dtems(nmx), dttrp(nmx), dtdif(nmx), dtcnv(nmx)
    REAL(CHEM_KIND_R8) :: dtwet(nmx), dtdry(nmx), dtstl(nmx)
    REAL(CHEM_KIND_R8) :: dtems2(nmx), dttrp2(nmx), dtdif2(nmx), dtcnv2(nmx)
    REAL(CHEM_KIND_R8) :: dtwet2(nmx), dtdry2(nmx), dtstl2(nmx)

    ! detailed integrated budgets for individual emissions
    REAL(CHEM_KIND_R8), TARGET :: ems_an(imx,jmx,nmx),    ems_bb(imx,jmx,nmx), ems_tp(imx,jmx)
    REAL(CHEM_KIND_R8), TARGET :: ems_ac(imx,jmx,lmx,nmx)
    REAL(CHEM_KIND_R8), TARGET :: ems_co(imx,jmx,nmx)

    ! executable statements
! decrease seasalt emissions (Colarco et al. 2010)
!
    !fudge_fac= 1. !.5
    !fudge_fac= .5 !lzhang
    fudge_fac= .25 !lzhang
!
    DO n = 1,nmx
!    if(ipr.eq.1)write(0,*)'in seasalt',n,ipr,ilwi
       bems(:,:,n) = 0.0
       rho_d = den_seas(n)
       r0 = ra(n)*frh
       r1 = rb(n)*frh
       r = r0
       nr = INT((r1-r0)/dr+.001)
!    if(ipr.eq.1.and.n.eq.1)write(0,*)'in seasalt',nr,r1,r0,dr,rho_d
       DO ir = 1,nr
          r_w = r + dr*0.5
          r = r + dr
          IF (new) THEN
             a = 4.7*(1.0 + theta*r_w)**(-0.017*r_w**(-1.44))
             c0 = c_new
             b0 = b_new
          ELSE
             a = 3.0
             c0 = c_old
             b0 = b_old
          END IF
          !
          b = (b0(1) - LOG10(r_w))/b0(2)
          dfn = (c0(1)/r_w**a)*(1.0 + c0(3)*r_w**c0(4))* &
               10**(c0(5)*EXP(-(b**2)))
          
          r_d = r_w/frh*1.0D-6  ! um -> m
          dfm = 4.0/3.0*pi*r_d**3*rho_d*frh*dfn*dr*dt1 ! 3600 !dt1
          DO i = 1,imx
             DO j = 1,jmx
!              IF (water(i,j) > 0.0) THEN
                IF (ilwi(i,j) == 0) THEN
!                 src = dfm*dxy(j)*water(i,j)*w10m(i,j)**c0(2)
                   src = dfm*dxy(j)*w10m(i,j)**c0(2)
!                 src = ch_ss(n,dt(1)%mn)*dfm*dxy(j)*w10m(i,j)**c0(2)
                   tc(i,j,1,n) = tc(i,j,1,n) + fudge_fac*src/airmas(i,j,1)
!                if(ipr.eq.1)write(0,*)n,dfm,c0(2),dxy(j),w10m(i,j),src,airmas(i,j,1)
                ELSE
                   src = 0.0
                END IF
                bems(i,j,n) = bems(i,j,n) + src*fudge_fac/(dxy(j)*dt1) !kg/m2/s
             END DO  ! i
          END DO ! j
       END DO ! ir
    END DO ! n

  END SUBROUTINE source_ss

   
   subroutine SeasaltEmission ( imx,jmx,lmx,dt,tc,rLow, rUp, method,dxy, w10m, ustar, &
                                airmas, memissions,fsstemis,rc )

! !DESCRIPTION: Calculates the seasalt mass emission flux every timestep.
!  The particular method (algorithm) used for the calculation is based
!  on the value of "method" passed on input.  Mostly these algorithms are
!  a function of wind speed and particle size (nominally at 80% RH).
!  Routine is called once for each size bin, passing in the edge radii
!  "rLow" and "rUp" (in dry radius, units of um).  Returned in the emission
!  mass flux [kg m-2 s-1].  A sub-bin assumption is made to break (possibly)
!  large size bins into a smaller space.
!
! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   real(CHEM_KIND_R8), intent(in)                  :: rLow, rUp,fsstemis   ! Dry particle bin edge radii [um]
   INTEGER, INTENT(IN)                             :: imx,jmx,lmx     !1d 
   real(CHEM_KIND_R8), intent(in)                  :: dxy(jmx)     !grid cell[m-2] 
   real(CHEM_KIND_R8), intent(in)          :: w10m(imx,jmx)   ! 10-m wind speed [m s-1]
   real(CHEM_KIND_R8), intent(in)          :: ustar(imx,jmx)  ! friction velocity [m s-1]
   real(CHEM_KIND_R8), intent(in)          :: airmas(imx,jmx,lmx)  ! 
   REAL(CHEM_KIND_R8),    INTENT(INOUT) :: tc(imx,jmx,lmx)
   REAL(CHEM_KIND_R8),    INTENT(INOUT) :: memissions(imx,jmx)! Mass EmissionsFlux [kg m-2 s-1]
   integer, intent(in)               :: method      ! Algorithm to use
   real, intent(in)                        :: dt                   !timestep[s]

! !OUTPUT PARAMETERS:

   real  nemissions      ! Number Emissions Flux [m-2 s-1]
   integer, intent(out)          :: rc              ! Error return code:
                                                    !  0 - all is well
                                                    !  1 - 
! !Local Variables
   integer       :: ir,i,j
   real          :: w                          ! Intermediary wind speed [ms-1]
   real          :: r, dr                           ! sub-bin radius spacing(dry, um)
   real          :: rwet, drwet                     ! sub-bin radius spacing(rh=80%, um)
   real          :: aFac, bFac, scalefac, rpow, exppow, wpow

   integer, parameter :: nr = 10                    ! Number of (linear)sub-size bins
   real, parameter :: emission_scale= 0.875      ! emission scale in NGAC namelist

!   character(len=*), parameter :: myname = 'SeasaltEmission'

!  Define the sub-bins (still in dry radius)
   !dr = (rUp - rLow)/nr
   dr = (rUp - rLow)/nr
   r  = rLow + 0.5*dr

    do i=1,imx
      do j=1,jmx
!  Loop over size bins
   nemissions = 0.
   memissions = 0.

   do ir = 1, nr
    rwet  = r80fac * r
    drwet = r80fac * dr
    select case(method)

     case(1)  ! Gong 2003
      aFac     = 4.7*(1.+30.*rwet)**(-0.017*rwet**(-1.44))
      bFac     = (0.433-log10(rwet))/0.433
      scalefac = 1.
      rpow     = 3.45
      exppow   = 1.607
      wpow     = 3.41
      w        = w10m(i,j)

     case(2)  ! Gong 1997
      aFac     = 3.
      bFac     = (0.380-log10(rwet))/0.650
      scalefac = 1.
      rpow     = 1.05
      exppow   = 1.19
      wpow     = 3.41
      w        = w10m (i,j)

     case(3)  ! GEOS5 2012
      aFac     = 4.7*(1.+30.*rwet)**(-0.017*rwet**(-1.44))
      bFac     = (0.433-log10(rwet))/0.433
      scalefac = 33.0e3
      rpow     = 3.45
      exppow   = 1.607
      wpow     = 3.41 - 1.
      w        = ustar (i,j)

     case default
      !if(MAPL_AM_I_ROOT()) print *, 'SeasaltEmission missing algorithm method'
      rc = 1
      return

    end select


!   Number emissions flux (# m-2 s-1)
    nemissions = nemissions + SeasaltEmissionGong( rwet, drwet, w, scalefac,aFac, bFac, rpow, exppow, wpow )
!   Mass emissions flux (kg m-2 s-1)
    scalefac = scalefac * 4./3.*pi*rhop*r**3.*1.e-18
    memissions(i,j) = memissions(i,j) + emission_scale*fsstemis*SeasaltEmissionGong( rwet, drwet, w, scalefac,aFac, bFac, rpow, exppow, wpow )

    tc(i,j,1)=memissions(i,j)*dxy(j)*dt/airmas(i,j,1)  !kg/kg
    r = r + dr

   end do

    enddo !j
    enddo !i
   rc = 0

  end subroutine SeasaltEmission

! Function to compute sea salt emissions following the Gong style
! parameterization.  Functional form is from Gong 2003:
!  dN/dr = scalefac * 1.373 * (w^wpow) * (r^-aFac) * (1+0.057*r^rpow) *
!  10^(exppow*exp(-bFac^2))
! where r is the particle radius at 80% RH, dr is the size bin width at 80% RH,
! and w is the wind speed

  function SeasaltEmissionGong ( r, dr, w, scalefac, aFac, bFac, rpow, exppow, wpow )

   real, intent(in)    :: r, dr     ! Wet particle radius, bin width [um]
   real, intent(in)    :: w    ! Grid box mean wind speed [m s-1](10-m or ustar wind)
   real, intent(in)    :: scalefac, aFac, bFac, rpow, exppow, wpow
   real                :: SeasaltEmissionGong

!  Initialize
   SeasaltEmissionGong = 0.

!  Particle size distribution function
   SeasaltEmissionGong = scalefac * 1.373*r**(-aFac)*(1.+0.057*r**rpow) &
                         *10**(exppow*exp(-bFac**2.))*dr
!  Apply wind speed function
   SeasaltEmissionGong = w**wpow * SeasaltEmissionGong

  end function SeasaltEmissionGong


end module seas_mod
