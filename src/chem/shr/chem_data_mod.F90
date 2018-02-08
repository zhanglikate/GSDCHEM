module chem_data_mod

  use chem_types_mod

  implicit none

  type chem_data_type
    ! -- input
    real(CHEM_KIND_R4), dimension(:),     allocatable :: p_gocart          ! GOCART pressure levels
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: clayfrac          ! clay fraction (AFWA dust scheme)
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: dm0               ! dms reference emissions
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: emiss_ab          ! emissions for all available species
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: emiss_abu         ! emissions for all available species
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: emiss_ash_dt      ! ash emissions
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: emiss_ash_height  ! ash emissions
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: emiss_ash_mass    ! ash emissions
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: emiss_tr_dt      ! ash emissions
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: emiss_tr_height  ! ash emissions
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: emiss_tr_mass    ! ash emissions
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: ero1              ! dust erosion factor
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: ero2              ! dust erosion factor
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: ero3              ! dust erosion factor
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: h2o2_backgd       ! H2O2 background for GOCART
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: no3_backgd        ! NO3 background for GOCART
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: oh_backgd         ! OH background for GOCART
    real(CHEM_KIND_R4), dimension(:,:,:), allocatable :: plumestuff        ! fire info
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: sandfrac          ! sand fraction (AFWA dust scheme)
    real(CHEM_KIND_R4), dimension(:,:),   allocatable :: th_pvsrf
    ! -- output
    real(CHEM_KIND_R4), dimension(:,:),     allocatable :: emi_d1
    real(CHEM_KIND_R4), dimension(:,:),     allocatable :: emi_d2
    real(CHEM_KIND_R4), dimension(:,:),     allocatable :: emi_d3
    real(CHEM_KIND_R4), dimension(:,:),     allocatable :: emi_d4
    real(CHEM_KIND_R4), dimension(:,:),     allocatable :: emi_d5
    real(CHEM_KIND_R4), dimension(:,:),     allocatable :: aod2d
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: pm10
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: pm25
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: ebu_oc
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: oh_bg
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: h2o2_bg
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: no3_bg
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: wet_dep
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: ext_cof
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: sscal
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: asymp
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: tr3d
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: trdp
  end type chem_data_type

contains

end module chem_data_mod
