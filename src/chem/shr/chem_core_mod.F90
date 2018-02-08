module chem_core_mod

  use chem_rc_mod
  use chem_types_mod

  implicit none

  type chem_core_type
    real(CHEM_KIND_R4), dimension(:,:),     allocatable :: firesize_agef
    real(CHEM_KIND_R4), dimension(:,:),     allocatable :: firesize_aggr
    real(CHEM_KIND_R4), dimension(:,:),     allocatable :: firesize_agsv
    real(CHEM_KIND_R4), dimension(:,:),     allocatable :: firesize_agtf
    real(CHEM_KIND_R4), dimension(:,:),     allocatable :: mean_fct_agef
    real(CHEM_KIND_R4), dimension(:,:),     allocatable :: mean_fct_aggr
    real(CHEM_KIND_R4), dimension(:,:),     allocatable :: mean_fct_agsv
    real(CHEM_KIND_R4), dimension(:,:),     allocatable :: mean_fct_agtf
    real(CHEM_KIND_R4), dimension(:,:),     allocatable :: th_pvsrf
    real(CHEM_KIND_R4), dimension(:,:),     allocatable :: u10
    real(CHEM_KIND_R4), dimension(:,:),     allocatable :: v10

    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: convfac
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: dz8w
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: ebu_in
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: e_bio
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: erod
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: exch_h
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: backg_h2o2
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: backg_no3
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: backg_oh
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: p_phy
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: p8w
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: relhum
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: rho_phy
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: rri
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: smois
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: t_phy
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: t8w
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: u_phy
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: v_phy
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: vvel
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: z_at_w
    real(CHEM_KIND_R4), dimension(:,:,:),   allocatable :: zmid

    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: chem
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: emis_ant
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: emis_vol
    real(CHEM_KIND_R4), dimension(:,:,:,:), allocatable :: moist
  end type chem_core_type

  private

  public :: chem_core_type
  public :: chem_core_init

contains

  subroutine chem_core_init(core, is, ie, js, je, ks, ke, &
    kemit, ne_area, num_emis_ant, num_emis_dust, num_emis_seas, num_emis_vol, &
    num_chem, rc)

    type(chem_core_type), intent(inout) :: core
    integer,              intent(in)    :: is, ie, js, je, ks, ke
    integer,              intent(in)    :: kemit, ne_area, num_emis_ant, num_emis_dust, &
                                           num_emis_seas, num_emis_vol
    integer,              intent(in)    :: num_chem
    integer, optional,    intent(out)   :: rc

    ! -- local variables
    integer :: localrc
    character(len=*), parameter :: errmsg = "Failed to allocate memory for core arrays"

    ! -- begin
    if (present(rc)) rc = CHEM_RC_SUCCESS

    allocate(core % e_bio(is:ie, js:je, ne_area), stat=localrc)
    if (chem_rc_test((localrc /= 0), msg=errmsg, file=__FILE__, line=__LINE__, rc=rc)) return
    core % e_bio = 0._CHEM_KIND_R4

    allocate(core % backg_h2o2(is:ie, ks:ke, js:je), stat=localrc)
    if (chem_rc_test((localrc /= 0), msg=errmsg, file=__FILE__, line=__LINE__, rc=rc)) return
    core % backg_h2o2 = 0._CHEM_KIND_R4

    allocate(core % backg_no3(is:ie, ks:ke, js:je), stat=localrc)
    if (chem_rc_test((localrc /= 0), msg=errmsg, file=__FILE__, line=__LINE__, rc=rc)) return
    core % backg_h2o2 = 0._CHEM_KIND_R4

    allocate(core % backg_oh(is:ie, ks:ke, js:je), stat=localrc)
    if (chem_rc_test((localrc /= 0), msg=errmsg, file=__FILE__, line=__LINE__, rc=rc)) return
    core % backg_h2o2 = 0._CHEM_KIND_R4

    allocate(core % chem(is:ie, js:je, ks:ke, num_chem), stat=localrc)
    if (chem_rc_test((localrc /= 0), msg=errmsg, file=__FILE__, line=__LINE__, rc=rc)) return
    core % chem = 0._CHEM_KIND_R4

    allocate(core % emis_ant(is:ie, 1:kemit, js:je, num_emis_ant), stat=localrc)
    if (chem_rc_test((localrc /= 0), msg=errmsg, file=__FILE__, line=__LINE__, rc=rc)) return
    core % emis_ant = 0._CHEM_KIND_R4

    allocate(core % emis_vol(is:ie, 1:kemit, js:je, num_emis_vol), stat=localrc)
    if (chem_rc_test((localrc /= 0), msg=errmsg, file=__FILE__, line=__LINE__, rc=rc)) return
    core % emis_vol = 0._CHEM_KIND_R4

    allocate(core % z_at_w(is:ie, ks:ke, js:je), stat=localrc)
    if (chem_rc_test((localrc /= 0), msg=errmsg, file=__FILE__, line=__LINE__, rc=rc)) return
    core % z_at_w = 0._CHEM_KIND_R4

    allocate(core % dz8w(is:ie, ks:ke, js:je), stat=localrc)
    if (chem_rc_test((localrc /= 0), msg=errmsg, file=__FILE__, line=__LINE__, rc=rc)) return
    core % dz8w = 0._CHEM_KIND_R4


  end subroutine chem_core_init

end module chem_core_mod
