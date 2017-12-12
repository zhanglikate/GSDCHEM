module chem_domain_mod

  use mpp_mod,         only : mpp_npes
  use mpp_domains_mod, only : domain2d, mpp_define_domains, mpp_get_compute_domain
  use mpp_domains_mod, only : domain1d

  implicit none

  ! -- atmospheric grid parameters
  integer :: nip                 =  0      ! # of icosaedral cells
  integer :: nvl                 =  0      ! Number of vertical native levels
  integer :: nvlp                =  0      ! # of isobaric vertical levels - ex.  1000-25 hPa
  integer :: nvlp1               =  0      ! # of vertical levels ( = layers+1)

  integer, parameter :: nvl_chem =  55     ! Number of vertical native levels

  ! -- domain decomposition
  type(domain2d) :: chem_domain

  integer :: ids, ide, jds, jde, kds, kde
  integer :: ims, ime, jms, jme, kms, kme
  integer :: its, ite, jts, jte, kts, kte

  ! -- PET list
  integer, dimension(:), allocatable :: chem_pelist

  integer, allocatable :: indx(:), map(:)

  public

contains

  subroutine chem_domain_setup(domain)
    type(domain2d), optional, intent(in) :: domain

    type(domain1d) :: domain_1d
    ! -- begin

    if (present(domain)) then
      chem_domain = domain
    else
      call mpp_define_domains( (/ 1, nvl, 1, nip /), (/ 1, mpp_npes() /), chem_domain)
    end if

    ! set up trivial indices
    call mpp_get_compute_domain(chem_domain, kts, kte, jts, jte)

    ids =   1 ; ide =   1 ; jds =   1 ; jde = nip ; kds =   1 ; kde = nvl + 1
    ims =   1 ; ime =   1 ; jms = jts ; jme = jte ; kms =   1 ; kme = nvl + 1
    its =   1 ; ite =   1 

  end subroutine chem_domain_setup

end module chem_domain_mod
