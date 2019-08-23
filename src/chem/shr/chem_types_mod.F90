module chem_types_mod

  implicit none

  integer, parameter :: CHEM_KIND_I4 = kind(1)
  integer, parameter :: CHEM_KIND_R4 = kind(1.0)
  integer, parameter :: CHEM_KIND_R8 = kind(1.d0)
  integer, parameter :: CHEM_KIND_C4 = kind(cmplx(1.,1.,kind=CHEM_KIND_R4))
  integer, parameter :: CHEM_KIND_C8 = kind(cmplx(1.,1.,kind=CHEM_KIND_R8))

  integer, parameter :: CHEM_MAXSTR  = 256

  public

end module chem_types_mod
