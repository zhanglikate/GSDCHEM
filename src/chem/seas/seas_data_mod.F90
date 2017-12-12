module seas_data_mod

  use chem_types_mod, only : CHEM_KIND_R8

  real(CHEM_KIND_R8), DIMENSION (4), PARAMETER :: ra        = (/   1.d-1,   5.d-1,   1.5d0,  5.0d0 /)
  real(CHEM_KIND_R8), DIMENSION (4), PARAMETER :: rb        = (/   5.d-1,   1.5d0,    5.d0,   1.d1 /)
  real(CHEM_KIND_R8), DIMENSION (4), PARAMETER :: den_seas  = (/   2.2d3,   2.2d3,   2.2d3,  2.2d3 /)
  real(CHEM_KIND_R8), DIMENSION (4), PARAMETER :: reff_seas = (/ 0.30D-6, 1.00D-6, 3.25D-6,7.50D-6 /)
  REAL(CHEM_KIND_R8), dimension(4,12)          :: ch_ss     = 1.0_CHEM_KIND_R8

end module seas_data_mod
