module chem_state_mod

  use chem_types_mod

  implicit none

  type chem_state_type

    real(CHEM_KIND_R8), dimension(:,:),     pointer :: area     => null()
    real(CHEM_KIND_R8), dimension(:,:),     pointer :: hf2d     => null()
    real(CHEM_KIND_R8), dimension(:,:),     pointer :: pb2d     => null()
    real(CHEM_KIND_R8), dimension(:,:),     pointer :: rc2d     => null()
    real(CHEM_KIND_R8), dimension(:,:),     pointer :: rn2d     => null()
    real(CHEM_KIND_R8), dimension(:,:),     pointer :: rsds     => null()
    real(CHEM_KIND_R8), dimension(:,:),     pointer :: slmsk2d  => null()
    real(CHEM_KIND_R8), dimension(:,:),     pointer :: snwdph2d => null()
    real(CHEM_KIND_R8), dimension(:,:),     pointer :: stype2d  => null()
    real(CHEM_KIND_R8), dimension(:,:),     pointer :: ts2d     => null()
    real(CHEM_KIND_R8), dimension(:,:),     pointer :: us2d     => null()
    real(CHEM_KIND_R8), dimension(:,:),     pointer :: vtype2d  => null()
    real(CHEM_KIND_R8), dimension(:,:),     pointer :: vfrac2d  => null()
    real(CHEM_KIND_R8), dimension(:,:),     pointer :: zorl2d   => null()

    real(CHEM_KIND_R8), dimension(:,:,:),   pointer :: exch     => null()
    real(CHEM_KIND_R8), dimension(:,:,:),   pointer :: ph3d     => null()
    real(CHEM_KIND_R8), dimension(:,:,:),   pointer :: phl3d    => null()
    real(CHEM_KIND_R8), dimension(:,:,:),   pointer :: pr3d     => null()
    real(CHEM_KIND_R8), dimension(:,:,:),   pointer :: prl3d    => null()
    real(CHEM_KIND_R8), dimension(:,:,:),   pointer :: sm3d     => null()
    real(CHEM_KIND_R8), dimension(:,:,:),   pointer :: tk3d     => null()
    real(CHEM_KIND_R8), dimension(:,:,:),   pointer :: us3d     => null()
    real(CHEM_KIND_R8), dimension(:,:,:),   pointer :: vs3d     => null()
    real(CHEM_KIND_R8), dimension(:,:,:),   pointer :: ws3d     => null()

    real(CHEM_KIND_R8), dimension(:,:,:,:), pointer :: tr3d     => null()

  end type chem_state_type

  public

end module chem_state_mod
