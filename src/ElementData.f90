Module ElementData 
  use ModPrecision, only : fp_kind

  public :: bc_type, bc_vel, bc_disp
  enum, bind(c)
    enumerator :: bc_type = 0
    enumerator :: bc_vel = 157839
    enumerator :: bc_disp = 23097
  end enum
  !Usage 
  ! integer(kind(bc_type)) :: bc
  ! bc = bc_disp  

Type Element
  Integer, Dimension(:), Allocatable :: ID
  !GENERAL
  real(fp_kind), dimension(:), Allocatable :: h, t, cs, rho, m, rho_0, drhodt !influence radius, temp

  !THERMAL
  real(fp_kind), dimension(:), allocatable :: cp_t, k_t
  real(fp_kind), dimension(:,:,:), allocatable :: sigma, str_rate, rot_rate, shear_stress
  real(fp_kind), dimension(:), allocatable:: pressure, strain, mat_g
  
  real(fp_kind), dimension(:), allocatable :: sigma_eq !ONLY CALCULATED AT OUTPUT
  
  !Matrices --assembles or by gauss point...
  !Updated lagrangian formulation
  real(fp_kind), dimension(:,:,:), allocatable :: BL,BNL, jacob, x2,dHxy
  !! STIFFNESS MATRICES
  real(fp_kind), dimension(:,:,:), allocatable :: matKL, matKNL
  
  Integer, Dimension(:,:), Allocatable :: elnod !Connectivity
  real(fp_kind), dimension(:), allocatable :: detj
  
  Integer solver_type
 
  
End Type

!contains


End Module ElementData
