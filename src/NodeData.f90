!STRUCTURE OF ARRAYS TYPE
Module NodeData 
  use ModPrecision, only : fp_kind

Type Node
  Integer, Dimension(:), Allocatable :: ID
  real(fp_kind), dimension(:,:), Allocatable :: x, x_prev
  !GENERAL
  real(fp_kind), dimension(:), Allocatable :: h, t, cs, rho, m, rho_0, drhodt !influence radius, temp

  !THERMAL
  real(fp_kind), dimension(:), allocatable :: cp_t, k_t
  !Mechanical
  !real(fp_kind), dimension(:), allocatable :: cs
  real(fp_kind), dimension(:,:), allocatable :: u, v, a, disp, f_hour, u_inc ! u_inc is from last step
  real(fp_kind), dimension(:,:,:), allocatable :: sigma, str_rate, rot_rate, shear_stress
  real(fp_kind), dimension(:), allocatable:: pressure, strain, mat_g
  
  !!! FOR GRAPHIAL AVERAGING PURPOSE(IN REDUCED INTEGRATION)
  integer, dimension(:), Allocatable :: elxnod !!!Elements shared by each node
  integer, dimension(:,:), Allocatable :: nodel !!! element node 

  logical, dimension(:,:), Allocatable :: is_bcv, is_fix !0, none, else, dimension
  real(fp_kind), dimension(:,:), Allocatable :: bcv !Node and dim
  
  real(fp_kind), dimension(:), allocatable :: sigma_eq !ONLY CALCULATED AT OUTPUT
  
  Integer solver_type
 
  
End Type

End Module NodeData
