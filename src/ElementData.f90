!STRUCTURE OF ARRAYS TYPE
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
  Integer, Dimension(:), Allocatable :: ID, gausspc
  !GENERAL
  !QUESTION: IS IT NECESARY TO STORE CS ON EACH ERLEMENT???
  real(fp_kind), dimension(:), Allocatable :: h, t, cs, mass, vol, drhodt, vol_0, vol_inc !influence radius, temp
  
  real(fp_kind), dimension(:,:), Allocatable :: rho_0,pressure,rho   !!!!elcount, gp
  !THERMAL
  real(fp_kind), dimension(:), allocatable :: cp_t, k_t
  real(fp_kind), dimension(:,:,:,:), allocatable :: sigma, str_rate, rot_rate , shear_stress,strain, str_inc,def_grad !tau is Cauchy Stress (do not confuse with shear)
  real(fp_kind), dimension(:,:,:,:), allocatable :: Umat,Rmat !!POLAR DECOMPOSITION, ONLY FOR GREEN NAGHDI
  real(fp_kind), dimension(:), allocatable:: mat_g
  
  real(fp_kind), dimension(:,:), allocatable :: c_s, p_visc !WAVE SPEED AND SHOCK VISCOSITY PRESSURE
  real(fp_kind), dimension(:), allocatable :: e_length
  real(fp_kind), dimension(:,:), allocatable :: sigma_eq !ONLY CALCULATED AT OUTPUT
  
  !Matrices --assembles or by gauss point...
  !Updated lagrangian formulation
  real(fp_kind), dimension(:,:,:,:), allocatable :: BL,BNL, jacob, dHxy,dHxy_detJ, dHxy0,math, dHrs !!!DIM: e,gp,,:,:, is it necesary to store dHrs??? is only because is used twice, at J and dHxy
  
  real(fp_kind), dimension(:,:,:), allocatable :: x2 !(rearranged nodes elem, x,y)
  
  real(fp_kind), dimension(:,:,:), allocatable :: hourg_nodf,f_int, f_ext !!!!HOURGLASS NODAL FORCES, Elem, node, dim
  

  !! STIFFNESS AND MASS MATRICES, ARE INTEGRATED MATRICES (NOT ON EACH GAUSS POINT)
  !!!!! MATMxDIM IS THE COMPLETE MATRIX (nodecount x dim ); USED FOR IMPLICIT PROBLEMS
  real(fp_kind), dimension(:,:,:), allocatable :: matKL, matKNL, matm, matmxdim,  uele,vele
  
  Integer, Dimension(:,:), Allocatable :: elnod,dof !Connectivity
  real(fp_kind), dimension(:,:), allocatable :: detj !(GAUSS POINT)
  
  Integer solver_type
 
  
End Type




End Module ElementData
