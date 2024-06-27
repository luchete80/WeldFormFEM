
module Domain
use NodeData
use ElementData
use ModPrecision, only : fp_kind
!use Neighbor 

implicit none 

REAL, PARAMETER :: PI = 3.14159265358979323846264338327950288419716939937510

type Dom_type
integer, dimension (:), allocatable ::slavenod

real(fp_kind):: mat_cs0, mat_K, mat_E, mat_G, mat_nu, mat_Et !IF ET IS CONSTANT (BILINEAR)
end Type Dom_type

integer :: Dim, node_count, elem_count, Nproc, dof 
!!!!! THIS SHOULD BE INSIDE OF DOMAIN
type(Node)	::nod
type(Element)	::elem
real(fp_kind):: tot_mass
real(fp_kind), dimension(:,:), Allocatable :: mat_C !TODO: CHANGE TO SEVERAL MATERIALS
real(fp_kind), dimension(:,:), Allocatable :: kglob, uglob, m_glob

integer, dimension(:,:), Allocatable :: cont_nodes(:) !!! ELEMENT AND SEGMENT

integer :: bind_dom_type  !1 plane stress, 2 plain strain,. 3 axisymm
logical :: axisymm_vol_weight !! PREFERABLE AREA WEIGHT 
! enum, bind(C) :: plane_mode_en
   ! enumerator :: pl_stress, pl_strain, axi_sym
! end type
! type(plane_mode_en) :: plane_mode


!THESE ARE VECTORS, NOT MATRICES (AND ARE NOT MULTIPLIED)
real(fp_kind), dimension(:,:), Allocatable :: rint_glob, fext_glob !Accelerations and internal forces
integer :: nodxelem !TODO: SET TO ELEMEENT VAR


real(fp_kind), dimension(3):: dommax, dommin

real(fp_kind)::mat_G, mat_K, mat_cs0, mat_E!TODO: change to material
real(fp_kind):: time
logical :: calc_m
contains 

  

  subroutine DomInit(proc)
    integer, intent(in)::proc
    nproc = proc
    Dim = 3
    DomMax (:) = -100000000000.0;
    DomMin (:) = 100000000000.0;
  end subroutine DomInit
  
  !After initialization of nodes and elements
  subroutine AllocateDomain()
    dof = node_count*dim
    
    allocate (kglob(node_count*dim,node_count*dim))
    allocate (uglob(node_count*dim,1))
    allocate (m_glob(node_count,node_count))
    
    allocate (rint_glob(node_count,dim))
    allocate (fext_glob(node_count,dim))
    
  end subroutine AllocateDomain
  
  subroutine AllocateNodes(pt_count)
    integer, intent(in):: pt_count
    
    node_count = pt_count
    !!!GENERAL 
    allocate (nod%x(node_count,dim))
    allocate (nod%x_prev(node_count,dim))
    allocate (nod%u_inc(node_count,dim))
    ! allocate (nod%rho(node_count))
    ! allocate (nod%drhodt(node_count))
    ! allocate (nod%h(node_count))
    allocate (nod%m(node_count))
    allocate (nod%id(node_count))
    
    allocate (nod%elxnod(node_count))
    if (dim .eq. 2) then 
      allocate (nod%nodel(node_count, 4))
    else 
      allocate (nod%nodel(node_count, 8))
    end if
    allocate (nod%rho(node_count))
    
    
    !! THERMAL
    ! allocate (nod%cp_t(node_count))
    ! allocate (nod%t(node_count))
    ! allocate (nod%k_t(node_count))
    
    !MECHANICAL PROPERTIES
    !if (solver_type==1) then 
    allocate (nod%u(node_count,dim))
    allocate (nod%v(node_count,dim))
    allocate (nod%a(node_count,dim))
    allocate (nod%f_hour(node_count,dim))
    allocate (nod%disp(node_count,dim))
    allocate (nod%sigma_eq(node_count))

    
    allocate (nod%sigma(node_count,3,3))

    allocate (nod%pl_strain(node_count))
    ! allocate (nod%str_rate(node_count,3,3))
    ! allocate (nod%rot_rate(node_count,3,3))
    ! allocate (nod%shear_stress(node_count,3,3))
    ! allocate (nod%strain(node_count))
    ! allocate (nod%pressure(node_count))
    ! allocate (nod%cs(node_count))
    

    

    !!!!! BOUNDARY CONDITIONS
    allocate (nod%is_bcv(node_count,dim))
    allocate (nod%is_fix(node_count,dim))
    allocate (nod%bcv(node_count,dim))
    !end if  
  end subroutine
  
  subroutine AllocateElements(el_count, gp) !gauss point count
    integer, intent(in):: el_count, gp
    
    print *, "Creating ", el_count, " elements"
    elem_count = el_count
    print *, "Node count per element: ", nodxelem
    allocate (elem%elnod(el_count,nodxelem))
    allocate (elem%gausspc(el_count))
    allocate (elem%dof(el_count,dim*nodxelem))
    allocate (elem%vol(el_count))
    allocate (elem%vol_inc(el_count))
    allocate (elem%vol_0(el_count))
    allocate (elem%x2(el_count,nodxelem,dim))
    allocate (elem%jacob(el_count,gp,dim,dim))
    allocate (elem%detj(el_count,gp))
    allocate (elem%sigma_eq(el_count,gp)) !But is constant??
    allocate (elem%sigma_y(el_count,gp)) !But is constant??
    allocate (elem%pl_strain(el_count,gp)) !But is constant??
    elem%pl_strain(:,:) = 0.0d0 
    !IF NOT YIELDING THING
    elem%sigma_y(:,:) = 1.0e20 !But is constant??
        
    allocate (elem%dHxy(el_count,gp,dim,nodxelem))
    allocate (elem%dHxy_detJ(el_count,gp,dim,nodxelem)) !!!! STORE LIKE THIS TO SAVE CALCULATION TIME (THIS IS USED  TO CALC FORCES INTEGRATING IT )
    allocate (elem%dHxy0(el_count,gp,dim,nodxelem)) !!!USED FOR DEFORMATION GRADIENT ONLY FOR FULL INTEGRATION ELEMENTS 
    allocate (elem%dHrs(el_count,gp,dim,nodxelem))

    
    allocate (elem%B_ax(el_count,gp,dim,nodxelem))

    allocate (elem%uele (el_count,dim*nodxelem,1)) 

    allocate (elem%vele (el_count,dim*nodxelem,1)) 
    
    allocate (elem%mass(el_count)) !Mass matrix    
    
    allocate (elem%c_s(el_count,gp))
    allocate (elem%p_visc(el_count,gp))
    allocate (elem%e_length(el_count))

    allocate (elem%matm(el_count,nodxelem,nodxelem)) !Mass matrix
    allocate (elem%math(el_count,gp,1,nodxelem)) !Mass matrix
    
    allocate (elem%hourg_nodf(el_count,nodxelem,dim)) !AS 1 COLUMN OR NOT????? Mass matrix
    
    allocate (elem%f_int(el_count,nodxelem,dim))
    allocate (elem%f_ext(el_count,nodxelem,dim))
    
    allocate (elem%rho(el_count,gp)) !AT FIRST ONLY ONE POINT
    allocate (elem%rho_0(el_count,gp))
    allocate (elem%pressure(el_count,gp))
 
    !!---
    allocate (elem%radius(el_count,gp)) !!ONLY IN AXISYMM
    allocate (elem%sigma_tg(el_count,gp)) !!!TANGENTIAL (HOOP) VALUES
    allocate (elem%str_tg(el_count,gp))
    !! --AXISYMM
    
    allocate (elem%cs(el_count))
    !!!! ORIGINAL DIMENSIONS WAS DIM
    allocate (elem%shear_stress(el_count,gp, 3,3))
    allocate (elem%str_rate(el_count,gp, 3,3))
    allocate (elem%str_inc(el_count,gp, 3,3))
    allocate (elem%str_tot(el_count,gp, 3,3))
    allocate (elem%rot_rate(el_count,gp, 3,3))
    allocate (elem%sigma(el_count,gp,3,3))  !!!THIS IS A DIMxDIM SYMMETRIC TENSOR
    
    allocate (elem%def_grad(el_count,gp, dim,dim))
    
    !!! ONLY IF POLAR DECOMP
    allocate (elem%umat(el_count,gp, dim,dim)) !TODO; CHANGE TO SYMM
    allocate (elem%rmat(el_count,gp, dim,dim))
    
    if (Dim .eq. 2) then
      allocate (elem%bl (el_count,gp,3,dim*nodxelem))
      allocate (elem%bnl(el_count,gp, 4,dim*nodxelem))
      allocate (elem%strain(el_count,gp, 4,1))
      !allocate (elem%str_rate(el_count,gp, 4,1))
      !allocate (elem%rot_rate(el_count,gp, 4,1))
    else 
      allocate (elem%bl (el_count,gp,6,dim*nodxelem)) 
      allocate (elem%strain(el_count,gp, 6,1)) !!VECTORIZED 
      !allocate (elem%str_rate(el_count,gp, 6,1))
      !allocate (elem%rot_rate(el_count,gp, 6,1))
    end if 
    
    elem%gausspc(:) = gp
    
    bind_dom_type = 1
    axisymm_vol_weight = .False. !!!AREA WEIGHTED


  end subroutine

  !!! THIS IS TO AVERAGE DATA
  !!! ALLOCATE THE FOLLOWING
  ! ! real(fp_kind), dimension(:) :: elxnod !!!Elements shared by each node
  ! ! real(fp_kind), dimension(:,:) :: nodel !!! element node 
  subroutine SearchNodelem
    implicit none 
    integer :: e, n
    nod%elxnod = 0
    nod%nodel = 0 !Unnecessary
    do e=1,elem_count
      do n=1,nodxelem
        nod%elxnod(elem%elnod(e,n)) = nod%elxnod(elem%elnod(e,n)) + 1
        nod%nodel(elem%elnod(e,n),nod%elxnod(elem%elnod(e,n))) = e
      end do
    end do
  end subroutine
  
  !Average element constant data to elements
  subroutine AverageData(el_data, nod_data)
    !implicit none 
    integer :: e, n
    real(fp_kind),dimension(node_count), intent(out) :: nod_data  
    real(fp_kind),dimension(elem_count), intent(in) :: el_data  
    nod_data (:)= 0.0d0
    do n=1,node_count
      do e=1,nod%elxnod(n)
        !nod%elxnod(elem%elnod(e,n)) = nod%elxnod(elem%elnod(e,n)) + 1
        !nod%nodel(elem%elnod(e,n),nod%elxnod(elem%elnod(e,n))) = e
        !print *, "nod data n e", nod_data(n), n, e
        nod_data(n) = nod_data(n) + el_data(nod%nodel(n,e))
      end do
      nod_data(n) = nod_data(n) / nod%elxnod(n)
    end do
  end subroutine

  !Average element Nodal Data data to elements
  subroutine AssembleElNodData(el_data, nod_data)
    !implicit none 
    integer :: e, n
    real(fp_kind),dimension(node_count), intent(out) :: nod_data  
    real(fp_kind),dimension(elem_count,nodxelem), intent(in) :: el_data  
    nod_data (:)= 0.0d0
    do e=1,elem_count !Element per node
      do n=1,nodxelem
        !nod%elxnod(elem%elnod(e,n)) = nod%elxnod(elem%elnod(e,n)) + 1
        !nod%nodel(elem%elnod(e,n),nod%elxnod(elem%elnod(e,n))) = e
        !print *, "nod data n e", nod_data(n), n, e
        nod_data(elem%elnod(e,n)) = nod_data(elem%elnod(e,n)) + el_data(e,n)
      end do
    end do
  end subroutine
  
  subroutine MeshCSVreader()
    implicit none
    ! real(fp_kind),dimension(node_count), intent(out) :: nod_data  
    ! real(fp_kind),dimension(node_count), intent(out) :: ele_data  
    real :: x, y, z
    INTEGER :: m, n
    CHARACTER first*30
    logical :: readblock
    integer :: i 
    
    dim = 3

    OPEN(UNIT = 7, FILE = "mesh.csv")
    ! READ(7,*) x, y, z

    ! READ(7,*) m, n, first
    READ(7,*) first    
    if (first == '*Nodes') then
      print *, "*Nodes found."
      readblock = .true.
    end if
    readblock = .True.
    node_count = 0
    do while (readblock .eqv. .true.)
      READ(7,*) first
      if (first == '*Elements') then
        readblock = .false.
      else
        node_count = node_count + 1
      end if
    end do 
    readblock = .true.
    print * , "Node count ", node_count
    elem_count = 0
    do while (readblock .eqv. .true.)
      READ(7,*) first
      if (first == '*End') then
        readblock = .false.
      else
        elem_count = elem_count + 1
      end if
    end do 
    print * , "Element count ", elem_count
    rewind (7)
    READ(7,*) first    
    print *, "Allocating nodes... "
    call AllocateNodes(node_count)
    print *, "Reading nodes "
    do i=1, node_count
      READ(7,*) nod%x(i,1), nod%x(i,2), nod%x(i,3)
    end do
  
    print *, "Allocating elements... "
    nodxelem = 8
    call AllocateElements(elem_count, 1)
    READ(7,*) first    
    print *, first  
    !print *, elem%elnod(1,1), elem%elnod(1,2)
    do i=1, elem_count
      !print *, i
      
      READ(7,*) elem%elnod(i,1), elem%elnod(i,2), elem%elnod(i,3), elem%elnod(i,4), &
              & elem%elnod(i,5), elem%elnod(i,6), elem%elnod(i,7), elem%elnod(i,8)
      ! print *, elem%elnod(i,1), elem%elnod(i,2), elem%elnod(i,3), elem%elnod(i,4), &
              ! & elem%elnod(i,5), elem%elnod(i,6), elem%elnod(i,7), elem%elnod(i,8)
    end do
    
    close (7)
    
    
    call AllocateDomain()
    i = 1
    do while ( i <= node_count)
      nod%is_bcv(i,:) = .false.
      i = i + 1
    end do
  
    ! nod%m(:)   = Density * Lx * Ly * Lz / node_count
    ! nod%rho(:)   = Density
    
    !print *, "Particle mass ", nod%m(2)
    
    !nod%id(:) = tag
    
    fext_glob (:,:) = 0.0d0
    
    call SearchNodelem
    
    print *, "Done. "
  end subroutine MeshCSVreader
  
  !!!!! NODE DISTRIBUTION ARE LIKE IN FLANAGAN (1981), GOUDREAU, AND BENSON (1992)
  subroutine AddBoxLength(tag, V, Lx, Ly, Lz, r, Density,  h, redint)			
    integer, intent(in):: tag
    logical, intent(in) :: redint
    !real(fp_kind), intent(in), allocatable :: V
    real(fp_kind), dimension(1:3), intent(in)  :: V ! input
    real(fp_kind), intent(in):: r, Lx, Ly, Lz, Density, h  
    real(fp_kind), dimension (1:3) :: Xp
    integer :: i, j, k, p, ex, ey, ez, nnodz, gp
      
    integer, dimension(1:3) :: nel 
    
    nel(1) = nint(Lx/(2.0*r)) 
    nel(2) = nint(Ly/(2.0*r)) 
    if (Dim .eq. 2) then
      nel(3) = 0
      nodxelem = 4
    else
      nel(3) = nint(Lz/(2.0*r)) 
      nodxelem = 8
    end if
    
    Xp(3) = V(3) 
    

    write (*,*) "Creating Mesh ...", "Elements ", nel(1), ", ",nel(2)
    
    call AllocateNodes((nel(1) +1)* (nel(2)+1) * (nel(3)+1))
    print *, "Element count in XYZ: ", nel(:)
    write (*,*) "Box Node count ", node_count
    
    
    write (*,*) "xp ", Xp(:)    
    
    if (dim .eq. 2) then
    !write(*,*) "Box Particle Count is ", node_count
    p = 1
    !do while (Xp(3) <= (V(3)+Lz))
      j = 1;         Xp(2) = V(2)
      do while (j <= (nel(2) +1))
        i = 1
        Xp(1) = V(1)
        do while (i <= (nel(1) +1))
          nod%x(p,:) = Xp(:)
          print *,"node ",p , "X: ",Xp(:)
          p = p + 1
          Xp(1) = Xp(1) + 2 * r
          i = i +1
        end do
        Xp(2) = Xp(2) + 2 * r
        j = j +1
      end do 
      ! Xp(3) = Xp(3) + 2 * r
    !end do
    
    else 
      p = 1
      k = 1; Xp(3) = V(3)
      do while (k <= (nel(3) +1))
        j = 1;         Xp(2) = V(2)
        do while (j <= (nel(2) +1))
          i = 1
          Xp(1) = V(1)
          do while (i <= (nel(1) +1))
            nod%x(p,:) = Xp(:)
            print *,"node ",p , "X: ",Xp(:)
            p = p + 1
            Xp(1) = Xp(1) + 2 * r
            i = i +1
          end do
          Xp(2) = Xp(2) + 2 * r
          j = j +1
        end do 
        Xp(3) = Xp(3) + 2 * r
        k = k + 1
      end do    
    end if
    
    !! ALLOCATE ELEMENTS
    !! DIMENSION = 2
    gp = 1
    if (dim .eq. 2) then
      if (redint .eqv. .False.) then
        gp = 4
      end if 
      call AllocateElements(nel(1) * nel(2),gp) !!!!REDUCED INTEGRATION
    else 
      if (redint .eqv. .False.) then
        gp = 8
      end if 
      call AllocateElements(nel(1) * nel(2)*nel(3),gp) 
    end if
    
    if (dim .eq. 2) then
      ey = 0
      i = 1
      do while ( ey < nel(2))
          ex = 0
          do while (ex < nel(1)) 
              elem%elnod(i,:)=[(nel(1)+1)*ey + ex+1,(nel(1)+1)*ey + ex+2,(nel(1)+1)*(ey+1)+ex+2,(nel(1)+1)*(ey+1)+ex+1]         
              print *, "Element ", i , "Elnod", elem%elnod(i,:) 
              i=i+1
            ex = ex + 1
          end do
        ey = ey + 1
      end do  
    else 
      ez = 0
      i = 1
      nnodz = (nel(1)+1)*(nel(2)+1)
      print *, "Element Nodes at z ", nnodz
      do while ( ez < nel(3))
        ey = 0    
        do while ( ey < nel(2))
            ex = 0
            do while (ex < nel(1)) 
                !elem%elnod(i,:)=[(nel(1)+1)*(ey+1)+ex+2,(nel(1)+1)*(ey+1)+ex+1,(nel(1)+1)*ey + ex+1,(nel(1)+1)*ey + ex+2]       
                elem%elnod(i,:) = [ nnodz*ez + (nel(1)+1)*ey + ex+1,nnodz*ez + (nel(1)+1)*ey + ex+2, &
                                    nnodz*ez + (nel(1)+1)*(ey+1)+ex+2,nnodz*ez + (nel(1)+1)*(ey+1)+ex+1, &
                                    nnodz*(ez + 1) + (nel(1)+1)*ey + ex+1,nnodz*(ez + 1) + (nel(1)+1)*ey + ex+2, &
                                    nnodz*(ez + 1) + (nel(1)+1)*(ey+1)+ex+2,nnodz*(ez + 1)+ (nel(1)+1)*(ey+1)+ex+1]
                print *, "Element ", i , "Elnod", elem%elnod(i,:) 
                i=i+1
              ex = ex + 1
            end do
          ey = ey + 1
        end do 
        ez = ez + 1
      end do !el z
    end if !dim
    
    call AllocateDomain()
    i = 1
    do while ( i <= node_count)
      nod%is_bcv(i,:) = .false.
      i = i + 1
    end do
  
    ! nod%m(:)   = Density * Lx * Ly * Lz / node_count
    ! nod%rho(:)   = Density
    elem%rho_0(:,:) = Density
    !print *, "Particle mass ", nod%m(2)
    
    !nod%id(:) = tag
    
    fext_glob (:,:) = 0.0d0
    
    elem%e_length(:) = Lx !TODO: CHANGE!
    
    tot_mass = Density * Lx * Ly * Lz
    if (dim == 2) then !!!assuming plain strain
      tot_mass = tot_mass / Lz
    end if
    print *, "Box Total Mass: ", tot_mass
    
    call SearchNodelem
  end subroutine AddBoxLength
  
  
!THIS IS DONE IN ORDER TO USE THE SAME ASSEMBLY SUBROUTINE FOR 2D AND 3D
subroutine set_edof_from_elnod()
  integer :: e,dof,d,n

  do e=1,elem_count
    do n=1,nodxelem
      do d=1,dim

        dof = dim * (elem%elnod(e,n) - 1 ) + d
        !print *, "elem ", e, " dof loc", dim*(n-1)+d, "glob ", dof
        elem%dof (e,dim*(n-1)+d) = dof
      end do
    end do
    !print *, "elem ", e, " dof ", elem%dof(e,:)
  end do

end subroutine

End Module Domain 
