
module Domain
use NodeData
use ElementData
use ModPrecision, only : fp_kind
!use Neighbor 

implicit none 

integer :: Dim, node_count, elem_count, Nproc
type(Node)	::nod
type(Element)	::elem


integer :: nodxelem

real(fp_kind), dimension(3):: dommax, dommin

real(fp_kind)::mat_G, mat_E!TODO: change to material
real(fp_kind):: time

contains 

  subroutine DomInit(proc)
    integer, intent(in)::proc
    nproc = proc
    Dim = 3
    DomMax (:) = -100000000000.0;
    DomMin (:) = 100000000000.0;
  end subroutine DomInit
  
  subroutine AllocateNodes(pt_count)
    integer, intent(in):: pt_count
    
    node_count = pt_count
    !!!GENERAL 
    allocate (nod%x(node_count,dim))
    ! allocate (nod%rho(node_count))
    ! allocate (nod%drhodt(node_count))
    ! allocate (nod%h(node_count))
    ! allocate (nod%m(node_count))
    allocate (nod%id(node_count))
    
    !! THERMAL
    ! allocate (nod%cp_t(node_count))
    ! allocate (nod%t(node_count))
    ! allocate (nod%k_t(node_count))
    
    !MECHANICAL PROPERTIES
    !if (solver_type==1) then 
    
    allocate (nod%v(node_count,dim))
    allocate (nod%a(node_count,dim))
    allocate (nod%disp(node_count,dim))
    
    ! allocate (nod%sigma(node_count,3,3))
    ! allocate (nod%str_rate(node_count,3,3))
    ! allocate (nod%rot_rate(node_count,3,3))
    ! allocate (nod%shear_stress(node_count,3,3))
    ! allocate (nod%strain(node_count))
    ! allocate (nod%pressure(node_count))
    ! allocate (nod%cs(node_count))
    
    ! allocate (nod%sigma_eq(node_count))
    
    ! allocate (nod%rho_0(node_count))
    
    !end if  
  end subroutine
  
  subroutine AllocateElements(el_count)
    integer, intent(in):: el_count
    print *, "Creating ", el_count, " elements"
    elem_count = el_count
    allocate (elem%elnod(el_count,nodxelem))
    allocate (elem%x2(el_count,nodxelem,dim))
    allocate (elem%jacob(el_count,dim,dim))
    allocate (elem%detj(el_count))
    allocate (elem%sigma_eq(el_count))
    allocate (elem%dHxy(el_count,dim,nodxelem))
    
    if (Dim .eq. 2) then
      allocate (elem%bl (el_count,3,16))
      allocate (elem%bnl(el_count,4,16))
    else 
    end if 
  end subroutine

  subroutine AddBoxLength(tag, V, Lx, Ly, Lz, r, Density,  h)			
    integer, intent(in):: tag
    !real(fp_kind), intent(in), allocatable :: V
    real(fp_kind), dimension(1:3), intent(in)  :: V ! input
    real(fp_kind), intent(in):: r, Lx, Ly, Lz, Density, h  
    real(fp_kind), dimension (1:3) :: Xp
    integer :: i, j, p, ex, ey
      
    integer, dimension(1:3) :: nel 
    
    nel(1) = nint(Lx/(2.0*r)) 
    nel(2) = nint(Ly/(2.0*r)) 
    if (Dim .eq. 2) then
      nel(3) = 0
      nodxelem = 4
    else
      nel(2) = nint(Lz/(2.0*r)) 
      nodxelem = 8
    end if
    
    Xp(3) = V(3) 
    

    write (*,*) "Creating Mesh ...", "Elements ", nel(1), ", ",nel(2)
    
    call AllocateNodes((nel(1) +1)* (nel(2)+1) * (nel(3)+1))
    
    write (*,*) "Box Node count ", node_count
    
    
    write (*,*) "xp ", Xp(:)    
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
      Xp(3) = Xp(3) + 2 * r
    !end do
    
    !! ALLOCATE ELEMENTS
    !! DIMENSION = 2
    if (dim .eq. 2) then
      call AllocateElements(nel(1) * nel(2))
    else 
      call AllocateElements(nel(1) * nel(2)*nel(3))
    end if
    ey = 0
    i = 1
    do while ( ey < nel(2))
        ex = 0
        do while (ex < nel(1)) 
            elem%elnod(i,:)=[(nel(1)+1)*(ey+1)+ex+2,(nel(1)+1)*(ey+1)+ex+1,(nel(1)+1)*ey + ex+1,(nel(1)+1)*ey + ex+2]         
            print *, "Element ", i , "Elnod", elem%elnod(i,:) 
            i=i+1
          ex = ex + 1
        end do
      ey = ey + 1
    end do
  
    ! nod%m(:)   = Density * Lx * Ly * Lz / node_count
    ! nod%rho(:)   = Density
    ! nod%rho_0(:) = Density
    !print *, "Particle mass ", nod%m(2)
    
    !nod%id(:) = tag
    
    
  end subroutine AddBoxLength
End Module Domain 
