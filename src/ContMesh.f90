!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WeldFormFEM - A C++/Fortran library to simulate Mechanical Solids                !
!               using Explicit Finite Element Method                               !
! Copyright (C) 2023 Luciano Buglioni                                              !
!                                                                                  !
!                                                                                  !
! This is free software; you can redistribute it andor modify it under the         !
! terms of the GNU General Public License as published by the Free Software        !
! Foundation; either version 3 of the License, or (at your option) any later       !
! version.                                                                         !
!                                                                                  !
! This program is distributed in the hope that it will be useful, but WITHOUT ANY  !
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  !
! PARTICULAR PURPOSE. See the GNU General Public License for more details.         !
!                                                                                  !
! You should have received a copy of the GNU General Public License along with     !
! PersianSPH; if not, see <http:www.gnu.orglicenses>                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! https://fortranwiki.org/fortran/show/Object-oriented+progrsamming
!!!! AIM TO GIVE LIKE STRUCT OF ARRAY (TO CONVERT IT EASILY TO CUDA)

module class_ContMesh
use ModPrecision, only : fp_kind
use NodeData
use ElementData
use Domain

implicit none
public :: Mesh, circle_area, circle_print

  real :: pi = 3.1415926535897931d0 ! Class-wide private constant

!!!!!!!! CONTACT MESH
!!!!!! THIS TYPE SHOULD BE THE SAME FOR ALL, RIGID AND SOLID
type Mesh !!! THIS SHOULD BE RIGID MESH THING
  real :: radius
  real(fp_kind), dimension(:,:), Allocatable:: x,v !POSITION AND VEL
  integer, Dimension(:,:), Allocatable :: elnod
  integer :: elem_count, node_count
  
  !!!type(Node)	::nod
  !!!type(Element)	::elem

end type Mesh

contains
 function circle_area(this) result(area)
    type(Mesh), intent(in) :: this
    real :: area
    area = pi * this%radius**2
  end function circle_area


  !AxisPlaneMesh(const int &axis, bool positaxisorent, const Vec3_t p1, const Vec3_t p2,  const int &dens
  !!! TODO: PASS DIM AS AN ARGUMENT
  subroutine AxisPlaneMesh(this, axis, positaxisorent, p1, p2,  dens)
    implicit none
    type(Mesh), intent(out) :: this
    integer, intent(in)::axis
    logical, intent(in):: positaxisorent
    real(fp_kind), intent(in),dimension(3) :: p1,p2
    integer, intent(in) :: dens
    real :: area
    
    integer, dimension(4) :: n
    integer, dimension(3) :: dir 
    integer :: i, j, e, test, elcon(2,3), k
    real(fp_kind) ::x1,x2,x3,dl,v(3)
    real(fp_kind), dimension(3) :: p
    
    this%elem_count = 2 * dens * dens
    this%node_count = (dens +1) * (dens + 1)
    
    print *, "Contact  Mesh with ", this%node_count, " nodes and ", this%elem_count, " elements was created."
    allocate (this%x(this%node_count,dim))
    allocate (this%v(this%node_count,dim))
    allocate (this%elnod(this%elem_count,dim))
        
    if (dim .eq. 2) then 
      elem_count = dens 
    end if
  
	! double x1,x2,x3;
	! double l1,l2;
	p = p2-p1;
  print *, "Point Length: ", p
	! int dir[3];
  select case (axis )	
    case (1) 
    dir(1) = 2; dir(2) = 3; 
    case (2) 
      dir(1) = 1; dir(2) = 3;
    case(3) 
      dir(1) = 1; dir(2) = 2;
	end select
	

  dir (3) = axis; !dir2 is which remains constant
  print *, "Direction: ", dir
	
	x3 = p1(dir(3));
	x2 = p1(dir(2)); 
  dl = p(dir(1))/dens;	!!!!Could be allowed 2 diff densities
  print *, "dl: ", dl
	! //cout <<"dens: "<<dens<<endl;
	! //Plane is in 0 and 1 dirs
	test =dens + 1 
  if (dim .eq. 2) then
    test = 1
  end if 
  
  !!! CREATING NODES
  k=1
  do j=1, test !for (int j=0; j<test; j++) {
    x1 = p1(dir(1))
    do i=1,dens+1
      v(dir(1))=x1;v(dir(2))=x2;v(dir(3))=x3;
        this%x(k,:)=v
      print *,  "xyz: ", x1 ,", ",x2,", ",x3
      x1 = x1 + dl;
      k=k+1
    end do
    x2 = x2 + dl;
  end do
  ! cout << "Created "<<node.Size()<< " nodes "<<endl;

	! int n[4];
	! int el =0;
	! int i;
  ! cout << "Creating elements "<<endl;
  print *, "Creating contact mesh elements.. "

  if (dim .eq. 3) then
    e = 1
    do j=1,dens
      do i=1,dens    
        n(1) = (dens + 1)* (j-1) + i; 		n(2) = n(1) + 1;   !!!! 3--4
        n(3) = (dens + 1)* j + i;         n(4) = n(3) + 1;   !!!! 1--2
        !print * , "Node ", n
        if (positaxisorent .eqv. .true.) then
          elcon(1,1) = n(1);elcon(1,2) = n(2);elcon(1,3) = n(3); !! 1, 2, 3
          elcon(2,1) = n(2);elcon(2,2) = n(4);elcon(2,3) = n(3);!!! 2,4,3
        else 
          elcon(1,1) = n(1);elcon(1,2) = n(3);elcon(1,3) = n(2);
          elcon(2,1) = n(2);elcon(2,2) = n(3);elcon(2,3) = n(4);          
        endif
        this%elnod(e,:)   = elcon(1,:)
        this%elnod(e+1,:) = elcon(2,:)
        !print *, "Element ", e,this%elnod(e,:) 
        !print *, "Element ", e+1,this%elnod(e+1,:)
        e = e + 2
      end do
    end do 

  else ! dim = 2
    ! for ( i = 0; i < dens; i++ ){
          ! n[0] = i; 		n[1] = n[0] + 1; 
        ! //cout <<" jj" << jj<<endl;
        ! int elcon[2];	// TODO: check x, y and z normals and node direction 
                          ! // For all plane orientations
        ! if (positaxisorent) {
          ! elcon[0] = n[0];elcon[1] = n[1];
        ! } else {
          ! elcon[0] = n[1];elcon[1] = n[0];		
        ! }

        ! element.Push(new Element(elcon[0],elcon[1],0));		
        ! //cout << "Element "<< el <<": ";
        ! // for (int en = 0 ; en<3; en++) cout << elcon[e][en]<<", ";
        ! // cout <<endl;
        
        ! Vec3_t v = ( *node[elcon[0]] + *node[elcon[1]] ) / 2. ;
        ! element[el] -> centroid = v; 
        ! //cout << "Centroid" << element[el] -> centroid << endl;
        ! el++;                
    ! }  
  end if !DIMENSION  
    print *, 'Circle: r = ', this%radius, ' area = ', area
  end subroutine AxisPlaneMesh
  
  subroutine circle_print(this)
    type(Mesh), intent(in) :: this
    real :: area
    area = circle_area(this)  ! Call the circle_area function
    print *, 'Circle: r = ', this%radius, ' area = ', area
  end subroutine circle_print
  
end module

!!!!! If wanted to implement
! program circle_test
  ! use class_Circle
  ! implicit none

  ! type(Circle) :: c     ! Declare a variable of type Circle.
  ! c = Circle(1.5)       ! Use the implicit constructor, radius = 1.5.
  ! call circle_print(c)  ! Call a class subroutine
! end program circle_test