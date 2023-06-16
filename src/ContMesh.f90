!!!!!! https://fortranwiki.org/fortran/show/Object-oriented+progrsamming
module class_ContMesh
use ModPrecision, only : fp_kind

implicit none
public :: Mesh, circle_area, circle_print

  real :: pi = 3.1415926535897931d0 ! Class-wide private constant

!!!!!!!! CONTACT MESH
!!!!!! THIS TYPE SHOULD BE THE SAME FOR ALL, RIGID AND SOLID
type Mesh
  real :: radius
  real(fp_kind), dimension(:):: x,v !POSITION AND VEL
  integer :: elemcount  
end type Mesh

contains
 function circle_area(this) result(area)
    type(Mesh), intent(in) :: this
    real :: area
    area = pi * this%radius**2
  end function circle_area

  subroutine AxisPlaneMesh(this, AxisPlaneMesh(const int &axis, bool positaxisorent, const Vec3_t p1, const Vec3_t p2,  const int &dens)
    type(Mesh), intent(in) :: this
    real :: area
	int elemcount = dens * dens;
  
  if (dim == 2) then 
    elemcount = dens; 
	end if
  
	! double x1,x2,x3;
	! double l1,l2;
	! Vec3_t p = p2-p1;
	! int dir[3];
	! if 			(axis == 0 )	{dir[0] = 1; dir[1] = 2;}
	! else if (axis == 1 )	{dir[0] = 0; dir[1] = 2;}
	! else									{dir[0] = 0; dir[1] = 1;}
	
	! dir [2] = axis; //dir2 is which remains constant
	
	! x3 = p1(dir[2]);

	! x2=p1(dir[1]); 
	! double dl = p(dir[0])/dens;	//Could be allowed 2 diff densities
	! //cout <<"dens: "<<dens<<endl;
	! //Plane is in 0 and 1 dirs
	! int test =dens+1;
  
  ! if (dimension == 2) test = 1;
  
	! for (int j=0; j<test; j++) {
		! x1 = p1(dir[0]);
		! for (int i=0; i<dens+1; i++){
			! Vec3_t v;
			! v(dir[0])=x1;v(dir[1])=x2;v(dir[2])=x3;
			! //cout << "i,j" << i << ", " << j<<endl; 
			! //node.Push(new Vec3_t(x1,x2,x3));
			! node.Push(new Vec3_t(v(0),v(1),v(2)));
			! node_v.Push(new Vec3_t(0.,0.,0.));
			! //cout << "xyz: "<<x1 << ", "<<x2<<", "<<x3<<endl;
			! x1+=dl;
		! }
		! x2+=dl;
	! }
  ! cout << "Created "<<node.Size()<< " nodes "<<endl;

	! int n[4];
	! int el =0;
	! int i;
  ! cout << "Creating elements "<<endl;
  ! if (dimension == 3){
    ! for (size_t j = 0 ;j  < dens; j++ ) {
          ! // cout <<"j, dens" <<j<<", "<<dens<<endl;
          ! // cout <<"j<dens"<< (j  < dens)<<endl;
      ! for ( i = 0; i < dens; i++ ){
          ! // cout <<"i, dens" <<i<<", "<<dens<<endl;
          ! // cout <<"i <dens"<< (i  < dens)<<endl;
          ! n[0] = (dens + 1)* j + i; 		n[1] = n[0] + 1; 
          ! n[2] = (dens + 1)* (j+1) + i; n[3] = n[2] + 1;
        ! //cout <<" jj" << jj<<endl;
        ! int elcon[2][3];	// TODO: check x, y and z normals and node direction 
                          ! // For all plane orientations
        ! //If connectivity  is anticlockwise normal is outwards
        ! if (positaxisorent) {
          ! elcon[0][0] = n[0];elcon[0][1] = n[1];elcon[0][2] = n[2];
          ! elcon[1][0] = n[1];elcon[1][1] = n[3];elcon[1][2] = n[2];
        ! } else {
          ! elcon[0][0] = n[0];elcon[0][1] = n[2];elcon[0][2] = n[1];
          ! elcon[1][0] = n[1];elcon[1][1] = n[2];elcon[1][2] = n[3];				
        ! }
        ! //cout << "elnodes"<<endl;
        ! for ( int e= 0; e<2;e++) { // 2 triangles
          ! element.Push(new Element(elcon[e][0],elcon[e][1],elcon[e][2]));		
          ! //cout << "Element "<< el <<": ";
          ! // for (int en = 0 ; en<3; en++) cout << elcon[e][en]<<", ";
          ! // cout <<endl;
          
          ! Vec3_t v = ( *node[elcon[e][0]] + *node[elcon[e][1]] + *node[elcon[e][2]] ) / 3. ;
          ! element[el] -> centroid = v; 
          ! //cout << "Centroid" << element[el] -> centroid << endl;
          ! el++;
        ! }
      ! }// i for
      
    ! }
  ! } else {
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
  ! }    
    print *, 'Circle: r = ', this%radius, ' area = ', area
  end subroutine circle_print
  
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