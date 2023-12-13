!https://gcc.gnu.org/onlinedocs/gfortran/C_005fF_005fPOINTER.html
! In case of array of character https://stackoverflow.com/questions/59247212/how-to-receive-a-string-from-a-c-function-called-by-fortran-by-iso-c-binding

 program main
   use iso_c_binding
   implicit none

   interface
      function func (a) bind (C, name="func")
         import
         integer(c_int):: func
         real(c_double), dimension(1:4), intent(in):: a
      end function func
      
      !!function ReadNastranTriMesh(fName, node, elcon)
      
      SUBROUTINE c_func(xval, s) BIND(C, name="c_func") ! ****
      import :: c_ptr

      TYPE(c_ptr) , intent(out):: xval ! ****
      integer,  intent(in) :: s! ****
      END SUBROUTINE c_func

      SUBROUTINE reader(fname, node, elnod) BIND(C, name="ReadNastranTriMesh") ! ****
      import :: c_ptr, c_char
      character(C_CHAR), intent(in)  :: fname(*)

      TYPE(c_ptr) , intent(out):: node, elnod ! ****

      END SUBROUTINE reader
      
   end interface

   real(c_double), dimension(1:4):: a = [ 2.3, 3.4, 4.5, 5.6 ]
   integer(c_int):: result
   
   
   !INTEGER, ALLOCATABLE :: X(:)
   INTEGER, pointer :: X(:)

   type(C_PTR) :: pX ! ****
   type(C_PTR) :: node,elnod ! ****
   
   CALL c_func(pX,100)
   CALL C_F_POINTER(pX, X, [100])
   
   print *, x(1), x(  2), x(3)
  
   result = func(a)
   
   call reader('tool_metal_cut_mm.nas', node, elnod)
   CALL C_F_POINTER(pX, X, [100])
   ! write *, "Node: ",node
   write (*,*) result
end program main