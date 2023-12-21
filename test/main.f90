!https://gcc.gnu.org/onlinedocs/gfortran/C_005fF_005fPOINTER.html
! In case of array of character https://stackoverflow.com/questions/59247212/how-to-receive-a-string-from-a-c-function-called-by-fortran-by-iso-c-binding

! https://stackoverflow.com/questions/16385372/allocating-memory-in-c-for-a-fortran-array
! https://www.ibm.com/docs/es/openxl-fortran-aix/17.1.0?topic=calls-mixing-fortran-c  

    ! FORTRAN TO C STING
!https://community.intel.com/t5/Intel-Fortran-Compiler/passing-string-from-Fortran-to-C/td-p/1138766
   
 program main
   use iso_c_binding
   !use String_mod
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

      SUBROUTINE reader(fname, node, elnod ,nodecount) BIND(C, name="ReadNastranTriMesh") ! ****
      import :: c_ptr, c_char
      character(C_CHAR), intent(in)  :: fname(*)

      TYPE(c_ptr) , intent(out):: node, elnod ! ****
      integer :: nodecount

      END SUBROUTINE reader

      SUBROUTINE lsdynareadnodes(fname, node, nodecount) BIND(C, name="readNodes") ! ****
      import :: c_ptr, c_char
      character(C_CHAR), intent(in)  :: fname(*)

      TYPE(c_ptr) , intent(out):: node ! ****
      integer :: nodecount

      END SUBROUTINE lsdynareadnodes
    end interface

   real(c_double), dimension(1:4):: a = [ 2.3, 3.4, 4.5, 5.6 ]
   integer(c_int):: result
   
   !!! ALLOCATION
   !INTEGER, ALLOCATABLE :: X(:)
   INTEGER, pointer:: X(:)
   real(8), pointer :: Node(:)
   integer :: length, nodecount, i

   type(C_PTR) :: pX ! ****
   type(C_PTR) :: pnode,elnod ! ****

  integer :: arg_no
  character(1024) :: cmd_arg
  do arg_no = 0, 4
  call getarg(arg_no, cmd_arg)
  print *, len_trim(cmd_arg), trim(cmd_arg)
   end do
   
   CALL c_func(pX,100)
   CALL C_F_POINTER(pX, X, [100])
   print *, "X values "
   print *, x(1), x(  2), x(3)
   print *, "size of X: ", size(x)
   ! print *, "size of pX: ", size(px)
   result = func(a)
   
   ! !! THIS WORKS !!
   ! call reader("tool_metal_cut_mm.nas", pnode, elnod, nodecount)
   ! length = 3*nodecount
   ! CALL C_F_POINTER(pNode, node, [length])
  
   ! do i = 1, length
    ! print *, node(i)
   ! end do
   ! !!!
   
   
   call lsdynareadnodes('sphere-plate.k', pnode, nodecount)
   ! length = 3*nodecount
   ! CALL C_F_POINTER(pNode, node, [length])
   
end program main