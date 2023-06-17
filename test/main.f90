 program main
   use iso_c_binding
   implicit none

   interface
      function func (a) bind (C, name="func")
         import
         integer(c_int):: func
         real(c_double), dimension(1:4), intent(in):: a
      end function func
   end interface

   real(c_double), dimension(1:4):: a = [ 2.3, 3.4, 4.5, 5.6 ]
   integer(c_int):: result
   result = func(a)
   write (*,*) result
end program main