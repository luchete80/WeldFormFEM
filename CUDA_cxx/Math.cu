#include "Math.cuh"
__device__  void adj(a) {
  real(fp_kind), dimension(dim,dim), intent (in) :: a 
  // real(fp_kind), dimension(dim,dim) :: cofactor,adj
  
  // if (dim .eq. 2) then
    // adj(1,:) = [ a(2,2),-a(1,2)]
    // adj(2,:) = [-a(2,1), a(1,1)]
  // else
    // cofactor(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
    // cofactor(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    // cofactor(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    // cofactor(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
    // cofactor(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
    // cofactor(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
    // cofactor(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
    // cofactor(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
    // cofactor(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))
    
    // adj = TRANSPOSE(COFACTOR)
  // end if
}
