module mymath 
use ModPrecision, only : fp_kind

contains

function trace(a) result(j)
    real(fp_kind), intent (in) :: a(3,3) ! input
    real(fp_kind)              :: j ! output

    j = a(1,1)+a(2,2)+a(3,3)
end function

!!!!! ATENTION! THIS EIGENVALUE SOLVER IS FOR SYMMETRIC MATRICES!!!
!FROM https://github.com/awvwgk/diag3x3/blob/master/diag3x3.f90
pure subroutine eigval3x3(a, w)

   !> The symmetric input matrix
   real(fp_kind), intent(in) :: a(3, 3)

   !> Contains eigenvalues on exit
   real(fp_kind), intent(out) :: w(3)

   real(fp_kind) :: q, p, r

   r = a(1, 2) * a(1, 2) + a(1, 3) * a(1, 3) + a(2, 3) * a(2, 3)
   q = (a(1, 1) + a(2, 2) + a(3, 3)) / 3.0_fp_kind
   w(1) = a(1, 1) - q
   w(2) = a(2, 2) - q
   w(3) = a(3, 3) - q
   p = sqrt((w(1) * w(1) + w(2) * w(2) + w(3) * w(3) + 2*r) / 6.0_fp_kind)
   r = (w(1) * (w(2) * w(3) - a(2, 3) * a(2, 3)) &
      & - a(1, 2) * (a(1, 2) * w(3) - a(2, 3) * a(1, 3)) &
      & + a(1, 3) * (a(1, 2) * a(2, 3) - w(2) * a(1, 3))) / (p*p*p) * 0.5_fp_kind

   if (r <= -1.0_fp_kind) then
      r = 0.5_fp_kind * twothirdpi
   else if (r >= 1.0_fp_kind) then
      r = 0.0_fp_kind
   else
      r = acos(r) / 3.0_fp_kind
   end if

   w(3) = q + 2 * p * cos(r)
   w(1) = q + 2 * p * cos(r + twothirdpi)
   w(2) = 3 * q - w(1) - w(3)

end subroutine eigval3x3

!> Calculates eigenvector using an analytical method based on vector cross
!  products.
pure subroutine eigvec3x3(a, w, q)

   !> The symmetric input matrix, destroyed while solving
   real(fp_kind), intent(inout) :: a(3,3)

   !> Contains eigenvalues on exit
   real(fp_kind), intent(out) :: w(3)

   !> Contains eigenvectors on exit
   real(fp_kind), intent(out) :: q(3,3)

   !> Local variables
   real(fp_kind) :: norm, n1, n2, n3, precon
   integer :: i

   w(1) = max(abs(a(1, 1)), abs(a(1, 2)))
   w(2) = max(abs(a(1, 3)), abs(a(2, 2)))
   w(3) = max(abs(a(2, 3)), abs(a(3, 3)))
   precon = max(w(1), max(w(2), w(3)))

   ! null matrix
   if (precon < eps) then
      w(1) = 0.0_fp_kind
      w(2) = 0.0_fp_kind
      w(3) = 0.0_fp_kind
      q(1, 1) = 1.0_fp_kind
      q(2, 2) = 1.0_fp_kind
      q(3, 3) = 1.0_fp_kind
      q(1, 2) = 0.0_fp_kind
      q(1, 3) = 0.0_fp_kind
      q(2, 3) = 0.0_fp_kind
      q(2, 1) = 0.0_fp_kind
      q(3, 1) = 0.0_fp_kind
      q(3, 2) = 0.0_fp_kind
      return
   end if

   norm = 1.0_fp_kind / precon

   a(1, 1) = a(1, 1) * norm
   a(1, 2) = a(1, 2) * norm
   a(2, 2) = a(2, 2) * norm
   a(1, 3) = a(1, 3) * norm
   a(2, 3) = a(2, 3) * norm
   a(3, 3) = a(3, 3) * norm

   ! Calculate eigenvalues
   call eigval3x3(a, w)

   ! Compute first eigenvector
   a(1, 1) = a(1, 1) - w(1)
   a(2, 2) = a(2, 2) - w(1)
   a(3, 3) = a(3, 3) - w(1)

   q(1, 1) = a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2)
   q(2, 1) = a(1, 3) * a(1, 2) - a(1, 1) * a(2, 3)
   q(3, 1) = a(1, 1) * a(2, 2) - a(1, 2) * a(1, 2)
   q(1, 2) = a(1, 2) * a(3, 3) - a(1, 3) * a(2, 3)
   q(2, 2) = a(1, 3) * a(1, 3) - a(1, 1) * a(3, 3)
   q(3, 2) = a(1, 1) * a(2, 3) - a(1, 2) * a(1, 3)
   q(1, 3) = a(2, 2) * a(3, 3) - a(2, 3) * a(2, 3)
   q(2, 3) = a(2, 3) * a(1, 3) - a(1, 2) * a(3, 3)
   q(3, 3) = a(1, 2) * a(2, 3) - a(2, 2) * a(1, 3)
   n1 = q(1, 1) * q(1, 1) + q(2, 1) * q(2, 1) + q(3, 1) * q(3, 1)
   n2 = q(1, 2) * q(1, 2) + q(2, 2) * q(2, 2) + q(3, 2) * q(3, 2)
   n3 = q(1, 3) * q(1, 3) + q(2, 3) * q(2, 3) + q(3, 3) * q(3, 3)

   norm = n1
   i = 1
   if (n2 > norm) then
      i = 2
      norm = n1
   end if
   if (n3 > norm) then
      i = 3
   end if

   if (i == 1) then
      norm = sqrt(1.0_fp_kind / n1)
      q(1, 1) = q(1, 1) * norm
      q(2, 1) = q(2, 1) * norm
      q(3, 1) = q(3, 1) * norm
   else if (i == 2) then
      norm = sqrt(1.0_fp_kind / n2)
      q(1, 1) = q(1, 2) * norm
      q(2, 1) = q(2, 2) * norm
      q(3, 1) = q(3, 2) * norm
   else
      norm = sqrt(1.0_fp_kind / n3)
      q(1, 1) = q(1, 3) * norm
      q(2, 1) = q(2, 3) * norm
      q(3, 1) = q(3, 3) * norm
   end if

   ! Robustly compute a right-hand orthonormal set (ev1, u, v)
   if (abs(q(1, 1)) > abs(q(2, 1))) then
      norm = sqrt(1.0_fp_kind / (q(1, 1) * q(1, 1) + q(3, 1) * q(3, 1)))
      q(1, 2) = -q(3, 1) * norm
      q(2, 2) = 0.0_fp_kind
      q(3, 2) = +q(1, 1) * norm
   else
      norm = sqrt(1.0_fp_kind / (q(2, 1) * q(2, 1) + q(3, 1) * q(3, 1)))
      q(1, 2) = 0.0_fp_kind
      q(2, 2) = +q(3, 1) * norm
      q(3, 2) = -q(2, 1) * norm
   end if
   q(1, 3) = q(2, 1) * q(3, 2) - q(3, 1) * q(2, 2)
   q(2, 3) = q(3, 1) * q(1, 2) - q(1, 1) * q(3, 2)
   q(3, 3) = q(1, 1) * q(2, 2) - q(2, 1) * q(1, 2)

   ! Reset A
   a(1, 1) = a(1, 1) + w(1)
   a(2, 2) = a(2, 2) + w(1)
   a(3, 3) = a(3, 3) + w(1)

   ! A*U
   n1 = a(1, 1) * q(1, 2) + a(1, 2) * q(2, 2) + a(1, 3) * q(3, 2)
   n2 = a(1, 2) * q(1, 2) + a(2, 2) * q(2, 2) + a(2, 3) * q(3, 2)
   n3 = a(1, 3) * q(1, 2) + a(2, 3) * q(2, 2) + a(3, 3) * q(3, 2)

   ! A*V, note out of order computation
   a(3, 3) = a(1, 3) * q(1, 3) + a(2, 3) * q(2, 3) + a(3, 3) * q(3, 3)
   a(1, 3) = a(1, 1) * q(1, 3) + a(1, 2) * q(2, 3) + a(1, 3) * q(3, 3)
   a(2, 3) = a(1, 2) * q(1, 3) + a(2, 2) * q(2, 3) + a(2, 3) * q(3, 3)

   ! UT*(A*U) - l2*E
   n1 = q(1, 2) * n1 + q(2, 2) * n2 + q(3, 2) * n3 - w(2)
   ! UT*(A*V)
   n2 = q(1, 2) * a(1, 3) + q(2, 2) * a(2, 3) + q(3, 2) * a(3, 3)
   ! VT*(A*V) - l2*E
   n3 = q(1, 3) * a(1, 3) + q(2, 3) * a(2, 3) + q(3, 3) * a(3, 3) - w(2)

   if (abs(n1) >= abs(n3)) then
      norm = max(abs(n1), abs(n2))
      if (norm > eps) then
         if (abs(n1) >= abs(n2)) then
            n2 = n2 / n1
            n1 = sqrt(1.0_fp_kind / (1.0_fp_kind + n2 * n2))
            n2 = n2 * n1
         else
            n1 = n1 / n2
            n2 = sqrt(1.0_fp_kind / (1.0_fp_kind + n1 * n1))
            n1 = n1 * n2
         end if
         q(1, 2) = n2 * q(1, 2) - n1 * q(1, 3)
         q(2, 2) = n2 * q(2, 2) - n1 * q(2, 3)
         q(3, 2) = n2 * q(3, 2) - n1 * q(3, 3)
      end if
   else
      norm = max(abs(n3), abs(n2))
      if (norm > eps) then
         if (abs(n3) >= abs(n2)) then
            n2 = n2 / n3
            n3 = sqrt(1.0_fp_kind / (1.0_fp_kind + n2 * n2))
            n2 = n2 * n3
         else
            n3 = n3 / n2
            n2 = sqrt(1.0_fp_kind / (1.0_fp_kind + n3 * n3))
            n3 = n3 * n2
         end if
         q(1, 2) = n3 * q(1, 2) - n2 * q(1, 3)
         q(2, 2) = n3 * q(2, 2) - n2 * q(2, 3)
         q(3, 2) = n3 * q(3, 2) - n2 * q(3, 3)
      end if
   end if

   ! Calculate third eigenvector from cross product
   q(1, 3) = q(2, 1) * q(3, 2) - q(3, 1) * q(2, 2)
   q(2, 3) = q(3, 1) * q(1, 2) - q(1, 1) * q(3, 2)
   q(3, 3) = q(1, 1) * q(2, 2) - q(2, 1) * q(1, 2)

   w(1) = w(1) * precon
   w(2) = w(2) * precon
   w(3) = w(3) * precon

end subroutine eigvec3x3


! Concerning the internal storage of data, the SymTensor2 data is stored in vector of $6$ components named \textsf{\_data} using the following storage scheme:
! \begin{equation*}
! \T=\left[\begin{array}{ccc}
    ! T_{0} & T_{1} & T_{2}\\
    ! T_{1} & T_{3} & T_{4}\\
    ! T_{2} & T_{4} & T_{5}
    ! \end{array}\right]
! \end{equation*}
! //-----------------------------------------------------------------------------
! void Tensor2::buildFTF(double FTF[3][3]) const
! //-----------------------------------------------------------------------------
! {
  ! FTF[0][0] = dnlSquare(_data[0]) + dnlSquare(_data[3]) + dnlSquare(_data[6]);
  ! FTF[0][1] = _data[0] * _data[1] + _data[3] * _data[4] + _data[6] * _data[7];
  ! FTF[0][2] = _data[0] * _data[2] + _data[3] * _data[5] + _data[6] * _data[8];
  ! FTF[1][0] = FTF[0][1];
  ! FTF[1][1] = dnlSquare(_data[1]) + dnlSquare(_data[4]) + dnlSquare(_data[7]);
  ! FTF[1][2] = _data[1] * _data[2] + _data[4] * _data[5] + _data[7] * _data[8];
  ! FTF[2][0] = FTF[0][2];
  ! FTF[2][1] = FTF[1][2];
  ! FTF[2][2] = dnlSquare(_data[2]) + dnlSquare(_data[5]) + dnlSquare(_data[8]);
! }
! \begin{equation*}
! \T=\left[\begin{array}{ccc}
    ! T_{0} & T_{1} & T_{2}\\
    ! T_{3} & T_{4} & T_{5}\\
    ! T_{6} & T_{7} & T_{8}
    ! \end{array}\right]
! \end{equation*}
! @END
! */
!!!! FI IS NOT SYMMETRIC
subroutine buildFTF (tin,FTF)
  real(fp_kind), intent(in) :: tin(9)
  real(fp_kind), intent(out) :: FTF(3,3)
  
  FTF(1,1) = tin(1) * tin(1) + tin(4) * tin(4) + tin(7) * tin(7)
  FTF(1,2) = tin(1) * tin(2) + tin(4) * tin(5) + tin(7) * tin(8)
  FTF(1,3) = tin(1) * tin(3) + tin(4) * tin(6) + tin(7) * tin(9)
  FTF(2,1) = FTF(1,2)
  FTF(2,2) = tin(2) * tin(2) + tin(5) * tin(5) + tin(8) * tin(8)
  FTF(2,3) = tin(2) * tin(3) + tin(4) * tin(6) + tin(8) * tin(9)
  FTF(3,1) = FTF(1,3)
  FTF(3,2) = FTF(2,3)
  FTF(3,3) = tin(3)*tin(3) + tin(6)*tin(6) + tin(9)*tin(9)
end subroutine buildFTF

! //-----------------------------------------------------------------------------
! void Tensor2::polarExtract(double eigenVectors[3][3], double eigenValues[3], SymTensor2 &U, Tensor2 &R) const
! //-----------------------------------------------------------------------------
! {
  subroutine polarExtract ( tin, eigenVectors, eigenValues, U, R)
  real(fp_kind), intent(in) :: tin(9), eigenVectors(3,3),eigenValues(3)
  real(fp_kind), intent(out) :: U(6), R(3,3)
  ! double sq[3];

  ! // eigenVectors 1
  ! double U0[6];
  real(fp_kind) ::U0(6), U1(6), U2(6), Um1(6), sq(3), deter
  real(fp_kind) :: t1,t2,t3,t4,t5,t6
  ! U0[0] = dnlSquare(eigenVectors[0][0]);
  ! U0[1] = eigenVectors[0][0] * eigenVectors[1][0];
  ! U0[2] = eigenVectors[0][0] * eigenVectors[2][0];
  ! U0[3] = dnlSquare(eigenVectors[1][0]);
  ! U0[4] = eigenVectors[1][0] * eigenVectors[2][0];
  ! U0[5] = dnlSquare(eigenVectors[2][0]);
  U0(1) = eigenVectors(1,1)*eigenVectors(1,1)
  U0(2) = eigenVectors(1,1) * eigenVectors(2,1)
  U0(3) = eigenVectors(1,1) * eigenVectors(3,1)
  U0(4) = eigenVectors(2,1) * eigenVectors(2,1)
  U0(5) = eigenVectors(2,1) * eigenVectors(3,1)
  U0(6) = eigenVectors(3,1) * eigenVectors(3,1)
  ! // eigenVectors 2
  ! double U1[6];
  ! U1[0] = dnlSquare(eigenVectors[0][1]);
  ! U1[1] = eigenVectors[0][1] * eigenVectors[1][1];
  ! U1[2] = eigenVectors[0][1] * eigenVectors[2][1];
  ! U1[3] = dnlSquare(eigenVectors[1][1]);
  ! U1[4] = eigenVectors[1][1] * eigenVectors[2][1];
  ! U1[5] = dnlSquare(eigenVectors[2][1]);
  U1(1) = eigenVectors(1,2) * eigenVectors(1,2)
  U1(2) = eigenVectors(1,2) * eigenVectors(2,2)
  U1(3) = eigenVectors(1,2) * eigenVectors(3,2)
  U1(4) = eigenVectors(2,2) * eigenVectors(2,2)
  U1(5) = eigenVectors(2,2) * eigenVectors(3,2)
  U1(6) = eigenVectors(3,2) * eigenVectors(3,2)
  ! // eigenVectors 3
  ! double U2[6];
  ! U2[0] = dnlSquare(eigenVectors[0][2]);
  ! U2[1] = eigenVectors[0][2] * eigenVectors[1][2];
  ! U2[2] = eigenVectors[0][2] * eigenVectors[2][2];
  ! U2[3] = dnlSquare(eigenVectors[1][2]);
  ! U2[4] = eigenVectors[1][2] * eigenVectors[2][2];
  ! U2[5] = dnlSquare(eigenVectors[2][2]);
  U2(1) = eigenVectors(1,3) * eigenVectors(1,3)
  U2(2) = eigenVectors(1,3) * eigenVectors(2,3)
  U2(3) = eigenVectors(1,3) * eigenVectors(3,3)
  U2(4) = eigenVectors(2,3) * eigenVectors(2,3)
  U2(5) = eigenVectors(2,3) * eigenVectors(3,3)
  U2(6) = eigenVectors(3,3) * eigenVectors(3,3)
  
  sq(:) = eigenValues(:)
  ! sq[0] = sqrt(eigenValues[0]);
  ! sq[1] = sqrt(eigenValues[1]);
  ! sq[2] = sqrt(eigenValues[2]);
  ! U._data[0] = sq[0] * U0[0] + sq[1] * U1[0] + sq[2] * U2[0];
  ! U._data[1] = sq[0] * U0[1] + sq[1] * U1[1] + sq[2] * U2[1];
  ! U._data[2] = sq[0] * U0[2] + sq[1] * U1[2] + sq[2] * U2[2];
  ! U._data[3] = sq[0] * U0[3] + sq[1] * U1[3] + sq[2] * U2[3];
  ! U._data[4] = sq[0] * U0[4] + sq[1] * U1[4] + sq[2] * U2[4];
  ! U._data[5] = sq[0] * U0[5] + sq[1] * U1[5] + sq[2] * U2[5];
  U(1) = sq(1) * U0(1) + sq(2) * U1(1) + sq(3) * U2(1);
  U(2) = sq(1) * U0(2) + sq(2) * U1(2) + sq(3) * U2(2);
  U(3) = sq(1) * U0(3) + sq(2) * U1(3) + sq(3) * U2(3);
  U(4) = sq(1) * U0(4) + sq(2) * U1(4) + sq(3) * U2(4);
  U(5) = sq(1) * U0(5) + sq(2) * U1(5) + sq(3) * U2(5);
  U(6) = sq(1) * U0(6) + sq(2) * U1(6) + sq(3) * U2(6);

  ! double Um1[6];
  ! double t1 = U._data[3] * U._data[5];
  ! double t2 = U._data[2] * U._data[4];
  ! double t4 = U._data[4] * U._data[4];
  ! double t5 = U._data[1] * U._data[5];
  ! double t6 = U._data[2] * U._data[3];

  t1 = U(4) * U(6); t2 = U(3) * U(5);
  t4 = U(5) * U(5); t5 = U(2) * U(6);
  t6 = U(3) * U(4);

  ! double deter = U._data[0] * t1 + 2.0 * U._data[1] * t2 - U._data[0] * t4 - U._data[1] * t5 - U._data[2] * t6;
  deter = U(1) * t1 + 2.0 * U(2) * t2 - U(1) * t4 - U(2) * t5 - U(3) * t6;

  ! Um1[0] = (t1 - t4) / deter;
  ! Um1[1] = (t2 - t5) / deter;
  ! Um1[2] = (U._data[1] * U._data[4] - t6) / deter;
  ! Um1[3] = (U._data[0] * U._data[5] - U._data[2] * U._data[2]) / deter;
  ! Um1[4] = (U._data[2] * U._data[1] - U._data[0] * U._data[4]) / deter;
  ! Um1[5] = (U._data[0] * U._data[3] - U._data[1] * U._data[1]) / deter;


  Um1(1) = (t1 - t4) / deter;
  Um1(2) = (t2 - t5) / deter;
  Um1(3) = (U(2) * U(5) - t6) / deter;
  Um1(4) = (U(1) * U(6) - U(3) * U(3)) / deter;
  Um1(5) = (U(3) * U(2) - U(1) * U(5)) / deter;
  Um1(6) = (U(1) * U(4) - U(2) * U(2)) / deter;

  R(1,1) = tin(1) * Um1(1) + tin(2) * Um1(2) + tin(3) * Um1(3)
  R(1,2) = tin(1) * Um1(2) + tin(2) * Um1(4) + tin(3) * Um1(5)
  R(1,3) = tin(1) * Um1(3) + tin(2) * Um1(5) + tin(3) * Um1(6);
  
  R(2,1) = tin(4) * Um1(1) + tin(5) * Um1(2) + tin(6) * Um1(3);
  R(2,2) = tin(4) * Um1(2) + tin(5) * Um1(4) + tin(6) * Um1(5);
  R(2,3) = tin(4) * Um1(3) + tin(5) * Um1(5) + tin(6) * Um1(6);

  R(3,1) = tin(7) * Um1(1) + tin(8) * Um1(2) + tin(9) * Um1(3);
  R(3,2) = tin(7) * Um1(2) + tin(8) * Um1(4) + tin(9) * Um1(5);
  R(3,3) = tin(7) * Um1(3) + tin(8) * Um1(5) + tin(9) * Um1(6);
  
  ! ! R._data[0] = _data[0] * Um1[0] + _data[1] * Um1[1] + _data[2] * Um1[2];
  ! ! R._data[1] = _data[0] * Um1[1] + _data[1] * Um1[3] + _data[2] * Um1[4];
  ! ! R._data[2] = _data[0] * Um1[2] + _data[1] * Um1[4] + _data[2] * Um1[5];
  
  ! ! R._data[3] = _data[3] * Um1[0] + _data[4] * Um1[1] + _data[5] * Um1[2];
  ! ! R._data[4] = _data[3] * Um1[1] + _data[4] * Um1[3] + _data[5] * Um1[4];
  ! ! R._data[5] = _data[3] * Um1[2] + _data[4] * Um1[4] + _data[5] * Um1[5];
  
  ! ! R._data[6] = _data[6] * Um1[0] + _data[7] * Um1[1] + _data[8] * Um1[2];
  ! ! R._data[7] = _data[6] * Um1[1] + _data[7] * Um1[3] + _data[8] * Um1[4];
  ! ! R._data[8] = _data[6] * Um1[2] + _data[7] * Um1[4] + _data[8] * Um1[5];
end subroutine polarExtract

  ! // This method computes the polar decomposition of a second order tensor with computation of the \f$ ln[U] \f$ and \f$ R \f$ tensors as the returning arguments.
  ! // The logarithm of a symmetric tensor is givent by the following formulation:
  ! // \f[ \ln U =\sum _{i=1}^{3}\ln \lambda _{i}(u_{i}\otimes u_{i}) \f]
  ! // \param U Return second order tensor containing \f$ ln[U] \f$
  ! // \param R Return second order tensor containing \f$ R \f$

! //-----------------------------------------------------------------------------
! void SymTensor2::polarCuppen(SymTensor2 &U, Tensor2 &R) const
! //-----------------------------------------------------------------------------
subroutine polarCuppen(tin, U, R)
  real(fp_kind), intent(in) :: tin(3,3)
  real(fp_kind), intent(out) :: U(6), R(3,3)
  real (fp_kind) :: tin_plane(9), FTF(3,3), eigenVectors(3,3),eigenValues(3)  
  ! double FTF[3][3];
  ! double eigenVectors[3][3];
  ! double eigenValues[3];

  !Build the F(T).F symmetric matrix
  tin_plane(1) = tin(1,1);tin_plane(2) = tin(1,2);tin_plane(3) = tin(1,3);
  tin_plane(4) = tin(2,1);tin_plane(5) = tin(2,2);tin_plane(6) = tin(2,3);
  tin_plane(7) = tin(3,1);tin_plane(8) = tin(3,2);tin_plane(9) = tin(3,3);  
  call buildFTF(tin_plane,FTF)

  !Compute the eigenvalues and eigenvectors
  !dsyevd3(FTF, eigenVectors, eigenValues); !!!// Cuppen
  call eigvec3x3(FTF,eigenValues,eigenVectors) !!!FTF is destroyed!!!!!

  !Extract the tensors for U and R
  call polarExtract(tin_plane, eigenVectors, eigenValues, U, R);
  print *, "U ", U
  print *, "R ", R
end subroutine 


end module mymath

! ! /*
! ! @LABEL:SymTensor2::polarCuppen(SymTensor2 U, Tensor2 R)
! ! @SHORT:Polar decomposition of a second order tensor using the Cuppen’s Divide and Conquer algorithm.
! ! @RETURN:SymTensor2 and Tensor2
! ! @ARG:SymTensor2&U&Symmetric tensor $\log[\U]$
! ! @ARG:Tensor2&R&Rotation tensor $\R$
! ! This method computes the polar decomposition of a second order tensor $\F$ and returns the symmetric tensor $\R$ and the tensor $\U$ so that:
! ! \begin{equation*}
! ! \F = \R \cdot \U
! ! \end{equation*}
! ! where $\F$ is a second order tensor defined by the object itself.
! ! @END
! //-----------------------------------------------------------------------------
! void SymTensor2::polarCuppen(SymTensor2 &U, Tensor2 &R) const
! //-----------------------------------------------------------------------------
! {
  ! double FTF[3][3];
  ! double eigenVectors[3][3];
  ! double eigenValues[3];

  ! // Build the F(T).F symmetric matrix
  ! buildFTF(FTF);

  ! // Compute the eigenvalues and eigenvectors
  ! dsyevd3(FTF, eigenVectors, eigenValues); // Cuppen

  ! // Extract the tensors for U and R
  ! polarExtract(eigenVectors, eigenValues, U, R);
! }

! //-----------------------------------------------------------------------------
! void Tensor2::buildFTF(double FTF[3][3]) const
! //-----------------------------------------------------------------------------
! {
  ! FTF[0][0] = dnlSquare(_data[0]) + dnlSquare(_data[3]) + dnlSquare(_data[6]);
  ! FTF[0][1] = _data[0] * _data[1] + _data[3] * _data[4] + _data[6] * _data[7];
  ! FTF[0][2] = _data[0] * _data[2] + _data[3] * _data[5] + _data[6] * _data[8];
  ! FTF[1][0] = FTF[0][1];
  ! FTF[1][1] = dnlSquare(_data[1]) + dnlSquare(_data[4]) + dnlSquare(_data[7]);
  ! FTF[1][2] = _data[1] * _data[2] + _data[4] * _data[5] + _data[7] * _data[8];
  ! FTF[2][0] = FTF[0][2];
  ! FTF[2][1] = FTF[1][2];
  ! FTF[2][2] = dnlSquare(_data[2]) + dnlSquare(_data[5]) + dnlSquare(_data[8]);
! }


! //-----------------------------------------------------------------------------
! void Tensor2::polarExtractLnU(double eigenVectors[3][3], double eigenValues[3], SymTensor2 &U, Tensor2 &R) const
! //-----------------------------------------------------------------------------
! {
  ! double sq[3];

  ! // eigenVectors 1
  ! double U0[6];
  ! U0[0] = dnlSquare(eigenVectors[0][0]);
  ! U0[1] = eigenVectors[0][0] * eigenVectors[1][0];
  ! U0[2] = eigenVectors[0][0] * eigenVectors[2][0];
  ! U0[3] = dnlSquare(eigenVectors[1][0]);
  ! U0[4] = eigenVectors[1][0] * eigenVectors[2][0];
  ! U0[5] = dnlSquare(eigenVectors[2][0]);
  ! // eigenVectors 2
  ! double U1[6];
  ! U1[0] = dnlSquare(eigenVectors[0][1]);
  ! U1[1] = eigenVectors[0][1] * eigenVectors[1][1];
  ! U1[2] = eigenVectors[0][1] * eigenVectors[2][1];
  ! U1[3] = dnlSquare(eigenVectors[1][1]);
  ! U1[4] = eigenVectors[1][1] * eigenVectors[2][1];
  ! U1[5] = dnlSquare(eigenVectors[2][1]);
  ! // eigenVectors 3
  ! double U2[6];
  ! U2[0] = dnlSquare(eigenVectors[0][2]);
  ! U2[1] = eigenVectors[0][2] * eigenVectors[1][2];
  ! U2[2] = eigenVectors[0][2] * eigenVectors[2][2];
  ! U2[3] = dnlSquare(eigenVectors[1][2]);
  ! U2[4] = eigenVectors[1][2] * eigenVectors[2][2];
  ! U2[5] = dnlSquare(eigenVectors[2][2]);

  ! sq[0] = sqrt(eigenValues[0]);
  ! sq[1] = sqrt(eigenValues[1]);
  ! sq[2] = sqrt(eigenValues[2]);
  ! U._data[0] = sq[0] * U0[0] + sq[1] * U1[0] + sq[2] * U2[0];
  ! U._data[1] = sq[0] * U0[1] + sq[1] * U1[1] + sq[2] * U2[1];
  ! U._data[2] = sq[0] * U0[2] + sq[1] * U1[2] + sq[2] * U2[2];
  ! U._data[3] = sq[0] * U0[3] + sq[1] * U1[3] + sq[2] * U2[3];
  ! U._data[4] = sq[0] * U0[4] + sq[1] * U1[4] + sq[2] * U2[4];
  ! U._data[5] = sq[0] * U0[5] + sq[1] * U1[5] + sq[2] * U2[5];

  ! double Um1[6];
  ! double t1 = U._data[3] * U._data[5];
  ! double t2 = U._data[2] * U._data[4];
  ! double t4 = U._data[4] * U._data[4];
  ! double t5 = U._data[1] * U._data[5];
  ! double t6 = U._data[2] * U._data[3];

  ! double deter = U._data[0] * t1 + 2.0 * U._data[1] * t2 - U._data[0] * t4 - U._data[1] * t5 - U._data[2] * t6;

  ! Um1[0] = (t1 - t4) / deter;
  ! Um1[1] = (t2 - t5) / deter;
  ! Um1[2] = (U._data[1] * U._data[4] - t6) / deter;
  ! Um1[3] = (U._data[0] * U._data[5] - U._data[2] * U._data[2]) / deter;
  ! Um1[4] = (U._data[2] * U._data[1] - U._data[0] * U._data[4]) / deter;
  ! Um1[5] = (U._data[0] * U._data[3] - U._data[1] * U._data[1]) / deter;

  ! R._data[0] = _data[0] * Um1[0] + _data[1] * Um1[1] + _data[2] * Um1[2];
  ! R._data[1] = _data[0] * Um1[1] + _data[1] * Um1[3] + _data[2] * Um1[4];
  ! R._data[2] = _data[0] * Um1[2] + _data[1] * Um1[4] + _data[2] * Um1[5];
  ! R._data[3] = _data[3] * Um1[0] + _data[4] * Um1[1] + _data[5] * Um1[2];
  ! R._data[4] = _data[3] * Um1[1] + _data[4] * Um1[3] + _data[5] * Um1[4];
  ! R._data[5] = _data[3] * Um1[2] + _data[4] * Um1[4] + _data[5] * Um1[5];
  ! R._data[6] = _data[6] * Um1[0] + _data[7] * Um1[1] + _data[8] * Um1[2];
  ! R._data[7] = _data[6] * Um1[1] + _data[7] * Um1[3] + _data[8] * Um1[4];
  ! R._data[8] = _data[6] * Um1[2] + _data[7] * Um1[4] + _data[8] * Um1[5];

  ! sq[0] = log(eigenValues[0]) / 2;
  ! sq[1] = log(eigenValues[1]) / 2;
  ! sq[2] = log(eigenValues[2]) / 2;
  ! U._data[0] = sq[0] * U0[0] + sq[1] * U1[0] + sq[2] * U2[0];
  ! U._data[1] = sq[0] * U0[1] + sq[1] * U1[1] + sq[2] * U2[1];
  ! U._data[2] = sq[0] * U0[2] + sq[1] * U1[2] + sq[2] * U2[2];
  ! U._data[3] = sq[0] * U0[3] + sq[1] * U1[3] + sq[2] * U2[3];
  ! U._data[4] = sq[0] * U0[4] + sq[1] * U1[4] + sq[2] * U2[4];
  ! U._data[5] = sq[0] * U0[5] + sq[1] * U1[5] + sq[2] * U2[5];
! }