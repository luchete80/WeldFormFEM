module mymath 
use ModPrecision, only : fp_kind

contains

function trace(a) result(j)
  real(fp_kind), intent (in) :: a(3,3) ! input
  real(fp_kind)              :: j ! output
  j = a(1,1)+a(2,2)+a(3,3)
end function

! function trace(a) result(j)
  ! real(fp_kind), intent (in) :: a(dim,dim) ! input
  ! real(fp_kind)              :: j ! output
  ! do i =1, dim
    ! j = j + a(dim,dim)
  ! end do 
! end function

function deviator(a) result(j)
  real(fp_kind), intent (in) :: a(3,3) ! input
  real(fp_kind)              :: j(3,3)! output

  real(fp_kind) :: ident(3,3)
  ident = 0.0d0
  ident (1,1) = 1.0d0; ident (2,2) = 1.0d0; ident (3,3) = 1.0d0
  
  j = a - 1.0d0/3.0d0 * trace(a) * ident
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
  FTF(2,3) = tin(2) * tin(3) + tin(5) * tin(6) + tin(8) * tin(9)
  FTF(3,1) = FTF(1,3)
  FTF(3,2) = FTF(2,3)
  FTF(3,3) = tin(3)*tin(3) + tin(6)*tin(6) + tin(9)*tin(9)
end subroutine buildFTF

! ! //-----------------------------------------------------------------------------
! ! void Tensor2::polarExtractLnU(double eigenVectors[3][3], double eigenValues[3], SymTensor2 &U, Tensor2 &R) const
! ! //-----------------------------------------------------------------------------
  subroutine polarExtractLnU ( tin, eigenVectors, eigenValues, U, R)
  real(fp_kind), intent(in) :: tin(9), eigenVectors(3,3),eigenValues(3)
  real(fp_kind), intent(out) :: U(6), R(3,3)
  !real(fp_kind):: U(6)
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
  U0(1) = eigenVectors(1,1) * eigenVectors(1,1)
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
  
  sq(:) = sqrt(eigenValues(:))
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

  ! U(1,1) = sq(1) * U0(1) + sq(2) * U1(1) + sq(3) * U2(1);
  ! U(1,2) = sq(1) * U0(2) + sq(2) * U1(2) + sq(3) * U2(2);
  ! U(1,3) = sq(1) * U0(3) + sq(2) * U1(3) + sq(3) * U2(3);
  ! U(2,2) = sq(1) * U0(4) + sq(2) * U1(4) + sq(3) * U2(4);
  ! U(2,3) = sq(1) * U0(5) + sq(2) * U1(5) + sq(3) * U2(5);
  ! U(3,3) = sq(1) * U0(6) + sq(2) * U1(6) + sq(3) * U2(6);
  ! U(2,1) = U(1,2);U(3,1) = U(1,3);U(3,2) = U(2,3);
  ! double Um1[6];
  ! double t1 = U._data[3] * U._data[5];
  ! double t2 = U._data[2] * U._data[4];
  ! double t4 = U._data[4] * U._data[4];
  ! double t5 = U._data[1] * U._data[5];
  ! double t6 = U._data[2] * U._data[3];

  t1 = U(4) * U(6); t2 = U(3) * U(5);
  t4 = U(5) * U(5); t5 = U(2) * U(6);
  t6 = U(3) * U(4);
  print *, "t1 a 6 ", t1,t2,t4,t5,t6
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
  ! Um(1,1) = U(1);    Um(1,2) = U(2); Um(1,2) = U(3);
  ! Um(2,1) = Um(1,2); Um(2,2) = U(4); Um(2,3) = U(5);
  ! Um(3,3) = U(6);
  
  sq(:) = log(eigenValues(:))/2.0d0
  
  U(1) = sq(1) * U0(1) + sq(2) * U1(1) + sq(3) * U2(1);
  U(2) = sq(1) * U0(2) + sq(2) * U1(2) + sq(3) * U2(2);
  U(3) = sq(1) * U0(3) + sq(2) * U1(3) + sq(3) * U2(3);
  U(4) = sq(1) * U0(4) + sq(2) * U1(4) + sq(3) * U2(4);
  U(5) = sq(1) * U0(5) + sq(2) * U1(5) + sq(3) * U2(5);
  U(6) = sq(1) * U0(6) + sq(2) * U1(6) + sq(3) * U2(6);  
end subroutine polarExtractLnU

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
  real (fp_kind) :: xx(3,3)
  real (fp_kind) :: abserr
  integer :: i
  
  abserr = 1.0e-15
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
  !call eigvec3x3(FTF,eigenValues,eigenVectors) !!!FTF is destroyed!!!!!
  !print *, "eigenvalues" , eigenValues
  
  call Jacobi(FTF,xx,abserr,3)
  !print *, "eigenvalues" , FTF(1,1),FTF(2,2),FTF(3,3)
  eigenValues(1) = FTF(1,1);eigenValues(2) = FTF(2,2);eigenValues(3) = FTF(3,3);
  do i=1,3
    eigenVectors(:,i) = xx(:,i)
  end do

  !Extract the tensors for U and R
  call polarExtractLnU(tin_plane, eigenVectors, eigenValues, U, R);

end subroutine 

!!!https://github.com/minar09/parallel-computing/blob/master/OMP_Fortran.f90
subroutine Jacobi(a,x,abserr,n)
!===========================================================
! Evaluate eigenvalues and eigenvectors
! of a real symmetric matrix a(n,n): a*x = lambda*x 
! method: Jacoby method for symmetric matrices 
! Alex G. (December 2009)
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - number of equations
! abserr - abs tolerance [sum of (off-diagonal elements)^2]
! output ...
! a(i,i) - eigenvalues
! x(i,j) - eigenvectors
! comments ...
!===========================================================
implicit none
integer i, j, k, n
real (fp_kind) :: a(n,n),x(n,n)
real (fp_kind) :: abserr, b2, bar
real (fp_kind) :: beta, coeff, c, s, cs, sc

! initialize x(i,j)=0, x(i,i)=1
! *** the array operation x=0.0 is specific for Fortran 90/95
x = 0.0
do i=1,n
  x(i,i) = 1.0
end do

! find the sum of all off-diagonal elements (squared)
b2 = 0.0
do i=1,n
  do j=1,n
    if (i.ne.j) b2 = b2 + a(i,j)**2
  end do
end do

if (b2 <= abserr) return

! average for off-diagonal elements /2
bar = 0.5*b2/float(n*n)

do while (b2.gt.abserr)
  do i=1,n-1
    do j=i+1,n
      if (a(j,i)**2 <= bar) cycle  ! do not touch small elements
      b2 = b2 - 2.0*a(j,i)**2
      bar = 0.5*b2/float(n*n)
! calculate coefficient c and s for Givens matrix
      beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
      coeff = 0.5*beta/sqrt(1.0+beta**2)
      s = sqrt(max(0.5+coeff,0.0))
      c = sqrt(max(0.5-coeff,0.0))
! recalculate rows i and j
      do k=1,n
        cs =  c*a(i,k)+s*a(j,k)
        sc = -s*a(i,k)+c*a(j,k)
        a(i,k) = cs
        a(j,k) = sc
      end do
! new matrix a_{k+1} from a_{k}, and eigenvectors 
      do k=1,n
        cs =  c*a(k,i)+s*a(k,j)
        sc = -s*a(k,i)+c*a(k,j)
        a(k,i) = cs
        a(k,j) = sc
        cs =  c*x(k,i)+s*x(k,j)
        sc = -s*x(k,i)+c*x(k,j)
        x(k,i) = cs
        x(k,j) = sc
      end do
    end do
  end do
end do
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