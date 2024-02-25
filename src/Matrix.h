/*  Copyright (c) 2013-2015 INGV, EDF, UniCT, JHU

    Istituto Nazionale di Geofisica e Vulcanologia, Sezione di Catania, Italy
    Électricité de France, Paris, France
    Università di Catania, Catania, Italy
    Johns Hopkins University, Baltimore (MD), USA

    This file is part of GPUSPH. Project founders:
        Alexis Hérault, Giuseppe Bilotta, Robert A. Dalrymple,
        Eugenio Rustico, Ciro Del Negro
    For a full list of authors and project partners, consult the logs
    and the project website <https://www.gpusph.org>

    GPUSPH is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GPUSPH is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GPUSPH.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _MATRIX_H_
#define _MATRIX_H_

#ifdef CUDA_BUILD

#include <cuda.h>
#define __spec __device__ __inline__
#else 
#define __spec __inline__
#endif

class Matrix {
public: 
  __spec Matrix(){}
  __spec Matrix(const int &row, const int &col);
  
  __spec double & getVal(int a, int b);
  __spec double & operator()(int a, int b);
  __spec Matrix  operator*(const double &f);
  __spec Matrix & Mul(const double &f);
  __spec void Set(const int &r, const int &c, const double &d);
	__spec void Print();
  __spec Matrix & Transpose();
  
  __spec ~Matrix(){/*cudaFree (m_data);*/
                       free(m_data);}
	
	__spec double calcDet();
  
  __spec void ToFlatSymPtr (double *flat, int initial);
  __spec void FromFlatSymPtr(double *flat, int initial);
  

	double *m_data;
  int m_row, m_col, m_dim;

};


__spec Matrix::Matrix(const int &row, const int &col) {
  m_row = row;
  m_col = col;
	if (m_row == m_col) m_dim = m_row;
  //cudaMalloc((void**)&m_data, row * col * sizeof(double)); //CRASHES
  m_data = (double*)malloc(row * col);
  //for (int i=0;i<row*col;i++) m_data[i] = 0.0;
}

//// OPERATION WITH MATRIX CREATION COULD PRODUCE MEM LEAKING
__spec Matrix MatMul(Matrix &A, Matrix &B){
  Matrix ret(A.m_row,B.m_col);
  for (int i = 0; i<A.m_row; i++)
    for (int j = 0; j<A.m_col; j++)
      for (int k = 0; k<A.m_col; k++)
        ret.m_data[i * A.m_row + j] += A.m_data[i * A.m_row + k] * B.m_data[k * B.m_row + j ];
  
  
  
  return ret;
}

__spec void MatMul(Matrix &A, Matrix &B, Matrix *ret){
  for (int i = 0; i<A.m_row; i++)
    for (int j = 0; j<A.m_col; j++)
      for (int k = 0; k<A.m_col; k++)
        ret->m_data[i * A.m_row + j] += A.m_data[i * A.m_row + k] * B.m_data[k * B.m_row + j ];

}

  __spec double & Matrix::getVal(int a, int b){
    return m_data[m_row*a+b];
  }
  
  __spec double & Matrix::operator()(int a, int b){
    return m_data[m_row*a+b];
  }

  __spec void Matrix::Set(const int &r, const int &c, const double &d){
    m_data[m_row*r+c] = d;
  }
  
  __spec Matrix & Matrix::Transpose(){
    for (int i=0;i<m_row*m_col;i++)
      for (int j=0;j<m_col*m_col;j++)
        this->Set(j,i,this->getVal(i,j) );
    return *this;
  }
  
	__spec Matrix operator*(const double &c, Matrix &A) {
	Matrix ret;
  for (int i=0;i<A.m_row*A.m_col;i++) ret.m_data[i] = A.m_data[i] * c;
	return ret;
}

	__spec Matrix Matrix::operator*(const double &f) {
  for (int i=0;i<m_row*m_col;i++) m_data[i] = f* m_data[i] ;
  return *this;
}

	__spec Matrix & Matrix::Mul(const double &f) {
  for (int i=0;i<m_row*m_col;i++) m_data[i] *= f;
  return *this;
}

	__spec void Matrix::Print() {
	//printf("%lf ",m_data[0]);
  for (int i=0;i<m_row;i++) {
		for (int j=0;j<m_col;j++) 
			printf("%lf ", getVal(i,j)) ;
		printf("\n");
	}

}



// subroutine M33INV (A, AINV, OK_FLAG)

  // IMPLICIT NONE

  // DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
  // DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
  // LOGICAL, INTENT(OUT) :: OK_FLAG

  // DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
  // DOUBLE PRECISION :: DETE
  // DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


  // ! DET =   A(1,1)*A(2,2)*A(3,3)  &
        // ! - A(1,1)*A(2,3)*A(3,2)  &
        // ! - A(1,2)*A(2,1)*A(3,3)  &
        // ! + A(1,2)*A(2,3)*A(3,1)  &
        // ! + A(1,3)*A(2,1)*A(3,2)  &
        // ! - A(1,3)*A(2,2)*A(3,1)
  
  // dete = det(A)

  // IF (ABS(DETE) .LE. EPS) THEN
     // AINV = 0.0D0
     // OK_FLAG = .FALSE.
     // RETURN
  // END IF

  // COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
  // COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
  // COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
  // COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
  // COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
  // COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
  // COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
  // COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
  // COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

  // AINV = TRANSPOSE(COFACTOR) / DETE

  // OK_FLAG = .TRUE.

// RETURN

// end subroutine M33INV

__spec double Matrix::calcDet (){
	double ret = 0.0;
  // real(fp_kind), dimension(dim,dim), intent (in) :: a 
  // real(fp_kind) :: det
  // if (dim .eq. 2) then
    // det = a(1,1)*a(2,2)-a(1,2)*a(2,1)
  // else 
  // DET =   A(1,1)*A(2,2)*A(3,3)  &
        // - A(1,1)*A(2,3)*A(3,2)  &
        // - A(1,2)*A(2,1)*A(3,3)  &
        // + A(1,2)*A(2,3)*A(3,1)  &
        // + A(1,3)*A(2,1)*A(3,2)  &
        // - A(1,3)*A(2,2)*A(3,1)  
  // end if
	if (m_dim == 2) {
		ret = getVal(0,0) * getVal(1,1) - getVal(0,1) * getVal(1,0);
	} else if (m_dim ==3) {
		ret =   getVal(0,0) * getVal(1,1) * getVal(2,2)
          - getVal(0,0) * getVal(1,2) * getVal(2,1)
					- getVal(0,1) * getVal(1,0) * getVal(2,2)
					+ getVal(0,1) * getVal(1,2) * getVal(2,0)
          + getVal(0,2) * getVal(1,0) * getVal(2,1)
					- getVal(0,2) * getVal(1,1) * getVal(2,0);


	}
	return ret;
}

__spec void InvMat(Matrix &A, Matrix *invA){
  if (A.m_dim ==2){
    
    
  } else if (A.m_dim == 3) {
    Matrix *cofactor = new Matrix(3,3);
    cofactor->Set(0,0, (A(1,1)*A(2,2)-A(1,2)*A(2,1)) );
    cofactor->Set(0,1,-(A(1,0)*A(2,2)-A(1,2)*A(2,0)) );
    cofactor->Set(0,2, (A(2,0)*A(2,1)-A(1,1)*A(2,0)) );
    cofactor->Set(1,0,-(A(0,1)*A(2,2)-A(0,2)*A(2,1)) );
    cofactor->Set(1,1, (A(0,0)*A(2,2)-A(0,2)*A(2,0)) );
    cofactor->Set(1,2,-(A(0,0)*A(2,1)-A(0,1)*A(2,0)) );
    cofactor->Set(2,0, (A(0,1)*A(1,2)-A(0,2)*A(1,1)) );
    cofactor->Set(2,1,-(A(0,0)*A(1,2)-A(0,2)*A(1,0)) );
    cofactor->Set(2,2, (A(0,0)*A(1,1)-A(0,1)*A(1,0)) );    
    
    cofactor->Transpose();
    *invA = cofactor->Mul(1.0/A.calcDet());
    
    delete cofactor;
  }
  
}

///// ONLY FOR 3X3 MATRICES
__spec void Matrix::ToFlatSymPtr(double *flat, int initial){
	flat [initial + 0] = getVal(0,0); 
  flat [initial + 1] = getVal(1,1); 
  flat [initial + 2] = getVal(2,2);
  flat [initial + 3] = getVal(0,1);
  flat [initial + 4] = getVal(1,2);
  flat [initial + 5] = getVal(0,2);
}
__spec void Matrix::FromFlatSymPtr(double *flat, int initial){
	Set(0,0,flat [initial + 0]);
  Set(0,1,flat [initial + 1]);
  Set(0,2,flat [initial + 2]);
  Set(1,1,flat [initial + 3]);
  Set(1,2,flat [initial + 4]);
  Set(2,2,flat [initial + 5]);
}



// function invmat (a)
  // real(fp_kind), dimension(dim,dim), intent (in) :: a 
  // real(fp_kind), dimension(dim,dim) :: invmat 
  // !if (dim .eq. 2) then
  // invmat(1,:) = 1.0d0/(det(a))*[ a(2,2),-a(1,2)]
  // invmat(2,:) = 1.0d0/(det(a))*[-a(2,1), a(1,1)]
  // !end if
// end function

// function adj (a)
  // real(fp_kind), dimension(dim,dim), intent (in) :: a 
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
// end function

#endif
