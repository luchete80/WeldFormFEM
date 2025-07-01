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

#include <stdarg.h>
#include "defs.h"

#ifdef CUDA_BUILD

#include <cuda.h>
#define __spec __device__ inline
#else 
#define __spec inline
#include <stdio.h>
#endif


//template <int ROWS, int COLS>
class Matrix {
public: 
  __spec Matrix() : m_data(nullptr), m_row(0), m_col(0), m_dim(0) {}
  __spec Matrix(const int row, const int col) {
  m_row = row;
  m_col = col;
	if (m_row == m_col) m_dim = m_row;
  //cudaMalloc((void**)&m_data, row * col * sizeof(double)); //CRASHES
  //m_data = (double*)malloc(row * col * sizeof (double));
  malloc_dev_t(m_data,double,m_col*m_row);
        if (m_data == nullptr) {
            printf("Memory allocation failed!\n");
            //exit(EXIT_FAILURE);  // Handle allocation failure
        }
        
  for (int i=0;i<row*col;i++) m_data[i] = 0.0;
}
  //~ void setIdentity(){
    //~ for (int r=0;r<m_row;r++){
        //~ if (r<m_col)Set(r,r,1.0);
    //~ }
  
  inline double & getVal(int a, int b);
  double& at(int i, int j)const;          // non-const accessor for assignment
  __spec double & operator()(int a, int b);
  __spec Matrix  operator*(const double &f);
  //__spec Matrix  operator=(const Matrix &m);

// Copy constructor
    __spec Matrix(const Matrix& other) : m_row(other.m_row), m_col(other.m_col) {
        //m_data = (double *)malloc(m_row * m_col * sizeof(double));
        malloc_dev_t(m_data,double,m_col*m_row);
        if (m_data == nullptr) {
            printf("Memory allocation failed!n");
            //exit(EXIT_FAILURE);
        }
        memcpy(m_data, other.m_data, m_row * m_col * sizeof(double));
    }

    // Assignment operator
    __spec Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            // Free existing memory
            if (m_data != nullptr) {
                free_dev_t(m_data);
                m_data = nullptr;
            }
            
            // Allocate new memory
            m_row = other.m_row;
            m_col = other.m_col;
            //m_data = (double *)malloc(m_row * m_col * sizeof(double));
            malloc_dev_t(m_data,double,m_col*m_row);
            if (m_data == nullptr) {
                printf("Memory allocation failed!n");
                //exit(EXIT_FAILURE);
            }

            // Copy data
            //memcpy(m_data, other.m_data, m_row * m_col * sizeof(double));
            memcpy_dev_t(m_data, other.m_data, m_row * m_col * sizeof(double));
        }
        return *this;
    }

  // Move Constructor
  __spec Matrix(Matrix&& other) noexcept
      : m_data(other.m_data), m_row(other.m_row), m_col(other.m_col), m_dim(other.m_dim) {
      other.m_data = nullptr;
      other.m_row = other.m_col = other.m_dim = 0;
  }

  // Move Assignment Operator
  __spec Matrix& operator=(Matrix&& other) noexcept {
      if (this != &other) {
          free_dev_t(m_data);  // Free current resources
          m_data = other.m_data;
          m_row = other.m_row;
          m_col = other.m_col;
          m_dim = other.m_dim;

          other.m_data = nullptr;
          other.m_row = other.m_col = other.m_dim = 0;
      }
      return *this;
  }

  __spec Matrix & Mul(const double &f);
  __spec void Set(const int r, const int c, const double d);
  __spec void SetZero();
	__spec void Print();
  __spec Matrix & Transpose();
  __spec Matrix getTranspose();
  __spec Matrix Inv();
  void Free(){free_dev_t(m_data);m_data=nullptr;}
  __spec ~Matrix(){/*cudaFree (m_data);*/
      if (m_data != nullptr) {
      //printf("deleting\n");
      free_dev_t(m_data); 
    }
  }
	
	__spec double calcDet();
  
  __spec void ToFlatSymPtr (double *flat, int initial);
  __spec void FromFlatSymPtr(double *flat, int initial);
  

	double *m_data;
  int m_row, m_col, m_dim;


__spec Matrix& operator+=(Matrix& other) {
    // Check dimensions match
    if (m_row != other.m_row || m_col != other.m_col) {
        printf("Matrix dimension mismatch in operator+=\n");
        return *this; // or handle error properly
    }
    for (int i = 0; i < m_row * m_col; ++i) {
        m_data[i] += other.m_data[i];
    }
    return *this;
}

__spec Matrix& operator+=(Matrix other) {
    // Check dimensions match
    if (m_row != other.m_row || m_col != other.m_col) {
        printf("Matrix dimension mismatch in operator+=\n");
        return *this; // or handle error properly
    }
    for (int i = 0; i < m_row * m_col; ++i) {
        m_data[i] += other.m_data[i];
    }
    return *this;
}

__spec Matrix& operator-(Matrix other) {
    Matrix ret(*this);
    // Check dimensions match
    if (m_row != other.m_row || m_col != other.m_col) {
        printf("Matrix dimension mismatch in operator+=\n");
        return *this; // or handle error properly
    }
    for (int i = 0; i < m_row * m_col; ++i) {
        ret.m_data[i] -= other.m_data[i];
    }
    return ret;
}


// Matrix& Matrix::operator+=(const Matrix& other) {
    // // your implementation here
    // return *this;
// }


}; //MATRIX


__spec Matrix Identity(const int &d){
    Matrix ret (d,d);
    for (int r=0;r<ret.m_row;r++){
        ret.Set(r,r,1.0);
    }
    return ret;
}
//// OPERATION WITH MATRIX CREATION COULD PRODUCE MEM LEAKING
__spec Matrix MatMul(const Matrix &A, const Matrix &B){
  Matrix ret(A.m_row,B.m_col);
  for (int i = 0; i<A.m_row; i++)
    for (int j = 0; j<B.m_col; j++)
      for (int k = 0; k<A.m_col; k++)
        //ret.m_data[i * A.m_ + j] += A.m_data[i * A.m_row + k] * B.m_data[k * B.m_row + j ];
        ret.getVal(i,j) += A.at(i,k) * B.at(k,j); //AT IS CONST
  
  
  
  return ret;
}

// __spec Matrix Matrix::operator=(const Matrix &A){
  // for (int i=0;i<A.m_row*A.m_col;i++) this->m_data[i] = A.m_data[i];
  // this->m_row=A.m_row;this->m_col=A.m_col;
  // return *this;
// }

__spec void MatMul(Matrix A, Matrix B, Matrix *ret){
  for (int i=0;i<A.m_row*B.m_col;i++) ret->m_data[i] = 0.0;
  for (int i = 0; i<A.m_row; i++)
    for (int j = 0; j<B.m_col; j++)
      for (int k = 0; k<A.m_col; k++)
        //ret->m_data[i * A.m_row + j] += A.m_data[i * A.m_row + k] * B.m_data[k * B.m_row + j ];
        ret->getVal(i,j) += A.getVal(i,k) * B.getVal(k,j); 

}

  __spec double & Matrix::getVal(int a, int b){
    return m_data[m_col*a+b];
  }

  __spec double & Matrix::at(int a, int b)const{
    return m_data[m_col*a+b];
  }  
  __spec double & Matrix::operator()(int a, int b){
    return m_data[m_col*a+b];
  }

  __spec void Matrix::Set(const int r, const int c, const double d){
    if (r>m_row) printf("ERROR, trying to set row %d, max row %d",r, m_row);
    m_data[m_col*r+c] = d;
  }

  __spec void Matrix::SetZero(){
      for (int i=0;i<m_row;i++)
        for (int j=i;j<m_col;j++)
          m_data[m_col*j+i] = 0.0;
  }
  
  __spec Matrix & Matrix::Transpose(){
    double t;
    for (int i=0;i<m_row;i++)
      for (int j=i+1;j<m_col;j++){
        t= this->getVal(i,j); this->Set(i,j,this->getVal(j,i)); this->Set(j,i,t);}
    int temp = this->m_col;
    this->m_col = this->m_row;
    this->m_row = temp;
    return *this;
  }

  __spec Matrix Matrix::getTranspose(){
    Matrix ret(m_col,m_row);
    double t;
    for (int i=0;i<m_row;i++)
      for (int j=0;j<m_col;j++){
        ret.Set(j,i,this->getVal(i,j)); 
      }
    return ret;
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
			printf("%.6e ", getVal(i,j)) ;
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


__spec Matrix Matrix::Inv(){
  Matrix invA(m_row,m_col);
  if (this->m_dim ==2){
    
    
  } else if (this->m_dim == 3) {
    
    Matrix *cofactor = new Matrix(3,3);
    
    cofactor->Set(0,0, (this->getVal(1,1)*this->getVal(2,2)-this->getVal(1,2)*this->getVal(2,1)) );
    cofactor->Set(0,1,-(this->getVal(1,0)*this->getVal(2,2)-this->getVal(1,2)*this->getVal(2,0)) );
    cofactor->Set(0,2, (this->getVal(1,0)*this->getVal(2,1)-this->getVal(1,1)*this->getVal(2,0)) );
    cofactor->Set(1,0,-(this->getVal(0,1)*this->getVal(2,2)-this->getVal(0,2)*this->getVal(2,1)) );
    cofactor->Set(1,1, (this->getVal(0,0)*this->getVal(2,2)-this->getVal(0,2)*this->getVal(2,0)) );
    cofactor->Set(1,2,-(this->getVal(0,0)*this->getVal(2,1)-this->getVal(0,1)*this->getVal(2,0)) );
    cofactor->Set(2,0, (this->getVal(0,1)*this->getVal(1,2)-this->getVal(0,2)*this->getVal(1,1)) );
    cofactor->Set(2,1,-(this->getVal(0,0)*this->getVal(1,2)-this->getVal(0,2)*this->getVal(1,0)) );
    cofactor->Set(2,2, (this->getVal(0,0)*this->getVal(1,1)-this->getVal(0,1)*this->getVal(1,0)) );    
    
    cofactor->Transpose();
    //printf ("COFACTOR: \n");
    //cofactor->Print();
    double f = 1.0/this->calcDet();
    invA = cofactor->Mul(f);
    //printf("INVA\n");
    //invA.Print();
    //for (int i=0;i<cofactor->m_row*cofactor->m_col;i++) invA->m_data[i] = cofactor->Mul(1.0/A.calcDet()).m_data[i];
    delete cofactor;
  }
  return invA;
}



__spec Matrix InvMat(Matrix &A, Matrix *invA){
  //Matrix invA(A.m_row,A.m_col);
  if (A.m_dim ==2){
    
    
  } else if (A.m_dim == 3) {
    
    Matrix *cofactor = new Matrix(3,3);
    
    cofactor->Set(0,0, (A(1,1)*A(2,2)-A(1,2)*A(2,1)) );
    cofactor->Set(0,1,-(A(1,0)*A(2,2)-A(1,2)*A(2,0)) );
    cofactor->Set(0,2, (A(1,0)*A(2,1)-A(1,1)*A(2,0)) );
    cofactor->Set(1,0,-(A(0,1)*A(2,2)-A(0,2)*A(2,1)) );
    cofactor->Set(1,1, (A(0,0)*A(2,2)-A(0,2)*A(2,0)) );
    cofactor->Set(1,2,-(A(0,0)*A(2,1)-A(0,1)*A(2,0)) );
    cofactor->Set(2,0, (A(0,1)*A(1,2)-A(0,2)*A(1,1)) );
    cofactor->Set(2,1,-(A(0,0)*A(1,2)-A(0,2)*A(1,0)) );
    cofactor->Set(2,2, (A(0,0)*A(1,1)-A(0,1)*A(1,0)) );    
    
    cofactor->Transpose();
    //printf ("COFACTOR: \n");
    cofactor->Print();
    Matrix *temp = new Matrix(3,3);
    //*temp=cofactor->Mul(1.0/A.calcDet());
    *invA = cofactor->Mul(1.0/A.calcDet());
    //THE ONE WHICH WORKS, OPERATOR= DOES NOT WORK
    //for (int i=0;i<A.m_row*A.m_col;i++){invA->m_data[i]=temp->m_data[i];}
    //printf("INVA\n");
    invA->Print();
    //for (int i=0;i<cofactor->m_row*cofactor->m_col;i++) invA->m_data[i] = cofactor->Mul(1.0/A.calcDet()).m_data[i];
    delete cofactor, temp;
  }
  return *invA;
}



__spec Matrix AdjMat(Matrix &A, Matrix *invA){
  //Matrix invA(A.m_row,A.m_col);
  if (A.m_dim ==2){
    invA->Set(0,0,A.getVal(1,1)); invA->Set(0,1,A.getVal(1,0));
    invA->Set(1,0,A.getVal(0,1)); invA->Set(1,1,A.getVal(1,1));
    
  } else if (A.m_dim == 3) {
    
    Matrix *cofactor = new Matrix(3,3);
    
    cofactor->Set(0,0, (A(1,1)*A(2,2)-A(1,2)*A(2,1)) );
    cofactor->Set(0,1,-(A(1,0)*A(2,2)-A(1,2)*A(2,0)) );
    cofactor->Set(0,2, (A(1,0)*A(2,1)-A(1,1)*A(2,0)) );
    cofactor->Set(1,0,-(A(0,1)*A(2,2)-A(0,2)*A(2,1)) );
    cofactor->Set(1,1, (A(0,0)*A(2,2)-A(0,2)*A(2,0)) );
    cofactor->Set(1,2,-(A(0,0)*A(2,1)-A(0,1)*A(2,0)) );
    cofactor->Set(2,0, (A(0,1)*A(1,2)-A(0,2)*A(1,1)) );
    cofactor->Set(2,1,-(A(0,0)*A(1,2)-A(0,2)*A(1,0)) );
    cofactor->Set(2,2, (A(0,0)*A(1,1)-A(0,1)*A(1,0)) );    
    
    cofactor->Transpose();
    //printf ("COFACTOR: \n");
    //cofactor->Print();

    *invA = *cofactor;
    //THE ONE WHICH WORKS, OPERATOR= DOES NOT WORK
    //for (int i=0;i<A.m_row*A.m_col;i++){invA->m_data[i]=cofactor->m_data[i];}
    //printf("ADJA\n");
    //invA->Print();
    //for (int i=0;i<cofactor->m_row*cofactor->m_col;i++) invA->m_data[i] = cofactor->Mul(1.0/A.calcDet()).m_data[i];
    delete cofactor;
  }
  return *invA;
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

__spec int SetMatVals(Matrix *mat, int argcount, ...){

  int total = 0;

  /* Declare a variable of type va_list. */
  va_list argptr;

  /* Initialize that variable.. */
  va_start (argptr, argcount);
  printf("arg count %d",argcount);
  for (int counter = 0; counter < argcount; counter++)
  {
    printf("pos %d %lf",counter, va_arg (argptr, double));
    //mat->m_data[counter] = va_arg (argptr, double);
  }

  /* End use of the argptr variable. */
  va_end (argptr);

  return total;

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

///// ONLY FOR SMALL PROBLEMS AND PROTOTYPNG. TOO SLOW
__spec Matrix SolveLinearSystem(Matrix A, Matrix b) {
    if (A.m_row != A.m_col || b.m_col != 1 || A.m_row != b.m_row) {
        printf("Invalid dimensions for linear system!\n");
        // Return an empty matrix or handle error
        return Matrix();
    }

    int n = A.m_row;
    Matrix Ab(n, n + 1); // Augmented matrix [A | b]

    // Copy A and b into Ab
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Ab(i, j) = A(i, j);
        }
        Ab(i, n) = b(i, 0);
    }

    // Gaussian elimination
    for (int i = 0; i < n; ++i) {
        // Pivot
        double maxEl = fabs(Ab(i, i));
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (fabs(Ab(k, i)) > maxEl) {
                maxEl = fabs(Ab(k, i));
                maxRow = k;
            }
        }

        // Swap rows
        for (int k = i; k < n + 1; ++k) {
            double tmp = Ab(maxRow, k);
            Ab(maxRow, k) = Ab(i, k);
            Ab(i, k) = tmp;
        }

        // Eliminate below
        for (int k = i + 1; k < n; ++k) {
            double factor = Ab(k, i) / Ab(i, i);
            for (int j = i; j < n + 1; ++j) {
                Ab(k, j) -= factor * Ab(i, j);
            }
        }
    }

    // Back substitution
    Matrix x(n, 1);
    for (int i = n - 1; i >= 0; --i) {
        x(i, 0) = Ab(i, n);
        for (int j = i + 1; j < n; ++j) {
            x(i, 0) -= Ab(i, j) * x(j, 0);
        }
        x(i, 0) /= Ab(i, i);
    }

    return x;
}

#endif
