/**********************************************************************************
 *                                                                                *
 *  DynELA Finite Element Code v.4.0                                              *
 *  by Olivier PANTALE                                                            *
 *  Olivier.Pantale@enit.fr                                                       *
 *                                                                                *
 *********************************************************************************/
//@!CODEFILE = DynELA-C-file
//@!BEGIN = PRIVATE

#include <lapacke.h>
#include <cblas.h>

#include <Matrix.h>
#include <NumpyInterface.h>
#include <MatrixDiag.h>
#include <Tensor2.h>

/*
@LABEL:Matrix::Matrix(long r, long c, double m)
@SHORT:Constructor of the Matrix class with initialization.
@RETURN:Matrix & The initialized matrix.
@ARG:long & r & Number of rows of the matrix to create.
@ARG:long & c & Number of cols of the matrix to create.
@ARG:double & m & Value to give to each element of the new matrix.
This method creates a new matrix of size $r \times c$ where all values are initialized to the scalar value $m$.
\begin{equation*}
\M=\left[\begin{array}{cccc}
  M_{11}=m & M_{12}=m & \hdots & M_{1c}=m\\
  M_{21}=m & M_{22}=m & \hdots & M_{2c}=m\\
  \vdots & \vdots & \hdots & \vdots\\
  M_{r1}=m & M_{r1}=m & \hdots & M_{rc}=m
  \end{array}\right]
\end{equation*}
where the $\M$ is the object itself and $m$ is scalar value defined by parameter m.
@END
*/
//-----------------------------------------------------------------------------
Matrix::Matrix(const long rows, const long cols, const double value)
//-----------------------------------------------------------------------------
{
  // allocation memoire;
  allocate(rows, cols);

  // remplissage initial
  setValue(value);
}

// Copy constructor
//-----------------------------------------------------------------------------
Matrix::Matrix(const Matrix &mat)
//-----------------------------------------------------------------------------
{
  // memory allocation
  allocate(mat._rows, mat._cols);

  memcpy(_data, mat._data, _dataLength * sizeof(double));
}

/*
@LABEL:Matrix::Matrix(int r, int c, double x1, double x2, ...)
@SHORT:Constructor of the matrix class with initialization.
@RETURN:Matrix & The initialized matrix.
@ARG:int & r & Number of rows of the matrix to create.
@ARG:int & c & Number of cols of the matrix to create.
@ARG:double & x1 & first component of the matrix to create.
@ARG:double & x2 & second component of the matrix to create.
@ARG:double & ... & $x^{th}$ component of the matrix to create.
All components have to be defined in row order.
\begin{equation*}
\M=\left[\begin{array}{cccc}
  M_{11}=x_1 & M_{12}=x_2 & \hdots & M_{1c}=x_c\\
  M_{21}=x_{c+1} & M_{22}=x_{c+2} & \hdots & M_{2c}=x_{2c}\\
  \vdots & \vdots & \hdots & \vdots\\
  M_{r1}=x_{(r-1)c+1} & M_{r1}=x_{(r-1)c+2} & \hdots & M_{rc}=x_{rc}
  \end{array}\right]
\end{equation*}
@END
*/
//-----------------------------------------------------------------------------
Matrix::Matrix(int rows, int cols, double firstValue, double secondValue, ...)
//-----------------------------------------------------------------------------
{
  // allocation memoire;
  allocate(rows, cols);

  // traitement de la premiere valeur
  _data[0] = firstValue;
  _data[1] = secondValue;

  // A place to store the list of arguments
  va_list arguments;

  // Store arguments values
  va_start(arguments, secondValue);

  long skip = 0;
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
    {
      if (skip >= 2)
        _data[skip] = va_arg(arguments, double);
      skip++;
    }

  // Cleans up the list
  va_end(arguments);
}

// Destructor
//-----------------------------------------------------------------------------
Matrix::~Matrix()
//-----------------------------------------------------------------------------
{
  // desallocation de memoire
  deallocate();
}

// Memory allocation
//-----------------------------------------------------------------------------
void Matrix::allocate(const long rows, const long cols)
//-----------------------------------------------------------------------------
{
  _rows = rows;
  _cols = cols;
  _dataLength = _rows * _cols;
  _data = new double[_dataLength];

#ifdef VERIF_alloc
  if (!_data)
  {
    internalFatalErrorLine("Matrix::allocate()",
                           "new double[%d] for a [%d,%d] Matrix allocation Error\n"
                           "Seems to have an overflow memory error\n"
                           "Check your memory size, and memory consumption first\n",
                           _dataLength, _rows, _cols);
  }
#endif
}

// Memory deallocation
//-----------------------------------------------------------------------------
void Matrix::deallocate()
//-----------------------------------------------------------------------------
{
  if (_dataLength >= 0)
    delete[] _data;

  _rows = 0;
  _cols = 0;
  _dataLength = 0;
}

/*
@LABEL:Matrix::redim(long r, long c)
@SHORT:Change the allocation size of a matrix.
@ARG:long & r & Number of rows of the matrix.
@ARG:long & c & Number of cols of the matrix.
@WARNING:This method cleans the content of the matrix.
This method changes the size of a matrix.
If the new size is the same as the actual size, this method does nothing.
@END
*/
//-----------------------------------------------------------------------------
void Matrix::redim(const long rows, const long cols)
//-----------------------------------------------------------------------------
{
  if ((rows == _rows) && (cols == _cols))
    return;

  deallocate();
  allocate(rows, cols);
}

// Send the content of a matrix to the output flux for display
//-----------------------------------------------------------------------------
std::ostream &operator<<(std::ostream &os, const Matrix &mat)
//-----------------------------------------------------------------------------
{
  mat.print(os);
  return os;
}

// Print the content of a matrix to the output flux for display
//-----------------------------------------------------------------------------
void Matrix::print(std::ostream &os) const
//-----------------------------------------------------------------------------
{
  long i, j;
  os << "matrix " << _rows << "x" << _cols << " ={{";
  for (i = 0; i < _rows; i++)
  {
    for (j = 0; j < _cols; j++)
    {
      if (j != 0)
        os << ",";
      os << _data[Ind(i, j, _rows, _cols)];
    }
    if (i != _rows - 1)
      os << "},{";
  }
  os << "}}";
}

/*
@LABEL:Matrix::setValue(double v)
@SHORT:Fill a matrix with a scalar value.
@ARG:double & v & Value to use for the operation.
This method is a surdefinition of the = operator for the matrix class.
@END
*/
//-----------------------------------------------------------------------------
void Matrix::setValue(double val)
//-----------------------------------------------------------------------------
{
  long i;

  for (i = 0; i < _dataLength; i++)
    _data[i] = val;
}

/*
@LABEL:Matrix::setToUnity()
@SHORT:Unity matrix.
@SMOD
@WARNING:This only works on a square matrix $n \times n$.
This method transforms the current matrix to a unity matrix.
\begin{equation*}
\M=\left[\begin{array}{cccc}
  M_{11}=1 & M_{12}=0 & \hdots & M_{1n}=0\\
  M_{21}=0 & M_{22}=1 & \hdots & M_{2n}=0\\
  \vdots & \vdots & \hdots & \vdots\\
  M_{n1}=0 & M_{n1}=0 & \hdots & M_{nn}=1
  \end{array}\right]
\end{equation*}
where $\M$ is a square matrix defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
void Matrix::setToUnity()
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  if (_rows != _cols)
    fatalError(" Matrix::setToUnity()",
               "Your matrix is not a square matrix, it's a [%d,%d] matrix",
               _rows, _cols);
#endif

  long i, j;

  for (i = 0; i < _rows; i++)
    for (j = 0; j < _cols; j++)
      _data[Ind(i, j, _rows, _cols)] = (double)(i == j);
}

/*
@LABEL:Matrix::maxVal()
@SHORT:Maximum component in a matrix.
@RETURN:double & The maximum component of the matrix.
This method returns the maximum component in a matrix.
@END
*/
//-----------------------------------------------------------------------------
double Matrix::maxVal() const
//-----------------------------------------------------------------------------
{
  double max = _data[0];
  for (int i = 1; i < _dataLength; i++)
  {
    if (_data[i] > max)
      max = _data[i];
  }
  return max;
}

/*
@LABEL:Matrix::minVal()
@SHORT:Minimum component in a matrix.
@RETURN:double & The minimum component of the matrix.
This method returns the minimum component in a matrix.
@END
*/
//-----------------------------------------------------------------------------
double Matrix::minVal() const
//-----------------------------------------------------------------------------
{
  double min = _data[0];
  for (int i = 1; i < _dataLength; i++)
  {
    if (_data[i] < min)
      min = _data[i];
  }
  return min;
}

/*
@LABEL:Matrix::maxAbs()
@SHORT:Maximum absolute component in a matrix.
@RETURN:double & The maximum component of the matrix.
This method returns the maximum absolute component in a matrix.
@END
*/
//-----------------------------------------------------------------------------
double Matrix::maxAbs() const
//-----------------------------------------------------------------------------
{
  double max = dnlAbs(_data[0]);
  for (int i = 1; i < _dataLength; i++)
  {
    if (dnlAbs(_data[i]) > max)
      max = dnlAbs(_data[i]);
  }
  return max;
}

/*
@LABEL:Matrix::minAbs()
@SHORT:Minimum absolute component in a matrix.
@RETURN:double & The minimum component of the matrix.
This method returns the minimum absolute component in a matrix.
@END
*/
//-----------------------------------------------------------------------------
double Matrix::minAbs() const
//-----------------------------------------------------------------------------
{
  double min = dnlAbs(_data[0]);
  for (int i = 1; i < _dataLength; i++)
  {
    if (dnlAbs(_data[i]) < min)
      min = dnlAbs(_data[i]);
  }
  return min;
}

/*
@LABEL:Matrix::operator=(double m)
@SHORT:Fill a matrix with a scalar value.
@RETURN:Matrix&
@ARG:double & m & Value to use for the operation.
This method is a surdefinition of the = operator for the matrix class.
\begin{equation*}
\M=\left[\begin{array}{cccc}
  M_{11}=m & M_{12}=m & \hdots & M_{1c}=m\\
  M_{21}=m & M_{22}=m & \hdots & M_{2c}=m\\
  \vdots & \vdots & \hdots & \vdots\\
  M_{r1}=m & M_{r1}=m & \hdots & M_{rc}=m
  \end{array}\right]
\end{equation*}
where $\M$ is a matrix defined by the object itself and $m$ is the scalar value defined by parameter m.
@END
*/
//-----------------------------------------------------------------------------
Matrix &Matrix::operator=(const double &val)
//-----------------------------------------------------------------------------
{
  setValue(val);
  return *this;
}

// Copy the content of a matrix into a new one
//-----------------------------------------------------------------------------
Matrix &Matrix::operator=(const Matrix &mat)
//-----------------------------------------------------------------------------
{
  // redimension de la matrice resultat
  redim(mat._rows, mat._cols);

  // recopie des valeurs
  memcpy(_data, mat._data, _dataLength * sizeof(double));

  return *this;
}

/*
@LABEL:Matrix::operator=(Tensor2 T)
@SHORT:Copy the content of a Tensor2 into a Matrix.
@RETURN:Matrix&
@ARG:Tensor2 & T & Second order tensor to copy.
The result of this operation is a matrix as a copy of a second order tensor where $\T$ is a matrix defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Matrix &Matrix::operator=(const Tensor2 &tens)
//-----------------------------------------------------------------------------
{
  // redimension de la matrice resultat
  redim(3, 3);

  // recopie des valeurs
  memcpy(_data, tens._data, _dataLength * sizeof(double));

  return *this;
}

/*
@LABEL:Matrix::operator+(Matrix B)
@SHORT:Addition of 2 matrices.
@ARG:Matrix & B & Matrix to add to the current one.
@RETURN:Matrix & Result of the addition operation.
This method defines the addition of 2 matrices.
The result of this operation is also a matrix defined by:
\begin{equation*}
\T = \A + \B
\end{equation*}
where $\A$ is a matrix defined by the object itself and $\B$ is the matrix defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
Matrix Matrix::operator+(const Matrix &mat) const
//-----------------------------------------------------------------------------
{
  long i;

#ifdef VERIF_maths
  if ((_rows != mat._rows) && (_cols != mat._cols))
    fatalError("Matrix::operator -=",
               "matrix sizes incompatible\n"
               "You're about to substract a [%d,%d] and a [%d,%d] matrix",
               _rows, _cols, mat._rows, mat._cols);
#endif

  // creation d'une nouvelle matrice
  Matrix resu(mat._rows, mat._cols);

  // calcul de la somme
  for (i = 0; i < _dataLength; i++)
    resu._data[i] = _data[i] + mat._data[i];

  // renvoi de l'objet
  return resu;
}

/*
@LABEL:Matrix::operator-(Matrix B)
@SHORT:Difference of 2 matrices.
@ARG:Matrix & B & Matrix to subtract to the current one.
@RETURN:Matrix & Result of the difference operation.
This method defines the difference of 2 matrices.
The result of this operation is also a matrix defined by:
\begin{equation*}
\T = \A - \B
\end{equation*}
where $\A$ is a matrix defined by the object itself and $\B$ is the matrix defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
Matrix Matrix::operator-(const Matrix &mat) const
//-----------------------------------------------------------------------------
{
  long i;

#ifdef VERIF_maths
  if ((_rows != mat._rows) && (_cols != mat._cols))
    fatalError("Matrix::operator -=",
               "matrix sizes incompatible\n"
               "You're about to substract a [%d,%d] and a [%d,%d] matrix",
               _rows, _cols, mat._rows, mat._cols);
#endif

  // creation d'une nouvelle matrice
  Matrix resu(_rows, _cols);

  // calcul de la somme
  for (i = 0; i < _dataLength; i++)
    resu._data[i] = _data[i] - mat._data[i];

  // renvoi de l'objet
  return resu;
}

/*
@LABEL:Matrix::operator-()
@SHORT:Opposite of a matrix.
@RETURN:Matrix & The opposite matrix.
This method defines the opposite of a matrix.
The result of this operation is also a matrix defined by:
\begin{equation*}
\T = - \A
\end{equation*}
where $\A$ is a matrix defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Matrix Matrix::operator-() const
//-----------------------------------------------------------------------------
{
  long i;

  // creation d'une nouvelle matrice
  Matrix resu(_rows, _cols);

  // calcul de la somme
  for (i = 0; i < _dataLength; i++)
    resu._data[i] = -_data[i];

  // renvoi de l'objet
  return resu;
}

/*
@LABEL:Matrix::operator+=(Matrix B)
@SHORT:Addition of 2 matrices.
@ARG:Matrix & B & Second matrix to add to the current one.
@RETURN:Matrix & Result of the addition operation.
This method defines the addition of 2 matrices.
The result of this operation is also a matrix defined by:
\begin{equation*}
\A += \B
\end{equation*}
where $\A$ is a matrix defined by the object itself and $\B$ is the matrix defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
void Matrix::operator+=(const Matrix &mat)
//-----------------------------------------------------------------------------
{
  long i;

#ifdef VERIF_maths
  if ((_rows != mat._rows) && (_cols != mat._cols))
    fatalError("Matrix::operator -=",
               "matrix sizes incompatible\n"
               "You're about to substract a [%d,%d] and a [%d,%d] matrix",
               _rows, _cols, mat._rows, mat._cols);
#endif

  // calcul de la somme
  for (i = 0; i < _dataLength; i++)
    _data[i] += mat._data[i];
}

/*
@LABEL:Matrix::operator-=(Matrix B)
@SHORT:Difference of 2 matrices.
@ARG:Matrix & B & Second matrix to add to the current one.
@RETURN:Matrix & Result of the difference operation.
This method defines the difference of 2 matrices.
The result of this operation is also a matrix defined by:
\begin{equation*}
\A -= \B
\end{equation*}
where $\A$ is a matrix defined by the object itself and $\B$ is the matrix defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
void Matrix::operator-=(const Matrix &mat)
//-----------------------------------------------------------------------------
{
  long i;

#ifdef VERIF_maths
  if ((_rows != mat._rows) && (_cols != mat._cols))
    fatalError("Matrix::operator -=",
               "matrix sizes incompatible\n"
               "You're about to substract a [%d,%d] and a [%d,%d] matrix",
               _rows, _cols, mat._rows, mat._cols);
#endif

  // calcul de la somme
  for (i = 0; i < _dataLength; i++)
    _data[i] -= mat._data[i];
}

/*
@LABEL:Matrix::operator*(double l)
@SHORT:Multiplication of a matrix by a scalar.
@ARG:double & l & Scalar value to use for the operation.
@RETURN:Matrix & Result of the multiplication operation.
This method defines the multiplication of a matrix by a scalar value.
The result of this operation is also a matrix defined by:
\begin{equation*}
\T = \lambda \A
\end{equation*}
where $\A$ is a matrix defined by the object itself and $\lambda$ is the scalar value defined by parameter l.
@END
*/
//-----------------------------------------------------------------------------
Matrix Matrix::operator*(const double &lambda) const
//-----------------------------------------------------------------------------
{
  long i;
  Matrix resu(_rows, _cols);

  for (i = 0; i < _dataLength; i++)
    resu._data[i] = lambda * _data[i];
  return resu;
}

/*
@LABEL:Matrix::operator*=(double l)
@SHORT:Multiplication of a matrix by a scalar.
@ARG:double & l & Scalar value to use for the operation.
@RETURN:Matrix & Result of the multiplication operation.
This method defines the multiplication of a matrix by a scalar value.
The result of this operation is also a matrix defined by:
\begin{equation*}
\lambda \A
\end{equation*}
where $\A$ is a matrix defined by the object itself and $\lambda$ is the scalar value defined by parameter l.
@END
*/
//-----------------------------------------------------------------------------
void Matrix::operator*=(const double &lambda)
//-----------------------------------------------------------------------------
{
  long i;

  for (i = 0; i < _dataLength; i++)
    _data[i] *= lambda;
}

/*
@LABEL:Matrix::operator/=(double l)
@SHORT:Division of a matrix by a scalar.
@ARG:double & l & Scalar value to use for the operation.
@RETURN:Matrix & Result of the division operation.
This method defines the division of a matrix by a scalar value.
The result of this operation is also a matrix defined by:
\begin{equation*}
\frac{1}{\lambda} \A
\end{equation*}
where $\A$ is a matrix defined by the object itself and $\lambda$ is the scalar value defined by parameter l.
@END
*/
//-----------------------------------------------------------------------------
void Matrix::operator/=(const double &lambda)
//-----------------------------------------------------------------------------
{
  long i;

  for (i = 0; i < _dataLength; i++)
    _data[i] /= lambda;
}

/*
@LABEL:Matrix::operator/(double l)
@SHORT:Division of a matrix by a scalar.
@ARG:double & l & Scalar value to use for the operation.
@RETURN:Matrix & Result of the division operation.
This method defines the division of a matrix by a scalar value.
The result of this operation is also a matrix defined by:
\begin{equation*}
\T = \frac{1}{\lambda} \A
\end{equation*}
where $\A$ is a matrix defined by the object itself and $\lambda$ is the scalar value defined by parameter l.
@END
*/
//-----------------------------------------------------------------------------
Matrix Matrix::operator/(const double &lambda) const
//-----------------------------------------------------------------------------
{
  long i;
  Matrix resu(_rows, _cols);

  for (i = 0; i < _dataLength; i++)
    resu._data[i] = _data[i] / lambda;
  return resu;
}

/*
@LABEL:operator*(double l, Matrix A)
@SHORT:Multiplication of a matrix by a scalar.
@ARG:double & l & Scalar value to use for the operation.
@ARG:Matrix & A & Matrix to use for the operation.
@RETURN:Matrix & Result of the multiplication operation.
This method defines the multiplication of a matrix by a scalar value.
The result of this operation is also a matrix defined by:
\begin{equation*}
\T = \lambda \A
\end{equation*}
where $\A$ is a matrix and $\lambda$ is the scalar value defined by parameter l.
@END
*/
//-----------------------------------------------------------------------------
Matrix operator*(const double &lambda, const Matrix &mat)
//-----------------------------------------------------------------------------
{
  long i;
  Matrix resu(mat._rows, mat._cols);

  for (i = 0; i < mat._dataLength; i++)
    resu._data[i] = lambda * mat._data[i];
  return resu;
}

/*
@LABEL:Matrix::operator*(Matrix B)
@SHORT:Single contracted product of two matrices.
@RETURN:Matrix & Result of the multiplication operation.
@ARG:Matrix & B & Second matrix for the multiplication operation.
This method defines a single contracted product of two matrices.
The result of this operation is also a matrix defined by:
\begin{equation*}
\T = \A \cdot \B
\end{equation*}
where $\A$ is a matrix defined by the object itself and $\B$ is the matrix defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
Matrix Matrix::operator*(const Matrix &mat) const
//-----------------------------------------------------------------------------
{
  return dot(mat);
}

/*
@LABEL:Matrix::dot()
@SHORT:Single contracted product of a matrix by itself.
@RETURN:Matrix & Result of the multiplication operation.
This method defines a single contracted product of a matrix by itself.
The result of this operation is also a matrix defined by:
\begin{equation*}
\T = \A \cdot \A
\end{equation*}
where $\A$ is a matrix defined by the object itself.
This method uses the Blas \textsf{dgemm} Fortran subroutine to perform the operation.
@END
*/
//-----------------------------------------------------------------------------
Matrix Matrix::dot() const
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  if (_rows != _cols)
    fatalError(" Matrix::dot()",
               "Your matrix is not a square matrix, it's a [%d,%d] matrix",
               _rows, _cols);
#endif

  Matrix resu(_rows, _rows);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, _rows, _rows, _rows, 1.0, _data, _rows, _data, _rows, 0.0, resu._data, _rows);

  return resu;
}

/*
@LABEL:Matrix::dotTxN()
@SHORT:Single contracted product of a matrix by its transpose.
@RETURN:Matrix & Result of the multiplication operation.
This method defines a single contracted product of a matrix by its transpose.
The result of this operation is also a matrix defined by:
\begin{equation*}
\T = \A^T\cdot \A
\end{equation*}
where $\A$ is a matrix defined by the object itself. Result is a symmetric matrix.
This method uses the Blas \textsf{dgemm} Fortran subroutine to perform the operation.
@END
*/
//-----------------------------------------------------------------------------
Matrix Matrix::dotTxN() const
//-----------------------------------------------------------------------------
{
  Matrix resu(_cols, _cols);

  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, _cols, _cols, _rows, 1.0, _data, _cols, _data, _cols, 0.0, resu._data, _cols);

  return resu;
}

/*
@LABEL:Matrix::dotNxT()
@SHORT:Single contracted product of a matrix by its transpose.
@RETURN:Matrix & Result of the multiplication operation.
This method defines a single contracted product of a matrix by its transpose.
The result of this operation is also a matrix defined by:
\begin{equation*}
\T = \A \cdot \A^T
\end{equation*}
where $\A$ is a matrix defined by the object itself. Result is a symmetric matrix.
This method uses the Blas \textsf{dgemm} Fortran subroutine to perform the operation.
@END
*/
//-----------------------------------------------------------------------------
Matrix Matrix::dotNxT() const
//-----------------------------------------------------------------------------
{
  Matrix resu(_rows, _rows);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, _rows, _rows, _cols, 1.0, _data, _cols, _data, _cols, 0.0, resu._data, _rows);

  return resu;
}

/*
@LABEL:Matrix::dot(Matrix B)
@SHORT:Single contracted product of two matrices.
@RETURN:Matrix & Result of the multiplication operation.
@ARG:Matrix & B & Second matrix for the multiplication operation.
This method defines a single contracted product of two matrices.
The result of this operation is also a matrix defined by:
\begin{equation*}
\T = \A \cdot \B
\end{equation*}
where $\A$ is a matrix defined by the object itself and $\B$ is the matrix defined by parameter B.
This method uses the Blas \textsf{dgemm} Fortran subroutine to perform the operation.
@END
*/
//-----------------------------------------------------------------------------
Matrix Matrix::dot(const Matrix mat) const
//-----------------------------------------------------------------------------
{
  Matrix resu(_rows, mat._cols);

#ifdef VERIF_maths
  if (_cols != mat.rows())
    fatalError("Matrix::squareMultiplyBy (const MatrixDiag & mat)",
               "matrix sizes incompatible\n"
               "You're about to multiply a [%d,%d] and a [%d,%d] matrix",
               _rows, _cols, mat.rows(), mat.cols());
#endif

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, _rows, mat._cols, _cols, 1.0, _data, _cols, mat._data, mat._cols, 0.0, resu._data, mat._cols);

  return resu;
}

/*
@LABEL:Matrix::doubleDot()
@SHORT:Double contracted product of a matrix by itself.
@RETURN:double & Result of the multiplication operation.
This method defines a double contracted product of a matrix by itself.
The result of this operation is a scalar $s$ defined by:
\begin{equation*}
s = \A : \A = \sum_{i=1}^{n} \sum_{j=1}^{m} A_{ij}\times A_{ij}
\end{equation*}
where $\A$ is a matrix defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
double Matrix::doubleDot() const
//-----------------------------------------------------------------------------
{
  double sum = 0.0;
  for (long i = 0; i < _rows; i++)
    for (long j = 0; j < _cols; j++)
      sum += _data[Ind(i, j, _rows, _cols)] * _data[Ind(i, j, _rows, _cols)];
  return (sum);
}

/*
@LABEL:Matrix::doubleDot(Matrix B)
@SHORT:Double contracted product of 2 matrices.
@RETURN:double & Result of the multiplication operation.
@ARG:Matrix & B & Second matrix for the multiplication operation.
This method defines a double contracted product of two matrices.
The result of this operation is a scalar $s$ defined by:
\begin{equation*}
s = \A : \B = \sum_{i=1}^{n} \sum_{j=1}^{m} A_{ij}\times B_{ij}
\end{equation*}
where $\A$ is a matrix defined by the object itself and $\B$ is a matrix defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
double Matrix::doubleDot(const Matrix mat) const
//-----------------------------------------------------------------------------
{
  double sum = 0.0;
  for (long i = 0; i < _rows; i++)
    for (long j = 0; j < _cols; j++)
      sum += _data[Ind(i, j, _rows, _cols)] * mat._data[Ind(i, j, _rows, _cols)];
  return (sum);
}

// multiplication de deux matrices
/*
  Cette methode permet de surdefinir l'operation de multiplication des matrices et d'ecrire simplement la multiplication de deux matrice sous la forme donnee en exemple. Cette operation correspond à la notion de produit contracte de deux matrices.

  Exemple :
  \code
  Matrix t1,t2,t3;
  t3=t1*t2; // produit contracte
  \endcode
*/
//-----------------------------------------------------------------------------
Matrix Matrix::operator*(const MatrixDiag &mat) const
//-----------------------------------------------------------------------------
{
  long i, k;

  Matrix resu(_rows, mat.cols());

  for (i = 0; i < _rows; i++)
  {
    for (k = 0; k < mat.cols(); k++)
    {
      resu._data[Ind(i, k, _rows, mat.cols())] = _data[Ind(i, k, _rows, _cols)] * mat(k);
    }
  }

  return resu;
}

// multiplication d'une matrices par une matrice diagonale
/*
  Cette methode permet de surdefinir l'operation de multiplication d'une matrice par une matrice diagonale et de stocker le resultat dans la matrice sur laquelle on applique la methode. Les deux matrices doivent etre de meme dimension pour que cette methode soit utilisable.
  Exemple :
  \code
  Matrix t1;
  MatrixDiag t2;
  t1.squareMultiplyBy(t2); // produit
  \endcode
*/
//-----------------------------------------------------------------------------
void Matrix::squareMultiplyBy(const MatrixDiag &mat)
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  if (_cols != mat.rows())
    fatalError("Matrix::squareMultiplyBy (const MatrixDiag & mat)",
               "matrix sizes incompatible\n"
               "You're about to multiply a [%d,%d] and a [%d,%d] matrix",
               _rows, _cols, mat.rows(), mat.cols());
  if (_cols != _rows)
    fatalError("Matrix::squareMultiplyBy (const MatrixDiag & mat)",
               "matrix sizes incompatible\n"
               "Matrix must be a square one to do this operation");
#endif

  long k = 0;
  long i, j;
  for (i = 0; i < _rows; i++)
  {
    double m = mat._data[i];
    for (j = 0; j < _rows; j++)
    {
      _data[k++] *= m;
    }
  }
}

// multiplication d'une matrices par l'inverse d'une matrice diagonale
/*
  Cette methode permet de surdefinir l'operation de multiplication d'une matrice par l'inverse d'une matrice diagonale et de stocker le resultat dans la matrice sur laquelle on applique la methode. Les deux matrices doivent etre de meme dimension pour que cette methode soit utilisable.
  Exemple :
  \code
  Matrix t1;
  MatrixDiag t2;
  t1.squareDivideBy(t2); // produit
  \endcode
*/
//-----------------------------------------------------------------------------
void Matrix::squareDivideBy(const MatrixDiag &mat)
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  if (_cols != mat.rows())
    fatalError("Matrix::squareMultiplyBy (const MatrixDiag & mat)",
               "matrix sizes incompatible\n"
               "You're about to multiply a [%d,%d] and a [%d,%d] matrix",
               _rows, _cols, mat.rows(), mat.cols());
  if (_cols != _rows)
    fatalError("Matrix::squareMultiplyBy (const MatrixDiag & mat)",
               "matrix sizes incompatible\n"
               "Matrix must be a square one to do this operation");
#endif

  long k = _rows * _cols - 1;
  long i, j;
  for (i = _rows - 1; i >= 0; i--)
  {
    double m = 1.0 / mat._data[i];
    for (j = _rows - 1; j >= 0; j--)
    {
      _data[k--] *= m;
    }
  }
}

/*
@LABEL:Matrix::dotTxN(Vector V)
@SHORT:Multiplication of a matrix by a vector.
@RETURN:Vector & Result of the multiplication operation.
@ARG:Vector & V & Vector to use for the multiplication operation.
This method defines the product of a matrix by a vector.
The result of this operation is also a vector defined by:
\begin{equation*}
\overrightarrow{y} = \A^T \cdot \overrightarrow{x}
\end{equation*}
where $\A$ is a matrix defined by the object itself and $\overrightarrow{x}$ is a Vector defined by parameter V.
This method uses the Blas \textsf{dgemv} Fortran subroutine to perform the operation.
@END
*/
//-----------------------------------------------------------------------------
Vector Matrix::dotTxN(const Vector &vec) const
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  if (_rows != vec.size())
    fatalError("Matrix::operator *",
               "matrix and vector sizes incompatible\n"
               "You're about to multiply a [%d,%d] matrix and a [%d] vector",
               _rows, _cols, vec.size());
#endif

  Vector resu(_cols);

  cblas_dgemv(CblasRowMajor, CblasTrans, _rows, _cols, 1, _data, _cols, vec._data, 1, 0, resu._data, 1);

  return resu;
}

/*
@LABEL:Matrix::operator*(Vector V)
@SHORT:Multiplication of a matrix by a vector.
@RETURN:Vector & Result of the multiplication operation.
@ARG:Vector & V & Vector to use for the multiplication operation.
This method defines the product of a matrix by a vector.
The result of this operation is also a vector defined by:
\begin{equation*}
\overrightarrow{y} = \A \cdot \overrightarrow{x}
\end{equation*}
where $\A$ is a matrix defined by the object itself and $\overrightarrow{x}$ is a Vector defined by parameter V.
This method uses the Blas \textsf{dgemv} Fortran subroutine to perform the operation.
@END
*/
//-----------------------------------------------------------------------------
Vector Matrix::operator*(const Vector &vec) const
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  if (_cols != vec.size())
    fatalError("Matrix::operator *",
               "matrix and vector sizes incompatible\n"
               "You're about to multiply a [%d,%d] matrix and a [%d] vector",
               _rows, _cols, vec.size());
#endif

  Vector resu(_rows);

  cblas_dgemv(CblasRowMajor, CblasNoTrans, _rows, _cols, 1, _data, _cols, vec._data, 1, 0, resu._data, 1);

  return resu;
}

/*
@LABEL:Matrix::dot(Vector V)
@SHORT:Multiplication of a matrix by a vector.
@ARG:Vector & V & Vector to use for the multiplication operation.
@WARNING: The result of the operation is the parameter V itself.
This method defines the product of a matrix by a vector.
The result of this operation is also a vector defined by:
\begin{equation*}
\overrightarrow{y} = \A \cdot \overrightarrow{x}
\end{equation*}
where $\A$ is a matrix defined by the object itself and $\overrightarrow{x}$ is a Vector defined by parameter V.
This method uses the Blas \textsf{dgemv} Fortran subroutine to perform the operation.
@END
*/
//-----------------------------------------------------------------------------
void Matrix::dot(Vector &resu) const
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  if (_rows != _cols)
    fatalError(" Matrix::dot()",
               "Your matrix is not a square matrix, it's a [%d,%d] matrix",
               _rows, _cols);

  if (_cols != resu.size())
    fatalError("Matrix::operator *",
               "matrix and vector sizes incompatible\n"
               "You're about to multiply a [%d,%d] matrix and a [%d] vector",
               _rows, _cols, resu.size());
#endif

  Vector vec = resu;

  cblas_dgemv(CblasRowMajor, CblasNoTrans, _rows, _cols, 1, _data, _rows, vec._data, 1, 0, resu._data, 1);
}

/*
@LABEL:Matrix::trace()
@SHORT:Returns the trace of a matrix.
@RETURN:double & The trace of the matrix.
This method returns the trace of a matrix, i.e. the sum $s$ of all the terms of the diagonal:
\begin{equation*}
s = \tr[\M] = M_{11}+M_{22}+\hdots+M_{nn}
\end{equation*}
where $\M$ is a matrix defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
double Matrix::trace() const
//-----------------------------------------------------------------------------
{
  long i;

#ifdef VERIF_maths
  if (_rows != _cols)
    fatalError("Matrix::trace()",
               "Your matrix is not a square matrix, it's a [%d,%d] matrix",
               _rows, _cols);
#endif

  double trace = 0.0;
  for (i = 0; i < _rows; i++)
    trace += _data[Ind(i, i, _rows, _cols)];
  return trace;
}

/*
@LABEL:Matrix::averageTrace()
@SHORT:Returns the average trace of a matrix.
@RETURN:double & The average trace of the matrix.
This method returns the average trace of a matrix, i.e. the sum $s$ of all the terms of the diagonal divided by the number of terms on the diagonal:
\begin{equation*}
s = \frac{1}{n} \tr[\M] = \frac{1}{n} M_{11}+M_{22}+\hdots+M_{nn}
\end{equation*}
where $\M$ is a matrix defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
double Matrix::averageTrace() const
//-----------------------------------------------------------------------------
{
  long i;

#ifdef VERIF_maths
  if (_rows != _cols)
    fatalError("Matrix::averageTrace()",
               "Your matrix is not a square matrix, it's a [%d,%d] matrix",
               _rows, _cols);
#endif

  double trace = 0.0;
  for (i = 0; i < _rows; i++)
    trace += _data[Ind(i, i, _rows, _cols)];
  return (trace / (double)_rows);
}

/*
@LABEL:Matrix::transpose()
@SHORT:Transpose of a matrix.
@RETURN:Matrix & The transpose of the matrix.
This method defines the transpose of a matrix.
@END
*/
//-----------------------------------------------------------------------------
Matrix Matrix::transpose() const
//-----------------------------------------------------------------------------
{
  Matrix resu(_cols, _rows);
  long i, j;

  for (i = 0; i < _rows; i++)
    for (j = 0; j < _cols; j++)
      resu._data[Ind(j, i, _cols, _rows)] = _data[Ind(i, j, _rows, _cols)];
  return resu;
}

// Sum of rows of a matrix
//-----------------------------------------------------------------------------
Vector Matrix::rowSum() const
//-----------------------------------------------------------------------------
{
  Vector resu(_rows);

  for (long i = 0; i < _rows; i++)
    for (long j = 0; j < _cols; j++)
      resu(i) += _data[Ind(i, j, _rows, _cols)];

  return resu;
}

// Sum of cols of a matrix
//-----------------------------------------------------------------------------
Vector Matrix::colSum() const
//-----------------------------------------------------------------------------
{
  Vector resu(_cols);

  for (long i = 0; i < _rows; i++)
    for (long j = 0; j < _cols; j++)
      resu(j) += _data[Ind(i, j, _rows, _cols)];

  return resu;
}

/*
@LABEL:Matrix::symmetric()
@SHORT:Symmetric part of a matrix.
@RETURN:Matrix & The symmetric part of the matrix.
This method defines the symmetric part of a matrix with the following:
\begin{equation*}
\B = \frac{1}{2} \left( \A + \A^T \right)
\end{equation*}
where $\M$ is a matrix defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Matrix Matrix::symmetric() const
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  if (_rows != _cols)
    fatalError("Matrix::symmetric()",
               "Your matrix is not a square matrix, it's a [%d,%d] matrix",
               _rows, _cols);
#endif

  Matrix result(_rows, _rows);

  for (long j = 0; j < _rows; j++)
    for (long i = 0; i < _rows; i++)
      result._data[Ind(i, j, _rows, _rows)] = (_data[Ind(i, j, _rows, _rows)] + _data[Ind(j, i, _rows, _rows)]) / 2;

  return result;
}

/*
@LABEL:Matrix::skewSymmetric()
@SHORT:Skew-symmetric part of a matrix.
@RETURN:Matrix & The skew-symmetric part of the matrix.
This method defines the skew-symmetric part of a matrix with the following:
\begin{equation*}
\B = \frac{1}{2} \left( \A - \A^T \right)
\end{equation*}
where $\M$ is a matrix defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Matrix Matrix::skewSymmetric() const
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  if (_rows != _cols)
    fatalError("Matrix::skewSymmetric()",
               "Your matrix is not a square matrix, it's a [%d,%d] matrix",
               _rows, _cols);
#endif

  Matrix result(_rows, _rows);

  for (long j = 0; j < _rows; j++)
    for (long i = 0; i < _rows; i++)
      result._data[Ind(i, j, _rows, _rows)] = (_data[Ind(i, j, _rows, _rows)] - _data[Ind(j, i, _rows, _rows)]) / 2;

  return result;
}

/*
@LABEL:Matrix::row(short r)
@SHORT:Extraction of a row from a matrix.
@RETURN:Vec3D & The extracted row.
@ARG:short & r & Row to extract
This method returns a vector as part of a matrix.
The result of this operation with the argument r is a vector defined by:
\begin{equation*}
v_{i} = T_{ri}
\end{equation*}
where $\T$ is a matrix defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Vector Matrix::row(long row) const
//-----------------------------------------------------------------------------
{
  indexOK(row, 0);

  Vector resu(_cols);

  for (long j = 0; j < _cols; j++)
    resu(j) = _data[Ind(row, j, _rows, _cols)];

  return resu;
}

/*
@LABEL:Matrix::col(short c)
@SHORT:Extraction of a column from a matrix.
@RETURN:Vec3D & The extracted col.
@ARG:short & c & Column to extract
This method returns a vector as part of a matrix.
The result of this operation with the argument c is a vector defined by:
\begin{equation*}
v_{i} = T_{ic}
\end{equation*}
where $\T$ is a matrix defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Vector Matrix::col(long col) const
//-----------------------------------------------------------------------------
{
  indexOK(0, col);

  Vector resu(_rows);

  for (long i = 0; i < _rows; i++)
    resu(i) = _data[Ind(i, col, _rows, _cols)];

  return resu;
}

//  Test the equality of two matrices
//-----------------------------------------------------------------------------
bool Matrix::operator==(const Matrix &mat) const
//-----------------------------------------------------------------------------
{
  long i;

#ifdef VERIF_maths
  if ((_rows != mat._rows) && (_cols != mat._cols))
    fatalError("Matrix::operator ==",
               "The dimensions of the two matrix are'nt corresponding\n"
               "they are respectively a [%d,%d] and a [%d,%d] matrix",
               _rows, _cols, mat._rows, mat._cols);
#endif

  for (i = 0; i < _dataLength; i++)
    if (_data[i] != mat._data[i])
      return false;
  return true;
}

//  Test the inequality of two matrices
//-----------------------------------------------------------------------------
bool Matrix::operator!=(const Matrix &mat) const
//-----------------------------------------------------------------------------
{
  return !(*this == mat);
}

//  Writes a matrix in a binary flux for storage
//-----------------------------------------------------------------------------
void Matrix::write(std::ofstream &ofs) const
//-----------------------------------------------------------------------------
{
  ofs.write((char *)&_rows, sizeof(long));
  ofs.write((char *)&_cols, sizeof(long));
  ofs.write((char *)_data, _dataLength * sizeof(double));
}

//  Reads a matrix in a binary flux from storage
//-----------------------------------------------------------------------------
void Matrix::read(std::ifstream &ifs)
//-----------------------------------------------------------------------------
{
  long rows, cols;

  ifs.read((char *)&rows, sizeof(long));
  ifs.read((char *)&cols, sizeof(long));
  redim(rows, cols);
  ifs.read((char *)_data, _dataLength * sizeof(double));
}

//  Writes a matrix in a binary flux for storage
//-----------------------------------------------------------------------------
std::ofstream &operator<<(std::ofstream &os, const Matrix &mat)
//-----------------------------------------------------------------------------
{
  mat.write(os);
  return os;
}

//  Reads a matrix from a binary flux for storage
//-----------------------------------------------------------------------------
std::ifstream &operator>>(std::ifstream &is, Matrix &mat)
//-----------------------------------------------------------------------------
{
  mat.read(is);
  return is;
}

//-----------------------------------------------------------------------------
long Matrix::Memory() const
//-----------------------------------------------------------------------------
{
  return (_dataLength * sizeof(double));
}

/*
@LABEL:Matrix::det2()
@SHORT:Determinant of a $2 \times 2$ matrix.
@RETURN:double & The determinant.
This method returns the determinant of a $2 \times 2$ matrix.
The result of this operation is a scalar value defined by:
\begin{equation*}
d = T_{11} T_{22} - T_{21} T_{12}
\end{equation*}
where $\A$ is a matrix defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
double Matrix::det2() const
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  if ((_rows != 2) || (_cols != 2))
    fatalError("Matrix::det",
               "Your matrix is not a 2x2 square matrix, it's a [%d,%d] matrix",
               _rows, _cols);
#endif
  return _data[0] * _data[3] - _data[2] * _data[1];
}

/*
@LABEL:Matrix::det3()
@SHORT:Determinant of a $3 \times 3$ matrix.
@RETURN:double & The determinant.
This method returns the determinant of a $3 \times 3$ matrix.
The result of this operation is a scalar value defined by:
\begin{equation*}
d = A_{11} A_{22} A_{33} + A_{21} A_{32} A_{13} + A_{31} A_{12} A_{23} - A_{31} A_{22} A_{13} - A_{11} A_{32} A_{23} - A_{21} A_{12} A_{33}
\end{equation*}
where $\A$ is a matrix defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
double Matrix::det3() const
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  if ((_rows != 3) || (_cols != 3))
    fatalError("Matrix::det",
               "Your matrix is not a 3x3 square matrix, it's a [%d,%d] matrix",
               _rows, _cols);
#endif
  return _data[1] * _data[5] * _data[6] -
         _data[2] * _data[4] * _data[6] +
         _data[2] * _data[3] * _data[7] -
         _data[0] * _data[5] * _data[7] -
         _data[1] * _data[3] * _data[8] +
         _data[0] * _data[4] * _data[8];
}

/*
@LABEL:Matrix::det()
@SHORT:Determinant of a matrix.
@RETURN:double & The determinant.
This method returns the determinant of a matrix.
This method uses the Lapack \textsf{dgetrf} Fortran subroutine to perform the operation.
@END
*/
//-----------------------------------------------------------------------------
double Matrix::det() const
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  if (_rows != _cols)
    fatalError("Matrix::det",
               "Your matrix is not a square matrix, it's a [%d,%d] matrix",
               _rows, _cols);
#endif

  lapack_int n, m, info, lda;
  double Det = 1.0;

  n = _rows;
  m = _rows;
  lda = _rows;
  lapack_int *ipiv = new lapack_int[_rows];

  // init a matrix
  double *A = new double[_dataLength];
  memcpy(A, _data, _dataLength * sizeof(double));

  info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, A, lda, ipiv);
  if (info != 0)
  {
    if (info < 0)
    {
      internalFatalErrorLine("Matrix::det",
                             "parameter %d of function dgetrf_ has illegal value",
                             -info);
    }
    if (info > 0)
    {
      fatalError("Matrix::det",
                 "pivot (%d) is exactly zero\n"
                 "The factorization has been completed, but this factor is exactly singular, and division by zero will occur if it is used to solve a system of equations",
                 info);
    }
  }

  for (int i = 0; i < _rows; i++)
  {
    if (ipiv[i] == i + 1)
      Det *= A[Ind(i, i, _rows, _cols)];
    else
      Det *= -A[Ind(i, i, _rows, _cols)];
  }

  delete[] ipiv;
  delete[] A;

  return Det;
}

/*
@LABEL:Matrix::cofactors()
@SHORT:Cofactors of a matrix.
@RETURN:Matrix & The cofactor of the matrix.
This method returns the cofactor of a matrix defined by the following equation:
\begin{equation*}
\C = \det[\A] \cdot {\A^{-1}}^T
\end{equation*}
where $\T$ is a matrix defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Matrix Matrix::cofactors() const
//-----------------------------------------------------------------------------
{
  Matrix result = inverse();
  return (result.transpose() * det());
}

/*
@LABEL:Matrix::computeInverse2x2(double d, Matrix B)
@SHORT:Inverse of a matrix.
@ARG:double & d & The determinant of the matrix.
@ARG:@>Matrix & B & The inverse of the matrix.
This method returns the inverse of a $2 \times 2$ matrix.
The result of this operation is a $2 \times 2$ matrix defined by:
\begin{equation*}
\B = \A^{-1} = \frac {1}{d} \left[\begin{array}{cc}
  A_{22} & - A_{12}\\
  -A_{21} & A_{11}
  \end{array}
  \right]
\end{equation*}
where $\A$ is a matrix defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
void Matrix::computeInverse2x2(double det, Matrix &inverse) const
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  if ((_rows != 2) || (_cols != 2))
    fatalError("Matrix::det",
               "Your matrix is not a 2x2 square matrix, it's a [%d,%d] matrix",
               _rows, _cols);
  if ((inverse._rows != 2) || (inverse._cols != 2))
    fatalError("Matrix::det",
               "The inverse matrix is not a 2x2 square matrix, it's a [%d,%d] matrix",
               inverse._rows, inverse._cols);
#endif

  inverse._data[0] = _data[3] / det;
  inverse._data[1] = -_data[1] / det;
  inverse._data[2] = -_data[2] / det;
  inverse._data[3] = _data[0] / det;
}

/*
@LABEL:Matrix::computeInverse3x3(double d, Matrix B)
@SHORT:Inverse of a matrix.
@ARG:double & d & The determinant of the matrix.
@ARG:@>Matrix & B & The inverse of the matrix.
This method returns the inverse of a $3 \times 3$ matrix.
The result of this operation is a $3 \times 3$ matrix defined by:
\begin{equation*}
\B = \A^{-1} = \frac {1}{d} \left[\begin{array}{ccc}
  A_{22}A_{33}-A_{23}A_{32}&A_{13}A_{32}-A_{12}A_{33}&A_{12}A_{23}-A_{13}A_{22}\\
  A_{23}A_{31}-A_{21}A_{33}&A_{11}A_{33}-A_{13}A_{31}&A_{13}A_{21}-A_{11}A_{23}\\
  A_{21}A_{32}-A_{22}A_{31}&A_{12}A_{31}-A_{11}A_{32}&A_{11}A_{22}-A_{12}A_{21}
  \end{array}
  \right]
\end{equation*}
where $\A$ is a matrix defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
void Matrix::computeInverse3x3(double det, Matrix &inverse) const
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  if ((_rows != 3) || (_cols != 3))
    fatalError("Matrix::det",
               "Your matrix is not a 3x3 square matrix, it's a [%d,%d] matrix",
               _rows, _cols);
  if ((inverse._rows != 3) || (inverse._cols != 3))
    fatalError("Matrix::det",
               "The inverse matrix is not a 3x3 square matrix, it's a [%d,%d] matrix",
               inverse._rows, inverse._cols);

#endif

  inverse._data[0] = (_data[4] * _data[8] - _data[5] * _data[7]) / det;
  inverse._data[1] = (_data[2] * _data[7] - _data[1] * _data[8]) / det;
  inverse._data[2] = (_data[1] * _data[5] - _data[2] * _data[4]) / det;
  inverse._data[3] = (_data[5] * _data[6] - _data[3] * _data[8]) / det;
  inverse._data[4] = (_data[0] * _data[8] - _data[2] * _data[6]) / det;
  inverse._data[5] = (_data[2] * _data[3] - _data[0] * _data[5]) / det;
  inverse._data[6] = (_data[3] * _data[7] - _data[4] * _data[6]) / det;
  inverse._data[7] = (_data[1] * _data[6] - _data[0] * _data[7]) / det;
  inverse._data[8] = (_data[0] * _data[4] - _data[1] * _data[3]) / det;
}

/*
@LABEL:Matrix::inverse()
@SHORT:Inverse of a matrix.
@RETURN:Matrix & The inverse of the matrix.
This method returns the inverse of a matrix.
This method uses the Lapack \textsf{dgetrf} and \textsf{dgetri} Fortran subroutines to perform the operation.
@END
*/
//-----------------------------------------------------------------------------
Matrix Matrix::inverse() const
//-----------------------------------------------------------------------------
{
  Matrix inv(_rows, _cols);
  lapack_int n, m, info, lda;
  lapack_int *ipiv;
  n = _rows;
  m = _rows;
  lda = _rows;
  ipiv = new lapack_int[_rows];

  // copie de la matrice
  inv = *this;

  info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, inv._data, lda, ipiv);
  if (info != 0)
  {
    if (info < 0)
    {
      internalFatalErrorLine("Matrix::inverse",
                             "parameter (%d) of function dgetrf_ has illegal value",
                             -info);
    }
    if (info > 0)
    {
      fatalError("Matrix::inverse",
                 "pivot (%d) is exactly zero\n"
                 "The factorization has been completed, but this factor is exactly singular, and division by zero will occur if it is used to solve a system of equations",
                 info);
    }
  }

  info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, inv._data, lda, ipiv);
  if (info != 0)
  {
    if (info < 0)
    {
      internalFatalErrorLine("Matrix::inverse",
                             "parameter %d of function dgetri_ has illegal value",
                             -info);
    }
  }

  delete[] ipiv;

  // and return the result
  return inv;
}

/*
@LABEL:Matrix::getSolve(Vector x)
@SHORT:Solves a linear system $\A\cdot \overrightarrow{x} = \overrightarrow{b}$.
@RETURN:Vector & The solution of the linear system.
This method returns the solution of a linear system with the following form:
\begin{equation*}
\overrightarrow{y} = \A \cdot \overrightarrow{x}
\end{equation*}
where $\A$ is a matrix defined by the object itself and $\overrightarrow{x}$ is a vector defined by parameter x.
This method uses the Lapack \textsf{dgesv} Fortran subroutine to perform the operation.
@END
*/
//-----------------------------------------------------------------------------
Vector Matrix::getSolve(const Vector &vect) const
//-----------------------------------------------------------------------------
{
  Vector result(_rows);
  lapack_int info;
  lapack_int *ipiv;

  ipiv = new lapack_int[_rows];
  result = vect;

  // init a matrix
  double *A = new double[_dataLength];
  memcpy(A, _data, _dataLength * sizeof(double));

  info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, _rows, 1, A, _rows, ipiv, result._data, 1);
  if (info != 0)
  {
    if (info < 0)
    {
      internalFatalErrorLine("Matrix::getSolve",
                             "parameter %d of function dgesv_ has illegal value",
                             -info);
    }
    if (info > 0)
    {
      fatalError("Matrix::getSolve",
                 "pivot (%d) is exactly zero\n"
                 "The factorization has been completed, but this factor is exactly singular, and division by zero will occur if it is used to solve a system of equations",
                 info);
    }
  }

  delete[] ipiv;
  delete[] A;

  return result;
}

/*
@LABEL:Matrix::solve(Vector x)
@SHORT:Solves a linear system $\A\cdot \overrightarrow{x} = \overrightarrow{b}$.
@WARNING: The result of the operation is the parameter x itself.
This method returns the solution of a small linear system with the following form:
\begin{equation*}
\overrightarrow{y} = \A \cdot \overrightarrow{x}
\end{equation*}
where $\A$ is a matrix defined by the object itself and $\overrightarrow{x}$ is a vector defined by parameter x.
This method uses the Lapack \textsf{dgesv} Fortran subroutine to perform the operation.
@END
*/
//-----------------------------------------------------------------------------
void Matrix::solve(Vector &b)
//-----------------------------------------------------------------------------
{
  lapack_int info;
  lapack_int *ipiv = new lapack_int[_rows];

  info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, _rows, 1, _data, _rows, ipiv, b._data, 1);

  if (info != 0)
  {
    if (info < 0)
    {
      internalFatalErrorLine("Matrix::solveLinearPivoting",
                             "parameter %d of function dgesv_ has illegal value",
                             -info);
    }
    if (info > 0)
    {
      fatalError("Matrix::solveLinearPivoting",
                 "pivot (%d) is exactly zero\n"
                 "The factorization has been completed, but this factor is exactly singular, and division by zero will occur if it is used to solve a system of equations",
                 info);
    }
  }

  delete[] ipiv;
}

// evaluation du Null Space d'une matrice
/*
  Cette methode calcule le Null Space d'une matrice \f$ \textbf{A} \f$. Le null-space de  \f$ \textbf{A}=N.S.(\textbf{A})\subset\mathcal{W} \f$  est l'ensemble des elements de \f$ \mathcal{V} \f$ tels que \f$ \textbf{A}(\overrightarrow{_data})=\overrightarrow{0} \f$. Cette methode utilise la routine DGESVD de la librairie Lapack.

  Cette methode utilise des copies de vecteurs et matrices donnes en argument et ne modifie pas les valeurs contenues dans ceux-ci au cours de l'appel.
  - relative Parametre booleen indiquant si la tolerance de recherche est relative ou non. Dans le cas d'une tolerance relative, la tolerance de recherche de la plus patite valeur preopres est recalculee par rapport à la plus grande valeur propre.
  - tol tolerance de recherche du Null-Space par defaut 1e-10.
  Return : matrice contenant les vecteurs du Null-Space en cols (taille N x nNs) avce nNs nombre de "modes rigides"
*/
//-----------------------------------------------------------------------------
Matrix Matrix::getNullSpace(bool relative, double tol)
//-----------------------------------------------------------------------------
{
  Vector S(_rows);
  lapack_int n, m, info, lda;
  lapack_int ldu, ldvt;
  double *U = NULL;
  double *VT = NULL;
  double *superb;

  n = _rows;
  m = _rows;
  lda = _rows;
  //  ldb = _rows;
  ldu = _rows;
  ldvt = _rows;

  // init a matrix
  Matrix A = *this;

  //  DGESVD computes the singular value decomposition (SVD) of a double
  //  M-by-N matrix A, optionally computing the left and/or right singular
  //  vectors. The SVD is written

  //       A = U * SIGMA * transpose(V)

  //  where SIGMA is an M-by-N matrix which is zero except for its
  //  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
  //  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
  //  are the singular values of A; they are double and non-negative, and
  //  are returned in descending order.  The first min(m,n) cols of
  //  U and V are the left and right singular vectors of A.

  //  Note that the routine returns V**T, not V.

  // computes the work length vector
  // lwork = -1;
  // dgesvd_("N", "O", &m, &n, &A(0,0), &lda, &S(0), U, &ldu, VT, &ldvt, &worksize, &lwork, &info);
  // lwork = (int)worksize;

  // allocate memory for the work buffer
  superb = new double[m];

  // compute eigenvalues and eigenvectors
  // dgesvd_("N", "O", &m, &n, &A(0,0), &lda, &S(0), U, &ldu, VT, &ldvt, work, &lwork, &info);
  info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'N', 'O', m, n, &A(0, 0), lda, &S(0), U, ldu, VT, ldvt, superb);

  if (info != 0)
  {
    if (info < 0)
    {
      internalFatalErrorLine("Matrix::getNullSpace",
                             "parameter %d of function dgesv_ has illegal value",
                             -info);
    }
    if (info > 0)
    {
      fatalError("Matrix::getNullSpace",
                 "Eigenvalues decomposition incorrect",
                 info);
    }
  }

  // on recalcule la tolerance relative en la multipliant par la plus grande valeur propre
  if (relative)
    tol *= S(0);

  // on compte le nombre de null cols
  long ic = 0;
  for (long i = n; i > 0;)
    if (S(--i) <= tol)
      ic++;

  // dimension de la matrice de sortie
  Matrix ret(n, ic);

  // remplissage du Null-Space
  ic = 0;
  for (long i = 0; i < n; i++)
  {
    if (S(i) < tol)
    {
      // copy the last rows
      for (long j = 0; j < n; j++)
        ret(j, ic) = A(i, j);
      ic++;
    }
  }

  delete[] superb;

  return ret;
}

//-----------------------------------------------------------------------------
Matrix Matrix::getNullSpace2(bool relative, double tol)
//-----------------------------------------------------------------------------
{
  Vector S(_rows);
  lapack_int n, m, info, lda;
  lapack_int ldu, ldvt;
  Matrix U;
  Matrix VT;

  n = _rows;
  m = _rows;
  lda = _rows;
  //  ldb = _rows;
  ldu = _rows;
  ldvt = _rows;

  // init a matrix
  Matrix A = *this;

  U.redim(ldu, m);
  VT.redim(ldvt, n);

  //  DGESVD computes the singular value decomposition (SVD) of a double
  //  M-by-N matrix A, optionally computing the left and/or right singular
  //  vectors. The SVD is written

  //       A = U * SIGMA * transpose(V)

  //  where SIGMA is an M-by-N matrix which is zero except for its
  //  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
  //  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
  //  are the singular values of A; they are double and non-negative, and
  //  are returned in descending order.  The first min(m,n) cols of
  //  U and V are the left and right singular vectors of A.

  //  Note that the routine returns V**T, not V.

  // computes the work length vector
  // iwork = new int[8 * m];

  // lwork = -1;
  // dgesdd_("A", &m, &n, &A(0, 0), &lda, &S(0), &U(0, 0), &ldu, &VT(0, 0), &ldvt, &worksize, &lwork, iwork, &info);
  // lwork = (int)worksize;

  // allocate memory for the work buffer
  // work = new double[lwork];

  // compute eigenvalues and eigenvectors
  // dgesdd_("A", &m, &n, &A(0, 0), &lda, &S(0), &U(0, 0), &ldu, &VT(0, 0), &ldvt, work, &lwork, iwork, &info);
  info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'A', m, n, &A(0, 0), lda, &S(0), &U(0, 0), ldu, &VT(0, 0), ldvt);

  if (info != 0)
  {
    if (info < 0)
    {
      internalFatalErrorLine("Matrix::getNullSpace",
                             "parameter %d of function dgesdd_ has illegal value",
                             -info);
    }
    if (info > 0)
    {
      fatalError("Matrix::getNullSpace",
                 "Eigenvalues decomposition incorrect",
                 info);
    }
  }

  // on recalcule la tolerance relative en la multipliant par la plus grande valeur propre
  if (relative)
    tol *= S(0);

  // on compte le nombre de null cols
  long ic = 0;
  for (long i = n; i > 0;)
    if (S(--i) <= tol)
      ic++;

  // dimension de la matrice de sortie
  Matrix ret(n, ic);

  // remplissage du Null-Space
  ic = 0;
  for (long i = 0; i < n; i++)
  {
    if (S(i) < tol)
    {
      // copy the last rows
      for (long j = 0; j < n; j++)
        ret(j, ic) = VT(i, j);
      ic++;
    }
  }

  return ret;
}

/*
@LABEL:Matrix::computeSVD(Vector w, Matrix L, Matrix R)
@SHORT:Eigenvalues and eigenvactors of a $n \times m$ matrix.
@ARG:Vector & w & Eigenvalues of the matrix
@ARG:Matrix & L & Left eigenvector of the matrix
@ARG:Matrix & T & Right eigenvector of the matrix
This method computes the eigenvalues and the left and right eigenvectors of a $n \times m$ matrix
This method uses the Lapack \textsf{dgesvd} Fortran subroutine to perform the operation.
@END
*/
//-----------------------------------------------------------------------------
void Matrix::computeSVD(Vector &eigenValues, Matrix &leftEigenVectors, Matrix &rightEigenVectors)
//-----------------------------------------------------------------------------
{
  lapack_int info, lwork;
  double *superb;

  Matrix A = *this;
  leftEigenVectors.redim(_rows, _rows);
  rightEigenVectors.redim(_cols, _cols);
  eigenValues.redim(_cols);

  superb = new double[_rows];

  info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', _rows, _cols, &A(0, 0), _cols, &eigenValues(0), &leftEigenVectors(0, 0), _rows, &rightEigenVectors(0, 0), _cols, superb);

  if (info != 0)
  {
    if (info < 0)
    {
      internalFatalErrorLine("Matrix::computeSVD",
                             "parameter %d of function LAPACKE_dgesvd has illegal value",
                             -info);
    }
    if (info > 0)
    {
      fatalError("Matrix::computeSVD",
                 "SVD incorrect",
                 info);
    }
  }

  delete[] superb;
}

/* //-----------------------------------------------------------------------------
void Matrix::computeEigenVectors2(Vector &eigenValues, Matrix &eigenVectors)
//-----------------------------------------------------------------------------
{
  lapack_int n, m, info, lda;
  lapack_int ldu, ldvt, lwork;

  n = _rows;
  m = _rows;
  lda = _rows;
  ldu = _rows;
  ldvt = _rows;

  // init a matrix
  eigenVectors = *this;

  // redim the matrices and vectors
  //  leftEigenVectors.redim(n, n);
  //  rightEigenVectors.redim(n, n);
  eigenValues.redim(n);

  //  DGESVD computes the singular value decomposition (SVD) of a double
  //  M-by-N matrix A, optionally computing the left and/or right singular
  //  vectors. The SVD is written

  //       A = U * SIGMA * transpose(V)

  //  where SIGMA is an M-by-N matrix which is zero except for its
  //  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
  //  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
  //  are the singular values of A; they are double and non-negative, and
  //  are returned in descending order.  The first min(m,n) cols of
  //  U and V are the left and right singular vectors of A.

  //  Note that the routine returns V**T, not V.

  // computes the work length vector
  // lwork = -1;
  // dgesvd_("A", "A", &m, &n, &A(0, 0), &lda, &eigenValues(0), &leftEigenVectors(0, 0), &ldu, &rightEigenVectors(0, 0), &ldvt, &worksize, &lwork, &info);
  // lwork = (int)worksize;

  // allocate memory for the work buffer
  //  work = new double[lwork];
  // superb = new double[m];

  // compute eigenvalues and eigenvectors
  // dgesvd_("A", "A", &m, &n, &A(0, 0), &lda, &eigenValues(0), &leftEigenVectors(0, 0), &ldu, &rightEigenVectors(0, 0), &ldvt, work, &lwork, &info);
  // info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', m, n, &A(0, 0), lda, &eigenValues(0), &leftEigenVectors(0, 0), ldu, &rightEigenVectors(0, 0), ldvt, superb);
  info = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', n, &eigenVectors(0, 0), lda, &eigenValues(0));

  if (info != 0)
  {
    if (info < 0)
    {
      internalFatalErrorLine("Matrix::computeSVD",
                             "parameter %d of function dgesv_ has illegal value",
                             -info);
    }
    if (info > 0)
    {
      fatalError("Matrix::computeSVD",
                 "Eigenvalues decomposition incorrect",
                 info);
    }
  }
} */

/* //-----------------------------------------------------------------------------
void Matrix::computeEigenVectors2(Vector &eigenValues)
//-----------------------------------------------------------------------------
{
  lapack_int n, m, info, lda;
  lapack_int ldu, ldvt, lwork;

  n = _rows;
  m = _rows;
  lda = _rows;
  ldu = _rows;
  ldvt = _rows;

  // init a matrix
  Matrix A = *this;

  // redim the matrices and vectors
  //  leftEigenVectors.redim(n, n);
  //  rightEigenVectors.redim(n, n);
  eigenValues.redim(n);

  //  DGESVD computes the singular value decomposition (SVD) of a double
  //  M-by-N matrix A, optionally computing the left and/or right singular
  //  vectors. The SVD is written

  //       A = U * SIGMA * transpose(V)

  //  where SIGMA is an M-by-N matrix which is zero except for its
  //  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
  //  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
  //  are the singular values of A; they are double and non-negative, and
  //  are returned in descending order.  The first min(m,n) cols of
  //  U and V are the left and right singular vectors of A.

  //  Note that the routine returns V**T, not V.

  // computes the work length vector
  // lwork = -1;
  // dgesvd_("A", "A", &m, &n, &A(0, 0), &lda, &eigenValues(0), &leftEigenVectors(0, 0), &ldu, &rightEigenVectors(0, 0), &ldvt, &worksize, &lwork, &info);
  // lwork = (int)worksize;

  // allocate memory for the work buffer
  //  work = new double[lwork];
  // superb = new double[m];

  // compute eigenvalues and eigenvectors
  // dgesvd_("A", "A", &m, &n, &A(0, 0), &lda, &eigenValues(0), &leftEigenVectors(0, 0), &ldu, &rightEigenVectors(0, 0), &ldvt, work, &lwork, &info);
  // info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', m, n, &A(0, 0), lda, &eigenValues(0), &leftEigenVectors(0, 0), ldu, &rightEigenVectors(0, 0), ldvt, superb);
  info = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', n, &A(0, 0), lda, &eigenValues(0));

  if (info != 0)
  {
    if (info < 0)
    {
      internalFatalErrorLine("Matrix::computeSVD",
                             "parameter %d of function dgesv_ has illegal value",
                             -info);
    }
    if (info > 0)
    {
      fatalError("Matrix::computeSVD",
                 "Eigenvalues decomposition incorrect",
                 info);
    }
  }
} */

// evaluation de la pseudo inverse d'une matrice
/*
  Cette methode calcule la pseudoInverse d'une matrice \f$ \textbf{A} \f$. Cette methode utilise la routine DGESVD de la librairie Lapack.

  La pseudo inverse d'une matrice, egalement appelee inverse de Moore-Penrose est generalement notee par \f$ A^{+} \f$

  Cette methode utilise des copies de vecteurs et matrices donnes en argument et ne modifie pas les valeurs contenues dans ceux-ci au cours de l'appel. De plus, les matrices et vecteurs de retour donnes en argument n'on pas besoin d'etre dimensionnes à un taille correcte avant appel (ceci est fait en interne dans la routine).
  \warning Pour plus d'efficacite, cette methode est (Me semble-t-il !!!) restreinte au calcul des matrices carrees symetriques uniquement. Dans le cas d'une matrice non symetrique, il conviendra de verifier le comportement.
  - relative Parametre booleen indiquant si la tolerance de recherche est relative ou non. Dans le cas d'une tolerance relative, la tolerance de recherche de la plus patite valeur preopres est recalculee par rapport à la plus grande valeur propre.
  - tol tolerance de recherche du Null-Space par defaut 1e-10.
  Return : pseudo inverse de la matrice donnee en argument
*/
//-----------------------------------------------------------------------------
Matrix Matrix::getPseudoInverse(bool relative, double tol)
//-----------------------------------------------------------------------------
{
  Vector ev;
  Matrix U;
  Matrix V;
  Matrix result(_cols, _rows);
  long i;

  // compute all eigenvalues end eigenvectors
  computeSVD(ev, U, V);

  // on recalcule la tolerance relative en la multipliant par la plus grande valeur propre
  if (relative)
    tol *= ev(0);

  // calcul des inverses des valeurs propres (si dans la tolerance)
  for (i = 0; i < ev.size(); i++)
  {
    result(i, i) = (ev(i) >= tol ? 1.0 / ev(i) : 0.0);
  }

  // compute the pseudo inverse
  result = (V.transpose() * result * U.transpose());

  // retour
  return result;
}

// evaluation de la pseudo inverse d'une matrice et du Null-Space
/*
  Cette methode calcule la pseudoInverse et le Null-Space d'une matrice \f$ \textbf{A} \f$. Cette methode utilise la routine DGESVD de la librairie Lapack.

  La pseudo inverse d'une matrice, egalement appelee inverse de Moore-Penrose est generalement notee par \f$ A^{+} \f$.

  Cette methode utilise des copies de vecteurs et matrices donnes en argument et ne modifie pas les valeurs contenues dans ceux-ci au cours de l'appel. De plus, les matrices et vecteurs de retour donnes en argument n'on pas besoin d'etre dimensionnes à un taille correcte avant appel (ceci est fait en interne dans la routine).
  \warning Pour plus d'efficacite, cette methode est (Me semble-t-il !!!) restreinte au calcul des matrices carrees symetriques uniquement. Dans le cas d'une matrice non symetrique, il conviendra de verifier le comportement.
  - Kplus Matrice de retour pour la pseudo-inverse de la matrice donnee en argument
  - NS Matrice de retour pour le Null-Space de la matrice donnee en argument
  - relative Parametre booleen indiquant si la tolerance de recherche est relative ou non. Dans le cas d'une tolerance relative, la tolerance de recherche de la plus patite valeur preopres est recalculee par rapport à la plus grande valeur propre.
  - tol tolerance de recherche du Null-Space par defaut 1e-10.
*/
//-----------------------------------------------------------------------------
void Matrix::computePseudoInverse(Matrix &Kplus, Matrix &NS, bool relative, double tol)
//-----------------------------------------------------------------------------
{
  Vector ev;
  Matrix U;
  Matrix V;
  long nFloats;
  long i;

  // compute all ev end eigenvectors
  computeSVD(ev, U, V);

  // on recalcule la tolerance relative en la multipliant par la plus grande valeur propre
  if (relative)
    tol *= ev(0);

  // dimension de Kplus
  Kplus.redim(_cols, _rows);

  // init de nFloats
  nFloats = 0;

  // calcul des inverses des valeurs propres (si dans la tolerance)
  for (i = 0; i < ev.size(); i++)
  {
    if (ev(i) > tol)
    {
      Kplus(i, i) = 1.0 / ev(i);
    }
    else
    {
      Kplus(i, i) = 0.0;
      nFloats++;
    }
  }

  // compute the pseudo inverse
  Kplus = (V.transpose() * Kplus * U.transpose());

  // dimension de la matrice de Null-Space
  NS.redim(_rows, nFloats);

  // remplissage du Null-Space
  long ic = 0;
  for (long i = 0; i < ev.size(); i++)
  {
    if (ev(i) <= tol)
    {
      // copy the last rows
      for (long j = 0; j < _rows; j++)
        NS(j, ic) = V(i, j);
      ic++;
    }
  }
}

// Assemblage d'une matrice dans une autre
/*
  Cette methode effectue la phase d'assemblage d'une matrice dans une autre conformement à ce qui se fait classiquement en Elements Finis.

  - M Matrice elementaire à assembler dans la grande matrice
  - int tableau d'long correspondant au nombre de NOEUDS de l'element et contenant les Id des noeuds que l'on utilise
  - numberOfDim nombre de dimensions du probleme
  \code
  Matrix stiffnessMatrix;
  long numberOfDim=...; // nombre de dimensions du probleme
  long *ind=new long[loop_I]; // allocation memoire pour ind
  for (long Iloop=0; Iloop < loop_I; Iloop++) ind[Iloop]=pel->nodes (Iloop)->Id; // affectation des indices
  K.gatherFrom(stiffnessMatrix, ind, numberOfDim);
  delete [] ind;
  \endcode

  \version 1.1.0
*/
//-----------------------------------------------------------------------------
void Matrix::gatherFrom(const Matrix &M, long *ind0, int numberOfDim)
//-----------------------------------------------------------------------------
{
  //   long mr=M.rows();
  //   long mc=M.cols();
  long mr = M._rows;
  long mc = M._cols;
  long loop_I = mc / numberOfDim;
  long loop_J = mr / numberOfDim;
  long *ind = new long[loop_I];

  // preparation
  long *pind = ind;
  long *pind0 = ind0;
  for (long Iloop = loop_I; Iloop > 0; Iloop--)
    *pind++ = *pind0++ * numberOfDim;

  // boucle principale
  for (long Iloop = loop_I - 1; Iloop >= 0; Iloop--)
  {
    for (long Jloop = loop_J - 1; Jloop >= 0; Jloop--)
    {
      switch (numberOfDim)
      {
      case 1:
      {
        double *gl = &_data[Ind(ind[Iloop], ind[Jloop], _rows, _cols)];
        double *lo = &M._data[Ind(Iloop * numberOfDim, Jloop * numberOfDim, mr, mc)];
        *gl += *lo;
      }
      break;
      case 2:
      {
        double *gl = &_data[Ind(ind[Iloop], ind[Jloop], _rows, _cols)];
        double *lo = &M._data[Ind(Iloop * numberOfDim, Jloop * numberOfDim, mr, mc)];
        *gl++ += *lo++;
        *gl += *lo;
        gl += (_rows);
        lo += (mr);
        *gl-- += *lo--;
        *gl += *lo;
      }
      break;
      case 3:
      {
        double *gl = &_data[Ind(ind[Iloop], ind[Jloop], _rows, _cols)];
        double *lo = &M._data[Ind(Iloop * numberOfDim, Jloop * numberOfDim, mr, mc)];
        *gl++ += *lo++;
        *gl++ += *lo++;
        *gl += *lo;
        gl += (_rows);
        lo += (mr);
        *gl-- += *lo--;
        *gl-- += *lo--;
        *gl += *lo;
        gl += (_rows);
        lo += (mr);
        *gl++ += *lo++;
        *gl++ += *lo++;
        *gl += *lo;
      }
      break;

        // cas general
      default:
      {
        // boucle sur les ddl
        long i = Iloop * numberOfDim;
        long j = Jloop * numberOfDim;
        long globI = ind[Iloop];
        long globJ = ind[Jloop];
        for (long Ii = 0; Ii < numberOfDim; Ii++)
          for (long Jj = 0; Jj < numberOfDim; Jj++)
            (*this)(globI + Ii, globJ + Jj) += M(i + Ii, j + Jj);
      }
      }
    }
  }
  delete[] ind;

#ifdef VERIF_assert
  assert(ind == NULL);
#endif
}

//-----------------------------------------------------------------------------
void Matrix::scatterFrom(const Matrix &M, long *ind0, int numberOfDim)
//-----------------------------------------------------------------------------
{
  //   long mr=M.rows();
  //   long mc=M.cols();
  long mr = _rows;
  long mc = _cols;
  long loop_I = mc / numberOfDim;
  long loop_J = mr / numberOfDim;
  long *ind = new long[loop_I];

  // preparation
  long *pind = ind;
  long *pind0 = ind0;
  for (long Iloop = loop_I; Iloop > 0; Iloop--)
    *pind++ = *pind0++ * numberOfDim;

  // boucle principale
  for (long Iloop = loop_I - 1; Iloop >= 0; Iloop--)
  {
    for (long Jloop = loop_J - 1; Jloop >= 0; Jloop--)
    {
      switch (numberOfDim)
      {
      case 1:
      {
        double *lo = &M._data[Ind(ind[Iloop], ind[Jloop], _rows, _cols)];
        double *gl = &_data[Ind(Iloop * numberOfDim, Jloop * numberOfDim, mr, mc)];
        *gl += *lo;
      }
      break;
      case 2:
      {
        double *lo = &M._data[Ind(ind[Iloop], ind[Jloop], _rows, _cols)];
        double *gl = &_data[Ind(Iloop * numberOfDim, Jloop * numberOfDim, mr, mc)];
        *gl++ += *lo++;
        *gl += *lo;
        gl += (_rows);
        lo += (mr);
        *gl-- += *lo--;
        *gl += *lo;
      }
      break;
      case 3:
      {
        double *lo = &M._data[Ind(ind[Iloop], ind[Jloop], _rows, _cols)];
        double *gl = &_data[Ind(Iloop * numberOfDim, Jloop * numberOfDim, mr, mc)];
        *gl++ += *lo++;
        *gl++ += *lo++;
        *gl += *lo;
        gl += (_rows);
        lo += (mr);
        *gl-- += *lo--;
        *gl-- += *lo--;
        *gl += *lo;
        gl += (_rows);
        lo += (mr);
        *gl++ += *lo++;
        *gl++ += *lo++;
        *gl += *lo;
      }
      break;

        // cas general
      default:
      {
        // boucle sur les ddl
        long i = Iloop * numberOfDim;
        long j = Jloop * numberOfDim;
        long globI = ind[Iloop];
        long globJ = ind[Jloop];
        for (long Ii = 0; Ii < numberOfDim; Ii++)
          for (long Jj = 0; Jj < numberOfDim; Jj++)
            (*this)(i + Ii, j + Jj) += M(globI + Ii, globJ + Jj);
      }
      }
    }
  }
  delete[] ind;

#ifdef VERIF_assert
  assert(ind == NULL);
#endif
}

// Prints the content of the Matrix
//-----------------------------------------------------------------------------
void Matrix::printOut()
//-----------------------------------------------------------------------------
{
  switch (outType)
  {
  case outMatrixTypeNormal:
  {
    std::cout << "matrix " << _rows << "x" << _cols << " ={{";
    for (long i = 0; i < _rows; i++)
    {
      for (long j = 0; j < _cols; j++)
      {
        if (j != 0)
          std::cout << ",";
        std::cout << (*this)(i, j);
      }
      if (i != _rows - 1)
        std::cout << "},{";
    }
    std::cout << "}}\n";
  }
  break;
  case outMatrixTypeMupad:
  {
    std::cout << "matrix([[";
    for (long i = 0; i < _rows; i++)
    {
      for (long j = 0; j < _cols; j++)
      {
        if (j != 0)
          std::cout << ",";
        std::cout << _data[Ind(i, j, _rows, _cols)];
      }
      if (i != _rows - 1)
        std::cout << "],[";
    }
    std::cout << "]])\n";
  }
  break;
  default:
    fatalError("Matrices::printOut()", "unknown type %d\n", outType);
  }
}

//  Saves the content of a Matrix into a NumPy file
//-----------------------------------------------------------------------------
void Matrix::numpyWrite(std::string filename, bool initialize) const
//-----------------------------------------------------------------------------
{
  std::string mode = "a";
  if (initialize)
    mode = "w";
  NumpyInterface::npySave(filename, &_data[0], {_rows, _cols}, mode);
}

//  Saves the content of a Matrix into a NumPyZ file
//-----------------------------------------------------------------------------
void Matrix::numpyWriteZ(std::string filename, std::string name, bool initialize) const
//-----------------------------------------------------------------------------
{
  std::string mode = "a";
  if (initialize)
    mode = "w";
  NumpyInterface::npzSave(filename, name, &_data[0], {_rows, _cols}, mode);
}

//  Read the content of a Matrix from a NumPy file
//-----------------------------------------------------------------------------
void Matrix::numpyRead(std::string filename)
//-----------------------------------------------------------------------------
{
  NumpyInterface::NumpyArray arr = NumpyInterface::npyLoad(filename);

  redim(arr.shape[0], arr.shape[1]);

  if (arr.num_vals != _dataLength)
  {
    std::cout << "ERROR\n";
  }
  memcpy(_data, arr.data<double *>(), arr.num_vals * arr.word_size);
}

//  Read the content of a Matrix from a NumPyZ file
//-----------------------------------------------------------------------------
void Matrix::numpyReadZ(std::string filename, std::string name)
//-----------------------------------------------------------------------------
{
  NumpyInterface::NumpyArray arr = NumpyInterface::npzLoad(filename, name);

  redim(arr.shape[0], arr.shape[1]);

  if (arr.num_vals != _dataLength)
  {
    std::cout << "ERROR\n";
  }
  memcpy(_data, arr.data<double *>(), arr.num_vals * arr.word_size);
}
