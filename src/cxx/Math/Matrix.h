/**********************************************************************************
 *                                                                                *
 *  DynELA Finite Element Code v.4.0                                              *
 *  by Olivier PANTALE                                                            *
 *  Olivier.Pantale@enit.fr                                                       *
 *                                                                                *
 *********************************************************************************/
//@!CODEFILE = DynELA-H-file
//@!BEGIN = PRIVATE

#ifndef __dnlMaths_Matrix_h__
#define __dnlMaths_Matrix_h__

#include <Matrices.h>
#include <Vector.h>

class MatrixDiag;
class Tensor2;

/*
@LABEL:Matrix::Matrix
@SHORT:Matrix class.
This class is used to store information for matrices.
A Matrix class is a two dimensional object with size $n\times m$ with the following form:
\begin{equation*}
\M=\left[\begin{array}{cccc}
  M_{11} & M_{12} & \hdots & M_{1m}\\
  M_{21} & M_{22} & \hdots & M_{2m}\\
  \vdots & \vdots & \hdots & \vdots\\
  M_{n1} & M_{n1} & \hdots & M_{nm}
  \end{array}\right]
\end{equation*}
@END
*/
class Matrix : public Matrices
{
  friend class Vector;

  unsigned long _dataLength;
  double *_data;

  void allocate(const long, const long);
  void deallocate();

public:
  Matrix(const long = 3, const long = 3, const double = 0.0);
  Matrix(const Matrix &);
  ~Matrix();

  // Interface methods excluded from SWIG
#ifndef SWIG
  Matrix(int, int, double, double, ...);
  double &operator()(long, long);
  friend Matrix operator*(const double &, const Matrix &);
  friend std::ifstream &operator>>(std::ifstream &, Matrix &);
  friend std::ofstream &operator<<(std::ofstream &, const Matrix &);
  friend std::ostream &operator<<(std::ostream &, const Matrix &);
  Matrix &operator=(const double &);
  Matrix &operator=(const Matrix &);
  Matrix &operator=(const Tensor2 &);
  void operator-=(const Matrix &);
  void operator*=(const double &);
  void operator/=(const double &);
  void operator+=(const Matrix &);
  void print(std::ostream &) const;
  void read(std::ifstream &);
  void write(std::ofstream &) const;
#endif

  // Interface methods excluded from basic SWIG support
#if !defined(SWIG) || defined(CSWIG)
#endif

  bool indexOK(long, long) const;
  bool operator!=(const Matrix &) const;
  bool operator==(const Matrix &) const;
  double doubleDot() const;
  double doubleDot(const Matrix) const;
  double averageTrace() const;
  double det() const;
  double det2() const;
  double det3() const;
  double trace() const;
  double maxAbs() const;
  double maxVal() const;
  double getMean() const;
  double minAbs() const;
  double minVal() const;
  double operator()(long, long) const;
  double getSum() const;
  long Memory() const;
  Matrix cofactors() const;
  Matrix inverse() const;
  Matrix getNullSpace(bool = false, double = 1e-10);
  Matrix getNullSpace2(bool = false, double = 1e-10);
  Matrix getPseudoInverse(bool = false, double = 1e-10);
  Matrix skewSymmetric() const;
  Matrix symmetric() const;
  Matrix transpose() const;
  Matrix operator-() const;
  Matrix operator-(const Matrix &) const;
  Matrix operator*(const double &) const;
  Matrix operator*(const Matrix &) const;
  Matrix operator*(const MatrixDiag &) const;
  Matrix operator/(const double &) const;
  Matrix operator+(const Matrix &) const;
  Matrix dot() const;
  Matrix dot(const Matrix) const;
  Matrix dotNxT() const;
  Matrix dotTxN() const;
  Vector colSum() const;
  Vector col(long) const;
  Vector row(long) const;
  Vector getSolve(const Vector &) const;
  Vector operator*(const Vector &) const;
  Vector rowSum() const;
  Vector dotTxN(const Vector &) const;
  // void computeEigenVectors2(Vector &, Matrix &);
  // void computeEigenVectors2(Vector &);
  void computeInverse2x2(double, Matrix &) const;
  void computeInverse3x3(double, Matrix &) const;
  void computePseudoInverse(Matrix &, Matrix &, bool = false, double = 1e-10);
  void computeSVD(Vector &, Matrix &, Matrix &);
  void gatherFrom(const Matrix &, long *, int);
  void numpyRead(std::string);
  void numpyReadZ(std::string, std::string);
  void numpyWrite(std::string, bool = false) const;
  void numpyWriteZ(std::string, std::string, bool = false) const;
  void outStdout() { std::cout << *this << std::endl; }
  void printOut();
  void dot(Vector &) const;
  void redim(const long, const long);
  void scatterFrom(const Matrix &, long *, int);
  void setToUnity();
  void setValue(double);
  void solve(Vector &);
  void squareDivideBy(const MatrixDiag &);
  void squareMultiplyBy(const MatrixDiag &);
};

//------inline functions-------------------------------------------------------

//  tests if the index is ok
//-----------------------------------------------------------------------------
inline bool Matrix::indexOK(long i, long j) const
//-----------------------------------------------------------------------------
{
  if ((i >= 0) && (i < _rows) && (j >= 0) && (j < _cols))
    return (true);

  printf("Matrix::indexOK\nindice(s) [%ld,%ld] out of allowed range [0-%ld,0-%ld]", i, j, _rows - 1, _cols - 1);
  exit(-1);

  // fonction inutile, mais necessaire pour le compilateur
  return (false);
}

// Access to the values _data[i,j] of matrix (Read only method)
//-----------------------------------------------------------------------------
inline double Matrix::operator()(long i, long j) const
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  indexOK(i, j);
#endif
  return _data[Ind(i, j, _rows, _cols)];
}

// Access to the values _data[i,j] of matrix
//-----------------------------------------------------------------------------
inline double &Matrix::operator()(long i, long j)
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  indexOK(i, j);
#endif
  return _data[Ind(i, j, _rows, _cols)];
}

// Sum of all terms of the matrix
//-----------------------------------------------------------------------------
inline double Matrix::getSum() const
//-----------------------------------------------------------------------------
{
  double val = 0;
  for (long i = 0; i < _dataLength; i++)
    val += _data[i];
  return val;
}

// Mean of all terms of the matrix
//-----------------------------------------------------------------------------
inline double Matrix::getMean() const
//-----------------------------------------------------------------------------
{
  double val = 0;
  for (long i = 0; i < _dataLength; i++)
    val += _data[i];
  return val / (_dataLength);
}

#endif
