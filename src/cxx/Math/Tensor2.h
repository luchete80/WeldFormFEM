/**********************************************************************************
 *                                                                                *
 *  DynELA Finite Element Code v.4.0                                              *
 *  by Olivier PANTALE                                                            *
 *  Olivier.Pantale@enit.fr                                                       *
 *                                                                                *
 *********************************************************************************/
//@!CODEFILE = DynELA-H-file
//@!BEGIN = PRIVATE

#ifndef __dnlMaths_Tensor2_h__
#define __dnlMaths_Tensor2_h__

#include <string>
#include <string.h>
// #include <dnlKernel.h>

class Vec3D;
class SymTensor2;
class Matrix;

/*
@LABEL:Tensor2::Tensor2
@SHORT:Second order tensor class.
The Tensor2 class is used to store and manage data relative to a second order tensor which can be seen as a $3\times3$ matrix with the following form:
\begin{equation*}
\T=\left[\begin{array}{ccc}
  T_{11} & T_{12} & T_{13}\\
  T_{21} & T_{22} & T_{23}\\
  T_{31} & T_{32} & T_{33}
  \end{array}\right]
\end{equation*}
where $T_{ij}$ are the components of the second order tensor $\T$.
Concerning the internal storage of data, the Tensor2 data is stored in a vector of $9$ components named \textsf{\_data} using the following storage scheme:
\begin{equation*}
\T=\left[\begin{array}{ccc}
    T_{0} & T_{1} & T_{2}\\
    T_{3} & T_{4} & T_{5}\\
    T_{6} & T_{7} & T_{8}
    \end{array}\right]
\end{equation*}
@END
*/
class Tensor2
{
  static const Tensor2Index _internalIndexes;
  friend class Vec3D;
  friend class SymTensor2;
  friend class Matrix;

  double _data[9]; // Data storage for 9 double

private:
  bool indexOK(short, short) const;
  void buildFTF(double[3][3]) const;
  void polarExtract(double[3][3], double[3], SymTensor2 &, Tensor2 &) const;
  void polarExtractLnU(double[3][3], double[3], SymTensor2 &, Tensor2 &) const;

public:
  Tensor2();
  Tensor2(double, double, double, double, double, double, double, double, double);
  Tensor2(const Tensor2 &);
  ~Tensor2();

  // Interface methods excluded from SWIG
#ifndef SWIG
  double &operator()(short, short);
  friend std::ifstream &operator>>(std::ifstream &, Tensor2 &);
  friend std::ofstream &operator<<(std::ofstream &, const Tensor2 &);
  friend std::ostream &operator<<(std::ostream &, const Tensor2 &);
  friend Tensor2 operator*(const double &, const Tensor2 &);
  Tensor2 &operator=(const double &);
  Tensor2 &operator=(const SymTensor2 &);
  Tensor2 &operator=(const Tensor2 &);
  void operator-=(const Tensor2 &);
  void operator*=(const double);
  void operator/=(const double);
  void operator+=(const Tensor2 &);
  void print(std::ostream &) const;
  void read(std::ifstream &);
  void write(std::ofstream &) const;
#endif

  // Interface methods excluded from basic SWIG support
#if !defined(SWIG) || defined(CSWIG)
#endif

  bool operator!=(const Tensor2 &) const;
  bool operator==(const Tensor2 &) const;
  double doubleDot() const;
  double doubleDot(const Tensor2) const;
  double det() const;
  double J2() const;
  double norm() const;
  double thirdTrace() const;
  double trace() const;
  double maxAbs();
  double maxVal();
  double minAbs();
  double minVal();
  double operator()(short, short) const;
  SymTensor2 dotNxT() const;
  SymTensor2 dotTxN() const;
  Tensor2 cofactors() const;
  Tensor2 deviator() const;
  Tensor2 inverse() const;
  Tensor2 minor() const;
  Tensor2 skewSymmetric() const;
  Tensor2 symmetric() const;
  Tensor2 transpose() const;
  Tensor2 operator-() const;
  Tensor2 operator-(const Tensor2 &) const;
  Tensor2 operator*(const double) const;
  Tensor2 operator*(const SymTensor2 &) const;
  Tensor2 operator*(const Tensor2 &) const;
  Tensor2 operator/(const double) const;
  Tensor2 operator+(const Tensor2 &) const;
  Tensor2 dot() const;
  Tensor2 dot(const Tensor2) const;
  Vec3D colSum() const;
  Vec3D col(short) const;
  Vec3D row(short) const;
  Vec3D operator*(const Vec3D &) const;
  Vec3D rowSum() const;
  Vec3D solve(const Vec3D &) const;
  void numpyRead(std::string);
  void numpyReadZ(std::string, std::string);
  void numpyWrite(std::string, bool = false) const;
  void numpyWriteZ(std::string, std::string, bool = false) const;
  void polar(SymTensor2 &, Tensor2 &) const;
  void polarQL(SymTensor2 &, Tensor2 &) const;
  void polarJacobi(SymTensor2 &, Tensor2 &) const;
  void polarCuppen(SymTensor2 &, Tensor2 &) const;
  void polarLapack(SymTensor2 &, Tensor2 &) const;
  void polarLnU(SymTensor2 &, Tensor2 &) const;
  void polarQLLnU(SymTensor2 &, Tensor2 &) const;
  void polarJacobiLnU(SymTensor2 &, Tensor2 &) const;
  void polarCuppenLnU(SymTensor2 &, Tensor2 &) const;
  void polarLapackLnU(SymTensor2 &, Tensor2 &) const;
  void setToUnity();
  void setValue(const double);
  void setToZero();
};

//------inline functions-------------------------------------------------------

//  tests if the couple of indexes is ok
//-----------------------------------------------------------------------------
inline bool Tensor2::indexOK(short i, short j) const
//-----------------------------------------------------------------------------
{
  if ((i >= 0) && (i < 3) && (j >= 0) && (j < 3))
  {
    return (true);
  }

  printf("Tensor2::indexOK\nindice(s) [%d,%d] out of allowed range [0-2,0-2]", i, j);
  exit(-1);

  // no sense
  return (false);
}

//  Access to the values T[i,j] of a second order tensor
//-----------------------------------------------------------------------------
inline double &Tensor2::operator()(short i, short j)
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  indexOK(i, j);
#endif
  return _data[_internalIndexes.index[i][j]];
}

//  Access to the values T[i,j] of a second order tensor (Read only method)
//-----------------------------------------------------------------------------
inline double Tensor2::operator()(short i, short j) const
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  indexOK(i, j);
#endif
  return _data[_internalIndexes.index[i][j]];
}

//  Fill a second order tensor with a scalar value
//-----------------------------------------------------------------------------
inline void Tensor2::setValue(const double val)
//-----------------------------------------------------------------------------
{
  _data[0] = val;
  _data[1] = val;
  _data[2] = val;
  _data[3] = val;
  _data[4] = val;
  _data[5] = val;
  _data[6] = val;
  _data[7] = val;
  _data[8] = val;
}

/*
@LABEL:Tensor2::operator+=(Tensor2 B)
@SHORT:Addition of 2 second order tensors.
@ARG:Tensor2 & B & Second order tensor to add to the current one.
@RETURN:Tensor2 & Result of the addition operation.
This method defines the addition of 2 second order tensors.
The result of this operation is also a second order tensor defined by:
\begin{equation*}
\A += \B
\end{equation*}
where $\A$ is a second order tensor defined by the object itself and $\B$ is the second order tensor defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
inline void Tensor2::operator+=(const Tensor2 &T)
//-----------------------------------------------------------------------------
{
  // calcul de la somme
  _data[0] += T._data[0];
  _data[1] += T._data[1];
  _data[2] += T._data[2];
  _data[3] += T._data[3];
  _data[4] += T._data[4];
  _data[5] += T._data[5];
  _data[6] += T._data[6];
  _data[7] += T._data[7];
  _data[8] += T._data[8];
}

/*
@LABEL:Tensor2::operator-=(Tensor2 B)
@SHORT:Difference of 2 second order tensors.
@ARG:Tensor2 & B & Second order tensor to add to the current one.
@RETURN:Tensor2 & Result of the difference operation.
This method defines the difference of 2 second order tensors.
The result of this operation is also a second order tensor defined by:
\begin{equation*}
\A -= \B
\end{equation*}
where $\A$ is a second order tensor defined by the object itself and $\B$ is the second order tensor defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
inline void Tensor2::operator-=(const Tensor2 &T)
//-----------------------------------------------------------------------------
{
  // calcul de la difference
  _data[0] -= T._data[0];
  _data[1] -= T._data[1];
  _data[2] -= T._data[2];
  _data[3] -= T._data[3];
  _data[4] -= T._data[4];
  _data[5] -= T._data[5];
  _data[6] -= T._data[6];
  _data[7] -= T._data[7];
  _data[8] -= T._data[8];
}

/*
@LABEL:Tensor2::operator*=(double l)
@SHORT:Multiplication of a second order tensor by a scalar.
@ARG:double & l & Scalar value to use for the operation.
@RETURN:Tensor2 & Result of the multiplication operation.
This method defines the multiplication of a second order tensor by a scalar value.
The result of this operation is also a second order tensor defined by:
\begin{equation*}
\lambda \A
\end{equation*}
where $\A$ is a second order tensor defined by the object itself and $\lambda$ is the scalar value defined by parameter l.
@END
*/
//-----------------------------------------------------------------------------
inline void Tensor2::operator*=(const double val)
//-----------------------------------------------------------------------------
{
  _data[0] *= val;
  _data[1] *= val;
  _data[2] *= val;
  _data[3] *= val;
  _data[4] *= val;
  _data[5] *= val;
  _data[6] *= val;
  _data[7] *= val;
  _data[8] *= val;
}

/*
@LABEL:Tensor2::operator/=(double l)
@SHORT:Division of a second order tensor by a scalar.
@ARG:double & l & Scalar value to use for the operation.
@RETURN:Tensor2 & Result of the division operation.
This method defines the division of a second order tensor by a scalar value.
The result of this operation is also a second order tensor defined by:
\begin{equation*}
\frac{1}{\lambda} \A
\end{equation*}
where $\A$ is a second order tensor defined by the object itself and $\lambda$ is the scalar value defined by parameter l.
@END
*/
//-----------------------------------------------------------------------------
inline void Tensor2::operator/=(const double val)
//-----------------------------------------------------------------------------
{
  _data[0] /= val;
  _data[1] /= val;
  _data[2] /= val;
  _data[3] /= val;
  _data[4] /= val;
  _data[5] /= val;
  _data[6] /= val;
  _data[7] /= val;
  _data[8] /= val;
}

/*
@LABEL:Tensor2::trace()
@SHORT:Returns the trace of a second order tensor.
@RETURN:double & The trace of the second order tensor.
This method returns the trace of a second order tensor, i.e. the sum $s$ of all the terms of the diagonal:
\begin{equation*}
s = \tr[\T] = T_{11}+T_{22}+T_{33}
\end{equation*}
where $\T$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
inline double Tensor2::trace() const
//-----------------------------------------------------------------------------
{
  return (_data[0] + _data[4] + _data[8]);
}

/*
@LABEL:Tensor2::thirdTrace()
@SHORT:Returns the average value of the trace of a second order tensor.
@RETURN:double & The third of the trace of the second order tensor.
This method returns average value of the trace of a second order tensor, i.e. the sum $s$ of all the terms of the diagonal divided by 3:
\begin{equation*}
s = \frac{1}{3} \tr[\T] =  \frac{1}{3} \left( T_{11}+T_{22}+T_{33} \right)
\end{equation*}
where $\T$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
inline double Tensor2::thirdTrace() const
//-----------------------------------------------------------------------------
{
  return ((_data[0] + _data[4] + _data[8]) / 3.0);
}

/*
  Returns the von-Mises stress of a second order tensor

  This method returns the von-Mises stress of a second order tensor defined by:
  \f[ \overline{\sigma} = \frac {1}{\sqrt{2}}\sqrt{(s_{11}-s_{22})^2+(s_{22}-s_{33})^2+(s_{33}-s_{11})^2+6(s_{12}^2+s_{23}^2+s_{31}^2)}\f]
  Return : von-Mises stress of a second order tensor
*/
/* //-----------------------------------------------------------------------------
inline double Tensor2::vonMises() const
//-----------------------------------------------------------------------------
{
  exit(-1); // No sense, since this is not a symmetric tensor2
  return (dnlVonMises(_data[0], _data[4], _data[8], _data[1], _data[2], _data[5]));
} */

/*
@LABEL:Tensor2::norm()
@SHORT:Norm of a second order tensor.
@RETURN:double & The norm of the second order tensor.
This method returns the norm $s$ of a second order tensor $\T$ defined by:\begin{equation*}
s = \left\Vert \T \right\Vert  = \sqrt {T_{ij}:T_{ij}}
\end{equation*}
where $\T$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
inline double Tensor2::norm() const
//-----------------------------------------------------------------------------
{
  return sqrt(dnlSquare(_data[0]) + dnlSquare(_data[1]) + dnlSquare(_data[2]) +
              dnlSquare(_data[3]) + dnlSquare(_data[4]) + dnlSquare(_data[5]) +
              dnlSquare(_data[6]) + dnlSquare(_data[7]) + dnlSquare(_data[8]));
}

/*
@LABEL:Tensor2::J2()
@SHORT:J2 norm of a second order tensor.
@RETURN:double & The J2 of the second order tensor.
This method returns the J2 norm $s$ of a second order tensor $\T$ defined by:
\begin{equation*}
s = \sqrt {\frac{3}{2}} \left\Vert \T \right\Vert  = \sqrt {\frac{3}{2} T_{ij}:T_{ij}}
\end{equation*}
where $\T$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
inline double Tensor2::J2() const
//-----------------------------------------------------------------------------
{
  return dnlSqrt32 * sqrt(dnlSquare(_data[0]) + dnlSquare(_data[1]) + dnlSquare(_data[2]) +
                          dnlSquare(_data[3]) + dnlSquare(_data[4]) + dnlSquare(_data[5]) +
                          dnlSquare(_data[6]) + dnlSquare(_data[7]) + dnlSquare(_data[8]));
}

/*
@LABEL:Tensor2::det()
@SHORT:Determinant of a second order tensor.
@RETURN:double & The determinant.
This method returns the determinant of a second order tensor.
The result of this operation is a scalar value defined by:
\begin{equation*}
d = T_{11} T_{22} T_{33} + T_{21} T_{32} T_{13} + T_{31} T_{12} T_{23} - T_{31} T_{22} T_{13} - T_{11} T_{32} T_{23} - T_{21} T_{12} T_{33}
\end{equation*}
where $\T$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
inline double Tensor2::det() const
//-----------------------------------------------------------------------------
{
  return (_data[0] * _data[4] * _data[8] +
          _data[2] * _data[3] * _data[7] +
          _data[1] * _data[5] * _data[6] -
          _data[2] * _data[4] * _data[6] -
          _data[0] * _data[5] * _data[7] -
          _data[1] * _data[3] * _data[8]);
}

#endif
