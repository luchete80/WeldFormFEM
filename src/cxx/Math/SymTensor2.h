/**********************************************************************************
 *                                                                                *
 *  DynELA Finite Element Code v.4.0                                              *
 *  by Olivier PANTALE                                                            *
 *  Olivier.Pantale@enit.fr                                                       *
 *                                                                                *
 *********************************************************************************/
//@!CODEFILE = DynELA-H-file
//@!BEGIN = PRIVATE

#ifndef __dnlMaths_SymTensor2_h__
#define __dnlMaths_SymTensor2_h__

#include <string>
#include <string.h>
#include <dnlKernel.h>
//#include <Macros.h>

class Vec3D;
class Tensor2;

/*
@LABEL:SymTensor2::SymTensor2
@SHORT:Second order tensor class.
The SymTensor2 class is used to store and manage data relative to a symmetric second order tensor which can be seen as a $3\times3$ matrix with the following form:
\begin{equation*}
\T=\left[\begin{array}{ccc}
  T_{11} & T_{12} & T_{13}\\
  T_{12} & T_{22} & T_{23}\\
  T_{13} & T_{23} & T_{33}
  \end{array}\right]
\end{equation*}
where $T_{ij}$ are the components of the symmetric second order tensor $\T$.
Concerning the internal storage of data, the SymTensor2 data is stored in vector of $6$ components named \textsf{\_data} using the following storage scheme:
\begin{equation*}
\T=\left[\begin{array}{ccc}
    T_{0} & T_{1} & T_{2}\\
    T_{1} & T_{3} & T_{4}\\
    T_{2} & T_{4} & T_{5}
    \end{array}\right]
\end{equation*}
@END
*/
class SymTensor2
{
  static const Tensor2Index _internalIndexes;
  friend class Vec3D;
  friend class Tensor2;

  double _data[6]; // Data storage for 6 double

private:
  bool indexOK(short, short) const;
  void buildFTF(double[3][3]) const;
  void polarExtract(double[3][3], double[3], SymTensor2 &, Tensor2 &) const;
  void polarExtractLnU(double[3][3], double[3], SymTensor2 &, Tensor2 &) const;

public:
  SymTensor2();
  SymTensor2(double, double, double, double, double, double);
  SymTensor2(const SymTensor2 &);
  ~SymTensor2();

  // Interface methods excluded from SWIG
#ifndef SWIG
  double &operator()(short, short);
  friend std::ifstream &operator>>(std::ifstream &, SymTensor2 &);
  friend std::ofstream &operator<<(std::ofstream &, const SymTensor2 &);
  friend std::ostream &operator<<(std::ostream &, const SymTensor2 &);
  friend SymTensor2 operator*(const double &, const SymTensor2 &);
  SymTensor2 &operator=(const double &);
  SymTensor2 &operator=(const SymTensor2 &);
  void operator-=(const SymTensor2 &);
  void operator*=(const double);
  void operator/=(const double);
  void operator+=(const SymTensor2 &);
  void print(std::ostream &) const;
  void read(std::ifstream &);
  void write(std::ofstream &) const;
#endif

  // Interface methods excluded from basic SWIG support
#if !defined(SWIG) || defined(CSWIG)
#endif

  bool operator!=(const SymTensor2 &) const;
  bool operator==(const SymTensor2 &) const;
  double doubleDot() const;
  double doubleDot(const SymTensor2) const;
  double det() const;
  double J2() const;
  double vonMises() const;
  double norm() const;
  double thirdTrace() const;
  double trace() const;
  double maxAbs();
  double maxVal();
  double minAbs();
  double minVal();
  double operator()(short, short) const;
  SymTensor2 cofactors() const;
  SymTensor2 deviator() const;
  SymTensor2 inverse() const;
  SymTensor2 minor() const;
  SymTensor2 operator-() const;
  SymTensor2 operator-(const SymTensor2 &) const;
  SymTensor2 operator*(const double) const;
  SymTensor2 operator/(const double) const;
  SymTensor2 operator+(const SymTensor2 &) const;
  SymTensor2 dotRTxR(const Tensor2) const;
  SymTensor2 dotRxRT(const Tensor2) const;
  SymTensor2 dot() const;
  Tensor2 operator*(const SymTensor2 &) const;
  Tensor2 operator*(const Tensor2 &) const;
  Tensor2 dot(const SymTensor2) const;
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
inline bool SymTensor2::indexOK(short i, short j) const
//-----------------------------------------------------------------------------
{
  if ((i >= 0) && (i < 3) && (j >= 0) && (j < 3))
  {
    return (true);
  }

  printf("SymTensor2::indexOK\nindice(s) [%d,%d] out of allowed range [0-2,0-2]", i, j);
  exit(-1);

  // no sense
  return (false);
}

//  Access to the values T[i,j] of a symmetric second order tensor
//-----------------------------------------------------------------------------
inline double &SymTensor2::operator()(short i, short j)
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  indexOK(i, j);
#endif
  return _data[_internalIndexes.index[i][j]];
}

//  Access to the values T[i,j] of a symmetric second order tensor (Read only method)
//-----------------------------------------------------------------------------
inline double SymTensor2::operator()(short i, short j) const
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  indexOK(i, j);
#endif
  return _data[_internalIndexes.index[i][j]];
}

//  Fill a symmetric second order tensor with a scalar value
//-----------------------------------------------------------------------------
inline void SymTensor2::setValue(const double val)
//-----------------------------------------------------------------------------
{
  _data[0] = val;
  _data[1] = val;
  _data[2] = val;
  _data[3] = val;
  _data[4] = val;
  _data[5] = val;
}

/*
@LABEL:SymTensor2::operator+=(SymTensor2 B)
@SHORT:Addition of 2 second order tensors.
@ARG:SymTensor2 & B & Second order tensor to add to the current one.
@RETURN:SymTensor2 & Result of the addition operation.
This method defines the addition of 2 second order tensors.
The result of this operation is also a symmetric second order tensor defined by:
\begin{equation*}
\A += \B
\end{equation*}
where $\A$ is a symmetric second order tensor defined by the object itself and $\B$ is the second order tensor value defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
inline void SymTensor2::operator+=(const SymTensor2 &tensor)
//-----------------------------------------------------------------------------
{
  // calcul de la somme
  _data[0] += tensor._data[0];
  _data[1] += tensor._data[1];
  _data[2] += tensor._data[2];
  _data[3] += tensor._data[3];
  _data[4] += tensor._data[4];
  _data[5] += tensor._data[5];
}

/*
@LABEL:SymTensor2::operator-=(SymTensor2 B)
@SHORT:Difference of 2 second order tensors.
@ARG:SymTensor2 & B & Second order tensor to add to the current one.
@RETURN:SymTensor2 & Result of the difference operation.
This method defines the difference of 2 second order tensors.
The result of this operation is also a symmetric second order tensor defined by:
\begin{equation*}
\A -= \B
\end{equation*}
where $\A$ is a symmetric second order tensor defined by the object itself and $\B$ is the second order tensor value defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
inline void SymTensor2::operator-=(const SymTensor2 &tensor)
//-----------------------------------------------------------------------------
{
  // calcul de la difference
  _data[0] -= tensor._data[0];
  _data[1] -= tensor._data[1];
  _data[2] -= tensor._data[2];
  _data[3] -= tensor._data[3];
  _data[4] -= tensor._data[4];
  _data[5] -= tensor._data[5];
}

/*
@LABEL:SymTensor2::operator*=(double l)
@SHORT:Multiplication of a symmetric second order tensor by a scalar.
@ARG:double & l & Scalar value to use for the operation.
@RETURN:SymTensor2 & Result of the multiplication operation.
This method defines the multiplication of a symmetric second order tensor by a scalar value.
The result of this operation is also a symmetric second order tensor defined by:
\begin{equation*}
\lambda \A
\end{equation*}
where $\A$ is a symmetric second order tensor defined by the object itself and $\lambda$ is the scalar value defined by parameter l.
@END
*/
//-----------------------------------------------------------------------------
inline void SymTensor2::operator*=(const double val)
//-----------------------------------------------------------------------------
{
  _data[0] *= val;
  _data[1] *= val;
  _data[2] *= val;
  _data[3] *= val;
  _data[4] *= val;
  _data[5] *= val;
}

/*
@LABEL:SymTensor2::operator/=(double l)
@SHORT:Division of a symmetric second order tensor by a scalar.
@ARG:double & l & Scalar value to use for the operation.
@RETURN:SymTensor2 & Result of the division operation.
This method defines the division of a symmetric second order tensor by a scalar value.
The result of this operation is also a symmetric second order tensor defined by:
\begin{equation*}
\frac{1}{\lambda} \A
\end{equation*}
where $\A$ is a symmetric second order tensor defined by the object itself and $\lambda$ is the scalar value defined by parameter l.
@END
*/
//-----------------------------------------------------------------------------
inline void SymTensor2::operator/=(const double val)
//-----------------------------------------------------------------------------
{
  _data[0] /= val;
  _data[1] /= val;
  _data[2] /= val;
  _data[3] /= val;
  _data[4] /= val;
  _data[5] /= val;
}

/*
@LABEL:SymTensor2::trace()
@SHORT:Returns the trace of a symmetric second order tensor.
@RETURN:double & The trace of the second order tensor.
This method returns the trace of a symmetric second order tensor, i.e. the sum $s$ of all the terms of the diagonal:
\begin{equation*}
s = \tr[\T] = T_{11}+T_{22}+T_{33}
\end{equation*}
where $\T$ is a symmetric second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
inline double SymTensor2::trace() const
//-----------------------------------------------------------------------------
{
  return (_data[0] + _data[3] + _data[5]);
}

/*
@LABEL:SymTensor2::thirdTrace()
@SHORT:Returns the average value of the trace of a symmetric second order tensor.
@RETURN:double & The third of the trace of the second order tensor.
This method returns average value of the trace of a symmetric second order tensor, i.e. the sum $s$ of all the terms of the diagonal divided by 3:
\begin{equation*}
s = \frac{1}{3} \tr[\T] =  \frac{1}{3} \left( T_{11}+T_{22}+T_{33} \right)
\end{equation*}
where $\T$ is a symmetric second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
inline double SymTensor2::thirdTrace() const
//-----------------------------------------------------------------------------
{
  return (_data[0] + _data[3] + _data[5]) / 3.0;
}

/*
@LABEL:SymTensor2::vonMises()
@SHORT:Returns the von Mises stress of a symmetric second order tensor.
@RETURN:double & the von Mises equivalent stress.
This method returns the von Mises stress of a symmetric second order tensor defined by:
\begin{equation*}
\overline{\sigma} = \frac {1}{\sqrt{2}}\sqrt{(s_{11}-s_{22})^2+(s_{22}-s_{33})^2+(s_{33}-s_{11})^2+6(s_{12}^2+s_{23}^2+s_{31}^2)}
\end{equation*}
where $\Sig$ is a symmetric second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
inline double SymTensor2::vonMises() const
//-----------------------------------------------------------------------------
{
  return (dnlVonMises(_data[0], _data[3], _data[5], _data[1], _data[2], _data[4]));
}

/*
@LABEL:SymTensor2::norm()
@SHORT:Norm of a symmetric second order tensor.
@RETURN:double & The norm of the second order tensor.
This method returns the norm $s$ of a symmetric second order tensor $\T$ defined by:\begin{equation*}
s = \left\Vert \T \right\Vert  = \sqrt {T_{ij}:T_{ij}}
\end{equation*}
where $\T$ is a symmetric second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
inline double SymTensor2::norm() const
//-----------------------------------------------------------------------------
{
  return sqrt(dnlSquare(_data[0]) + dnlSquare(_data[3]) + dnlSquare(_data[5]) +
              2.0 * (dnlSquare(_data[1]) + dnlSquare(_data[2]) + dnlSquare(_data[4])));
}

/*
@LABEL:SymTensor2::J2()
@SHORT:J2 norm of a symmetric second order tensor.
@RETURN:double
This method returns the J2 norm of a symmetric second order tensor defined by:
\begin{equation*}
\sqrt {\frac{3}{2}} \left\Vert s \right\Vert  = \sqrt {\frac{3}{2} s_{ij}:s_{ij}}
\end{equation*}
@END
*/
//-----------------------------------------------------------------------------
inline double SymTensor2::J2() const
//-----------------------------------------------------------------------------
{
  return dnlSqrt32 * sqrt(dnlSquare(_data[0]) + dnlSquare(_data[3]) + dnlSquare(_data[5]) +
                          2.0 * (dnlSquare(_data[1]) + dnlSquare(_data[2]) + dnlSquare(_data[4])));
}

/*
@LABEL:SymTensor2::det()
@SHORT:Determinant of a symmetric second order tensor.
@RETURN:double & The determinant.
This method returns the determinant of a symmetric second second order tensor.
The result of this operation is a scalar value defined by:
\begin{equation*}
D = T_{11} T_{22} T_{33} + 2 T_{12} T_{23} T_{13} - T_{13} T_{22} T_{13} - T_{11} T_{23} T_{23} - T_{12} T_{12} T_{33}
\end{equation*}
where $\T$ is a symmetric second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
inline double SymTensor2::det() const
//-----------------------------------------------------------------------------
{
  return (_data[0] * _data[3] * _data[5] +
          2.0 * _data[1] * _data[4] * _data[2] -
          _data[2] * _data[3] * _data[2] -
          _data[0] * _data[4] * _data[4] -
          _data[1] * _data[1] * _data[5]);
}

#endif
