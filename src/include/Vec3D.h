/**********************************************************************************
 *                                                                                *
 *  DynELA Finite Element Code v.4.0                                              *
 *  by Olivier PANTALE                                                            *
 *  Olivier.Pantale@enit.fr                                                       *
 *                                                                                *
 *********************************************************************************/
//@!CODEFILE = DynELA-H-file
//@!BEGIN = PRIVATE

#ifndef __dnlMaths_Vec3D_h__
#define __dnlMaths_Vec3D_h__

#include <string>
#include <string.h>
//#include <dnlKernel.h>
#include <cmath> //sqrt
#include "Macros.h"

class Tensor2;
class SymTensor2;

/*
@LABEL:Vec3D::Vec3D
@SHORT:Vec3D class.
This class is used to store information for coordinate vectors.
This file is the declaration file for the 3D vector class.
A 3D vector class is a vector with the following form:
\begin{equation*}
\overrightarrow{x}=\left[\begin{array}{c}
  x_{1}\\
  x_{2}\\
  x_{3}
  \end{array}\right]
\end{equation*}
where $x_{i}$ are the components of the vector $\overrightarrow{x}$.
@END
*/
class Vec3D
{
  friend class Tensor2;    // allows a direct access to private data for class Tensor2
  friend class SymTensor2; // allows a direct access to private data for class Tensor2

private:
  double _data[3]; // staticaly allocated vector storage
  bool indexOK(int) const;

public:
  Vec3D(double = 0, double = 0, double = 0);
  Vec3D(const Vec3D &);
  ~Vec3D();

  // Interface methods excluded from SWIG
#ifndef SWIG
  double &operator()(int);
  friend std::ifstream &operator>>(std::ifstream &, Vec3D &);
  friend std::ofstream &operator<<(std::ofstream &, const Vec3D &);
  friend std::ostream &operator<<(std::ostream &, const Vec3D &);
  friend Vec3D operator*(const double, const Vec3D &);
  Vec3D &operator=(const double *);
  Vec3D &operator=(const Vec3D &);
  Vec3D &operator=(double);
  Vec3D &read(std::ifstream &);
  void operator-=(const Vec3D &);
  void operator*=(const double);
  void operator/=(const double);
  void operator+=(const Vec3D &);
  void print(std::ostream &) const;
  void write(std::ofstream &) const;
#endif

  // Interface methods excluded from basic SWIG support
#if !defined(SWIG) || defined(CSWIG)
#endif

  bool isInsideBox(const Vec3D &, const Vec3D &) const;
  bool operator!=(const Vec3D &) const;
  bool operator==(const Vec3D &) const;
  double distance(const Vec3D &) const;
  double dot();
  double dot(const Vec3D &) const;
  double maxAbs();
  double maxVal();
  double minAbs();
  double minVal();
  double norm();
  double operator()(int) const;
  double squareDistance(const Vec3D &) const;
  double squareNorm();
  int size() const;
  SymTensor2 dyadic() const;
  Tensor2 dyadic(const Vec3D &) const;
  Tensor2 ewProduct(const Tensor2 &) const;
  Vec3D ewProduct(const Vec3D &) const;
  Vec3D getNormalized();
  Vec3D operator-() const;
  Vec3D operator-(const Vec3D &) const;
  Vec3D operator*(const double) const;
  Vec3D operator/(const double) const;
  Vec3D operator+(const Vec3D &) const;
  Vec3D vectorProduct(const Vec3D &) const;
  void normalize();
  void numpyRead(std::string);
  void numpyReadZ(std::string, std::string);
  void numpyWrite(std::string, bool = false) const;
  void numpyWriteZ(std::string, std::string, bool = false) const;
  void setNegativeValuesToZero();
  void setValue(double = 0);
  void setValue(double, double, double);
};

//------inline functions-------------------------------------------------------

//  tests if the index is ok
//-----------------------------------------------------------------------------
inline bool Vec3D::indexOK(int i) const
//-----------------------------------------------------------------------------
{
  // test
  if ((i >= 0) && (i < 3))
    return (true);

  // here, we detected an error in the index
  //fatalError("Vec3D::indexOK", "indice %d out of allowed range [0:3]", i);

  // stupid, because we'll never pass here because of the fatalError but it's to avoid a remark from the compiler
  return (false);
}

// Access to the values _data[i] of a 3D vector
//-----------------------------------------------------------------------------
inline double &Vec3D::operator()(int i)
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  indexOK(i);
#endif
  return _data[i];
}

// Access to the values _data[i] of 3D vector (Read only method)
//-----------------------------------------------------------------------------
inline double Vec3D::operator()(int i) const
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
  indexOK(i);
#endif
  return _data[i];
}

// Fill a 3D vector with a scalar value
//-----------------------------------------------------------------------------
inline void Vec3D::setValue(double val)
//-----------------------------------------------------------------------------
{
  _data[0] = _data[1] = _data[2] = val;
}

/*
@LABEL:Vec3D::setValue(double x, double y, double z)
@SHORT:Fill a vector with a scalar value.
@RETURN:Vec3D & The new Vec3D object created by the constructor.
@ARG:double & x & Component $x_{1}$ of the vector to create.
@ARG:double & y & Component $x_{2}$ of the vector to create.
@ARG:double & z & Component $x_{3}$ of the vector to create.
Fill a vector with a scalar value for the $3$ components of the vector.
\begin{equation*}
\overrightarrow{x}=\left[\begin{array}{c}
  x_{1} = x\\
  x_{2} = y\\
  x_{3} = z
  \end{array}\right]
\end{equation*}
@END
*/
//-----------------------------------------------------------------------------
inline void Vec3D::setValue(double x, double y, double z)
//-----------------------------------------------------------------------------
{
  _data[0] = x;
  _data[1] = y;
  _data[2] = z;
}

// Size of the 3D vector
//-----------------------------------------------------------------------------
inline int Vec3D::size() const
//-----------------------------------------------------------------------------
{
  return 3;
}

// Copy the content of a 3D vector into a new one
//-----------------------------------------------------------------------------
inline Vec3D &Vec3D::operator=(const Vec3D &vec)
//-----------------------------------------------------------------------------
{
  memcpy(_data, vec._data, 3 * sizeof(double));
  return *this;
}

/*
@LABEL:Vec3D::operator=(double v)
@SHORT:Fill a vector with a scalar value.
@RETURN:Vec3D
@ARG:double & v & Value to use for the operation.
This method is a surdefinition of the = operator for the vector class.
\begin{equation*}
\overrightarrow{x}=\left[\begin{array}{c}
v \\
v \\
v
\end{array}\right]
\end{equation*}
where $\overrightarrow{x}$ is a vector defined by the object itself and $v$ is the scalar value defined by parameter v.
@END
*/
//-----------------------------------------------------------------------------
inline Vec3D &Vec3D::operator=(double val)
//-----------------------------------------------------------------------------
{
  // recopie des valeurs
  _data[0] = _data[1] = _data[2] = val;

  // retour du vecteur
  return *this;
}

// Fill a 3D vector with a table of values
//-----------------------------------------------------------------------------
inline Vec3D &Vec3D::operator=(const double *vals)
//-----------------------------------------------------------------------------
{
  memcpy(_data, vals, 3 * sizeof(double));
  return *this;
}

/*
@LABEL:Vec3D::operator+(Vec3D y)
@SHORT:Addition of 2 vectors.
@ARG:Vec3D & y & Vector to add to the current one.
@RETURN:Vec3D & Result of the addition operation.
This method defines the addition of 2 vectors.
The result of this operation is also a vector defined by:
\begin{equation*}
\overrightarrow{z} = \overrightarrow{x} + \overrightarrow{y}
\end{equation*}
where $\overrightarrow{x}$ is a vector defined by the object itself and $\overrightarrow{y}$ is the vector value defined by parameter y.
@END
*/
//-----------------------------------------------------------------------------
inline Vec3D Vec3D::operator+(const Vec3D &vect) const
//-----------------------------------------------------------------------------
{
  return Vec3D(_data[0] + vect._data[0], _data[1] + vect._data[1], _data[2] + vect._data[2]);
}

/*
@LABEL:Vec3D::operator-(Vec3D y)
@SHORT:Difference of 2 vectors.
@ARG:Vec3D & y & Vector to add to the current one.
@RETURN:Vec3D & Result of the difference operation.
This method defines the difference of 2 vectors.
The result of this operation is also a vector defined by:
\begin{equation*}
\overrightarrow{z} = \overrightarrow{x} - \overrightarrow{y}
\end{equation*}
where $\overrightarrow{x}$ is a vector defined by the object itself and $\overrightarrow{y}$ is the vector value defined by parameter y.
@END
*/
//-----------------------------------------------------------------------------
inline Vec3D Vec3D::operator-(const Vec3D &vect) const
//-----------------------------------------------------------------------------
{
  return Vec3D(_data[0] - vect._data[0], _data[1] - vect._data[1], _data[2] - vect._data[2]);
}

/*
@LABEL:Vec3D::operator-(Vec3D x)
@SHORT:Opposite of a vector.
@RETURN:Vec3D & Result of the opposite operation.
This method defines the opposite of a vector.
The result of this operation is also a vector defined by:
\begin{equation*}
\overrightarrow{y} = -\overrightarrow{x}
\end{equation*}
where $\overrightarrow{x}$ is a vector defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
inline Vec3D Vec3D::operator-() const
//-----------------------------------------------------------------------------
{
  return Vec3D(-_data[0], -_data[1], -_data[2]);
}

/*
@LABEL:Vec3D::operator+=(Vec3D y)
@SHORT:Addition of 2 vectors.
@ARG:Vec3D & y & Vector to add to the current one.
@RETURN:Vec3D & Result of the addition operation.
This method defines the addition of 2 vectors.
The result of this operation is also a vector defined by:
\begin{equation*}
\overrightarrow{x} += \overrightarrow{y}
\end{equation*}
where $\overrightarrow{x}$ is a vector defined by the object itself and $\overrightarrow{y}$ is the vector value defined by parameter y.
@END
*/
//-----------------------------------------------------------------------------
inline void Vec3D::operator+=(const Vec3D &vect)
//-----------------------------------------------------------------------------
{
  _data[0] += vect._data[0];
  _data[1] += vect._data[1];
  _data[2] += vect._data[2];
}

/*
@LABEL:Vec3D::operator-=(Vec3D y)
@SHORT:Difference of 2 vectors.
@ARG:Vec3D & y & Vector to add to the current one.
@RETURN:Vec3D & Result of the difference operation.
This method defines the difference of 2 vectors.
The result of this operation is also a vector defined by:
\begin{equation*}
\overrightarrow{x} -= \overrightarrow{y}
\end{equation*}
where $\overrightarrow{x}$ is a vector defined by the object itself and $\overrightarrow{y}$ is the vector value defined by parameter y.
@END
*/
//-----------------------------------------------------------------------------
inline void Vec3D::operator-=(const Vec3D &vect)
//-----------------------------------------------------------------------------
{
  _data[0] -= vect._data[0];
  _data[1] -= vect._data[1];
  _data[2] -= vect._data[2];
}

/*
@LABEL:Vec3D::operator*=(double l)
@SHORT:Multiplication of a vector by a scalar.
@ARG:double & l & Scalar value to use for the operation.
@RETURN:Vec3D & Result of the multiplication operation.
This method defines the multiplication of a vector by a scalar value.
The result of this operation is also a vector defined by:
\begin{equation*}
\overrightarrow{x} *= \lambda
\end{equation*}
where $\overrightarrow{x}$ is a vector defined by the object itself and $\lambda$ is the vector value defined by parameter l.
@END
*/
//-----------------------------------------------------------------------------
inline void Vec3D::operator*=(const double lambda)
//-----------------------------------------------------------------------------
{
  _data[0] *= lambda;
  _data[1] *= lambda;
  _data[2] *= lambda;
}

/*
@LABEL:Vec3D::operator/=(double l)
@SHORT:Division of a vector by a scalar.
@ARG:double & l & Scalar value to use for the operation.
@RETURN:Vec3D & Result of the division operation.
This method defines the division of a vector by a scalar value.
The result of this operation is also a vector defined by:
\begin{equation*}
\overrightarrow{x} /= \lambda
\end{equation*}
where $\overrightarrow{x}$ is a vector defined by the object itself and $\lambda$ is the vector value defined by parameter l.
@END
*/
//-----------------------------------------------------------------------------
inline void Vec3D::operator/=(const double lambda)
//-----------------------------------------------------------------------------
{
  if (lambda == 0.0)
  {
//    fatalError("Vec3D:: operator /=", "divide by zero");
  }
  _data[0] /= lambda;
  _data[1] /= lambda;
  _data[2] /= lambda;
}

/*
@LABEL:Vec3D::J2()
@SHORT:J2 norm of a vector.
@RETURN:double & The J2 of the vector.
This method returns the J2 norm $s$ of a vector $\overrightarrow{x}$ defined by:
\begin{equation*}
s = \left\Vert \overrightarrow{x} \right\Vert  = \sqrt {x_{1}^2 + x_{2}^2 + x_{3}^2}
\end{equation*}
where $\overrightarrow{x}$ is a vector defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
inline double Vec3D::norm()
//-----------------------------------------------------------------------------
{
  return sqrt(dnlSquare(_data[0]) + dnlSquare(_data[1]) + dnlSquare(_data[2]));
}

/*
@LABEL:Vec3D::squareNorm()
@SHORT:Square norm of a vector.
@RETURN:double & The square of the vector.
This method returns the square norm $s$ of a vector $\overrightarrow{x}$ defined by:
\begin{equation*}
s = \left\Vert \overrightarrow{x} \right\Vert^2  = {x_{1}^2 + x_{2}^2 + x_{3}^2}
\end{equation*}
where $\overrightarrow{x}$ is a vector defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
inline double Vec3D::squareNorm()
//-----------------------------------------------------------------------------
{
  return dnlSquare(_data[0]) + dnlSquare(_data[1]) + dnlSquare(_data[2]);
}

/*
@LABEL:Vec3D::dot()
@SHORT:Returns the dot product of a Vec3D by itself.
@RETURN:double & The dot product.
This method returns dot product of a Vec3D by itself defined by:
\begin{equation*}
\left\Vert \overrightarrow{x} \right\Vert^2  = x_{1}^2 + x_{2}^2 + x_{3}^2
\end{equation*}
where $\overrightarrow{x}$ is a vector defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
inline double Vec3D::dot()
//-----------------------------------------------------------------------------
{
  return dnlSquare(_data[0]) + dnlSquare(_data[1]) + dnlSquare(_data[2]);
}

/*
@LABEL:Vec3D::normalize()
@SHORT:Normalize the Vec3D.
This method modifies the given vector and makes its norm equal to $1$.
@END
*/
//-----------------------------------------------------------------------------
inline void Vec3D::normalize()
//-----------------------------------------------------------------------------
{
  // calcul de la norme
  double _norm = norm();

  // recalcul des composantes
  _data[0] /= _norm;
  _data[1] /= _norm;
  _data[2] /= _norm;
}

/*
@LABEL:Vec3D::getNormalized()
@SHORT:Get the normalized Vec3D.
This method returns a colinear Vec3D with a norm equal to $1$.
@END
*/
//-----------------------------------------------------------------------------
inline Vec3D Vec3D::getNormalized()
//-----------------------------------------------------------------------------
{
  // calcul de la norme
  double _norm = norm();

  return Vec3D(_data[0] / _norm, _data[1] / _norm, _data[2] / _norm);
}

// Sets all negative components to zero
//-----------------------------------------------------------------------------
inline void Vec3D::setNegativeValuesToZero()
//-----------------------------------------------------------------------------
{
  if (_data[0] < 0)
    _data[0] = 0;
  if (_data[1] < 0)
    _data[1] = 0;
  if (_data[2] < 0)
    _data[2] = 0;
}

/*
@LABEL:Vec3D::distance(Vec3D y)
@SHORT:Distance between two points.
@RETURN:double & The distance.
This method computes the distance between two points using an Euclidian norm.
\begin{equation*}
d = \left\Vert \overrightarrow{y} - \overrightarrow{x} \right\Vert
\end{equation*}
where the $\overrightarrow{x}$ is the object itself and $\overrightarrow{y}$ is the second vector of the operation.
@END
*/
//-----------------------------------------------------------------------------
inline double Vec3D::distance(const Vec3D &vect) const
//-----------------------------------------------------------------------------
{
  Vec3D x = vect - *this;
  return x.norm();
}

/*
@LABEL:Vec3D::squareDistance(Vec3D y)
@SHORT:Square of distance between two points.
@RETURN:double & The distance.
This method computes the square of the distance between two points using an Euclidian norm.
\begin{equation*}
d = {\left\Vert \overrightarrow{y} - \overrightarrow{x} \right\Vert}^2
\end{equation*}
where the $\overrightarrow{x}$ is the object itself and $\overrightarrow{y}$ is the second vector of the operation.
@END
*/
//-----------------------------------------------------------------------------
inline double Vec3D::squareDistance(const Vec3D &vect) const
//-----------------------------------------------------------------------------
{
  Vec3D x = vect - *this;
  return x.dot();
}

#endif
