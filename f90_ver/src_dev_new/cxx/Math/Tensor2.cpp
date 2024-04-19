/**********************************************************************************
 *                                                                                *
 *  DynELA Finite Element Code v.4.0                                              *
 *  by Olivier PANTALE                                                            *
 *  Olivier.Pantale@enit.fr                                                       *
 *                                                                                *
 *********************************************************************************/
//@!CODEFILE = DynELA-C-file
//@!BEGIN = PRIVATE

#include <fstream>
#include "Tensor2.h"
#include "SymTensor2.h"
#include "Vec3D.h"

// #include <Eigen3x3.h>
// #include <lapacke.h>

const Tensor2Index Tensor2::_internalIndexes = {0, 1, 2, 3, 4, 5, 6, 7, 8};

/*
@LABEL:Tensor2::Tensor2()
@SHORT:Default constructor of the Tensor2 class.
@RETURN:Tensor2 & The new Tensor2 object created by the constructor.
This is the default constructor of the Tensor2 class, where all components are initialized to zero by default.
\begin{equation*}
\T=\left[\begin{array}{ccc}
0 & 0 & 0\\
0 & 0 & 0\\
0 & 0 & 0
\end{array}\right]
\end{equation*}
@END
*/
//-----------------------------------------------------------------------------
Tensor2::Tensor2()
//-----------------------------------------------------------------------------
{
  setValue(0);
}

// Copy constructor
//-----------------------------------------------------------------------------
Tensor2::Tensor2(const Tensor2 &T)
//-----------------------------------------------------------------------------
{
  memcpy(_data, T._data, 9 * sizeof(double));
}

/*
@LABEL:Tensor2::Tensor2(double t11, double t12, ...)
@SHORT:Constructor of the Tensor2 class.
@RETURN:Tensor2 & The new Tensor2 object created by the constructor.
@ARG:double & t11 & Component $T_{11}$ of the second order tensor.
@ARG:double & t12 & Component $T_{12}$ of the second order tensor.
@ARG:double & t13 & Component $T_{13}$ of the second order tensor.
@ARG:double & t21 & Component $T_{21}$ of the second order tensor.
@ARG:double & t22 & Component $T_{22}$ of the second order tensor.
@ARG:double & t23 & Component $T_{23}$ of the second order tensor.
@ARG:double & t31 & Component $T_{31}$ of the second order tensor.
@ARG:double & t32 & Component $T_{32}$ of the second order tensor.
@ARG:double & t33 & Component $T_{33}$ of the second order tensor.
Constructor of a second order tensor with explicit initialization of the $9$ components of the tensor.
\begin{equation*}
\T=\left[\begin{array}{ccc}
  T_{11} & T_{12} & T_{13}\\
  T_{21} & T_{22} & T_{23}\\
  T_{31} & T_{32} & T_{33}
  \end{array}\right]
\end{equation*}
@END
*/
//-----------------------------------------------------------------------------
Tensor2::Tensor2(double t11, double t12, double t13, double t21, double t22, double t23, double t31, double t32, double t33)
//-----------------------------------------------------------------------------
{
  _data[0] = t11;
  _data[1] = t12;
  _data[2] = t13;
  _data[3] = t21;
  _data[4] = t22;
  _data[5] = t23;
  _data[6] = t31;
  _data[7] = t32;
  _data[8] = t33;
}

// Destructor
//-----------------------------------------------------------------------------
Tensor2::~Tensor2()
//-----------------------------------------------------------------------------
{
}

// Send the content of a second order tensor to the output flux for display
//-----------------------------------------------------------------------------
std::ostream &operator<<(std::ostream &outputStream, const Tensor2 &tensor)
//-----------------------------------------------------------------------------
{
  tensor.print(outputStream);
  return outputStream;
}

// Print the content of a second order tensor to the output flux for display
//-----------------------------------------------------------------------------
void Tensor2::print(std::ostream &outputStream) const
//-----------------------------------------------------------------------------
{
  short i, j;
  outputStream << "tensor 3x3 ={{";
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      if (j != 0)
      {
        outputStream << ",";
      }
      outputStream << _data[_internalIndexes.index[i][j]];
    }
    if (i != 2)
    {
      outputStream << "},{";
    }
  }
  outputStream << "}}";
}

/*
@LABEL:Tensor2::setToUnity()
@SHORT:Unity second order tensor.
@SMOD
This method transforms the current tensor to a unity tensor.
\begin{equation*}
\T=\left[\begin{array}{ccc}
1&0&0\\
0&1&0\\
0&0&1
\end{array}\right]
\end{equation*}
where $\T$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
void Tensor2::setToUnity()
//-----------------------------------------------------------------------------
{
  _data[0] = 1;
  _data[1] = 0;
  _data[2] = 0;
  _data[3] = 0;
  _data[4] = 1;
  _data[5] = 0;
  _data[6] = 0;
  _data[7] = 0;
  _data[8] = 1;
}

/*
@LABEL:Tensor2::setToZero()
@SHORT:Zero second order tensor.
@SMOD
This method transforms the current second order tensor to a zero tensor.
\begin{equation*}
\T=\left[\begin{array}{ccc}
0 & 0 & 0\\
0 & 0 & 0\\
0 & 0 & 0
\end{array}\right]
\end{equation*}
where $\T$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
void Tensor2::setToZero()
//-----------------------------------------------------------------------------
{
  setValue(0);
}

/*
@LABEL:Tensor2::operator=(double v)
@SHORT:Fill a second order tensor with a scalar value.
@RETURN:Tensor2
@ARG:double & v & Value to use for the operation.
This method is a surdefinition of the = operator for the second order tensor class.
\begin{equation*}
\T=\left[\begin{array}{ccc}
v & v & v\\
v & v & v\\
v & v & v
\end{array}\right]
\end{equation*}
where $\T$ is a second order tensor defined by the object itself and $v$ is the scalar value defined by parameter v.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 &Tensor2::operator=(const double &val)
//-----------------------------------------------------------------------------
{
  setValue(val);
  return *this;
}

// Copy the content of a second order tensor into a new one
//-----------------------------------------------------------------------------
Tensor2 &Tensor2::operator=(const Tensor2 &T)
//-----------------------------------------------------------------------------
{
  memcpy(_data, T._data, 9 * sizeof(double));
  return *this;
}

/*
@LABEL:Tensor2::operator=(SymTensor2 T)
@SHORT:Copy the content of a SymTensor2 into a Tensor2.
@RETURN:Tensor2
@ARG:SymTensor2 & T & Symmetric second order tensor to copy.
The result of this operation is a second order tensor as a copy of a symmetric second order tensor where $\T$ is a symmetric second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 &Tensor2::operator=(const SymTensor2 &T)
//-----------------------------------------------------------------------------
{
  _data[0] = T._data[0];
  _data[1] = T._data[1];
  _data[2] = T._data[2];
  _data[3] = T._data[1];
  _data[4] = T._data[3];
  _data[5] = T._data[4];
  _data[6] = T._data[2];
  _data[7] = T._data[4];
  _data[8] = T._data[5];
  return *this;
}

/*
@LABEL:Tensor2::operator+(Tensor2 B)
@SHORT:Addition of 2 second order tensors.
@ARG:Tensor2 & B & Second order tensor to add to the current one.
@RETURN:Tensor2 & Result of the addition operation.
This method defines the addition of 2 second order tensors.
The result of this operation is also a second order tensor defined by:
\begin{equation*}
\T = \A + \B
\end{equation*}
where $\A$ is a second order tensor defined by the object itself and $\B$ is the second order tensor defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 Tensor2::operator+(const Tensor2 &T) const
//-----------------------------------------------------------------------------
{
  // creation d'un nouveau tenseur
  Tensor2 result;

  // calcul de la somme
  result._data[0] = _data[0] + T._data[0];
  result._data[1] = _data[1] + T._data[1];
  result._data[2] = _data[2] + T._data[2];
  result._data[3] = _data[3] + T._data[3];
  result._data[4] = _data[4] + T._data[4];
  result._data[5] = _data[5] + T._data[5];
  result._data[6] = _data[6] + T._data[6];
  result._data[7] = _data[7] + T._data[7];
  result._data[8] = _data[8] + T._data[8];

  // renvoi de l'objet
  return result;
}

/*
@LABEL:Tensor2::operator-(Tensor2 B)
@SHORT:Difference of 2 second order tensors.
@ARG:Tensor2 & B & Second order tensor to subtract to the current one.
@RETURN:Tensor2 & Result of the difference operation.
This method defines the subtraction of 2 second order tensors.
The result of this operation is also a second order tensor defined by:
\begin{equation*}
\T = \A - \B
\end{equation*}
where $\A$ is a second order tensor defined by the object itself and $\B$ is the second order tensor defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 Tensor2::operator-(const Tensor2 &T) const
//-----------------------------------------------------------------------------
{
  // creation d'un nouveau tenseur
  Tensor2 result;

  // calcul de la somme
  result._data[0] = _data[0] - T._data[0];
  result._data[1] = _data[1] - T._data[1];
  result._data[2] = _data[2] - T._data[2];
  result._data[3] = _data[3] - T._data[3];
  result._data[4] = _data[4] - T._data[4];
  result._data[5] = _data[5] - T._data[5];
  result._data[6] = _data[6] - T._data[6];
  result._data[7] = _data[7] - T._data[7];
  result._data[8] = _data[8] - T._data[8];

  // renvoi de l'objet
  return result;
}

/*
@LABEL:Tensor2::operator-()
@SHORT:Opposite of a second order tensor.
@RETURN:Tensor2 & The opposite second order tensor.
This method defines the opposite of a second order tensor.
The result of this operation is also a second order tensor defined by:
\begin{equation*}
\T = - \A
\end{equation*}
where $\A$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 Tensor2::operator-() const
//-----------------------------------------------------------------------------
{
  Tensor2 result;

  result._data[0] = -_data[0];
  result._data[1] = -_data[1];
  result._data[2] = -_data[2];
  result._data[3] = -_data[3];
  result._data[4] = -_data[4];
  result._data[5] = -_data[5];
  result._data[6] = -_data[6];
  result._data[7] = -_data[7];
  result._data[8] = -_data[8];

  return result;
}

/*
@LABEL:Tensor2::operator*(double l)
@SHORT:Multiplication of a second order tensor by a scalar.
@ARG:double & l & Scalar value to use for the operation.
@RETURN:Tensor2 & Result of the multiplication operation.
This method defines the multiplication of a second order tensor by a scalar value.
The result of this operation is also a second order tensor defined by:
\begin{equation*}
\T = \lambda \A
\end{equation*}
where $\A$ is a second order tensor defined by the object itself and $\lambda$ is the scalar value defined by parameter l.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 Tensor2::operator*(const double l) const
//-----------------------------------------------------------------------------
{
  Tensor2 result;

  result._data[0] = l * _data[0];
  result._data[1] = l * _data[1];
  result._data[2] = l * _data[2];
  result._data[3] = l * _data[3];
  result._data[4] = l * _data[4];
  result._data[5] = l * _data[5];
  result._data[6] = l * _data[6];
  result._data[7] = l * _data[7];
  result._data[8] = l * _data[8];

  return result;
}

/*
@LABEL:Tensor2::operator/(double l)
@SHORT:Division of a second order tensor by a scalar.
@ARG:double & l & Scalar value to use for the operation.
@RETURN:Tensor2 & Result of the division operation.
This method defines the division of a second order tensor by a scalar value.
The result of this operation is also a second order tensor defined by:
\begin{equation*}
\T = \frac{1}{\lambda} \A
\end{equation*}
where $\A$ is a second order tensor defined by the object itself and $\lambda$ is the scalar value defined by parameter l.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 Tensor2::operator/(const double l) const
//-----------------------------------------------------------------------------
{
  Tensor2 result;

#ifdef VERIF_maths
  if (l == 0)
  {
    fatalError("Tensor2:: operator /", "divide by zero");
  }
#endif

  result._data[0] = _data[0] / l;
  result._data[1] = _data[1] / l;
  result._data[2] = _data[2] / l;
  result._data[3] = _data[3] / l;
  result._data[4] = _data[4] / l;
  result._data[5] = _data[5] / l;
  result._data[6] = _data[6] / l;
  result._data[7] = _data[7] / l;
  result._data[8] = _data[8] / l;

  return result;
}

/*
@LABEL:operator*(double l, Tensor2 A)
@SHORT:Multiplication of a second order tensor by a scalar.
@ARG:double & l & Scalar value to use for the operation.
@ARG:Tensor2 & A & Second order tensor to use for the operation.
@RETURN:Tensor2 & Result of the multiplication operation.
This method defines the multiplication of a second order tensor by a scalar value.
The result of this operation is also a second order tensor defined by:
\begin{equation*}
\T = \lambda \A
\end{equation*}
where $\A$ is a second order tensor and $\lambda$ is the scalar value defined by parameter l.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 operator*(const double &l, const Tensor2 &T)
//-----------------------------------------------------------------------------
{
  Tensor2 result;

  result._data[0] = l * T._data[0];
  result._data[1] = l * T._data[1];
  result._data[2] = l * T._data[2];
  result._data[3] = l * T._data[3];
  result._data[4] = l * T._data[4];
  result._data[5] = l * T._data[5];
  result._data[6] = l * T._data[6];
  result._data[7] = l * T._data[7];
  result._data[8] = l * T._data[8];

  return result;
}

/*
@LABEL:Tensor2::dot()
@SHORT:Single contracted product of a second order tensor by itself.
@RETURN:Tensor2 & Result of the multiplication operation.
This method defines a single contracted product of a second order tensor by itself.
The result of this operation is also a second order tensor defined by:
\begin{equation*}
\T = \A \cdot \A
\end{equation*}
where $\A$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 Tensor2::dot() const
//-----------------------------------------------------------------------------
{
  Tensor2 result;

  result._data[0] = _data[0] * _data[0] + _data[1] * _data[3] + _data[2] * _data[6];
  result._data[1] = _data[0] * _data[1] + _data[1] * _data[4] + _data[2] * _data[7];
  result._data[2] = _data[0] * _data[2] + _data[1] * _data[5] + _data[2] * _data[8];
  result._data[3] = _data[3] * _data[0] + _data[4] * _data[3] + _data[5] * _data[6];
  result._data[4] = _data[3] * _data[1] + _data[4] * _data[4] + _data[5] * _data[7];
  result._data[5] = _data[3] * _data[2] + _data[4] * _data[5] + _data[5] * _data[8];
  result._data[6] = _data[6] * _data[0] + _data[7] * _data[3] + _data[8] * _data[6];
  result._data[7] = _data[6] * _data[1] + _data[7] * _data[4] + _data[8] * _data[7];
  result._data[8] = _data[6] * _data[2] + _data[7] * _data[5] + _data[8] * _data[8];

  return result;
}

/*
@LABEL:Tensor2::dot(Tensor2 B)
@SHORT:Single contracted product of two second order tensors.
@RETURN:Tensor2 & Result of the multiplication operation.
@ARG:Tensor2 & B & Second tensor for the multiplication operation.
This method defines a single contracted product of two second order tensors.
The result of this operation is also a second order tensor defined by:
\begin{equation*}
\T = \A \cdot \B
\end{equation*}
where $\A$ is a second order tensor defined by the object itself and $\B$ is the second order tensor defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 Tensor2::dot(const Tensor2 T) const
//-----------------------------------------------------------------------------
{
  //  return (*this) * T;
  Tensor2 result;

  result._data[0] = _data[0] * T._data[0] + _data[1] * T._data[3] + _data[2] * T._data[6];
  result._data[1] = _data[0] * T._data[1] + _data[1] * T._data[4] + _data[2] * T._data[7];
  result._data[2] = _data[0] * T._data[2] + _data[1] * T._data[5] + _data[2] * T._data[8];
  result._data[3] = _data[3] * T._data[0] + _data[4] * T._data[3] + _data[5] * T._data[6];
  result._data[4] = _data[3] * T._data[1] + _data[4] * T._data[4] + _data[5] * T._data[7];
  result._data[5] = _data[3] * T._data[2] + _data[4] * T._data[5] + _data[5] * T._data[8];
  result._data[6] = _data[6] * T._data[0] + _data[7] * T._data[3] + _data[8] * T._data[6];
  result._data[7] = _data[6] * T._data[1] + _data[7] * T._data[4] + _data[8] * T._data[7];
  result._data[8] = _data[6] * T._data[2] + _data[7] * T._data[5] + _data[8] * T._data[8];

  return result;
}

/*
@LABEL:Tensor2::dotTxN()
@SHORT:Single contracted product of a second order tensor by its transpose.
@RETURN:SymTensor2 & Result of the multiplication operation.
This method defines a single contracted product of a second order tensor by its transpose.
The result of this operation is also a second order tensor defined by:
\begin{equation*}
\T = \A^T\cdot \A
\end{equation*}
where $\A$ is a second order tensor defined by the object itself. Result is a symmetric second order tensor.
@END
*/
//-----------------------------------------------------------------------------
SymTensor2 Tensor2::dotTxN() const
//-----------------------------------------------------------------------------
{
  SymTensor2 result;

  result._data[0] = _data[0] * _data[0] + _data[3] * _data[3] + _data[6] * _data[6];
  result._data[1] = _data[0] * _data[1] + _data[3] * _data[4] + _data[6] * _data[7];
  result._data[2] = _data[0] * _data[2] + _data[3] * _data[5] + _data[6] * _data[8];
  result._data[3] = _data[1] * _data[1] + _data[4] * _data[4] + _data[7] * _data[7];
  result._data[4] = _data[1] * _data[2] + _data[4] * _data[5] + _data[7] * _data[8];
  result._data[5] = _data[2] * _data[2] + _data[5] * _data[5] + _data[8] * _data[8];

  return result;
}

/*
@LABEL:Tensor2::dotNxT()
@SHORT:Single contracted product of a second order tensor by its transpose.
@RETURN:SymTensor2 & Result of the multiplication operation.
This method defines a single contracted product of a second order tensor by its transpose.
The result of this operation is also a second order tensor defined by:
\begin{equation*}
\T = \A \cdot \A^T
\end{equation*}
where $\A$ is a second order tensor defined by the object itself. Result is a symmetric second order tensor.
@END
*/
//-----------------------------------------------------------------------------
SymTensor2 Tensor2::dotNxT() const
//-----------------------------------------------------------------------------
{
  SymTensor2 result;

  result._data[0] = _data[0] * _data[0] + _data[1] * _data[1] + _data[2] * _data[2];
  result._data[1] = _data[0] * _data[3] + _data[1] * _data[4] + _data[2] * _data[5];
  result._data[2] = _data[0] * _data[6] + _data[1] * _data[7] + _data[2] * _data[8];
  result._data[3] = _data[3] * _data[3] + _data[4] * _data[4] + _data[5] * _data[5];
  result._data[4] = _data[3] * _data[6] + _data[4] * _data[7] + _data[5] * _data[8];
  result._data[5] = _data[6] * _data[6] + _data[7] * _data[7] + _data[8] * _data[8];

  return result;
}

/*
@LABEL:Tensor2::operator*(Tensor2 B)
@SHORT:Single contracted product of two second order tensors.
@RETURN:Tensor2 & Result of the multiplication operation.
@ARG:Tensor2 & B & Second tensor for the multiplication operation.
This method defines a single contracted product of two second order tensors.
The result of this operation is also a second order tensor defined by:
\begin{equation*}
\T = \A \cdot \B
\end{equation*}
where $\A$ is a second order tensor defined by the object itself and $\B$ is the second order tensor defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 Tensor2::operator*(const Tensor2 &T) const
//-----------------------------------------------------------------------------
{
  Tensor2 result;

  result._data[0] = _data[0] * T._data[0] + _data[1] * T._data[3] + _data[2] * T._data[6];
  result._data[1] = _data[0] * T._data[1] + _data[1] * T._data[4] + _data[2] * T._data[7];
  result._data[2] = _data[0] * T._data[2] + _data[1] * T._data[5] + _data[2] * T._data[8];
  result._data[3] = _data[3] * T._data[0] + _data[4] * T._data[3] + _data[5] * T._data[6];
  result._data[4] = _data[3] * T._data[1] + _data[4] * T._data[4] + _data[5] * T._data[7];
  result._data[5] = _data[3] * T._data[2] + _data[4] * T._data[5] + _data[5] * T._data[8];
  result._data[6] = _data[6] * T._data[0] + _data[7] * T._data[3] + _data[8] * T._data[6];
  result._data[7] = _data[6] * T._data[1] + _data[7] * T._data[4] + _data[8] * T._data[7];
  result._data[8] = _data[6] * T._data[2] + _data[7] * T._data[5] + _data[8] * T._data[8];

  return result;
}

/*
@LABEL:Tensor2::operator*(SymTensor2 B)
@SHORT:Single contracted product of a second order tensor and a symmetric second order tensor.
@RETURN:Tensor2 & Result of the multiplication operation.
@ARG:Tensor2 & B & Second tensor for the multiplication operation.
This method defines a single contracted product of a second order tensor and a symmetric second order tensor.
The result of this operation is also a second order tensor defined by:
\begin{equation*}
\T = \A \cdot \B
\end{equation*}
where $\A$ is a second order tensor defined by the object itself and $\B$ is a symmetric second order tensor defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 Tensor2::operator*(const SymTensor2 &T) const
//-----------------------------------------------------------------------------
{
  Tensor2 result;

  result._data[0] = _data[0] * T._data[0] + _data[1] * T._data[1] + _data[2] * T._data[2];
  result._data[1] = _data[0] * T._data[1] + _data[1] * T._data[3] + _data[2] * T._data[4];
  result._data[2] = _data[0] * T._data[2] + _data[1] * T._data[4] + _data[2] * T._data[5];
  result._data[3] = _data[3] * T._data[0] + _data[4] * T._data[1] + _data[5] * T._data[2];
  result._data[4] = _data[3] * T._data[1] + _data[4] * T._data[3] + _data[5] * T._data[4];
  result._data[5] = _data[3] * T._data[2] + _data[4] * T._data[4] + _data[5] * T._data[5];
  result._data[6] = _data[6] * T._data[0] + _data[7] * T._data[1] + _data[8] * T._data[2];
  result._data[7] = _data[6] * T._data[1] + _data[7] * T._data[3] + _data[8] * T._data[4];
  result._data[8] = _data[6] * T._data[2] + _data[7] * T._data[4] + _data[8] * T._data[5];

  return result;
}

/*
@LABEL:Tensor2::operator*(Vec3D V)
@SHORT:Multiplication of a second order tensor by a vector.
@RETURN:Vec3D & Result of the multiplication operation.
@ARG:Vec3D & V & Vec3D to use for the multiplication operation.
This method defines the product of a second order tensor by a vector.
The result of this operation is also a vector defined by:
\begin{equation*}
\overrightarrow{y} = \A \cdot \overrightarrow{x}
\end{equation*}
where $\A$ is a second order tensor defined by the object itself and $\overrightarrow{x}$ is a Vec3D defined by parameter V.
@END
*/
//-----------------------------------------------------------------------------
Vec3D Tensor2::operator*(const Vec3D &V) const
//-----------------------------------------------------------------------------
{
  Vec3D result;

  result._data[0] = _data[0] * V._data[0] + _data[1] * V._data[1] + _data[2] * V._data[2];
  result._data[1] = _data[3] * V._data[0] + _data[4] * V._data[1] + _data[5] * V._data[2];
  result._data[2] = _data[6] * V._data[0] + _data[7] * V._data[1] + _data[8] * V._data[2];

  return result;
}

/*
@LABEL:Tensor2::doubleDot()
@SHORT:Double contracted product of a second order tensor by itself.
@RETURN:double & Result of the multiplication operation.
This method defines a double contracted product of a second order tensor by itself.
The result of this operation is a scalar $s$ defined by:
\begin{equation*}
s = \A : \A = \sum_{i=1}^{3} \sum_{j=1}^{3} A_{ij}\times A_{ij}
\end{equation*}
where $\A$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
double Tensor2::doubleDot() const
//-----------------------------------------------------------------------------
{
  return (_data[0] * _data[0] + _data[1] * _data[1] + _data[2] * _data[2] +
          _data[3] * _data[3] + _data[4] * _data[4] + _data[5] * _data[5] +
          _data[6] * _data[6] + _data[7] * _data[7] + _data[8] * _data[8]);
}

/*
@LABEL:Tensor2::doubleDot(Tensor2 B)
@SHORT:Double contracted product of 2 second order tensors.
@RETURN:double & Result of the multiplication operation.
@ARG:Tensor2 & B & Second tensor for the multiplication operation.
This method defines a double contracted product of two second order tensors.
The result of this operation is a scalar $s$ defined by:
\begin{equation*}
s = \A : \B = \sum_{i=1}^{3} \sum_{j=1}^{3} A_{ij}\times B_{ij}
\end{equation*}
where $\A$ is a second order tensor defined by the object itself and $\B$ is a second order tensor defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
double Tensor2::doubleDot(const Tensor2 T) const
//-----------------------------------------------------------------------------
{
  return (_data[0] * T._data[0] + _data[1] * T._data[1] + _data[2] * T._data[2] +
          _data[3] * T._data[3] + _data[4] * T._data[4] + _data[5] * T._data[5] +
          _data[6] * T._data[6] + _data[7] * T._data[7] + _data[8] * T._data[8]);
}

/*
@LABEL:Tensor2::deviator()
@SHORT:Deviatoric part of a second order tensor.
@RETURN:Tensor2 & The deviatoric part of the second order tensor.
This method defines the deviatoric part of a second order tensor.
The result of this operation is a second order tensor defined by the following equation:
\begin{equation*}
\A^d=\A-\frac{1}{3}\tr[\A].\Id
\end{equation*}
where $\A$ is a second order tensor defined by the object itself and $\Id$ is the unit tensor.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 Tensor2::deviator() const
//-----------------------------------------------------------------------------
{
  Tensor2 result(*this);
  double pressure = thirdTrace();

  result._data[0] -= pressure;
  result._data[4] -= pressure;
  result._data[8] -= pressure;

  return result;
}

/*
@LABEL:Tensor2::transpose()
@SHORT:Transpose of a second order tensor.
@RETURN:Tensor2 & The transpose of the second order tensor.
This method defines the transpose of a second order tensor.
The result of this operation is a second order tensor defined by the following equation:
\begin{equation*}
\B=\A^T =\left[\begin{array}{ccc}
  A_{11} & A_{21} & A_{31}\\
  A_{12} & A_{22} & A_{32}\\
  A_{13} & A_{23} & A_{33}
  \end{array}\right]
\end{equation*}
where $\A$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 Tensor2::transpose() const
//-----------------------------------------------------------------------------
{
  return Tensor2(_data[0], _data[3], _data[6], _data[1], _data[4], _data[7], _data[2], _data[5], _data[8]);
}

/*
@LABEL:Tensor2::rowSum()
@SHORT:Sum of the rows of a second order tensor.
@RETURN:Vec3D & The sums of the rows of the second order tensor.
This method returns a vector by computing the sum of the components on all rows of a second order tensor.
The result of this operation is a vector $\overrightarrow{v}$ defined by:
\begin{equation*}
v_{i}=\sum_{j=1}^{3} T_{ji}
\end{equation*}
where $\T$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Vec3D Tensor2::rowSum() const
//-----------------------------------------------------------------------------
{
  return Vec3D(_data[0] + _data[1] + _data[2],
               _data[3] + _data[4] + _data[5],
               _data[6] + _data[7] + _data[8]);
}

/*
@LABEL:Tensor2::colSum()
@SHORT:Sum of the columns of a second order tensor.
@RETURN:Vec3D & The sums of the columns of the second order tensor.
This method returns a vector by computing the sum of the components on all columns of a second order tensor.
The result of this operation is a vector $\overrightarrow{v}$ defined by:
\begin{equation*}
v_{i}=\sum_{j=1}^{3}T_{ij}
\end{equation*}
where $\T$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Vec3D Tensor2::colSum() const
//-----------------------------------------------------------------------------
{
  return Vec3D(_data[0] + _data[3] + _data[6],
               _data[1] + _data[4] + _data[7],
               _data[2] + _data[5] + _data[8]);
}

/*
@LABEL:Tensor2::symmetric()
@SHORT:Symmetric part of a second order tensor.
@RETURN:Tensor2 & Symmetric part of the second order tensor.
This method returns the symmetric part of a second order tensor.
The result of this operation is a second order tensor defined by:
\begin{equation*}
\B = \left[\begin{array}{ccc}
 A_{11} & \frac{A_{12} + A_{21}}{2} & \frac{A_{13} + A_{31}}{2}\\
 \frac{A_{12} + A_{21}}{2} & A_{22} & \frac {A_{23} + A_{32}}{2}\\
 \frac{A_{13} + A_{31}}{2} & \frac {A_{23} + A_{32}}{2} & A_{33}\end{array}
\right]
\end{equation*}
where $\A$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 Tensor2::symmetric() const
//-----------------------------------------------------------------------------
{
  return Tensor2(_data[0],
                 (_data[1] + _data[3]) / 2,
                 (_data[2] + _data[6]) / 2,
                 (_data[1] + _data[3]) / 2,
                 _data[4],
                 (_data[5] + _data[7]) / 2,
                 (_data[2] + _data[6]) / 2,
                 (_data[5] + _data[7]) / 2,
                 _data[8]);
}

/*
@LABEL:Tensor2::skewSymmetric()
@SHORT:Skew-symmetric part of a second order tensor.
@RETURN:Tensor2 & Symmetric part of the second order tensor.
This method returns the skew-symmetric part of a second order tensor.
The result of this operation is a second order tensor defined by:
\begin{equation*}
\B = \left[\begin{array}{ccc}
 A_{11} & \frac{A_{12} - A_{21}}{2} & \frac{A_{13} - A_{31}}{2}\\
 -\frac{A_{12} -  A_{21}}{2} & A_{22} & \frac {A_{23} - A_{32}}{2}\\
 -\frac{A_{13} - A_{31}}{2} & -\frac {A_{23} - A_{32}}{2} & A_{33}\end{array}
\right]
\end{equation*}
where $\A$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 Tensor2::skewSymmetric() const
//-----------------------------------------------------------------------------
{
  return Tensor2(0,
                 (_data[1] - _data[3]) / 2,
                 (_data[2] - _data[6]) / 2,
                 (_data[3] - _data[1]) / 2,
                 0,
                 (_data[5] - _data[7]) / 2,
                 (_data[6] - _data[2]) / 2,
                 (_data[7] - _data[5]) / 2,
                 0);
}

/*
@LABEL:Tensor2::row(short r)
@SHORT:Extraction of a row from a second order tensor.
@RETURN:Vec3D & The extracted row.
@ARG:short & r & Row to extract
This method returns a vector as part of a second order tensor.
The result of this operation with the argument r is a vector defined by:
\begin{equation*}
v_{i} = T_{ri}
\end{equation*}
where $\T$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Vec3D Tensor2::row(short row) const
//-----------------------------------------------------------------------------
{
  Vec3D result;

  for (short i = 0; i < 3; i++)
  {
    result(i) = _data[_internalIndexes.index[row][i]];
  }

  return result;
}

/*
@LABEL:Tensor2::col(short c)
@SHORT:Extraction of a column from a second order tensor.
@RETURN:Vec3D & The extracted col.
@ARG:short & c & Column to extract
This method returns a vector as part of a second order tensor.
The result of this operation with the argument c is a vector defined by:
\begin{equation*}
v_{i} = T_{ic}
\end{equation*}
where $\T$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Vec3D Tensor2::col(short col) const
//-----------------------------------------------------------------------------
{
  Vec3D result;

  for (short i = 0; i < 3; i++)
  {
    result(i) = _data[_internalIndexes.index[i][col]];
  }

  return result;
}

//  Test the equality of two second order tensors
//-----------------------------------------------------------------------------
bool Tensor2::operator==(const Tensor2 &tensor) const
//-----------------------------------------------------------------------------
{
  short i;

  for (i = 0; i < 9; i++)
    if (_data[i] != tensor._data[i])
    {
      return false;
    }
  return true;
}

//  Test the inequality of two second order tensors
//-----------------------------------------------------------------------------
bool Tensor2::operator!=(const Tensor2 &tensor) const
//-----------------------------------------------------------------------------
{
  return !(*this == tensor);
}

//  Writes a second order tensor in a binary flux for storage
//-----------------------------------------------------------------------------
void Tensor2::write(std::ofstream &ofs) const
//-----------------------------------------------------------------------------
{
  ofs.write((char *)_data, 9 * sizeof(double));
}

//  Reads a second order tensor in a binary flux from storage
//-----------------------------------------------------------------------------
void Tensor2::read(std::ifstream &ifs)
//-----------------------------------------------------------------------------
{
  ifs.read((char *)_data, 9 * sizeof(double));
}

//  Writes a second order tensor in a binary flux for storage
//-----------------------------------------------------------------------------
std::ofstream &operator<<(std::ofstream &os, const Tensor2 &tensor)
//-----------------------------------------------------------------------------
{
  tensor.write(os);
  return os;
}

//  Reads a second order tensor from a binary flux for storage
//-----------------------------------------------------------------------------
std::ifstream &operator>>(std::ifstream &is, Tensor2 &tensor)
//-----------------------------------------------------------------------------
{
  tensor.read(is);
  return is;
}

/*
@LABEL:Tensor2::maxVal()
@SHORT:Maximum component in a second order tensor.
@RETURN:double & The maximum component of the second order tensor.
This method returns the maximum component in a second order tensor.
@END
*/
//-----------------------------------------------------------------------------
double Tensor2::maxVal()
//-----------------------------------------------------------------------------
{
  double max = _data[0];
  for (short i = 1; i < 9; i++)
  {
    if (_data[i] > max)
      max = _data[i];
  }
  return max;
}

/*
@LABEL:Tensor2::minVal()
@SHORT:Minimum component in a second order tensor.
@RETURN:double & The minimum component of the second order tensor.
This method returns the minimum component in a second order tensor.
@END
*/
//-----------------------------------------------------------------------------
double Tensor2::minVal()
//-----------------------------------------------------------------------------
{
  double min = _data[0];
  for (short i = 1; i < 9; i++)
  {
    if (_data[i] < min)
      min = _data[i];
  }
  return min;
}

/*
@LABEL:Tensor2::maxAbs()
@SHORT:Maximum absolute component in a second order tensor.
@RETURN:double & The maximum component of the second order tensor.
This method returns the maximum absolute component in a second order tensor.
@END
*/
//-----------------------------------------------------------------------------
double Tensor2::maxAbs()
//-----------------------------------------------------------------------------
{
  double max = dnlAbs(_data[0]);
  for (short i = 1; i < 9; i++)
  {
    if (dnlAbs(_data[i]) > max)
      max = dnlAbs(_data[i]);
  }
  return max;
}

/*
@LABEL:Tensor2::minAbs()
@SHORT:Minimum absolute component in a second order tensor.
@RETURN:double & The minimum component of the second order tensor.
This method returns the minimum absolute component in a second order tensor.
@END
*/
//-----------------------------------------------------------------------------
double Tensor2::minAbs()
//-----------------------------------------------------------------------------
{
  double min = dnlAbs(_data[0]);
  for (short i = 1; i < 9; i++)
  {
    if (dnlAbs(_data[i]) < min)
      min = dnlAbs(_data[i]);
  }
  return min;
}

/*
@LABEL:Tensor2::inverse()
@SHORT:Inverse of a second order tensor.
@RETURN:Tensor2 & The inverse of the second order tensor.
This method returns the inverse of a second order tensor.
The result of this operation is a second order tensor defined by:
\begin{equation*}
d = T_{11} T_{22} T_{33} + T_{21} T_{32} T_{13} + T_{31} T_{12} T_{23} - T_{31} T_{22} T_{13} - T_{11} T_{32} T_{23} - T_{21} T_{12} T_{33}
\end{equation*}
\begin{equation*}
T^{-1} = \frac {1}{d} \left[\begin{array}{ccc}
  T_{22}T_{33}-T_{23}T_{32}&T_{13}T_{32}-T_{12}T_{33}&T_{12}T_{23}-T_{13}T_{22}\\
  T_{23}T_{31}-T_{21}T_{33}&T_{11}T_{33}-T_{13}T_{31}&T_{13}T_{21}-T_{11}T_{23}\\
  T_{21}T_{32}-T_{22}T_{31}&T_{12}T_{31}-T_{11}T_{32}&T_{11}T_{22}-T_{12}T_{21}
  \end{array}
  \right]
\end{equation*}
where $\T$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 Tensor2::inverse() const
//-----------------------------------------------------------------------------
{
  double t1 = _data[4] * _data[8];
  double t2 = _data[2] * _data[7];
  double t3 = _data[1] * _data[5];
  double t4 = _data[5] * _data[7];
  double t5 = _data[1] * _data[8];
  double t6 = _data[2] * _data[4];

  double unSurDeter = 1 / (_data[0] * t1 +
                           _data[3] * t2 +
                           _data[6] * t3 -
                           _data[0] * t4 -
                           _data[3] * t5 -
                           _data[6] * t6);
  return Tensor2((t1 - t4) * unSurDeter,
                 (t2 - t5) * unSurDeter,
                 (t3 - t6) * unSurDeter,
                 (_data[5] * _data[6] - _data[3] * _data[8]) * unSurDeter,
                 (_data[0] * _data[8] - _data[2] * _data[6]) * unSurDeter,
                 (_data[2] * _data[3] - _data[0] * _data[5]) * unSurDeter,
                 (_data[3] * _data[7] - _data[4] * _data[6]) * unSurDeter,
                 (_data[1] * _data[6] - _data[0] * _data[7]) * unSurDeter,
                 (_data[0] * _data[4] - _data[1] * _data[3]) * unSurDeter);
}

/*
@LABEL:Tensor2::minor()
@SHORT:Minor of a second order tensor.
@RETURN:Tensor2 & The minor of the second order tensor.
This method returns the minor of a second order tensor.
\begin{equation*}
T^{minor} = \left[\begin{array}{ccc}
T_{22}T_{33}-T_{32}T_{23} & T_{33}T_{21}-T_{23}T_{31} & T_{21}T_{32}-T_{31}T_{22}\\
T_{12}T_{33}-T_{13}T_{32} & T_{33}T_{11}-T_{13}T_{31} & T_{11}T_{32}-T_{31}T_{12}\\
T_{12}T_{23}-T_{22}T_{13} & T_{23}T_{11}-T_{13}T_{21} & T_{11}T_{22}-T_{21}T_{12}
\end{array}
\right]
\end{equation*}
where $\T$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 Tensor2::minor() const
//-----------------------------------------------------------------------------
{
  return Tensor2(_data[4] * _data[8] - _data[7] * _data[5],
                 _data[8] * _data[3] - _data[5] * _data[6],
                 _data[3] * _data[7] - _data[6] * _data[4],
                 _data[1] * _data[8] - _data[2] * _data[7],
                 _data[8] * _data[0] - _data[2] * _data[6],
                 _data[0] * _data[7] - _data[6] * _data[1],
                 _data[1] * _data[5] - _data[4] * _data[2],
                 _data[5] * _data[0] - _data[2] * _data[3],
                 _data[0] * _data[4] - _data[3] * _data[1]);
}

/*
@LABEL:Tensor2::cofactors()
@SHORT:Cofactors of a second order tensor.
@RETURN:Tensor2 & The cofactor of the second order tensor.
This method returns the cofactor of a second order tensor.
\begin{equation*}
T^{cof} = \left[\begin{array}{ccc}
T_{22}T_{33}-T_{32}T_{23} & T_{23}T_{31}-T_{33}T_{21} & T_{21}T_{32}-T_{31}T_{22}\\
T_{13}T_{32}-T_{12}T_{33} & T_{33}T_{11}-T_{13}T_{31} & T_{31}T_{12}-T_{11}T_{32}\\
T_{12}T_{23}-T_{22}T_{13} & T_{13}T_{21}-T_{23}T_{11} & T_{11}T_{22}-T_{21}T_{12}
\end{array}
\right]
\end{equation*}
where $\T$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 Tensor2::cofactors() const
//-----------------------------------------------------------------------------
{
  return Tensor2(_data[4] * _data[8] - _data[7] * _data[5],
                 _data[5] * _data[6] - _data[8] * _data[3],
                 _data[3] * _data[7] - _data[6] * _data[4],
                 _data[2] * _data[7] - _data[1] * _data[8],
                 _data[8] * _data[0] - _data[2] * _data[6],
                 _data[6] * _data[1] - _data[0] * _data[7],
                 _data[1] * _data[5] - _data[4] * _data[2],
                 _data[2] * _data[3] - _data[5] * _data[0],
                 _data[0] * _data[4] - _data[3] * _data[1]);
}

/*
@LABEL:Tensor2::solve(Vec3D x)
@SHORT:Solves a small linear system $\A\cdot \overrightarrow{x} = \overrightarrow{b}$.
@RETURN:Vec3D & The solution of the linear system.
This method returns the solution of a small linear system with the following form:
\begin{equation*}
\overrightarrow{y} = \A \cdot \overrightarrow{x}
\end{equation*}
where $\A$ is a second order tensor defined by the object itself and $\overrightarrow{x}$ is a vector defined by parameter x.
@END
*/
//-----------------------------------------------------------------------------
Vec3D Tensor2::solve(const Vec3D &b) const
//-----------------------------------------------------------------------------
{
  return Vec3D(inverse() * b);
}
