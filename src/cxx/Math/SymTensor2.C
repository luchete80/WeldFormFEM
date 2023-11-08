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
#include <SymTensor2.h>
#include <NumpyInterface.h>
#include <Tensor2.h>
#include <Vec3D.h>
#include <Eigen3x3.h>
#include <lapacke.h>

const Tensor2Index SymTensor2::_internalIndexes = {0, 1, 2, 1, 3, 4, 2, 4, 5};

/*
@LABEL:SymTensor2::SymTensor2()
@SHORT:Default constructor of the SymTensor2 class.
@RETURN:SymTensor2 & The new SymTensor2 object created by the constructor.
This is the default constructor of the SymTensor2 class, where all components are initialized to zero by default.
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
SymTensor2::SymTensor2()
//-----------------------------------------------------------------------------
{
  setValue(0);
}

// Copy constructor
//-----------------------------------------------------------------------------
SymTensor2::SymTensor2(const SymTensor2 &T)
//-----------------------------------------------------------------------------
{
  memcpy(_data, T._data, 6 * sizeof(double));
}

/*
@LABEL:SymTensor2::SymTensor2(double t11, double t12, ...)
@SHORT:Constructor of the SymTensor2 class.
@RETURN:SymTensor2 & The new SymTensor2 object created by the constructor.
@ARG:double & t11 & Component $T_{11}$ of the symmetric second order tensor.
@ARG:double & t12 & Component $T_{12}$ of the symmetric second order tensor.
@ARG:double & t13 & Component $T_{13}$ of the symmetric second order tensor.
@ARG:double & t22 & Component $T_{22}$ of the symmetric second order tensor.
@ARG:double & t23 & Component $T_{23}$ of the symmetric second order tensor.
@ARG:double & t33 & Component $T_{33}$ of the symmetric second order tensor.
Constructor of a second order symmetric tensor with explicit initialization of the $6$ components of the tensor.
\begin{equation*}
\T=\left[\begin{array}{ccc}
  T_{11} & T_{12} & T_{13}\\
  T_{12} & T_{22} & T_{23}\\
  T_{13} & T_{23} & T_{33}
  \end{array}\right]
\end{equation*}
@END
*/
//-----------------------------------------------------------------------------
SymTensor2::SymTensor2(double t11, double t12, double t13, double t22, double t23, double t33)
//-----------------------------------------------------------------------------
{
  _data[0] = t11;
  _data[1] = t12;
  _data[2] = t13;
  _data[3] = t22;
  _data[4] = t23;
  _data[5] = t33;
}

// Destructor
//-----------------------------------------------------------------------------
SymTensor2::~SymTensor2()
//-----------------------------------------------------------------------------
{
}

// Send the content of a second order tensor to the output flux for display
//-----------------------------------------------------------------------------
std::ostream &operator<<(std::ostream &outputStream, const SymTensor2 &tensor)
//-----------------------------------------------------------------------------
{
  tensor.print(outputStream);
  return outputStream;
}

// Print the content of a second order tensor to the output flux for display
//-----------------------------------------------------------------------------
void SymTensor2::print(std::ostream &outputStream) const
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
@LABEL:SymTensor2::setToUnity()
@SHORT:Unity symmetric second order tensor.
@SMOD
This method transforms the current symmetric second order tensor to a unity tensor.
\begin{equation*}
\T=\left[\begin{array}{ccc}
1&0&0\\
0&1&0\\
0&0&1
\end{array}\right]
\end{equation*}
where $\T$ is a symmetric second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
void SymTensor2::setToUnity()
//-----------------------------------------------------------------------------
{
  _data[0] = 1;
  _data[1] = 0;
  _data[2] = 0;
  _data[3] = 1;
  _data[4] = 0;
  _data[5] = 1;
}

/*
@LABEL:SymTensor2::setToZero()
@SHORT:Zero symmetric second order tensor.
@SMOD
This method transforms the current symmetric second order tensor to a zero tensor.
\begin{equation*}
\T=\left[\begin{array}{ccc}
0 & 0 & 0\\
0 & 0 & 0\\
0 & 0 & 0
\end{array}\right]
\end{equation*}
where $\T$ is a symmetric second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
void SymTensor2::setToZero()
//-----------------------------------------------------------------------------
{
  setValue(0);
}

/*
@LABEL:SymTensor2::operator=(double v)
@SHORT:Fill a symmetric second order tensor with a scalar value.
@RETURN:SymTensor2
@ARG:double & v & Value to use for the operation.
This method is a surdefinition of the = operator for the symmetric second order tensor class.
\begin{equation*}
\T=\left[\begin{array}{ccc}
v & v & v\\
v & v & v\\
v & v & v
\end{array}\right]
\end{equation*}
where $\T$ is a symmetric second order tensor defined by the object itself and $v$ is the scalar value defined by parameter v.
@END
*/
//-----------------------------------------------------------------------------
SymTensor2 &SymTensor2::operator=(const double &val)
//-----------------------------------------------------------------------------
{
  setValue(val);
  return *this;
}

// Copy the content of a symmetric second order tensor into a new one
//-----------------------------------------------------------------------------
SymTensor2 &SymTensor2::operator=(const SymTensor2 &T)
//-----------------------------------------------------------------------------
{
  memcpy(_data, T._data, 6 * sizeof(double));
  return *this;
}

/*
@LABEL:SymTensor2::operator+(SymTensor2 B)
@SHORT:Addition of 2 symmetric second order tensors.
@ARG:SymTensor2 & B & Symmetric second order tensor to add to the current one.
@RETURN:SymTensor2 & Result of the addition operation.
This method defines the addition of 2 symmetric second order tensors.
The result of this operation is also a symmetric second order tensor defined by:
\begin{equation*}
\T = \A + \B
\end{equation*}
where $\A$ is a symmetric second order tensor defined by the object itself and $\B$ is the symmetric second order tensor value defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
SymTensor2 SymTensor2::operator+(const SymTensor2 &T) const
//-----------------------------------------------------------------------------
{
  // creation d'un nouveau tenseur
  SymTensor2 result;

  // calcul de la somme
  result._data[0] = _data[0] + T._data[0];
  result._data[1] = _data[1] + T._data[1];
  result._data[2] = _data[2] + T._data[2];
  result._data[3] = _data[3] + T._data[3];
  result._data[4] = _data[4] + T._data[4];
  result._data[5] = _data[5] + T._data[5];

  // renvoi de l'objet
  return result;
}

/*
@LABEL:SymTensor2::operator-(SymTensor2 B)
@SHORT:Subtraction of 2 symmetric second order tensors.
@ARG:SymTensor2 & B & Symmetric second order tensor to subtract to the current one.
@RETURN:SymTensor2 & Result of the difference operation.
This method defines the subtraction of 2 symmetric second order tensors.
The result of this operation is also a symmetric second order tensor defined by:
\begin{equation*}
\T = \A - \B
\end{equation*}
where $\A$ is a symmetric second order tensor defined by the object itself and $\B$ is the symmetric second order tensor value defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
SymTensor2 SymTensor2::operator-(const SymTensor2 &T) const
//-----------------------------------------------------------------------------
{
  // creation d'un nouveau tenseur
  SymTensor2 result;

  // calcul de la somme
  result._data[0] = _data[0] - T._data[0];
  result._data[1] = _data[1] - T._data[1];
  result._data[2] = _data[2] - T._data[2];
  result._data[3] = _data[3] - T._data[3];
  result._data[4] = _data[4] - T._data[4];
  result._data[5] = _data[5] - T._data[5];

  // renvoi de l'objet
  return result;
}

/*
@LABEL:SymTensor2::operator-()
@SHORT:Opposite of a symmetric second order tensor.
@RETURN:SymTensor2 & The opposite symmetric second order tensor.
This method defines the opposite of a symmetric second order tensor.
The result of this operation is also a symmetric second order tensor defined by:
\begin{equation*}
\T = - \A
\end{equation*}
where $\A$ is a symmetric second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
SymTensor2 SymTensor2::operator-() const
//-----------------------------------------------------------------------------
{
  SymTensor2 result;

  result._data[0] = -_data[0];
  result._data[1] = -_data[1];
  result._data[2] = -_data[2];
  result._data[3] = -_data[3];
  result._data[4] = -_data[4];
  result._data[5] = -_data[5];

  return result;
}

/*
@LABEL:SymTensor2::operator*(double l)
@SHORT:Multiplication of a symmetric second order tensor by a scalar.
@ARG:double & l & Scalar value to use for the operation.
@RETURN:SymTensor2 & Result of the multiplication operation.
This method defines the multiplication of a symmetric second order tensor by a scalar value.
The result of this operation is also a symmetric second order tensor defined by:
\begin{equation*}
\T = \lambda \A
\end{equation*}
where $\A$ is a symmetric second order tensor defined by the object itself and $\lambda$ is the scalar value defined by parameter l.
@END
*/
//-----------------------------------------------------------------------------
SymTensor2 SymTensor2::operator*(const double l) const
//-----------------------------------------------------------------------------
{
  SymTensor2 result;

  result._data[0] = l * _data[0];
  result._data[1] = l * _data[1];
  result._data[2] = l * _data[2];
  result._data[3] = l * _data[3];
  result._data[4] = l * _data[4];
  result._data[5] = l * _data[5];

  return result;
}

/*
@LABEL:SymTensor2::operator/(double l)
@SHORT:Division of a symmetric second order tensor by a scalar.
@ARG:double & l & Scalar value to use for the operation.
@RETURN:SymTensor2 & Result of the division operation.
This method defines the division of a symmetric second order tensor by a scalar value.
The result of this operation is also a symmetric second order tensor defined by:
\begin{equation*}
\T = \frac{1}{\lambda} \A
\end{equation*}
where $\A$ is a symmetric second order tensor defined by the object itself and $\lambda$ is the scalar value defined by parameter l.
@END
*/
//-----------------------------------------------------------------------------
SymTensor2 SymTensor2::operator/(const double l) const
//-----------------------------------------------------------------------------
{
  SymTensor2 result;

#ifdef VERIF_maths
  if (l == 0)
  {
    fatalError("SymTensor2:: operator /", "divide by zero");
  }
#endif

  result._data[0] = _data[0] / l;
  result._data[1] = _data[1] / l;
  result._data[2] = _data[2] / l;
  result._data[3] = _data[3] / l;
  result._data[4] = _data[4] / l;
  result._data[5] = _data[5] / l;

  return result;
}

/*
@LABEL:operator*(double l, SymTensor2 A)
@SHORT:Multiplication of a symmetric second order tensor by a scalar.
@ARG:double & l & Scalar value to use for the operation.
@ARG:SymTensor2 & A & Symmetric second order tensor to use for the operation.
@RETURN:SymTensor2 & Result of the multiplication operation.
This method defines the multiplication of a symmetric second order tensor by a scalar value.
The result of this operation is also a symmetric second order tensor defined by:
\begin{equation*}
\T = \lambda \A
\end{equation*}
where $\A$ is a symmetric second order tensor and $\lambda$ is the scalar value defined by parameter l.
@END
*/
//-----------------------------------------------------------------------------
SymTensor2 operator*(const double &l, const SymTensor2 &T)
//-----------------------------------------------------------------------------
{
  SymTensor2 result;

  result._data[0] = l * T._data[0];
  result._data[1] = l * T._data[1];
  result._data[2] = l * T._data[2];
  result._data[3] = l * T._data[3];
  result._data[4] = l * T._data[4];
  result._data[5] = l * T._data[5];

  return result;
}

/*
@LABEL:SymTensor2::dot()
@SHORT:Single contracted product of a symmetric second order tensor by itself.
@RETURN:SymTensor2 & Result of the multiplication operation.
This method defines a single contracted product of a symmetric second order tensor by itself.
The result of this operation is also a symmetric second order tensor defined by:
\begin{equation*}
\T = \A \cdot \A
\end{equation*}
where $\A$ is a symmetric second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
SymTensor2 SymTensor2::dot() const
//-----------------------------------------------------------------------------
{
  SymTensor2 result;

  result._data[0] = _data[0] * _data[0] + _data[1] * _data[1] + _data[2] * _data[2];
  result._data[1] = _data[0] * _data[1] + _data[1] * _data[3] + _data[2] * _data[4];
  result._data[2] = _data[0] * _data[2] + _data[1] * _data[4] + _data[2] * _data[5];
  result._data[3] = _data[1] * _data[1] + _data[3] * _data[3] + _data[4] * _data[4];
  result._data[4] = _data[1] * _data[2] + _data[3] * _data[4] + _data[4] * _data[5];
  result._data[5] = _data[2] * _data[2] + _data[4] * _data[4] + _data[5] * _data[5];

  return result;
}

/*
@LABEL:SymTensor2::operator*(SymTensor2 B)
@SHORT:Single contracted product of two symmetric second order tensors.
@RETURN:Tensor2 & Result of the multiplication operation.
@ARG:SymTensor2 & B & Symmetric second tensor for the multiplication operation.
@WARNING: The result of this operation is a non symmetric second order tensor.
This method defines a single contracted product of two symmetric second order tensors.
The result of this operation is also a symmetric second order tensor defined by:
\begin{equation*}
\T = \A \cdot \B
\end{equation*}
where $\A$ is a symmetric second order tensor defined by the object itself and $\B$ is the symmetric second order tensor defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 SymTensor2::operator*(const SymTensor2 &T) const
//-----------------------------------------------------------------------------
{
  Tensor2 result;

  result._data[0] = _data[0] * T._data[0] + _data[1] * T._data[1] + _data[2] * T._data[2];
  result._data[1] = _data[0] * T._data[1] + _data[1] * T._data[3] + _data[2] * T._data[4];
  result._data[2] = _data[0] * T._data[2] + _data[1] * T._data[4] + _data[2] * T._data[5];
  result._data[3] = _data[1] * T._data[0] + _data[3] * T._data[1] + _data[4] * T._data[2];
  result._data[4] = _data[1] * T._data[1] + _data[3] * T._data[3] + _data[4] * T._data[4];
  result._data[5] = _data[1] * T._data[2] + _data[3] * T._data[4] + _data[4] * T._data[5];
  result._data[6] = _data[2] * T._data[0] + _data[4] * T._data[1] + _data[5] * T._data[2];
  result._data[7] = _data[2] * T._data[1] + _data[4] * T._data[3] + _data[5] * T._data[4];
  result._data[8] = _data[2] * T._data[2] + _data[4] * T._data[4] + _data[5] * T._data[5];

  return result;
}

/*
@LABEL:SymTensor2::operator*(Tensor2 B)
@SHORT:Single contracted product of two symmetric second order tensors.
@RETURN:Tensor2 & Result of the multiplication operation.
@ARG:Tensor2 & B & Symmetric second tensor for the multiplication operation.
@WARNING: The result of this operation is a non symmetric second order tensor.
This method defines a single contracted product of two symmetric second order tensors.
The result of this operation is also a symmetric second order tensor defined by:
\begin{equation*}
\T = \A \cdot \B
\end{equation*}
where $\A$ is a symmetric second order tensor defined by the object itself and $\B$ is the symmetric second order tensor defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 SymTensor2::operator*(const Tensor2 &T) const
//-----------------------------------------------------------------------------
{
  Tensor2 result;

  result._data[0] = _data[0] * T._data[0] + _data[1] * T._data[3] + _data[2] * T._data[6];
  result._data[1] = _data[0] * T._data[1] + _data[1] * T._data[4] + _data[2] * T._data[7];
  result._data[2] = _data[0] * T._data[2] + _data[1] * T._data[5] + _data[2] * T._data[8];
  result._data[3] = _data[1] * T._data[0] + _data[3] * T._data[3] + _data[4] * T._data[6];
  result._data[4] = _data[1] * T._data[1] + _data[3] * T._data[4] + _data[4] * T._data[7];
  result._data[5] = _data[1] * T._data[2] + _data[3] * T._data[5] + _data[4] * T._data[8];
  result._data[6] = _data[2] * T._data[0] + _data[5] * T._data[3] + _data[5] * T._data[6];
  result._data[7] = _data[2] * T._data[1] + _data[4] * T._data[4] + _data[5] * T._data[7];
  result._data[8] = _data[2] * T._data[2] + _data[4] * T._data[5] + _data[5] * T._data[8];

  return result;
}

/*
@LABEL:SymTensor2::dotRxRT(Tensor2 R)
@SHORT:Special combination for a multiplication of 3 symmetric second order tensors.
@RETURN:SymTensor2 & Result of the multiplication operation.
@ARG:Tensor2 & R & Symmetric second Tensor for the multiplication operation.
This method defines the product of a symmetric tensor by two rotations defined by the following equation:
\begin{equation*}
\R = \Q \cdot \A \cdot \Q^T =\left[\begin{array}{ccc}
R_0 & R_1 & R_2\\
& R_3 & R_4\\
sym &  & R_5
\end{array}\right]
\end{equation*}
with:
\begin{align*}
R_0 & = A_0 Q_0^2 + 2 (A_1 Q_1 + A_2 Q_2)Q_0 + A_3 Q_1^2 + 2 A_4 Q_1 Q_2 + A_5 Q_2^2\\
R_1 & = (A_0 Q_3 + A_1 Q_4 + A_2 Q_5)Q_0 + (A_1 Q_3 + A_3 Q_4 + A_4 Q_5)Q_1 + (A_2 Q_3 + A_4 Q_4 + A_5 Q_5)Q_2\\
R_2 & = (A_0 Q_6 + A_1 Q_7 + A_2 Q_8)Q_0 + (A_1 Q_6 + A_3 Q_7 + A_4 Q_8)Q_1 + (A_2 Q_6 + A_4 Q_7 + A_5 Q_8)Q_2\\
R_3 & = A_0 Q_3^2 + 2 (A_1 Q_4 + A_2 Q_5)Q_3+ A_3 Q_4^2 + 2 A_4 Q_4 Q_5 + A_5 Q_5^2\\
R_4 & = (A_0 Q_6 + A_1 Q_7 + A_2 Q_8)Q_3 + (A_1 Q_6 + A_3 Q_7 + A_4 Q_8)Q_4 + (A_2 Q_6 + A_4 Q_7 + A_5 Q_8)Q_5\\
R_5 & = A_0 Q_6^2 + 2 (A_1 Q_7 + A_2 Q_8)Q_6+ A_3 Q_7^2 + 2 A_4 Q_7 Q_8 + A_5 Q_8^2
\end{align*}
where $\A$ is a symmetric second order tensor defined by the object itself and $\Q$ an orthogonal tensor.
@END
*/
//-----------------------------------------------------------------------------
SymTensor2 SymTensor2::dotRxRT(const Tensor2 R) const
//-----------------------------------------------------------------------------
{
  SymTensor2 result;

  result._data[0] = R._data[0] * R._data[0] * _data[0] + 2 * R._data[0] * (R._data[1] * _data[1] + R._data[2] * _data[2]) + R._data[1] * R._data[1] * _data[3] + 2 * R._data[1] * R._data[2] * _data[4] + R._data[2] * R._data[2] * _data[5];
  result._data[1] = R._data[0] * (R._data[3] * _data[0] + R._data[4] * _data[1] + R._data[5] * _data[2]) + R._data[1] * (R._data[3] * _data[1] + R._data[4] * _data[3] + R._data[5] * _data[4]) + R._data[2] * (R._data[3] * _data[2] + R._data[4] * _data[4] + R._data[5] * _data[5]);
  result._data[2] = R._data[0] * (R._data[6] * _data[0] + R._data[7] * _data[1] + R._data[8] * _data[2]) + R._data[1] * (R._data[6] * _data[1] + R._data[7] * _data[3] + R._data[8] * _data[4]) + R._data[2] * (R._data[6] * _data[2] + R._data[7] * _data[4] + R._data[8] * _data[5]);
  result._data[3] = R._data[3] * R._data[3] * _data[0] + 2 * R._data[3] * (R._data[4] * _data[1] + R._data[5] * _data[2]) + R._data[4] * R._data[4] * _data[3] + 2 * R._data[4] * R._data[5] * _data[4] + R._data[5] * R._data[5] * _data[5];
  result._data[4] = R._data[3] * (R._data[6] * _data[0] + R._data[7] * _data[1] + R._data[8] * _data[2]) + R._data[4] * (R._data[6] * _data[1] + R._data[7] * _data[3] + R._data[8] * _data[4]) + R._data[5] * (R._data[6] * _data[2] + R._data[7] * _data[4] + R._data[8] * _data[5]);
  result._data[5] = R._data[6] * R._data[6] * _data[0] + 2 * R._data[6] * (R._data[7] * _data[1] + R._data[8] * _data[2]) + R._data[7] * R._data[7] * _data[3] + 2 * R._data[7] * R._data[8] * _data[4] + R._data[8] * R._data[8] * _data[5];

  return result;
}

/*
@LABEL:SymTensor2::dotRTxR(Tensor2 R)
@SHORT:Special combination for a multiplication of 3 second order tensors.
@RETURN:SymTensor2 & Result of the multiplication operation.
@ARG:Tensor2&R&Second Tensor for the multiplication operation.
This method defines the product of a symmetric tensor by two rotations defined by the following equation:
\begin{equation*}
\R = \Q^T \cdot \A \cdot \Q =\left[\begin{array}{ccc}
R_0 & R_1 & R_2\\
& R_3 & R_4\\
sym &  & R_5
\end{array}\right]
\end{equation*}
with:
\begin{align*}
R_0 &= A_0 Q_0^2 + 2 (A_1 Q_3 + A_2 Q_6) Q_0 + A_3 Q_3^2 + 2 A_4 Q_3 Q_6 + A_5 Q_6^2 \\
R_1 &= A_1 Q_1 Q_3 +A_2 Q_1 Q_6 + (A_0 Q_1 + A_1 Q_4 + A_2 Q_7) Q_0+ A_3 Q_3 Q_4 + A_4 Q_4 Q_6 + A_4 Q_3 Q_7 + A_5 Q_6 Q_7\\
R_2 &= A_1 Q_2 Q_3 + A_2 Q_2 Q_6 + Q_0 (A_0 Q_2 + A_1 Q_5 + A_2 Q_8) + A_3 Q_3 Q_5 + A_4 Q_5 Q_6 + A_4 Q_3 Q_8 + A_5 Q_6 Q_8\\
R_3 &= A_0 Q_1^2 + 2 (A_1 Q_4 + A_2 Q_7)Q_1 + A_3 Q_4^2 + 2 A_4 Q_4 Q_7 + A_5 Q_7^2\\
R_4 &= A_1 Q_2 Q_4 + A_2 Q_2 Q_7 + (A_0 Q_2 + A_1 Q_5 + A_2 Q_8)Q_1 + A_3 Q_4 Q_5 + A_4 Q_5 Q_7 + A_4 Q_4 Q_8 + A_5 Q_7 Q_8\\
R_5 &= A_0 Q_2^2 + 2 (A_1 Q_5 + A_2 Q_8) Q_2 + A_3 Q_5^2 + 2 A_4 Q_5 Q_8 + A_5 Q_8^2
\end{align*}
where $\A$ is a symmetric second order tensor defined by the object itself and $\Q$ an orthogonal tensor.
@END
*/
//-----------------------------------------------------------------------------
SymTensor2 SymTensor2::dotRTxR(const Tensor2 R) const
//-----------------------------------------------------------------------------
{
  SymTensor2 result;

  result._data[0] = R._data[0] * R._data[0] * _data[0] + 2 * R._data[0] * (R._data[3] * _data[1] + R._data[6] * _data[2]) + R._data[3] * R._data[3] * _data[3] + 2 * R._data[3] * R._data[6] * _data[4] + R._data[6] * R._data[6] * _data[5];
  result._data[1] = R._data[1] * R._data[3] * _data[1] + R._data[1] * R._data[6] * _data[2] + R._data[0] * (R._data[1] * _data[0] + R._data[4] * _data[1] + R._data[7] * _data[2]) + R._data[3] * R._data[4] * _data[3] + R._data[4] * R._data[6] * _data[4] + R._data[3] * R._data[7] * _data[4] + R._data[6] * R._data[7] * _data[5];
  result._data[2] = R._data[2] * R._data[3] * _data[1] + R._data[2] * R._data[6] * _data[2] + R._data[0] * (R._data[2] * _data[0] + R._data[5] * _data[1] + R._data[8] * _data[2]) + R._data[3] * R._data[5] * _data[3] + R._data[5] * R._data[6] * _data[4] + R._data[3] * R._data[8] * _data[4] + R._data[6] * R._data[8] * _data[5];
  result._data[3] = R._data[1] * R._data[1] * _data[0] + 2 * R._data[1] * (R._data[4] * _data[1] + R._data[7] * _data[2]) + R._data[4] * R._data[4] * _data[3] + 2 * R._data[4] * R._data[7] * _data[4] + R._data[7] * R._data[7] * _data[5];
  result._data[4] = R._data[2] * R._data[4] * _data[1] + R._data[2] * R._data[7] * _data[2] + R._data[1] * (R._data[2] * _data[0] + R._data[5] * _data[1] + R._data[8] * _data[2]) + R._data[4] * R._data[5] * _data[3] + R._data[5] * R._data[7] * _data[4] + R._data[4] * R._data[8] * _data[4] + R._data[7] * R._data[8] * _data[5];
  result._data[5] = R._data[2] * R._data[2] * _data[0] + 2 * R._data[2] * (R._data[5] * _data[1] + R._data[8] * _data[2]) + R._data[5] * R._data[5] * _data[3] + 2 * R._data[5] * R._data[8] * _data[4] + R._data[8] * R._data[8] * _data[5];

  return result;
}

/*
@LABEL:SymTensor2::operator*(Vec3D V)
@SHORT:Multiplication of a symmetric second order tensor by a vector.
@RETURN:Vec3D & Result of the multiplication operation.
@ARG:Vec3D & V & Vec3D to use for the multiplication operation.
This method defines the product of a symmetric second order tensor by a vector.
The result of this operation is also a vector defined by:
\begin{equation*}
\overrightarrow{y} = \A \cdot \overrightarrow{x}
\end{equation*}
where $\A$ is a symmetric second order tensor defined by the object itself and $\overrightarrow{x}$ is a Vec3D defined by parameter V.
@END
*/
//-----------------------------------------------------------------------------
Vec3D SymTensor2::operator*(const Vec3D &V) const
//-----------------------------------------------------------------------------
{
  Vec3D result;

  result._data[0] = _data[0] * V._data[0] + _data[1] * V._data[1] + _data[2] * V._data[2];
  result._data[1] = _data[1] * V._data[0] + _data[3] * V._data[1] + _data[4] * V._data[2];
  result._data[2] = _data[2] * V._data[0] + _data[4] * V._data[1] + _data[5] * V._data[2];

  return result;
}

/*
@LABEL:SymTensor2::dot(SymTensor2 B)
@SHORT:Single contracted product of two symmetric second order tensors.
@RETURN:Tensor2 & Result of the multiplication operation.
@ARG:SymTensor2 & B & Second tensor for the multiplication operation.
This method defines a single contracted product of two symmetric second order tensors.
The result of this operation is also a symmetric second order tensor defined by:
\begin{equation*}
\T = \A \cdot \B
\end{equation*}
where $\A$ is a symmetric second order tensor defined by the object itself and $\B$ is the symmetric second order tensor defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 SymTensor2::dot(const SymTensor2 T) const
//-----------------------------------------------------------------------------
{
  Tensor2 result;

  result._data[0] = _data[0] * T._data[0] + _data[1] * T._data[1] + _data[2] * T._data[2];
  result._data[1] = _data[0] * T._data[1] + _data[1] * T._data[3] + _data[2] * T._data[4];
  result._data[2] = _data[0] * T._data[2] + _data[1] * T._data[4] + _data[2] * T._data[5];
  result._data[3] = _data[1] * T._data[0] + _data[3] * T._data[1] + _data[4] * T._data[2];
  result._data[4] = _data[1] * T._data[1] + _data[3] * T._data[3] + _data[4] * T._data[4];
  result._data[5] = _data[1] * T._data[2] + _data[3] * T._data[4] + _data[4] * T._data[5];
  result._data[6] = _data[2] * T._data[0] + _data[4] * T._data[1] + _data[5] * T._data[2];
  result._data[7] = _data[2] * T._data[1] + _data[4] * T._data[3] + _data[5] * T._data[4];
  result._data[8] = _data[2] * T._data[2] + _data[4] * T._data[4] + _data[5] * T._data[5];

  return result;
}

/*
@LABEL:SymTensor2::doubleDot(SymTensor2 B)
@SHORT:Double contracted product of 2 symmetric second order tensors.
@RETURN:double & Result of the multiplication operation.
@ARG:SymTensor2 & B & Second tensor for the multiplication operation.
This method defines a double contracted product of two symmetric second order tensors.
The result of this operation is a scalar $s$ defined by:
\begin{equation*}
s = \A : \B = \sum_{i=1}^{3} \sum_{j=1}^{3} A_{ij}\times B_{ij}
\end{equation*}
where $\A$ is a symmetric second order tensor defined by the object itself and $\B$ is a symmetric second order tensor defined by parameter B.
@END
*/
//-----------------------------------------------------------------------------
double SymTensor2::doubleDot(const SymTensor2 T) const
//-----------------------------------------------------------------------------
{
  return (_data[0] * T._data[0] +
          _data[3] * T._data[3] +
          _data[5] * T._data[5] +
          2.0 * (_data[1] * T._data[1] +
                 _data[2] * T._data[2] +
                 _data[4] * T._data[4]));
}

/*
@LABEL:SymTensor2::doubleDot()
@SHORT:Double contracted product of a symmetric second order tensor by itself.
@RETURN:double & Result of the multiplication operation.
This method defines a double contracted product of a symmetric second order tensor by itself.
The result of this operation is a scalar $s$ defined by:
\begin{equation*}
s = \A : \A = \sum_{i=1}^{3} \sum_{j=1}^{3} A_{ij}\times A_{ij}
\end{equation*}
where $\A$ is a symmetric second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
double SymTensor2::doubleDot() const
//-----------------------------------------------------------------------------
{
  return (_data[0] * _data[0] +
          _data[3] * _data[3] +
          _data[5] * _data[5] +
          2.0 * (_data[1] * _data[1] +
                 _data[2] * _data[2] +
                 _data[4] * _data[4]));
}

/*
@LABEL:SymTensor2::deviator()
@SHORT:Deviatoric part of a symmetric second order tensor.
@RETURN:SymTensor2 & The deviatoric part of the symmetric second order tensor.
This method defines the deviatoric part of a symmetric second order tensor.
The result of this operation is a symmetric second order tensor defined by the following equation:
\begin{equation*}
\A^d=\A-\frac{1}{3}\tr[\A].\Id
\end{equation*}
where $\A$ is a symmetric second order tensor defined by the object itself and $\Id$ is the unit tensor.
@END
*/
//-----------------------------------------------------------------------------
SymTensor2 SymTensor2::deviator() const
//-----------------------------------------------------------------------------
{
  SymTensor2 result(*this);
  double pressure = thirdTrace();

  result._data[0] -= pressure;
  result._data[3] -= pressure;
  result._data[5] -= pressure;

  return result;
}

/*
@LABEL:SymTensor2::rowSum()
@SHORT:Sum of the rows of a symmetric second order tensor.
@RETURN:Vec3D & The sums of the rows of the symmetric second order tensor.
This method returns a vector by computing the sum of the components on all rows of a symmetric second order tensor.
The result of this operation is a vector $\overrightarrow{v}$ defined by:
\begin{equation*}
v_{i}=\sum_{j=1}^{3} T_{ji}
\end{equation*}
where $\T$ is a symmetric second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Vec3D SymTensor2::rowSum() const
//-----------------------------------------------------------------------------
{
  return Vec3D(_data[0] + _data[1] + _data[2],
               _data[1] + _data[3] + _data[4],
               _data[2] + _data[4] + _data[5]);
}

/*
@LABEL:SymTensor2::colSum()
@SHORT:Sum of the columns of a symmetric second order tensor.
@RETURN:Vec3D & The sums of the columns of the symmetric second order tensor.
This method returns a vector by computing the sum of the components on all columns of a symmetric second order tensor.
The result of this operation is a vector $\overrightarrow{v}$ defined by:
\begin{equation*}
v_{i}=\sum_{j=1}^{3}T_{ij}
\end{equation*}
where $\T$ is a symmetric second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Vec3D SymTensor2::colSum() const
//-----------------------------------------------------------------------------
{
  return Vec3D(_data[0] + _data[1] + _data[2],
               _data[1] + _data[3] + _data[4],
               _data[2] + _data[4] + _data[5]);
}

/*
@LABEL:SymTensor2::row(short r)
@SHORT:Extraction of a row from a symmetric second order tensor.
@RETURN:Vec3D & The extracted row.
@ARG:short & r & Row to extract
This method returns a vector as part of a symmetric second order tensor.
The result of this operation with the argument r is a vector defined by:
\begin{equation*}
v_{i} = T_{ri}
\end{equation*}
where $\T$ is a symmetric second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Vec3D SymTensor2::row(short row) const
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
@LABEL:SymTensor2::col(short c)
@SHORT:Extraction of a column from a symmetric second order tensor.
@RETURN:Vec3D & The extracted col.
@ARG:short & c & Column to extract
This method returns a vector as part of a symmetric second order tensor.
The result of this operation with the argument c is a vector defined by:
\begin{equation*}
v_{i} = T_{ic}
\end{equation*}
where $\T$ is a symmetric second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
Vec3D SymTensor2::col(short col) const
//-----------------------------------------------------------------------------
{
  Vec3D result;

  for (short i = 0; i < 3; i++)
  {
    result(i) = _data[_internalIndexes.index[i][col]];
  }

  return result;
}

//  Test the equality of two symmetric second order tensors
//-----------------------------------------------------------------------------
bool SymTensor2::operator==(const SymTensor2 &T) const
//-----------------------------------------------------------------------------
{
  short i;

  for (i = 0; i < 6; i++)
    if (_data[i] != T._data[i])
    {
      return false;
    }
  return true;
}

//  Test the inequality of two symmetric second order tensors
//-----------------------------------------------------------------------------
bool SymTensor2::operator!=(const SymTensor2 &tensor) const
//-----------------------------------------------------------------------------
{
  return !(*this == tensor);
}

//  Writes a symmetric second order tensor in a binary flux for storage
//-----------------------------------------------------------------------------
void SymTensor2::write(std::ofstream &ofs) const
//-----------------------------------------------------------------------------
{
  ofs.write((char *)_data, 6 * sizeof(double));
}

//  Reads a symmetric second order tensor in a binary flux from storage
//-----------------------------------------------------------------------------
void SymTensor2::read(std::ifstream &ifs)
//-----------------------------------------------------------------------------
{
  ifs.read((char *)_data, 6 * sizeof(double));
}

//  Writes a symmetric second order tensor in a binary flux for storage
//-----------------------------------------------------------------------------
std::ofstream &operator<<(std::ofstream &os, const SymTensor2 &tensor)
//-----------------------------------------------------------------------------
{
  tensor.write(os);
  return os;
}

//  Reads a symmetric second order tensor from a binary flux for storage
//-----------------------------------------------------------------------------
std::ifstream &operator>>(std::ifstream &is, SymTensor2 &tensor)
//-----------------------------------------------------------------------------
{
  tensor.read(is);
  return is;
}

/*
@LABEL:SymTensor2::maxVal()
@SHORT:Maximum component in a symmetric second order tensor.
@RETURN:double & The maximum component of the symmetric second order tensor.
This method returns the maximum component in a symmetric second order tensor.
@END
*/
//-----------------------------------------------------------------------------
double SymTensor2::maxVal()
//-----------------------------------------------------------------------------
{
  double max = _data[0];
  for (short i = 1; i < 6; i++)
  {
    if (_data[i] > max)
      max = _data[i];
  }
  return max;
}

/*
@LABEL:SymTensor2::minVal()
@SHORT:Minimum component in a symmetric second order tensor.
@RETURN:double & The minimum component of the symmetric second order tensor.
This method returns the minimum component in a symmetric second order tensor.
@END
*/
//-----------------------------------------------------------------------------
double SymTensor2::minVal()
//-----------------------------------------------------------------------------
{
  double min = _data[0];
  for (short i = 1; i < 6; i++)
  {
    if (_data[i] < min)
      min = _data[i];
  }
  return min;
}

/*
@LABEL:SymTensor2::maxAbs()
@SHORT:Maximum absolute component in a symmetric second order tensor.
@RETURN:double & The maximum component of the symmetric second order tensor.
This method returns the maximum absolute component in a symmetric second order tensor.
@END
*/
//-----------------------------------------------------------------------------
double SymTensor2::maxAbs()
//-----------------------------------------------------------------------------
{
  double max = dnlAbs(_data[0]);
  for (short i = 1; i < 6; i++)
  {
    if (dnlAbs(_data[i]) > max)
      max = dnlAbs(_data[i]);
  }
  return max;
}

/*
@LABEL:SymTensor2::minAbs()
@SHORT:Minimum absolute component in a symmetric second order tensor.
@RETURN:double & The minimum component of the symmetric second order tensor.
This method returns the minimum absolute component in a symmetric second order tensor.
@END
*/
//-----------------------------------------------------------------------------
double SymTensor2::minAbs()
//-----------------------------------------------------------------------------
{
  double min = dnlAbs(_data[0]);
  for (short i = 1; i < 6; i++)
  {
    if (dnlAbs(_data[i]) < min)
      min = dnlAbs(_data[i]);
  }
  return min;
}

/*
@LABEL:SymTensor2::inverse()
@SHORT:Inverse of a symmetric second order tensor.
@RETURN:SymTensor2 & The inverse of the symmetric second order tensor.
This method returns the inverse of a symmetric second order tensor.
The result of this operation is a symmetric second order tensor defined by:
\begin{equation*}
D = T_{11} T_{22} T_{33} + 2 T_{12} T_{23} T_{13} - T_{22} T_{13}^2 - T_{11} T_{23}^2 - T_{33} T_{12}^2
\end{equation*}
\begin{equation*}
T^{-1} = \frac {1}{D} \left[\begin{array}{ccc}
  T_{22}T_{33}-T_{23}^2&T_{13}T_{23}-T_{12}T_{33}&T_{12}T_{23}-T_{13}T_{22}\\
  T_{13}T_{23}-T_{12}T_{33}&T_{11}T_{33}-T_{13}^2&T_{12}T_{13}-T_{11}T_{23}\\
  T_{12}T_{23}-T_{13}T_{22}&T_{12}T_{13}-T_{11}T_{23}&T_{11}T_{22}-T_{12}^2
  \end{array}
  \right]
\end{equation*}
where $\T$ is a symmetric second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
SymTensor2 SymTensor2::inverse() const
//-----------------------------------------------------------------------------
{

  double t1 = _data[3] * _data[5];
  double t2 = _data[2] * _data[4];
  double t4 = _data[4] * _data[4];
  double t5 = _data[1] * _data[5];
  double t6 = _data[2] * _data[3];

  double unSurDeter = 1 / (_data[0] * t1 + 2.0 * _data[1] * t2 - _data[0] * t4 - _data[1] * t5 - _data[2] * t6);
  return SymTensor2((t1 - t4) * unSurDeter,
                    (t2 - t5) * unSurDeter,
                    (_data[1] * _data[4] - t6) * unSurDeter,
                    (_data[0] * _data[5] - _data[2] * _data[2]) * unSurDeter,
                    (_data[2] * _data[1] - _data[0] * _data[4]) * unSurDeter,
                    (_data[0] * _data[3] - _data[1] * _data[1]) * unSurDeter);
}

/*
@LABEL:SymTensor2::minor()
@SHORT:Minor of a symmetric second order tensor.
@RETURN:SymTensor2 & The minor of the symmetric second order tensor.
This method returns the minor of a symmetric second order tensor.
\begin{equation*}
T^{minor} = \left[\begin{array}{ccc}
T_{22}T_{33}-T_{23}T_{23} & T_{33}T_{12}-T_{23}T_{13} & T_{12}T_{23}-T_{13}T_{22}\\
T_{12}T_{33}-T_{13}T_{23} & T_{33}T_{11}-T_{13}T_{13} & T_{11}T_{23}-T_{13}T_{12}\\
T_{12}T_{23}-T_{22}T_{13} & T_{23}T_{11}-T_{13}T_{12} & T_{11}T_{22}-T_{12}T_{12}
\end{array}
\right]
\end{equation*}
where $\T$ is a symmetric second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
SymTensor2 SymTensor2::minor() const
//-----------------------------------------------------------------------------
{
  return SymTensor2(_data[3] * _data[5] - _data[4] * _data[4],
                    _data[5] * _data[1] - _data[4] * _data[2],
                    _data[1] * _data[4] - _data[2] * _data[3],
                    _data[5] * _data[0] - _data[2] * _data[2],
                    _data[0] * _data[4] - _data[2] * _data[1],
                    _data[0] * _data[3] - _data[1] * _data[1]);
}

/*
@LABEL:SymTensor2::cofactors()
@SHORT:Cofactors of a symmetric second order tensor.
@RETURN:SymTensor2 & The cofactor of the symmetric second order tensor.
This method returns the cofactor of a symmetric second order tensor.
\begin{equation*}
T^{cof} = \left[\begin{array}{ccc}
T_{22}T_{33}-T_{23}T_{23} & T_{23}T_{13}-T_{33}T_{12} & T_{12}T_{23}-T_{13}T_{22}\\
T_{13}T_{23}-T_{12}T_{33} & T_{33}T_{11}-T_{13}T_{13} & T_{13}T_{12}-T_{11}T_{23}\\
T_{12}T_{23}-T_{22}T_{13} & T_{13}T_{12}-T_{23}T_{11} & T_{11}T_{22}-T_{12}T_{12}
\end{array}
\right]
\end{equation*}
where $\T$ is a symmetric second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
SymTensor2 SymTensor2::cofactors() const
//-----------------------------------------------------------------------------
{
  return SymTensor2(_data[3] * _data[5] - _data[4] * _data[4],
                    _data[4] * _data[2] - _data[5] * _data[1],
                    _data[1] * _data[4] - _data[2] * _data[3],
                    _data[5] * _data[0] - _data[2] * _data[2],
                    _data[2] * _data[1] - _data[0] * _data[4],
                    _data[0] * _data[3] - _data[1] * _data[1]);
}

/*
@LABEL:SymTensor2::solve(Vec3D x)
@SHORT:Solves a small linear system A.x = b.
@RETURN:Vec3D & The solution of the linear system.
This method returns the solution of a small linear system with the following form:
\begin{equation*}
\overrightarrow{y} = \A \cdot \overrightarrow{x}
\end{equation*}
where $\A$ is a symmetric second order tensor defined by the object itself and $\overrightarrow{x}$ is a vector defined by parameter x.
@END
*/
//-----------------------------------------------------------------------------
Vec3D SymTensor2::solve(const Vec3D &b) const
//-----------------------------------------------------------------------------
{
  return Vec3D(inverse() * b);
}

//-----------------------------------------------------------------------------
void SymTensor2::polarExtract(double eigenVectors[3][3], double eigenValues[3], SymTensor2 &U, Tensor2 &R) const
//-----------------------------------------------------------------------------
{
  double sq[3];

  // eigenVectors 1
  double U0[6];
  U0[0] = dnlSquare(eigenVectors[0][0]);
  U0[1] = eigenVectors[0][0] * eigenVectors[1][0];
  U0[2] = eigenVectors[0][0] * eigenVectors[2][0];
  U0[3] = dnlSquare(eigenVectors[1][0]);
  U0[4] = eigenVectors[1][0] * eigenVectors[2][0];
  U0[5] = dnlSquare(eigenVectors[2][0]);
  // eigenVectors 2
  double U1[6];
  U1[0] = dnlSquare(eigenVectors[0][1]);
  U1[1] = eigenVectors[0][1] * eigenVectors[1][1];
  U1[2] = eigenVectors[0][1] * eigenVectors[2][1];
  U1[3] = dnlSquare(eigenVectors[1][1]);
  U1[4] = eigenVectors[1][1] * eigenVectors[2][1];
  U1[5] = dnlSquare(eigenVectors[2][1]);
  // eigenVectors 3
  double U2[6];
  U2[0] = dnlSquare(eigenVectors[0][2]);
  U2[1] = eigenVectors[0][2] * eigenVectors[1][2];
  U2[2] = eigenVectors[0][2] * eigenVectors[2][2];
  U2[3] = dnlSquare(eigenVectors[1][2]);
  U2[4] = eigenVectors[1][2] * eigenVectors[2][2];
  U2[5] = dnlSquare(eigenVectors[2][2]);

  sq[0] = sqrt(eigenValues[0]);
  sq[1] = sqrt(eigenValues[1]);
  sq[2] = sqrt(eigenValues[2]);
  U._data[0] = sq[0] * U0[0] + sq[1] * U1[0] + sq[2] * U2[0];
  U._data[1] = sq[0] * U0[1] + sq[1] * U1[1] + sq[2] * U2[1];
  U._data[2] = sq[0] * U0[2] + sq[1] * U1[2] + sq[2] * U2[2];
  U._data[3] = sq[0] * U0[3] + sq[1] * U1[3] + sq[2] * U2[3];
  U._data[4] = sq[0] * U0[4] + sq[1] * U1[4] + sq[2] * U2[4];
  U._data[5] = sq[0] * U0[5] + sq[1] * U1[5] + sq[2] * U2[5];

  double Um1[6];
  double t1 = U._data[3] * U._data[5];
  double t2 = U._data[2] * U._data[4];
  double t4 = U._data[4] * U._data[4];
  double t5 = U._data[1] * U._data[5];
  double t6 = U._data[2] * U._data[3];

  double deter = U._data[0] * t1 + 2.0 * U._data[1] * t2 - U._data[0] * t4 - U._data[1] * t5 - U._data[2] * t6;

  Um1[0] = (t1 - t4) / deter;
  Um1[1] = (t2 - t5) / deter;
  Um1[2] = (U._data[1] * U._data[4] - t6) / deter;
  Um1[3] = (U._data[0] * U._data[5] - U._data[2] * U._data[2]) / deter;
  Um1[4] = (U._data[2] * U._data[1] - U._data[0] * U._data[4]) / deter;
  Um1[5] = (U._data[0] * U._data[3] - U._data[1] * U._data[1]) / deter;

  R._data[0] = _data[0] * Um1[0] + _data[1] * Um1[1] + _data[2] * Um1[2];
  R._data[1] = _data[0] * Um1[1] + _data[1] * Um1[3] + _data[2] * Um1[4];
  R._data[2] = _data[0] * Um1[2] + _data[1] * Um1[4] + _data[2] * Um1[5];
  R._data[3] = _data[3] * Um1[0] + _data[4] * Um1[1] + _data[5] * Um1[2];
  R._data[4] = _data[3] * Um1[1] + _data[4] * Um1[3] + _data[5] * Um1[4];
  R._data[5] = _data[3] * Um1[2] + _data[4] * Um1[4] + _data[5] * Um1[5];
  R._data[6] = _data[6] * Um1[0] + _data[7] * Um1[1] + _data[8] * Um1[2];
  R._data[7] = _data[6] * Um1[1] + _data[7] * Um1[3] + _data[8] * Um1[4];
  R._data[8] = _data[6] * Um1[2] + _data[7] * Um1[4] + _data[8] * Um1[5];
}

//-----------------------------------------------------------------------------
void SymTensor2::polarExtractLnU(double eigenVectors[3][3], double eigenValues[3], SymTensor2 &U, Tensor2 &R) const
//-----------------------------------------------------------------------------
{
  double sq[3];

  // eigenVectors 1
  double U0[6];
  U0[0] = dnlSquare(eigenVectors[0][0]);
  U0[1] = eigenVectors[0][0] * eigenVectors[1][0];
  U0[2] = eigenVectors[0][0] * eigenVectors[2][0];
  U0[3] = dnlSquare(eigenVectors[1][0]);
  U0[4] = eigenVectors[1][0] * eigenVectors[2][0];
  U0[5] = dnlSquare(eigenVectors[2][0]);
  // eigenVectors 2
  double U1[6];
  U1[0] = dnlSquare(eigenVectors[0][1]);
  U1[1] = eigenVectors[0][1] * eigenVectors[1][1];
  U1[2] = eigenVectors[0][1] * eigenVectors[2][1];
  U1[3] = dnlSquare(eigenVectors[1][1]);
  U1[4] = eigenVectors[1][1] * eigenVectors[2][1];
  U1[5] = dnlSquare(eigenVectors[2][1]);
  // eigenVectors 3
  double U2[6];
  U2[0] = dnlSquare(eigenVectors[0][2]);
  U2[1] = eigenVectors[0][2] * eigenVectors[1][2];
  U2[2] = eigenVectors[0][2] * eigenVectors[2][2];
  U2[3] = dnlSquare(eigenVectors[1][2]);
  U2[4] = eigenVectors[1][2] * eigenVectors[2][2];
  U2[5] = dnlSquare(eigenVectors[2][2]);

  sq[0] = sqrt(eigenValues[0]);
  sq[1] = sqrt(eigenValues[1]);
  sq[2] = sqrt(eigenValues[2]);
  U._data[0] = sq[0] * U0[0] + sq[1] * U1[0] + sq[2] * U2[0];
  U._data[1] = sq[0] * U0[1] + sq[1] * U1[1] + sq[2] * U2[1];
  U._data[2] = sq[0] * U0[2] + sq[1] * U1[2] + sq[2] * U2[2];
  U._data[3] = sq[0] * U0[3] + sq[1] * U1[3] + sq[2] * U2[3];
  U._data[4] = sq[0] * U0[4] + sq[1] * U1[4] + sq[2] * U2[4];
  U._data[5] = sq[0] * U0[5] + sq[1] * U1[5] + sq[2] * U2[5];

  double Um1[6];
  double t1 = U._data[3] * U._data[5];
  double t2 = U._data[2] * U._data[4];
  double t4 = U._data[4] * U._data[4];
  double t5 = U._data[1] * U._data[5];
  double t6 = U._data[2] * U._data[3];

  double deter = U._data[0] * t1 + 2.0 * U._data[1] * t2 - U._data[0] * t4 - U._data[1] * t5 - U._data[2] * t6;

  Um1[0] = (t1 - t4) / deter;
  Um1[1] = (t2 - t5) / deter;
  Um1[2] = (U._data[1] * U._data[4] - t6) / deter;
  Um1[3] = (U._data[0] * U._data[5] - U._data[2] * U._data[2]) / deter;
  Um1[4] = (U._data[2] * U._data[1] - U._data[0] * U._data[4]) / deter;
  Um1[5] = (U._data[0] * U._data[3] - U._data[1] * U._data[1]) / deter;

  R._data[0] = _data[0] * Um1[0] + _data[1] * Um1[1] + _data[2] * Um1[2];
  R._data[1] = _data[0] * Um1[1] + _data[1] * Um1[3] + _data[2] * Um1[4];
  R._data[2] = _data[0] * Um1[2] + _data[1] * Um1[4] + _data[2] * Um1[5];
  R._data[3] = _data[3] * Um1[0] + _data[4] * Um1[1] + _data[5] * Um1[2];
  R._data[4] = _data[3] * Um1[1] + _data[4] * Um1[3] + _data[5] * Um1[4];
  R._data[5] = _data[3] * Um1[2] + _data[4] * Um1[4] + _data[5] * Um1[5];
  R._data[6] = _data[6] * Um1[0] + _data[7] * Um1[1] + _data[8] * Um1[2];
  R._data[7] = _data[6] * Um1[1] + _data[7] * Um1[3] + _data[8] * Um1[4];
  R._data[8] = _data[6] * Um1[2] + _data[7] * Um1[4] + _data[8] * Um1[5];

  sq[0] = log(eigenValues[0]) / 2;
  sq[1] = log(eigenValues[1]) / 2;
  sq[2] = log(eigenValues[2]) / 2;
  U._data[0] = sq[0] * U0[0] + sq[1] * U1[0] + sq[2] * U2[0];
  U._data[1] = sq[0] * U0[1] + sq[1] * U1[1] + sq[2] * U2[1];
  U._data[2] = sq[0] * U0[2] + sq[1] * U1[2] + sq[2] * U2[2];
  U._data[3] = sq[0] * U0[3] + sq[1] * U1[3] + sq[2] * U2[3];
  U._data[4] = sq[0] * U0[4] + sq[1] * U1[4] + sq[2] * U2[4];
  U._data[5] = sq[0] * U0[5] + sq[1] * U1[5] + sq[2] * U2[5];
}

//-----------------------------------------------------------------------------
void SymTensor2::buildFTF(double FTF[3][3]) const
//-----------------------------------------------------------------------------
{
  FTF[0][0] = dnlSquare(_data[0]) + dnlSquare(_data[1]) + dnlSquare(_data[2]);
  FTF[0][1] = _data[0] * _data[1] + _data[1] * _data[3] + _data[2] * _data[4];
  FTF[0][2] = _data[0] * _data[2] + _data[1] * _data[4] + _data[2] * _data[5];
  FTF[1][0] = FTF[0][1];
  FTF[1][1] = dnlSquare(_data[1]) + dnlSquare(_data[3]) + dnlSquare(_data[4]);
  FTF[1][2] = _data[1] * _data[2] + _data[3] * _data[4] + _data[4] * _data[5];
  FTF[2][0] = FTF[0][2];
  FTF[2][1] = FTF[1][2];
  FTF[2][2] = dnlSquare(_data[2]) + dnlSquare(_data[4]) + dnlSquare(_data[5]);
}

/*
@LABEL:SymTensor2::polarQL(SymTensor2 U, Tensor2 R)
@SHORT:Polar decomposition of a second order tensor using the QL with implicit shifts algorithm.
@RETURN:SymTensor2 and Tensor2
@ARG:SymTensor2&U&Symmetric tensor $\U$
@ARG:Tensor2&R&Rotation tensor $\R$
This method computes the polar decomposition of a second order tensor $\F$ and returns the symmetric tensor $\R$ and the tensor $\U$ so that:
\begin{equation*}
\F = \R \cdot \U
\end{equation*}
where $\F$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
void SymTensor2::polarQL(SymTensor2 &U, Tensor2 &R) const
//-----------------------------------------------------------------------------
{
  double FTF[3][3];
  double eigenVectors[3][3];
  double eigenValues[3];

  // Build the F(T).F symmetric matrix
  buildFTF(FTF);

  // Compute the eigenvalues and eigenvectors
  dsyevq3(FTF, eigenVectors, eigenValues); // QL with implicit shifts

  // Extract the tensors for U and R
  polarExtract(eigenVectors, eigenValues, U, R);
}

/*
@LABEL:SymTensor2::polarQLLnU(SymTensor2 U, Tensor2 R)
@SHORT:Polar decomposition of a second order tensor using the QL with implicit shifts algorithm.
@RETURN:SymTensor2 and Tensor2
@ARG:SymTensor2&U&Symmetric tensor $\log[\U]$
@ARG:Tensor2&R&Rotation tensor $\R$
This method computes the polar decomposition of a second order tensor $\F$ and returns the symmetric tensor $\R$ and the tensor $\log[\U]$ so that:
\begin{equation*}
\F = \R \cdot \U
\end{equation*}
\begin{equation*}
\log [\U] =\sum _{i=1}^{3}\log[\lambda_{i}](\overrightarrow{u}_{i}\otimes \overrightarrow{u}_{i})
\end{equation*}
where $\F$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
void SymTensor2::polarQLLnU(SymTensor2 &U, Tensor2 &R) const
//-----------------------------------------------------------------------------
{
  double FTF[3][3];
  double eigenVectors[3][3];
  double eigenValues[3];

  // Build the F(T).F symmetric matrix
  buildFTF(FTF);

  // Compute the eigenvalues and eigenvectors
  dsyevq3(FTF, eigenVectors, eigenValues); // QL with implicit shifts

  // Extract the tensors for U and R
  polarExtractLnU(eigenVectors, eigenValues, U, R);
}

/*
@LABEL:SymTensor2::polarCuppen(SymTensor2 U, Tensor2 R)
@SHORT:Polar decomposition of a second order tensor using the Cuppen’s Divide and Conquer algorithm.
@RETURN:SymTensor2 and Tensor2
@ARG:SymTensor2&U&Symmetric tensor $\log[\U]$
@ARG:Tensor2&R&Rotation tensor $\R$
This method computes the polar decomposition of a second order tensor $\F$ and returns the symmetric tensor $\R$ and the tensor $\U$ so that:
\begin{equation*}
\F = \R \cdot \U
\end{equation*}
where $\F$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
void SymTensor2::polarCuppen(SymTensor2 &U, Tensor2 &R) const
//-----------------------------------------------------------------------------
{
  double FTF[3][3];
  double eigenVectors[3][3];
  double eigenValues[3];

  // Build the F(T).F symmetric matrix
  buildFTF(FTF);

  // Compute the eigenvalues and eigenvectors
  dsyevd3(FTF, eigenVectors, eigenValues); // Cuppen

  // Extract the tensors for U and R
  polarExtract(eigenVectors, eigenValues, U, R);
}

/*
@LABEL:SymTensor2::polarCuppenLnU(SymTensor2 U, Tensor2 R)
@SHORT:Polar decomposition of a second order tensor using the Cuppen’s Divide and Conquer algorithm.
@RETURN:SymTensor2 and Tensor2
@ARG:SymTensor2&U&Symmetric tensor $\log[\U]$
@ARG:Tensor2&R&Rotation tensor $\R$
This method computes the polar decomposition of a second order tensor $\F$ and returns the symmetric tensor $\R$ and the tensor $\log[\U]$ so that:
\begin{equation*}
\F = \R \cdot \U
\end{equation*}
\begin{equation*}
\log [\U] =\sum _{i=1}^{3}\log[\lambda_{i}](\overrightarrow{u}_{i}\otimes \overrightarrow{u}_{i})
\end{equation*}
where $\F$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
void SymTensor2::polarCuppenLnU(SymTensor2 &U, Tensor2 &R) const
//-----------------------------------------------------------------------------
{
  double FTF[3][3];
  double eigenVectors[3][3];
  double eigenValues[3];

  // Build the F(T).F symmetric matrix
  buildFTF(FTF);

  // Compute the eigenvalues and eigenvectors
  dsyevd3(FTF, eigenVectors, eigenValues); // Cuppen

  // Extract the tensors for U and R
  polarExtractLnU(eigenVectors, eigenValues, U, R);
}

/*
@LABEL:SymTensor2::polarJacobi(SymTensor2 U, Tensor2 R)
@SHORT:Polar decomposition of a second order tensor using the Jacobi algorithm.
@RETURN:SymTensor2 and Tensor2
@ARG:SymTensor2&U&Symmetric tensor $\log[\U]$
@ARG:Tensor2&R&Rotation tensor $\R$
This method computes the polar decomposition of a second order tensor $\F$ and returns the symmetric tensor $\R$ and the tensor $\U$ so that:
\begin{equation*}
\F = \R \cdot \U
\end{equation*}
where $\F$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
void SymTensor2::polarJacobi(SymTensor2 &U, Tensor2 &R) const
//-----------------------------------------------------------------------------
{
  double FTF[3][3];
  double eigenVectors[3][3];
  double eigenValues[3];

  // Build the F(T).F symmetric matrix
  buildFTF(FTF);

  // Compute the eigenvalues and eigenvectors
  dsyevj3(FTF, eigenVectors, eigenValues); // Jacobi

  // Extract the tensors for U and R
  polarExtract(eigenVectors, eigenValues, U, R);
}

/*
@LABEL:SymTensor2::polarJacobiLnU(SymTensor2 U, Tensor2 R)
@SHORT:Polar decomposition of a second order tensor using the Jacobi algorithm.
@RETURN:SymTensor2 and Tensor2
@ARG:SymTensor2&U&Symmetric tensor $\log[\U]$
@ARG:Tensor2&R&Rotation tensor $\R$
This method computes the polar decomposition of a second order tensor $\F$ and returns the symmetric tensor $\R$ and the tensor $\log[\U]$ so that:
\begin{equation*}
\F = \R \cdot \U
\end{equation*}
\begin{equation*}
\log [\U] =\sum _{i=1}^{3}\log[\lambda_{i}](\overrightarrow{u}_{i}\otimes \overrightarrow{u}_{i})
\end{equation*}
where $\F$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
void SymTensor2::polarJacobiLnU(SymTensor2 &U, Tensor2 &R) const
//-----------------------------------------------------------------------------
{
  double FTF[3][3];
  double eigenVectors[3][3];
  double eigenValues[3];

  // Build the F(T).F symmetric matrix
  buildFTF(FTF);

  // Compute the eigenvalues and eigenvectors
  dsyevj3(FTF, eigenVectors, eigenValues); // Jacobi

  // Extract the tensors for U and R
  polarExtractLnU(eigenVectors, eigenValues, U, R);
}

/*
@LABEL:SymTensor2::polarLapack(SymTensor2 U, Tensor2 R)
@SHORT:Polar decomposition of a second order tensor using the Jacobi algorithm.
@RETURN:SymTensor2 and Tensor2
@ARG:SymTensor2&U&Symmetric tensor $\log[\U]$
@ARG:Tensor2&R&Rotation tensor $\R$
This method computes the polar decomposition of a second order tensor $\F$ and returns the symmetric tensor $\R$ and the tensor $\U$ so that:
\begin{equation*}
\F = \R \cdot \U
\end{equation*}
where $\F$ is a second order tensor defined by the object itself.
It uses the LAPACKE\_dgeev function of the Lapack library which is far from efficient for a trivial 3x3 matrix. So this method is very slow.
@END
*/
//-----------------------------------------------------------------------------
void SymTensor2::polarLapack(SymTensor2 &U, Tensor2 &R) const
//-----------------------------------------------------------------------------
{
  double FTF[9];

  FTF[0] = dnlSquare(_data[0]) + dnlSquare(_data[1]) + dnlSquare(_data[2]);
  FTF[1] = _data[0] * _data[1] + _data[1] * _data[3] + _data[2] * _data[4];
  FTF[2] = _data[0] * _data[2] + _data[1] * _data[4] + _data[2] * _data[5];
  FTF[3] = FTF[1];
  FTF[4] = dnlSquare(_data[1]) + dnlSquare(_data[3]) + dnlSquare(_data[4]);
  FTF[5] = _data[1] * _data[2] + _data[3] * _data[4] + _data[4] * _data[5];
  FTF[6] = FTF[2];
  FTF[7] = FTF[5];
  FTF[8] = dnlSquare(_data[2]) + dnlSquare(_data[4]) + dnlSquare(_data[5]);

  double eigenVectors[9];
  double tmpMatrix[9];
  double eigenValues[3];
  double tmpVector[3];
  double sq[3];

  int info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', 3, FTF, 3, eigenValues, tmpVector, tmpMatrix, 3, eigenVectors, 3);

  // eigenVectors 1
  double U0[6];
  U0[0] = dnlSquare(eigenVectors[0]);
  U0[1] = eigenVectors[0] * eigenVectors[3];
  U0[2] = eigenVectors[0] * eigenVectors[6];
  U0[3] = dnlSquare(eigenVectors[3]);
  U0[4] = eigenVectors[3] * eigenVectors[6];
  U0[5] = dnlSquare(eigenVectors[6]);
  // eigenVectors 2
  double U1[6];
  U1[0] = dnlSquare(eigenVectors[1]);
  U1[1] = eigenVectors[1] * eigenVectors[4];
  U1[2] = eigenVectors[1] * eigenVectors[7];
  U1[3] = dnlSquare(eigenVectors[4]);
  U1[4] = eigenVectors[4] * eigenVectors[7];
  U1[5] = dnlSquare(eigenVectors[7]);
  // eigenVectors 3
  double U2[6];
  U2[0] = dnlSquare(eigenVectors[2]);
  U2[1] = eigenVectors[2] * eigenVectors[5];
  U2[2] = eigenVectors[2] * eigenVectors[8];
  U2[3] = dnlSquare(eigenVectors[5]);
  U2[4] = eigenVectors[5] * eigenVectors[8];
  U2[5] = dnlSquare(eigenVectors[8]);

  sq[0] = sqrt(eigenValues[0]);
  sq[1] = sqrt(eigenValues[1]);
  sq[2] = sqrt(eigenValues[2]);
  U._data[0] = sq[0] * U0[0] + sq[1] * U1[0] + sq[2] * U2[0];
  U._data[1] = sq[0] * U0[1] + sq[1] * U1[1] + sq[2] * U2[1];
  U._data[2] = sq[0] * U0[2] + sq[1] * U1[2] + sq[2] * U2[2];
  U._data[3] = sq[0] * U0[3] + sq[1] * U1[3] + sq[2] * U2[3];
  U._data[4] = sq[0] * U0[4] + sq[1] * U1[4] + sq[2] * U2[4];
  U._data[5] = sq[0] * U0[5] + sq[1] * U1[5] + sq[2] * U2[5];

  double Um1[6];
  double t1 = U._data[3] * U._data[5];
  double t2 = U._data[2] * U._data[4];
  double t4 = U._data[4] * U._data[4];
  double t5 = U._data[1] * U._data[5];
  double t6 = U._data[2] * U._data[3];

  double deter = U._data[0] * t1 + 2.0 * U._data[1] * t2 - U._data[0] * t4 - U._data[1] * t5 - U._data[2] * t6;
  Um1[0] = (t1 - t4) / deter;
  Um1[1] = (t2 - t5) / deter;
  Um1[2] = (U._data[1] * U._data[4] - t6) / deter;
  Um1[3] = (U._data[0] * U._data[5] - U._data[2] * U._data[2]) / deter;
  Um1[4] = (U._data[2] * U._data[1] - U._data[0] * U._data[4]) / deter;
  Um1[5] = (U._data[0] * U._data[3] - U._data[1] * U._data[1]) / deter;

  R._data[0] = _data[0] * Um1[0] + _data[1] * Um1[1] + _data[2] * Um1[2];
  R._data[1] = _data[0] * Um1[1] + _data[1] * Um1[3] + _data[2] * Um1[4];
  R._data[2] = _data[0] * Um1[2] + _data[1] * Um1[4] + _data[2] * Um1[5];
  R._data[3] = _data[3] * Um1[0] + _data[4] * Um1[1] + _data[5] * Um1[2];
  R._data[4] = _data[3] * Um1[1] + _data[4] * Um1[3] + _data[5] * Um1[4];
  R._data[5] = _data[3] * Um1[2] + _data[4] * Um1[4] + _data[5] * Um1[5];
  R._data[6] = _data[6] * Um1[0] + _data[7] * Um1[1] + _data[8] * Um1[2];
  R._data[7] = _data[6] * Um1[1] + _data[7] * Um1[3] + _data[8] * Um1[4];
  R._data[8] = _data[6] * Um1[2] + _data[7] * Um1[4] + _data[8] * Um1[5];
}

/*
@LABEL:SymTensor2::polarLapackLnU(SymTensor2 U, Tensor2 R)
@SHORT:Polar decomposition of a second order tensor using the Jacobi algorithm.
@RETURN:SymTensor2 and Tensor2
@ARG:SymTensor2&U&Symmetric tensor $\log[\U]$
@ARG:Tensor2&R&Rotation tensor $\R$
This method computes the polar decomposition of a second order tensor $\F$ and returns the symmetric tensor $\R$ and the tensor $\log[\U]$ so that:
\begin{equation*}
\F = \R \cdot \U
\end{equation*}
\begin{equation*}
\log [\U] =\sum _{i=1}^{3}\log[\lambda_{i}](\overrightarrow{u}_{i}\otimes \overrightarrow{u}_{i})
\end{equation*}
where $\F$ is a second order tensor defined by the object itself.
It uses the LAPACKE\_dgeev function of the Lapack library which is far from efficient for a trivial 3x3 matrix. So this method is very slow.
@END
*/
//-----------------------------------------------------------------------------
void SymTensor2::polarLapackLnU(SymTensor2 &U, Tensor2 &R) const
//-----------------------------------------------------------------------------
{
  double FTF[9];

  FTF[0] = dnlSquare(_data[0]) + dnlSquare(_data[1]) + dnlSquare(_data[2]);
  FTF[1] = _data[0] * _data[1] + _data[1] * _data[3] + _data[2] * _data[4];
  FTF[2] = _data[0] * _data[2] + _data[1] * _data[4] + _data[2] * _data[5];
  FTF[3] = FTF[1];
  FTF[4] = dnlSquare(_data[1]) + dnlSquare(_data[3]) + dnlSquare(_data[4]);
  FTF[5] = _data[1] * _data[2] + _data[3] * _data[4] + _data[4] * _data[5];
  FTF[6] = FTF[2];
  FTF[7] = FTF[5];
  FTF[8] = dnlSquare(_data[2]) + dnlSquare(_data[4]) + dnlSquare(_data[5]);

  double eigenVectors[9];
  double tmpMatrix[9];
  double eigenValues[3];
  double tmpVector[3];
  double sq[3];

  int info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', 3, FTF, 3, eigenValues, tmpVector, tmpMatrix, 3, eigenVectors, 3);

  // eigenVectors 1
  double U0[6];
  U0[0] = dnlSquare(eigenVectors[0]);
  U0[1] = eigenVectors[0] * eigenVectors[3];
  U0[2] = eigenVectors[0] * eigenVectors[6];
  U0[3] = dnlSquare(eigenVectors[3]);
  U0[4] = eigenVectors[3] * eigenVectors[6];
  U0[5] = dnlSquare(eigenVectors[6]);
  // eigenVectors 2
  double U1[6];
  U1[0] = dnlSquare(eigenVectors[1]);
  U1[1] = eigenVectors[1] * eigenVectors[4];
  U1[2] = eigenVectors[1] * eigenVectors[7];
  U1[3] = dnlSquare(eigenVectors[4]);
  U1[4] = eigenVectors[4] * eigenVectors[7];
  U1[5] = dnlSquare(eigenVectors[7]);
  // eigenVectors 3
  double U2[6];
  U2[0] = dnlSquare(eigenVectors[2]);
  U2[1] = eigenVectors[2] * eigenVectors[5];
  U2[2] = eigenVectors[2] * eigenVectors[8];
  U2[3] = dnlSquare(eigenVectors[5]);
  U2[4] = eigenVectors[5] * eigenVectors[8];
  U2[5] = dnlSquare(eigenVectors[8]);

  sq[0] = sqrt(eigenValues[0]);
  sq[1] = sqrt(eigenValues[1]);
  sq[2] = sqrt(eigenValues[2]);
  U._data[0] = sq[0] * U0[0] + sq[1] * U1[0] + sq[2] * U2[0];
  U._data[1] = sq[0] * U0[1] + sq[1] * U1[1] + sq[2] * U2[1];
  U._data[2] = sq[0] * U0[2] + sq[1] * U1[2] + sq[2] * U2[2];
  U._data[3] = sq[0] * U0[3] + sq[1] * U1[3] + sq[2] * U2[3];
  U._data[4] = sq[0] * U0[4] + sq[1] * U1[4] + sq[2] * U2[4];
  U._data[5] = sq[0] * U0[5] + sq[1] * U1[5] + sq[2] * U2[5];

  double Um1[6];
  double t1 = U._data[3] * U._data[5];
  double t2 = U._data[2] * U._data[4];
  double t4 = U._data[4] * U._data[4];
  double t5 = U._data[1] * U._data[5];
  double t6 = U._data[2] * U._data[3];

  double deter = U._data[0] * t1 + 2.0 * U._data[1] * t2 - U._data[0] * t4 - U._data[1] * t5 - U._data[2] * t6;
  Um1[0] = (t1 - t4) / deter;
  Um1[1] = (t2 - t5) / deter;
  Um1[2] = (U._data[1] * U._data[4] - t6) / deter;
  Um1[3] = (U._data[0] * U._data[5] - U._data[2] * U._data[2]) / deter;
  Um1[4] = (U._data[2] * U._data[1] - U._data[0] * U._data[4]) / deter;
  Um1[5] = (U._data[0] * U._data[3] - U._data[1] * U._data[1]) / deter;

  R._data[0] = _data[0] * Um1[0] + _data[1] * Um1[1] + _data[2] * Um1[2];
  R._data[1] = _data[0] * Um1[1] + _data[1] * Um1[3] + _data[2] * Um1[4];
  R._data[2] = _data[0] * Um1[2] + _data[1] * Um1[4] + _data[2] * Um1[5];
  R._data[3] = _data[3] * Um1[0] + _data[4] * Um1[1] + _data[5] * Um1[2];
  R._data[4] = _data[3] * Um1[1] + _data[4] * Um1[3] + _data[5] * Um1[4];
  R._data[5] = _data[3] * Um1[2] + _data[4] * Um1[4] + _data[5] * Um1[5];
  R._data[6] = _data[6] * Um1[0] + _data[7] * Um1[1] + _data[8] * Um1[2];
  R._data[7] = _data[6] * Um1[1] + _data[7] * Um1[3] + _data[8] * Um1[4];
  R._data[8] = _data[6] * Um1[2] + _data[7] * Um1[4] + _data[8] * Um1[5];

  sq[0] = log(eigenValues[0]) / 2;
  sq[1] = log(eigenValues[1]) / 2;
  sq[2] = log(eigenValues[2]) / 2;
  U._data[0] = sq[0] * U0[0] + sq[1] * U1[0] + sq[2] * U2[0];
  U._data[1] = sq[0] * U0[1] + sq[1] * U1[1] + sq[2] * U2[1];
  U._data[2] = sq[0] * U0[2] + sq[1] * U1[2] + sq[2] * U2[2];
  U._data[3] = sq[0] * U0[3] + sq[1] * U1[3] + sq[2] * U2[3];
  U._data[4] = sq[0] * U0[4] + sq[1] * U1[4] + sq[2] * U2[4];
  U._data[5] = sq[0] * U0[5] + sq[1] * U1[5] + sq[2] * U2[5];
}

/*
@LABEL:SymTensor2::polarLnU(SymTensor2 U, Tensor2 R)
@SHORT:Polar decomposition of a second order tensor using the old \DynELA algorithm.
@RETURN:SymTensor2 and Tensor2
@ARG:SymTensor2&U&Symmetric tensor $\log[\U]$
@ARG:Tensor2&R&Rotation tensor $\R$
This method computes the polar decomposition of a second order tensor $\F$ and returns the symmetric tensor $\R$ and the tensor $\log[\U]$ so that:
\begin{equation*}
\F = \R \cdot \U
\end{equation*}
\begin{equation*}
\log [\U] =\sum _{i=1}^{3}\log[\lambda_{i}](\overrightarrow{u}_{i}\otimes \overrightarrow{u}_{i})
\end{equation*}
where $\F$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
void SymTensor2::polarLnU(SymTensor2 &LnU, Tensor2 &R) const
//-----------------------------------------------------------------------------
{
  short p, q;
  double alpha, t, Ctheta, Stheta;
  double _tmp1, _tmp2, _tmp3, _tmp4, _tmp5, _tmp6, _tmp7, _tmp8;
  double Rot11, Rot12, Rot13, Rot21, Rot22, Rot23, Rot31, Rot32, Rot33;
  SymTensor2 FF;
  short it = 1;

  // calcul de FF=F(T).F
  FF._data[0] = dnlSquare(_data[0]) + dnlSquare(_data[1]) + dnlSquare(_data[2]);
  FF._data[1] = _data[0] * _data[1] + _data[1] * _data[3] + _data[2] * _data[4];
  FF._data[2] = _data[0] * _data[2] + _data[1] * _data[4] + _data[2] * _data[5];
  FF._data[3] = dnlSquare(_data[1]) + dnlSquare(_data[3]) + dnlSquare(_data[4]);
  FF._data[4] = _data[1] * _data[2] + _data[3] * _data[4] + _data[4] * _data[5];
  FF._data[5] = dnlSquare(_data[2]) + dnlSquare(_data[4]) + dnlSquare(_data[5]);

  // initialisation de Rot
  Rot11 = 1;
  Rot22 = 1;
  Rot33 = 1;
  Rot12 = 0;
  Rot13 = 0;
  Rot21 = 0;
  Rot23 = 0;
  Rot31 = 0;
  Rot32 = 0;

  while (it < 30)
  {
    // recherche du terme maxi
    q = (dnlAbs(FF._data[1]) > dnlAbs(FF._data[2]) ? 1 : 2);
    p = (dnlAbs(FF._data[q]) > dnlAbs(FF._data[4]) ? 0 : 1);
    if (p == 1)
      q = 2;

    // test de convergence
    if (dnlAbs(FF._data[_internalIndexes.index[p][q]]) < precisionPolarDecompose)
    {
      _tmp1 = 0.5 * log(FF._data[0]);
      _tmp2 = 0.5 * log(FF._data[3]);
      _tmp3 = 0.5 * log(FF._data[5]);
      LnU._data[0] = dnlSquare(Rot11) * _tmp1 + dnlSquare(Rot12) * _tmp2 + dnlSquare(Rot13) * _tmp3;
      LnU._data[3] = dnlSquare(Rot21) * _tmp1 + dnlSquare(Rot22) * _tmp2 + dnlSquare(Rot23) * _tmp3;
      LnU._data[5] = dnlSquare(Rot31) * _tmp1 + dnlSquare(Rot32) * _tmp2 + dnlSquare(Rot33) * _tmp3;
      LnU._data[1] = (Rot11 * Rot21 * _tmp1 + Rot12 * Rot22 * _tmp2 + Rot13 * Rot23 * _tmp3);
      LnU._data[2] = (Rot11 * Rot31 * _tmp1 + Rot12 * Rot32 * _tmp2 + Rot13 * Rot33 * _tmp3);
      LnU._data[4] = (Rot21 * Rot31 * _tmp1 + Rot22 * Rot32 * _tmp2 + Rot23 * Rot33 * _tmp3);

      _tmp1 = 1 / sqrt(FF._data[0]);
      _tmp2 = 1 / sqrt(FF._data[3]);
      _tmp3 = 1 / sqrt(FF._data[5]);
      _tmp4 = (Rot13 * _data[0] + Rot23 * _data[1] + Rot33 * _data[2]) * _tmp3;
      _tmp5 = (Rot12 * _data[0] + Rot22 * _data[1] + Rot32 * _data[2]) * _tmp2;
      _tmp6 = (Rot11 * _data[0] + Rot21 * _data[1] + Rot31 * _data[2]) * _tmp1;
      R._data[0] = Rot13 * _tmp4 + Rot12 * _tmp5 + Rot11 * _tmp6;
      R._data[1] = Rot23 * _tmp4 + Rot22 * _tmp5 + Rot21 * _tmp6;
      R._data[2] = Rot33 * _tmp4 + Rot32 * _tmp5 + Rot31 * _tmp6;
      _tmp4 = (Rot13 * _data[1] + Rot23 * _data[3] + Rot33 * _data[4]) * _tmp3;
      _tmp5 = (Rot12 * _data[1] + Rot22 * _data[3] + Rot32 * _data[4]) * _tmp2;
      _tmp6 = (Rot11 * _data[1] + Rot21 * _data[3] + Rot31 * _data[4]) * _tmp1;
      R._data[3] = Rot13 * _tmp4 + Rot12 * _tmp5 + Rot11 * _tmp6;
      R._data[4] = Rot23 * _tmp4 + Rot22 * _tmp5 + Rot21 * _tmp6;
      R._data[5] = Rot33 * _tmp4 + Rot32 * _tmp5 + Rot31 * _tmp6;
      _tmp4 = (Rot13 * _data[2] + Rot23 * _data[4] + Rot33 * _data[5]) * _tmp3;
      _tmp5 = (Rot12 * _data[2] + Rot22 * _data[4] + Rot32 * _data[5]) * _tmp2;
      _tmp6 = (Rot11 * _data[2] + Rot21 * _data[4] + Rot31 * _data[5]) * _tmp1;
      R._data[6] = Rot13 * _tmp4 + Rot12 * _tmp5 + Rot11 * _tmp6;
      R._data[7] = Rot23 * _tmp4 + Rot22 * _tmp5 + Rot21 * _tmp6;
      R._data[8] = Rot33 * _tmp4 + Rot32 * _tmp5 + Rot31 * _tmp6;
      return;
    }

    // calcul de l'angle
    alpha = (FF._data[_internalIndexes.index[q][q]] - FF._data[_internalIndexes.index[p][p]]) / (2. * FF._data[_internalIndexes.index[p][q]]);
    t = (alpha != 0 ? 1 / (alpha + sqrt(alpha * alpha + 1) * (alpha > 0 ? 1 : -1)) : 1);

    // calcul des cosinus
    Ctheta = 1 / (sqrt(1 + t * t));
    Stheta = t * Ctheta;

    // calcul du nouveau tenseur
    R.setToUnity();
    R(p, q) = Stheta;
    R(q, p) = -Stheta;
    R(p, p) = R(q, q) = Ctheta;

    // produit R(T).FF.R avec FF symétrique
    _tmp1 = FF._data[0] * R._data[0] + FF._data[1] * R._data[3] + FF._data[2] * R._data[6];
    _tmp2 = FF._data[1] * R._data[0] + FF._data[3] * R._data[3] + FF._data[4] * R._data[6];
    _tmp3 = FF._data[2] * R._data[0] + FF._data[4] * R._data[3] + FF._data[5] * R._data[6];
    _tmp4 = FF._data[0] * R._data[1] + FF._data[1] * R._data[4] + FF._data[2] * R._data[7];
    _tmp5 = FF._data[1] * R._data[1] + FF._data[3] * R._data[4] + FF._data[4] * R._data[7];
    _tmp6 = FF._data[2] * R._data[1] + FF._data[4] * R._data[4] + FF._data[5] * R._data[7];
    FF._data[5] = R._data[2] * (FF._data[0] * R._data[2] + FF._data[1] * R._data[5] + FF._data[2] * R._data[8]) +
                  R._data[5] * (FF._data[1] * R._data[2] + FF._data[3] * R._data[5] + FF._data[4] * R._data[8]) +
                  R._data[8] * (FF._data[2] * R._data[2] + FF._data[4] * R._data[5] + FF._data[5] * R._data[8]);
    FF._data[0] = R._data[0] * _tmp1 + R._data[3] * _tmp2 + R._data[6] * _tmp3;
    FF._data[1] = R._data[1] * _tmp1 + R._data[4] * _tmp2 + R._data[7] * _tmp3;
    FF._data[2] = R._data[2] * _tmp1 + R._data[5] * _tmp2 + R._data[8] * _tmp3;
    FF._data[3] = R._data[1] * _tmp4 + R._data[4] * _tmp5 + R._data[7] * _tmp6;
    FF._data[4] = R._data[2] * _tmp4 + R._data[5] * _tmp5 + R._data[8] * _tmp6;

    // cumulate of R
    _tmp1 = R._data[0] * Rot11 + R._data[3] * Rot12 + R._data[6] * Rot13;
    _tmp2 = R._data[1] * Rot11 + R._data[4] * Rot12 + R._data[7] * Rot13;
    Rot13 = R._data[2] * Rot11 + R._data[5] * Rot12 + R._data[8] * Rot13;
    _tmp4 = R._data[0] * Rot21 + R._data[3] * Rot22 + R._data[6] * Rot23;
    _tmp5 = R._data[1] * Rot21 + R._data[4] * Rot22 + R._data[7] * Rot23;
    Rot23 = R._data[2] * Rot21 + R._data[5] * Rot22 + R._data[8] * Rot23;
    _tmp7 = R._data[0] * Rot31 + R._data[3] * Rot32 + R._data[6] * Rot33;
    _tmp8 = R._data[1] * Rot31 + R._data[4] * Rot32 + R._data[7] * Rot33;
    Rot33 = R._data[2] * Rot31 + R._data[5] * Rot32 + R._data[8] * Rot33;
    Rot11 = _tmp1;
    Rot12 = _tmp2;
    Rot21 = _tmp4;
    Rot22 = _tmp5;
    Rot31 = _tmp7;
    Rot32 = _tmp8;
    it++;
  }

  std::cout << (*this) << std::endl;
  fatalError("SymTensor2::polarLnU", "No convergence");
}

/*
@LABEL:SymTensor2::polar(SymTensor2 U, Tensor2 R)
@SHORT:Polar decomposition of a second order tensor using the old \DynELA algorithm.
@RETURN:SymTensor2 and Tensor2
@ARG:SymTensor2&U&Symmetric tensor $\log[\U]$
@ARG:Tensor2&R&Rotation tensor $\R$
This method computes the polar decomposition of a second order tensor $\F$ and returns the symmetric tensor $\R$ and the tensor $\U$ so that:
\begin{equation*}
\F = \R \cdot \U
\end{equation*}
where $\F$ is a second order tensor defined by the object itself.
@END
*/
//-----------------------------------------------------------------------------
void SymTensor2::polar(SymTensor2 &U, Tensor2 &R) const
//-----------------------------------------------------------------------------
{
  short p, q;
  double alpha, t, Ctheta, Stheta;
  double _tmp1, _tmp2, _tmp3, _tmp4, _tmp5, _tmp6, _tmp7, _tmp8;
  double Rot11, Rot12, Rot13, Rot21, Rot22, Rot23, Rot31, Rot32, Rot33;
  SymTensor2 FF;
  short it = 1;

  // calcul de FF=F(T).F
  FF._data[0] = dnlSquare(_data[0]) + dnlSquare(_data[1]) + dnlSquare(_data[2]);
  FF._data[1] = _data[0] * _data[1] + _data[1] * _data[3] + _data[2] * _data[4];
  FF._data[2] = _data[0] * _data[2] + _data[1] * _data[4] + _data[2] * _data[5];
  FF._data[3] = dnlSquare(_data[1]) + dnlSquare(_data[3]) + dnlSquare(_data[4]);
  FF._data[4] = _data[1] * _data[2] + _data[3] * _data[4] + _data[4] * _data[5];
  FF._data[5] = dnlSquare(_data[2]) + dnlSquare(_data[4]) + dnlSquare(_data[5]);

  // initialisation de Rot
  Rot11 = 1;
  Rot22 = 1;
  Rot33 = 1;
  Rot12 = 0;
  Rot13 = 0;
  Rot21 = 0;
  Rot23 = 0;
  Rot31 = 0;
  Rot32 = 0;

  while (it < 30)
  {
    // recherche du terme maxi
    q = (dnlAbs(FF._data[1]) > dnlAbs(FF._data[2]) ? 1 : 2);
    p = (dnlAbs(FF._data[q]) > dnlAbs(FF._data[4]) ? 0 : 1);
    if (p == 1)
      q = 2;

    // test de convergence
    if (dnlAbs(FF._data[_internalIndexes.index[p][q]]) < precisionPolarDecompose)
    {
      _tmp1 = sqrt(FF._data[0]);
      _tmp2 = sqrt(FF._data[3]);
      _tmp3 = sqrt(FF._data[5]);
      U._data[0] = (dnlSquare(Rot11) * _tmp1 + dnlSquare(Rot12) * _tmp2 + dnlSquare(Rot13) * _tmp3);
      U._data[3] = (dnlSquare(Rot21) * _tmp1 + dnlSquare(Rot22) * _tmp2 + dnlSquare(Rot23) * _tmp3);
      U._data[5] = (dnlSquare(Rot31) * _tmp1 + dnlSquare(Rot32) * _tmp2 + dnlSquare(Rot33) * _tmp3);
      U._data[1] = (Rot11 * Rot21 * _tmp1 + Rot12 * Rot22 * _tmp2 + Rot13 * Rot23 * _tmp3);
      U._data[2] = (Rot11 * Rot31 * _tmp1 + Rot12 * Rot32 * _tmp2 + Rot13 * Rot33 * _tmp3);
      U._data[4] = (Rot21 * Rot31 * _tmp1 + Rot22 * Rot32 * _tmp2 + Rot23 * Rot33 * _tmp3);

      _tmp1 = 1 / sqrt(FF._data[0]);
      _tmp2 = 1 / sqrt(FF._data[3]);
      _tmp3 = 1 / sqrt(FF._data[5]);
      _tmp4 = (Rot13 * _data[0] + Rot23 * _data[1] + Rot33 * _data[2]) * _tmp3;
      _tmp5 = (Rot12 * _data[0] + Rot22 * _data[1] + Rot32 * _data[2]) * _tmp2;
      _tmp6 = (Rot11 * _data[0] + Rot21 * _data[1] + Rot31 * _data[2]) * _tmp1;
      R._data[0] = Rot13 * _tmp4 + Rot12 * _tmp5 + Rot11 * _tmp6;
      R._data[1] = Rot23 * _tmp4 + Rot22 * _tmp5 + Rot21 * _tmp6;
      R._data[2] = Rot33 * _tmp4 + Rot32 * _tmp5 + Rot31 * _tmp6;
      _tmp4 = (Rot13 * _data[1] + Rot23 * _data[3] + Rot33 * _data[4]) * _tmp3;
      _tmp5 = (Rot12 * _data[1] + Rot22 * _data[3] + Rot32 * _data[4]) * _tmp2;
      _tmp6 = (Rot11 * _data[1] + Rot21 * _data[3] + Rot31 * _data[4]) * _tmp1;
      R._data[3] = Rot13 * _tmp4 + Rot12 * _tmp5 + Rot11 * _tmp6;
      R._data[4] = Rot23 * _tmp4 + Rot22 * _tmp5 + Rot21 * _tmp6;
      R._data[5] = Rot33 * _tmp4 + Rot32 * _tmp5 + Rot31 * _tmp6;
      _tmp4 = (Rot13 * _data[2] + Rot23 * _data[4] + Rot33 * _data[5]) * _tmp3;
      _tmp5 = (Rot12 * _data[2] + Rot22 * _data[4] + Rot32 * _data[5]) * _tmp2;
      _tmp6 = (Rot11 * _data[2] + Rot21 * _data[4] + Rot31 * _data[5]) * _tmp1;
      R._data[6] = Rot13 * _tmp4 + Rot12 * _tmp5 + Rot11 * _tmp6;
      R._data[7] = Rot23 * _tmp4 + Rot22 * _tmp5 + Rot21 * _tmp6;
      R._data[8] = Rot33 * _tmp4 + Rot32 * _tmp5 + Rot31 * _tmp6;
      return;
    }

    // calcul de l'angle
    alpha = (FF._data[_internalIndexes.index[q][q]] - FF._data[_internalIndexes.index[p][p]]) / (2. * FF._data[_internalIndexes.index[p][q]]);
    t = (alpha != 0 ? 1 / (alpha + sqrt(alpha * alpha + 1) * (alpha > 0 ? 1 : -1)) : 1);

    // calcul des cosinus
    Ctheta = 1 / (sqrt(1 + t * t));
    Stheta = t * Ctheta;

    // calcul du nouveau tenseur
    R.setToUnity();
    R(p, q) = Stheta;
    R(q, p) = -Stheta;
    R(p, p) = R(q, q) = Ctheta;

    // produit R(T).FF.R avec FF symétrique
    _tmp1 = FF._data[0] * R._data[0] + FF._data[1] * R._data[3] + FF._data[2] * R._data[6];
    _tmp2 = FF._data[1] * R._data[0] + FF._data[3] * R._data[3] + FF._data[4] * R._data[6];
    _tmp3 = FF._data[2] * R._data[0] + FF._data[4] * R._data[3] + FF._data[5] * R._data[6];
    _tmp4 = FF._data[0] * R._data[1] + FF._data[1] * R._data[4] + FF._data[2] * R._data[7];
    _tmp5 = FF._data[1] * R._data[1] + FF._data[3] * R._data[4] + FF._data[4] * R._data[7];
    _tmp6 = FF._data[2] * R._data[1] + FF._data[4] * R._data[4] + FF._data[5] * R._data[7];
    FF._data[5] = R._data[2] * (FF._data[0] * R._data[2] + FF._data[1] * R._data[5] + FF._data[2] * R._data[8]) +
                  R._data[5] * (FF._data[1] * R._data[2] + FF._data[3] * R._data[5] + FF._data[4] * R._data[8]) +
                  R._data[8] * (FF._data[2] * R._data[2] + FF._data[4] * R._data[5] + FF._data[5] * R._data[8]);
    FF._data[0] = R._data[0] * _tmp1 + R._data[3] * _tmp2 + R._data[6] * _tmp3;
    FF._data[1] = R._data[1] * _tmp1 + R._data[4] * _tmp2 + R._data[7] * _tmp3;
    FF._data[2] = R._data[2] * _tmp1 + R._data[5] * _tmp2 + R._data[8] * _tmp3;
    FF._data[3] = R._data[1] * _tmp4 + R._data[4] * _tmp5 + R._data[7] * _tmp6;
    FF._data[4] = R._data[2] * _tmp4 + R._data[5] * _tmp5 + R._data[8] * _tmp6;

    // cumulate of R
    _tmp1 = R._data[0] * Rot11 + R._data[3] * Rot12 + R._data[6] * Rot13;
    _tmp2 = R._data[1] * Rot11 + R._data[4] * Rot12 + R._data[7] * Rot13;
    Rot13 = R._data[2] * Rot11 + R._data[5] * Rot12 + R._data[8] * Rot13;
    _tmp4 = R._data[0] * Rot21 + R._data[3] * Rot22 + R._data[6] * Rot23;
    _tmp5 = R._data[1] * Rot21 + R._data[4] * Rot22 + R._data[7] * Rot23;
    Rot23 = R._data[2] * Rot21 + R._data[5] * Rot22 + R._data[8] * Rot23;
    _tmp7 = R._data[0] * Rot31 + R._data[3] * Rot32 + R._data[6] * Rot33;
    _tmp8 = R._data[1] * Rot31 + R._data[4] * Rot32 + R._data[7] * Rot33;
    Rot33 = R._data[2] * Rot31 + R._data[5] * Rot32 + R._data[8] * Rot33;
    Rot11 = _tmp1;
    Rot12 = _tmp2;
    Rot21 = _tmp4;
    Rot22 = _tmp5;
    Rot31 = _tmp7;
    Rot32 = _tmp8;

    it++;
  }

  std::cout << (*this) << std::endl;
  fatalError("SymTensor2::polar", "No convergence");
}

//  Saves the content of a SymTensor2 into a NumPy file
//-----------------------------------------------------------------------------
void SymTensor2::numpyWrite(std::string filename, bool initialize) const
//-----------------------------------------------------------------------------
{
  double nst[9];
  nst[0] = _data[0];
  nst[1] = _data[1];
  nst[2] = _data[2];
  nst[3] = _data[1];
  nst[4] = _data[3];
  nst[5] = _data[4];
  nst[6] = _data[2];
  nst[7] = _data[4];
  nst[8] = _data[5];
  std::string mode = "a";
  if (initialize)
    mode = "w";
  NumpyInterface::npySave(filename, &nst[0], {3, 3}, mode);
}

//  Saves the content of a SymTensor2 into a NumPyZ file
//-----------------------------------------------------------------------------
void SymTensor2::numpyWriteZ(std::string filename, std::string name, bool initialize) const
//-----------------------------------------------------------------------------
{
  double nst[9];
  nst[0] = _data[0];
  nst[1] = _data[1];
  nst[2] = _data[2];
  nst[3] = _data[1];
  nst[4] = _data[3];
  nst[5] = _data[4];
  nst[6] = _data[2];
  nst[7] = _data[4];
  nst[8] = _data[5];
  std::string mode = "a";
  if (initialize)
    mode = "w";
  NumpyInterface::npzSave(filename, name, &nst[0], {3, 3}, mode);
}

//  Read the content of a SymTensor2 from a NumPy file
//-----------------------------------------------------------------------------
void SymTensor2::numpyRead(std::string filename)
//-----------------------------------------------------------------------------
{
  NumpyInterface::NumpyArray arr = NumpyInterface::npyLoad(filename);
  if (arr.num_vals != 9)
  {
    std::cout << "ERROR\n";
  }
  _data[0] = arr.data<double>()[0];
  _data[1] = arr.data<double>()[1];
  _data[2] = arr.data<double>()[2];
  _data[3] = arr.data<double>()[4];
  _data[4] = arr.data<double>()[5];
  _data[5] = arr.data<double>()[8];

  // memcpy(_data, arr.data<double *>(), arr.num_vals * arr.word_size);
}

//  Read the content of a SymTensor2 from a NumPyZ file
//-----------------------------------------------------------------------------
void SymTensor2::numpyReadZ(std::string filename, std::string name)
//-----------------------------------------------------------------------------
{
  NumpyInterface::NumpyArray arr = NumpyInterface::npzLoad(filename, name);
  if (arr.num_vals != 9)
  {
    std::cout << "ERROR\n";
  }
  _data[0] = arr.data<double>()[0];
  _data[1] = arr.data<double>()[1];
  _data[2] = arr.data<double>()[2];
  _data[3] = arr.data<double>()[4];
  _data[4] = arr.data<double>()[5];
  _data[5] = arr.data<double>()[8];
}
