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
#include <Vec3D.h>
#include <NumpyInterface.h>
#include <Tensor2.h>
#include <SymTensor2.h>

/*
@LABEL:Vec3D::Vec3D()
@SHORT:Default constructor of the Vec3D class with initialization.
@RETURN:Vec3D
By default this constructor initializes all components of the vector to zero.
@END
*/
/*
@LABEL:Vec3D::Vec3D(double x, double y, double z)
@SHORT:Constructor of the Vec3D class with initialization.
@RETURN:Vec3D
@ARG:double & x & first component of the Vec3D to create.
@ARG:double & y & second component of the Vec3D to create.
@ARG:double & z & third component of the Vec3D to create.
@END
*/
//-----------------------------------------------------------------------------
Vec3D::Vec3D(double x, double y, double z)
//-----------------------------------------------------------------------------
{
    _data[0] = x;
    _data[1] = y;
    _data[2] = z;
}

// Copy constructor
//-----------------------------------------------------------------------------
Vec3D::Vec3D(const Vec3D &V)
//-----------------------------------------------------------------------------
{
    memcpy(_data, V._data, 3 * sizeof(double));
}

// Destructor
//-----------------------------------------------------------------------------
Vec3D::~Vec3D()
//-----------------------------------------------------------------------------
{
}

/*
@LABEL:Vec3D::operator*(double val, Vec3D x)
@SHORT:Multiplication of a vector by a scalar value.
@ARG:double & val & Scalar value to use for the multiplication operation.
@ARG:Vec3D & x & Vec3D to use for the multiplication operation.
This method defines the multiplication of a vector by a scalar value.
The result of this operation is also a vector defined by:
\begin{equation*}
\overrightarrow{y} = \lambda \overrightarrow{x}
\end{equation*}
where the $\overrightarrow{x}$ is the object itself and $\lambda$ is scalar value defined by parameter val.
@END
*/
//-----------------------------------------------------------------------------
Vec3D operator*(const double val, const Vec3D &V)
//-----------------------------------------------------------------------------
{
    Vec3D resu;

    resu._data[0] = val * V._data[0];
    resu._data[1] = val * V._data[1];
    resu._data[2] = val * V._data[2];

    return resu;
}

/*
@LABEL:Vec3D::operator*(double val)
@SHORT:Multiplication of a vector by a scalar value.
@ARG:double&val&Scalar value to use for the multiplication operation.
This method defines the multiplication of a vector by a scalar value.
The result of this operation is also a vector defined by:
\begin{equation*}
\overrightarrow{y} = \lambda \overrightarrow{x}
\end{equation*}
@END
*/
//-----------------------------------------------------------------------------
Vec3D Vec3D::operator*(const double val) const
//-----------------------------------------------------------------------------
{
    return Vec3D(val * _data[0], val * _data[1], val * _data[2]);
}

/*
@LABEL:Vec3D::operator/(double val)
@SHORT:Division of a vector by a scalar value.
@ARG:double&val&Scalar value to use for the division operation.
This method defines the division of a vector by a scalar value.
The result of this operation is also a vector defined by:
\begin{equation*}
\overrightarrow{y} = \frac{1}{\lambda}\overrightarrow{x}
\end{equation*}
@END
*/
//-----------------------------------------------------------------------------
Vec3D Vec3D::operator/(const double val) const
//-----------------------------------------------------------------------------
{
#ifdef VERIF_maths
    if (val == 0)
    {
        fatalError("Vec3D:: operator /", "divide by zero");
    }
#endif

    return Vec3D(_data[0] / val, _data[1] / val, _data[2] / val);
}

// Test the equality of two vectors
//-----------------------------------------------------------------------------
bool Vec3D::operator==(const Vec3D &vect) const
//-----------------------------------------------------------------------------
{
    if (_data[0] != vect._data[0])
        return false;
    if (_data[1] != vect._data[1])
        return false;
    if (_data[2] != vect._data[2])
        return false;

    return true;
}

// Test the inequality of two vectors
//-----------------------------------------------------------------------------
bool Vec3D::operator!=(const Vec3D &vect) const
//-----------------------------------------------------------------------------
{
    return !(*this == vect);
}

/*
@LABEL:Vec3D::maxVal()
@SHORT:Maximum value in a Vec3D.
@RETURN:double
This method returns the maximum value in a Vec3D.
@END
*/
//-----------------------------------------------------------------------------
double Vec3D::maxVal()
//-----------------------------------------------------------------------------
{
    if (_data[0] >= _data[1])
    {
        if (_data[0] >= _data[2])
            return _data[0];
        else
            return _data[2];
    }
    else
    {
        if (_data[1] >= _data[2])
            return _data[1];
        else
            return _data[2];
    }
}

/*
@LABEL:Vec3D::minVal()
@SHORT:Minumum value in a Vec3D.
@RETURN:double
This method returns the minimum value in a Vec3D.
@END
*/
//-----------------------------------------------------------------------------
double Vec3D::minVal()
//-----------------------------------------------------------------------------
{
    if (_data[0] <= _data[1])
    {
        if (_data[0] <= _data[2])
            return _data[0];
        else
            return _data[2];
    }
    else
    {
        if (_data[1] <= _data[2])
            return _data[1];
        else
            return _data[2];
    }
}

/*
@LABEL:Vec3D::maxAbs()
@SHORT:Maximum absolute value in a Vec3D.
@RETURN:double
This method returns the maximum absolute value in a Vec3D.
@END
*/
//-----------------------------------------------------------------------------
double Vec3D::maxAbs()
//-----------------------------------------------------------------------------
{
    if (dnlAbs(_data[0]) >= dnlAbs(_data[1]))
    {
        if (dnlAbs(_data[0]) >= dnlAbs(_data[2]))
            return dnlAbs(_data[0]);
        else
            return dnlAbs(_data[2]);
    }
    else
    {
        if (dnlAbs(_data[1]) >= dnlAbs(_data[2]))
            return dnlAbs(_data[1]);
        else
            return dnlAbs(_data[2]);
    }
}

/*
@LABEL:Vec3D::minAbs()
@SHORT:Minumum absolute value in a Vec3D.
@RETURN:double
This method returns the minimum absolute value in a Vec3D.
@END
*/
//-----------------------------------------------------------------------------
double Vec3D::minAbs()
//-----------------------------------------------------------------------------
{
    if (dnlAbs(_data[0]) <= dnlAbs(_data[1]))
    {
        if (dnlAbs(_data[0]) <= dnlAbs(_data[2]))
            return dnlAbs(_data[0]);
        else
            return dnlAbs(_data[2]);
    }
    else
    {
        if (dnlAbs(_data[1]) <= dnlAbs(_data[2]))
            return dnlAbs(_data[1]);
        else
            return dnlAbs(_data[2]);
    }
}

// Test if a given point is inside a given box defined by two opposite points
//-----------------------------------------------------------------------------
bool Vec3D::isInsideBox(const Vec3D &vect_m, const Vec3D &vect_M) const
//-----------------------------------------------------------------------------
{
    if (_data[0] < vect_m._data[0] || vect_M._data[0] < _data[0])
        return false;
    if (_data[1] < vect_m._data[1] || vect_M._data[1] < _data[1])
        return false;
    if (_data[2] < vect_m._data[2] || vect_M._data[2] < _data[2])
        return false;
    return true;
}

/*
@LABEL:Vec3D::dot(Vec3D y)
@SHORT:Dot product of two Vec3D.
@RETURN:double & The dot product.
@ARG:Vec3D & y & Vec3D $\overrightarrow{y}$ to use for the dot product operation.
This method returns the dot product of two Vec3D defined by the following equation:
\begin{equation*}
m = \overrightarrow{x}\cdot\overrightarrow{y},
\end{equation*}
where the $\overrightarrow{x}$ is the object itself and $\overrightarrow{y}$ is given by the parameter y.
@END
*/
//-----------------------------------------------------------------------------
double Vec3D::dot(const Vec3D &V) const
//-----------------------------------------------------------------------------
{
    return (_data[0] * V._data[0] + _data[1] * V._data[1] + _data[2] * V._data[2]);
}

/*
@LABEL:Vec3D::vectorProduct(Vec3D y)
@SHORT:Vector product of two Vec3D.
@RETURN:Vec3D & Result of the operation.
@ARG:Vec3D & y & Vec3D $\overrightarrow{y}$ to use for the vector product operation.
This method returns the vector product of two Vec3D defined by the following equation:
\begin{equation*}
\overrightarrow{w} = \overrightarrow{x}\land\overrightarrow{y},
\end{equation*}
where the $\overrightarrow{x}$ is the object itself and $\overrightarrow{y}$ is given by the parameter y.
@END
*/
//-----------------------------------------------------------------------------
Vec3D Vec3D::vectorProduct(const Vec3D &V) const
//-----------------------------------------------------------------------------
{
    return Vec3D(_data[1] * V._data[2] - _data[2] * V._data[1],
                 _data[2] * V._data[0] - _data[0] * V._data[2],
                 _data[0] * V._data[1] - _data[1] * V._data[0]);
}

/*
@LABEL:Vec3D::ewProduct(Vec3D y)
@SHORT:Element-wise product of two Vec3D.
@RETURN:Vec3D
@ARG:Vec3D & y & Vec3D $\overrightarrow{y}$ to use for the element-wise product operation.
This method returns the vector product of two Vec3D defined by the following equation:
\begin{equation*}
w_i = x_i y_i,
\end{equation*}
where the $\overrightarrow{x}$ is the object itself and $\overrightarrow{y}$ is given by the parameter y.
@END
*/
//-----------------------------------------------------------------------------
Vec3D Vec3D::ewProduct(const Vec3D &V) const
//-----------------------------------------------------------------------------
{
    return Vec3D(_data[0] * V._data[0], _data[1] * V._data[1], _data[2] * V._data[2]);
}

/*
@LABEL:Vec3D::ewProduct(Tensor2 B)
@SHORT:Element-wise product of two Vec3D.
@RETURN:Tensor2
@ARG:Tensor2 & B & Tensor2 $\B$ to use for the element-wise product operation.
This method returns the vector product of a Vec3D and a Tensor2 defined by the following equation:
\begin{equation*}
A_{ij} = x_i T_{ij},
\end{equation*}
where the $\overrightarrow{x}$ is the object itself and $\B$ is given by the parameter B.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 Vec3D::ewProduct(const Tensor2 &T) const
//-----------------------------------------------------------------------------
{
    Tensor2 result;

    result._data[0] = _data[0] * T._data[0];
    result._data[1] = _data[0] * T._data[1];
    result._data[2] = _data[0] * T._data[2];
    result._data[3] = _data[1] * T._data[3];
    result._data[4] = _data[1] * T._data[4];
    result._data[5] = _data[1] * T._data[5];
    result._data[6] = _data[2] * T._data[6];
    result._data[7] = _data[2] * T._data[7];
    result._data[8] = _data[2] * T._data[8];

    return result;
}

/*
@LABEL:Vec3D::dyadic(Vec3D y)
@SHORT:Dyadic product of two Vec3D.
@RETURN:Tensor2
@ARG:Vec3D&y&Vec3D $\overrightarrow{y}$ to use for the dyadic product operation.
This method returns the dyadic product of two Vec3D defined by the following equation:
\begin{equation*}
\T = \overrightarrow{x}\otimes\overrightarrow{y},
\end{equation*}
where the $\overrightarrow{x}$ is the object itself.
@END
*/
//-----------------------------------------------------------------------------
Tensor2 Vec3D::dyadic(const Vec3D &V) const
//-----------------------------------------------------------------------------
{
    Tensor2 result;

    result._data[0] = _data[0] * V._data[0];
    result._data[1] = _data[0] * V._data[1];
    result._data[2] = _data[0] * V._data[2];
    result._data[3] = _data[1] * V._data[0];
    result._data[4] = _data[1] * V._data[1];
    result._data[5] = _data[1] * V._data[2];
    result._data[6] = _data[2] * V._data[0];
    result._data[7] = _data[2] * V._data[1];
    result._data[8] = _data[2] * V._data[2];

    return result;
}

/*
@LABEL:Vec3D::dyadic()
@SHORT:Dyadic product of a Vec3D by itself.
@RETURN:SymTensor2
This method returns the dyadic product of two Vec3D defined by the following equation:
\begin{equation*}
\T = \overrightarrow{x}\otimes\overrightarrow{x},
\end{equation*}
where the $\overrightarrow{x}$ is the object itself. The result of this operation is a symmetric second order tensor.
@END
*/
//-----------------------------------------------------------------------------
SymTensor2 Vec3D::dyadic() const
//-----------------------------------------------------------------------------
{
    SymTensor2 result;

    result._data[0] = _data[0] * _data[0];
    result._data[1] = _data[0] * _data[1];
    result._data[2] = _data[0] * _data[2];
    result._data[3] = _data[1] * _data[1];
    result._data[4] = _data[1] * _data[2];
    result._data[5] = _data[2] * _data[2];

    return result;
}

// Writes a vector in a binary flux for storage
//-----------------------------------------------------------------------------
void Vec3D::write(std::ofstream &ofs) const
//-----------------------------------------------------------------------------
{
    ofs.write((char *)_data, 3 * sizeof(double));
}

// Reads a vector from a binary flux for storage
//-----------------------------------------------------------------------------
Vec3D &Vec3D::read(std::ifstream &ifs)
//-----------------------------------------------------------------------------
{
    ifs.read((char *)_data, 3 * sizeof(double));
    return *this;
}

//  Writes a vector to a binary flux for storage
//-----------------------------------------------------------------------------
std::ofstream &operator<<(std::ofstream &os, const Vec3D &vect)
//-----------------------------------------------------------------------------
{
    vect.write(os);
    return os;
}

// Reads a vector from a binary flux for storage
//-----------------------------------------------------------------------------
std::ifstream &operator>>(std::ifstream &is, Vec3D &vect)
//-----------------------------------------------------------------------------
{
    vect.read(is);
    return is;
}

// Display the content of a Vec3D
//-----------------------------------------------------------------------------
std::ostream &operator<<(std::ostream &os, const Vec3D &vector)
//-----------------------------------------------------------------------------
{
    vector.print(os);
    return os;
}

// Display the content of a Vec3D
//-----------------------------------------------------------------------------
void Vec3D::print(std::ostream &os) const
//-----------------------------------------------------------------------------
{
    os << "Vec3D coords={";
    os << _data[0];
    for (int i = 1; i < 3; i++)
    {
        os << "," << _data[i];
    }
    os << "}";
}

// Saves the content of a Vec3D into a NumPy file
//-----------------------------------------------------------------------------
void Vec3D::numpyWrite(std::string filename, bool initialize) const
//-----------------------------------------------------------------------------
{
    std::string mode = "a";
    if (initialize)
        mode = "w";
    NumpyInterface::npySave(filename, &_data[0], {3}, mode);
}

// Saves the content of a Vec3D into a NumPyZ file
//-----------------------------------------------------------------------------
void Vec3D::numpyWriteZ(std::string filename, std::string name, bool initialize) const
//-----------------------------------------------------------------------------
{
    std::string mode = "a";
    if (initialize)
        mode = "w";
    NumpyInterface::npzSave(filename, name, &_data[0], {3}, mode);
}

// Read the content of a Vec3D from a NumPy file
//-----------------------------------------------------------------------------
void Vec3D::numpyRead(std::string filename)
//-----------------------------------------------------------------------------
{
    NumpyInterface::NumpyArray arr = NumpyInterface::npyLoad(filename);
    if (arr.num_vals != 3)
    {
        std::cout << "ERROR\n";
    }
    memcpy(_data, arr.data<double *>(), arr.num_vals * arr.word_size);
}

// Read the content of a Vec3D from a NumPyZ file
//-----------------------------------------------------------------------------
void Vec3D::numpyReadZ(std::string filename, std::string name)
//-----------------------------------------------------------------------------
{
    NumpyInterface::NumpyArray arr = NumpyInterface::npzLoad(filename, name);
    if (arr.num_vals != 3)
    {
        std::cout << "ERROR\n";
    }
    memcpy(_data, arr.data<double *>(), arr.num_vals * arr.word_size);
}
