/**********************************************************************************
 *                                                                                *
 *  DynELA Finite Element Code v.4.0                                              *
 *  by Olivier PANTALE                                                            *
 *  Olivier.Pantale@enit.fr                                                       *
 *                                                                                *
 *********************************************************************************/
//@!CODEFILE = DynELA-H-file
//@!BEGIN = PRIVATE

/*
  \file Macros.h
  Definition of the basic mathematical operators and constants.

  This file defines all the mathematical constants, tools and basic operators.
  \ingroup dnlKernel
*/

#ifndef __dnlKernel_Macros_h__
#define __dnlKernel_Macros_h__

#include <cmath>

/*
  Structure to decode indexes
*/
struct Tensor2Index
{
  short index[3][3]; // Set of indexes
};

/*
 * Predefined constants
 */
// If PI is not defined on this system the we define it
#ifndef PI
#define PI 3.14159265358979324 // Numerical value of the constant \f$ \pi \f$ if not already defined on the system
#endif

// Usefull constants to be defined
#define dnl2PI 6.28318530717958648        // Numerical value of the constant \f$ 2\pi \f$
#define dnlDegToRad 0.0174532925199432958 // Numerical value of the constant used to convert degree to radians \f$ \frac{\pi}{180} \f$
#define dnlInvSqrt2 0.707106781186547524  // Numerical value of the constant \f$ \frac{1}{\sqrt{2}} \f$
#define dnlInvSqrt3 0.577350269189625764  // Numerical value of the constant \f$ \frac{1}{\sqrt{3}} \f$
#define dnlRadToDeg 57.2957795130823209   // Numerical value of the constant used to convert radians to degree \f$ \frac{180}{\pi} \f$
#define dnlSqrt2 1.41421356237309505      // Numerical value of the constant \f$ \sqrt{2} \f$
#define dnlSqrt23 0.816496580927726034    // Numerical value of the constant \f$ \frac{2}{3} \f$
#define dnlSqrt3 1.73205080756887729      // Numerical value of the constant \f$ \sqrt{3} \f$
#define dnlSqrt32 1.22474487139158894     // Numerical value of the constant \f$ \frac{3}{2} \f$
#define dnlSqrt3m1 0.732050807568877294   // Numerical value of the constant \f$ \sqrt{3} - 1\f$
#define dnlSqrt3m2 -0.267949192431122706  // Numerical value of the constant \f$ \sqrt{3} - 2\f$
#define dnlSqrt3m3 -1.26794919243112271   // Numerical value of the constant \f$ \sqrt{3} - 3\f$
#define dnlSqrt3p1 2.73205080756887729    // Numerical value of the constant \f$ \sqrt{3} + 1\f$
#define dnlSqrt3p2 3.73205080756887729    // Numerical value of the constant \f$ \sqrt{3} + 2\f$
#define dnlSqrt3p3 4.73205080756887729    // Numerical value of the constant \f$ \sqrt{3} + 3\f$
#define dnlEpsilon 2.2204460492503131e-16 // Numerical value of the constant Epsilon

/*
 * Predefined mathematical functions
 */
#define dnlMax(i, j) ((i) > (j) ? (i) : (j))                                                                                                       // Returns the maximum of 2 values
#define dnlMin(i, j) ((i) < (j) ? (i) : (j))                                                                                                       // Returns the minimum of 2 values
#define dnlMax3(i, j, k) ((i) > (j) ? ((i) > (k) ? (i) : (k)) : ((j) > (k) ? (j) : (k)))                                                           // Returns the maximum of 3 values.
#define dnlMin3(i, j, k) ((i) < (j) ? ((i) < (k) ? (i) : (k)) : ((j) < (k) ? (j) : (k)))                                                           // Returns the minimum of 3 values.
#define dnlAbs(i) ((i) < 0 ? -(i) : (i))                                                                                                           // Returns the absolute value of the given argument.
#define dnlSign(i) ((i) >= 0 ? 1 : -1)                                                                                                             // Returns the sign of a given value (+1) if i>0, (-1) if i<0.
#define dnlMacAuley(x) ((x) <= 0 ? 0 : (x))                                                                                                        // Returns the MacAuley operator defined by: \f[ \left\langle x\right\rangle =\frac{1}{2}\left[x+\left|x\right|\right]\f]
#define dnlNormVec3D(x1, x2, x3) (sqrt(dnlSquare(x1) + dnlSquare(x2) + dnlSquare(x3)))                                                             // Returns the norm of a 3D vector defined by: \f[ \left\Vert x\right\Vert=\sqrt{x_1^2+x_2^2+x_3^2} \f]
#define dnlSurfTriangle(a, b, c) ((sqrt(((a) + (b) + (c)) * ((a) - (b) + (c)) * ((a) + (b) - (c)) * (-(a) + (b) + (c)))) / 4)                      // Returns the surface of a triangle function of the 3 edges lengths from the Heron formula: \f[ \frac{1}{4} \sqrt{-(a-b-c) (a+b-c) (a-b+c) (a+b+c)}\f]
#define dnlSurfTrianglePts(xa, ya, xb, yb, xc, yc) (dnlAbs((xb) * (ya) - (xc) * (ya) - (xa) * (yb) + (xc) * (yb) + (xa) * (yc) - (xb) * (yc)) / 2) // Returns the surface of a triangle function of the coords of the 3 corners from the Heron formula
#define dnlVonMises(s1, s2, s3, s4, s5, s6) \
  (sqrt((dnlSquare((s1) - (s2)) + dnlSquare((s1) - (s3)) + dnlSquare((s2) - (s3)) + 6 * (dnlSquare(s4) + dnlSquare(s5) + dnlSquare(s6))) / 2)) // Returns the J2 von-Mises stress with (s1,s2,s3) the normal stresses and (s4,s5,s6) the shear stresses : \f[ s = \frac {1}{\sqrt{2}}\sqrt{(s_{11}-s_{22})^2+(s_{22}-s_{33})^2+(s_{33}-s_{11})^2+6(s_{12}^2+s_{23}^2+s_{31}^2)}\f]
#define dnlLinearInterp(x1, x2, dx) ((dx) * (x2) + (1 - (dx)) * (x1))                                                                          // Returns the interpolated value between 2 points from the following equation: \f[ x=\lambda x_{1}+(1-\lambda)x_{2} \f]
#define dnlSquare(a) ((a) * (a))                                                                                                               // Returns the square value of the given argument from the following \f[ x^2 =x*x \f]
#define dnlBetweenValues(a, b, c) ((a) < (b) ? false : ((a) > (c) ? false : true))                                                             // Returns \ref true if the given parameter a is within the range [b,c], false if not
#define dnlBetweenValuesExclude(a, b, c) ((a) <= (b) ? false : ((a) >= (c) ? false : true))                                                    // Returns \ref true if the given parameter a is within the range ]b,c[, false if not
#define dnlTensor2Ind(i, j, m) ((i) * (m) + (j))                                                                                               // Returns the value of T(i,j) in a second order tensor of dimension m
#define dnlTensor3Ind(i, j, k, m) ((i) * (m * m) + (j) * (m) + k)                                                                              // Returns the value of T(i,j,k) in a third order tensor of dimension m
#define dnlTensor4Ind(i, j, k, l, m) ((i) * (m * m * m) + (j) * (m * m) + k * (m) + l)                                                         // Returns the value of T(i,j,k,l) in a fourth order tensor of dimension m
#define dnlKronecker(i, j) ((i) == (j) ? 1 : 0)                                                                                                // Returns the Kronecker: (1) if i=j, (0) in all other cases
/*
 * Standard mathematical defines functions
 */
#define precisionPolarDecompose 10.0E-10

#endif
