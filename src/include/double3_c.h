#ifndef _DOUBLE3_C_H
#define _DOUBLE3_C_H

#include <math.h>
////NEXT WILL BE 
struct vector_t_{
  vector_t_(int m_dim){x = new double[m_dim];}
  double *x;
  ~vector_t_(){delete[] x;}
};


struct double2 {
  double x, y;
  double2(){x=y=0.0;}
  //double3(double x_, double y_, double z_){x=x_;y=y_;z=z_;}
};

struct double3 {
  double x, y, z;
  double3(){x=y=z=0.0;}
  //double3(double x_, double y_, double z_){x=x_;y=y_;z=z_;}
};

inline double3 make_double3(double x_,double y_,double z_){double3 ret; ret.x=x_;
                                                                        ret.y=y_;
                                                                        ret.z=z_;
                                                          return ret;}

inline double3 operator+(double3 v1, double3 v2){return make_double3(v1.x+v2.x,v1.y+v2.y,v1.z+v2.z);}
inline double3 operator-(double3 v1, double3 v2){return make_double3(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z);}
//inline double3 operator*(double  s, double3 v  ){return make_double3(v.x*s,v.y*s,v.z*s);}


inline double2 make_double2(double x_,double y_){double2 ret; ret.x=x_;
                                                                        ret.y=y_;
                                                          return ret;}

inline double2 operator+(double2 v1, double2 v2){return make_double2(v1.x+v2.x,v1.y+v2.y);}
inline double2 operator-(double2 v1, double2 v2){return make_double2(v1.x-v2.x,v1.y-v2.y);}
inline double2 operator*(double  s, double2 v  ){return make_double2(v.x*s,v.y*s);}


static double3 operator*(const double3 &a, const double3 &b)
{
	return make_double3(a.x * b.x, a.y * b.y, a.z * b.z);
}




static double3 operator*(const double3 &a, const double &s)
{
	return make_double3(a.x * s, a.y * s, a.z * s);
}

inline double3 operator*(const double &s, const double3 &a)
{
	return make_double3(a.x * s, a.y * s, a.z * s);
}
static void operator*=(double3 &a, const double &s)
{
	a.x *= s; a.y *= s; a.z *= s;
}

inline double3 operator/(const double3 &a, const double &s)
{
	double inv = 1.0 / s;
	return a * inv;
}

inline double3 operator/(const double &s, const double3 &a)
{
	double inv = 1.0 / s;
	return a * inv;
}

inline double norm2(double3 &v){
  return v.x*v.x + v.y*v.y + v.z*v.z ;
}

inline double sqlength(double3 &v){
  return norm2(v);
}


inline double norm(double3 &v){
  return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

inline double3 cross(const double3 &a, const double3 &b)
{
	return make_double3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

// dot product
inline double dot(const double3 &a, const double3 &b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}
// squared length
inline double sqlength(const double3 &v)
{
	return dot(v, v);
}

// length
inline double length(const double3 &v)
{
	return sqrt(sqlength(v));
}

#endif
